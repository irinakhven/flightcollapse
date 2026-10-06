"""The mono-exonic track.

Runs entirely separately from the spliced pipeline.  The empty intron chain is
a vacuous suffix of every chain, so allowing mono-exonic reads anywhere near
the suffix logic lets a single-exon stub absorb unbounded full-length read mass
-- 9.2% of reads are genuinely unspliced but ``isoseq collapse`` put 57.9% of
read mass on single-exon models.  Here they only ever group with each other.

Three cases, following the biology of an oligo-dT protocol:

1. **Last exon / 3'UTR of a known gene.**  Trusted.  With polyA capture, a long
   3'UTR that the RT did not fully traverse leaves exactly this: UTR plus part
   of the last exon, with a correct 3' end.  Emitted as a model, flagged as a
   fragment, and named after its parent transcript.
2. **Internal to a gene body / a long internal exon.**  Only credible if there
   is an alternative polyA site there.  Requires positive polyA evidence
   (signal hexamer or a catalogued site) and absence of internal priming.
3. **Intergenic.**  Same requirement as (2) with a higher support floor.

Note the scale trap that cost a full analysis round: ``perc_A_downstream`` is a
*percent* (0-100).  Comparing it to 0.6 flags essentially everything as
internally primed.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import MonoexonParams
from .ends import cluster_3p_ends, five_prime_end, tail_evidence
from .genome import Genome, PolyASiteAtlas
from .intervals import Exon, overlap
from .model import TranscriptModel, polya_evidence
from .molecules import NO_CELL, NO_UMI, n_cells, n_molecules
from .reads import ContigReads
from .reference import ReferenceIndex


def _split_by_overlap(
    starts: np.ndarray, ends: np.ndarray, idx: np.ndarray, min_frac: float
) -> List[np.ndarray]:
    """Within a 3' peak, split reads that do not mutually overlap enough.

    Single-linkage on "overlap covers at least ``min_frac`` of the shorter
    read", walking left to right.
    """
    if idx.size <= 1:
        return [idx]
    order = idx[np.argsort(starts[idx], kind="stable")]
    out: List[List[int]] = [[int(order[0])]]
    cur_s, cur_e = int(starts[order[0]]), int(ends[order[0]])
    for j in order[1:]:
        s, e = int(starts[j]), int(ends[j])
        ov = overlap((cur_s, cur_e), (s, e))
        if ov >= min_frac * min(cur_e - cur_s, e - s):
            out[-1].append(int(j))
            cur_s, cur_e = min(cur_s, s), max(cur_e, e)
        else:
            out.append([int(j)])
            cur_s, cur_e = s, e
    return [np.array(g, np.int64) for g in out]


def collapse_monoexonic(
    reads: ContigReads,
    contig: str,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    params: MonoexonParams,
    read_ok: Optional[np.ndarray] = None,
    five_prime_percentile: float = 90.0,
    umi_hamming: int = 1,
    spliced_3p: Optional[Dict[str, np.ndarray]] = None,
    rejections: Optional[List[Dict[str, object]]] = None,
    end_params: Optional[object] = None,
    demoted: Optional[List[Tuple[str, int, np.ndarray]]] = None,
    spliced_spans: Optional[Dict[str, np.ndarray]] = None,
    terminal_params: Optional[object] = None,
    cov: Optional[object] = None,
    cov_params: Optional[object] = None,
) -> Tuple[List[TranscriptModel], Dict[str, int]]:
    stats = {
        "monoexon_reads": 0,
        "clusters": 0,
        "models_utr3": 0,
        "models_internal": 0,
        "models_intergenic": 0,
        "rejected_low_support": 0,
        "rejected_no_polya": 0,
        "rejected_internal_priming": 0,
        "rejected_demoted_fragment": 0,
        "rejected_sr_continuation": 0,
        "reads_rejected_sr_continuation": 0,
        "reads_demoted_fragment": 0,
        "rescued_by_tail": 0,
        "reads_rejected_low_support": 0,
        "reads_rejected_no_polya": 0,
        "reads_rejected_internal_priming": 0,
        "reads_in_models": 0,
    }
    models: List[TranscriptModel] = []
    if not params.enabled or reads.n == 0:
        return models, stats

    mono = np.diff(reads.joff) == 0
    if read_ok is not None:
        mono = mono & read_ok
    idx_all = np.flatnonzero(mono)
    stats["monoexon_reads"] = int(idx_all.size)
    if idx_all.size == 0:
        return models, stats

    counter = defaultdict(int)
    for schar, sint in (("+", 1), ("-", -1)):
        idx = idx_all[reads.strand[idx_all] == sint]
        if idx.size == 0:
            continue
        labels, peaks = cluster_3p_ends(reads.tts[idx], params.max_3p_diff)
        for li in range(peaks.size):
            sel = idx[labels == li]
            for grp in _split_by_overlap(reads.start, reads.end, sel, params.min_overlap_frac):
                stats["clusters"] += 1
                m = _make_monoexon_model(
                    reads, contig, schar, grp, ref, genome, atlas, params,
                    five_prime_percentile, umi_hamming, counter, stats,
                    spliced_3p, rejections, end_params, demoted,
                    spliced_spans, terminal_params, cov, cov_params,
                )
                if m is not None:
                    models.append(m)
                    stats["reads_in_models"] += m.n_reads
    return models, stats


def _make_monoexon_model(
    reads: ContigReads,
    contig: str,
    strand: str,
    grp: np.ndarray,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    params: MonoexonParams,
    five_prime_percentile: float,
    umi_hamming: int,
    counter: Dict[str, int],
    stats: Dict[str, int],
    spliced_3p: Optional[Dict[str, np.ndarray]] = None,
    rejections: Optional[List[Dict[str, object]]] = None,
    end_params: Optional[object] = None,
    demoted: Optional[List[Tuple[str, int, np.ndarray]]] = None,
    spliced_spans: Optional[Dict[str, np.ndarray]] = None,
    terminal_params: Optional[object] = None,
    cov: Optional[object] = None,
    cov_params: Optional[object] = None,
) -> Optional[TranscriptModel]:
    if grp.size == 0:
        return None
    pairs = [(int(reads.cell[j]), int(reads.umi[j])) for j in grp]
    has_umi = any(c != int(NO_CELL) and u != int(NO_UMI) for c, u in pairs)
    n_mols = n_molecules(pairs, umi_hamming) if has_umi else 0
    n_cell = n_cells(pairs) if has_umi else 0
    n_reads = int(grp.size)

    three = int(np.median(reads.tts[grp]))
    five = five_prime_end(reads.tss[grp], three, strand, five_prime_percentile)
    exon: Exon = (min(five, three), max(five, three))
    if exon[1] - exon[0] < 1:
        return None

    gene_id = ref.assign_gene(contig, [exon], strand)
    ev = polya_evidence(
        contig, three, strand, genome, atlas,
        motif_window=params.polya_motif_window,
        perc_a_window=params.perc_a_window,
        max_perc_a=params.max_perc_a_downstream,
        motif_min_dist=params.polya_motif_min_dist,
        motif_max_dist=params.polya_motif_max_dist,
        motif_strict=params.polya_motif_strict,
        require_atlas_when_available=params.require_atlas_when_available,
    )

    flags: List[str] = []
    parent_tx: Optional[str] = None
    if gene_id is None:
        category = "monoexon_intergenic"
        min_reads = params.intergenic_min_reads
        need_polya = params.intergenic_require_polya
        name_stub = f"NOVELMONO_{contig}"
    else:
        frac_term = ref.in_terminal_exon(gene_id, exon)
        frac_utr3 = ref.in_utr3(gene_id, exon)
        if max(frac_term, frac_utr3) >= 0.5:
            category = "monoexon_3UTR"
            min_reads = params.utr3_min_reads
            need_polya = False
            parent_tx = _nearest_parent_tx(ref, gene_id, three)
            name_stub = parent_tx or ref.genes[gene_id].gene_name
            d, _site = ref.nearest_tts(gene_id, three)
            if abs(d) <= params.max_3p_diff:
                flags.append("rt_dropoff_fragment_of_annotated_3p_end")
        else:
            category = "monoexon_internal"
            min_reads = params.min_reads
            need_polya = params.internal_require_polya
            name_stub = ref.genes[gene_id].gene_name or gene_id

    def _reject(reason: str, detail: str = "") -> None:
        """Record WHY, not just that. The SIRV audit could not answer this."""
        stats[f"rejected_{reason}"] = stats.get(f"rejected_{reason}", 0) + 1
        stats[f"reads_rejected_{reason}"] = (
            stats.get(f"reads_rejected_{reason}", 0) + n_reads
        )
        if rejections is not None:
            rejections.append({
                "contig": contig, "strand": strand,
                "start": exon[0], "end": exon[1], "length": exon[1] - exon[0],
                "category": category, "gene_id": gene_id or "",
                "n_reads": n_reads, "n_molecules": n_mols, "n_cells": n_cell,
                "perc_a_downstream": ev.get("perc_a_downstream", ""),
                "internal_priming": int(bool(ev.get("internal_priming"))),
                "polya_motif_found": ev.get("polya_motif_found", ""),
                "tail_molecule_frac": round(float(tail["tail_molecule_frac"]), 4),
                "median_tail_len": tail["median_tail_len"],
                "reason": reason, "detail": detail,
            })

    tail = tail_evidence(reads, grp, np.ones(grp.size, float), end_params) \
        if end_params is not None else \
        {"tail_molecule_frac": 0.0, "n_tail_reads": 0.0, "median_tail_len": 0.0}

    # -- 0.3.0 demotion: is this a fragment of a spliced model already called?
    # -- 0.4.0 Rule M: for a 3'UTR cluster, demotion is the DEFAULT ---------
    # The evidence that lets one survive, computed here because the demotion
    # decisions below need it and the model may never be built.
    _disp = float(np.std(reads.tts[grp])) if grp.size else 0.0
    _max_disp = 12.0 if terminal_params is None \
        else float(terminal_params.max_dispersion_promote)
    _d_site = ev.get("dist_to_polya_site")
    _at_site = _d_site is not None and abs(_d_site) <= params.max_3p_diff
    _at_ann = False
    if gene_id is not None and gene_id in ref.genes:
        _d_tts, _ = ref.nearest_tts(gene_id, three)
        _at_ann = abs(_d_tts) <= params.max_3p_diff
    # A PAS hexamer plus a peak tight enough to be a cleavage event is the
    # locus-level stand-in where no atlas exists. Either alone is not enough:
    # hexamers are six bases and occur everywhere, and a tight peak on its own
    # is what a reproducible RT drop-off also looks like.
    _evidenced = params.keep_fragment_if_evidenced and (
        _at_ann or _at_site or (bool(ev.get("polya_motif_found")) and _disp <= _max_disp)
    )

    if (category == "monoexon_3UTR" and not _evidenced
            and params.demote_contained_fragments and spliced_spans
            and gene_id in spliced_spans):
        # Containment, not coincidence. The 0.3.x rule looks for a fragment
        # that stops WHERE ITS PARENT STOPS; on SIRV not one of the 186
        # surviving mono-exon models does that, because RT drop-off stops
        # anywhere. Lying wholly inside an emitted spliced model is what a
        # fragment actually looks like.
        sp = spliced_spans[gene_id]
        inside = sp[(sp[:, 0] <= exon[0]) & (sp[:, 1] >= exon[1])]
        if inside.size:
            stats["reads_demoted_fragment"] += n_reads
            _reject("demoted_fragment",
                    "inside the span of an emitted spliced model of this gene")
            if demoted is not None:
                demoted.append((gene_id, int(inside[0, 2]), grp))
            return None

    if (category == "monoexon_3UTR" and not _evidenced
            and params.demote_on_sr_continuation and cov is not None
            and cov_params is not None and cov_params.enabled):
        # Short reads showing NO drop across this 3' end mean the transcript
        # demonstrably continues past it. Only a POSITIVE measurement counts:
        # `supported is None` is no coverage or coverage below the floor, and
        # silence must never be read as evidence of truncation.
        step = cov.three_prime_step(contig, three, strand, cov_params)
        if step.get("sr_3p_step_supported") is False:
            # Demote only when there is a model to fold the reads INTO. With no
            # containing spliced model the reads would have to be discarded,
            # and a discarded read needs its own accounting key or the
            # unspliced-read invariant silently stops adding up. Short reads
            # disagreeing with an end we cannot re-home is a thing to record,
            # not a licence to drop the evidence -- so it is flagged on the
            # model instead and the model is still emitted.
            # OVERLAP, not containment. A containing model would already have
            # been demoted by the clause above and returned, so searching for
            # one here could only ever come up empty -- which is exactly what
            # happened: `monoexon_models_demoted_by_short_reads` was 0 on
            # BD176c, not because short reads never disagreed but because this
            # clause was structurally unable to act on it. The models that
            # reach here are the ones NO spliced model contains, so the fold
            # target has to be one that overlaps.
            tgt = None
            if spliced_spans and gene_id in spliced_spans:
                sp = spliced_spans[gene_id]
                ov = np.minimum(sp[:, 1], exon[1]) - np.maximum(sp[:, 0], exon[0])
                if ov.size and ov.max() >= 0.5 * (exon[1] - exon[0]):
                    tgt = int(sp[int(np.argmax(ov)), 2])
            if tgt is not None:
                stats["reads_demoted_fragment"] += n_reads
                _reject("sr_continuation",
                        f"short-read coverage does not drop here "
                        f"(ratio {step.get('sr_3p_step_ratio')})")
                if demoted is not None:
                    demoted.append((gene_id, tgt, grp))
                return None
            flags.append("sr_coverage_continues_past_3p")

    # A cluster inside a terminal exon whose 3' peak coincides with an emitted
    # spliced model's 3' end is that transcript's last exon, not a transcript.
    if (params.demote_terminal_exon_fragments and category == "monoexon_3UTR"
            and spliced_3p is not None and gene_id in spliced_3p):
        ends = spliced_3p[gene_id]
        if ends.size and np.min(np.abs(ends - three)) <= params.max_3p_diff:
            # `_reject` owns the model and read counters for every reason,
            # including this one -- incrementing rejected_demoted_fragment here
            # as well double-counted it, and the QC reported exactly 2x the
            # true number of demoted fragments (the per-candidate rejection
            # table, written from the same call, had it right).
            #
            # reads_demoted_fragment is a SEPARATE key and is not what _reject
            # writes (that is reads_rejected_demoted_fragment). It feeds
            # monoexon_reads_folded_into_spliced, which is a different claim:
            # these reads were re-homed onto a spliced model, not discarded
            # like the reads of a low-support or internally-primed candidate.
            stats["reads_demoted_fragment"] += n_reads
            _reject("demoted_fragment",
                    f"3p within {params.max_3p_diff}bp of an emitted spliced model")
            # The reads are FOLDED IN, not discarded. Demoting a fragment is a
            # statement about which model owns this read mass, not a reason to
            # lose it -- and the read-accounting invariant would (correctly)
            # fail if it were dropped.
            if demoted is not None:
                demoted.append(
                    (gene_id, int(ends[int(np.argmin(np.abs(ends - three)))]), grp)
                )
            return None

    if n_reads < min_reads or (has_umi and n_mols < params.min_umis):
        _reject("low_support", f"{n_reads} reads / {n_mols} molecules")
        return None
    if ev["internal_priming"]:
        _reject("internal_priming",
                f"perc_a_downstream={ev.get('perc_a_downstream')}")
        return None
    if need_polya and not ev["polya_supported"]:
        # -- 0.3.0 rescue: the genomic gate is a claim about the locus; the
        # read's own terminal A-run is a claim about the molecule. The spliced
        # path has accepted that substitution since 0.1.18.
        rescued = (
            params.tail_rescues_polya
            and tail["tail_molecule_frac"] >= params.min_tail_molecule_frac
            and n_cell >= params.rescue_min_cells
        )
        if not rescued:
            _reject("no_polya",
                    f"tail_frac={tail['tail_molecule_frac']:.2f} "
                    f"cells={n_cell}")
            return None
        stats["rescued_by_tail"] += 1
        flags.append("polya_from_read_tail")

    # Naming is deliberately NOT done here. This counter was per *contig* and
    # keyed on a gene NAME, and gene names repeat: `Y_RNA` exists on chr9, chr10
    # and chr16 as three different gene_ids, so all three got `Y_RNA|APA1`.
    # That is the same collision class as DNAJC9-AS1, for the third time, and it
    # survived the 0.1.10 fix only because this function bypassed IdAssigner by
    # setting model_id itself. It no longer does: the assigner counts on the
    # exact string the name is built from, globally, and cannot repeat.
    counter[name_stub] += 1
    stats[
        {"monoexon_3UTR": "models_utr3",
         "monoexon_internal": "models_internal",
         "monoexon_intergenic": "models_intergenic"}[category]
    ] += 1

    m = TranscriptModel(
        model_id="",
        contig=contig,
        strand=strand,
        exons=[exon],
        chain=(),
        category=category,
        reads=grp,
        n_reads=n_reads,
        n_mols=n_mols,
        n_cells=n_cell,
        parent_tx=parent_tx,
        gene_id=gene_id,
        gene_name=ref.genes[gene_id].gene_name if gene_id in ref.genes else "",
        flags=flags,
    )
    m.evidence.update(ev)
    m.evidence["three_prime_dispersion"] = float(np.std(reads.tts[grp]))
    m.evidence.update(tail)

    # -- 0.3.1: the terminal-score fields. --------------------------------
    # These were set only on the spliced path, so every mono-exonic model
    # reached score_terminals() with end_is_annotated = end_at_polya_site = 0.
    # `pos = ann | atlas` was therefore FALSE for all of them by construction:
    # a mono-exonic model could be a decoy or an unlabelled candidate, never a
    # positive. On BD176c that produced a unanimous verdict -- all 7,566
    # monoexon_3UTR, all 12,382 monoexon_internal and all 1,383
    # monoexon_intergenic models sat above any ceiling, 100.0% of each. A
    # filter that removes a whole category is not judging it.
    #
    # The distances behind both flags were already measured above by
    # polya_evidence() and ref.nearest_tts(); only the flags were missing.
    d_site = ev.get("dist_to_polya_site")
    m.evidence["end_at_polya_site"] = int(
        d_site is not None and abs(d_site) <= params.max_3p_diff
    )
    m.evidence["dist_polya_site"] = "" if d_site is None else int(abs(d_site))
    if gene_id is not None and gene_id in ref.genes:
        d_tts, _ = ref.nearest_tts(gene_id, three)
        m.evidence["end_is_annotated"] = int(abs(d_tts) <= params.max_3p_diff)
        m.evidence["dist_ann_tts"] = int(abs(d_tts))
    else:
        # intergenic: there is no annotated end to be near, which is a real
        # statement about the model, not a missing measurement
        m.evidence["end_is_annotated"] = 0
        m.evidence["dist_ann_tts"] = ""
    # A mono-exonic model has no intron chain, so it IS the whole of its own
    # group: its chain share is 1 and it does not compete with a sibling peak.
    # 1.0 is the neutral value for both, not a filler.
    m.evidence["chain_share"] = 1.0
    m.evidence["ratio_to_dominant_peak"] = 1.0
    m.evidence["n_3p_peaks"] = 1

    # -- 0.4.0: tier the mono-exon track too -------------------------------
    # Rule M can only DEMOTE a fragment when there is a spliced model to fold
    # its reads into. That is not a detail of the implementation, it is the
    # binding constraint: a mono-exon pile in a gene for which no spliced model
    # was emitted has nowhere to send its reads, so removing it would discard
    # them and the unspliced-read invariant would stop adding up.
    #
    # On SIRV that is 77 of the 89 models the rule cannot touch, and all 89
    # have wrong 3' ends. They are emitted because their reads must go
    # somewhere, not because the end is believed -- so they are tiered, and a
    # downstream filter on `terminal_tier` recovers the precision that
    # demotion cannot safely take. Filtering `unsupported` out of the SIRV fair
    # arm gives 110 models at 75.5% 3'-end precision against 199 at 41.7%,
    # with the same 83 correct ends either way.
    if _at_ann:
        m.evidence["terminal_tier"] = "annotated"
    elif _at_site:
        m.evidence["terminal_tier"] = "high"
    elif bool(ev.get("polya_motif_found")) and _disp <= _max_disp:
        m.evidence["terminal_tier"] = "reported"
    else:
        m.evidence["terminal_tier"] = "unsupported"
        if not (spliced_spans and gene_id in spliced_spans):
            flags.append("no_spliced_model_to_fold_into")
            m.flags = sorted(set(m.flags) | {"no_spliced_model_to_fold_into"})
    m.evidence.setdefault("candidate_share", "")
    return m


def _nearest_parent_tx(ref: ReferenceIndex, gene_id: str, pos: int) -> Optional[str]:
    """Transcript of the gene whose annotated 3' end is nearest, canonical first."""
    g = ref.genes.get(gene_id)
    if g is None or not g.tx_ids:
        return None
    return min(
        g.tx_ids,
        key=lambda t: (abs(ref.tx[t].tts - pos), ref.tx[t].rank, -ref.tx[t].tx_len),
    )
