"""Terminal-end resolution within a chain group.

Two asymmetries drive the design, both specific to oligo-dT chemistry:

**The 3' end is the anchored end.**  Priming happens at the polyA tail, so the
3' ends of a group are tight (median offset ~4 bp in the isoseq groups) and a
genuinely different 3' end means a genuinely different cleavage site.  3' ends
are therefore clustered into peaks, and a peak that matches no annotated end
becomes its own model -- but only with positive polyA evidence, so that
internal priming does not become an isoform.

**The 5' end is the degraded end.**  57.8% of spliced read mass is 5'-truncated.
A model's 5' end is a high percentile of its members' 5' ends, never the
minimum -- taking the minimum is the single line that produced the isoseq
artifact.  A 5'-truncated pile is folded back onto its parent's annotated start
and recorded as ``n_5p_truncated_reads`` rather than emitted, because otherwise
every well-expressed gene grows a ladder of truncation models, which is the
diversity inflation this tool exists to avoid.

Naming, since 0.1.7
-------------------
The categories are stated as **terminal-end differences**, not as UTR
differences.  An identical intron chain guarantees only that the difference sits
in the terminal exons; whether it is *untranslated* depends on where the CDS
ends, and is reported separately as ``is_utr_only``.

===================  ====================================================
``FSM``              chain and both ends match one annotated transcript
``end3_annotated``   3' end is the annotated end of another isoform
``end3_novel``       3' end matches nothing annotated; polyA-supported
``end5_alt``         5' end sits on an annotated TSS of another isoform
``end5_extended``    5' end reaches beyond every annotated start
``IR``               chain = a reference chain minus introns read through
``READTHROUGH``      3' end has run past the gene into the next one
===================  ====================================================

Three of these are decided against **every transcript sharing the chain**, not
against the canonical one.  Comparing a read pile only to the MANE transcript
turns every other isoform's genuine 3' end into a discovery: measured on this
data, 35% of 3' ends fall within 10 bp of the associated transcript's end but
72% fall within 10 bp of the nearest annotated end of *any* transcript of the
gene.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import EndParams, MonoexonParams
from .genome import Genome, PolyASiteAtlas
from .intervals import Chain, Exon, build_exons, signed_5p_offset
from .model import TranscriptModel, polya_evidence
from .molecules import NO_CELL, NO_UMI, n_cells, n_molecules
from .reads import ContigReads
from .reference import ReferenceIndex
from .similarity import nearest_transcript, readthrough_gene, retained_intron_rank

#: emission order when two 3' peaks resolve onto the same final end
RANK = {
    "FSM": 0,
    "end3_annotated": 1,
    "end3_novel": 2,
    "READTHROUGH": 3,
    "NOVEL": 4,
    # lowest priority: if a refused-fallback peak lands on the same coordinate
    # as a peak that resolved properly, the resolved statement is the true one
    "end3_unresolved": 5,
}


def molecule_weights(
    reads: ContigReads, read_idx: np.ndarray, umi_hamming: int = 0
) -> np.ndarray:
    """1/(family size) per read, so one molecule contributes one unit of mass.

    New in 0.1.17.  Peak discovery, the dominant-peak choice and the peak
    coordinate were all driven by raw read counts, while the support *gates*
    counted molecules -- so PCR duplicates of a single molecule could decide
    which cleavage site a chain group was named after.

    Families are keyed on exact ``(cell, UMI)``.  With ``umi_hamming=0`` that is
    the molecule definition the rest of the tool uses and the weighting is
    exact; with a larger hamming radius it under-merges, which makes this an
    approximation -- but a strictly better one than counting reads.  Reads with
    no cell or no UMI keep weight 1: unknown provenance is not evidence of
    duplication.
    """
    if read_idx.size == 0:
        return np.zeros(0, float)
    cells = reads.cell[read_idx]
    umis = reads.umi[read_idx]
    known = (cells != np.uint32(NO_CELL)) & (umis != np.uint64(NO_UMI))
    w = np.ones(read_idx.size, float)
    if not known.any():
        return w
    keys = np.stack([cells.astype(np.int64), umis.astype(np.int64)], axis=1)[known]
    _, inv, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    w[known] = 1.0 / counts[inv]
    return w


def weighted_median(values: np.ndarray, weights: np.ndarray) -> int:
    if values.size == 0:
        return 0
    o = np.argsort(values, kind="stable")
    v, cw = values[o], np.cumsum(weights[o])
    if cw[-1] <= 0:
        return int(np.median(values))
    return int(v[min(int(np.searchsorted(cw, cw[-1] / 2.0)), v.size - 1)])


def cluster_3p_ends(
    positions: np.ndarray, window: int, weights: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Greedy peak clustering of 3' ends.

    Repeatedly take the position with the most support inside +/- ``window``,
    claim those reads, and continue.  This gives well-separated peaks; naive
    single-linkage would chain a whole 3'UTR into one cluster.

    Returns ``(label_per_position, peak_positions)``.
    """
    n = positions.size
    labels = np.full(n, -1, np.int64)
    if n == 0:
        return labels, np.zeros(0, np.int64)
    w = np.ones(n) if weights is None else weights.astype(float)
    order = np.argsort(positions, kind="stable")
    pos_s = positions[order]
    w_s = w[order]

    remaining = np.ones(n, bool)
    peaks: List[int] = []
    label = 0
    while remaining.any():
        idx_rem = np.flatnonzero(remaining)
        pr = pos_s[idx_rem]
        wr = w_s[idx_rem]
        cs = np.concatenate([[0.0], np.cumsum(wr)])
        lo = np.searchsorted(pr, pr - window, side="left")
        hi = np.searchsorted(pr, pr + window, side="right")
        mass = cs[hi] - cs[lo]
        k = int(np.argmax(mass))
        centre = int(pr[k])
        sel = idx_rem[(pr >= centre - window) & (pr <= centre + window)]
        # peak position = weighted median of the claimed reads
        claimed = pos_s[sel]
        cw = w_s[sel]
        srt = np.argsort(claimed)
        cc = np.cumsum(cw[srt])
        peak = int(claimed[srt][np.searchsorted(cc, cc[-1] / 2.0)])
        labels[order[sel]] = label
        peaks.append(peak)
        remaining[sel] = False
        label += 1
    return labels, np.array(peaks, np.int64)


def five_prime_end(
    tss: np.ndarray, tts_ref: int, strand: str, percentile: float
) -> int:
    """Model 5' end = the ``percentile``-th percentile of member extension.

    Extension is measured as exonic-agnostic genomic distance upstream of the
    group's 3' end, so the statistic is monotone in "how full-length is this
    read" regardless of strand.
    """
    if tss.size == 0:
        return tts_ref
    ext = np.abs(tss.astype(np.int64) - tts_ref)
    d = float(np.percentile(ext, percentile))
    return int(tts_ref + d) if strand == "-" else int(tts_ref - d)


def _nearest_by_tts(ref: ReferenceIndex, tids: Sequence[str], pos: int) -> Optional[str]:
    if not tids:
        return None
    return min(tids, key=lambda t: (abs(ref.tx[t].tts - pos), ref.tx[t].rank,
                                    -ref.tx[t].tx_len, t))


def _tts_owners(ref: ReferenceIndex, gene_id: Optional[str], site: int) -> str:
    """Which transcripts of the gene end exactly at ``site``."""
    if not gene_id or gene_id not in ref.genes:
        return ""
    hits = [t for t in ref.genes[gene_id].tx_ids if ref.tx[t].tts == site]
    return ",".join(hits[:4])


# ---------------------------------------------------------------------- #
def build_models_for_group(
    reads: ContigReads,
    contig: str,
    strand: str,
    chain: Chain,
    read_idx: np.ndarray,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    params: EndParams,
    mono_params: MonoexonParams,
    chain_tx: Sequence[str],
    gene_id: Optional[str],
    umi_hamming: int = 1,
) -> List[TranscriptModel]:
    """Split one chain group into 3'-end models and name them against the reference.

    ``chain_tx`` is **every** annotated transcript whose intron chain equals this
    group's chain, not the canonical one among them.  Those transcripts differ
    only in their ends, which is exactly the thing being decided here.

    Three stages, and the order matters.  A 3' peak that gets *rejected* is
    folded back onto the parent's annotated end, which means several peaks can
    land on the same final 3' end.  Those are the same transcript and must be
    merged into one model before the 5' end is computed -- emitting them
    separately produced several models sharing one ENST name, split that
    transcript's read mass across them, and gave each a different 5' end.
    """
    if read_idx.size == 0:
        return []
    # a bare transcript id is iterable, and iterating it yields characters --
    # accept the older single-parent form rather than failing obscurely
    if isinstance(chain_tx, str):
        chain_tx = [chain_tx]
    chain_tx = [t for t in (chain_tx or ()) if t in ref.tx]
    tts = reads.tts[read_idx]
    tss = reads.tss[read_idx]
    # Separable, but by default the same: measured on BD144, the narrow window
    # consolidates rather than fragments (see EndParams.peak_window).
    window = params.max_3p_diff if params.peak_window is None else params.peak_window
    # 0.1.17: peak GEOMETRY is molecule-weighted; support GATES still count
    # whole reads and whole molecules, so `min_end_reads` keeps its meaning.
    wts = (molecule_weights(reads, read_idx, umi_hamming)
           if getattr(params, "weight_ends_by_molecule", False)
           else np.ones(read_idx.size, float))
    labels, peaks = cluster_3p_ends(tts, window, weights=wts)
    sizes = np.bincount(labels, minlength=peaks.size)              # reads
    mass = np.bincount(labels, weights=wts, minlength=peaks.size)  # molecules
    order = np.argsort(-mass)

    # group default parent: the transcript matching the *dominant* peak's end
    dominant = int(peaks[int(order[0])]) if peaks.size else None
    default_tx = _nearest_by_tts(ref, chain_tx, dominant) if dominant is not None else None
    dtx = ref.tx.get(default_tx) if default_tx else None

    def support(sel: np.ndarray):
        pairs = [(int(reads.cell[j]), int(reads.umi[j])) for j in sel]
        has = any(c != int(NO_CELL) and u != int(NO_UMI) for c, u in pairs)
        return (n_molecules(pairs, umi_hamming) if has else 0,
                n_cells(pairs) if has else 0, has)

    # keep peaks with enough support; fold the rest into the nearest kept peak.
    # A peak sitting on the annotated end of ANY transcript with this chain is
    # kept whatever its support -- it is not a discovery, it is an annotation.
    ann_ends = {ref.tx[t].tts for t in chain_tx}
    keep = []
    for li in order:
        sel = read_idx[labels == li]
        nm, _nc, has = support(sel)
        strong = sizes[li] >= params.min_end_reads and (
            not has or nm >= params.min_end_umis
        )
        if any(abs(int(peaks[li]) - e) <= params.max_3p_diff for e in ann_ends):
            strong = True
        if strong:
            keep.append(int(li))
    if not keep:
        keep = [int(order[0])]
    keep_set = set(keep)
    kept_peaks = np.array([peaks[li] for li in keep], np.int64)
    remap = labels.copy()
    for li in range(peaks.size):
        if li in keep_set:
            continue
        j = int(np.argmin(np.abs(kept_peaks - peaks[li])))
        remap[labels == li] = keep[j]

    group_reads = int(read_idx.size)
    group_mols = 0
    for li in keep:
        nm, _nc, has = support(read_idx[remap == li])
        if has:
            group_mols += nm

    def evidence_at(pos: int):
        return polya_evidence(
            contig, pos, strand, genome, atlas,
            motif_window=mono_params.polya_motif_window,
            perc_a_window=mono_params.perc_a_window,
            max_perc_a=mono_params.max_perc_a_downstream,
            motif_min_dist=mono_params.polya_motif_min_dist,
            motif_max_dist=mono_params.polya_motif_max_dist,
            motif_strict=mono_params.polya_motif_strict,
            require_atlas_when_available=mono_params.require_atlas_when_available,
        )

    # -- stage A: decide the final 3' end of every kept peak ----------------
    decided: List[Tuple[int, str, str, str, List[str], np.ndarray]] = []
    for li in keep:
        sel_mask = remap == li
        sel = read_idx[sel_mask]
        if sel.size == 0:
            continue
        nm, _nc, _has = support(sel)
        three = weighted_median(tts[sel_mask], wts[sel_mask])
        observed_three = three
        flags: List[str] = []
        cat = "NOVEL"
        end3_from = ""
        # 1. an annotated end of some transcript carrying THIS chain
        near = _nearest_by_tts(ref, chain_tx, three)
        ntx = ref.tx.get(near) if near else None
        if ntx is not None and abs(signed_5p_offset(three, ntx.tts, strand)) <= params.max_3p_diff:
            three, cat, end3_from = ntx.tts, "FSM", near
            decided.append((three, cat, near, end3_from, flags, sel_mask))
            continue
        parent = default_tx
        if dtx is None:
            decided.append((three, cat, parent or "", end3_from, flags, sel_mask))
            continue
        # 2. an annotated end of any transcript of the gene -- the chain is that
        #    of one isoform, the cleavage site that of another.  Both parts are
        #    annotated; the pairing is not, and the flag says so.
        snapped = False
        if params.snap_to_gene_level_3p_ends and gene_id in ref.genes:
            d_gene, site = ref.nearest_tts(gene_id, three)
            if abs(d_gene) <= params.max_3p_diff and site >= 0:
                three, cat, snapped = site, "end3_annotated", True
                end3_from = _tts_owners(ref, gene_id, site)
                flags.append("annotated_alt_3p_end")
                flags.append("novel_pairing_of_annotated_ends")
        if not snapped:
            share = ((nm / group_mols) if group_mols
                     else (sel.size / max(group_reads, 1)))
            rt = readthrough_gene(ref, contig, strand, gene_id, three)
            tail = tail_evidence(reads, sel, wts[sel_mask], params)
            polya_ok = evidence_at(three)["polya_supported"]
            if (not polya_ok and params.polya_tail_gates_novel_end
                    and tail["tail_molecule_frac"] >= params.min_tail_molecule_frac):
                polya_ok = True
                flags.append("polya_from_read_tail")
            if share < params.min_end_share:
                three, cat = _fallback(
                    observed_three, dtx, strand, params, flags,
                    "3p_variant_rejected_low_share")
            elif params.require_polya_evidence_for_utr_variant and not polya_ok:
                three, cat = _fallback(
                    observed_three, dtx, strand, params, flags,
                    "3p_variant_rejected_no_polya")
            elif rt is not None:
                cat = "READTHROUGH"
                flags.append("3p_shifted")
                flags.append(f"readthrough_into={rt}")
            else:
                cat = "end3_novel"
                flags.append("3p_shifted")
        decided.append((three, cat, parent or "", end3_from, flags, sel_mask))

    # -- stage B: peaks that landed on the same 3' end are one transcript ---
    merged: Dict[int, List] = {}
    for three, cat, parent, end3_from, flags, sel_mask in decided:
        e = merged.get(three)
        if e is None:
            merged[three] = [cat, parent, end3_from, list(flags), sel_mask.copy()]
        else:
            if RANK.get(cat, 9) < RANK.get(e[0], 9):
                e[0], e[1], e[2] = cat, parent, end3_from
            e[3] = sorted(set(e[3]) | set(flags))
            e[4] |= sel_mask

    # -- stage C: 5' end, category refinement, emit ------------------------
    ordered = sorted(merged.items(), key=lambda kv: -float(wts[kv[1][4]].sum()))
    models: List[TranscriptModel] = []
    variant_rank = 0
    for three, (cat, parent, end3_from, flags, sel_mask) in ordered:
        sel = read_idx[sel_mask]
        if sel.size == 0:
            continue
        nm, nc, _has = support(sel)
        five = five_prime_end(tss[sel_mask], three, strand, params.five_prime_percentile)
        ptx = ref.tx.get(parent) if parent else None
        d3 = signed_5p_offset(three, ptx.tts, strand) if ptx else None
        d5 = signed_5p_offset(five, ptx.tss, strand) if ptx else None
        n_trunc = 0
        if ptx is not None:
            if abs(d5) <= params.max_5p_diff:
                five = ptx.tss
            elif d5 < 0:
                cat, flags = _classify_extension(
                    ref, gene_id, cat, flags, five, ptx.tss, strand, params
                )
            else:
                # 5' truncation.  Absorbed by default: emitting it would add a
                # model per degradation peak, and degradation peaks are not
                # isoforms.  The mass is kept and counted, not thrown away.
                if params.allow_5p_truncation_models:
                    cat = "end5_truncated"
                    flags = sorted(set(flags) | {"5p_truncation"})
                else:
                    five = ptx.tss
                    flags = sorted(set(flags) | {"5p_truncation_absorbed"})
        # how much of this model's own mass starts downstream of its 5' end.
        # The truncated reads are counted here instead of becoming models.
        n_trunc = int(
            np.sum(
                np.array(
                    [signed_5p_offset(int(t), five, strand) for t in tss[sel_mask]]
                )
                > params.max_5p_diff
            )
        )
        if cat == params.unresolved_3p_category and not params.emit_unresolved_3p:
            continue
        if cat in ("end3_novel", params.unresolved_3p_category):
            # an unresolved peak counts against the same per-parent budget as a
            # novel end -- it is one more terminal-end model on this chain, and
            # letting it in for free is how a group grows a ladder of them
            variant_rank += 1
            if variant_rank > params.max_utr_variants_per_parent:
                continue
        ev = evidence_at(three)
        exons = build_exons(chain, strand, five, three)
        m = TranscriptModel(
            model_id="",
            contig=contig,
            strand=strand,
            exons=exons,
            chain=chain,
            category=cat,
            reads=sel,
            n_reads=int(sel.size),
            n_mols=nm,
            n_cells=nc,
            parent_tx=parent or None,
            gene_id=gene_id,
            gene_name=ref.genes[gene_id].gene_name if gene_id in ref.genes else "",
            flags=list(flags),
        )
        m.evidence.update(ev)
        m.evidence["dist_to_ref_tts"] = d3
        m.evidence["dist_to_ref_tss"] = d5
        m.evidence["end3_from"] = end3_from
        m.evidence["utr_variant_rank"] = variant_rank
        m.evidence["n_3p_peaks"] = len(ordered)
        m.evidence["three_prime_dispersion"] = float(np.std(tts[sel_mask]))
        m.evidence["n_mol_equivalents"] = round(float(wts[sel_mask].sum()), 3)
        m.evidence.update(tail_evidence(reads, sel, wts[sel_mask], params))
        # how far this model's 3' end sits from where its own reads end.
        # Non-zero means annotation moved it; it is the number the NSL1 class of
        # failure shows up in, and it is now in every row rather than inferable
        # only from flags.
        m.evidence["tts_shift_from_reads"] = int(
            signed_5p_offset(three, weighted_median(tts[sel_mask], wts[sel_mask]), strand)
        )
        m.evidence["n_5p_truncated_reads"] = n_trunc
        offs = np.array(
            [signed_5p_offset(int(t), three, strand) for t in tss[sel_mask]], float
        )
        for p in (10, 50, 90):
            m.evidence[f"five_prime_p{p}"] = float(np.percentile(offs, p))
        models.append(m)

    return _dedupe(models)


def tail_evidence(
    reads: ContigReads, sel: np.ndarray, w: np.ndarray, params: EndParams
) -> Dict[str, float]:
    """Direct polyA-tail evidence for one 3' peak, from the reads' own clips.

    New in 0.1.17, and RECORDED rather than enforced by default.  Everything the
    0.1.16 3' path called "polyA evidence" was genomic -- a PAS motif, downstream
    A-richness, an optional atlas -- which is a statement about the locus, not
    about the molecules that were primed there.  Those cannot distinguish a real
    cleavage site from an A-rich stretch that primes internally in every library.
    The soft-clipped tail can.

    Returns zeros when ``AlignmentGates.measure_polya_tail`` was off, which is
    indistinguishable here from "no read had a tail" -- so a caller that acts on
    these numbers must confirm the flag first.  ``EndParams`` gates that do act
    on them are off by default for exactly this reason.
    """
    if sel.size == 0 or reads.tail_len.size == 0:
        return {"tail_molecule_frac": 0.0, "n_tail_reads": 0.0, "median_tail_len": 0.0}
    tl = reads.tail_len[sel].astype(np.int64)
    tf = reads.tail_frac[sel].astype(float)
    # 0.1.17 required BOTH a terminal A-run and >=75% A over the clip's first 60
    # bases. Those are not independent: a genuine 12 bp tail inside a 60 bp clip
    # that also holds adapter or TSO sequence scores frac ~0.2 and failed, so the
    # frac term silently demanded a ~45 bp tail and biased every tail number
    # downward. The run is the evidence -- a template-free A-run outside the
    # alignment is the tail -- so it alone decides, and frac is kept as an
    # optional extra check, off by default.
    ok = tl >= params.min_polya_tail_len
    if params.min_polya_tail_frac > 0:
        ok = ok & (tf >= params.min_polya_tail_frac)
    tot = float(w.sum())
    return {
        "tail_molecule_frac": float(w[ok].sum() / tot) if tot > 0 else 0.0,
        "n_tail_reads": float(int(ok.sum())),
        "median_tail_len": float(np.median(tl)) if tl.size else 0.0,
    }


def _fallback(
    observed: int,
    dtx,
    strand: str,
    params: EndParams,
    flags: List[str],
    reason: str,
) -> Tuple[int, str]:
    """Where a rejected 3' peak's model actually ends.  The 0.1.17 fix.

    Rejecting a peak means "there is not enough evidence to call a new cleavage
    site here".  It does not mean "the molecules really ended somewhere else",
    and through 0.1.16 the code treated it as if it did: the reads were
    re-labelled with the default parent's annotated end at any distance.  At
    NSL1 that moved the model 11.6 kb past every read supporting it.

    Snapping remains correct when the annotated end is close -- that is
    annotation refining an observed coordinate.  Beyond ``max_3p_fallback_dist``
    it would be annotation overriding the observation, so the observed
    coordinate stands and the model says it is unresolved.
    """
    d = abs(signed_5p_offset(observed, dtx.tts, strand))
    flags.append(reason)
    if d <= params.max_3p_fallback_dist:
        return dtx.tts, "FSM"
    flags.append(f"3p_fallback_refused_dist={d}")
    # By construction this coordinate matched no annotated end of this chain or
    # of the gene -- both snap tests ran first and failed. So it IS a novel 3'
    # end; what it lacks is polyA corroboration, and the flag says which.
    flags.append("novel_3p_end_no_polya_support")
    return observed, params.unresolved_3p_category


def _classify_extension(
    ref: ReferenceIndex,
    gene_id: Optional[str],
    cat: str,
    flags: List[str],
    five: int,
    ref_tss: int,
    strand: str,
    params: EndParams,
) -> Tuple[str, List[str]]:
    """A 5' end reaching upstream of the parent's start is one of two things.

    The discriminator is whether the extension **crosses an annotated junction**.

    * It does not -- some isoform simply has a longer first exon, and the model
      is an ordinary alternative or extended start.
    * It does -- the longer isoform splices there and this molecule read
      straight through.  That is a retained first intron, not a new start site,
      and calling it a novel TSS would be wrong in the specific way that
      inflates diversity: the same molecule counted as a new promoter.
    """
    lo, hi = (five, ref_tss) if five < ref_tss else (ref_tss, five)
    crossed = ref.gene_junctions_within(gene_id, lo, hi)
    f = set(flags)
    if crossed:
        # a structural statement, so it outranks whatever the 3' end decided
        f.add("5p_extension_crosses_annotated_junction")
        f.add(f"unspliced_introns={len(crossed)}")
        return "IR", sorted(f)
    f.add("5p_extension")
    if cat != "FSM":
        # the 3' end already found something more specific to say, and it is the
        # informative end under oligo-dT; do not overwrite it with the 5' one
        return cat, sorted(f)
    if gene_id and gene_id in ref.genes:
        d_tss, site = ref.nearest_tss(gene_id, five)
        if abs(d_tss) <= params.max_5p_diff and site >= 0:
            f.add("annotated_alt_5p_end")
            return "end5_alt", sorted(f)
    return "end5_extended", sorted(f)


# ---------------------------------------------------------------------- #
def annotate_similarity(
    models: Sequence[TranscriptModel],
    ref: ReferenceIndex,
    ir_min_flank: int = 10,
) -> None:
    """Attach the nearest reference transcript and the structural diff.

    Annotation only.  Nothing here changes which models exist or how much read
    mass they carry; it changes what a row of ``models.tsv`` tells you.  The one
    exception is intron retention, which is a *classification* the chain lookup
    cannot express: a chain that is a reference chain minus introns the molecule
    read through is not a novel combination of junctions, and calling it NIC
    hides a common and well-understood event behind a novelty label.
    """
    for m in models:
        if m.n_exons < 2:
            # the mono-exonic track already chooses its parent by nearest
            # annotated 3' end, which *is* a similarity search; surface it in
            # the same column so the table has no holes
            m.associated_tx = m.associated_tx or m.parent_tx
            continue
        match = nearest_transcript(
            ref, m.gene_id, m.chain, m.exons, m.strand, ir_min_flank,
            require_retention=(m.category == "IR"),
        )
        if match is None:
            continue
        m.associated_tx = match.tx_id
        m.structural_diff = match.describe()
        m.n_equally_close = match.n_equally_close
        m.is_utr_only = match.is_utr_only
        m.evidence["n_shared_junctions"] = match.n_shared
        m.evidence["n_novel_junctions_vs_ref"] = match.n_model_only
        m.evidence["n_retained_introns"] = match.n_retained
        if match.retained:
            m.evidence["ir_introns"] = ";".join(f"{a}-{b}" for a, b in match.retained)
            # Retention deep in the 3' end is a confident observation: the
            # molecule was captured by its polyA tail and sequenced through it.
            # Retention of the 5'-most intron is the one confounded with
            # incomplete reverse transcription, so it is marked rather than
            # trusted silently.
            rank = retained_intron_rank(ref.tx[match.tx_id].chain, match.retained)
            m.evidence["ir_intron_rank"] = rank
            if rank == 0:
                m.flags = sorted(set(m.flags) | {"ir_at_5p_most_intron"})
        if match.n_equally_close > 1:
            m.flags = sorted(set(m.flags) | {f"n_equally_close={match.n_equally_close}"})
        if match.is_pure_retention and m.category in ("NIC", "NNC", "NOVEL"):
            m.category = "IR"
            m.flags = sorted(set(m.flags) | {"intron_retention"})
        elif match.retained and m.category == "FSM":
            # the chain matched a `retained_intron` GENCODE model exactly, so it
            # is an FSM and an IR event at once.  Both facts are true; recording
            # only the first would hide a large, prep-sensitive class.
            m.flags = sorted(set(m.flags) | {"ir_event"})


def _dedupe(models: List[TranscriptModel]) -> List[TranscriptModel]:
    by_struct: Dict[Tuple, TranscriptModel] = {}
    for m in models:
        key = (m.strand, tuple(m.exons))
        cur = by_struct.get(key)
        if cur is None:
            by_struct[key] = m
        else:
            cur.reads = np.concatenate([cur.reads, m.reads])
            cur.n_reads += m.n_reads
            cur.n_mols += m.n_mols
            cur.n_cells = max(cur.n_cells, m.n_cells)
            cur.flags = sorted(set(cur.flags) | set(m.flags))
    return list(by_struct.values())
