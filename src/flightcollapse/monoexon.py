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
from .ends import cluster_3p_ends, five_prime_end
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

    if n_reads < min_reads or (has_umi and n_mols < params.min_umis):
        stats["rejected_low_support"] += 1
        stats["reads_rejected_low_support"] += n_reads
        return None
    if ev["internal_priming"]:
        stats["rejected_internal_priming"] += 1
        stats["reads_rejected_internal_priming"] += n_reads
        return None
    if need_polya and not ev["polya_supported"]:
        stats["rejected_no_polya"] += 1
        stats["reads_rejected_no_polya"] += n_reads
        return None

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
