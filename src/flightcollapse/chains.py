"""Intron-chain grouping and the suffix resolution that fixes the isoseq bug.

Background measurement (BD144_kin, chr1+chr12+chrX): **57.8% of spliced read
mass sits on chains that are a proper suffix of some longer observed chain.**
Merging on the suffix relation alone is therefore not an edge case, it is the
majority of the data -- and the direction of the relation is not what a naive
rule assumes:

===========================  =========
class (threshold 1 read)     read mass
===========================  =========
debris      (S <= P/3)          3.7%
dominant    (S >= 3P)          53.0%
comparable                      1.0%
===========================  =========

So most suffix chains carry *more* reads than their own "parent": the suffix is
the real isoform and the parent is a rare 5'-extended variant.  Merging those
away is exactly the isoseq failure.  Only the small `debris` class is genuine
5'-truncation.

Two hard rules, both enforced here:

1. **The representative is the maximal chain, never the minimal.**  Merges only
   ever go from a shorter chain into a longer one, so the surviving model is by
   construction the longest member of its group.
2. **The empty chain never enters this logic.**  Mono-exonic reads are routed to
   :mod:`flightcollapse.monoexon`.  The empty chain is a vacuous suffix of every
   chain, which is how a mono-exonic stub absorbed 41% of read mass.

The damage metric to watch is not the debris chain's own read mass but the
*parent* mass dragged onto a truncated model: a 20-read suffix of a 5,000-read
parent contributes 20 reads to `debris` but mis-models 5,020 reads if merged the
wrong way.  ``reads_restructured_by_merge`` in the QC report is that number.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from .config import ChainParams
from .intervals import Chain, is_suffix
from .molecules import NO_CELL, NO_UMI, n_cells, n_molecules
from .reads import ContigReads
from .reference import ReferenceIndex


@dataclass
class ChainGroup:
    gid: int
    contig: str
    strand: str
    chain: Chain
    reads: np.ndarray                       # read indices into ContigReads
    n_reads: int = 0
    n_mols: int = 0
    n_cells: int = 0
    tx_ids: List[str] = field(default_factory=list)   # exact annotated matches
    gene_id: Optional[str] = None
    #: gene the positional assignment would have chosen, when it disagrees with
    #: the one the exact chain match implies (paralogues, pseudogenes, overlaps)
    gene_conflict: Optional[str] = None
    n_novel_junctions: int = 0
    n_decoy_junctions: int = 0
    min_junction_posterior: float = 1.0
    frac_junctions_annotated: float = 1.0
    # resolution outcome
    merged_into: Optional[int] = None
    merge_class: str = ""                   # debris | annotated_parent | ambiguous | ...
    ambiguous: bool = False
    n_parents: int = 0
    eff_n_parents: float = 0.0
    ratio_to_top_parent: float = float("nan")
    #: model 5' end (high percentile of member reads) sits at an annotated TSS
    #: of the same gene -- evidence for a real alternative start rather than
    #: 5' degradation
    tss_evidence: bool = False

    @property
    def annotated(self) -> bool:
        return bool(self.tx_ids)

    @property
    def support(self) -> int:
        return self.n_mols if self.n_mols > 0 else self.n_reads

    def __len__(self) -> int:
        return len(self.chain)


# ---------------------------------------------------------------------- #
def build_chain_groups(
    reads: ContigReads,
    contig: str,
    ref: ReferenceIndex,
    read_ok: Optional[np.ndarray] = None,
    umi_hamming: int = 1,
) -> Tuple[List[ChainGroup], np.ndarray]:
    """Group spliced reads by exact (post-snapping) intron chain.

    Returns the groups and a per-read array of group id (-1 for mono-exonic or
    excluded reads).
    """
    gid_of_read = np.full(reads.n, -1, np.int64)
    if reads.n == 0:
        return [], gid_of_read

    counts = np.diff(reads.joff)
    spliced = counts > 0
    if read_ok is not None:
        spliced = spliced & read_ok

    buckets: Dict[bytes, List[int]] = defaultdict(list)
    donor, acceptor, joff, strand = reads.donor, reads.acceptor, reads.joff, reads.strand
    for i in np.flatnonzero(spliced):
        a, b = int(joff[i]), int(joff[i + 1])
        key = bytes([1 if strand[i] > 0 else 0]) + donor[a:b].tobytes() + acceptor[a:b].tobytes()
        buckets[key].append(int(i))

    groups: List[ChainGroup] = []
    for gid, (_key, idx) in enumerate(buckets.items()):
        arr = np.array(idx, np.int64)
        i0 = int(arr[0])
        st = reads.strand_char(i0)
        ch = reads.chain(i0)
        pairs = [(int(reads.cell[j]), int(reads.umi[j])) for j in arr]
        has_umi = any(c != int(NO_CELL) and u != int(NO_UMI) for c, u in pairs)
        g = ChainGroup(
            gid=gid,
            contig=contig,
            strand=st,
            chain=ch,
            reads=arr,
            n_reads=int(arr.size),
            n_mols=n_molecules(pairs, umi_hamming) if has_umi else 0,
            n_cells=n_cells(pairs) if has_umi else 0,
        )
        key = (contig, st)
        g.tx_ids = list(ref.transcripts_for_chain(key, ch))
        groups.append(g)
        gid_of_read[arr] = gid
    return groups, gid_of_read


def annotate_groups_with_junctions(
    groups: Sequence[ChainGroup],
    jinfo: Dict[Tuple[str, int, int], Dict[str, float]],
) -> None:
    """Attach per-junction novelty statistics to each chain group."""
    for g in groups:
        n_novel = n_decoy = 0
        post = 1.0
        n_ann = 0
        for d, a in g.chain:
            rec = jinfo.get((g.strand, int(d), int(a)))
            if rec is None:
                n_novel += 1
                post = min(post, 0.0)
                continue
            if rec["annotated"]:
                n_ann += 1
            else:
                n_novel += 1
                p = rec.get("posterior_real", float("nan"))
                post = min(post, 1.0 if np.isnan(p) else float(p))
                if rec.get("decoy", 0):
                    n_decoy += 1
        g.n_novel_junctions = n_novel
        g.n_decoy_junctions = n_decoy
        g.min_junction_posterior = post
        g.frac_junctions_annotated = n_ann / max(len(g.chain), 1)


def assign_genes(groups: Sequence[ChainGroup], ref: ReferenceIndex) -> Dict[str, int]:
    """Attach a gene to each chain group.

    An exact chain match decides, always -- that is the reference-anchored rule.
    But recent paralogues, processed pseudogenes and overlapping annotations
    make the chain match and the positional assignment disagree, and letting the
    chain win *silently* means nobody ever finds out how often.  The
    disagreement is recorded on the group and counted here.
    """
    stats = {"gene_conflicts": 0}
    for g in groups:
        gj = g.chain if g.strand == "+" else g.chain[::-1]
        exon_proxy = [(int(gj[0][0]) - 1, int(gj[0][0]))] + [
            (int(gj[k][1]), int(gj[k + 1][0])) for k in range(len(gj) - 1)
        ] + [(int(gj[-1][1]), int(gj[-1][1]) + 1)]
        by_overlap = ref.assign_gene(g.contig, exon_proxy, g.strand)
        if g.tx_ids:
            g.gene_id = ref.tx[g.tx_ids[0]].gene_id
            if by_overlap and by_overlap != g.gene_id:
                g.gene_conflict = by_overlap
                stats["gene_conflicts"] += 1
            continue
        g.gene_id = by_overlap
    return stats


# ---------------------------------------------------------------------- #
# suffix resolution
# ---------------------------------------------------------------------- #
class _DSU:
    def __init__(self, n: int) -> None:
        self.p = list(range(n))

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union_into(self, child: int, parent: int) -> None:
        self.p[self.find(child)] = self.find(parent)


def find_parents(
    groups: Sequence[ChainGroup], candidate: np.ndarray
) -> Dict[int, List[int]]:
    """``gid -> [maximal parent gids]`` among candidates, per (strand, 3'-most intron).

    Only chains sharing the 3'-most intron can be in a suffix relation, which
    keeps this near-linear.  Parents are then reduced to *maximal* ones -- the
    raw parent counts in the earlier census were inflated because parents are
    themselves nested under the same relation.
    """
    by_last: Dict[Tuple[str, Tuple[int, int]], List[int]] = defaultdict(list)
    for g in groups:
        if candidate[g.gid] and g.chain:
            by_last[(g.strand, g.chain[-1])].append(g.gid)
    for k in by_last:
        by_last[k].sort(key=lambda i: len(groups[i].chain))

    out: Dict[int, List[int]] = {}
    for g in groups:
        if not candidate[g.gid] or not g.chain:
            continue
        cands = by_last.get((g.strand, g.chain[-1]), ())
        parents = [
            p for p in cands
            if p != g.gid and is_suffix(g.chain, groups[p].chain)
        ]
        if not parents:
            continue
        maximal = [
            p for p in parents
            if not any(
                q != p and is_suffix(groups[p].chain, groups[q].chain) for q in parents
            )
        ]
        out[g.gid] = maximal or parents
    return out


def resolve_suffixes(
    groups: List[ChainGroup],
    candidate: np.ndarray,
    params: ChainParams,
) -> Dict[str, float]:
    """Decide, for every suffix relation, whether to merge and in which direction.

    Mutates ``groups`` (``merged_into``, ``merge_class``, ``ambiguous``, ...) and
    returns diagnostics including the read mass whose structure a merge changed.
    """
    parents = find_parents(groups, candidate)
    dsu = _DSU(len(groups))
    stats = {
        "n_suffix_chains": 0,
        "n_merged": 0,
        "reads_merged": 0,
        "reads_restructured_by_merge": 0,
        "n_dominant_protected": 0,
        "n_ambiguous": 0,
    }

    # process shortest first so that transitive merges land on maximal chains
    order = sorted(parents, key=lambda gid: len(groups[gid].chain))
    for gid in order:
        g = groups[gid]
        ps = parents[gid]
        stats["n_suffix_chains"] += 1
        sup = np.array([max(groups[p].support, 1) for p in ps], float)
        g.n_parents = len(ps)
        g.eff_n_parents = float(1.0 / np.sum((sup / sup.sum()) ** 2))
        top = ps[int(np.argmax(sup))]
        top_sup = max(groups[top].support, 1)
        ratio = g.support / top_sup
        g.ratio_to_top_parent = ratio
        frac_top = float(sup.max() / sup.sum())

        ann_parents = [p for p in ps if groups[p].annotated]

        target: Optional[int] = None
        why = ""

        # 1. a chain that IS an annotated transcript is never merged away.
        if g.annotated:
            g.merge_class = "annotated_kept"
            continue

        # 2. reference-anchored rule: an unannotated suffix of an ANNOTATED
        #    chain is 5'-truncation of that transcript, whatever the support
        #    ratio says -- unless its 5' end sits at an annotated TSS, which is
        #    independent evidence of a real alternative start.
        if params.merge_unannotated_suffix_into_annotated_parent and ann_parents:
            if params.require_tss_evidence_to_keep_unannotated_suffix and g.tss_evidence:
                g.merge_class = "alt_tss_kept"
                stats["n_dominant_protected"] += 1
                continue
            ann_sup = np.array([max(groups[p].support, 1) for p in ann_parents], float)
            best_ann = ann_parents[int(np.argmax(ann_sup))]
            target = best_ann
            if len(ann_parents) == 1:
                why = "truncation_of_annotated_parent"
            else:
                # Several annotated transcripts share this 3' intron block, so
                # the truncated read is consistent with all of them.  It is
                # still a truncation of a known transcript, not a new isoform:
                # emitting it as one is exactly the over-calling this tool
                # exists to avoid.  Assign to the best-supported parent and
                # flag it, so the ambiguity stays visible in group_chains.tsv.
                why = "truncation_of_annotated_parents_ambiguous"
                g.ambiguous = True
                stats["n_ambiguous_resolved_to_annotated"] = (
                    stats.get("n_ambiguous_resolved_to_annotated", 0) + 1
                )

        # 3. neither side annotated: decide by support direction only.
        if target is None:
            if ratio >= params.dominant_ratio:
                g.merge_class = "dominant_kept"
                stats["n_dominant_protected"] += 1
                continue
            if ratio <= params.debris_ratio and (len(ps) == 1 or frac_top >= 0.9):
                target, why = top, "debris"
            elif params.annotation_breaks_ties and len(ann_parents) == 1:
                target, why = ann_parents[0], "annotated_parent_tiebreak"

        if target is None:
            g.ambiguous = True
            g.merge_class = "ambiguous"
            stats["n_ambiguous"] += 1
            continue

        g.merged_into = target
        g.merge_class = why
        dsu.union_into(gid, target)
        stats["n_merged"] += 1
        stats["reads_merged"] += g.n_reads
        stats["reads_restructured_by_merge"] += g.n_reads

    # resolve transitive targets
    for g in groups:
        if g.merged_into is not None:
            g.merged_into = dsu.find(g.merged_into)
    return stats


def attach_subthreshold(
    groups: List[ChainGroup],
    candidate: np.ndarray,
    params: ChainParams,
) -> Dict[str, int]:
    """Re-home chains that failed the support threshold.

    A sub-threshold chain is merged into a surviving chain that it is a suffix
    of, when that destination is unambiguous.  Otherwise its reads are marked
    unassigned rather than being promoted into a novel isoform -- silently
    inventing a model out of three reads is how transcriptional diversity gets
    overestimated.
    """
    surviving = np.array(
        [candidate[g.gid] and g.merged_into is None for g in groups], bool
    )
    by_last: Dict[Tuple[str, Tuple[int, int]], List[int]] = defaultdict(list)
    for g in groups:
        if surviving[g.gid] and g.chain:
            by_last[(g.strand, g.chain[-1])].append(g.gid)

    stats = {"n_rehomed": 0, "reads_rehomed": 0, "n_dropped": 0, "reads_dropped": 0}
    for g in groups:
        if candidate[g.gid] or g.merged_into is not None or not g.chain:
            continue
        cands = [
            p for p in by_last.get((g.strand, g.chain[-1]), ())
            if is_suffix(g.chain, groups[p].chain)
        ]
        if not cands:
            g.merge_class = "dropped_no_parent"
            stats["n_dropped"] += 1
            stats["reads_dropped"] += g.n_reads
            continue
        sup = np.array([max(groups[p].support, 1) for p in cands], float)
        top = cands[int(np.argmax(sup))]
        ann = [p for p in cands if groups[p].annotated]
        if len(cands) == 1 or sup.max() / sup.sum() >= 0.9:
            g.merged_into, g.merge_class = top, "subthreshold_merge"
        elif len(ann) == 1:
            g.merged_into, g.merge_class = ann[0], "subthreshold_annotated_parent"
        else:
            g.merge_class = "dropped_ambiguous"
            stats["n_dropped"] += 1
            stats["reads_dropped"] += g.n_reads
            continue
        stats["n_rehomed"] += 1
        stats["reads_rehomed"] += g.n_reads
    return stats


def mark_candidates(
    groups: Sequence[ChainGroup],
    params: ChainParams,
    gene_support: Optional[Dict[str, int]] = None,
    lfdr: Optional[Dict[int, float]] = None,
) -> np.ndarray:
    """Which chains are allowed to become models in their own right.

    Annotated chains pass on annotation alone (this is a reference-anchored
    tool).  Novel chains must clear absolute support, depth-relative support
    within their gene, and -- when the calibrated model is available -- the
    local-FDR ceiling.
    """
    out = np.zeros(len(groups), bool)
    use_mols = any(g.n_mols > 0 for g in groups)
    for g in groups:
        if g.annotated:
            out[g.gid] = True
            continue
        if use_mols and g.n_mols < params.min_chain_umis:
            continue
        if g.n_reads < params.min_chain_reads:
            continue
        if gene_support and g.gene_id:
            tot = gene_support.get(g.gene_id, 0)
            if tot > 0 and (g.support / tot) < params.min_novel_chain_share:
                continue
        if lfdr is not None and params.max_chain_lfdr is not None:
            v = lfdr.get(g.gid, float("nan"))
            if not np.isnan(v) and v > params.max_chain_lfdr:
                continue
        out[g.gid] = True
    return out
