"""Splice-junction catalogue: census, snapping, features, curation.

Rationale (from the chain census on BD144_kin, chr1+chr12+chrX):

* 83.5% of *distinct* junctions are novel but they carry only ~1.2% of junction
  read mass.  The splice-site catalogue underlying the data is small, clean and
  almost entirely annotated.
* Fuzzy tolerance is not a lever -- 0 -> 20 bp removes only 15% of chains.
  Junction *support* filtering is the lever; junction *snapping* is not.
* Therefore: filter at the junction level first, rebuild chains from the
  surviving catalogue, and only then threshold chains.  A chain with 40 reads
  containing one spurious junction survives a chain threshold of 10; it does
  not survive a junction threshold.

Snapping is **annotation-first**: a read junction within the fuzzy tolerance of
an annotated splice site is moved onto the annotated site, before any
read-derived clustering happens.  That converts a large share of "novel"
junctions into known ones without inventing anything.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import JunctionParams
from .genome import Genome
from .molecules import NO_CELL, NO_UMI
from .reads import ContigReads
from .reference import ReferenceIndex


def _pair_hash(cell: np.ndarray, umi: np.ndarray) -> np.ndarray:
    """Stable 64-bit hash of (cell, umi); sentinel-preserving."""
    c = cell.astype(np.uint64)
    u = umi.astype(np.uint64)
    h = (c * np.uint64(0x9E3779B97F4A7C15)) ^ (u + np.uint64(0x165667B19E3779F9))
    h ^= h >> np.uint64(29)
    h *= np.uint64(0xBF58476D1CE4E5B9)
    h ^= h >> np.uint64(32)
    bad = (cell == NO_CELL) | (umi == NO_UMI)
    return np.where(bad, np.uint64(0), h)


def _distinct_per_group(gid: np.ndarray, val: np.ndarray, n_groups: int) -> np.ndarray:
    """Number of distinct ``val`` per group id, ignoring val == 0."""
    out = np.zeros(n_groups, np.int64)
    if gid.size == 0:
        return out
    keep = val != 0
    if not keep.any():
        return out
    g, v = gid[keep], val[keep]
    order = np.lexsort((v, g))
    g, v = g[order], v[order]
    new = np.empty(g.size, bool)
    new[0] = True
    new[1:] = (g[1:] != g[:-1]) | (v[1:] != v[:-1])
    np.add.at(out, g[new], 1)
    return out


class JunctionCensus:
    """Unique junctions on one contig with read / molecule / cell support."""

    def __init__(
        self,
        strand: np.ndarray,
        donor: np.ndarray,
        acceptor: np.ndarray,
        n_reads: np.ndarray,
        n_mols: np.ndarray,
        n_cells: np.ndarray,
        jid_of_record: np.ndarray,
    ) -> None:
        self.strand = strand
        self.donor = donor
        self.acceptor = acceptor
        self.n_reads = n_reads
        self.n_mols = n_mols
        self.n_cells = n_cells
        self.jid_of_record = jid_of_record

    @property
    def size(self) -> int:
        return int(self.donor.size)

    @classmethod
    def from_reads(cls, reads: ContigReads) -> "JunctionCensus":
        if reads.donor.size == 0:
            z = np.zeros(0, np.int64)
            return cls(np.zeros(0, np.int8), z, z, z, z, z, np.zeros(0, np.int64))
        counts = np.diff(reads.joff)
        jstrand = np.repeat(reads.strand, counts).astype(np.int8)
        key = np.stack([jstrand.astype(np.int64), reads.donor, reads.acceptor], axis=1)
        uniq, jid = np.unique(key, axis=0, return_inverse=True)
        jid = jid.astype(np.int64)
        n = uniq.shape[0]
        n_reads = np.bincount(jid, minlength=n).astype(np.int64)
        molh = np.repeat(_pair_hash(reads.cell, reads.umi), counts)
        cellh = np.repeat(
            np.where(reads.cell == NO_CELL, 0, reads.cell.astype(np.uint64) + 1), counts
        )
        n_mols = _distinct_per_group(jid, molh, n)
        n_cells = _distinct_per_group(jid, cellh, n)
        return cls(
            uniq[:, 0].astype(np.int8), uniq[:, 1], uniq[:, 2],
            n_reads, n_mols, n_cells, jid,
        )

    def strand_chars(self) -> np.ndarray:
        return np.where(self.strand > 0, "+", "-")


# ---------------------------------------------------------------------- #
# exon-internal deletions
# ---------------------------------------------------------------------- #
def promote_internal_gaps(
    reads: ContigReads,
    contig: str,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    params: JunctionParams,
) -> Dict[str, int]:
    """Turn short introns the aligner spelled as ``D`` into real junctions.

    Runs *before* the junction census, so a promoted gap is snapped, featurised,
    scored and curated exactly like every other junction -- it is not waved
    through because of how it was spelled.

    Two ways in, and nothing else:

    * the gap is an annotated junction of the reference.  Annotation decides;
      no motif check is needed or wanted.
    * the gap recurs across reads **and** has a canonical motif at its
      boundaries.  A one-off gap with no motif is alignment noise, and
      promoting it would manufacture a junction that curation then rejects,
      which would throw the whole read away.

    Everything not promoted stays recorded on the read, so a model built from
    reads carrying unexplained internal deletions can say so.
    """
    stats = {
        "gap_records": int(reads.gap_start.size),
        "gap_distinct": 0,
        "gap_promoted_annotated": 0,
        "gap_promoted_motif": 0,
        "gap_records_promoted": 0,
    }
    if not params.promote_internal_deletions or reads.gap_start.size == 0:
        return stats
    census = reads.gap_census()
    stats["gap_distinct"] = len(census)
    promote = set()
    for (sc, a, b), n in census.items():
        if (b - a) > params.max_internal_deletion_len:
            continue
        if ref.is_annotated_junction((contig, sc), (a, b)):
            promote.add((sc, a, b))
            stats["gap_promoted_annotated"] += 1
            continue
        if genome is None or n < params.min_internal_deletion_reads:
            continue
        if genome.motif_class(genome.splice_motif(contig, a, b, sc)) == "canonical":
            promote.add((sc, a, b))
            stats["gap_promoted_motif"] += 1
    stats["gap_records_promoted"] = reads.promote_gaps(promote)
    return stats


def reads_with_unexplained_gap(reads: ContigReads, promoted: bool = True) -> np.ndarray:
    """Per read: does it carry an exon-internal deletion that is not a junction?"""
    out = np.zeros(reads.n, bool)
    if reads.gap_start.size == 0:
        return out
    counts = np.diff(reads.goff)
    out[counts > 0] = True
    if not promoted:
        return out
    # a gap that became a junction is now explained; drop those reads again
    jset = set(zip(reads.donor.tolist(), reads.acceptor.tolist()))
    for i in np.flatnonzero(counts > 0):
        a, b = int(reads.goff[i]), int(reads.goff[i + 1])
        if all(
            (int(x), int(y)) in jset
            for x, y in zip(reads.gap_start[a:b], reads.gap_end[a:b])
        ):
            out[i] = False
    return out


# ---------------------------------------------------------------------- #
# snapping
# ---------------------------------------------------------------------- #
def build_snap_maps(
    census: JunctionCensus,
    contig: str,
    ref: ReferenceIndex,
    params: JunctionParams,
) -> Tuple[Dict[Tuple[str, int], int], Dict[Tuple[str, int], int]]:
    """Return ``{(strand, raw_pos): canonical_pos}`` for donors and acceptors.

    Order of preference, per site, taking sites most-supported first:

    1. an annotated splice site within ``fuzzy_tolerance``  (annotation wins)
    2. an already-accepted seed within ``fuzzy_tolerance``   (dominant read site)
    3. itself
    """
    tol = params.fuzzy_tolerance
    maps: List[Dict[Tuple[str, int], int]] = [{}, {}]
    if tol < 0 or census.size == 0:
        return maps[0], maps[1]

    schars = census.strand_chars()
    for which, positions in (("donor", census.donor), ("acceptor", census.acceptor)):
        support: Dict[Tuple[str, int], int] = defaultdict(int)
        for s, p, n in zip(schars, positions, census.n_reads):
            support[(str(s), int(p))] += int(n)
        seeds: Dict[str, List[int]] = defaultdict(list)
        out = maps[0] if which == "donor" else maps[1]
        for (s, p), _n in sorted(support.items(), key=lambda kv: (-kv[1], kv[0])):
            key = (contig, s)
            target = None
            if params.snap_to_annotation_first and tol >= 0:
                target = ref.snap_site(key, p, tol, which)
            if target is None and tol > 0:
                arr = seeds[s]
                i = int(np.searchsorted(arr, p))
                for j in (i - 1, i):
                    if 0 <= j < len(arr) and abs(arr[j] - p) <= tol:
                        target = arr[j]
                        break
            if target is None:
                target = p
                arr = seeds[s]
                arr.insert(int(np.searchsorted(arr, p)), p)
            out[(s, p)] = target
    return maps[0], maps[1]


# ---------------------------------------------------------------------- #
# features
# ---------------------------------------------------------------------- #
def junction_features(
    census: JunctionCensus,
    contig: str,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    params: JunctionParams,
    sj: Optional["object"] = None,
) -> pd.DataFrame:
    """One row per unique junction, with everything the novelty model needs."""
    n = census.size
    schars = census.strand_chars()
    donor = census.donor
    acceptor = census.acceptor

    annotated = np.zeros(n, bool)
    d_ann = np.zeros(n, bool)
    a_ann = np.zeros(n, bool)
    d_dist = np.zeros(n, np.int32)
    a_dist = np.zeros(n, np.int32)
    for i in range(n):
        key = (contig, str(schars[i]))
        j = (int(donor[i]), int(acceptor[i]))
        annotated[i] = ref.is_annotated_junction(key, j)
        d_ann[i] = ref.is_annotated_donor(key, j[0])
        a_ann[i] = ref.is_annotated_acceptor(key, j[1])
        d_dist[i] = min(ref.site_distance(key, j[0], "donor"), 1 << 20)
        a_dist[i] = min(ref.site_distance(key, j[1], "acceptor"), 1 << 20)

    motif = np.full(n, "NN-NN", object)
    motif_cls = np.full(n, "unknown", object)
    repeat = np.zeros(n, np.int16)
    if genome is not None:
        for i in range(n):
            m = genome.splice_motif(contig, int(donor[i]), int(acceptor[i]), str(schars[i]))
            motif[i] = m
            motif_cls[i] = genome.motif_class(m)
            repeat[i] = genome.direct_repeat_len(contig, int(donor[i]), int(acceptor[i]))

    # site usage shares: how dominant is this junction at its own donor/acceptor?
    d_tot = defaultdict(int)
    a_tot = defaultdict(int)
    for i in range(n):
        d_tot[(str(schars[i]), int(donor[i]))] += int(census.n_reads[i])
        a_tot[(str(schars[i]), int(acceptor[i]))] += int(census.n_reads[i])
    d_share = np.array(
        [census.n_reads[i] / max(d_tot[(str(schars[i]), int(donor[i]))], 1) for i in range(n)]
    )
    a_share = np.array(
        [census.n_reads[i] / max(a_tot[(str(schars[i]), int(acceptor[i]))], 1) for i in range(n)]
    )

    df = pd.DataFrame(
        {
            "contig": contig,
            "strand": schars,
            "donor": donor,
            "acceptor": acceptor,
            "intron_len": acceptor - donor,
            "n_reads": census.n_reads,
            "n_mols": census.n_mols,
            "n_cells": census.n_cells,
            "annotated": annotated,
            "donor_annotated": d_ann,
            "acceptor_annotated": a_ann,
            "n_sites_annotated": d_ann.astype(int) + a_ann.astype(int),
            "dist_ann_donor": d_dist,
            "dist_ann_acceptor": a_dist,
            "motif": motif,
            "motif_class": motif_cls,
            "direct_repeat": repeat,
            "donor_share": d_share,
            "acceptor_share": a_share,
            "min_site_share": np.minimum(d_share, a_share),
        }
    )
    df["support"] = np.where(df["n_mols"] > 0, df["n_mols"], df["n_reads"])
    _add_short_read_support(df, contig, schars, donor, acceptor, params, sj)
    return df


def _add_short_read_support(
    df: pd.DataFrame,
    contig: str,
    schars: np.ndarray,
    donor: np.ndarray,
    acceptor: np.ndarray,
    params: JunctionParams,
    sj: Optional["object"],
) -> None:
    """Attach STAR support, the supported flag, and the anchor/holdout split.

    ``sr_anchor`` and ``sr_holdout`` partition the short-read-supported *novel*
    junctions.  The anchor half joins the positive training set; the holdout
    half deliberately does not, so that the local FDR it receives is a
    prediction rather than a memory.  The split is a hash of the junction's own
    coordinates, so it is identical in every sample and stays valid when two
    samples are compared.
    """
    n = len(df)
    df["sr_rescue_enabled"] = bool(params.short_read_rescue)
    if sj is None:
        for c in ("sr_uniq", "sr_multi", "sr_overhang", "sr_n_libs"):
            df[c] = np.zeros(n, np.int64)
        for c in ("sr_in_sjdb", "sr_found", "sr_supported", "sr_anchor", "sr_holdout"):
            df[c] = np.zeros(n, bool)
        return
    arrays = sj.support_arrays(contig, schars, donor, acceptor)
    for k, v in arrays.items():
        df[k] = v
    # clamped, so supplying one library does not silently support nothing
    need_libs = max(1, min(params.min_sr_replicates, sj.n_libraries))
    supported = (
        (df["sr_uniq"].to_numpy() >= params.min_sr_unique_reads)
        & (df["sr_overhang"].to_numpy() >= params.min_sr_overhang)
        & (df["sr_n_libs"].to_numpy() >= need_libs)
    )
    df["sr_supported"] = supported
    novel_supported = supported & ~df["annotated"].to_numpy()
    hold = sj.holdout_mask(schars, donor, acceptor, params.sr_holdout_frac)
    df["sr_holdout"] = novel_supported & hold
    df["sr_anchor"] = novel_supported & ~hold


# ---------------------------------------------------------------------- #
# curation
# ---------------------------------------------------------------------- #
def curate(
    feat: pd.DataFrame,
    params: JunctionParams,
    lfdr: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Add ``keep`` and ``drop_reason`` columns.

    Annotated junctions are always kept -- the tool is reference-anchored, and
    dropping an annotated junction because it is rare in one sample would
    fabricate novel chains out of well-known transcripts.
    """
    n = len(feat)
    keep = np.ones(n, bool)
    reason = np.full(n, "", object)
    ann = feat["annotated"].to_numpy()

    def _reject(mask: np.ndarray, why: str) -> None:
        m = mask & keep & ~ann
        keep[m] = False
        reason[m] = why

    use_mols = bool((feat["n_mols"] > 0).any())
    if use_mols:
        _reject(feat["n_mols"].to_numpy() < params.min_novel_junction_umis, "low_molecules")
    _reject(feat["n_reads"].to_numpy() < params.min_novel_junction_reads, "low_reads")
    _reject(
        feat["min_site_share"].to_numpy() < params.min_novel_junction_share,
        "minor_at_both_sites",
    )
    if params.require_canonical_motif:
        cls = feat["motif_class"]
        ok = cls.isin(params.novel_motif_classes).to_numpy()
        # 0.1.16: a conditional class (GC-AG, AT-AC) is admitted on independent
        # evidence rather than on its motif. Short-read confirmation first;
        # failing that, a high long-read floor. Without either it is rejected
        # exactly as before, so the artefact population this rule was built for
        # is unaffected.
        cond = cls.isin(getattr(params, "conditional_motif_classes", ())).to_numpy()
        if cond.any():
            sr_reps = (feat["sr_n_libs"].to_numpy()
                       if "sr_n_libs" in feat.columns else np.zeros(n))
            rescued = cond & (
                (sr_reps >= params.conditional_motif_min_sr_replicates)
                | (feat["n_reads"].to_numpy() >= params.conditional_motif_min_reads)
            )
            ok = ok | rescued
        _reject(~ok, "motif_class_not_allowed_for_novel")
    if params.require_annotated_site_for_novel_junction and "n_sites_annotated" in feat:
        _reject(feat["n_sites_annotated"].to_numpy() < 1, "no_annotated_splice_site")
    if params.rt_switch_repeat_len > 0:
        _reject(
            feat["direct_repeat"].to_numpy() >= params.rt_switch_repeat_len,
            "rt_switch_repeat",
        )
    if lfdr is not None and params.max_junction_lfdr is not None:
        _reject(lfdr > params.max_junction_lfdr, "lfdr")

    # Rescue, opt-in.  Note the asymmetry, which is the whole point: short-read
    # support can *save* a junction, and its absence can never condemn one.
    # With 24% unique mapping and a 3'-biased long-read library, "the short
    # reads did not see it" is not evidence that it is not there.
    if params.short_read_rescue and "sr_supported" in feat.columns:
        rescued = feat["sr_supported"].to_numpy().astype(bool) & ~keep
        keep[rescued] = True
        reason[rescued] = "rescued_by_short_reads"

    out = feat.copy()
    out["keep"] = keep
    out["drop_reason"] = reason
    if lfdr is not None:
        out["lfdr"] = lfdr
    return out


def kept_junction_set(cur: pd.DataFrame) -> set:
    sub = cur.loc[cur["keep"]]
    return set(zip(sub["strand"], sub["donor"].astype(int), sub["acceptor"].astype(int)))


def read_junction_status(
    reads: ContigReads, census: JunctionCensus, keep_mask: np.ndarray
) -> np.ndarray:
    """Per read: True if *every* junction it carries survived curation."""
    ok = np.ones(reads.n, bool)
    if reads.donor.size == 0:
        return ok
    bad_rec = ~keep_mask[census.jid_of_record]
    if not bad_rec.any():
        return ok
    counts = np.diff(reads.joff)
    read_of_record = np.repeat(np.arange(reads.n), counts)
    ok[np.unique(read_of_record[bad_rec])] = False
    return ok
