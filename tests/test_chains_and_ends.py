import numpy as np
import pytest

from flightcollapse.chains import ChainGroup, find_parents, resolve_suffixes
from flightcollapse.config import ChainParams
from flightcollapse.ends import cluster_3p_ends, five_prime_end


def _g(gid, chain, n_reads, annotated=False, tss_evidence=False):
    g = ChainGroup(
        gid=gid, contig="chr1", strand="+", chain=chain,
        reads=np.arange(n_reads), n_reads=n_reads, n_mols=n_reads, n_cells=n_reads,
    )
    if annotated:
        g.tx_ids = [f"ENST{gid}"]
    g.tss_evidence = tss_evidence
    return g


A = (10, 20)
B = (30, 40)
C = (50, 60)
FULL = (A, B, C)      # transcript orientation, C is 3'-most
TRUNC = (B, C)
TRUNC2 = (C,)


def test_parents_are_reduced_to_maximal():
    groups = [_g(0, FULL, 100), _g(1, TRUNC, 50), _g(2, TRUNC2, 10)]
    cand = np.ones(3, bool)
    parents = find_parents(groups, cand)
    assert parents[1] == [0]
    # TRUNC2 is a suffix of both FULL and TRUNC, but only FULL is maximal
    assert parents[2] == [0]


def test_truncation_merges_into_annotated_parent_even_when_it_dominates():
    """5' degradation is pervasive, so the truncated form routinely outnumbers
    the intact one.  Support direction alone must not rescue it."""
    groups = [_g(0, FULL, 20, annotated=True), _g(1, TRUNC, 500)]
    cand = np.ones(2, bool)
    resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into == 0
    assert groups[1].merge_class == "truncation_of_annotated_parent"


def test_alt_tss_evidence_rescues_a_suffix():
    groups = [_g(0, FULL, 200, annotated=True), _g(1, TRUNC, 300, tss_evidence=True)]
    cand = np.ones(2, bool)
    resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into is None
    assert groups[1].merge_class == "alt_tss_kept"


def test_annotated_chain_is_never_merged_away():
    """The short transcript really exists in GENCODE; a rare 5'-extended
    variant must not swallow it."""
    groups = [_g(0, FULL, 5), _g(1, TRUNC, 900, annotated=True)]
    cand = np.ones(2, bool)
    resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into is None
    assert groups[1].merge_class == "annotated_kept"


def test_dominant_unannotated_suffix_is_protected():
    groups = [_g(0, FULL, 5), _g(1, TRUNC, 900)]
    cand = np.ones(2, bool)
    stats = resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into is None
    assert stats["n_dominant_protected"] == 1


def test_debris_merges_upward_when_neither_is_annotated():
    groups = [_g(0, FULL, 900), _g(1, TRUNC, 10)]
    cand = np.ones(2, bool)
    resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into == 0
    assert groups[1].merge_class == "debris"


def test_merges_are_transitive_to_the_maximal_chain():
    groups = [_g(0, FULL, 900, annotated=True), _g(1, TRUNC, 40), _g(2, TRUNC2, 5)]
    cand = np.ones(3, bool)
    resolve_suffixes(groups, cand, ChainParams())
    assert groups[1].merged_into == 0
    assert groups[2].merged_into == 0     # never left pointing at the middle chain


# ---------------------------------------------------------------------- #
def test_cluster_3p_ends_finds_separate_peaks():
    pos = np.concatenate([
        np.full(100, 1000) + np.random.default_rng(0).integers(-20, 20, 100),
        np.full(60, 5000) + np.random.default_rng(1).integers(-20, 20, 60),
    ])
    labels, peaks = cluster_3p_ends(pos, window=100)
    assert peaks.size == 2
    assert abs(int(peaks[0]) - 1000) < 50
    assert abs(int(peaks[1]) - 5000) < 50
    assert set(np.bincount(labels)) == {100, 60}


def test_cluster_does_not_chain_a_whole_utr_into_one_peak():
    pos = np.arange(0, 3000, 10)          # a uniform smear, no real peak
    _labels, peaks = cluster_3p_ends(pos, window=50)
    assert peaks.size > 1


def test_five_prime_end_is_not_the_shortest_member():
    """The single line that produced the isoseq artifact."""
    tss = np.array([1000, 1010, 1020, 4800, 4900])   # + strand, 3' end at 5000
    five = five_prime_end(tss, 5000, "+", percentile=90.0)
    assert five < 2000            # near the long reads, not the truncated ones
    assert five != int(tss.max())  # not the most 5'-truncated read


def test_five_prime_end_minus_strand():
    tss = np.array([9000, 8990, 8980, 5200, 5100])   # - strand, 3' end at 5000
    five = five_prime_end(tss, 5000, "-", percentile=90.0)
    assert five > 8000
