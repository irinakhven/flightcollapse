"""The 3' clip must not be able to reject a read for being short.

Through 0.4.0 coverage was aligned-query over FULL query length, so the polyA
tail plus barcode plus UMI sat in the denominator. That clip is near-constant in
absolute terms while aligned length spans an order of magnitude, which turned
`min_aln_coverage` into a length filter: on BD176c it rejected 8,406,872 of
24,254,151 scanned reads, 88% of them in the 500-800 bp band, and 85% of what it
rejected was multi-exon.

These tests pin the arithmetic, not the biology: a perfectly aligned read with a
normal adapter passes at 0.90 regardless of its length, and a read that is
genuinely half unaligned at its 5' end still fails.
"""
from __future__ import annotations

import pytest

from flightcollapse.config import AlignmentGates
from flightcollapse.reads import alignment_stats

CMATCH, CREF_SKIP, CSOFT = 0, 3, 4


class FakeRead:
    """The four things alignment_stats actually reads."""

    def __init__(self, cigar, is_reverse=False, nm=0):
        self.cigartuples = list(cigar)
        self.is_reverse = is_reverse
        self._nm = nm

    def infer_read_length(self):
        return sum(ln for op, ln in self.cigartuples if op != CREF_SKIP)

    @property
    def query_length(self):
        return self.infer_read_length()

    def get_tag(self, tag):
        if tag == "NM":
            return self._nm
        raise KeyError(tag)


def read_with(aligned: int, clip3: int, clip5: int = 0, is_reverse: bool = False):
    """One read: clip5 | aligned | clip3, in TRANSCRIPT orientation."""
    cig = []
    if clip5:
        cig.append((CSOFT, clip5))
    cig.append((CMATCH, aligned))
    if clip3:
        cig.append((CSOFT, clip3))
    if is_reverse:
        cig = list(reversed(cig))
    return FakeRead(cig, is_reverse=is_reverse)


# The measured BD144 3' clip: tail + 27 nt barcode + 8 nt UMI.
CLIP = 99


@pytest.mark.parametrize("aligned", [200, 400, 600, 800, 1200, 2000, 4000])
@pytest.mark.parametrize("is_reverse", [False, True])
def test_adapter_never_rejects_a_clean_alignment(aligned, is_reverse):
    """Every length passes 0.90 once the clip leaves the denominator."""
    r = read_with(aligned, CLIP, is_reverse=is_reverse)
    cov, _, _, clip3 = alignment_stats(r, exclude_3p_clip=True)
    assert clip3 == CLIP
    assert cov == pytest.approx(1.0)
    assert cov >= AlignmentGates().min_aln_coverage


@pytest.mark.parametrize(
    "aligned,passes_old",
    [(200, False), (400, False), (600, False), (891, True), (2000, True)],
)
def test_old_statistic_was_a_length_filter(aligned, passes_old):
    """0.4.0's arithmetic, pinned: the break-even is clip*t/(1-t) = 891 bp."""
    r = read_with(aligned, CLIP)
    cov, _, _, _ = alignment_stats(r, exclude_3p_clip=False)
    assert (cov >= 0.90) is passes_old


def test_the_two_disagree_exactly_where_predicted():
    """A 600 bp read with a median adapter: rejected before, kept now."""
    r = read_with(600, CLIP)
    old, _, _, _ = alignment_stats(r, exclude_3p_clip=False)
    new, _, _, _ = alignment_stats(r, exclude_3p_clip=True)
    assert old < 0.90 < new
    assert old == pytest.approx(600 / 699, abs=1e-6)
    assert new == pytest.approx(1.0)


def test_five_prime_clip_still_counts():
    """The 5' clip is unexplained sequence and keeps its say in coverage.

    Excluding the 3' clip must not become 'ignore every clip': a read whose 5'
    half did not align is a real alignment problem, and max_5p_softclip gates
    the length separately but only coverage sees the ratio.
    """
    r = read_with(500, CLIP, clip5=500)
    cov, _, clip5, _ = alignment_stats(r, exclude_3p_clip=True)
    assert clip5 == 500
    assert cov == pytest.approx(0.5)
    assert cov < AlignmentGates().min_aln_coverage


def test_reverse_strand_picks_the_right_end():
    """clip3 is transcript-3', so on a reverse read it is the LEFT clip.

    Getting this backwards would exclude the 5' clip from the denominator and
    keep the 3' one -- the exact inverse of the fix, and silent.
    """
    fwd = read_with(500, 99, clip5=20, is_reverse=False)
    rev = read_with(500, 99, clip5=20, is_reverse=True)
    for r in (fwd, rev):
        cov, _, clip5, clip3 = alignment_stats(r, exclude_3p_clip=True)
        assert (clip5, clip3) == (20, 99)
        assert cov == pytest.approx(500 / 520)


def test_all_clip_read_does_not_divide_by_zero():
    r = FakeRead([(CSOFT, 150)])
    cov, _, _, _ = alignment_stats(r, exclude_3p_clip=True)
    assert cov == 0.0


def test_coverage_is_capped_at_one():
    """An insertion-heavy alignment must not report coverage above 1."""
    r = FakeRead([(CMATCH, 300), (1, 200), (CSOFT, 400)])  # 200 bp of I
    cov, _, _, _ = alignment_stats(r, exclude_3p_clip=True)
    assert cov <= 1.0


def test_gate_default_is_on_and_reversible():
    assert AlignmentGates().coverage_excludes_3p_clip is True
    assert AlignmentGates(coverage_excludes_3p_clip=False).coverage_excludes_3p_clip is False


def test_function_default_preserves_old_behaviour():
    """alignment_stats' own default stays 0.4.0, so any external caller that
    does not pass the flag is unchanged; only the gate opts in."""
    r = read_with(600, CLIP)
    assert alignment_stats(r)[0] == pytest.approx(600 / 699, abs=1e-6)
