"""Short-read coverage as terminal evidence (0.3.0, T6).

Built against a real bigWig rather than a mock, because the two things most
likely to be wrong here are both file-level: which way round the flanking
windows go on the minus strand, and whether an Ensembl-named track (`1`) is
found from a UCSC-named reference (`chr1`). Both fail silently -- the first
inverts the ratio for half the genome, the second returns "no coverage"
everywhere, which is indistinguishable from a locus with no expression.
"""

import numpy as np
import pytest

from flightcollapse.config import CoverageParams
from flightcollapse.coverage import (
    HAVE_PYBIGWIG,
    CoverageTrack,
    _flanks,
    blank_coverage_evidence,
)

pytestmark = pytest.mark.skipif(not HAVE_PYBIGWIG, reason="pyBigWig not installed")

P = CoverageParams()


@pytest.fixture(scope="module")
def bw(tmp_path_factory):
    """A track with a sharp drop at 5000 on contig `1` (Ensembl naming)."""
    import pyBigWig

    p = str(tmp_path_factory.mktemp("cov") / "t.bw")
    f = pyBigWig.open(p, "w")
    f.addHeader([("1", 20000)])
    starts = list(range(0, 20000, 100))
    values = [100.0 if s < 5000 else 2.0 for s in starts]
    f.addEntries(["1"] * len(starts), starts,
                 ends=[s + 100 for s in starts], values=values)
    f.close()
    return p


# ----------------------------------------------------------------- windows --
def test_flanks_are_transcript_oriented():
    assert _flanks(1000, "+", 100, True) == ((900, 1000), (1000, 1100))
    assert _flanks(1000, "-", 100, True) == ((1000, 1100), (900, 1000))


# ------------------------------------------------------------------ lookup --
def test_missing_track_is_no_track_not_an_error():
    assert CoverageTrack.open([]) is None
    assert CoverageTrack.open(["/nowhere/none.bw"], verbose=False) is None


def test_contig_naming_is_normalised_both_ways(bw):
    t = CoverageTrack.open([bw], verbose=False)
    assert t.mean("1", 0, 1000) == pytest.approx(100.0)
    assert t.mean("chr1", 0, 1000) == pytest.approx(100.0)
    assert t.mean("chrZ", 0, 1000) is None
    t.close()


def test_three_prime_step_finds_a_real_cleavage_site(bw):
    t = CoverageTrack.open([bw], verbose=False)
    ev = t.three_prime_step("chr1", 5000, "+", P)
    assert ev["sr_cov_upstream"] == pytest.approx(100.0)
    assert ev["sr_cov_downstream"] == pytest.approx(2.0)
    assert ev["sr_3p_step_ratio"] == pytest.approx(0.02)
    assert ev["sr_3p_step_supported"] is True
    t.close()


def test_no_step_where_the_transcript_continues(bw):
    """An internal-priming site has coverage on both sides -- that is the whole
    point of measuring the step rather than the depth."""
    t = CoverageTrack.open([bw], verbose=False)
    ev = t.three_prime_step("chr1", 2000, "+", P)
    assert ev["sr_3p_step_ratio"] == pytest.approx(1.0)
    assert ev["sr_3p_step_supported"] is False
    t.close()


def test_the_step_is_strand_aware(bw):
    """On the minus strand the same coordinate is a RISE, not a drop.

    Uses a lower min_coverage than the default: on the minus strand the
    UPSTREAM window here is the low side, and the default floor would
    (correctly) withhold an opinion rather than report a ratio off 2x.
    """
    q = CoverageParams(min_coverage=1.0)
    t = CoverageTrack.open([bw], verbose=False)
    plus = t.three_prime_step("chr1", 5000, "+", q)
    minus = t.three_prime_step("chr1", 5000, "-", q)
    assert plus["sr_3p_step_supported"] is True
    assert minus["sr_3p_step_supported"] is False
    assert minus["sr_3p_step_ratio"] == pytest.approx(50.0)
    t.close()


def test_low_coverage_means_no_opinion(bw):
    """Below min_coverage the ratio is noise, so it is withheld rather than
    reported -- absent and negative are different statements."""
    t = CoverageTrack.open([bw], verbose=False)
    ev = t.three_prime_step("chr1", 10000, "+", P)
    assert ev["sr_3p_step_ratio"] is None
    assert ev["sr_3p_step_supported"] is None
    t.close()


def test_tss_ratio(bw):
    q = CoverageParams(min_coverage=1.0)
    t = CoverageTrack.open([bw], verbose=False)
    # on the minus strand, "downstream of the start" is the high side
    ev = t.tss_ratio("chr1", 5000, "-", q)
    assert ev["sr_tss_ratio"] == pytest.approx(50.0)
    assert ev["sr_tss_supported"] is True
    ev = t.tss_ratio("chr1", 5000, "+", q)
    assert ev["sr_tss_supported"] is False
    t.close()


def test_continuity_is_length_weighted(bw):
    t = CoverageTrack.open([bw], verbose=False)
    v = t.continuity("chr1", [(0, 1000), (6000, 7000)])
    assert v == pytest.approx(51.0)
    assert t.continuity("chr1", []) is None
    t.close()


def test_blank_evidence_keeps_every_key():
    """A run without coverage must still emit the columns, blank, so a table
    from a no-coverage run lines up with one from a coverage run."""
    b = blank_coverage_evidence()
    assert set(b) == {
        "sr_cov_upstream", "sr_cov_downstream", "sr_3p_step_ratio",
        "sr_3p_step_supported", "sr_tss_ratio", "sr_tss_supported",
        "sr_cov_continuity"}
    assert all(v == "" for v in b.values())
