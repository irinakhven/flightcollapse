import math

import pytest

from flightcollapse.genome import Genome, POLYA_MOTIFS, revcomp


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    import pysam

    root = tmp_path_factory.mktemp("g")
    fa = root / "t.fa"
    #            0         1         2         3         4
    #            0123456789012345678901234567890123456789012345
    seq = "CCCCGTCCCCCCCCCCCCCCCCAGCCCCAATAAACCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAA"
    with open(fa, "w") as fh:
        fh.write(">chrT\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(fa))
    return Genome(str(fa))


def test_splice_motif_plus_strand(tiny):
    # intron starts at 4 ("GT") and the two bases before 24 are "AG"
    assert tiny.splice_motif("chrT", 4, 24, "+") == "GT-AG"


def test_splice_motif_minus_strand_is_reverse_complemented(tiny):
    m = tiny.splice_motif("chrT", 4, 24, "-")
    assert m == "CT-AC"
    assert tiny.motif_class(m) == "antisense"


def test_motif_classes(tiny):
    assert tiny.motif_class("GT-AG") == "canonical"
    assert tiny.motif_class("GC-AG") == "semi_canonical"
    assert tiny.motif_class("AT-AC") == "semi_canonical"
    assert tiny.motif_class("TT-TT") == "non_canonical"


def test_direct_repeat(tiny):
    # identical flanks -> non-zero repeat; the exact length depends on sequence
    r = tiny.direct_repeat_len("chrT", 10, 20)
    assert r >= 1


def test_polya_motif_found_upstream(tiny):
    found, motif, dist = tiny.polya_motif("chrT", 55, "+", window=50)
    assert found and motif == "AATAAA"
    assert 10 <= dist <= 50


def test_perc_a_is_on_a_0_100_scale(tiny):
    """The 0-100 trap: comparing this to 0.6 flags everything as internal priming."""
    pa = tiny.perc_a_downstream("chrT", 54, "+", window=20)
    assert pa == pytest.approx(100.0)
    assert pa > 1.5   # i.e. clearly not a fraction


def test_revcomp():
    assert revcomp("AATAAA") == "TTTATT"


def test_polya_motif_list_starts_with_the_two_dominant_signals():
    assert POLYA_MOTIFS[0] == "AATAAA"
    assert POLYA_MOTIFS[1] == "ATTAAA"
