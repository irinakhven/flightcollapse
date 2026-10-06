"""TSSAtlas loader and lookup mechanics.

The failure this file exists to prevent is the silent one: a coordinate
convention slip or a contig-naming mismatch does not raise, it returns "no
evidence" for every 5' end, which is indistinguishable from an atlas that
genuinely supports nothing.
"""

import gzip
import os

import pytest

from flightcollapse.genome import TSSAtlas


def _write(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


# BED9, FANTOM5 shape: chrom start end name score strand thickStart thickEnd rgb
BED9 = [
    ("chr1", 1000, 1060, "chr1:1000..1060,+", 42, "+", 1020, 1021, "0,0,0"),
    ("chr1", 5000, 5040, "chr1:5000..5040,-", 7, "-", 5030, 5031, "0,0,0"),
    ("chr2", 200, 260, "chr2:200..260,+", 99, "+", 210, 211, "0,0,0"),
]


def test_bed9_reads_the_representative_tss_from_the_thick_fields(tmp_path):
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    assert len(a) == 3
    # the repTSS is 1020, not the peak start
    assert a.dist_to_reptss("chr1", 1020, "+") == 0
    assert a.dist_to_reptss("chr1", 1000, "+") == 20


def test_containment_is_half_open(tmp_path):
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    assert a.in_peak("chr1", 1000, "+")        # start is inside
    assert a.in_peak("chr1", 1059, "+")
    assert not a.in_peak("chr1", 1060, "+")    # end is not
    assert not a.in_peak("chr1", 999, "+")


def test_strand_matters(tmp_path):
    """A peak on the other strand is a different promoter, not this one."""
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    assert a.in_peak("chr1", 1020, "+")
    assert not a.in_peak("chr1", 1020, "-")


def test_distance_to_a_peak_outside_it(tmp_path):
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    ev = a.evidence("chr1", 900, "+")
    assert ev["in_cage_peak"] is False
    assert ev["dist_to_cage_peak"] == 100
    ev = a.evidence("chr1", 1200, "+")
    assert ev["dist_to_cage_peak"] == 141   # 1200 - 1060 + 1


def test_unknown_contig_and_strand_return_blanks_not_errors(tmp_path):
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    ev = a.evidence("chrZ", 1000, "+")
    assert ev["in_cage_peak"] is False
    assert ev["dist_to_cage_peak"] is None
    assert ev["cage_peak_id"] == ""


def test_contig_prefix_is_normalised_both_ways(tmp_path):
    rows = [("1", 1000, 1060, "p1", 5, "+", 1020, 1021, "0,0,0")]
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), rows), verbose=False)
    assert a.in_peak("chr1", 1020, "+")
    assert a.in_peak("1", 1020, "+")


def test_unsorted_input_is_sorted_on_load(tmp_path):
    """CrossMap does not preserve input order."""
    rows = list(reversed(BED9))
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), rows), verbose=False)
    assert a.in_peak("chr1", 1020, "+")
    assert a.in_peak("chr2", 210, "+")


def test_gzip_is_transparent(tmp_path):
    p = str(tmp_path / "p.bed.gz")
    with gzip.open(p, "wt") as fh:
        for r in BED9:
            fh.write("\t".join(str(x) for x in r) + "\n")
    assert len(TSSAtlas.from_bed(p, verbose=False)) == 3


def test_track_and_comment_lines_are_skipped(tmp_path):
    p = str(tmp_path / "p.bed")
    with open(p, "w") as fh:
        fh.write("track name=cage\n# a comment\n")
        for r in BED9:
            fh.write("\t".join(str(x) for x in r) + "\n")
    assert len(TSSAtlas.from_bed(p, verbose=False)) == 3


def test_min_score_filters(tmp_path):
    p = _write(str(tmp_path / "p.bed"), BED9)
    assert len(TSSAtlas.from_bed(p, min_score=50, verbose=False)) == 1


def test_bed6_falls_back_to_the_strand_aware_five_prime_edge(tmp_path):
    rows = [
        ("chr1", 1000, 1060, "p+", 1, "+"),
        ("chr1", 5000, 5040, "p-", 1, "-"),
    ]
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), rows), verbose=False)
    assert a.dist_to_reptss("chr1", 1000, "+") == 0     # + : peak start
    assert a.dist_to_reptss("chr1", 5039, "-") == 0     # - : peak end - 1


def test_reptss_side_table_joins_on_peak_name(tmp_path):
    """The lifted-mouse case: peaks and repTSS arrive as two BED6 files."""
    peaks = [("chr1", 1000, 1060, "pA", 1, "+"), ("chr1", 5000, 5040, "pB", 1, "-")]
    reps = [("chr1", 1033, 1034, "pA", 1, "+"), ("chr1", 5011, 5012, "pB", 1, "-")]
    a = TSSAtlas.from_bed(
        _write(str(tmp_path / "p.bed"), peaks),
        _write(str(tmp_path / "r.bed"), reps),
        verbose=False,
    )
    assert a.n_reptss_joined == 2
    assert a.dist_to_reptss("chr1", 1033, "+") == 0
    assert a.dist_to_reptss("chr1", 5011, "-") == 0


def test_split_peaks_keep_both_intervals_and_join_the_right_one(tmp_path):
    """liftOver split 6 of the mouse peaks in two, so names are not unique."""
    peaks = [
        ("chr1", 1000, 1030, "pSplit", 1, "+"),
        ("chr1", 9000, 9030, "pSplit", 1, "+"),
    ]
    reps = [("chr1", 9010, 9011, "pSplit", 1, "+")]
    a = TSSAtlas.from_bed(
        _write(str(tmp_path / "p.bed"), peaks),
        _write(str(tmp_path / "r.bed"), reps),
        verbose=False,
    )
    assert len(a) == 2
    assert a.n_reptss_joined == 1
    # the joined repTSS attached to the interval that contains it ...
    assert a.dist_to_reptss("chr1", 9010, "+") == 0
    # ... and the other interval fell back to its own 5' edge rather than
    # borrowing a coordinate 8 kb away
    assert a.evidence("chr1", 1000, "+")["dist_to_cage_reptss"] == 0


def test_overlapping_peaks_are_both_findable(tmp_path):
    rows = [
        ("chr1", 1000, 1200, "wide", 1, "+"),
        ("chr1", 1050, 1080, "narrow", 1, "+"),
    ]
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), rows), verbose=False)
    assert a.in_peak("chr1", 1060, "+")
    assert a.in_peak("chr1", 1190, "+")


def test_peak_slack_widens_containment(tmp_path):
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    assert not a.in_peak("chr1", 1075, "+")
    assert a.in_peak("chr1", 1075, "+", slack=20)


def test_flip_strands_moves_every_peak(tmp_path):
    """The negative control used to estimate the rescue false-positive rate."""
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), BED9), verbose=False)
    b = a.flip_strands()
    assert a.in_peak("chr1", 1020, "+") and not a.in_peak("chr1", 1020, "-")
    assert b.in_peak("chr1", 1020, "-") and not b.in_peak("chr1", 1020, "+")
    assert len(b) == len(a)


def test_coverage_check_flags_a_wrong_assembly(tmp_path, sim):
    """An atlas whose peaks land nowhere near the annotation must say so."""
    from flightcollapse.reference import ReferenceIndex

    ref = ReferenceIndex.from_gtf(sim.gtf, verbose=False)
    rows = [("chrNOWHERE", 10, 70, "p", 1, "+", 20, 21, "0,0,0")]
    a = TSSAtlas.from_bed(_write(str(tmp_path / "p.bed"), rows), verbose=False)
    rep = a.coverage_of_annotated_tss(ref)
    assert rep["n_genes_tested"] > 0
    assert rep["frac_in_peak"] == 0.0
    assert rep["ok"] is False
    assert "GRCm39" in rep["diagnosis"]
