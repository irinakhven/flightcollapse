"""STAR short-read junctions as an anchor.

The coordinate convention is the thing to guard hardest. A one-base error does
not raise: it returns "no short-read support" for every junction, which looks
exactly like an orthogonal library that confirms nothing. Verified against the
real BD144b table this agreement is 99.99%; these tests make sure a regression
that breaks it is loud.
"""

import gzip

import numpy as np
import pandas as pd
import pytest

from flightcollapse.config import JunctionParams, ScoringParams
from flightcollapse.junctions import _add_short_read_support, curate
from flightcollapse.reference import ReferenceIndex
from flightcollapse.scoring import junction_anchors
from flightcollapse.shortreads import SpliceJunctionCatalogue, norm_contig


#: intron 1201-2000 (1-based inclusive) == (1200, 2000) half-open here
STAR_ROWS = [
    # contig start end strand motif sjdb uniq multi overhang
    ("1", 1201, 2000, 1, 1, 1, 40, 3, 30),      # annotated, well supported
    ("1", 2201, 3000, 1, 1, 1, 25, 0, 28),      # annotated
    ("1", 4001, 4500, 1, 1, 0, 12, 1, 25),      # NOVEL, well supported
    ("1", 5001, 5500, 1, 0, 0, 1, 0, 6),        # NOVEL, one read, short overhang
]


def _write(tmp_path, rows, name="SJ.out.tab", gz=False):
    text = "\n".join("\t".join(str(x) for x in r) for r in rows) + "\n"
    p = tmp_path / (name + (".gz" if gz else ""))
    if gz:
        p.write_bytes(gzip.compress(text.encode()))
    else:
        p.write_text(text)
    return str(p)


def _gtf(tmp_path):
    """A transcript whose introns are exactly the two annotated STAR rows."""
    a = ('gene_id "G1"; transcript_id "ENST1"; gene_name "G1"; '
         'transcript_type "protein_coding"; tag "basic";')
    ex = [(1000, 1200), (2000, 2200), (3000, 3500)]
    rows = [f"chr1\tT\ttranscript\t1001\t3500\t.\t+\t.\t{a}"]
    rows += [f"chr1\tT\texon\t{s + 1}\t{e}\t.\t+\t.\t{a}" for s, e in ex]
    p = tmp_path / "g.gtf"
    p.write_text("\n".join(rows) + "\n")
    return ReferenceIndex.from_gtf(str(p), verbose=False)


# ---------------------------------------------------------------------- #
def test_star_coordinates_convert_to_half_open_donor_acceptor(tmp_path):
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    assert cat.lookup("chr1", "+", 1200, 2000)[0] == 40
    assert cat.lookup("chr1", "+", 1201, 2001) is None      # the off-by-one
    assert cat.lookup("chr1", "+", 1199, 1999) is None


def test_contig_naming_is_normalised_on_both_sides(tmp_path):
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    # the table says "1", the BAM says "chr1"; both must resolve
    assert cat.lookup("1", "+", 1200, 2000) is not None
    assert cat.lookup("chr1", "+", 1200, 2000) is not None
    assert norm_contig("chrM") == norm_contig("MT") == "MT"


def test_replicates_are_summed(tmp_path):
    a = _write(tmp_path, STAR_ROWS, "a.tab")
    b = _write(tmp_path, STAR_ROWS, "b.tab")
    cat = SpliceJunctionCatalogue.from_star([a, b], verbose=False)
    hit = cat.lookup("chr1", "+", 1200, 2000)
    assert hit[0] == 80 and hit[1] == 6          # reads summed
    assert hit[2] == 30                          # overhang is a max, not a sum
    assert len(cat) == len(STAR_ROWS)


def test_gzipped_input_works(tmp_path):
    cat = SpliceJunctionCatalogue.from_star(
        [_write(tmp_path, STAR_ROWS, "c.tab", gz=True)], verbose=False
    )
    assert len(cat) == len(STAR_ROWS)


# ---------------------------------------------------------------------- #
def test_verification_passes_when_coordinates_agree(tmp_path):
    ref = _gtf(tmp_path)
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    v = cat.verify_against(ref)
    assert v["sjdb_junctions"] == 2
    assert v["sjdb_found_in_reference"] == 2
    assert v["agreement"] == 1.0


def test_verification_names_an_off_by_one(tmp_path):
    ref = _gtf(tmp_path)
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    broken = SpliceJunctionCatalogue()
    broken.j = {(c, s, d + 1, a + 1): v for (c, s, d, a), v in cat.j.items()}
    v = broken.verify_against(ref)
    assert v["agreement"] == 0.0
    assert v["systematic_offset_hits"].get(-1) == 2
    assert "coordinate-convention mismatch" in v["diagnosis"]


def test_verification_names_a_contig_mismatch(tmp_path):
    ref = _gtf(tmp_path)
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    wrong = SpliceJunctionCatalogue()
    wrong.j = {("scaffold_" + c, s, d, a): v for (c, s, d, a), v in cat.j.items()}
    v = wrong.verify_against(ref)
    assert not v["ok"]
    assert "contigs do not line up" in v["diagnosis"]


# ---------------------------------------------------------------------- #
def _feat(tmp_path, **over):
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    df = pd.DataFrame({
        # 0,1 annotated | 2 novel+supported | 3 novel, seen but under the floor
        # 4 novel non-canonical, unseen (a decoy) | 5 novel canonical, unseen but
        # well supported in the long reads -- this one must survive
        "strand": ["+"] * 6,
        "donor": [1200, 2200, 4000, 5000, 9000, 12000],
        "acceptor": [2000, 3000, 4500, 5500, 9500, 12500],
        "annotated": [True, True, False, False, False, False],
        "motif_class": ["canonical"] * 3 + ["non_canonical", "non_canonical",
                                            "canonical"],
        "direct_repeat": [0] * 6,
        "n_reads": [100, 100, 30, 4, 4, 30],
        "n_mols": [50, 50, 15, 2, 2, 15],
        "min_site_share": [1.0] * 6,
    })
    params = JunctionParams(**over)
    _add_short_read_support(
        df, "chr1", df.strand.to_numpy(), df.donor.to_numpy(),
        df.acceptor.to_numpy(), params, cat,
    )
    return df, params


def test_support_flags_respect_the_overhang_floor(tmp_path):
    df, _ = _feat(tmp_path)
    # row 2 is the well-supported novel junction; row 3 has 1 read / 6 bp overhang
    assert list(df.sr_supported) == [True, True, True, False, False, False]
    assert list(df.sr_found) == [True, True, True, True, False, False]


def test_supported_novel_junctions_split_into_anchor_and_holdout(tmp_path):
    df, _ = _feat(tmp_path)
    novel_supported = (~df.annotated) & df.sr_supported
    # every supported novel junction is in exactly one of the two halves
    assert ((df.sr_anchor | df.sr_holdout) == novel_supported).all()
    assert not (df.sr_anchor & df.sr_holdout).any()
    # annotated junctions are in neither -- they are already positives
    assert not df.loc[df.annotated, ["sr_anchor", "sr_holdout"]].to_numpy().any()


def test_the_holdout_split_is_deterministic_and_position_derived(tmp_path):
    a, _ = _feat(tmp_path)
    b, _ = _feat(tmp_path)
    assert list(a.sr_holdout) == list(b.sr_holdout)
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, STAR_ROWS)], verbose=False)
    # ...and depends only on the junction, not on the order it was seen in
    fwd = cat.holdout_mask(["+", "-"], [1200, 4000], [2000, 4500], 0.5)
    rev = cat.holdout_mask(["-", "+"], [4000, 1200], [4500, 2000], 0.5)
    assert list(fwd) == list(rev)[::-1]


def test_short_reads_add_positives_and_shrink_the_decoy_set(tmp_path):
    df, _ = _feat(tmp_path)
    pos, decoy, exempt = junction_anchors(df)
    # a novel junction the short reads confirm is a true positive that
    # annotation could never have supplied
    assert pos.sum() >= df.annotated.sum()
    # ...and it is not exempt: it still earns its local FDR from long reads
    assert list(exempt) == list(df.annotated)
    # the unsupported non-canonical junction is a decoy; a supported one is not
    assert decoy[4]
    supported_noncanon = df.sr_supported.to_numpy() & (df.motif_class == "non_canonical")
    assert not decoy[supported_noncanon].any()


def test_rescue_is_opt_in_and_absence_never_drops(tmp_path):
    df, _ = _feat(tmp_path)
    df["lfdr"] = 0.0

    # by default a junction failing curation stays failed even if short reads saw it
    off = curate(df, JunctionParams(min_novel_junction_reads=50))
    assert not off.loc[2, "keep"]

    on = curate(df, JunctionParams(min_novel_junction_reads=50, short_read_rescue=True))
    assert on.loc[2, "keep"]
    assert on.loc[2, "drop_reason"] == "rescued_by_short_reads"

    # a junction with NO short-read support is never dropped for that alone
    plain = curate(df, JunctionParams())
    assert plain.loc[5, "keep"]
    assert not plain.sr_found[5]
    # nothing in the drop reasons mentions short reads, in either direction
    assert not any("short" in r or "sr_" in r for r in plain.drop_reason if r)


def test_features_are_inert_when_no_short_reads_are_supplied():
    df = pd.DataFrame({
        "strand": ["+"], "donor": [1200], "acceptor": [2000],
        "annotated": [False], "motif_class": ["canonical"], "direct_repeat": [0],
        "n_reads": [10], "n_mols": [5],
    })
    _add_short_read_support(
        df, "chr1", df.strand.to_numpy(), df.donor.to_numpy(),
        df.acceptor.to_numpy(), JunctionParams(), None,
    )
    assert not df.sr_supported.any() and not df.sr_anchor.any()
    pos, _decoy, exempt = junction_anchors(df)
    assert list(pos) == list(exempt) == [False]


# ---------------------------------------------------------------------- #
# replicate concordance
# ---------------------------------------------------------------------- #
#: the same novel junction, but reaching 2 reads two different ways
ROWS_A = [("1", 4001, 4500, 1, 1, 0, 2, 0, 25),    # both reads in library A
          ("1", 7001, 7500, 1, 1, 0, 1, 0, 25)]    # one read here...
ROWS_B = [("1", 7001, 7500, 1, 1, 0, 1, 0, 25)]    # ...and one in library B


def test_a_junction_confined_to_one_replicate_is_not_supported(tmp_path):
    """Both are at 2 summed unique reads; only one of them reproduced.

    On the real BD144 a/b/c data 4,302 novel junctions reach 2 summed reads in a
    single replicate where sampling predicts 1,812 -- so most of that class is
    library-specific, and an anchor set wants precision far more than recall.
    """
    cat = SpliceJunctionCatalogue.from_star(
        [_write(tmp_path, ROWS_A, "a.tab"), _write(tmp_path, ROWS_B, "b.tab")],
        verbose=False,
    )
    assert cat.lookup("chr1", "+", 4000, 4500)[0] == 2
    assert cat.lookup("chr1", "+", 7000, 7500)[0] == 2
    assert cat.n_libs(("1", "+", 4000, 4500)) == 1
    assert cat.n_libs(("1", "+", 7000, 7500)) == 2

    df = pd.DataFrame({
        "strand": ["+", "+"], "donor": [4000, 7000], "acceptor": [4500, 7500],
        "annotated": [False, False], "motif_class": ["canonical"] * 2,
        "direct_repeat": [0, 0], "n_reads": [10, 10], "n_mols": [5, 5],
        "min_site_share": [1.0, 1.0],
    })
    _add_short_read_support(
        df, "chr1", df.strand.to_numpy(), df.donor.to_numpy(),
        df.acceptor.to_numpy(), JunctionParams(), cat,
    )
    assert list(df.sr_supported) == [False, True]


def test_the_replicate_floor_is_clamped_to_the_libraries_supplied(tmp_path):
    """One SJ file must not silently support nothing."""
    cat = SpliceJunctionCatalogue.from_star([_write(tmp_path, ROWS_A)], verbose=False)
    df = pd.DataFrame({
        "strand": ["+"], "donor": [4000], "acceptor": [4500],
        "annotated": [False], "motif_class": ["canonical"], "direct_repeat": [0],
        "n_reads": [10], "n_mols": [5], "min_site_share": [1.0],
    })
    _add_short_read_support(
        df, "chr1", df.strand.to_numpy(), df.donor.to_numpy(),
        df.acceptor.to_numpy(), JunctionParams(min_sr_replicates=3), cat,
    )
    assert bool(df.sr_supported[0])


def test_concordance_compares_against_the_sampling_expectation(tmp_path):
    cat = SpliceJunctionCatalogue.from_star(
        [_write(tmp_path, ROWS_A, "a.tab"), _write(tmp_path, ROWS_B, "b.tab")],
        verbose=False,
    )
    c = cat.concordance(min_unique=2)
    assert c["n_libraries"] == 2
    assert c["by_n_replicates"] == {1: 1, 2: 1}
    assert c["seen_in_one_replicate"] == 1
    # a single library reports nothing rather than a meaningless number
    solo = SpliceJunctionCatalogue.from_star([_write(tmp_path, ROWS_A)], verbose=False)
    assert solo.concordance() == {}
