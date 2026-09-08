"""0.1.16: the two gates that were losing real junctions.

Both rules were measured against junctions *known* to be real, by hiding 7,224
GENCODE transcripts from a BD144_HM run and asking which of their junctions the
caller then failed to emit. Of 230 curation rejections, 174 came from the share
rule and 44 from the motif rule -- 95% between them -- while the lFDR rejected
2 and the support floors 9.

The junctions the share rule rejected carried 51 reads and 49.5 cells at the
median, against 13 and 12 for the novel junctions it kept, at lfdr 0.000. Both
quantities the two rules cut on (share, motif) are already features of the
calibrated scorer, so the hard thresholds were overruling a model on its own
inputs.

These tests pin the new behaviour AND pin the escape hatch back to 0.1.15, so a
future reader can see exactly what changed and undo it in one argument.
"""

import pandas as pd

from flightcollapse.config import JunctionParams
from flightcollapse.junctions import curate


def _feat(**over):
    f = pd.DataFrame({
        "annotated":         [True,             False,            False],
        "motif_class":       ["semi_canonical", "semi_canonical", "canonical"],
        "motif":             ["GC-AG",          "GC-AG",          "GT-AG"],
        "n_sites_annotated": [2, 2, 2],
        "direct_repeat":     [0, 0, 0],
        "n_reads":           [100, 20, 20],
        "n_mols":            [50, 10, 10],
        "min_site_share":    [1.0, 1.0, 1.0],
    })
    for k, v in over.items():
        f[k] = v
    return f


# ---------------------------------------------------------------- motif ---- #
def test_semi_canonical_is_admitted_on_short_read_confirmation():
    """SOD1's junction: GC-AG, 1,005 unique short reads in 3 of 3 replicates,
    both sites annotated, 10,419 long reads on the chain -- rejected by
    0.1.15 for its motif alone."""
    f = _feat(sr_n_libs=[0, 3, 0])
    out = curate(f, JunctionParams())
    assert out.loc[1, "keep"]
    assert out.loc[1, "drop_reason"] == ""


def test_semi_canonical_is_admitted_on_a_high_long_read_floor():
    f = _feat(n_reads=[100, 80, 20])
    assert curate(f, JunctionParams()).loc[1, "keep"]


def test_semi_canonical_without_evidence_is_still_rejected():
    """The artefact population the rule was built for is untouched: novel GC-AG
    with no external support and few reads confirms at 0.13%."""
    out = curate(_feat(sr_n_libs=[0, 1, 0]), JunctionParams())
    assert not out.loc[1, "keep"]
    assert out.loc[1, "drop_reason"] == "motif_class_not_allowed_for_novel"


def test_one_short_read_replicate_is_not_enough():
    """Single-replicate short-read support is within sampling noise: 4,302
    single-replicate novel junctions were seen against 1,812 expected."""
    assert not curate(_feat(sr_n_libs=[0, 1, 0]), JunctionParams()).loc[1, "keep"]


def test_the_conditional_rule_can_be_switched_off_entirely():
    """`conditional_motif_classes=()` is exactly 0.1.15."""
    p = JunctionParams(conditional_motif_classes=())
    assert not curate(_feat(sr_n_libs=[0, 3, 0]), p).loc[1, "keep"]


def test_an_annotated_semi_canonical_junction_needs_no_evidence():
    assert curate(_feat(), JunctionParams()).loc[0, "keep"]


def test_a_missing_short_read_column_falls_back_to_the_read_floor():
    """No short reads supplied: the run must not crash, and the long-read floor
    is then the only route in."""
    f = _feat(n_reads=[100, 80, 20])
    assert "sr_n_libs" not in f.columns
    assert curate(f, JunctionParams()).loc[1, "keep"]
    assert not curate(_feat(), JunctionParams()).loc[1, "keep"]


# ---------------------------------------------------------------- share ---- #
def test_a_minority_isoform_junction_survives_the_new_share_floor():
    """0.007 is the median share of the real junctions 0.1.15 threw away. A
    minority isoform at a well-expressed gene has a low share by definition."""
    f = _feat(min_site_share=[1.0, 0.007, 0.007], sr_n_libs=[0, 3, 0])
    out = curate(f, JunctionParams())
    assert out.loc[2, "keep"], "novel GT-AG at 0.7% share must now survive"
    assert curate(f, JunctionParams(min_novel_junction_share=0.01)) \
        .loc[2, "drop_reason"] == "minor_at_both_sites"


def test_the_floor_still_catches_genuine_scatter():
    """Lowered, not removed: misalignment scatter around a dominant junction
    does sit near zero, and the lFDR is a population-level control that one
    pathological locus can slip past."""
    f = _feat(min_site_share=[1.0, 1.0, 0.0001])
    out = curate(f, JunctionParams())
    assert not out.loc[2, "keep"]
    assert out.loc[2, "drop_reason"] == "minor_at_both_sites"


def test_annotated_junctions_are_never_touched_by_either_rule():
    f = _feat(min_site_share=[0.0001, 1.0, 1.0])
    assert curate(f, JunctionParams()).loc[0, "keep"]
