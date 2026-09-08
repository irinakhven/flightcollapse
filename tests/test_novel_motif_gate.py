"""Novel junctions must be GT-AG, and semi-canonical ones are decoys.

Both rules come from short reads rather than from first principles. Of the novel
junctions the 0.1.6 caller kept on BD144, the GC-AG/AT-AC class was confirmed by
an orthogonal library at 0.13% while annotated semi-canonical junctions were
confirmed at 33.4% -- so the motif class is visible to short reads and the novel
members of it are simply not real. The same numbers reproduced in all three
libraries (0.13 / 0.12 / 0.15%).
"""

import numpy as np
import pandas as pd
import pytest

from flightcollapse.config import JunctionParams
from flightcollapse.junctions import curate
from flightcollapse.scoring import junction_anchors


def _feat():
    return pd.DataFrame({
        "annotated":         [True,             True,   False,       False,    False],
        "motif_class":       ["semi_canonical", "canonical",
                              "canonical",      "semi_canonical", "canonical"],
        "motif":             ["GC-AG", "GT-AG", "GT-AG", "GC-AG", "GT-AG"],
        "n_sites_annotated": [2, 2, 1, 0, 0],
        "direct_repeat":     [0, 0, 0, 0, 0],
        "n_reads":           [100, 100, 40, 40, 40],
        "n_mols":            [50, 50, 20, 20, 20],
        "min_site_share":    [1.0, 1.0, 1.0, 1.0, 1.0],
    })


def test_a_novel_semi_canonical_junction_is_rejected():
    out = curate(_feat(), JunctionParams())
    assert list(out.keep) == [True, True, True, False, True]
    assert out.loc[3, "drop_reason"] == "motif_class_not_allowed_for_novel"


def test_an_annotated_semi_canonical_junction_is_untouched():
    """A real GC-AG intron that GENCODE knows about must survive."""
    out = curate(_feat(), JunctionParams())
    assert out.loc[0, "keep"]
    assert out.loc[0, "drop_reason"] == ""


def test_the_motif_gate_can_be_widened_again():
    p = JunctionParams(novel_motif_classes=("canonical", "semi_canonical", "unknown"))
    assert curate(_feat(), p).loc[3, "keep"]


def test_unknown_motifs_survive_when_no_genome_was_supplied():
    f = _feat()
    f["motif_class"] = ["unknown"] * len(f)
    assert curate(f, JunctionParams()).keep.all()


def test_requiring_an_annotated_site_is_opt_in():
    f = _feat()
    assert curate(f, JunctionParams()).loc[4, "keep"]
    strict = curate(f, JunctionParams(require_annotated_site_for_novel_junction=True))
    assert not strict.loc[4, "keep"]
    assert strict.loc[4, "drop_reason"] == "no_annotated_splice_site"
    assert strict.loc[2, "keep"]           # reuses one annotated site


def test_novel_semi_canonical_junctions_are_decoys_for_the_calibration():
    """Without this the model saw annotated GC-AG in the positive set, learned
    the motif was fine, and scored 97.1% of the novel ones at lfdr <= 0.05."""
    f = _feat()
    _pos, decoy, _exempt = junction_anchors(f)
    assert decoy[3]                        # novel semi-canonical
    assert not decoy[0]                    # annotated semi-canonical
    assert not decoy[2] and not decoy[4]   # novel canonical


def test_short_read_support_still_rescues_a_semi_canonical_junction_from_the_decoys():
    f = _feat()
    f["sr_supported"] = [False, False, False, True, False]
    _pos, decoy, _exempt = junction_anchors(f)
    assert not decoy[3]
