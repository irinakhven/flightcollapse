"""A contig with no spliced read must not fail.

The regression, and the reason it needed two attempts.

``_chain_features`` builds its frame with ``pd.DataFrame(rows)``. On a contig
whose admitted reads are all unspliced there are no chain groups, ``rows`` is
empty, and the result has **no columns at all** -- so every downstream column
access raises ``KeyError`` and the runner reports the contig as FAILED.

The first fix guarded ``score_chains``'s access to ``n_novel_junctions``. That
moved the failure exactly one column along, to ``cfeat["gid"]`` in the caller:
same six contigs, same run, new KeyError. Guarding consumers one at a time
cannot terminate, because each guard only reveals the next consumer. The frame
itself has to carry its columns.

Observed on chrY (3 reads, 3 unspliced), five unplaced scaffolds, and any contig
where every admitted read is mono-exonic. Cosmetic at that scale, but it is a
whole contig silently dropped, and on a sample where chrY is real it would not
be cosmetic at all.
"""

import re
import inspect

import pandas as pd
import pytest

from flightcollapse import pipeline
from flightcollapse.config import ScoringParams
from flightcollapse.pipeline import _CHAIN_FEATURE_COLS, _chain_features
from flightcollapse.scoring import score_chains


def _empty():
    # groups=[] short-circuits the loop, so the other arguments are never read
    return _chain_features([], None, None, None, None, None, {})


def test_no_chain_groups_still_yields_a_usable_frame():
    f = _empty()
    assert len(f) == 0
    assert list(f.columns) == [c for c, _ in _CHAIN_FEATURE_COLS]


def test_the_two_columns_that_actually_broke():
    f = _empty()
    f["gid"].to_numpy()                 # the second failure: pipeline.py:419
    f["n_novel_junctions"].to_numpy()   # the first: scoring.py score_chains


def test_scoring_an_empty_frame_does_not_raise():
    score, model = score_chains(_empty(), ScoringParams())
    assert len(score) == 0
    assert not model.fitted


def test_every_consumed_column_is_declared():
    """The drift guard the source comment promises.

    Reading the builder's own dict literal is deliberate: a column added to the
    populated path but not to the empty one reintroduces this bug in a form that
    only appears on contigs the test suite has no reason to simulate.
    """
    src = inspect.getsource(_chain_features)
    body = src.split("rows.append(")[1]
    built = set(re.findall(r'"([a-z_0-9]+)":', body))
    declared = {c for c, _ in _CHAIN_FEATURE_COLS}
    assert built == declared, (
        f"only built: {sorted(built - declared)}   "
        f"only declared: {sorted(declared - built)}"
    )


def test_lfdr_lookup_over_an_empty_frame_is_empty_not_an_error():
    """The exact expression at pipeline.py:418-424."""
    cfeat = _empty()
    cscore, _ = score_chains(cfeat, ScoringParams())
    by_gid = {
        int(g): float(v)
        for g, v in zip(cfeat["gid"].to_numpy(), cscore["lfdr"].to_numpy())
    }
    assert by_gid == {}
