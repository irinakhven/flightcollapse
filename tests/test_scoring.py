"""The calibrated novelty model.

These tests build a junction feature table with a known truth: 'real' novel
junctions drawn with high support / canonical motif / high site share, and
'artefact' junctions drawn with low support, non-canonical motifs and long
direct repeats.  The model never sees the truth labels -- it anchors on
*annotated* junctions as positives and on the decoy definition as negatives --
so recovering the truth is a real test of the calibration, not a tautology.
"""

import numpy as np
import pandas as pd
import pytest

from flightcollapse.config import ScoringParams
from flightcollapse.scoring import (
    fit_logistic,
    junction_anchors,
    junction_design,
    local_fdr,
    score_junctions,
)


def _junctions(n_ann=400, n_real=200, n_art=400, seed=0):
    rng = np.random.default_rng(seed)

    def block(n, annotated, real, support_lo, support_hi, canon_p, repeat_lo, repeat_hi):
        return pd.DataFrame(
            {
                "n_reads": rng.integers(support_lo, support_hi, n),
                "n_mols": rng.integers(support_lo, support_hi, n),
                "n_cells": rng.integers(max(support_lo // 2, 1), support_hi, n),
                "annotated": annotated,
                "truth_real": real,
                "donor_share": rng.beta(6, 2, n) if real else rng.beta(1, 12, n),
                "acceptor_share": rng.beta(6, 2, n) if real else rng.beta(1, 12, n),
                "intron_len": rng.integers(200, 40_000, n),
                "direct_repeat": rng.integers(repeat_lo, repeat_hi, n),
                "motif_class": np.where(
                    rng.random(n) < canon_p, "canonical", "non_canonical"
                ),
                "dist_ann_donor": 0 if annotated else rng.integers(30, 5000, n),
                "dist_ann_acceptor": 0 if annotated else rng.integers(30, 5000, n),
                "n_sites_annotated": 2 if annotated else rng.integers(0, 3, n),
            }
        )

    df = pd.concat(
        [
            block(n_ann, True, True, 20, 4000, 0.99, 0, 3),
            block(n_real, False, True, 15, 600, 0.97, 0, 4),
            block(n_art, False, False, 1, 6, 0.15, 6, 18),
        ],
        ignore_index=True,
    )
    return df


def test_logistic_separates_a_linearly_separable_problem():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(-2, 1, 300), rng.normal(2, 1, 300)])
    y = np.concatenate([np.zeros(300), np.ones(300)])
    X = np.column_stack([np.ones(600), x])
    w = fit_logistic(X, y, l2=1.0)
    assert w[1] > 0
    pred = (X @ w > 0).astype(float)
    assert (pred == y).mean() > 0.9


def test_local_fdr_is_bounded_and_monotone():
    rng = np.random.default_rng(1)
    decoy = rng.normal(0, 1, 3000)
    cand = np.concatenate([rng.normal(0, 1, 2000), rng.normal(4, 1, 1000)])
    lf, pi0 = local_fdr(cand, decoy, n_bins=50)
    assert np.all((lf >= 0) & (lf <= 1))
    assert 0.0 <= pi0 <= 1.0
    order = np.argsort(cand)
    assert np.all(np.diff(lf[order]) <= 1e-9)          # non-increasing in score
    assert lf[cand > 5].mean() < 0.2                    # clear signal is confident
    assert lf[cand < -1].mean() > 0.6                   # noise region is not


def test_model_separates_real_from_artefact_novel_junctions():
    df = _junctions()
    res, model = score_junctions(df, ScoringParams())
    assert model.fitted, model.note
    novel = ~df["annotated"].to_numpy()
    real = df["truth_real"].to_numpy()[novel]
    lf = res["lfdr"].to_numpy()[novel]
    assert np.isfinite(lf).all()
    # a 5% local-FDR ceiling should keep most real ones and almost no artefacts
    kept = lf <= 0.05
    recall = kept[real].mean()
    fpr = kept[~real].mean()
    assert recall > 0.7, recall
    assert fpr < 0.05, fpr


def test_annotated_junctions_are_never_discoveries():
    df = _junctions()
    res, _ = score_junctions(df, ScoringParams())
    assert (res["lfdr"].to_numpy()[df["annotated"].to_numpy()] == 0).all()


def test_falls_back_cleanly_when_anchors_are_too_small():
    df = _junctions(n_ann=5, n_real=5, n_art=5)
    res, model = score_junctions(df, ScoringParams())
    assert not model.fitted
    assert "too few anchors" in model.note
    assert res["lfdr"].isna().all()


def test_disabled_scoring_returns_nan():
    df = _junctions()
    res, model = score_junctions(df, ScoringParams(enabled=False))
    assert res["lfdr"].isna().all()
    assert not model.fitted


def test_decoy_definition_targets_artefacts():
    df = _junctions()
    pos, decoy, _exempt = junction_anchors(df)
    assert pos.sum() == int(df["annotated"].sum())
    # decoys should be overwhelmingly drawn from the artefact block
    assert df.loc[decoy, "truth_real"].mean() < 0.15


def test_design_matrix_has_no_nans():
    X = junction_design(_junctions())
    assert np.isfinite(X).all()
