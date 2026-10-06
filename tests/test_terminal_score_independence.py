"""The terminal score must not read the columns that define its own labels.

Both bugs pinned here shipped in 0.3.0 and survived 268 passing tests, because
every existing assertion asked whether the score EXISTED, not whether it meant
anything. They were found by running the thing on real data and looking at the
output, which is the reason these tests are written the way they are: each one
fails on a score that is merely present.

Measured on BD176c (155,192 models, flightcollapse 0.3.0):

* a ceiling of 0.05 agreed with the bare "this end is not catalogued" flag on
  93.0% of models -- the score was very nearly a rename of that flag;
* median ``three_prime_dispersion`` across lFDR bands ran
  4.0 -> 7.5 -> 16.5 -> 8.0 -> 11.0, NOT monotonic, so the discriminator the
  SIRV benchmark identified was contributing nothing, while
  ``perc_a_downstream`` ran a clean 15 -> 20 -> 25 -> 30 -> 40 because it *was*
  the decoy definition;
* every one of the 21,331 mono-exonic models sat above any ceiling -- 100.0% of
  monoexon_3UTR, monoexon_internal and monoexon_intergenic alike.
"""

import numpy as np
import pandas as pd
import pytest

from flightcollapse.config import ScoringParams
from flightcollapse.scoring import (
    TERMINAL_FEATURES,
    score_terminals,
    terminal_anchors,
    terminal_design,
)

#: Columns ``terminal_anchors`` thresholds to build ``pos`` and ``decoy``.
#: Anything derived from one of these is the label wearing a feature's name.
ANCHOR_COLUMNS = ("dist_ann_tts", "dist_polya_site", "perc_a_downstream",
                  "end_is_annotated", "end_at_polya_site")


def test_no_feature_is_derived_from_an_anchor_column():
    """The structural version of the bug, stated once.

    ``log_dist_ann_tts`` and ``log_dist_polya_site`` are ``pos`` after a
    threshold; ``perc_a_downstream`` is ``decoy`` after a threshold. A feature
    that defines the label dominates the fit and separates almost perfectly,
    and the local FDR stops being a statement about evidence.
    """
    for feat in TERMINAL_FEATURES:
        stem = feat.removeprefix("log_").removeprefix("logit_")
        assert stem not in ANCHOR_COLUMNS, (
            f"{feat} is derived from {stem}, which terminal_anchors uses to "
            f"build the labels. Scoring on it makes the lFDR circular.")


def test_the_design_matrix_refuses_a_feature_it_cannot_build():
    """A renamed feature must fail loudly rather than be silently all-zero."""
    import flightcollapse.scoring as sc

    old = sc.TERMINAL_FEATURES
    sc.TERMINAL_FEATURES = old + ("a_feature_that_does_not_exist",)
    try:
        with pytest.raises(KeyError, match="a_feature_that_does_not_exist"):
            terminal_design(pd.DataFrame({"n_reads": [1, 2, 3]}))
    finally:
        sc.TERMINAL_FEATURES = old


# --------------------------------------------------------------------------
def _frame(n=4000, seed=0, informative=False):
    """Ends whose catalogue status is real but whose EVIDENCE is pure noise.

    Half the rows are at a catalogued site, half are not. Every scoring feature
    is drawn from the same distribution for both halves, so there is nothing
    legitimate to separate them with. A score that nonetheless calls them apart
    is reading the labels.

    With ``informative=True`` the dispersion genuinely differs between the two,
    which is the case the score is supposed to detect.
    """
    rng = np.random.default_rng(seed)
    atlas = (np.arange(n) % 2).astype(int)
    disp = (rng.gamma(2.0, 4.0, n) if not informative
            else np.where(atlas == 1, rng.gamma(2.0, 2.0, n), rng.gamma(4.0, 9.0, n)))
    return pd.DataFrame({
        "n_reads": rng.poisson(30, n) + 1,
        "n_mols": rng.poisson(20, n) + 1,
        "n_cells": rng.poisson(8, n) + 1,
        "three_prime_dispersion": disp,
        "chain_share": rng.uniform(0.2, 1.0, n),
        "ratio_to_dominant_peak": rng.uniform(0.1, 1.0, n),
        "tail_molecule_frac": rng.uniform(0.5, 1.0, n),
        "median_tail_len": rng.integers(5, 40, n),
        "polya_motif_found": rng.integers(0, 2, n),
        "sr_3p_step_ratio": rng.uniform(0, 1, n),
        "n_3p_peaks": rng.integers(1, 5, n),
        # the label-defining columns, perfectly informative about the label
        "end_is_annotated": np.zeros(n, int),
        "end_at_polya_site": atlas,
        "dist_polya_site": np.where(atlas == 1, rng.integers(0, 40, n),
                                    rng.integers(5_000, 200_000, n)),
        "dist_ann_tts": rng.integers(5_000, 200_000, n),
        "perc_a_downstream": np.where(atlas == 1, rng.uniform(5, 25, n),
                                      rng.uniform(42, 70, n)),
    })


def _agreement(feat):
    """How often a 0.05 ceiling reproduces the bare 'not catalogued' flag."""
    res, _ = score_terminals(feat, ScoringParams(), decoy_perc_a=40.0)
    lf = pd.to_numeric(res["lfdr"], errors="coerce")
    pos, _, _ = terminal_anchors(feat, 40.0)
    ok = lf.notna().to_numpy()
    return float(((lf.to_numpy()[ok] > 0.05) == (~pos[ok])).mean())


def test_the_score_does_not_reproduce_the_catalogue_flag():
    """With no real evidence to go on, the score must NOT separate the classes.

    This is the test that fails on the 0.3.0 feature set: with the distances
    and %A in the design, agreement is essentially total, because the model can
    read the answer off the label definition.
    """
    assert _agreement(_frame()) < 0.90


def test_the_score_still_moves_when_the_evidence_is_real():
    """...and removing those features must not leave it inert.

    Dispersion now differs between catalogued and uncatalogued ends, which is a
    real signal the remaining features can see. The score should track it
    rather than sitting at one value.
    """
    feat = _frame(informative=True)
    res, _ = score_terminals(feat, ScoringParams(), decoy_perc_a=40.0)
    lf = pd.to_numeric(res["lfdr"], errors="coerce").dropna()
    assert lf.nunique() > 20, "the score collapsed to a constant"
    d = pd.DataFrame({"lfdr": lf, "disp": feat["three_prime_dispersion"][lf.index]})
    lo = d.loc[d["lfdr"] <= d["lfdr"].median(), "disp"].median()
    hi = d.loc[d["lfdr"] > d["lfdr"].median(), "disp"].median()
    assert hi > lo, ("tight peaks should score better than diffuse ones; "
                     f"got {lo:.1f} bp for the low-lFDR half against {hi:.1f}")


# --------------------------------------------------------------------------
def test_a_monoexonic_model_can_be_a_terminal_positive(sim, tmp_path_factory):
    """The second bug: mono-exonic models never carried the anchor flags.

    ``pos = end_is_annotated | end_at_polya_site`` was therefore False for every
    mono-exonic model by construction -- each one could be a decoy or an
    unlabelled candidate, never a positive -- and all 21,331 of them in BD176c
    landed above every ceiling. A category that is rejected unanimously is not
    being judged.
    """
    from flightcollapse import Config
    from flightcollapse.pipeline import run

    # Demotion is left ON everywhere else, so the default run emits no
    # mono-exonic model at all -- the one 3'UTR fragment is folded into its
    # spliced parent. Switch it off here so there IS a mono-exonic model whose
    # terminal fields can be inspected; the fields are set before demotion is
    # considered, so this tests the same code path.
    out = tmp_path_factory.mktemp("monoterm")
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.molecules.barcode_umi_tsv = sim.barcode_umi
    cfg.output.outdir, cfg.output.prefix = str(out), "mt"
    cfg.verbose, cfg.strict_invariants = False, False
    cfg.monoexon.demote_terminal_exon_fragments = False
    # 0.5.0: this test is about the FIRST pass -- whether the mono-exon
    # track emits and tiers a fragment model. The second pass would then
    # remove that model, correctly: the simulated genes are multi-exon only,
    # so a single-exon model in one of them is class C/D. Pinning the second
    # pass off here keeps this test measuring what it was written to measure;
    # tests/test_secondpass.py covers the removal itself.
    cfg.secondpass.enabled = False
    # 0.6.0: and the final consolidation off for the same reason -- its
    # sub-chain rule folds a mono-exon model lying inside a spliced model's
    # exon into it, which is exactly the fragment this test needs to survive.
    cfg.posthoc.enabled = False
    run(cfg)

    tbl = pd.read_csv(out / "mt.models.tsv", sep="\t", low_memory=False)
    mono = tbl[tbl["category"].str.startswith("monoexon")].copy()
    assert len(mono), "the simulation emitted no mono-exonic model to check"

    for key in ("end_is_annotated", "end_at_polya_site",
                "chain_share", "ratio_to_dominant_peak", "n_3p_peaks"):
        assert key in mono.columns, f"models.tsv has no {key}"
        assert mono[key].notna().all() and (mono[key].astype(str) != "").all(), (
            f"{key} is blank on a mono-exonic model -- it is emitted by the "
            f"spliced path only, which is the bug")

    # and the flags must be capable of being SET, not merely present: the
    # simulation puts a mono-exon cluster on an annotated 3' end on purpose
    for c in ("end_is_annotated", "end_at_polya_site", "perc_a_downstream"):
        mono[c] = pd.to_numeric(mono[c], errors="coerce").fillna(0)
    pos, _, _ = terminal_anchors(mono, 40.0)
    assert pos.any(), (
        "no mono-exonic model qualified as a terminal positive -- the flags are "
        "present but never set, which is the same bug one layer down")
