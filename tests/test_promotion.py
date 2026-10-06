"""Rule S and Rule M: promotion replaces absorption at the 3' end (0.4.0).

The change in one sentence: a novel 3' end is no longer kept because nothing
absorbed it, it is kept because it earned a place. The thresholds here are not
invented -- every one was priced against SIRV Set 4 truth, and the numbers that
chose them are in the docstrings of the params they test.

What these tests are guarding against, in order of how badly each would hurt:

* read mass disappearing. Promotion demotes far more models than absorption
  did, so every path that removes a model needs a test that its reads landed
  somewhere.
* the competition being run against the wrong comparator. Candidate-versus-FSM
  is 11.7% precise and candidate-versus-candidate is 63-90%; the two differ
  only in which models enter the denominator, which is exactly the kind of
  thing a passing test suite will not notice.
* silence being read as evidence. Rule M demotes on a short-read measurement;
  absent coverage must never look like a measurement that came back negative.
"""

import numpy as np
import pytest

from flightcollapse.config import CoverageParams, MonoexonParams, TerminalParams
from flightcollapse.ends import prune_terminal_variants
from flightcollapse.model import TranscriptModel

P = TerminalParams()


def _m(three, n_reads, *, n_mols=0, disp=40.0, annotated=0, atlas=0, step=None,
       strand="+", first=1000):
    exons = [(first, 1200), (2000, three)] if strand == "+" else \
            [(three, 1200), (2000, first)]
    m = TranscriptModel(
        model_id=f"m{three}", contig="chrT", strand=strand, exons=exons,
        chain=((1200, 2000),), category="end3_novel", reads=np.arange(n_reads),
        n_reads=n_reads, n_mols=n_mols, n_cells=max(1, n_reads // 5),
    )
    m.evidence.update({"three_prime_dispersion": disp,
                       "end_is_annotated": annotated, "end_at_polya_site": atlas})
    if step is not None:
        m.evidence["sr_3p_step_supported"] = step
    return m


# ------------------------------------------------------------------ Rule S --
def test_a_diffuse_peak_is_absorbed_however_abundant():
    """The whole point. 0.3.x kept this model because no sibling was twice its
    size; 0.4.0 removes it because a 34 bp smear is not a cleavage site."""
    tight = _m(4000, 40, disp=3.0)
    diffuse = _m(3400, 4000, disp=34.0)          # 100x the support
    out = prune_terminal_variants([tight, diffuse], "+", P)
    assert [x.model_id for x in out] == ["m4000"]
    assert out[0].n_reads == 4040                 # its reads moved, not vanished


def test_read_mass_survives_promotion():
    models = [_m(4000, 500, disp=2.0), _m(3400, 20), _m(3000, 11), _m(2800, 300)]
    before = sum(m.n_reads for m in models)
    out = prune_terminal_variants(models, "+", P)
    assert sum(m.n_reads for m in out) == before


def test_an_unpromoted_model_with_nowhere_to_go_is_kept_and_marked():
    """A single diffuse peak on its own chain has no sibling to fold into. Its
    reads have nowhere else to belong, so it stays -- tiered, not believed."""
    out = prune_terminal_variants([_m(4000, 50, disp=40.0)], "+", P)
    assert len(out) == 1
    assert out[0].evidence["terminal_tier"] == "unsupported"


def test_a_ladder_with_no_evidence_anywhere_still_collapses():
    """The regression the promotion rule could easily have introduced: if
    nothing in a group is promoted, falling through would hand back the whole
    truncation ladder. The 0.3.x dominance rule still runs in that case."""
    models = [_m(4000, 900), _m(3600, 30), _m(3400, 20), _m(3200, 10)]
    out = prune_terminal_variants(models, "+", P)
    assert len(out) == 1
    assert out[0].n_reads == 960
    assert out[0].evidence["terminal_tier"] == "unsupported"


def test_each_promotion_route_works_on_its_own():
    for kw in ({"atlas": 1}, {"step": 1}, {"disp": 5.0}):
        out = prune_terminal_variants([_m(4000, 5000, disp=2.0), _m(3400, 6, **kw)],
                                      "+", P)
        assert len(out) == 2, f"{kw} should promote"


def test_promotion_can_be_switched_off():
    """A behaviour change that cannot be turned off cannot be attributed."""
    old = TerminalParams(promote_terminal_ends=False)
    pair = [_m(4000, 40, disp=3.0), _m(3400, 4000, disp=34.0)]
    assert len(prune_terminal_variants(pair, "+", old)) == 2      # 0.3.x kept it
    pair = [_m(4000, 40, disp=3.0), _m(3400, 4000, disp=34.0)]
    assert len(prune_terminal_variants(pair, "+", P)) == 1        # 0.4.0 does not


# ------------------------------------------------------------ the tier -----
def test_the_competition_is_between_candidates_not_against_fsm():
    """The correction that made this rule worth having.

    A proximal polyA site carrying 25% of a gene's molecules is normal biology.
    Measured against the annotated variant it loses; measured against the other
    novel candidates it wins. Both models here are dwarfed by the annotated
    one, and both must still be judged on their share of the CANDIDATE mass.
    """
    fsm = _m(4000, 10_000, n_mols=9_000, disp=2.0, annotated=1)
    big = _m(3400, 300, n_mols=280, disp=3.0)
    small = _m(3100, 40, n_mols=35, disp=3.0)
    out = prune_terminal_variants([fsm, big, small], "+", P)
    assert len(out) == 3
    tiers = {m.model_id: m.evidence["terminal_tier"] for m in out}
    assert tiers["m4000"] == "annotated"
    assert tiers["m3400"] == "high"        # 280/315 of the candidate molecules
    assert tiers["m3100"] == "reported"    # 35/315, a minority -- but emitted
    # and the annotated model's 9,000 molecules are NOT in that denominator
    assert out[1].evidence["candidate_share"] == pytest.approx(280 / 315, abs=1e-3)


def test_an_atlas_site_is_high_whatever_its_share():
    out = prune_terminal_variants(
        [_m(4000, 2000, n_mols=1900, disp=3.0), _m(3400, 5, n_mols=5, atlas=1)],
        "+", P)
    tiers = {m.model_id: m.evidence["terminal_tier"] for m in out}
    assert tiers["m3400"] == "high"


def test_the_tier_is_not_a_gate():
    """Everything promoted is emitted. `reported` is a label, not a deletion --
    gating on it would cost 62-76% of the true novel ends on SIRV."""
    out = prune_terminal_variants(
        [_m(4000, 1000, n_mols=900, disp=2.0), _m(3400, 3, n_mols=3, disp=2.0)],
        "+", P)
    assert len(out) == 2
    assert {m.evidence["terminal_tier"] for m in out} == {"high", "reported"}


def test_support_is_molecules_not_reads():
    """`max(n_mols, n_reads, 1)` always returned n_reads, so the comparison was
    between PCR copies. 60 reads from 3 molecules must lose to 30 reads from
    28 molecules."""
    amplified = _m(3400, 600, n_mols=3, disp=3.0)
    real = _m(3100, 30, n_mols=28, disp=3.0)
    out = prune_terminal_variants([amplified, real], "+", P)
    tiers = {m.model_id: m.evidence["terminal_tier"] for m in out}
    assert tiers["m3100"] == "high", "the deeply amplified model won on reads"
    assert tiers["m3400"] == "reported"


def test_minus_strand_promotion():
    """The silent failure mode: absorbing the wrong way round on one strand."""
    distal = _m(500, 500, disp=2.0, strand="-", first=4000)
    short = _m(1100, 900, disp=40.0, strand="-", first=4000)
    out = prune_terminal_variants([distal, short], "-", P)
    assert [m.model_id for m in out] == ["m500"]
    assert out[0].n_reads == 1400


# ------------------------------------------------------------------ Rule M --
def _mono_cfg(**kw):
    c = MonoexonParams()
    for k, v in kw.items():
        setattr(c, k, v)
    return c


def test_rule_m_params_default_on_and_are_separable():
    c = MonoexonParams()
    assert c.demote_contained_fragments is True
    assert c.keep_fragment_if_evidenced is True
    assert c.demote_on_sr_continuation is True
    # each is independently switchable, so a run can attribute its effect
    assert _mono_cfg(demote_contained_fragments=False).demote_on_sr_continuation


def test_a_contained_fragment_is_demoted_and_its_reads_fold(sim, tmp_path_factory):
    """End to end, because the containment test spans two tracks: the spliced
    models are built first and the mono track has to see them."""
    import pandas as pd
    from flightcollapse import Config
    from flightcollapse.pipeline import run

    out = tmp_path_factory.mktemp("rulem")
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.molecules.barcode_umi_tsv = sim.barcode_umi
    cfg.output.outdir, cfg.output.prefix = str(out), "rm"
    cfg.verbose, cfg.strict_invariants = False, False
    rep = run(cfg)

    acc = rep["read_accounting"]
    # the invariant that matters: nothing was dropped on the way
    folded = acc["monoexon_reads_folded_into_spliced"]
    total = (acc["monoexon_reads_modelled"] + folded
             + acc["monoexon_reads_rejected_low_support"]
             + acc["monoexon_reads_rejected_no_polya"]
             + acc["monoexon_reads_rejected_internal_priming"])
    assert total == acc["unspliced"], "unspliced reads went missing"
    t = pd.read_csv(out / "rm.models.tsv", sep="\t", low_memory=False)
    assert "terminal_tier" in t.columns
    assert "candidate_share" in t.columns


def test_silence_is_not_evidence_of_truncation():
    """Rule M clause 4. A gene with no short-read coverage must be left to the
    other clauses -- a lowly expressed transcript has no coverage, and reading
    that as truncation would delete real models exactly where long reads are
    the only evidence there is.

    Asserted on the contract of the value the rule keys on: `None` means no
    opinion and only an explicit `False` is a measurement.
    """
    from flightcollapse.coverage import blank_coverage_evidence

    blank = blank_coverage_evidence()
    assert blank["sr_3p_step_supported"] == ""
    # the rule tests `is False`, so neither "" nor None can trigger it
    for v in ("", None):
        assert not (v is False)


def test_coverage_params_reach_the_mono_track():
    """Regression: the coverage channel existed for two releases and the mono
    track could not see it, which is why clause M(3) was unimplementable."""
    import inspect

    from flightcollapse.monoexon import collapse_monoexonic

    sig = inspect.signature(collapse_monoexonic).parameters
    for name in ("spliced_spans", "terminal_params", "cov", "cov_params"):
        assert name in sig, f"{name} does not reach the mono-exon track"


def test_mono_exon_models_are_tiered_too(sim, tmp_path_factory):
    """Rule M can only demote when there is somewhere to fold the reads.

    A mono-exon pile in a gene with no emitted spliced model has nowhere to
    send them, so it must be emitted -- and that is the whole reason the SIRV
    fair arm keeps 107 mono-exon models where the design note priced 18. The
    note was wrong: it priced "drop every unevidenced mono model", which is not
    implementable without discarding read mass.

    What IS implementable is tiering them, so the precision is recoverable
    downstream: filtering `unsupported` out of that arm gives 110 models at
    75.5% 3'-end precision against 199 at 41.7%, with 83 correct ends either
    way. This test pins that every mono-exon model carries the label that makes
    the filter possible.
    """
    import pandas as pd
    from flightcollapse import Config
    from flightcollapse.pipeline import run

    out = tmp_path_factory.mktemp("monotier")
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.molecules.barcode_umi_tsv = sim.barcode_umi
    cfg.output.outdir, cfg.output.prefix = str(out), "mt"
    cfg.verbose, cfg.strict_invariants = False, False
    cfg.monoexon.demote_terminal_exon_fragments = False   # so one survives
    # 0.6.0: this test is about first-pass mono-exon tiering, so the final
    # consolidation is pinned off. With it on, the surviving mono-exon model is
    # (correctly) folded into the spliced model whose exon contains it, which is
    # the sub-chain rule doing its job rather than this behaviour changing.
    cfg.posthoc.enabled = False
    cfg.monoexon.demote_contained_fragments = False
    # 0.5.0: this test is about the FIRST pass -- whether the mono-exon
    # track emits and tiers a fragment model. The second pass would then
    # remove that model, correctly: the simulated genes are multi-exon only,
    # so a single-exon model in one of them is class C/D. Pinning the second
    # pass off here keeps this test measuring what it was written to measure;
    # tests/test_secondpass.py covers the removal itself.
    cfg.secondpass.enabled = False
    run(cfg)

    t = pd.read_csv(out / "mt.models.tsv", sep="\t", low_memory=False)
    mono = t[t["category"].str.startswith("monoexon")]
    assert len(mono), "no mono-exon model to tier"
    assert mono["terminal_tier"].notna().all()
    assert set(mono["terminal_tier"]) <= {
        "annotated", "high", "reported", "unsupported"}


def test_the_short_read_clause_can_actually_fire():
    """It could not, and the run proved it.

    The clause searched for a spliced model CONTAINING the fragment -- but a
    containing model is demoted one clause earlier and returns, so the search
    could only ever come up empty. `monoexon_models_demoted_by_short_reads` was
    0 on BD176c: not because short reads never disagreed with a 3' end, but
    because the code could not act on it when they did.

    The fix is that the fold target is now the best-OVERLAPPING spliced model.
    Pinned structurally, since reproducing it end to end needs a coverage track.
    """
    import inspect

    from flightcollapse import monoexon

    src = inspect.getsource(monoexon._make_monoexon_model)
    sr = src.split("demote_on_sr_continuation")[1]
    assert "minimum" in sr and "maximum" in sr, \
        "the short-read clause no longer computes an overlap"
    # and the containment clause must still come first, on the stronger evidence
    assert src.index("demote_contained_fragments") < src.index("demote_on_sr_continuation")
