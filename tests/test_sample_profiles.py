"""Sample-type presets and the evidence-channel report (0.3.0).

The available external evidence is a property of the SAMPLE, not the algorithm,
and it changes between runs of the same pipeline: the retinal organoids have
matched short reads, the mouse-brain nuclei and HEK benchmarks do not, and the
SIRV spike-ins additionally have no CAGE atlas and no polyA-site atlas because
they are synthetic.

Two rules this file pins:

* a preset never changes the INTRINSIC algorithm. Terminal pruning and the
  terminal score compare a model against its own siblings, need no external
  input, and were motivated by the SIRV benchmark -- so they are on for every
  profile including the one meant for SIRVs.
* a declared input that is missing degrades the run and is REPORTED. It does
  not fail, and it does not silently pretend the evidence was there.
"""

import pytest

from flightcollapse.config import SAMPLE_PROFILES, Config


def _cfg(**kw):
    c = Config(bam="b.bam", reference_gtf="r.gtf", genome_fasta="g.fa")
    c.molecules.write_matrix = False      # no barcode table in these fixtures
    for k, v in kw.items():
        setattr(c, k, v)
    return c


def test_the_three_shipped_profiles_exist():
    assert set(SAMPLE_PROFILES) == {"retinal_organoid", "generic", "minimal"}
    for name, prof in SAMPLE_PROFILES.items():
        assert prof["description"]
        assert isinstance(prof["expects"], list)
        assert isinstance(prof["settings"], dict)


def test_organoids_expect_everything_including_short_reads():
    c = _cfg()
    info = c.apply_sample_profile("retinal_organoid")
    assert "short_read_sj" in info["expects"]
    assert "short_read_coverage" in info["expects"]
    assert c.coverage.enabled is True
    assert c.tss.enabled is True


def test_other_cell_types_get_everything_except_short_reads():
    """The user's rule: atlases yes, short reads no, for anything not organoid."""
    c = _cfg()
    info = c.apply_sample_profile("generic")
    assert "cage_peak_bed" in info["expects"]
    assert "polya_site_bed" in info["expects"]
    assert "short_read_sj" not in info["expects"]
    assert "short_read_coverage" not in info["expects"]
    assert c.coverage.enabled is False
    assert c.tss.enabled is True


def test_minimal_expects_nothing():
    c = _cfg()
    info = c.apply_sample_profile("minimal")
    assert info["expects"] == []
    assert info["missing_inputs"] == []
    assert c.tss.enabled is False
    assert c.coverage.enabled is False


@pytest.mark.parametrize("name", sorted(SAMPLE_PROFILES))
def test_no_profile_disables_the_intrinsic_fixes(name):
    """Pruning and the terminal score need no external evidence, so no sample
    type may switch them off -- least of all the SIRV-facing one, which is the
    benchmark that motivated them."""
    c = _cfg()
    c.apply_sample_profile(name)
    assert c.terminal.prune_terminal_variants is True
    assert c.terminal.score_terminals is True


def test_a_missing_declared_input_degrades_rather_than_raising():
    c = _cfg()
    info = c.apply_sample_profile("retinal_organoid")
    assert set(info["missing_inputs"]) == {
        "cage_peak_bed", "polya_site_bed", "short_read_sj", "short_read_coverage"}
    c.validate()          # must not raise: a sample without short reads is normal
    assert c.evidence_channels()["short_read_junctions"] is False


def test_supplying_the_inputs_clears_the_missing_list():
    c = _cfg(cage_peak_bed="c.bed", polya_site_bed="p.bed",
             short_read_sj=["sj.tab"], short_read_coverage=["cov.bw"])
    info = c.apply_sample_profile("retinal_organoid")
    assert info["missing_inputs"] == []
    ch = c.evidence_channels()
    assert ch["cage_atlas"] and ch["polya_atlas"]
    assert ch["short_read_junctions"] and ch["short_read_coverage"]


def test_evidence_channels_reports_what_was_actually_active():
    """Two catalogues that differ only in their inputs are not the same
    experiment, and nothing else in the outputs says so."""
    c = _cfg(short_read_sj=["sj.tab"])
    c.apply_sample_profile("generic")
    ch = c.evidence_channels()
    assert ch["sample_type"] == "generic"
    assert ch["short_read_junctions"] is True     # supplied even if not expected
    assert ch["cage_atlas"] is False              # declared, not supplied
    assert ch["short_read_coverage"] is False


def test_an_unknown_sample_type_is_refused():
    with pytest.raises(ValueError, match="unknown sample_type"):
        _cfg().apply_sample_profile("retina_ish")
    c = _cfg()
    c.sample_type = "nonsense"
    with pytest.raises(ValueError, match="sample_type"):
        c.validate()


def test_explicit_settings_survive_the_preset():
    """The CLI applies the profile before --set, so a hand-set value wins."""
    c = _cfg()
    c.apply_sample_profile("minimal")
    c.set_path("terminal.sibling_support_ratio", "5")
    assert c.terminal.sibling_support_ratio == 5.0


def test_coverage_set_but_disabled_is_refused():
    c = _cfg(short_read_coverage=["cov.bw"])
    c.apply_sample_profile("generic")     # generic turns coverage off
    with pytest.raises(ValueError, match="loaded and ignored"):
        c.validate()


def test_a_terminal_ceiling_without_scoring_is_refused():
    c = _cfg()
    c.terminal.score_terminals = False
    c.terminal.max_terminal_lfdr = 0.05
    with pytest.raises(ValueError, match="max_terminal_lfdr"):
        c.validate()
