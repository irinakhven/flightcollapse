import json

import pytest

from flightcollapse.config import Config


def _minimal(**kw) -> Config:
    cfg = Config(bam="a.bam", reference_gtf="a.gtf", genome_fasta="a.fa", **kw)
    cfg.molecules.write_matrix = False
    return cfg


def test_minimal_config_validates():
    _minimal().validate()


def test_missing_genome_is_rejected_when_it_would_be_used():
    cfg = _minimal()
    cfg.genome_fasta = None
    with pytest.raises(ValueError, match="genome_fasta is required"):
        cfg.validate()


def test_no_genome_is_fine_once_the_dependent_features_are_off():
    cfg = _minimal()
    cfg.genome_fasta = None
    cfg.junctions.require_canonical_motif = False
    cfg.junctions.rt_switch_repeat_len = 0
    cfg.monoexon.enabled = False
    cfg.ends.require_polya_evidence_for_utr_variant = False
    cfg.scoring.enabled = False
    cfg.validate()


def test_matrix_without_a_barcode_table_is_rejected():
    """A flag that is accepted but silently ignored is worse than one that is
    rejected -- this is the --max-5p-diff lesson."""
    cfg = _minimal()
    cfg.molecules.write_matrix = True
    with pytest.raises(ValueError, match="barcode_umi_tsv"):
        cfg.validate()


def test_five_prime_percentile_guardrail():
    cfg = _minimal()
    cfg.ends.five_prime_percentile = 10
    with pytest.raises(ValueError, match="isoseq collapse bug"):
        cfg.validate()


def test_ratio_guardrails():
    cfg = _minimal()
    cfg.chains.debris_ratio = 1.5
    with pytest.raises(ValueError):
        cfg.validate()
    cfg = _minimal()
    cfg.chains.dominant_ratio = 0.5
    with pytest.raises(ValueError):
        cfg.validate()


def test_set_path_coerces_types():
    cfg = _minimal()
    cfg.set_path("ends.max_3p_diff", "50")
    cfg.set_path("scoring.enabled", "false")
    cfg.set_path("junctions.max_junction_lfdr", "0.01")
    cfg.set_path("junctions.max_junction_lfdr", "none")
    assert cfg.ends.max_3p_diff == 50
    assert cfg.scoring.enabled is False
    assert cfg.junctions.max_junction_lfdr is None


def test_set_path_rejects_unknown_keys():
    with pytest.raises(ValueError, match="unknown config key"):
        _minimal().set_path("ends.no_such_thing", "1")


def test_json_roundtrip(tmp_path):
    cfg = _minimal()
    cfg.ends.max_3p_diff = 42
    p = tmp_path / "c.json"
    cfg.dump(str(p))
    back = Config.load(str(p))
    assert back.ends.max_3p_diff == 42
    assert back.bam == cfg.bam


def test_cli_config_subcommand(capsys):
    from flightcollapse.cli import main

    assert main(["config"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert "junctions" in d and "chains" in d
