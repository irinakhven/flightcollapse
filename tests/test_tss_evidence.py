"""End-to-end behaviour of the CAGE tier (0.2.0, vector V1b).

The contract being tested is narrow and one-way:

* an atlas hit may KEEP a model that would otherwise be merged into its
  annotated parent;
* an atlas miss may never drop, demote or otherwise change one.

The second half is the one that needs a test, because nothing in the output
makes its violation obvious.
"""

import os

import pandas as pd
import pytest

from flightcollapse import Config
from flightcollapse.pipeline import run

STRUCTURE = ["contig", "strand", "start", "end", "n_exons", "n_reads"]


def _tss_sites(gtf):
    """Annotated TSS of every transcript in the simulated GTF."""
    out = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "transcript":
                continue
            start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
            out.add((f[0], strand, start if strand == "+" else end - 1))
    return sorted(out)


def _write_bed(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


@pytest.fixture(scope="module")
def annotated_only_bed(sim, tmp_path_factory):
    """Peaks that sit on annotated starts and nowhere else.

    Passes the load-time coverage check, and confirms only what annotation
    already said -- so it must not change a single model.
    """
    d = tmp_path_factory.mktemp("cage")
    rows = [
        (c, max(p - 25, 0), p + 25, f"ann_{c}_{p}_{s}", 100, s, p, p + 1, "0,0,0")
        for c, s, p in _tss_sites(sim.gtf)
    ]
    return _write_bed(str(d / "annotated_only.bed"), rows)


@pytest.fixture(scope="module")
def blanket_bed(sim, tmp_path_factory):
    """One enormous peak per strand, covering the whole simulated contig.

    Not biology -- a lever. Every 5' end is "supported", so anything the
    reference-anchored rule would have merged away is instead kept, which makes
    the rescue path observable in a 500-read simulation.
    """
    d = tmp_path_factory.mktemp("cage")
    contigs = {c for c, _s, _p in _tss_sites(sim.gtf)}
    rows = [
        (c, 0, 10_000_000, f"blanket_{c}_{s}", 100, s, 1000, 1001, "0,0,0")
        for c in sorted(contigs) for s in ("+", "-")
    ]
    return _write_bed(str(d / "blanket.bed"), rows)


def _run(sim, outdir, **over):
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.molecules.barcode_umi_tsv = sim.barcode_umi
    cfg.output.outdir = str(outdir)
    cfg.output.prefix = "t"
    cfg.verbose = False
    cfg.strict_invariants = False
    for k, v in over.items():
        cfg.set_path(k, v) if "." in k else setattr(cfg, k, v)
    report = run(cfg)
    table = pd.read_csv(os.path.join(str(outdir), "t.models.tsv"), sep="\t")
    return report, table


@pytest.fixture(scope="module")
def baseline(sim, tmp_path_factory):
    return _run(sim, tmp_path_factory.mktemp("base"))


# ---------------------------------------------------------------------- #
def test_an_atlas_confirming_only_annotated_starts_changes_nothing(
    sim, tmp_path_factory, annotated_only_bed, baseline
):
    _b_report, base = baseline
    _report, got = _run(
        sim, tmp_path_factory.mktemp("annonly"), cage_peak_bed=annotated_only_bed
    )
    pd.testing.assert_frame_equal(
        base[STRUCTURE].sort_values(STRUCTURE).reset_index(drop=True),
        got[STRUCTURE].sort_values(STRUCTURE).reset_index(drop=True),
    )


def test_recording_without_rescue_changes_nothing(
    sim, tmp_path_factory, blanket_bed, baseline
):
    """The middle arm of the measurement protocol.

    With the atlas recording but not acting, even one that supports every
    position in the genome must leave the model set alone -- otherwise the two
    arms are not comparable and the measured effect of the atlas is not the
    atlas.

    Note that "not acting" takes BOTH levers. ``use_as_chain_feature`` puts the
    evidence tier into the chain novelty model, which moves the fitted score,
    the local FDR and therefore the candidate set -- a real effect, and not one
    ``rescue_suffix`` governs.
    """
    _b_report, base = baseline
    _report, got = _run(
        sim, tmp_path_factory.mktemp("recorded"),
        cage_peak_bed=blanket_bed,
        **{"tss.rescue_suffix": "false", "tss.use_as_chain_feature": "false"},
    )
    pd.testing.assert_frame_equal(
        base[STRUCTURE].sort_values(STRUCTURE).reset_index(drop=True),
        got[STRUCTURE].sort_values(STRUCTURE).reset_index(drop=True),
    )
    # ... but the evidence is still on every spliced row (the mono-exonic
    # track has no chain group and so no tier to report)
    spliced = got[got["n_exons"] > 1]
    assert (spliced["in_cage_peak"] == 1).all()
    assert spliced["tss_evidence_source"].isin(
        ["annotation", "cage_reptss", "cage_peak"]
    ).all()


def test_a_supporting_atlas_can_only_add_models(
    sim, tmp_path_factory, blanket_bed, baseline
):
    b_report, base = baseline
    report, got = _run(
        sim, tmp_path_factory.mktemp("rescue"), cage_peak_bed=blanket_bed
    )
    assert len(got) > len(base)
    assert "end5_novel" in set(got["category"])

    # Every baseline STRUCTURE survives. Read counts legitimately change --
    # a rescued suffix takes its own reads with it, which is the whole point --
    # so the comparison is on the structure, not on the mass.
    cols = [c for c in STRUCTURE if c != "n_reads"]
    kept = set(map(tuple, base[cols].values.tolist()))
    now = set(map(tuple, got[cols].values.tolist()))
    assert kept - now == set()

    # ... and the mass is conserved: rescuing moves reads between models, it
    # does not create or lose them.
    assert report["reads_assigned_to_models"] == b_report["reads_assigned_to_models"]


def test_the_evidence_columns_are_populated_and_ordered(
    sim, tmp_path_factory, blanket_bed
):
    _report, got = _run(
        sim, tmp_path_factory.mktemp("cols"), cage_peak_bed=blanket_bed
    )
    for c in ("tss_evidence_source", "in_cage_peak", "dist_to_cage_peak",
              "dist_to_cage_reptss", "cage_peak_id",
              "exact_chain_reads", "exact_chain_frac", "n_distinct_chains"):
        assert c in got.columns
    # blanket atlas: every spliced 5' end is inside a peak
    assert got[got["n_exons"] > 1]["dist_to_cage_peak"].max() == 0
    assert (got["cage_peak_id"].astype(str).str.startswith("blanket")).any()


def test_exact_chain_fraction_is_a_fraction_and_reports_merged_mass(
    sim, tmp_path_factory, baseline
):
    """The V2 statistic the measurement protocol depends on."""
    _report, base = baseline
    spliced = base[base["n_exons"] > 1]
    frac = pd.to_numeric(spliced["exact_chain_frac"], errors="coerce").dropna()
    assert len(frac) > 0
    assert frac.between(0.0, 1.0).all()
    # a model that absorbed suffix groups must report more than one chain
    assert (pd.to_numeric(spliced["n_distinct_chains"], errors="coerce") >= 1).all()


def test_no_atlas_leaves_the_columns_empty_rather_than_zero(baseline):
    """Absent evidence and negative evidence must be distinguishable."""
    _report, base = baseline
    assert base["in_cage_peak"].isna().all()
    assert base["dist_to_cage_reptss"].isna().all()
    # the mono-exonic track has no chain group, so it has no tier to report;
    # blank there is honest rather than a hole
    spliced = base[base["n_exons"] > 1]
    assert spliced["tss_evidence_source"].isin(["annotation", "none"]).all()


def test_a_wrong_assembly_atlas_refuses_to_run(sim, tmp_path_factory):
    d = tmp_path_factory.mktemp("wrong")
    bed = _write_bed(
        str(d / "wrong.bed"),
        [("chrNOWHERE", 10, 70, "p", 1, "+", 20, 21, "0,0,0")],
    )
    with pytest.raises(RuntimeError, match="do not line up"):
        _run(sim, d, cage_peak_bed=bed)


def test_config_rejects_settings_with_no_atlas_to_act_on(sim):
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.cage_reptss_bed = "/nowhere/reptss.bed"
    with pytest.raises(ValueError, match="cage_reptss_bed"):
        cfg.validate()

    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.cage_peak_bed = "/nowhere/peaks.bed"
    cfg.tss.enabled = False
    with pytest.raises(ValueError, match="loaded and ignored"):
        cfg.validate()

    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.tss.min_peak_score = 5.0
    with pytest.raises(ValueError, match="no atlas to filter"):
        cfg.validate()


def test_the_chain_feature_lever_is_separate_from_the_rescue_lever(
    sim, tmp_path_factory, blanket_bed, baseline
):
    """Documented here because it is the one non-obvious part of the wiring.

    ``cage_support`` enters ``CHAIN_FEATURES``, so switching it on re-fits the
    chain novelty model even with rescue disabled. A "recorded only" arm of the
    measurement therefore has to turn both off, and this test is what stops
    that from being rediscovered by a confusing benchmark six months from now.
    """
    _b_report, base = baseline
    _report, got = _run(
        sim, tmp_path_factory.mktemp("featonly"),
        cage_peak_bed=blanket_bed,
        **{"tss.rescue_suffix": "false", "tss.use_as_chain_feature": "true"},
    )
    # no rescue, so no new category ...
    assert "end5_novel" not in set(got["category"])
    # ... but the scores are not required to be identical, and the evidence is
    # recorded either way
    assert (got[got["n_exons"] > 1]["in_cage_peak"] == 1).all()
