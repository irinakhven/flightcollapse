"""End-to-end tests against the simulated dataset with ground truth.

The simulation deliberately reproduces the conditions that broke
``isoseq collapse``: 5'-truncated reads outnumbering intact ones almost 3:1,
mono-exonic reads sitting inside the last exon of a spliced gene, and singleton
novel junctions scattered across loci.
"""

import gzip
import os

import numpy as np
import pandas as pd
import pytest


def _read_stat(outdir):
    p = os.path.join(outdir, "sim.read_stat.txt.gz")
    if not os.path.exists(p):
        p = os.path.join(outdir, "sim.read_stat.txt")
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return pd.read_csv(fh, sep="\t")


# ---------------------------------------------------------------------- #
def test_all_annotated_transcripts_are_recovered(collapsed, models_table):
    truth = collapsed["sim"].transcripts
    fsm = set(models_table.loc[models_table["category"] == "FSM", "model_id"])
    missing = set(truth) - fsm
    assert not missing, f"annotated transcripts not recovered: {sorted(missing)}"


def test_truncated_reads_land_on_the_full_length_model(models_table):
    """25 intact + 70 truncated reads must end up as one model with 95 reads,
    carrying the full 4-exon structure -- not as a 2-exon stub with 70."""
    row = models_table.set_index("model_id").loc["ENSTSIM00011"]
    assert row["n_exons"] == 4
    assert row["n_reads"] == 95


def test_no_truncation_derived_novel_models(models_table):
    """Every emitted model must be either an annotated transcript, a real novel
    isoform, a UTR variant or a mono-exon model. A 5'-truncated fragment of a
    known transcript is none of those."""
    novel = models_table[models_table["category"].isin(["NIC", "NNC"])]
    assert len(novel) == 1, novel[["model_id", "category", "n_reads"]].to_string()
    assert novel.iloc[0]["n_reads"] >= 30


def test_model_exon_counts_are_never_below_their_reads(collapsed):
    inv = {r["invariant"]: r for r in collapsed["report"]["invariants"]}
    assert inv["model_is_not_always_the_shortest_member"]["passed"]
    assert inv["model_is_not_always_the_shortest_member"]["observed"] < 0.95


def test_step0_invariant_holds(collapsed):
    """The number that exposed the original bug: mono-exon model read mass must
    track the genuinely unspliced read fraction."""
    rep = collapsed["report"]
    inv = {r["invariant"]: r for r in rep["invariants"]}
    r = inv["monoexon_read_mass_matches_unspliced_reads"]
    assert r["passed"]
    assert r["observed"] <= 2 * rep["frac_reads_unspliced"] + 1e-9


def test_all_invariants_pass(collapsed):
    bad = [r for r in collapsed["report"]["invariants"] if not r["passed"]]
    assert not bad, bad


# ---------------------------------------------------------------------- #
def test_genuine_novel_skipping_isoform_is_called(models_table):
    novel = models_table[models_table["category"].isin(["NIC", "NNC"])]
    assert (novel["gene_name"] == "Sim2").any()
    assert novel.iloc[0]["n_exons"] == 3


def test_singleton_spurious_junctions_are_rejected(collapsed):
    """Six one-read non-canonical junctions were simulated across six loci."""
    j = pd.read_csv(os.path.join(collapsed["outdir"], "sim.junctions.tsv.gz"), sep="\t")
    novel = j[~j["annotated"]]
    assert (novel["n_reads"] == 1).any()
    assert not novel.loc[novel["n_reads"] == 1, "keep"].any()


def test_pcr_duplicated_novel_junction_is_rejected(collapsed):
    """30 reads but a single molecule. Read counts say 'real', molecules say no."""
    j = pd.read_csv(os.path.join(collapsed["outdir"], "sim.junctions.tsv.gz"), sep="\t")
    dup = j[(~j["annotated"]) & (j["n_reads"] >= 25) & (j["n_mols"] <= 2)]
    assert len(dup) == 1, dup.to_string()
    assert not dup.iloc[0]["keep"]
    assert dup.iloc[0]["drop_reason"] == "low_molecules"


def test_alternative_polya_site_becomes_a_novel_3p_end(models_table):
    v = models_table[models_table["category"] == "end3_novel"]
    assert len(v) == 1
    assert v.iloc[0]["model_id"] == "ENSTSIM00011|novel3end1"
    assert v.iloc[0]["parent_transcript"] == "ENSTSIM00011"
    assert bool(v.iloc[0]["polya_motif_found"])
    assert v.iloc[0]["n_reads"] >= 25


def test_monoexon_3utr_fragment_is_kept(models_table):
    m = models_table[models_table["category"] == "monoexon_3UTR"]
    assert len(m) == 1
    assert m.iloc[0]["n_exons"] == 1
    assert m.iloc[0]["n_reads"] >= 30
    assert "rt_dropoff_fragment" in str(m.iloc[0]["flags"])


def test_monoexon_gene_body_fragment_without_polya_is_rejected(collapsed, models_table):
    assert (models_table["category"] == "monoexon_internal").sum() == 0
    mono = collapsed["report"]["per_contig"]["chr1"]["monoexon"]
    assert mono["rejected_no_polya"] >= 1


# ---------------------------------------------------------------------- #
def test_outputs_exist(collapsed):
    out = collapsed["outdir"]
    for f in (
        "sim.gtf", "sim.gff", "sim.models.tsv", "sim.abundance.txt",
        "sim.flnc_count.txt", "sim.junctions.tsv.gz", "sim.group_chains.tsv",
        "sim.qc.json", "sim.qc.md", "sim.config.json",
        "sim_umi_corrected_isoform_matrix.mtx", "sim_cells.txt",
        "sim_isoforms.txt", "sim_isoform_metadata.txt",
    ):
        assert os.path.exists(os.path.join(out, f)), f


def test_read_stat_is_read_name_keyed_and_unique(collapsed, models_table):
    rs = _read_stat(collapsed["outdir"])
    assert list(rs.columns[:5]) == ["id", "length", "is_fl", "stat", "pbid"]
    assert rs["id"].is_unique
    assert rs["pbid"].isin(set(models_table["model_id"])).all()
    # every model's read count must match its read_stat rows exactly
    counts = rs["pbid"].value_counts()
    for r in models_table.itertuples(index=False):
        assert counts.get(r.model_id, 0) == r.n_reads, r.model_id


def test_group_file_matches_read_stat(collapsed):
    p = os.path.join(collapsed["outdir"], "sim.group.txt.gz")
    if not os.path.exists(p):
        p = os.path.join(collapsed["outdir"], "sim.group.txt")
    op = gzip.open if p.endswith(".gz") else open
    total = 0
    with op(p, "rt") as fh:
        for line in fh:
            _mid, names = line.rstrip("\n").split("\t")
            total += len(names.split(","))
    assert total == len(_read_stat(collapsed["outdir"]))


def test_gtf_is_parseable_and_sorted(collapsed):
    rows = []
    with open(os.path.join(collapsed["outdir"], "sim.gtf")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            assert len(f) == 9
            if f[2] == "transcript":
                rows.append((f[0], int(f[3])))
    assert rows == sorted(rows)
    assert len(rows) > 0


def test_molecule_matrix_dimensions(collapsed, models_table):
    out = collapsed["outdir"]
    with open(os.path.join(out, "sim_umi_corrected_isoform_matrix.mtx")) as fh:
        assert fh.readline().startswith("%%MatrixMarket")
        n_cells, n_iso, nnz = (int(x) for x in fh.readline().split())
        entries = [tuple(int(v) for v in l.split()) for l in fh]
    cells = open(os.path.join(out, "sim_cells.txt")).read().split()
    isos = open(os.path.join(out, "sim_isoforms.txt")).read().split()
    assert len(cells) == n_cells
    assert len(isos) == n_iso == len(models_table)
    assert len(entries) == nnz
    assert max(r for r, _c, _v in entries) <= n_cells
    assert max(c for _r, c, _v in entries) <= n_iso
    # molecules per model must not exceed reads per model
    per_iso = {}
    for _r, c, v in entries:
        per_iso[c] = per_iso.get(c, 0) + v
    by_id = models_table.set_index("model_id")
    for c, tot in per_iso.items():
        assert tot <= int(by_id.loc[isos[c - 1], "n_reads"])


def test_config_roundtrip(collapsed):
    from flightcollapse.config import Config

    cfg = Config.load(os.path.join(collapsed["outdir"], "sim.config.json"))
    assert cfg.ends.five_prime_percentile == collapsed["cfg"].ends.five_prime_percentile


# ---------------------------------------------------------------------- #
# regression: a late failure must not discard completed work
# ---------------------------------------------------------------------- #
def test_tiny_contig_does_not_trip_the_gate_guard(sim, tmp_path):
    """chrY in a female line has a few hundred repeat-mapping alignments.
    Rejecting 99% of those is the gates working, not a misconfiguration —
    and it must not kill a run that has already done 23 chromosomes."""
    from flightcollapse import Config
    from flightcollapse.pipeline import run

    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.output.outdir = str(tmp_path / "out")
    cfg.output.prefix = "t"
    cfg.verbose = False
    cfg.strict_invariants = False
    cfg.molecules.write_matrix = False
    # a gate nothing can pass, on a read count below the guard's floor
    cfg.gates.min_aln_coverage = 1.01
    report = run(cfg)
    assert report["failed_contigs"] == {}
    assert os.path.exists(os.path.join(cfg.output.outdir, "t.qc.json"))


def test_a_failing_contig_still_leaves_the_others_written(sim, tmp_path, monkeypatch):
    from flightcollapse import Config
    from flightcollapse import pipeline as pl

    real = pl._run_contig

    def boom(cfg, bam, contig, *a, **kw):
        if contig == "chr1":
            raise RuntimeError("synthetic failure")
        return real(cfg, bam, contig, *a, **kw)

    monkeypatch.setattr(pl, "_run_contig", boom)
    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.output.outdir = str(tmp_path / "out2")
    cfg.output.prefix = "t"
    cfg.verbose = False
    cfg.strict_invariants = False
    cfg.molecules.write_matrix = False
    report = run_quiet(cfg)
    assert "chr1" in report["failed_contigs"]
    assert not report["all_passed"]
    # the run still produced its outputs rather than losing everything
    for f in ("t.qc.json", "t.qc.md", "t.gtf", "t.models.tsv", "t.config.json"):
        assert os.path.exists(os.path.join(cfg.output.outdir, f)), f


def run_quiet(cfg):
    from flightcollapse.pipeline import run

    return run(cfg)


def test_gate_guard_still_fires_on_a_real_contig(sim, tmp_path):
    """The guard must not have been defanged: above the read floor it still stops."""
    from flightcollapse import Config
    from flightcollapse.pipeline import run

    cfg = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    cfg.output.outdir = str(tmp_path / "out3")
    cfg.output.prefix = "t"
    cfg.verbose = False
    cfg.strict_invariants = False
    cfg.molecules.write_matrix = False
    cfg.gates.min_aln_coverage = 1.01
    cfg.gates.min_reads_for_gate_check = 10     # below the sim's read count
    report = run(cfg)
    assert "chr1" in report["failed_contigs"]
    assert "GateRejectionError" in report["failed_contigs"]["chr1"]


def test_model_ids_are_unique(models_table, collapsed):
    """A duplicate transcript_id silently corrupts read_stat, group.txt and the
    count matrix. Several 3' peaks folding back onto the same annotated end used
    to produce exactly that."""
    dup = models_table[models_table.duplicated("model_id", keep=False)]
    assert dup.empty, dup[["model_id", "category", "start", "end", "flags"]].to_string()
    inv = {r["invariant"]: r for r in collapsed["report"]["invariants"]}
    assert inv["model_ids_are_unique"]["passed"]
    assert collapsed["report"]["model_id_collisions"] == 0


def test_one_model_per_final_three_prime_end(models_table):
    """Two models of the same transcript may differ at the 3' end, never at the
    5' end alone — the 5' end is a percentile statistic, not a distinguishing
    feature."""
    key = models_table.groupby(["contig", "strand", "parent_transcript", "n_exons"])
    for (_c, strand, _p, _n), grp in key:
        three = grp["end"] if strand == "+" else grp["start"]
        assert three.is_unique, grp[["model_id", "start", "end", "flags"]].to_string()
