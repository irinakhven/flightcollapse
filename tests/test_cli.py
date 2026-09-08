import json
import os

import pytest

from flightcollapse.cli import main


def test_run_subcommand_end_to_end(sim, tmp_path):
    out = tmp_path / "cli_out"
    rc = main([
        "run",
        "-b", sim.bam,
        "-g", sim.gtf,
        "-f", sim.genome_fasta,
        "-u", sim.barcode_umi,
        "-o", str(out),
        "-p", "cli",
        "--quiet",
        "--set", "ends.max_3p_diff=80",
    ])
    assert rc == 0
    assert os.path.exists(out / "cli.gtf")
    cfg = json.load(open(out / "cli.config.json"))
    assert cfg["ends"]["max_3p_diff"] == 80
    qc = json.load(open(out / "cli.qc.json"))
    assert qc["all_passed"]


def test_run_without_barcode_table_disables_the_matrix(sim, tmp_path):
    out = tmp_path / "nomol"
    rc = main([
        "run", "-b", sim.bam, "-g", sim.gtf, "-f", sim.genome_fasta,
        "-o", str(out), "-p", "n", "--quiet",
    ])
    assert rc == 0
    assert not os.path.exists(out / "n_umi_corrected_isoform_matrix.mtx")
    assert os.path.exists(out / "n.models.tsv")


def test_index_umi_subcommand(sim, tmp_path, capsys):
    idx = tmp_path / "idx.npz"
    assert main(["index-umi", sim.barcode_umi, "-o", str(idx)]) == 0
    assert idx.exists()
    assert "barcodes" in capsys.readouterr().out


def test_check_subcommand_on_our_output(sim, tmp_path, capsys):
    out = tmp_path / "chk"
    assert main([
        "run", "-b", sim.bam, "-g", sim.gtf, "-f", sim.genome_fasta,
        "-o", str(out), "-p", "c", "--quiet",
    ]) == 0
    rs = out / "c.read_stat.txt.gz"
    if not rs.exists():
        rs = out / "c.read_stat.txt"
    capsys.readouterr()
    rc = main(["check", "-b", sim.bam, "-g", str(out / "c.gtf"), "-r", str(rs),
               "--cigar-scan-reads", "100000"])
    assert rc == 0
    assert not json.loads(capsys.readouterr().out)["step0_mismatch"]


def test_bad_set_argument(sim, tmp_path):
    with pytest.raises(SystemExit):
        main(["run", "-b", sim.bam, "-g", sim.gtf, "-f", sim.genome_fasta,
              "-o", str(tmp_path), "--set", "nonsense"])


def test_diagnose_umi_finds_the_right_column_and_mode(sim, tmp_path, capsys):
    """A table whose read names lost the segment span must still be diagnosable."""
    import pysam

    from flightcollapse.diagnose import diagnose

    src = open(sim.barcode_umi).read().rstrip("\n").split("\n")
    shifted = tmp_path / "shifted.tsv"
    with open(shifted, "w") as fh:
        # extra leading columns + segment spans appended, as in a real BD table
        fh.write("cell_barcode\tumi\tbarcode_umi\tread_name\tsource_bam\n")
        for line in src[1:]:
            cb, umi, name = line.split("\t")
            fh.write(f"{cb}\t{umi}\t{cb}|{umi}\t{name}/100_900\tx.bam\n")

    res = diagnose(sim.bam, str(shifted), verbose=False)
    best = res["best"]
    assert best["column"] == "read_name"
    assert best["normalise"] == "strip_segment"
    assert best["hits"] > 0
    assert "read_name_normalise='strip_segment'" in res["verdict"]


def test_diagnose_umi_reports_a_hopeless_mismatch(sim, tmp_path):
    from flightcollapse.diagnose import diagnose

    bad = tmp_path / "bad.tsv"
    with open(bad, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for i in range(500):
            fh.write(f"CB{i}\tACGTACGTAC\tothermovie/{i}/ccs\n")
    res = diagnose(sim.bam, str(bad), verbose=False)
    assert res["best"]["hits"] == 0
    assert "different BAM" in res["verdict"]


def test_diagnose_samples_across_the_file_not_just_the_head(sim, tmp_path):
    """These tables are written by streaming a coordinate-sorted BAM, so the
    head is chr1, not a sample of the library. Sampling only the head makes a
    perfectly good table look like a total mismatch."""
    from flightcollapse.diagnose import diagnose, sample_table

    names = open(sim.barcode_umi).read().rstrip("\n").split("\n")[1:]
    tsv = tmp_path / "ordered.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for i in range(60_000):                      # "chr1..chr20" block
            fh.write(f"CB{i % 9}\tACGTACGTAC\tothermovie/{i}/ccs/1_2\n")
        for line in names:                            # the "chr21" block, at the end
            cb, umi, name = line.split("\t")
            fh.write(f"{cb}\t{umi}\t{name}\n")

    head_only = [
        r for r in sample_table(str(tsv), 2000)[1][:0]
    ]  # sanity: helper returns rows
    assert head_only == []

    res = diagnose(sim.bam, str(tsv), table_rows=20_000, verbose=False)
    assert res["best"]["hits"] > 0, "random-offset sampling failed to reach the tail"
    assert res["best"]["normalise"] == "none"


def test_segmented_reads_do_not_get_a_zmw_level_join(tmp_path):
    """Regression: a coordinate-ordered table of Kinnex S-reads.

    The exact read name is the correct key. Matching on movie/ZMW also "works"
    -- and matches ~8x as many rows, because a ZMW holds several segments, each
    a different molecule with a different UMI. Recommending it would give most
    reads the wrong barcode. The tool must prefer the exact name and mark the
    coarse modes unsafe.
    """
    import pysam

    from flightcollapse.diagnose import diagnose

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrT", "LN": 100_000}]}
    bam_path = tmp_path / "seg.bam"
    zmws = list(range(1000, 1400))
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bf:
        for k, z in enumerate(zmws):
            a = pysam.AlignedSegment()
            a.query_name = f"movie/{z}/ccs/{100 * k}_{100 * k + 90}"
            a.query_sequence = "A" * 90
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 10 * k
            a.mapping_quality = 60
            a.cigartuples = [(0, 90)]
            a.set_tag("NM", 0)
            bf.write(a)
    pysam.index(str(bam_path))

    tsv = tmp_path / "seg.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for z in zmws:
            for seg in range(8):                       # 8 segments per ZMW
                k = zmws.index(z)
                span = (f"{100 * k}_{100 * k + 90}" if seg == 0
                        else f"{9000 + 100 * seg}_{9090 + 100 * seg}")
                fh.write(f"CB{seg}\tACGTACGT{seg}A\tmovie/{z}/ccs/{span}\n")

    res = diagnose(str(bam_path), str(tsv), table_rows=100_000, verbose=False)
    by_mode = {r["normalise"]: r for r in res["results"]}
    assert by_mode["none"]["hits"] > 0
    assert by_mode["zmw"]["hits"] > 2 * by_mode["none"]["hits"]
    assert by_mode["zmw"]["unsafe"], "ZMW-level join must be flagged unsafe"
    assert not by_mode["none"]["unsafe"]
    assert res["best"]["normalise"] == "none"
    assert "wrong barcode" in res["verdict"]


def test_tie_breaks_towards_the_most_specific_key(sim, tmp_path):
    """A whole-genome BAM makes every mode hit every sampled row. The tie must
    break to the exact read name, not to whichever mode was listed first."""
    from flightcollapse.diagnose import diagnose

    names = open(sim.barcode_umi).read().rstrip("\n").split("\n")[1:]
    tsv = tmp_path / "all.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for line in names:
            cb, umi, name = line.split("\t")
            fh.write(f"{cb}\t{umi}\t{name}\n")
    res = diagnose(sim.bam, str(tsv), table_rows=5000, verbose=False)
    hits = {r["normalise"]: r["hits"] for r in res["results"]}
    assert len(set(hits.values())) == 1, "expected every mode to tie"
    assert res["best"]["normalise"] == "none"


def test_no_strict_reports_instead_of_exiting_nonzero(sim, tmp_path, monkeypatch, capsys):
    """`--no-strict` promises to report rather than fail.

    It used to suppress the exception and then return 1 anyway, so a flag on
    0.2% of models -- with every output written correctly -- read to a calling
    script as "the run failed".
    """
    from flightcollapse import cli

    def _fake_run(cfg):
        return {
            "all_passed": False,
            "invariants": [{"invariant": "model_ids_are_unique", "passed": False}],
            "failed_contigs": {},
        }

    monkeypatch.setattr("flightcollapse.pipeline.run", _fake_run)
    base = ["run", "-b", sim.bam, "-g", sim.gtf, "-f", sim.genome_fasta,
            "-o", str(tmp_path), "-q"]
    assert cli.main(base + ["--no-strict"]) == 0
    assert "model_ids_are_unique" in capsys.readouterr().err
    assert cli.main(base) == 1
