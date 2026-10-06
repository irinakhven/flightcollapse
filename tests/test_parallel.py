"""Parallel and serial runs must produce the same files, byte for byte.

This is the acceptance test for 0.5.0's parallelism, and it is written as an
equality assertion rather than a smoke test on purpose.  A parallel collapse
that is merely "close" to the serial one is useless: every downstream
comparison in this project -- three-sample reproducibility, SIRV precision,
arm-versus-arm catalogues -- assumes that two runs differing only in machine
settings are the same run.  If ``-j 4`` changed a single model id, none of
those comparisons could be trusted again.

The things most likely to break are all exercised here: model ids and the
``PB.N`` gene numbering (assigned globally, in contig order), the read_stat and
group files (stitched from per-contig shards), and the molecule matrix
(accumulated from count shards).
"""

from __future__ import annotations

import os
import shutil
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import simulate  # noqa: E402

from flightcollapse import Config  # noqa: E402
from flightcollapse.pipeline import run  # noqa: E402


def _merge_sims(root: str, sims):
    """One BAM / GTF / FASTA over several contigs, from several simulations.

    Read names are prefixed with the contig: the simulator numbers reads from 1
    in every build, so an unprefixed merge would give two different molecules
    the same name and the barcode table would silently resolve both to whichever
    row it saw last.
    """
    os.makedirs(root, exist_ok=True)
    fa = os.path.join(root, "merged.fa")
    with open(fa, "w") as out:
        for s in sims:
            with open(s.genome_fasta) as fh:
                shutil.copyfileobj(fh, out)
    pysam.faidx(fa)

    gtf = os.path.join(root, "merged.gtf")
    with open(gtf, "w") as out:
        for s in sims:
            with open(s.gtf) as fh:
                shutil.copyfileobj(fh, out)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            dict(SN=s.contig, LN=pysam.FastaFile(s.genome_fasta).get_reference_length(s.contig))
            for s in sims
        ],
    }
    bam = os.path.join(root, "merged.bam")
    unsorted = bam + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for s in sims:
            with pysam.AlignmentFile(s.bam, "rb") as src:
                for rec in src.fetch(until_eof=True):
                    # Rebuilt against the merged header rather than re-pointed:
                    # reference_id is an index into the header the record came
                    # from, and assigning across headers is what pysam refuses.
                    d = rec.to_dict()
                    d["name"] = f"{s.contig}|{d['name']}"
                    d["ref_name"] = s.contig
                    out.write(pysam.AlignedSegment.from_dict(d, out.header))
    pysam.sort("-o", bam, unsorted)
    pysam.index(bam)
    os.remove(unsorted)

    bc = os.path.join(root, "merged.barcode_umi.tsv")
    with open(bc, "w") as out:
        out.write("cell_barcode\tumi\tread_name\n")
        for s in sims:
            with open(s.barcode_umi) as fh:
                next(fh)
                for line in fh:
                    cb, umi, name = line.rstrip("\n").split("\t")
                    out.write(f"{cb}\t{umi}\t{s.contig}|{name}\n")
    return fa, gtf, bam, bc


@pytest.fixture(scope="module")
def multi(tmp_path_factory):
    root = tmp_path_factory.mktemp("multi")
    sims = [
        simulate.build(str(root / f"s{i}"), seed=i, contig=c, n_genes=6)
        for i, c in enumerate(("chr1", "chr2", "chr3", "chr4"))
    ]
    fa, gtf, bam, bc = _merge_sims(str(root / "merged"), sims)
    return {"fa": fa, "gtf": gtf, "bam": bam, "bc": bc,
            "contigs": [s.contig for s in sims]}


def _run(multi, outdir: str, workers: int) -> dict:
    cfg = Config(bam=multi["bam"], reference_gtf=multi["gtf"],
                 genome_fasta=multi["fa"])
    cfg.molecules.barcode_umi_tsv = multi["bc"]
    cfg.output.outdir, cfg.output.prefix = outdir, "m"
    cfg.verbose, cfg.strict_invariants = False, False
    cfg.parallel.workers = workers
    # The clamp is about the machine, not about correctness, and a CI box with
    # little free memory would silently turn the parallel arm into a serial one
    # -- which would make this test pass by not testing anything.
    cfg.parallel.clamp_by_memory = False
    return {"report": run(cfg), "outdir": outdir, "cfg": cfg}


#: Every file the run writes that must not depend on how many workers ran.
COMPARED = (
    ".gtf", ".gff", ".models.tsv", ".abundance.txt", ".flnc_count.txt",
    ".read_stat.txt", ".group.txt", ".group_chains.tsv", ".secondpass.tsv",
    "_umi_corrected_isoform_matrix.mtx", "_cells.txt", "_isoforms.txt",
    "_isoform_metadata.txt",
)


@pytest.fixture(scope="module")
def two_runs(multi, tmp_path_factory):
    a = str(tmp_path_factory.mktemp("serial"))
    b = str(tmp_path_factory.mktemp("parallel"))
    return _run(multi, a, workers=1), _run(multi, b, workers=4)


def test_the_parallel_run_used_more_than_one_worker(two_runs):
    """Guard the guard: a clamped-to-1 parallel arm would make every other
    assertion in this file vacuously true."""
    _serial, par = two_runs
    assert par["report"]["parallel"]["workers"] == 4
    assert par["report"]["parallel"]["start_method"] == "fork"
    assert _serial["report"]["parallel"]["workers"] == 1


def _read(outdir: str, suffix: str) -> bytes:
    """Contents of an output, transparently un-gzipping the big tables.

    Compared decompressed rather than byte-for-byte on disk because gzip
    records a modification time in its header, so two identical tables written
    a second apart differ as files while being the same table.
    """
    path = os.path.join(outdir, "m" + suffix)
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return fh.read()
    if os.path.exists(path + ".gz"):
        import gzip

        with gzip.open(path + ".gz", "rb") as fh:
            return fh.read()
    raise AssertionError(f"neither {path} nor {path}.gz exists")


@pytest.mark.parametrize("suffix", COMPARED)
def test_outputs_are_identical(two_runs, suffix):
    ser, par = two_runs
    a, b = _read(ser["outdir"], suffix), _read(par["outdir"], suffix)
    assert a == b, f"{suffix} differs between a 1-worker and a 4-worker run"


def test_the_qc_verdict_is_the_same(two_runs):
    """The invariants are accumulated per contig and merged; merging must be
    exact, not approximately exact."""
    ser, par = two_runs
    for key in ("n_models", "read_accounting", "invariants", "all_passed",
                "novel_model_burden", "model_id_collisions",
                "read_stat_rows", "reads_unassigned"):
        assert ser["report"].get(key) == par["report"].get(key), key
    # second_pass carries the audit table's path, which is the one thing that
    # legitimately differs: the two runs write to different directories.
    sa = dict(ser["report"]["second_pass"])
    sb = dict(par["report"]["second_pass"])
    sa.pop("audit_table"), sb.pop("audit_table")
    assert sa == sb


def test_every_contig_was_actually_processed(two_runs, multi):
    ser, par = two_runs
    for r in (ser, par):
        assert not r["report"]["failed_contigs"], r["report"]["failed_contigs"]
        assert set(r["report"]["per_contig"]) == set(multi["contigs"])


def test_shards_are_cleaned_up(two_runs):
    """A run that leaves its scratch behind fills the disk on the fourth sample."""
    for r in two_runs:
        stray = [p for p in os.listdir(r["outdir"]) if p.startswith(".m_shards")]
        assert not stray, stray


def test_a_failing_contig_does_not_take_the_run_with_it(multi, tmp_path_factory,
                                                        monkeypatch):
    """The serial path has always finished the other contigs after one failed.
    Moving that across a process boundary is exactly where it would get lost."""
    from flightcollapse import pipeline as pipemod

    real = pipemod._run_contig

    def boom(shared, contig):
        if contig == "chr2":
            raise RuntimeError("synthetic failure")
        return real(shared, contig)

    monkeypatch.setattr(pipemod, "_run_contig", boom)
    out = str(tmp_path_factory.mktemp("failing"))
    # serial: the monkeypatch is visible without needing it to survive a fork
    res = _run(multi, out, workers=1)
    rep = res["report"]
    assert "chr2" in rep["failed_contigs"]
    assert "synthetic failure" in rep["failed_contigs"]["chr2"]
    assert rep["all_passed"] is False
    assert set(rep["per_contig"]) == {"chr1", "chr3", "chr4"}
    assert rep["n_models"] > 0, "the surviving contigs were still written"


# --------------------------------------------------------------------- #
# A killed worker must not hang the run.
#
# BD144_kin sat at 23 of 24 contigs for two and a half hours with an idle pool
# and a silent log, because one worker had been SIGKILLed and
# ``Pool.imap_unordered`` was still waiting for a result that no longer had a
# process behind it.  This is that failure, reproduced deterministically.

_PARENT_PID = None


def _suicidal_worker(contig: str):
    """Dies on one contig, but only in a child -- so the parent's retry works.

    Standing in for the OOM killer: the point is not how the process died but
    that it died without putting anything on the result queue.
    """
    import signal

    from flightcollapse.parallel import ContigResult

    if contig == "chrKill" and os.getpid() != _PARENT_PID:
        os.kill(os.getpid(), signal.SIGKILL)
    return ContigResult(contig=contig, ok=True, seconds=0.5)


@pytest.mark.skipif("fork" not in __import__("multiprocessing").get_all_start_methods(),
                    reason="fork-only")
def test_a_killed_worker_is_noticed_and_its_contig_recomputed():
    from flightcollapse import parallel as pmod

    global _PARENT_PID
    _PARENT_PID = os.getpid()
    contigs = ["chrA", "chrB", "chrKill", "chrC"]
    lines = []

    out = pmod.map_contigs(_suicidal_worker, contigs, workers=3,
                           log=lines.append, poll=0.1, grace=0.5)

    assert set(out) == set(contigs), "the run must not lose a contig"
    assert all(r.ok for r in out.values())
    assert any("disappeared without returning a result" in ln for ln in lines)
    assert any("re-running serially in the parent" in ln for ln in lines)
