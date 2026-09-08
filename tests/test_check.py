"""The standalone validation harness.

Two directions matter: it must pass this tool's own output, and it must *fail*
a deliberately isoseq-shaped output where the representative model is the
minimal member of its read group.  A harness that only ever passes is worthless.
"""

import os

import pysam
import pytest

from flightcollapse.check import inspect, parse_models


def test_check_passes_our_own_output(collapsed):
    out = collapsed["outdir"]
    rs = os.path.join(out, "sim.read_stat.txt.gz")
    if not os.path.exists(rs):
        rs = os.path.join(out, "sim.read_stat.txt")
    res = inspect(
        collapsed["sim"].bam, os.path.join(out, "sim.gtf"), rs,
        cigar_scan_reads=None, verbose=False,
    )
    s = res["summary"]
    assert not s["step0_mismatch"]
    assert s["read_ids_in_multiple_models"] == 0
    assert s["read_stat_rows"] == s["read_stat_unique_ids"]
    multi = s["by_class"].get("multi-exon")
    assert multi is not None
    assert multi["frac_model_is_min"] < 0.95
    assert abs(multi["median_log2_span_ratio"]) < 1.0


@pytest.fixture(scope="module")
def broken(sim, tmp_path_factory):
    """An isoseq-shaped failure: every model is a mono-exonic stub inside the
    last exon of its gene, and every read of the gene is assigned to it."""
    root = tmp_path_factory.mktemp("broken")
    gtf = root / "broken.gtf"
    rst = root / "broken.read_stat.txt"

    genes = {}
    for t in sim.transcripts.values():
        if not t.tx_id.endswith("1"):
            continue
        last = t.exons[-1] if t.strand == "+" else t.exons[0]
        stub = (last[1] - 300, last[1]) if t.strand == "+" else (last[0], last[0] + 300)
        genes[t.gene_id] = (t, stub)

    with open(gtf, "w") as fh:
        for i, (gid, (t, stub)) in enumerate(sorted(genes.items()), 1):
            a = f'gene_id "{gid}"; transcript_id "PB.{i}.1";'
            fh.write(
                f"chr1\tBAD\ttranscript\t{stub[0] + 1}\t{stub[1]}\t.\t{t.strand}\t.\t{a}\n"
            )
            fh.write(f"chr1\tBAD\texon\t{stub[0] + 1}\t{stub[1]}\t.\t{t.strand}\t.\t{a}\n")

    bam = pysam.AlignmentFile(sim.bam)
    with open(rst, "w") as fh:
        fh.write("id\tlength\tis_fl\tstat\tpbid\n")
        for i, (gid, (t, _stub)) in enumerate(sorted(genes.items()), 1):
            lo = min(e[0] for e in t.exons)
            hi = max(e[1] for e in t.exons)
            seen = set()
            for r in bam.fetch("chr1", lo, hi):
                if r.is_unmapped or r.query_name in seen:
                    continue
                seen.add(r.query_name)
                fh.write(
                    f"{r.query_name}\t{r.reference_length}\tY\tunique\tPB.{i}.1\n"
                )
    bam.close()
    return {"gtf": str(gtf), "read_stat": str(rst)}


def test_check_detects_the_isoseq_failure_mode(sim, broken):
    res = inspect(sim.bam, broken["gtf"], broken["read_stat"],
                  cigar_scan_reads=None, verbose=False)
    s = res["summary"]
    assert s["step0_mismatch"], s
    assert s["frac_read_mass_on_monoexon_models"] > 0.5
    assert s["frac_reads_unspliced"] < 0.15
    mono = s["by_class"]["mono-exon"]
    assert mono["median_read_exons"] > 1        # the reads are spliced
    assert mono["frac_model_is_min"] == 1.0     # the model is the shortest member
    assert mono["median_log2_span_ratio"] > 1.0  # reads are far longer than the model


def test_parse_models_reads_gtf_exons(collapsed):
    m = parse_models(os.path.join(collapsed["outdir"], "sim.gtf"))
    assert m
    for v in m.values():
        assert v["exons"] == sorted(v["exons"])
        assert v["n_exons"] >= 1
        assert v["strand"] in "+-"
