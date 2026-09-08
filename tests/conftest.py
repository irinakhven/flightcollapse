import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import simulate  # noqa: E402

from flightcollapse import Config  # noqa: E402
from flightcollapse.pipeline import run  # noqa: E402


@pytest.fixture(scope="session")
def sim(tmp_path_factory):
    root = tmp_path_factory.mktemp("sim")
    return simulate.build(str(root), seed=0)


@pytest.fixture(scope="session")
def collapsed(sim, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    cfg = Config(
        bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta
    )
    cfg.molecules.barcode_umi_tsv = sim.barcode_umi
    cfg.output.outdir = str(out)
    cfg.output.prefix = "sim"
    cfg.verbose = False
    cfg.strict_invariants = False
    report = run(cfg)
    return {"cfg": cfg, "report": report, "outdir": str(out), "sim": sim}


@pytest.fixture(scope="session")
def models_table(collapsed):
    import pandas as pd

    return pd.read_csv(os.path.join(collapsed["outdir"], "sim.models.tsv"), sep="\t")
