"""flightcollapse -- reference-anchored long-read transcript-model collapse.

A replacement for ``isoseq collapse`` for oligo-dT single-cell Iso-Seq
(FLIGHT-seq) data, built around two rules the original violates:

* the representative model of a read group is the **maximal** intron chain,
  never the minimal one;
* the **empty** intron chain is not a suffix of anything, so mono-exonic reads
  can never absorb full-length spliced read mass.

Entry points::

    from flightcollapse import Config, run
    cfg = Config(bam="...", reference_gtf="...", genome_fasta="...")
    report = run(cfg)

or from the shell::

    flightcollapse run -b in.bam -g gencode.v44.annotation.gtf -f GRCh38.fa \\
                       -u sample_barcode_umi_strict.tsv -o out -p sample
"""

from .config import (
    AlignmentGates,
    ChainParams,
    Config,
    EndParams,
    JunctionParams,
    MoleculeParams,
    MonoexonParams,
    OutputParams,
    ScoringParams,
)
from .model import TranscriptModel

__version__ = "0.1.18"
__all__ = [
    "Config",
    "AlignmentGates",
    "JunctionParams",
    "ChainParams",
    "EndParams",
    "MonoexonParams",
    "ScoringParams",
    "MoleculeParams",
    "OutputParams",
    "TranscriptModel",
    "run",
    "__version__",
]


def run(cfg: "Config"):
    """Run the full pipeline (imported lazily so ``import flightcollapse`` is cheap)."""
    from .pipeline import run as _run

    return _run(cfg)
