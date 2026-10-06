"""A fixture with soft-clipped polyA tails, for the 0.5.0 second pass.

``simulate.py`` emits reads with no soft clip at all, which is fine for
everything it was written to test but makes the templated-tail measurement
vacuous: with no clip every read's terminal A-run is 0, so a genuine cleavage
site and an internal-priming site look identical.  Under oligo-dT the tail is
the physical thing that was captured, so a fixture without one cannot exercise
the discriminator this release is built on.

Six loci, each a specific verdict the second pass has to reach:

``MonoOnly``    gene with one single-exon transcript -> class A, kept
``Both``        gene with a spliced and a single-exon transcript -> class B, kept
``MultiData``   multi-exon gene; a mono pile spanning intron 1 -> class C,
                replaced by the gene's own best-supported spliced model
``MultiRef``    multi-exon gene with NO spliced reads; a mono pile inside one
                exon -> class D, replaced by the reference canonical, INFERRED
``Intronic``    multi-exon gene (minus strand); a mono pile wholly inside an
                intron -> class E, replaced
``Apa``         a genuine alternative 3' end: no genomic A-run under it, reads
                carrying a 25 nt NON-templated tail -> end3_novel, kept
``IntPrime``    an internal-priming decoy: a 12 nt genomic A-run, reads whose
                12 nt "tail" is exactly that run -> end3_novel in the first
                pass (the 0.1.18 tail rescue fires, because an A-run is an
                A-run), removed by the second pass

``Apa`` and ``IntPrime`` are the pair that matters.  They are indistinguishable
to every 0.4.6 measurement -- same read counts, same tail fraction, same tight
peak -- and differ only in what the genome says underneath the tail.
"""

from __future__ import annotations

import os
import random
from dataclasses import dataclass
from typing import Dict, List, Tuple

import pysam

BASES = "ACGT"
CONTIG = "chr1"
CONTIG_LEN = 80_000


@dataclass
class TailSim:
    root: str
    contig: str
    genome_fasta: str
    gtf: str
    bam: str
    barcode_umi: str
    expect: Dict[str, str]


def _patch(seq: List[str], pos: int, text: str) -> None:
    for k, ch in enumerate(text):
        if 0 <= pos + k < len(seq):
            seq[pos + k] = ch


def build(root: str, seed: int = 7) -> TailSim:
    rng = random.Random(seed)
    os.makedirs(root, exist_ok=True)
    seq = [rng.choice(BASES) for _ in range(CONTIG_LEN)]

    # transcript_id -> (gene_id, gene_name, strand, exons, tags)
    tx: Dict[str, Tuple[str, str, str, List[Tuple[int, int]], Tuple[str, ...]]] = {}

    tx["ENSTMONO1"] = ("GMONO", "MonoOnly", "+", [(5_000, 6_500)],
                       ("basic", "Ensembl_canonical"))
    tx["ENSTBOTH1"] = ("GBOTH", "Both", "+",
                       [(12_000, 12_400), (13_000, 13_400), (14_000, 15_000)],
                       ("basic", "MANE_Select"))
    tx["ENSTBOTH2"] = ("GBOTH", "Both", "+", [(14_000, 15_000)], ("basic",))
    tx["ENSTMD1"] = ("GMD", "MultiData", "+",
                     [(20_000, 20_400), (21_000, 21_400), (22_000, 23_000)],
                     ("basic", "MANE_Select"))
    tx["ENSTMR1"] = ("GMR", "MultiRef", "+",
                     [(30_000, 30_400), (31_000, 31_400), (32_000, 33_000)],
                     ("basic", "MANE_Select"))
    tx["ENSTMR2"] = ("GMR", "MultiRef", "+",
                     [(30_000, 30_400), (32_000, 33_000)], ("basic",))
    tx["ENSTINT1"] = ("GINT", "Intronic", "-",
                      [(36_000, 37_000), (38_000, 38_400), (39_000, 39_400)],
                      ("basic", "MANE_Select"))
    tx["ENSTAPA1"] = ("GAPA", "Apa", "+",
                      [(40_000, 40_400), (41_000, 41_400), (42_000, 44_000)],
                      ("basic", "MANE_Select"))
    tx["ENSTIP1"] = ("GIP", "IntPrime", "+",
                     [(50_000, 50_400), (51_000, 51_400), (52_000, 54_000)],
                     ("basic", "MANE_Select"))

    # -- sequence patches ------------------------------------------------
    # A polyA signal ~20 bp upstream of every annotated 3' end, as in real
    # 3'UTRs, so the first pass has its usual genomic evidence.
    for tid, (_g, _n, strand, exons, _t) in tx.items():
        tts = exons[-1][1] if strand == "+" else exons[0][0]
        if strand == "+":
            _patch(seq, tts - 26, "AATAAA")
        else:
            _patch(seq, tts + 20, "TTTATT")  # revcomp of AATAAA

    # Apa: a real alternative cleavage site at 43,000. Nothing A-rich under it
    # -- a genuine tail is non-templated, so the genome must NOT supply one.
    _patch(seq, 42_974, "AATAAA")
    _patch(seq, 43_000, "GTGTCTGTTCGTTGCTGTTCGGTCTGTTCG")

    # IntPrime: a 12 nt genomic A-run at 53,000 and a hexamer-free upstream
    # region, so the ONLY thing that makes the first pass call a 3' end here is
    # the reads' own A-run -- which is the genome's, not the molecule's.
    _patch(seq, 52_950, "GCGCGCTGCTGCGGCTGCGCTGCGGCTGCTGCGGCTGCGCTGCGGCTGCT")
    _patch(seq, 53_000, "A" * 12)
    _patch(seq, 53_012, "CAGCTGCAGCTGCAGCTGCAGCTGCAGCTG")

    fasta = os.path.join(root, "genome.fa")
    with open(fasta, "w") as fh:
        s = "".join(seq)
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")
    pysam.faidx(fasta)

    # -- GTF --------------------------------------------------------------
    gtf = os.path.join(root, "ref.gtf")
    rows: List[Tuple[int, str]] = []
    genes: Dict[str, Tuple[str, str, int, int]] = {}
    for tid, (gid, gname, strand, exons, tags) in tx.items():
        lo, hi = exons[0][0], exons[-1][1]
        if gid in genes:
            g = genes[gid]
            genes[gid] = (gname, strand, min(g[2], lo), max(g[3], hi))
        else:
            genes[gid] = (gname, strand, lo, hi)
        attr = (f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}"; '
                f'transcript_type "protein_coding"; '
                + " ".join(f'tag "{t}";' for t in tags))
        rows.append((lo, f"{CONTIG}\tsim\ttranscript\t{lo + 1}\t{hi}\t.\t{strand}\t.\t{attr}"))
        for s, e in exons:
            rows.append((s, f"{CONTIG}\tsim\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{attr}"))
    for gid, (gname, strand, lo, hi) in genes.items():
        attr = f'gene_id "{gid}"; gene_name "{gname}";'
        rows.append((lo, f"{CONTIG}\tsim\tgene\t{lo + 1}\t{hi}\t.\t{strand}\t.\t{attr}"))
    rows.sort(key=lambda r: r[0])
    with open(gtf, "w") as fh:
        for _p, line in rows:
            fh.write(line + "\n")

    # -- BAM ---------------------------------------------------------------
    records: List[Tuple[int, str, List[Tuple[int, int]], str, int]] = []
    bc_rows: List[Tuple[str, str, str]] = []
    counter = [0]

    def add(exons, strand: str, tail: int = 0, n: int = 1) -> None:
        for _ in range(n):
            counter[0] += 1
            name = f"read/{counter[0]}/ccs"
            records.append((exons[0][0], name, list(exons), strand, tail))
            bc_rows.append((f"CELL{counter[0] % 40:03d}",
                            "".join(rng.choices(BASES, k=10)), name))

    T = 25  # a realistic non-templated tail length for this library

    add([(5_100, 6_500)], "+", tail=T, n=20)                       # A
    # B carries no spliced reads on purpose. With one, the mono track's
    # containment rule would demote the pile into the spliced model before the
    # second pass ever saw it -- correctly, but then the class-B branch would
    # be untested, which is the opposite of what this locus is for.
    add([(14_200, 15_000)], "+", tail=T, n=20)                     # B
    add([(20_000, 20_400), (21_000, 21_400), (22_000, 23_000)], "+", tail=T, n=30)
    add([(20_300, 21_100)], "+", tail=T, n=20)                     # C: spans intron 1
    add([(32_200, 33_000)], "+", tail=T, n=20)                     # D: inside exon 3
    add([(37_300, 37_800)], "-", tail=T, n=20)                     # E: inside an intron
    add([(40_000, 40_400), (41_000, 41_400), (42_000, 44_000)], "+", tail=T, n=40)
    add([(40_000, 40_400), (41_000, 41_400), (42_000, 43_000)], "+", tail=T, n=25)
    add([(50_000, 50_400), (51_000, 51_400), (52_000, 54_000)], "+", tail=T, n=40)
    # The decoy: the alignment stops where the aligner ran out of A's to consume
    # and the "tail" it leaves behind is the remainder of the genomic run.
    add([(50_000, 50_400), (51_000, 51_400), (52_000, 53_000)], "+", tail=12, n=25)

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    bam = os.path.join(root, "reads.bam")
    unsorted = bam + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bf:
        for _s, name, exons, strand, tail in sorted(records, key=lambda r: r[0]):
            bf.write(_segment(name, exons, strand, tail, rng))
    pysam.sort("-o", bam, unsorted)
    pysam.index(bam)
    os.remove(unsorted)

    bcp = os.path.join(root, "barcode_umi.tsv")
    with open(bcp, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for cb, umi, name in bc_rows:
            fh.write(f"{cb}\t{umi}\t{name}\n")

    return TailSim(
        root=root, contig=CONTIG, genome_fasta=fasta, gtf=gtf, bam=bam,
        barcode_umi=bcp,
        expect={
            "GMONO": "A", "GBOTH": "B", "GMD": "C", "GMR": "D", "GINT": "E",
        },
    )


def _segment(name: str, exons, strand: str, tail: int, rng) -> pysam.AlignedSegment:
    """One alignment, with the polyA tail as a 3' soft clip.

    Orientation is the part worth being careful about: the stored read is in
    reference orientation, so on a minus-strand alignment the transcript's 3'
    end is at the START of the sequence and the tail reads as poly-T.  That is
    the convention ``reads.polya_tail_stats`` decodes, and a fixture that got it
    backwards would make the minus-strand half of the test meaningless.
    """
    a = pysam.AlignedSegment()
    a.query_name = name
    aligned = sum(e - s for s, e in exons)
    # a few non-A bases after the tail, standing in for the untrimmed UMI and
    # barcode this library carries between the polyA and the 3' primer
    handle = "".join(rng.choices("CGT", k=10)) if tail else ""
    body = "A" * aligned
    cig: List[Tuple[int, int]] = []
    for i, (s, e) in enumerate(exons):
        if i:
            cig.append((3, s - exons[i - 1][1]))
        cig.append((0, e - s))
    clip = tail + len(handle)
    if strand == "-":
        a.flag = 16
        a.query_sequence = ("T" * tail + handle)[::-1] + body if clip else body
        if clip:
            cig = [(4, clip)] + cig
    else:
        a.flag = 0
        a.query_sequence = body + "A" * tail + handle
        if clip:
            cig = cig + [(4, clip)]
    a.reference_id = 0
    a.reference_start = exons[0][0]
    a.mapping_quality = 60
    a.cigartuples = cig
    a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
    a.set_tag("NM", 0)
    return a
