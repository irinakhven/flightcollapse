"""Synthetic dataset with ground truth.

Builds a small genome, a GENCODE-shaped GTF, an aligned BAM and a barcode/UMI
table, with the specific structures the collapse has to get right:

* pervasive 5' truncation, deliberately at a level where the truncated form
  **outnumbers** the intact one (as it does in the real data) -- a tool that
  picks the representative by support or by minimality gets this wrong
* mono-exonic reads inside 3'UTRs with a proper polyA signal (real) and
  mono-exonic reads in the middle of a gene body with no signal (artefact)
* a genuine novel exon-skipping isoform with real support
* one-off spurious junctions (a single read each, non-canonical)
* an alternative 3' end with its own polyA signal
* PCR duplicates: many reads, few molecules, for one novel junction

Everything is deterministic given ``seed``.
"""

from __future__ import annotations

import os
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

BASES = "ACGT"


@dataclass
class SimTranscript:
    tx_id: str
    gene_id: str
    gene_name: str
    strand: str
    exons: List[Tuple[int, int]]
    tags: Tuple[str, ...] = ("basic",)

    @property
    def tss(self) -> int:
        return self.exons[0][0] if self.strand == "+" else self.exons[-1][1]

    @property
    def tts(self) -> int:
        return self.exons[-1][1] if self.strand == "+" else self.exons[0][0]

    def chain(self) -> Tuple[Tuple[int, int], ...]:
        ch = tuple((self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1))
        return ch[::-1] if self.strand == "-" else ch


@dataclass
class SimData:
    root: str
    contig: str
    genome_fasta: str
    gtf: str
    bam: str
    barcode_umi: str
    transcripts: Dict[str, SimTranscript] = field(default_factory=dict)
    truth: Dict[str, object] = field(default_factory=dict)


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _random_seq(n: int, rng: random.Random) -> List[str]:
    # C/G-leaning so that stray AATAAA hexamers and A-rich stretches are rare
    return rng.choices(BASES, weights=[0.2, 0.3, 0.3, 0.2], k=n)


def build(
    root: str,
    seed: int = 0,
    n_genes: int = 12,
    contig: str = "chr1",
    contig_len: int = 400_000,
    n_cells: int = 40,
) -> SimData:
    rng = random.Random(seed)
    os.makedirs(root, exist_ok=True)
    seq = _random_seq(contig_len, rng)

    transcripts: Dict[str, SimTranscript] = {}
    genes: List[Tuple[str, str, str]] = []   # gene_id, gene_name, strand
    layout: Dict[str, List[Tuple[int, int]]] = {}

    pos = 5_000
    for gi in range(1, n_genes + 1):
        strand = "+" if gi % 2 else "-"
        gene_id = f"SIMG{gi:04d}"
        gene_name = f"Sim{gi}"
        exons: List[Tuple[int, int]] = []
        p = pos
        for k in range(4):
            L = rng.randint(200, 500) if k < 3 else rng.randint(900, 1800)  # long last exon
            exons.append((p, p + L))
            p += L + rng.randint(1200, 2500)
        pos = p + 6_000
        layout[gene_id] = exons
        genes.append((gene_id, gene_name, strand))

        full = SimTranscript(
            f"ENSTSIM{gi:04d}1", gene_id, gene_name, strand, list(exons),
            tags=("basic", "MANE_Select"),
        )
        transcripts[full.tx_id] = full
        # a second annotated isoform for a third of the genes: skips exon 3
        if gi % 3 == 0:
            alt_ex = [exons[0], exons[1], exons[3]]
            alt = SimTranscript(
                f"ENSTSIM{gi:04d}2", gene_id, gene_name, strand, alt_ex, tags=("basic",)
            )
            transcripts[alt.tx_id] = alt

    # --- write real splice motifs and polyA signals into the sequence -----
    for t in transcripts.values():
        for i in range(len(t.exons) - 1):
            d, a = t.exons[i][1], t.exons[i + 1][0]
            if t.strand == "+":
                seq[d], seq[d + 1] = "G", "T"
                seq[a - 2], seq[a - 1] = "A", "G"
            else:
                seq[d], seq[d + 1] = "C", "T"
                seq[a - 2], seq[a - 1] = "A", "C"
        _place_polya(seq, t.tts, t.strand)

    # a real alternative polyA site inside the long last exon of gene 1
    g1 = transcripts["ENSTSIM00011"]
    apa_pos = g1.exons[-1][1] - 600 if g1.strand == "+" else g1.exons[0][0] + 600
    _place_polya(seq, apa_pos, g1.strand)

    genome_fasta = os.path.join(root, "sim.fa")
    with open(genome_fasta, "w") as fh:
        s = "".join(seq)
        fh.write(f">{contig}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")
    pysam.faidx(genome_fasta)

    gtf = os.path.join(root, "sim.gtf")
    _write_gtf(gtf, contig, transcripts, genes)

    bam = os.path.join(root, "sim.bam")
    bc = os.path.join(root, "sim_barcode_umi.tsv")
    truth = _write_reads(bam, bc, contig, contig_len, transcripts, apa_pos, rng, n_cells)

    return SimData(root, contig, genome_fasta, gtf, bam, bc, transcripts, truth)


def _place_polya(seq: List[str], tts: int, strand: str, offset: int = 20) -> None:
    motif = "AATAAA"
    if strand == "+":
        start = tts - offset - len(motif)
        if start < 0:
            return
        for k, ch in enumerate(motif):
            seq[start + k] = ch
    else:
        start = tts + offset
        if start + len(motif) >= len(seq):
            return
        rc = _revcomp(motif)
        for k, ch in enumerate(rc):
            seq[start + k] = ch


def _write_gtf(path: str, contig: str, transcripts, genes) -> None:
    rows: List[Tuple[int, str]] = []
    by_gene: Dict[str, List[SimTranscript]] = {}
    for t in transcripts.values():
        by_gene.setdefault(t.gene_id, []).append(t)
    for gene_id, gene_name, strand in genes:
        ts = by_gene[gene_id]
        gs = min(t.exons[0][0] for t in ts)
        ge = max(t.exons[-1][1] for t in ts)
        attr = f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "protein_coding";'
        rows.append((gs, f"{contig}\tSIM\tgene\t{gs + 1}\t{ge}\t.\t{strand}\t.\t{attr}"))
        for t in ts:
            a = (
                f'gene_id "{gene_id}"; transcript_id "{t.tx_id}"; gene_name "{gene_name}"; '
                f'transcript_type "protein_coding"; '
                + "".join(f'tag "{tg}"; ' for tg in t.tags)
            )
            rows.append(
                (t.exons[0][0],
                 f"{contig}\tSIM\ttranscript\t{t.exons[0][0] + 1}\t{t.exons[-1][1]}"
                 f"\t.\t{t.strand}\t.\t{a}")
            )
            for k, (s, e) in enumerate(t.exons, 1):
                rows.append(
                    (s, f"{contig}\tSIM\texon\t{s + 1}\t{e}\t.\t{t.strand}\t.\t"
                        f'{a}exon_number "{k}";')
                )
            # a CDS that stops 400 nt before the 3' end, so the 3'UTR is defined
            cds = _cds_exons(t, 400)
            for s, e in cds:
                rows.append((s, f"{contig}\tSIM\tCDS\t{s + 1}\t{e}\t.\t{t.strand}\t0\t{a}"))
            if cds:
                sp = (cds[-1][1] - 3, cds[-1][1]) if t.strand == "+" else (cds[0][0], cds[0][0] + 3)
                rows.append(
                    (sp[0], f"{contig}\tSIM\tstop_codon\t{sp[0] + 1}\t{sp[1]}"
                            f"\t.\t{t.strand}\t0\t{a}")
                )
    rows.sort(key=lambda r: r[0])
    with open(path, "w") as fh:
        fh.write("##description: simulated annotation\n")
        for _p, line in rows:
            fh.write(line + "\n")


def _cds_exons(t: SimTranscript, utr3_len: int) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    remaining = utr3_len
    ex = t.exons if t.strand == "+" else t.exons[::-1]
    trimmed = []
    for s, e in reversed(ex):
        L = e - s
        if remaining >= L:
            remaining -= L
            continue
        if remaining > 0:
            if t.strand == "+":
                trimmed.append((s, e - remaining))
            else:
                trimmed.append((s + remaining, e))
            remaining = 0
        else:
            trimmed.append((s, e))
    out = sorted(trimmed)
    return [x for x in out if x[1] > x[0]]


# ---------------------------------------------------------------------- #
def _cigar_for(exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    cig: List[Tuple[int, int]] = []
    for i, (s, e) in enumerate(exons):
        if i:
            cig.append((3, s - exons[i - 1][1]))
        cig.append((0, e - s))
    return cig


def _emit(bam, contig_tid: int, name: str, exons, strand: str) -> None:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "A" * sum(e - s for s, e in exons)
    a.flag = 16 if strand == "-" else 0
    a.reference_id = contig_tid
    a.reference_start = exons[0][0]
    a.mapping_quality = 60
    a.cigartuples = _cigar_for(exons)
    a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
    a.set_tag("NM", 0)
    bam.write(a)


def _truncate(t: SimTranscript, k: int) -> List[Tuple[int, int]]:
    """Drop ``k`` exons from the 5' end (in transcript orientation)."""
    ex = list(t.exons)
    if t.strand == "+":
        return ex[k:]
    return ex[: len(ex) - k]


def _clip_3p(exons, strand: str, new_tts: int):
    out = []
    for s, e in exons:
        if strand == "+":
            if s >= new_tts:
                continue
            out.append((s, min(e, new_tts)))
        else:
            if e <= new_tts:
                continue
            out.append((max(s, new_tts), e))
    return [x for x in out if x[1] > x[0]]


def _write_reads(
    bam_path: str,
    bc_path: str,
    contig: str,
    contig_len: int,
    transcripts: Dict[str, SimTranscript],
    apa_pos: int,
    rng: random.Random,
    n_cells: int,
) -> Dict[str, object]:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": contig, "LN": contig_len}]}
    records: List[Tuple[int, str, List[Tuple[int, int]], str]] = []
    rows: List[Tuple[str, str, str]] = []
    truth: Dict[str, object] = {
        "full_length_reads": 0,
        "truncated_reads": 0,
        "monoexon_utr3_reads": 0,
        "monoexon_internal_reads": 0,
        "spurious_junction_reads": 0,
        "novel_skip_reads": 0,
        "apa_reads": 0,
        "expect_models": [],
        "expect_absent": [],
    }
    counter = [0]
    cells = [f"CELL{i:03d}" for i in range(n_cells)]

    def add(exons, strand, tag, molecule: Optional[str] = None) -> None:
        counter[0] += 1
        name = f"read/{counter[0]}/ccs"
        records.append((exons[0][0], name, list(exons), strand))
        cb = molecule.split(":")[0] if molecule else rng.choice(cells)
        umi = (
            molecule.split(":")[1] if molecule
            else "".join(rng.choices(BASES, k=10))
        )
        rows.append((cb, umi, name))

    tx_list = sorted(transcripts.values(), key=lambda t: t.tx_id)
    for t in tx_list:
        gi = int(t.tx_id[7:11])
        is_primary = t.tx_id.endswith("1")

        # ---- full-length + heavily 5'-truncated ------------------------
        n_full = 25 if is_primary else 12
        for _ in range(n_full):
            add(t.exons, t.strand, "full")
            truth["full_length_reads"] += 1
        if is_primary:
            # the truncated form OUTNUMBERS the intact one, as in real data
            for _ in range(70):
                ex = _truncate(t, rng.choice([1, 2]))
                if len(ex) >= 2:
                    add(ex, t.strand, "trunc")
                    truth["truncated_reads"] += 1
            truth["expect_models"].append(t.tx_id)

    # ---- mono-exonic 3'UTR fragments of gene 1 (real: polyA capture) ----
    g1 = transcripts["ENSTSIM00011"]
    last = g1.exons[-1] if g1.strand == "+" else g1.exons[0]
    for _ in range(40):
        if g1.strand == "+":
            s = rng.randint(last[0] + 100, last[1] - 300)
            ex = [(s, last[1])]
        else:
            e = rng.randint(last[0] + 300, last[1] - 100)
            ex = [(last[0], e)]
        add(ex, g1.strand, "mono_utr3")
        truth["monoexon_utr3_reads"] += 1

    # ---- mono-exonic fragments in the middle of gene 3, no polyA signal --
    g3 = transcripts["ENSTSIM00031"]
    mid = g3.exons[1]
    for _ in range(30):
        s = rng.randint(mid[0], max(mid[0] + 1, mid[1] - 120))
        add([(s, min(mid[1], s + 120))], g3.strand, "mono_internal")
        truth["monoexon_internal_reads"] += 1
    truth["expect_absent"].append("monoexon_internal_gene3")

    # ---- a genuine novel exon-skipping isoform of gene 2 -----------------
    g2 = transcripts["ENSTSIM00021"]
    skip = [g2.exons[0], g2.exons[1], g2.exons[3]]
    for _ in range(35):
        add(skip, g2.strand, "novel_skip")
        truth["novel_skip_reads"] += 1
    truth["expect_models"].append("novel_skip_gene2")

    # ---- one-off spurious junctions (singletons, non-canonical) ---------
    for gi, t in enumerate(tx_list[:6]):
        ex = list(t.exons)
        if len(ex) < 3:
            continue
        bad = [ex[0], (ex[1][0] + 37, ex[1][1]), ex[2]]
        add(bad, t.strand, "spurious")
        truth["spurious_junction_reads"] += 1
    truth["expect_absent"].append("spurious_singletons")

    # ---- PCR-duplicate novel junction: many reads, one molecule ---------
    t = tx_list[4]
    ex = list(t.exons)
    if len(ex) >= 3:
        dup = [ex[0], (ex[1][0] + 61, ex[1][1]), ex[2]]
        for _ in range(30):
            add(dup, t.strand, "pcr_dup", molecule="CELL000:AAAACCCCGG")
        truth["expect_absent"].append("pcr_duplicate_junction")

    # ---- alternative 3' end on gene 1 (real APA, own polyA signal) ------
    for _ in range(30):
        ex = _clip_3p(g1.exons, g1.strand, apa_pos)
        if len(ex) >= 2:
            add(ex, g1.strand, "apa")
            truth["apa_reads"] += 1
    truth["expect_models"].append("ENSTSIM00011|3UTR_mod1")

    records.sort(key=lambda r: r[0])
    tmp = bam_path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as bf:
        for _s, name, exons, strand in records:
            _emit(bf, 0, name, exons, strand)
    pysam.sort("-o", bam_path, tmp)
    pysam.index(bam_path)
    os.remove(tmp)

    with open(bc_path, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for cb, umi, name in rows:
            fh.write(f"{cb}\t{umi}\t{name}\n")

    truth["total_reads"] = len(records)
    return truth
