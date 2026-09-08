"""Reference annotation index.

Parses a GENCODE/Ensembl GTF once into the structures the collapse needs:

* ``chains[(contig, strand)][chain] -> [tx_id, ...]``  exact intron chains
* ``donors / acceptors / junctions`` per ``(contig, strand)``  -- the annotated
  splice-site catalogue that read junctions are snapped onto **first**
* per-transcript exon structure, TSS/TTS, stop codon, 3'UTR length, tags
* per-gene: transcript list, annotated 3' ends, terminal exons, exon union,
  and a binned positional index for "which gene is this read in?"

Everything is 0-based half-open (see :mod:`flightcollapse.intervals`).
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from .intervals import (
    Chain,
    Exon,
    exonic_bases_downstream,
    exonic_length,
    exons_to_chain,
    gtf_to_exon,
    total_overlap,
    tss_of,
    tts_of,
)

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')
BIN = 100_000


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_attrs(field_str: str) -> Dict[str, str]:
    return dict(_ATTR.findall(field_str))


@dataclass
class Transcript:
    tx_id: str
    gene_id: str
    gene_name: str
    contig: str
    strand: str
    tx_type: str = ""
    tags: Tuple[str, ...] = ()
    exons: List[Exon] = field(default_factory=list)
    cds: List[Exon] = field(default_factory=list)
    stop: List[Exon] = field(default_factory=list)
    # filled by finalise()
    chain: Chain = ()
    tss: int = 0
    tts: int = 0
    tx_len: int = 0
    n_exons: int = 0
    stop_pos: Optional[int] = None
    start_pos: Optional[int] = None
    utr3_len: int = 0
    rank: int = 99  # lower = more canonical

    def finalise(self) -> None:
        self.exons.sort()
        self.cds.sort()
        self.stop.sort()
        self.n_exons = len(self.exons)
        self.tx_len = exonic_length(self.exons)
        self.tss = tss_of(self.exons, self.strand)
        self.tts = tts_of(self.exons, self.strand)
        self.chain = exons_to_chain(self.exons, self.strand)
        self.stop_pos = _stop_position(self)
        self.start_pos = _start_position(self)
        self.utr3_len = (
            exonic_bases_downstream(self.exons, self.strand, self.stop_pos)
            if self.stop_pos is not None
            else 0
        )
        self.rank = _canonical_rank(self)

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def terminal_exon(self) -> Exon:
        """Last exon in transcript orientation."""
        return self.exons[-1] if self.strand == "+" else self.exons[0]

    @property
    def is_readthrough(self) -> bool:
        return "readthrough_transcript" in self.tags


def _stop_position(t: Transcript) -> Optional[int]:
    """Genomic coordinate immediately 3' of the last coding base."""
    if t.stop:
        return max(e for _, e in t.stop) if t.strand == "+" else min(s for s, _ in t.stop)
    if t.cds:
        return max(e for _, e in t.cds) if t.strand == "+" else min(s for s, _ in t.cds)
    return None


def _start_position(t: Transcript) -> Optional[int]:
    """Genomic coordinate of the first coding base.

    Needed to answer "is this difference actually in the UTR?".  A model whose
    5' end sits downstream of this coordinate does not have a shorter 5'UTR, it
    has a different coding sequence.
    """
    if not t.cds:
        return None
    return min(s for s, _ in t.cds) if t.strand == "+" else max(e for _, e in t.cds)


def _canonical_rank(t: Transcript) -> int:
    tags = t.tags
    if "MANE_Select" in tags:
        return 0
    if "Ensembl_canonical" in tags:
        return 1
    if t.tx_type == "protein_coding" and "basic" in tags:
        return 2
    if "basic" in tags:
        return 3
    if t.tx_type in ("retained_intron", "processed_transcript", "nonsense_mediated_decay"):
        return 8
    return 5


@dataclass
class Gene:
    gene_id: str
    gene_name: str
    contig: str
    strand: str
    tx_ids: List[str] = field(default_factory=list)
    start: int = 0
    end: int = 0
    exon_union: List[Exon] = field(default_factory=list)
    tts_sites: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    tss_sites: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    terminal_exons: List[Exon] = field(default_factory=list)
    utr3_regions: List[Exon] = field(default_factory=list)
    representative: str = ""


def _merge(intervals: Iterable[Exon]) -> List[Exon]:
    out: List[Exon] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


class ReferenceIndex:
    """Everything the collapse needs to know about the annotation."""

    def __init__(self) -> None:
        self.tx: Dict[str, Transcript] = {}
        self.genes: Dict[str, Gene] = {}
        self.chains: Dict[Tuple[str, str], Dict[Chain, List[str]]] = defaultdict(dict)
        self.junctions: Dict[Tuple[str, str], Set[Tuple[int, int]]] = defaultdict(set)
        self.donors: Dict[Tuple[str, str], np.ndarray] = {}
        self.acceptors: Dict[Tuple[str, str], np.ndarray] = {}
        self._donor_sets: Dict[Tuple[str, str], Set[int]] = defaultdict(set)
        self._acceptor_sets: Dict[Tuple[str, str], Set[int]] = defaultdict(set)
        self._gene_bins: Dict[Tuple[str, str], Dict[int, List[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        self._gene_bins_any: Dict[str, Dict[int, List[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        self.name_to_gene: Dict[str, str] = {}
        self._chain_sets: Dict[str, frozenset] = {}

    def chain_set(self, tid: str) -> frozenset:
        """Cached ``set(transcript.chain)`` -- the similarity search is quadratic
        in the transcripts of a gene and rebuilds these constantly."""
        s = self._chain_sets.get(tid)
        if s is None:
            s = frozenset(self.tx[tid].chain)
            self._chain_sets[tid] = s
        return s

    # ------------------------------------------------------------------ #
    @classmethod
    def from_gtf(
        cls,
        path: str,
        contigs: Optional[Sequence[str]] = None,
        verbose: bool = True,
    ) -> "ReferenceIndex":
        self = cls()
        want = set(contigs) if contigs else None
        keep = {"exon", "CDS", "stop_codon"}
        n_lines = 0
        with _open(path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in keep:
                    continue
                if want is not None and f[0] not in want:
                    continue
                a = parse_attrs(f[8])
                tid = a.get("transcript_id")
                if tid is None:
                    continue
                t = self.tx.get(tid)
                if t is None:
                    t = Transcript(
                        tx_id=tid,
                        gene_id=a.get("gene_id", ""),
                        gene_name=a.get("gene_name", a.get("gene_id", "")),
                        contig=f[0],
                        strand=f[6],
                        tx_type=a.get("transcript_type", a.get("transcript_biotype", "")),
                        tags=tuple(v for k, v in _ATTR.findall(f[8]) if k == "tag"),
                    )
                    self.tx[tid] = t
                iv = gtf_to_exon(int(f[3]), int(f[4]))
                if f[2] == "exon":
                    t.exons.append(iv)
                elif f[2] == "CDS":
                    t.cds.append(iv)
                else:
                    t.stop.append(iv)
                n_lines += 1
        self._finalise(verbose=verbose)
        if verbose:
            print(
                f"[ref] {len(self.tx):,} transcripts / {len(self.genes):,} genes / "
                f"{sum(len(v) for v in self.junctions.values()):,} junctions "
                f"from {n_lines:,} GTF records"
            )
        return self

    # ------------------------------------------------------------------ #
    def _finalise(self, verbose: bool = True) -> None:
        drop = [tid for tid, t in self.tx.items() if not t.exons]
        for tid in drop:
            del self.tx[tid]

        by_gene: Dict[str, List[str]] = defaultdict(list)
        for tid, t in self.tx.items():
            t.finalise()
            by_gene[t.gene_id].append(tid)
            key = (t.contig, t.strand)
            if t.chain:
                self.chains[key].setdefault(t.chain, []).append(tid)
                for j in t.chain:
                    self.junctions[key].add(j)
                    self._donor_sets[key].add(j[0])
                    self._acceptor_sets[key].add(j[1])

        for key in self.junctions:
            self.donors[key] = np.array(sorted(self._donor_sets[key]), dtype=np.int64)
            self.acceptors[key] = np.array(sorted(self._acceptor_sets[key]), dtype=np.int64)

        for gid, tids in by_gene.items():
            first = self.tx[tids[0]]
            g = Gene(
                gene_id=gid,
                gene_name=first.gene_name,
                contig=first.contig,
                strand=first.strand,
                tx_ids=sorted(tids, key=lambda x: (self.tx[x].rank, -self.tx[x].tx_len)),
            )
            g.start = min(self.tx[t].start for t in tids)
            g.end = max(self.tx[t].end for t in tids)
            g.exon_union = _merge(e for t in tids for e in self.tx[t].exons)
            g.tts_sites = np.array(sorted({self.tx[t].tts for t in tids}), dtype=np.int64)
            g.tss_sites = np.array(sorted({self.tx[t].tss for t in tids}), dtype=np.int64)
            g.terminal_exons = _merge(self.tx[t].terminal_exon for t in tids)
            g.representative = g.tx_ids[0]
            regions = []
            for t in tids:
                if self.tx[t].stop_pos is None:
                    continue
                sp = self.tx[t].stop_pos
                for s, e in self.tx[t].exons:
                    if first.strand == "+":
                        if e > sp:
                            regions.append((max(s, sp), e))
                    else:
                        if s < sp:
                            regions.append((s, min(e, sp)))
            g.utr3_regions = _merge(r for r in regions if r[1] > r[0])
            self.genes[gid] = g
            self.name_to_gene.setdefault(g.gene_name, gid)
            for b in range(g.start // BIN, g.end // BIN + 1):
                self._gene_bins[(g.contig, g.strand)][b].append(gid)
                self._gene_bins_any[g.contig][b].append(gid)

    # ------------------------------------------------------------------ #
    # queries
    # ------------------------------------------------------------------ #
    def snap_site(self, key: Tuple[str, str], pos: int, tol: int, which: str) -> Optional[int]:
        """Nearest annotated donor/acceptor within ``tol``, else ``None``."""
        arr = (self.donors if which == "donor" else self.acceptors).get(key)
        if arr is None or arr.size == 0:
            return None
        i = int(np.searchsorted(arr, pos))
        best, bd = None, tol + 1
        for j in (i - 1, i):
            if 0 <= j < arr.size:
                d = abs(int(arr[j]) - pos)
                if d < bd:
                    best, bd = int(arr[j]), d
        return best

    def site_distance(self, key: Tuple[str, str], pos: int, which: str) -> int:
        """Distance to the nearest annotated donor/acceptor (unbounded)."""
        arr = (self.donors if which == "donor" else self.acceptors).get(key)
        if arr is None or arr.size == 0:
            return 1 << 20
        i = int(np.searchsorted(arr, pos))
        best = 1 << 20
        for j in (i - 1, i):
            if 0 <= j < arr.size:
                best = min(best, abs(int(arr[j]) - pos))
        return best

    def is_annotated_junction(self, key: Tuple[str, str], j: Tuple[int, int]) -> bool:
        return j in self.junctions.get(key, ())

    def is_annotated_donor(self, key: Tuple[str, str], p: int) -> bool:
        return p in self._donor_sets.get(key, ())

    def is_annotated_acceptor(self, key: Tuple[str, str], p: int) -> bool:
        return p in self._acceptor_sets.get(key, ())

    def transcripts_for_chain(self, key: Tuple[str, str], chain: Chain) -> List[str]:
        return self.chains.get(key, {}).get(chain, [])

    def best_transcript(self, tids: Sequence[str]) -> Optional[str]:
        """Most canonical of a set of transcripts, ignoring the observed ends.

        Kept for callers that have no ends to compare (gene-level defaults).
        Anywhere a model's own ends are known, use
        :meth:`best_transcript_for_ends` instead -- picking MANE there is what
        turned every non-MANE isoform's genuine 3' end into a discovery.
        """
        if not tids:
            return None
        return min(tids, key=lambda t: (self.tx[t].rank, -self.tx[t].tx_len, t))

    def best_transcript_for_ends(
        self,
        tids: Sequence[str],
        five: Optional[int],
        three: Optional[int],
        strand: str,
    ) -> Optional[str]:
        """Of the transcripts in ``tids``, the one whose ends match the model's.

        Several GENCODE transcripts routinely share one intron chain and differ
        only in where they start and stop.  Choosing the MANE one among them and
        then measuring the model's ends against it manufactures "alternative
        3'UTR" calls for every read pile that is in fact an exact match to a
        different annotated isoform.

        The 3' end decides first: with oligo-dT capture it is the anchored end,
        while the 5' end is a percentile over a degraded pile.  The canonical
        rank is consulted only to break a genuine tie, so that two
        indistinguishable transcripts are not named differently in two samples
        of the same cDNA.
        """
        if not tids:
            return None
        if three is None and five is None:
            return self.best_transcript(tids)

        def key(t: str):
            tx = self.tx[t]
            d3 = abs(three - tx.tts) if three is not None else 0
            d5 = abs(five - tx.tss) if five is not None else 0
            return (d3, d5, tx.rank, -tx.tx_len, t)

        return min(tids, key=key)

    def annotated_superchains(self, gene_id: Optional[str], chain: Chain) -> List[str]:
        """Annotated transcripts of the gene that have ``chain`` as a proper suffix.

        A model whose chain is a 5'-truncated form of a longer annotated
        transcript is still called FSM when GENCODE happens to contain the
        short form too (a ``processed_transcript`` fragment, say).  That is a
        deliberate choice, but it moves read mass off the long isoform and onto
        the fragment in proportion to how degraded the library is, so the count
        is reported rather than left implicit.
        """
        if not chain or not gene_id or gene_id not in self.genes:
            return []
        n = len(chain)
        out = []
        for tid in self.genes[gene_id].tx_ids:
            c = self.tx[tid].chain
            if len(c) > n and c[-n:] == tuple(chain):
                out.append(tid)
        return out

    def genes_overlapping(
        self, contig: str, start: int, end: int, strand: Optional[str] = None
    ) -> List[str]:
        bins = self._gene_bins[(contig, strand)] if strand else self._gene_bins_any[contig]
        out: Set[str] = set()
        for b in range(start // BIN, end // BIN + 1):
            for gid in bins.get(b, ()):
                g = self.genes[gid]
                if g.start < end and start < g.end:
                    out.add(gid)
        return sorted(out)

    def assign_gene(
        self, contig: str, exons: Sequence[Exon], strand: Optional[str] = None
    ) -> Optional[str]:
        """Gene with the largest exonic overlap; falls back to span overlap."""
        cands = self.genes_overlapping(contig, exons[0][0], exons[-1][1], strand)
        if not cands:
            return None
        best, best_ov = None, 0
        for gid in cands:
            ov = total_overlap(exons, self.genes[gid].exon_union)
            if ov > best_ov:
                best, best_ov = gid, ov
        if best is not None:
            return best
        # intronic / no exonic overlap: fall back to the enclosing gene span
        return min(cands, key=lambda g: (self.genes[g].end - self.genes[g].start))

    def gene_junctions_within(
        self, gene_id: Optional[str], lo: int, hi: int
    ) -> List[Tuple[int, int]]:
        """Annotated junctions of ``gene_id`` lying entirely inside ``[lo, hi)``.

        This is the discriminator between a 5' extension and a retained first
        intron.  If a model's first exon reaches upstream into a region that
        some isoform splices out, the model did not find a new start site -- it
        read through an intron.
        """
        if not gene_id or gene_id not in self.genes or hi <= lo:
            return []
        out = set()
        for tid in self.genes[gene_id].tx_ids:
            for d, a in self.tx[tid].chain:
                if d >= lo and a <= hi:
                    out.add((int(d), int(a)))
        return sorted(out)

    def rank_genes_by_overlap(
        self, contig: str, exons: Sequence[Exon], strand: Optional[str] = None
    ) -> List[Tuple[str, int]]:
        """``[(gene_id, exonic overlap), ...]`` best first.

        Used to notice when the chain match and the positional assignment
        disagree -- recent paralogues and processed pseudogenes make that
        happen, and the chain match is allowed to win, but silently.
        """
        cands = self.genes_overlapping(contig, exons[0][0], exons[-1][1], strand)
        scored = [(g, total_overlap(exons, self.genes[g].exon_union)) for g in cands]
        scored = [s for s in scored if s[1] > 0]
        scored.sort(key=lambda kv: (-kv[1], kv[0]))
        return scored

    def nearest_tts(self, gene_id: str, pos: int) -> Tuple[int, int]:
        """(signed distance in transcript orientation, nearest annotated 3' end)."""
        g = self.genes[gene_id]
        arr = g.tts_sites
        if arr.size == 0:
            return (1 << 30, -1)
        i = int(np.argmin(np.abs(arr - pos)))
        raw = int(arr[i]) - pos
        return (raw if g.strand == "+" else -raw, int(arr[i]))

    def nearest_tss(self, gene_id: str, pos: int) -> Tuple[int, int]:
        g = self.genes[gene_id]
        arr = g.tss_sites
        if arr.size == 0:
            return (1 << 30, -1)
        i = int(np.argmin(np.abs(arr - pos)))
        raw = int(arr[i]) - pos
        return (raw if g.strand == "+" else -raw, int(arr[i]))

    def in_terminal_exon(self, gene_id: str, exon: Exon) -> float:
        g = self.genes[gene_id]
        if not g.terminal_exons:
            return 0.0
        ov = total_overlap([exon], g.terminal_exons)
        return ov / max(exon[1] - exon[0], 1)

    def in_utr3(self, gene_id: str, exon: Exon) -> float:
        g = self.genes[gene_id]
        if not g.utr3_regions:
            return 0.0
        ov = total_overlap([exon], g.utr3_regions)
        return ov / max(exon[1] - exon[0], 1)

    # ------------------------------------------------------------------ #
    def contigs(self) -> List[str]:
        return sorted({g.contig for g in self.genes.values()})
