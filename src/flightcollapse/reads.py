"""BAM scanning.

One contig at a time, into compact parallel numpy arrays rather than Python
objects -- whole-genome collapse OOM-killed the VM on the merged-organoid
reference, and a naive ``list[ReadRecord]`` for chr1 costs ~5 GB where the
array form costs ~1 GB.

Read *names* are deliberately never held in memory.  They are needed only when
writing ``read_stat.txt`` / ``group.txt``, and the alignment iteration order is
deterministic, so the writer simply re-walks the BAM and pairs the n-th
admitted read with row ``n`` of the arrays.  :func:`iter_alignments` is the
single definition of "admitted read" that both passes share; if you change a
gate, both passes change together.
"""

from __future__ import annotations

from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from .config import AlignmentGates
from .intervals import Chain, Exon

# CIGAR opcodes
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP = 0, 1, 2, 3
BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF = 4, 5, 6, 7, 8


def resolve_contig(bam, contig: str) -> Optional[str]:
    refs = set(bam.references)
    for c in (contig, contig[3:] if contig.startswith("chr") else "chr" + contig):
        if c in refs:
            return c
    return None


def alignment_stats(read) -> Tuple[float, float, int, int]:
    """``(coverage, identity, clip5, clip3)`` for a long-read alignment.

    * coverage = aligned query bases / full query length (hard clips included)
    * identity = 1 - NM / (aligned columns, introns excluded)
    * ``clip5`` / ``clip3`` are the terminal soft+hard clips **in transcript
      orientation**, so ``clip3`` is where an untrimmed polyA tail lives.

    Reporting the clips separately matters: on an untrimmed FLNC BAM the polyA
    tail is soft-clipped, which drags coverage down for reasons that have
    nothing to do with alignment quality.  A ``--min-aln-coverage 0.99`` gate
    applied to such a BAM keeps ~6% of reads and *enriches* for short,
    mono-exonic ones.
    """
    q_total = read.infer_read_length() or read.query_length or 0
    aligned = 0
    columns = 0
    cig = read.cigartuples or []
    for op, ln in cig:
        if op in (BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF, BAM_CINS):
            aligned += ln
            if op != BAM_CINS:
                columns += ln
            else:
                columns += ln
        elif op == BAM_CDEL:
            columns += ln
    cov = aligned / q_total if q_total else 0.0
    try:
        nm = read.get_tag("NM")
    except KeyError:
        nm = 0
    ident = 1.0 - (nm / columns) if columns else 0.0

    def _clip(ops):
        n = 0
        for op, ln in ops:
            if op in (BAM_CSOFT_CLIP, BAM_CHARD_CLIP):
                n += ln
            else:
                break
        return n

    left, right = _clip(cig), _clip(reversed(cig))
    clip5, clip3 = (right, left) if read.is_reverse else (left, right)
    return cov, ident, clip5, clip3


def polya_tail_stats(read, scan: int = 60) -> Tuple[int, int, float]:
    """``(soft_clip3, terminal A-run, A-fraction)`` of the 3' clip.

    New in 0.1.17.  ``alignment_stats`` already measures how *long* the 3' clip
    is; this reads what is actually in it.  Under oligo-dT priming the tail is
    the physical thing that was captured, so a run of A's immediately outside
    the aligned block -- in transcript orientation -- is the one piece of direct
    per-molecule evidence that this coordinate is a cleavage site rather than an
    internal-priming or degradation artefact.

    Only SOFT clips are counted, because hard-clipped bases are not in
    ``query_sequence``; a hard-clipped BAM therefore reports 0 rather than
    silently indexing into the wrong bases.  On a minus-strand alignment the
    stored read is in reference orientation, so the tail sits at the START of
    the sequence and reads as poly-T; the slice is reversed so that index 0 is
    always the base adjacent to the cleavage site.
    """
    cig = read.cigartuples or []
    if not cig:
        return 0, 0, 0.0

    def _soft(ops):
        n = 0
        for op, ln in ops:
            if op == BAM_CSOFT_CLIP:
                n += ln
            elif op == BAM_CHARD_CLIP:
                continue
            else:
                break
        return n

    clip3 = _soft(cig) if read.is_reverse else _soft(reversed(cig))
    if clip3 <= 0:
        return 0, 0, 0.0
    seq = read.query_sequence
    if not seq:
        return clip3, 0, 0.0
    clip3 = min(clip3, len(seq))
    if read.is_reverse:
        tail = seq[:clip3][::-1].upper()
        base = "T"
    else:
        tail = seq[len(seq) - clip3:].upper()
        base = "A"
    run = 0
    for ch in tail:
        if ch != base:
            break
        run += 1
    window = tail[:scan]
    frac = window.count(base) / len(window) if window else 0.0
    return clip3, run, frac


def introns_of(read, min_len: int, max_len: int) -> List[Tuple[int, int]]:
    """``N`` operations as 0-based half-open ``(donor, acceptor)``, genomic order."""
    out: List[Tuple[int, int]] = []
    ref = read.reference_start
    for op, ln in read.cigartuples or ():
        if op in (BAM_CMATCH, BAM_CDEL, BAM_CEQUAL, BAM_CDIFF):
            ref += ln
        elif op == BAM_CREF_SKIP:
            if min_len <= ln <= max_len:
                out.append((ref, ref + ln))
            ref += ln
    return out


def internal_gaps(read, min_len: int) -> List[Tuple[int, int]]:
    """Long ``D`` operations as 0-based half-open reference intervals.

    An exact intron-chain match guarantees that every *internal* exon boundary
    is identical to the reference, because those boundaries are defined by the
    flanking junctions.  It does not guarantee the exon *contents* are: a
    deletion is a ``D``, not an ``N``, so two structures with the same chain can
    still differ inside an exon and the comparison never notices.

    That bucket holds two very different things.  A ~20-60 bp gap with GT..AG at
    its boundaries, recurring across molecules, is a short intron the aligner
    declined to call as ``N``; leaving it as a name suffix would bury a real
    splicing event inside an "FSM".  Everything else is alignment noise.  Both
    are collected here; :mod:`flightcollapse.junctions` decides which is which.
    """
    out: List[Tuple[int, int]] = []
    ref = read.reference_start
    for op, ln in read.cigartuples or ():
        if op in (BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF):
            ref += ln
        elif op == BAM_CDEL:
            if ln >= min_len:
                out.append((ref, ref + ln))
            ref += ln
        elif op == BAM_CREF_SKIP:
            ref += ln
    return out


def exon_blocks(read) -> List[Exon]:
    blocks: List[Exon] = []
    ref = read.reference_start
    s = e = ref
    for op, ln in read.cigartuples or ():
        if op in (BAM_CMATCH, BAM_CDEL, BAM_CEQUAL, BAM_CDIFF):
            e = ref + ln
            ref += ln
        elif op == BAM_CREF_SKIP:
            if e > s:
                blocks.append((s, e))
            ref += ln
            s = e = ref
    if e > s:
        blocks.append((s, e))
    return blocks


def iter_alignments(
    bam, contig: str, gates: AlignmentGates, stats: Optional[Dict[str, int]] = None
) -> Iterator:
    """The one definition of an admitted read.  Deterministic order.

    Both the scanning pass and the output pass go through here, so a gate
    change can never desynchronise them.  ``stats`` accumulates a rejection
    breakdown, which is reported per contig -- a gate that quietly discards
    most of the data is the kind of thing that should be loud.
    """
    c = resolve_contig(bam, contig)
    if c is None:
        return

    def _bump(k: str) -> None:
        if stats is not None:
            stats[k] = stats.get(k, 0) + 1

    for r in bam.fetch(c):
        if r.is_unmapped:
            continue
        _bump("scanned")
        if gates.drop_secondary and r.is_secondary:
            _bump("secondary")
            continue
        if gates.drop_supplementary and r.is_supplementary:
            _bump("supplementary")
            continue
        if r.mapping_quality < gates.min_mapq:
            _bump("low_mapq")
            continue
        cov, ident, clip5, clip3 = alignment_stats(r)
        if cov < gates.min_aln_coverage:
            _bump("low_coverage")
            continue
        if ident < gates.min_aln_identity:
            _bump("low_identity")
            continue
        if gates.max_5p_softclip is not None and clip5 > gates.max_5p_softclip:
            _bump("long_5p_clip")
            continue
        if gates.max_3p_softclip is not None and clip3 > gates.max_3p_softclip:
            _bump("long_3p_clip")
            continue
        _bump("admitted")
        yield r


def gate_diagnostics(bam, contig: str, gates: AlignmentGates, max_reads: int = 200_000):
    """Coverage / identity / clip quantiles over the raw alignments.

    Used to explain a gate rejection rather than just reporting one.
    """
    c = resolve_contig(bam, contig)
    if c is None:
        return {}
    cov: List[float] = []
    ident: List[float] = []
    c5: List[int] = []
    c3: List[int] = []
    for r in bam.fetch(c):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        a, b, x, y = alignment_stats(r)
        cov.append(a)
        ident.append(b)
        c5.append(x)
        c3.append(y)
        if len(cov) >= max_reads:
            break
    if not cov:
        return {}
    q = [1, 5, 10, 25, 50]
    return {
        "n": len(cov),
        "coverage_percentiles": {
            f"p{p}": round(float(np.percentile(cov, p)), 4) for p in q
        },
        "identity_percentiles": {
            f"p{p}": round(float(np.percentile(ident, p)), 4) for p in q
        },
        "median_5p_clip": int(np.median(c5)),
        "median_3p_clip": int(np.median(c3)),
        "p90_5p_clip": int(np.percentile(c5, 90)),
        "p90_3p_clip": int(np.percentile(c3, 90)),
    }


class ContigReads:
    """Compact per-contig read table.

    Attributes are parallel arrays of length ``n``:

    ``strand``   int8, +1 / -1
    ``tss``      int64, genomic coordinate of the read 5' end
    ``tts``      int64, genomic coordinate of the read 3' end
    ``start``    int64, leftmost aligned base
    ``end``      int64, rightmost aligned base (exclusive)
    ``joff``     int64, offset into ``donor``/``acceptor`` (length n+1)
    ``cell``     uint32, sentinel 0xFFFFFFFF when unknown
    ``umi``      uint64, sentinel 0xFFFFFFFFFFFFFFFF when unknown
    ``clip3``    uint16, 3' SOFT clip length, transcript orientation (0.1.17)
    ``tail_len`` uint16, terminal A-run inside that clip (0.1.17)
    ``tail_frac`` float32, A-fraction over the clip's first 60 bases (0.1.17)

    The three tail arrays are zero-filled when ``AlignmentGates.measure_polya_tail``
    is off, so downstream code may read them unconditionally; it must not read a
    zero as "no tail" without checking that flag.
    """

    __slots__ = (
        "contig", "n", "strand", "tss", "tts", "start", "end",
        "joff", "donor", "acceptor", "cell", "umi", "n_scanned", "n_admitted",
        "goff", "gap_start", "gap_end", "clip3", "tail_len", "tail_frac",
    )

    def __init__(self, contig: str) -> None:
        self.contig = contig
        self.n = 0
        self.n_scanned = 0
        self.n_admitted = 0
        self.strand = np.zeros(0, np.int8)
        self.tss = np.zeros(0, np.int64)
        self.tts = np.zeros(0, np.int64)
        self.start = np.zeros(0, np.int64)
        self.end = np.zeros(0, np.int64)
        self.joff = np.zeros(1, np.int64)
        self.donor = np.zeros(0, np.int64)
        self.acceptor = np.zeros(0, np.int64)
        self.cell = np.zeros(0, np.uint32)
        self.umi = np.zeros(0, np.uint64)
        #: long exon-internal deletions, same offset-array layout as junctions
        self.goff = np.zeros(1, np.int64)
        self.gap_start = np.zeros(0, np.int64)
        self.gap_end = np.zeros(0, np.int64)
        self.clip3 = np.zeros(0, np.uint16)
        self.tail_len = np.zeros(0, np.uint16)
        self.tail_frac = np.zeros(0, np.float32)

    # ------------------------------------------------------------------ #
    @classmethod
    def scan(
        cls,
        bam,
        contig: str,
        gates: AlignmentGates,
        molecule_index=None,
        name_batch: int = 200_000,
        gate_stats: Optional[Dict[str, int]] = None,
    ) -> "ContigReads":
        self = cls(contig)
        strand: List[int] = []
        tss: List[int] = []
        tts: List[int] = []
        start: List[int] = []
        end: List[int] = []
        joff: List[int] = [0]
        don: List[int] = []
        acc: List[int] = []
        goff: List[int] = [0]
        gs: List[int] = []
        ge: List[int] = []
        c3: List[int] = []
        tl: List[int] = []
        tf: List[float] = []
        names: List[str] = []
        cells: List[np.ndarray] = []
        umis: List[np.ndarray] = []

        def _flush() -> None:
            if not names:
                return
            if molecule_index is None:
                from .molecules import NO_CELL, NO_UMI

                cells.append(np.full(len(names), NO_CELL, np.uint32))
                umis.append(np.full(len(names), NO_UMI, np.uint64))
            else:
                c, u = molecule_index.lookup(names)
                cells.append(c)
                umis.append(u)
            names.clear()

        n_admitted = 0
        for r in iter_alignments(bam, contig, gates, gate_stats):
            n_admitted += 1
            iv = introns_of(r, gates.min_intron_len, gates.max_intron_len)
            s = -1 if r.is_reverse else 1
            strand.append(s)
            start.append(r.reference_start)
            end.append(r.reference_end)
            if s == 1:
                tss.append(r.reference_start)
                tts.append(r.reference_end)
            else:
                tss.append(r.reference_end)
                tts.append(r.reference_start)
            for d, a in iv:
                don.append(d)
                acc.append(a)
            joff.append(len(don))
            for a, b in internal_gaps(r, gates.min_intron_len):
                gs.append(a)
                ge.append(b)
            goff.append(len(gs))
            if gates.measure_polya_tail:
                a, b, c = polya_tail_stats(r)
                c3.append(min(a, 65535))
                tl.append(min(b, 65535))
                tf.append(c)
            names.append(r.query_name)
            if len(names) >= name_batch:
                _flush()
        _flush()

        self.n = len(strand)
        self.n_admitted = n_admitted
        self.strand = np.array(strand, np.int8)
        self.tss = np.array(tss, np.int64)
        self.tts = np.array(tts, np.int64)
        self.start = np.array(start, np.int64)
        self.end = np.array(end, np.int64)
        self.joff = np.array(joff, np.int64)
        self.donor = np.array(don, np.int64)
        self.acceptor = np.array(acc, np.int64)
        self.goff = np.array(goff, np.int64)
        self.gap_start = np.array(gs, np.int64)
        self.gap_end = np.array(ge, np.int64)
        if gates.measure_polya_tail and len(c3) == self.n:
            self.clip3 = np.array(c3, np.uint16)
            self.tail_len = np.array(tl, np.uint16)
            self.tail_frac = np.array(tf, np.float32)
        else:
            self.clip3 = np.zeros(self.n, np.uint16)
            self.tail_len = np.zeros(self.n, np.uint16)
            self.tail_frac = np.zeros(self.n, np.float32)
        from .molecules import NO_CELL, NO_UMI

        self.cell = (
            np.concatenate(cells) if cells else np.zeros(self.n, np.uint32) + NO_CELL
        )
        self.umi = np.concatenate(umis) if umis else np.zeros(self.n, np.uint64) + NO_UMI
        return self

    # ------------------------------------------------------------------ #
    def n_introns(self, i: int) -> int:
        return int(self.joff[i + 1] - self.joff[i])

    def is_monoexonic(self, i: int) -> bool:
        return self.joff[i + 1] == self.joff[i]

    def chain(self, i: int) -> Chain:
        """Intron chain in transcript orientation."""
        a, b = int(self.joff[i]), int(self.joff[i + 1])
        if a == b:
            return ()
        ch = tuple(zip(self.donor[a:b].tolist(), self.acceptor[a:b].tolist()))
        return ch[::-1] if self.strand[i] < 0 else ch

    def strand_char(self, i: int) -> str:
        return "+" if self.strand[i] > 0 else "-"

    def n_gaps(self, i: int) -> int:
        return int(self.goff[i + 1] - self.goff[i])

    def has_internal_gap(self, i: int) -> bool:
        return self.goff[i + 1] > self.goff[i]

    def gap_census(self) -> Dict[Tuple[str, int, int], int]:
        """``(strand, start, end) -> reads carrying this exon-internal deletion``."""
        out: Dict[Tuple[str, int, int], int] = {}
        if self.gap_start.size == 0:
            return out
        counts = np.diff(self.goff)
        gstrand = np.repeat(self.strand, counts)
        for s, a, b in zip(gstrand, self.gap_start, self.gap_end):
            k = ("+" if s > 0 else "-", int(a), int(b))
            out[k] = out.get(k, 0) + 1
        return out

    def promote_gaps(self, promote: "set") -> int:
        """Move selected exon-internal deletions into the intron chain.

        ``promote`` holds ``(strand_char, start, end)`` keys.  Reads keep their
        deletion records either way -- a non-promoted gap still flags the read
        as carrying an unexplained internal difference -- but a promoted gap
        becomes an ordinary junction and goes through snapping, feature
        extraction and curation like any other, which is the point: a real short
        intron should be judged on the same evidence as every other junction,
        not accepted because of how the aligner happened to spell it.

        Returns the number of junction records added.
        """
        if not promote or self.gap_start.size == 0:
            return 0
        counts_g = np.diff(self.goff)
        don: List[int] = []
        acc: List[int] = []
        joff: List[int] = [0]
        added = 0
        for i in range(self.n):
            a, b = int(self.joff[i]), int(self.joff[i + 1])
            merged = list(zip(self.donor[a:b].tolist(), self.acceptor[a:b].tolist()))
            if counts_g[i]:
                sc = "+" if self.strand[i] > 0 else "-"
                ga, gb = int(self.goff[i]), int(self.goff[i + 1])
                extra = [
                    (int(x), int(y))
                    for x, y in zip(self.gap_start[ga:gb], self.gap_end[ga:gb])
                    if (sc, int(x), int(y)) in promote
                ]
                if extra:
                    merged = sorted(merged + extra)
                    added += len(extra)
            for d, ac in merged:
                don.append(d)
                acc.append(ac)
            joff.append(len(don))
        if added:
            self.donor = np.array(don, np.int64)
            self.acceptor = np.array(acc, np.int64)
            self.joff = np.array(joff, np.int64)
        return added

    def pair(self, i: int) -> Tuple[int, int]:
        return (int(self.cell[i]), int(self.umi[i]))

    def mono_exon(self, i: int) -> Exon:
        return (int(self.start[i]), int(self.end[i]))

    # ------------------------------------------------------------------ #
    def apply_snap(
        self,
        donor_map: Dict[Tuple[str, int], int],
        acceptor_map: Dict[Tuple[str, int], int],
    ) -> None:
        """Re-key every junction onto the curated splice-site catalogue.

        Maps are keyed by ``(strand_char, position)`` so that a coordinate can
        snap differently on the two strands (they have different annotated
        catalogues).
        """
        if self.n == 0 or self.donor.size == 0:
            return
        counts = np.diff(self.joff)
        jstrand = np.repeat(self.strand, counts)  # int8, +1 / -1
        for arr, mapping in ((self.donor, donor_map), (self.acceptor, acceptor_map)):
            for schar, sint in (("+", 1), ("-", -1)):
                items = [(k[1], v) for k, v in mapping.items() if k[0] == schar]
                if not items:
                    continue
                items.sort()
                keys = np.fromiter((k for k, _ in items), np.int64, len(items))
                vals = np.fromiter((v for _, v in items), np.int64, len(items))
                sel = np.flatnonzero(jstrand == sint)
                if sel.size == 0:
                    continue
                sub = arr[sel]
                idx = np.searchsorted(keys, sub)
                idx_c = np.clip(idx, 0, keys.size - 1)
                hit = keys[idx_c] == sub
                if hit.any():
                    arr[sel[hit]] = vals[idx_c[hit]]

    # ------------------------------------------------------------------ #
    def summary(self) -> Dict[str, int]:
        mono = int(np.sum(np.diff(self.joff) == 0))
        return {
            "reads": self.n,
            "mono_exonic_reads": mono,
            "spliced_reads": self.n - mono,
            "junction_records": int(self.donor.size),
        }
