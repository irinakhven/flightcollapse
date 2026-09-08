"""Small interval / coordinate helpers.

Coordinate convention used everywhere in this package
-----------------------------------------------------
* Exons are 0-based, half-open ``[start, end)``  (BED / pysam style).
* An intron derived from a CIGAR ``N`` operation is ``(donor, acceptor)`` where
  ``donor`` is the first intronic base and ``acceptor`` is the first exonic base
  after the intron -- i.e. also 0-based half-open.  This is exactly what
  ``exon_end, next_exon_start`` gives you, on either strand, so read-derived and
  GTF-derived junctions are directly comparable without any +/-1 bookkeeping.
* GTF is 1-based inclusive; :func:`gtf_to_exon` converts.
* "Transcript orientation" means the chain tuple is reversed on the minus strand,
  so element ``[-1]`` is always the 3'-most intron and a 5'-truncated read is
  always a *suffix* of the chain it came from.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

Exon = Tuple[int, int]
Junction = Tuple[int, int]
Chain = Tuple[Junction, ...]


def gtf_to_exon(start: int, end: int) -> Exon:
    """GTF 1-based inclusive -> 0-based half-open."""
    return (start - 1, end)


def exons_to_chain(exons: Sequence[Exon], strand: str) -> Chain:
    """Ordered exons (genomic, ascending) -> intron chain in transcript orientation."""
    if len(exons) < 2:
        return ()
    ch = tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))
    return ch[::-1] if strand == "-" else ch


def chain_to_genomic(chain: Chain, strand: str) -> Chain:
    """Transcript-orientation chain -> ascending genomic order."""
    return chain[::-1] if strand == "-" else chain


def exonic_length(exons: Iterable[Exon]) -> int:
    return sum(e - s for s, e in exons)


def span(exons: Sequence[Exon]) -> int:
    return exons[-1][1] - exons[0][0]


def overlap(a: Exon, b: Exon) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def total_overlap(a: Sequence[Exon], b: Sequence[Exon]) -> int:
    """Exonic overlap between two ascending, non-overlapping exon lists."""
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        tot += overlap(a[i], b[j])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def tss_of(exons: Sequence[Exon], strand: str) -> int:
    """Genomic coordinate of the transcript 5' end."""
    return exons[0][0] if strand == "+" else exons[-1][1]


def tts_of(exons: Sequence[Exon], strand: str) -> int:
    """Genomic coordinate of the transcript 3' end."""
    return exons[-1][1] if strand == "+" else exons[0][0]


def is_suffix(short: Chain, long: Chain) -> bool:
    """True if ``short`` is a *proper* suffix of ``long``.

    Both must be in transcript orientation.  The empty chain is deliberately
    NOT a suffix of anything: the vacuous-suffix case is precisely the bug that
    let mono-exonic reads absorb full-length spliced read mass in
    ``isoseq collapse``.
    """
    n = len(short)
    if n == 0 or n >= len(long):
        return False
    return long[-n:] == short


def build_exons(chain: Chain, strand: str, tss: int, tts: int) -> List[Exon]:
    """Rebuild an exon structure from an intron chain plus the two ends.

    ``tss``/``tts`` are genomic coordinates of the 5'/3' ends.  The result is in
    ascending genomic order, 0-based half-open.
    """
    g = chain_to_genomic(chain, strand)
    left, right = (tss, tts) if strand == "+" else (tts, tss)
    if not g:
        return [(left, right)]
    exons = [(left, g[0][0])]
    for i in range(len(g) - 1):
        exons.append((g[i][1], g[i + 1][0]))
    exons.append((g[-1][1], right))
    # guard against ends that fall inside the first/last intron after clustering
    exons[0] = (min(exons[0][0], exons[0][1] - 1), exons[0][1])
    exons[-1] = (exons[-1][0], max(exons[-1][1], exons[-1][0] + 1))
    return exons


def exonic_bases_downstream(exons: Sequence[Exon], strand: str, pos: int) -> int:
    """Exonic bases of ``exons`` strictly 3' of genomic coordinate ``pos``."""
    tot = 0
    if strand == "+":
        for s, e in exons:
            if e <= pos:
                continue
            tot += e - max(s, pos)
    else:
        for s, e in exons:
            if s >= pos:
                continue
            tot += min(e, pos) - s
    return max(tot, 0)


def signed_5p_offset(a: int, b: int, strand: str) -> int:
    """``a - b`` in transcript orientation; negative = ``a`` is upstream of ``b``."""
    return (a - b) if strand == "+" else (b - a)
