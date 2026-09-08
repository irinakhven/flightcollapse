"""Genome-sequence features: splice motifs, RT-switching repeats, polyA evidence.

These are the features that separate a real novel splice junction or a real
alternative polyA site from an alignment artefact, and they are what the
calibrated novelty model in :mod:`flightcollapse.scoring` is built on.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

#: The 12 canonical human polyA signal hexamers, in descending frequency
#: (Beaudoing et al. 2000).  The first two account for ~80% of real sites.
POLYA_MOTIFS: Tuple[str, ...] = (
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
    "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA",
)

#: The subset used when calling a *novel* cleavage site without an atlas.
#: The first two account for ~80% of real sites; the tail members are common
#: enough by chance in an AT-rich 3'UTR to be worthless as evidence.
POLYA_MOTIFS_STRICT: Tuple[str, ...] = (
    "AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "CATAAA", "GATAAA",
)

CANONICAL = "GT-AG"
SEMI_CANONICAL = ("GC-AG", "AT-AC")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


class Genome:
    """Thin, contig-cached wrapper over a pysam FastaFile.

    Holds at most ``cache_contigs`` chromosome sequences in memory (chr1 is
    ~250 MB as a Python str, so the default of 1 is deliberate -- the pipeline
    streams one contig at a time anyway).
    """

    def __init__(self, path: str, cache_contigs: int = 1) -> None:
        import pysam

        self.path = path
        self._fa = pysam.FastaFile(path)
        self._cache: Dict[str, str] = {}
        self._order: List[str] = []
        self._cache_n = max(1, cache_contigs)
        self.references = set(self._fa.references)

    def close(self) -> None:
        self._fa.close()

    def __enter__(self) -> "Genome":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ------------------------------------------------------------------ #
    def resolve(self, contig: str) -> Optional[str]:
        if contig in self.references:
            return contig
        alt = contig[3:] if contig.startswith("chr") else "chr" + contig
        return alt if alt in self.references else None

    def _seq(self, contig: str) -> str:
        c = self.resolve(contig)
        if c is None:
            raise KeyError(f"contig {contig!r} not in {self.path}")
        s = self._cache.get(c)
        if s is None:
            s = self._fa.fetch(c).upper()
            self._cache[c] = s
            self._order.append(c)
            while len(self._order) > self._cache_n:
                self._cache.pop(self._order.pop(0), None)
        return s

    def fetch(self, contig: str, start: int, end: int) -> str:
        """0-based half-open, clipped to the contig, upper case."""
        s = self._seq(contig)
        return s[max(0, start):max(0, min(end, len(s)))]

    def length(self, contig: str) -> int:
        c = self.resolve(contig)
        return len(self._seq(c)) if c else 0

    # ------------------------------------------------------------------ #
    # splice junctions
    # ------------------------------------------------------------------ #
    def splice_motif(self, contig: str, donor: int, acceptor: int, strand: str) -> str:
        """``"GT-AG"``-style dinucleotide pair, read on the transcript strand.

        ``donor``/``acceptor`` are the 0-based half-open intron bounds, i.e.
        the intron is ``[donor, acceptor)`` in genomic coordinates.
        """
        if acceptor - donor < 4:
            return "NN-NN"
        five = self.fetch(contig, donor, donor + 2)
        three = self.fetch(contig, acceptor - 2, acceptor)
        if strand == "-":
            five, three = revcomp(three), revcomp(five)
        if len(five) < 2 or len(three) < 2:
            return "NN-NN"
        return f"{five}-{three}"

    def motif_class(self, motif: str) -> str:
        if motif == CANONICAL:
            return "canonical"
        if motif in SEMI_CANONICAL:
            return "semi_canonical"
        # a junction whose motif is the reverse complement of a canonical one is
        # a strand-assignment error; under a stranded protocol these are almost
        # always artefacts and make an excellent decoy set.
        if motif in ("CT-AC", "CT-GC", "GT-AT"):
            return "antisense"
        return "non_canonical"

    def direct_repeat_len(
        self, contig: str, donor: int, acceptor: int, max_len: int = 20
    ) -> int:
        """Longest direct repeat flanking the intron (RT template-switching signature).

        Reverse transcriptase can jump between two positions sharing a short
        direct repeat, producing a fake intron.  The signature is that the
        sequence ending at the donor site equals the sequence ending at the
        acceptor site (SQANTI's RT-switching test).
        """
        k = min(max_len, donor, acceptor - donor)
        if k <= 0:
            return 0
        up_d = self.fetch(contig, donor - k, donor)
        up_a = self.fetch(contig, acceptor - k, acceptor)
        n = 0
        for i in range(1, min(len(up_d), len(up_a)) + 1):
            if up_d[-i] != up_a[-i]:
                break
            n = i
        return n

    # ------------------------------------------------------------------ #
    # 3' ends
    # ------------------------------------------------------------------ #
    def perc_a_downstream(
        self, contig: str, pos: int, strand: str, window: int = 20
    ) -> float:
        """Percent genomic A in the ``window`` bases downstream of a 3' end.

        Scale is 0-100, matching pigeon's ``perc_A_downstream_TTS``.  The
        SQANTI internal-priming cutoff is 60, not 0.6 -- getting this wrong
        flags essentially every isoform as internally primed.
        """
        if strand == "+":
            s = self.fetch(contig, pos, pos + window)
        else:
            s = revcomp(self.fetch(contig, max(0, pos - window), pos))
        if not s:
            return float("nan")
        return 100.0 * s.count("A") / len(s)

    def polya_motif(
        self,
        contig: str,
        pos: int,
        strand: str,
        window: int = 50,
        min_dist: int = 0,
        max_dist: int = 10 ** 6,
        strict: bool = False,
    ) -> Tuple[bool, str, int]:
        """Search upstream of a 3' end for a polyA signal hexamer.

        Returns ``(found, motif, distance)`` where ``distance`` is the number of
        bases from the hexamer start to the cleavage site (typically 10-30).
        The highest-priority (most frequent) motif found wins.
        """
        if strand == "+":
            s = self.fetch(contig, max(0, pos - window), pos)
            offset_from_end = lambda i: len(s) - i  # noqa: E731
        else:
            s = revcomp(self.fetch(contig, pos, pos + window))
            offset_from_end = lambda i: len(s) - i  # noqa: E731
        if not s:
            return (False, "", -1)
        motifs = POLYA_MOTIFS_STRICT if strict else POLYA_MOTIFS
        best: Tuple[bool, str, int] = (False, "", -1)
        for m in motifs:
            start = 0
            while True:
                i = s.find(m, start)
                if i < 0:
                    break
                d = offset_from_end(i)
                if min_dist <= d <= max_dist:
                    return (True, m, d)
                start = i + 1
        return best


class PolyASiteAtlas:
    """Optional external polyA-site catalogue (PolyASite 2.0 / PolyA_DB v4 BED).

    Turns "is this 3' end a real cleavage site?" from a motif guess into a
    lookup.  Strongly recommended for the mono-exonic track, where the whole
    question is whether an internal 3' end is a genuine alternative polyA site.
    """

    def __init__(self) -> None:
        self._sites: Dict[Tuple[str, str], np.ndarray] = {}

    @classmethod
    def from_bed(cls, path: str, verbose: bool = True) -> "PolyASiteAtlas":
        import gzip as _gz
        from collections import defaultdict

        self = cls()
        acc: Dict[Tuple[str, str], List[int]] = defaultdict(list)
        opener = _gz.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line[0] in "#t":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                contig = f[0] if f[0].startswith("chr") else "chr" + f[0]
                strand = f[5] if len(f) > 5 and f[5] in "+-" else "+"
                start, end = int(f[1]), int(f[2])
                pos = end if strand == "+" else start
                acc[(contig, strand)].append(pos)
        for k, v in acc.items():
            self._sites[k] = np.array(sorted(set(v)), dtype=np.int64)
        if verbose:
            print(f"[polyA] {sum(a.size for a in self._sites.values()):,} sites from {path}")
        return self

    def distance(self, contig: str, pos: int, strand: str) -> int:
        arr = self._sites.get((contig, strand))
        if arr is None:
            arr = self._sites.get(
                (contig[3:] if contig.startswith("chr") else "chr" + contig, strand)
            )
        if arr is None or arr.size == 0:
            return 1 << 30
        i = int(np.searchsorted(arr, pos))
        best = 1 << 30
        for j in (i - 1, i):
            if 0 <= j < arr.size:
                best = min(best, abs(int(arr[j]) - pos))
        return best
