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


# ---------------------------------------------------------------------- #
# CAGE / TSS peak atlas
# ---------------------------------------------------------------------- #
class TSSAtlas:
    """Optional external CAGE peak catalogue (FANTOM5 phase1&2, refTSS, ...).

    The 5' counterpart of :class:`PolyASiteAtlas`, and deliberately *not* the
    same shape.  A polyA site is a point, so the polyA atlas reduces every
    record to one coordinate.  A CAGE peak is an interval 20-100 bp wide that
    additionally carries a representative TSS in the BED9 thick fields, and the
    two facts answer different questions:

    * ``in_peak``   -- the 5' end lies inside a region where transcription is
      known to initiate.  Weaker, because a peak is wide.
    * ``dist_to_reptss`` -- the 5' end is close to the *modal* start of that
      peak.  Stronger, and directly comparable to ``EndParams.max_5p_diff``.

    :mod:`flightcollapse.pipeline` tiers them in that order.

    The atlas can only ever *support* a start site.  Absence of a peak is not
    evidence of absence -- FANTOM5 has no retinal organoids and no brain nuclei
    -- so nothing in this package may reject a model for missing one.  This is
    the same asymmetry :mod:`flightcollapse.shortreads` states for junctions.

    Input formats
    -------------
    * BED9 (native FANTOM5): columns 7/8 are the representative TSS.
    * BED6: no representative TSS.  Either supply one via ``reptss_path``
      (a BED6 joined on the peak name in column 4 -- which is what you have
      after lifting a mouse file to GRCm39, where liftOver drops the thick
      fields) or fall back to the peak's 5'-most base by strand.
    """

    __slots__ = ("_peaks", "n_peaks", "n_reptss_joined", "max_width",
                 "sources", "min_score")

    def __init__(self) -> None:
        #: (contig, strand) -> dict of parallel arrays sorted by start
        self._peaks: Dict[Tuple[str, str], Dict[str, object]] = {}
        self.n_peaks = 0
        self.n_reptss_joined = 0
        self.max_width = 0
        self.sources: List[str] = []
        self.min_score = 0.0

    # ------------------------------------------------------------------ #
    @staticmethod
    def _open(path: str):
        import gzip as _gz

        return _gz.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

    @staticmethod
    def _norm(contig: str) -> str:
        return contig if contig.startswith("chr") else "chr" + contig

    @classmethod
    def _read_bed(cls, path: str, min_score: float):
        """Yield ``(contig, strand, start, end, name, reptss_or_None)``."""
        with cls._open(path) as fh:
            for line in fh:
                if not line or line[0] == "#" or line.startswith(("track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                strand = f[5]
                if strand not in ("+", "-"):
                    continue
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                if min_score:
                    try:
                        if float(f[4]) < min_score:
                            continue
                    except ValueError:
                        pass
                rep = None
                if len(f) >= 8:
                    try:
                        r = int(f[6])
                        # a BED9 thickStart of 0 on a non-zero peak is "unset",
                        # not a coordinate on chr1
                        if start <= r < end:
                            rep = r
                    except ValueError:
                        rep = None
                yield cls._norm(f[0]), strand, start, end, f[3], rep

    @classmethod
    def from_bed(
        cls,
        path: str,
        reptss_path: Optional[str] = None,
        min_score: float = 0.0,
        verbose: bool = True,
    ) -> "TSSAtlas":
        from collections import defaultdict as _dd

        self = cls()
        self.min_score = float(min_score)
        self.sources = [path] + ([reptss_path] if reptss_path else [])

        # optional representative-TSS side table, keyed on peak name.  liftOver
        # can split one peak into two intervals, so a name is NOT unique and the
        # join has to pick the position that falls inside the interval.
        side: Dict[str, List[Tuple[str, str, int]]] = _dd(list)
        if reptss_path:
            for c, s, st, en, name, _r in cls._read_bed(reptss_path, 0.0):
                side[name].append((c, s, st))

        acc: Dict[Tuple[str, str], List[Tuple[int, int, int, str]]] = _dd(list)
        for contig, strand, start, end, name, rep in cls._read_bed(path, self.min_score):
            if rep is None and name in side:
                for c, s, pos in side[name]:
                    if c == contig and s == strand and start <= pos < end:
                        rep = pos
                        self.n_reptss_joined += 1
                        break
            if rep is None:
                rep = start if strand == "+" else end - 1
            acc[(contig, strand)].append((start, end, rep, name))

        for key, rows in acc.items():
            rows.sort()
            starts = np.array([r[0] for r in rows], np.int64)
            ends = np.array([r[1] for r in rows], np.int64)
            reps = np.array([r[2] for r in rows], np.int64)
            self._peaks[key] = {
                "start": starts,
                "end": ends,
                "reptss": reps,
                "name": [r[3] for r in rows],
                # nearest-representative lookups need their own sorted view
                "reptss_sorted": np.sort(reps),
            }
            self.n_peaks += len(rows)
            w = int((ends - starts).max()) if len(rows) else 0
            self.max_width = max(self.max_width, w)

        if verbose:
            joined = (f", {self.n_reptss_joined:,} representative TSS joined"
                      if reptss_path else "")
            print(f"[cage] {self.n_peaks:,} peaks from {path}"
                  f" (max width {self.max_width:,} bp{joined})")
        return self

    # ------------------------------------------------------------------ #
    def _get(self, contig: str, strand: str):
        d = self._peaks.get((contig, strand))
        if d is None:
            alt = contig[3:] if contig.startswith("chr") else "chr" + contig
            d = self._peaks.get((alt, strand))
        return d

    def __len__(self) -> int:
        return self.n_peaks

    # ------------------------------------------------------------------ #
    def evidence(
        self,
        contig: str,
        pos: int,
        strand: str,
        peak_slack: int = 0,
        reptss_window: int = 100,
    ) -> Dict[str, object]:
        """Everything the tiering in the pipeline needs, in one lookup."""
        out: Dict[str, object] = {
            "in_cage_peak": False,
            "dist_to_cage_peak": None,
            "dist_to_cage_reptss": None,
            "cage_peak_id": "",
        }
        d = self._get(contig, strand)
        if d is None:
            return out
        starts = d["start"]
        ends = d["end"]
        if starts.size == 0:
            return out

        # Peaks may overlap after a liftover, so containment cannot assume a
        # disjoint set.  Every peak that could contain `pos` starts no earlier
        # than pos - max_width and no later than pos, which bounds the scan to
        # a handful of rows.
        lo = int(np.searchsorted(starts, pos - self.max_width - peak_slack, "left"))
        hi = int(np.searchsorted(starts, pos + peak_slack, "right"))
        best_d = 1 << 30
        best_i = -1
        for i in range(lo, hi):
            s, e = int(starts[i]), int(ends[i])
            if s - peak_slack <= pos < e + peak_slack:
                out["in_cage_peak"] = True
                out["dist_to_cage_peak"] = 0
                out["cage_peak_id"] = d["name"][i]
                best_i = i
                best_d = 0
                break
            gap = s - pos if pos < s else pos - e + 1
            if gap < best_d:
                best_d, best_i = gap, i
        if best_i < 0:
            # nothing to the left within max_width; look right for the nearest
            j = int(np.searchsorted(starts, pos, "left"))
            for k in (j - 1, j):
                if 0 <= k < starts.size:
                    s, e = int(starts[k]), int(ends[k])
                    gap = s - pos if pos < s else max(pos - e + 1, 0)
                    if gap < best_d:
                        best_d, best_i = gap, k
        if best_i >= 0 and not out["in_cage_peak"]:
            out["dist_to_cage_peak"] = int(best_d)
            out["cage_peak_id"] = d["name"][best_i]

        rs = d["reptss_sorted"]
        i = int(np.searchsorted(rs, pos))
        best = 1 << 30
        for j in (i - 1, i):
            if 0 <= j < rs.size:
                best = min(best, abs(int(rs[j]) - pos))
        out["dist_to_cage_reptss"] = int(best)
        return out

    def in_peak(self, contig: str, pos: int, strand: str, slack: int = 0) -> bool:
        return bool(self.evidence(contig, pos, strand, slack)["in_cage_peak"])

    def dist_to_reptss(self, contig: str, pos: int, strand: str) -> int:
        v = self.evidence(contig, pos, strand)["dist_to_cage_reptss"]
        return 1 << 30 if v is None else int(v)

    # ------------------------------------------------------------------ #
    def flip_strands(self) -> "TSSAtlas":
        """Return a copy with every peak's strand inverted.

        The negative control for the rescue rule: any 5' end still "supported"
        after this is supported by position alone, so the number of models
        rescued under the flip is a direct false-positive estimate.
        """
        other = TSSAtlas()
        other.min_score = self.min_score
        other.max_width = self.max_width
        other.n_peaks = self.n_peaks
        other.n_reptss_joined = self.n_reptss_joined
        other.sources = list(self.sources) + ["<strand-flipped>"]
        for (contig, strand), d in self._peaks.items():
            other._peaks[(contig, "-" if strand == "+" else "+")] = d
        return other

    # ------------------------------------------------------------------ #
    def coverage_of_annotated_tss(self, ref, max_genes: int = 0) -> Dict[str, object]:
        """What share of annotated TSSs land in a peak.

        The mirror of ``SpliceJunctionCatalogue.verify_against``, and it exists
        for the same reason: a half-open slip or a contig-naming mismatch does
        not fail here, it silently returns "no evidence" for every model, which
        is indistinguishable from a genuine result.  A large majority of
        annotated starts should sit in a CAGE peak.  Near zero means the
        coordinates are wrong, not that the biology is surprising.
        """
        n = hit = near = 0
        for gid, g in ref.genes.items():
            if not g.tx_ids:
                continue
            tx = ref.tx[g.tx_ids[0]]
            ev = self.evidence(g.contig, tx.tss, g.strand)
            n += 1
            hit += bool(ev["in_cage_peak"])
            d = ev["dist_to_cage_reptss"]
            near += bool(d is not None and d <= 100)
            if max_genes and n >= max_genes:
                break
        frac = hit / n if n else 0.0
        return {
            "n_genes_tested": n,
            "annotated_tss_in_peak": hit,
            "frac_in_peak": round(frac, 4),
            "frac_within_100bp_of_reptss": round(near / n, 4) if n else 0.0,
            "ok": bool(n == 0 or frac >= 0.10),
            "diagnosis": (
                "" if n == 0 or frac >= 0.10 else
                "fewer than 10% of annotated TSSs fall inside a CAGE peak. "
                "Check that the atlas assembly matches the reference GTF "
                "(FANTOM5 ships mm10/GRCm38; GENCODE vM26+ is GRCm39) and that "
                "contig names agree."
            ),
        }
