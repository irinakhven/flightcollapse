"""Short-read coverage as terminal and length evidence.  New in 0.3.0.

Junctions -- which :mod:`flightcollapse.shortreads` already reads -- say which
introns exist.  Coverage says something no junction table can: how the signal
*behaves* across a terminal boundary, and whether a long transcript's distal
exons are expressed at all.

Three measurements, and it is worth being precise about what each can and
cannot establish.

``three_prime_step``
    A real cleavage site has a sharp coverage drop across it: transcription
    products end there.  An internal-priming site does not, because the
    transcript continues and the short reads keep covering it.  This is the
    only signal in the package that distinguishes those two cases by
    *mechanism* rather than by genomic context.

``tss_ratio``
    Downstream over upstream coverage at a putative start -- SQANTI3's
    ``ratio_TSS``, and the measurement the original review asked for.

``continuity``
    Mean coverage across a set of exons.  Where long reads stop short of an
    annotated long transcript but coverage runs on at comparable depth, the
    longer form is probably real.

The boundary that matters for the third one: coverage establishes that distal
exons are *expressed*.  It cannot establish **linkage** -- that those exons sit
on the same molecule as the proximal ones.  Connectivity is the one thing only
a long read provides.  So continuity may RESCUE an annotated long transcript
that the long reads under-support; it must never build a chain no read spans.
This module therefore reports numbers and never proposes structure.

Strandedness
------------
Both short-read libraries measured for this project are **unstranded** (the
bulk organoid alignment and the FLASH-seq plates both give STAR forward and
reverse counts at ~0.57 of unstranded).  At a locus overlapping an antisense
gene, unstranded coverage cannot be attributed, and a bidirectional promoter
makes the TSS ratio meaningless specifically.  Those loci are flagged and left
unscored rather than scored wrongly -- see ``CoverageParams.skip_antisense_overlap``.

Dependency
----------
``pyBigWig`` is optional.  Without it, or without a track, every method returns
``None`` and every downstream feature is absent rather than zero: "no opinion"
and "no signal" are different statements and the tables keep them apart.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import CoverageParams

#: Set once at import so a missing dependency is reported, not raised per call.
try:  # pragma: no cover - environment dependent
    import pyBigWig  # type: ignore

    HAVE_PYBIGWIG = True
except ImportError:  # pragma: no cover - environment dependent
    pyBigWig = None
    HAVE_PYBIGWIG = False


class CoverageTrack:
    """One or more bigWig tracks, summed, with strand-agnostic lookups.

    Several tracks are summed rather than kept apart because depth is the point
    here: these are replicate libraries of the same material, and one
    under-powered opinion per replicate is worse than their sum.  That is the
    opposite of the decision made for *libraries of different kinds*, which
    stay separate channels -- a bulk library and a plate library answer
    different questions and must not be pooled.
    """

    __slots__ = ("_bw", "sources", "_contigs", "unstranded", "n_failed")

    def __init__(self) -> None:
        self._bw: List[object] = []
        self.sources: List[str] = []
        self._contigs: Dict[str, str] = {}
        self.unstranded = True
        self.n_failed = 0

    # ------------------------------------------------------------------ #
    @classmethod
    def open(
        cls,
        paths: Sequence[str],
        unstranded: bool = True,
        verbose: bool = True,
    ) -> Optional["CoverageTrack"]:
        """Return a track, or ``None`` when coverage is unavailable.

        ``None`` rather than an exception: a run without coverage is an
        ordinary run with one fewer evidence channel, and the QC report says
        which channels were active.
        """
        if not paths:
            return None
        if not HAVE_PYBIGWIG:
            if verbose:
                print(
                    "[coverage] pyBigWig is not installed, so the coverage "
                    "channel is inactive. `pip install pyBigWig` to enable "
                    "the 3' step, TSS ratio and continuity features."
                )
            return None
        self = cls()
        self.unstranded = unstranded
        for p in paths:
            if not os.path.exists(p):
                if verbose:
                    print(f"[coverage] missing track, skipping: {p}")
                self.n_failed += 1
                continue
            try:
                bw = pyBigWig.open(p)
            except Exception as exc:  # pragma: no cover - file dependent
                if verbose:
                    print(f"[coverage] could not open {p}: {exc}")
                self.n_failed += 1
                continue
            self._bw.append(bw)
            self.sources.append(p)
            for c in bw.chroms():
                self._contigs.setdefault(_norm(c), c)
        if not self._bw:
            return None
        if verbose:
            print(
                f"[coverage] {len(self._bw)} track(s), "
                f"{len(self._contigs):,} contigs, "
                f"{'unstranded' if unstranded else 'stranded'}"
            )
        return self

    def close(self) -> None:
        for bw in self._bw:
            try:
                bw.close()
            except Exception:  # pragma: no cover
                pass

    def __len__(self) -> int:
        return len(self._bw)

    # ------------------------------------------------------------------ #
    def _resolve(self, contig: str) -> Optional[str]:
        """Map a reference contig name onto the track's own naming.

        Both short-read libraries here are Ensembl-style (``1``) while the
        GENCODE references are UCSC-style (``chr1``), so this is not a
        hypothetical: without it every lookup silently returns nothing, which
        is indistinguishable from a locus with no coverage.
        """
        return self._contigs.get(_norm(contig))

    def mean(self, contig: str, start: int, end: int) -> Optional[float]:
        """Mean coverage over ``[start, end)``; ``None`` if out of range."""
        name = self._resolve(contig)
        if name is None or end <= start:
            return None
        total = 0.0
        seen = False
        for bw in self._bw:
            try:
                length = bw.chroms().get(name)
                if length is None:
                    continue
                lo, hi = max(0, int(start)), min(int(length), int(end))
                if hi <= lo:
                    continue
                v = bw.stats(name, lo, hi, type="mean")[0]
            except Exception:  # pragma: no cover - file dependent
                continue
            if v is not None:
                total += float(v)
                seen = True
        return total if seen else None

    # ------------------------------------------------------------------ #
    def three_prime_step(
        self, contig: str, pos: int, strand: str, params: CoverageParams
    ) -> Dict[str, object]:
        """Coverage either side of a putative cleavage site.

        A real 3' end: high upstream, low downstream.  An internal-priming
        site: similar on both sides, because the transcript did not end there.
        """
        w = params.window
        up, down = _flanks(pos, strand, w, three_prime=True)
        cu = self.mean(contig, *up)
        cd = self.mean(contig, *down)
        out: Dict[str, object] = {
            "sr_cov_upstream": None if cu is None else round(cu, 3),
            "sr_cov_downstream": None if cd is None else round(cd, 3),
            "sr_3p_step_ratio": None,
            "sr_3p_step_supported": None,
        }
        if cu is None or cd is None or cu < params.min_coverage:
            return out
        ratio = cd / cu if cu > 0 else None
        if ratio is None:
            return out
        out["sr_3p_step_ratio"] = round(ratio, 4)
        out["sr_3p_step_supported"] = bool(ratio <= params.max_downstream_ratio)
        return out

    def tss_ratio(
        self, contig: str, pos: int, strand: str, params: CoverageParams
    ) -> Dict[str, object]:
        """SQANTI3-style ratio: signal downstream of a start over upstream."""
        w = params.window
        up, down = _flanks(pos, strand, w, three_prime=False)
        cu = self.mean(contig, *up)
        cd = self.mean(contig, *down)
        out: Dict[str, object] = {
            "sr_tss_ratio": None,
            "sr_tss_supported": None,
        }
        if cu is None or cd is None or cd < params.min_coverage:
            return out
        ratio = cd / max(cu, 1e-6)
        out["sr_tss_ratio"] = round(float(min(ratio, 999.0)), 4)
        out["sr_tss_supported"] = bool(ratio >= params.min_tss_ratio)
        return out

    def continuity(
        self, contig: str, exons: Sequence[Tuple[int, int]]
    ) -> Optional[float]:
        """Mean coverage across a model's exons, length-weighted.

        Used to compare a model against a longer annotated transcript at the
        same locus.  It answers "are the distal exons expressed", never "are
        they on this molecule".
        """
        tot = 0.0
        length = 0
        for s, e in exons:
            v = self.mean(contig, s, e)
            if v is None:
                continue
            tot += v * (e - s)
            length += e - s
        return (tot / length) if length else None


# ---------------------------------------------------------------------- #
def _norm(contig: str) -> str:
    c = str(contig)
    if c.startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "MT") else c


def _flanks(
    pos: int, strand: str, window: int, three_prime: bool
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """``(upstream, downstream)`` intervals in TRANSCRIPT orientation.

    Getting this wrong is the silent failure mode of the whole module: the
    windows would simply be swapped on one strand and the ratio inverted for
    half the genome, with no error anywhere. Hence one function, used by both
    measurements, and a test that pins both strands.
    """
    pos = int(pos)
    if strand == "+":
        before = (max(0, pos - window), pos)
        after = (pos, pos + window)
    else:
        before = (pos, pos + window)
        after = (max(0, pos - window), pos)
    return (before, after)


def blank_coverage_evidence() -> Dict[str, object]:
    """The keys every model carries, so absent and zero stay distinguishable."""
    return {
        "sr_cov_upstream": "",
        "sr_cov_downstream": "",
        "sr_3p_step_ratio": "",
        "sr_3p_step_supported": "",
        "sr_tss_ratio": "",
        "sr_tss_supported": "",
        "sr_cov_continuity": "",
    }
