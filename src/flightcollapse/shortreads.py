"""Short-read splice junctions from STAR, as an independent anchor.

Why this is not a filter
------------------------
The obvious use of short reads -- drop any long-read junction they do not
confirm -- is the wrong one here, for two reasons.

First, absence of support is weak evidence in this data.  The BD144b library
maps 24% uniquely with 70% of reads discarded as "too short", and an oligo-dT
long-read library is 3'-biased in a way a short-read library is not.  A junction
missing from the short reads may simply be in a region the short reads did not
cover.

Second, and more importantly: a set used to filter cannot also be used to
validate.  The calibration in :mod:`flightcollapse.scoring` currently anchors
its positives on *annotation*, which means the model can only ever learn "looks
like something already in GENCODE" -- precisely the wrong prior for finding
genuine novelty.  Short reads supply what annotation cannot: **17,337 junctions
in this sample that are supported by an orthogonal library and are not in the
reference.**  Those are true positives that are not annotated, and they are what
makes the fitted score something other than an annotation detector.

So the default use is:

* short-read-supported novel junctions **join the positive anchor set**
* they **leave the decoy set** -- a non-canonical motif with fifty short reads
  behind it is a real minor-class junction, not an alignment error
* **half of them are held out** by a deterministic hash, so the recall measured
  on that half is an honest estimate rather than the model grading its own work
* nothing is ever *dropped* for lacking short-read support

Optional rescue (``junctions.short_read_rescue``, off by default) additionally
exempts a short-read-supported junction from curation.  That one *does* change
which models are called, which is why it is opt-in.

The STAR format
---------------
``SJ.out.tab`` columns, all 1-based::

    1 contig                     5 intron motif (0 non-canonical, 1 GT/AG, ...)
    2 first base of the intron   6 0 novel / 1 annotated in the sjdb
    3 last base of the intron    7 uniquely-mapping reads crossing the junction
    4 strand (0 undef, 1 +, 2 -) 8 multi-mapping reads crossing the junction
                                 9 maximum spliced alignment overhang

This package stores junctions as ``(donor, acceptor)`` 0-based half-open, where
``donor`` is the first intronic base and ``acceptor`` the first exonic base
after the intron.  So ``donor = col2 - 1`` and ``acceptor = col3``.  That
off-by-one is the classic way to get a silent zero-support result, so
:meth:`SpliceJunctionCatalogue.verify_against` checks it rather than trusting
it: STAR's own column 6 says which junctions are annotated, and if the
conversion or the contig naming is wrong the agreement with the reference index
collapses to near zero instead of quietly returning "no support".
"""

from __future__ import annotations

import gzip
import hashlib
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

#: STAR intron-motif codes -> the motif as this package spells it
STAR_MOTIF = {
    0: "non-canonical", 1: "GT-AG", 2: "GT-AG", 3: "GC-AG",
    4: "GC-AG", 5: "AT-AC", 6: "AT-AC",
}
STAR_STRAND = {0: ".", 1: "+", 2: "-"}


def norm_contig(c: str) -> str:
    """``chr1`` and ``1`` are the same contig; ``chrM`` and ``MT`` are too.

    The user has hit this before: a reference FASTA and a barcode table that
    disagreed on the prefix. Normalising on both sides costs nothing and turns a
    silent zero-overlap into a non-issue.
    """
    c = str(c)
    if c.startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "MT") else c


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _holdout(strand: str, donor: int, acceptor: int, frac: float) -> bool:
    """Deterministic, position-derived split -- stable across contigs and runs.

    Must not depend on iteration order or on a random seed: the same junction
    has to land in the same half in every sample, or the held-out set stops
    being a held-out set the moment two samples are compared.
    """
    h = hashlib.blake2b(f"{strand}:{donor}:{acceptor}".encode(), digest_size=8).digest()
    return (int.from_bytes(h, "big") % 10_000) < int(frac * 10_000)


@dataclass
class SpliceJunctionCatalogue:
    """Summed STAR junctions from one or more libraries."""

    #: (norm_contig, strand, donor, acceptor) -> [uniq, multi, overhang, motif, sjdb]
    j: Dict[Tuple[str, str, int, int], List] = field(default_factory=dict)
    #: strand-agnostic fallback for STAR's strand-0 (undefined) junctions
    j_any: Dict[Tuple[str, int, int], List] = field(default_factory=dict)
    #: per-input-library unique-read counts, so that replicate concordance can be
    #: asked about separately from the sum
    lib_uniq: List[Dict[Tuple, int]] = field(default_factory=list)
    sources: List[str] = field(default_factory=list)
    n_rows: int = 0

    def __len__(self) -> int:
        return len(self.j) + len(self.j_any)

    @property
    def n_libraries(self) -> int:
        return max(len(self.lib_uniq), 1)

    def n_libs(self, key) -> int:
        """Replicates in which this junction was seen with a unique read.

        Summing replicates is right for *evidence*, but it destroys the
        information that matters most at low counts: a junction that reaches two
        reads only because both came from one library is a different object from
        one that appeared once in each of two.  These are technical replicates
        of the same cells, so the second reproduces and the first may not.
        """
        if not self.lib_uniq:
            return 1
        return sum(1 for d in self.lib_uniq if d.get(key, 0) > 0)

    # ------------------------------------------------------------------ #
    @classmethod
    def from_star(
        cls,
        paths: Sequence[str],
        min_overhang: int = 0,
        verbose: bool = True,
    ) -> "SpliceJunctionCatalogue":
        """Load and **sum** several ``SJ.out.tab`` files.

        Technical replicates of the same cells are summed: they are the same
        molecules sampled twice, so pooling them is the right way to get one
        junction's evidence rather than three under-powered opinions of it.
        """
        self = cls()
        for p in paths:
            if not os.path.exists(p):
                raise FileNotFoundError(f"short-read SJ file not found: {p}")
            self.sources.append(p)
            this_lib: Dict[Tuple, int] = {}
            self.lib_uniq.append(this_lib)
            n = 0
            with _open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 9:
                        continue
                    contig = norm_contig(f[0])
                    donor = int(f[1]) - 1          # 1-based first intronic base
                    acceptor = int(f[2])           # 1-based last intronic base
                    strand = STAR_STRAND.get(int(f[3]), ".")
                    motif = STAR_MOTIF.get(int(f[4]), "non-canonical")
                    sjdb = int(f[5])
                    uniq, multi, over = int(f[6]), int(f[7]), int(f[8])
                    if over < min_overhang:
                        continue
                    n += 1
                    if strand == ".":
                        d = self.j_any
                        k = (contig, donor, acceptor)
                    else:
                        d = self.j
                        k = (contig, strand, donor, acceptor)
                    cur = d.get(k)
                    if cur is None:
                        d[k] = [uniq, multi, over, motif, sjdb]
                    else:
                        cur[0] += uniq
                        cur[1] += multi
                        cur[2] = max(cur[2], over)
                        cur[4] = max(cur[4], sjdb)
                    if uniq > 0:
                        this_lib[k] = this_lib.get(k, 0) + uniq
            self.n_rows += n
            if verbose:
                print(f"[sj] {n:,} junctions from {os.path.basename(p)}")
        if verbose and len(paths) > 1:
            print(f"[sj] {len(self):,} distinct junctions after summing "
                  f"{len(paths)} libraries")
        return self

    # ------------------------------------------------------------------ #
    def lookup(self, contig: str, strand: str, donor: int, acceptor: int):
        """``[uniq, multi, overhang, motif, sjdb]`` or ``None``."""
        c = norm_contig(contig)
        hit = self.j.get((c, strand, int(donor), int(acceptor)))
        if hit is not None:
            return hit
        return self.j_any.get((c, int(donor), int(acceptor)))

    def support_arrays(
        self,
        contig: str,
        strands: Sequence[str],
        donors: Sequence[int],
        acceptors: Sequence[int],
    ) -> Dict[str, np.ndarray]:
        """Vectorised per-junction lookup for a whole contig's census."""
        n = len(donors)
        uniq = np.zeros(n, np.int64)
        multi = np.zeros(n, np.int64)
        over = np.zeros(n, np.int64)
        libs = np.zeros(n, np.int64)
        sjdb = np.zeros(n, bool)
        found = np.zeros(n, bool)
        c = norm_contig(contig)
        for i in range(n):
            key = (c, str(strands[i]), int(donors[i]), int(acceptors[i]))
            hit = self.j.get(key)
            if hit is None:
                hit = self.j_any.get((c, int(donors[i]), int(acceptors[i])))
                if hit is None:
                    continue
            else:
                libs[i] = self.n_libs(key)
            uniq[i], multi[i], over[i] = hit[0], hit[1], hit[2]
            sjdb[i] = bool(hit[4])
            found[i] = True
        return {
            "sr_uniq": uniq, "sr_multi": multi, "sr_overhang": over,
            "sr_n_libs": libs, "sr_in_sjdb": sjdb, "sr_found": found,
        }

    # ------------------------------------------------------------------ #
    def concordance(self, min_unique: int = 2) -> Dict[str, object]:
        """Are the single-replicate novel junctions just sampling, or artefacts?

        Technical replicates of the same cells should not disagree about which
        junctions exist -- except by sampling, which at two or three reads is a
        large effect.  So the observed count of junctions seen in exactly one
        replicate is compared against the multinomial expectation given each
        junction's own total and the libraries' relative spliced-read depth.

        The *excess* over that expectation is the part sampling cannot explain,
        and on this data it sits entirely at 2-5 reads and vanishes by 10.

        One caveat this cannot address: replicates catch stochastic noise, not
        systematic misalignment.  A junction produced by a pseudogene or a
        repeat reproduces perfectly in all three libraries.
        """
        if len(self.lib_uniq) < 2:
            return {}
        keys = [k for k, v in self.j.items() if not v[4] and v[0] >= min_unique]
        if not keys:
            return {}
        depth = np.array([sum(d.values()) for d in self.lib_uniq], float)
        p = depth / max(depth.sum(), 1.0)

        obs_one = 0
        exp_one = 0.0
        by_libs = {i: 0 for i in range(1, len(self.lib_uniq) + 1)}
        for k in keys:
            counts = np.array([d.get(k, 0) for d in self.lib_uniq], float)
            seen = int((counts > 0).sum())
            by_libs[seen] = by_libs.get(seen, 0) + 1
            if seen == 1:
                obs_one += 1
            exp_one += float((p ** counts.sum()).sum())
        excess = obs_one - exp_one
        return {
            "n_libraries": len(self.lib_uniq),
            "library_depth_share": [round(float(x), 3) for x in p],
            "novel_junctions_considered": len(keys),
            "by_n_replicates": by_libs,
            "seen_in_one_replicate": obs_one,
            "expected_by_sampling": round(exp_one, 1),
            "excess_over_sampling": round(excess, 1),
            "pct_of_set_unexplained": round(100 * excess / len(keys), 1),
        }

    def holdout_mask(
        self, strands: Sequence[str], donors: Sequence[int],
        acceptors: Sequence[int], frac: float,
    ) -> np.ndarray:
        return np.array(
            [_holdout(str(s), int(d), int(a), frac)
             for s, d, a in zip(strands, donors, acceptors)], bool
        )

    # ------------------------------------------------------------------ #
    def verify_against(self, ref, contigs: Optional[Sequence[str]] = None) -> Dict[str, object]:
        """Does the coordinate conversion actually land on annotated junctions?

        STAR's column 6 marks the junctions that were in the splice-junction
        database it was built with.  Those should also be annotated junctions in
        the reference index here.  If the off-by-one is wrong, or the contigs are
        named differently, this agreement goes to roughly zero -- and a zero here
        is the difference between "the short reads confirm nothing" and "the
        lookup never matched anything".  Worth one pass over 130k rows.
        """
        ref_by_norm: Dict[Tuple[str, str], set] = {}
        for (c, s), js in ref.junctions.items():
            ref_by_norm.setdefault((norm_contig(c), s), set()).update(js)
        want = {norm_contig(c) for c in contigs} if contigs else None

        n_sjdb = n_hit = 0
        n_novel = n_novel_hit = 0
        shifted: Dict[int, int] = {}
        for (c, s, d, a), v in self.j.items():
            if want is not None and c not in want:
                continue
            js = ref_by_norm.get((c, s))
            if js is None:
                continue
            hit = (d, a) in js
            if v[4]:
                n_sjdb += 1
                n_hit += hit
                if not hit:
                    # if we are off by one, the miss is systematic; count the
                    # offset that WOULD have hit so the error names itself
                    for off in (-1, 1):
                        if (d + off, a + off) in js:
                            shifted[off] = shifted.get(off, 0) + 1
                            break
            else:
                n_novel += 1
                n_novel_hit += hit

        frac = n_hit / max(n_sjdb, 1)
        out = {
            "sjdb_junctions": n_sjdb,
            "sjdb_found_in_reference": n_hit,
            "agreement": round(frac, 4),
            "novel_junctions": n_novel,
            "novel_also_in_reference": n_novel_hit,
            "systematic_offset_hits": shifted,
            "ok": n_sjdb > 0 and frac >= 0.8,
        }
        if n_sjdb < 100:
            out["note"] = (
                f"only {n_sjdb} STAR-annotated junctions were comparable; that is "
                "too few to conclude much either way"
            )
        if not out["ok"]:
            out["diagnosis"] = _diagnose(out, self, ref)
        return out


def _diagnose(res: Dict[str, object], sj: "SpliceJunctionCatalogue", ref) -> str:
    # Most specific first. A systematic shift is diagnostic on any number of
    # junctions, so it must be tested before the "nothing compared" branch --
    # otherwise a small run reports a contig problem it does not have.
    off = res.get("systematic_offset_hits") or {}
    if off:
        k = max(off, key=off.get)
        return (
            f"{off[k]:,} STAR-annotated junctions match the reference only after "
            f"shifting by {k:+d}. That is a coordinate-convention mismatch, not a "
            f"biological disagreement -- do not use these results."
        )
    sj_contigs = {k[0] for k in sj.j}
    ref_contigs = {norm_contig(c) for c, _ in ref.junctions}
    if not (sj_contigs & ref_contigs) or res["sjdb_junctions"] == 0:
        return (
            "almost no STAR junctions could even be compared. The contigs do not "
            f"line up: SJ file has {sorted(sj_contigs)[:5]}, reference has "
            f"{sorted(ref_contigs)[:5]}. Both are normalised by stripping 'chr', "
            "so this means a genuinely different assembly or reference build."
        )
    return (
        f"only {res['agreement']:.1%} of STAR's own annotated junctions are "
        "annotated in this reference GTF. Most likely the STAR index was built "
        "from a different annotation than the one passed to --reference-gtf. "
        "Short-read support will still work, but the sjdb flag is not comparable."
    )


# ---------------------------------------------------------------------- #
def holdout_summary(frames: Sequence["object"]) -> Dict[str, object]:
    """Genome-wide recall on junctions the model was never shown.

    This is the only number in the package that is not, in some form, the
    long-read data grading its own homework.  A held-out short-read-supported
    novel junction is known-real by evidence the calibration had no access to,
    so the fraction clearing an lFDR threshold is an honest recall.  Beside it
    sits the fraction of *unsupported* novel junctions clearing the same
    threshold: the gap between the two is what the score is actually buying.

    Pooled across contigs rather than averaged, so a contig with forty
    junctions does not count as much as chr1.
    """
    import pandas as pd

    if not len(frames):
        return {}
    df = pd.concat(list(frames), ignore_index=True)
    if "sr_holdout" not in df.columns or not df["sr_holdout"].any():
        return {}
    hold = df["sr_holdout"].to_numpy().astype(bool)
    ann = df["annotated"].to_numpy().astype(bool)
    sup = df["sr_supported"].to_numpy().astype(bool)
    unsup = ~ann & ~sup
    lf = df["lfdr"].to_numpy(float) if "lfdr" in df.columns else np.full(len(df), np.nan)

    out: Dict[str, object] = {
        "n_holdout_junctions": int(hold.sum()),
        "n_anchor_junctions": int(df.get("sr_anchor", pd.Series(dtype=bool)).sum()),
        "n_novel_unsupported": int(unsup.sum()),
        "holdout_kept_by_curation": round(float(df.loc[hold, "keep"].mean()), 4),
        "unsupported_novel_kept_by_curation": round(
            float(df.loc[unsup, "keep"].mean()), 4
        ) if unsup.any() else None,
    }
    h, u = lf[hold], lf[unsup]
    h, u = h[~np.isnan(h)], u[~np.isnan(u)]
    for t in (0.01, 0.05, 0.10):
        out[f"recall_at_lfdr_{t:g}"] = (
            round(float(np.mean(h <= t)), 4) if h.size else None
        )
        out[f"unsupported_kept_at_lfdr_{t:g}"] = (
            round(float(np.mean(u <= t)), 4) if u.size else None
        )
    return out


def evaluate_run(
    junctions_tsv: str,
    sj: SpliceJunctionCatalogue,
    min_sr_unique: int = 2,
) -> Tuple["object", "object"]:
    """Score an existing run's junction table against the short reads.

    Runs on ``{prefix}.junctions.tsv.gz`` from any previous collapse, so a
    finished run can be evaluated without collapsing it again.

    The number that matters is the last column: of the **novel** junctions the
    caller decided to keep, how many an orthogonal library also saw.  That is a
    precision estimate which nothing in the long-read pipeline could have
    manufactured.
    """
    import pandas as pd

    df = pd.read_csv(junctions_tsv, sep="\t")
    need = {"contig", "strand", "donor", "acceptor", "annotated", "keep"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{junctions_tsv} is missing columns: {sorted(missing)}")

    uniq = np.zeros(len(df), np.int64)
    for i, (c, s, d, a) in enumerate(
        zip(df.contig, df.strand, df.donor.astype(int), df.acceptor.astype(int))
    ):
        hit = sj.lookup(str(c), str(s), int(d), int(a))
        if hit is not None:
            uniq[i] = hit[0]
    df["sr_uniq"] = uniq
    df["sr_supported"] = uniq >= min_sr_unique

    rows = []
    for (ann, keep), sub in df.groupby(["annotated", "keep"]):
        rows.append({
            "annotated": bool(ann),
            "kept": bool(keep),
            "junctions": len(sub),
            "sr_supported": int(sub.sr_supported.sum()),
            "pct_sr_supported": round(100 * sub.sr_supported.mean(), 1),
            "median_sr_uniq": float(sub.sr_uniq.median()),
        })
    summary = pd.DataFrame(rows).sort_values(["annotated", "kept"], ascending=False)

    novel = df[~df.annotated.astype(bool)].copy()
    bands = []
    if len(novel):
        support = np.where(novel.get("n_mols", 0) > 0, novel.get("n_mols", 0), novel.n_reads)
        novel["lr_support"] = support
        for lo, hi in ((1, 2), (3, 4), (5, 9), (10, 19), (20, 49), (50, 10 ** 9)):
            sub = novel[(novel.lr_support >= lo) & (novel.lr_support <= hi)]
            if not len(sub):
                continue
            bands.append({
                "long_read_support": f"{lo}-{hi}" if hi < 10 ** 9 else f"{lo}+",
                "novel_junctions": len(sub),
                "kept": int(sub.keep.sum()),
                "sr_supported": int(sub.sr_supported.sum()),
                "pct_sr_supported": round(100 * sub.sr_supported.mean(), 1),
                "pct_sr_supported_among_kept": (
                    round(100 * sub[sub.keep].sr_supported.mean(), 1)
                    if sub.keep.any() else float("nan")
                ),
            })
    return summary, pd.DataFrame(bands)
