"""Cross-sample model concordance.

Three libraries built from the *same cDNA* are the strongest validation
available without a ground truth: a transcript model called in one and not the
others is either an artifact that a preparation step removed, or a false
positive of the caller. Reproducibility across them is therefore a usable stand
-in for precision, and — unlike anything internal to a single run — it is not
circular.

**Models are matched on structure, never on name.** ``GAPDH|NIC_3`` is a
per-gene counter assigned in emission order, so the same name in two samples is
two different transcripts. The key here is the exon structure itself:

    (contig, strand, intron chain within `fuzzy`, 3' end within `max_3p_diff`)

The 5' end is deliberately excluded. It is a percentile statistic over a
5'-degraded read pile, so it drifts by hundreds of bases between libraries of
the same molecules; including it would report real reproducibility as noise.

Reads any tool's GTF/GFF, so ``isoseq collapse`` output goes in on the same
footing and the comparison is like for like.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .intervals import Chain, Exon, exons_to_chain


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _attr(field: str, key: str) -> Optional[str]:
    i = field.find(f'{key} "')
    if i < 0:
        return None
    j = i + len(key) + 2
    k = field.find('"', j)
    return field[j:k] if k > 0 else None


#: columns a pigeon / SQANTI3 classification table is recognised by
_PIGEON_ID = ("isoform", "transcript_id", "pbid", "id")
_PIGEON_CAT = ("structural_category", "category", "class")


def load_categories(path: str, verbose: bool = True) -> Dict[str, Tuple[str, float]]:
    """``transcript_id -> (category, support)`` from another tool's own table.

    A collapse GFF carries structure and nothing else -- ``isoseq collapse``
    writes no category and no counts, because the classification happens later
    in ``pigeon classify``.  This reads that classification back so the same
    comparison can be drawn for any pipeline, using **its own** vocabulary
    rather than one imposed on it.

    Accepts a pigeon/SQANTI3 ``*_classification.txt`` (recognised by an
    ``isoform`` and a ``structural_category`` column) or any TSV whose first
    column is the id and second the category, optionally with a third for
    support.
    """
    # A pigeon classification for an unfiltered collapse runs to hundreds of MB
    # across fifty-odd columns. Read the header first, decide which three are
    # wanted, and load only those -- otherwise three samples of this is tens of
    # gigabytes of strings nobody looks at.
    head = pd.read_csv(path, sep="\t", nrows=0)
    cols = {c.lower(): c for c in head.columns}
    idc = next((cols[c] for c in _PIGEON_ID if c in cols), head.columns[0])
    catc = next((cols[c] for c in _PIGEON_CAT if c in cols), head.columns[1])
    # pigeon writes the full-length count as FL.<sample>; SQANTI3 as FL
    supc = next(
        (c for c in head.columns
         if c == "FL" or (c.startswith("FL.") and not c.startswith("FL_TPM"))),
        None,
    )
    if supc is None and len(head.columns) > 2 and catc != head.columns[2]:
        supc = head.columns[2]
    want = [idc, catc] + ([supc] if supc else [])
    df = pd.read_csv(path, sep="\t", usecols=want, low_memory=False)
    if supc is not None and not pd.api.types.is_numeric_dtype(df[supc]):
        df = df.drop(columns=[supc])
        supc = None
    sup = df[supc].astype(float) if supc is not None else pd.Series(0.0, index=df.index)
    out = {str(i): (str(c), float(s))
           for i, c, s in zip(df[idc], df[catc], sup)}
    if verbose:
        print(f"[compare] {len(out):,} categories from {path}"
              f"{f' (support from {supc})' if supc else ' (no support column)'}")
    return out


def load_models(
    path: str,
    verbose: bool = True,
    categories: Optional[Dict[str, Tuple[str, float]]] = None,
    min_support: float = 0.0,
    support_attr: Optional[str] = None,
) -> pd.DataFrame:
    """One row per transcript in a GTF/GFF, with its exon structure.

    ``min_support`` is applied **while parsing** whenever the support is already
    known -- from ``categories``, or from a GTF attribute. An unfiltered
    ``isoseq collapse`` GFF holds 2.3M transcripts for one sample; building the
    exon table for all of them and then discarding 94% costs several GB per
    sample for nothing.

    ``support_attr`` picks which GTF attribute counts as support. It matters
    when comparing across pipelines: ``isoseq`` reports FL *reads*, so a
    like-for-like threshold has to read ``n_reads`` here rather than the
    UMI-collapsed ``n_molecules``.
    """
    attrs = ((support_attr,) if support_attr else ("n_molecules", "n_reads"))
    ex: Dict[str, List[Exon]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str, str, str, str]] = {}
    rejected: set = set()
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = _attr(f[8], "transcript_id")
            if tid is None:
                continue
            if min_support > 0:
                if tid in rejected:
                    continue
                if tid not in ex:          # first exon line: decide once
                    if categories is not None:
                        sup = categories.get(tid, ("", 0.0))[1]
                    else:
                        sup = float(next(
                            (v for v in (_attr(f[8], a) for a in attrs) if v), 0))
                    if sup < min_support:
                        rejected.add(tid)
                        continue
            ex[tid].append((int(f[3]) - 1, int(f[4])))
            if tid not in meta:
                meta[tid] = (
                    f[0], f[6],
                    _attr(f[8], "category") or "",
                    _attr(f[8], "gene_name") or _attr(f[8], "gene_id") or "",
                    next((v for v in (_attr(f[8], a) for a in attrs) if v), "0"),
                )
    rows = []
    for tid, e in ex.items():
        e.sort()
        contig, strand, cat, gene, sup = meta[tid]
        rows.append(
            dict(model_id=tid, contig=contig, strand=strand, category=cat, gene=gene,
                 support=float(sup), n_exons=len(e),
                 start=e[0][0], end=e[-1][1],
                 three=e[-1][1] if strand == "+" else e[0][0],
                 chain=exons_to_chain(e, strand))
        )
    df = pd.DataFrame(rows)
    if verbose and rejected:
        print(f"[compare] {len(rejected):,} models below support {min_support:g} "
              f"were skipped while parsing")
    if categories:
        hit = df.model_id.map(lambda t: categories.get(t))
        found = hit.notna()
        df.loc[found, "category"] = [v[0] for v in hit[found]]
        df.loc[found, "support"] = [v[1] for v in hit[found]]
        df.loc[~found, "category"] = df.loc[~found, "category"].replace("", "unclassified")
        if verbose:
            print(f"[compare] {int(found.sum()):,}/{len(df):,} models matched a "
                  f"category ({100 * found.mean():.1f}%)")
    if verbose:
        print(f"[compare] {len(df):,} models from {path}")
    return df


def _chain_key(chain: Chain, fuzzy: int) -> Tuple:
    """Round junction coordinates so that near-identical chains collide."""
    if not chain:
        return ()
    if fuzzy <= 0:
        return chain
    return tuple((d // fuzzy, a // fuzzy) for d, a in chain)


def _match_block(three: np.ndarray, sample: np.ndarray, tol: int) -> np.ndarray:
    """Assign one structural block to groups, at most one model per sample.

    This replaces a single-linkage chain plus a patch, and the patch was the
    problem. Previously the block was cut wherever consecutive 3' ends were more
    than ``tol`` apart, and any sample appearing twice in the resulting group had
    its extra models *forced into unique singleton groups*. That makes a looser
    tolerance produce MORE structures and fewer matches -- the opposite of what
    loosening should do -- because merging two neighbouring sites converts one
    honest pair of groups into one group plus a handful of manufactured
    singletons. Any comparison between two callers that fragment 3' ends to
    different degrees was reading that artefact as a difference in the data.

    Instead: repeatedly take the position covering the most distinct samples
    within +/- ``tol``, claim the nearest model from each, and continue. Same
    greedy-peak idea as the 3'-end clustering in ``ends.py``, for the same
    reason -- a cleavage site should not be split by where a linkage chain
    happens to break.
    """
    n = three.size
    grp = np.full(n, -1, np.int64)
    if n == 1:
        grp[0] = 0
        return grp
    # ``three`` is sorted within the block, so each sample's own positions are
    # sorted too and the nearest neighbour is one searchsorted away. Scoring
    # every candidate centre against every sample this way is O(S n log n) per
    # round; the obvious nested-loop version was O(n^2 k) and did not finish on
    # a whole-genome comparison at a wide tolerance.
    samples = list(dict.fromkeys(sample.tolist()))
    where = {s: np.flatnonzero(sample == s) for s in samples}
    alive = np.ones(n, bool)
    g = 0
    while alive.any():
        cand = np.flatnonzero(alive)
        pos = three[cand]
        cover = np.zeros(cand.size, np.int64)
        cost = np.zeros(cand.size, float)
        nearest: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for s in samples:
            av = where[s][alive[where[s]]]
            if av.size == 0:
                continue
            sp = three[av]
            j = np.searchsorted(sp, pos)
            jl = np.clip(j - 1, 0, sp.size - 1)
            jr = np.clip(j, 0, sp.size - 1)
            dl, dr = np.abs(sp[jl] - pos), np.abs(sp[jr] - pos)
            left = dl <= dr
            best = np.where(left, av[jl], av[jr])
            dist = np.where(left, dl, dr)
            ok = dist <= tol
            cover += ok
            cost += np.where(ok, dist, 0)
            nearest[s] = (best, ok)
        # most samples covered, then tightest, then leftmost for determinism
        k = int(np.lexsort((pos, cost, -cover))[0])
        for s, (best, ok) in nearest.items():
            if ok[k]:
                i = int(best[k])
                grp[i] = g
                alive[i] = False
        g += 1
    return grp


def _match_block_sweep(three: np.ndarray, sample: np.ndarray, tol: int) -> np.ndarray:
    """Left-to-right variant for very large blocks.

    Same guarantee -- at most one model per sample per group -- but seeded on
    the leftmost unassigned model instead of the densest position, which makes
    it linear in the block size. Slightly worse pairing, vastly cheaper, and
    only used where the densest-first search would not finish.
    """
    n = three.size
    grp = np.full(n, -1, np.int64)
    used = np.zeros(n, bool)
    samples = list(dict.fromkeys(sample.tolist()))
    where = {s: np.flatnonzero(sample == s) for s in samples}
    pos = {s: three[where[s]] for s in samples}
    g = 0
    for i in range(n):
        if used[i]:
            continue
        centre = three[i]
        used[i] = True
        grp[i] = g
        for s in samples:
            if s == sample[i]:
                continue
            arr, sp = where[s], pos[s]
            lo = int(np.searchsorted(sp, centre - tol, side="left"))
            hi = int(np.searchsorted(sp, centre + tol, side="right"))
            best, bd = -1, tol + 1
            for j in range(lo, hi):
                k = int(arr[j])
                if used[k]:
                    continue
                d = abs(int(three[k]) - int(centre))
                if d < bd:
                    bd, best = d, k
            if best >= 0:
                used[best] = True
                grp[best] = g
        g += 1
    return grp


def _match_all(A: pd.DataFrame, tol: int, verbose: bool = True) -> np.ndarray:
    """Run :func:`_match_block` over every structural block of a sorted frame."""
    keys = A["_k"].values
    three = A["three"].to_numpy(np.int64)
    sample = A["sample"].to_numpy(object)
    out = np.empty(len(A), np.int64)
    start = 0
    g = 0
    big = 0
    for i in range(1, len(A) + 1):
        if i < len(A) and keys[i] == keys[start]:
            continue
        sl = slice(start, i)
        n = i - start
        if n > 2000:
            # Mono-exonic models share the empty chain, so every one of them on
            # a contig lands in a single block -- tens of thousands of rows.
            # The densest-first search is quadratic in that; sweep instead.
            # (This branch used to assign the whole block to ONE group, which
            # silently merged thousands of unrelated models into one
            # "reproduced" structure. It fired twice on the 0.1.10 run.)
            big += 1
            sub = _match_block_sweep(three[sl], sample[sl], tol)
        else:
            sub = _match_block(three[sl], sample[sl], tol)
        out[sl] = sub + g
        g += int(sub.max()) + 1
        start = i
    if verbose and big:
        print(f"[compare] {big} large block(s) matched by sweep rather than "
              f"densest-first")
    return out


def concordance(
    paths: Dict[str, str],
    fuzzy: int = 5,
    max_3p_diff: int = 100,
    verbose: bool = True,
    categories: Optional[Dict[str, str]] = None,
    min_support: float = 0.0,
    support_attr: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Match models across samples; return (per-structure table, summary).

    ``categories`` optionally maps a sample name to a classification table, so
    another pipeline's own categories and counts can be attached to its GFF.
    """
    cats = {k: load_categories(v, verbose) for k, v in (categories or {}).items()}
    frames = []
    for name, p in paths.items():
        d = load_models(p, verbose, categories=cats.get(name),
                        min_support=min_support, support_attr=support_attr)
        d["sample"] = name
        frames.append(d)
    A = pd.concat(frames, ignore_index=True)
    A["ckey"] = [_chain_key(c, fuzzy) for c in A["chain"]]

    A["_k"] = list(zip(A.contig, A.strand, A.ckey, A.n_exons))
    A = A.sort_values(["_k", "three"], kind="mergesort").reset_index(drop=True)
    A["grp"] = _match_all(A, max_3p_diff, verbose)

    g = A.groupby("grp")
    C = pd.DataFrame({
        "n_samples": g["sample"].nunique(),
        "samples": g["sample"].apply(lambda s: ",".join(sorted(s))),
        "category": g["category"].agg(lambda s: s.mode().iat[0] if len(s) else ""),
        "gene": g["gene"].first(),
        "model_id": g["model_id"].first(),
        "contig": g["contig"].first(),
        "strand": g["strand"].first(),
        "three": g["three"].first(),
        "n_exons": g["n_exons"].first(),
        "mean_support": g["support"].mean(),
        "min_support": g["support"].min(),
        #: how far apart the samples put the same structure's 3' end. This is
        #: the number that separates "the caller resolved a real second
        #: cleavage site" from "the caller split one site at a different offset
        #: in each library": a real site reproduces to within the cleavage
        #: heterogeneity (~20-30 bp), a split does not.
        "three_spread": g["three"].max() - g["three"].min(),
    }).reset_index(drop=True)

    n = len(paths)

    def _row(label, sub):
        d = {"category": label, "structures": len(sub)}
        for k in range(n, 0, -1):
            d[f"in_{k}"] = int((sub.n_samples == k).sum())
        d["pct_reproduced"] = round(100 * (sub.n_samples == n).mean(), 1)
        full = sub[sub.n_samples == n]
        d["median_3p_spread"] = (
            round(float(full.three_spread.median()), 1) if len(full) else float("nan")
        )
        return d

    rows = [_row(cat, sub) for cat, sub in C.groupby("category")]
    rows.append(_row("ALL", C))
    return C, pd.DataFrame(rows)


def reproducibility_curve(
    C: pd.DataFrame, n_samples: int, col: str = "mean_support",
    thresholds: Sequence[float] = (0, 3, 5, 10, 20, 50, 100, 200, 500),
) -> pd.DataFrame:
    """How reproducibility trades against yield as a support threshold rises.

    This is the curve that sets the threshold empirically: replication stands in
    for truth, so the knee is where extra models stop being reproducible.
    """
    rows = []
    for t in thresholds:
        sub = C[C[col] >= t]
        if sub.empty:
            continue
        rows.append({
            "min_support": t,
            "structures": len(sub),
            "reproduced_all": int((sub.n_samples == n_samples).sum()),
            "pct_reproduced": round(100 * (sub.n_samples == n_samples).mean(), 1),
            "singletons": int((sub.n_samples == 1).sum()),
            "pct_singleton": round(100 * (sub.n_samples == 1).mean(), 1),
        })
    return pd.DataFrame(rows)
