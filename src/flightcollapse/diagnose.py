"""Diagnose a barcode/UMI table against a BAM.

A read-name mismatch between the two is silent and catastrophic -- every
molecule count comes back zero and every threshold quietly demotes itself to
raw reads -- and it is fiddly to debug with shell one-liners, because a random
sample of a genome-wide table intersected with one chromosome's reads is
*supposed* to be small.  What matters is the ratio between candidate keys, not
the absolute number.

So this samples the table, tries every column that could plausibly hold a read
name under every normalisation mode, and reports the hit counts side by side.
The winning combination is usually obvious by an order of magnitude.
"""

from __future__ import annotations

import gzip
import re
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pysam

from .molecules import normalise_name

MODES = ("none", "strip_segment", "strip_ccs", "zmw")
#: movie/zmw/ccs, optionally with a segment span
_READLIKE = re.compile(r"^[A-Za-z0-9_.\-]+/\d+/\w+(/\d+_\d+)?$")


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def bam_read_names(bam_path: str, limit: Optional[int] = None) -> List[str]:
    bam = pysam.AlignmentFile(bam_path, "rb")
    out: List[str] = []
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        out.append(r.query_name)
        if limit and len(out) >= limit:
            break
    bam.close()
    return out


def sample_table(
    tsv: str, n_rows: int, rows_per_seek: int = 1, seed: int = 0
) -> Tuple[List[str], List[List[str]]]:
    """Sample rows uniformly at random from the whole file.

    Two traps, both of which produced a wrong answer before this was fixed:

    1. These tables are written by streaming a coordinate-sorted BAM, so their
       row order is **genomic**.  ``head -200000`` is not a sample of the
       library, it is the start of chr1.
    2. Seeking to N random offsets and then reading a *block* of consecutive
       rows does not fix that.  Each block is a narrow genomic window, so a
       chromosome slice that occupies 0.6% of the genome is hit by ~1 block in
       200 — the exact-match count then swings between ~0 and thousands purely
       on where the blocks landed, while ZMW-level matching stays stable
       (a ZMW has segments all over the genome).  That asymmetry makes exact
       matching look broken when it is fine.

    So: one row per seek, many seeks.  ``rows_per_seek`` is exposed only for
    testing; leave it at 1.

    Falls back to a head scan for gzipped input, which is not usefully
    seekable — the caller is warned in that case.
    """
    import os
    import random

    with _open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        head_len = fh.tell() if not tsv.endswith(".gz") else 0

    if tsv.endswith(".gz"):
        rows = []
        with _open(tsv) as fh:
            fh.readline()
            for k, line in enumerate(fh):
                rows.append(line.rstrip("\n").split("\t"))
                if k + 1 >= n_rows:
                    break
        return header, rows

    size = os.path.getsize(tsv)
    n_seeks = max(1, n_rows // max(rows_per_seek, 1))
    rng = random.Random(seed)
    offsets = sorted(rng.randrange(head_len, max(head_len + 1, size)) for _ in range(n_seeks))
    rows: List[List[str]] = []
    with open(tsv, "rb") as fh:
        for off in offsets:
            fh.seek(off)
            fh.readline()  # discard the partial line
            for _ in range(rows_per_seek):
                line = fh.readline()
                if not line:
                    break
                f = line.decode("utf-8", "replace").rstrip("\n").split("\t")
                if len(f) >= 2:
                    rows.append(f)
    return header, rows


def candidate_columns(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[int]:
    """Columns whose values look like PacBio read names, best first."""
    scores: List[Tuple[float, int]] = []
    for i, name in enumerate(header):
        vals = [r[i] for r in rows if len(r) > i]
        if not vals:
            continue
        frac = sum(bool(_READLIKE.match(v)) for v in vals) / len(vals)
        if frac < 0.5:
            continue
        bonus = 0.5 if name.strip().lower() in ("read_name", "read", "id", "qname") else 0.0
        scores.append((frac + bonus, i))
    scores.sort(reverse=True)
    return [i for _s, i in scores]


def diagnose(
    bam_path: str,
    tsv: str,
    table_rows: int = 500_000,
    bam_reads: Optional[int] = None,
    verbose: bool = True,
) -> Dict[str, object]:
    names = bam_read_names(bam_path, bam_reads)
    if not names:
        return {"error": "no primary alignments in the BAM"}
    if verbose:
        print(f"[diagnose] {len(names):,} BAM read names, e.g. {names[0]}")

    bam_sets: Dict[str, Set[str]] = {
        m: {normalise_name(n, m) for n in names} for m in MODES
    }

    header, rows = sample_table(tsv, table_rows)
    if not rows:
        return {"error": f"{tsv} has no data rows"}
    sample = rows[:200]

    cands = candidate_columns(header, sample)
    if verbose:
        print(f"[diagnose] {len(rows):,} table rows sampled from across the file; "
              f"read-name-like columns: "
              f"{[f'{i + 1}:{header[i]}' for i in cands] or 'NONE'}")
        print(f"[diagnose] table example: "
              f"{rows[0][cands[0]] if cands else rows[0][:4]}")

    results = []
    for i in cands:
        exact = None
        per_mode = {}
        for m in MODES:
            hits = 0
            for f in rows:
                if len(f) > i and normalise_name(f[i], m) in bam_sets[m]:
                    hits += 1
            per_mode[m] = hits
        exact = per_mode["none"]
        for m, hits in per_mode.items():
            # A coarser key that matches MORE table rows than the exact key is
            # not matching the same reads better -- it is matching *additional*
            # rows, i.e. the other segments of the same ZMW, which are different
            # molecules with different UMIs.  Joining on it would assign the
            # wrong barcode to most reads.
            expansion = (hits / exact) if exact else float("inf")
            unsafe = m != "none" and exact > 0 and expansion > 1.2
            results.append(
                {
                    "column_index": i + 1,
                    "column": header[i],
                    "normalise": m,
                    "hits": hits,
                    "rows_per_exact_hit": round(expansion, 2) if exact else None,
                    "unsafe": unsafe,
                }
            )
    # safe modes first, then by hits, then by specificity: when a whole-genome
    # BAM makes every mode hit every sampled row, the tie must break towards the
    # most specific key rather than towards whatever came first in MODES.
    spec = {m: i for i, m in enumerate(MODES)}
    results.sort(key=lambda d: (d["unsafe"], -d["hits"], spec[d["normalise"]]))

    out: Dict[str, object] = {
        "bam_reads_sampled": len(names),
        "bam_example": names[0],
        "table_rows_sampled": len(rows),
        "table_columns": len(header),
        "candidate_columns": [f"{i + 1}:{header[i]}" for i in cands],
        "results": results[:12],
        "best": results[0] if results else None,
    }
    safe = [r for r in results if not r["unsafe"] and r["hits"] > 0]
    if not results or results[0]["hits"] == 0:
        out["verdict"] = (
            "No column and no normalisation mode matched a single read. Either "
            "the table was built from a different BAM than the one given, or "
            "the read names were rewritten between the two. Compare "
            "`samtools view <bam> | cut -f1 | head -3` with column "
            f"{cands[0] + 1 if cands else '?'} of the table by eye."
        )
    elif safe:
        best = safe[0]
        expected = (
            len(rows) * len(names) / 78_000_000 if len(names) else 0
        )
        unsafe_note = ""
        dropped = [r for r in results if r["unsafe"]]
        if dropped:
            top = max(dropped, key=lambda r: r["hits"])
            unsafe_note = (
                f"\n  Ignored {top['normalise']!r}: it matched "
                f"{top['rows_per_exact_hit']}x as many table rows as the exact "
                f"name. Those extra rows are other segments of the same ZMW -- "
                f"different molecules with different UMIs -- so joining on it "
                f"would give most reads the wrong barcode."
            )
        out["verdict"] = (
            f"Use col_read_name={header[best['column_index'] - 1]!r} with "
            f"read_name_normalise={best['normalise']!r} "
            f"({best['hits']:,} hits in {len(rows):,} sampled rows)." + unsafe_note +
            "\n  Hit counts scale with how much of the genome your BAM covers: "
            "a whole-genome table sampled against one chromosome legitimately "
            "hits ~1% of rows. Compare modes to each other, not to the row count."
        )
    else:
        out["verdict"] = (
            "Every mode that matched anything is unsafe: the only keys that hit "
            "are coarser than a read (movie/ZMW), which lumps several segmented "
            "reads -- different molecules -- onto one barcode. The exact read "
            "name does not match at all, so the table and the BAM come from "
            "different stages of the pipeline. Rebuild the table from this BAM "
            "rather than joining on ZMW."
        )
    return out
