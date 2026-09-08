"""Standalone validation harness for *any* collapse output.

This is the diagnostic that settled the isoseq question, generalised: give it a
BAM, a collapsed GTF/GFF and a ``read_stat.txt``, and it reports what the reads
assigned to each model actually look like.  It makes no assumption about which
tool produced the GTF, so it works equally on ``isoseq collapse``,
``flightcollapse``, bambu or StringTie output and is the right way to compare
them.

Reported per model class:

* median exons in the model vs median exons among its assigned reads
* fraction of assigned reads that are mono-exonic
* median ``log2(read span / model span)``
* ``model_is_min`` -- is the model the shortest member of its own read group?
* median 5'/3' offsets in transcript orientation

Plus the collapse-independent step-0 number: the fraction of aligned reads that
are genuinely unspliced, straight from the CIGARs.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam

from .intervals import Exon
from .reads import exon_blocks, introns_of, resolve_contig


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_models(path: str) -> Dict[str, Dict[str, object]]:
    ex: Dict[str, List[Exon]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str]] = {}
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = f[8]
            i = a.find('transcript_id "')
            if i < 0:
                continue
            tid = a[i + 15: a.find('"', i + 15)]
            ex[tid].append((int(f[3]) - 1, int(f[4])))
            meta.setdefault(tid, (f[0], f[6]))
    out = {}
    for tid, e in ex.items():
        e.sort()
        contig, strand = meta[tid]
        out[tid] = {
            "contig": contig,
            "strand": strand,
            "exons": e,
            "n_exons": len(e),
            "start": e[0][0],
            "end": e[-1][1],
            "span": e[-1][1] - e[0][0],
        }
    return out


def _detect_cols(header: Sequence[str], first: Sequence[str]) -> Tuple[int, int]:
    low = [h.strip().lower() for h in header]
    id_i = pb_i = None
    for i, h in enumerate(low):
        if h in ("id", "read_id", "read", "name", "qname"):
            id_i = i
        if h in ("pbid", "isoform", "transcript_id", "pb_id", "model_id"):
            pb_i = i
    if pb_i is None:
        for i, v in enumerate(first):
            if str(v).startswith("PB.") and str(v).count(".") >= 2:
                pb_i = i
                break
    if pb_i is None:
        pb_i = len(header) - 1
    if id_i is None:
        id_i = 0 if pb_i != 0 else 1
    return id_i, pb_i


def read_assignments(
    path: str, wanted: Optional[set] = None, cap_per_model: int = 300
) -> Tuple[Dict[str, List[str]], Counter, int, int, int]:
    out: Dict[str, List[str]] = defaultdict(list)
    counts: Counter = Counter()
    seen: Dict[str, str] = {}
    dup_multi = 0
    n_rows = 0
    with _open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        first = next(rdr, None)
        id_i, pb_i = _detect_cols(header, first or header)

        def _add(row):
            nonlocal dup_multi, n_rows
            if len(row) <= max(id_i, pb_i):
                return
            n_rows += 1
            rid, pb = row[id_i], row[pb_i]
            counts[pb] += 1
            prev = seen.get(rid)
            if prev is None:
                seen[rid] = pb
            elif prev != pb:
                dup_multi += 1
            if (wanted is None or pb in wanted) and len(out[pb]) < cap_per_model:
                out[pb].append(rid)

        if first is not None:
            _add(first)
        for row in rdr:
            _add(row)
    return out, counts, n_rows, len(seen), dup_multi


def cigar_scan(
    bam_path: str,
    contigs: Optional[Sequence[str]] = None,
    max_reads: Optional[int] = 2_000_000,
    min_intron_len: int = 20,
) -> Tuple[Counter, int]:
    """Exons per aligned read straight from the CIGAR -- the step-0 invariant."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    if contigs:
        it: Iterable = (
            r for c in contigs if resolve_contig(bam, c) for r in bam.fetch(resolve_contig(bam, c))
        )
    else:
        it = bam.fetch(until_eof=True)
    hist: Counter = Counter()
    n = 0
    for r in it:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        n_ex = 1 + len(introns_of(r, min_intron_len, 10 ** 9))
        hist[min(n_ex, 20)] += 1
        n += 1
        if max_reads and n >= max_reads:
            break
    bam.close()
    return hist, n


def inspect(
    bam_path: str,
    gtf_path: str,
    read_stat_path: str,
    n_per_class: int = 25,
    max_reads_per_model: int = 300,
    window_pad: int = 5000,
    contigs: Optional[Sequence[str]] = None,
    cigar_scan_reads: Optional[int] = 2_000_000,
    verbose: bool = True,
) -> Dict[str, object]:
    models = parse_models(gtf_path)
    if verbose:
        # stderr, so that the JSON summary on stdout stays machine-readable
        print(f"[check] {len(models):,} models in {gtf_path}", file=sys.stderr)

    hist, n_scanned = cigar_scan(bam_path, contigs, cigar_scan_reads)
    mono_reads = hist.get(1, 0)
    frac_unspliced = mono_reads / max(n_scanned, 1)

    _all, counts, n_rows, n_uniq, dup_multi = read_assignments(read_stat_path, None, 1)
    tot_assigned = sum(counts.values())
    mono_mass = sum(n for pb, n in counts.items() if models.get(pb, {}).get("n_exons") == 1)
    frac_mono_mass = mono_mass / max(tot_assigned, 1)

    # pick the top models per exon-count class
    by_class: Dict[str, List[str]] = defaultdict(list)
    for pb, n in counts.most_common():
        m = models.get(pb)
        if m is None:
            continue
        cls = "mono-exon" if m["n_exons"] == 1 else "multi-exon"
        if len(by_class[cls]) < n_per_class:
            by_class[cls].append(pb)
    wanted = {pb for v in by_class.values() for pb in v}
    assign, _c, _r, _u, _d = read_assignments(read_stat_path, wanted, max_reads_per_model)

    bam = pysam.AlignmentFile(bam_path, "rb")
    rows = []
    for cls, pbs in by_class.items():
        for pb in pbs:
            m = models[pb]
            want = set(assign.get(pb, ()))
            if not want:
                continue
            contig = resolve_contig(bam, m["contig"])
            if contig is None:
                continue
            recs = []
            seen = set()
            for r in bam.fetch(contig, max(0, m["start"] - window_pad), m["end"] + window_pad):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if r.query_name not in want or r.query_name in seen:
                    continue
                seen.add(r.query_name)
                b = exon_blocks(r)
                if not b:
                    continue
                rs, re_ = b[0][0], b[-1][1]
                if m["strand"] == "+":
                    d5, d3 = rs - m["start"], re_ - m["end"]
                else:
                    d5, d3 = m["end"] - re_, m["start"] - rs
                recs.append((len(b), re_ - rs, d5, d3))
            if not recs:
                continue
            arr = np.array(recs, float)
            rows.append(
                {
                    "model": pb,
                    "class": cls,
                    "model_exons": m["n_exons"],
                    "model_span": m["span"],
                    "n_reads_checked": len(recs),
                    "fl_read_stat": counts.get(pb, 0),
                    "frac_reads_monoexon": float(np.mean(arr[:, 0] == 1)),
                    "median_read_exons": float(np.median(arr[:, 0])),
                    "median_log2_span_ratio": float(
                        np.median(np.log2((arr[:, 1] + 1) / (m["span"] + 1)))
                    ),
                    "median_d5": float(np.median(arr[:, 2])),
                    "median_d3": float(np.median(arr[:, 3])),
                    "frac_d3_within_100": float(np.mean(np.abs(arr[:, 3]) <= 100)),
                    "model_is_min": bool(m["n_exons"] <= arr[:, 0].min()),
                    "model_is_max": bool(m["n_exons"] >= arr[:, 0].max()),
                }
            )
    bam.close()

    df = pd.DataFrame(rows)
    summary = {
        "aligned_reads_scanned": n_scanned,
        "frac_reads_unspliced": round(frac_unspliced, 5),
        "read_stat_rows": n_rows,
        "read_stat_unique_ids": n_uniq,
        "read_ids_in_multiple_models": dup_multi,
        "frac_read_mass_on_monoexon_models": round(frac_mono_mass, 5),
        "step0_mismatch": bool(frac_mono_mass > max(2 * frac_unspliced, frac_unspliced + 0.02)),
    }
    if not df.empty:
        g = df.groupby("class").agg(
            n_models=("model", "size"),
            median_model_exons=("model_exons", "median"),
            median_read_exons=("median_read_exons", "median"),
            mean_frac_reads_monoexon=("frac_reads_monoexon", "mean"),
            median_log2_span_ratio=("median_log2_span_ratio", "median"),
            median_d5=("median_d5", "median"),
            median_d3=("median_d3", "median"),
            frac_model_is_min=("model_is_min", "mean"),
            frac_model_is_max=("model_is_max", "mean"),
        )
        summary["by_class"] = g.round(4).to_dict(orient="index")
    return {"summary": summary, "per_model": df}
