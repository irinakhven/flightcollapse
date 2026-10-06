"""Command line interface."""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
from typing import List, Optional, Sequence

from .config import Config


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("-c", "--config", help="JSON (or YAML) config file")
    p.add_argument(
        "--set", dest="overrides", action="append", default=[], metavar="KEY=VALUE",
        help="override any config field, e.g. --set ends.max_3p_diff=50 "
             "(repeatable, dotted paths)",
    )


def _build_config(args) -> Config:
    cfg = Config.load(args.config) if getattr(args, "config", None) else Config()
    for k in ("bam", "reference_gtf", "genome_fasta", "polya_site_bed",
              "cage_peak_bed", "cage_reptss_bed"):
        v = getattr(args, k, None)
        if v:
            setattr(cfg, k, v)
    if getattr(args, "short_read_sj", None):
        cfg.short_read_sj = list(args.short_read_sj)
    if getattr(args, "short_read_coverage", None):
        cfg.short_read_coverage = list(args.short_read_coverage)
    if getattr(args, "barcode_umi", None):
        cfg.molecules.barcode_umi_tsv = args.barcode_umi
    if getattr(args, "read_name_normalise", None):
        cfg.molecules.read_name_normalise = args.read_name_normalise
    if getattr(args, "outdir", None):
        cfg.output.outdir = args.outdir
    if getattr(args, "prefix", None):
        cfg.output.prefix = args.prefix
    if getattr(args, "contigs", None):
        cfg.contigs = list(args.contigs)
    if getattr(args, "workers", None) is not None:
        cfg.parallel.workers = int(args.workers)
    if getattr(args, "no_second_pass", False):
        cfg.secondpass.enabled = False
    if getattr(args, "no_posthoc", False):
        cfg.posthoc.enabled = False
    if getattr(args, "end3_tolerance", None) is not None:
        cfg.posthoc.end3_tolerance = int(args.end3_tolerance)
    if getattr(args, "end5_window", None) is not None:
        cfg.posthoc.end5_window = int(args.end5_window)
    if getattr(args, "no_subchain_filter", False):
        cfg.posthoc.subchain_filter = False
    if getattr(args, "subchain_ratio", None) is not None:
        cfg.posthoc.subchain_ratio = float(args.subchain_ratio)
    if getattr(args, "quiet", False):
        cfg.verbose = False
    if getattr(args, "no_strict", False):
        cfg.strict_invariants = False
    if not cfg.molecules.barcode_umi_tsv:
        cfg.molecules.write_matrix = False
    if getattr(args, "sample_type", None):
        cfg.sample_type = args.sample_type
    if getattr(args, "cmd", None) == "run":
        info = cfg.apply_sample_profile()
        if info["missing_inputs"] and not getattr(args, "quiet", False):
            print(f"[profile] {info['sample_type']}: {info['description']}")
            print("[profile] declared but not supplied: "
                  + ", ".join(info["missing_inputs"])
                  + " -- those channels are inactive for this run "
                    "(see evidence_channels in the QC report)")
    for ov in getattr(args, "overrides", []):
        if "=" not in ov:
            raise SystemExit(f"--set expects KEY=VALUE, got {ov!r}")
        k, v = ov.split("=", 1)
        cfg.set_path(k.strip(), v.strip())
    return cfg


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="flightcollapse",
        description="Reference-anchored collapse of aligned long reads into "
                    "transcript models, for oligo-dT single-cell Iso-Seq data.",
    )
    from . import __version__

    ap.add_argument("--version", action="version", version=f"flightcollapse {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # -- run ------------------------------------------------------------- #
    r = sub.add_parser("run", help="collapse a BAM into transcript models")
    r.add_argument("-b", "--bam")
    r.add_argument("-g", "--reference-gtf", dest="reference_gtf")
    r.add_argument("-f", "--genome-fasta", dest="genome_fasta")
    r.add_argument("--polya-bed", dest="polya_site_bed",
                   help="PolyASite / PolyA_DB BED of catalogued cleavage sites")
    r.add_argument("--cage-bed", dest="cage_peak_bed",
                   help="CAGE / TSS peak atlas (FANTOM5, refTSS) as BED6 or BED9. "
                        "Supports an alternative 5' end; its absence never drops "
                        "a model")
    r.add_argument("--cage-reptss", dest="cage_reptss_bed",
                   help="representative TSS as a separate BED6, joined on the "
                        "peak name. Needed when the peak file lost its BED9 thick "
                        "fields, e.g. after lifting FANTOM5 mouse to GRCm39")
    r.add_argument("--sample-type", dest="sample_type",
                   help="which SAMPLE_PROFILES preset to apply: retinal_organoid "
                        "(matched short reads), generic (atlases only, the "
                        "default), minimal (no external evidence). A declared "
                        "input that is missing warns and is reported as an "
                        "inactive channel; it never fails the run")
    r.add_argument("--short-read-coverage", dest="short_read_coverage", nargs="+",
                   help="short-read coverage bigWig(s): the 3' cleavage step, "
                        "the TSS ratio and long-transcript continuity. Needs "
                        "pyBigWig; without it the channel is simply inactive")
    r.add_argument("--short-read-sj", dest="short_read_sj", nargs="+",
                   help="STAR SJ.out.tab file(s); several are summed. Used as an "
                        "independent positive anchor for the calibration, never "
                        "as a filter")
    r.add_argument("-u", "--barcode-umi", help="read-name-keyed barcode/UMI TSV")
    r.add_argument("--umi-name-normalise", dest="read_name_normalise",
                   choices=["none", "strip_segment", "strip_ccs", "zmw"],
                   help="reconcile read names between the BAM and that table")
    r.add_argument("-o", "--outdir")
    r.add_argument("-p", "--prefix")
    r.add_argument("--contigs", nargs="+")
    r.add_argument("-j", "--workers", type=int, metavar="N",
                   help="run N contigs in parallel (0 = all cores). Output is "
                        "byte-identical to a serial run; the count is clamped "
                        "by available memory, since each worker holds one "
                        "contig's reads plus its own chromosome sequence")
    r.add_argument("--no-second-pass", action="store_true",
                   help="emit the first-pass catalogue unfiltered: no mono-exon "
                        "A-E classification, no end3_novel gate. The second-pass "
                        "evidence columns are still measured and written")
    r.add_argument("--end3-tolerance", dest="end3_tolerance", type=int, metavar="BP",
                   help="post-hoc 3' tolerance (bp): models with an identical "
                        "intron chain, 5' ends within posthoc.end5_window and "
                        "3' ends within BP are merged onto the annotated / most "
                        "supported end (default 500, the largest tolerance "
                        "that cost no SIRV recall; 0 = identical ends only, "
                        "-1 disables the 3' rule)")
    r.add_argument("--end5-window", dest="end5_window", type=int, metavar="BP",
                   help="guard on the 3' rule: a merge is blocked when two "
                        "same-chain models' 5' ends differ by more than BP "
                        "(default 100). Negative removes the guard, so the 3' "
                        "end is the whole rule and 5' truncation is left to the "
                        "sub-chain filter")
    r.add_argument("--no-subchain-filter", action="store_true",
                   help="keep models whose intron chain is a contiguous "
                        "sub-chain of a better-supported model (5'-truncation "
                        "fragments); by default they are folded into it")
    r.add_argument("--subchain-ratio", dest="subchain_ratio", type=float, metavar="X",
                   help="container must carry >= X times the fragment's support "
                        "for the fragment to be removed (default 2)")
    r.add_argument("--no-posthoc", action="store_true",
                   help="skip the final consolidation entirely (0.5.x behaviour)")
    r.add_argument("--no-strict", action="store_true",
                   help="report failed invariants instead of exiting non-zero")
    r.add_argument("-q", "--quiet", action="store_true")
    _common(r)

    # -- index-umi ------------------------------------------------------- #
    i = sub.add_parser("index-umi", help="pre-build the compact barcode/UMI index")
    i.add_argument("tsv")
    i.add_argument("-o", "--out",
                   help="index path (default: <tsv>.fcidx.<normalise>.npz, which is "
                        "exactly where `run` looks for it)")
    i.add_argument("--col-cell", default="cell_barcode")
    i.add_argument("--col-umi", default="umi")
    i.add_argument("--col-read", default="read_name")
    i.add_argument("--normalise", default="none",
                   choices=["none", "strip_segment", "strip_ccs", "zmw"],
                   help="reconcile read names between BAM and table")

    # -- diagnose-umi ---------------------------------------------------- #
    d = sub.add_parser(
        "diagnose-umi",
        help="work out which column and name form links a barcode/UMI table to a BAM",
    )
    d.add_argument("-b", "--bam", required=True)
    d.add_argument("-u", "--barcode-umi", required=True)
    d.add_argument("--rows", type=int, default=500_000,
                   help="table rows to sample (default 500k)")
    d.add_argument("--bam-reads", type=int, default=None,
                   help="cap BAM reads read into memory (default: all)")

    # -- check ----------------------------------------------------------- #
    k = sub.add_parser(
        "check",
        help="validate ANY collapse output (isoseq, bambu, stringtie, this tool)",
    )
    k.add_argument("-b", "--bam", required=True)
    k.add_argument("-g", "--gtf", required=True)
    k.add_argument("-r", "--read-stat", required=True)
    k.add_argument("--contigs", nargs="+")
    k.add_argument("--n-per-class", type=int, default=25)
    k.add_argument("--cigar-scan-reads", type=int, default=2_000_000)
    k.add_argument("-o", "--out", help="write the per-model table here")

    # -- compare --------------------------------------------------------- #
    m = sub.add_parser(
        "compare",
        help="cross-sample model concordance; works on any tool's GTF/GFF")
    m.add_argument("gtf", nargs="+", metavar="NAME=FILE.gtf",
                   help="one per sample, e.g. kin=out/kin/BD144_kin.gtf")
    m.add_argument("--fuzzy", type=int, default=5,
                   help="junction tolerance when matching intron chains")
    m.add_argument("--max-3p-diff", type=int, default=100)
    m.add_argument("--categories", nargs="+", default=[], metavar="NAME=FILE",
                   help="pigeon/SQANTI classification (or any id/category TSV) "
                        "supplying that sample's own categories and FL counts; "
                        "a collapse GFF carries neither")
    m.add_argument("-o", "--out", help="prefix for the concordance tables")

    # -- sj -------------------------------------------------------------- #
    s = sub.add_parser(
        "sj",
        help="check STAR short-read junctions and score a finished run against them",
    )
    s.add_argument("sj", nargs="+", metavar="SJ.out.tab",
                   help="one or more STAR SJ.out.tab files; replicates are summed")
    s.add_argument("-g", "--reference-gtf", dest="reference_gtf",
                   help="verify the coordinate convention against this annotation")
    s.add_argument("-j", "--junctions", dest="junctions_tsv",
                   help="a previous run's {prefix}.junctions.tsv.gz to evaluate")
    s.add_argument("--min-sr-unique", type=int, default=2)
    s.add_argument("-o", "--out", help="write the evaluation tables here (prefix)")

    # -- config ---------------------------------------------------------- #
    c = sub.add_parser("config", help="print a fully-populated default config")
    c.add_argument("-o", "--out")

    args = ap.parse_args(argv)

    if args.cmd == "run":
        from .pipeline import run as run_pipeline

        cfg = _build_config(args)
        try:
            report = run_pipeline(cfg)
        except RuntimeError as exc:
            print(f"\n{exc}", file=sys.stderr)
            return 2
        if report.get("all_passed", True):
            return 0
        # --no-strict promises to *report* rather than exit non-zero. It used to
        # suppress the exception and then return 1 anyway, which read as "the
        # run failed" for what is often a flag on 0.2% of models with every
        # output written correctly.
        failed = [r["invariant"] for r in report.get("invariants", []) if not r["passed"]]
        if report.get("failed_contigs"):
            failed.append(f"failed_contigs={list(report['failed_contigs'])}")
        if not cfg.strict_invariants:
            print(
                "\n[run] completed with " + str(len(failed)) + " reported problem(s): "
                + ", ".join(failed)
                + "\n[run] all outputs were written; see the qc report for detail.",
                file=sys.stderr,
            )
            return 0
        return 1

    if args.cmd == "index-umi":
        from .molecules import MoleculeIndex, default_index_path

        idx = MoleculeIndex.build(
            args.tsv, col_cell=args.col_cell, col_umi=args.col_umi,
            col_read=args.col_read, normalise=args.normalise,
        )
        # must match what MoleculeIndex.open_or_build looks for, otherwise the
        # pre-built index is silently ignored and `run` rebuilds it
        out = args.out or default_index_path(args.tsv, args.normalise)
        idx.save(out)
        print(f"wrote {out} ({len(idx):,} reads, {len(idx.barcodes):,} barcodes)")
        return 0

    if args.cmd == "diagnose-umi":
        from .diagnose import diagnose

        res = diagnose(args.bam, args.barcode_umi, table_rows=args.rows,
                       bam_reads=args.bam_reads)
        print()
        for r in res.get("results", []):
            flag = "  UNSAFE (collapses molecules)" if r.get("unsafe") else ""
            print(f"  col {r['column_index']:>2} {r['column']:<16} "
                  f"normalise={r['normalise']:<14} hits={r['hits']:>8,}{flag}")
        print()
        print(res.get("verdict") or res.get("error", ""))
        return 0 if (res.get("best") or {}).get("hits", 0) else 1

    if args.cmd == "check":
        from .check import inspect

        res = inspect(
            args.bam, args.gtf, args.read_stat,
            n_per_class=args.n_per_class,
            contigs=args.contigs,
            cigar_scan_reads=args.cigar_scan_reads,
        )
        print(json.dumps(res["summary"], indent=2, default=str))
        if args.out:
            res["per_model"].to_csv(args.out, sep="\t", index=False)
            print(f"\nper-model table -> {args.out}")
        if res["summary"].get("step0_mismatch"):
            print(
                "\nSTEP-0 MISMATCH: far fewer reads are unspliced than are assigned "
                "to mono-exon models.\nThe mono-exon class is being created "
                "downstream of alignment, i.e. by the collapse.",
                file=sys.stderr,
            )
            return 1
        return 0

    if args.cmd == "compare":
        from .compare import concordance, reproducibility_curve

        paths = {}
        for spec in args.gtf:
            if "=" not in spec:
                raise SystemExit(f"expected NAME=FILE, got {spec!r}")
            k, v = spec.split("=", 1)
            paths[k] = v
        cats = {}
        for spec in args.categories:
            if "=" not in spec:
                raise SystemExit(f"expected NAME=FILE, got {spec!r}")
            k, v = spec.split("=", 1)
            cats[k] = v
        C, summary = concordance(paths, fuzzy=args.fuzzy,
                                 max_3p_diff=args.max_3p_diff,
                                 categories=cats or None)
        print("\nSTRUCTURES, by how many samples call them")
        print(summary.to_string(index=False))
        curve = reproducibility_curve(C, len(paths))
        print("\nREPRODUCIBILITY vs SUPPORT THRESHOLD")
        print(curve.to_string(index=False))
        if args.out:
            C.to_csv(f"{args.out}.concordance.tsv.gz", sep="\t", index=False)
            summary.to_csv(f"{args.out}.summary.tsv", sep="\t", index=False)
            curve.to_csv(f"{args.out}.curve.tsv", sep="\t", index=False)
            print(f"\ntables -> {args.out}.{{concordance.tsv.gz,summary.tsv,curve.tsv}}")
        return 0

    if args.cmd == "sj":
        from .shortreads import SpliceJunctionCatalogue, evaluate_run

        cat = SpliceJunctionCatalogue.from_star(args.sj)
        n_novel = sum(1 for v in cat.j.values() if not v[4])
        print(f"\n{len(cat):,} distinct junctions, {n_novel:,} not in STAR's sjdb")

        conc = cat.concordance(min_unique=args.min_sr_unique)
        if conc:
            print("\nREPLICATE CONCORDANCE  (novel junctions with "
                  f">={args.min_sr_unique} summed unique reads)")
            for k in sorted(conc["by_n_replicates"], reverse=True):
                v = conc["by_n_replicates"][k]
                if v:
                    print(f"  seen in {k}/{conc['n_libraries']} replicates : {v:>7,}")
            print(f"  single-replicate expected by sampling : "
                  f"{conc['expected_by_sampling']:>7,.0f}")
            print(f"  excess sampling cannot explain        : "
                  f"{conc['excess_over_sampling']:>7,.0f} "
                  f"({conc['pct_of_set_unexplained']}% of the set)")
            print("  -> junctions.min_sr_replicates=2 removes these from the anchor set.")

        if args.reference_gtf:
            from .reference import ReferenceIndex

            ref = ReferenceIndex.from_gtf(args.reference_gtf, verbose=True)
            v = cat.verify_against(ref)
            print("\nCOORDINATE CHECK")
            print(f"  STAR-annotated junctions      : {v['sjdb_junctions']:,}")
            print(f"  ...also annotated here        : {v['sjdb_found_in_reference']:,} "
                  f"({v['agreement']:.1%})")
            print(f"  STAR-novel junctions          : {v['novel_junctions']:,}")
            print(f"  ...but annotated here         : {v['novel_also_in_reference']:,}")
            if not v["ok"]:
                print(f"\n  FAILED: {v.get('diagnosis', '')}", file=sys.stderr)
                return 1
            print("  coordinates line up.")

        if args.junctions_tsv:
            summary, bands = evaluate_run(
                args.junctions_tsv, cat, min_sr_unique=args.min_sr_unique
            )
            print("\nSHORT-READ SUPPORT BY CLASS")
            print(summary.to_string(index=False))
            if len(bands):
                print("\nNOVEL JUNCTIONS BY LONG-READ SUPPORT")
                print(bands.to_string(index=False))
            if args.out:
                summary.to_csv(f"{args.out}.sr_summary.tsv", sep="\t", index=False)
                bands.to_csv(f"{args.out}.sr_bands.tsv", sep="\t", index=False)
                print(f"\ntables -> {args.out}.sr_{{summary,bands}}.tsv")
        return 0

    if args.cmd == "config":
        cfg = Config()
        text = json.dumps(cfg.to_dict(), indent=2, default=list)
        if args.out:
            with open(args.out, "w") as fh:
                fh.write(text + "\n")
            print(f"wrote {args.out}")
        else:
            # tolerate `flightcollapse config | head`
            with contextlib.suppress(BrokenPipeError):
                print(text)
                sys.stdout.flush()
        return 0

    return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
