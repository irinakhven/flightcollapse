"""Orchestration.

Stages, in order, per contig (whole-genome collapse OOM-killed the VM, so
contigs are streamed and the read arrays are released before the next one):

1. **scan**        BAM -> compact per-read arrays, alignment gates applied once
2. **junctions**   census -> annotation-first snapping -> features -> calibrated
                   local FDR -> curation.  Reads carrying a rejected junction are
                   set aside, not silently promoted into novel isoforms.
3. **chains**      exact grouping on the curated chain; census; gene assignment
4. **novelty**     chain-level features -> calibrated local FDR -> candidate set
5. **suffixes**    directed resolution; the representative is always the maximal
                   chain
6. **ends**        3'-anchored peak clustering; 5' end from a high percentile;
                   naming against the reference, with UTR variants
7. **mono-exon**   a separate track that never touches the suffix logic
8. **second pass** (0.5.0) the catalogue is filtered, never rebuilt: mono-exon
                   models are judged against their gene's exon structure and
                   ``end3_novel`` models against per-molecule tail evidence
9. **outputs**     GTF/GFF, read_stat, group, abundance, molecule matrix
10. **validation** the invariants, checked before anything is declared done

Parallelism (0.5.0)
-------------------
Contigs are independent -- including the three calibrated models, which are
fitted per contig -- so ``parallel.workers`` distributes them over processes.
Everything genuinely global (model ids, the read-assignment files, the molecule
matrix) is assembled in the parent in the configured contig order, so a
parallel run is byte-identical to a serial one.  See :mod:`flightcollapse.parallel`.
"""

from __future__ import annotations

import os
import time
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam

from . import chains as chainmod
from . import parallel as pmod
from . import secondpass as spmod
from . import posthoc as phmod
from . import substitute as submod
from . import junctions as jmod
from . import output as omod
from .config import Config
from .coverage import CoverageTrack
from .ends import annotate_similarity, build_models_for_group, five_prime_end
from .genome import Genome, PolyASiteAtlas, TSSAtlas
from .model import TranscriptModel, polya_evidence
from .molecules import NO_CELL, MoleculeIndex, collapse_umis
from .monoexon import collapse_monoexonic
from .parallel import ContigResult
from .reads import ContigReads, gate_diagnostics, iter_alignments, resolve_contig
from .reference import ReferenceIndex
from .scoring import score_chains, score_junctions, score_terminals
from .validate import Validator, render_markdown


class GateRejectionError(RuntimeError):
    """The alignment gates discarded most of a contig's reads."""


def _log(cfg: Config, *a) -> None:
    if cfg.verbose:
        print(*a, flush=True)


def _check_name_match(cfg: Config, mol_index: MoleculeIndex, probe: int = 2000) -> None:
    """Fail fast if BAM read names do not match the barcode/UMI table.

    A silent mismatch here is invisible and catastrophic: every molecule count
    would be zero and every threshold would quietly fall back to raw reads.
    """
    bam = pysam.AlignmentFile(cfg.bam, "rb")
    names = []
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        names.append(r.query_name)
        if len(names) >= probe:
            break
    bam.close()
    if not names:
        return
    cid, _umi = mol_index.lookup(names)
    hit = float(np.mean(cid != NO_CELL))
    _log(cfg, f"[umi] {hit:.1%} of the first {len(names):,} BAM reads matched the table")
    if hit < 0.5:
        raise SystemExit(
            f"only {hit:.1%} of BAM read names were found in "
            f"{cfg.molecules.barcode_umi_tsv}.\n"
            f"  BAM example   : {names[0]}\n"
            f"  Try molecules.read_name_normalise = strip_segment | strip_ccs | zmw,\n"
            f"  or check molecules.col_read_name (currently "
            f"{cfg.molecules.col_read_name!r})."
        )


def select_contigs(cfg: Config, bam) -> List[str]:
    if cfg.contigs:
        return list(cfg.contigs)
    excl = set(cfg.exclude_contigs)
    out = []
    for c in bam.references:
        base = c if c.startswith("chr") else "chr" + c
        if c in excl or base in excl:
            continue
        out.append(c)
    return out


# ---------------------------------------------------------------------- #
class _Shared:
    """Read-only context every contig needs, plus the handles that must not fork.

    ``ReferenceIndex``, the atlases and the short-read junction catalogue are
    plain Python and numpy, so ``fork`` gives every worker a copy-on-write view
    of them for free -- which is the entire reason the pool uses ``fork``.

    The pysam and pyBigWig handles are different: a file descriptor and its
    buffered state shared across a fork produce interleaved reads and silently
    wrong data, not an error.  They are therefore opened lazily and re-opened
    the moment the process id changes, so a forked child never touches a handle
    it inherited.
    """

    def __init__(
        self,
        cfg: Config,
        ref: ReferenceIndex,
        atlas: Optional[PolyASiteAtlas],
        cage: Optional[TSSAtlas],
        sj,
        mol_index: Optional[MoleculeIndex],
        shard_dir: str,
    ) -> None:
        self.cfg = cfg
        self.ref = ref
        self.atlas = atlas
        self.cage = cage
        self.sj = sj
        self.mol_index = mol_index
        self.shard_dir = shard_dir
        self._pid: Optional[int] = None
        self._bam = None
        self._genome: Optional[Genome] = None
        self._cov: Optional[CoverageTrack] = None

    def handles(self):
        pid = os.getpid()
        if self._pid != pid:
            self._pid = pid
            cfg = self.cfg
            self._bam = pysam.AlignmentFile(cfg.bam, "rb")
            self._genome = Genome(cfg.genome_fasta) if cfg.genome_fasta else None
            self._cov = None
            if cfg.short_read_coverage and cfg.coverage.enabled:
                self._cov = CoverageTrack.open(
                    cfg.short_read_coverage, cfg.coverage.unstranded, verbose=False
                )
        return self._bam, self._genome, self._cov

    def close(self) -> None:
        if self._bam is not None:
            self._bam.close()
        if self._genome is not None:
            self._genome.close()
        if self._cov is not None:
            self._cov.close()
        self._pid = None
        self._bam = self._genome = self._cov = None


#: Set in the parent before the pool is created so ``fork`` inherits it.  A
#: module global rather than a closure because ``Pool`` pickles the worker by
#: qualified name and would otherwise have to pickle the reference index once
#: per contig.
_SHARED: Optional[_Shared] = None


def _contig_worker(contig: str) -> ContigResult:
    """Pool entry point.  Never raises: a failed contig is a reported result.

    Losing a completed genome's work to one late failure is a worse outcome
    than an incomplete result that says which contig failed and why -- the same
    rule the serial path has always followed, moved across the process
    boundary.
    """
    assert _SHARED is not None
    t0 = time.time()
    try:
        res = _run_contig(_SHARED, contig)
    except Exception as exc:  # noqa: BLE001 - deliberately broad, see docstring
        res = ContigResult(contig=contig, ok=False,
                           error=f"{type(exc).__name__}: {exc}")
    res.seconds = round(time.time() - t0, 1)
    return res


# ---------------------------------------------------------------------- #
def run(cfg: Config) -> Dict[str, object]:
    cfg.validate()
    t0 = time.time()
    os.makedirs(cfg.output.outdir, exist_ok=True)

    mol_index: Optional[MoleculeIndex] = None
    if cfg.molecules.barcode_umi_tsv:
        mol_index = MoleculeIndex.open_or_build(
            cfg.molecules.barcode_umi_tsv,
            cfg.molecules.index_path,
            verbose=cfg.verbose,
            col_cell=cfg.molecules.col_cell_barcode,
            col_umi=cfg.molecules.col_umi,
            col_read=cfg.molecules.col_read_name,
            normalise=cfg.molecules.read_name_normalise,
        )
        _check_name_match(cfg, mol_index)

    ref = ReferenceIndex.from_gtf(cfg.reference_gtf, cfg.contigs, verbose=cfg.verbose)
    genome = Genome(cfg.genome_fasta) if cfg.genome_fasta else None
    atlas = PolyASiteAtlas.from_bed(cfg.polya_site_bed, cfg.verbose) if cfg.polya_site_bed else None

    cage = None
    cage_check: Dict[str, object] = {}
    if cfg.cage_peak_bed and cfg.tss.enabled:
        cage = TSSAtlas.from_bed(
            cfg.cage_peak_bed, cfg.cage_reptss_bed,
            min_score=cfg.tss.min_peak_score, verbose=cfg.verbose,
        )
        # Same trap as the short-read junctions: a half-open slip or an
        # assembly mismatch does not fail, it returns "no evidence" for every
        # model, which looks exactly like an atlas that supports nothing.
        cage_check = cage.coverage_of_annotated_tss(ref)
        _log(cfg, f"[cage] {cage_check['annotated_tss_in_peak']:,}/"
                  f"{cage_check['n_genes_tested']:,} annotated TSSs fall inside a "
                  f"CAGE peak ({cage_check['frac_in_peak']:.1%}); "
                  f"{cage_check['frac_within_100bp_of_reptss']:.1%} within 100 bp "
                  f"of a representative TSS")
        if not cage_check["ok"]:
            raise RuntimeError(
                "CAGE peak coordinates do not line up with the reference:\n"
                f"  {cage_check.get('diagnosis', '')}\n"
                "  Refusing to run: every 5' end would come back unsupported, "
                "which is indistinguishable from a real result."
            )

    cov = None
    if cfg.short_read_coverage and cfg.coverage.enabled:
        cov = CoverageTrack.open(
            cfg.short_read_coverage, cfg.coverage.unstranded, cfg.verbose
        )
        if cov is None:
            _log(cfg, "[coverage] no usable track; the coverage channel is "
                      "inactive for this run (see evidence_channels in the QC)")

    sj = None
    sj_check: Dict[str, object] = {}
    if cfg.short_read_sj:
        from .shortreads import SpliceJunctionCatalogue

        sj = SpliceJunctionCatalogue.from_star(
            cfg.short_read_sj, min_overhang=0, verbose=cfg.verbose
        )
        # An off-by-one or a chr-prefix mismatch does not fail loudly here, it
        # returns "no support" for everything -- which looks exactly like a
        # library that confirms nothing. Check rather than trust.
        sj_check = sj.verify_against(ref, cfg.contigs)
        _log(cfg, f"[sj] {sj_check['sjdb_found_in_reference']:,}/"
                  f"{sj_check['sjdb_junctions']:,} of STAR's annotated junctions are "
                  f"annotated here ({sj_check['agreement']:.1%}); "
                  f"{sj_check['novel_junctions']:,} short-read junctions are novel")
        if not sj_check["ok"]:
            raise RuntimeError(
                "short-read junction coordinates do not line up with the "
                f"reference:\n  {sj_check.get('diagnosis', '')}\n"
                "  Refusing to run: every junction would silently come back "
                "unsupported, which is indistinguishable from a real result."
            )

    bam = pysam.AlignmentFile(cfg.bam, "rb")
    contigs = [c for c in select_contigs(cfg, bam) if resolve_contig(bam, c) is not None]
    lengths = dict(zip(bam.references, bam.lengths))
    bam.close()
    _log(cfg, f"[run] {len(contigs)} contigs: {', '.join(contigs[:8])}"
              f"{' ...' if len(contigs) > 8 else ''}")

    # The parent's own copies of the fork-unsafe handles are no longer needed:
    # every consumer now goes through _Shared.handles(), which opens one set per
    # process. Closing them here rather than leaving them open is what makes
    # "inherited but never touched" true rather than merely intended.
    if genome is not None:
        genome.close()
        genome = None
    if cov is not None:
        cov.close()
        cov = None

    workers = cfg.worker_count()
    requested = cfg.parallel.workers if cfg.parallel.workers != 1 else cfg.threads
    # Say so whenever the clamp bit -- ESPECIALLY when it bit all the way down
    # to 1. A `-j 12` that silently runs serially is the version of this that
    # wastes an afternoon before anyone thinks to check the QC file.
    if requested not in (0, 1) and requested != workers:
        _log(cfg, f"[run] {requested} workers requested, running {workers}: "
                  f"clamped by available memory "
                  f"({cfg.parallel.gb_per_worker:g} GB budgeted per worker, "
                  f"{cfg.parallel.reserve_gb:g} GB reserved for the parent). "
                  f"Lower parallel.gb_per_worker if your contigs are smaller "
                  f"than that, or set parallel.clamp_by_memory=false to "
                  f"override.")
    shard_dir = cfg.parallel.shard_dir or os.path.join(
        cfg.output.outdir, f".{cfg.output.prefix}_shards"
    )
    os.makedirs(shard_dir, exist_ok=True)

    global _SHARED
    _SHARED = _Shared(cfg, ref, atlas, cage, sj, mol_index, shard_dir)

    _log(cfg, f"[run] {workers} worker{'s' if workers != 1 else ''} over "
              f"{len(contigs)} contigs")
    order = pmod.order_contigs(contigs, lengths, cfg.parallel.longest_first)
    results = pmod.map_contigs(
        _contig_worker, order, workers, cfg.parallel.start_method,
        log=(lambda s: _log(cfg, s)) if cfg.verbose else None,
    )
    _SHARED.close()

    # -- merge, in the CONFIGURED contig order, never the completion order ---
    validator = Validator()
    ids = omod.IdAssigner(cfg.output.id_style)
    all_models: List[TranscriptModel] = []
    cell_counts: Dict[str, Dict[int, int]] = {}
    group_chain_rows: List[Tuple[str, str, int, int, bool]] = []
    substitution_rows: List[Tuple] = []
    junction_frames: List[pd.DataFrame] = []
    model_reports: List[Dict[str, object]] = []
    mono_reject_rows: List[Dict[str, object]] = []
    secondpass_rows: List[Dict[str, object]] = []
    posthoc_rows: List[Dict[str, object]] = []
    per_contig: Dict[str, Dict[str, object]] = {}
    failed: Dict[str, str] = {}
    stitch_plan: List[Tuple[str, List[str]]] = []
    n_unassigned = 0

    for contig in contigs:
        res = results.get(contig)
        if res is None:
            continue
        if not res.ok:
            failed[contig] = res.error
            _log(cfg, f"[{contig}] FAILED, continuing with the remaining contigs\n"
                      f"        {res.error}")
            continue
        if not res.summary:
            continue

        # Ids are assigned here, contig by contig in order, which is exactly
        # what the serial path did -- so the names, the PB.N numbering and the
        # collision counter all come out identical.
        ids.assign(res.models)
        contig_ids = [m.model_id for m in res.models]
        base = len(all_models)
        all_models.extend(res.models)

        for anchor, ch, nr, nm, is_model in res.group_chain_rows:
            group_chain_rows.append((contig_ids[anchor], ch, nr, nm, is_model))
        substitution_rows.extend(res.substitution_rows)
        if res.junction_frame is not None and len(res.junction_frame):
            junction_frames.append(res.junction_frame)
        model_reports.extend(res.model_reports)
        mono_reject_rows.extend(res.mono_reject_rows)
        for row in res.posthoc_audit:
            mi = int(row.get("into_model_index", -1))
            row["into_model_id"] = contig_ids[mi] if 0 <= mi < len(contig_ids) else ""
            posthoc_rows.append(row)
        for row in res.secondpass_audit:
            mi = int(row.get("model_index", -1))
            row["model_id"] = contig_ids[mi] if 0 <= mi < len(contig_ids) else ""
            secondpass_rows.append(row)
        validator.merge(res.validator)
        per_contig[contig] = res.summary
        n_unassigned += int(res.n_unassigned)
        if res.shard_stem:
            stitch_plan.append((res.shard_stem, contig_ids))
            if cfg.molecules.write_matrix and mol_index is not None:
                cell_counts.update(pmod.load_counts(res.shard_stem, contig_ids))
        _ = base  # ids are contig-local in the shards; base is informational

    stitched = pmod.stitch_shards(
        stitch_plan,
        cfg.outpath(".read_stat.txt"),
        cfg.outpath(".group.txt"),
        gz=cfg.output.gzip_big_tables,
    )
    pmod.cleanup(shard_dir, cfg.parallel.keep_shards)

    if secondpass_rows and cfg.secondpass.write_audit:
        spmod.write_audit(secondpass_rows, cfg.outpath(".secondpass.tsv"))
    if cfg.posthoc.enabled and cfg.posthoc.write_audit:
        phmod.write_audit(posthoc_rows, cfg.outpath(".posthoc.tsv"))

    _write_all(cfg, all_models, cell_counts, group_chain_rows,
               substitution_rows, junction_frames,
               mol_index, ref)

    validator.extra["model_id_collisions"] = ids.collisions
    report = validator.report()
    report["runtime_seconds"] = round(time.time() - t0, 1)
    report["n_models"] = len(all_models)
    report["models_by_category"] = dict(
        pd.Series([m.category for m in all_models]).value_counts()
    ) if all_models else {}
    report["novelty_models"] = model_reports
    report["novel_model_burden"] = _novel_burden(all_models)
    report["three_prime"] = _three_prime_report(all_models, cfg)
    report["model_id_collisions"] = ids.collisions
    report["per_contig"] = per_contig
    report["evidence_channels"] = cfg.evidence_channels()
    report["parallel"] = {
        "workers": workers,
        "requested": requested,
        "start_method": cfg.parallel.start_method if workers > 1 else "serial",
        "contigs": len(contigs),
        "longest_first": bool(cfg.parallel.longest_first),
        "slowest_contig": max(
            ((c, r.seconds) for c, r in results.items()),
            key=lambda kv: kv[1], default=("", 0.0),
        )[0],
        "contig_seconds": {c: r.seconds for c, r in sorted(results.items())},
    }
    report["read_stat_rows"] = int(stitched.get("read_stat_rows", 0))
    report["reads_unassigned"] = int(n_unassigned)
    if cfg.posthoc.enabled:
        _ph = [r.summary.get("posthoc", {}) for r in results.values() if r.ok and r.summary]
        report["posthoc"] = {
            "end3_tolerance": int(cfg.posthoc.end3_tolerance),
            "end5_window": int(cfg.posthoc.end5_window),
            "subchain_filter": bool(cfg.posthoc.subchain_filter),
            "subchain_ratio": float(cfg.posthoc.subchain_ratio),
            "models_before": int(sum(d.get("models_before", 0) for d in _ph)),
            "models_after": int(sum(d.get("models_after", 0) for d in _ph)),
            "end3_models_absorbed": int(sum(d.get("end3_models_absorbed", 0) for d in _ph)),
            "subchain_removed": int(sum(d.get("subchain_removed", 0) for d in _ph)),
            "subchain_reads_folded": int(sum(d.get("subchain_reads_folded", 0) for d in _ph)),
            "audit": cfg.outpath(".posthoc.tsv") if cfg.posthoc.write_audit else "",
        }
    if cfg.secondpass.enabled:
        report["second_pass"] = _second_pass_report(
            per_contig, secondpass_rows, cfg
        )
    if cfg.short_read_coverage and cfg.coverage.enabled:
        report["short_read_coverage"] = {
            "sources": list(cfg.short_read_coverage),
            "unstranded": bool(cfg.coverage.unstranded),
        }
    if mono_reject_rows:
        import pandas as _pd
        _mr = _pd.DataFrame(mono_reject_rows)
        _mr.to_csv(cfg.outpath(".monoexon_rejected.tsv"), sep="\t", index=False)
        report["monoexon_rejections"] = {
            "n": int(len(_mr)),
            "by_reason": {str(k): int(v) for k, v in
                          _mr["reason"].value_counts().items()},
            "reads_by_reason": {str(k): int(v) for k, v in
                                _mr.groupby("reason").n_reads.sum().items()},
            "table": cfg.outpath(".monoexon_rejected.tsv"),
        }
    if cage_check:
        report["cage_atlas"] = {
            "sources": list(cage.sources) if cage else [],
            "n_peaks": int(len(cage)) if cage else 0,
            "rescue_suffix": bool(cfg.tss.rescue_suffix),
            **cage_check,
        }
    if sj_check:
        from .shortreads import holdout_summary

        report["short_read_junctions"] = {
            "sources": list(cfg.short_read_sj or []),
            **sj_check,
            "validation": holdout_summary(junction_frames),
        }
    if junction_frames:
        import pandas as _pd
        _jf = _pd.concat(junction_frames, ignore_index=True)
        _novel = _jf[~_jf["annotated"].astype(bool)]
        _rej = _novel[~_novel["keep"].astype(bool)]
        report["junction_curation"] = {
            "novel_junctions": int(len(_novel)),
            "novel_kept": int(_novel["keep"].astype(bool).sum()),
            "rejected_by_reason": {
                str(k): int(v) for k, v in _rej["drop_reason"].value_counts().items()
            },
            # a rule discarding junctions the calibrated model scores as real is
            # the signature this release was written to remove; make it visible
            # without needing an experiment to find it
            "median_lfdr_of_rejected": (
                float(_rej["lfdr"].median()) if "lfdr" in _rej and len(_rej) else None
            ),
            "median_reads_of_rejected": (
                float(_rej["n_reads"].median()) if len(_rej) else None
            ),
            "median_reads_of_kept": (
                float(_novel.loc[_novel["keep"].astype(bool), "n_reads"].median())
                if _novel["keep"].astype(bool).any() else None
            ),
        }
    # Only meaningful once there is enough data for the floor to mean anything:
    # a 500-read simulation legitimately has no decoys, and failing the run for
    # that would make the tripwire useless by making it noisy.
    _acc = report.get("read_accounting", {}) or {}
    _ndecoy = int(_acc.get("decoy_chain_groups", 0))
    _nspliced = int(_acc.get("spliced", 0))
    if _nspliced >= 1_000_000 and _ndecoy < cfg.junctions.min_decoy_chain_groups:
        _log(cfg, f"[run] WARNING: only {_ndecoy:,} decoy chain groups from "
                  f"{_nspliced:,} spliced reads (floor "
                  f"{cfg.junctions.min_decoy_chain_groups:,}). The chain-level "
                  "novelty model may have fallen back to flat thresholds -- "
                  "see junction_curation in the QC report.")
        report["all_passed"] = False
    report["failed_contigs"] = failed
    if failed:
        report["all_passed"] = False
    report["config"] = cfg.to_dict()

    omod.write_json(report, cfg.outpath(".qc.json"))
    with open(cfg.outpath(".qc.md"), "w") as fh:
        fh.write(render_markdown(report))
    _log(cfg, f"[run] {len(all_models):,} models in {report['runtime_seconds']}s "
              f"-> {cfg.output.outdir}")

    _SHARED = None
    if failed:
        _log(cfg, "\n[run] contigs that failed: "
                  + ", ".join(f"{k} ({v})" for k, v in failed.items())
                  + "\n[run] everything else was written normally.")
    if cfg.strict_invariants:
        validator.raise_on_failure()
    return report


# ---------------------------------------------------------------------- #
def _run_contig(shared: "_Shared", contig: str) -> ContigResult:
    """One contig, end to end, with no reference to any other contig's state.

    Returns everything the parent needs rather than mutating shared
    accumulators.  That is the whole refactor: the stage logic below is byte
    for byte what 0.4.6 did, but a function that only reads its inputs and
    returns its outputs is one that can run in another process.
    """
    cfg = shared.cfg
    ref = shared.ref
    atlas = shared.atlas
    cage = shared.cage
    sj = shared.sj
    mol_index = shared.mol_index
    bam, genome, cov = shared.handles()

    out = ContigResult(contig=contig)
    validator = Validator()
    model_reports: List[Dict[str, object]] = []
    mono_reject_rows: List[Dict[str, object]] = []
    substitution_rows: List[Tuple] = []

    t0 = time.time()
    gate_stats: Dict[str, int] = {}
    reads = ContigReads.scan(bam, contig, cfg.gates, mol_index, gate_stats=gate_stats)
    scanned = gate_stats.get("scanned", 0)
    if scanned:
        rejected = 1.0 - gate_stats.get("admitted", 0) / scanned
        _log(cfg, f"[{contig}] alignment gates: {gate_stats.get('admitted', 0):,}/"
                  f"{scanned:,} admitted ({rejected:.1%} rejected: "
                  + ", ".join(f"{k}={v:,}" for k, v in sorted(gate_stats.items())
                              if k not in ("scanned", "admitted")) + ")")
        if (rejected > cfg.gates.max_gate_rejection_frac
                and scanned >= cfg.gates.min_reads_for_gate_check):
            diag = gate_diagnostics(bam, contig, cfg.gates)
            raise GateRejectionError(
                f"[{contig}] the alignment gates rejected {rejected:.1%} of aligned "
                f"reads. That is a configuration error, not a filter.\n"
                f"  observed: {diag}\n"
                f"  If this is an untrimmed FLNC BAM the polyA tail is soft-clipped, "
                f"so coverage against the full query length is low by construction.\n"
                f"  Lower gates.min_aln_coverage (currently "
                f"{cfg.gates.min_aln_coverage}), run on the polyA-trimmed BAM, or "
                f"raise gates.max_gate_rejection_frac to accept this deliberately."
            )
        if rejected > cfg.gates.max_gate_rejection_frac:
            _log(cfg, f"[{contig}] high rejection on only {scanned:,} reads -- "
                      f"below gates.min_reads_for_gate_check, not treated as an error")
    if reads.n == 0:
        _log(cfg, f"[{contig}] no admitted reads")
        return out
    summ = reads.summary()
    _log(cfg, f"[{contig}] {summ['reads']:,} reads "
              f"({summ['mono_exonic_reads']:,} unspliced, "
              f"{summ['mono_exonic_reads'] / summ['reads']:.1%})")

    # -- 2. junction catalogue -------------------------------------------
    # An exon-internal deletion never touches the intron chain, so a short
    # intron the aligner spelled as `D` would hide inside an otherwise perfect
    # FSM.  Promote the ones that are annotated or recurrent-and-canonical
    # BEFORE the census, so they are curated like every other junction.
    gap_stats = jmod.promote_internal_gaps(reads, contig, ref, genome, cfg.junctions)
    if gap_stats.get("gap_distinct"):
        _log(cfg, f"[{contig}] internal deletions: {gap_stats['gap_distinct']:,} distinct, "
                  f"{gap_stats['gap_promoted_annotated']:,} annotated + "
                  f"{gap_stats['gap_promoted_motif']:,} canonical promoted to junctions "
                  f"({gap_stats['gap_records_promoted']:,} records)")
    unexplained_gap = jmod.reads_with_unexplained_gap(reads)

    census = jmod.JunctionCensus.from_reads(reads)
    dmap, amap = jmod.build_snap_maps(census, contig, ref, cfg.junctions)
    reads.apply_snap(dmap, amap)
    census = jmod.JunctionCensus.from_reads(reads)
    feat = jmod.junction_features(census, contig, ref, genome, cfg.junctions, sj)
    jscore, jmodel = score_junctions(feat, cfg.scoring)
    feat = pd.concat([feat, jscore], axis=1)
    cur = jmod.curate(feat, cfg.junctions, feat["lfdr"].to_numpy() if "lfdr" in feat else None)
    cur["contig"] = contig
    out.junction_frame = cur
    if jmodel.fitted or jmodel.note:
        rep = jmodel.report()
        rep["contig"] = contig
        model_reports.append(rep)
    _log(cfg, f"[{contig}] junctions: {len(cur):,} unique, "
              f"{int(cur['annotated'].sum()):,} annotated, "
              f"{int((~cur['annotated'] & cur['keep']).sum()):,} novel kept, "
              f"{int((~cur['keep']).sum()):,} rejected")

    keep_mask = cur["keep"].to_numpy()
    read_ok = jmod.read_junction_status(reads, census, keep_mask)

    jinfo: Dict[Tuple[str, int, int], Dict[str, float]] = {}
    for row in cur.itertuples(index=False):
        if not row.keep:
            continue
        jinfo[(row.strand, int(row.donor), int(row.acceptor))] = {
            "annotated": bool(row.annotated),
            # Rule J needs the raw count, and needs ABSENT and ZERO to stay
            # apart: a junction with no SJ-table entry is unmeasured, not
            # unsupported, and must never be substituted away.
            "sr_uniq": (int(row.sr_uniq) if hasattr(row, "sr_uniq")
                        and row.sr_uniq == row.sr_uniq else None),
            "posterior_real": float(getattr(row, "posterior_real", float("nan"))),
            "decoy": int(
                (not row.annotated)
                and (
                    row.motif_class in ("antisense", "non_canonical")
                    or row.direct_repeat >= 10
                )
            ),
        }

    # -- 3./4. chains -----------------------------------------------------
    groups, gid_of_read = chainmod.build_chain_groups(
        reads, contig, ref, read_ok, cfg.molecules.umi_hamming
    )
    chainmod.annotate_groups_with_junctions(groups, jinfo)
    gene_stats = chainmod.assign_genes(groups, ref)

    gene_support: Dict[str, int] = defaultdict(int)
    for g in groups:
        if g.gene_id:
            gene_support[g.gene_id] += g.support

    # Chains built from reads that carry a REJECTED junction are, by
    # construction, artefactual chains.  They never become models, but they are
    # exactly the negative anchor the chain-level novelty model needs -- without
    # them the decoy set is empty and the calibration silently falls back to
    # flat thresholds.
    decoy_groups: List["chainmod.ChainGroup"] = []
    if (~read_ok).any():
        decoy_groups, _ = chainmod.build_chain_groups(
            reads, contig, ref, ~read_ok, cfg.molecules.umi_hamming
        )
        chainmod.assign_genes(decoy_groups, ref)
        for dg in decoy_groups:
            dg.n_decoy_junctions = 1
            dg.min_junction_posterior = 0.0
            dg.frac_junctions_annotated = 0.0
            dg.tx_ids = []

    cfeat = _chain_features(groups, reads, ref, genome, atlas, cage, cfg, gene_support)
    if decoy_groups:
        dfeat = _chain_features(
            decoy_groups, reads, ref, genome, atlas, cage, cfg, gene_support
        )
        combined = pd.concat([cfeat, dfeat], ignore_index=True)
        cscore_all, cmodel = score_chains(combined, cfg.scoring)
        cscore = cscore_all.iloc[: len(cfeat)].reset_index(drop=True)
    else:
        cscore, cmodel = score_chains(cfeat, cfg.scoring)
    if cmodel.fitted or cmodel.note:
        rep = cmodel.report()
        rep["contig"] = contig
        model_reports.append(rep)
    lfdr_by_gid = {
        int(g): float(v)
        for g, v in zip(cfeat["gid"].to_numpy(), cscore["lfdr"].to_numpy())
    }
    post_by_gid = {
        int(g): float(v)
        for g, v in zip(cfeat["gid"].to_numpy(), cscore["posterior_real"].to_numpy())
    }

    candidate = chainmod.mark_candidates(groups, cfg.chains, gene_support, lfdr_by_gid)

    # -- Rule J (0.4.3): novel junctions the short reads say did not happen.
    # Before suffix resolution, so a substituted group is already pointed at its
    # annotated destination and cannot be re-homed somewhere else on the way.
    sub_stats: Dict[str, int] = {}
    if cfg.substitution.enabled:
        sr_counts = {
            k: int(v["sr_uniq"]) for k, v in jinfo.items()
            if v.get("sr_uniq") is not None
        }
        sub_plans, sub_stats = submod.plan_substitutions(
            groups, candidate, ref, sr_counts, cfg.substitution,
            cfg.junctions.fuzzy_tolerance,
        )
        submod.apply_substitutions(groups, sub_plans)
        if sub_plans and cfg.substitution.write_plan:
            substitution_rows.extend(
                (groups[g].contig, groups[g].strand, g, t, groups[g].n_reads, sp)
                for g, t, plan in sub_plans for sp in plan
            )
    merge_stats = chainmod.resolve_suffixes(groups, candidate, cfg.chains)
    merge_stats.update(chainmod.attach_subthreshold(groups, candidate, cfg.chains))
    validator.add_merge_stats(merge_stats)
    _log(cfg, f"[{contig}] chains: {len(groups):,} groups, "
              f"{int(candidate.sum()):,} candidates, "
              f"{merge_stats['n_merged']:,} merged "
              f"({merge_stats['n_dominant_protected']:,} dominant protected, "
              f"{merge_stats['n_ambiguous']:,} ambiguous)")

    # -- 5./6. models -----------------------------------------------------
    children: Dict[int, List[int]] = defaultdict(list)
    for g in groups:
        if g.merged_into is not None:
            children[g.merged_into].append(g.gid)

    contig_models: List[TranscriptModel] = []
    #: ``(anchor_model, chain_repr, n_reads, n_molecules, is_model_chain)``.
    #: 0.5.0 anchors on the MODEL OBJECT rather than its index: both the
    #: terminal-lFDR ceiling and the second pass can remove models, and an
    #: integer index into a list that is about to be filtered silently re-points
    #: a group's chain table at whichever model slid into that slot.
    pending_chain_rows: List[Tuple[TranscriptModel, str, int, int, bool]] = []

    for g in groups:
        if g.merged_into is not None or not candidate[g.gid]:
            continue
        if g.ambiguous and not cfg.chains.keep_ambiguous_as_models:
            continue
        member_gids = [g.gid] + children.get(g.gid, [])
        idx = np.concatenate([groups[k].reads for k in member_gids])
        # every transcript sharing this chain, not the canonical one among them
        ms = build_models_for_group(
            reads, contig, g.strand, g.chain, idx, ref, genome, atlas,
            cfg.ends, cfg.monoexon, g.tx_ids, g.gene_id, cfg.molecules.umi_hamming,
            cage=cage, tss_params=cfg.tss, terminal_params=cfg.terminal,
            cov=cov, cov_params=cfg.coverage,
        )
        n_super = len(ref.annotated_superchains(g.gene_id, g.chain))
        n_gap = int(np.sum(unexplained_gap[idx])) if unexplained_gap.size else 0
        for m in ms:
            m.lfdr = lfdr_by_gid.get(g.gid, float("nan"))
            m.posterior_real = post_by_gid.get(g.gid, float("nan"))
            if not g.annotated:
                m.category = _novel_category(g)
                if g.merge_class == "alt_tss_kept":
                    # A chain that is a proper suffix of an annotated one, kept
                    # because its start site has independent support. Calling it
                    # NIC is misleading -- it is not a novel COMBINATION of
                    # junctions, it is a known chain with a novel start, and the
                    # distinction is the whole point of keeping it.
                    m.category = "end5_novel"
            m.evidence["tss_evidence_source"] = g.tss_evidence_source
            # How much of this model's own read mass actually carries the
            # model's full intron chain. Members merged in from suffix groups do
            # not, by construction. Computed here because `gid_of_read` is the
            # only place the pre-merge chain identity still exists -- and it is
            # the statistic that separates a well-supported parent from one
            # inferred from a handful of reads.
            gids = gid_of_read[m.reads]
            n_exact = int(np.count_nonzero(gids == g.gid))
            m.evidence["exact_chain_reads"] = n_exact
            m.evidence["exact_chain_frac"] = (
                round(n_exact / len(m.reads), 6) if len(m.reads) else ""
            )
            m.evidence["n_distinct_chains"] = int(np.unique(gids).size)
            if g.ambiguous:
                m.flags.append("ambiguous_suffix_parent")
            if g.n_novel_junctions:
                m.flags.append(f"novel_junctions={g.n_novel_junctions}")
            m.evidence["n_annotated_superchains"] = n_super
            m.evidence["n_internal_indel_reads"] = n_gap
            if n_super:
                # FSM to a chain that is a 5'-truncated form of a longer
                # annotated transcript.  Deliberate, but it moves read mass off
                # the long isoform onto the fragment in proportion to how
                # degraded the library is, so it is never implicit.
                m.flags.append(f"suffix_of_annotated={n_super}")
            if g.gene_conflict:
                m.evidence["gene_conflict"] = g.gene_conflict
                m.flags.append("gene_assignment_conflict")
        if not ms:
            continue
        anchor = ms[0]
        contig_models.extend(ms)
        for k in member_gids:
            chg = groups[k]
            pending_chain_rows.append(
                (anchor, omod.chain_repr(chg.chain), chg.n_reads, chg.n_mols, k == g.gid)
            )

    # -- 7. mono-exonic track ---------------------------------------------
    # 0.3.0: the mono track now sees what the spliced track already emitted, so
    # a cluster sitting on an existing model's 3' end can be recognised as that
    # transcript's last exon rather than becoming a transcript of its own.
    spliced_3p: Dict[str, List[int]] = defaultdict(list)
    for m in contig_models:
        if m.gene_id and m.n_exons > 1:
            spliced_3p[m.gene_id].append(
                m.exons[-1][1] if m.strand == "+" else m.exons[0][0]
            )
    spliced_3p_arr = {k: np.array(sorted(v), np.int64) for k, v in spliced_3p.items()}
    # 0.4.0 Rule M needs SPANS, not just ends: a fragment is recognised by
    # lying wholly inside an emitted spliced model, not by stopping where it
    # stops. Columns are (start, end, three_prime) so the fold target can be
    # named by the same (gene, 3' end) key the demotion list already uses.
    spliced_span: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for m in contig_models:
        if m.gene_id and m.n_exons > 1:
            three = m.exons[-1][1] if m.strand == "+" else m.exons[0][0]
            spliced_span[m.gene_id].append(
                (int(m.exons[0][0]), int(m.exons[-1][1]), int(three))
            )
    spliced_span_arr = {k: np.array(v, np.int64) for k, v in spliced_span.items()}
    mono_rejections: List[Dict[str, object]] = []
    mono_demoted: List[Tuple[str, int, np.ndarray]] = []
    mono_models, mono_stats = collapse_monoexonic(
        reads, contig, ref, genome, atlas, cfg.monoexon, read_ok,
        cfg.ends.five_prime_percentile, cfg.molecules.umi_hamming,
        spliced_3p=spliced_3p_arr,
        rejections=mono_rejections if cfg.monoexon.write_rejection_reasons else None,
        end_params=cfg.ends,
        demoted=mono_demoted,
        spliced_spans=spliced_span_arr,
        terminal_params=cfg.terminal,
        cov=cov,
        cov_params=cfg.coverage,
    )
    if mono_rejections:
        mono_reject_rows.extend(mono_rejections)
    # A demoted fragment's reads belong to the spliced model whose 3' end it
    # sits on. Attaching them here keeps the read accounting exact: the mono
    # track gives up a model, not the evidence behind it.
    if mono_demoted:
        by_key: Dict[Tuple[str, int], TranscriptModel] = {}
        for m in contig_models:
            if m.gene_id and m.n_exons > 1:
                e = m.exons[-1][1] if m.strand == "+" else m.exons[0][0]
                by_key.setdefault((m.gene_id, int(e)), m)
        n_folded = 0
        for gene_id, end, grp in mono_demoted:
            tgt = by_key.get((gene_id, int(end)))
            if tgt is None:
                continue
            tgt.reads = np.concatenate([tgt.reads, grp])
            tgt.n_reads += int(grp.size)
            tgt.evidence["n_3p_fragment_reads"] = int(
                tgt.evidence.get("n_3p_fragment_reads", 0) or 0
            ) + int(grp.size)
            tgt.flags = sorted(set(tgt.flags) | {"absorbed_monoexon_fragment"})
            n_folded += int(grp.size)
        if n_folded:
            _log(cfg, f"[{contig}] mono-exon: {n_folded:,} fragment reads folded "
                      f"into the spliced models they belong to")

    contig_models.extend(mono_models)
    _log(cfg, f"[{contig}] mono-exon: {mono_stats['monoexon_reads']:,} reads -> "
              f"{len(mono_models):,} models "
              f"(3'UTR {mono_stats['models_utr3']}, "
              f"internal {mono_stats['models_internal']}, "
              f"intergenic {mono_stats['models_intergenic']}; "
              f"{mono_stats.get('rejected_demoted_fragment', 0):,} demoted as "
              f"fragments, {mono_stats.get('rescued_by_tail', 0):,} rescued on "
              f"tail evidence)")

    # -- 7b. nearest reference transcript ----------------------------------
    # Annotation only: this cannot change which models exist, only what each
    # row says about itself.  Its one classification effect is intron
    # retention, which the exact-chain lookup has no way to express.
    annotate_similarity(contig_models, ref)
    n_ir = sum(1 for m in contig_models if m.category == "IR")
    n_rt = sum(1 for m in contig_models if m.category == "READTHROUGH")
    n_named = sum(1 for m in contig_models if m.associated_tx)
    _log(cfg, f"[{contig}] associated a reference transcript with "
              f"{n_named:,}/{len(contig_models):,} models "
              f"({n_ir:,} intron retention, {n_rt:,} readthrough)")

    # -- 7b-bis. terminal evidence the second pass and the score both read ---
    # 0.5.0. Deliberately BEFORE score_terminals, so the calibrated model sees
    # the templated-tail and composition features too rather than the second
    # pass being their only consumer. It writes into m.evidence and decides
    # nothing -- `m.reads` still holds contig-local indices at this point,
    # which is what makes the per-read tail comparison possible at all.
    spmod.measure_terminal_evidence(
        contig_models, contig, genome, reads, cfg.secondpass,
        measure_tails=cfg.gates.measure_polya_tail,
    )

    # -- 7c. calibrated terminal score --------------------------------------
    if cfg.terminal.score_terminals:
        tfeat = _terminal_features(contig_models)
        if len(tfeat):
            tscore, tmodel = score_terminals(
                tfeat, cfg.scoring, cfg.terminal.decoy_perc_a
            )
            for m, lf, sc in zip(contig_models, tscore["lfdr"], tscore["score"]):
                m.evidence["terminal_lfdr"] = "" if lf != lf else round(float(lf), 6)
                m.evidence["terminal_score"] = "" if sc != sc else round(float(sc), 4)
            if tmodel.fitted or tmodel.note:
                rep = tmodel.report()
                rep["contig"] = contig
                model_reports.append(rep)
            ceiling = cfg.terminal.max_terminal_lfdr
            if ceiling is not None:
                keep, dropped = [], 0
                for m in contig_models:
                    v = m.evidence.get("terminal_lfdr", "")
                    novel = not int(m.evidence.get("end_is_annotated", 0) or 0)
                    if novel and v != "" and float(v) > ceiling:
                        dropped += 1
                        continue
                    keep.append(m)
                if dropped:
                    _log(cfg, f"[{contig}] terminal lFDR > {ceiling}: "
                              f"{dropped:,} novel 3' ends dropped")
                contig_models = keep

    # -- 8. second pass -----------------------------------------------------
    # The catalogue is filtered here, never rebuilt. Running it inside the
    # contig is not an optimisation: it needs `m.reads`, `reads` and the genome
    # handle, all of which exist only here, and gene-level questions ("the most
    # supported model of this gene") are contig-local because genes do not
    # cross contigs. It runs BEFORE bookkeeping so that model_of_read, the
    # read_stat shard and every invariant describe the catalogue that was
    # actually emitted.
    n_before = len(contig_models)
    sp_res = spmod.run_second_pass(contig_models, contig, ref, reads, cfg)
    contig_models = sp_res.models
    out.secondpass_audit = sp_res.audit
    sp_stats = sp_res.stats
    if cfg.secondpass.enabled and n_before != len(contig_models):
        _log(cfg, f"[{contig}] second pass: {n_before:,} -> "
                  f"{len(contig_models):,} models "
                  f"({sp_stats.get('monoexon_replaced', 0):,} mono-exon replaced, "
                  f"{sp_stats.get('inferred_models_created', 0):,} inferred, "
                  f"{sp_stats.get('end3_dropped', 0):,} end3_novel dropped)")

    # -- 8b. final consolidation (0.6.0) -----------------------------------
    # Same-chain 3' collapse, then the reference-free sub-chain filter. Like
    # the second pass it runs BEFORE bookkeeping, so every read table and the
    # per-cell matrix describe the consolidated catalogue.
    n_before_posthoc = len(contig_models)
    ph_res = phmod.run_posthoc(contig_models, contig, reads, cfg)
    contig_models = ph_res.models
    out.posthoc_audit = ph_res.audit
    ph_stats = ph_res.stats
    if cfg.posthoc.enabled and n_before_posthoc != len(contig_models):
        _log(cfg, f"[{contig}] post-hoc: {n_before_posthoc:,} -> "
                  f"{len(contig_models):,} models "
                  f"({ph_stats.get('end3_models_absorbed', 0):,} merged at 3' "
                  f"tolerance {cfg.posthoc.end3_tolerance} bp, "
                  f"{ph_stats.get('subchain_removed', 0):,} sub-chain fragments removed)")

    # -- 9. bookkeeping ----------------------------------------------------
    # Local model indices, not ids: the parent assigns ids once every contig
    # has reported, and substitutes them into the shards while stitching.
    model_of_read = np.full(reads.n, -1, np.int64)
    shard = pmod.ContigShard(os.path.join(shared.shard_dir, contig))
    for j, m in enumerate(contig_models):
        model_of_read[m.reads] = j
        if mol_index is not None:
            pairs = [(int(reads.cell[k]), int(reads.umi[k])) for k in m.reads]
            counts = collapse_umis(pairs, cfg.molecules.umi_hamming)
            per_cell: Dict[int, int] = defaultdict(int)
            for (cid, _u) in counts:
                per_cell[cid] += 1
            if per_cell:
                shard.add_counts(j, dict(per_cell))
    # Resolve the chain-table anchors against the FINAL catalogue. A group whose
    # anchor model did not survive has no model to report its chains under, so
    # the rows go with it rather than being re-pointed at a stranger.
    where = {id(m): j for j, m in enumerate(contig_models)}
    out.group_chain_rows = [
        (where[id(a)], ch, nr, nm, is_model)
        for a, ch, nr, nm, is_model in pending_chain_rows
        if id(a) in where
    ]

    validator.observe_contig(reads, contig_models)
    _write_read_shard(shard, bam, contig, cfg, reads, model_of_read)
    shard.close()
    out.shard_stem = shard.stem
    out.n_written = shard.n_written
    out.n_unassigned = shard.n_unassigned
    # read indices are contig-local; release them once bookkeeping is done
    for m in contig_models:
        m.reads = np.zeros(0, np.int64)

    n_spliced = int(np.sum(np.diff(reads.joff) > 0))
    validator.add_accounting({
        "admitted": reads.n,
        "spliced": n_spliced,
        "unspliced": reads.n - n_spliced,
        "reads_with_rejected_junction": int(np.sum(~read_ok & (np.diff(reads.joff) > 0))),
        # Rule J. reads_substituted MOVE onto an annotated model; they are not
        # lost, which is the difference between this and dropping a junction.
        "groups_substituted": sub_stats.get("groups_substituted", 0),
        "reads_substituted": sub_stats.get("reads_substituted", 0),
        "junctions_substituted": sub_stats.get("junctions_substituted", 0),
        "junctions_swallowing_exons": sub_stats.get("junctions_swallowing_exons", 0),
        "substitution_no_destination": sub_stats.get("no_destination", 0),
        # 0.1.16: the decoy set is the chain model's negative anchor. Relaxing
        # the junction gates shrinks it, and an empty one degrades calibration
        # silently -- so it is counted rather than assumed.
        "decoy_chain_groups": len(decoy_groups),
        "monoexon_reads_modelled": mono_stats["reads_in_models"],
        "monoexon_reads_rejected_low_support": mono_stats["reads_rejected_low_support"],
        "monoexon_reads_rejected_no_polya": mono_stats["reads_rejected_no_polya"],
        "monoexon_reads_rejected_internal_priming":
            mono_stats["reads_rejected_internal_priming"],
        "monoexon_models_demoted_by_short_reads":
            mono_stats.get("rejected_sr_continuation", 0),
        "monoexon_reads_folded_into_spliced":
            mono_stats.get("reads_demoted_fragment", 0),
        "monoexon_models_demoted_as_fragments":
            mono_stats.get("rejected_demoted_fragment", 0),
        "monoexon_models_rescued_by_tail": mono_stats.get("rescued_by_tail", 0),
        # 0.5.0: the second pass moves read mass too, so it is accounted for
        # beside every other rule that does.
        "secondpass_monoexon_replaced_reads": sp_stats.get(
            "monoexon_replaced_reads", 0),
        "secondpass_end3_dropped_reads": sp_stats.get("end3_dropped_reads", 0),
        "secondpass_end3_reads_folded": sp_stats.get(
            "end3_dropped_reads_folded", 0),
        "posthoc_end3_models_absorbed": ph_stats.get("end3_models_absorbed", 0),
        "posthoc_end3_reads_moved": ph_stats.get("end3_reads_moved", 0),
        "posthoc_subchain_models_removed": ph_stats.get("subchain_removed", 0),
        "posthoc_subchain_reads_folded": ph_stats.get("subchain_reads_folded", 0),
        "posthoc_subchain_reads_unassigned": ph_stats.get("subchain_reads_unassigned", 0),
        "secondpass_reads_unassigned": (
            sp_stats.get("end3_dropped_reads_unassigned", 0)
            + sp_stats.get("monoexon_removed_no_target_reads", 0)
            + sp_stats.get("monoexon_intergenic_removed_reads", 0)
        ),
    })

    out.models = contig_models
    out.model_reports = model_reports
    out.mono_reject_rows = mono_reject_rows
    out.substitution_rows = substitution_rows
    out.validator = validator
    out.merge_stats = dict(merge_stats)
    out.summary = {
        "reads": summ["reads"],
        "unspliced_reads": summ["mono_exonic_reads"],
        "unique_junctions": int(len(cur)),
        "junctions_kept": int(cur["keep"].sum()),
        "chain_groups": len(groups),
        "candidate_chains": int(candidate.sum()),
        "models": len(contig_models),
        "models_before_second_pass": n_before,
        "models_before_posthoc": n_before_posthoc,
        "posthoc": ph_stats,
        "models_intron_retention": n_ir,
        "models_readthrough": n_rt,
        "models_with_associated_transcript": n_named,
        "monoexon": mono_stats,
        "merge": merge_stats,
        "second_pass": sp_stats,
        "gates": gate_stats,
        "internal_deletions": gap_stats,
        "gene_conflicts": gene_stats.get("gene_conflicts", 0),
        "seconds": round(time.time() - t0, 1),
    }
    return out


def _write_read_shard(
    shard: pmod.ContigShard,
    bam,
    contig: str,
    cfg: Config,
    reads: ContigReads,
    model_of_read: np.ndarray,
) -> None:
    """One contig's read assignments, keyed on the contig-local model index.

    Read names are never held in memory: the BAM is re-walked with the same
    admission rule used during scanning and the n-th admitted read is paired
    with row ``n`` of the arrays -- exactly as ``ReadAssignmentWriter`` has
    always done.  The only change is that the model is named by an integer the
    parent resolves later, because ids do not exist yet.
    """
    groups: Dict[int, List[str]] = defaultdict(list)
    for i, r in enumerate(iter_alignments(bam, contig, cfg.gates)):
        if i >= reads.n:
            break
        mi = int(model_of_read[i])
        if mi < 0:
            shard.n_unassigned += 1
            continue
        shard.write_read(r.query_name, int(reads.end[i] - reads.start[i]), mi)
        groups[mi].append(r.query_name)
    for mi, names in groups.items():
        shard.write_group(mi, names)


def _second_pass_report(
    per_contig: Dict[str, Dict[str, object]],
    audit: Sequence[Dict[str, object]],
    cfg: Config,
) -> Dict[str, object]:
    """What the second pass removed, replaced and kept, and on what evidence.

    Written next to the invariants rather than left to be reconstructed from
    the audit table, because the number that matters -- how much of the
    catalogue a filter deleted -- should never need a script to find.
    """
    totals: Dict[str, int] = defaultdict(int)
    for summ in per_contig.values():
        for k, v in (summ.get("second_pass") or {}).items():
            try:
                totals[k] += int(v)
            except (TypeError, ValueError):
                continue
    classes = {
        c: int(totals.get(f"monoexon_class_{c}", 0))
        for c in spmod.MONOEXON_CLASSES
    }
    class_reads = {
        c: int(totals.get(f"monoexon_class_{c}_reads", 0))
        for c in spmod.MONOEXON_CLASSES
    }
    tiers = {
        t: int(totals.get(f"end3_kept_{t}", 0)) for t in spmod.END3_TIERS
    }
    kept = int(totals.get("end3_kept", 0))
    dropped = int(totals.get("end3_dropped", 0))
    return {
        "monoexon_classes": classes,
        "monoexon_class_reads": class_reads,
        "monoexon_classes_removed": sorted(cfg.secondpass.remove_classes),
        "monoexon_models_replaced": int(totals.get("monoexon_replaced", 0)),
        "monoexon_replaced_from_data": int(totals.get("monoexon_replaced_data", 0)),
        "monoexon_replaced_from_reference": int(
            totals.get("monoexon_replaced_reference", 0)),
        "monoexon_reads_replaced": int(totals.get("monoexon_replaced_reads", 0)),
        "inferred_models_created": int(totals.get("inferred_models_created", 0)),
        "monoexon_removed_no_inference_target": int(
            totals.get("monoexon_removed_no_inference_target", 0)),
        "end3_models_judged": kept + dropped,
        "end3_models_kept": kept,
        "end3_models_dropped": dropped,
        "end3_frac_dropped": round(dropped / max(kept + dropped, 1), 4),
        "end3_kept_by_tier": tiers,
        # The share kept only by the backstop is the honest measure of how much
        # of this category rests on "nothing contradicted it" rather than on
        # evidence. A high number here is not a pass -- it is the thing to look
        # at next.
        "end3_frac_kept_on_gene_dominance_alone": round(
            tiers.get("gene_dominant", 0) / max(kept, 1), 4),
        "end3_reads_dropped": int(totals.get("end3_dropped_reads", 0)),
        "end3_reads_folded": int(totals.get("end3_dropped_reads_folded", 0)),
        "end3_reads_unassigned": int(totals.get("end3_dropped_reads_unassigned", 0)),
        "audit_rows": len(audit),
        "audit_table": cfg.outpath(".secondpass.tsv") if audit else "",
    }


def _three_prime_report(models: Sequence[TranscriptModel], cfg) -> Dict[str, object]:
    """What the 3' end path actually did, in the run's own QC.

    New in 0.1.17.  The NSL1 failure -- a model extended 11.6 kb past every read
    supporting it -- was invisible in 0.1.16's report; finding it took a masking
    experiment and a per-locus hunt.  Two numbers make that class of failure
    self-reporting:

    ``max_abs_tts_shift``  the largest distance annotation moved any 3' end away
                           from its own reads.  Bounded by ``max_3p_fallback_dist``
                           by construction, so a value at the bound means the
                           bound is doing work.
    ``tail_molecule_frac`` how many of a peak's molecules carry a polyA tail in
                           their soft clip.  Recorded, not enforced -- this is the
                           distribution to look at before choosing a threshold.
    """
    if not models:
        return {}
    shifts = [abs(int(m.evidence.get("tts_shift_from_reads", 0) or 0)) for m in models]
    moved = [s for s in shifts if s > 0]
    unresolved = [m for m in models
                  if m.category == cfg.ends.unresolved_3p_category]
    refused = [m for m in models
               if any(str(f).startswith("3p_fallback_refused") for f in m.flags)]
    tails = [float(m.evidence.get("tail_molecule_frac", 0.0) or 0.0)
             for m in models if m.n_exons >= 1]
    tailed = [m for m in models
              if float(m.evidence.get("tail_molecule_frac", 0.0) or 0.0)
              >= cfg.ends.min_tail_molecule_frac]
    out: Dict[str, object] = {
        "max_3p_fallback_dist": cfg.ends.max_3p_fallback_dist,
        "models_moved_by_annotation": len(moved),
        "median_abs_tts_shift": float(np.median(moved)) if moved else 0.0,
        "max_abs_tts_shift": int(max(shifts)) if shifts else 0,
        "unresolved_3p_models": len(unresolved),
        "fallback_refused": len(refused),
        "polya_tail_measured": bool(cfg.gates.measure_polya_tail),
        "models_with_tail_support": len(tailed),
        "median_tail_molecule_frac": float(np.median(tails)) if tails else 0.0,
    }
    if shifts and out["max_abs_tts_shift"] > cfg.ends.max_3p_fallback_dist:
        # a shift past the bound means something moved an end outside _fallback
        out["WARNING"] = (
            f"a 3' end moved {out['max_abs_tts_shift']} bp from its reads, past "
            f"the {cfg.ends.max_3p_fallback_dist} bp bound -- not via the "
            f"annotation fallback, so it is a different path and a bug"
        )
    return out


def _novel_burden(models: Sequence[TranscriptModel]) -> Dict[str, object]:
    """How much novelty was called, and where.

    Over-calling transcriptional diversity is the failure mode this tool exists
    to avoid, so the burden is reported next to the invariants rather than left
    for the user to compute.
    """
    novel = [m for m in models if m.category in ("NIC", "NNC")]
    if not models:
        return {}
    by_cat: Dict[str, int] = defaultdict(int)
    for m in models:
        by_cat[m.category] += 1
    per_gene: Dict[str, int] = defaultdict(int)
    for m in novel:
        per_gene[m.gene_name or m.gene_id or "?"] += 1
    tot_reads = sum(m.n_reads for m in models) or 1
    lf = [m.lfdr for m in novel if m.lfdr == m.lfdr]
    return {
        "n_novel_models": len(novel),
        "frac_models_novel": round(len(novel) / len(models), 4),
        "frac_read_mass_novel": round(sum(m.n_reads for m in novel) / tot_reads, 4),
        "genes_with_novel_models": len(per_gene),
        "max_novel_models_in_one_gene": max(per_gene.values()) if per_gene else 0,
        "top_genes_by_novel_models": dict(
            sorted(per_gene.items(), key=lambda kv: -kv[1])[:15]
        ),
        "novel_models_at_lfdr_0.05": int(sum(1 for v in lf if v <= 0.05)),
        "novel_models_at_lfdr_0.01": int(sum(1 for v in lf if v <= 0.01)),
        # non-reference but well-understood classes, kept out of the novelty
        # burden because calling them "novel isoforms" is the overstatement
        "n_intron_retention": by_cat.get("IR", 0),
        "n_readthrough": by_cat.get("READTHROUGH", 0),
        "n_end3_annotated": by_cat.get("end3_annotated", 0),
        "n_end3_novel": by_cat.get("end3_novel", 0),
        "n_end5_alt": by_cat.get("end5_alt", 0),
        "n_end5_extended": by_cat.get("end5_extended", 0),
        "models_without_associated_transcript": int(
            sum(1 for m in models if not m.associated_tx)
        ),
        "models_with_ambiguous_association": int(
            sum(1 for m in models if (m.n_equally_close or 0) > 1)
        ),
    }


def _novel_category(g: "chainmod.ChainGroup") -> str:
    if g.n_novel_junctions == 0:
        return "NIC"          # novel combination of annotated junctions
    return "NNC"              # contains at least one novel junction


# ---------------------------------------------------------------------- #
def _chain_features(
    groups: Sequence["chainmod.ChainGroup"],
    reads: ContigReads,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    cage: Optional[TSSAtlas],
    cfg: Config,
    gene_support: Dict[str, int],
) -> pd.DataFrame:
    rows = []
    for g in groups:
        idx = g.reads
        three = int(np.median(reads.tts[idx]))
        five = five_prime_end(
            reads.tss[idx], three, g.strand, cfg.ends.five_prime_percentile
        )
        d_tts = d_tss = 1 << 16
        if g.gene_id and g.gene_id in ref.genes:
            d_tts = ref.nearest_tts(g.gene_id, three)[0]
            d_tss = ref.nearest_tss(g.gene_id, five)[0]
        g.tss_evidence_source = _tss_tier(
            g.contig, five, g.strand, d_tss, cage, cfg
        )
        # Annotation has always been allowed to keep a suffix; the atlas tiers
        # only do so when rescue is enabled.  With it off they are still
        # recorded on every model and no merge changes -- which is what makes
        # the effect of turning it on measurable against the same run.
        g.tss_evidence = g.tss_evidence_source == "annotation" or (
            g.tss_evidence_source != "none" and cfg.tss.rescue_suffix
        )
        pa = 0.0
        if genome is not None or atlas is not None:
            ev = polya_evidence(
                g.contig, three, g.strand, genome, atlas,
                motif_window=cfg.monoexon.polya_motif_window,
                perc_a_window=cfg.monoexon.perc_a_window,
                max_perc_a=cfg.monoexon.max_perc_a_downstream,
            )
            pa = 1.0 if ev["polya_supported"] else 0.0
        tot = gene_support.get(g.gene_id, 0) if g.gene_id else 0
        rows.append(
            {
                "gid": g.gid,
                "n_reads": g.n_reads,
                "n_mols": g.n_mols,
                "n_cells": g.n_cells,
                "gene_share": (g.support / tot) if tot else 1.0,
                "frac_junctions_annotated": g.frac_junctions_annotated,
                "min_junction_posterior": g.min_junction_posterior,
                "n_novel_junctions": g.n_novel_junctions,
                "n_decoy_junctions": g.n_decoy_junctions,
                "n_introns": len(g.chain),
                "dist_ann_tts": d_tts,
                "dist_ann_tss": d_tss,
                "polya_evidence": pa,
                "cage_support": (
                    float(TSS_TIER_ORDINAL.get(g.tss_evidence_source, 0))
                    if cfg.tss.use_as_chain_feature else 0.0
                ),
                "chain_annotated": g.annotated,
            }
        )
    # A contig whose admitted reads are ALL unspliced produces no chain groups,
    # and `pd.DataFrame([])` has no columns at all -- so every downstream column
    # access raises KeyError and the runner reports the contig as FAILED. It
    # surfaces on chrY (3 reads here), unplaced scaffolds, and any contig with
    # no spliced read. Pre-existing; fixed at the source rather than by guarding
    # each consumer, because the consumers are several and the next one only
    # shows up after the previous is patched.
    if not rows:
        return pd.DataFrame({c: pd.Series(dtype=t) for c, t in _CHAIN_FEATURE_COLS})
    return pd.DataFrame(rows)


#: columns of the chain feature frame, for the empty case. Kept beside the
#: builder so the two cannot drift apart unnoticed; the test asserts they match.
_CHAIN_FEATURE_COLS = (
    ("gid", "int64"), ("n_reads", "int64"), ("n_mols", "int64"),
    ("n_cells", "int64"), ("gene_share", "float64"),
    ("frac_junctions_annotated", "float64"), ("min_junction_posterior", "float64"),
    ("n_novel_junctions", "int64"), ("n_decoy_junctions", "int64"),
    ("n_introns", "int64"), ("dist_ann_tts", "float64"), ("dist_ann_tss", "float64"),
    ("polya_evidence", "float64"), ("cage_support", "float64"),
    ("chain_annotated", "bool"),
)


#: Evidence tiers for a model 5' end, strongest first.  The ordinal is what the
#: chain-level novelty model sees; the string is what the table reports.
TSS_TIER_ORDINAL: Dict[str, int] = {
    "none": 0, "cage_peak": 1, "cage_reptss": 2, "annotation": 3,
}


def _tss_tier(
    contig: str,
    five: int,
    strand: str,
    d_tss: float,
    cage: Optional[TSSAtlas],
    cfg: Config,
) -> str:
    """Which tier of evidence supports this 5' end being a real start site.

    First hit wins, strongest first, so the reported source is unambiguous.
    Annotation outranks the atlas because a GENCODE TSS is a curated claim
    about this gene, while a CAGE peak is a claim about the locus pooled over
    ~1,800 samples that do not include yours.
    """
    if abs(d_tss) <= cfg.ends.max_5p_diff:
        return "annotation"
    if cage is None or not cfg.tss.enabled:
        return "none"
    ev = cage.evidence(
        contig, int(five), strand,
        peak_slack=cfg.tss.peak_slack,
        reptss_window=cfg.tss.max_dist_to_reptss,
    )
    d_rep = ev["dist_to_cage_reptss"]
    if d_rep is not None and d_rep <= cfg.tss.max_dist_to_reptss:
        return "cage_reptss"
    if ev["in_cage_peak"]:
        return "cage_peak"
    return "none"


# ---------------------------------------------------------------------- #
def _write_all(
    cfg: Config,
    models: List[TranscriptModel],
    cell_counts: Dict[str, Dict[int, int]],
    group_chain_rows: List[Tuple[str, str, int, int, bool]],
    substitution_rows: List[Tuple],
    junction_frames: List[pd.DataFrame],
    mol_index: Optional[MoleculeIndex],
    ref: ReferenceIndex,
) -> None:
    o = cfg.output
    # 0.4.2: the SQANTI3 labelling, written into m.evidence only. It runs here,
    # after every model is final, so it can never influence what was emitted --
    # which is the point: the two labellings have to sit on the same catalogue
    # for the difference between them to mean anything.
    if getattr(cfg, "sqanti_classification", True):
        from . import structural as smod
        smod.annotate_models(models, ref, cfg)
    if o.write_gtf:
        omod.write_gtf(models, cfg.outpath(".gtf"))
    if o.write_gff:
        omod.write_gff(models, cfg.outpath(".gff"))
    if o.write_models_table:
        omod.write_models_table(models, cfg.outpath(".models.tsv"))
    if o.write_abundance:
        omod.write_abundance(models, cfg.outpath(".abundance.txt"),
                             ["generated by flightcollapse"])
        omod.write_flnc_count(models, cfg.outpath(".flnc_count.txt"))
    if junction_frames:
        jf = pd.concat(junction_frames, ignore_index=True)
        jf.to_csv(cfg.outpath(".junctions.tsv.gz"), sep="\t", index=False)
    if o.write_group_chains and group_chain_rows:
        omod.write_group_chains(cfg.outpath(".group_chains.tsv"), group_chain_rows)
    if substitution_rows:
        # One row per moved junction. A read move should be auditable from a
        # file rather than reconstructed by diffing two catalogues.
        with open(cfg.outpath(".junction_substitutions.tsv"), "w") as fh:
            fh.write("contig\tstrand\tgid\ttarget_gid\tn_reads\tdonor\tacceptor\t"
                     "shift\texons_swallowed\tnovel_sr\tcompetitor_sr\treplaced_by\n")
            for contig, strand, gid, tgt, nr, sp in substitution_rows:
                fh.write(f"{contig}\t{strand}\t{gid}\t{tgt}\t{nr}\t"
                         f"{sp.junction[0]}\t{sp.junction[1]}\t{sp.shift}\t"
                         f"{sp.swallowed}\t{sp.novel_sr}\t{sp.competitor_sr}\t"
                         + ";".join(f"{x}-{y}" for x, y in sp.replaces) + "\n")
    if cfg.molecules.write_matrix and mol_index is not None:
        omod.write_molecule_matrix(models, cell_counts, mol_index.barcodes, cfg)
    cfg.dump(cfg.outpath(".config.json"))


# ---------------------------------------------------------------------- #
def _terminal_features(models: Sequence[TranscriptModel]) -> pd.DataFrame:
    """One row per model, carrying everything the terminal score reads.

    Every field here is already computed during model building; this only
    gathers them. That is deliberate -- a scoring level that needed new
    measurements would also need new failure modes.
    """
    rows = []
    for m in models:
        ev = m.evidence
        rows.append({
            "n_reads": m.n_reads,
            "n_mols": m.n_mols,
            "n_cells": m.n_cells,
            "three_prime_dispersion": ev.get("three_prime_dispersion", 0.0),
            "chain_share": ev.get("chain_share", 1.0),
            "ratio_to_dominant_peak": ev.get("ratio_to_dominant_peak", 1.0),
            "n_3p_peaks": ev.get("n_3p_peaks", 1),
            "tail_molecule_frac": ev.get("tail_molecule_frac", 0.0),
            "median_tail_len": ev.get("median_tail_len", 0.0),
            "perc_a_downstream": ev.get("perc_a_downstream", 0.0),
            # 0.5.0. The two measurements that actually separate internal
            # priming from cleavage, fed to the calibrated model rather than
            # used only by the second pass -- the terminal score already has
            # the right positives (annotated TTS / atlas hits) and the right
            # decoys (A-rich downstream, nothing annotated), and these are the
            # features that null was missing.
            "median_tail_excess": _f(ev.get("median_tail_excess"), 0.0),
            "frac_reads_nontemplated": _f(ev.get("frac_reads_nontemplated"), 0.0),
            "genomic_a_run_3p": _f(ev.get("genomic_a_run_3p"), 0.0),
            "dse_gu_frac": _f(ev.get("dse_gu_frac"), 0.0),
            "dse_gu_minus_a": _f(ev.get("dse_gu_minus_a"), 0.0),
            "upstream_u_frac": _f(ev.get("upstream_u_frac"), 0.0),
            "polya_motif_found": int(bool(ev.get("polya_motif_found", 0))),
            "dist_ann_tts": ev.get("dist_ann_tts", 1e6),
            "dist_polya_site": ev.get("dist_polya_site", 1e6),
            "sr_3p_step_ratio": ev.get("sr_3p_step_ratio", 0.0),
            "end_is_annotated": bool(int(ev.get("end_is_annotated", 0) or 0)),
            "end_at_polya_site": bool(int(ev.get("end_at_polya_site", 0) or 0)),
        })
    return pd.DataFrame(rows)


def _f(v: object, default: float) -> float:
    """``evidence`` mixes "", None and numbers; the feature frame needs floats."""
    if v is None or v == "":
        return default
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return default if f != f else f
