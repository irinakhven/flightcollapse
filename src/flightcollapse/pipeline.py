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
8. **outputs**     GTF/GFF, read_stat, group, abundance, molecule matrix
9. **validation**  the invariants, checked before anything is declared done
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
from . import junctions as jmod
from . import output as omod
from .config import Config
from .ends import annotate_similarity, build_models_for_group, five_prime_end
from .genome import Genome, PolyASiteAtlas
from .model import TranscriptModel, polya_evidence
from .molecules import NO_CELL, NO_UMI, MoleculeIndex, collapse_umis
from .monoexon import collapse_monoexonic
from .reads import ContigReads, gate_diagnostics, resolve_contig
from .reference import ReferenceIndex
from .scoring import score_chains, score_junctions
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
    contigs = select_contigs(cfg, bam)
    _log(cfg, f"[run] {len(contigs)} contigs: {', '.join(contigs[:8])}"
              f"{' ...' if len(contigs) > 8 else ''}")

    validator = Validator()
    writer = omod.ReadAssignmentWriter(cfg)
    ids = omod.IdAssigner(cfg.output.id_style)
    all_models: List[TranscriptModel] = []
    cell_counts: Dict[str, Dict[int, int]] = {}
    group_chain_rows: List[Tuple[str, str, int, int, bool]] = []
    junction_frames: List[pd.DataFrame] = []
    model_reports: List[Dict[str, object]] = []
    per_contig: Dict[str, Dict[str, object]] = {}
    failed: Dict[str, str] = {}

    for contig in contigs:
        if resolve_contig(bam, contig) is None:
            _log(cfg, f"[{contig}] not in BAM, skipping")
            continue
        try:
            res = _run_contig(
                cfg, bam, contig, ref, genome, atlas, mol_index, validator,
                writer, ids, all_models, cell_counts, group_chain_rows,
                junction_frames, model_reports, sj,
            )
        except Exception as exc:
            # Whatever went wrong on this contig, the other 23 are still worth
            # writing out. Losing hours of completed work to a late failure is
            # a worse outcome than an incomplete result that says so.
            failed[contig] = f"{type(exc).__name__}: {exc}"
            _log(cfg, f"[{contig}] FAILED, continuing with the remaining contigs\n"
                      f"        {exc}")
            continue
        if res is not None:
            per_contig[contig] = res

    writer.close()
    bam.close()

    _write_all(cfg, all_models, cell_counts, group_chain_rows, junction_frames,
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

    if genome is not None:
        genome.close()
    if failed:
        _log(cfg, "\n[run] contigs that failed: "
                  + ", ".join(f"{k} ({v})" for k, v in failed.items())
                  + "\n[run] everything else was written normally.")
    if cfg.strict_invariants:
        validator.raise_on_failure()
    return report


# ---------------------------------------------------------------------- #
def _run_contig(
    cfg: Config,
    bam,
    contig: str,
    ref: ReferenceIndex,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    mol_index: Optional[MoleculeIndex],
    validator: Validator,
    writer: omod.ReadAssignmentWriter,
    ids: omod.IdAssigner,
    all_models: List[TranscriptModel],
    cell_counts: Dict[str, Dict[int, int]],
    group_chain_rows: List[Tuple[str, str, int, int, bool]],
    junction_frames: List[pd.DataFrame],
    model_reports: List[Dict[str, object]],
    sj=None,
) -> Optional[Dict[str, object]]:
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
        return None
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
    junction_frames.append(cur)
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

    cfeat = _chain_features(groups, reads, ref, genome, atlas, cfg, gene_support)
    if decoy_groups:
        dfeat = _chain_features(
            decoy_groups, reads, ref, genome, atlas, cfg, gene_support
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
    model_of_read = np.full(reads.n, -1, np.int64)
    pending_chain_rows: List[Tuple[int, str, int, int, bool]] = []

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
        )
        n_super = len(ref.annotated_superchains(g.gene_id, g.chain))
        n_gap = int(np.sum(unexplained_gap[idx])) if unexplained_gap.size else 0
        for m in ms:
            m.lfdr = lfdr_by_gid.get(g.gid, float("nan"))
            m.posterior_real = post_by_gid.get(g.gid, float("nan"))
            if not g.annotated:
                m.category = _novel_category(g)
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
        anchor = len(contig_models)
        contig_models.extend(ms)
        for k in member_gids:
            chg = groups[k]
            pending_chain_rows.append(
                (anchor, omod.chain_repr(chg.chain), chg.n_reads, chg.n_mols, k == g.gid)
            )

    # -- 7. mono-exonic track ---------------------------------------------
    mono_models, mono_stats = collapse_monoexonic(
        reads, contig, ref, genome, atlas, cfg.monoexon, read_ok,
        cfg.ends.five_prime_percentile, cfg.molecules.umi_hamming,
    )
    contig_models.extend(mono_models)
    _log(cfg, f"[{contig}] mono-exon: {mono_stats['monoexon_reads']:,} reads -> "
              f"{len(mono_models):,} models "
              f"(3'UTR {mono_stats['models_utr3']}, "
              f"internal {mono_stats['models_internal']}, "
              f"intergenic {mono_stats['models_intergenic']})")

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

    # -- 8. bookkeeping ----------------------------------------------------
    base = len(all_models)
    ids.assign(contig_models)
    for j, m in enumerate(contig_models):
        model_of_read[m.reads] = base + j
        if mol_index is not None:
            pairs = [(int(reads.cell[k]), int(reads.umi[k])) for k in m.reads]
            counts = collapse_umis(pairs, cfg.molecules.umi_hamming)
            per_cell: Dict[int, int] = defaultdict(int)
            for (cid, _u) in counts:
                per_cell[cid] += 1
            cell_counts[m.model_id] = dict(per_cell)
    for anchor, ch, nr, nm, is_model in pending_chain_rows:
        group_chain_rows.append((contig_models[anchor].model_id, ch, nr, nm, is_model))
    all_models.extend(contig_models)

    validator.observe_contig(reads, contig_models)
    writer.write_contig(
        bam, contig, cfg, reads, model_of_read,
        [m.model_id for m in all_models],
    )
    # read indices are contig-local; release them once bookkeeping is done
    for m in contig_models:
        m.reads = np.zeros(0, np.int64)

    n_spliced = int(np.sum(np.diff(reads.joff) > 0))
    validator.add_accounting({
        "admitted": reads.n,
        "spliced": n_spliced,
        "unspliced": reads.n - n_spliced,
        "reads_with_rejected_junction": int(np.sum(~read_ok & (np.diff(reads.joff) > 0))),
        # 0.1.16: the decoy set is the chain model's negative anchor. Relaxing
        # the junction gates shrinks it, and an empty one degrades calibration
        # silently -- so it is counted rather than assumed.
        "decoy_chain_groups": len(decoy_groups),
        "monoexon_reads_modelled": mono_stats["reads_in_models"],
        "monoexon_reads_rejected_low_support": mono_stats["reads_rejected_low_support"],
        "monoexon_reads_rejected_no_polya": mono_stats["reads_rejected_no_polya"],
        "monoexon_reads_rejected_internal_priming":
            mono_stats["reads_rejected_internal_priming"],
    })

    out = {
        "reads": summ["reads"],
        "unspliced_reads": summ["mono_exonic_reads"],
        "unique_junctions": int(len(cur)),
        "junctions_kept": int(cur["keep"].sum()),
        "chain_groups": len(groups),
        "candidate_chains": int(candidate.sum()),
        "models": len(contig_models),
        "models_intron_retention": n_ir,
        "models_readthrough": n_rt,
        "models_with_associated_transcript": n_named,
        "monoexon": mono_stats,
        "merge": merge_stats,
        "gates": gate_stats,
        "internal_deletions": gap_stats,
        "gene_conflicts": gene_stats.get("gene_conflicts", 0),
        "seconds": round(time.time() - t0, 1),
    }
    return out


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
            g.tss_evidence = abs(d_tss) <= cfg.ends.max_5p_diff
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
    ("polya_evidence", "float64"), ("chain_annotated", "bool"),
)


# ---------------------------------------------------------------------- #
def _write_all(
    cfg: Config,
    models: List[TranscriptModel],
    cell_counts: Dict[str, Dict[int, int]],
    group_chain_rows: List[Tuple[str, str, int, int, bool]],
    junction_frames: List[pd.DataFrame],
    mol_index: Optional[MoleculeIndex],
    ref: ReferenceIndex,
) -> None:
    o = cfg.output
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
    if cfg.molecules.write_matrix and mol_index is not None:
        omod.write_molecule_matrix(models, cell_counts, mol_index.barcodes, cfg)
    cfg.dump(cfg.outpath(".config.json"))
