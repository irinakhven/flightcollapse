"""Built-in validation invariants.

These are the checks that would have caught the ``isoseq collapse``
representative-selection bug on the day it happened.  They are cheap, they run
on every build, and by default a violation fails the run rather than being
written into a log nobody reads.

1. **Step-0 invariant.** The fraction of read mass on mono-exon models must
   match the fraction of genuinely unspliced aligned reads.  This number comes
   straight from the CIGARs and is independent of every downstream stage, so a
   gap can only have been created by the collapse.  Observed on BD144_kin:
   9.2% of reads unspliced, 57.9% of read mass on single-exon models.
2. **``model_is_min`` must not be ~1.0.**  If the emitted model is always the
   shortest member of its own read group, the representative rule is inverted.
3. **median log2(read span / model span) ~ 0** within groups.
4. **Reads restructured by merge**, reported explicitly.  The damage of a bad
   merge is not the debris chain's own reads, it is the parent read mass
   dragged onto a truncated model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .model import TranscriptModel
from .reads import ContigReads


@dataclass
class InvariantResult:
    name: str
    passed: bool
    observed: float
    expected: Optional[float]
    message: str

    def to_dict(self) -> Dict[str, object]:
        return {
            "invariant": self.name,
            "passed": self.passed,
            "observed": None if self.observed != self.observed else round(self.observed, 5),
            "expected": self.expected,
            "message": self.message,
        }


class Validator:
    """Accumulates per-contig evidence and renders a verdict at the end."""

    def __init__(self) -> None:
        self.total_reads = 0
        self.unspliced_reads = 0
        self.read_mass_on_monoexon_models = 0
        self.read_mass_on_models = 0
        self.model_is_min = 0
        self.model_is_max = 0
        self.n_multiexon_models = 0
        self.span_ratios: List[float] = []
        self.exon_deltas: List[float] = []
        self.merge_stats: Dict[str, int] = {}
        self.accounting: Dict[str, int] = {}
        self.extra: Dict[str, object] = {}

    # ------------------------------------------------------------------ #
    def observe_contig(
        self, reads: ContigReads, models: Sequence[TranscriptModel]
    ) -> None:
        n_intr = np.diff(reads.joff)
        self.total_reads += reads.n
        self.unspliced_reads += int(np.sum(n_intr == 0))

        read_span = (reads.end - reads.start).astype(float)
        for m in models:
            if m.reads.size == 0:
                continue
            self.read_mass_on_models += m.n_reads
            if m.n_exons == 1:
                self.read_mass_on_monoexon_models += m.n_reads
                continue
            self.n_multiexon_models += 1
            rex = n_intr[m.reads] + 1
            if m.n_exons <= int(rex.min()):
                self.model_is_min += 1
            if m.n_exons >= int(rex.max()):
                self.model_is_max += 1
            self.exon_deltas.append(float(m.n_exons - np.median(rex)))
            self.span_ratios.append(
                float(np.median(np.log2((read_span[m.reads] + 1.0) / (m.span + 1.0))))
            )

    def add_accounting(self, stats: Dict[str, int]) -> None:
        for k, v in stats.items():
            self.accounting[k] = self.accounting.get(k, 0) + int(v)

    def add_merge_stats(self, stats: Dict[str, int]) -> None:
        for k, v in stats.items():
            self.merge_stats[k] = self.merge_stats.get(k, 0) + int(v)

    # ------------------------------------------------------------------ #
    def results(self) -> List[InvariantResult]:
        out: List[InvariantResult] = []

        frac_unspliced = self.unspliced_reads / max(self.total_reads, 1)
        frac_mono_mass = self.read_mass_on_monoexon_models / max(self.read_mass_on_models, 1)
        ok = frac_mono_mass <= max(2.0 * frac_unspliced, frac_unspliced + 0.02)
        out.append(
            InvariantResult(
                "monoexon_read_mass_matches_unspliced_reads",
                ok,
                frac_mono_mass,
                round(frac_unspliced, 5),
                (
                    f"{frac_mono_mass:.1%} of assigned read mass sits on single-exon "
                    f"models but only {frac_unspliced:.1%} of aligned reads are "
                    f"genuinely unspliced. A gap this size can only be created "
                    f"between alignment and model assignment."
                    if not ok
                    else f"{frac_mono_mass:.1%} vs {frac_unspliced:.1%} unspliced reads"
                ),
            )
        )

        n = max(self.n_multiexon_models, 1)
        f_min = self.model_is_min / n
        out.append(
            InvariantResult(
                "model_is_not_always_the_shortest_member",
                f_min < 0.95 or self.n_multiexon_models == 0,
                f_min,
                "< 0.95",
                (
                    "the emitted model is the shortest member of its own read group "
                    "in essentially every case -- the representative rule is inverted"
                    if f_min >= 0.95
                    else f"model_is_min = {f_min:.3f}, model_is_max = {self.model_is_max / n:.3f}"
                ),
            )
        )

        med_span = float(np.median(self.span_ratios)) if self.span_ratios else 0.0
        out.append(
            InvariantResult(
                "median_log2_span_ratio_near_zero",
                abs(med_span) < 1.0,
                med_span,
                "|.| < 1.0",
                (
                    f"median log2(read span / model span) = {med_span:.2f}; reads are "
                    f"systematically {'longer' if med_span > 0 else 'shorter'} than "
                    f"the models they were assigned to"
                    if abs(med_span) >= 1.0
                    else f"median log2(read span / model span) = {med_span:.2f}"
                ),
            )
        )

        # the same comparison in the other direction: mono-exon models should
        # not capture *far less* than the unspliced read mass either. Reads that
        # vanish are as much a modelling failure as reads that pile onto a stub,
        # they are just a quieter one.
        mono_modelled = self.accounting.get("monoexon_reads_modelled", 0)
        n_unspliced = self.accounting.get("unspliced", self.unspliced_reads)
        cap = mono_modelled / max(n_unspliced, 1)
        out.append(
            InvariantResult(
                "unspliced_reads_are_accounted_for",
                cap >= 0.25 or n_unspliced == 0,
                cap,
                ">= 0.25",
                (
                    f"only {cap:.1%} of genuinely unspliced reads ended on a "
                    f"mono-exon model; the rest were discarded. Check the "
                    f"rejection reasons in read_accounting -- if they are "
                    f"dominated by no_polya in repeat-rich regions that is "
                    f"expected, otherwise the mono-exon track is too strict."
                    if cap < 0.25
                    else f"{cap:.1%} of unspliced reads modelled"
                ),
            )
        )

        coll = int(self.extra.get("model_id_collisions", 0))
        out.append(
            InvariantResult(
                "model_ids_are_unique",
                coll == 0,
                float(coll),
                "0",
                (f"{coll:,} models were emitted with a transcript_id already in use. "
                 f"A duplicate id silently corrupts read_stat, group.txt and the "
                 f"count matrix; they have been suffixed .dupN and flagged."
                 if coll else "no duplicate transcript ids"),
            )
        )

        med_exon = float(np.median(self.exon_deltas)) if self.exon_deltas else 0.0
        out.append(
            InvariantResult(
                "model_exon_count_matches_reads",
                med_exon >= -1.0,
                med_exon,
                ">= -1",
                f"median(model exons - median read exons) = {med_exon:.2f}",
            )
        )
        return out

    # ------------------------------------------------------------------ #
    def report(self) -> Dict[str, object]:
        res = self.results()
        return {
            "aligned_reads": self.total_reads,
            "unspliced_reads": self.unspliced_reads,
            "frac_reads_unspliced": round(self.unspliced_reads / max(self.total_reads, 1), 5),
            "reads_assigned_to_models": self.read_mass_on_models,
            "reads_on_monoexon_models": self.read_mass_on_monoexon_models,
            "multiexon_models": self.n_multiexon_models,
            "frac_model_is_min": round(self.model_is_min / max(self.n_multiexon_models, 1), 5),
            "frac_model_is_max": round(self.model_is_max / max(self.n_multiexon_models, 1), 5),
            "merge": self.merge_stats,
            "read_accounting": self.accounting,
            "invariants": [r.to_dict() for r in res],
            "all_passed": all(r.passed for r in res),
            **self.extra,
        }

    def raise_on_failure(self) -> None:
        bad = [r for r in self.results() if not r.passed]
        if bad:
            msg = "\n".join(f"  - {r.name}: {r.message}" for r in bad)
            raise RuntimeError(f"validation invariants failed:\n{msg}")


def render_markdown(report: Dict[str, object]) -> str:
    lines = ["# flightcollapse QC report", ""]
    lines.append(f"- aligned reads: **{report.get('aligned_reads', 0):,}**")
    lines.append(
        f"- genuinely unspliced (no `N` in CIGAR): "
        f"**{report.get('frac_reads_unspliced', 0):.2%}**"
    )
    lines.append(f"- reads assigned to models: **{report.get('reads_assigned_to_models', 0):,}**")
    lines.append(
        f"- multi-exon models where the model is the shortest / longest member of "
        f"its own group: **{report.get('frac_model_is_min', 0):.2f}** / "
        f"**{report.get('frac_model_is_max', 0):.2f}**"
    )
    lines.append(
        f"- reads on single-exon models: **{report.get('reads_on_monoexon_models', 0):,}**"
    )
    acc = report.get("read_accounting", {}) or {}
    if acc:
        lines += ["", "## Read accounting", ""]
        for k, v in acc.items():
            lines.append(f"- `{k}`: {v:,}")
    jc = report.get("junction_curation", {}) or {}
    if jc:
        lines += ["", "## Junction curation", ""]
        lines.append(f"- novel junctions: **{jc.get('novel_junctions', 0):,}**, "
                     f"kept **{jc.get('novel_kept', 0):,}**")
        mr, mk = jc.get("median_reads_of_rejected"), jc.get("median_reads_of_kept")
        if mr is not None and mk is not None:
            lines.append(f"- median reads, rejected vs kept: **{mr:.0f}** vs **{mk:.0f}**"
                         + ("  <- rejected junctions carry MORE evidence than kept "
                            "ones; a gate is mis-set" if mr > mk else ""))
        ml = jc.get("median_lfdr_of_rejected")
        if ml is not None:
            lines.append(f"- median lFDR of rejected: **{ml:.3f}**"
                         + ("  <- the calibrated model considers these real"
                            if ml < 0.05 else ""))
        by = jc.get("rejected_by_reason", {}) or {}
        if by:
            tot = sum(by.values()) or 1
            lines += ["", "| rule | junctions rejected | share |", "|---|---:|---:|"]
            for k, v in sorted(by.items(), key=lambda x: -x[1]):
                lines.append(f"| `{k}` | {v:,} | {100 * v / tot:.0f}% |")
    tp = report.get("three_prime", {}) or {}
    if tp:
        lines += ["", "## 3' ends", ""]
        lines.append(
            f"- annotation moved **{tp.get('models_moved_by_annotation', 0):,}** "
            f"model ends off their own reads; median "
            f"**{tp.get('median_abs_tts_shift', 0):.0f} bp**, max "
            f"**{tp.get('max_abs_tts_shift', 0):,} bp** "
            f"(bound {tp.get('max_3p_fallback_dist', 0):,} bp)"
        )
        lines.append(
            f"- **{tp.get('fallback_refused', 0):,}** peaks kept their observed "
            f"coordinate because the nearest annotated end was past that bound, "
            f"emitted as **{tp.get('unresolved_3p_models', 0):,}** unresolved models"
        )
        if tp.get("polya_tail_measured"):
            lines.append(
                f"- polyA tail read off the 3' soft clip: median "
                f"**{tp.get('median_tail_molecule_frac', 0):.2f}** of a model's "
                f"molecules carry one; **{tp.get('models_with_tail_support', 0):,}** "
                f"models clear the (currently unenforced) tail threshold"
            )
        else:
            lines.append("- polyA tail measurement was **off** "
                         "(`gates.measure_polya_tail=false`)")
        if tp.get("WARNING"):
            lines.append(f"- **WARNING**: {tp['WARNING']}")
    merge = report.get("merge", {}) or {}
    if merge:
        lines += ["", "## Suffix resolution", ""]
        for k, v in merge.items():
            lines.append(f"- `{k}`: {v:,}")
    nb = report.get("novel_model_burden") or {}
    if nb:
        lines += ["", "## Novelty burden", ""]
        lines.append(f"- novel (NIC/NNC) models: **{nb.get('n_novel_models', 0):,}** "
                     f"({nb.get('frac_models_novel', 0):.1%} of models, "
                     f"{nb.get('frac_read_mass_novel', 0):.1%} of read mass)")
        lines.append(f"- genes with novel models: {nb.get('genes_with_novel_models', 0):,}"
                     f" (max {nb.get('max_novel_models_in_one_gene', 0)} in one gene)")
        lines.append(f"- surviving a chain local-FDR ceiling of 0.05 / 0.01: "
                     f"{nb.get('novel_models_at_lfdr_0.05', 0):,} / "
                     f"{nb.get('novel_models_at_lfdr_0.01', 0):,}")
    sr = report.get("short_read_junctions") or {}
    if sr:
        lines += ["", "## Short-read validation", ""]
        lines.append(
            f"- coordinate check: **{sr.get('sjdb_found_in_reference', 0):,}/"
            f"{sr.get('sjdb_junctions', 0):,}** of STAR's annotated junctions are "
            f"annotated here ({sr.get('agreement', 0):.1%})"
        )
        lines.append(
            f"- short-read junctions not in the annotation: "
            f"**{sr.get('novel_junctions', 0):,}**"
        )
        v = sr.get("validation") or {}
        if v:
            lines.append(
                f"- anchors used / held out: {v.get('n_anchor_junctions', 0):,} / "
                f"**{v.get('n_holdout_junctions', 0):,}**"
            )
            lines.append("")
            lines.append("| lFDR | recall on held-out real junctions | novel junctions with no short-read support kept |")
            lines.append("|---|---|---|")
            for t in ("0.01", "0.05", "0.1"):
                r = v.get(f"recall_at_lfdr_{t}")
                u = v.get(f"unsupported_kept_at_lfdr_{t}")
                if r is None:
                    continue
                lines.append(
                    f"| {t} | {r:.1%} | {'-' if u is None else f'{u:.1%}'} |"
                )
            lines.append("")
            lines.append(
                "Recall is measured on junctions an orthogonal library confirmed "
                "and the calibration never saw, so it is not circular. The right "
                "column is the same threshold applied to novel junctions the "
                "short reads did not confirm; the gap between the columns is what "
                "the score is buying."
            )

    lines += ["", "## Invariants", "", "| invariant | ok | observed | expected |", "|---|---|---|---|"]
    for r in report.get("invariants", []):
        lines.append(
            f"| {r['invariant']} | {'PASS' if r['passed'] else 'FAIL'} | "
            f"{r['observed']} | {r['expected']} |"
        )
    lines += ["", "## Notes", ""]
    for r in report.get("invariants", []):
        if not r["passed"]:
            lines.append(f"- **{r['invariant']}**: {r['message']}")
    return "\n".join(lines) + "\n"
