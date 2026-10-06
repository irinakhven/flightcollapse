"""The second pass.  New in 0.5.0.

The first pass is unchanged: the same reads, the same junction curation, the
same suffix resolution, the same 3'-peak clustering, the same models.  This
stage runs afterwards, on the catalogue the first pass produced, and it may
only do three things -- keep a model, remove it, or replace it with an
explicitly inferred one.  It never builds a structure from scratch.

That separation is the point.  A filter folded into model construction cannot
be turned off, cannot be compared against its own absence, and cannot report
what it cost.  This one writes ``<prefix>.secondpass.tsv`` with one row per
decision, and the QC report carries the class and tier histograms next to the
invariants.

--------------------------------------------------------------------------
Rule 1 -- mono-exonic models judged against their gene's exon structure
--------------------------------------------------------------------------

The mono-exon track is deliberately blind to the gene it lands in beyond
"terminal exon or not".  That is the right primitive for an oligo-dT library,
where a long 3'UTR the RT did not traverse leaves exactly a single-exon pile
with a correct 3' end -- but it also means a single-exon model in a gene with
no single-exon isoform is emitted with no one asking whether such a transcript
could exist.

Five classes, decided from the annotation alone:

======  ==================================================  ========
class   gene / model relationship                           verdict
======  ==================================================  ========
``A``   gene is mono-exonic only                            keep
``B``   gene has both mono- and multi-exon isoforms         keep
``C``   multi-exon gene; the model spans an intron          replace
``D``   multi-exon gene; the model lies inside one exon     replace
``E``   multi-exon gene; anything else (e.g. wholly         replace
        intronic, or running off the gene)
``X``   no gene (intergenic) -- out of scope for A-E        policy
======  ==================================================  ========

"Spans an intron" uses *constitutively* intronic sequence: inside the gene's
span and exonic in no isoform.  Using any one transcript's introns instead
would call an ordinary alternative-exon model a read-through, because one
isoform's exon is another's intron all over the genome.

A replaced model does not simply vanish.  The gene's most-supported isoform is
recorded in its place and marked ``INFERRED`` -- the category, a flag, and an
``inferred=1`` evidence column all say the same thing three times, because the
one failure this must not have is an inferred model being read as an observed
one.  Its reads move onto it, so the read accounting stays exact.

--------------------------------------------------------------------------
Rule 2 -- ``end3_novel`` models must earn their 3' end
--------------------------------------------------------------------------

On SIRV Set 4 this category was wrong 12/12 in the reference arm and 11/11 in
the masked arm.  The mechanism is internal priming: oligo-dT anneals to a
genomic A-stretch inside the transcript, and the resulting pile accumulates
support exactly like a real cleavage site -- reproducibly, because it is driven
by genomic sequence.  Three consequences follow, and the third is the one that
catches people out:

* a read floor cannot separate them (false ends reach p90 548 reads while true
  ones go down to 5);
* ``tail_molecule_frac`` cannot separate them (0.96 vs 1.00 on BD176c), because
  it scores an A-run without asking what the genome says underneath it;
* **cross-sample reproducibility cannot separate them either.**  An internal
  priming site recurs in every library made from the same genome.  A model
  called in all three HEK samples is not thereby real.

Six tiers, strongest first.  The first one that speaks decides, and the tier is
written on the model so a downstream filter can pick its own operating point.

``atlas``          the end sits on a catalogued polyA site
``nontemplated``   the reads' own tails are longer than the genome's A-run
``composition``    hexamer + a U/GU-rich downstream element + a tight peak
``short_read``     short-read coverage demonstrably drops across the end
``calibrated``     the terminal local FDR clears a ceiling (off by default)
``gene_dominant``  the backstop: sole or best-supported model of its gene
------------------ ---------------------------------------------------------

``nontemplated`` is the only per-*molecule* test in the list, and the only one
that is new information rather than a rearrangement of what 0.4.6 already had.
A genuine polyA tail is not in the genome: the A-run in a read's 3' soft clip
runs far past the genomic A-run at that coordinate.  An internally primed read
has no such excess, because the A's it was primed on are templated and the
aligner consumes them.  Everything needed was already being measured --
``reads.tail_len`` per read since 0.1.17, ``genome.fetch`` since the beginning
-- and nothing compared the two.

``gene_dominant`` is last on purpose.  Used first it selects against the
biology: a genuine proximal polyA site usually carries a minority of a gene's
molecules, while a strong internal-priming site in a well-expressed gene is
exactly the thing that becomes dominant.  It exists to rescue models the
evidence tiers cannot speak to, not to adjudicate the ones they can.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import Config, SecondPassParams
from .genome import Genome, revcomp
from .intervals import Exon, total_overlap
from .model import TranscriptModel
from .molecules import NO_CELL, NO_UMI, n_cells, n_molecules
from .reads import ContigReads
from .reference import ReferenceIndex

#: What a mono-exon model's gene context makes it.  Order matters for reporting
#: only; the verdict comes from ``SecondPassParams.remove_classes``.
MONOEXON_CLASSES = ("A", "B", "C", "D", "E", "X")

#: Tiers of the end3 gate, strongest first.  Also the report order.
END3_TIERS = (
    "atlas",
    "nontemplated",
    "composition",
    "short_read",
    "calibrated",
    "gene_dominant",
)


@dataclass
class SecondPassResult:
    """Everything the parent needs to report and audit one contig's second pass."""

    models: List[TranscriptModel] = field(default_factory=list)
    audit: List[Dict[str, object]] = field(default_factory=list)
    stats: Dict[str, int] = field(default_factory=dict)


# ====================================================================== #
# terminal evidence
# ====================================================================== #
def three_prime_position(m: TranscriptModel) -> int:
    """The model's 3' coordinate, in the convention ``polya_evidence`` uses."""
    return int(m.exons[-1][1]) if m.strand == "+" else int(m.exons[0][0])


def _downstream(genome: Genome, contig: str, pos: int, strand: str, n: int) -> str:
    """``n`` genomic bases immediately 3' of ``pos``, in transcript orientation."""
    if strand == "+":
        return genome.fetch(contig, pos, pos + n)
    return revcomp(genome.fetch(contig, max(0, pos - n), pos))


def _upstream(genome: Genome, contig: str, pos: int, strand: str, n: int) -> str:
    """``n`` genomic bases immediately 5' of ``pos``, in transcript orientation."""
    if strand == "+":
        return genome.fetch(contig, max(0, pos - n), pos)
    return revcomp(genome.fetch(contig, pos, pos + n))


def _leading_a_run(s: str) -> int:
    n = 0
    for ch in s:
        if ch != "A":
            break
        n += 1
    return n


def _frac(s: str, bases: str) -> float:
    if not s:
        return float("nan")
    return sum(s.count(b) for b in bases) / len(s)


class _ARunCache:
    """Genomic A-run length at a 3' coordinate, memoised.

    The per-read test needs the run length at each *read's own* 3' end, not at
    the model's, because a peak is 4-35 bp wide and the genomic context can
    change inside it.  A peak has at most a few dozen distinct end positions,
    so memoising on the position makes the exact test as cheap as the
    approximate one.
    """

    def __init__(self, genome: Genome, contig: str, scan: int) -> None:
        self._g = genome
        self._contig = contig
        self._scan = scan
        self._cache: Dict[Tuple[int, str], int] = {}

    def at(self, pos: int, strand: str) -> int:
        key = (int(pos), strand)
        v = self._cache.get(key)
        if v is None:
            v = _leading_a_run(
                _downstream(self._g, self._contig, int(pos), strand, self._scan)
            )
            self._cache[key] = v
        return v


def measure_terminal_evidence(
    models: Sequence[TranscriptModel],
    contig: str,
    genome: Optional[Genome],
    reads: Optional[ContigReads],
    params: SecondPassParams,
    measure_tails: bool = True,
) -> None:
    """Fill in the 3'-end evidence the second pass and the terminal score read.

    Called after model building and *before* ``score_terminals``, so the
    calibrated model sees these features too rather than the second pass being
    the only consumer.  Writes only into ``m.evidence``; it decides nothing.

    Every key is always written, blank where it could not be measured, so that
    "no genome supplied" and "measured and found nothing" stay distinguishable
    in ``models.tsv`` -- the same rule the CAGE and coverage channels follow.
    """
    blank = {
        "genomic_a_run_3p": "",
        "median_tail_excess": "",
        "frac_reads_nontemplated": "",
        "n_nontemplated_reads": "",
        "dse_gu_frac": "",
        "dse_a_frac": "",
        "dse_gu_minus_a": "",
        "upstream_u_frac": "",
        "max_a_run_downstream": "",
    }
    if genome is None:
        for m in models:
            for k, v in blank.items():
                m.evidence.setdefault(k, v)
        return

    cache = _ARunCache(genome, contig, params.a_run_scan)
    have_tails = bool(
        measure_tails
        and reads is not None
        and getattr(reads, "tail_len", None) is not None
        and reads.tail_len.size
    )

    for m in models:
        for k, v in blank.items():
            m.evidence.setdefault(k, v)
        pos = three_prime_position(m)
        down = _downstream(genome, contig, pos, m.strand, params.a_run_scan)
        a_run = _leading_a_run(down)
        m.evidence["genomic_a_run_3p"] = int(a_run)
        # longest A-run anywhere in the first `a_run_scan` bases, which is what
        # catches a priming site a few bases downstream of where the aligner
        # happened to stop
        best, cur = 0, 0
        for ch in down:
            cur = cur + 1 if ch == "A" else 0
            best = max(best, cur)
        m.evidence["max_a_run_downstream"] = int(best)

        if params.use_composition:
            dse = _downstream(
                genome, contig, pos, m.strand, params.dse_end
            )[params.dse_start:]
            gu, a = _frac(dse, "GT"), _frac(dse, "A")
            m.evidence["dse_gu_frac"] = _round(gu)
            m.evidence["dse_a_frac"] = _round(a)
            m.evidence["dse_gu_minus_a"] = _round(
                gu - a if (gu == gu and a == a) else float("nan")
            )
            up = _upstream(genome, contig, pos, m.strand, params.upstream_u_window)
            m.evidence["upstream_u_frac"] = _round(_frac(up, "T"))

        if have_tails and m.reads.size:
            sel = m.reads
            tl = reads.tail_len[sel].astype(np.int64)
            # the run length the GENOME already accounts for, at each read's own
            # end -- this is the whole test
            gen = np.fromiter(
                (cache.at(int(p), m.strand) for p in reads.tts[sel]),
                dtype=np.int64,
                count=sel.size,
            )
            excess = tl - gen
            ok = excess >= params.min_tail_excess
            m.evidence["median_tail_excess"] = int(np.median(excess))
            m.evidence["n_nontemplated_reads"] = int(ok.sum())
            m.evidence["frac_reads_nontemplated"] = _round(
                float(ok.sum()) / float(sel.size)
            )


def _round(v: float) -> object:
    if v is None or v != v:
        return ""
    return round(float(v), 4)


# ====================================================================== #
# rule 1 -- mono-exon gene-context classes
# ====================================================================== #
class _GeneStructure:
    """Per-gene exon/intron facts the A-E classifier needs, memoised.

    Built lazily because a contig touches a small fraction of its genes and the
    constitutive-intron computation is not free.
    """

    def __init__(self, ref: ReferenceIndex) -> None:
        self._ref = ref
        self._cache: Dict[str, Tuple[int, int, List[Exon], List[Exon]]] = {}

    def get(self, gene_id: str) -> Tuple[int, int, List[Exon], List[Exon]]:
        """``(n_mono_tx, n_multi_tx, distinct exons, constitutive introns)``."""
        v = self._cache.get(gene_id)
        if v is not None:
            return v
        g = self._ref.genes[gene_id]
        n_mono = n_multi = 0
        exons: set = set()
        for tid in g.tx_ids:
            t = self._ref.tx[tid]
            if t.n_exons <= 1:
                n_mono += 1
            else:
                n_multi += 1
            exons.update((int(s), int(e)) for s, e in t.exons)
        # Constitutive introns: inside the gene span, exonic in no isoform.
        introns: List[Exon] = []
        prev = g.start
        for s, e in g.exon_union:
            if s > prev:
                introns.append((prev, s))
            prev = max(prev, e)
        if g.end > prev:
            introns.append((prev, g.end))
        v = (n_mono, n_multi, sorted(exons), introns)
        self._cache[gene_id] = v
        return v


def classify_monoexon(
    m: TranscriptModel,
    ref: ReferenceIndex,
    struct: _GeneStructure,
    params: SecondPassParams,
) -> str:
    """A-E (or X) for one single-exon model.  Annotation only; no read evidence."""
    if m.n_exons != 1:
        return ""
    gid = m.gene_id
    if not gid or gid not in ref.genes:
        return "X"
    n_mono, n_multi, exons, introns = struct.get(gid)
    if n_multi == 0:
        return "A"
    if n_mono > 0:
        return "B"
    exon: Exon = (int(m.exons[0][0]), int(m.exons[0][1]))
    # D before C: a model contained in a single annotated exon cannot, by
    # construction, cover constitutively intronic sequence, so the two are
    # mutually exclusive -- but checking containment first makes that explicit
    # rather than a consequence the reader has to derive.
    for s, e in exons:
        if s <= exon[0] and exon[1] <= e:
            return "D"
    # "Spans an intron" means exonic AND intronic sequence, not merely intronic:
    # a pile lying wholly inside an intron has not spanned anything, it is an
    # unspliced precursor or a mis-assignment, and calling that C would put two
    # quite different failures under one label. It is E.
    intronic = total_overlap([exon], introns) if introns else 0
    exonic = total_overlap([exon], self_exon_union(ref, gid))
    if intronic >= params.min_intron_overlap and exonic >= params.min_intron_overlap:
        return "C"
    return "E"


def self_exon_union(ref: ReferenceIndex, gene_id: str) -> List[Exon]:
    return ref.genes[gene_id].exon_union


def _support_of(m: TranscriptModel) -> int:
    return m.n_mols if m.n_mols > 0 else m.n_reads


def _recount(
    m: TranscriptModel, reads: Optional[ContigReads], umi_hamming: int
) -> None:
    """Recompute a model's support after reads were folded onto it."""
    m.n_reads = int(m.reads.size)
    if reads is None or m.reads.size == 0:
        return
    pairs = [(int(reads.cell[j]), int(reads.umi[j])) for j in m.reads]
    if any(c != int(NO_CELL) and u != int(NO_UMI) for c, u in pairs):
        m.n_mols = n_molecules(pairs, umi_hamming)
        m.n_cells = n_cells(pairs)


def _inferred_model(
    tid: str,
    ref: ReferenceIndex,
    contig: str,
    params: SecondPassParams,
) -> TranscriptModel:
    """A model standing in for a reference transcript nothing observed intact.

    Structure is the reference transcript's, verbatim -- this stage does not
    invent coordinates.  Three independent markers say it was inferred
    (category, flag, evidence column) because the single failure mode that
    matters here is an inferred row being read as an observation.
    """
    t = ref.tx[tid]
    m = TranscriptModel(
        model_id="",
        contig=contig,
        strand=t.strand,
        exons=[(int(s), int(e)) for s, e in t.exons],
        chain=t.chain,
        category=params.inferred_category,
        reads=np.zeros(0, np.int64),
        n_reads=0,
        n_mols=0,
        n_cells=0,
        parent_tx=tid,
        gene_id=t.gene_id,
        gene_name=t.gene_name,
        flags=["INFERRED", "HYPOTHETICAL", "no_observed_full_structure"],
    )
    m.associated_tx = tid
    m.evidence["inferred"] = 1
    m.evidence["inferred_source"] = "reference"
    m.evidence["inferred_transcript"] = tid
    m.evidence["end_is_annotated"] = 1
    m.evidence["three_prime_dispersion"] = ""
    return m


def apply_monoexon_classes(
    models: List[TranscriptModel],
    contig: str,
    ref: ReferenceIndex,
    reads: Optional[ContigReads],
    params: SecondPassParams,
    umi_hamming: int = 1,
    audit: Optional[List[Dict[str, object]]] = None,
) -> Tuple[List[TranscriptModel], Dict[str, int]]:
    """Rule 1.  Returns the surviving catalogue and what it did."""
    stats: Dict[str, int] = defaultdict(int)
    if not params.monoexon_classes:
        return models, dict(stats)

    struct = _GeneStructure(ref)
    remove = set(params.remove_classes)

    # The replacement target is chosen from the catalogue BEFORE anything is
    # removed, so a gene's verdict does not depend on the order its mono-exon
    # models happen to be visited in.
    best_spliced: Dict[str, TranscriptModel] = {}
    for m in models:
        if m.n_exons > 1 and m.gene_id:
            cur = best_spliced.get(m.gene_id)
            if cur is None or _support_of(m) > _support_of(cur):
                best_spliced[m.gene_id] = m

    keep: List[TranscriptModel] = []
    new_inferred: Dict[Tuple[str, str], TranscriptModel] = {}
    folded: List[Tuple[TranscriptModel, np.ndarray]] = []

    for m in models:
        cls = classify_monoexon(m, ref, struct, params)
        if not cls:
            keep.append(m)
            continue
        m.evidence["monoexon_gene_class"] = cls
        stats[f"monoexon_class_{cls}"] += 1
        stats[f"monoexon_class_{cls}_reads"] += int(m.n_reads)

        if cls == "X":
            if params.intergenic_policy == "remove":
                stats["monoexon_intergenic_removed"] += 1
                stats["monoexon_intergenic_removed_reads"] += int(m.n_reads)
                _audit(audit, m, cls, "removed", "intergenic_policy=remove", "")
                continue
            keep.append(m)
            _audit(audit, m, cls, "kept", "intergenic_out_of_scope", "")
            continue

        if cls not in remove:
            keep.append(m)
            _audit(audit, m, cls, "kept", "gene_admits_monoexonic_isoforms", "")
            continue

        # -- replace ------------------------------------------------------
        target: Optional[TranscriptModel] = None
        source = ""
        tid = ""
        if params.infer_from == "data_then_reference":
            target = best_spliced.get(m.gene_id or "")
            if target is not None:
                source = "data"
                tid = target.parent_tx or target.associated_tx or target.model_id
        if target is None:
            g = ref.genes.get(m.gene_id or "")
            # tx_ids is pre-sorted by (canonical rank, -length), so [0] is
            # MANE_Select where one exists, then Ensembl_canonical, then basic,
            # then the longest of whatever is left.
            tid = g.tx_ids[0] if (g and g.tx_ids) else ""
            if not tid:
                # A gene with no usable transcript cannot be inferred from. The
                # model is removed and the reads are counted as unassigned,
                # because inventing a structure here is the one thing this
                # stage must never do.
                stats["monoexon_removed_no_inference_target"] += 1
                stats["monoexon_removed_no_target_reads"] += int(m.n_reads)
                _audit(audit, m, cls, "removed", "no_reference_transcript", "")
                continue
            key = (m.gene_id or "", tid)
            target = new_inferred.get(key)
            if target is None:
                target = _inferred_model(tid, ref, contig, params)
                new_inferred[key] = target
            source = "reference"

        if params.fold_reads_into_inferred and m.reads.size:
            folded.append((target, m.reads))
        target.evidence["inferred_from_monoexon"] = int(
            target.evidence.get("inferred_from_monoexon", 0) or 0
        ) + 1
        target.evidence.setdefault("inferred_transcript", tid)
        target.flags = sorted(set(target.flags) | {"absorbed_monoexon_class_" + cls})
        stats["monoexon_replaced"] += 1
        stats[f"monoexon_replaced_{source}"] += 1
        stats["monoexon_replaced_reads"] += int(m.n_reads)
        _audit(audit, m, cls, "replaced", f"infer_from={source}", tid)

    for target, idx in folded:
        target.reads = np.concatenate([target.reads, idx])
    for target in {id(t): t for t, _ in folded}.values():
        _recount(target, reads, umi_hamming)

    for key in sorted(new_inferred):
        m = new_inferred[key]
        # An inferred model with no reads at all is a claim with nothing behind
        # it; that only happens with fold_reads_into_inferred off, and then it
        # is the point.
        keep.append(m)
        stats["inferred_models_created"] += 1

    return keep, dict(stats)


def _audit(
    audit: Optional[List[Dict[str, object]]],
    m: TranscriptModel,
    cls: str,
    verdict: str,
    reason: str,
    target: str,
) -> None:
    if audit is None:
        return
    audit.append({
        # Model ids do not exist yet -- they are assigned in the parent, after
        # the second pass has settled which models there are. `_obj` is patched
        # into `model_index` by run_second_pass and the id filled in later, so
        # a surviving row can be joined back to models.tsv; a removed row keeps
        # -1, because a removed model never gets an id to join on.
        "_obj": m,
        "rule": "monoexon_class",
        "model_index": -1,
        "model_id": "",
        "contig": m.contig,
        "strand": m.strand,
        "start": int(m.start) + 1,
        "end": int(m.end),
        "category": m.category,
        "gene_id": m.gene_id or "",
        "gene_name": m.gene_name,
        "n_reads": int(m.n_reads),
        "n_molecules": int(m.n_mols),
        "monoexon_gene_class": cls,
        "verdict": verdict,
        "reason": reason,
        "target": target,
        "tier": "",
    })


# ====================================================================== #
# rule 2 -- the end3_novel gate
# ====================================================================== #
def _truthy(v: object) -> bool:
    """`evidence` mixes "", None, 0/1 and bools; only a positive counts."""
    if v is None or v == "":
        return False
    if isinstance(v, bool):
        return v
    try:
        return float(v) > 0
    except (TypeError, ValueError):
        return False


def _num(v: object, default: float = float("nan")) -> float:
    if v is None or v == "":
        return default
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return f


def end3_tier(
    m: TranscriptModel,
    params: SecondPassParams,
    is_sole: bool,
    is_dominant: bool,
) -> Optional[str]:
    """Which tier of evidence, if any, lets this 3' end survive.

    First hit wins, strongest first, so the recorded tier is unambiguous and a
    downstream filter can reproduce the decision from ``models.tsv`` alone.
    """
    ev = m.evidence

    if params.use_atlas and _truthy(ev.get("end_at_polya_site")):
        return "atlas"

    if params.use_templated_tail:
        frac = _num(ev.get("frac_reads_nontemplated"))
        excess = _num(ev.get("median_tail_excess"))
        if (frac == frac and frac >= params.min_nontemplated_frac
                and excess == excess and excess >= params.min_tail_excess
                and not _truthy(ev.get("internal_priming"))):
            return "nontemplated"

    if params.use_composition:
        # Positive evidence (a hexamer) AND a peak tight enough to be a cleavage
        # event AND no internal-priming signature downstream. All three: a
        # hexamer is six bases and 3'UTRs supply one by chance, a tight peak on
        # its own is also what a reproducible RT drop-off looks like, and
        # A-rich downstream is the artefact this whole tier exists to exclude.
        disp = _num(ev.get("three_prime_dispersion"), float("inf"))
        contrast = _num(ev.get("dse_gu_minus_a"))
        if (_truthy(ev.get("polya_motif_found"))
                and contrast == contrast
                and contrast >= params.min_dse_gu_minus_a
                and disp <= params.max_dispersion
                and not _truthy(ev.get("internal_priming"))):
            return "composition"

    if params.use_short_read_step and ev.get("sr_3p_step_supported") is True:
        return "short_read"

    if params.max_terminal_lfdr is not None:
        lf = _num(ev.get("terminal_lfdr"))
        if lf == lf and lf <= params.max_terminal_lfdr:
            return "calibrated"

    if params.gene_dominance_backstop and (is_sole or is_dominant):
        return "gene_dominant"

    return None


def apply_end3_filter(
    models: List[TranscriptModel],
    reads: Optional[ContigReads],
    params: SecondPassParams,
    umi_hamming: int = 1,
    audit: Optional[List[Dict[str, object]]] = None,
) -> Tuple[List[TranscriptModel], Dict[str, int]]:
    """Rule 2.  Returns the surviving catalogue and what it did."""
    stats: Dict[str, int] = defaultdict(int)
    if not params.end3_filter:
        return models, dict(stats)

    judged = set(params.end3_categories)
    by_gene: Dict[str, List[TranscriptModel]] = defaultdict(list)
    for m in models:
        if m.gene_id:
            by_gene[m.gene_id].append(m)

    # Gene dominance is read off the catalogue as it stands BEFORE any removal,
    # so a model is never promoted to "dominant" by the removal of a rival in
    # the same pass -- which would make the result depend on iteration order.
    #
    # A model with no gene gets neither flag. "Sole model of its gene" is a
    # claim about competition within a locus; a model that belongs to no locus
    # has not won that competition, it was never in it, and defaulting such a
    # model to `sole=True` would hand every intergenic novel 3' end a free pass
    # through the backstop.
    sole: Dict[int, bool] = defaultdict(bool)
    dominant: Dict[int, bool] = defaultdict(bool)
    for gid, group in by_gene.items():
        sup = sorted((_support_of(x) for x in group), reverse=True)
        top = sup[0] if sup else 0
        runner = sup[1] if len(sup) > 1 else 0
        for x in group:
            sole[id(x)] = len(group) == 1
            dominant[id(x)] = (
                _support_of(x) == top
                and top > 0
                and (runner == 0 or top >= params.dominance_ratio * runner)
            )

    keep: List[TranscriptModel] = []
    dropped: List[TranscriptModel] = []
    for m in models:
        if m.category not in judged:
            keep.append(m)
            continue
        tier = end3_tier(m, params, sole[id(m)], dominant[id(m)])
        m.evidence["end3_second_pass_tier"] = tier or "none"
        if tier is None:
            dropped.append(m)
            stats["end3_dropped"] += 1
            stats["end3_dropped_reads"] += int(m.n_reads)
            _audit_end3(audit, m, "removed", "no_surviving_tier", "", "")
            continue
        m.flags = sorted(set(m.flags) | {f"end3_kept_by_{tier}"})
        stats["end3_kept"] += 1
        stats[f"end3_kept_{tier}"] += 1
        stats["end3_kept_reads"] += int(m.n_reads)
        keep.append(m)
        _audit_end3(audit, m, "kept", f"tier={tier}", tier, "")

    if not dropped:
        return keep, dict(stats)

    # -- re-home the reads of every dropped model -------------------------
    # Same principle the 5' path has always used: rejecting a 3' peak is a
    # statement about which model owns this read mass, not a licence to lose it.
    if params.fold_dropped_reads:
        by_chain: Dict[Tuple[str, object], List[TranscriptModel]] = defaultdict(list)
        surviving_by_gene: Dict[str, List[TranscriptModel]] = defaultdict(list)
        for x in keep:
            by_chain[(x.strand, x.chain)].append(x)
            if x.gene_id:
                surviving_by_gene[x.gene_id].append(x)
        folds: List[Tuple[TranscriptModel, np.ndarray]] = []
        for m in dropped:
            sibs = by_chain.get((m.strand, m.chain)) or []
            if not sibs and m.gene_id:
                sibs = surviving_by_gene.get(m.gene_id) or []
            if not sibs or not m.reads.size:
                stats["end3_dropped_reads_unassigned"] += int(m.n_reads)
                continue
            tgt = max(sibs, key=_support_of)
            folds.append((tgt, m.reads))
            tgt.flags = sorted(set(tgt.flags) | {"absorbed_end3_novel"})
            tgt.evidence["n_3p_absorbed_reads"] = int(
                tgt.evidence.get("n_3p_absorbed_reads", 0) or 0
            ) + int(m.n_reads)
            stats["end3_dropped_reads_folded"] += int(m.n_reads)
        for tgt, idx in folds:
            tgt.reads = np.concatenate([tgt.reads, idx])
        for tgt in {id(t): t for t, _ in folds}.values():
            _recount(tgt, reads, umi_hamming)
    else:
        for m in dropped:
            stats["end3_dropped_reads_unassigned"] += int(m.n_reads)

    return keep, dict(stats)


def _audit_end3(
    audit: Optional[List[Dict[str, object]]],
    m: TranscriptModel,
    verdict: str,
    reason: str,
    tier: str,
    target: str,
) -> None:
    if audit is None:
        return
    audit.append({
        "_obj": m,
        "rule": "end3_novel",
        "model_index": -1,
        "model_id": "",
        "contig": m.contig,
        "strand": m.strand,
        "start": int(m.start) + 1,
        "end": int(m.end),
        "category": m.category,
        "gene_id": m.gene_id or "",
        "gene_name": m.gene_name,
        "n_reads": int(m.n_reads),
        "n_molecules": int(m.n_mols),
        "monoexon_gene_class": "",
        "verdict": verdict,
        "reason": reason,
        "target": target,
        "tier": tier,
    })


# ====================================================================== #
def run_second_pass(
    models: List[TranscriptModel],
    contig: str,
    ref: ReferenceIndex,
    reads: Optional[ContigReads],
    cfg: Config,
) -> SecondPassResult:
    """Both rules, in order, on one contig's catalogue.

    Mono-exon classes run first: they can create an inferred multi-exon model,
    and that model is then a legitimate sibling for the end3 rule to fold a
    dropped peak into.  The reverse order would make the two rules interact
    through whichever happened to run first.
    """
    params = cfg.secondpass
    res = SecondPassResult(models=models)
    if not params.enabled:
        return res

    audit: Optional[List[Dict[str, object]]] = [] if params.write_audit else None
    stats: Dict[str, int] = {}

    models, s1 = apply_monoexon_classes(
        models, contig, ref, reads, params, cfg.molecules.umi_hamming, audit
    )
    stats.update(s1)
    models, s2 = apply_end3_filter(
        models, reads, params, cfg.molecules.umi_hamming, audit
    )
    stats.update(s2)

    # Resolve the surviving rows to their position in the final catalogue, and
    # drop the object references so the result can cross a process boundary.
    if audit is not None:
        where = {id(m): i for i, m in enumerate(models)}
        for row in audit:
            obj = row.pop("_obj", None)
            if obj is not None:
                row["model_index"] = where.get(id(obj), -1)

    res.models = models
    res.audit = audit or []
    res.stats = stats
    return res


AUDIT_COLUMNS = (
    "rule", "model_id", "contig", "strand", "start", "end", "category",
    "gene_id", "gene_name", "n_reads", "n_molecules", "monoexon_gene_class",
    "verdict", "reason", "target", "tier",
)


def write_audit(rows: Sequence[Dict[str, object]], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(AUDIT_COLUMNS) + "\n")
        for r in rows:
            fh.write(
                "\t".join(str(r.get(c, "")) for c in AUDIT_COLUMNS) + "\n"
            )
