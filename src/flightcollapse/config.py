"""Configuration.

Every tunable lives here with a default and a one-line justification.  The CLI
can override any field from a JSON/YAML file or ``--set key=value``.

Design rule inherited from the isoseq post-mortem: **no parameter may be
silently inactive**.  :meth:`Config.validate` raises on combinations where a
setting would have no effect, rather than accepting and ignoring it.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, fields, asdict
from typing import Any, Dict, List, Optional


@dataclass
class AlignmentGates:
    """Read-level admission criteria.

    NOTE on ``min_aln_coverage``.  The isoseq run being replaced used 0.99, but
    on an **untrimmed** FLNC BAM the polyA tail (and any TSO remnant) is
    soft-clipped, so coverage measured against the full query length is ~0.93
    at the median.  A 0.99 gate there keeps ~6% of reads *and* enriches for
    short mono-exonic ones, which is the opposite of what you want.  The
    default is therefore 0.90, with the terminal clips gated separately and
    explicitly.  Set 0.99 when running on a polyA-trimmed BAM.
    """

    min_aln_coverage: float = 0.90      # aligned query bases / query length
    #: 0.4.1.  Take the transcript-3' terminal clip OUT of the coverage
    #: denominator.  The docstring above already said the terminal clips should
    #: be "gated separately and explicitly" -- but the 3' clip was still counted
    #: as unaligned query, so ``min_aln_coverage`` was partly a length filter:
    #: a near-constant ~99 bp clip against an aligned length spanning an order
    #: of magnitude.  On BD176c the 0.90 gate dropped 8,406,872 of 24,254,151
    #: scanned reads -- 88% of those 500-800 bp long, and 85% of the rejected
    #: reads MULTI-exon.  With the clip excluded the same threshold drops 0.3%.
    #: False reproduces 0.2.0-0.4.0 read admission exactly.
    coverage_excludes_3p_clip: bool = True
    min_aln_identity: float = 0.95      # 1 - NM / (aligned columns, introns excluded)
    #: terminal clips in TRANSCRIPT orientation; None disables.  The 3' clip is
    #: where an untrimmed polyA tail lives, so it is left open by default.
    max_5p_softclip: Optional[int] = 200
    max_3p_softclip: Optional[int] = None
    min_mapq: int = 0                   # minimap2 long-read MAPQ is not well calibrated
    min_intron_len: int = 20            # shorter N ops are alignment artefacts, not introns
    max_intron_len: int = 500_000
    drop_supplementary: bool = True
    drop_secondary: bool = True
    #: refuse to run silently when the gates discard more than this share of
    #: aligned reads -- a gate that throws away most of the data is a
    #: configuration error, not a filter
    max_gate_rejection_frac: float = 0.5
    #: ...but only once there are enough reads for that fraction to mean
    #: anything.  A contig with 715 alignments -- chrY in a female line, an
    #: unplaced scaffold -- is repeat-mapping junk, and rejecting 99% of it is
    #: the gates working, not a misconfiguration.
    min_reads_for_gate_check: int = 50_000
    #: Read the 3' soft-clipped SEQUENCE, not just its length, and record the
    #: terminal homopolymer run in transcript orientation.  This is the only
    #: direct observation of the polyA tail the tool has: everything else the
    #: 3' path calls "polyA evidence" is genomic (a PAS motif, downstream
    #: A-richness, an atlas), i.e. a statement about the locus rather than about
    #: the molecule that was primed there.  Costs one query_sequence decode per
    #: admitted read; set False to skip it.  New in 0.1.17, measurement only --
    #: see EndParams.polya_tail_gates_novel_end.
    measure_polya_tail: bool = True


@dataclass
class JunctionParams:
    fuzzy_tolerance: int = 5            # bp; census showed 0->20 bp changes chain count by 15%
    snap_to_annotation_first: bool = True
    #: absolute floors, applied on molecules (UMIs) when available, else reads
    min_novel_junction_umis: int = 3
    min_novel_junction_reads: int = 5
    #: A novel junction must hold this share of the reads crossing its donor.
    #:
    #: 0.1.16: lowered from 0.01, which was the single largest source of lost
    #: recall. Masking 7,224 real GENCODE transcripts from BD144_HM and asking
    #: which of their junctions the caller failed to emit: this rule rejected
    #: 174 of 230 (76%), and the junctions it rejected carried MORE absolute
    #: evidence than the novel junctions it kept --
    #:
    #:                     rejected (real)   kept (novel)
    #:   min_site_share          0.007          0.095
    #:   n_reads                    51             13
    #:   n_cells                  49.5             12
    #:   lfdr                    0.000          0.000
    #:
    #: A minority isoform at a well-expressed gene has a low share by
    #: definition; that is what alternative splicing looks like, not what noise
    #: looks like. And the information is not lost by lowering this: the scorer
    #: already consumes ``logit_donor_share`` and ``logit_acceptor_share`` as
    #: features and rated these junctions confidently real. A hard cut here was
    #: overruling a calibrated estimate on one of its own inputs.
    #:
    #: A floor is kept rather than removed: genuine misalignment scatter around
    #: a dominant junction does sit near zero share, and the lFDR is a
    #: population-level control that a single pathological locus can slip past.
    min_novel_junction_share: float = 0.001
    #: local-FDR ceiling from the calibrated model (None -> use thresholds only)
    max_junction_lfdr: Optional[float] = 0.05
    #: motif requirement for novel junctions
    require_canonical_motif: bool = True
    canonical_motifs: tuple = ("GT-AG", "GC-AG", "AT-AC")
    #: Motif classes a **novel** junction may have.  Annotated junctions are
    #: unaffected -- a real GC-AG intron that GENCODE knows about stays.
    #:
    #: This used to include ``semi_canonical`` (GC-AG, AT-AC), and the short
    #: reads showed that was wrong.  Measured on BD144 kin/hm/nd, of the novel
    #: junctions the caller kept:
    #:
    #:   GT-AG   42,118 / 42,605 / 39,014   short-read confirmed  5.8 / 5.8 / 6.2 %
    #:   GC-AG    4,562 (+137 AT-AC)        short-read confirmed  0.13 / 0.12 / 0.15 %
    #:
    #: while *annotated* semi-canonical junctions confirm at 33.4%, so the short
    #: reads see this motif class perfectly well -- it is not a detection bias.
    #: Nor are they shifted copies of real junctions: only 12% have an annotated
    #: splice site within 5 bp, against 51% for novel GT-AG.  Dropping the whole
    #: class removes 4,699 junctions and costs 6 confirmed ones.
    novel_motif_classes: tuple = ("canonical", "unknown")
    #: 0.1.16. The evidence above is about junctions that are novel *and* absent
    #: from the annotation -- overwhelmingly artefacts, and the rule is right
    #: about them. Masking says something different about the other population:
    #: junctions that are real but that the annotation happens not to know.
    #: Of 525 real junctions the caller failed to emit, 44 GC-AG and 6 AT-AC
    #: were rejected by this rule, 8-10% of the total, and 27 of the 44 were
    #: short-read confirmed in >=2 replicates.
    #:
    #: The case that settles the design: SOD1's canonical transcript. 10,419
    #: long reads on the chain, 1,005 unique short reads in 3 of 3 replicates,
    #: both splice sites annotated -- rejected for being GC-AG.
    #:
    #: So the class is not admitted outright (that would re-admit the 4,699
    #: artefacts per sample) but conditionally, on independent evidence. Set
    #: ``conditional_motif_classes = ()`` to restore 0.1.15 behaviour exactly.
    conditional_motif_classes: tuple = ("semi_canonical",)
    #: short-read replicates that must carry a conditional-motif junction...
    conditional_motif_min_sr_replicates: int = 2
    #: ...or, with no short reads available, the long-read floor it must clear.
    #: Deliberately high: without external evidence this is the only thing
    #: standing between a GC-AG call and the artefact population above.
    conditional_motif_min_reads: int = 50
    #: Additionally require a novel junction to reuse at least one annotated
    #: splice site.  Off by default: it is a real trade-off rather than a free
    #: one.  Novel GT-AG junctions with two invented sites confirm at 0.08-0.20
    #: of the annotated rate versus 0.26-0.35 for those reusing a site, so the
    #: class is much weaker -- but it is also the only place a genuinely new
    #: splice site can appear.  Turning this on removes a further 9,160
    #: junctions and costs 182 confirmed ones.
    require_annotated_site_for_novel_junction: bool = False
    #: RT template-switching: reject if the direct repeat at the intron boundaries
    #: is at least this long (SQANTI-style)
    rt_switch_repeat_len: int = 8
    #: Reads carrying a rejected junction are diverted into the decoy set, which
    #: is the negative anchor the chain-level novelty model trains on. Relaxing
    #: the gates above shrinks that set as a side effect, and an empty decoy set
    #: makes the chain calibration fall back to flat thresholds *silently*.
    #: Below this many decoy chain groups genome-wide, say so loudly.
    #:
    #: The bulk of decoys come from non-canonical and antisense motifs, neither
    #: of which 0.1.16 relaxes, so this is a tripwire rather than an expectation.
    min_decoy_chain_groups: int = 2_000
    #: An identical intron chain fixes every *internal* exon boundary, but not
    #: the exon contents: a short intron spelled as a CIGAR ``D`` rather than an
    #: ``N`` leaves the chain untouched, so a real splicing event can hide
    #: inside an "FSM".  Promote those to junctions before the census, then let
    #: ordinary curation judge them.
    promote_internal_deletions: bool = True
    #: recurrence needed to promote a gap that is not already annotated
    min_internal_deletion_reads: int = 3
    #: above this, a deletion is a structural variant or a misalignment, not an
    #: unrecognised intron
    max_internal_deletion_len: int = 5_000

    # -- short-read support (STAR SJ.out.tab) ---------------------------- #
    #: uniquely-mapping short reads needed to call a junction short-read
    #: supported.  Multi-mapping reads are recorded but not counted here.
    min_sr_unique_reads: int = 2
    #: STAR's maximum spliced overhang; a junction seen only through 5 bp of
    #: flanking sequence is an alignment guess
    min_sr_overhang: int = 10
    #: Replicates a junction must appear in, clamped to the number of SJ files
    #: supplied (so a single library is unaffected).
    #:
    #: Summing replicates is right for evidence but destroys the information
    #: that matters most at low counts.  Measured on BD144 a/b/c: of the 10,007
    #: novel junctions reaching 2 summed unique reads, 4,302 appear in only one
    #: replicate where sampling predicts 1,812 -- an excess of ~2,490 that is
    #: concentrated entirely at 2-5 reads and gone by 10.  An anchor set wants
    #: precision far more than recall, so requiring two replicates is nearly
    #: free: it costs ~1,800 genuine junctions out of ~7,500 and removes ~2,490
    #: library-specific ones.
    min_sr_replicates: int = 2
    #: Use short-read-supported NOVEL junctions as positive anchors for the
    #: calibration.  This is the point of having short reads: the positive set
    #: was previously "annotated", so the model could only learn to recognise
    #: things already in GENCODE.  Supported novel junctions are true positives
    #: that annotation cannot supply.
    use_short_reads_as_anchor: bool = True
    #: Fraction of short-read-supported novel junctions withheld from the
    #: anchor set, so that recall measured on them is not the model grading its
    #: own work.  Split deterministically on position, so the same junction is
    #: held out in every sample.
    sr_holdout_frac: float = 0.5
    #: Exempt a short-read-supported junction from curation entirely.  This one
    #: changes which models are called, so it is opt-in.  Absence of short-read
    #: support is NEVER a reason to drop a junction, in either setting.
    short_read_rescue: bool = False


@dataclass
class ChainParams:
    #: Candidate chains below this are not emitted as their own model.
    #:
    #: These are deliberately *low*.  Since 0.1.7 the tool emits at a low floor
    #: and carries the calibrated local FDR, the category and every flag on each
    #: row, so the threshold is chosen downstream on the finished table instead
    #: of being baked into a 3.5-hour run.  Nothing is discarded silently: a
    #: chain below the floor is re-homed onto a surviving chain it is a suffix
    #: of, or counted in ``reads_dropped``.
    min_chain_umis: int = 3
    min_chain_reads: int = 5
    #: chain must hold this share of its gene's molecules to be emitted as novel
    min_novel_chain_share: float = 0.01
    #: Advisory by default (None).  Unlike the junction level, the chain level
    #: has no clean null -- there is no chain-shaped analogue of a wrong-strand
    #: splice motif -- so the decoy set has to lean on low support and on
    #: "carries a junction that failed curation".  The score is written to
    #: models.tsv for every model; set a value here (0.05 is a reasonable start)
    #: once you have looked at that distribution for your data.
    max_chain_lfdr: Optional[float] = None
    #: direction of the suffix relation, by support ratio  reads(S)/reads(top parent)
    debris_ratio: float = 0.33          # S << P  -> S is 5'-truncation debris, merge up
    dominant_ratio: float = 3.0         # S >> P  -> S is real, P is a rare 5'-extension
    #: annotation overrides the ratio in the undecidable band
    annotation_breaks_ties: bool = True
    #: A non-annotated chain that is a suffix of an ANNOTATED chain is
    #: 5'-truncation of that transcript unless it has independent evidence of
    #: its own start site.  Support ratio does not rescue it: 5' degradation is
    #: pervasive here (57.8% of spliced read mass), so a truncated form
    #: routinely outnumbers the intact one.  This is the core reference-anchored
    #: rule -- with it off, the tool falls back to support direction alone.
    merge_unannotated_suffix_into_annotated_parent: bool = True
    #: ...and the escape hatch: keep it as its own model if its 5' end sits at an
    #: annotated TSS of the same gene (a real alternative start), within
    #: ``EndParams.max_5p_diff``.
    require_tss_evidence_to_keep_unannotated_suffix: bool = True
    #: keep ambiguous chains as their own flagged models rather than guessing
    keep_ambiguous_as_models: bool = True


@dataclass
class EndParams:
    #: Cleavage is heterogeneous over roughly 20-30 bp, so this is a biological
    #: constant, not a fudge factor: two 3' ends further apart than this are
    #: different cleavage sites.  It decides **naming**: FSM vs an alternative
    #: end.
    max_3p_diff: int = 30
    #: Width of the *peak-clustering* window. ``None`` means "use max_3p_diff",
    #: which is the measured-best setting; the knob exists because clustering
    #: ("is this one cleavage site?") and naming ("is this the annotated end?")
    #: are different questions and should be separable.
    #:
    #: History, because the reasoning here was wrong once. 0.1.13 set this to
    #: 100 on the theory that the narrower window introduced in 0.1.7 was
    #: fragmenting cleavage sites: three-way reproducibility at >=20 molecules
    #: appeared to fall from 88.0% to 73.5%, and pairs of models sharing a
    #: parent transcript with 3' ends 31-100 bp apart rose from 160 to 3,639.
    #:
    #: That measurement was made with a broken matcher (fixed in 0.1.14 -- it
    #: forced any sample appearing twice in a group into its own singleton).
    #: Re-measured properly, the narrow window is **better at every support
    #: threshold from 5 molecules up**: at >=10, 86.9% vs 78.3%; at >=20, 94.4%
    #: vs 91.5% on *fewer* structures (64,959 vs 67,178). It consolidates rather
    #: than fragments. The extra 3' peaks it resolves are real: for structures
    #: called in all three libraries the median 3'-end spread is 0 bp, including
    #: for ``end3_novel`` models, which are not snapped to annotation at all.
    #: Cleavage in this data really is that precise.
    peak_window: Optional[int] = None
    #: The 5' tolerance is doing a completely different job -- it says how much
    #: degradation to forgive before a start counts as different -- so it has no
    #: business sharing a value with the 3' one.
    max_5p_diff: int = 100
    #: 5' end of a model = this percentile of member reads' 5' ends, in
    #: transcript orientation.  Never the minimum -- that is the isoseq bug.
    five_prime_percentile: float = 90.0
    three_prime_statistic: str = "mode"  # "mode" | "median"
    #: support needed to split a chain group into a separate 3'-end model
    min_end_umis: int = 3
    min_end_reads: int = 5
    #: ...and, for an UNANNOTATED 3' end, this share of the chain group's
    #: molecules.  The absolute floor alone is far too permissive at deep loci:
    #: 10 reads out of 200,000 is noise, 10 out of 60 is the dominant site.
    min_end_share: float = 0.05
    #: emit a 3'UTR variant only with positive polyA evidence
    require_polya_evidence_for_utr_variant: bool = True
    #: A 3' peak that misses the parent transcript's annotated end may still sit
    #: on the annotated end of ANOTHER transcript of the same gene.  Checking
    #: only the parent badly undercounts: measured on this data, 35% of 3' ends
    #: fall within 10 bp of the associated transcript's end but 72% fall within
    #: 10 bp of the nearest annotated end of *any* transcript of the gene.
    #: Those get their own category rather than being called novel.
    snap_to_gene_level_3p_ends: bool = True
    #: a 5' *extension* past the annotated TSS is credible; a 5' *truncation* is
    #: usually degradation, so it needs an independent TSS signature.
    allow_5p_truncation_models: bool = False
    max_utr_variants_per_parent: int = 3

    # -- 0.1.17: the annotation fallback is bounded ------------------------
    #: How far a REJECTED 3' peak may be moved onto the default parent's
    #: annotated end.
    #:
    #: Through 0.1.16 this was unbounded.  When a peak failed the novel-end
    #: evidence rules its reads were re-labelled with ``dtx.tts`` whatever the
    #: distance, so a model could be extended past every read that supports it.
    #: NSL1 is the worked example: the read pile ends near chr1:212,738,000 and
    #: the emitted model ended at chr1:212,726,153, 11.6 kb away, because the
    #: empirical peak carried no polyA evidence.
    #:
    #: That inverts the chemistry the tool is built on.  Under oligo-dT the 3'
    #: end is the OBSERVED end; annotation may refine a read-supported cleavage
    #: coordinate, it must not manufacture a distant extension.  Beyond this
    #: distance the peak keeps its observed coordinate and is emitted as
    #: ``unresolved_3p_category`` instead.
    #:
    #: 300 bp is deliberately generous -- wide enough to absorb ordinary 3'UTR
    #: heterogeneity and short annotation disagreements, narrow enough that a
    #: kilobase-scale move can no longer happen silently.
    max_3p_fallback_dist: int = 300
    #: Category for a peak whose fallback was refused.  It is neither FSM (the
    #: annotated end is not where the molecules end) nor end3_novel (it failed
    #: the evidence rules for calling a new site), and collapsing it into either
    #: would state something the data does not support.
    unresolved_3p_category: str = "end3_unresolved"
    #: Emit those models (True) or drop them and count the reads (False).
    #:
    #: 0.3.0 flips this to False on benchmark evidence. `end3_unresolved` exists
    #: for a 3' peak whose annotation fallback was REFUSED -- the code has
    #: already concluded the evidence is insufficient, and then emitted a model
    #: anyway. On SIRV Set 4 the category was wrong 20/20 in the reference arm
    #: and 23/23 in the masked arm: not one correct call, in any arm. Dropping
    #: it costs zero recall and takes masked precision from 42.0% to 49.6%.
    #: The reads are still counted -- they fold into the surviving peak -- and
    #: the knob remains for anyone who wants the diagnostic back.
    emit_unresolved_3p: bool = False

    # -- 0.1.17: peaks are molecule-weighted -------------------------------
    #: Weight each read by 1/(reads sharing its cell+UMI) when discovering
    #: peaks, choosing the dominant peak, and placing the peak coordinate.
    #: Support GATES still count whole molecules and whole reads as before --
    #: only the mass that decides peak geometry changes.
    #:
    #: NOTE this is a no-op on an already-deduplicated BAM run with
    #: ``molecules.umi_hamming=0``, which is how the BD144 production runs and
    #: the masking experiment were done: there is one read per family there
    #: already.  It matters for runs over a raw FLNC BAM.
    weight_ends_by_molecule: bool = True

    # -- 0.1.17: direct polyA tail evidence, RECORDED not enforced ----------
    #: A read counts as tail-bearing when its 3' soft clip opens with at least
    #: this many A's (transcript orientation) and the clip is at least this
    #: A-rich over its first 60 bases.
    min_polya_tail_len: int = 8
    #: Extra A-richness check over the clip's first 60 bases. 0 disables it, and
    #: it is disabled by default: combined with the run requirement it demanded
    #: a tail nearly as long as the clip, which is wrong whenever the clip also
    #: carries adapter or TSO sequence. Kept as a knob, not a default.
    min_polya_tail_frac: float = 0.0
    #: Let that direct evidence satisfy ``require_polya_evidence_for_utr_variant``.
    #:
    #: ON since 0.1.18, on evidence rather than assumption. Measured on BD144_HM:
    #: 12,541 of 19,145 refused peaks (66%) failed on the GENOMIC polyA test, and
    #: the soft clips are not trimmed on this BAM, so the tail is both the
    #: binding constraint and available. Under oligo-dT the molecule's own tail
    #: outranks a hexamer at the locus: the hexamer describes the genome, the
    #: tail describes the molecule that was actually primed.
    #:
    #: The distinction it does NOT collapse is internal priming, which produces
    #: a read pile with no clipped tail -- it primes on genomic A's, which align.
    #: A tailed peak becomes ``end3_novel``; an untailed one stays
    #: ``unresolved_3p_category`` and is flagged. Set
    #: ``unresolved_3p_category=end3_novel`` to merge them anyway.
    polya_tail_gates_novel_end: bool = True
    #: ...and if it does gate, the share of a peak's molecules that must be
    #: tail-bearing.
    min_tail_molecule_frac: float = 0.30


@dataclass
class TerminalParams:
    """Competitive 3'-end pruning and the calibrated terminal score. New in 0.3.0.

    The problem this block exists for, measured on BD176c (159,068 models):
    34.2% of the catalogue is a terminal-3' variant, carrying 38% of read mass,
    and 31.4% of those sit on a parent transcript that already has an FSM model.
    One chain can carry up to twelve 3' models. The 3' end PLACEMENT is
    excellent -- 13.9 bp dispersion, zero annotation shift -- so this is not a
    coordinate problem. It is a multiplicity problem.

    The rule that fixes it is relative, not absolute. On SIRV Set 4, all 58
    right-chain-wrong-end models sat on a chain that already had an exact
    model, so absorbing every one of them costs no recall at all; whereas a
    flat read floor at 100 reads still keeps 13 false models and loses 2 true
    ones, because RT drop-off is reproducible and accumulates support exactly
    like a real cleavage site.

    This mirrors what the 5' path has always done: a 5'-truncated pile is
    absorbed into its parent and counted in ``n_5p_truncated_reads`` rather
    than emitted, "because otherwise every well-expressed gene grows a ladder
    of truncation models". The 3' path had no equivalent because under oligo-dT
    the 3' end is the trustworthy one -- and it is, for placement, but not for
    how many of them there are.
    """

    prune_terminal_variants: bool = True
    #: A shorter 3' variant is absorbed into a compatible sibling only when the
    #: sibling is at least this many times better supported. 1.0 would absorb
    #: on a tie; the default demands the sibling clearly dominates.
    sibling_support_ratio: float = 2.0
    #: ...and only when the variant sits at most this far inside the sibling's
    #: 3' end. Beyond it the two are different enough that "truncation of" stops
    #: being the simpler explanation.
    max_truncation_dist: int = 5_000
    #: Evidence that makes a shorter variant survive regardless of the ratio:
    #: its 3' end is an annotated end, or a catalogued polyA site, or the peak
    #: is tight enough to be a real cleavage site.
    keep_if_annotated_end: bool = True
    keep_if_atlas_site: bool = True
    #: Peak dispersion below which a 3' end is credible on its own.
    #:
    #: The single most discriminating number in the whole terminal problem, and
    #: it was computed and written to models.tsv all along without any decision
    #: reading it. BD176c medians: FSM 4.1 bp, end3_annotated 6.4, mono-exon
    #: 3.2 -- against end3_novel 34.6 and end3_unresolved 14.5. Cleavage is
    #: precise, so a real site is tight and a drop-off pile is diffuse. Read
    #: support (10 vs 11) and tail fraction (0.96 vs 1.00) separate nothing.
    max_dispersion_for_credible_end: float = 12.0

    # -- 0.4.0: promotion rather than absorption ---------------------------
    #: Flip the 3' default from "keep unless absorbed" to "do not emit unless
    #: promoted". This is the paradigm change, and the reason for it is that
    #: the opt-out form leaks: on SIRV Set 4 the 0.3.1 catalogue is 279 models
    #: against 55 observable truths, and only 29.7% of its 3' ends land within
    #: 100 bp of a real transcript end.
    #:
    #: A model with a NON-ANNOTATED 3' end is promoted when it has a catalogued
    #: polyA site, a short-read coverage step, or a peak tight enough to be a
    #: cleavage event. Everything else is absorbed into the most distal
    #: promoted sibling of its chain group -- its reads move, they are not
    #: dropped -- and a candidate with no sibling to absorb it is kept and
    #: tiered ``unsupported``, because reads with nowhere else to go are not
    #: evidence of nothing.
    #:
    #: Measured on SIRV with the atlas and coverage clauses inert, so this is
    #: the weakest form of the rule: 279 -> 103 models, 3'-end precision
    #: 29.7% -> 78.6%, two correct ends lost of 83. Masked arm 297 -> 115,
    #: 32.0% -> 80.0%, three lost of 95.
    promote_terminal_ends: bool = True
    #: Peak dispersion at or below which a non-annotated 3' end is promoted on
    #: its own. Distinct from ``max_dispersion_for_credible_end`` only so the
    #: promotion threshold can be tuned without changing 0.3.x absorption
    #: behaviour; both default to 12 bp.
    #:
    #: 12 bp is measured on SIRV, where cleavage is clean and coverage deep:
    #: median dispersion is 0.0-2.0 bp for correct 3' ends against 21-26 bp for
    #: wrong ones, AT EVERY READ DEPTH INCLUDING <= 5 READS. That last point is
    #: why this rule is not an abundance rule: what the reads agree on is the
    #: signal, not how many there are.
    max_dispersion_promote: float = 12.0

    #: Share of a chain group's CANDIDATE molecules above which a promoted
    #: novel 3' end is marked ``terminal_tier = high``.
    #:
    #: The competition is candidate-versus-candidate, never candidate-versus-
    #: FSM. Measured on SIRV: beating the annotated variant is 11.7%/11.8%
    #: precise -- barely better than keeping everything, because a genuine
    #: proximal polyA site usually carries a MINORITY of a gene's molecules, so
    #: the test selects against the biology. Ranking candidates against each
    #: other instead: 63.6%/90.0% precise at this threshold, and 75%/100% in
    #: combination with the dispersion test.
    #:
    #: It is a TIER, not a gate, because gating on it would cost 62-76% of the
    #: true novel ends. Everything promoted is emitted; this only decides what
    #: is marked, so a downstream filter can choose the precision it wants.
    min_candidate_share: float = 0.30

    # -- the calibrated terminal score (third scoring level) ---------------
    #: Score 3' ends with the same decoy-anchored mixture used for junctions.
    #:
    #: This level has what the CHAIN level lacks: a clean, orthogonal null.
    #: Positives are ends matching an annotated TTS or a catalogued polyA site;
    #: decoys are A-rich-downstream ends with neither. That is a real
    #: internal-priming population rather than "low support", so unlike
    #: ``chains.max_chain_lfdr`` this one can enforce rather than advise.
    score_terminals: bool = True
    #: Local-FDR ceiling for a NOVEL 3' end. None = advisory (score recorded,
    #: nothing dropped). Start advisory on a new dataset, look at the
    #: distribution, then set it.
    max_terminal_lfdr: Optional[float] = None
    #: %A downstream above which an unannotated, un-catalogued end joins the
    #: decoy set. Deliberately below monoexon.max_perc_a_downstream (60): the
    #: decoy set wants the internal-priming POPULATION, not only the cases
    #: extreme enough to reject outright. BD176c: %A runs 20 (FSM) -> 25
    #: (end3_novel) -> 30 (end3_unresolved), so the signal sits well under 60.
    decoy_perc_a: float = 40.0


@dataclass
class CoverageParams:
    """Short-read coverage as terminal and length evidence. New in 0.3.0.

    Junctions (``short_read_sj``) say which introns exist. Coverage says
    something no junction table can: how the signal BEHAVES across a terminal
    boundary, and whether a long transcript's distal exons are expressed at all.

    Three uses, all optional and all degrading to "no opinion" when the track
    is absent:

    * **3' step** -- a real cleavage site has a sharp coverage drop across it;
      an internal-priming site does not, because transcription continues past.
    * **TSS ratio** -- the SQANTI3 signal the original review asked for.
    * **continuity** -- where long reads stop short of an annotated long
      transcript but coverage runs on at comparable depth, the longer form is
      probably real. This RESCUES a reference model; it never builds a chain no
      read spans, because coverage cannot establish that two distal exons sit
      on the same molecule.
    """

    enabled: bool = True
    #: bp either side of a terminal boundary used for the step/ratio.
    window: int = 100
    #: A 3' end is "step-supported" when downstream mean coverage falls to at
    #: most this fraction of upstream mean.
    max_downstream_ratio: float = 0.30
    #: A 5' end is "ratio-supported" when downstream/upstream exceeds this,
    #: i.e. signal starts here. SQANTI3's ratio_TSS uses a similar cut.
    min_tss_ratio: float = 1.5
    #: Minimum upstream coverage for the step/ratio to mean anything at all.
    min_coverage: float = 3.0
    #: These libraries are UNSTRANDED (measured: both the bulk organoid
    #: alignment and the FLASH-seq plates give forward/reverse ~0.57 each), so
    #: at a locus overlapping an antisense gene the signal cannot be attributed.
    #: Such loci are FLAGGED and left unscored rather than scored wrongly.
    skip_antisense_overlap: bool = True
    #: Treat the track as unstranded. Set False only for a genuinely stranded
    #: pair of tracks.
    unstranded: bool = True


@dataclass
class MonoexonParams:
    """The mono-exonic track. Runs entirely separately; an empty intron chain
    never enters the suffix logic."""

    enabled: bool = True
    min_umis: int = 3
    min_reads: int = 5
    #: single-linkage clustering of unspliced reads is 3'-anchored
    max_3p_diff: int = 100
    min_overlap_frac: float = 0.5
    #: case 1 (last exon / 3'UTR of a known gene): trusted for oligo-dT
    utr3_min_reads: int = 5
    #: case 2/3 (gene body-internal, intergenic): require a polyA site
    internal_require_polya: bool = True
    intergenic_require_polya: bool = True
    intergenic_min_reads: int = 10
    #: %A in the genomic window downstream of the 3' end, on a 0-100 scale.
    #: NOTE: pigeon writes this as a PERCENT.  A 0.6 cutoff flags everything.
    max_perc_a_downstream: float = 60.0
    perc_a_window: int = 20
    polya_motif_window: int = 50
    #: A polyA hexamer sits 10-35 nt upstream of the cleavage site.  Accepting
    #: one anywhere in a 50 nt window roughly doubles the false-positive rate,
    #: and 3'UTRs are AT-rich enough that a stray AATAAA is common.
    polya_motif_min_dist: int = 10
    polya_motif_max_dist: int = 40
    #: When no polyA atlas is supplied, restrict novel calls to the strong
    #: signals.  AAAAAG / AAGAAA are essentially noise in an A-rich 3'UTR.
    polya_motif_strict: bool = True
    #: With an atlas supplied, require a catalogued site -- a motif alone is
    #: not enough to invent a cleavage site.
    require_atlas_when_available: bool = True

    # -- 0.3.0: the track over-called and under-called at the same time -----
    #: On the SIRV masked arm this track emitted 25 models of which 2 were
    #: correct, and in the SAME run rejected the four genuinely novel
    #: mono-exon truths (SIRV205/206/512/617, 642 exact reads between them),
    #: leaving their reads unassigned. One cause: the track decides on genomic
    #: and terminal features plus a flat floor, and never compares a candidate
    #: against the SPLICED models at the same locus. A fragment of a
    #: well-expressed transcript's last exon and a real single-exon transcript
    #: look identical to that code.
    #:
    #: Demotion is the over-calling half: a cluster inside a terminal exon
    #: whose 3' peak coincides with an emitted spliced model's 3' end is a
    #: fragment OF that model, so its reads fold in instead of becoming a
    #: model. The flag this keys on, `rt_dropoff_fragment_of_annotated_3p_end`,
    #: was already being computed and then ignored.
    demote_terminal_exon_fragments: bool = True

    # -- 0.4.0: demote by default, survive on evidence ---------------------
    #: The 0.3.x demotion looks for a fragment that stops WHERE ITS PARENT
    #: STOPS. The SIRV benchmark says that is the wrong shape of rule: of the
    #: 186 mono-exon models surviving in the fair arm, ZERO have a 3' end
    #: within 100 bp of any spliced model's 3' end, so the rule could not fire
    #: on a single one of them -- and only 12 of the 186 have a 3' end within
    #: 100 bp of a real transcript end at all. RT drop-off stops ANYWHERE.
    #:
    #: So containment replaces coincidence: a mono-exonic model lying wholly
    #: inside the span of an emitted spliced model of the same gene is a
    #: fragment of it. On SIRV that is 81 models of which 81 have wrong 3'
    #: ends -- zero correct lost.
    #:
    #: Scope note, and it matters: SIRV Set 4 contains NO true mono-exonic
    #: transcript with a novel 3' end, so this bounds the rule's precision and
    #: says nothing about its recall cost on real tissue, where single-exon
    #: genes and genuine alternative 3'UTR ends both exist. Clause (2) below is
    #: what protects those, and it has not been measured.
    demote_contained_fragments: bool = True
    #: A contained fragment survives demotion on external evidence its position
    #: cannot supply: an annotated 3' end, a catalogued polyA site, or a PAS
    #: hexamer together with a peak tight enough to be a cleavage event
    #: (``TerminalParams.max_dispersion_promote``).
    keep_fragment_if_evidenced: bool = True
    #: Short reads showing NO drop across the putative 3' end mean the
    #: transcript demonstrably continues past it, so this is a truncation.
    #:
    #: Only ever fires on a positive measurement. Absent coverage, or coverage
    #: below ``CoverageParams.min_coverage``, leaves the model to the rules
    #: above -- a lowly expressed gene has no coverage, and treating silence as
    #: evidence of truncation would delete real transcripts in exactly the
    #: genes where long reads matter most.
    demote_on_sr_continuation: bool = True
    #: Rescue is the under-calling half. The polyA gate here is purely genomic
    #: -- a PAS hexamer, downstream A-richness, an atlas hit -- which is a
    #: statement about the LOCUS. The spliced path learned in 0.1.18 that the
    #: read's own soft-clipped tail is a statement about the MOLECULE and can
    #: stand in; the mono track never got that. SIRVs carry a defined 30 nt
    #: tail and were rejected anyway.
    tail_rescues_polya: bool = True
    min_tail_molecule_frac: float = 0.30
    #: ...and a novel mono-exon rescued on tail evidence must also be seen in
    #: this many distinct cells, so one deeply-sequenced artefact cannot do it.
    rescue_min_cells: int = 3
    #: Write one row per REJECTED candidate with the rule that rejected it.
    #: The SIRV audit could not say why those four were dropped, because
    #: nothing recorded it.
    write_rejection_reasons: bool = True


@dataclass
class TssParams:
    """External 5'-end evidence (CAGE / TSS atlas).  New in 0.2.0.

    The 3' path has three independent tiers of evidence -- annotation, a polyA
    site atlas, and the read's own soft-clipped tail.  Through 0.1.18 the 5'
    path had exactly one: proximity to an already-annotated TSS.  That made a
    genuinely novel promoter impossible to keep, because the only escape from
    `chains.merge_unannotated_suffix_into_annotated_parent` was evidence that
    GENCODE already knew the start site.

    An atlas cannot fix that on its own -- a peak says the *locus* initiates
    transcription, not that *this molecule* started at the cap -- but it is the
    first external opinion the 5' end has ever had.

    Direction of effect is one-way and enforced by test: an atlas hit can keep
    a model that would otherwise be merged away; a miss can never drop or
    demote one.
    """

    enabled: bool = True
    #: bp outside a peak that still counts as inside.  0 = strict containment;
    #: the peak is already 20-100 bp wide, so this is for annotation-vs-atlas
    #: coordinate slop, not for widening the evidence.
    peak_slack: int = 0
    #: Distance to a peak's representative TSS that counts as a start-site
    #: match.  Deliberately the same default as ``EndParams.max_5p_diff``: it is
    #: answering the same question ("is this the same start?") against a
    #: different catalogue, and two different numbers would be arbitrary.
    max_dist_to_reptss: int = 100
    #: FANTOM5 column 5 peak score floor.  0 keeps everything; the knob exists
    #: so "only strong peaks" is a measurement rather than an argument.
    min_peak_score: float = 0.0
    #: May atlas evidence alone keep an unannotated suffix from being merged
    #: into its annotated parent?  With this False no merge changes.
    #:
    #: NOTE for the recorded-but-not-acting arm of a measurement: this lever is
    #: not sufficient on its own.  ``use_as_chain_feature`` below re-fits the
    #: chain novelty model, which moves the local FDR and therefore the
    #: candidate set even with no rescue.  A run that must be structurally
    #: identical to one without an atlas needs BOTH off.
    rescue_suffix: bool = True
    #: Feed the evidence tier into the chain-level novelty model as an ordinal
    #: (none 0 / cage_peak 1 / cage_reptss 2 / annotation 3).  See the note
    #: above: this changes model selection independently of ``rescue_suffix``.
    use_as_chain_feature: bool = True


@dataclass
class ScoringParams:
    """Calibrated novelty model."""

    enabled: bool = True
    #: anchors for the supervised step
    positive_anchor: str = "annotated"          # annotated junctions / chains
    negative_anchor: str = "decoy"              # wrong-strand + non-canonical novel
    min_anchor_size: int = 50                   # below this, fall back to thresholds
    l2: float = 1.0                             # ridge penalty on the logistic fit
    max_iter: int = 50
    #: local-FDR estimation
    n_score_bins: int = 60
    pi0_method: str = "decoy"                   # "decoy" | "fixed"
    pi0_fixed: float = 0.8
    smooth_window: int = 5


@dataclass
class MoleculeParams:
    """Barcode/UMI handling.

    ``isoseq tag`` is unusable with variable-length BD Rhapsody barcodes and
    ``groupdedup`` was never run, so the BAM carries no CB/UB tags and the
    mapping arrives as an external read-name-keyed table.
    """

    barcode_umi_tsv: Optional[str] = None
    #: column names in the TSV header
    col_cell_barcode: str = "cell_barcode"
    col_umi: str = "umi"
    col_read_name: str = "read_name"
    #: reconcile read names between BAM and table when a segmentation step
    #: renamed them: none | strip_segment | strip_ccs | zmw
    read_name_normalise: str = "none"
    #: cached compact index (sorted uint64 read-name hashes) written next to the TSV
    index_path: Optional[str] = None
    #: collapse UMIs differing by <= this Hamming distance within a
    #: (cell, transcript) group.  0 disables UMI error correction.
    umi_hamming: int = 1
    #: emit the cell x transcript matrix
    write_matrix: bool = True


@dataclass
class OutputParams:
    outdir: str = "flightcollapse_out"
    prefix: str = "sample"
    #: "ref" -> ENST00000229239.10, ENST...|3UTR_mod1, NIC_GAPDH_1
    #: "pb"  -> PB.1.1 (with the reference-anchored name kept as an attribute)
    id_style: str = "ref"
    write_gtf: bool = True
    write_gff: bool = True
    write_read_stat: bool = True
    write_group: bool = True
    write_abundance: bool = True
    write_models_table: bool = True
    write_group_chains: bool = True
    write_qc: bool = True
    gzip_big_tables: bool = True


#: Sample-type presets. New in 0.3.0.
#:
#: These exist because the external evidence available is a property of the
#: SAMPLE, not of the algorithm, and it changes between runs of the same
#: pipeline. Retinal organoids have matched short reads; the mouse-brain
#: nuclei, HEK benchmarks and SIRV spike-ins do not, and SIRVs additionally
#: have no CAGE atlas and no polyA-site atlas because they are synthetic.
#: Encoding that once, by name, beats remembering six flags per run.
#:
#: What a profile does NOT govern: the intrinsic fixes. Competitive terminal
#: pruning and the terminal score need no external input -- they compare a
#: model against its own siblings -- and the SIRV benchmark is precisely what
#: motivated them, so they are on everywhere. A profile only decides which
#: EXTERNAL channels are expected.
#:
#: ``expects`` is advisory. A declared channel with no file supplied warns and
#: is reported as inactive in the QC; it never fails the run and never silently
#: pretends the evidence was there.
SAMPLE_PROFILES: Dict[str, Dict[str, Any]] = {
    "retinal_organoid": {
        "description": "matched short reads available: junctions and coverage",
        "expects": ["cage_peak_bed", "polya_site_bed", "short_read_sj",
                    "short_read_coverage"],
        "settings": {
            "tss.enabled": True,
            "tss.rescue_suffix": True,
            "coverage.enabled": True,
            "terminal.prune_terminal_variants": True,
            "terminal.score_terminals": True,
            "junctions.use_short_reads_as_anchor": True,
        },
    },
    "generic": {
        "description": "atlases but no short reads (mouse brain, HEK, any new tissue)",
        "expects": ["cage_peak_bed", "polya_site_bed"],
        "settings": {
            "tss.enabled": True,
            "tss.rescue_suffix": True,
            "coverage.enabled": False,
            "terminal.prune_terminal_variants": True,
            "terminal.score_terminals": True,
        },
    },
    "minimal": {
        "description": "no external evidence at all (SIRVs, benchmark parity runs)",
        "expects": [],
        "settings": {
            "tss.enabled": False,
            "coverage.enabled": False,
            "terminal.prune_terminal_variants": True,
            "terminal.score_terminals": True,
        },
    },
}


@dataclass
class SubstitutionParams:
    """Rule J, new in 0.4.3: a novel junction the short reads say did not happen.

    Substitution, never deletion. A chain group whose only difference from an
    annotated transcript is junctions that (a) have no short-read support and
    (b) displace annotated junctions that DO, merges into that transcript's
    group and its reads move with it.

    The short-read conditions are both required, and that is the whole design.
    ``sr_uniq == 0`` alone covers 46.7% of BD67's novel-junction read mass --
    zero short reads is the normal state of a novel junction when the SR library
    is a different, shallower batch. Only the conjunction with a well-covered
    annotated competitor turns absence of evidence into evidence of absence.

    Inert with no SJ tables: ``min_competitor_sr`` cannot be met, so HEK, mouse
    brain and SIRV runs are untouched.
    """

    enabled: bool = True
    #: how far a displaced site may sit inside the intron. 30 bp is the
    #: microexon ceiling; ARR3 is +12 over a 10 bp exon, BSG +19 over 25 bp.
    max_shift: int = 30
    #: the displaced junction must have at most this many unique short reads
    max_novel_sr: int = 0
    #: each annotated junction it displaces must have at least this many, which
    #: is what proves the locus was visible to the short reads at all
    min_competitor_sr: int = 10
    #: write <prefix>.junction_substitutions.tsv -- one row per moved junction,
    #: because a read move should be auditable rather than inferred from a diff
    write_plan: bool = True


@dataclass
class ParallelParams:
    """Contig-level parallelism.  New in 0.5.0.

    ``run()`` already streams one contig at a time and releases its read arrays
    before the next, and -- this is the part that makes parallelism safe here --
    all three calibrated models (``score_junctions``, ``score_chains``,
    ``score_terminals``) are fitted **per contig**.  Nothing in a contig's
    result depends on another contig, so the only shared state is the
    accumulators, the ``IdAssigner`` and the read-assignment writer.

    Determinism is a requirement, not a bonus: workers return results, the
    parent merges them in the configured contig order, ids are assigned in the
    parent, and ``read_stat``/``group`` are stitched from per-contig shards in
    the same order.  A parallel run is therefore byte-identical to a serial one.

    Memory is the binding constraint, not cores.  Each worker holds one contig's
    read arrays plus its own cached chromosome sequence (~250 MB as a Python str
    for chr1), so peak RSS is roughly ``workers x per-contig peak``.  The
    default clamps the worker count to what the machine can actually hold.
    """

    #: 0 = os.cpu_count(); 1 = the serial path (identical code, no fork)
    workers: int = 1
    #: Only "fork" shares the parsed ReferenceIndex for free.  On a platform
    #: without it the run falls back to serial with a warning rather than
    #: re-parsing the GTF in every worker.
    start_method: str = "fork"
    #: Clamp ``workers`` by available memory as well as by cores.
    clamp_by_memory: bool = True
    #: Memory budgeted per worker, GB.  One contig's reads plus one cached
    #: chromosome; chr1 of a deep FLNC BAM is the worst case.
    gb_per_worker: float = 8.0
    #: Leave this much for the parent, which holds the reference index, every
    #: contig's models and the molecule matrix.
    reserve_gb: float = 8.0
    #: Dispatch the longest contigs first.  Wall-clock is floored by the
    #: slowest single contig, so starting chr1 last wastes the whole tail.
    longest_first: bool = True
    #: Where per-contig shards live.  None = a temp dir under output.outdir,
    #: removed when the run finishes.
    shard_dir: Optional[str] = None
    #: Keep the shards after stitching.  Debug only.
    keep_shards: bool = False


@dataclass
class SecondPassParams:
    """Filtering of the first-pass catalogue.  New in 0.5.0.

    The model-*building* algorithm is deliberately untouched.  This stage only
    accepts, rewrites or removes models that the first pass already produced,
    and every decision lands in ``<prefix>.secondpass.tsv`` with its reason.

    Two rules, in order:

    **1. Mono-exonic models are judged against their gene's exon structure.**
    A single-exon model is credible when the gene it sits in is itself
    single-exon (class A) or has both single- and multi-exon isoforms (class
    B).  In a gene that is multi-exon only, a single-exon model is a fragment
    or a mis-assignment however good its 3' end looks -- it spans an intron
    (C), sits inside one exon (D), or is something else entirely (E).  Those
    are removed, and the gene's most-supported isoform is recorded in their
    place, marked INFERRED so that it is never mistaken for an observation.

    **2. ``end3_novel`` models must earn their 3' end.**  The SIRV Set 4
    benchmark put this category at 0/12 correct in the reference arm and 0/11
    in the masked arm: internal priming on genomic A-stretches produces read
    piles that accumulate support exactly like real cleavage sites.  The first
    pass already checks the polyA atlas and a hexamer; this stage adds the two
    measurements that actually separate the two populations, and falls back to
    a gene-dominance rule when neither can speak.
    """

    enabled: bool = True
    #: One row per second-pass decision, kept or removed.
    write_audit: bool = True

    # -- rule 1: mono-exon gene-context classes ----------------------------
    monoexon_classes: bool = True
    #: Classes removed and replaced by an inferred reference isoform.
    #: A = gene is mono-exonic only; B = gene has mono + multi isoforms;
    #: C = multi-exon gene, model spans an intron; D = multi-exon gene, model
    #: lies inside one exon; E = multi-exon gene, anything else.
    remove_classes: List[str] = field(default_factory=lambda: ["C", "D", "E"])
    #: Bases of *constitutively* intronic sequence -- inside the gene span and
    #: exonic in no isoform -- that a model must cover to count as class C.
    #: Using the constitutive definition rather than any one transcript's
    #: introns is what stops an ordinary alternative-exon model being called a
    #: read-through.
    min_intron_overlap: int = 10
    #: A mono-exon model with no gene has no exon structure to be judged
    #: against, so it is out of scope for A-E.  "keep" leaves the existing
    #: monoexon_intergenic handling alone; "remove" drops the class entirely.
    intergenic_policy: str = "keep"
    #: Where the replacement comes from.  "data_then_reference" takes the
    #: gene's best-supported *emitted* multi-exon model -- which is what this
    #: sample actually shows -- and falls back to the reference canonical
    #: (MANE_Select > Ensembl_canonical > basic > longest) when the gene has no
    #: emitted spliced model.  "reference" always takes the canonical.
    infer_from: str = "data_then_reference"
    #: Category written on a model that exists only by inference.
    inferred_category: str = "INFERRED"
    #: Move the removed model's reads onto the inferred transcript.  With this
    #: off the reads become unassigned and are counted, which makes the
    #: inference obviously evidence-free but loses the read mass.
    fold_reads_into_inferred: bool = True

    # -- rule 2: the end3_novel terminal gate ------------------------------
    end3_filter: bool = True
    #: Categories this gate judges.  ``end3_unresolved`` is not emitted by
    #: default (ends.emit_unresolved_3p = False); add it here if you turn it on.
    end3_categories: List[str] = field(default_factory=lambda: ["end3_novel"])

    # tier 1: a catalogued cleavage site
    use_atlas: bool = True

    # tier 2: the templated-tail test -- the one per-MOLECULE discriminator.
    #
    # A genuine polyA tail is not in the genome: the A-run in the read's 3'
    # soft clip greatly exceeds the genomic A-run at that coordinate.  Internal
    # priming is templated -- the aligner consumes the genomic A's and the
    # excess collapses to nothing.  `tail_molecule_frac` cannot see this
    # because it scores an A-run without asking what the genome says underneath.
    use_templated_tail: bool = True
    #: Bases by which a read's terminal A-run must exceed the genomic A-run at
    #: its own 3' end before that read counts as non-templated.
    min_tail_excess: int = 8
    #: ...and the share of the model's reads that must clear it.
    min_nontemplated_frac: float = 0.50
    #: How far downstream to scan for the genomic A-run.  Longer than any
    #: plausible tail would need, so the run length is never truncated by the
    #: window.
    a_run_scan: int = 60

    # tier 3: positional composition around the cleavage site.
    #
    # A real 3' end has a U/GU-rich downstream element at roughly +5..+30 and a
    # U-rich stretch upstream; an internal-priming site has A-rich downstream
    # and no DSE.  A single `perc_a_downstream` number cannot express that --
    # measured on BD176c it moves only 20 -> 25 between FSM and end3_novel.
    #
    # The tier is deliberately built as positive evidence (a hexamer) plus a
    # tight peak plus the ABSENCE of the internal-priming signature, rather
    # than as a claim that a real site must look U-rich.  That matters for the
    # threshold below.
    use_composition: bool = True
    dse_start: int = 5
    dse_end: int = 30
    #: (G+T) minus A fraction in the downstream element.
    #:
    #: A contrast, not an absolute G+U cut, and the reason is measurement
    #: noise: the window is ~25 bases, so an absolute threshold at 0.45 is one
    #: base away from 0.40 and neutral sequence crosses it by chance about half
    #: the time.  The contrast asks the question the tier actually cares about
    #: -- is this downstream sequence A-rich? -- and neutral sequence sits at
    #: about +0.25 while an internal-priming site goes sharply negative, so the
    #: two populations are separated by far more than the sampling error.
    #:
    #: 0.10 is a starting point, not a measurement.  Look at the distribution
    #: of `dse_gu_frac` and `dse_a_frac` over your FSM models -- which are
    #: written to models.tsv for exactly this purpose -- before trusting it.
    min_dse_gu_minus_a: float = 0.10
    upstream_u_window: int = 40
    #: Peak dispersion at or below which the end is tight enough to be a
    #: cleavage event.  BD176c: FSM 4.1 bp against end3_novel 34.6 bp.
    max_dispersion: float = 12.0

    # tier 4: short-read coverage must actually drop across the end
    use_short_read_step: bool = True

    # tier 5: the calibrated terminal score.  None = not used.  Left off by
    # default because the 0.4.6 SIRV-derived rescoring did not beat ranking on
    # raw molecule count; turn it on once you have looked at the distribution.
    max_terminal_lfdr: Optional[float] = None

    # tier 6: the hard backstop.  Survive if this is the ONLY model of its gene,
    # or the best-supported one.  Deliberately last: on its own it is permissive
    # in the wrong direction -- a highly expressed gene with a strong internal
    # priming site makes that artefact the dominant model -- so it only ever
    # rescues a model the evidence tiers could not speak to.
    gene_dominance_backstop: bool = True
    #: How far ahead of the runner-up "most supported" has to be.  1.0 = a
    #: strict maximum; raise it to demand a clear winner.
    dominance_ratio: float = 1.0
    #: Move a dropped model's reads onto the best compatible surviving sibling
    #: (same intron chain first, then the gene's best model).  Off = the reads
    #: become unassigned and are counted as such.
    fold_dropped_reads: bool = True


@dataclass
class PosthocParams:
    """Final catalogue consolidation.  New in 0.6.0; see ``posthoc.py``.

    Two independent benchmarks set these, and they disagreed about one number.

    SIRV Set 4 (BD144_kin, 69 conventional SIRVs; exact match = junctions
    +/- 5 bp, ends +/- 50 bp) set the first defaults.  Against the same run with
    both rules off (0.28 / 0.84 / 0.42 precision / recall / F1 without
    annotation), 500 bp gave 0.82 / 0.84 / 0.83 without annotation and
    0.86 / 0.87 / 0.86 with it, recall unchanged in both arms.  The fragment
    filter alone, at end3_tolerance = 0, already takes annotation-free precision
    from 0.28 to 0.54 at no cost in recall; that result is unchanged and is why
    the filter is on by default.

    0.6.0 shipped 500 on BD144_kin alone, with a note that HM and ND would
    confirm it before release.  They did not.  Three-way reproducibility across
    BD144_kin / HM / ND (202,742 pooled transcript structures; reproducible read
    mass as the measure, since a share always falls when a consolidation merges
    two reproducible structures into one) costs

        200 bp   -16.3% models, reproducible read mass -0.32 pp, +0.11 pp
                 at >= 50 reads and >= 100 barcodes, positive in 34 of 81 cells
        500 bp   -21.5% models, reproducible read mass -1.11 pp, -0.25 pp
                 at the same cutoff, positive in 17 of 81 cells

    so the extra 300 bp buys consolidation with reproducible mass.  At the
    splice-chain level the two are indistinguishable (122,270 against 122,174
    chains), which is what says the tolerance is doing 3' work only.  The
    default is therefore 200 as of 0.6.3; 500 remains reasonable where
    precision on a known panel matters more than replicate agreement.

    The 69-transcript SIRV panel has ground truth and no 3'-end diversity worth
    the name; the HEK libraries have 418,097 models and no ground truth.  They
    measure different things, and both numbers are kept here on purpose.
    """

    enabled: bool = True
    write_audit: bool = True
    #: 3' end tolerance (bp) for merging models with an IDENTICAL intron chain.
    #: 200 is the knee of the three-library reproducibility curve (see above);
    #: 500 is the largest tolerance that cost no exact-match recall in either
    #: SIRV arm.  A negative value switches the 3' rule off, 0 merges only
    #: identical 3' ends.
    end3_tolerance: int = 200
    #: 5' window (bp) that two same-chain models must also share to merge. A
    #: guard on rule 1, not a second tolerance. NEGATIVE removes it, which makes
    #: the 3' end the whole rule and leaves 5' truncation to the sub-chain
    #: filter -- at the price of merging two identical chains that differ only
    #: in where transcription started, since the sub-chain filter never
    #: compares two models that share a chain exactly.  Off by default as of
    #: 0.6.3: across the three HEK libraries removing it freed 24 to 5,664
    #: merges per arm for a read-mass difference of at most 0.03 pp, and it
    #: makes the two ends' treatment explicit -- the 3' end gets a tolerance,
    #: the 5' end gets the structural sub-chain filter and nothing else.
    end5_window: int = -1
    #: Reference-free removal of models whose chain is a contiguous sub-chain
    #: of a better-supported model (5'-truncation / fragment class).
    subchain_filter: bool = True
    #: Container support must be >= ratio x the fragment's support.
    subchain_ratio: float = 2.0
    #: Junction tolerance (bp) when matching a sub-chain.
    subchain_fuzzy: int = 5
    #: A fragment's terminal exons must lie inside the container's matching
    #: exons, +/- this many bp. Protects alternative first/last exons that
    #: start inside an upstream intron, which are not truncations.
    subchain_terminal_slack: int = 50
    #: Categories never removed by the sub-chain rule (an exact annotated
    #: intron chain is evidence in itself).
    protected_categories: List[str] = field(default_factory=lambda: ["FSM"])
    #: Move a removed fragment's reads onto its container (True) or leave
    #: them unassigned and count them (False).
    fold_reads: bool = True


@dataclass
class Config:
    bam: str = ""
    reference_gtf: str = ""
    genome_fasta: Optional[str] = None     # required for motif/%A features
    polya_site_bed: Optional[str] = None   # PolyASite / PolyA_DB, optional but recommended
    #: CAGE / TSS peak atlas (FANTOM5 phase1&2, refTSS).  BED6 or BED9; a BED9
    #: carries its representative TSS in columns 7/8.
    cage_peak_bed: Optional[str] = None
    #: Representative TSS as a separate BED6, joined on the peak name.  Needed
    #: only when the peak file lost its thick fields -- which is what happens
    #: when a FANTOM5 mouse file is lifted from mm10 to GRCm39.
    cage_reptss_bed: Optional[str] = None
    #: STAR ``SJ.out.tab`` files. Several are summed -- technical replicates are
    #: the same molecules sampled twice, so pooling them gives one junction's
    #: evidence rather than three under-powered opinions of it.
    #:
    #: For plate-based data, pool the per-well tables into one file PER PLATE
    #: and pass the plates. ``min_sr_replicates`` then means "seen in >= N
    #: plates", which is a real replicate structure. Passing every well
    #: separately would make it mean "seen in >= N cells" -- a stronger
    #: independence claim, but the loader keeps a dict per library and 1,500 of
    #: them is not worth it.
    short_read_sj: Optional[List[str]] = None
    #: Short-read COVERAGE tracks (bigWig). New in 0.3.0; needs pyBigWig.
    #: Junctions say which introns exist; coverage says how signal behaves
    #: across a terminal boundary, which no junction table can.
    short_read_coverage: Optional[List[str]] = None
    #: Which sample this is, naming one of SAMPLE_PROFILES. Decides which
    #: external evidence channels are expected; never changes the intrinsic
    #: algorithm. Applied by `apply_sample_profile()`, which the CLI calls
    #: before any explicit --set override, so --set always wins.
    #: 0.4.2. Write the SQANTI3 structural_category / subcategory / terminal_class
    #: columns beside `category`. Labels only -- it runs after every model is
    #: final and cannot change which models exist. False omits the columns.
    sqanti_classification: bool = True
    sample_type: str = "generic"

    contigs: Optional[List[str]] = None    # None -> all except `exclude_contigs`
    exclude_contigs: List[str] = field(default_factory=lambda: ["chrM", "MT", "M"])

    gates: AlignmentGates = field(default_factory=AlignmentGates)
    junctions: JunctionParams = field(default_factory=JunctionParams)
    substitution: SubstitutionParams = field(default_factory=SubstitutionParams)
    chains: ChainParams = field(default_factory=ChainParams)
    ends: EndParams = field(default_factory=EndParams)
    tss: TssParams = field(default_factory=TssParams)
    terminal: TerminalParams = field(default_factory=TerminalParams)
    coverage: CoverageParams = field(default_factory=CoverageParams)
    monoexon: MonoexonParams = field(default_factory=MonoexonParams)
    scoring: ScoringParams = field(default_factory=ScoringParams)
    molecules: MoleculeParams = field(default_factory=MoleculeParams)
    output: OutputParams = field(default_factory=OutputParams)
    parallel: ParallelParams = field(default_factory=ParallelParams)
    secondpass: SecondPassParams = field(default_factory=SecondPassParams)
    posthoc: PosthocParams = field(default_factory=PosthocParams)

    #: Deprecated alias for ``parallel.workers``, kept because it was already a
    #: config key.  Used only when ``parallel.workers`` is left at 1.
    threads: int = 1
    seed: int = 0
    verbose: bool = True
    #: fail the run if a validation invariant is violated
    strict_invariants: bool = True

    # ------------------------------------------------------------------ #
    def apply_sample_profile(self, name: Optional[str] = None) -> Dict[str, Any]:
        """Apply a SAMPLE_PROFILES preset. Returns what it did, for the report.

        Called before any explicit ``--set`` override, so a hand-set value
        always wins over the preset. Missing expected inputs are reported, not
        raised: a run on a sample that happens to lack short reads is a normal
        run, it just has fewer channels, and the report says which.
        """
        name = name or self.sample_type
        if name not in SAMPLE_PROFILES:
            raise ValueError(
                f"unknown sample_type {name!r}; known: "
                + ", ".join(sorted(SAMPLE_PROFILES))
            )
        self.sample_type = name
        prof = SAMPLE_PROFILES[name]
        for key, value in prof["settings"].items():
            self.set_value(key, value)
        missing = [k for k in prof["expects"] if not getattr(self, k, None)]
        return {
            "sample_type": name,
            "description": prof["description"],
            "expects": list(prof["expects"]),
            "missing_inputs": missing,
        }

    def evidence_channels(self) -> Dict[str, Any]:
        """Which external evidence is actually active in this run.

        Written into the QC report so that a catalogue can never be compared
        against another one without knowing what each had to work with. The
        0.2.0 lesson: two runs that differ only in their inputs are not the
        same experiment, and nothing else in the outputs says so.
        """
        return {
            "sample_type": self.sample_type,
            "cage_atlas": bool(self.cage_peak_bed) and self.tss.enabled,
            "polya_atlas": bool(self.polya_site_bed),
            "genome_fasta": bool(self.genome_fasta),
            "short_read_junctions": bool(self.short_read_sj),
            "short_read_coverage": bool(self.short_read_coverage) and self.coverage.enabled,
            "barcode_umi": bool(self.molecules.barcode_umi_tsv),
            "terminal_pruning": self.terminal.prune_terminal_variants,
            "terminal_scoring": self.terminal.score_terminals,
            "second_pass": self.secondpass.enabled,
            "second_pass_monoexon_classes": (
                self.secondpass.enabled and self.secondpass.monoexon_classes
            ),
            "second_pass_end3_filter": (
                self.secondpass.enabled and self.secondpass.end3_filter
            ),
            "templated_tail_test": (
                self.secondpass.enabled
                and self.secondpass.end3_filter
                and self.secondpass.use_templated_tail
                and bool(self.genome_fasta)
                and self.gates.measure_polya_tail
            ),
            "workers": self.worker_count(),
            "live_end3_tiers": self.live_end3_tiers(),
        }

    # ------------------------------------------------------------------ #
    def live_end3_tiers(self) -> List[str]:
        """Which tiers of the end3 gate actually have their inputs in this run.

        A tier without its input has no opinion; it is not a rejection.  The
        distinction matters because the tiers are tried strongest-first and a
        model that falls through to the backstop because the atlas is missing
        is a very different statement from one that falls through because the
        atlas was consulted and said no.  Reported so the two can be told
        apart from the QC file alone.
        """
        sp = self.secondpass
        live: List[str] = []
        if sp.use_atlas and self.polya_site_bed:
            live.append("atlas")
        if (sp.use_templated_tail and self.genome_fasta
                and self.gates.measure_polya_tail):
            live.append("nontemplated")
        if sp.use_composition and self.genome_fasta:
            live.append("composition")
        if (sp.use_short_read_step and self.short_read_coverage
                and self.coverage.enabled):
            live.append("short_read")
        if sp.max_terminal_lfdr is not None and self.terminal.score_terminals:
            live.append("calibrated")
        return live

    # ------------------------------------------------------------------ #
    def worker_count(self) -> int:
        """How many contig workers this run will actually use.

        Resolves ``0`` to the core count, honours the deprecated ``threads``
        alias, and -- because memory is the binding constraint here, not cores
        -- clamps to what the machine can hold.  Returning the resolved number
        rather than the requested one means the QC report says what happened,
        not what was asked for.
        """
        want = self.parallel.workers
        if want == 1 and self.threads > 1:
            want = self.threads
        if want == 0:
            want = os.cpu_count() or 1
        want = max(1, int(want))
        if want > 1 and self.parallel.clamp_by_memory:
            free = _available_gb()
            if free is not None:
                budget = max(0.0, free - self.parallel.reserve_gb)
                by_mem = int(budget // self.parallel.gb_per_worker)
                want = max(1, min(want, by_mem))
        return want

    def set_value(self, dotted: str, value: Any) -> None:
        """``cfg.set_value("tss.enabled", True)`` -- typed, unlike set_path."""
        obj: Any = self
        parts = dotted.split(".")
        for part in parts[:-1]:
            obj = getattr(obj, part)
        if not hasattr(obj, parts[-1]):
            raise ValueError(f"unknown config key: {dotted}")
        setattr(obj, parts[-1], value)

    # ------------------------------------------------------------------ #
    def validate(self) -> None:
        """Reject configurations where a setting would be silently ignored."""
        problems = []
        if not self.bam:
            problems.append("config.bam is required")
        if not self.reference_gtf:
            problems.append("config.reference_gtf is required")

        needs_genome = (
            self.junctions.require_canonical_motif
            or self.junctions.rt_switch_repeat_len > 0
            or self.monoexon.enabled
            or self.ends.require_polya_evidence_for_utr_variant
            or self.scoring.enabled
        )
        if needs_genome and not self.genome_fasta:
            problems.append(
                "genome_fasta is required by: canonical-motif checks, RT-switch "
                "detection, polyA evidence and the scoring model. Either supply "
                "it or turn those off explicitly "
                "(junctions.require_canonical_motif=false, "
                "junctions.rt_switch_repeat_len=0, monoexon.enabled=false, "
                "ends.require_polya_evidence_for_utr_variant=false, "
                "scoring.enabled=false)."
            )
        if self.molecules.write_matrix and not self.molecules.barcode_umi_tsv:
            problems.append(
                "molecules.write_matrix=true but molecules.barcode_umi_tsv is unset; "
                "set the table or molecules.write_matrix=false"
            )
        if self.ends.five_prime_percentile <= 50:
            problems.append(
                "ends.five_prime_percentile <= 50 pushes models towards the most "
                "truncated members. That is the isoseq collapse bug; use >= 75."
            )
        if self.chains.debris_ratio >= 1.0:
            problems.append("chains.debris_ratio must be < 1")
        if self.chains.dominant_ratio <= 1.0:
            problems.append("chains.dominant_ratio must be > 1")
        # -- 0.2.0: the CAGE atlas, under the same no-silent-no-op rule -------
        if self.cage_reptss_bed and not self.cage_peak_bed:
            problems.append(
                "cage_reptss_bed is set but cage_peak_bed is not. The "
                "representative-TSS table is joined onto the peaks by name; on "
                "its own it is never read."
            )
        if self.cage_peak_bed and not self.tss.enabled:
            problems.append(
                "cage_peak_bed is set but tss.enabled=false, so the atlas would "
                "be loaded and ignored. Drop the file or enable the tier."
            )
        if not self.cage_peak_bed:
            if self.tss.min_peak_score:
                problems.append(
                    "tss.min_peak_score has no atlas to filter; set cage_peak_bed "
                    "or leave the score at 0"
                )
            if self.tss.peak_slack:
                problems.append(
                    "tss.peak_slack has no atlas to widen; set cage_peak_bed or "
                    "leave the slack at 0"
                )
        if self.tss.max_dist_to_reptss < 0:
            problems.append("tss.max_dist_to_reptss must be >= 0")
        # -- 0.3.0: terminal pruning and the coverage tracks -------------------
        if self.sample_type not in SAMPLE_PROFILES:
            problems.append(
                f"sample_type {self.sample_type!r} is not one of: "
                + ", ".join(sorted(SAMPLE_PROFILES))
            )
        if self.terminal.sibling_support_ratio < 1.0:
            problems.append(
                "terminal.sibling_support_ratio < 1 would absorb a BETTER "
                "supported model into a worse one; use >= 1"
            )
        if self.terminal.max_terminal_lfdr is not None and not self.terminal.score_terminals:
            problems.append(
                "terminal.max_terminal_lfdr is set but terminal.score_terminals "
                "is false, so nothing would compute the score it filters on"
            )
        if self.short_read_coverage and not self.coverage.enabled:
            problems.append(
                "short_read_coverage is set but coverage.enabled is false, so "
                "the tracks would be loaded and ignored"
            )
        if self.output.id_style not in ("ref", "pb"):
            problems.append("output.id_style must be 'ref' or 'pb'")
        # -- 0.5.0: parallelism and the second pass, same no-silent-no-op rule -
        if self.parallel.workers < 0:
            problems.append("parallel.workers must be >= 0 (0 = all cores)")
        if self.parallel.gb_per_worker <= 0:
            problems.append("parallel.gb_per_worker must be > 0")
        sp = self.secondpass
        if sp.enabled:
            bad = [c for c in sp.remove_classes if c not in ("A", "B", "C", "D", "E")]
            if bad:
                problems.append(
                    f"secondpass.remove_classes contains unknown classes {bad}; "
                    "valid classes are A B C D E"
                )
            if sp.intergenic_policy not in ("keep", "remove"):
                problems.append(
                    "secondpass.intergenic_policy must be 'keep' or 'remove'"
                )
            if sp.infer_from not in ("data_then_reference", "reference"):
                problems.append(
                    "secondpass.infer_from must be 'data_then_reference' or "
                    "'reference'"
                )
            if sp.dominance_ratio < 1.0:
                problems.append(
                    "secondpass.dominance_ratio < 1 would let a WORSE supported "
                    "model count as the gene's dominant one; use >= 1"
                )
            if sp.dse_end <= sp.dse_start:
                problems.append("secondpass.dse_end must be > secondpass.dse_start")
            if sp.min_nontemplated_frac > 1.0 or sp.min_nontemplated_frac < 0.0:
                problems.append(
                    "secondpass.min_nontemplated_frac is a fraction of a model's "
                    "reads and must lie in [0, 1]"
                )
            if sp.max_terminal_lfdr is not None and not self.terminal.score_terminals:
                problems.append(
                    "secondpass.max_terminal_lfdr is set but "
                    "terminal.score_terminals is false, so nothing would compute "
                    "the score it filters on"
                )
            # A tier whose input is absent does not fail -- a run without a
            # polyA atlas or without short reads is a normal run, and each
            # channel simply has no opinion (``live_end3_tiers`` in the QC
            # report says which ones spoke). What must never happen silently is
            # a gate with NO live tier and no backstop: that deletes every
            # model in the category while looking like a filter.
            if sp.end3_filter and not self.live_end3_tiers():
                if not sp.gene_dominance_backstop:
                    problems.append(
                        "secondpass.end3_filter is on, no evidence tier has its "
                        "inputs (needs genome_fasta with "
                        "gates.measure_polya_tail, a polya_site_bed, short-read "
                        "coverage, or terminal.max_terminal_lfdr) and "
                        "secondpass.gene_dominance_backstop is off -- every "
                        f"model in {sp.end3_categories} would be removed "
                        "unconditionally. Supply an input, enable the backstop, "
                        "or set secondpass.end3_filter=false."
                    )
        ph = self.posthoc
        if ph.end5_window < 0 and not ph.subchain_filter:
            problems.append(
                "posthoc.end5_window is negative, which removes the 5' guard "
                "from the 3' rule, but posthoc.subchain_filter is off -- then "
                "nothing is looking at the 5' end at all")
        if ph.subchain_ratio < 1.0:
            problems.append("posthoc.subchain_ratio must be >= 1 (a container "
                            "less supported than its fragment is not a container)")
        if ph.subchain_fuzzy < 0 or ph.subchain_terminal_slack < 0:
            problems.append("posthoc.subchain_fuzzy and subchain_terminal_slack must be >= 0")
        if problems:
            raise ValueError("invalid configuration:\n  - " + "\n  - ".join(problems))

    # ------------------------------------------------------------------ #
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def dump(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2, default=list)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Config":
        sub = {
            "gates": AlignmentGates,
            "junctions": JunctionParams,
            "substitution": SubstitutionParams,
            "chains": ChainParams,
            "ends": EndParams,
            "tss": TssParams,
            "terminal": TerminalParams,
            "coverage": CoverageParams,
            "monoexon": MonoexonParams,
            "scoring": ScoringParams,
            "molecules": MoleculeParams,
            "output": OutputParams,
            "parallel": ParallelParams,
            "secondpass": SecondPassParams,
            "posthoc": PosthocParams,
        }
        kw: Dict[str, Any] = {}
        known = {f.name for f in fields(cls)}
        for k, v in d.items():
            if k not in known:
                raise ValueError(f"unknown config key: {k}")
            kw[k] = sub[k](**v) if k in sub and isinstance(v, dict) else v
        return cls(**kw)

    @classmethod
    def load(cls, path: str) -> "Config":
        with open(path) as fh:
            text = fh.read()
        if path.endswith((".yaml", ".yml")):
            try:
                import yaml  # type: ignore
            except ImportError as exc:  # pragma: no cover
                raise SystemExit("PyYAML is needed for YAML configs; use JSON instead") from exc
            d = yaml.safe_load(text)
        else:
            d = json.loads(text)
        return cls.from_dict(d)

    # ------------------------------------------------------------------ #
    def set_path(self, dotted: str, value: str) -> None:
        """``--set ends.max_3p_diff=50`` style override with type coercion."""
        obj: Any = self
        parts = dotted.split(".")
        for p in parts[:-1]:
            obj = getattr(obj, p)
        name = parts[-1]
        if not hasattr(obj, name):
            raise ValueError(f"unknown config key: {dotted}")
        cur = getattr(obj, name)
        setattr(obj, name, _coerce(value, cur))

    def outpath(self, suffix: str) -> str:
        os.makedirs(self.output.outdir, exist_ok=True)
        return os.path.join(self.output.outdir, f"{self.output.prefix}{suffix}")


def _available_gb() -> Optional[float]:
    """Memory this machine can actually give a worker, GB.

    Reads ``MemAvailable`` rather than ``MemFree`` -- the kernel's own estimate
    of what is obtainable without swapping, which is the question being asked.
    Returns None where it cannot be answered, and the caller then trusts the
    requested worker count rather than guessing.
    """
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / (1024.0 * 1024.0)
    except OSError:
        pass
    try:  # pragma: no cover - platform dependent
        return (os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
                / (1024.0 ** 3))
    except (ValueError, OSError, AttributeError):
        return None


def _coerce(value: str, template: Any) -> Any:
    if value.lower() in ("none", "null"):
        return None
    if isinstance(template, bool) or value.lower() in ("true", "false"):
        return value.lower() == "true"
    if isinstance(template, int) and not isinstance(template, bool):
        return int(value)
    if isinstance(template, float):
        return float(value)
    if isinstance(template, (list, tuple)):
        return type(template)(v.strip() for v in value.split(","))
    return value
