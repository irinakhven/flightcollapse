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

    min_aln_coverage: float = 0.90      # aligned query bases / full query length
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
    emit_unresolved_3p: bool = True

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


@dataclass
class Config:
    bam: str = ""
    reference_gtf: str = ""
    genome_fasta: Optional[str] = None     # required for motif/%A features
    polya_site_bed: Optional[str] = None   # PolyASite / PolyA_DB, optional but recommended
    #: STAR ``SJ.out.tab`` files. Several are summed -- technical replicates are
    #: the same molecules sampled twice, so pooling them gives one junction's
    #: evidence rather than three under-powered opinions of it.
    short_read_sj: Optional[List[str]] = None

    contigs: Optional[List[str]] = None    # None -> all except `exclude_contigs`
    exclude_contigs: List[str] = field(default_factory=lambda: ["chrM", "MT", "M"])

    gates: AlignmentGates = field(default_factory=AlignmentGates)
    junctions: JunctionParams = field(default_factory=JunctionParams)
    chains: ChainParams = field(default_factory=ChainParams)
    ends: EndParams = field(default_factory=EndParams)
    monoexon: MonoexonParams = field(default_factory=MonoexonParams)
    scoring: ScoringParams = field(default_factory=ScoringParams)
    molecules: MoleculeParams = field(default_factory=MoleculeParams)
    output: OutputParams = field(default_factory=OutputParams)

    threads: int = 1
    seed: int = 0
    verbose: bool = True
    #: fail the run if a validation invariant is violated
    strict_invariants: bool = True

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
        if self.output.id_style not in ("ref", "pb"):
            problems.append("output.id_style must be 'ref' or 'pb'")
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
            "chains": ChainParams,
            "ends": EndParams,
            "monoexon": MonoexonParams,
            "scoring": ScoringParams,
            "molecules": MoleculeParams,
            "output": OutputParams,
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
