# Technical notes from the v0.1.18 source archive

This is the detailed README shipped inside the original v0.1.18 source
archive. It is retained for implementation context. The top-level README is the
current installation and evaluation guide; its test count and validation note
supersede historical counts and unavailable internal-document links below.

Reference-anchored collapse of aligned long reads into transcript models, for
oligo-dT single-cell Iso-Seq (FLIGHT-seq) data. A replacement for
`isoseq collapse` that does not invent transcriptional diversity.

## Why

`isoseq collapse` in single-cell mode merges isoforms that differ only by extra
5' exons — a *suffix* relation on the intron chain — and then writes the
**shortest** member of each group as the representative. Two consequences:

* the empty intron chain of a mono-exonic read is vacuously a suffix of every
  chain, so a single-exon stub can absorb unbounded full-length read mass.
  Measured on BD144_kin: **9.2%** of aligned reads are genuinely unspliced,
  but **57.9%** of read mass landed on single-exon models;
* `--max-5p-diff` is silently inactive unless `--do-not-collapse-extra-5exons`
  is set, so a passed value is accepted and ignored.

This package fixes both, and adds the things that were missing when the problem
had to be diagnosed after the fact: read-name-keyed outputs, per-group chain
heterogeneity, and validation invariants that fail the build.

## Install

```bash
pip install -e .            # needs pysam, numpy, pandas
pytest -q                   # 80 tests, ~1 s, all on simulated data with ground truth
```

## Use

```bash
flightcollapse run \
    -b BD144_kin.fltnc.aligned.bam \
    -g /data/references/gencode.v44.annotation.gtf \
    -f /data/references/GRCh38.primary_assembly.genome.fa \
    -u BD144_kin_barcode_umi_strict.tsv \
    --polya-bed polyasite_2.0_hg38.bed \
    -o out -p BD144_kin
```

Everything is configurable and nothing is silently ignored:

```bash
flightcollapse config > my.json           # fully populated defaults
flightcollapse run -c my.json --set ends.max_3p_diff=50 --set chains.min_chain_umis=5
```

Validate **any** collapse output — this tool's, isoseq's, bambu's, StringTie's:

```bash
flightcollapse check -b reads.bam -g collapsed.gff -r read_stat.txt
```

It exits non-zero when the fraction of read mass on mono-exon models does not
match the fraction of genuinely unspliced reads.

## If molecule counts come back empty

The barcode/UMI table joins to the BAM by read name, and a mismatch there is
silent — every count is zero and every threshold quietly falls back to raw
reads. Rather than debugging it with shell one-liners, ask:

```bash
flightcollapse diagnose-umi -b reads.bam -u sample_barcode_umi_strict.tsv
```

It samples the table **from random byte offsets across the whole file**, tries
every column that looks like a read name under every normalisation mode, and
prints the hit counts side by side.

One trap it exists to avoid: these tables are usually written by streaming a
*coordinate-sorted* BAM, so `head -200000` is not a sample of the library, it is
the start of chr1. Intersect that with a chr21 slice and you get single-digit
hits from a table that is completely fine.

`run` refuses to start if fewer than half the first 2,000 BAM reads are found
in the table, so this cannot pass unnoticed.

## Two defaults worth knowing about

**`gates.min_aln_coverage` is 0.90, not 0.99.** On an untrimmed FLNC BAM the
polyA tail is soft-clipped, so median coverage against the full query length is
0.93 and a 0.99 gate keeps ~6% of reads *while enriching for mono-exonic ones*.
Set 0.99 if you run on the polyA-trimmed BAM. A run that rejects more than half
its reads at the gates stops and prints the observed percentiles.

**`chains.max_chain_lfdr` is `None`.** The chain-level novelty score is computed
and written to `models.tsv`, but it does not remove anything by default — unlike
the junction level, it has no clean statistical null (see
the algorithm validation note). Look at the distribution, then set
it with `--set chains.max_chain_lfdr=0.05`.

## Outputs

| file | contents |
|---|---|
| `{prefix}.gtf` / `.gff` | transcript models, sorted, GENCODE-shaped attributes |
| `{prefix}.read_stat.txt` | one row per read, **keyed by read name** |
| `{prefix}.group.txt` | model → comma-separated read names |
| `{prefix}.abundance.txt`, `.flnc_count.txt` | pigeon/cupcake-compatible counts, reads *and* molecules |
| `{prefix}.models.tsv` | per model: category, parent transcript, support, posterior, polyA evidence, flags |
| `{prefix}.junctions.tsv.gz` | every junction with support, motif, novelty score, keep/drop reason |
| `{prefix}.group_chains.tsv` | the distinct chains inside each group, with counts |
| `{prefix}_umi_corrected_isoform_matrix.mtx` + `_cells.txt` / `_isoforms.txt` / `_isoform_metadata.txt` | cell × transcript molecule counts |
| `{prefix}.qc.json` / `.qc.md` | the validation invariants |

## Naming

Models are named against the reference, so identity survives downstream:

```
ENST00000229239.10                 exact intron chain + matching ends
ENST00000229239.10|alt3end1        same chain, 3' end of a DIFFERENT isoform
ENST00000229239.10|novel3end1      same chain, an unannotated (supported) polyA site
ENST00000229239.10|alt5end1        same chain, 5' end at another isoform's TSS
ENST00000229239.10|ext5end1        same chain, a start beyond every annotated one
ENST00000229239.10|IR1             that transcript with an intron read through
ENST00000229239.10|readthrough1    3' end has run past the gene into the next
ENST00000229239.10|NIC_1           novel combination of annotated junctions
ENST00000229239.10|NNC_1           contains a novel junction
ENST00000229239.10|3UTRfrag1       mono-exonic 3'UTR fragment of that transcript
GAPDH|APA1                         mono-exonic model at an internal polyA site
NOVELMONO_chr1_3                   intergenic mono-exonic model
```

`alt` and `novel` are kept apart deliberately. An end that coincides with the
annotated end of another isoform of the same gene, and an end that matches
nothing annotated and rests only on read mass plus polyA evidence, are separated
by a large confidence gap — 35% vs 72% coincidence on this data — and one word
for both would hide it.

Novel models are named against their **closest** reference transcript, not just
their gene. `GAPDH|NIC_3` says which locus; `ENST00000229239.10|NIC_3` plus the
`structural_diff` column (`+1jxn;-2jxn`) says what the model actually is.

## Terminal ends, not UTRs

An identical intron chain guarantees that a difference sits in the **terminal
exons**, because every internal boundary is defined by the flanking junctions.
It does not guarantee the difference is *untranslated*: a model starting
downstream of the start codon has a different coding sequence, not a shorter
5'UTR, and a non-coding reference transcript has no UTR at all. The categories
are therefore stated as terminal-end differences, with `is_utr_only` reported
separately.

Two further things an identical chain does not guarantee:

* **Exon contents.** A short intron the aligner spelled as a CIGAR `D` rather
  than an `N` leaves the chain untouched. Those gaps are collected during the
  scan; ones that are annotated, or that recur with a canonical motif, are
  promoted to real junctions *before* the census so they face the same curation
  as everything else. The rest flag their reads and are counted.
* **Direction.** The four end cases need four different evidence bars. A 5'
  extension is interesting and a 5' truncation is degradation; a 3' extension is
  usually a real distal site and a 3' *truncation* is the internal-priming
  signature that needs the polyA test.

A `PB.x.y` id is always emitted alongside (`pb_id` attribute), or as the primary
id with `--set output.id_style=pb`.

## What it does

See [the algorithm scheme and validation](SCHEME_VALIDATION.md) for the audited
overview. In short:

1. **Junction catalogue first.** 83.5% of distinct junctions in this data are
   novel but carry only ~1.2% of junction read mass. Junctions are snapped onto
   the annotated catalogue *first*, then curated on molecule support, splice
   motif, RT-switch direct repeats and a calibrated local FDR. Reads carrying a
   rejected junction are set aside, not promoted into novel isoforms.
2. **Chains second.** Grouped on the curated chain. A chain that matches an
   annotated intron chain gets that transcript's identity.
3. **Suffix resolution.** A non-annotated chain that is a suffix of an annotated
   one is 5' truncation of that transcript — regardless of support ratio, since
   the truncated form routinely outnumbers the intact one — unless its 5' end
   sits at an annotated TSS. Neither side annotated: decided by support
   direction. **The representative is always the maximal chain.**
4. **Ends.** 3' ends are peak-clustered (the anchored end under oligo-dT) and
   compared against **every transcript sharing the chain**, then against every
   annotated end of the gene — not against the MANE transcript, which would turn
   every other isoform's genuine end into a discovery. A peak matching nothing
   annotated becomes a `novel3end` model only with polyA evidence. The model's
   5' end is a high percentile of its members, never the minimum; a 5'-truncated
   pile is folded back onto its parent and counted in `n_5p_truncated_reads`
   rather than emitted, and a 5' *extension* that crosses an annotated junction
   is intron retention, not a new start site.
4b. **Nearest reference transcript.** Every model gets the transcript of its gene
   it is structurally closest to, plus a structural diff and `n_equally_close`.
   This is annotation only: it never feeds back into merging, thresholding or
   the novelty calibration, which is fitted against annotation and would
   otherwise become circular. Its one classification effect is intron retention,
   which an exact-chain lookup has no way to express.
5. **Mono-exonic reads run in a separate track** and never touch the suffix
   logic. Last-exon/3'UTR fragments are trusted; gene-body-internal and
   intergenic ones need an alternative polyA site.
6. **Molecules, not reads.** With the barcode/UMI table, every threshold is on
   UMI-collapsed molecules and cell diversity. 30 reads from one molecule is a
   PCR artefact, not support.

## Calibrated novelty

Instead of "a novel event must appear at least N times" — which is far too
permissive at GAPDH-depth loci and far too strict at shallow ones — novel
junctions and chains get a feature vector, a ridge-logistic score anchored on
annotated events (positives) versus a decoy set (wrong-strand and non-canonical
motifs, RT-switch repeats), and a two-component mixture that converts the score
into a **local false discovery rate**. `--set junctions.max_junction_lfdr=0.01`
is a calibrated dial; `N=10` is not.

The model degrades gracefully: too few anchors, or `scoring.enabled=false`, and
it returns `NaN` and the explicit thresholds take over. It never silently
disables a threshold you set.

### Short reads as an anchor, never a filter

```bash
flightcollapse run -b long.bam -g gencode.v44.annotation.gtf -f GRCh38.fa \
  --short-read-sj repA_SJ.out.tab repB_SJ.out.tab repC_SJ.out.tab
```

Anchoring positives on *annotation* means the score can only learn to recognise
things already in GENCODE — the wrong prior for finding novelty. STAR's `sjdb`
flag supplies what annotation cannot: junctions an orthogonal library saw and the
reference does not contain. Those join the positive set and leave the decoy set;
**half are held out** so the recall measured on them is honest. Replicate files
are summed. Nothing is ever dropped for lacking short-read support —
`junctions.short_read_rescue` can make support *save* a junction, but its absence
never condemns one.

The coordinate convention is checked rather than trusted. STAR writes 1-based
inclusive intron bounds; a one-base error returns "unsupported" for everything,
which looks exactly like a library that confirms nothing. `run` compares STAR's
own annotated flag against the reference index and refuses to start if they
disagree, naming the offset or the contig mismatch.

To score a run you have already finished:

```bash
flightcollapse sj repA_SJ.out.tab repB_SJ.out.tab \
  -g gencode.v44.annotation.gtf -j out/sample.junctions.tsv.gz
```

## Validation

Every run asserts, and by default fails on:

1. read mass on mono-exon models ≈ fraction of genuinely unspliced aligned reads
2. `model_is_min` is not ≈ 1.0 (the model is not always the shortest member)
3. median `log2(read span / model span)` ≈ 0
4. reads whose structure a merge changed, reported explicitly

Use `--no-strict` to report instead of exiting non-zero.

## v0.5.0: contig parallelism

`run()` distributes contigs over processes. The reason it is safe is worth
stating rather than assuming: a contig's result depends on no other contig, and
all three calibrated models — `score_junctions`, `score_chains`,
`score_terminals` — were already fitted **per contig**, so nothing about the
calibration changes when contigs move between processes.

What *is* global is small, and all of it is done in the parent:

* **Model ids.** Genes do not cross contigs, so the per-gene counters would be
  safe in a worker — but the collision set and the `PB.N` gene numbering are
  global, and `gene_name` repeats across contigs (`Y_RNA` exists on chr9, chr10
  and chr16 as three different gene ids; this is the collision class that has
  bitten this codebase three times). Ids are assigned in the parent, contig by
  contig in the configured order, which is exactly what the serial path did.
* **`read_stat` / `group`.** Workers write per-contig shards keyed on a
  *contig-local model index*; the parent substitutes ids while stitching, in
  contig order. Row order within a contig is the BAM's, as before. This costs
  one extra pass over a text file.
* **The molecule matrix**, accumulated from per-contig count shards (`.npz`
  rather than pickled dicts, which is the difference between seconds and
  minutes on a real catalogue).

`fork` is required and is not incidental: it gives every worker a
copy-on-write view of the parsed `ReferenceIndex`. Under `spawn` each worker
would re-read the GTF, which costs more than the parallelism saves, so a
platform without `fork` falls back to serial with a warning rather than
silently paying that cost.

pysam and pyBigWig handles are opened lazily per process and re-opened whenever
the pid changes. A file descriptor shared across a fork does not raise — it
interleaves reads and returns wrong data — so `_Shared.handles()` never hands a
child a handle it inherited.

### Memory, not cores

Peak RSS is roughly `workers × per-contig peak`: one contig's read arrays plus
that worker's own cached chromosome sequence (~250 MB for chr1 as a Python
string). `Config.worker_count()` clamps the requested count against
`MemAvailable`, and the run logs it whenever the clamp bites — including, and
especially, when it bites all the way down to 1. A `-j 12` that silently ran
serially would waste an afternoon before anyone thought to check.

Wall-clock is floored by the slowest single contig. With ~24 contigs the
speedup saturates well below the core count; contigs are dispatched
longest-first (`parallel.longest_first`) so the tail is not a whole chr1 long.

### Determinism

`tests/test_parallel.py` runs a four-contig fixture at 1 and 4 workers and
asserts that the GTF, GFF, models table, abundance, flnc counts, read_stat,
group, group_chains, second-pass audit and the whole molecule matrix are
identical, and that the QC verdict matches. This is an equality assertion, not
a smoke test, because every comparison in this project — three-sample
reproducibility, SIRV precision, arm-versus-arm catalogues — assumes two runs
differing only in machine settings are the same run.

## v0.5.0: the second pass

Runs after model building, inside the contig (it needs `m.reads`, `reads` and
the genome handle, all of which exist only there, and gene-level questions are
contig-local because genes do not cross contigs). It may keep, remove or
replace a model; it never builds a structure.

### Terminal evidence

`secondpass.measure_terminal_evidence` runs *before* `score_terminals`, so the
calibrated terminal model sees the new features too rather than the second pass
being their only consumer — that scoring level already has the right positives
(annotated TTS / atlas hits) and the right decoys (A-rich downstream, nothing
annotated), and these are the features its null was missing.

New columns in `models.tsv`:

| column | meaning |
|--------|---------|
| `genomic_a_run_3p` | genomic A-run starting at the model's 3' end |
| `max_a_run_downstream` | longest A-run in the first `a_run_scan` bases |
| `median_tail_excess` | median over member reads of (read's terminal A-run − genomic A-run **at that read's own 3' end**) |
| `n_nontemplated_reads`, `frac_reads_nontemplated` | how many cleared `min_tail_excess` |
| `dse_gu_frac`, `dse_a_frac`, `dse_gu_minus_a` | composition of the downstream element |
| `upstream_u_frac` | U-richness upstream of the cleavage site |

The genomic A-run is memoised per *position*, not per model: a peak is 4–35 bp
wide and the genomic context can change inside it, so the exact test costs the
same as the approximate one.

### Thresholds are starting points

`min_tail_excess = 8`, `min_nontemplated_frac = 0.5`, `min_dse_gu_minus_a =
0.10` and `max_dispersion = 12` are defensible defaults, not measurements on
your data. The columns above are written for every model, including FSM ones,
precisely so the distributions can be looked at first. The recommended
sequence on a new dataset:

1. run with `--no-second-pass` and plot `median_tail_excess` and
   `frac_reads_nontemplated` for FSM and `end3_annotated` models — those are
   the closest thing to a positive control the sample contains;
2. do the same for `dse_gu_minus_a`;
3. set the thresholds where the two populations separate, then run with the
   second pass on and compare the catalogues.

`min_dse_gu_minus_a` is a *contrast* rather than an absolute G+U cut for a
concrete reason: the window is ~25 bases, so an absolute threshold at 0.45 sits
one base away from 0.40 and neutral sequence crosses it by chance about half
the time. The contrast asks the question the tier actually cares about — is the
downstream A-rich? — where neutral sequence sits near +0.25 and an
internal-priming site goes sharply negative.

### What the second pass cannot do

* It cannot rescue a model the first pass never emitted. If a genuine
  alternative 3' end failed `require_polya_evidence_for_utr_variant`, it is
  already an `end3_unresolved` (not emitted by default) before this stage runs.
* The `nontemplated` tier needs untrimmed 3' soft clips. On a polyA-trimmed BAM
  every excess is 0 and the tier has no opinion — `live_end3_tiers` in the QC
  report says so, and `gates.measure_polya_tail=false` is rejected at config
  time rather than silently rejecting every model.
* Class D removal has a real cost: a short last-exon-only model in a multi-exon
  gene *can* be a genuine proximal APA event. `secondpass.remove_classes` is a
  list for that reason — `["C", "E"]` keeps class D while still removing the
  intron-spanning and out-of-structure cases.

## v0.6.0: the final consolidation

Two rules run per contig as step 8b, after the second pass and before any read
bookkeeping. The placement is the whole point: `group.txt`, `read_stat`, the
abundance tables and the molecule matrix are all produced in step 9, so a model
removed here was never in them, and nothing downstream has to be reconciled
afterwards. The order inside the step is fixed — 3' consolidation first, then
the fragment filter — so the support ratio the fragment filter tests is a
consolidated count rather than a count that is about to change.

Support means `n_mols` where it is positive and `n_reads` otherwise, in both
rules and in the audit.

### Why the representative is not the longest model

The obvious choice for a 3' collapse is to keep the longest end, on the theory
that it is the least truncated. On SIRV Set 4 that choice loses exact-match
recall as the tolerance grows, because a genuine proximal cleavage site gets
absorbed into a longer neighbour and the model that matched the reference stops
existing. Ordering by annotated end, then support, then length holds recall at
its N = 0 value in both the annotated and the annotation-free arm while still
removing the duplicates. The annotated flag comes from `end_is_annotated`,
which the terminal-ends stage writes only when a reference was supplied; with
no reference it is absent everywhere and the rule degrades cleanly to support
order.

### Why the clusters are anchored

Within a chain, a model is absorbed when both its ends are within tolerance *of
the representative* — never of whichever model was absorbed last. Chaining
pairwise lets a dense run of 3' ends 400 bp apart merge things kilobases apart
at a 500 bp setting, and the span of a cluster then grows with how many models
the locus happens to carry rather than with the parameter. This is the same
reason `compare._match_block` rejected single linkage. The test fixture pins
it: ends at 0, 400 and 800 with a 500 bp tolerance give two clusters, not one.

### Why the 5' window is not a parameter

`end5_window` is a guard, fixed at 100 bp. Widening it lets the rule merge
models that differ by an alternative TSS, which is a different transcript and
not terminal-end redundancy. The user-facing tolerance is on the 3' end only
because that is the end an oligo-dT protocol actually measures; the 5' end is a
percentile over a degraded read pile and has no business driving a merge.

### Why the fragment filter needs the terminal check

The contiguous sub-chain test alone would remove alternative-first-exon
isoforms, which share every junction of a longer model and differ only in where
they start. The terminal check asks where the first exon *begins*: a 5'
truncation starts inside the container's matching exon, while an alternative
first exon starts inside the upstream intron, before that exon. Dropping the
check costs recall in the annotation-free arm. `subchain_terminal_slack`
(50 bp) is the tolerance on "inside", not a licence to ignore the question.

Decisions are taken in one pass against the pre-filter supports. Iterating
would let a model removed early change the ratio for a model considered later,
which makes the result depend on iteration order. When a container is itself
removed, the fragment resolves to the final surviving container, so reads never
land on a model that is not in the catalogue.

### Scale

Both rules are indexed. The fragment search looks a candidate's first intron up
in a positional index of every intron on the contig and checks only the models
that share that coordinate; mono-exonic candidates go through an exon index
with a running-maximum end so the scan can stop early. The 3' rule bisects a
position-sorted view of each chain partition rather than scanning it.

The all-pairs form of the fragment filter is what a direct reading of the rule
produces, and it is 88x slower at 1,800 models on one synthetic contig — a gap
that widens quadratically. At 9,000 models per contig the indexed form takes
0.1 s; all-pairs takes about 25 s, which across a whole genome is 20 minutes
per sample spent on a stage that is supposed to tidy the output. The indexed
form is checked against the literal definition on random catalogues in
`tests/test_posthoc.py::test_indexed_filter_matches_the_brute_force_definition`,
because an optimisation that quietly decides something different from the rule
it implements is worse than the slow version.

### What the consolidation cannot do

* It cannot merge two different splice structures. Partitioning is on the exact
  intron chain, so the 3' rule has no way to reach across one.
* It cannot recover a model the earlier stages dropped; it only removes.
* It overlaps with `monoexon.demote_contained_fragments` and with the second
  pass's class C/D removal, which reach the same mono-exonic fragments earlier.
  On a default run there is usually little left for it to do there, and its
  work is on the spliced 5'-truncation class instead. Tests that are about the
  first-pass mono-exon track therefore pin `posthoc.enabled = False`, the same
  way they already pin the second pass off.
* `--no-posthoc` reproduces 0.5.x output exactly. Verified by running 0.5.1 and
  0.6.0 `--no-posthoc` over the same simulated dataset: all fifteen output
  files match, twelve byte for byte and the three gzipped ones byte for byte
  after decompression (a `.gz` header carries the wall clock of the run).
