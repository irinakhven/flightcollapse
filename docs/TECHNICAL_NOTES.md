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
