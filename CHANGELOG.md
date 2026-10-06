# Changelog

Versions 0.2.0 through 0.4.6 were developed and run internally but never
published, and no changelog was kept for them. What they changed is described,
version by version, in the "New in" sections of the README — external 5'-end
evidence (0.2.0), the terminal end and sample-type profiles (0.3.0), and the
gate rework behind 0.4.6. The entries below are the ones that were written at
the time.

---

# flightcollapse 0.6.3

## Changed

- **`posthoc.end3_tolerance` 500 -> 200** and **`posthoc.end5_window` 100 -> -1
  (off)**. 0.6.0 shipped 500 on the strength of SIRV Set 4 on BD144_kin alone,
  with a note in `PosthocParams` that BD144_HM and BD144_ND would confirm it
  before release. They did not.

  Three-way reproducibility across the three HEK libraries, measured as the
  share of assigned reads sitting on a transcript structure found in all three
  — the only measure a consolidation step cannot bias, because merging two
  reproducible structures into one always lowers a count-based share:

  | 3' tolerance | models | reproducible read mass | at >=50 reads, >=100 barcodes | positive cells |
  |---|---|---|---|---|
  | **200 bp** | 499,553 -> 418,097 (-16.3%) | **-0.32 pp** | **+0.11 pp** | 34 of 81 |
  | 500 bp | 499,553 -> 392,017 (-21.5%) | -1.11 pp | -0.25 pp | 17 of 81 |

  The extra 300 bp buys its further 26,080 merges with reproducible read mass.
  At the splice-chain level the two arms are indistinguishable (122,270 chains
  against 122,174), which is what says the 3' tolerance is doing 3' work only
  and that the chain-level cost belongs entirely to the fragment filter.

  Turning the 5' guard off costs almost nothing (24-5,664 further merges per
  arm, at most 0.03 pp of mass) and makes the asymmetry explicit: the 3' end
  gets a tolerance, the 5' end gets the structural sub-chain filter and nothing
  else.

  `--end3-tolerance 500 --set posthoc.end5_window=100` restores 0.6.0-0.6.2
  behaviour exactly.

- **`--no-subchain-filter` on its own is now rejected.** With the 5' guard off
  by default, turning the filter off as well would leave nothing looking at the
  5' end at all. `validate()` says so and names the fix. Pass
  `--set posthoc.end5_window=100` alongside it.

- SIRV Set 4 is not retracted, and both benchmarks are kept in the
  `PosthocParams` docstring. The fragment filter is still worth 26 points of
  annotation-free precision at no cost in recall, and 500 still loses no
  exact-match recall on that panel. The panel has 69 transcripts with ground
  truth and almost no 3'-end diversity; the HEK libraries have 418,097 models
  and no ground truth. They measure different things.

- Two tests asserted the default's *value* where they meant to test something
  else. `test_pipeline_writes_the_posthoc_audit` now compares against
  `Config().posthoc.end3_tolerance`, and the rule-ordering test pins the 500 bp
  its fixture needs rather than inheriting whatever the default happens to be.

## Added

- **The logic sheet is current again.** `docs/flightcollapse_logic_v0.6.3.html`
  (and the PDF beside it) carries the 0.5.1 sheet forward with a third band for
  the consolidation: one card per decision, each drawn as the exon geometry the
  decision is about — same chain with ends 150 bp apart against 400 bp apart,
  anchored absorption against single-linkage chaining, a 5' truncation against
  an alternative first exon that starts inside an intron, FSM protection, the
  2x ratio held at 1.6x, and what turning the 5' guard off costs and buys. Built
  by `docs/make_logic_sheet.py` from the 0.5.1 sheet, which is kept beside it.
- `docs/make_schematic.py` and `flightcollapse_algorithm_distilled.svg` — the
  same pipeline in five steps, for a figure panel. The README's previous diagram
  was three versions out of date; generating these from scripts is what stops
  that recurring.

## Independent check

pigeon (SQANTI3) classifies every model without knowing anything about these
rules, and the classes the sub-chain filter empties are exactly the fragment
classes. In BD144_kin: ISM 3' fragments -60%, ISM 5' fragments -44%, against
-4% for exact reference matches — while the 3' tolerance's own signature shows
up as -17% in *FSM with an alternative 3' end* and only -7% in *FSM with an
alternative 5' end*. A reference-free rule sorting itself along a reference it
cannot see is the result worth having.

---

# flightcollapse 0.6.2

## Added

- `posthoc.end5_window` accepts a **negative** value, which removes the 5' guard
  from the 3' rule entirely: same contig, same strand, same exact intron chain
  and 3' ends within `end3_tolerance` becomes the whole of rule 1, and the 5'
  end is left to the sub-chain filter. CLI: `--end5-window BP`.

  This is a real choice rather than a loosening, and the docstring says so. The
  sub-chain filter only ever compares a model whose chain is a STRICT sub-chain
  of another, so it never sees two models that share a chain exactly — with the
  guard off, two identical chains differing only in where transcription started
  merge, alternative TSS included. With the guard on, a 5'-truncated copy of a
  model never merges either, because it shares that model's chain exactly and so
  is invisible to the filter as well. The guard decides which of those two you
  would rather have.

  `validate()` rejects a negative `end5_window` together with
  `subchain_filter = false`: nothing would be looking at the 5' end at all.

Defaults are unchanged — `end5_window` is still 100.

Suite: 434 passed, 10 skipped.

---

# flightcollapse 0.6.1

## Fixed

- `posthoc.filter_subchains` indexed candidate containers on `(strand, position
  bucket)` without the contig. `_is_fragment_of` rejects a cross-contig pair, so
  the definition was always right; the index was not. The pipeline calls the
  filter once per contig, so **no 0.6.0 run is affected** — every model it saw
  was already on one contig. It bites the moment the function is called on a
  whole catalogue at once, which is exactly what the 3' tolerance sweep does:
  two contigs share coordinates, and a chr2 model would have been offered as the
  container for an identically placed chr1 one. Now keyed on
  `(contig, strand, bucket)`, with the exon index likewise.
  `tests/test_posthoc.py::test_a_fragment_never_crosses_a_contig` pins it, and
  checks the same pair on one contig is still removed so the guard is not a
  blanket refusal.

Suite: 432 passed, 10 skipped.

## Not changed

Rules, defaults, outputs and the audit are identical to 0.6.0. A 0.6.0 run and a
0.6.1 run of the same sample produce the same catalogue.

---

# flightcollapse 0.6.0

Two post-hoc rules, from the SIRV Set 4 benchmark on BD144_kin. Nothing about
how models are *built* changed: same reads, same junction curation, same suffix
resolution, same 3'-peak clustering, same second pass. `--no-posthoc` reproduces
0.5.1 output exactly.

## Added

`src/flightcollapse/posthoc.py` — final per-contig consolidation, run as step 8b
after the second pass and **before** read bookkeeping, so `group.txt`,
`read_stat`, the abundance tables and the molecule matrix all describe the
catalogue that is actually written. The order inside the step is fixed: 3'
consolidation, then the fragment filter, so the ratio in the second rule is
tested on consolidated counts. Support means `n_mols` where positive, else
`n_reads`.

**Rule 1 — same-chain 3' consolidation.** Models are partitioned by contig,
strand and *exact* intron chain; nothing is ever compared across a partition, so
the rule cannot merge two splice structures. Within a partition the
representative is chosen by annotated 3' end > support > transcript length >
start, and absorbs every model whose 5' end is within `end5_window` and whose 3'
end is within `end3_tolerance` **of it**. Anchoring to the representative rather
than chaining pairwise bounds every cluster at the tolerance. Absorbed models'
reads move to the representative and molecules are recounted from them.
INFERRED models are untouched.

**Rule 2 — reference-free sub-chain (fragment) filter.** A model is removed when
its ordered intron chain is a contiguous sub-chain of a longer model on the same
contig and strand (junctions within `subchain_fuzzy`), its terminal exons lie
inside the matching exons of that model (within `subchain_terminal_slack`, so an
alternative first exon starting in an upstream intron is kept), and that model
carries at least `subchain_ratio` x its support. A mono-exonic model is a
fragment when it lies inside one exon of such a model. `protected_categories`
(FSM) is never removed; INFERRED is neither removed nor used as a container.
Decisions are taken in one pass on pre-filter supports; a container that is
itself removed resolves to the final survivor. Reads fold into the surviving
container, or become unassigned and are counted when `fold_reads` is off.

Also: `{prefix}.posthoc.tsv` audit (one row per removal, with the index *and*
id of the model it went into), a `posthoc` block in `qc.json`, posthoc counters
in `read_accounting`, and `posthoc_end3_absorbed` / `posthoc_fragments_absorbed`
in surviving models' evidence.

## Parameters (config section `posthoc`)

| key | default | notes |
|---|---|---|
| `enabled` | `true` | `--no-posthoc` turns both rules off |
| `end3_tolerance` | **500** | `--end3-tolerance BP`; 0 = identical 3' ends only; < 0 disables rule 1 |
| `end5_window` | 100 | a guard, not a knob; must be >= 0 |
| `subchain_filter` | `true` | `--no-subchain-filter` |
| `subchain_ratio` | 2.0 | `--subchain-ratio X`; validated >= 1 |
| `subchain_fuzzy` | 5 | validated >= 0 |
| `subchain_terminal_slack` | 50 | validated >= 0 |
| `protected_categories` | `["FSM"]` | |
| `fold_reads` | `true` | |
| `write_audit` | `true` | |

## Evidence

SIRV Set 4, BD144_kin, 69 conventional SIRVs. Exact match = junctions +/- 5 bp,
ends +/- 50 bp. Values are precision / recall / F1, both rules applied.

| arm | read filter | N = 0 | N = 200 | **N = 500** | N = 1000 |
|---|---|---|---|---|---|
| no annotation | none | 0.54 / 0.84 / 0.66 | 0.73 / 0.84 / 0.78 | **0.82 / 0.84 / 0.83** | 0.86 / 0.83 / 0.84 |
| no annotation | >= 50 reads | 0.76 / 0.84 / 0.80 | 0.83 / 0.84 / 0.83 | **0.88 / 0.84 / 0.86** | 0.89 / 0.83 / 0.86 |
| reference | none | 0.74 / 0.87 / 0.80 | 0.80 / 0.87 / 0.83 | **0.86 / 0.87 / 0.86** | 0.92 / 0.87 / 0.90 |
| reference | >= 50 reads | 0.87 / 0.87 / 0.87 | 0.92 / 0.87 / 0.90 | **0.94 / 0.87 / 0.90** | 0.94 / 0.87 / 0.90 |

The same run with neither rule scores 0.28 / 0.84 / 0.42 without annotation, so
the fragment filter alone (N = 0) is worth 26 points of precision at no cost in
recall. N = 500 is the largest tolerance that loses no recall in either arm; at
N = 1000 the annotation-free arm gives up one true model.

Other annotation-free callers on the same data: IsoQuant 0.89 / 0.67 / 0.76,
FLAIR transcriptome mode 0.53 / 0.86 / 0.65 (0.69 / 0.86 / 0.77 at >= 50 reads),
Bambu 1.00 / 0.25 / 0.40.

**BD144_HM and BD144_ND have not confirmed these defaults yet.** The release
ships 500 on BD144_kin alone; the SIRV benchmark reports N in {0, 50, 100, 200,
300, 500, 750, 1000, 2000} for all three, and the default is a one-line change
in `PosthocParams` if the other two disagree.

## Changed from the handed-over reference implementation

- **`end3_tolerance` default is 500, not 200**, with the evidence table above in
  the config docstring, the CLI help, the README and the technical notes.
- **Both rules are indexed instead of quadratic.** The fragment search looked
  every model up against every other model on the contig, which is 88x slower at
  1,800 models on one synthetic contig and widens quadratically: about 25 s per
  contig at 9,000 models, roughly 20 minutes per sample across a genome, on a
  stage meant to tidy the output. It now looks a candidate's first intron up in
  a positional index; mono-exonic candidates go through an exon index with a
  running-maximum end. The 3' rule bisects a position-sorted view of each
  partition instead of scanning it. The per-survivor fragment tally was also
  quadratic and is now a single pass.
- **Read tallies survive a recount.** `_move_reads` returns the count it folded
  in from a source with no read-index array, because `_recount` sets `n_reads`
  from the index array and would otherwise discard it.
- **The audit's `into_model_id` is always present**, defaulted at build time
  rather than left to the writer.
- **Container choice is explicitly deterministic** (best support, then leftmost,
  then longest) rather than depending on list order.
- `_genomic_introns` is computed once per model per rule instead of per
  comparison, and its docstring says why `TranscriptModel.chain` is not used:
  the chain runs in transcript orientation and so is reversed on the minus
  strand, which would silently break the sub-chain test there.

## Changed

- Three first-pass tests (`test_pipeline`, `test_promotion`,
  `test_terminal_score_independence`) pin `posthoc.enabled = False`, as they
  already pinned the second pass. They test first-pass mono-exon behaviour, and
  the new filter correctly folds that simulated fragment into its container.
  Only the first of the three was pinned in the handed-over tree, which is why
  it reported 419 passing with 2 failures outstanding.

## Tests

`tests/test_posthoc.py` — 29 tests. Beyond the handoff list:

- the indexed filter is checked against the literal all-pairs definition of the
  rule on 40 randomised catalogues built from truncations, alternative first
  exons, exon skips and mono-exons;
- the two rules are shown to run in the documented order, by a fixture where the
  fragment survives if they do not;
- `--no-posthoc` is shown byte-identical to physically bypassing step 8b,
  comparing decompressed content so a `.gz` header timestamp cannot make it pass
  or fail by accident;
- every downstream table is shown to name only surviving models, and every audit
  row to point at a model that is in the catalogue.

Suite: 431 passed, 10 skipped.

Cross-version check, outside the suite: 0.5.1 and 0.6.0 `--no-posthoc` run over
the same simulated dataset produce all fifteen output files identical — twelve
byte for byte, and the three gzipped ones byte for byte after decompression,
the containers differing only in the wall clock gzip writes into its header.

---

# flightcollapse 0.5.1 — a killed worker no longer hangs the run

`multiprocessing.Pool` has no abrupt-death detection. When a worker is SIGKILLed — the
OOM killer is how that happens here — the pool quietly forks a replacement, but the task
it was holding never produces a result and `imap_unordered` waits for it forever.

Observed in production: BD144_kin sat at 23 of 24 contigs for two and a half hours with
every worker idle and nothing in the log, because one worker running the smallest contig
in the genome had been killed. A run that hangs silently after 95% of the work is the
worst failure mode this tool has.

The parent now announces which pid owns which contig, notices when one disappears, and
re-runs that contig itself — serially, in the parent, after the pool has closed, because
if the cause was memory then re-dispatching into a full pool would just kill it again.
The loss and the recovery are both logged by name.

`tests/test_parallel.py` gained a test that SIGKILLs a worker deterministically and
asserts no contig is lost. Full suite: 402 passed, 10 skipped.

---

---

# flightcollapse 0.5.0 — the second pass

A filtering stage that runs **after** the catalogue is built and may only do three
things: keep a model, remove it, or replace it with an explicitly inferred one. It
never builds a structure from scratch.

The separation is deliberate. A filter folded into model construction cannot be turned
off, cannot be compared against its own absence, and cannot report what it cost.

#### Rule 1 — mono-exonic models judged against their gene's exon structure

The mono-exon track asks "terminal exon or not". It never asks whether a single-exon
transcript could exist in that gene at all. Five classes, decided from the annotation
alone:

| class | gene / model relationship | verdict |
|---|---|---|
| A | gene is mono-exonic only | keep |
| B | gene has both mono- and multi-exon isoforms | keep |
| C | multi-exon gene; model spans an intron | replace |
| D | multi-exon gene; model lies inside one exon | replace |
| E | multi-exon gene; anything else | replace |
| X | no gene (intergenic) — out of scope for A–E | policy, default keep |

"Spans an intron" uses **constitutively** intronic sequence — inside the gene span and
exonic in no isoform. Using any one transcript's introns would call an ordinary
alternative-exon model a read-through, because one isoform's exon is another's intron
all over the genome. Default `min_intron_overlap = 10` bp.

A replaced model does not vanish. The gene's most-supported **emitted** multi-exon model
is recorded in its place, falling back to the reference canonical (MANE_Select >
Ensembl_canonical > basic > longest) when the gene has no emitted spliced model. The
removed model's reads move onto it, so read accounting stays exact. The replacement is
marked `INFERRED` as a category, as a flag, and as an `inferred=1` column — three times,
because the one failure this must not have is an inferred model read as an observed one.

Measured across the nine runs: 167,467 mono-exon models replaced, 93.7% from the
sample's own data and 6.3% from the reference canonical, creating 7,715 distinct
INFERRED models. On BD144_kin the class split is A 3,006 · B 2,283 · C 10,734 ·
D 4,959 · E 5,548 · X 4,816 — roughly four in five in-gene mono-exon models are C/D/E.

#### Rule 2 — `end3_novel` models must earn their 3' end

On SIRV Set 4 this category was wrong 12/12 in the reference arm and 11/11 in the masked
arm. The mechanism is internal priming: oligo-dT anneals to a genomic A-stretch inside
the transcript and the resulting pile accumulates support exactly like a real cleavage
site.

Three things that **cannot** separate the two populations, all checked:

- a read floor — false ends reach p90 548 reads while true ones go down to 5;
- `tail_molecule_frac` — 0.96 against 1.00 on BD176c, because it scores an A-run without
  asking what the genome says underneath it;
- **cross-sample reproducibility** — an internal-priming site recurs in every library
  made from the same genome, so a model called in all three HEK samples is not thereby
  real. This one matters for the concordance work: reproducibility is a stand-in for
  precision everywhere except here.

Six tiers, strongest first; the first that speaks decides, and the tier is written on
every judged model so a downstream filter can pick its own operating point without
re-running anything.

| tier | test | share of survivors |
|---|---|---|
| `atlas` | end sits on a catalogued polyA site | 63.5% |
| `nontemplated` | reads' own tails exceed the genome's A-run | 35.2% |
| `composition` | hexamer + U/GU-rich DSE + tight peak | 0.1% |
| `short_read` | short-read coverage demonstrably drops across the end | 0.0% |
| `calibrated` | terminal local FDR clears a ceiling | off by default |
| `gene_dominant` | sole or best-supported model of its gene | 1.2% |

**`nontemplated` is the only per-molecule test and the only genuinely new information.**
A real poly(A) tail is not in the genome: the A-run in a read's 3' soft clip runs far
past the genomic A-run at that coordinate, while an internally primed read has no excess
because the aligner consumed the templated A's. Everything needed was already measured —
`reads.tail_len` per read since 0.1.17, `genome.fetch` since the beginning — and nothing
compared the two. Defaults: excess ≥ 8 bp on ≥ 50% of the model's reads, scanning 60 bp
downstream.

This tier is why the collapse must read the **untrimmed** dedup BAM. A trimmed read has
no terminal A-run and `gates.measure_polya_tail` would return zero for every read, so the
tier would silently have no opinion.

`composition` scores a **contrast**, `(G+U) − A ≥ 0.10`, not an absolute G+U cut: the
window is ~25 bases, so an absolute 0.45 threshold is one base from 0.40 and neutral
sequence crosses it by chance about half the time. Neutral sequence sits near +0.25;
an internal-priming site goes sharply negative. A single `perc_a_downstream` number
cannot express this — on BD176c it moves only 20 → 25 between FSM and end3_novel.

`gene_dominant` is last on purpose. Used first it selects against the biology: a genuine
proximal polyA site usually carries a minority of a gene's molecules, while a strong
internal-priming site in a well-expressed gene is exactly what becomes dominant.

Measured: 239,695 end3_novel models judged across the nine runs, 219,368 kept, 7.3–10.0%
dropped per sample. Dropped models' reads fold onto the best compatible survivor (same
intron chain first, then the gene's best model); zero reads were left unassigned.

**Net effect of the second pass on catalogue size: −8.5% to −13.0% per sample.**

#### Contig-parallel execution

`-j / --workers`. A contig's result depends on no other contig, and all three calibrated
models (`score_junctions`, `score_chains`, `score_terminals`) were already fitted per
contig — so distributing contigs across processes changes no calibration, no threshold
and no decision.

What is genuinely shared is small and all of it lives in the parent: model IDs and the
`PB.N` gene numbering (assigned in contig order), `read_stat` and `group` (stitched from
per-contig shards), and the molecule matrix (accumulated from count shards). `fork` only
— under `spawn` every worker would re-read the GTF, which costs more than the
parallelism saves.

Measured on BD144_kin: 2.24 h of summed contig work, 0.65 h wall at `-j 8` — **3.5×**,
with the floor set by chr1 at 997 s, 43% of wall. Longest-first dispatch. Memory, not
cores, is the binding constraint; `worker_count()` clamps against MemAvailable.

**Byte-identity is an acceptance test, not a nice-to-have.** `tests/test_parallel.py`
asserts a parallel run is byte-identical to a serial one across model IDs, `PB.N`
numbering, `read_stat`, `group` and the molecule matrix. If `-j 4` changed one model ID,
no cross-sample comparison in this project could be trusted again.

#### New outputs and columns

- `<prefix>.secondpass.tsv` — one row per second-pass decision, with its reason.
- `<prefix>.monoexon_rejected.tsv` — mono-exon models the first pass refused, by reason.
- `models.tsv` gains `monoexon_gene_class`, `inferred`, `inferred_source`,
  `inferred_transcript`, `inferred_from_monoexon`, `end3_second_pass_tier`,
  `genomic_a_run_3p`, `max_a_run_downstream`, `median_tail_excess`,
  `n_nontemplated_reads`, `frac_reads_nontemplated`, `dse_gu_frac`, `dse_a_frac`,
  `dse_gu_minus_a`, `upstream_u_frac`, `n_3p_absorbed_reads`.
- `qc.json` gains `second_pass` (class and tier histograms, reads moved),
  `parallel` (workers, contig seconds, slowest contig), and two invariants.
- CLI: `-j/--workers`, `--no-second-pass`.

---

---

## 0.1.15

**0.1.13 was wrong. The regression it "fixed" was not real, and reverting is the
fix.** Re-measured with the 0.1.14 matcher, on the same whole-genome output:

| min molecules | 0.1.6 reproduced | 0.1.10 reproduced | delta |
|---|---|---|---|
| 0 | 58.2% | 54.9% | -3.3 |
| 3 | 64.6% | 60.0% | -4.6 |
| **5** | 66.6% | **69.7%** | **+3.1** |
| **10** | 78.3% | **86.9%** | **+8.6** |
| **20** | 91.5% | **94.4%** | **+2.9** |
| 50 | 96.8% | 97.5% | +0.7 |
| 500 | 99.2% | 99.7% | +0.5 |

0.1.10 is better at every threshold from 5 molecules up, and at >=20 it reaches
94.4% on **fewer** structures than 0.1.6 (64,959 vs 67,178) with a third the
singletons (1,835 vs 3,056). It consolidates; it does not fragment. The
"88.0% -> 73.5%" that prompted 0.1.13 was the broken matcher.

The narrow clustering window is therefore correct and `peak_window` reverts to
inheriting `max_3p_diff` (30). The knob stays, because the two questions really
are distinct -- the default is now set by measurement rather than by argument.

The extra 3' peaks are real, and `median_3p_spread` says so: for structures
called in all three libraries it is **0 bp in every category, in both versions**
-- including `end3_novel`, which is not snapped to annotation at all. Cleavage
in this data reproduces to the base across independent library preparations.

The only place 0.1.10 loses is unfiltered (58.2% -> 54.9%), which is the lowered
support floors doing exactly what they were asked to do: emit low, carry the
lFDR, threshold downstream. The unfiltered percentage is the wrong headline for
a tool configured that way; the curve is the headline.

### Large-block matching

`compare` assigned any structural block over 4,000 models to a **single group**,
silently merging thousands of unrelated models into one "reproduced" structure.
It fired twice on the 0.1.10 run. Mono-exonic models all share the empty chain,
so every one of them on a contig lands in one block -- this was not a rare path.
Replaced with a left-to-right sweep that keeps the one-model-per-sample
guarantee and is linear in the block size.

## 0.1.14

**The confirmatory test I proposed for the 0.1.13 regression was invalid, and
the flaw was in `compare`.** Loosening `--max-3p-diff` from 100 to 300 on the
0.1.10 output produced *more* structures (283,662 -> 308,889) and lower
reproducibility (47.5% -> 40.5%). Merging neighbouring sites cannot create
structures; that was the tool.

The matcher cut each structural block wherever consecutive 3' ends were more
than the tolerance apart, then patched the result by forcing any sample
appearing twice in a group into its **own unique singleton group**. Widening the
window merges two honest pairs into one group and manufactures singletons out of
the surplus. Worse, it means any comparison between two callers that fragment 3'
ends to different degrees was reading that artefact as a difference in the data
-- so the 0.1.6-vs-0.1.10 gap was measured with a contaminated instrument.

Replaced with proper one-per-sample matching: repeatedly take the position
covering the most distinct samples within the tolerance, claim the nearest model
from each, continue. Same greedy-peak idea as the 3'-end clustering in `ends.py`.
Vectorised per sample (`O(S n log n)` a round); the obvious nested-loop version
did not finish a whole-genome comparison at a wide tolerance.

Monotonicity now holds, on the real 0.1.6 data:

| tolerance | structures | in all three | at >=20 molecules |
|---|---|---|---|
| 100 | 192,690 | 58.2% | 91.5% |
| 300 | 188,313 | 61.2% | 92.5% |

Fewer structures and higher agreement as the window widens, which is the only
direction that makes sense. The 0.1.6 baseline moves from 57.0% to 58.2% (88.0%
-> 91.5% at >=20 molecules) purely from fixing the instrument.

### `three_spread`, the number that settles the argument

New per-structure column, and `median_3p_spread` in the summary: how far apart
the samples put the same structure's 3' end.

This separates "the caller resolved a genuinely distinct cleavage site" from
"the caller split one site at a different offset in each library". A real site
reproduces to within cleavage heterogeneity (~20-30 bp); a split does not.

Measured on 0.1.6, for structures called in all three: **median 0 bp, 97.7%
within 30 bp, 99.9% within 100 bp**. A median of exactly zero is the signature
of 3'-end snapping -- 0.1.6's agreement is substantially *produced* by snapping
to annotation rather than observed. That is a real caveat on the 0.1.6 numbers
and it was not visible before.

Run the same on 0.1.10 and the spread says whether its extra 3'-end resolution
is signal or noise.

## 0.1.13

A regression fix. **0.1.10 made three-way reproducibility worse, and the cause
was mine.**

### What the whole-genome 0.1.10 run showed

| at >=20 molecules | 0.1.6 | 0.1.10 |
|---|---|---|
| structures | 68,903 | 77,732 |
| reproduced in all three | **88.0%** | **73.5%** |
| plateau at >=500 | 91.6% | 71.6% |

More structures at the same support threshold, reproducing worse. That cannot
be the lowered support floors -- those admit *low*-support models, which this
threshold removes. It is `max_3p_diff` 100 -> 30 in 0.1.7: the peak-clustering
window read the same field as the FSM naming tolerance, so a single cleavage
site's read pile started splitting into two or three peaks, at offsets that
differ between libraries.

Measured on BD144 kin, model pairs sharing a parent transcript whose 3' ends sit
31-100 bp apart:

| | 0.1.6 | 0.1.10 |
|---|---|---|
| all such pairs | 160 | **3,639** |
| among models with >=20 molecules | 39 | **1,416** |

`EndParams.peak_window` (default 100) is now separate from `max_3p_diff`
(default 30). Clustering asks "is this one cleavage site?" and wants to be
generous; naming asks "is this the annotated end?" and wants to be strict. They
were never the same question.

### Mono-exon naming, the same collision for the third time

`Y_RNA` is three different gene_ids on chr9, chr10 and chr16. The mono-exon
track named its own models with a counter that was **per contig** and keyed on
the gene *name*, so all three came out `Y_RNA|APA1`. This is the DNAJC9-AS1 bug
again, and it survived the 0.1.10 fix only because `monoexon.py` bypassed
`IdAssigner` by setting `model_id` itself.

It no longer sets it. All naming now happens in one place, on a counter keyed on
the exact string the name is built from, globally. The class is gone rather than
the instance.

### Not the cause, checked

Internal-deletion promotion was a candidate and is not: of 39,139 distinct gaps
in kin, **12** were promoted on motif and 0 on annotation, 46 records in total.

## 0.1.12

`flightcollapse compare --categories NAME=FILE`.

A collapse GFF carries structure and nothing else -- `isoseq collapse` writes no
category and no counts, because classification happens later in `pigeon
classify`. This reads that table back, so the same cross-sample comparison can
be drawn for any pipeline in **its own** vocabulary rather than one imposed on
it. Accepts a pigeon/SQANTI3 `*_classification.txt` (recognised by `isoform` +
`structural_category`, with support taken from `FL.*`) or any id/category TSV.

Matching remains structural in every case. `PB.7.2` is a per-run counter emitted
in collapse order, so the same string in two files is two unrelated transcripts;
an id-based overlap would be noise. There is now a test asserting exactly that,
in both directions: identical ids with different structures must not match, and
different ids with the same structure must.

## 0.1.11

Two CLI fixes, both from the three-way run.

* **`--no-strict` now returns 0.** Its help said "report failed invariants
  instead of exiting non-zero", but it only suppressed the exception and then
  returned 1 anyway. A calling script therefore reported "at least one run
  failed" for what was, in that run, a flag on 0.2% of models with every output
  written correctly. Failed invariants are now printed to stderr and the exit
  code is 0; without `--no-strict` it still exits 1.
* `flightcollapse --version`.

## 0.1.10

Two id-collision fixes, and the first change to what gets called that came from
evidence outside the long-read data.

### Novel junctions must be GT-AG

Scoring the finished BD144 runs against the summed short reads
(`flightcollapse sj -j ...`), matched on long-read support so that depth is not
doing the work:

| kept novel junctions | n (kin) | short-read confirmed |
|---|---|---|
| GT-AG | 42,118 | 5.8% |
| GC-AG + AT-AC | 4,699 | **0.13%** |

...while *annotated* semi-canonical junctions are confirmed at 33.4%. The short
reads see the motif class perfectly well; the novel members of it are simply not
there. Nor are they shifted copies of real junctions -- only 12% have an
annotated splice site within 5 bp, against 51% for novel GT-AG -- so widening
the snap window would not have helped. The numbers reproduce across all three
libraries (0.13 / 0.12 / 0.15%).

`junctions.novel_motif_classes` therefore defaults to canonical only. It removes
4,699 junctions per sample and costs 6 confirmed ones.

Novel semi-canonical junctions also join the **decoy set**. Their absence from
it is why the calibration could not see the problem: annotated GC-AG junctions
are real and sit in the positive anchor, so the model learned that the motif was
fine and applied that to the novel ones, scoring 97.1% of them at lfdr <= 0.05.
Tightening `max_junction_lfdr` from 0.05 to 0.001 moved precision 5.2% -> 6.0%
while discarding 8,600 junctions; the motif rule alone does better than that for
almost nothing.

`junctions.require_annotated_site_for_novel_junction` is new and **off**. Novel
GT-AG junctions with two invented splice sites confirm at 0.08-0.20 of the
annotated rate against 0.26-0.35 for those reusing a site, so the class is much
weaker -- but it is also the only place a genuinely new splice site can appear.
Turning it on removes a further 9,160 junctions and costs 182 confirmed ones.

### Model-id collisions

The whole-genome runs emitted 333 / 323 / 325 duplicate transcript ids out of
~151,000 models. Nothing was corrupted -- the `.dupN` suffix is applied before
anything is written -- but both causes were naming bugs:

* the 3'-end suffix index came from `utr_variant_rank`, a counter that advances
  only for *novel* 3' ends. Once a group emitted one, every later
  annotated-end model in it reused the index: RBBP4 produced three genuinely
  different transcripts at three different annotated 3' ends, all named
  `|3UTR_alt1`.
* the per-gene counter was keyed on `gene_id` while the name was built from
  `gene_name`, so two gene_ids sharing a name each started counting at 1 and
  collided (DNAJC9-AS1, two models 6 kb apart).

Suffix indices now count on `(stem, kind)` -- exactly the two pieces the name is
made of -- which makes both impossible by construction.

Separately, `IdAssigner` tracked "already assigned" in a set of `id(model)`.
That is a memory address: once a contig's models were released the addresses
were reused and the next contig's models would be silently skipped, leaving them
unnamed. It never fired in the pipeline, where every model is retained in
`all_models`, but it fired immediately in a test that frees each batch. The flag
now lives on the model.

## 0.1.9

Replicate concordance for the short-read anchor set, measured on BD144 a/b/c.

Summing technical replicates is right for *evidence* but destroys the
information that matters most at low counts: a junction reaching two reads
because both came from one library is a different object from one that appeared
once in each of two. These are replicates of the same cells, so the second
reproduces and the first may not.

Measured, of the 35,974 novel short-read junctions:

| | junctions |
|---|---|
| reaching 2 summed unique reads | 10,007 |
| seen in 3/3 replicates | 2,047 |
| seen in 2/3 | 3,658 |
| seen in 1/3 | 4,302 |
| ...expected in 1/3 by sampling alone | 1,812 |
| **excess sampling cannot explain** | **~2,490** |

The excess sits entirely at 2-5 summed reads (32% at n=2, 33% at n=3, gone by
n=10), which is the signature of library-specific noise rather than biology.

An anchor set wants precision far more than recall -- 5,705 positives is already
two orders of magnitude past what the calibration needs -- so
`junctions.min_sr_replicates` now defaults to **2**, clamped to the number of SJ
files supplied so a single library is unaffected. It costs roughly 1,800 genuine
junctions and removes roughly 2,490 artefacts.

The caveat this cannot address, and it matters: replicates catch *stochastic*
noise, not systematic misalignment. A junction produced by a pseudogene or a
repeat reproduces perfectly in all three libraries. Motif class and overhang are
the only guards against that class.

`flightcollapse sj` prints the concordance table and the sampling expectation
whenever more than one SJ file is given, and `sr_n_libs` is now a column in
`junctions.tsv.gz`.

## 0.1.8

Short-read splice junctions from STAR, used as an independent **anchor** for the
calibration rather than as a filter.

### Why an anchor and not a filter

The calibration's positive set was "annotated", which means the fitted score
could only ever learn to recognise things already in GENCODE -- the wrong prior
for finding real novelty. STAR's own `sjdb` flag separates its junctions into
annotated and not, and in `1317_BD_144b_HEK` that gives **17,280 junctions an
orthogonal library saw and the annotation does not contain**. Those are true
positives that annotation cannot supply, and they are what stops the model from
being an annotation detector.

Filtering on short reads would forfeit that, and would also make validation
impossible: a set used to filter cannot also be used to check the result. With
this library mapping 24% uniquely and 70.5% of reads discarded as "too short",
and with the long-read library 3'-biased in a way the short-read one is not,
absence of short-read support is in any case weak evidence.

So:

* short-read-supported novel junctions **join the positive anchor set**;
* they **leave the decoy set** -- a non-canonical motif with fifty
  uniquely-mapping short reads behind it is a real minor-class junction, and
  leaving it in the null biases the whole calibration;
* **half are held out** on a hash of their own coordinates, so recall measured
  on them is a prediction rather than a memory, and the split is identical in
  every sample;
* nothing is ever dropped for lacking short-read support.

`junctions.short_read_rescue` (default **off**) additionally exempts a supported
junction from curation. That one changes which models are called, which is why
it is opt-in. The asymmetry is deliberate and permanent: short reads can save a
junction, never condemn one.

### The coordinate check

STAR writes 1-based inclusive intron bounds; this package stores 0-based
half-open `(donor, acceptor)`, so `donor = col2 - 1` and `acceptor = col3`. Get
that wrong by one and nothing raises -- every junction comes back unsupported,
which is indistinguishable from an orthogonal library that confirms nothing. So
it is checked, not trusted: STAR's column 6 says which junctions were in its
splice-junction database, and those must be annotated here too.

Measured on the real data (chr1, 2, 17, 19, X): **37,403 / 37,407 = 99.99%**.
A deliberately shifted table reports `agreement 0.0` and names the offset; a
mismatched assembly names the contigs. `run` refuses to start on either.

Contig naming is normalised on both sides, so the SJ file's `1` and the BAM's
`chr1` are the same contig.

### New

* `--short-read-sj a.tab b.tab c.tab` on `run`; several files are **summed**,
  because technical replicates of the same cells are the same molecules sampled
  twice.
* `flightcollapse sj` -- check the coordinate convention and score a **finished**
  run's `junctions.tsv.gz` against the short reads, without collapsing anything
  again.
* `sr_uniq`, `sr_multi`, `sr_overhang`, `sr_in_sjdb`, `sr_supported`,
  `sr_anchor`, `sr_holdout` columns in `junctions.tsv.gz`.
* A short-read validation section in `qc.md`: recall on the held-out real
  junctions at each lFDR threshold, next to the fraction of *unsupported* novel
  junctions clearing the same threshold. The gap between those two columns is
  what the score is actually buying.

### Defaults

`min_sr_unique_reads = 2`, `min_sr_overhang = 10` (STAR's maximum spliced
overhang; a junction seen only through 5 bp of flank is an alignment guess),
`sr_holdout_frac = 0.5`.

## 0.1.7

Terminal-end language, a nearest-reference-transcript assignment for every
model, and the edge cases an exact intron-chain match does not cover.

### Changes that affect which models are called

* **Ends are compared against every transcript sharing the chain**, not the
  MANE one among them. Several GENCODE transcripts routinely share an intron
  chain and differ only in where they start and stop; measuring a read pile's
  ends against the canonical one manufactured an "alternative 3'UTR" for every
  pile that was in fact an exact match to a different annotated isoform.
  `ReferenceIndex.best_transcript_for_ends` ranks by 3' distance, then 5', with
  the canonical rank only as a deterministic tiebreak so that two
  indistinguishable transcripts are not named differently in two samples.
* **Intron retention is its own category.** A chain that is a reference chain
  minus introns the molecule read through was previously NIC. Retention
  requires exonic flank on both sides of the missing junction — without that
  test a 5'-truncated read looks identical to an unspliced one.
* **A 5' extension that crosses an annotated junction is retention, not a new
  TSS.** The junction-crossing test is the discriminator; the old code called
  both an extension.
* **`READTHROUGH`** for a 3' end that runs past the gene into the next one's
  exons on the same strand. Still requires polyA evidence, otherwise it cannot
  be told from internal priming in the downstream gene body.
* **Exon-internal deletions are no longer invisible.** Long CIGAR `D`
  operations are captured during the scan. Ones that are annotated junctions,
  or that recur with a canonical motif, are promoted to junctions *before* the
  census and curated like any other. The rest flag their reads
  (`n_internal_indel_reads`) and are never promoted — an invented junction
  would fail curation and take the whole read with it.
* **`max_3p_diff` 100 → 30.** Cleavage heterogeneity is ~20-30 bp wide; the 5'
  tolerance is a statement about how much degradation to forgive and now has
  its own value (`max_5p_diff`, 100).
* **Support floors lowered** (chains/ends/mono-exon: 10 reads → 5). The tool
  emits at a low floor and carries the calibrated local FDR, the category and
  every flag on each row, so the threshold is chosen on the finished table
  instead of being baked into a 3.5-hour run.

### Categories renamed

| 0.1.6 | 0.1.7 |
|---|---|
| `UTR3_annotated` | `end3_annotated` |
| `UTR3_variant` | `end3_novel` |
| `UTR5_extension` | `end5_extended`, or `end5_alt` at an annotated TSS |
| `UTR5_truncation` | `end5_truncated` |
| — | `IR`, `READTHROUGH` |

An identical intron chain confines a difference to the terminal exons; it does
not make that difference untranslated. `is_utr_only` reports the UTR question
separately, and is empty for a non-coding reference transcript where it has no
answer.

Name suffixes follow: `|alt3end1`, `|novel3end1`, `|alt5end1`, `|ext5end1`,
`|IR1`, `|readthrough1`. Novel models are named against their closest reference
transcript rather than their gene, so `GAPDH|NIC_3` becomes
`ENST00000229239.10|NIC_3`.

### Annotation added, no effect on calls

* `associated_transcript`, `structural_diff`, `n_equally_close` on every model.
  For a gene with forty isoforms "most similar" is often not well defined, and
  saying so beats pretending otherwise.
* `n_annotated_superchains` — how many annotated transcripts of the gene have
  this model's chain as a proper suffix. FSM to an annotated fragment stays
  FSM, but it moves read mass off the long isoform in proportion to how
  degraded the library is, and that is now visible rather than implicit.
* `ir_introns`, `ir_intron_rank`, and the `ir_event` flag on an FSM whose chain
  matches a `retained_intron` GENCODE model exactly. Retention at the 3' end is
  a confident observation under oligo-dT; retention of the 5'-most intron is
  confounded with incomplete reverse transcription and is flagged.
* `gene_conflict` — the gene the positional assignment would have chosen when
  it disagrees with the chain match. The chain match still wins; it just no
  longer wins silently.
* `n_5p_truncated_reads` and `five_prime_p10/p50/p90`. 5'-truncated piles merge
  into the parent rather than becoming models, and their mass is counted here.

## 0.1.6

Fragmentation fix in `build_models_for_group` (69,169 of 455,536 rows were
duplicate-named fragments of one transcript differing only at the 5' end), the
`model_ids_are_unique` invariant, and `flightcollapse compare` for cross-sample
structural concordance.
