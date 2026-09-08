# Validation of the v0.1.18 algorithm scheme

Validation date: 2026-09-08

## Scope

This was a static audit of the packaged source plus automated tests. It checks
that the diagram describes the implementation; it is not an independent
biological replication of the full BD67, BD70, or BD144 benchmarks.

The audited source archive was:

```text
flightcollapse-0.1.18.tar.gz
SHA-256 e193dd056e7643a511fbce52fe40c80b90807370e5443e6a708e6aba27e93dd9
```

Every runtime `.py` file in the archive matched the installed v0.1.18 package
byte-for-byte. The diagram PDF supplied for review matched the copy retained in
the stable evaluation directory:

```text
flightcollapse_logic_v0.1.18.pdf
SHA-256 a017bcd5290a2d2bbab2c30119c4342be4fb3c9a03e94afbe36e14e045521e64
```

## Findings

| Diagram claim | Implementation evidence | Result |
|---|---|---|
| Shared read-admission gates | `ContigReads.scan` applies coverage, identity, terminal-clip, mapping, and alignment-status gates once; the admitted object feeds census and output assignment. | Confirmed |
| Annotation-first junction snapping | `junctions.build_snap_maps` prefers a reference donor/acceptor within the fuzzy tolerance before a read-derived seed. | Confirmed |
| Internal deletions promoted before ordinary junction curation | The pipeline promotes recurrent/annotated CIGAR deletions and then subjects them to the same census and curation path. | Confirmed |
| Junctions curated before chains | The pipeline builds features, scores and curates junctions, excludes reads carrying a rejected junction, and only then calls `build_chain_groups`. | Confirmed |
| Annotated junctions are not rejected | `junctions.curate` applies each rejection mask only to non-annotated junctions. | Confirmed |
| Novel-junction defaults shown in the figure | Defaults are >=3 molecules when available, >=5 reads, >=0.1% share at both splice sites, direct repeat <8 bp, and local FDR <=0.05. GC-AG/AT-AC novel junctions require short-read or high long-read support. | Confirmed |
| Short reads are an anchor, not a negative filter | Short-read-supported novel junctions can enter the positive anchor set. Absence of support never rejects a junction; rescue is opt-in. | Confirmed |
| Exact-chain grouping and maximal suffix representative | `build_chain_groups`, `find_parents`, and `resolve_suffixes` implement exact curated-chain groups and directed suffix merging into maximal parents. Annotated chains are retained; an annotated TSS can protect an alternative start. | Confirmed |
| Molecule-weighted 3' peak geometry | `ends.split_ends` uses `molecule_weights` for clustering, ranking and placement while support gates retain whole-read/whole-molecule counts. | Confirmed |
| Direct read-tail evidence in v0.1.18 | `reads.polya_tail_stats`, `ends.tail_evidence`, and the default `polya_tail_gates_novel_end=True` implement the clipped-tail evidence described in the figure. | Confirmed |
| Annotation cannot manufacture a distant 3' end | `_fallback` snaps only within `max_3p_fallback_dist=300`; farther peaks keep the observed coordinate and become `end3_unresolved`. | Confirmed |
| 5' truncation is absorbed; 5' extension may be retained or classified as IR | Stage C of `ends.split_ends` uses the 90th-percentile start, absorbs truncation by default, and checks whether an extension crosses an annotated junction. | Confirmed |
| Mono-exonic reads use a separate track | `collapse_monoexonic` selects reads with zero introns and groups only those reads; last-exon/3'UTR fragments, internal models, and intergenic models follow the three policies shown. | Confirmed |
| Per-run invariants | `validate.Validator` checks mono-exon mass, representative minimality, span agreement, unspliced-read accounting, ID uniqueness, and exon-count agreement; strict mode fails the run when an invariant fails. | Confirmed |

## Important qualifications

- Molecule-based support requires a working barcode/UMI table. When no molecule
  identifiers are available, the implementation deliberately falls back to
  read support.
- Motif claims assume a matching genome FASTA was supplied. Without sequence,
  motif class is `unknown`, which the default configuration permits rather than
  discarding all otherwise evaluable junctions.
- A read carrying a rejected junction is excluded from transcript construction
  and contributes to the negative/decoy evidence used by the novelty model.
- `parent_transcript` has the strict exact-chain meaning for multi-exon terminal
  models. `associated_transcript` is a separate nearest-reference annotation for
  every model. Mono-exon models use a nearest-end parent assignment.
- The diagram's output label `junctions.tsv` corresponds to the compressed file
  actually written as `.junctions.tsv.gz`.

## Test result and source-archive issue

The original v0.1.18 source archive contains the test modules but accidentally
omits `tests/conftest.py` and `tests/simulate.py`, leaving 37 integration tests
without fixtures. The two support files were recovered from the adjacent
development tree and added only to this GitHub-ready copy. They do not modify
runtime code. With them present:

```text
203 passed, 4 warnings in 60.38s
```

The warnings are NumPy `Mean of empty slice` warnings in simulated end-to-end
cases. No test failed.
