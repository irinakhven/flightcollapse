"""Writers.

Everything ``isoseq collapse`` emits, plus the things whose absence made the
original bug hard to find:

* ``read_stat.txt`` keyed by **read name** -- keeping this is what made the
  whole diagnosis possible in the first place; do not drop it.
* ``group_chains.tsv`` -- per emitted model, the distinct intron chains among
  its member reads with support counts, so a heterogeneous group is visible as
  a heterogeneous group instead of being collapsed to one line.
* ``models.tsv`` -- one row per model with category, parent transcript, support
  in reads / molecules / cells, calibrated posterior, and every flag raised.

Read names are never held in memory: the writer re-walks the BAM with the same
admission rule used during scanning and pairs the n-th admitted read with row
``n`` of the arrays.
"""

from __future__ import annotations

import gzip
import json
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .config import Config
from .model import TranscriptModel
from .molecules import NO_CELL, NO_UMI, collapse_umis, decode_umi
from .reads import ContigReads, iter_alignments

GTF_SOURCE = "flightcollapse"


def _open_w(path: str, gz: bool = False):
    return gzip.open(path + ".gz", "wt") if gz else open(path, "w")


# ---------------------------------------------------------------------- #
class IdAssigner:
    """Stable ids across contigs.

    Reference-anchored names carry meaning (``ENST00000229239.10|3UTR_mod1``);
    PB ids keep cupcake-era downstream tools happy.  Both are always present:
    in ``ref`` style ``pb_id`` is a GTF attribute, in ``pb`` style the two are
    swapped.  The assigner is stateful so that ids stay unique when models are
    emitted contig by contig, and it is applied exactly once per model.
    """

    def __init__(self, id_style: str = "ref") -> None:
        self.id_style = id_style
        self.gene_index: Dict[str, int] = {}
        self.per_gene: Dict[str, int] = defaultdict(int)
        #: Suffix counter keyed on **exactly the string the name is built from**,
        #: which is the fix for two collision classes seen on the real data:
        #:
        #: * the 3'-end suffix index used to come from ``utr_variant_rank``, a
        #:   counter incremented only for *novel* 3' ends. Once a group emitted
        #:   one, every later annotated-end model in that group reused the same
        #:   index -- RBBP4 produced three genuinely different transcripts, at
        #:   three different annotated 3' ends, all called ``|3UTR_alt1``.
        #: * the per-gene counter was keyed on ``gene_id`` while the name was
        #:   built from ``gene_name``. Two gene_ids sharing a name each started
        #:   their own count at 1 and collided (DNAJC9-AS1, 6 kb apart).
        #:
        #: Counting on the name key makes both impossible by construction.
        self.per_name: Dict[Tuple[str, str], int] = defaultdict(int)
        self._names: set = set()
        #: A duplicate transcript_id silently corrupts read_stat, group.txt and
        #: the count matrix (two columns with one name). Kept as the net that
        #: makes a regression loud instead of invisible.
        self.collisions: int = 0

    def assign(self, models: Sequence[TranscriptModel]) -> None:
        for m in sorted(models, key=lambda x: x.sort_key()):
            # "already assigned" is marked on the model, not by object identity.
            # This used to be a set of id(m), which is a memory address: once a
            # contig's models were released the addresses got reused and the
            # next contig's models were silently skipped, leaving them unnamed.
            # It never fired in the pipeline (every model is retained in
            # all_models) but it is exactly the kind of thing that starts firing
            # the day someone frees them.
            if m.pb_id:
                continue
            gkey = m.gene_id or f"novel_{m.contig}"
            if gkey not in self.gene_index:
                self.gene_index[gkey] = len(self.gene_index) + 1
            self.per_gene[gkey] += 1
            n = self.per_gene[gkey]
            m.pb_id = f"PB.{self.gene_index[gkey]}.{n}"
            if not m.model_id:
                m.model_id = self._name(m)
            if m.model_id in self._names:
                self.collisions += 1
                k = 2
                while f"{m.model_id}.dup{k}" in self._names:
                    k += 1
                m.model_id = f"{m.model_id}.dup{k}"
                m.flags = sorted(set(m.flags) | {"duplicate_model_id"})
            self._names.add(m.model_id)
            if self.id_style == "pb":
                m.model_id, m.pb_id = m.pb_id, m.model_id

    def _next(self, stem: str, kind: str) -> int:
        self.per_name[(stem, kind)] += 1
        return self.per_name[(stem, kind)]

    def _name(self, m: TranscriptModel) -> str:
        """Reference-anchored name whose index cannot repeat.

        The index counts within ``(stem, kind)`` -- the two pieces the name is
        actually made of -- so the n-th ``alt3end`` of a transcript is numbered
        by how many ``alt3end`` models that transcript already has, and nothing
        else can share it.
        """
        if m.category == "FSM" and m.parent_tx:
            return m.parent_tx
        if m.category == "monoexon_intergenic":
            return f"NOVELMONO_{m.contig}_{self._next(m.contig, m.category)}"
        # For a novel model the exact-chain lookup found nothing, so the handle
        # is the *most similar* reference transcript. `GAPDH|NIC_3` says only
        # which locus; `ENST00000229239.10|NIC_3` says which transcript it is a
        # variation on, which is the question anyone reading the table has.
        # Falling back to the gene keeps the readable form when there is no
        # transcript to point at -- an intergenic mono-exon pile, say.
        stem = m.parent_tx or m.associated_tx or m.gene_name or m.gene_id or m.contig
        kind = _SUFFIX.get(m.category)
        if kind:
            return f"{stem}|{kind}{self._next(stem, kind)}"
        return f"{stem}|{m.category}_{self._next(stem, m.category)}"


def assign_ids(models: List[TranscriptModel], id_style: str = "ref") -> None:
    IdAssigner(id_style).assign(models)


#: category -> name suffix.  The distinction between ``alt`` (the end coincides
#: with an annotated end of another isoform of the same gene) and ``novel`` (the
#: end matches nothing annotated and rests on read mass plus polyA evidence) is
#: a large confidence gap -- roughly the 35% vs 72% gap measured on this data --
#: so the two must not share a word.
_SUFFIX = {
    "end3_annotated": "alt3end",
    "end3_novel": "novel3end",
    "end5_alt": "alt5end",
    "end5_extended": "ext5end",
    "end5_truncated": "trunc5end",
    "IR": "IR",
    "READTHROUGH": "readthrough",
    "monoexon_3UTR": "3UTRfrag",
    "monoexon_internal": "APA",
}


# ---------------------------------------------------------------------- #
def write_gtf(models: Sequence[TranscriptModel], path: str, feature: str = "exon") -> None:
    with open(path, "w") as fh:
        fh.write("##description: flightcollapse transcript models\n")
        fh.write("##format: gtf\n")
        for m in sorted(models, key=lambda x: x.sort_key()):
            attrs = _attrs(m)
            fh.write(
                "\t".join(
                    [m.contig, GTF_SOURCE, "transcript", str(m.start + 1), str(m.end),
                     ".", m.strand, ".", attrs]
                )
                + "\n"
            )
            exons = m.exons if m.strand == "+" else m.exons
            for k, (s, e) in enumerate(exons, 1):
                fh.write(
                    "\t".join(
                        [m.contig, GTF_SOURCE, feature, str(s + 1), str(e),
                         ".", m.strand, ".", attrs + f' exon_number "{k}";']
                    )
                    + "\n"
                )


def write_gff(models: Sequence[TranscriptModel], path: str) -> None:
    """cupcake/pigeon-style collapsed GFF (transcript + exon rows, GTF attributes)."""
    write_gtf(models, path, feature="exon")


def _attrs(m: TranscriptModel) -> str:
    gene = m.gene_id or (m.gene_name or "novelGene")
    parts = [
        f'gene_id "{gene}";',
        f'transcript_id "{m.model_id}";',
        f'gene_name "{m.gene_name or gene}";',
        f'pb_id "{m.pb_id}";',
        f'category "{m.category}";',
        f'n_reads "{m.n_reads}";',
        f'n_molecules "{m.n_mols}";',
        f'n_cells "{m.n_cells}";',
    ]
    if m.parent_tx:
        parts.append(f'parent_transcript "{m.parent_tx}";')
    if m.associated_tx:
        parts.append(f'associated_transcript "{m.associated_tx}";')
    if m.structural_diff:
        parts.append(f'structural_diff "{m.structural_diff}";')
    if m.is_utr_only is not None:
        parts.append(f'is_utr_only "{int(m.is_utr_only)}";')
    if m.posterior_real == m.posterior_real:
        parts.append(f'posterior_real "{m.posterior_real:.4f}";')
    if m.flags:
        parts.append(f'flags "{",".join(m.flags)}";')
    return " ".join(parts)


# ---------------------------------------------------------------------- #
def write_models_table(models: Sequence[TranscriptModel], path: str) -> None:
    cols = [
        "model_id", "pb_id", "contig", "strand", "start", "end", "n_exons",
        "length", "category", "parent_transcript", "associated_transcript",
        "structural_diff", "n_equally_close", "is_utr_only",
        "gene_id", "gene_name", "gene_conflict",
        "n_reads", "n_molecules", "n_cells", "posterior_real", "lfdr",
        "polya_motif_found", "polya_motif", "dist_to_polya_site",
        "perc_a_downstream", "internal_priming", "dist_to_ref_tts",
        "dist_to_ref_tss", "end3_from", "three_prime_dispersion",
        "n_5p_truncated_reads", "five_prime_p10", "five_prime_p50", "five_prime_p90",
        "n_annotated_superchains", "n_retained_introns", "ir_introns",
        "ir_intron_rank", "n_internal_indel_reads",
        # 0.1.17: the 3'-end diagnostics. tts_shift_from_reads is how far
        # annotation moved this model's end off its own read pile -- the number
        # the NSL1 class of failure lives in, and unreadable from any other
        # column. The tail columns are the direct polyA evidence; they are
        # recorded so a threshold can be chosen from the distribution rather
        # than guessed, so they must survive into the table.
        "tts_shift_from_reads", "n_mol_equivalents",
        "tail_molecule_frac", "n_tail_reads", "median_tail_len",
        "flags",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in sorted(models, key=lambda x: x.sort_key()):
            ev = m.evidence
            row = [
                m.model_id, m.pb_id, m.contig, m.strand, m.start + 1, m.end,
                m.n_exons, m.tx_len, m.category, m.parent_tx or "",
                m.associated_tx or "", m.structural_diff,
                m.n_equally_close or "",
                "" if m.is_utr_only is None else int(m.is_utr_only),
                m.gene_id or "", m.gene_name, ev.get("gene_conflict", ""),
                m.n_reads, m.n_mols, m.n_cells,
                _fmt(m.posterior_real), _fmt(m.lfdr),
                ev.get("polya_motif_found", ""), ev.get("polya_motif", ""),
                ev.get("dist_to_polya_site", ""), _fmt(ev.get("perc_a_downstream")),
                ev.get("internal_priming", ""), ev.get("dist_to_ref_tts", ""),
                ev.get("dist_to_ref_tss", ""), ev.get("end3_from", ""),
                _fmt(ev.get("three_prime_dispersion")),
                ev.get("n_5p_truncated_reads", ""),
                _fmt(ev.get("five_prime_p10")), _fmt(ev.get("five_prime_p50")),
                _fmt(ev.get("five_prime_p90")),
                ev.get("n_annotated_superchains", ""),
                ev.get("n_retained_introns", ""), ev.get("ir_introns", ""),
                ev.get("ir_intron_rank", ""), ev.get("n_internal_indel_reads", ""),
                ev.get("tts_shift_from_reads", ""),
                _fmt(ev.get("n_mol_equivalents")),
                _fmt(ev.get("tail_molecule_frac")),
                ev.get("n_tail_reads", ""), _fmt(ev.get("median_tail_len")),
                ",".join(m.flags),
            ]
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def _fmt(v) -> str:
    if v is None:
        return ""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return "" if f != f else f"{f:.4g}"


def write_abundance(models: Sequence[TranscriptModel], path: str, prefix_lines: Sequence[str] = ()) -> None:
    tot_reads = sum(m.n_reads for m in models) or 1
    tot_mols = sum(m.n_mols for m in models) or 1
    with open(path, "w") as fh:
        for line in prefix_lines:
            fh.write(f"#{line}\n")
        fh.write("pbid\tcount_fl\tnorm_fl\tcount_umi\tnorm_umi\n")
        for m in sorted(models, key=lambda x: x.sort_key()):
            fh.write(
                f"{m.model_id}\t{m.n_reads}\t{m.n_reads / tot_reads:.6e}\t"
                f"{m.n_mols}\t{m.n_mols / tot_mols:.6e}\n"
            )


def write_flnc_count(models: Sequence[TranscriptModel], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("id\tcount_fl\n")
        for m in sorted(models, key=lambda x: x.sort_key()):
            fh.write(f"{m.model_id}\t{m.n_reads}\n")


# ---------------------------------------------------------------------- #
class ReadAssignmentWriter:
    """Streams ``read_stat.txt`` and ``group.txt`` one contig at a time."""

    def __init__(self, cfg: Config) -> None:
        gz = cfg.output.gzip_big_tables
        self.read_stat_path = cfg.outpath(".read_stat.txt")
        self.group_path = cfg.outpath(".group.txt")
        self._rs = _open_w(self.read_stat_path, gz)
        self._rs.write("id\tlength\tis_fl\tstat\tpbid\n")
        self._grp = _open_w(self.group_path, gz)
        self.n_written = 0
        self.n_unassigned = 0

    def write_contig(
        self,
        bam,
        contig: str,
        cfg: Config,
        reads: ContigReads,
        model_of_read: np.ndarray,
        model_ids: Sequence[str],
    ) -> None:
        groups: Dict[int, List[str]] = defaultdict(list)
        for i, r in enumerate(iter_alignments(bam, contig, cfg.gates)):
            if i >= reads.n:
                break
            mi = int(model_of_read[i])
            if mi < 0:
                self.n_unassigned += 1
                continue
            mid = model_ids[mi]
            length = int(reads.end[i] - reads.start[i])
            self._rs.write(f"{r.query_name}\t{length}\tY\tunique\t{mid}\n")
            groups[mi].append(r.query_name)
            self.n_written += 1
        for mi, names in groups.items():
            self._grp.write(f"{model_ids[mi]}\t{','.join(names)}\n")

    def close(self) -> None:
        self._rs.close()
        self._grp.close()


# ---------------------------------------------------------------------- #
def write_group_chains(
    path: str,
    rows: Sequence[Tuple[str, str, int, int, bool]],
) -> None:
    """``model_id, chain_repr, n_reads, n_molecules, is_model_chain``."""
    with open(path, "w") as fh:
        fh.write("model_id\tintron_chain\tn_reads\tn_molecules\tis_model_chain\n")
        for mid, ch, nr, nm, isc in rows:
            fh.write(f"{mid}\t{ch}\t{nr}\t{nm}\t{int(isc)}\n")


def chain_repr(chain) -> str:
    return ";".join(f"{d}-{a}" for d, a in chain) if chain else "."


# ---------------------------------------------------------------------- #
def write_molecule_matrix(
    models: Sequence[TranscriptModel],
    cell_counts: Dict[str, Dict[int, int]],
    barcodes: Sequence[str],
    cfg: Config,
) -> Optional[str]:
    """Cell x transcript UMI matrix in MatrixMarket form.

    ``cell_counts`` maps ``model_id -> {cell_id: n_molecules}`` and is built
    contig by contig so the read arrays can be released as we go.

    File naming follows the convention already used in the project:
    ``{prefix}_umi_corrected_isoform_matrix.mtx`` / ``_cells.txt`` /
    ``_isoforms.txt`` / ``_isoform_metadata.txt``.
    """
    ordered = sorted(models, key=lambda x: x.sort_key())
    used_cells: Dict[int, int] = {}
    entries: List[Tuple[int, int, int]] = []

    for col, m in enumerate(ordered, 1):
        for cid, n in (cell_counts.get(m.model_id) or {}).items():
            row = used_cells.get(cid)
            if row is None:
                row = used_cells[cid] = len(used_cells) + 1
            entries.append((row, col, n))
    if not entries:
        return None

    mtx = cfg.outpath("_umi_corrected_isoform_matrix.mtx")
    with open(mtx, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{len(used_cells)} {len(ordered)} {len(entries)}\n")
        for r, c, v in entries:
            fh.write(f"{r} {c} {v}\n")

    inv = [""] * len(used_cells)
    for cid, row in used_cells.items():
        inv[row - 1] = barcodes[cid] if 0 <= cid < len(barcodes) else f"cell{cid}"
    with open(cfg.outpath("_cells.txt"), "w") as fh:
        fh.write("\n".join(inv) + "\n")
    with open(cfg.outpath("_isoforms.txt"), "w") as fh:
        fh.write("\n".join(m.model_id for m in ordered) + "\n")
    with open(cfg.outpath("_isoform_metadata.txt"), "w") as fh:
        fh.write(
            "isoform\tgene_id\tgene_name\tcategory\tparent_transcript\t"
            "associated_transcript\tstructural_diff\tn_exons\tlength\n"
        )
        for m in ordered:
            fh.write(
                f"{m.model_id}\t{m.gene_id or ''}\t{m.gene_name}\t{m.category}\t"
                f"{m.parent_tx or ''}\t{m.associated_tx or ''}\t"
                f"{m.structural_diff}\t{m.n_exons}\t{m.tx_len}\n"
            )
    return mtx


def write_json(obj, path: str) -> None:
    def _default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.ndarray,)):
            return o.tolist()
        if isinstance(o, (set, frozenset)):
            return sorted(o)
        return str(o)

    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_default)
