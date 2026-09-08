"""Which reference transcript is a model *most like*, and how does it differ?

Until 0.1.7 a model carried a ``parent_transcript`` only when its intron chain
matched a GENCODE chain **exactly** -- the chain was a dictionary key, and the
lookup either hit or missed.  Everything novel therefore came out with a gene
and nothing else: ``GAPDH|NIC_3`` tells you the locus and hides whether the
model is "the canonical transcript minus exon 4" or "the canonical transcript
with a retained intron".  The pipeline knows which; it was throwing it away.

This module supplies the missing half: for any model, the transcript of its
gene with the most similar exon structure, plus a description of the
difference.  It is deliberately **annotation only**.  The nearest transcript is
never allowed to feed back into merging, thresholding or the novelty
calibration -- the calibration is fitted against annotation, so letting an
annotation-derived similarity score influence what gets called would make it
circular.

Three things are computed, all confined to the transcripts of one gene (a few
dozen, so the quadratic comparison is free):

``shared`` / ``model_only`` / ``ref_only`` junction sets
    The structural diff.  ``model_only`` are junctions the model has and the
    reference does not; ``ref_only`` the reverse.

retained introns
    A ``ref_only`` junction whose whole span sits inside a single model exon,
    with real exonic flank on **both** sides, is a retained intron.  The flank
    requirement is what separates retention from truncation: if the model's 5'
    end lies inside the intron there is no upstream flank, and the molecule is
    a 5'-truncated form of the long isoform rather than an unspliced one.

``is_utr_only``
    "UTR" is a coding-frame concept, not a structural one.  An identical intron
    chain guarantees only that the difference is confined to the terminal
    exons; whether that difference is *untranslated* depends on where the CDS
    ends.  A model whose 5' end lies downstream of the start codon has a
    different coding sequence, and calling that a 5'UTR difference overclaims.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .intervals import Chain, Exon, Junction, signed_5p_offset

#: exonic bases required on each side of a ref-only junction before it counts as
#: a retained intron rather than an end that happens to fall inside an intron
IR_MIN_FLANK = 10


@dataclass
class RefMatch:
    """The most similar reference transcript, and how the model differs from it."""

    tx_id: Optional[str] = None
    n_shared: int = 0
    model_only: Tuple[Junction, ...] = ()
    ref_only: Tuple[Junction, ...] = ()
    retained: Tuple[Junction, ...] = ()          # ref_only, contiguously covered
    d5: int = 0                                  # transcript orientation, <0 = upstream
    d3: int = 0                                  # transcript orientation, >0 = downstream
    chain_equal: bool = False
    chain_suffix: bool = False                   # model chain is a proper suffix
    is_utr_only: Optional[bool] = None           # None = no CDS, question undefined
    n_equally_close: int = 1
    score: float = 0.0

    @property
    def n_model_only(self) -> int:
        return len(self.model_only)

    @property
    def n_ref_only(self) -> int:
        return len(self.ref_only)

    @property
    def n_retained(self) -> int:
        return len(self.retained)

    @property
    def is_pure_retention(self) -> bool:
        """Model = reference minus one or more introns it simply read through."""
        return (
            self.tx_id is not None
            and not self.model_only
            and bool(self.retained)
            and len(self.retained) == len(self.ref_only)
        )

    def describe(self) -> str:
        """Compact structural diff, e.g. ``+1jxn;-2IR;d3=+412;d5=-88``."""
        bits: List[str] = []
        if self.model_only:
            bits.append(f"+{len(self.model_only)}jxn")
        n_lost = len(self.ref_only) - len(self.retained)
        if self.retained:
            bits.append(f"-{len(self.retained)}IR")
        if n_lost > 0:
            bits.append(f"-{n_lost}jxn")
        if self.d3:
            bits.append(f"d3={self.d3:+d}")
        if self.d5:
            bits.append(f"d5={self.d5:+d}")
        if not bits:
            return "identical"
        if self.chain_equal:
            bits.insert(0, "chain=")
        return ";".join(bits)


# ---------------------------------------------------------------------- #
def _covers(exons: Sequence[Exon], j: Junction, flank: int) -> bool:
    """Is junction ``j`` spanned by a single exon, with ``flank`` bases either side?

    This is the retained-intron test.  ``j`` is ``(donor, acceptor)`` 0-based
    half-open, so the intronic bases are ``[donor, acceptor)`` and the exon must
    contain ``[donor - flank, acceptor + flank)``.
    """
    lo, hi = int(j[0]), int(j[1])
    for s, e in exons:
        if s <= lo - flank and e >= hi + flank:
            return True
    return False


def compare_to_transcript(
    chain: Chain,
    exons: Sequence[Exon],
    strand: str,
    tx,
    ir_min_flank: int = IR_MIN_FLANK,
    ref_chain_set: Optional[frozenset] = None,
) -> RefMatch:
    """Structural comparison of one model against one reference transcript."""
    M = set(chain)
    R = ref_chain_set if ref_chain_set is not None else set(tx.chain)
    shared = M & R
    model_only = tuple(sorted(M - R))
    ref_only = tuple(sorted(R - M))
    retained = tuple(j for j in ref_only if _covers(exons, j, ir_min_flank))

    five = exons[0][0] if strand == "+" else exons[-1][1]
    three = exons[-1][1] if strand == "+" else exons[0][0]
    d5 = signed_5p_offset(five, tx.tss, strand)
    d3 = signed_5p_offset(three, tx.tts, strand)

    n = len(chain)
    m = RefMatch(
        tx_id=tx.tx_id,
        n_shared=len(shared),
        model_only=model_only,
        ref_only=ref_only,
        retained=retained,
        d5=int(d5),
        d3=int(d3),
        chain_equal=(n > 0 and tuple(chain) == tuple(tx.chain)),
        chain_suffix=(
            0 < n < len(tx.chain) and tuple(tx.chain)[-n:] == tuple(chain)
        ),
    )
    m.is_utr_only = _utr_only(m, five, three, strand, tx)
    # A junction the model invents is a much bigger difference than one it
    # merely fails to show, because 5' degradation removes junctions for free.
    m.score = (
        2.0 * len(shared)
        - 3.0 * len(model_only)
        - 1.0 * (len(ref_only) - len(retained))
        - 0.5 * len(retained)
    )
    return m


def _utr_only(m: RefMatch, five: int, three: int, strand: str, tx) -> Optional[bool]:
    """Does the difference lie entirely outside the reference CDS?"""
    if m.model_only or m.ref_only:
        return False                      # a structural difference, not an end one
    start_pos = getattr(tx, "start_pos", None)
    stop_pos = getattr(tx, "stop_pos", None)
    if start_pos is None or stop_pos is None:
        return None                       # non-coding: the question has no answer
    if m.d5 and signed_5p_offset(five, start_pos, strand) > 0:
        return False                      # model starts inside the CDS
    if m.d3 and signed_5p_offset(three, stop_pos, strand) < 0:
        return False                      # model ends inside the CDS
    return True


# ---------------------------------------------------------------------- #
def nearest_transcript(
    ref,
    gene_id: Optional[str],
    chain: Chain,
    exons: Sequence[Exon],
    strand: str,
    ir_min_flank: int = IR_MIN_FLANK,
    candidates: Optional[Sequence[str]] = None,
    require_retention: bool = False,
) -> Optional[RefMatch]:
    """The transcript of ``gene_id`` whose structure is closest to the model.

    Ranking, in order: structural score, then |3' distance|, then |5' distance|,
    then the canonical rank, then length, then the id.  The 3' end outranks the
    5' end because with oligo-dT capture it is the anchored, informative one;
    the canonical rank appears only as a deterministic last resort, so that two
    genuinely indistinguishable transcripts do not get named differently in two
    samples of the same cDNA.

    ``n_equally_close`` reports how many transcripts tied on the first three
    keys -- for a gene with forty annotated isoforms "most similar" is often not
    well defined, and saying so is better than pretending otherwise.

    ``require_retention`` restricts the search to transcripts the model has
    actually read an intron through.  It exists for one case: a model whose
    chain equals a short isoform's chain exactly, but whose first exon runs back
    across a *longer* isoform's junction.  Scored plainly, the exact chain match
    wins and the row would say "5' extension of the short isoform" when the
    informative statement is "unspliced first intron of the long one".
    """
    tids = list(candidates) if candidates is not None else None
    if tids is None:
        if not gene_id or gene_id not in ref.genes:
            return None
        tids = ref.genes[gene_id].tx_ids
    if not tids:
        return None

    scored: List[Tuple[Tuple, RefMatch]] = []
    for tid in tids:
        tx = ref.tx.get(tid)
        if tx is None:
            continue
        m = compare_to_transcript(
            chain, exons, strand, tx, ir_min_flank, ref.chain_set(tid)
        )
        key = (-m.score, abs(m.d3), abs(m.d5), tx.rank, -tx.tx_len, tid)
        scored.append((key, m))
    if not scored:
        return None
    if require_retention:
        with_ir = [s for s in scored if s[1].retained]
        if with_ir:
            scored = with_ir
    scored.sort(key=lambda kv: kv[0])
    best_key, best = scored[0]
    head = best_key[:3]
    best.n_equally_close = sum(1 for k, _ in scored if k[:3] == head)
    return best


# ---------------------------------------------------------------------- #
def retained_intron_rank(chain_ref: Chain, retained: Sequence[Junction]) -> int:
    """Position of the most 5' retained intron in the reference chain.

    ``0`` means the 5'-most intron of the transcript was the one read through,
    which is the case confounded with 5' degradation and incomplete reverse
    transcription; a large index means retention deep in the 3' end of the
    molecule, which for an oligo-dT library is a confident observation because
    the molecule was captured by its polyA tail and sequenced through.
    """
    if not retained:
        return -1
    idx = {j: i for i, j in enumerate(chain_ref)}
    # chain is in transcript orientation, so index 0 is the 5'-most intron
    return min((idx.get(j, len(chain_ref)) for j in retained), default=-1)


def readthrough_gene(
    ref,
    contig: str,
    strand: str,
    gene_id: Optional[str],
    three: int,
) -> Optional[str]:
    """Gene whose exons the model's 3' end has run into, if any.

    A 3' end past the parent gene's annotated extent that lands in another
    gene's exons on the same strand is transcriptional readthrough, not an
    alternative cleavage site.  It is worth its own category: readthrough is
    stress- and preparation-responsive, so it is a plausible source of
    disagreement between libraries built from the same cDNA.
    """
    if not gene_id or gene_id not in ref.genes:
        return None
    g = ref.genes[gene_id]
    if g.start <= three < g.end:
        return None
    for other in ref.genes_overlapping(contig, three, three + 1, strand):
        if other == gene_id:
            continue
        og = ref.genes[other]
        if any(s <= three < e for s, e in og.exon_union):
            return other
    return None
