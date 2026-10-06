"""SQANTI3-compatible structural classification.

New in 0.4.2.  Emits three orthogonal fields beside the existing ``category``,
which is left untouched so every downstream consumer keeps working:

``structural_category``  full-splice_match / incomplete-splice_match /
                         novel_in_catalog / novel_not_in_catalog / fusion
``subcategory``          SQANTI3's own subcategory strings
``terminal_class``       where the 3' end sits, independent of the chain

WHY THIS EXISTS.  ``_novel_category`` was two lines -- no novel junction means
NIC, otherwise NNC -- so flightcollapse had no ISM concept at all and every
5'-truncated fragment of a reference transcript landed in NIC.  Measured on
BD67: about 78% of NIC UMI mass is on models SQANTI3 calls ISM, and 25.7% of
NNC UMI is on models whose every non-annotated junction uses two ANNOTATED
sites, which is NIC by definition.  At RTN4 that produced twelve separate NIC
models that are all the same single-intron fragment of RTN4-202.

The second over-call is in the names.  ``end3_annotated``, ``end3_novel``,
``end5_alt`` and ``end5_extended`` all have an intron chain identical to an
annotated transcript (see the ends.py docstring), so under SQANTI3 they are FSM
with an alternative-end subcategory -- not four separate kinds of novelty.  On
BD176c that is the difference between reporting 16.3% FSM and 48.9%, and
between 41.6% of read mass and 80.5%.

NOTHING IS DROPPED HERE.  This module only labels.  A 3'-truncated fragment
gets ``incomplete-splice_match`` / ``3prime_fragment`` and its internal-priming
evidence columns, and is emitted; SQANTI3 reports those too.  Whether to filter
them is a decision for the caller, made on a label rather than inside a
classifier.

IR AND FUSION ARE REAL CATEGORIES.  Intron retention is an ISM *subcategory* in
SQANTI3 ("ISM with Intron Retention", for a reference SJ lost to an IR event),
and fusion is a top-level category.  Neither is an artifact class, and neither
is filtered here.  The IR flag is computed for every model regardless of
category, because a retained intron is a fact about the model, not a verdict on
it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .intervals import Chain, Exon

# ---------------------------------------------------------------------- #
# SQANTI3's own strings, so a downstream tool that already speaks SQANTI3
# needs no translation table.
FSM = "full-splice_match"
ISM = "incomplete-splice_match"
NIC = "novel_in_catalog"
NNC = "novel_not_in_catalog"
FUSION = "fusion"
GENIC = "genic"
INTERGENIC = "intergenic"

#: The short forms flightcollapse's own QC and figures use.
SHORT = {FSM: "FSM", ISM: "ISM", NIC: "NIC", NNC: "NNC", FUSION: "fusion",
         GENIC: "genic", INTERGENIC: "intergenic"}


@dataclass
class Structural:
    structural_category: str = ""
    subcategory: str = ""
    #: reference transcripts satisfying the FSM/ISM test, best first
    matches: Tuple[str, ...] = ()
    primary: Optional[str] = None
    #: ISM only, in TRANSCRIPT orientation
    miss5: int = 0
    miss3: int = 0
    #: annotated introns lying entirely inside one of the model's exons
    retained_introns: Tuple[Tuple[int, int], ...] = ()
    #: junctions whose donor or acceptor is not an annotated site
    novel_sites: int = 0
    #: junctions that are unannotated but use two annotated sites
    novel_combination_junctions: int = 0

    def as_columns(self) -> Dict[str, object]:
        return {
            "structural_category": self.structural_category,
            "subcategory": self.subcategory,
            "sqanti_primary_ref": self.primary or "",
            "miss5_introns": self.miss5,
            "miss3_introns": self.miss3,
            "n_retained_introns_any": len(self.retained_introns),
            "n_novel_splice_sites": self.novel_sites,
            "n_novel_combination_junctions": self.novel_combination_junctions,
        }


# ---------------------------------------------------------------------- #
def _near(a: int, b: int, tol: int) -> bool:
    return abs(a - b) <= tol


def _reciprocal(a, b, frac: float) -> bool:
    """Reciprocal overlap, for the one case with no junctions to compare."""
    ov = min(a[1], b[1]) - max(a[0], b[0])
    if ov <= 0:
        return False
    return ov >= frac * (a[1] - a[0]) and ov >= frac * (b[1] - b[0])


def genomic(chain: Chain) -> Chain:
    """Any chain, in ascending genomic order -- whatever order it arrived in.

    ORIENTATION IS THE TRAP HERE, and it has now bitten twice in opposite
    directions, so this function exists to make the question unaskable.

    ``ContigReads.chain`` returns ``ch[::-1]`` on the minus strand and
    ``Transcript.finalise`` calls ``exons_to_chain``, which does the same: both
    a model chain and a reference chain are stored in TRANSCRIPT orientation,
    which is why the collapser's own FSM lookup has always been strand-balanced
    (BD144 kin: parent_transcript present on 56.5% of plus-strand multi-exon
    models and 55.7% of minus).  0.4.2 compared them as they came and was right.

    What was wrong in 0.4.2 was ``miss5``/``miss3``: it flipped them by strand
    on top of a chain that was already transcript-oriented, so the two were
    swapped on the minus strand.  And 0.4.4 "fixed" a strand skew that was never
    in flightcollapse at all -- it was in a post-hoc script that re-sorted chains
    into genomic order before calling in -- by reversing the reference, which
    would have broken the minus strand inside the pipeline instead.

    Junctions in a chain are disjoint, so sorting is unambiguous and idempotent.
    Everything below therefore works in genomic order, ``left``/``right`` mean
    genomic left and right, and the strand flip that turns them into 5'/3' is
    correct exactly once.
    """
    return tuple(sorted(chain))


def chains_equal(c1: Chain, c2: Chain, tol: int) -> bool:
    return len(c1) == len(c2) and all(
        _near(x[0], y[0], tol) and _near(x[1], y[1], tol) for x, y in zip(c1, c2)
    )


def subchain_offsets(c: Chain, ref: Chain, tol: int) -> List[int]:
    """Offsets k where ``c`` sits inside ``ref`` as a CONTIGUOUS proper subchain.

    Proper: ``len(c) < len(ref)``.  Contiguous: no reference junction may be
    skipped in the middle, because skipping one is an exon-skip event and that
    is NIC, not ISM.  Both chains are in genomic order here; the caller flips
    miss5/miss3 for strand.
    """
    n, m = len(c), len(ref)
    if n == 0 or n >= m:
        return []
    out = []
    for k in range(m - n + 1):
        if all(
            _near(c[t][0], ref[k + t][0], tol) and _near(c[t][1], ref[k + t][1], tol)
            for t in range(n)
        ):
            out.append(k)
    return out


def subsequence_gaps(c: Chain, ref: Chain, tol: int):
    """``c`` as a possibly NON-contiguous subsequence of ``ref``.

    Returns ``(skipped_internal, miss_left, miss_right)`` or None.

    This exists because of intron retention, and getting it wrong is how IR
    ends up mislabelled.  SQANTI3 defines "ISM with Intron Retention" as a
    transcript "for which the loss of a SJ in the reference is due to an IR
    event" -- so an IR model has FEWER junctions than its reference, and the
    missing one sits in the MIDDLE.  Its chain is therefore not a contiguous
    subchain, and a classifier that only tests contiguity calls it NIC.

    The distinction that matters is what happened to the skipped reference
    junction.  If it lies inside one of the model's exons, the model read
    through it: intron retention, and ISM.  If the model splices straight over
    it, that is an exon-skip event and a genuinely novel combination: NIC.  The
    caller makes that check; this function only reports which junctions were
    skipped and where.
    """
    n, m = len(c), len(ref)
    if n == 0 or n >= m:
        return None
    i = 0
    matched: List[int] = []
    skipped: List[Tuple[int, Tuple[int, int]]] = []
    for r_idx, rj in enumerate(ref):
        if i < n and _near(c[i][0], rj[0], tol) and _near(c[i][1], rj[1], tol):
            matched.append(r_idx)
            i += 1
        else:
            skipped.append((r_idx, rj))
    if i != n or not matched:
        return None
    lo, hi = matched[0], matched[-1]
    left = sum(1 for r_idx, _ in skipped if r_idx < lo)
    right = sum(1 for r_idx, _ in skipped if r_idx > hi)
    internal = [j for r_idx, j in skipped if lo < r_idx < hi]
    return internal, left, right


def retained_introns(exons: Sequence[Exon], ref_chain: Chain) -> List[Tuple[int, int]]:
    """Annotated introns lying entirely inside one of the model's exons.

    Reported for every category.  An NNC with a retained intron stays NNC --
    this is a property of the model, not a competing classification.
    """
    out = []
    for d, a in ref_chain:
        for s, e in exons:
            if s <= d and a <= e:
                out.append((d, a))
                break
    return out


# ---------------------------------------------------------------------- #
def classify(
    chain: Chain,
    exons: Sequence[Exon],
    strand: str,
    ref,
    gene_ids: Sequence[str],
    tol: int,
    *,
    multi_gene: bool = False,
) -> Structural:
    """A2: FSM > ISM > NIC > NNC, evaluated in that order.

    ``gene_ids`` are the loci whose transcripts form the reference catalogue R.
    More than one **and** ``multi_gene`` means the model spans genes, which
    SQANTI3 calls ``fusion`` -- kept as a real category, not an artifact.

    The ISM test is per reference transcript, never over the union: a chain
    assembled from R1's 5' half and R2's 3' half is NIC, and calling it ISM
    would invent a transcript nobody annotated.
    """
    st = Structural()
    chain = genomic(chain)
    tids: List[str] = []
    for gid in gene_ids:
        g = ref.genes.get(gid)
        if g is not None:
            tids.extend(g.tx_ids)

    if not chain:
        # A single exon has no splice match to make, so "the gene exists" is not
        # grounds for full-splice_match: doing that put 22,127 mono-exon
        # fragments of BD144 kin -- 13% of its catalogue -- into FSM. SQANTI3
        # keeps mono-exon transcripts in the same taxonomy with a `mono-exon`
        # subcategory, but the category still has to be earned.
        st.subcategory = "mono-exon"
        st.primary = _primary(tids, ref)
        if not tids or not exons:
            st.structural_category = INTERGENIC if not tids else GENIC
            return st
        e = (exons[0][0], exons[-1][1])
        if any(not ref.tx[t].chain and _reciprocal(e, (ref.tx[t].exons[0][0],
                                                       ref.tx[t].exons[-1][1]), 0.5)
               for t in tids if ref.tx[t].exons):
            st.structural_category = FSM       # matches a mono-exon transcript
        elif any(s <= e[0] and e[1] <= x
                 for t in tids for s, x in ref.tx[t].exons):
            st.structural_category = ISM       # sits inside one reference exon
        else:
            st.structural_category = GENIC     # overlaps the gene, matches nothing
        return st

    # retained introns, for every category, against the whole catalogue
    seen = set()
    for tid in tids:
        for j in retained_introns(exons, genomic(ref.tx[tid].chain)):
            seen.add(j)
    st.retained_introns = tuple(sorted(seen))

    # --- 1. FSM
    fsm = [t for t in tids
           if chains_equal(chain, genomic(ref.tx[t].chain), tol)]
    if fsm:
        st.structural_category = FSM
        st.matches = tuple(fsm)
        st.primary = _primary(fsm, ref)
        return st

    # --- 2. ISM, per reference transcript
    ism: List[Tuple[str, int, int]] = []
    for t in tids:
        rc = genomic(ref.tx[t].chain)
        for k in subchain_offsets(chain, rc, tol):
            left = k                               # genomic-left introns missing
            right = len(rc) - len(chain) - k
            m5, m3 = (left, right) if strand == "+" else (right, left)
            ism.append((t, m5, m3))
    if ism:
        st.structural_category = ISM
        st.matches = tuple(t for t, _, _ in ism)
        st.primary = _primary([t for t, _, _ in ism], ref)
        best = min(ism, key=lambda x: (x[1] + x[2], x[1]))
        st.miss5, st.miss3 = best[1], best[2]
        st.subcategory = _ism_subcategory(st.miss5, st.miss3, False)
        return st

    # --- 2b. ISM with intron retention: a non-contiguous subsequence whose
    # every skipped INTERNAL reference junction is retained inside a model
    # exon. A skipped junction the model splices over instead is an exon skip,
    # which falls through to NIC below -- that is the whole difference between
    # "read through an intron" and "invented a new junction".
    exon_list = list(exons)
    ir_hits: List[Tuple[str, int, int]] = []
    for t in tids:
        got = subsequence_gaps(chain, genomic(ref.tx[t].chain), tol)
        if got is None:
            continue
        internal, left, right = got
        if not internal:
            continue                       # contiguous; step 2 already saw it
        if not all(retained_introns(exon_list, (j,)) for j in internal):
            continue                       # spliced over, not read through
        m5, m3 = (left, right) if strand == "+" else (right, left)
        ir_hits.append((t, m5, m3))
    if ir_hits:
        st.structural_category = ISM
        st.matches = tuple(t for t, _, _ in ir_hits)
        st.primary = _primary([t for t, _, _ in ir_hits], ref)
        best = min(ir_hits, key=lambda x: (x[1] + x[2], x[1]))
        st.miss5, st.miss3 = best[1], best[2]
        st.subcategory = "intron_retention"
        return st

    # --- 3/4. NIC vs NNC, by splice-site membership
    key = (ref.tx[tids[0]].contig, strand) if tids else None
    novel_sites = 0
    novel_combo = 0
    for d, a in chain:
        d_ok = key is not None and ref.is_annotated_donor(key, d)
        a_ok = key is not None and ref.is_annotated_acceptor(key, a)
        if not (d_ok and a_ok):
            novel_sites += 1
        elif key is not None and not ref.is_annotated_junction(key, (d, a)):
            novel_combo += 1
    st.novel_sites = novel_sites
    st.novel_combination_junctions = novel_combo

    if multi_gene and len({ref.tx[t].gene_id for t in tids}) > 1:
        st.structural_category = FUSION
        st.subcategory = "multi-exon"
        st.primary = _primary(tids, ref)
        return st

    if novel_sites == 0:
        st.structural_category = NIC
        # RP1 lives here: one unannotated junction, both of its sites annotated.
        # Today that call is NNC, which is what made a plain exon skip between
        # two known sites read as a novel splice site.
        st.subcategory = ("combination_of_known_junctions" if novel_combo == 0
                          else "combination_of_known_splicesites")
    else:
        st.structural_category = NNC
        st.subcategory = "at_least_one_novel_splicesite"
    st.primary = _primary(tids, ref)
    return st


def _ism_subcategory(miss5: int, miss3: int, has_ir: bool) -> str:
    """SQANTI3's ISM subcategories.

    ``has_ir`` is passed False from the contiguous branch on purpose. A model
    can carry a retained intron somewhere and still be a plain 3' fragment of
    its reference; ``intron_retention`` is for the case where the retention is
    what REMOVED the reference junction, which only the subsequence branch can
    establish. Keying the subcategory on "any retained intron anywhere" instead
    made every 5' fragment whose first exon happens to span an upstream intron
    read as IR.
    """
    if has_ir:
        return "intron_retention"
    if miss5 > 0 and miss3 > 0:
        return "internal_fragment"
    if miss3 > 0:
        return "3prime_fragment"
    return "5prime_fragment"


def _primary(tids: Sequence[str], ref) -> Optional[str]:
    """MANE > Ensembl canonical > basic > rest, via the existing rank."""
    if not tids:
        return None
    return min(tids, key=lambda t: (ref.tx[t].rank, t))


# ---------------------------------------------------------------------- #
def terminal_class(category: str, evidence: Dict[str, object]) -> str:
    """A3/A4's orthogonal 3'-end field, from columns that already exist.

    Deliberately NOT derived from `structural_category`: the whole point of
    splitting them is that "which chain is this" and "where does it end" are
    different questions, and overloading one category string with both is what
    produced thirteen flightcollapse categories for four SQANTI3 ones.
    """
    def _i(k: str) -> int:
        try:
            return int(float(evidence.get(k, 0) or 0))
        except (TypeError, ValueError):
            return 0

    if _i("end_is_annotated"):
        return "reference_end" if category == "FSM" else "alt_end_annotated"
    if _i("end_at_polya_site"):
        return "alt_end_atlas"
    if _i("internal_priming"):
        return "truncated_exonic"
    if category in ("monoexon_3UTR", "monoexon_internal", "monoexon_intergenic"):
        return "unresolved"
    return "unresolved"


def fsm_subcategory(category: str) -> str:
    """The legacy terminal categories, restated as SQANTI3 FSM subcategories.

    ``end5_novel`` is deliberately absent: pipeline.py assigns it to a chain
    that is a proper SUFFIX of an annotated one, which is ISM/5prime_fragment,
    not FSM. It is the one legacy name that does not simply move.
    """
    return {
        "FSM": "reference_match",
        "end3_annotated": "alternative_3end",
        "end3_novel": "alternative_3end",
        "end5_alt": "alternative_5end",
        "end5_extended": "alternative_5end",
    }.get(category, "")


# ---------------------------------------------------------------------- #
def annotate_models(models, ref, cfg) -> Dict[str, int]:
    """Set the SQANTI3 fields on every model. A pure post-pass.

    Runs AFTER everything else and writes only into ``m.evidence``: it never
    touches ``m.category``, ``m.reads`` or the model set. That is deliberate --
    the classification is a second opinion recorded beside the first, so the
    0.4.1 catalogue and the SQANTI3 one are the same models with two labels and
    the difference between them is measurable rather than asserted.

    An FSM whose legacy category already names its terminal difference keeps
    that information in ``subcategory``; ``end5_novel`` is the one legacy name
    that moves category, because pipeline.py assigns it to a proper suffix of an
    annotated chain, which is ISM.
    """
    tol = cfg.junctions.fuzzy_tolerance
    counts: Dict[str, int] = {}
    for m in models:
        # ``m.gene_id`` is the single best positional/chain assignment chosen
        # during collapse.  It cannot, by construction, describe a model whose
        # exons overlap more than one reference gene, which made the fusion
        # branch in ``classify`` unreachable through the production pipeline.
        # Resolve all same-strand genes with real exonic overlap for this
        # annotation-only pass.  ``classify`` still evaluates FSM and ISM before
        # fusion, so an ordinary model inside overlapping annotations keeps its
        # more specific match.
        ranked = ref.rank_genes_by_overlap(m.contig, m.exons, m.strand)
        gids = [gid for gid, overlap in ranked if overlap > 0]
        if m.gene_id and m.gene_id not in gids:
            gids.append(m.gene_id)
        multi_gene = len(gids) > 1
        try:
            st = classify(
                tuple(m.chain), m.exons, m.strand, ref, gids, tol,
                multi_gene=multi_gene,
            )
        except Exception:                      # annotation gaps must not fail a run
            continue
        if st.structural_category == FSM and not st.subcategory:
            st.subcategory = fsm_subcategory(m.category) or "reference_match"
        ev = m.evidence
        ev.update(st.as_columns())
        ev["structural_gene_ids"] = ",".join(gids)
        ev["structural_gene_count"] = len(gids)
        ev["terminal_class"] = terminal_class(m.category, ev)
        counts[st.structural_category] = counts.get(st.structural_category, 0) + 1
    return counts
