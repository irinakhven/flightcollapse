"""Part A of the NIC/ISM spec: the A6 fixtures, one per rule.

The regression that matters is #7 (RP1): one unannotated junction whose two
sites are both annotated. flightcollapse called that NNC, because
`_novel_category` asked only "is there a junction I have not seen" and never
"are its sites known". On BD67 that mislabelled 25.7% of NNC UMI mass.

#3 (RTN4) is the other one: twelve NIC models at that locus are the same
single-intron fragment of RTN4-202. Here it must come back ISM, and it must be
EMITTED with a label rather than dropped -- the classifier does not filter.
"""
from __future__ import annotations

import pytest

from flightcollapse.structural import (FSM, FUSION, ISM, NIC, NNC, annotate_models, classify,
                                       chains_equal, retained_introns,
                                       subchain_offsets, terminal_class,
                                       fsm_subcategory)

TOL = 5


class FakeTx:
    def __init__(self, tid, chain, strand="+", rank=0, gene_id="G", contig="chr1"):
        self.tx_id = tid
        # ReferenceIndex stores the chain in TRANSCRIPT orientation --
        # exons_to_chain reverses it on the minus strand -- so the fake has to
        # as well, or the strand tests pass against a fixture the real
        # reference never produces.
        self.chain = tuple(chain)[::-1] if strand == "-" else tuple(chain)
        self.strand = strand
        self.rank = rank
        self.gene_id = gene_id
        self.contig = contig
        # exon list implied by the chain, so the mono-exon path has spans to test
        _g = tuple(sorted(self.chain))
        self.exons = ([] if not _g else
                      [(_g[0][0] - 100, _g[0][0])]
                      + [(a, d) for (_, a), (d, _) in zip(_g[:-1], _g[1:])]
                      + [(_g[-1][1], _g[-1][1] + 100)])


class FakeGene:
    def __init__(self, tx_ids):
        self.tx_ids = list(tx_ids)


class FakeRef:
    """Only the four things `classify` touches."""

    def __init__(self, txs, genes=None):
        self.tx = {t.tx_id: t for t in txs}
        gid = txs[0].gene_id if txs else "G"
        self.genes = genes or {gid: FakeGene([t.tx_id for t in txs])}
        self._d = {p for t in txs for p, _ in t.chain}
        self._a = {p for t in txs for _, p in t.chain}
        self._j = {j for t in txs for j in t.chain}

    def is_annotated_donor(self, key, p):
        return any(abs(p - x) <= TOL for x in self._d)

    def is_annotated_acceptor(self, key, p):
        return any(abs(p - x) <= TOL for x in self._a)

    def is_annotated_junction(self, key, j):
        return any(abs(j[0] - d) <= TOL and abs(j[1] - a) <= TOL for d, a in self._j)


# A reference transcript with 4 introns, and a second isoform sharing its ends.
R1 = [(100, 200), (300, 400), (500, 600), (700, 800)]
R2 = [(100, 200), (300, 400), (505, 600), (700, 800)]   # alt donor, still annotated


def ref_std(strand="+"):
    return FakeRef([FakeTx("R1", R1, strand, rank=0),
                    FakeTx("R2", R2, strand, rank=1)])


def ex(chain, pad=40):
    """Exons implied by a chain, with realistic terminal exons.

    The first exon starts just before the first donor rather than at a fixed
    coordinate. That matters: a first exon anchored far upstream would span the
    reference introns the fragment is missing, and those would then read as
    retained -- which is exactly the false IR call this fixture must not make.
    """
    out, prev = [], chain[0][0] - pad
    for d, a in chain:
        out.append((prev, d))
        prev = a
    out.append((prev, chain[-1][1] + pad))
    return out


# --------------------------------------------------------------------- 1
def test_fsm_exact_and_at_tolerance():
    r = ref_std()
    s = classify(tuple(R1), ex(R1), "+", r, ["G"], TOL)
    assert s.structural_category == FSM
    jit = [(100, 200), (303, 400), (500, 600), (700, 800)]     # +3, inside TOL
    assert classify(tuple(jit), ex(jit), "+", r, ["G"], TOL).structural_category == FSM


def test_one_bp_past_tolerance_is_not_fsm():
    r = ref_std()
    off = [(100, 200), (306, 400), (500, 600), (700, 800)]     # +6, outside TOL
    s = classify(tuple(off), ex(off), "+", r, ["G"], TOL)
    assert s.structural_category != FSM
    # both sites are still annotated (R2 has a 505 donor, 306 is not near it),
    # so this is a genuinely novel donor
    assert s.structural_category == NNC


# --------------------------------------------------------------------- 2,3
def test_ism_5prime_fragment():
    """First introns dropped, 3' side intact -- RT drop-off, the expected shape."""
    r = ref_std()
    c = tuple(R1[2:])
    s = classify(c, ex(c), "+", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert (s.miss5, s.miss3) == (2, 0)
    assert s.subcategory == "5prime_fragment"


def test_ism_3prime_fragment_is_emitted_not_dropped():
    """RTN4's shape. The classifier LABELS it; it never filters."""
    r = ref_std()
    c = tuple(R1[:2])
    s = classify(c, ex(c), "+", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert (s.miss5, s.miss3) == (0, 2)
    assert s.subcategory == "3prime_fragment"
    assert s.primary in ("R1", "R2")


def test_ism_internal_fragment():
    r = ref_std()
    c = (R1[1], R1[2])
    s = classify(c, ex(c), "+", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert (s.miss5, s.miss3) == (1, 1)
    assert s.subcategory == "internal_fragment"


# --------------------------------------------------------------------- 6
def test_union_chain_is_nic_never_ism():
    """R1's 5' half plus R2's 3' half. The ISM test is per transcript."""
    r = ref_std()
    c = (R1[0], R1[1], R2[2], R2[3])
    # that IS R2 exactly, so build a real union instead
    R3 = [(100, 200), (300, 400), (500, 600), (700, 800)]
    r2 = FakeRef([FakeTx("A", [(100, 200), (300, 400)], rank=0),
                  FakeTx("B", [(500, 600), (700, 800)], rank=1)])
    u = tuple(R3)
    s = classify(u, ex(u), "+", r2, ["G"], TOL)
    assert s.structural_category == NIC
    assert s.subcategory == "combination_of_known_junctions"


# --------------------------------------------------------------------- 7
def test_rp1_exon_skip_is_nic_not_nnc():
    """THE regression. An unannotated junction from two annotated sites."""
    r = ref_std()
    skip = (R1[0], (300, 600), R1[3])       # donor 300 known, acceptor 600 known
    s = classify(skip, ex(skip), "+", r, ["G"], TOL)
    assert s.structural_category == NIC, "this is today's NNC call"
    assert s.subcategory == "combination_of_known_splicesites"
    assert s.novel_sites == 0
    assert s.novel_combination_junctions == 1


def test_genuinely_novel_site_is_still_nnc():
    r = ref_std()
    c = (R1[0], (317, 400), R1[2], R1[3])   # donor 317 is near nothing annotated
    s = classify(c, ex(c), "+", r, ["G"], TOL)
    assert s.structural_category == NNC
    assert s.novel_sites == 1


# --------------------------------------------------------------------- 8
def test_single_intron_inside_a_longer_reference_is_ism():
    r = ref_std()
    c = (R1[2],)
    s = classify(c, ex(c), "+", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert (s.miss5, s.miss3) == (2, 1)


# --------------------------------------------------------------------- 9
def test_minus_strand_flips_miss5_and_miss3():
    """The same genomic chain, read the other way, must swap the two."""
    r = ref_std("-")
    c = tuple(R1[2:])
    plus = classify(c, ex(c), "+", ref_std("+"), ["G"], TOL)
    minus = classify(c, ex(c), "-", r, ["G"], TOL)
    assert (plus.miss5, plus.miss3) == (2, 0)
    assert (minus.miss5, minus.miss3) == (0, 2)
    assert plus.subcategory == "5prime_fragment"
    assert minus.subcategory == "3prime_fragment"


# --------------------------------------------------------------------- 10
def test_intron_retention_is_an_ism_subcategory_not_a_drop():
    """FSM chain minus one intron, that intron retained inside an exon.

    SQANTI3 calls this ISM / intron_retention. It stays a real model with a
    real category -- IR is never an artifact verdict here.
    """
    r = ref_std()
    c = (R1[0], R1[1], R1[3])               # (500,600) read through
    exons = [(50, 100), (200, 300), (400, 700), (800, 900)]
    s = classify(c, exons, "+", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert s.subcategory == "intron_retention"
    assert (500, 600) in s.retained_introns


def test_retained_introns_reported_for_nnc_too():
    """An NNC with a retained intron stays NNC; the flag is orthogonal."""
    r = ref_std()
    c = ((100, 200), (317, 400), (700, 800))
    exons = [(50, 100), (200, 317), (400, 700), (800, 900)]
    s = classify(c, exons, "+", r, ["G"], TOL)
    assert s.structural_category == NNC
    assert (500, 600) in s.retained_introns


# --------------------------------------------------------------------- fusion
def test_multi_gene_is_fusion_and_is_a_real_category():
    a = FakeTx("A", [(100, 200)], gene_id="G1")
    b = FakeTx("B", [(700, 800)], gene_id="G2")
    r = FakeRef([a, b], genes={"G1": FakeGene(["A"]), "G2": FakeGene(["B"])})
    c = ((100, 200), (700, 800))
    s = classify(c, ex(c), "+", r, ["G1", "G2"], TOL, multi_gene=True)
    assert s.structural_category == FUSION


def test_production_annotation_supplies_all_overlapping_genes_for_fusion():
    """Regression for 0.4.5: annotate_models passed exactly one gene id, so
    classify's fusion branch was correct but unreachable in the pipeline.
    """
    from types import SimpleNamespace

    a = FakeTx("A", [(100, 200)], gene_id="G1")
    b = FakeTx("B", [(700, 800)], gene_id="G2")
    r = FakeRef([a, b], genes={"G1": FakeGene(["A"]), "G2": FakeGene(["B"])})
    r.rank_genes_by_overlap = lambda contig, exons, strand: [("G1", 100), ("G2", 100)]
    c = ((100, 200), (700, 800))
    model = SimpleNamespace(
        contig="chr1", strand="+", chain=c, exons=ex(c), gene_id="G1",
        category="NNC", evidence={},
    )
    cfg = SimpleNamespace(junctions=SimpleNamespace(fuzzy_tolerance=TOL))

    counts = annotate_models([model], r, cfg)

    assert model.evidence["structural_category"] == FUSION
    assert model.evidence["structural_gene_ids"] == "G1,G2"
    assert model.evidence["structural_gene_count"] == 2
    assert counts == {FUSION: 1}


# --------------------------------------------------------------------- helpers
def test_subchain_offsets_requires_contiguity():
    ref = tuple(R1)
    assert subchain_offsets(((100, 200), (300, 400)), ref, TOL) == [0]
    # skipping ref[1] is an exon-skip event, which is NIC, not ISM
    assert subchain_offsets(((100, 200), (500, 600)), ref, TOL) == []


def test_subchain_is_proper():
    ref = tuple(R1)
    assert subchain_offsets(ref, ref, TOL) == [], "equal chains are FSM, not ISM"


def test_chains_equal_needs_same_length():
    assert not chains_equal(tuple(R1[:2]), tuple(R1), TOL)


def test_retained_introns_needs_full_containment():
    assert retained_introns([(400, 700)], ((500, 600),)) == [(500, 600)]
    assert retained_introns([(400, 550)], ((500, 600),)) == []


def test_legacy_terminal_categories_are_all_fsm_subcategories():
    """The second over-call: four category names for one SQANTI3 category."""
    for c in ("FSM", "end3_annotated", "end3_novel", "end5_alt", "end5_extended"):
        assert fsm_subcategory(c) != ""
    # end5_novel is a proper SUFFIX of an annotated chain -> ISM, not FSM
    assert fsm_subcategory("end5_novel") == ""


def test_terminal_class_is_independent_of_the_chain():
    assert terminal_class("FSM", {"end_is_annotated": 1}) == "reference_end"
    assert terminal_class("NIC", {"end_is_annotated": 1}) == "alt_end_annotated"
    assert terminal_class("NNC", {"end_at_polya_site": 1}) == "alt_end_atlas"
    assert terminal_class("NIC", {"internal_priming": 1}) == "truncated_exonic"
    assert terminal_class("NIC", {}) == "unresolved"


def test_same_chain_ir_vs_exon_skip_depends_on_the_exons():
    """The crux, and the bug this fixture caught during development.

    Chain (100,200),(300,400),(700,800) is missing R1's (500,600). Whether that
    is ISM/intron_retention or NIC is not a property of the chain at all -- it
    is whether the model's exons cover the skipped junction. A classifier that
    only tests contiguous subchains calls both NIC.
    """
    r = ref_std()
    c = ((100, 200), (300, 400), (700, 800))
    read_through = [(60, 100), (200, 300), (400, 700), (800, 840)]
    spliced_over = [(60, 100), (200, 300), (400, 700), (800, 840)]
    s_ir = classify(c, read_through, "+", r, ["G"], TOL)
    assert s_ir.structural_category == ISM
    assert s_ir.subcategory == "intron_retention"
    # now make the model splice over (500,600) instead: the exon must not cover
    # it, which means a real skip junction, i.e. a different chain
    skip = ((100, 200), (300, 400), (500, 800))      # acceptor 800, donor 500
    s_sk = classify(skip, ex(skip), "+", r, ["G"], TOL)
    assert s_sk.structural_category == NIC
    assert s_sk.subcategory == "combination_of_known_splicesites"


# ------------------------------------------------------------- mono-exon
def test_monoexon_does_not_become_fsm_just_because_the_gene_exists():
    """The bug this caught: 22,127 single-exon fragments of BD144 kin -- 13% of
    the catalogue -- were landing in full-splice_match because their gene had
    transcripts. A single exon has no splice match to make."""
    from flightcollapse.structural import GENIC, INTERGENIC
    r = ref_std()
    inside = classify((), [(210, 290)], "+", r, ["G"], TOL)   # inside exon 200-300
    assert inside.structural_category == ISM
    assert inside.subcategory == "mono-exon"
    # spanning an intron, matching no mono-exon transcript and inside no exon
    across = classify((), [(150, 350)], "+", r, ["G"], TOL)
    assert across.structural_category == GENIC
    assert classify((), [(10, 20)], "+", r, [], TOL).structural_category == INTERGENIC


def test_monoexon_matching_a_monoexon_transcript_is_fsm():
    mono = FakeTx("M", [], rank=0)
    mono.exons = [(1000, 2000)]
    r = FakeRef([mono])
    r.tx["M"].chain = ()
    s = classify((), [(1010, 1990)], "+", r, ["G"], TOL)
    assert s.structural_category == FSM
    assert s.subcategory == "mono-exon"


# ------------------------------------------------- orientation, both ways
@pytest.mark.parametrize("orient", ["transcript", "genomic"])
def test_minus_strand_is_correct_whichever_order_the_chain_arrives_in(orient):
    """The trap that bit twice, in opposite directions.

    ContigReads.chain reverses on the minus strand and so does exons_to_chain,
    so the pipeline hands classify() TRANSCRIPT-oriented chains on both sides --
    which is why the collapser's own FSM lookup was always strand-balanced
    (BD144 kin: 56.5% of plus-strand multi-exon models carry a parent_transcript
    against 55.7% of minus). A post-hoc caller that re-sorts hands it GENOMIC
    ones. Both must give the same answer.
    """
    r = ref_std("-")
    g = tuple(R1)                     # genomic
    t = g[::-1]                       # transcript orientation, minus strand
    c = t if orient == "transcript" else g
    assert classify(c, ex(g), "-", r, ["G"], TOL).structural_category == FSM


@pytest.mark.parametrize("orient", ["transcript", "genomic"])
def test_minus_strand_miss5_miss3_are_not_double_flipped(orient):
    """0.4.2's real minus-strand bug: miss5/miss3 swapped.

    It flipped by strand on top of an already transcript-oriented chain. The
    genomic-left introns are missing here; on the minus strand those are the
    3' side, so miss3 must be 2 and miss5 zero.
    """
    r = ref_std("-")
    g = tuple(R1[2:])
    c = g[::-1] if orient == "transcript" else g
    s = classify(c, ex(g), "-", r, ["G"], TOL)
    assert s.structural_category == ISM
    assert (s.miss5, s.miss3) == (0, 2)
    assert s.subcategory == "3prime_fragment"


def test_genomic_is_idempotent():
    from flightcollapse.structural import genomic
    g = tuple(R1)
    assert genomic(g) == g == genomic(g[::-1]) == genomic(genomic(g[::-1]))
