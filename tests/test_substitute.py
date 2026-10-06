"""Rule J: ARR3 and BSG geometry, and the things that must NOT move.

The two real cases, as coordinates:

  ARR3  donor +12 bp into the intron, jumping a 10 bp MANE microexon.
        novel jx 35,244 LR / 0 SR; the two annotated jx 111 and 109 SR.
  BSG   donor +19 bp, jumping a 25 bp microexon.
        novel jx 30,684 LR / 0 SR; annotated pair 2,804 and 2,962 SR.

The negative cases matter more than the positive ones. `sr_uniq == 0` alone
covers 46.7% of BD67's novel-junction read mass, so every test here that
asserts "nothing happens" is guarding against the version of this rule that
deletes half the catalogue.
"""
from __future__ import annotations

import pytest

from flightcollapse.substitute import (Substitution, apply_substitutions,
                                       derive_from_reference, plan_for_chain,
                                       plan_substitutions)

TOL, SHIFT = 5, 30

# A reference transcript with a microexon: exons ... 200-300, [310-320], 400-...
# so its chain is (100,200) (300,310) (320,400). The 10 bp exon is 310..320.
REF_MICRO = ((100, 200), (300, 310), (320, 400))
# the aligner's version: one junction from 300+12 straight to 400
ARR3_LIKE = ((100, 200), (312, 400))

REF_PLAIN = ((100, 200), (300, 400), (500, 600))


def sr_table(**kw):
    """(strand, donor, acceptor) -> unique short reads."""
    return {("+",) + tuple(int(x) for x in k.split("_")): v for k, v in kw.items()}


class P:
    max_shift = SHIFT
    max_novel_sr = 0
    min_competitor_sr = 10


# --------------------------------------------------------------- geometry
def test_arr3_geometry_is_derivable():
    ops = derive_from_reference(ARR3_LIKE, REF_MICRO, TOL, SHIFT)
    assert ops is not None
    assert ops[0] is None, "first junction matches exactly"
    lo, hi, shift = ops[1]
    assert (lo, hi) == (1, 3), "one junction spans both reference junctions"
    assert shift == 12


def test_shift_beyond_max_is_not_derivable():
    far = ((100, 200), (340, 400))          # +40, past max_shift
    assert derive_from_reference(far, REF_MICRO, TOL, SHIFT) is None


def test_plain_site_shift_onto_one_reference_junction():
    c = ((100, 200), (309, 400), (500, 600))
    ops = derive_from_reference(c, REF_PLAIN, TOL, SHIFT)
    assert ops is not None
    assert ops[1] == (1, 2, 9), "one reference junction, donor displaced 9 bp"


def test_exact_chain_derives_with_no_displacement():
    ops = derive_from_reference(REF_PLAIN, REF_PLAIN, TOL, SHIFT)
    assert ops == [None, None, None]


def test_longer_than_reference_never_derives():
    assert derive_from_reference(REF_PLAIN, ((100, 200),), TOL, SHIFT) is None


# --------------------------------------------------------------- the rule
def test_arr3_substitutes_when_the_competitor_is_covered():
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 111, "320_400": 109})
    plan = plan_for_chain(ARR3_LIKE, REF_MICRO, "+", sr, TOL, SHIFT, 0, 10)
    assert plan is not None and len(plan) == 1
    s = plan[0]
    assert s.junction == (312, 400)
    assert s.replaces == ((300, 310), (320, 400))
    assert (s.shift, s.swallowed, s.novel_sr, s.competitor_sr) == (12, 1, 0, 109)


def test_bsg_geometry_with_a_25bp_microexon():
    ref = ((100, 200), (300, 325), (350, 400))      # 25 bp exon 325..350
    obs = ((100, 200), (319, 400))                  # +19
    sr = sr_table(**{"100_200": 500, "319_400": 0, "300_325": 2804, "350_400": 2962})
    plan = plan_for_chain(obs, ref, "+", sr, TOL, SHIFT, 0, 10)
    assert plan is not None and plan[0].shift == 19


# ------------------------------------------------- what must NOT substitute
def test_no_substitution_when_the_competitor_has_no_short_reads():
    """The locus was never visible to the short reads, so 0 proves nothing.

    This is the guard against the rule that deletes 46.7% of novel read mass.
    """
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 0, "320_400": 0})
    assert plan_for_chain(ARR3_LIKE, REF_MICRO, "+", sr, TOL, SHIFT, 0, 10) is None


def test_no_substitution_when_the_competitor_is_thinly_covered():
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 3, "320_400": 2})
    assert plan_for_chain(ARR3_LIKE, REF_MICRO, "+", sr, TOL, SHIFT, 0, 10) is None


def test_no_substitution_when_the_novel_junction_has_short_reads():
    """RP1's protection: a novel junction the short reads DID see stays."""
    sr = sr_table(**{"100_200": 90, "312_400": 21, "300_310": 111, "320_400": 109})
    assert plan_for_chain(ARR3_LIKE, REF_MICRO, "+", sr, TOL, SHIFT, 0, 10) is None


def test_exon_skip_between_two_annotated_sites_is_not_derivable():
    """Both ends sit exactly on annotated sites, so nothing is displaced.

    270 junctions on BD67 are short-read-depleted exon skips, 100% canonical,
    and RP1 is one of that shape. Condition 1 has to exclude them structurally,
    not by a threshold.
    """
    skip = ((100, 200), (300, 600))         # donor 300 and acceptor 600 both real
    ops = derive_from_reference(skip, REF_PLAIN, TOL, SHIFT)
    assert ops is None or all(o is None or o[2] == 0 for o in ops)
    sr = sr_table(**{"100_200": 90, "300_600": 0, "300_400": 50, "500_600": 50})
    assert plan_for_chain(skip, REF_PLAIN, "+", sr, TOL, SHIFT, 0, 10) is None


def test_an_exact_match_yields_no_plan():
    sr = sr_table(**{"100_200": 90, "300_400": 90, "500_600": 90})
    assert plan_for_chain(REF_PLAIN, REF_PLAIN, "+", sr, TOL, SHIFT, 0, 10) is None


def test_missing_short_read_entry_is_not_treated_as_zero():
    """A junction absent from the SJ tables is unmeasured, not unsupported."""
    sr = sr_table(**{"100_200": 90, "300_310": 111, "320_400": 109})   # no 312_400
    assert plan_for_chain(ARR3_LIKE, REF_MICRO, "+", sr, TOL, SHIFT, 0, 10) is None


# ------------------------------------------------------------ group level
class G:
    def __init__(self, gid, chain, annotated=False, novel=0, reads=10,
                 gene_id="G", strand="+", contig="chr1"):
        self.gid = gid
        self.chain = tuple(chain)
        self.annotated = annotated
        self.n_novel_junctions = novel
        self.n_reads = reads
        self.gene_id = gene_id
        self.strand = strand
        self.contig = contig
        self.merged_into = None
        self.merge_class = ""


class Tx:
    def __init__(self, chain):
        self.chain = tuple(chain)


class Ref:
    def __init__(self, chains):
        self.tx = {f"T{i}": Tx(c) for i, c in enumerate(chains)}
        class Gene:
            tx_ids = list(self.tx)
        self.genes = {"G": Gene()}


def test_group_merges_into_the_annotated_group():
    groups = [G(0, ARR3_LIKE, novel=1, reads=35244),
              G(1, REF_MICRO, annotated=True, reads=338)]
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 111, "320_400": 109})
    plans, st = plan_substitutions(groups, [True, True], Ref([REF_MICRO]), sr, P, TOL)
    assert st["groups_substituted"] == 1
    assert st["reads_substituted"] == 35244
    assert st["junctions_swallowing_exons"] == 1
    apply_substitutions(groups, plans)
    assert groups[0].merged_into == 1
    assert groups[0].merge_class == "sr_substituted"
    assert groups[1].merged_into is None


def test_no_destination_leaves_the_model_alone():
    """Condition 5. Without a surviving annotated group there is nowhere
    defensible for the reads, so nothing moves and it is counted."""
    groups = [G(0, ARR3_LIKE, novel=1, reads=35244)]
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 111, "320_400": 109})
    plans, st = plan_substitutions(groups, [True], Ref([REF_MICRO]), sr, P, TOL)
    assert plans == []
    assert st["no_destination"] == 1
    assert groups[0].merged_into is None


def test_inert_without_short_reads():
    """HEK, mouse brain and SIRVs must be untouched."""
    groups = [G(0, ARR3_LIKE, novel=1), G(1, REF_MICRO, annotated=True)]
    plans, st = plan_substitutions(groups, [True, True], Ref([REF_MICRO]), {}, P, TOL)
    assert plans == [] and st["groups_substituted"] == 0
    assert groups[0].merged_into is None


def test_already_merged_groups_are_skipped():
    groups = [G(0, ARR3_LIKE, novel=1), G(1, REF_MICRO, annotated=True)]
    groups[0].merged_into = 1
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 111, "320_400": 109})
    plans, _ = plan_substitutions(groups, [True, True], Ref([REF_MICRO]), sr, P, TOL)
    assert plans == []


def test_plan_is_written_before_it_is_applied():
    """The decision and the edit are separable, so the move can be audited."""
    groups = [G(0, ARR3_LIKE, novel=1, reads=35244),
              G(1, REF_MICRO, annotated=True)]
    sr = sr_table(**{"100_200": 90, "312_400": 0, "300_310": 111, "320_400": 109})
    plans, _ = plan_substitutions(groups, [True, True], Ref([REF_MICRO]), sr, P, TOL)
    assert groups[0].merged_into is None, "planning must not mutate"
    apply_substitutions(groups, plans)
    assert groups[0].merged_into == 1
