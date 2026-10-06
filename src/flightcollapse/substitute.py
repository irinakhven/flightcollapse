"""Rule J: a novel junction the short reads say did not happen.

New in 0.4.3.  Substitution, never deletion.

THE CASE.  At ARR3 a single novel donor 12 bp into the intron carries 35,244
long reads and **zero** unique short reads, while the two annotated junctions it
displaces -- across a 10 bp MANE microexon -- carry 111 and 109.  At BSG the
same shape: +19 bp, 30,684 long reads, 0 short reads, against 2,804 and 2,962
on the annotated pair.  A spliced aligner cannot seed inside a 10-25 bp exon, so
it extends the flanking exon a few bases and jumps the microexon.  The result is
highly supported, sharply positioned, canonical and *reproducible across
molecules*, which is why no support-, motif- or lFDR-based gate can see it.

WHY NOT SIMPLY DROP JUNCTIONS WITH NO SHORT-READ SUPPORT.  Measured on BD67's
68,563 kept novel junctions, ``sr_uniq == 0`` covers 46,393 of them carrying
808,434 long reads -- **46.7% of all novel-junction read mass** -- and it is
that high in every class (64-76%), because the short-read library is a
different, shallower organoid batch.  Zero short reads is the normal state of a
novel junction here, not a verdict.  What makes ARR3 and BSG decisive is the
CONJUNCTION: zero where many were expected, AND a well-covered annotated
competitor proving the locus was visible to the short reads at all.

Dropping would also lose the reads rather than move them.  A read carrying a
rejected junction goes to the decoy set and never reaches a model -- 735,864
reads on BD176c already take that path.  The 31k ARR3 UMIs belong on MANE, so
the operation has to be a substitution.

THE RULE.  A chain group merges into an annotated group when, for one reference
transcript T:

  1. every junction of the group either matches a junction of T within TOL, or
     is a single junction spanning a contiguous RUN of T's junctions with one
     end anchored (within TOL) and the other displaced into the intron by at
     most ``max_shift``;
  2. every displaced junction has ``sr_uniq <= max_novel_sr``;
  3. every reference junction it displaces has ``sr_uniq >= min_competitor_sr``;
  4. at least one junction is actually displaced -- otherwise this is an
     ordinary FSM and nothing to do;
  5. a surviving group with T's exact chain exists to merge into.

(5) is the part that keeps this honest.  Without a destination there is no
defensible place for the reads, so the model is left alone and only flagged.

WHAT IT DELIBERATELY DOES NOT TOUCH.  An exon skip whose two splice sites are
both annotated is ordinary alternative splicing; 270 such junctions on BD67 are
short-read-depleted and 100% canonical, and RP1 -- a real 172 bp skip with 21
short reads -- lives in that class.  Condition 1 excludes them, because a skip
displaces no site: both of its ends sit exactly on annotated sites.

AND IT IS INERT WITHOUT SHORT READS.  Conditions 2 and 3 are unsatisfiable when
no SJ table was supplied, so on HEK, mouse brain or SIRVs this does nothing at
all.  That is the intended behaviour, not a gap: without a second library there
is no evidence that separates an aligner artifact from a real novel junction.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .intervals import Chain
from .structural import genomic

Junction = Tuple[int, int]


@dataclass
class Substitution:
    """One model junction replaced by a run of reference junctions."""
    model_index: int
    junction: Junction
    replaces: Tuple[Junction, ...]
    shift: int
    novel_sr: int
    competitor_sr: int

    @property
    def swallowed(self) -> int:
        """Reference exons the junction jumped over. 0 for a plain site shift."""
        return max(len(self.replaces) - 1, 0)


def _near(a: int, b: int, tol: int) -> bool:
    return abs(a - b) <= tol


def derive_from_reference(
    chain: Chain,
    ref_chain: Chain,
    tol: int,
    max_shift: int,
) -> Optional[List[Optional[Tuple[int, int, int]]]]:
    """Can ``chain`` be read as ``ref_chain`` with some junctions displaced?

    Returns one entry per model junction: ``None`` for an exact match, or
    ``(ref_lo, ref_hi, shift)`` for a junction that spans ``ref_chain[lo:hi]``
    with one end anchored and the other displaced by ``shift`` bases into the
    intron.  ``None`` overall when no such reading exists.

    Greedy left to right, which is sound because junctions within a chain are
    disjoint and sorted: a model junction can only ever consume a prefix of the
    remaining reference junctions, and consuming the shortest one that works
    never blocks a later match.
    """
    if not chain or not ref_chain or len(chain) > len(ref_chain):
        return None
    ops: List[Optional[Tuple[int, int, int]]] = []
    i = r = 0
    while i < len(chain) and r < len(ref_chain):
        md, ma = chain[i]
        rd, ra = ref_chain[r]
        if _near(md, rd, tol) and _near(ma, ra, tol):
            ops.append(None)
            i += 1
            r += 1
            continue
        hit = None
        for s in range(r, len(ref_chain)):
            sd, sa = ref_chain[s]
            # donor anchored, acceptor pulled back into the intron
            if _near(md, rd, tol) and 0 < (sa - ma) <= max_shift:
                hit = (r, s + 1, int(sa - ma))
            # acceptor anchored, donor pushed forward into the intron
            elif _near(ma, sa, tol) and 0 < (md - rd) <= max_shift:
                hit = (r, s + 1, int(md - rd))
            if hit is not None:
                # a run of one is a plain site shift; a longer run means the
                # junction jumped over s-r reference exons, which is the
                # microexon case
                ops.append(hit)
                i += 1
                r = s + 1
                break
        if hit is None:
            return None
    if i != len(chain) or r != len(ref_chain):
        return None
    return ops


def plan_for_chain(
    chain: Chain,
    ref_chain: Chain,
    strand: str,
    sr: Dict[Tuple[str, int, int], int],
    tol: int,
    max_shift: int,
    max_novel_sr: int,
    min_competitor_sr: int,
) -> Optional[List[Substitution]]:
    """The substitutions needed to turn ``chain`` into ``ref_chain``, or None.

    None when the chains do not relate that way, when nothing is displaced (an
    ordinary FSM), or when any displaced junction fails the short-read test.
    A partial plan is never returned: moving a model's reads onto a reference
    chain is only defensible if EVERY difference is explained.
    """
    ops = derive_from_reference(chain, ref_chain, tol, max_shift)
    if ops is None:
        return None
    out: List[Substitution] = []
    for idx, op in enumerate(ops):
        if op is None:
            continue
        lo, hi, shift = op
        j = (int(chain[idx][0]), int(chain[idx][1]))
        n_sr = sr.get((strand, j[0], j[1]))
        if n_sr is None or n_sr > max_novel_sr:
            return None
        replaced = tuple((int(d), int(a)) for d, a in ref_chain[lo:hi])
        comp = [sr.get((strand, d, a)) for d, a in replaced]
        if any(c is None or c < min_competitor_sr for c in comp):
            return None
        out.append(Substitution(idx, j, replaced, shift, int(n_sr), int(min(comp))))
    return out or None


# ---------------------------------------------------------------------- #
def plan_substitutions(
    groups: Sequence,
    candidate,
    ref,
    sr: Dict[Tuple[str, int, int], int],
    params,
    tol: int,
) -> Tuple[List[Tuple[int, int, List[Substitution]]], Dict[str, int]]:
    """Which groups should merge into which annotated group.

    Returns ``(plans, stats)`` without mutating anything -- the caller applies
    them, so the decision and the edit stay separable and the plan can be
    written out and checked.
    """
    stats = {
        "groups_considered": 0,
        "groups_substituted": 0,
        "reads_substituted": 0,
        "junctions_substituted": 0,
        "junctions_swallowing_exons": 0,
        "no_destination": 0,
        "failed_short_read_test": 0,
    }
    if not sr:
        return [], stats

    # surviving annotated groups, by exact chain
    dest: Dict[Tuple[str, Chain], int] = {}
    for g in groups:
        if g.merged_into is None and g.chain and getattr(g, "annotated", False):
            dest.setdefault((g.strand, tuple(g.chain)), g.gid)

    plans: List[Tuple[int, int, List[Substitution]]] = []
    for g in groups:
        if g.merged_into is not None or not g.chain or g.n_novel_junctions == 0:
            continue
        if not candidate[g.gid]:
            continue
        stats["groups_considered"] += 1
        tids: List[str] = []
        if g.gene_id and g.gene_id in ref.genes:
            tids = list(ref.genes[g.gene_id].tx_ids)
        best: Optional[Tuple[int, List[Substitution]]] = None
        saw_shape = False
        for tid in tids:
            rc = genomic(ref.tx[tid].chain)
            if derive_from_reference(genomic(g.chain), rc, tol, params.max_shift) is None:
                continue
            saw_shape = True
            target = dest.get((g.strand, tuple(rc)))
            if target is None or target == g.gid:
                continue
            plan = plan_for_chain(
                genomic(g.chain), rc, g.strand, sr, tol, params.max_shift,
                params.max_novel_sr, params.min_competitor_sr,
            )
            if plan is None:
                continue
            # prefer the reading that displaces fewest junctions
            if best is None or len(plan) < len(best[1]):
                best = (target, plan)
        if best is None:
            if saw_shape:
                stats["no_destination"] += 1
            continue
        target, plan = best
        plans.append((g.gid, target, plan))
        stats["groups_substituted"] += 1
        stats["reads_substituted"] += int(g.n_reads)
        stats["junctions_substituted"] += len(plan)
        stats["junctions_swallowing_exons"] += sum(1 for s in plan if s.swallowed)
    return plans, stats


def apply_substitutions(groups: Sequence, plans) -> None:
    """Point each planned group at its destination.

    Uses the existing ``merged_into`` / ``merge_class`` machinery, so read
    accounting, `exact_chain_frac` and the group-chains table all handle these
    reads the same way they handle a suffix merge -- including recording that
    the moved reads do NOT carry the destination's intron chain, which is true
    and must stay visible.
    """
    for gid, target, _plan in plans:
        groups[gid].merged_into = target
        groups[gid].merge_class = "sr_substituted"


def write_plan(plans, groups, path: str) -> None:
    """One row per substituted junction. The audit trail for a read move."""
    with open(path, "w") as fh:
        fh.write("contig\tstrand\tgid\ttarget_gid\tn_reads\tdonor\tacceptor\t"
                 "shift\texons_swallowed\tnovel_sr\tcompetitor_sr\treplaced_by\n")
        for gid, target, plan in plans:
            g = groups[gid]
            for s in plan:
                fh.write(
                    f"{g.contig}\t{g.strand}\t{gid}\t{target}\t{g.n_reads}\t"
                    f"{s.junction[0]}\t{s.junction[1]}\t{s.shift}\t{s.swallowed}\t"
                    f"{s.novel_sr}\t{s.competitor_sr}\t"
                    + ";".join(f"{d}-{a}" for d, a in s.replaces) + "\n"
                )
