"""Final consolidation of the emitted catalogue.  New in 0.6.0.

Two rules, applied per contig after the second pass and before read
bookkeeping, so the read assignments, ``group.txt``, ``read_stat`` and the
per-cell matrix all describe the catalogue that is actually written.  The order
is fixed -- 3' consolidation, then the fragment filter -- so the support ratio
in the second rule is tested on consolidated counts.

Support throughout is ``n_mols`` when it is positive, else ``n_reads``.

**1. Same-chain 3' consolidation** (``posthoc.end3_tolerance``, default 500 bp)

``end5_window`` is a guard on this rule, not a second tolerance: it blocks a
merge when two same-chain models' 5' ends differ by more than it. Set it
negative to remove the guard entirely, which makes the 3' end the whole of the
rule and leaves 5' truncation to the fragment filter below. That is a real
choice, not a loosening: the fragment filter only ever looks at a model whose
chain is a STRICT sub-chain of another, so it never sees two models that share
a chain exactly -- with the guard off, two identical chains differing only in
where transcription started are merged, alternative TSS included.

Models are partitioned by contig, strand and their exact intron chain.  No
model is ever compared across partitions, so this rule cannot merge two
different splice structures -- it removes terminal-end redundancy only.
Within a partition the representative is chosen as

    annotated 3' end  >  support  >  transcript length  >  start

and every model whose 5' end lies within ``end5_window`` and whose 3' end lies
within ``end3_tolerance`` **of that representative** is absorbed: its reads move
onto the representative and support is recounted from the reads.  The next
unabsorbed model starts the next cluster.  Anchoring to the representative
rather than chaining pairwise is the point: with 3' ends at 0, 400 and 800 and
a 500 bp tolerance, the model at 800 is *not* absorbed, so no cluster can span
more than the tolerance however dense the catalogue is.

The 5' window is a fixed guard, not a knob users are meant to turn: it stops
the rule from merging models that differ by an alternative TSS.

Choosing the *longest* end as representative loses exact-match recall on SIRVs
as the tolerance grows; choosing the annotated, then best-supported end holds
recall at its N = 0 value.  That is why the sort reads the way it does.

**2. Reference-free sub-chain (fragment) filter** (``posthoc.subchain_filter``)

A model is removed when its ordered intron chain is a contiguous sub-chain of a
longer model on the same contig and strand (every junction within
``subchain_fuzzy`` bp), its terminal exons lie inside the matching exons of that
model (allowing ``subchain_terminal_slack`` bp, so an alternative first exon
that starts inside an upstream intron is *not* treated as a fragment), and that
model carries at least ``subchain_ratio`` x its support.  A mono-exonic model is
a fragment when it lies inside one exon of such a model.  This is the
5'-truncation / internal-priming class that dominated FlightCollapse's
annotation-free false positives on SIRV Set 4.

Models whose category is in ``protected_categories`` (default: FSM, an exact
annotated chain) are never removed, and INFERRED models are neither removed nor
used as containers.  Decisions are taken in one pass on the pre-filter supports,
never iteratively; if a container is itself removed, the fragment resolves to
the final surviving container.  Removed models' reads are folded into the
most-supported surviving container when ``fold_reads`` is set, and otherwise
become unassigned and are counted.

Every decision is written to ``<prefix>.posthoc.tsv``.

SCALE
-----
Both rules are indexed rather than quadratic.  The fragment search looks a
candidate's first intron up in a positional index of every intron in the contig
instead of comparing all pairs, and the 3' rule bisects a position-sorted view
of each partition.  A naive all-pairs form of either is minutes per contig on a
whole-genome catalogue, which is how a rule meant to tidy the output turns into
the slowest stage of the run.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .model import TranscriptModel

AUDIT_COLUMNS = (
    "contig", "rule", "removed_strand", "removed_start", "removed_end",
    "removed_n_exons", "removed_category", "removed_n_reads", "removed_n_mols",
    "into_model_index", "into_model_id", "detail",
)

#: positional index granularity.  Any power of two comfortably larger than the
#: junction fuzz works; it only decides how many buckets a lookup touches.
_BUCKET = 1 << 10


@dataclass
class PosthocResult:
    models: List[TranscriptModel]
    audit: List[Dict[str, object]] = field(default_factory=list)
    stats: Dict[str, int] = field(default_factory=dict)


# ---------------------------------------------------------------------- #
def _support(m: TranscriptModel) -> int:
    return m.n_mols if m.n_mols > 0 else m.n_reads


def _five(m: TranscriptModel) -> int:
    return m.end if m.strand == "-" else m.start


def _three(m: TranscriptModel) -> int:
    return m.start if m.strand == "-" else m.end


def _end_annotated(m: TranscriptModel) -> bool:
    """Does the reference place a 3' end here?

    ``end_is_annotated`` is written by the terminal-ends stage when a reference
    was supplied.  Without one the flag is absent everywhere and the category
    fallback is also false, so the whole rule degrades to support order -- which
    is the intended reference-free behaviour, not an accident.
    """
    flag = m.evidence.get("end_is_annotated", None)
    if flag not in (None, ""):
        try:
            return bool(int(flag))
        except (TypeError, ValueError):
            pass
    return m.category in ("FSM", "end3_annotated")


def _is_inferred(m: TranscriptModel, inferred_category: str) -> bool:
    return m.category == inferred_category or "INFERRED" in (m.flags or [])


def _genomic_introns(m: TranscriptModel) -> List[Tuple[int, int]]:
    """Introns in genomic order, which is what both rules compare.

    ``TranscriptModel.chain`` is in transcript orientation and so runs backwards
    on the minus strand; using it here would make a minus-strand sub-chain test
    silently compare reversed sequences.
    """
    ex = m.exons
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]


def _recount(m: TranscriptModel, reads, umi_hamming: int) -> None:
    from .secondpass import _recount as _sp_recount

    _sp_recount(m, reads, umi_hamming)


def _move_reads(dst: TranscriptModel, src: TranscriptModel) -> int:
    """Move src's reads onto dst.  Returns a count that recounting cannot see.

    When a source carries read indices they are unioned in and the caller
    recounts from them.  When it carries only a tally (no index array -- the
    case when read retention is off) there is nothing to union, so its reads are
    added to ``n_reads`` directly and returned, because the recount sets
    ``n_reads`` from the index array and would otherwise discard them.
    """
    if src.reads.size:
        dst.reads = np.unique(np.concatenate([dst.reads, src.reads]))
        return 0
    dst.n_reads += int(src.n_reads)
    return int(src.n_reads)


def _audit_row(contig: str, rule: str, m: TranscriptModel,
               into: Optional[TranscriptModel], detail: str) -> Dict[str, object]:
    return {
        "contig": contig, "rule": rule, "removed_strand": m.strand,
        "removed_start": m.start, "removed_end": m.end, "removed_n_exons": m.n_exons,
        "removed_category": m.category, "removed_n_reads": m.n_reads,
        "removed_n_mols": m.n_mols, "_into": into, "detail": detail,
    }


# ------------------------------------------------------------- rule 1 -- #
def collapse_3p(
    models: List[TranscriptModel], end3: int, end5: int, reads, umi_hamming: int,
    inferred_category: str, contig: str, audit: Optional[List[Dict[str, object]]],
) -> Tuple[List[TranscriptModel], Dict[str, int]]:
    """Same contig + strand + exact intron chain; 5' within end5, 3' within end3."""
    stats = {"end3_groups": 0, "end3_models_absorbed": 0, "end3_reads_moved": 0}
    if end3 < 0:
        return models, stats

    parts: Dict[tuple, List[TranscriptModel]] = defaultdict(list)
    kept: List[TranscriptModel] = []
    for m in models:
        if _is_inferred(m, inferred_category):
            kept.append(m)                      # INFERRED models are untouched
            continue
        parts[(m.contig, m.strand, tuple(_genomic_introns(m)))].append(m)

    unlimited5 = end5 < 0
    for members in parts.values():
        if len(members) == 1:
            kept.append(members[0])
            continue
        three = [_three(m) for m in members]
        five = [_five(m) for m in members]
        # rank order decides who gets to be a representative; the position order
        # is only how candidates are found
        rank = sorted(range(len(members)),
                      key=lambda i: (-int(_end_annotated(members[i])),
                                     -_support(members[i]), -members[i].tx_len,
                                     members[i].start))
        by_three = sorted(range(len(members)), key=lambda i: three[i])
        sorted_three = [three[i] for i in by_three]
        alive = [True] * len(members)
        for ri in rank:
            if not alive[ri]:
                continue
            alive[ri] = False
            rep = members[ri]
            lo = bisect.bisect_left(sorted_three, three[ri] - end3)
            hi = bisect.bisect_right(sorted_three, three[ri] + end3)
            absorbed = 0
            phantom = 0
            for p in range(lo, hi):
                j = by_three[p]
                if not alive[j] or (not unlimited5
                                    and abs(five[j] - five[ri]) > end5):
                    continue
                alive[j] = False
                m = members[j]
                stats["end3_reads_moved"] += int(m.n_reads)
                phantom += _move_reads(rep, m)
                absorbed += 1
                if audit is not None:
                    audit.append(_audit_row(
                        contig, "end3_collapse", m, rep,
                        f"d5={abs(five[j] - five[ri])};d3={abs(three[j] - three[ri])}"))
            if absorbed:
                stats["end3_groups"] += 1
                stats["end3_models_absorbed"] += absorbed
                rep.evidence["posthoc_end3_absorbed"] = absorbed
                if rep.reads.size:
                    _recount(rep, reads, umi_hamming)
                    rep.n_reads += phantom
            kept.append(rep)
    kept.sort(key=lambda m: (m.start, m.end, m.strand))
    return kept, stats


# ------------------------------------------------------------- rule 2 -- #
def _contiguous_at(ai, bi, k: int, fuzzy: int) -> bool:
    return all(abs(ai[t][0] - bi[k + t][0]) <= fuzzy
               and abs(ai[t][1] - bi[k + t][1]) <= fuzzy
               for t in range(len(ai)))


def _terminals_inside(a: TranscriptModel, b: TranscriptModel, k: int, n: int,
                      slack: int) -> bool:
    """A's first and last exons must sit inside B's matching exons.

    This is what separates a 5' truncation from an alternative first exon: a
    model whose first exon begins inside B's upstream intron starts *before*
    B's matching exon, so it fails here and is kept.
    """
    first, last = b.exons[k], b.exons[k + n]
    return a.exons[0][0] >= first[0] - slack and a.exons[-1][1] <= last[1] + slack


def _is_fragment_of(a: TranscriptModel, b: TranscriptModel, fuzzy: int,
                    slack: int) -> bool:
    """Kept as the readable definition; the filter uses the index below."""
    if a is b or a.strand != b.strand or a.contig != b.contig:
        return False
    ai, bi = _genomic_introns(a), _genomic_introns(b)
    if len(ai) >= len(bi):
        return False
    if not ai:
        s, e = a.start, a.end
        return any(x - fuzzy <= s and e <= y + fuzzy for x, y in b.exons)
    n = len(ai)
    for k in range(len(bi) - n + 1):
        if _contiguous_at(ai, bi, k, fuzzy) and _terminals_inside(a, b, k, n, slack):
            return True
    return False


def filter_subchains(
    models: List[TranscriptModel], ratio: float, fuzzy: int, slack: int,
    protected: Sequence[str], fold_reads: bool, reads, umi_hamming: int,
    inferred_category: str, contig: str, audit: Optional[List[Dict[str, object]]],
) -> Tuple[List[TranscriptModel], Dict[str, int]]:
    """Decisions are taken on the pre-filter supports, in one pass."""
    stats = {"subchain_removed": 0, "subchain_reads_folded": 0,
             "subchain_reads_unassigned": 0}
    prot = set(protected)
    cand = [m for m in models if not _is_inferred(m, inferred_category)]
    if len(cand) < 2:
        return models, stats

    introns = {id(m): _genomic_introns(m) for m in cand}
    sup = {id(m): _support(m) for m in cand}

    # positional index: every intron of every possible container, and every
    # exon, so a fragment is matched against the handful of models that share a
    # coordinate instead of against the whole contig
    # Keyed on (contig, strand): the pipeline calls this once per contig, but the
    # function is public and _is_fragment_of rejects a cross-contig pair, so the
    # index has to as well -- two contigs share coordinates, and a bucket that
    # ignored the contig would offer chr2 models as containers for chr1 ones.
    by_donor: Dict[Tuple[str, str, int], List[Tuple[TranscriptModel, int]]] = defaultdict(list)
    exon_rows: Dict[Tuple[str, str], List[Tuple[int, int, TranscriptModel]]] = defaultdict(list)
    for b in cand:
        bi = introns[id(b)]
        for k, (d, _a) in enumerate(bi):
            by_donor[(b.contig, b.strand, d // _BUCKET)].append((b, k))
        if bi:                                   # only spliced models contain
            for (x, y) in b.exons:
                exon_rows[(b.contig, b.strand)].append((x, y, b))
    exon_idx: Dict[Tuple[str, str], Tuple[List[int], List[Tuple[int, int, TranscriptModel]], List[int]]] = {}
    for key, rows in exon_rows.items():
        rows.sort(key=lambda r: r[0])
        starts = [r[0] for r in rows]
        run_max = []
        hi = -1
        for r in rows:                           # running max end, so the scan
            hi = max(hi, r[1])                   # below can stop early
            run_max.append(hi)
        exon_idx[key] = (starts, rows, run_max)

    def _containers_for(a: TranscriptModel) -> List[TranscriptModel]:
        ai = introns[id(a)]
        need = ratio * max(sup[id(a)], 1)
        out: List[TranscriptModel] = []
        if not ai:
            ent = exon_idx.get((a.contig, a.strand))
            if ent is None:
                return out
            starts, rows, run_max = ent
            # every exon starting at or before a.start + fuzzy; walk back while
            # some exon in the prefix still reaches far enough right
            p = bisect.bisect_right(starts, a.start + fuzzy) - 1
            while p >= 0 and run_max[p] >= a.end - fuzzy:
                x, y, b = rows[p]
                if (b is not a and sup[id(b)] >= need
                        and x - fuzzy <= a.start and a.end <= y + fuzzy):
                    out.append(b)
                p -= 1
            return out
        n = len(ai)
        d0, a0 = ai[0]
        seen = set()
        for bucket in range((d0 - fuzzy) // _BUCKET, (d0 + fuzzy) // _BUCKET + 1):
            for b, k in by_donor.get((a.contig, a.strand, bucket), ()):
                if b is a or id(b) in seen:
                    continue
                bi = introns[id(b)]
                if len(bi) <= n or k + n > len(bi):
                    continue
                if abs(bi[k][0] - d0) > fuzzy or abs(bi[k][1] - a0) > fuzzy:
                    continue
                if sup[id(b)] < need:
                    continue
                if _contiguous_at(ai, bi, k, fuzzy) and _terminals_inside(a, b, k, n, slack):
                    seen.add(id(b))
                    out.append(b)
        return out

    container: Dict[int, TranscriptModel] = {}
    for a in cand:
        if a.category in prot:
            continue
        hits = _containers_for(a)
        if not hits:
            continue
        # deterministic: best supported, then leftmost, then longest
        container[id(a)] = min(hits, key=lambda b: (-sup[id(b)], b.start, -b.end))
    if not container:
        return models, stats

    def resolve(m: TranscriptModel) -> TranscriptModel:
        seen = set()
        while id(m) in container and id(m) not in seen:
            seen.add(id(m))
            m = container[id(m)]
        return m

    kept = [m for m in models if id(m) not in container]
    folded: Dict[int, int] = defaultdict(int)
    phantom: Dict[int, int] = defaultdict(int)
    for m in models:
        if id(m) not in container:
            continue
        direct = container[id(m)]
        dst = resolve(m)
        stats["subchain_removed"] += 1
        if fold_reads and id(dst) not in container:
            stats["subchain_reads_folded"] += int(m.n_reads)
            phantom[id(dst)] += _move_reads(dst, m)
            folded[id(dst)] += 1
            into = dst
        else:
            stats["subchain_reads_unassigned"] += int(m.n_reads)
            into = None
        detail = (f"container_support={sup[id(direct)]};own_support={sup[id(m)]};"
                  f"ratio={ratio}")
        if into is not None and into is not direct:
            detail += f";resolved_to={into.start}-{into.end}"
        if audit is not None:
            audit.append(_audit_row(contig, "subchain_fragment", m, into, detail))
    for m in kept:
        n = folded.get(id(m), 0)
        if not n:
            continue
        m.evidence["posthoc_fragments_absorbed"] = int(
            m.evidence.get("posthoc_fragments_absorbed", 0) or 0) + n
        if m.reads.size:
            _recount(m, reads, umi_hamming)
            m.n_reads += phantom.get(id(m), 0)
    return kept, stats


# ---------------------------------------------------------------------- #
def run_posthoc(models: List[TranscriptModel], contig: str, reads, cfg) -> PosthocResult:
    p = cfg.posthoc
    res = PosthocResult(models=models)
    if not p.enabled:
        return res
    audit: Optional[List[Dict[str, object]]] = [] if p.write_audit else None
    inferred = cfg.secondpass.inferred_category
    n0 = len(models)
    models, s1 = collapse_3p(models, int(p.end3_tolerance), int(p.end5_window), reads,
                             cfg.molecules.umi_hamming, inferred, contig, audit)
    res.stats.update(s1)
    if p.subchain_filter:
        models, s2 = filter_subchains(
            models, float(p.subchain_ratio), int(p.subchain_fuzzy),
            int(p.subchain_terminal_slack), list(p.protected_categories),
            bool(p.fold_reads), reads, cfg.molecules.umi_hamming, inferred, contig, audit)
        res.stats.update(s2)
    res.stats["models_before"] = n0
    res.stats["models_after"] = len(models)
    res.models = models
    if audit is not None:
        where = {id(m): j for j, m in enumerate(models)}
        for row in audit:
            into = row.pop("_into")
            row["into_model_index"] = where.get(id(into), -1) if into is not None else -1
            row.setdefault("into_model_id", "")
        res.audit = audit
    return res


def write_audit(rows: Sequence[Dict[str, object]], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(AUDIT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in AUDIT_COLUMNS) + "\n")
