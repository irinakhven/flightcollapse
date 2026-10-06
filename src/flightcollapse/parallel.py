"""Contig-level parallelism.  New in 0.5.0.

``run()`` has always streamed one contig at a time, because a whole-genome
collapse OOM-killed the VM.  That constraint turns out to be the thing that
makes parallelism safe: a contig's result depends on no other contig, and -- the
part worth stating explicitly -- all three calibrated models are fitted **per
contig** already, so distributing contigs across processes changes no
calibration, no threshold and no decision.

What is genuinely shared is small and all of it lives in the parent:

* ``IdAssigner``    model ids must be unique and stable across the whole run.
  Genes do not cross contigs, so per-gene counters would be safe in a worker,
  but the collision set and the ``PB.N`` gene numbering are global, and
  ``gene_name`` repeats across contigs (``Y_RNA`` exists three times).  Ids are
  therefore assigned in the parent, in the configured contig order, which makes
  them identical to a serial run's.
* ``read_stat`` / ``group``  are one file each, in contig order.  Workers write
  per-contig shards keyed on a contig-local model index; the parent stitches
  them in order once ids exist.  One extra pass over a text file, against hours
  of collapse.
* the molecule matrix, which is accumulated from per-contig count shards.

Determinism is the acceptance test, not a nice-to-have: a parallel run must be
byte-identical to a serial one, and ``tests/test_parallel.py`` asserts it.

Memory, not cores, is the binding constraint.  Each worker holds one contig's
read arrays plus its own cached chromosome sequence (chr1 is ~250 MB as a
Python str), so peak RSS is roughly ``workers x per-contig peak``.
``Config.worker_count()`` clamps for this; a run that swaps is slower than a
serial one.
"""

from __future__ import annotations

import os
import shutil
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------- #
@dataclass
class ContigResult:
    """One contig's entire contribution, in a form that can cross a pipe.

    ``models`` carry no read arrays by the time this is returned (they are
    released after the shard is written, exactly as in the serial path), so the
    pickle cost is the evidence dicts and nothing else.
    """

    contig: str
    ok: bool = True
    error: str = ""
    models: List = field(default_factory=list)
    #: ``(anchor_model_index, chain_repr, n_reads, n_molecules, is_model_chain)``
    group_chain_rows: List[Tuple[int, str, int, int, bool]] = field(default_factory=list)
    substitution_rows: List = field(default_factory=list)
    junction_frame: object = None
    model_reports: List[Dict[str, object]] = field(default_factory=list)
    mono_reject_rows: List[Dict[str, object]] = field(default_factory=list)
    secondpass_audit: List[Dict[str, object]] = field(default_factory=list)
    posthoc_audit: List[Dict[str, object]] = field(default_factory=list)
    accounting: Dict[str, int] = field(default_factory=dict)
    merge_stats: Dict[str, int] = field(default_factory=dict)
    summary: Dict[str, object] = field(default_factory=dict)
    validator: object = None
    #: shard stems written by this contig, or "" where nothing was written
    shard_stem: str = ""
    n_written: int = 0
    n_unassigned: int = 0
    seconds: float = 0.0


# ---------------------------------------------------------------------- #
class ContigShard:
    """Per-contig ``read_stat`` / ``group`` / molecule-count shards.

    The model is referred to by its **contig-local index**, because the worker
    cannot know the model's id -- that is assigned in the parent once every
    contig has reported.  The parent substitutes ids while stitching.
    """

    def __init__(self, stem: str) -> None:
        self.stem = stem
        os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
        self._rs = open(stem + ".read_stat", "w")
        self._grp = open(stem + ".group", "w")
        self.n_written = 0
        self.n_unassigned = 0
        self._counts: List[Tuple[int, int, int]] = []

    def write_read(self, name: str, length: int, local_idx: int) -> None:
        self._rs.write(f"{name}\t{length}\t{local_idx}\n")
        self.n_written += 1

    def write_group(self, local_idx: int, names: Sequence[str]) -> None:
        self._grp.write(f"{local_idx}\t{','.join(names)}\n")

    def add_counts(self, local_idx: int, per_cell: Dict[int, int]) -> None:
        for cid, n in per_cell.items():
            self._counts.append((local_idx, int(cid), int(n)))

    def close(self) -> None:
        self._rs.close()
        self._grp.close()
        if self._counts:
            arr = np.array(self._counts, dtype=np.int64)
            np.savez_compressed(
                self.stem + ".counts.npz",
                midx=arr[:, 0], cell=arr[:, 1], count=arr[:, 2],
            )
        self._counts = []


def _open_out(path: str, gz: bool):
    import gzip

    return gzip.open(path + ".gz", "wt") if gz else open(path, "w")


def stitch_shards(
    stems_and_ids: Sequence[Tuple[str, Sequence[str]]],
    read_stat_path: str,
    group_path: str,
    gz: bool = False,
) -> Dict[str, int]:
    """Concatenate per-contig shards into the final read_stat / group files.

    ``stems_and_ids`` is in the run's contig order, and each entry pairs a shard
    stem with that contig's model ids indexed by the local index the shard
    wrote.  Row order inside a contig is the shard's, which is the BAM's, so the
    stitched file is the same file a serial run produces.
    """
    n_written = 0
    with _open_out(read_stat_path, gz) as rs, _open_out(group_path, gz) as grp:
        rs.write("id\tlength\tis_fl\tstat\tpbid\n")
        for stem, ids in stems_and_ids:
            if not stem:
                continue
            rs_path = stem + ".read_stat"
            if os.path.exists(rs_path):
                with open(rs_path) as fh:
                    for line in fh:
                        name, length, idx = line.rstrip("\n").split("\t")
                        rs.write(f"{name}\t{length}\tY\tunique\t{ids[int(idx)]}\n")
                        n_written += 1
            g_path = stem + ".group"
            if os.path.exists(g_path):
                with open(g_path) as fh:
                    for line in fh:
                        idx, names = line.rstrip("\n").split("\t", 1)
                        grp.write(f"{ids[int(idx)]}\t{names}\n")
    return {"read_stat_rows": n_written}


def load_counts(stem: str, ids: Sequence[str]) -> Dict[str, Dict[int, int]]:
    """``model_id -> {cell_id: n_molecules}`` from one contig's count shard."""
    path = stem + ".counts.npz"
    out: Dict[str, Dict[int, int]] = {}
    if not stem or not os.path.exists(path):
        return out
    with np.load(path) as z:
        midx, cell, count = z["midx"], z["cell"], z["count"]
    for i in range(midx.size):
        mid = ids[int(midx[i])]
        out.setdefault(mid, {})[int(cell[i])] = int(count[i])
    return out


def cleanup(shard_dir: Optional[str], keep: bool) -> None:
    if shard_dir and not keep and os.path.isdir(shard_dir):
        shutil.rmtree(shard_dir, ignore_errors=True)


# ---------------------------------------------------------------------- #
def order_contigs(
    contigs: Sequence[str], lengths: Dict[str, int], longest_first: bool
) -> List[str]:
    """Dispatch order.  The return order of results is not this order.

    Wall-clock is floored by the slowest single contig, so a schedule that
    starts chr1 last wastes the entire tail: with 24 contigs and 12 workers,
    first-come order finishes roughly a whole chr1 later than longest-first.
    """
    if not longest_first:
        return list(contigs)
    return sorted(contigs, key=lambda c: (-int(lengths.get(c, 0)), c))


def start_method_available(name: str) -> bool:
    import multiprocessing as mp

    try:
        return name in mp.get_all_start_methods()
    except Exception:  # pragma: no cover - defensive
        return False


#: set in the parent before the pool is forked; workers inherit both.
_WORKER = None
_ANNOUNCE = None


def _announced(contig: str) -> "ContigResult":
    """Tell the parent which pid owns this contig, then do the work.

    The announcement is what makes a lost worker detectable.  Without it the
    parent knows only that a result has not arrived, which is indistinguishable
    from a contig that is still running.
    """
    try:
        _ANNOUNCE.put((contig, os.getpid()))
    except Exception:  # pragma: no cover - the work matters more than the note
        pass
    return _WORKER(contig)


def _pid_alive(pid: int) -> bool:
    if os.path.isdir("/proc"):
        return os.path.isdir(f"/proc/{pid}")
    try:  # pragma: no cover - non-Linux
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def map_contigs(
    worker,
    contigs: Sequence[str],
    workers: int,
    start_method: str = "fork",
    log=None,
    poll: float = 2.0,
    grace: float = 30.0,
) -> Dict[str, ContigResult]:
    """Run ``worker(contig)`` over ``contigs``, returning results by contig.

    ``worker`` must be a module-level function; everything it needs beyond the
    contig name is read from a module global set **before** the pool is created,
    so ``fork`` shares the parsed reference index without copying or re-parsing
    it.  That is why only ``fork`` is supported: under ``spawn`` every worker
    would re-read the GTF, which costs more than the parallelism saves.

    WHY THIS IS NOT ``imap_unordered``
    ----------------------------------
    ``multiprocessing.Pool`` has no abrupt-death detection.  When a worker is
    SIGKILLed -- and the OOM killer is the way that happens here -- the pool
    quietly forks a replacement, but the task that worker was holding is gone:
    no result is ever put on the result queue, and ``imap_unordered`` waits for
    it forever.  Observed in the wild on BD144_kin, which sat at 23 of 24
    contigs for two and a half hours with every worker idle and nothing in the
    log, because a single tiny contig's worker had been killed.  A run that
    hangs silently after doing 95% of the work is the worst failure mode this
    tool has, so the parent now tracks which pid owns which contig, notices
    when one disappears, and re-runs that contig itself.

    The retry is deliberately done **in the parent, serially, after the pool has
    closed**: if the cause was memory then re-dispatching into a full pool would
    just kill it again, whereas one contig alone in a quiet process is the
    cheapest configuration available.  ``grace`` covers the race between a child
    exiting normally and its result arriving down the pipe.
    """
    import multiprocessing as mp

    if workers <= 1 or len(contigs) <= 1:
        return {c: worker(c) for c in contigs}
    if not start_method_available(start_method):  # pragma: no cover - platform
        if log:
            log(f"[run] start method {start_method!r} is not available on this "
                f"platform; falling back to a serial run")
        return {c: worker(c) for c in contigs}

    def _say(msg: str) -> None:
        if log:
            log(msg)

    def _report(res: "ContigResult") -> None:
        _say(f"[{res.contig}] done in {res.seconds:.1f}s "
             f"({len(res.models):,} models)"
             if res.ok else f"[{res.contig}] FAILED: {res.error}")

    global _WORKER, _ANNOUNCE
    ctx = mp.get_context(start_method)
    _WORKER = worker
    _ANNOUNCE = ctx.SimpleQueue()

    out: Dict[str, ContigResult] = {}
    lost: List[str] = []
    owner: Dict[str, int] = {}
    vanished_at: Dict[str, float] = {}

    try:
        # maxtasksperchild=1 bounds a worker's RSS to a single contig.
        # Re-forking is cheap under copy-on-write and a worker that has already
        # held chr1 has a page table nobody wants to carry into chr22.
        with ctx.Pool(processes=workers, maxtasksperchild=1) as pool:
            pending = {c: pool.apply_async(_announced, (c,))
                       for c in list(contigs)}
            while pending:
                while not _ANNOUNCE.empty():
                    c, pid = _ANNOUNCE.get()
                    owner[c] = int(pid)

                for c in [c for c, r in pending.items() if r.ready()]:
                    handle = pending.pop(c)
                    vanished_at.pop(c, None)
                    try:
                        res = handle.get()
                    except Exception as exc:  # pragma: no cover - defensive
                        res = ContigResult(contig=c, ok=False, error=repr(exc))
                    out[res.contig] = res
                    _report(res)

                now = time.monotonic()
                for c in list(pending):
                    pid = owner.get(c)
                    # No pid yet means the task is still queued, not lost.
                    if pid is None or _pid_alive(pid):
                        vanished_at.pop(c, None)
                        continue
                    since = vanished_at.setdefault(c, now)
                    if now - since < grace:
                        continue
                    pending.pop(c)
                    lost.append(c)
                    _say(f"[{c}] worker pid {pid} disappeared without "
                         f"returning a result (killed, most likely by the OOM "
                         f"killer); will re-run this contig in the parent")

                if pending:
                    time.sleep(poll)
    finally:
        _WORKER = None
        _ANNOUNCE = None

    for c in lost:
        _say(f"[{c}] re-running serially in the parent process")
        try:
            res = worker(c)
        except Exception as exc:  # pragma: no cover - defensive
            res = ContigResult(contig=c, ok=False, error=repr(exc))
        out[c] = res
        _report(res)
    if lost:
        _say(f"[run] {len(lost)} contig(s) lost a worker and were recomputed "
             f"serially: {', '.join(lost)}. If this recurs, lower -j or run "
             f"fewer samples at once -- the usual cause is memory.")
    return out
