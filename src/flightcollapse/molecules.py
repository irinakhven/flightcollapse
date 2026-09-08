"""Barcode / UMI handling.

The BAMs carry no ``CB``/``UB`` tags -- ``isoseq tag`` cannot handle
variable-length BD Rhapsody barcodes and ``groupdedup`` was never run -- so the
mapping arrives as an external, read-name-keyed TSV that is tens of gigabytes.

Holding 85 M ``read_name -> (barcode, umi)`` entries in a Python dict would need
~10 GB.  Instead the table is compiled once into three parallel numpy arrays

    read_hash : uint64  (sorted; 64-bit blake2b of the read name)
    cell_id   : uint32  (index into a barcode vocabulary)
    umi_code  : uint64  (2-bit packed UMI, exact and reversible)

which is ~20 bytes per read (1.7 GB for 85 M reads) and supports vectorised
lookup by ``np.searchsorted``.  The index is cached next to the TSV so the
19 GB scan happens once.

**Why this matters for the collapse.**  Read counts are the wrong support
statistic for a PCR-amplified library: a novel junction seen in 40 reads from
2 molecules is an artefact, the same junction in 40 reads from 35 molecules in
20 cells is real.  Every threshold in this package is expressed on molecules
when a barcode/UMI table is available, and only falls back to reads when it is
not.
"""

from __future__ import annotations

import hashlib
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

_B2I = {"A": 0, "C": 1, "G": 2, "T": 3}
NO_UMI = np.uint64(0xFFFFFFFFFFFFFFFF)
NO_CELL = np.uint32(0xFFFFFFFF)
#: length field value marking a hashed (non-ACGT / over-long) UMI.  Any value
#: > 31 works; it makes the code a valid identity key that is deliberately not
#: decodable and not Hamming-correctable.
HASHED = 63


def hash_read_name(name: str) -> int:
    return int.from_bytes(
        hashlib.blake2b(name.encode("ascii", "ignore"), digest_size=8).digest(), "little"
    )


# ---------------------------------------------------------------------- #
# read-name normalisation
# ---------------------------------------------------------------------- #
_SEGMENT = re.compile(r"/\d+_\d+$")


def normalise_name(name: str, mode: str = "none") -> str:
    """Reconcile read names between the BAM and the barcode/UMI table.

    Kinnex / segmented-read pipelines rename reads at several points, so the
    two sides do not always agree:

    ``none``           use the name as-is
    ``strip_segment``  drop a trailing ``/<start>_<end>`` segment span, e.g.
                       ``movie/12345/ccs/2645_3379`` -> ``movie/12345/ccs``
    ``zmw``            keep only ``movie/zmw``, dropping everything after the
                       second field -- the coarsest match, and the right one if
                       the table predates read segmentation
    ``strip_ccs``      drop a trailing ``/ccs`` (with or without a segment span)

    Applied identically when the index is built and when it is queried, so the
    two sides can never drift apart.
    """
    if mode == "none":
        return name
    if mode == "strip_segment":
        return _SEGMENT.sub("", name)
    if mode == "strip_ccs":
        n = _SEGMENT.sub("", name)
        return n[:-4] if n.endswith("/ccs") else n
    if mode == "zmw":
        parts = name.split("/")
        return "/".join(parts[:2]) if len(parts) >= 2 else name
    raise ValueError(f"unknown read-name normalisation mode: {mode!r}")


def encode_umi(umi: str) -> int:
    """2-bit pack an ACGT UMI of length <= 31 into a uint64 (length in the top bits).

    Non-ACGT characters or over-long UMIs fall back to a 58-bit hash, which is
    still a valid identity key but cannot be Hamming-corrected.
    """
    n = len(umi)
    if n == 0 or n > 31:
        return (hash_read_name(umi) >> 6) | (HASHED << 58)
    v = 0
    for ch in umi:
        b = _B2I.get(ch)
        if b is None:
            return (hash_read_name(umi) >> 6) | (HASHED << 58)
        v = (v << 2) | b
    return v | (n << 58)


def decode_umi(code: int) -> Optional[str]:
    n = (code >> 58) & 0x3F
    if n == 0 or n > 31:
        return None
    bases = "ACGT"
    out = []
    v = code & ((1 << 58) - 1)
    for i in range(n):
        out.append(bases[(v >> (2 * (n - 1 - i))) & 3])
    return "".join(out)


def umi_length(code: int) -> int:
    return (code >> 58) & 0x3F


def hamming1_neighbours(code: int) -> Iterable[int]:
    """All codes at Hamming distance 1 from a 2-bit packed UMI."""
    n = umi_length(code)
    if n == 0 or n > 31:
        return ()
    body = code & ((1 << 58) - 1)
    head = n << 58
    out = []
    for i in range(n):
        shift = 2 * i
        cur = (body >> shift) & 3
        for b in range(4):
            if b == cur:
                continue
            out.append(((body & ~(3 << shift)) | (b << shift)) | head)
    return out


def default_index_path(tsv: str, normalise: str = "none") -> str:
    """The one place the cached-index filename is defined.

    ``flightcollapse index-umi`` and ``run`` must agree on this, or a
    pre-built index is silently ignored and rebuilt.
    """
    return f"{tsv}.fcidx.{normalise}.npz"


class MoleculeIndex:
    """Compact read-name -> (cell, UMI) index."""

    def __init__(
        self,
        read_hash: np.ndarray,
        cell_id: np.ndarray,
        umi_code: np.ndarray,
        barcodes: List[str],
        normalise: str = "none",
    ) -> None:
        self.read_hash = read_hash
        self.cell_id = cell_id
        self.umi_code = umi_code
        self.barcodes = barcodes
        self.normalise = normalise

    # ------------------------------------------------------------------ #
    @classmethod
    def build(
        cls,
        tsv: str,
        col_cell: str = "cell_barcode",
        col_umi: str = "umi",
        col_read: str = "read_name",
        normalise: str = "none",
        verbose: bool = True,
        progress_every: int = 10_000_000,
    ) -> "MoleculeIndex":
        import gzip

        opener = gzip.open if tsv.endswith(".gz") else open
        hashes: List[int] = []
        cells: List[int] = []
        umis: List[int] = []
        vocab: Dict[str, int] = {}
        with opener(tsv, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                ic = header.index(col_cell)
                iu = header.index(col_umi)
                ir = header.index(col_read)
            except ValueError as exc:
                raise SystemExit(
                    f"{tsv}: expected columns {col_cell}/{col_umi}/{col_read}, got {header}"
                ) from exc
            wide = max(ic, iu, ir)
            for n, line in enumerate(fh, 1):
                f = line.rstrip("\n").split("\t")
                if len(f) <= wide:
                    continue
                cb = f[ic]
                cid = vocab.get(cb)
                if cid is None:
                    cid = vocab[cb] = len(vocab)
                hashes.append(hash_read_name(normalise_name(f[ir], normalise)))
                cells.append(cid)
                umis.append(encode_umi(f[iu]))
                if verbose and n % progress_every == 0:
                    print(f"[umi] {n:,} rows, {len(vocab):,} barcodes")

        h = np.array(hashes, dtype=np.uint64)
        del hashes
        c = np.array(cells, dtype=np.uint32)
        del cells
        u = np.array(umis, dtype=np.uint64)
        del umis
        order = np.argsort(h, kind="stable")
        h, c, u = h[order], c[order], u[order]
        barcodes = [""] * len(vocab)
        for k, v in vocab.items():
            barcodes[v] = k
        if verbose:
            print(f"[umi] indexed {h.size:,} reads / {len(barcodes):,} cell barcodes")
        return cls(h, c, u, barcodes, normalise)

    # ------------------------------------------------------------------ #
    def save(self, path: str) -> None:
        np.savez(
            path,
            read_hash=self.read_hash,
            cell_id=self.cell_id,
            umi_code=self.umi_code,
            barcodes=np.array(self.barcodes, dtype=object),
            normalise=np.array([self.normalise], dtype=object),
        )

    @classmethod
    def load(cls, path: str) -> "MoleculeIndex":
        z = np.load(path, allow_pickle=True)
        mode = str(z["normalise"][0]) if "normalise" in z else "none"
        return cls(z["read_hash"], z["cell_id"], z["umi_code"], list(z["barcodes"]), mode)

    @classmethod
    def open_or_build(
        cls,
        tsv: str,
        index_path: Optional[str] = None,
        verbose: bool = True,
        **kw,
    ) -> "MoleculeIndex":
        mode = kw.get("normalise", "none")
        idx = index_path or default_index_path(tsv, mode)
        if os.path.exists(idx) and os.path.getmtime(idx) >= os.path.getmtime(tsv):
            if verbose:
                print(f"[umi] loading cached index {idx}")
            return cls.load(idx)
        obj = cls.build(tsv, verbose=verbose, **kw)
        obj.save(idx)
        if verbose:
            print(f"[umi] wrote {idx}")
        return obj

    # ------------------------------------------------------------------ #
    def lookup(self, names: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
        """Vectorised ``read_name -> (cell_id, umi_code)``; misses are sentinels."""
        if len(names) == 0:
            return (np.zeros(0, np.uint32), np.zeros(0, np.uint64))
        q = np.fromiter(
            (hash_read_name(normalise_name(n, self.normalise)) for n in names),
            dtype=np.uint64, count=len(names),
        )
        i = np.searchsorted(self.read_hash, q)
        i_clip = np.clip(i, 0, self.read_hash.size - 1)
        hit = self.read_hash[i_clip] == q
        cid = np.where(hit, self.cell_id[i_clip], NO_CELL)
        umi = np.where(hit, self.umi_code[i_clip], NO_UMI)
        return cid.astype(np.uint32), umi.astype(np.uint64)

    def __len__(self) -> int:
        return int(self.read_hash.size)


# ---------------------------------------------------------------------- #
# UMI collapsing
# ---------------------------------------------------------------------- #
def collapse_umis(
    pairs: Sequence[Tuple[int, int]], hamming: int = 1
) -> Dict[Tuple[int, int], int]:
    """Count molecules per (cell, UMI), with directional Hamming-1 correction.

    ``pairs`` is an iterable of ``(cell_id, umi_code)``.  Returns
    ``{(cell_id, representative_umi): n_reads}``.

    Correction follows the UMI-tools *directional* rule: a UMI is absorbed by a
    neighbour at Hamming distance 1 when the neighbour's count is at least
    ``2n - 1``.  This is conservative -- it will not merge two genuinely
    distinct molecules that happen to be sequencing-error neighbours at similar
    abundance.
    """
    per_cell: Dict[int, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for cid, umi in pairs:
        if cid == int(NO_CELL) or umi == int(NO_UMI):
            continue
        per_cell[cid][umi] += 1

    out: Dict[Tuple[int, int], int] = {}
    for cid, counts in per_cell.items():
        if hamming <= 0 or len(counts) == 1:
            for u, n in counts.items():
                out[(cid, u)] = n
            continue
        order = sorted(counts, key=lambda u: (-counts[u], u))
        parent: Dict[int, int] = {}
        for u in order:
            if u in parent:
                continue
            parent[u] = u
            nu = counts[u]
            for v in hamming1_neighbours(u):
                nv = counts.get(v)
                if nv is not None and v not in parent and nu >= 2 * nv - 1:
                    parent[v] = u
        merged: Dict[int, int] = defaultdict(int)
        for u, n in counts.items():
            merged[parent.get(u, u)] += n
        for u, n in merged.items():
            out[(cid, u)] = n
    return out


def n_molecules(pairs: Sequence[Tuple[int, int]], hamming: int = 1) -> int:
    return len(collapse_umis(pairs, hamming))


def n_cells(pairs: Sequence[Tuple[int, int]]) -> int:
    return len({c for c, _ in pairs if c != int(NO_CELL)})
