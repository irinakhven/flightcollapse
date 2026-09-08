"""Exon-internal deletions: the difference an identical intron chain hides.

An exact chain match fixes every internal exon *boundary*, because those are
defined by the flanking junctions.  It says nothing about the exon *contents*:
a short intron spelled as a CIGAR ``D`` instead of an ``N`` leaves the chain
untouched, so a real splicing event can sit inside a model labelled FSM.
"""

import numpy as np
import pytest

from flightcollapse.config import AlignmentGates, JunctionParams
from flightcollapse.junctions import promote_internal_gaps, reads_with_unexplained_gap
from flightcollapse.molecules import NO_CELL, NO_UMI
from flightcollapse.reads import ContigReads, internal_gaps
from flightcollapse.reference import ReferenceIndex

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP = 0, 1, 2, 3


class _FakeRead:
    def __init__(self, start, cigar):
        self.reference_start = start
        self.cigartuples = cigar


class _FakeGenome:
    """Reports a canonical motif for one nominated interval and nothing else."""

    def __init__(self, canonical=()):
        self.canonical = set(canonical)

    def splice_motif(self, contig, donor, acceptor, strand):
        return "GT-AG" if (donor, acceptor) in self.canonical else "AC-TT"

    def motif_class(self, motif):
        return "canonical" if motif == "GT-AG" else "non_canonical"


def _gtf(tmp_path):
    """One transcript with a genuine 40 bp intron at 1240-1280."""
    a = ('gene_id "G1"; transcript_id "ENST1"; gene_name "G1"; '
         'transcript_type "protein_coding"; tag "basic";')
    ex = [(1000, 1240), (1280, 1500), (2000, 2400)]
    rows = [f"chrT\tT\ttranscript\t1001\t2400\t.\t+\t.\t{a}"]
    rows += [f"chrT\tT\texon\t{s + 1}\t{e}\t.\t+\t.\t{a}" for s, e in ex]
    p = tmp_path / "g.gtf"
    p.write_text("\n".join(rows) + "\n")
    return ReferenceIndex.from_gtf(str(p), verbose=False)


def _reads_with_gap(n, gap, chain=((1500, 2000),)):
    r = ContigReads("chrT")
    r.n = n
    r.strand = np.ones(n, np.int8)
    r.tss = np.full(n, 1000, np.int64)
    r.tts = np.full(n, 2400, np.int64)
    r.start = np.full(n, 1000, np.int64)
    r.end = np.full(n, 2400, np.int64)
    k = len(chain)
    r.joff = np.arange(0, k * (n + 1), k, dtype=np.int64)
    r.donor = np.tile(np.array([d for d, _ in chain], np.int64), n)
    r.acceptor = np.tile(np.array([a for _, a in chain], np.int64), n)
    r.goff = np.arange(n + 1, dtype=np.int64)
    r.gap_start = np.full(n, gap[0], np.int64)
    r.gap_end = np.full(n, gap[1], np.int64)
    r.cell = np.full(n, NO_CELL, np.uint32)
    r.umi = np.full(n, NO_UMI, np.uint64)
    return r


# ---------------------------------------------------------------------- #
def test_internal_gaps_reads_long_deletions_from_the_cigar():
    r = _FakeRead(1000, [(BAM_CMATCH, 240), (BAM_CDEL, 40), (BAM_CMATCH, 220),
                         (BAM_CREF_SKIP, 500), (BAM_CMATCH, 400)])
    assert internal_gaps(r, 20) == [(1240, 1280)]
    # short indels are sequencing noise, not structure
    assert internal_gaps(r, 60) == []


def test_a_deletion_over_an_annotated_junction_is_promoted(tmp_path):
    ref = _gtf(tmp_path)
    reads = _reads_with_gap(10, (1240, 1280))
    stats = promote_internal_gaps(reads, "chrT", ref, None, JunctionParams())
    assert stats["gap_promoted_annotated"] == 1
    assert stats["gap_records_promoted"] == 10
    # the gap is now an ordinary junction and will be curated like one
    assert reads.chain(0) == ((1240, 1280), (1500, 2000))


def test_an_unannotated_gap_needs_recurrence_and_a_canonical_motif(tmp_path):
    ref = _gtf(tmp_path)
    genome = _FakeGenome(canonical=[(1600, 1650)])

    lonely = _reads_with_gap(2, (1600, 1650))     # canonical, but below recurrence
    assert promote_internal_gaps(
        lonely, "chrT", ref, genome, JunctionParams()
    )["gap_promoted_motif"] == 0

    recurrent = _reads_with_gap(10, (1600, 1650))
    assert promote_internal_gaps(
        recurrent, "chrT", ref, genome, JunctionParams()
    )["gap_promoted_motif"] == 1

    noncanonical = _reads_with_gap(10, (1700, 1750))
    assert promote_internal_gaps(
        noncanonical, "chrT", ref, genome, JunctionParams()
    )["gap_promoted_motif"] == 0


def test_a_gap_that_stays_a_gap_still_flags_its_reads(tmp_path):
    """Alignment noise must not be promoted -- curation would then reject the
    invented junction and throw the whole read away -- but it is not silent."""
    ref = _gtf(tmp_path)
    reads = _reads_with_gap(10, (1700, 1750))
    promote_internal_gaps(reads, "chrT", ref, _FakeGenome(), JunctionParams())
    assert reads.chain(0) == ((1500, 2000),)          # untouched
    assert reads_with_unexplained_gap(reads).all()


def test_promotion_can_be_switched_off(tmp_path):
    ref = _gtf(tmp_path)
    reads = _reads_with_gap(10, (1240, 1280))
    params = JunctionParams(promote_internal_deletions=False)
    assert promote_internal_gaps(
        reads, "chrT", ref, None, params
    )["gap_records_promoted"] == 0
    assert reads.chain(0) == ((1500, 2000),)


def test_oversized_deletions_are_never_promoted(tmp_path):
    ref = _gtf(tmp_path)
    reads = _reads_with_gap(10, (1240, 1280))
    params = JunctionParams(max_internal_deletion_len=10)
    assert promote_internal_gaps(
        reads, "chrT", ref, None, params
    )["gap_records_promoted"] == 0
