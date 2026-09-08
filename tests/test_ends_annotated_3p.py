"""Gene-level 3'-end snapping.

The handout's methodological trap, made into a test: testing 3'-end coincidence
against the single associated transcript badly undercounts (35% within 10 bp)
versus the nearest annotated end of *any* transcript of the gene (72%). Every
miss would otherwise be reported as a novel APA site.
"""

import numpy as np
import pysam
import pytest

from flightcollapse.config import EndParams, MonoexonParams
from flightcollapse.ends import build_models_for_group
from flightcollapse.reads import ContigReads
from flightcollapse.reference import ReferenceIndex

EXONS = [(1000, 1200), (2000, 2200), (3000, 4000)]
LONG_TTS = 4000
SHORT_TTS = 3400


def _gtf(tmp_path):
    """One gene, two transcripts, identical intron chain, different 3' ends."""
    p = tmp_path / "two_ends.gtf"
    rows = []
    for tid, last_end in (("ENSTLONG", LONG_TTS), ("ENSTSHORT", SHORT_TTS)):
        a = (f'gene_id "G1"; transcript_id "{tid}"; gene_name "Gene1"; '
             f'transcript_type "protein_coding"; tag "basic";')
        ex = EXONS[:-1] + [(EXONS[-1][0], last_end)]
        rows.append(f"chrT\tT\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t+\t.\t{a}")
        for s, e in ex:
            rows.append(f"chrT\tT\texon\t{s + 1}\t{e}\t.\t+\t.\t{a}")
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def _reads(n, three_prime):
    r = ContigReads("chrT")
    r.n = n
    r.strand = np.ones(n, np.int8)
    r.tss = np.full(n, 1000, np.int64)
    r.tts = np.full(n, three_prime, np.int64)
    r.start = np.full(n, 1000, np.int64)
    r.end = np.full(n, three_prime, np.int64)
    r.joff = np.arange(0, 2 * (n + 1), 2, dtype=np.int64)
    r.donor = np.tile(np.array([1200, 2200], np.int64), n)
    r.acceptor = np.tile(np.array([2000, 3000], np.int64), n)
    from flightcollapse.molecules import NO_CELL, NO_UMI
    r.cell = np.full(n, NO_CELL, np.uint32)
    r.umi = np.full(n, NO_UMI, np.uint64)
    return r


CHAIN = ((1200, 2000), (2200, 3000))


def _build(tmp_path, three_prime, **kw):
    ref = ReferenceIndex.from_gtf(_gtf(tmp_path), verbose=False)
    reads = _reads(60, three_prime)
    params = EndParams(**kw)
    return ref, build_models_for_group(
        reads, "chrT", "+", CHAIN, np.arange(60), ref, None, None,
        params, MonoexonParams(), "ENSTLONG", "G1",
    )


def test_end_matching_another_transcript_is_not_called_novel(tmp_path):
    _ref, models = _build(tmp_path, SHORT_TTS)
    assert len(models) == 1
    m = models[0]
    assert m.category == "end3_annotated"
    assert "annotated_alt_3p_end" in m.flags
    assert m.exons[-1][1] == SHORT_TTS      # snapped exactly onto the annotated end


def test_a_genuinely_novel_end_is_still_called_novel(tmp_path):
    _ref, models = _build(tmp_path, 3700)   # between the two annotated ends
    assert models[0].category in ("end3_novel", "FSM")
    assert "annotated_alt_3p_end" not in models[0].flags


def test_the_parents_own_end_is_still_FSM(tmp_path):
    _ref, models = _build(tmp_path, LONG_TTS)
    assert models[0].category == "FSM"
    assert models[0].exons[-1][1] == LONG_TTS


def test_gene_level_snapping_can_be_switched_off(tmp_path):
    _ref, models = _build(tmp_path, SHORT_TTS, snap_to_gene_level_3p_ends=False)
    assert models[0].category != "end3_annotated"
