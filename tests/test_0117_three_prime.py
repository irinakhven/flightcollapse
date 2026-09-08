"""0.1.17: the 3' end may be refined by annotation, never manufactured by it.

Three changes are covered here.

* **Bounded annotation fallback.**  A rejected 3' peak used to be re-labelled
  with the default parent's annotated end at any distance.  The NSL1 case moved
  a model 11.6 kb past every read supporting it.  Rejecting a peak means "not
  enough evidence to call a new cleavage site", not "the molecules ended
  somewhere else".
* **Molecule-weighted peak geometry.**  Peak discovery, the dominant peak and
  the peak coordinate were driven by raw reads while the gates counted
  molecules, so a PCR family could decide which site a chain group was named
  after.
* **Direct polyA tail evidence**, read off the 3' soft clip.  Recorded by
  default, enforced only when asked.
"""

import numpy as np
import pytest

from flightcollapse.config import AlignmentGates, EndParams, MonoexonParams
from flightcollapse.ends import (
    build_models_for_group,
    molecule_weights,
    tail_evidence,
    weighted_median,
)
from flightcollapse.molecules import NO_CELL, NO_UMI
from flightcollapse.reads import ContigReads, polya_tail_stats
from flightcollapse.reference import ReferenceIndex

EXONS = [(1000, 1200), (2000, 2200), (3000, 4000)]
CHAIN = ((1200, 2000), (2200, 3000))
LONG_TTS, SHORT_TTS = 4000, 3400


def _gtf(tmp_path):
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


def _reads(tts_per_read, cells=None, umis=None):
    n = len(tts_per_read)
    r = ContigReads("chrT")
    r.n = n
    r.strand = np.ones(n, np.int8)
    r.tss = np.full(n, 1000, np.int64)
    r.tts = np.array(tts_per_read, np.int64)
    r.start = np.full(n, 1000, np.int64)
    r.end = r.tts.copy()
    r.joff = np.arange(0, 2 * (n + 1), 2, dtype=np.int64)
    r.donor = np.tile(np.array([1200, 2200], np.int64), n)
    r.acceptor = np.tile(np.array([2000, 3000], np.int64), n)
    r.cell = (np.full(n, NO_CELL, np.uint32) if cells is None
              else np.array(cells, np.uint32))
    r.umi = (np.full(n, NO_UMI, np.uint64) if umis is None
             else np.array(umis, np.uint64))
    return r


def _build(tmp_path, reads, chain_tx=("ENSTLONG", "ENSTSHORT"), **kw):
    ref = ReferenceIndex.from_gtf(_gtf(tmp_path), verbose=False)
    return build_models_for_group(
        reads, "chrT", "+", CHAIN, np.arange(reads.n), ref, None, None,
        EndParams(**kw), MonoexonParams(), list(chain_tx), "G1", umi_hamming=0,
    )


# ------------------------------------------------- the NSL1 failure mode ----
def test_distant_annotated_end_no_longer_overrides_the_reads(tmp_path):
    """The regression. 8 kb from any annotated end, no polyA evidence."""
    models = _build(tmp_path, _reads([12000] * 40))
    assert len(models) == 1
    m = models[0]
    assert m.exons[-1][1] == 12000, "model was moved off its own read pile"
    assert m.category == "end3_unresolved"
    assert any(f.startswith("3p_fallback_refused_dist=") for f in m.flags)
    assert m.evidence["tts_shift_from_reads"] == 0


def test_a_near_annotated_end_still_wins(tmp_path):
    """Refinement is the correct behaviour and must survive the fix."""
    models = _build(tmp_path, _reads([4200] * 40))
    m = models[0]
    assert m.exons[-1][1] == LONG_TTS       # snapped: 200 bp <= max_3p_fallback_dist
    assert m.category == "FSM"
    assert "3p_variant_rejected_no_polya" in m.flags
    assert not any(f.startswith("3p_fallback_refused") for f in m.flags)


def test_the_boundary_is_the_configured_one(tmp_path):
    near = _build(tmp_path, _reads([4300] * 40), max_3p_fallback_dist=300)[0]
    far = _build(tmp_path, _reads([4301] * 40), max_3p_fallback_dist=300)[0]
    assert near.category == "FSM" and near.exons[-1][1] == LONG_TTS
    assert far.category == "end3_unresolved" and far.exons[-1][1] == 4301


def test_0116_behaviour_is_recoverable(tmp_path):
    """A huge bound reproduces the old collapse, so the change is attributable."""
    m = _build(tmp_path, _reads([12000] * 40), max_3p_fallback_dist=10 ** 9)[0]
    assert m.category == "FSM"
    assert m.exons[-1][1] == LONG_TTS


def test_unresolved_can_be_dropped_instead_of_emitted(tmp_path):
    assert _build(tmp_path, _reads([12000] * 40), emit_unresolved_3p=False) == []


def test_shift_from_reads_is_recorded_when_annotation_moves_the_end(tmp_path):
    m = _build(tmp_path, _reads([4200] * 40))[0]
    assert m.evidence["tts_shift_from_reads"] == LONG_TTS - 4200


# ------------------------------------------------- molecule weighting -------
def _pcr_heavy():
    """20 reads from 2 molecules at the long end; 6 reads from 6 at the short."""
    tts = [LONG_TTS] * 20 + [SHORT_TTS] * 6
    cells = [1] * 26
    umis = [100] * 10 + [101] * 10 + list(range(200, 206))
    return _reads(tts, cells, umis)


def test_a_pcr_family_no_longer_decides_the_dominant_peak(tmp_path):
    models = _build(tmp_path, _pcr_heavy(), weight_ends_by_molecule=True)
    ends = [m.exons[-1][1] for m in models]
    assert set(ends) == {LONG_TTS, SHORT_TTS}
    assert ends[0] == SHORT_TTS, "6 molecules should outrank 2 molecules"
    by_end = {m.exons[-1][1]: m for m in models}
    assert by_end[LONG_TTS].evidence["n_mol_equivalents"] == pytest.approx(2.0)
    assert by_end[SHORT_TTS].evidence["n_mol_equivalents"] == pytest.approx(6.0)
    assert by_end[LONG_TTS].n_reads == 20      # gates still see whole reads


def test_read_weighting_is_the_old_answer(tmp_path):
    models = _build(tmp_path, _pcr_heavy(), weight_ends_by_molecule=False)
    assert models[0].exons[-1][1] == LONG_TTS


def test_molecule_weights_are_exact_family_reciprocals():
    r = _reads([1, 2, 3, 4], cells=[1, 1, 1, 2], umis=[7, 7, 8, 7])
    w = molecule_weights(r, np.arange(4), umi_hamming=0)
    assert w.tolist() == [0.5, 0.5, 1.0, 1.0]
    assert w.sum() == pytest.approx(3.0)       # three distinct molecules


def test_unknown_provenance_is_not_treated_as_duplication():
    r = _reads([1, 2, 3])                      # all NO_CELL / NO_UMI
    assert molecule_weights(r, np.arange(3), 0).tolist() == [1.0, 1.0, 1.0]


def test_weighted_median_follows_the_weights():
    v = np.array([10, 20, 30], np.int64)
    assert weighted_median(v, np.array([1.0, 1.0, 1.0])) == 20
    assert weighted_median(v, np.array([9.0, 1.0, 1.0])) == 10
    assert weighted_median(v, np.array([1.0, 1.0, 9.0])) == 30


# ------------------------------------------------- polyA tail from clips ----
class _Aln:
    def __init__(self, cigar, seq, reverse=False):
        self.cigartuples, self.query_sequence, self.is_reverse = cigar, seq, reverse


def test_tail_is_read_off_the_forward_soft_clip():
    # 100 M then a 30 base clip that opens with 12 A's
    seq = "C" * 100 + "A" * 12 + "GCGCGCGCGCGCGCGCGC"
    clip, run, frac = polya_tail_stats(_Aln([(0, 100), (4, 30)], seq))
    assert (clip, run) == (30, 12)
    assert frac == pytest.approx(12 / 30)


def test_minus_strand_tail_reads_as_poly_t_from_the_cleavage_site_outward():
    # reference orientation: clip first, then the aligned block
    seq = "GCGCGCGCGCGCGCGCGC" + "T" * 12 + "C" * 100
    clip, run, frac = polya_tail_stats(_Aln([(4, 30), (0, 100)], seq, reverse=True))
    assert (clip, run) == (30, 12), "the run adjacent to the alignment is the tail"
    assert frac == pytest.approx(12 / 30)


def test_an_internal_priming_style_clip_scores_low():
    seq = "C" * 100 + "AAGACAGCTTTTGGCAGACT"
    _clip, run, frac = polya_tail_stats(_Aln([(0, 100), (4, 20)], seq))
    assert run == 2 and frac < 0.5


def test_hard_clipped_bases_are_not_counted_as_tail():
    """query_sequence lacks them, so counting them would index wrong bases."""
    seq = "C" * 100
    assert polya_tail_stats(_Aln([(0, 100), (5, 40)], seq)) == (0, 0, 0.0)


def test_tail_evidence_is_zero_when_measurement_was_off():
    r = _reads([100, 100, 100])                # tail arrays left empty
    ev = tail_evidence(r, np.arange(3), np.ones(3), EndParams())
    assert ev == {"tail_molecule_frac": 0.0, "n_tail_reads": 0.0,
                  "median_tail_len": 0.0}


def test_tail_evidence_is_molecule_weighted():
    r = _reads([100] * 4, cells=[1, 1, 1, 1], umis=[1, 1, 2, 3])
    r.tail_len = np.array([20, 20, 20, 0], np.uint16)
    r.tail_frac = np.array([1.0, 1.0, 1.0, 0.0], np.float32)
    w = molecule_weights(r, np.arange(4), 0)          # [.5, .5, 1, 1]
    ev = tail_evidence(r, np.arange(4), w, EndParams())
    assert ev["n_tail_reads"] == 3.0
    assert ev["tail_molecule_frac"] == pytest.approx(2.0 / 3.0)


def test_tail_evidence_can_rescue_a_novel_end_when_asked(tmp_path):
    """On by default since 0.1.18; off, the peak stays unresolved."""
    r = _reads([12000] * 40)
    r.tail_len = np.full(40, 25, np.uint16)
    r.tail_frac = np.full(40, 0.95, np.float32)
    off = _build(tmp_path, r, polya_tail_gates_novel_end=False)[0]
    on = _build(tmp_path, r)[0]
    assert off.category == "end3_unresolved"
    assert on.category == "end3_novel" and "polya_from_read_tail" in on.flags
    assert on.exons[-1][1] == 12000


# ------------------------------------------------- config surface -----------
def test_new_defaults_are_the_documented_ones():
    p = EndParams()
    assert p.max_3p_fallback_dist == 300
    assert p.unresolved_3p_category == "end3_unresolved"
    assert p.emit_unresolved_3p is True
    assert p.weight_ends_by_molecule is True
    assert p.polya_tail_gates_novel_end is True        # 0.1.18: the tail decides
    assert p.min_polya_tail_frac == 0.0                # the biased conjunction is off
    assert AlignmentGates().measure_polya_tail is True


def test_version():
    import flightcollapse
    assert flightcollapse.__version__ >= "0.1.18"
