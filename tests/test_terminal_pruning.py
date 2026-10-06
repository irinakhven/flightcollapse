"""Competitive terminal pruning and the terminal score (0.3.0, T1/T2/T3).

The contract:

* a shorter 3' variant explained as truncation of a stronger same-chain sibling
  is absorbed, and its reads go to the sibling rather than anywhere else;
* a variant with evidence the sibling lacks -- an annotated end, a catalogued
  polyA site, a tight peak, a short-read coverage step -- survives regardless
  of how much better supported the sibling is;
* read mass is conserved by every one of those decisions.

The third is the one worth a test of its own. Pruning that loses reads is a
worse failure than the ladder it replaces, and it would be invisible in a
model count.
"""

import numpy as np
import pytest

from flightcollapse.config import TerminalParams
from flightcollapse.ends import prune_terminal_variants
from flightcollapse.model import TranscriptModel


def _model(three, n_reads, *, chain=((1200, 2000),), strand="+", disp=40.0,
           annotated=0, atlas=0, step=None, first=1000):
    exons = [(first, 1200), (2000, three)] if strand == "+" else \
            [(three, 1200), (2000, first)]
    m = TranscriptModel(
        model_id=f"m{three}", contig="chrT", strand=strand, exons=exons,
        chain=chain, category="end3_novel", reads=np.arange(n_reads),
        n_reads=n_reads, n_mols=n_reads, n_cells=max(1, n_reads // 5),
    )
    m.evidence.update({
        "three_prime_dispersion": disp,
        "end_is_annotated": annotated,
        "end_at_polya_site": atlas,
    })
    if step is not None:
        m.evidence["sr_3p_step_supported"] = step
    return m


P = TerminalParams()


def test_a_weak_shorter_variant_is_absorbed():
    distal = _model(4000, 500)
    short = _model(3400, 20)
    out = prune_terminal_variants([distal, short], "+", P)
    assert len(out) == 1
    assert out[0] is distal
    assert out[0].n_reads == 520
    assert out[0].evidence["n_3p_truncated_reads"] == 20
    assert "absorbed_3p_variant" in out[0].flags


def test_read_mass_is_conserved():
    """Pruning moves reads; it must never lose them."""
    models = [_model(4000, 500), _model(3400, 20), _model(3000, 11)]
    before = sum(m.n_reads for m in models)
    out = prune_terminal_variants(models, "+", P)
    assert sum(m.n_reads for m in out) == before


def test_a_comparable_sibling_does_not_absorb():
    """The sibling must clearly dominate, not merely lead."""
    out = prune_terminal_variants([_model(4000, 100), _model(3400, 80)], "+", P)
    assert len(out) == 2


def test_an_annotated_end_survives_anything():
    out = prune_terminal_variants(
        [_model(4000, 5000), _model(3400, 6, annotated=1)], "+", P)
    assert len(out) == 2


def test_a_catalogued_polya_site_survives_anything():
    out = prune_terminal_variants(
        [_model(4000, 5000), _model(3400, 6, atlas=1)], "+", P)
    assert len(out) == 2


def test_a_tight_peak_survives():
    """The discriminator the benchmark found and nothing was using.

    BD176c medians: 4.1 bp for FSM against 34.6 for end3_novel. Cleavage is
    precise, so a tight peak is a real site even when a sibling is deeper.
    """
    out = prune_terminal_variants(
        [_model(4000, 5000), _model(3400, 6, disp=3.0)], "+", P)
    assert len(out) == 2
    # ...and the same model with a diffuse peak is absorbed
    out = prune_terminal_variants(
        [_model(4000, 5000), _model(3400, 6, disp=60.0)], "+", P)
    assert len(out) == 1


def test_a_short_read_coverage_step_survives():
    out = prune_terminal_variants(
        [_model(4000, 5000), _model(3400, 6, step=1)], "+", P)
    assert len(out) == 2


def test_distance_bounds_the_claim():
    """Beyond max_truncation_dist, 'truncation of' stops being the simpler story."""
    far = TerminalParams(max_truncation_dist=100)
    out = prune_terminal_variants([_model(4000, 500), _model(3400, 20)], "+", far)
    assert len(out) == 2


def test_minus_strand_orientation():
    """The silent failure mode: absorbing the wrong way round on one strand."""
    distal = _model(500, 500, strand="-", first=4000)
    short = _model(1100, 20, strand="-", first=4000)
    out = prune_terminal_variants([distal, short], "-", P)
    assert len(out) == 1
    assert out[0] is distal
    assert out[0].n_reads == 520


def test_a_ladder_collapses_onto_the_single_strongest():
    models = [_model(4000, 900), _model(3600, 30), _model(3400, 20), _model(3200, 10)]
    out = prune_terminal_variants(models, "+", P)
    assert len(out) == 1
    assert out[0].n_reads == 960
    assert out[0].evidence["n_3p_truncated_reads"] == 60


def test_pruning_can_be_switched_off():
    off = TerminalParams(prune_terminal_variants=False)
    assert off.prune_terminal_variants is False
    # the function itself is unconditional; the pipeline consults the flag,
    # so assert the flag exists and defaults on rather than duplicating logic
    assert P.prune_terminal_variants is True


def test_a_single_model_is_untouched():
    one = [_model(4000, 500)]
    assert prune_terminal_variants(one, "+", P) == one
