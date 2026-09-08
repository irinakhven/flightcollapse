import pytest

from flightcollapse.intervals import (
    build_exons,
    exonic_bases_downstream,
    exonic_length,
    exons_to_chain,
    is_suffix,
    total_overlap,
)


def test_chain_is_transcript_oriented():
    exons = [(100, 200), (300, 400), (500, 600)]
    plus = exons_to_chain(exons, "+")
    minus = exons_to_chain(exons, "-")
    assert plus == ((200, 300), (400, 500))
    assert minus == ((400, 500), (200, 300))
    # the 3'-most intron is always last, on both strands
    assert plus[-1] == (400, 500)
    assert minus[-1] == (200, 300)


def test_five_prime_truncation_is_a_suffix_on_both_strands():
    exons = [(100, 200), (300, 400), (500, 600), (700, 800)]
    for strand, trunc in (("+", exons[1:]), ("-", exons[:-1])):
        full = exons_to_chain(exons, strand)
        short = exons_to_chain(trunc, strand)
        assert is_suffix(short, full), strand


def test_empty_chain_is_never_a_suffix():
    """The vacuous-suffix case is the isoseq mono-exon bug. It must not hold."""
    full = ((200, 300), (400, 500))
    assert not is_suffix((), full)
    assert not is_suffix(full, full)


def test_build_exons_roundtrip():
    exons = [(100, 200), (300, 400), (500, 600)]
    for strand in "+-":
        ch = exons_to_chain(exons, strand)
        tss = exons[0][0] if strand == "+" else exons[-1][1]
        tts = exons[-1][1] if strand == "+" else exons[0][0]
        assert build_exons(ch, strand, tss, tts) == exons


def test_exonic_bases_downstream():
    exons = [(100, 200), (300, 400)]
    assert exonic_bases_downstream(exons, "+", 150) == 50 + 100
    assert exonic_bases_downstream(exons, "-", 350) == 50 + 100
    assert exonic_bases_downstream(exons, "+", 0) == exonic_length(exons)
    assert exonic_bases_downstream(exons, "+", 1000) == 0


def test_total_overlap():
    assert total_overlap([(0, 10), (20, 30)], [(5, 25)]) == 10
    assert total_overlap([(0, 10)], [(10, 20)]) == 0
