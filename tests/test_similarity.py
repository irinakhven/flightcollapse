"""The 0.1.7 rules: nearest reference transcript, retained introns, readthrough.

Each test here is one of the edge cases that motivated the release, written so
that a regression names itself.
"""

import numpy as np
import pytest

from flightcollapse.config import EndParams, MonoexonParams
from flightcollapse.ends import annotate_similarity, build_models_for_group
from flightcollapse.intervals import exons_to_chain
from flightcollapse.molecules import NO_CELL, NO_UMI
from flightcollapse.reads import ContigReads
from flightcollapse.reference import ReferenceIndex
from flightcollapse.similarity import (
    compare_to_transcript,
    nearest_transcript,
    readthrough_gene,
    retained_intron_rank,
)


# ---------------------------------------------------------------------- #
# a two-gene locus:  G1 with three isoforms, G2 immediately downstream
# ---------------------------------------------------------------------- #
def _rows(gene, tid, exons, strand="+", extra="", cds=None):
    a = (f'gene_id "{gene}"; transcript_id "{tid}"; gene_name "{gene}"; '
         f'transcript_type "protein_coding"; tag "basic";{extra}')
    out = [f"chrT\tT\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{a}"]
    for s, e in exons:
        out.append(f"chrT\tT\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{a}")
    for s, e in cds or ():
        out.append(f"chrT\tT\tCDS\t{s + 1}\t{e}\t.\t{strand}\t.\t{a}")
    return out


#: LONG   1000-1200 | 2000-2200 | 3000-4000      (MANE, long 3'UTR)
#: SHORT  1000-1200 | 2000-2200 | 3000-3400      same chain, earlier cleavage
#: FIVE    500-700  | 1000-1200 | 2000-2200 | 3000-4000   extra 5' exon
LONG = [(1000, 1200), (2000, 2200), (3000, 4000)]
SHORT = [(1000, 1200), (2000, 2200), (3000, 3400)]
FIVE = [(500, 700), (1000, 1200), (2000, 2200), (3000, 4000)]
G2 = [(5000, 5200), (6000, 6500)]


@pytest.fixture()
def ref(tmp_path):
    p = tmp_path / "locus.gtf"
    rows = []
    rows += _rows("G1", "ENSTLONG", LONG, extra=' tag "MANE_Select";',
                  cds=[(1100, 1200), (2000, 2200), (3000, 3300)])
    rows += _rows("G1", "ENSTSHORT", SHORT,
                  cds=[(1100, 1200), (2000, 2200), (3000, 3300)])
    rows += _rows("G1", "ENSTFIVE", FIVE,
                  cds=[(1100, 1200), (2000, 2200), (3000, 3300)])
    rows += _rows("G2", "ENSTG2", G2)
    p.write_text("\n".join(rows) + "\n")
    return ReferenceIndex.from_gtf(str(p), verbose=False)


def _reads(n, chain, tss, tts, strand=1):
    r = ContigReads("chrT")
    r.n = n
    k = len(chain)
    r.strand = np.full(n, strand, np.int8)
    r.tss = np.full(n, tss, np.int64)
    r.tts = np.full(n, tts, np.int64)
    r.start = np.full(n, min(tss, tts), np.int64)
    r.end = np.full(n, max(tss, tts), np.int64)
    r.joff = np.arange(0, k * (n + 1), k, dtype=np.int64)
    r.donor = np.tile(np.array([d for d, _ in chain], np.int64), n)
    r.acceptor = np.tile(np.array([a for _, a in chain], np.int64), n)
    r.goff = np.zeros(n + 1, np.int64)
    r.cell = np.full(n, NO_CELL, np.uint32)
    r.umi = np.full(n, NO_UMI, np.uint64)
    return r


CHAIN3 = ((1200, 2000), (2200, 3000))
CHAIN4 = ((700, 1000), (1200, 2000), (2200, 3000))


# ---------------------------------------------------------------------- #
# (3) the end comparison must use every isoform, not the MANE one
# ---------------------------------------------------------------------- #
def test_ends_are_matched_against_every_isoform_not_mane(ref):
    """A pile ending exactly at ENSTSHORT's 3' end is an FSM to ENSTSHORT.

    Comparing it to the MANE transcript instead would report an alternative
    3' end for a transcript that is annotated in full.
    """
    reads = _reads(60, CHAIN3, 1000, 3400)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        EndParams(), MonoexonParams(), ["ENSTLONG", "ENSTSHORT"], "G1",
    )
    assert len(models) == 1
    assert models[0].category == "FSM"
    assert models[0].parent_tx == "ENSTSHORT"


def test_best_transcript_for_ends_prefers_the_closer_end_over_the_rank(ref):
    got = ref.best_transcript_for_ends(
        ["ENSTLONG", "ENSTSHORT"], five=1000, three=3400, strand="+"
    )
    assert got == "ENSTSHORT"
    # ...and MANE still wins when the ends give no reason to prefer either
    assert ref.best_transcript(["ENSTSHORT", "ENSTLONG"]) == "ENSTLONG"


# ---------------------------------------------------------------------- #
# (6) a 5' extension that crosses a junction is retention, not a new TSS
# ---------------------------------------------------------------------- #
def test_five_prime_extension_across_a_junction_is_intron_retention(ref):
    """The pile starts at 500 -- inside ENSTFIVE's first exon -- but has no
    junction at 700-1000, so it read through ENSTFIVE's first intron."""
    reads = _reads(60, CHAIN3, 500, 4000)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        EndParams(), MonoexonParams(), ["ENSTLONG", "ENSTSHORT"], "G1",
    )
    m = models[0]
    assert m.category == "IR"
    assert "5p_extension_crosses_annotated_junction" in m.flags


def test_five_prime_extension_with_no_junction_crossed_is_a_start_difference(ref):
    """Same direction, but the extension stops short of ENSTFIVE's junction."""
    reads = _reads(60, CHAIN3, 850, 4000)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        EndParams(), MonoexonParams(), ["ENSTLONG"], "G1",
    )
    assert models[0].category in ("end5_extended", "end5_alt")
    assert "5p_extension_crosses_annotated_junction" not in models[0].flags


def test_five_prime_truncation_is_absorbed_and_counted(ref):
    reads = _reads(60, CHAIN3, 3000, 4000)     # starts in the last exon
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        EndParams(), MonoexonParams(), ["ENSTLONG"], "G1",
    )
    m = models[0]
    assert m.category == "FSM"
    assert "5p_truncation_absorbed" in m.flags
    assert m.exons[0][0] == 1000               # folded back onto the annotated start
    assert m.evidence["n_5p_truncated_reads"] == 60


# ---------------------------------------------------------------------- #
# retained introns
# ---------------------------------------------------------------------- #
def test_retention_needs_exonic_flank_on_both_sides(ref):
    tx = ref.tx["ENSTFIVE"]
    covered = [(400, 1300), (2000, 2200), (3000, 4000)]
    assert compare_to_transcript(CHAIN3, covered, "+", tx).n_retained == 1
    # the same missing junction, but the model starts inside the intron: that is
    # a 5'-truncated molecule, not an unspliced one
    truncated = [(800, 1300), (2000, 2200), (3000, 4000)]
    assert compare_to_transcript(CHAIN3, truncated, "+", tx).n_retained == 0


def test_retained_intron_rank_marks_the_five_prime_most_intron(ref):
    tx = ref.tx["ENSTFIVE"]
    m = compare_to_transcript(CHAIN3, [(400, 1300), (2000, 2200), (3000, 4000)], "+", tx)
    assert retained_intron_rank(tx.chain, m.retained) == 0


def test_pure_retention_is_reclassified_out_of_NIC(ref):
    """A chain that is a reference chain minus one intron it read through is a
    retention event, not a novel combination of junctions."""
    from flightcollapse.model import TranscriptModel

    m = TranscriptModel(
        model_id="", contig="chrT", strand="+",
        exons=[(1000, 1200), (2000, 4000)],       # second intron unspliced
        chain=((1200, 2000),), category="NIC", gene_id="G1",
    )
    annotate_similarity([m], ref)
    assert m.category == "IR"
    assert m.associated_tx == "ENSTLONG"
    assert "intron_retention" in m.flags
    assert m.evidence["ir_introns"] == "2200-3000"


def test_retention_reported_against_the_isoform_that_explains_it(ref):
    """Chain-equal to ENSTLONG, but the first exon crosses ENSTFIVE's junction.

    Scored plainly the exact chain match wins; the useful statement is the
    retention, so the IR path looks for a transcript that explains it.
    """
    reads = _reads(60, CHAIN3, 400, 4000)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        EndParams(), MonoexonParams(), ["ENSTLONG"], "G1",
    )
    annotate_similarity(models, ref)
    assert models[0].category == "IR"
    assert models[0].associated_tx == "ENSTFIVE"


# ---------------------------------------------------------------------- #
# nearest transcript, ties, and the UTR question
# ---------------------------------------------------------------------- #
def test_nearest_transcript_reports_ties(ref):
    """LONG and SHORT share a chain and differ only at the 3' end; a model
    landing midway between them is genuinely ambiguous and says so."""
    m = nearest_transcript(
        ref, "G1", CHAIN3, [(1000, 1200), (2000, 2200), (3000, 3700)], "+"
    )
    assert m.tx_id in ("ENSTLONG", "ENSTSHORT")
    assert m.n_shared == 2 and m.n_model_only == 0


def test_is_utr_only_is_false_when_the_model_starts_inside_the_cds(ref):
    # CDS starts at 1100; a model beginning at 1150 has a different protein,
    # not a shorter 5'UTR
    inside = nearest_transcript(
        ref, "G1", CHAIN3, [(1150, 1200), (2000, 2200), (3000, 4000)], "+"
    )
    assert inside.is_utr_only is False
    outside = nearest_transcript(
        ref, "G1", CHAIN3, [(1000, 1200), (2000, 2200), (3000, 3800)], "+"
    )
    assert outside.is_utr_only is True


def test_structural_diff_is_human_readable(ref):
    m = nearest_transcript(
        ref, "G1", CHAIN3, [(1000, 1200), (2000, 2200), (3000, 3800)], "+"
    )
    assert "d3=" in m.describe()


# ---------------------------------------------------------------------- #
# (1) fragment superchains, (5) readthrough, (7) gene conflicts
# ---------------------------------------------------------------------- #
def test_annotated_superchains_are_counted(ref):
    assert ref.annotated_superchains("G1", CHAIN3) == ["ENSTFIVE"]
    assert ref.annotated_superchains("G1", CHAIN4) == []


def test_readthrough_into_the_next_gene_is_detected(ref):
    assert readthrough_gene(ref, "chrT", "+", "G1", 5100) == "G2"
    assert readthrough_gene(ref, "chrT", "+", "G1", 3500) is None
    assert readthrough_gene(ref, "chrT", "+", "G1", 4600) is None   # intergenic


def test_readthrough_becomes_its_own_category(ref):
    reads = _reads(60, CHAIN3, 1000, 5100)
    params = EndParams(require_polya_evidence_for_utr_variant=False)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        params, MonoexonParams(), ["ENSTLONG"], "G1",
    )
    assert models[0].category == "READTHROUGH"
    assert any(f.startswith("readthrough_into=") for f in models[0].flags)


# ---------------------------------------------------------------------- #
# the 0.1.7 regression: clustering window vs naming tolerance
# ---------------------------------------------------------------------- #
def test_the_clustering_window_is_separable_and_defaults_to_the_tolerance(ref):
    """Clustering width and naming width are different questions.

    They are separable knobs; the default ties them, because that is what the
    data supports. 0.1.13 briefly defaulted `peak_window` to 100 on the theory
    that the narrow window fragmented cleavage sites -- a conclusion drawn from
    a comparison run through a broken matcher. Re-measured with the matcher
    fixed (0.1.14), the narrow window is better at every support threshold from
    5 molecules up: 86.9% vs 78.3% reproducible at >=10, and 94.4% vs 91.5% at
    >=20 on *fewer* structures. The sub-peaks it resolves are real -- structures
    called in all three libraries agree on the 3' end to a median of 0 bp,
    including unsnapped `end3_novel` models.
    """
    n = 90
    reads = _reads(n, CHAIN3, 1000, 4000)
    reads.tts = np.array(([3900] * 30) + ([3945] * 30) + ([3990] * 30), np.int64)
    reads.end = reads.tts.copy()
    P = dict(require_polya_evidence_for_utr_variant=False, min_end_share=0.0)

    # default: window == max_3p_diff == 30, so the three sub-peaks are resolved
    resolved = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(n), ref, None, None,
        EndParams(**P), MonoexonParams(), ["ENSTLONG"], "G1",
    )
    assert len(resolved) == 3

    # a wider window treats the pile as one site
    merged = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(n), ref, None, None,
        EndParams(peak_window=100, **P), MonoexonParams(), ["ENSTLONG"], "G1",
    )
    assert len(merged) == 1


def test_the_naming_tolerance_still_decides_FSM(ref):
    """Widening the clustering window must not widen what counts as FSM."""
    reads = _reads(60, CHAIN3, 1000, 3950)      # 50 bp short of ENSTLONG's end
    params = EndParams(require_polya_evidence_for_utr_variant=False,
                       min_end_share=0.0)
    models = build_models_for_group(
        reads, "chrT", "+", CHAIN3, np.arange(60), ref, None, None,
        params, MonoexonParams(), ["ENSTLONG"], "G1",
    )
    assert len(models) == 1
    assert models[0].category != "FSM"          # 50 > max_3p_diff of 30
