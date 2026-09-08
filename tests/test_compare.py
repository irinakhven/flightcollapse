"""Cross-sample concordance.

The key property under test: models are matched on structure, never on name.
`GAPDH|NIC_3` is a per-gene counter assigned in emission order, so the same name
in two samples is two different transcripts — matching on it would be wrong in
both directions at once.
"""
import os

import pytest

from flightcollapse.compare import concordance, load_models, reproducibility_curve


def _gtf(path, models):
    """models: (tid, strand, [(s,e)...], category, gene, support)"""
    with open(path, "w") as fh:
        for tid, strand, exons, cat, gene, sup in models:
            a = (f'gene_id "{gene}"; transcript_id "{tid}"; gene_name "{gene}"; '
                 f'category "{cat}"; n_molecules "{sup}";')
            for s, e in exons:
                fh.write(f"chrT\tfc\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{a}\n")
    return str(path)


A = [(100, 200), (300, 400), (500, 700)]
B = [(100, 200), (300, 400), (500, 900)]        # same chain, different 3' end
SKIP = [(100, 200), (500, 700)]                  # different chain


def test_same_structure_matches_under_different_names(tmp_path):
    p1 = _gtf(tmp_path / "s1.gtf", [("GENE|NIC_1", "+", A, "NIC", "GENE", 40)])
    p2 = _gtf(tmp_path / "s2.gtf", [("GENE|NIC_7", "+", A, "NIC", "GENE", 55)])
    C, summary = concordance({"a": p1, "b": p2}, verbose=False)
    assert len(C) == 1
    assert C.iloc[0].n_samples == 2
    assert summary.set_index("category").loc["ALL", "pct_reproduced"] == 100.0


def test_different_names_but_different_structures_do_not_match(tmp_path):
    p1 = _gtf(tmp_path / "s1.gtf", [("GENE|NIC_1", "+", A, "NIC", "GENE", 40)])
    p2 = _gtf(tmp_path / "s2.gtf", [("GENE|NIC_1", "+", SKIP, "NIC", "GENE", 40)])
    C, _s = concordance({"a": p1, "b": p2}, verbose=False)
    assert len(C) == 2
    assert set(C.n_samples) == {1}


def test_a_different_3prime_end_is_a_different_model(tmp_path):
    p1 = _gtf(tmp_path / "s1.gtf", [("t", "+", A, "FSM", "GENE", 40)])
    p2 = _gtf(tmp_path / "s2.gtf", [("t", "+", B, "FSM", "GENE", 40)])
    C, _s = concordance({"a": p1, "b": p2}, verbose=False)
    assert len(C) == 2


def test_the_5prime_end_is_not_part_of_the_key(tmp_path):
    """It is a percentile over a 5'-degraded pile and drifts between libraries."""
    shifted = [(40, 200)] + A[1:]
    p1 = _gtf(tmp_path / "s1.gtf", [("t", "+", A, "FSM", "GENE", 40)])
    p2 = _gtf(tmp_path / "s2.gtf", [("t", "+", shifted, "FSM", "GENE", 40)])
    C, _s = concordance({"a": p1, "b": p2}, verbose=False)
    assert len(C) == 1 and C.iloc[0].n_samples == 2


def test_minus_strand_uses_the_left_coordinate_as_the_3prime_end(tmp_path):
    """On the minus strand the leftmost coordinate is the 3' end, so a shift
    there is a different cleavage site — while the same shift on the plus strand
    is just 5' drift and must be ignored. Same chain in both cases."""
    base = [(1000, 1100), (1300, 1400), (1500, 1700)]
    shifted_left = [(500, 1100), (1300, 1400), (1500, 1700)]
    p1 = _gtf(tmp_path / "m1.gtf", [("t", "-", base, "FSM", "GENE", 40)])
    p2 = _gtf(tmp_path / "m2.gtf", [("t", "-", shifted_left, "FSM", "GENE", 40)])
    C, _s = concordance({"a": p1, "b": p2}, verbose=False)
    assert len(C) == 2, "a shifted minus-strand 3' end must not be merged"

    p3 = _gtf(tmp_path / "p1.gtf", [("t", "+", base, "FSM", "GENE", 40)])
    p4 = _gtf(tmp_path / "p2.gtf", [("t", "+", shifted_left, "FSM", "GENE", 40)])
    C2, _s2 = concordance({"a": p3, "b": p4}, verbose=False)
    assert len(C2) == 1, "the same shift on the plus strand is 5' drift"


def test_one_sample_cannot_fill_two_slots_of_a_structure(tmp_path):
    p1 = _gtf(tmp_path / "s1.gtf", [("t1", "+", A, "FSM", "GENE", 40),
                                    ("t2", "+", A, "FSM", "GENE", 12)])
    p2 = _gtf(tmp_path / "s2.gtf", [("t", "+", A, "FSM", "GENE", 40)])
    C, _s = concordance({"a": p1, "b": p2}, verbose=False)
    assert sorted(C.n_samples) == [1, 2]


def test_reproducibility_curve_is_monotone_in_yield(tmp_path):
    ms1 = [(f"t{i}", "+", [(1000 * i, 1000 * i + 100), (1000 * i + 300, 1000 * i + 400)],
            "NIC", "G", i) for i in range(1, 40)]
    ms2 = [m for m in ms1 if m[5] >= 10]
    C, _s = concordance({"a": _gtf(tmp_path / "a.gtf", ms1),
                         "b": _gtf(tmp_path / "b.gtf", ms2)}, verbose=False)
    curve = reproducibility_curve(C, 2, thresholds=(0, 5, 10, 20))
    assert curve.structures.is_monotonic_decreasing
    assert curve.pct_reproduced.iloc[-1] >= curve.pct_reproduced.iloc[0]


def test_load_models_reads_our_own_gtf(collapsed):
    d = load_models(os.path.join(collapsed["outdir"], "sim.gtf"), verbose=False)
    assert len(d) > 0
    assert d.model_id.is_unique
    assert set(d.strand) <= {"+", "-"}


def _pigeon(tmp_path):
    p = tmp_path / "cls_classification.txt"
    p.write_text(
        "isoform\tchrom\tstructural_category\tsubcategory\tFL.sample\n"
        "PB.1.1\tchr1\tfull-splice_match\treference_match\t120\n"
        "PB.1.2\tchr1\tincomplete-splice_match\t3prime_fragment\t45\n"
    )
    return str(p)


def test_categories_come_from_the_pipelines_own_classification(tmp_path):
    """A collapse GFF has no category and no counts -- classification is a
    separate step -- so the comparison reads them back from that table."""
    from flightcollapse.compare import load_categories

    got = load_categories(_pigeon(tmp_path), verbose=False)
    assert got["PB.1.1"] == ("full-splice_match", 120.0)
    assert got["PB.1.2"] == ("incomplete-splice_match", 45.0)


def test_a_plain_two_column_table_also_works(tmp_path):
    from flightcollapse.compare import load_categories

    p = tmp_path / "simple.tsv"
    p.write_text("id\tcategory\tsupport\nT1\tnovel\t7\n")
    got = load_categories(str(p), verbose=False)
    assert got["T1"] == ("novel", 7.0)


def test_categories_override_the_gtf_and_unmatched_models_are_marked(tmp_path):
    from flightcollapse.compare import load_categories, load_models

    gtf = tmp_path / "iso.gff"
    gtf.write_text(
        'chr1\tPacBio\texon\t100\t200\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1";\n'
        'chr1\tPacBio\texon\t400\t500\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1";\n'
        'chr1\tPacBio\texon\t100\t250\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.9.9";\n'
    )
    d = load_models(str(gtf), verbose=False,
                    categories=load_categories(_pigeon(tmp_path), verbose=False))
    by = d.set_index("model_id")
    assert by.loc["PB.1.1", "category"] == "full-splice_match"
    assert by.loc["PB.1.1", "support"] == 120.0
    # a model absent from the classification is named, not silently blank
    assert by.loc["PB.9.9", "category"] == "unclassified"


def test_transcript_ids_are_never_used_to_match_across_samples(tmp_path):
    """`PB.7.2` is a per-run counter: the same string in two files is two
    unrelated transcripts, so matching must be structural."""
    from flightcollapse.compare import concordance

    a = tmp_path / "a.gff"
    b = tmp_path / "b.gff"
    # same id, completely different structure
    a.write_text(
        'chr1\tPacBio\texon\t100\t200\t.\t+\t.\ttranscript_id "PB.1.1";\n'
        'chr1\tPacBio\texon\t400\t500\t.\t+\t.\ttranscript_id "PB.1.1";\n')
    b.write_text(
        'chr1\tPacBio\texon\t9000\t9100\t.\t+\t.\ttranscript_id "PB.1.1";\n'
        'chr1\tPacBio\texon\t9400\t9500\t.\t+\t.\ttranscript_id "PB.1.1";\n')
    C, _ = concordance({"a": str(a), "b": str(b)}, verbose=False)
    assert set(C.n_samples) == {1}          # shared id, but nothing in common

    # ...and different ids with the same structure DO match
    b.write_text(
        'chr1\tPacBio\texon\t100\t200\t.\t+\t.\ttranscript_id "PB.42.7";\n'
        'chr1\tPacBio\texon\t400\t500\t.\t+\t.\ttranscript_id "PB.42.7";\n')
    C, _ = concordance({"a": str(a), "b": str(b)}, verbose=False)
    assert list(C.n_samples) == [2]


def test_min_support_filters_while_parsing(tmp_path):
    """An unfiltered isoseq GFF holds millions of transcripts; building the exon
    table for all of them only to drop 94% costs gigabytes for nothing."""
    from flightcollapse.compare import load_models

    g = tmp_path / "x.gtf"
    g.write_text(
        'chr1\tX\texon\t100\t200\t.\t+\t.\ttranscript_id "A"; n_reads "50";\n'
        'chr1\tX\texon\t400\t500\t.\t+\t.\ttranscript_id "A"; n_reads "50";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\ttranscript_id "B"; n_reads "3";\n'
        'chr1\tX\texon\t700\t800\t.\t+\t.\ttranscript_id "B"; n_reads "3";\n')
    assert set(load_models(str(g), verbose=False).model_id) == {"A", "B"}
    kept = load_models(str(g), verbose=False, min_support=20)
    assert set(kept.model_id) == {"A"}
    assert kept.iloc[0].n_exons == 2          # both exons survived the filter


def test_support_attr_picks_reads_or_molecules(tmp_path):
    """isoseq reports FL reads, so a like-for-like threshold must compare reads
    with reads, not reads with UMI-collapsed molecules."""
    from flightcollapse.compare import load_models

    g = tmp_path / "y.gtf"
    g.write_text(
        'chr1\tX\texon\t100\t200\t.\t+\t.\ttranscript_id "A"; n_reads "50"; n_molecules "4";\n'
        'chr1\tX\texon\t400\t500\t.\t+\t.\ttranscript_id "A"; n_reads "50"; n_molecules "4";\n')
    assert load_models(str(g), verbose=False).iloc[0].support == 4.0
    assert load_models(str(g), verbose=False,
                       support_attr="n_reads").iloc[0].support == 50.0
    assert len(load_models(str(g), verbose=False, min_support=20)) == 0
    assert len(load_models(str(g), verbose=False, min_support=20,
                           support_attr="n_reads")) == 1


# ---------------------------------------------------------------------- #
# the matcher itself
# ---------------------------------------------------------------------- #
def _gtf_at(tmp_path, name, ends):
    """One two-exon transcript per 3' end given, all sharing an intron chain."""
    p = tmp_path / f"{name}.gtf"
    rows = []
    for i, e in enumerate(ends):
        t = f"{name}.{i}"
        rows.append(f'chr1\tX\texon\t100\t200\t.\t+\t.\ttranscript_id "{t}"; n_reads "50";')
        rows.append(f'chr1\tX\texon\t400\t{e}\t.\t+\t.\ttranscript_id "{t}"; n_reads "50";')
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def test_a_looser_tolerance_never_creates_more_structures(tmp_path):
    """The property that exposed the old matcher.

    Merging two neighbouring sites must reduce the structure count, never raise
    it. The old code cut the block on linkage gaps and then forced any repeated
    sample into a unique singleton group, so widening the window manufactured
    singletons: on the real BD144 0.1.10 output, going from 100 to 300 bp took
    283,662 structures to 308,889 and 'reproducibility' from 47.5% to 40.5%.
    """
    from flightcollapse.compare import concordance

    a = _gtf_at(tmp_path, "a", [1000, 1200])
    b = _gtf_at(tmp_path, "b", [1100, 1300])
    counts = []
    for tol in (20, 150, 300, 1000):
        C, _ = concordance({"a": a, "b": b}, max_3p_diff=tol, verbose=False)
        counts.append(len(C))
    assert counts == sorted(counts, reverse=True), counts
    assert counts[0] == 4      # tol 20: nothing is close enough to pair
    assert counts[-1] == 2     # wider: two pairs, and it cannot go below that


def test_each_sample_appears_at_most_once_per_structure(tmp_path):
    from flightcollapse.compare import concordance

    a = _gtf_at(tmp_path, "a", [1000, 1040, 1080])
    b = _gtf_at(tmp_path, "b", [1005])
    C, _ = concordance({"a": a, "b": b}, max_3p_diff=200, verbose=False)
    assert (C.n_samples <= 2).all()
    assert C.n_samples.max() == 2          # b pairs with exactly one a model
    assert len(C) == 3                     # the other two a models stand alone


def test_fragmentation_in_one_sample_does_not_destroy_the_match(tmp_path):
    """One sample splitting a site into three must still match the other's one.

    This is the shape the 0.1.7 clustering regression produced, and under the
    old matcher it cost the *pair* as well as adding singletons.
    """
    from flightcollapse.compare import concordance

    a = _gtf_at(tmp_path, "a", [1000, 1045, 1090])
    b = _gtf_at(tmp_path, "b", [1090])
    C, _ = concordance({"a": a, "b": b}, max_3p_diff=100, verbose=False)
    assert int((C.n_samples == 2).sum()) == 1
    paired = C[C.n_samples == 2].iloc[0]
    assert abs(paired.three - 1090) <= 100


def test_the_densest_position_wins_not_the_leftmost(tmp_path):
    """a@1000 alone, then a@1200+b@1205: the pair must survive."""
    from flightcollapse.compare import concordance

    a = _gtf_at(tmp_path, "a", [1000, 1200])
    b = _gtf_at(tmp_path, "b", [1205])
    C, _ = concordance({"a": a, "b": b}, max_3p_diff=250, verbose=False)
    assert int((C.n_samples == 2).sum()) == 1
    assert abs(C[C.n_samples == 2].iloc[0].three - 1200) <= 10
