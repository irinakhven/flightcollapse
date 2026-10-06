"""The 0.5.0 second pass: mono-exon gene classes, and the end3_novel gate.

Two kinds of test here.  The unit tests pin the decision functions against
hand-built inputs, so a threshold change shows up as a named failure rather
than as a count drifting somewhere downstream.  The end-to-end tests run the
whole pipeline over ``simulate_tails``, whose loci were built to produce one
specific verdict each -- including the pair that matters most, a genuine
alternative polyA site and an internal-priming decoy that are identical to
every 0.4.6 measurement and differ only in the genome under the tail.
"""

from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import simulate_tails  # noqa: E402

from flightcollapse import Config  # noqa: E402
from flightcollapse import secondpass as sp  # noqa: E402
from flightcollapse.config import SecondPassParams  # noqa: E402
from flightcollapse.model import TranscriptModel  # noqa: E402
from flightcollapse.pipeline import run  # noqa: E402
from flightcollapse.reference import ReferenceIndex  # noqa: E402


# ====================================================================== #
# fixtures
# ====================================================================== #
@pytest.fixture(scope="module")
def tails(tmp_path_factory):
    return simulate_tails.build(str(tmp_path_factory.mktemp("tails")))


def _run(tails, outdir, **overrides):
    cfg = Config(bam=tails.bam, reference_gtf=tails.gtf,
                 genome_fasta=tails.genome_fasta)
    cfg.molecules.barcode_umi_tsv = tails.barcode_umi
    cfg.output.outdir, cfg.output.prefix = str(outdir), "t"
    cfg.verbose, cfg.strict_invariants = False, False
    for k, v in overrides.items():
        cfg.set_value(k, v)
    report = run(cfg)
    tbl = pd.read_csv(os.path.join(outdir, "t.models.tsv"), sep="\t",
                      low_memory=False)
    audit_path = os.path.join(outdir, "t.secondpass.tsv")
    audit = (pd.read_csv(audit_path, sep="\t")
             if os.path.exists(audit_path) else pd.DataFrame())
    return {"report": report, "models": tbl, "audit": audit, "cfg": cfg}


@pytest.fixture(scope="module")
def filtered(tails, tmp_path_factory):
    return _run(tails, tmp_path_factory.mktemp("filtered"))


@pytest.fixture(scope="module")
def unfiltered(tails, tmp_path_factory):
    """The same run with the second pass off -- the control arm.

    Every claim about what the second pass removed is only meaningful against
    a run that differs in nothing else, which is why this is a fixture rather
    than a remembered number.
    """
    return _run(tails, tmp_path_factory.mktemp("unfiltered"),
                **{"secondpass.enabled": False})


# ====================================================================== #
# unit: the A-E classifier
# ====================================================================== #
@pytest.fixture(scope="module")
def ref(tails):
    return ReferenceIndex.from_gtf(tails.gtf, verbose=False)


def _mono(contig, start, end, strand, gene_id):
    return TranscriptModel(
        model_id="", contig=contig, strand=strand, exons=[(start, end)],
        chain=(), category="monoexon_internal", gene_id=gene_id,
        n_reads=10, n_mols=10,
    )


@pytest.mark.parametrize("gene_id,start,end,strand,expect", [
    ("GMONO", 5_100, 6_400, "+", "A"),    # gene has only a single-exon isoform
    ("GBOTH", 14_200, 15_000, "+", "B"),  # gene has both
    ("GMD", 20_300, 21_100, "+", "C"),    # exonic + intronic: spans intron 1
    ("GMD", 22_200, 22_900, "+", "D"),    # wholly inside the last exon
    ("GINT", 37_300, 37_800, "-", "E"),   # wholly inside an intron: not "spans"
])
def test_the_classifier_reads_the_gene_not_the_model(ref, gene_id, start, end,
                                                     strand, expect):
    params = SecondPassParams()
    struct = sp._GeneStructure(ref)
    m = _mono("chr1", start, end, strand, gene_id)
    assert sp.classify_monoexon(m, ref, struct, params) == expect


def test_a_model_with_no_gene_is_out_of_scope(ref):
    """A-E is a statement about a gene's exon structure. Intergenic piles have
    no such structure, and forcing them into a class would be an answer to a
    question nobody asked."""
    params = SecondPassParams()
    struct = sp._GeneStructure(ref)
    assert sp.classify_monoexon(_mono("chr1", 1_000, 2_000, "+", None),
                                ref, struct, params) == "X"


def test_a_spliced_model_is_never_classified(ref):
    params = SecondPassParams()
    struct = sp._GeneStructure(ref)
    m = TranscriptModel(model_id="", contig="chr1", strand="+",
                        exons=[(20_000, 20_400), (21_000, 21_400)],
                        chain=((20_400, 21_000),), category="FSM",
                        gene_id="GMD")
    assert sp.classify_monoexon(m, ref, struct, params) == ""


# ====================================================================== #
# unit: the end3 tier ladder
# ====================================================================== #
def _end3(**evidence):
    m = TranscriptModel(model_id="", contig="chr1", strand="+",
                        exons=[(100, 200)], chain=(), category="end3_novel",
                        gene_id="G", n_reads=10, n_mols=10)
    m.evidence.update(evidence)
    return m


def test_a_catalogued_polya_site_is_the_strongest_tier():
    m = _end3(end_at_polya_site=1)
    assert sp.end3_tier(m, SecondPassParams(), False, False) == "atlas"


def test_a_non_templated_tail_survives():
    """A genuine tail is longer than the genome's A-run at the same place."""
    m = _end3(frac_reads_nontemplated=0.9, median_tail_excess=22,
              genomic_a_run_3p=0)
    assert sp.end3_tier(m, SecondPassParams(), False, False) == "nontemplated"


def test_a_templated_tail_does_not():
    """The decoy: the same tail length, but the genome already accounts for it.

    This is the case 0.4.6 could not see. `tail_molecule_frac` is 1.0 for both
    reads piles; only the comparison with the genome separates them.
    """
    m = _end3(frac_reads_nontemplated=0.0, median_tail_excess=0,
              genomic_a_run_3p=12, polya_motif_found=False)
    assert sp.end3_tier(m, SecondPassParams(), False, False) is None


def test_the_composition_tier_needs_all_three_of_its_parts():
    good = dict(polya_motif_found=True, dse_gu_minus_a=0.3,
                three_prime_dispersion=2.0, frac_reads_nontemplated=0.0,
                median_tail_excess=0)
    p = SecondPassParams()
    assert sp.end3_tier(_end3(**good), p, False, False) == "composition"
    # a hexamer alone is six bases and 3'UTRs supply one by chance
    assert sp.end3_tier(_end3(**{**good, "dse_gu_minus_a": -0.4}), p,
                        False, False) is None
    # a diffuse peak is what RT drop-off looks like, hexamer or not
    assert sp.end3_tier(_end3(**{**good, "three_prime_dispersion": 40.0}), p,
                        False, False) is None
    # and an A-rich downstream is the artefact the tier exists to exclude
    assert sp.end3_tier(_end3(**{**good, "internal_priming": True}), p,
                        False, False) is None


def test_the_backstop_is_tried_last_and_can_be_turned_off():
    bare = _end3(frac_reads_nontemplated=0.0, median_tail_excess=0)
    p = SecondPassParams()
    assert sp.end3_tier(bare, p, True, False) == "gene_dominant"
    assert sp.end3_tier(bare, p, False, True) == "gene_dominant"
    assert sp.end3_tier(bare, p, False, False) is None
    off = SecondPassParams(gene_dominance_backstop=False)
    assert sp.end3_tier(bare, off, True, True) is None


def test_evidence_outranks_dominance():
    """A model that would pass on dominance anyway must still report WHICH
    evidence kept it, or the tier histogram overstates the backstop."""
    m = _end3(end_at_polya_site=1)
    assert sp.end3_tier(m, SecondPassParams(), True, True) == "atlas"


def test_an_absent_measurement_is_not_a_pass():
    """Blank means unmeasured. Reading it as 0, or as satisfied, is the bug
    class this codebase keeps the empty string for."""
    m = _end3(frac_reads_nontemplated="", median_tail_excess="",
              dse_gu_minus_a="", polya_motif_found="")
    assert sp.end3_tier(m, SecondPassParams(), False, False) is None


# ====================================================================== #
# unit: the templated-tail measurement itself
# ====================================================================== #
def test_the_a_run_is_read_in_transcript_orientation(tails):
    """On the minus strand "downstream" runs the other way along the genome.

    Getting this backwards would not raise -- it would quietly measure the
    wrong 30 bases for half the transcriptome.
    """
    from flightcollapse.genome import Genome

    g = Genome(tails.genome_fasta)
    try:
        # the decoy's 12 nt A-run sits at 53,000..53,012 on the plus strand
        assert sp._downstream(g, "chr1", 53_000, "+", 20).startswith("A" * 12)
        # read as a minus-strand 3' end, the same coordinate looks upstream and
        # onto the reverse complement, so the A-run is not there
        assert not sp._downstream(g, "chr1", 53_000, "-", 20).startswith("AA")
        # ...it is there when approached from the far side
        assert sp._downstream(g, "chr1", 53_012, "-", 20).startswith("T" * 12)
    finally:
        g.close()


# ====================================================================== #
# end to end
# ====================================================================== #
def test_every_monoexon_class_was_reached(filtered, tails):
    counts = filtered["report"]["second_pass"]["monoexon_classes"]
    for cls in ("A", "B", "C", "D", "E"):
        assert counts[cls] >= 1, f"the fixture produced no class {cls} model"


def test_classes_a_and_b_survive(filtered):
    audit = filtered["audit"]
    kept = audit[(audit["rule"] == "monoexon_class") & (audit["verdict"] == "kept")]
    assert set(kept["monoexon_gene_class"]) == {"A", "B"}


def test_classes_c_d_e_are_replaced(filtered):
    audit = filtered["audit"]
    rep = audit[(audit["rule"] == "monoexon_class")
                & (audit["verdict"] == "replaced")]
    assert set(rep["monoexon_gene_class"]) == {"C", "D", "E"}
    # and none of them is still in the catalogue
    models = filtered["models"]
    mono = models[models["n_exons"] == 1]
    assert set(mono["monoexon_gene_class"].dropna()) <= {"A", "B", "X"}


def test_a_replacement_prefers_the_gene_s_own_best_supported_model(filtered):
    """MultiData has an emitted spliced model, so the inference comes from the
    data rather than from GENCODE's idea of the canonical isoform."""
    audit = filtered["audit"]
    row = audit[audit["gene_name"] == "MultiData"].iloc[0]
    assert row["reason"] == "infer_from=data"
    assert row["target"] == "ENSTMD1"
    models = filtered["models"]
    tgt = models[models["model_id"] == "ENSTMD1"].iloc[0]
    assert "absorbed_monoexon_class_C" in tgt["flags"]
    # the reads moved with the model: 30 spliced + 20 folded
    assert tgt["n_reads"] == 50


def test_a_gene_with_no_spliced_model_falls_back_to_the_reference(filtered):
    audit = filtered["audit"]
    row = audit[audit["gene_name"] == "MultiRef"].iloc[0]
    assert row["reason"] == "infer_from=reference"
    # MANE_Select, not merely the first transcript of the gene
    assert row["target"] == "ENSTMR1"


def test_an_inferred_model_says_so_three_times(filtered):
    """Category, flag and evidence column. The one failure that must not
    happen is an inferred row being read as an observation."""
    models = filtered["models"]
    inf = models[models["category"] == "INFERRED"]
    assert len(inf) >= 1
    for _i, row in inf.iterrows():
        assert row["category"] == "INFERRED"
        assert "INFERRED" in row["flags"] and "HYPOTHETICAL" in row["flags"]
        assert int(row["inferred"]) == 1
        assert row["inferred_source"] == "reference"
        assert row["model_id"].endswith("|inferred1")
        # it carries the reference structure verbatim, not an invented one
        assert row["n_exons"] > 1


def test_the_genuine_alternative_polya_site_survives(filtered):
    models = filtered["models"]
    apa = models[(models["gene_name"] == "Apa")
                 & (models["category"] == "end3_novel")]
    assert len(apa) == 1, "the genuine APA model was removed"
    row = apa.iloc[0]
    assert row["end3_second_pass_tier"] == "nontemplated"
    assert row["genomic_a_run_3p"] == 0        # nothing templated underneath
    assert row["median_tail_excess"] >= 20     # the molecules supplied the tail
    assert row["frac_reads_nontemplated"] == 1.0


def test_the_internal_priming_decoy_is_removed(filtered):
    models = filtered["models"]
    assert models[(models["gene_name"] == "IntPrime")
                  & (models["category"] == "end3_novel")].empty
    audit = filtered["audit"]
    row = audit[(audit["rule"] == "end3_novel")
                & (audit["gene_name"] == "IntPrime")].iloc[0]
    assert row["verdict"] == "removed"
    assert row["reason"] == "no_surviving_tier"


def test_0_4_6_could_not_have_removed_it(unfiltered):
    """The point of the decoy: without the second pass it is emitted, and it is
    emitted for a *defensible* reason -- the 0.1.18 rule that a molecule's own
    A-run outranks a genomic hexamer. That rule is right; it just cannot tell
    whose A-run it is."""
    models = unfiltered["models"]
    decoy = models[(models["gene_name"] == "IntPrime")
                   & (models["category"] == "end3_novel")]
    assert len(decoy) == 1
    row = decoy.iloc[0]
    # every 0.4.6 signal says this end is fine
    assert row["tail_molecule_frac"] == 1.0
    assert row["three_prime_dispersion"] == 0.0
    assert row["n_reads"] >= 20
    # and the one 0.5.0 signal says it is not
    assert row["genomic_a_run_3p"] == 12
    assert row["median_tail_excess"] == 0


def test_dropped_reads_are_re_homed_not_lost(filtered, unfiltered):
    """Rejecting a 3' peak says which model owns this read mass, not that the
    mass should disappear."""
    rep = filtered["report"]["second_pass"]
    assert rep["end3_reads_dropped"] == rep["end3_reads_folded"]
    assert rep["end3_reads_unassigned"] == 0
    before = unfiltered["models"]
    after = filtered["models"]
    assert (after[after["model_id"] == "ENSTIP1"].iloc[0]["n_reads"]
            == before[before["model_id"] == "ENSTIP1"].iloc[0]["n_reads"]
            + rep["end3_reads_dropped"])


def test_total_read_mass_is_conserved(filtered, unfiltered):
    """Across the whole catalogue, not just at the one locus."""
    assert (filtered["models"]["n_reads"].sum()
            == unfiltered["models"]["n_reads"].sum())


def test_turning_the_second_pass_off_restores_the_first_pass_catalogue(
        filtered, unfiltered):
    """The first pass is unchanged by this release -- assert it rather than
    assume it, because 'only a filter' is exactly the claim that rots."""
    assert unfiltered["report"]["n_models"] > 0
    assert "second_pass" not in unfiltered["report"]
    assert len(unfiltered["models"]) > len(filtered["models"])
    # the evidence columns are still measured with the filter off: they are a
    # measurement, and turning off a decision must not turn off a measurement
    assert unfiltered["models"]["genomic_a_run_3p"].notna().any()
    assert unfiltered["models"]["median_tail_excess"].notna().any()


def test_the_audit_joins_back_to_the_models_table(filtered):
    audit, models = filtered["audit"], filtered["models"]
    kept = audit[audit["verdict"] == "kept"]
    assert len(kept)
    assert set(kept["model_id"]) <= set(models["model_id"])
    # a removed model has no id to join on, and inventing one would be worse
    removed = audit[audit["verdict"].isin(["removed", "replaced"])]
    assert removed["model_id"].isna().all() or (removed["model_id"] == "").all()


def test_the_report_says_how_much_rests_on_the_backstop_alone(filtered):
    """The number to look at first on a real dataset: a gate that keeps most of
    its models because nothing contradicted them is not evidence, and the QC
    should not make that comfortable to ignore."""
    rep = filtered["report"]["second_pass"]
    assert "end3_frac_kept_on_gene_dominance_alone" in rep
    assert rep["end3_kept_by_tier"]["nontemplated"] == 1
    assert rep["end3_frac_kept_on_gene_dominance_alone"] == 0.0


def test_removal_classes_are_configurable(tails, tmp_path_factory):
    """Keeping class D is a defensible choice -- a short last-exon-only model
    in a multi-exon gene can be a real proximal APA event -- so it must be a
    setting, not a decision baked into the code."""
    res = _run(tails, tmp_path_factory.mktemp("keepD"),
               **{"secondpass.remove_classes": ["C", "E"]})
    audit = res["audit"]
    d = audit[audit["monoexon_gene_class"] == "D"]
    assert len(d) and set(d["verdict"]) == {"kept"}


def test_a_gate_that_could_only_ever_delete_is_refused():
    """No live evidence tier and no backstop removes the whole category while
    looking like a filter. That is the one configuration worth refusing."""
    cfg = Config(bam="x.bam", reference_gtf="r.gtf")
    cfg.secondpass.gene_dominance_backstop = False
    cfg.monoexon.enabled = False
    cfg.scoring.enabled = False
    cfg.junctions.require_canonical_motif = False
    cfg.junctions.rt_switch_repeat_len = 0
    cfg.ends.require_polya_evidence_for_utr_variant = False
    with pytest.raises(ValueError, match="would be removed unconditionally"):
        cfg.validate()
