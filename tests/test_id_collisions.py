"""Model-id uniqueness, from the two collision classes seen on real data.

The whole-genome BD144 runs emitted 333 duplicate transcript ids out of 151,268
models. Nothing was corrupted -- the `.dupN` suffix is applied before anything
is written -- but a duplicate id is a real hazard (read_stat, group.txt and the
count matrix all key on it), and both causes were naming bugs rather than
modelling ones.
"""

import numpy as np
import pytest

from flightcollapse.model import TranscriptModel
from flightcollapse.output import IdAssigner


def _m(contig, start, end, category, **kw):
    return TranscriptModel(
        model_id="", contig=contig, strand="+",
        exons=[(start, start + 100), (end - 100, end)],
        chain=((start + 100, end - 100),), category=category, **kw
    )


def test_several_annotated_3p_ends_of_one_transcript_get_distinct_names():
    """RBBP4: three real transcripts at three annotated 3' ends of one gene.

    All three used to be named `|3UTR_alt1`, because the suffix index came from
    a counter that only advanced for *novel* 3' ends.
    """
    ms = [
        _m("chr1", 32651208, 32679908, "end3_annotated",
           parent_tx="ENST00000373493.10", gene_id="G", gene_name="RBBP4"),
        _m("chr1", 32651208, 32680053, "end3_annotated",
           parent_tx="ENST00000373493.10", gene_id="G", gene_name="RBBP4"),
        _m("chr1", 32651208, 32680204, "end3_annotated",
           parent_tx="ENST00000373493.10", gene_id="G", gene_name="RBBP4"),
    ]
    a = IdAssigner()
    a.assign(ms)
    assert [x.model_id for x in ms] == [
        "ENST00000373493.10|alt3end1",
        "ENST00000373493.10|alt3end2",
        "ENST00000373493.10|alt3end3",
    ]
    assert a.collisions == 0


def test_a_novel_3p_end_does_not_consume_the_annotated_index():
    """The two suffixes count independently -- that coupling was the bug."""
    ms = [
        _m("chr1", 1000, 5000, "end3_novel", parent_tx="ENST1", gene_id="G"),
        _m("chr1", 1000, 6000, "end3_annotated", parent_tx="ENST1", gene_id="G"),
        _m("chr1", 1000, 7000, "end3_annotated", parent_tx="ENST1", gene_id="G"),
    ]
    a = IdAssigner()
    a.assign(ms)
    assert [x.model_id for x in ms] == [
        "ENST1|novel3end1", "ENST1|alt3end1", "ENST1|alt3end2",
    ]
    assert a.collisions == 0


def test_two_gene_ids_sharing_a_gene_name_do_not_collide():
    """DNAJC9-AS1: the name came from gene_name, the counter from gene_id, so
    each id started its own count at 1 and produced the same string."""
    ms = [
        _m("chr10", 73246883, 73249358, "NIC", gene_id="ENSG_A",
           gene_name="DNAJC9-AS1"),
        _m("chr10", 73252747, 73254349, "NIC", gene_id="ENSG_B",
           gene_name="DNAJC9-AS1"),
    ]
    a = IdAssigner()
    a.assign(ms)
    assert ms[0].model_id != ms[1].model_id
    assert a.collisions == 0
    # ...while the PB ids still count per gene_id, which is what they mean
    assert ms[0].pb_id.endswith(".1") and ms[1].pb_id.endswith(".1")


def test_novel_models_are_named_against_their_closest_transcript():
    ms = [
        _m("chr1", 1000, 5000, "NIC", associated_tx="ENST9", gene_id="G",
           gene_name="GAPDH"),
        _m("chr1", 1000, 6000, "NIC", associated_tx="ENST9", gene_id="G",
           gene_name="GAPDH"),
    ]
    IdAssigner().assign(ms)
    assert [x.model_id for x in ms] == ["ENST9|NIC_1", "ENST9|NIC_2"]


def test_the_collision_net_still_works_if_a_name_is_forced():
    ms = [_m("chr1", 1000, 5000, "NIC", gene_id="G", gene_name="X"),
          _m("chr1", 2000, 6000, "NIC", gene_id="G", gene_name="X")]
    ms[0].model_id = "forced"
    ms[1].model_id = "forced"
    a = IdAssigner()
    a.assign(ms)
    assert a.collisions == 1
    assert ms[1].model_id == "forced.dup2"
    assert "duplicate_model_id" in ms[1].flags


def test_ids_stay_unique_across_contig_by_contig_assignment():
    a = IdAssigner()
    seen = set()
    for contig in ("chr1", "chr2", "chr3"):
        ms = [_m(contig, 1000 + 10 * i, 5000 + 10 * i, "end3_annotated",
                 parent_tx="ENST_SHARED", gene_id="G") for i in range(5)]
        a.assign(ms)
        for m in ms:
            assert m.model_id not in seen
            seen.add(m.model_id)
    assert a.collisions == 0


def test_repeated_gene_names_across_contigs_do_not_collide():
    """`Y_RNA` is three different gene_ids on chr9, chr10 and chr16.

    The mono-exon track used to name its own models with a counter that was
    per-contig and keyed on the gene *name*, so all three came out `Y_RNA|APA1`.
    Naming now happens only in IdAssigner, whose counter is global and keyed on
    the exact string the name is built from.
    """
    ms = [_m(c, 100, 5000, "monoexon_internal", gene_id=f"ENSG_{c}",
             gene_name="Y_RNA") for c in ("chr9", "chr10", "chr16")]
    for m in ms:
        m.exons = [(100, 5000)]
        m.chain = ()
    a = IdAssigner()
    a.assign(ms)
    assert len({m.model_id for m in ms}) == 3
    assert a.collisions == 0
    assert all(m.model_id.startswith("Y_RNA|APA") for m in ms)


def test_intergenic_monoexon_models_are_named_per_contig():
    ms = [_m("chr1", 100 + 10 * i, 5000, "monoexon_intergenic") for i in range(3)]
    for m in ms:
        m.exons = [(m.exons[0][0], 5000)]
        m.chain = ()
    IdAssigner().assign(ms)
    assert [m.model_id for m in ms] == [
        "NOVELMONO_chr1_1", "NOVELMONO_chr1_2", "NOVELMONO_chr1_3"]
