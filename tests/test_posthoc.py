"""0.6.0 final consolidation: same-chain 3' collapse + sub-chain fragment filter."""

from __future__ import annotations

import numpy as np
import pytest

from flightcollapse import Config
from flightcollapse.model import TranscriptModel
from flightcollapse.posthoc import collapse_3p, filter_subchains, run_posthoc

_NEXT = [0]


def mk(exons, strand="+", n=10, category="NIC", flags=None, end_annotated=None):
    """A model whose reads are n fresh, disjoint read indices."""
    exons = [tuple(e) for e in exons]
    introns = tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))
    chain = introns[::-1] if strand == "-" else introns
    reads = np.arange(_NEXT[0], _NEXT[0] + n, dtype=np.int64)
    _NEXT[0] += n
    m = TranscriptModel(model_id="", contig="chr1", strand=strand, exons=list(exons),
                        chain=chain, category=category, reads=reads, n_reads=n,
                        flags=list(flags or []))
    if end_annotated is not None:
        m.evidence["end_is_annotated"] = int(end_annotated)
    return m


def cfg(**ph):
    c = Config()
    for k, v in ph.items():
        setattr(c.posthoc, k, v)
    return c


FOUR = [(1000, 1200), (2000, 2100), (3000, 3100), (4000, 4100), (5000, 5600)]


# --------------------------------------------------------------------- 3' rule
def test_same_chain_3p_variants_merge_within_tolerance():
    a = mk(FOUR, n=50)
    b = mk(FOUR[:-1] + [(5000, 5750)], n=20)          # 3' end 150 bp further
    out, st = collapse_3p([a, b], 200, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1 and st["end3_models_absorbed"] == 1
    assert out[0] is a, "the most supported end is the representative"
    assert out[0].n_reads == 70 and out[0].reads.size == 70
    out, _ = collapse_3p([mk(FOUR, n=50), mk(FOUR[:-1] + [(5000, 5750)], n=20)],
                         100, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2, "150 bp apart must survive N = 100"


def test_annotated_end_wins_over_support():
    a = mk(FOUR, n=10, end_annotated=1)
    b = mk(FOUR[:-1] + [(5000, 5700)], n=90, end_annotated=0)
    out, _ = collapse_3p([b, a], 200, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1 and out[0] is a and out[0].n_reads == 100


def test_different_chains_never_merge():
    a = mk(FOUR, n=50)
    b = mk([(1000, 1200), (2000, 2130)] + FOUR[2:], n=50)  # one donor moved 30 bp
    out, st = collapse_3p([a, b], 2000, 2000, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2 and st["end3_models_absorbed"] == 0


def test_5p_window_is_respected():
    a = mk(FOUR, n=50)
    b = mk([(700, 1200)] + FOUR[1:], n=20)               # 5' end 300 bp upstream
    out, _ = collapse_3p([a, b], 200, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2


def test_minus_strand_uses_transcript_orientation():
    a = mk(FOUR, strand="-", n=50)
    b = mk([(850, 1200)] + FOUR[1:], strand="-", n=20)   # 3' (left) end 150 bp further
    out, _ = collapse_3p([a, b], 200, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1


# -------------------------------------------------------------- sub-chain rule
def _filter(models, **kw):
    p = dict(ratio=2.0, fuzzy=5, slack=50, protected=["FSM"], fold_reads=True,
             reads=None, umi_hamming=0, inferred_category="INFERRED", contig="chr1", audit=[])
    p.update(kw)
    return filter_subchains(models, **p)


def test_truncated_fragment_is_folded_into_its_container():
    full = mk(FOUR, n=100)
    frag = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=10)   # last two introns
    out, st = _filter([full, frag])
    assert out == [full] and st["subchain_removed"] == 1
    assert full.n_reads == 110


def test_well_supported_subchain_survives_the_ratio():
    full = mk(FOUR, n=100)
    alt = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=60)
    out, _ = _filter([full, alt])
    assert len(out) == 2


def test_alternative_first_exon_inside_an_intron_is_not_a_fragment():
    full = mk(FOUR, n=100)
    alt_tss = mk([(2600, 3100), (4000, 4100), (5000, 5600)], n=5)  # starts 400 bp into intron
    out, _ = _filter([full, alt_tss])
    assert len(out) == 2


def test_annotated_chain_is_protected():
    full = mk(FOUR, n=100)
    fsm = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=5, category="FSM")
    out, _ = _filter([full, fsm])
    assert len(out) == 2


def test_monoexon_inside_an_exon_is_a_fragment_but_not_across_an_intron():
    full = mk(FOUR, n=100)
    inside = mk([(5100, 5600)], n=5, category="monoexon_3UTR")
    across = mk([(4050, 5600)], n=5, category="monoexon_3UTR")
    out, _ = _filter([full, inside, across])
    assert inside not in out and across in out


def test_inferred_models_are_left_alone():
    full = mk(FOUR, n=100)
    inf = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=0, category="INFERRED",
             flags=["INFERRED"])
    out, _ = _filter([full, inf])
    assert inf in out


def test_fragment_chains_resolve_to_the_surviving_container():
    full = mk(FOUR, n=400)
    mid = mk([(2050, 2100), (3000, 3100), (4000, 4100), (5000, 5600)], n=100)
    tip = mk([(4050, 4100), (5000, 5600)], n=10)
    out, st = _filter([full, mid, tip])
    assert out == [full] and full.n_reads == 510 and st["subchain_reads_folded"] == 110


# ----------------------------------------------------------------- end to end
def test_run_posthoc_defaults_and_audit():
    full = mk(FOUR, n=100, end_annotated=1)
    tail = mk(FOUR[:-1] + [(5000, 5700)], n=30)
    frag = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=10)
    res = run_posthoc([full, tail, frag], "chr1", None, cfg())
    assert len(res.models) == 1 and res.models[0].n_reads == 140
    rules = sorted(r["rule"] for r in res.audit)
    assert rules == ["end3_collapse", "subchain_fragment"]
    assert all(r["into_model_index"] == 0 for r in res.audit)


def test_disabled_is_a_no_op():
    ms = [mk(FOUR, n=100), mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=10)]
    res = run_posthoc(ms, "chr1", None, cfg(enabled=False))
    assert res.models == ms


def test_negative_tolerance_switches_the_3p_rule_off():
    ms = [mk(FOUR, n=50), mk(FOUR[:-1] + [(5000, 5650)], n=20)]
    res = run_posthoc(ms, "chr1", None, cfg(end3_tolerance=-1, subchain_filter=False))
    assert len(res.models) == 2


def test_config_validation_and_cli_flags():
    c = Config()
    c.posthoc.subchain_ratio = 0.5
    with pytest.raises(ValueError, match="subchain_ratio"):
        c.validate()
    from flightcollapse.cli import _build_config
    import argparse
    ns = argparse.Namespace(end3_tolerance=500, no_subchain_filter=True, subchain_ratio=3.0,
                            no_posthoc=False, overrides=["posthoc.end5_window=50"])
    c = _build_config(ns)
    assert c.posthoc.end3_tolerance == 500 and not c.posthoc.subchain_filter
    assert c.posthoc.subchain_ratio == 3.0 and c.posthoc.end5_window == 50
    assert Config().posthoc.end3_tolerance == 200, "the documented default"
    assert Config().posthoc.end5_window == -1, "the 5' guard is off by default"
    assert Config().posthoc.subchain_filter, "the fragment filter is on by default"


def test_pipeline_writes_the_posthoc_audit(sim, tmp_path_factory):
    from flightcollapse.pipeline import run

    out = tmp_path_factory.mktemp("posthoc")
    c = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    c.molecules.barcode_umi_tsv = sim.barcode_umi
    c.output.outdir, c.output.prefix = str(out), "ph"
    c.verbose, c.strict_invariants = False, False
    rep = run(c)
    assert (out / "ph.posthoc.tsv").exists()
    ph = rep["posthoc"]
    assert ph["end3_tolerance"] == Config().posthoc.end3_tolerance
    assert ph["models_after"] == rep["n_models"]
    assert ph["models_before"] >= ph["models_after"]


# ------------------------------------------------- rule 1, the rest of the list
def test_opposite_strands_never_merge():
    a = mk(FOUR, strand="+", n=50)
    b = mk(FOUR, strand="-", n=20)
    out, st = collapse_3p([a, b], 2000, 2000, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2 and st["end3_models_absorbed"] == 0


def test_support_wins_over_length():
    """Neither end annotated, so support decides -- not the longer transcript."""
    short_rich = mk(FOUR, n=90)
    long_poor = mk(FOUR[:-1] + [(5000, 5800)], n=10)
    out, _ = collapse_3p([long_poor, short_rich], 500, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1 and out[0] is short_rich
    assert out[0].end == 5600, "the survivor keeps its own coordinates"


def test_zero_tolerance_merges_only_identical_3p_ends():
    same = [mk(FOUR, n=50), mk(FOUR, n=20)]
    out, _ = collapse_3p(same, 0, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1
    off_by_one = [mk(FOUR, n=50), mk(FOUR[:-1] + [(5000, 5601)], n=20)]
    out, _ = collapse_3p(off_by_one, 0, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2


def test_anchored_not_chained():
    """3' ends 0 / 400 / 800 at N = 500: the far one starts its own cluster.

    Single linkage would absorb all three, and the cluster would then span 800
    bp at a 500 bp tolerance. The span of a cluster has to be bounded by the
    parameter, not by how densely the ends happen to be packed.
    """
    rep = mk(FOUR, n=100)                                  # 3' end 5600
    near = mk(FOUR[:-1] + [(5000, 6000)], n=50)            # +400
    far = mk(FOUR[:-1] + [(5000, 6400)], n=40)             # +800
    out, st = collapse_3p([rep, near, far], 500, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2 and st["end3_models_absorbed"] == 1
    survivors = {m.end for m in out}
    assert survivors == {5600, 6400}
    assert rep.n_reads == 150 and far.n_reads == 40


def test_reads_move_and_molecules_are_recounted():
    """n_mols must be recomputed from the moved reads, not added up."""
    from flightcollapse.secondpass import NO_CELL, NO_UMI

    class _Reads:                                   # the two models share a cell
        cell = np.array([7, 7, 7, 7], dtype=np.int64)
        umi = np.array([1, 1, 2, 2], dtype=np.int64)

    a = TranscriptModel(model_id="", contig="chr1", strand="+", exons=list(FOUR),
                        chain=tuple((FOUR[i][1], FOUR[i + 1][0]) for i in range(4)),
                        category="NIC", reads=np.array([0, 1], dtype=np.int64),
                        n_reads=2, n_mols=1)
    b = TranscriptModel(model_id="", contig="chr1", strand="+",
                        exons=FOUR[:-1] + [(5000, 5700)],
                        chain=tuple((FOUR[i][1], FOUR[i + 1][0]) for i in range(4)),
                        category="NIC", reads=np.array([2, 3], dtype=np.int64),
                        n_reads=2, n_mols=1)
    out, _ = collapse_3p([a, b], 200, 100, _Reads(), 0, "INFERRED", "chr1", [])
    assert len(out) == 1 and out[0].n_reads == 4
    assert out[0].n_mols == 2, "two UMIs in one cell, not 1 + 1 summed blindly"
    assert out[0].evidence["posthoc_end3_absorbed"] == 1


# ------------------------------------------------- rule 2, the rest of the list
def test_non_contiguous_subchain_is_kept():
    """Sharing introns 1 and 3 but not 2 is exon skipping, not truncation."""
    full = mk(FOUR, n=200)
    skipper = mk([(1000, 1200), (2000, 2100), (4000, 4100), (5000, 5600)], n=5)
    out, _ = _filter([full, skipper])
    assert len(out) == 2


def test_ratio_boundary():
    frag = [(3050, 3100), (4000, 4100), (5000, 5600)]
    just_under, _ = _filter([mk(FOUR, n=19), mk(frag, n=10)])
    assert len(just_under) == 2, "1.9x is not a container"
    exactly, _ = _filter([mk(FOUR, n=20), mk(frag, n=10)])
    assert len(exactly) == 1, "2.0x is"


def test_unfolded_fragment_reads_are_counted_not_lost():
    full = mk(FOUR, n=100)
    frag = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=7)
    out, st = _filter([full, frag], fold_reads=False)
    assert out == [full] and full.n_reads == 100
    assert st["subchain_reads_unassigned"] == 7 and st["subchain_reads_folded"] == 0


def test_indexed_filter_matches_the_brute_force_definition():
    """The index is an optimisation; it must decide exactly what the definition does.

    The filter looks a candidate's first intron up in a positional index rather
    than comparing all pairs, because the all-pairs form is minutes per contig
    on a real catalogue. This checks the two agree on random catalogues built
    from the shapes the rule is about: truncations, alternative first exons,
    exon skips and mono-exons.
    """
    import random

    from flightcollapse.posthoc import _is_fragment_of, _support

    def brute(models, ratio, fuzzy, slack, protected):
        out = set()
        for a in models:
            if a.category in protected or "INFERRED" in a.flags:
                continue
            sa = _support(a)
            for b in models:
                if "INFERRED" in b.flags:
                    continue
                if _support(b) >= ratio * max(sa, 1) and _is_fragment_of(a, b, fuzzy, slack):
                    out.add(id(a))
                    break
        return out

    rng = random.Random(20260605)
    for trial in range(40):
        models = []
        for g in range(rng.randint(2, 6)):
            base = 10_000 + g * 40_000
            k = rng.randint(2, 6)
            ex = [(base, base + 300)]
            for i in range(k):
                s = base + 2_000 * (i + 1)
                ex.append((s, s + rng.randint(80, 400)))
            strand = rng.choice("+-")
            models.append(mk(ex, strand=strand, n=rng.randint(40, 200)))
            for _ in range(rng.randint(0, 3)):
                cut = rng.randint(1, max(1, len(ex) - 2))
                sub = [(ex[cut][0] + rng.randint(0, 60), ex[cut][1])] + ex[cut + 1:]
                if len(sub) >= 2:
                    models.append(mk(sub, strand=strand, n=rng.randint(1, 120),
                                     category=rng.choice(["NIC", "NNC", "FSM"])))
            if rng.random() < .5:                   # alternative first exon
                alt = [(ex[0][1] + 400, ex[1][1])] + ex[2:]
                if len(alt) >= 2:
                    models.append(mk(alt, strand=strand, n=rng.randint(1, 50)))
            if rng.random() < .6:                   # a mono-exon inside an exon
                x, y = ex[rng.randrange(len(ex))]
                if y - x > 40:
                    models.append(mk([(x + 10, y - 10)], strand=strand,
                                     n=rng.randint(1, 60), category="monoexon_3UTR"))
        rng.shuffle(models)
        want = brute(models, 2.0, 5, 50, {"FSM"})
        kept, _ = _filter(list(models), fold_reads=False)
        got = {id(m) for m in models} - {id(m) for m in kept}
        assert got == want, f"trial {trial}: {len(got)} removed vs {len(want)} expected"


# ------------------------------------------------------------------ integration
_DATA_OUTPUTS = (
    "x.gtf", "x.gff", "x.models.tsv", "x.abundance.txt", "x.flnc_count.txt",
    "x.group_chains.tsv", "x_umi_corrected_isoform_matrix.mtx", "x_cells.txt",
    "x_isoforms.txt", "x_isoform_metadata.txt",
)


def _run(sim, outdir, mono_demotion=True, **ph):
    from flightcollapse.pipeline import run

    c = Config(bam=sim.bam, reference_gtf=sim.gtf, genome_fasta=sim.genome_fasta)
    c.molecules.barcode_umi_tsv = sim.barcode_umi
    c.output.outdir, c.output.prefix = str(outdir), "x"
    c.verbose, c.strict_invariants = False, False
    if not mono_demotion:
        # The first pass and the second pass each already remove the
        # simulation's fragments, so with defaults step 8b has nothing left to
        # do and an end-to-end test of it would assert nothing. Turning both
        # off leaves the fragment standing, which is the only way to watch the
        # sub-chain rule act on a real pipeline run rather than on a fixture.
        c.monoexon.demote_terminal_exon_fragments = False
        c.monoexon.demote_contained_fragments = False
        c.secondpass.enabled = False
    for k, v in ph.items():
        setattr(c.posthoc, k, v)
    return run(c)


def test_no_posthoc_is_byte_identical_to_not_running_the_step(sim, tmp_path_factory,
                                                              monkeypatch):
    """`--no-posthoc` must reproduce 0.5.x output, not merely resemble it.

    The 0.5.x pipeline had no step 8b at all, so the check is against a run
    where the step is physically bypassed: if a disabled run and a bypassed run
    differ in any byte of any data file, the flag is not a true escape hatch.
    """
    import gzip
    import hashlib
    import os

    from flightcollapse import posthoc as phmod
    from flightcollapse import pipeline as plmod

    a = tmp_path_factory.mktemp("ph_off")
    _run(sim, a, enabled=False)

    b = tmp_path_factory.mktemp("ph_bypassed")
    monkeypatch.setattr(plmod.phmod, "run_posthoc",
                        lambda models, contig, reads, cfg: phmod.PosthocResult(models=models))
    _run(sim, b)                      # posthoc ENABLED, but the step does nothing

    def digest(d, name):
        # gzip writes the wall clock into its header, so a .gz container
        # differs between two runs even when the stream is identical. Compare
        # what the file says, not what second it was written in.
        p = os.path.join(d, name)
        assert os.path.exists(p), f"{name} missing from {d}"
        op = gzip.open if name.endswith(".gz") else open
        with op(p, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    for name in _DATA_OUTPUTS:
        assert digest(a, name) == digest(b, name), name
    for name in ("x.read_stat.txt.gz", "x.read_stat.txt",
                 "x.group.txt.gz", "x.group.txt", "x.junctions.tsv.gz"):
        if os.path.exists(os.path.join(a, name)):
            assert digest(a, name) == digest(b, name), name
    assert not os.path.exists(os.path.join(a, "x.posthoc.tsv")), \
        "a disabled run must not write an audit"


def test_outputs_describe_the_consolidated_catalogue(sim, tmp_path_factory):
    """Every downstream table must name only models that survived step 8b."""
    import csv as _csv
    import gzip as _gzip
    import os

    import pandas as pd

    out = tmp_path_factory.mktemp("ph_consistent")
    rep = _run(sim, out, mono_demotion=False)
    assert rep["posthoc"]["models_before"] > rep["posthoc"]["models_after"], \
        "this fixture exists to watch step 8b remove something"
    tbl = pd.read_csv(out / "x.models.tsv", sep="\t")
    live = set(tbl["model_id"])
    assert len(live) == rep["n_models"] == rep["posthoc"]["models_after"]

    p = out / "x.read_stat.txt.gz"
    op = _gzip.open if p.exists() else open
    rs = pd.read_csv(op(p if p.exists() else out / "x.read_stat.txt", "rt"), sep="\t")
    assert set(rs["pbid"]) <= live, "a read is assigned to a model that was removed"
    counts = rs["pbid"].value_counts()
    for r in tbl.itertuples(index=False):
        assert counts.get(r.model_id, 0) == r.n_reads, r.model_id

    isos = (out / "x_isoforms.txt").read_text().split()
    assert set(isos) <= live

    audit = out / "x.posthoc.tsv"
    assert audit.exists()
    with open(audit) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert rows, "something was removed, so the audit cannot be empty"
    n_removed = (rep["posthoc"]["end3_models_absorbed"]
                 + rep["posthoc"]["subchain_removed"])
    assert len(rows) == n_removed == rep["posthoc"]["models_before"] - len(live)
    from flightcollapse.posthoc import AUDIT_COLUMNS
    with open(audit) as fh:
        assert fh.readline().rstrip("\n").split("\t") == list(AUDIT_COLUMNS)
    for r in rows:
        assert r["rule"] in ("end3_collapse", "subchain_fragment")
        if int(r["into_model_index"]) >= 0:
            assert r["into_model_id"] in live, \
                "the audit points at a model that is not in the catalogue"


def test_the_two_rules_run_in_the_documented_order(sim, tmp_path_factory):
    """Rule 2's ratio has to see consolidated counts, so rule 1 goes first.

    Two same-chain models each at 10 reads consolidate to 20, which then clears
    the 2x bar against a 15-read fragment that would have survived against
    either half alone.
    """
    base = [(1000, 1200), (2000, 2100), (3000, 3100), (5000, 5600)]
    full_a = mk(base, n=10)
    full_b = mk(base[:-1] + [(5000, 5900)], n=10)   # same chain, 3' end 300 further
    frag = mk([(3050, 3100), (5000, 5600)], n=9)    # the last intron only
    # 500, not the default: the fixture's two 3' ends are 300 bp apart on
    # purpose, and this test is about the ORDER of the two rules, not about
    # what the default tolerance happens to be.
    res = run_posthoc([full_a, full_b, frag], "chr1", None, cfg(end3_tolerance=500))
    assert len(res.models) == 1
    assert res.models[0] is full_b, "equal support, so the longer end breaks the tie"
    assert res.models[0].n_reads == 29, "10 + 10 consolidated, then the fragment's 9"
    # with rule 1 off, neither half alone clears 2 x 9 and the fragment survives
    res = run_posthoc([mk(base, n=10), mk(base[:-1] + [(5000, 5900)], n=10),
                       mk([(3050, 3100), (5000, 5600)], n=9)],
                      "chr1", None, cfg(end3_tolerance=-1))
    assert len(res.models) == 3


def test_a_fragment_never_crosses_a_contig():
    """The pipeline calls the filter per contig; the function must not rely on it.

    Two contigs share coordinates. An index bucketed on position alone would
    offer a chr2 model as the container for an identically-placed chr1 one, and
    the 0.6.0 sweep -- which applies the rules to a whole catalogue at once --
    would quietly delete real models.
    """
    full = mk(FOUR, n=100)
    frag = mk([(3050, 3100), (4000, 4100), (5000, 5600)], n=10)
    frag.contig = "chr2"
    out, st = _filter([full, frag])
    assert len(out) == 2 and st["subchain_removed"] == 0
    mono = mk([(5100, 5600)], n=5, category="monoexon_3UTR")
    mono.contig = "chr2"
    out, st = _filter([full, mono])
    assert len(out) == 2 and st["subchain_removed"] == 0
    # and the same pair on one contig is still removed, so the guard is not a
    # blanket refusal
    frag.contig = "chr1"
    out, _ = _filter([full, frag])
    assert out == [full]


def test_a_negative_5p_window_removes_the_guard():
    """The 3' end becomes the whole rule, and 5' truncation is the filter's job.

    Two models that share a chain exactly can only differ at the 5' inside their
    first exon, so the fixture needs a long one. Such a pair is a truncation the
    sub-chain filter will never see -- it only compares a model whose chain is a
    STRICT sub-chain of another. Either the 3' rule merges them or nothing does.
    """
    LONG = [(1000, 6000), (8000, 8100), (9000, 9100), (10000, 10600)]
    rep = mk(LONG, n=100)
    deep5 = mk([(5000, 6000)] + LONG[1:], n=20)          # 5' end 4 kb later
    from flightcollapse.posthoc import _genomic_introns
    assert _genomic_introns(deep5) == _genomic_introns(rep), "same chain"
    out, _ = collapse_3p([rep, deep5], 200, 100, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2, "the 100 bp guard blocks it"
    rep = mk(LONG, n=100)
    deep5 = mk([(5000, 6000)] + LONG[1:], n=20)
    out, _ = collapse_3p([rep, deep5], 200, -1, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 1 and out[0] is rep and out[0].n_reads == 120
    # the 3' end still rules, however open the 5' is
    far = mk(LONG[:-1] + [(10000, 11600)], n=20)
    out, _ = collapse_3p([mk(LONG, n=100), far], 200, -1, None, 0, "INFERRED", "chr1", [])
    assert len(out) == 2


def test_config_rejects_an_open_5p_window_with_no_fragment_filter():
    """Nothing would be looking at the 5' end at all, which is not a setting."""
    def problems(**ph):
        c = Config(bam="x.bam", reference_gtf="x.gtf")
        for k, v in ph.items():
            setattr(c.posthoc, k, v)
        try:
            c.validate()
        except ValueError as e:
            return str(e)
        return ""

    assert "nothing is looking at the 5' end" in problems(
        end5_window=-1, subchain_filter=False)
    assert "nothing is looking at the 5' end" not in problems(
        end5_window=-1, subchain_filter=True)
    assert "nothing is looking at the 5' end" not in problems(
        end5_window=100, subchain_filter=False)
