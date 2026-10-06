#!/usr/bin/env python3
"""Bring the flightcollapse logic sheet from 0.5.1 to 0.6.3.

The 0.5.1 sheet is the detailed figure: one card per DECISION, each drawn as the
exon geometry the decision is about. This adds the third band — the 0.6.x
consolidation — in the same idiom, and updates the header, the legend and the
invariants.
"""
import pathlib, re

SRC = pathlib.Path(__file__).parent / "flightcollapse_logic_v0.5.1.html"
OUT = pathlib.Path(__file__).parent / "flightcollapse_logic_v0.6.3.html"

GREY, BLUE, ORNG, GREEN, PURP, RED, AMBR = (
    "#b9b6b0", "#2978d5", "#c2521f", "#0e8f62", "#4a3aa6", "#c7342e", "#a25900")


def ex(w, c=GREY):
    return f'<div class="ex" style="width:{w}px;background:{c}"></div>'


def iv(c=""):
    return f'<div class="iv {c}"></div>'


def gap(w):
    return f'<div style="width:{w}px;flex:none"></div>'


def trk(label, *parts):
    return f'<div class="trk"><div class="lbl">{label}</div>' + "".join(parts) + "</div>"


def dia(*rows, cap=None, capcls=""):
    c = f'<div class="cap {capcls}">{cap}</div>' if cap else ""
    # .cmp keeps the tracks compact: this band's columns are three times the
    # width of the first band's, and .iv is flex:1, so an unconstrained track
    # stretches its introns across the whole card and the exons vanish
    return '<div class="dia cmp">' + "".join(rows) + c + "</div>"


def card(title, *body, cls="", tag=None, tagcls="", sub=None, hcls=""):
    h = f'<div class="h {hcls}">{title}</div>' if title else ""
    t = f'<div class="s">{sub}</div>' if sub else ""
    g = f'<div class="tag {tagcls}">{tag}</div>' if tag else ""
    return f'<div class="p {cls}">{h}' + "".join(body) + t + g + "</div>"


def text(t):
    return f'<div class="t">{t}</div>'


def head(n, title, acc=True):
    a = " a" if acc else ""
    return (f'<div class="ph"><div class="num{a}">{n}</div>'
            f'<div class="ttl">{title}</div></div>')


# ---------------------------------------------------------------- panel 11
CHAIN = (ex(30), iv(), ex(40), iv(), ex(26))


def panel_11():
    cards = [
        f'<div class="p la">{head(11, "Rule 1 — same chain, 3′ ends within 200 bp")}'
        + text("Models are partitioned by contig, strand and <b>exact</b> intron chain, "
               "so the rule can never merge two splice structures. Within a partition "
               "the representative is the annotated 3′ end if there is one, then the "
               "best supported, then the longest. Everything within the tolerance "
               "<b>of that representative</b> is absorbed.")
        + '<div class="s">new default 200 bp in 0.6.3 · 500 bp in 0.6.0–0.6.2 · '
          '0 merges only identical ends · −1 switches the rule off</div></div>',

        card("identical chain, 3′ ends 150 bp apart &nbsp;→&nbsp; MERGE",
             dia(trk("GENCODE", *CHAIN, ex(18)),
                 trk("model A", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model B", ex(30, ORNG), iv("o"), ex(40, ORNG), iv("o"), ex(34, ORNG)),
                 cap="B's 3′ end is 150 bp past A's · A carries the annotated end"),
             text("B's reads move onto A and B's molecules are recounted from them. "
                  "Nothing is discarded; the catalogue is one model shorter."),
             cls="lr", tag="posthoc_end3_absorbed = 1", tagcls="o"),

        card("identical chain, 3′ ends 400 bp apart &nbsp;→&nbsp; KEPT",
             dia(trk("model A", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model B", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(58, BLUE)),
                 cap="beyond the tolerance — two 3′ ends, two models"),
             text("At 500 bp these would have merged. That is the whole of the "
                  "difference between the two defaults: 26,080 further merges across "
                  "the three HEK libraries, bought with 0.79 pp of reproducible read "
                  "mass."),
             cls="l"),

        card("anchored, not chained",
             dia(trk("A ·rep", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("B", ex(30, ORNG), iv("o"), ex(40, ORNG), iv("o"), ex(36, ORNG)),
                 trk("C", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(46, BLUE)),
                 cap="A→B 180 bp · B→C 180 bp · A→C 360 bp"),
             text("B is absorbed because it is within 200 bp <b>of A</b>. C is not, "
                  "even though it is within 200 bp of B. Single-linkage chaining would "
                  "let one cluster span the whole locus; anchoring to the "
                  "representative bounds every cluster at the tolerance."),
             cls="l", tag="C starts the next cluster"),

        card("one junction differs &nbsp;→&nbsp; NEVER COMPARED",
             dia(trk("model A", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model B", ex(30, PURP), iv("p"), ex(17, PURP), iv("p"),
                     ex(17, PURP), iv("p"), ex(26, PURP)),
                 cap="one extra junction — a different partition"),
             text("The partition key is the exact intron chain, so no tolerance, "
                  "however large, can merge two splice structures."),
             cls="lp"),

        card("INFERRED models are untouched",
             text("A model the second pass substituted for a rejected one is neither "
                  "absorbed nor used as a representative. It is not evidence about a 3′ "
                  "end — it is a stand-in for a model that failed."),
             cls="plain", tag="category = INFERRED", tagcls="a"),
    ]
    return '<div class="col w12">' + "".join(cards) + "</div>"


# ---------------------------------------------------------------- panel 12
def panel_12():
    cards = [
        f'<div class="p la">{head(12, "Rule 2 — reference-free sub-chain (fragment) filter")}'
        + text("A model is removed, and its reads folded into the container, when its "
               "ordered intron chain is a <b>contiguous</b> sub-chain of a longer "
               "model on the same contig and strand (junctions ±5 bp), its terminal "
               "exons lie inside the matching exons of that model (±50 bp), and the "
               "container carries ≥ 2× its support.")
        + '<div class="s">no annotation is consulted at any point · '
          'decisions are taken in one pass on pre-filter supports</div></div>',

        card("5′-truncated copy &nbsp;→&nbsp; FOLDED",
             dia(trk("container", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("fragment", gap(34), ex(40, ORNG), iv("o"), ex(26, ORNG)),
                 cap="last two exons, both terminal exons inside the container's"),
             text("The classic 5′ truncation: a shorter chain that is a suffix of a "
                  "longer one, with nothing of its own at either end."),
             cls="lr", tag="posthoc_fragments_absorbed", tagcls="o"),

        card("alternative first exon inside an intron &nbsp;→&nbsp; KEPT",
             dia(trk("container", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model", gap(16), ex(22, PURP), iv("p"), ex(40, PURP), iv("p"),
                     ex(26, PURP)),
                 cap="its first exon starts in the container's intron — outside the "
                     "±50 bp slack"),
             text("This is what the terminal-exon slack is for. A chain that is a "
                  "sub-chain but whose own first exon starts somewhere the container "
                  "has no exon is an alternative start, not a truncation."),
             cls="lp", tag="kept as its own model", tagcls="p"),

        card("mono-exon inside one exon &nbsp;→&nbsp; FOLDED",
             dia(trk("container", ex(30, BLUE), iv("b"), ex(54, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model", gap(44), ex(34, ORNG), gap(32)),
                 cap="wholly inside one exon of a 2×-better model"),
             text("A single-exon model has no chain, so containment is judged on the "
                  "exon alone."),
             cls="lr"),

        card("exact annotated chain (FSM) &nbsp;→&nbsp; PROTECTED",
             dia(trk("GENCODE", ex(30), iv(), ex(40), iv(), ex(26)),
                 trk("model", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 cap="a transcript the annotation asserts exists"),
             text("FSM is never removed, whatever contains it. A reference-free rule "
                  "deleting a reference transcript would be the filter overreaching "
                  "into territory it has no evidence about."),
             cls="l", tag="protected_categories = [FSM]"),

        card("container only 1.6× better &nbsp;→&nbsp; BOTH KEPT",
             dia(trk("container", ex(30, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model", gap(34), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 cap="50 molecules against 32 — below the 2× bar"),
             text("Comparable support means the shorter form is not obviously debris. "
                  "The ratio is the whole of the evidence the rule has, so it is held "
                  "to deliberately."),
             cls="l", tag="subchain_ratio = 2.0"),

        card("a container that is itself removed",
             text("Decisions are taken in one pass on pre-filter supports, then "
                  "resolved: a fragment whose container was itself folded away lands on "
                  "the final survivor, with cycle protection. Reads never stop at a "
                  "model that is not in the catalogue."),
             cls="plain", tag="into_model_id in {prefix}.posthoc.tsv"),
    ]
    return '<div class="col w12">' + "".join(cards) + "</div>"


# ---------------------------------------------------------------- panel 13
def panel_13():
    rows = [
        ("models", "499,553 → 418,097", "499,553 → 392,017"),
        ("", "−16.3%", "−21.5%"),
        ("reproducible read mass", "−0.32 pp", "−1.11 pp"),
        ("at ≥50 reads, ≥100 bc", "+0.11 pp", "−0.25 pp"),
        ("cells where mass rises", "34 of 81", "17 of 81"),
        ("splice chains", "122,270", "122,174"),
    ]
    tbl = "".join(
        f'<div class="kv{" hi" if i == 2 else ""}"><span>{k}</span>'
        f'<span>{a} &nbsp;·&nbsp; {b}</span></div>'
        for i, (k, a, b) in enumerate(rows))

    cards = [
        f'<div class="p la">{head(13, "Why the two ends are not treated alike")}'
        + text("The 3′ end is a <b>position</b> — the mode of a cleavage-site peak — "
               "and only a tolerance can say whether two of them are the same site. "
               "The 5′ end is a <b>containment relation</b> between two intron chains, "
               "which the sub-chain filter decides from structure alone, with no "
               "annotation and no distance threshold.")
        + '<div class="s">this asymmetry is the design, not an accident of tuning</div></div>',

        card("what the guard being off costs, and buys",
             dia(trk("model A", ex(22, BLUE), iv("b"), ex(40, BLUE), iv("b"), ex(26, BLUE)),
                 trk("model B", gap(10), ex(14, ORNG), iv("o"), ex(40, ORNG), iv("o"),
                     ex(26, ORNG)),
                 cap="same chain, same 3′ end, different TSS — these now merge"),
             text("With <span class=\"mono\">end5_window</span> on, that pair never "
                  "merged. But neither did a 5′-<b>truncated</b> copy, because it "
                  "shares the chain exactly and so is invisible to the sub-chain "
                  "filter too. 0.6.3 accepts the first case to fix the second: across "
                  "the three libraries it freed 24–5,664 merges per arm for at most "
                  "0.03 pp of mass."),
             cls="la", hcls="a",
             tag="posthoc.end5_window = −1  ·  =100 restores the guard", tagcls="a"),

        card("200 bp against 500 bp, on three libraries",
             f'<div class="dia">{tbl}</div>',
             text("Read mass, not the count-based share: merging two structures that "
                  "were both found in all three removes one from the numerator and one "
                  "from the denominator, so a share always falls. The chain-level row "
                  "is the one that settles it — the two tolerances are "
                  "indistinguishable there, which says the 3′ tolerance is doing 3′ "
                  "work only."),
             cls="l", sub="BD144_kin / HM / ND · structures matched on chain ±5 bp and "
                          "3′ end ±100 bp"),

        card("an outside check nobody designed for",
             '<div class="dia">'
             + '<div class="kv"><span>ISM, 3′ fragment</span><span>6,542 → 2,649 &nbsp;−60%</span></div>'
             + '<div class="kv"><span>ISM, 5′ fragment</span><span>20,590 → 11,616 &nbsp;−44%</span></div>'
             + '<div class="kv"><span>FSM, alt 3′ end</span><span>44,985 → 37,469 &nbsp;−17%</span></div>'
             + '<div class="kv"><span>FSM (reference match)</span><span>33,018 → 31,581 &nbsp;−4%</span></div>'
             + '<div class="cap">pigeon / SQANTI3 classes, BD144_kin</div></div>',
             text("The sub-chain filter never consults an annotation, and pigeon is "
                  "nothing but an annotation — yet the classes it empties are exactly "
                  "the fragment classes, while the 3′ tolerance's signature appears in "
                  "<i>FSM with an alternative 3′ end</i>. A rule cutting real "
                  "transcripts would not sort itself along a reference it cannot see."),
             cls="lg", tag="reference-free rule, reference-based verdict", tagcls="g"),

        card("where the stage sits",
             text("Step 8b: last on every contig, <b>after</b> the second pass and "
                  "<b>before</b> read bookkeeping — so "
                  "<span class=\"mono\">group.txt</span>, "
                  "<span class=\"mono\">read_stat</span>, the abundance tables and the "
                  "cell × isoform matrix all describe the catalogue that is actually "
                  "written. The two rules run 3′ first, so the ratio in rule 2 is "
                  "tested on consolidated counts."),
             cls="plain", tag="--no-posthoc reproduces 0.5.1 byte for byte"),
    ]
    return '<div class="col">' + "".join(cards) + "</div>"


BAND = f'''
<div class="band">
  <div class="bandhead b4">
    <span class="t">CONSOLIDATION — new in 0.6.0, retuned in 0.6.3</span>
    <span class="d">runs last on every contig · both rules only MOVE support between
    existing models, so the read totals are identical with the stage on and off</span>
  </div>
  <div class="cols">
    {panel_11()}
    {panel_12()}
    {panel_13()}
  </div>
</div>
'''


def main():
    s = SRC.read_text()

    s = s.replace("<title>flightcollapse 0.5.1 logic</title>",
                  "<title>flightcollapse 0.6.3 logic</title>")
    s = s.replace('<span class="ver">v0.5.1</span>', '<span class="ver">v0.6.3</span>')

    # a fourth band colour, matching the existing three
    s = s.replace(".b3{background:#f1faf5;border-left:4px solid #0e8f62}",
                  ".b3{background:#f1faf5;border-left:4px solid #0e8f62}\n"
                  ".b4{background:#f3f1fb;border-left:4px solid #4a3aa6}\n"
                  ".dia.cmp .trk{max-width:318px}\n"
                  ".dia.cmp .kv{max-width:none}")

    # legend: the 0.5.x pill becomes the 0.6.x one
    s = s.replace('<span><span class="sw" style="background:#a25900"></span>new in 0.5.x</span>',
                  '<span><span class="sw" style="background:#a25900"></span>'
                  'new or retuned in 0.6.x</span>')
    s = s.replace('<span class="pill">panels 0, 8, 9 are new</span>',
                  '<span class="pill">panels 11–13 are new</span>')

    # the orientation note
    old_note = ("0.5.x changes NOTHING in panels 1–6. It adds a second pass over the\n"
                "finished catalogue (panels 8–9) and splits the run across contigs (panel 0).")
    new_note = ("0.6.x changes NOTHING in panels 1–9. It adds a final consolidation of the\n"
                "finished catalogue (panels 11–13), which merges same-chain models at the 3′ end\n"
                "and folds 5′ fragments into their containers. "
                "<span class=\"mono\">--no-posthoc</span> reproduces 0.5.1 exactly.")
    if old_note in s:
        s = s.replace(old_note, new_note, 1)
    else:
        print("  ! orientation note not matched verbatim — check the header")

    # the conservation invariant belongs with the others
    inv = ('<li><b>mono-exon read mass = unspliced read fraction</b>'
           if '<li><b>mono-exon read mass' in s else None)
    marker = "Invariants — checked every run</div>"
    i = s.index(marker) + len(marker)
    j = s.index("<ul", i) + s[s.index("<ul", i):].index(">") + 1
    s = (s[:j] + '<li><b>consolidation moves support, never drops it</b> — '
         'reads assigned and reads unassigned are identical to the digit with the '
         'stage on and off, in all three HEK libraries. A difference is a bug in the '
         'fold, not a tuning question.</li>' + s[j:])

    # third band, before the footer
    k = s.index('<div class="foot">')
    s = s[:k] + BAND + "\n" + s[k:]

    s = s.replace("flightcollapse v0.5.1 · snapshot", "flightcollapse v0.6.3 · snapshot")
    s = s.replace("flightcollapse 0.5.1 · snapshot", "flightcollapse 0.6.3 · snapshot")

    OUT.write_text(s if s.endswith("\n") else s + "\n")
    print(f"wrote {OUT} ({len(s):,} bytes)")


if __name__ == "__main__":
    main()
