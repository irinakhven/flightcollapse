#!/usr/bin/env python3
"""Draw the two flightcollapse schematics.

Laid out on a grid in code rather than by hand, so the baselines and gaps are
exact and the figures can be regenerated when the algorithm changes.

  flightcollapse_algorithm_distilled.svg  the pipeline in five steps

The detailed figure is a different kind of drawing — one card per decision, each
drawn as the exon geometry the decision is about — and lives in
docs/flightcollapse_logic_v0.6.3.html, built by docs/make_logic_sheet.py.

Run it from the repository root after changing a stage or a default:

    python3 docs/make_schematic.py

Both inherit the page's foreground through currentColor and carry one accent
for the stage the current release changed.
"""
import pathlib, html

OUT = pathlib.Path(__file__).parent

ACC = "var(--acc, #2a78d6)"      # the consolidation stage
ACC2 = "var(--acc2, #eb6834)"    # the mono-exon track
MUT = "var(--mut, #898781)"


def esc(t):
    return html.escape(t, quote=False)


def box(x, y, w, h, label, sub=None, *, rx=5, stroke="currentColor", fill="none",
        sw=1.4, lsize=13.5, ssize=11, dash=None, lw=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    s = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
         f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>']
    cx = x + w / 2
    lines = wrap(label, lw) if lw else [label]
    sub_lines = []
    if sub:
        # authored line breaks are deliberate; re-wrapping them is what broke
        # "molecule support · splice motif · local / site share"
        sub_lines = sub.split("\n")
    n = len(lines) + len(sub_lines)
    total = len(lines) * (lsize + 3) + len(sub_lines) * (ssize + 2.5)
    ty = y + h / 2 - total / 2 + lsize
    for ln in lines:
        s.append(f'<text x="{cx:.1f}" y="{ty:.1f}" text-anchor="middle" '
                 f'font-size="{lsize}" font-weight="600" fill="currentColor">{esc(ln)}</text>')
        ty += lsize + 3
    ty += 1
    for ln in sub_lines:
        s.append(f'<text x="{cx:.1f}" y="{ty:.1f}" text-anchor="middle" '
                 f'font-size="{ssize}" fill="{MUT}">{esc(ln)}</text>')
        ty += ssize + 2.5
    return "\n".join(s)


def wrap(text, width):
    if not width:
        return [text]
    out, cur = [], ""
    for w in text.split():
        if cur and len(cur) + 1 + len(w) > width:
            out.append(cur); cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        out.append(cur)
    return out


def arrow(x1, y1, x2, y2, label=None, *, color="currentColor", dash=None,
          side="right", sw=1.4, fs=10.5, mid=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    s = [f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
         f'stroke="{color}" stroke-width="{sw}"{d} marker-end="url(#ah)"/>']
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        if mid:
            mx, my = mid
        anchor = "start" if side == "right" else "end" if side == "left" else "middle"
        dx = 7 if side == "right" else -7 if side == "left" else 0
        dy = -5 if side == "middle" else 3.5
        for i, ln in enumerate(label.split("\n")):
            s.append(f'<text x="{mx+dx:.1f}" y="{my+dy+i*12:.1f}" text-anchor="{anchor}" '
                     f'font-size="{fs}" fill="{MUT}">{esc(ln)}</text>')
    return "\n".join(s)


DEFS = '''<defs>
  <marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6"
          markerHeight="6" orient="auto-start-reverse">
    <path d="M 0 1 L 9 5 L 0 9 z" fill="context-stroke"/>
  </marker>
</defs>'''


# =============================================================== detailed
def detailed():
    W, H = 1020, 962
    IX, IW = 24, 170          # evidence column
    CH = 236                  # the channel every evidence arrow routes down
    BY = 272                  # the mono-exonic bypass line
    SX, SW = 296, 404         # spine
    NX = 724                  # notes column
    s = [DEFS]

    rows = [
        (62,  68, "1 · Admit reads",
         "one shared set of alignment gates\nMAPQ · clipping · indels · NM", None),
        (164, 68, "2 · Build the junction catalogue",
         "every intron in every admitted read,\nsnapped to annotated sites", None),
        (266, 82, "3 · Curate novel junctions",
         "molecule support · splice motif · local site share\n"
         "RT-switch evidence · calibrated local FDR", None),
        (380, 82, "4 · Group identical chains, resolve 5′ suffixes",
         "a chain that is a 5′ suffix of another folds into it\n"
         "unless its own TSS evidence says otherwise", None),
        (494, 82, "5 · Cluster 3′ ends",
         "molecule-weighted 3′ peaks + clipped polyA tails\n"
         "annotation refines a nearby end, never invents a distant one", None),
        (608, 82, "6 · Second pass  (0.5.0)",
         "mono-exon models classed A–E against their gene;\n"
         "C/D/E replaced · end3_novel held to a six-tier gate", None),
        (722, 94, "7 · Consolidate  (0.6.x)",
         "3′: same chain, ends within 200 bp → merge\n"
         "5′: contiguous sub-chain of a 2×-better model → fold\n"
         "reads move with the models; nothing is discarded", ACC),
        (848, 76, "8 · Write the catalogue",
         "models.gtf · read_stat · abundance\ncell × isoform matrix · posthoc audit", None),
    ]
    for y, h, t, sub, acc in rows:
        st = acc or "currentColor"
        s.append(box(SX, y, SW, h, t, sub, stroke=st, sw=2.2 if acc else 1.4, lw=46))
    for i in range(len(rows) - 1):
        s.append(arrow(SX + SW / 2, rows[i][0] + rows[i][1],
                       SX + SW / 2, rows[i + 1][0] - 2))

    ins = [
        (62,  "aligned BAM", "minimap2 / pbmm2", 0),
        (138, "reference GTF", "annotation", 1),
        (214, "genome FASTA", "splice motifs", 2),
        (290, "barcode + UMI table", "molecule counts", 2),
        (366, "CAGE / TSS atlas", "optional · 5′ tier 2", 3),
        (442, "short-read coverage", "optional", 3),
    ]
    for y, t, sub, target in ins:
        s.append(box(IX, y, IW, 52, t, sub, stroke=MUT, sw=1.1, rx=4, lsize=11.5,
                     ssize=9.5, lw=24))
        ty = rows[target][0] + rows[target][1] / 2
        s.append(f'<path d="M {IX+IW} {y+26} H {CH} V {ty:.1f} H {SX-2}" '
                 f'fill="none" stroke="{MUT}" stroke-width="1.1" '
                 f'stroke-dasharray="3 3" marker-end="url(#ah)"/>')

    # the mono-exonic track, as a bypass rather than a box: single-exon reads
    # skip the chain stages entirely and rejoin at the second pass
    y0 = rows[0][0] + rows[0][1] - 14
    y1 = rows[5][0] + rows[5][1] / 2
    s.append(f'<path d="M {SX} {y0} H {BY} V {y1:.1f} H {SX-2}" fill="none" '
             f'stroke="{ACC2}" stroke-width="1.5" stroke-dasharray="6 4" '
             f'marker-end="url(#ah)"/>')
    for i, ln in enumerate(["unspliced reads", "bypass stages 2–5"]):
        s.append(f'<text x="{BY-8}" y="{548+i*13}" text-anchor="end" '
                 f'font-size="10.5" fill="{ACC2}">{esc(ln)}</text>')

    notes = [
        (62,  "A gate applied here is applied to every\nread that reaches any later stage."),
        (164, "A junction snaps to an annotated site\nonly within ±5 bp — it is never moved\nonto a plausible-looking one."),
        (266, "Curation decides which junctions may\nexist. Nothing downstream can bring\nback a junction rejected here."),
        (380, "Suffix resolution moves a 5′ end.\nIt never alters an intron chain."),
        (494, "Annotation may refine an end the reads\nalready support; it cannot manufacture\none the reads never reached."),
        (608, "The second pass may keep, remove or\nreplace a model. It never builds a\nstructure from scratch."),
        (722, "Both rules only move support between\nexisting models, so the read totals are\nidentical with and without them.\n--no-posthoc reproduces 0.5.1 exactly."),
        (848, "Every removal is one row in\n{prefix}.posthoc.tsv, naming the\nmodel it was folded into."),
    ]
    for y, txt in notes:
        for i, ln in enumerate(txt.split("\n")):
            s.append(f'<text x="{NX}" y="{y+16+i*13.5}" font-size="10.5" '
                     f'fill="{MUT}">{esc(ln)}</text>')

    for x, lab, anc in ((IX, "EVIDENCE", "start"), (SX, "PIPELINE", "start"),
                        (NX, "WHAT EACH STAGE MAY NOT DO", "start")):
        s.append(f'<text x="{x}" y="36" font-size="10" fill="{MUT}" '
                 f'text-anchor="{anc}" letter-spacing="0.08em">{lab}</text>')

    body = "\n".join(s)
    return (f'<svg viewBox="0 0 {W} {H}" role="img" xmlns="http://www.w3.org/2000/svg" '
            f'aria-label="The flightcollapse pipeline: six evidence inputs feed eight '
            f'stages, from read admission through junction curation, chain grouping, '
            f'3-prime end clustering, a second filtering pass and the 0.6.x '
            f'consolidation. Single-exon reads bypass the chain stages and rejoin at '
            f'the second pass. Each stage is annotated with what it is not permitted '
            f'to do.">\n{body}\n</svg>')


# =============================================================== distilled
def distilled():
    W, H = 1000, 286
    s = [DEFS]
    y, h = 86, 104
    gap = 26
    n = 5
    bw = (W - 48 - gap * (n - 1)) / n
    steps = [
        ("Admit reads", "one set of\nalignment gates", None),
        ("Curate junctions", "support · motif\nRT-switch · FDR", None),
        ("Resolve 5′ ends", "suffix chains fold\ninto the full chain", None),
        ("Place 3′ ends", "molecule-weighted\npeaks + polyA", None),
        ("Consolidate", "3′ ≤ 200 bp merge\n5′ sub-chain filter", ACC),
    ]
    xs = []
    for i, (t, sub, acc) in enumerate(steps):
        x = 24 + i * (bw + gap)
        xs.append(x)
        st = acc or "currentColor"
        s.append(box(x, y, bw, h, t, sub, stroke=st, sw=2.2 if acc else 1.4,
                     lsize=13, ssize=10.5, lw=20))
        if i:
            s.append(arrow(xs[i - 1] + bw, y + h / 2, x - 2, y + h / 2))

    s.append(f'<text x="24" y="56" font-size="11.5" fill="{MUT}">'
             f'aligned long reads, one cell barcode and UMI each</text>')
    s.append(f'<text x="{W-24}" y="56" text-anchor="end" font-size="11.5" fill="{MUT}">'
             f'reference-anchored transcript models</text>')

    # the one claim the paper needs: the two ends are treated differently
    ya = y + h + 34
    x4, x5 = xs[3], xs[4]
    s.append(f'<path d="M {x4+bw/2} {y+h} V {ya-14}" fill="none" stroke="{MUT}" '
             f'stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(f'<path d="M {x5+bw/2} {y+h} V {ya-14}" fill="none" stroke="{ACC}" '
             f'stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(f'<text x="{W/2}" y="{ya}" text-anchor="middle" font-size="12" '
             f'font-weight="600" fill="currentColor">'
             f'The two ends are not treated alike</text>')
    s.append(f'<text x="{W/2}" y="{ya+17}" text-anchor="middle" font-size="11" fill="{MUT}">'
             f'the 3′ end is a peak position, so it gets a tolerance; '
             f'the 5′ end is a containment relation, so it gets a structural filter</text>')

    body = "\n".join(s)
    return (f'<svg viewBox="0 0 {W} {H}" role="img" xmlns="http://www.w3.org/2000/svg" '
            f'aria-label="flightcollapse in five steps: admit reads, curate junctions, '
            f'resolve 5-prime ends, place 3-prime ends, consolidate. The final step '
            f'treats the two ends differently — a 200 bp tolerance at the 3-prime end '
            f'and a structural sub-chain filter at the 5-prime end.">\n{body}\n</svg>')


if __name__ == "__main__":
    for name, fn in (("flightcollapse_algorithm_distilled", distilled),):
        svg = fn()
        (OUT / f"{name}.svg").write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            + svg.replace("var(--acc, #2a78d6)", "#2a78d6")
                 .replace("var(--acc2, #eb6834)", "#eb6834")
                 .replace("var(--mut, #898781)", "#6b6a66")
                 .replace("currentColor", "#141413") + "\n")
        print(f"{name}: {len(svg):,} bytes")
