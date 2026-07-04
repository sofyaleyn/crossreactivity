"""Representative-set scaling diagram (the plan.md productization path).

The demo proves the single-antibody core; this is the one-slide picture of how it
scales to a 100–1000-candidate panel WITHOUT running the full candidate x target
grid (demo-plan.md guardrail #1: represent plan.md as a diagram, not code).

The visual *is* the argument: the naive approach fills a candidate x target
matrix (100 x 2,896 = 289,600 cofolds). Clustering collapses the rows to K≈10
representatives; a discovery screen collapses the columns to a ~20-protein
shortlist. You only ever compute a thin **cross** — a few representative rows
across all columns (discovery) plus all rows across a few shortlist columns
(panel) — a small fraction of the square. Every number is quoted from plan.md.

Emitted as a self-contained, theme-aware SVG: solid light-fill panels + dark ink
via presentation attributes (readable even where <style> is stripped, e.g.
GitHub's markdown <img>), with a dark palette layered on via an embedded <style>
scoped to .pipeline-svg.
"""
from __future__ import annotations

from . import paths

W, H = 1180, 636

# accent colors read on both light and dark surfaces
BLUE = "#2a78d6"    # representatives / discovery
GREEN = "#0ca30c"   # shortlist / panel
TEAL = "#0f766e"    # overlap (rep × shortlist)
AMBER = "#c67f00"

# grid geometry (candidates = rows, targets = columns)
GX, GY, GW, GH = 430, 176, 470, 300
REP_FRAC = 0.13          # 10 reps of ~100 candidates → top band height
PANEL_X = GX + 0.60 * GW  # where the ~20-protein shortlist columns sit
PANEL_W = 30

_STYLE = """
<style>
.pipeline-svg .surface { fill:#fcfcfb; }
.pipeline-svg .box  { fill:#fbfbfa; }
.pipeline-svg .grid { fill:#f4f3ee; }
.pipeline-svg .ink  { fill:#17170b; }
.pipeline-svg .ink2 { fill:#52514e; }
.pipeline-svg .bar  { fill:#f1f0ea; }
.pipeline-svg .cell { stroke:#d7d6cd; }
.pipeline-svg .faint{ fill:#c9c8c0; }
.pipeline-svg .edge { stroke:#c3c2b7; }
@media (prefers-color-scheme: dark) {
  .pipeline-svg .surface { fill:#1a1a19; }
  .pipeline-svg .box  { fill:#1f1f1e; }
  .pipeline-svg .grid { fill:#141413; }
  .pipeline-svg .ink  { fill:#f5f5f0; }
  .pipeline-svg .ink2 { fill:#c3c2b7; }
  .pipeline-svg .bar  { fill:#262624; }
  .pipeline-svg .cell { stroke:#33322f; }
  .pipeline-svg .faint{ fill:#4a4945; }
  .pipeline-svg .edge { stroke:#3a3a37; }
}
:root[data-theme="dark"] .pipeline-svg .surface { fill:#1a1a19; }
:root[data-theme="dark"] .pipeline-svg .box  { fill:#1f1f1e; }
:root[data-theme="dark"] .pipeline-svg .grid { fill:#141413; }
:root[data-theme="dark"] .pipeline-svg .ink  { fill:#f5f5f0; }
:root[data-theme="dark"] .pipeline-svg .ink2 { fill:#c3c2b7; }
:root[data-theme="dark"] .pipeline-svg .bar  { fill:#262624; }
:root[data-theme="dark"] .pipeline-svg .cell { stroke:#33322f; }
:root[data-theme="dark"] .pipeline-svg .faint{ fill:#4a4945; }
:root[data-theme="dark"] .pipeline-svg .edge { stroke:#3a3a37; }
:root[data-theme="light"] .pipeline-svg .surface { fill:#fcfcfb; }
:root[data-theme="light"] .pipeline-svg .box  { fill:#fbfbfa; }
:root[data-theme="light"] .pipeline-svg .grid { fill:#f4f3ee; }
:root[data-theme="light"] .pipeline-svg .ink  { fill:#17170b; }
:root[data-theme="light"] .pipeline-svg .ink2 { fill:#52514e; }
:root[data-theme="light"] .pipeline-svg .bar  { fill:#f1f0ea; }
:root[data-theme="light"] .pipeline-svg .cell { stroke:#d7d6cd; }
:root[data-theme="light"] .pipeline-svg .faint{ fill:#c9c8c0; }
:root[data-theme="light"] .pipeline-svg .edge { stroke:#c3c2b7; }
</style>"""


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _t(x, y, s, size=12, cls="ink2", weight="normal", anchor="start", style=""):
    st = f' font-style="{style}"' if style else ""
    return (f'<text class="{cls}" x="{x}" y="{y}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}"{st}>{_esc(s)}</text>')


def _badge(cx, cy, n, color, r=13):
    return (f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}"/>'
            f'<text x="{cx}" y="{cy + r*0.34}" font-size="{r*1.1:.0f}" font-weight="700" '
            f'fill="#ffffff" text-anchor="middle">{n}</text>')


def _cluster_inset():
    """① candidate set narrows: a scatter of variants → K medoids (blue)."""
    px, py, pw, ph = 40, 150, 250, 300
    parts = [
        f'<rect class="box" x="{px}" y="{py}" width="{pw}" height="{ph}" rx="14" '
        f'stroke="{BLUE}" stroke-width="2"/>',
        _badge(px + 24, py + 26, "1", BLUE),
        _t(px + 46, py + 31, "Cluster candidates", 14.5, "ink", "700"),
        _t(px + 20, py + 58, "100–1000 Fv variants of one lead", 11, "ink2"),
    ]
    # five clusters: (medoid x,y) + satellite dot offsets
    clusters = [
        (px + 55, py + 108, [(-16, -8), (-9, 12), (13, -13), (10, 10)]),
        (px + 150, py + 96, [(-14, -10), (-4, 14), (15, 4)]),
        (px + 200, py + 150, [(-15, -6), (-2, 15), (14, 9), (6, -14)]),
        (px + 78, py + 178, [(-13, 10), (12, -11), (16, 8)]),
        (px + 150, py + 200, [(-15, 6), (-3, 16), (14, -8), (13, 12)]),
    ]
    for mx, my, sats in clusters:
        # faint dashed hull suggestion
        parts.append(f'<circle class="faint" cx="{mx}" cy="{my}" r="26" '
                     f'fill="none" stroke-dasharray="3 3" stroke-width="1" '
                     f'stroke="{BLUE}" opacity="0.35"/>')
        for dx, dy in sats:
            parts.append(f'<circle class="faint" cx="{mx+dx}" cy="{my+dy}" r="3.4"/>')
        parts.append(f'<circle cx="{mx}" cy="{my}" r="6.5" fill="{BLUE}"/>')
    parts.append(_t(px + pw/2, py + ph - 40, "cluster by paratope biochemistry",
                    10.5, "ink2", anchor="middle"))
    parts.append(_t(px + pw/2, py + ph - 20,
                    "→ K ≈ 10 representatives (medoids)", 12.5, "ink", "700",
                    anchor="middle"))
    # arrow from inset into the grid's discovery band
    ay = GY + GH * REP_FRAC / 2
    parts.append(f'<path class="edge" d="M {px+pw} {ay} H {GX-8}" stroke-width="2.5" '
                 f'fill="none" marker-end="url(#arrow)"/>')
    parts.append(_t((px+pw+GX)/2, ay - 8, "reps", 10, "ink2", anchor="middle"))
    return "\n".join(parts)


def _grid():
    """The candidate × target matrix with the computed 'cross' highlighted."""
    band_h = GH * REP_FRAC
    parts = [
        # grid background + cell texture
        f'<rect class="grid" x="{GX}" y="{GY}" width="{GW}" height="{GH}" rx="4"/>',
        f'<rect x="{GX}" y="{GY}" width="{GW}" height="{GH}" rx="4" '
        f'fill="url(#cells)"/>',
        # discovery band (rep rows × all columns)
        f'<rect x="{GX}" y="{GY}" width="{GW}" height="{band_h:.0f}" fill="{BLUE}" '
        f'opacity="0.85"/>',
        # panel band (all rows × shortlist columns)
        f'<rect x="{PANEL_X:.0f}" y="{GY}" width="{PANEL_W}" height="{GH}" '
        f'fill="{GREEN}" opacity="0.85"/>',
        # overlap (reps × shortlist)
        f'<rect x="{PANEL_X:.0f}" y="{GY}" width="{PANEL_W}" height="{band_h:.0f}" '
        f'fill="{TEAL}"/>',
        # grid border
        f'<rect class="edge" x="{GX}" y="{GY}" width="{GW}" height="{GH}" rx="4" '
        f'fill="none" stroke-width="1.5"/>',
        # "never computed" label in the empty quadrant
        _t(GX + GW*0.28, GY + GH*0.62, "the grid you never fill", 12.5, "ink2",
           "600", "middle", style="italic"),
        _t(GX + GW*0.28, GY + GH*0.62 + 18, "289,600 cofolds · ≈ $58k", 11, "ink2",
           anchor="middle"),

        # top axis: targets 2,896 → shortlist
        _t(GX, GY - 44, "SELF-PROTEIN TARGETS", 11, "ink", "700"),
        _t(GX, GY - 28, "2,896 surfaceome proteins", 10.5, "ink2"),
        f'<line class="edge" x1="{GX}" y1="{GY-10}" x2="{GX+GW}" y2="{GY-10}" '
        f'stroke-width="1.5"/>',
        # shortlist bracket over the panel columns
        f'<rect x="{PANEL_X:.0f}" y="{GY-16}" width="{PANEL_W}" height="8" '
        f'fill="{GREEN}"/>',
        _t(PANEL_X + PANEL_W/2, GY - 22, "shortlist ~20", 10.5, "ink", "700",
           "middle"),

        # left axis: candidates
        f'<text class="ink" transform="translate({GX-118},{GY+GH/2}) rotate(-90)" '
        f'font-size="11" font-weight="700" text-anchor="middle">CANDIDATES · 100</text>',
    ]
    return "\n".join(parts)


def _annotations():
    band_h = GH * REP_FRAC
    dy = GY + band_h / 2
    parts = [
        # ② discovery label (right of band)
        _badge(GX + GW + 24, dy, "2", BLUE),
        _t(GX + GW + 42, dy - 4, "Discovery screen", 13, "ink", "700"),
        _t(GX + GW + 42, dy + 12, "K reps × all 2,896 → hits", 10.5, "ink2"),
        # ③ panel label (below panel band)
        _badge(PANEL_X + PANEL_W/2, GY + GH + 24, "3", GREEN),
        _t(PANEL_X + PANEL_W/2, GY + GH + 51, "Panel screen", 13, "ink", "700",
           "middle"),
        _t(PANEL_X + PANEL_W/2, GY + GH + 67, "all candidates × shortlist only",
           10.5, "ink2", "middle"),
        # causal arrow: discovery hits DEFINE the shortlist columns
        f'<path class="edge" d="M {GX+GW+20} {dy+16} '
        f'C {GX+GW+30} {GY-40}, {PANEL_X+90} {GY-46}, {PANEL_X+PANEL_W+2} {GY-12}" '
        f'fill="none" stroke-width="2" stroke-dasharray="4 3" '
        f'marker-end="url(#arrow)"/>',
        _t(GX + GW - 4, GY - 54, "hits define the shortlist", 10, "ink2", "600",
           "end", style="italic"),
    ]
    return "\n".join(parts)


def _cost_bar():
    y = 560
    parts = [
        f'<rect class="bar" x="40" y="{y}" width="{W-80}" height="72" rx="12"/>',
        # legend swatches
        f'<rect x="66" y="{y+20}" width="16" height="16" fill="{BLUE}" opacity="0.85"/>',
        f'<rect x="66" y="{y+20}" width="8" height="16" fill="{GREEN}" opacity="0.85"/>',
        _t(92, y+32, "What you actually cofold — the cross:", 12, "ink", "700"),
        _t(92, y+50, "~31,000 cofolds  ≈  $4–7k per 100-candidate panel", 12, "ink"),
        f'<rect class="faint" x="470" y="{y+20}" width="16" height="16"/>',
        _t(496, y+32, "What you skip — the empty grid:", 12, "ink", "700"),
        _t(496, y+50, "258,600 cofolds you never run  ≈  the $52k you save", 12, "ink2"),
        _t(W-66, y+32, "Full 100 × 2,896 grid = 289,600 ≈ $58k", 11.5, "ink2", "600", "end"),
        _t(W-66, y+50, "wet-lab MPA screen $10–30k · scales with K, not panel×set",
           11, "ink2", anchor="end"),
    ]
    return "\n".join(parts)


def build_svg() -> str:
    parts = [
        f'<svg class="pipeline-svg" xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {W} {H}" width="100%" role="img" '
        f'aria-label="Representative-set screening never fills the candidate by target '
        f'grid: cluster candidates to K representatives, discovery-screen to a ~20 protein '
        f'shortlist, then screen all candidates against the shortlist — a thin cross of the '
        f'full grid.">',
        _STYLE,
        # defs: arrow marker + faint cell pattern
        f'<defs>'
        f'<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        f'markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M 0 0 L 10 5 L 0 10 z" class="faint"/></marker>'
        f'<pattern id="cells" width="13" height="13" patternUnits="userSpaceOnUse">'
        f'<path class="cell" d="M 13 0 L 0 0 0 13" fill="none" stroke-width="0.6"/>'
        f'</pattern>'
        f'</defs>',

        # opaque themed surface so every label has a readable backing everywhere
        f'<rect class="surface" x="0" y="0" width="{W}" height="{H}" rx="10"/>',

        # title + subtitle
        f'<rect class="bar" x="30" y="12" width="{W-60}" height="34" rx="8"/>',
        _t(W/2, 35, "How it scales: you never fill the candidate × target grid",
           19, "ink", "700", "middle"),
        f'<rect class="bar" x="30" y="52" width="{W-60}" height="30" rx="8"/>',
        _t(W/2, 72, "Variants of one lead share ~all their off-targets — so cluster the "
           "candidates, discover a shortlist, and cofold only a thin cross of the grid.",
           12.5, "ink2", "normal", "middle", style="italic"),

        _cluster_inset(),
        _grid(),
        _annotations(),
        _cost_bar(),
        "</svg>",
    ]
    return "\n".join(parts)


def build() -> str:
    svg = build_svg()
    paths.FIGURES.mkdir(parents=True, exist_ok=True)
    out = paths.FIGURES / "pipeline.svg"
    out.write_text(svg + "\n")
    return "pipeline.svg"
