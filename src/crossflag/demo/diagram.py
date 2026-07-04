"""Representative-set scaling diagram (the plan.md productization path).

The demo proves the single-antibody core; this diagram is the one-slide
representation of how it scales to a 100–1000-candidate panel WITHOUT running the
full candidate x target grid (demo-plan.md guardrail #1: represent plan.md as a
diagram, not code). Every number is quoted from plan.md's cost model.

Emitted as a self-contained, theme-aware SVG:
- Boxes carry a solid light fill + dark ink via presentation attributes, so the
  diagram is readable even where <style> is stripped (e.g. GitHub's markdown
  <img> rendering).
- An embedded <style> block flips to dark palette under prefers-color-scheme:dark
  or an ancestor :root[data-theme="dark"] (standalone file and inside the
  dashboard), scoped to .pipeline-svg so it never leaks onto a host page.
"""
from __future__ import annotations

from . import paths

# --- geometry ---
W, H = 1160, 452
CARD_Y, CARD_H = 122, 206
CARDS = [30, 420, 810]        # left x of the three stage cards
CARD_W = 320

# --- palette (light defaults as presentation attrs; dark via <style>) ---
BOX = "#fbfbfa"
INK = "#17170b"
INK2 = "#52514e"
LINE = "#898781"
BLUE = "#2a78d6"
AQUA = "#1baf7a"
GREEN = "#0ca30c"
AMBER = "#c67f00"

_STYLE = """
<style>
.pipeline-svg .box  { fill:#fbfbfa; }
.pipeline-svg .ink  { fill:#17170b; }
.pipeline-svg .ink2 { fill:#52514e; }
.pipeline-svg .line { stroke:#898781; }
.pipeline-svg .bar  { fill:#f1f0ea; }
@media (prefers-color-scheme: dark) {
  .pipeline-svg .box  { fill:#1f1f1e; }
  .pipeline-svg .ink  { fill:#f5f5f0; }
  .pipeline-svg .ink2 { fill:#c3c2b7; }
  .pipeline-svg .line { stroke:#a7a69f; }
  .pipeline-svg .bar  { fill:#262624; }
}
:root[data-theme="dark"] .pipeline-svg .box  { fill:#1f1f1e; }
:root[data-theme="dark"] .pipeline-svg .ink  { fill:#f5f5f0; }
:root[data-theme="dark"] .pipeline-svg .ink2 { fill:#c3c2b7; }
:root[data-theme="dark"] .pipeline-svg .line { stroke:#a7a69f; }
:root[data-theme="dark"] .pipeline-svg .bar  { fill:#262624; }
:root[data-theme="light"] .pipeline-svg .box  { fill:#fbfbfa; }
:root[data-theme="light"] .pipeline-svg .ink  { fill:#17170b; }
:root[data-theme="light"] .pipeline-svg .ink2 { fill:#52514e; }
:root[data-theme="light"] .pipeline-svg .line { stroke:#898781; }
:root[data-theme="light"] .pipeline-svg .bar  { fill:#f1f0ea; }
</style>"""


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _lines(x, y, rows, size=11.5, dy=17, cls="ink2", weight="normal"):
    out = []
    for i, txt in enumerate(rows):
        out.append(
            f'<text class="{cls}" x="{x}" y="{y + i * dy}" font-size="{size}" '
            f'font-weight="{weight}">{_esc(txt)}</text>'
        )
    return "\n".join(out)


def _badge(cx, cy, n, color):
    return (
        f'<circle cx="{cx}" cy="{cy}" r="13" fill="{color}"/>'
        f'<text x="{cx}" y="{cy + 4.5}" font-size="15" font-weight="700" '
        f'fill="#ffffff" text-anchor="middle">{n}</text>'
    )


def _card(x, accent):
    return (
        f'<rect class="box" x="{x}" y="{CARD_Y}" width="{CARD_W}" height="{CARD_H}" '
        f'rx="14" stroke="{accent}" stroke-width="2"/>'
    )


def _arrow(x1, x2, y):
    return (
        f'<line class="line" x1="{x1}" y1="{y}" x2="{x2 - 6}" y2="{y}" stroke-width="2.5"/>'
        f'<path d="M {x2} {y} l -11 -6 v 12 z" fill="{LINE}"/>'
    )


def build_svg() -> str:
    a, b, c = CARDS
    ct = CARD_Y  # card top
    parts = [
        f'<svg class="pipeline-svg" xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {W} {H}" width="100%" role="img" '
        f'aria-label="Representative-set off-target screening pipeline">',
        _STYLE,

        # title (on a themed bar so it reads in every render path)
        f'<rect class="bar" x="30" y="10" width="{W-60}" height="36" rx="8"/>',
        f'<text class="ink" x="{W/2}" y="34" font-size="19" font-weight="700" '
        f'text-anchor="middle">How it scales to a 100–1000-candidate panel '
        f'— without the full grid</text>',

        # insight banner
        f'<rect class="bar" x="30" y="52" width="{W-60}" height="34" rx="8"/>',
        f'<text class="ink2" x="{W/2}" y="74" font-size="12.5" font-style="italic" '
        f'text-anchor="middle">Variants of one lead share ~all their off-targets — so '
        f'cofold a few representatives to find the shortlist, then screen every candidate '
        f'against just that shortlist.</text>',

        # cards
        _card(a, BLUE), _card(b, AQUA), _card(c, GREEN),

        # --- card 1: cluster ---
        _badge(a + 26, ct + 28, "1", BLUE),
        _lines(a + 50, ct + 33, ["Cluster candidates"], size=15, cls="ink", weight="700"),
        _lines(a + 18, ct + 62, [
            "Input: 100–1000 Fv variants",
            "of one antibody lead.",
            "",
            "Fingerprint each paratope —",
            "CDR charge · hydrophobicity ·",
            "aromatics · length (weight H3/L3).",
            "k-medoids / sphere-exclusion.",
        ]),
        f'<text class="ink" x="{a + 18}" y="{ct + CARD_H - 16}" font-size="13" '
        f'font-weight="700" fill="{BLUE}">→ K ≈ 10 representatives</text>',

        # --- card 2: discover shortlist ---
        _badge(b + 26, ct + 28, "2", AQUA),
        _lines(b + 50, ct + 33, ["Discover off-target shortlist"], size=15, cls="ink", weight="700"),
        _lines(b + 18, ct + 62, [
            "Cofold the K reps ×",
            "2,896 self-proteins (surfaceome),",
            "scored vs the calibrated panel.",
            "",
            "Union of hits = SHORTLIST",
            "(~20 proteins).",
        ]),
        _badge(b + 26, ct + CARD_H - 24, "3", AMBER),
        f'<text class="ink2" x="{b + 50}" y="{ct + CARD_H - 20}" font-size="10.5">'
        f'saturation / Chao1 → add reps</text>',
        f'<text class="ink2" x="{b + 50}" y="{ct + CARD_H - 8}" font-size="10.5">'
        f'until recall ≥ 95%</text>',

        # --- card 3: panel screen ---
        _badge(c + 26, ct + 28, "4", GREEN),
        _lines(c + 50, ct + 33, ["Screen all candidates"], size=15, cls="ink", weight="700"),
        _lines(c + 18, ct + 62, [
            "Cofold every candidate ×",
            "the shortlist only —",
            "not the full grid.",
            "",
            "→ per-candidate off-target profile.",
        ]),
        _badge(c + 26, ct + CARD_H - 24, "5", GREEN),
        f'<text class="ink2" x="{c + 50}" y="{ct + CARD_H - 20}" font-size="10.5">'
        f'rank + name confirmation assay</text>',
        f'<text class="ink2" x="{c + 50}" y="{ct + CARD_H - 8}" font-size="10.5">'
        f'→ prioritized wet-lab plate</text>',

        # arrows between cards
        _arrow(a + CARD_W, b, CARD_Y + CARD_H / 2),
        _arrow(b + CARD_W, c, CARD_Y + CARD_H / 2),

        # cost strip
        f'<rect class="bar" x="30" y="346" width="{W-60}" height="56" rx="10"/>',
        f'<text class="ink" x="{W/2}" y="369" font-size="13" font-weight="700" '
        f'text-anchor="middle">Cost per ~100-candidate panel: '
        f'<tspan fill="{GREEN}">~$4–7k</tspan>  vs full 100×2,896 grid '
        f'<tspan fill="{AMBER}">~$58k</tspan>  vs wet-lab MPA screen $10–30k</text>',
        f'<text class="ink2" x="{W/2}" y="390" font-size="11.5" text-anchor="middle">'
        f'Scales with K (saturation-bounded), not with panel size × reference-set size.</text>',

        # footer (on a themed bar for the same reason)
        f'<rect class="bar" x="30" y="414" width="{W-60}" height="28" rx="8"/>',
        f'<text class="ink2" x="{W/2}" y="432" font-size="10.5" text-anchor="middle" '
        f'font-style="italic">Productization path (plan.md), Boltz-2 cofold confirmer '
        f'validated in this demo · numbers quoted from plan.md · not built here.</text>',

        "</svg>",
    ]
    return "\n".join(parts)


def build() -> str:
    svg = build_svg()
    paths.FIGURES.mkdir(parents=True, exist_ok=True)
    out = paths.FIGURES / "pipeline.svg"
    out.write_text(svg + "\n")
    return "pipeline.svg"
