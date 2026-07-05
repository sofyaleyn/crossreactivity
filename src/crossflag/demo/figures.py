"""Generate every demo figure as a deterministic PNG from committed data.

Two families:
  Dashboard charts (Phase 1): chart_a_panel, chart_b_family, chart_c_control.
  Interface structures (Phase 2): pae_heatmap_fzd5, epitope_map_fzd5.

Colors come from the validated data-viz reference palette (status good/warning
for confirmed/missed; a one-hue blue ramp for the sequential heatmap/epitope
count). Figures render on a fixed light card surface so they read the same in a
light or dark page. Output is byte-deterministic (Agg, pinned Software tag, no
timestamps).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from . import panel, paths, scoring

# --- validated palette slots (references/palette.md) ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
GOOD = "#0ca30c"      # confirmed
WARNING = "#fab219"   # missed (amber)
SERIES_1 = "#2a78d6"  # blue  — SHR-1210
SERIES_8 = "#eb6834"  # orange — pembrolizumab
NEUTRAL = "#b8b7b0"   # reference / decoy marks

# one-hue blue sequential ramp, light -> dark (palette.md § Sequential hue)
_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
BLUE_SEQ = LinearSegmentedColormap.from_list("blue_seq", _BLUE)
BLUE_SEQ_R = LinearSegmentedColormap.from_list("blue_seq_r", _BLUE[::-1])

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

def _save(fig, name: str, tight: bool = True) -> str:
    paths.FIGURES.mkdir(parents=True, exist_ok=True)
    out = paths.FIGURES / name
    kw = dict(dpi=140, metadata={"Software": "crossflag-demo"})
    if tight:
        kw["bbox_inches"] = "tight"
    fig.savefig(out, **kw)
    plt.close(fig)
    return name


# ------------------------------------------------------------------ Chart A
def chart_a_panel(rows) -> str:
    """Calibrated panel: reprod (x) vs PAE_IF (y, inverted so up = tighter)."""
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ceil, floor = rows[panel.CEILING_RUN], rows[panel.FLOOR_RUN]

    # reference bands: PD-1 ceiling (top) and lysozyme floor (bottom)
    ax.axhline(ceil.PAE_IF, color=NEUTRAL, ls="--", lw=1)
    ax.axhline(floor.PAE_IF, color=NEUTRAL, ls="--", lw=1)
    ax.text(0.02, ceil.PAE_IF - 0.15, "PD-1 ceiling (cognate)", ha="left",
            va="bottom", color=INK2, fontsize=9, transform=ax.get_yaxis_transform())
    # floor label on the empty bottom-right, clear of the lysozyme marker/label (left)
    ax.text(0.98, floor.PAE_IF + 0.15, "lysozyme floor (non-binder)", ha="right",
            va="top", color=INK2, fontsize=9, transform=ax.get_yaxis_transform())

    style = {
        "confirmed": (GOOD, "o", 150),
        "missed": (WARNING, "D", 120),
        "ceiling": (SERIES_1, "*", 340),
        "floor": (NEUTRAL, "s", 120),
    }
    # explicit per-point label offsets (points) to avoid collisions
    label_off = {
        "cofold-pd1": (0, 13), "cofold-fzd5": (0, 13), "cofold-ulbp2": (0, 13),
        "cofold-wt": (11, 10), "cofold-lyz": (11, -16),
    }
    for run in panel.PANEL_RUNS:
        r = rows[run]
        v = panel.verdict(run, r)
        color, marker, size = style[v]
        ax.scatter(r.epitope_reprod, r.PAE_IF, s=size, c=color, marker=marker,
                   edgecolors="white", linewidths=1.2, zorder=3)
        ax.annotate(r.display, (r.epitope_reprod, r.PAE_IF),
                    textcoords="offset points", xytext=label_off[run],
                    ha="center", fontsize=10, color=INK, weight="bold")
    ax.annotate("weakest off-target —\nhonestly out of reach in silico",
                (rows["cofold-wt"].epitope_reprod, rows["cofold-wt"].PAE_IF),
                textcoords="offset points", xytext=(60, 22), fontsize=8.5,
                color=WARNING, ha="left",
                arrowprops=dict(arrowstyle="->", color=WARNING, lw=1))

    ax.invert_yaxis()  # up = tighter interface
    # headroom so the top point labels clear the title and the lysozyme label clears the edge
    _paes = [rows[r].PAE_IF for r in panel.PANEL_RUNS]
    ax.set_ylim(max(_paes) + 0.9, min(_paes) - 1.05)
    ax.set_xlabel("epitope reproducibility  (Jaccard across 5 samples →)")
    ax.set_ylabel("PAE_IF  (Å; ↑ = tighter interface)")
    ax.set_title("SHR-1210 self-protein panel — from sequence alone",
                 fontsize=13, weight="bold", color=INK)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    # legend by status (identity never color-alone)
    from matplotlib.lines import Line2D
    leg = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GOOD, markersize=10, label="confirmed off-target"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=WARNING, markersize=9, label="missed"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor=SERIES_1, markersize=14, label="cognate ceiling"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=NEUTRAL, markersize=9, label="non-binder floor"),
    ]
    ax.legend(handles=leg, loc="center", bbox_to_anchor=(0.60, 0.42),
              frameon=False, fontsize=8.5)
    return _save(fig, "chart_a_panel.png")


# ------------------------------------------------------------------ Chart B
def chart_b_family(rows) -> str:
    """FZD5 vs 11 same-fold decoys, ranked by PAE_IF; FZD5 highlighted."""
    fam = panel.family_rows(rows)  # already PAE_IF ascending
    auroc, rank, n = panel.family_auroc(rows)
    labels = [r.antigen.replace(" CRD", "") for r in fam]
    vals = [r.PAE_IF for r in fam]
    colors = [GOOD if r.run == panel.FAMILY_POSITIVE else NEUTRAL for r in fam]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = range(len(fam))
    ax.barh(list(y), vals, color=colors, height=0.66, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()  # best (lowest PAE) on top
    ax.set_xlabel("PAE_IF  (Å; lower = tighter interface)")
    ax.set_title(f"Within-fold discrimination — FZD5 vs Frizzled-family decoys\n"
                 f"FZD5 rank {rank}/{n} · within-family AUROC {auroc:.3f}",
                 fontsize=12, weight="bold", color=INK)
    for i, r in enumerate(fam):
        if r.run == panel.FAMILY_POSITIVE:
            ax.text(vals[i] + 0.2, i, "  FZD5 — true off-target", va="center",
                    color=GOOD, fontsize=9, weight="bold")
    ax.grid(True, axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    return _save(fig, "chart_b_family.png")


# ------------------------------------------------------------------ Chart C
def chart_c_control(rows) -> str:
    """SHR-1210 vs pembrolizumab epitope-reprod on the shared antigen set."""
    pairs = panel.CONTROL_PAIRS
    labels = [p[0] for p in pairs]
    shr = [rows[p[1]].epitope_reprod for p in pairs]
    pem = [rows[p[2]].epitope_reprod for p in pairs]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = range(len(pairs))
    w = 0.38
    ax.bar([i - w / 2 for i in x], shr, w, label="SHR-1210 (camrelizumab)",
           color=SERIES_1, zorder=3)
    ax.bar([i + w / 2 for i in x], pem, w, label="pembrolizumab (control)",
           color=SERIES_8, zorder=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("epitope reproducibility (Jaccard)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Antibody-side specificity control\n"
                 "both bind PD-1; pembrolizumab collapses on SHR-1210's off-targets",
                 fontsize=12, weight="bold", color=INK)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.grid(True, axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for i, (s, p) in enumerate(zip(shr, pem)):
        ax.text(i - w / 2, s + 0.02, f"{s:.2f}", ha="center", fontsize=8, color=INK2)
        ax.text(i + w / 2, p + 0.02, f"{p:.2f}", ha="center", fontsize=8, color=INK2)
    return _save(fig, "chart_c_control.png")


# ------------------------------------------------ Phase 2: interface figures
def pae_heatmap_fzd5(rows) -> str:
    """Antibody×antigen PAE block, SHR-1210×FZD5 vs pembro×FZD5, shared scale."""
    a, _hA, _lA, _vA = scoring.interface_pae_block("cofold-fzd5")
    b, _hB, _lB, _vB = scoring.interface_pae_block("cofold-pembro-fzd5")
    vmin = min(a.min(), b.min())
    vmax = max(a.max(), b.max())

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.8))
    fig.subplots_adjust(top=0.78, bottom=0.13, left=0.08, right=0.88, wspace=0.28)
    for ax, blk, title, if_val in (
        (axes[0], a, "SHR-1210 × FZD5", rows["cofold-fzd5"].PAE_IF),
        (axes[1], b, "pembrolizumab × FZD5", rows["cofold-pembro-fzd5"].PAE_IF),
    ):
        im = ax.imshow(blk, cmap=BLUE_SEQ_R, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_title(f"{title}\nPAE_IF = {if_val:.2f} Å", fontsize=11, weight="bold",
                     color=(GOOD if "SHR" in title else INK2), pad=8)
        ax.set_xlabel("antigen (FZD5) residue")
        ax.set_ylabel("antibody H+L residue")
        ax.tick_params(labelsize=8)
    cbar = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.03)
    cbar.set_label("PAE (Å) — darker = more confident (tighter)", fontsize=9)
    fig.suptitle("Interface PAE: a tight low-error block vs a diffuse one",
                 fontsize=13, weight="bold", color=INK, y=0.98)
    return _save(fig, "pae_heatmap_fzd5.png", tight=False)


def epitope_map_fzd5(rows) -> str:
    """Per-residue contact count (0..5 samples): SHR-1210×FZD5 vs pembro×FZD5."""
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 4.2), sharex=False)
    fig.subplots_adjust(top=0.80, bottom=0.12, left=0.06, right=0.90, hspace=0.75)
    for ax, run, title in (
        (axes[0], "cofold-fzd5", "SHR-1210 × FZD5"),
        (axes[1], "cofold-pembro-fzd5", "pembrolizumab × FZD5"),
    ):
        resnums, counts = scoring.epitope_profile(run)
        reprod = rows[run].epitope_reprod
        if len(resnums):
            lo, hi = min(resnums), max(resnums)
            strip = [[0] * (hi - lo + 1)]
            for rn, c in zip(resnums, counts):
                strip[0][rn - lo] = c
            ax.imshow(strip, cmap=BLUE_SEQ, vmin=0, vmax=5, aspect="auto",
                      extent=[lo, hi, 0, 1])
            ax.set_xlim(lo, hi)
        ax.set_yticks([])
        ax.set_title(f"{title}   epitope_reprod = {reprod:.3f}",
                     fontsize=10.5, weight="bold",
                     color=(GOOD if "SHR" in title else INK2), loc="left")
        ax.set_xlabel("antigen residue", fontsize=9)
    # shared colorbar proxy
    sm = plt.cm.ScalarMappable(cmap=BLUE_SEQ, norm=plt.Normalize(0, 5))
    cbar = fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.02, ticks=[0, 1, 2, 3, 4, 5])
    cbar.set_label("# of 5 samples contacting", fontsize=9)
    fig.suptitle("Epitope reproducibility: a sharp reproduced patch vs scatter",
                 fontsize=13, weight="bold", color=INK, y=0.97)
    return _save(fig, "epitope_map_fzd5.png", tight=False)


ALL_FIGURES = [
    chart_a_panel, chart_b_family, chart_c_control,
    pae_heatmap_fzd5, epitope_map_fzd5,
]


def build_all() -> list[str]:
    rows = panel.load_rows()
    return [fn(rows) for fn in ALL_FIGURES]
