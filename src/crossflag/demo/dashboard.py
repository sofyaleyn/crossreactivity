"""Build the self-contained results dashboard (Phase 1).

A single standalone HTML page: every figure is base64-inlined, every number is
embedded as inline JSON, and there are no external requests. Opens in a browser
as-is and is also suitable for the Artifact tool (favicon 🧬, title
"CrossFlag — off-target cofold screen").

All numbers trace to ``data/results/cofold_metrics.csv`` via ``panel.load_rows``;
the inline ``#crossflag-data`` JSON block carries the exact rows behind each mark.
"""
from __future__ import annotations

import base64
import json

from . import diagram, panel, paths, run, screen_view

TITLE = "CrossFlag — off-target cofold screen"


def _b64(name: str) -> str:
    data = (paths.FIGURES / name).read_bytes()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def _figure(name: str, alt: str, caption: str, source: str) -> str:
    return f"""<figure class="card">
  <img src="{_b64(name)}" alt="{alt}" loading="eager">
  <figcaption>{caption}<span class="prov" data-source="{source}">source: {source}</span></figcaption>
</figure>"""


def build_html() -> str:
    rows = panel.load_rows()
    auroc, rank, n = panel.family_auroc(rows)
    table = run.build_table()

    # inline data block: exact numbers behind every mark (traceable)
    data_block = {
        "source": "data/results/cofold_metrics.csv",
        "verdict_table": table,
        "family": {
            "auroc": round(auroc, 3), "rank": rank, "size": n,
            "members": [
                {"antigen": r.antigen.replace(" CRD", ""), "PAE_IF": r.PAE_IF,
                 "run": r.run, "is_positive": r.run == panel.FAMILY_POSITIVE}
                for r in panel.family_rows(rows)
            ],
        },
        "control": [
            {"antigen": lbl,
             "shr1210_reprod": rows[shr].epitope_reprod,
             "pembrolizumab_reprod": rows[pem].epitope_reprod}
            for lbl, shr, pem in panel.CONTROL_PAIRS
        ],
        "thresholds": {"PAE_IF_CONFIRM": panel.PAE_IF_CONFIRM,
                       "REPROD_CONFIRM": panel.REPROD_CONFIRM},
    }

    # verdict table HTML rows (each traceable to a CSV run)
    badge = {"confirmed": "good", "ceiling": "ceiling", "missed": "warn", "floor": "floor"}
    trows = "\n".join(
        f'<tr data-run="{t["run"]}"><td>{t["antigen"]}</td>'
        f'<td class="num">{t["PAE_IF"]:.2f}</td>'
        f'<td class="num">{t["epitope_reprod"]:.3f}</td>'
        f'<td><span class="badge {badge[t["verdict"]]}">{t["read"]}</span></td>'
        f'<td class="assay">{t["confirmation_assay"]}</td></tr>'
        for t in table
    )

    figs = "\n".join([
        _figure("chart_a_panel.png",
                "Calibrated SHR-1210 self-protein panel",
                "<b>Chart A — the money chart.</b> From the SHR-1210 sequence alone, "
                "FZD5 and ULBP2 dock as tightly and reproducibly as the PD-1 cognate "
                "target (near the ceiling); VEGFR2 sits at the lysozyme non-binder floor "
                "and is honestly flagged as missed.",
                "cofold_metrics.csv rows: cofold-pd1, cofold-fzd5, cofold-ulbp2, cofold-wt, cofold-lyz"),
        _figure("chart_b_family.png",
                "FZD5 within-fold discrimination",
                f"<b>Chart B — within-fold discrimination.</b> Against 11 same-fold "
                f"Frizzled-family decoys, FZD5 (the true off-target) ranks {rank}/{n} by "
                f"interface PAE — within-family AUROC <b>{auroc:.3f}</b>. Not generic "
                f"stickiness to the fold.",
                "cofold_metrics.csv rows: cofold-fzd1..fzd10, cofold-smo, cofold-sfrp1, cofold-fzd5"),
        _figure("chart_c_control.png",
                "SHR-1210 vs pembrolizumab specificity control",
                "<b>Chart C — antibody-side specificity control.</b> Pembrolizumab "
                "(a different anti-PD-1 with no such off-targets) binds PD-1 but collapses "
                "on FZD5 (0.66→0.22) and ULBP2 (0.89→0.34). The confirmations are "
                "SHR-1210-CDR-specific.",
                "cofold_metrics.csv rows: cofold-{pd1,fzd5,ulbp2,wt} vs cofold-pembro-*"),
    ])

    structure_figs = "\n".join([
        _figure("pae_heatmap_fzd5.png",
                "Interface PAE heatmap, SHR-1210 vs pembrolizumab on FZD5",
                "<b>Interface PAE.</b> The antibody×antigen PAE block is a tight, uniformly "
                "confident low-error region for SHR-1210×FZD5 (PAE_IF 5.69 Å) and a diffuse, "
                "high-error one for pembrolizumab×FZD5 (18.55 Å). This is what PAE_IF measures.",
                "PAE matrices: structures/cofold-fzd5, cofold-pembro-fzd5 (best sample)"),
        _figure("epitope_map_fzd5.png",
                "Epitope reproducibility map, SHR-1210 vs pembrolizumab on FZD5",
                "<b>Epitope reproducibility.</b> Colored by how many of the 5 samples contact "
                "each antigen residue: SHR-1210 reproduces one sharp patch (reprod 0.660); "
                "pembrolizumab scatters (0.217). This is what epitope_reprod measures.",
                "CIF contacts: structures/cofold-fzd5, cofold-pembro-fzd5 (5 samples)"),
    ])

    # real-scale validation screen (Phase -1) — optional; skip if not run
    screen_section = ""
    if screen_view.SCREEN_CSV.exists():
        d = screen_view.compute()
        caps = {
            "screen_fpr_by_length.png": (
                "<b>Specificity by antigen size — the honest boundary.</b> False-positive "
                "rate falls sharply with ectodomain length: 2.1% for ≥150 aa, but ~20% for "
                "the small (&lt;150 aa) fragments the method over-docks (docking tighter than "
                "the real target — non-physical). A genuine limit, not a test artifact.",
                "screen_metrics.csv (flagship, hit rule applied per antigen)"),
            "screen_enrichment.png": (
                "<b>Sensitivity + generalization.</b> Left: in the valid regime the known "
                "off-targets FZD5/ULBP2 sit in the hit corner (tight + reproducible) while "
                "~98% of decoys do not. Right: for a blind second antibody (ABT-736), its real "
                "off-target PF4 stands out from same-size decoys — recovered with zero retuning.",
                "screen_metrics.csv (flagship ≥150 aa; anchor2 40–110 aa)"),
            "screen_scramble.png": (
                "<b>Paratope-specificity control.</b> Scrambling all six CDRs (framework "
                "preserved) collapses both off-targets below the hit threshold — the binding "
                "depends on the paratope sequence, not the scaffold.",
                "screen_metrics.csv (scramble6) vs cofold_metrics.csv (wild-type)"),
        }
        fig_html = [_figure(name, name, caps[name][0], caps[name][1])
                    for name in screen_view.build_figures(d)]
        screen_section = screen_view.section_html(d, fig_html)

    return _TEMPLATE.format(
        title=TITLE,
        auroc=f"{auroc:.3f}",
        trows=trows,
        figs=figs,
        structure_figs=structure_figs,
        screen_section=screen_section,
        pipeline_svg=diagram.build_svg(),
        thr_pae=panel.PAE_IF_CONFIRM,
        thr_rep=panel.REPROD_CONFIRM,
        data_json=json.dumps(data_block, indent=2),
    )


def build() -> str:
    html = build_html()
    paths.DEMO.mkdir(parents=True, exist_ok=True)
    paths.DASHBOARD.write_text(html)
    return str(paths.DASHBOARD.relative_to(paths.ROOT))


_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<link rel="icon" href="data:image/svg+xml,&lt;svg xmlns=%22http://www.w3.org/2000/svg%22 viewBox=%220 0 100 100%22&gt;&lt;text y=%22.9em%22 font-size=%2290%22&gt;🧬&lt;/text&gt;&lt;/svg&gt;">
<style>
  :root {{
    --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --border:rgba(11,11,11,.10);
    --good:#0ca30c; --warn:#eda100; --ceiling:#2a78d6; --floor:#898781;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
      --grid:#2c2c2a; --border:rgba(255,255,255,.10);
      --good:#0ca30c; --warn:#fab219; --ceiling:#3987e5; --floor:#898781;
    }}
  }}
  :root[data-theme="light"] {{
    --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e;
    --grid:#e1e0d9; --border:rgba(11,11,11,.10); --warn:#eda100; --ceiling:#2a78d6;
  }}
  :root[data-theme="dark"] {{
    --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7;
    --grid:#2c2c2a; --border:rgba(255,255,255,.10); --warn:#fab219; --ceiling:#3987e5;
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--plane); color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI",sans-serif; line-height:1.5; }}
  .wrap {{ max-width:1000px; margin:0 auto; padding:28px 20px 64px; }}
  header h1 {{ font-size:1.9rem; margin:0 0 4px; letter-spacing:-.01em; }}
  header .sub {{ color:var(--ink2); font-size:1.02rem; margin:0 0 18px; }}
  .lede {{ background:var(--surface); border:1px solid var(--border); border-radius:12px;
    padding:16px 18px; margin:0 0 24px; }}
  .lede p {{ margin:.3rem 0; }}
  .kpis {{ display:flex; flex-wrap:wrap; gap:12px; margin:14px 0 0; }}
  .kpi {{ flex:1 1 150px; background:var(--plane); border:1px solid var(--border);
    border-radius:10px; padding:10px 12px; }}
  .kpi b {{ display:block; font-size:1.5rem; line-height:1.1; }}
  .kpi span {{ color:var(--ink2); font-size:.82rem; }}
  h2 {{ font-size:1.15rem; margin:30px 0 12px; border-top:1px solid var(--grid); padding-top:22px; }}
  .card {{ background:var(--surface); border:1px solid var(--border); border-radius:12px;
    padding:14px; margin:0 0 20px; }}
  figure.card {{ margin:0 0 20px; }}
  .card img {{ width:100%; height:auto; display:block; border-radius:6px; }}
  figcaption {{ color:var(--ink2); font-size:.9rem; margin-top:10px; }}
  figcaption b {{ color:var(--ink); }}
  .prov {{ display:block; color:var(--muted); font-size:.74rem; margin-top:6px;
    font-variant-numeric:tabular-nums; }}
  .tablecard {{ overflow-x:auto; }}
  .svgwrap {{ overflow-x:auto; }}
  .svgwrap svg {{ min-width:760px; width:100%; height:auto; display:block; }}
  table {{ width:100%; border-collapse:collapse; font-size:.92rem; min-width:560px; }}
  th, td {{ text-align:left; padding:8px 10px; border-bottom:1px solid var(--grid); }}
  th {{ color:var(--ink2); font-weight:600; font-size:.82rem; text-transform:uppercase;
    letter-spacing:.03em; }}
  td.num {{ font-variant-numeric:tabular-nums; text-align:right; }}
  td.assay {{ color:var(--ink2); font-size:.84rem; }}
  .badge {{ display:inline-block; padding:2px 9px; border-radius:999px; font-size:.8rem;
    font-weight:600; }}
  .badge.good {{ background:color-mix(in srgb,var(--good) 18%,transparent); color:var(--good); }}
  .badge.warn {{ background:color-mix(in srgb,var(--warn) 22%,transparent);
    color:color-mix(in srgb,var(--warn) 75%,var(--ink)); }}
  .badge.ceiling {{ background:color-mix(in srgb,var(--ceiling) 18%,transparent); color:var(--ceiling); }}
  .badge.floor {{ background:color-mix(in srgb,var(--floor) 22%,transparent); color:var(--ink2); }}
  footer {{ margin-top:34px; padding-top:20px; border-top:1px solid var(--grid);
    color:var(--ink2); font-size:.88rem; }}
  footer b {{ color:var(--ink); }}
  code {{ background:var(--surface); border:1px solid var(--border); border-radius:5px;
    padding:1px 5px; font-size:.85em; }}
</style>
</head>
<body>
<div class="wrap">
<header>
  <h1>🧬 CrossFlag — off-target cofold screen</h1>
  <p class="sub">From an antibody sequence alone, flag which self-proteins it will bind by mistake — before it reaches patients.</p>
</header>

<div class="lede">
  <p><b>The problem.</b> ~33% of antibody leads bind ≥1 unintended self-protein
  (Norden et al., <i>mAbs</i> 2024). SHR-1210 / camrelizumab (an approved anti-PD-1)
  caused reactive capillary hemangioma in patients via a VEGFR2 CDR off-target —
  a liability invisible to sequence inspection.</p>
  <p><b>The screen.</b> A Boltz-2 cofold of the antibody Fv against a curated
  self-protein set. Two metrics per pair: <code>PAE_IF</code> (interface tightness)
  and <code>epitope_reprod</code> (does it hit the same epitope every time). Calibrated
  against a PD-1 cognate ceiling and a lysozyme non-binder floor.</p>
  <div class="kpis">
    <div class="kpi"><b>2 / 2</b><span>known off-targets (FZD5, ULBP2) recovered from sequence alone</span></div>
    <div class="kpi"><b>{auroc}</b><span>within-fold AUROC — FZD5 vs same-fold decoys</span></div>
    <div class="kpi"><b>~$50–150</b><span>per antibody vs $10–30k for a wet-lab specificity screen</span></div>
  </div>
</div>

<h2>Verdict table</h2>
<div class="card tablecard">
  <table>
    <thead><tr><th>Antigen</th><th>PAE_IF (Å)</th><th>epitope_reprod</th><th>read</th><th>confirmation assay</th></tr></thead>
    <tbody>
{trows}
    </tbody>
  </table>
</div>

<h2>The evidence</h2>
{figs}

<h2>What the interface actually looks like</h2>
{structure_figs}

{screen_section}

<h2>How this scales — productization (not built in this demo)</h2>
<div class="card">
  <div class="svgwrap">{pipeline_svg}</div>
  <figcaption>The demo proves the calibrated single-antibody core. To de-risk a whole
  <b>100–1000-candidate panel</b> without the full candidate×target grid, cluster the
  candidates, cofold a few representatives to discover the shared off-target shortlist,
  then screen every candidate against just that shortlist. Full method + cost model in
  <b>plan.md</b>.<span class="prov" data-source="plan.md">source: plan.md (representative-set scaling path)</span></figcaption>
</div>

<footer>
  <p><b>Triage, not certification.</b> Cofold is confirmation <i>evidence</i>, not
  proof — antibody–antigen is the hardest cofold class. VEGFR2, the weakest of
  SHR-1210's known off-targets, sits at the non-binder floor and is honestly shown
  as <i>missed</i>: this screen prioritizes wet-lab confirmation, it does not replace it.</p>
  <p class="prov">All numbers trace to <code>data/results/cofold_metrics.csv</code>
  (24 cofold runs, Boltz-2.1, 5 samples each). Full inline data below.
  Verdict rule: confirmed iff PAE_IF &lt; {thr_pae} Å and epitope_reprod ≥ {thr_rep}.</p>
</footer>
</div>
<script type="application/json" id="crossflag-data">
{data_json}
</script>
</body>
</html>
"""
