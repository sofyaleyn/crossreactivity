# CrossFlag — demo

**From an antibody sequence alone, flag the self-proteins it will bind by mistake — before it reaches patients.** This demo replays a validated Boltz-2 cofold screen on a drug that actually harmed people (SHR-1210 / camrelizumab) and shows the screen recovering its known off-targets — specifically, not as generic stickiness. **Triage, not certification.**

Everything here is regenerated from committed data with **no network and no GPU**:

```bash
make setup         # one-time: pip install -e . + matplotlib/pytest into crossflag-spike
make demo          # -> python -m crossflag.demo.build
# regenerates: demo/dashboard.html, demo/figures/*.png, demo/verdict_table.{md,json}
open demo/dashboard.html
make test-demo     # acceptance gate: pytest tests/test_demo.py
```

Prereq: the `crossflag-spike` conda env (Python 3.10; pinned deps in
[`data/results/scripts/requirements.lock.txt`](../data/results/scripts/requirements.lock.txt)).
`make setup` installs this package editable so `python -m crossflag.demo.build`
resolves. No network or GPU is used at build time.

Every number below traces to [`data/results/cofold_metrics.csv`](../data/results/cofold_metrics.csv) or a numbered experiment in [`findings.md`](../findings.md). Open `dashboard.html` for the visual version; its inline `#crossflag-data` block carries the exact rows behind each mark.

---

## The ~5-minute pitch

### 1. Problem (30s)
Antibody CDR-mediated off-target binding is a leading cause of biologic drug failure. **~33% of antibody leads bind ≥1 unintended self-protein** (Norden et al., *mAbs* 2024). It is driven by surface shape, not sequence — so it is invisible to sequence inspection. **SHR-1210 / camrelizumab**, an approved anti-PD-1, caused **reactive capillary hemangioma** in patients through a CDR off-target on **VEGFR2**. The question: could we have flagged its off-targets *before* the clinic, from sequence alone?

### 2. Result — the money chart (90s) → **Chart A** (`figures/chart_a_panel.png`)
Cofold the SHR-1210 Fv against a curated self-protein set. Two metrics per pair (semantics in [`data/results/README.md`](../data/results/README.md)):
- **`PAE_IF`** — interface tightness (Å; lower = tighter).
- **`epitope_reprod`** — does it hit the *same* epitope every time (Jaccard across 5 samples; higher = reproducible). **The primary discriminating metric.**

Calibrated against a **PD-1 cognate ceiling** and a **lysozyme non-binder floor** ([`findings.md`](../findings.md) Exp 4):

| antigen | PAE_IF (Å) | epitope_reprod | read |
|---|---|---|---|
| FZD5 | 5.69 | 0.660 | **confirmed** off-target |
| ULBP2 | 5.74 | 0.893 | **confirmed** off-target |
| PD-1 | 7.24 | 0.936 | ceiling (cognate) |
| VEGFR2 | 11.55 | 0.430 | **missed** (weakest) |
| lysozyme | 12.38 | 0.446 | floor (non-binder) |

**FZD5 and ULBP2 — both known SHR-1210 off-targets — dock as tightly and reproducibly as the cognate PD-1 target, recovered from sequence alone.**

### 3. Specificity — "not generic stickiness" (60s) → **Charts B + C**
- **Chart B — within-fold** (`figures/chart_b_family.png`, Exp 6): against **11 same-fold Frizzled-family decoys**, the true off-target FZD5 ranks **2/12** by interface PAE — **within-family AUROC 0.909**. The screen discriminates the real off-target from its own structural family, not just "binds Frizzled-shaped things."
- **Chart C — antibody-side control** (`figures/chart_c_control.png`, Exp 8): **pembrolizumab**, a different anti-PD-1 with no such off-targets, binds PD-1 but **collapses** on FZD5 (reprod 0.66→0.22) and ULBP2 (0.89→0.34). The confirmations are **SHR-1210-CDR-specific**, not an artifact of the antigens or the method.

### 4. What the interface looks like (60s) → **structure figures**
- **Interface PAE heatmap** (`figures/pae_heatmap_fzd5.png`): the antibody×antigen PAE block is a **tight, uniformly confident** region for SHR-1210×FZD5 (5.69 Å) and a **diffuse** one for pembrolizumab×FZD5 (18.55 Å). This *is* `PAE_IF`, visualized.
- **Epitope-reproducibility map** (`figures/epitope_map_fzd5.png`): SHR-1210 reproduces **one sharp epitope patch** across samples; pembrolizumab scatters. This *is* `epitope_reprod`, visualized.

### 5. Cost + the honest frontier (30s)
- **Cost** ([`findings.md`](../findings.md) §cost reframe): ~$0.20/cofold → **~$50–150/antibody** for a few-hundred-protein set, vs **$10–30k** for a wet-lab specificity screen.
- **Honesty.** **VEGFR2 — the weakest of SHR-1210's known off-targets — sits at the non-binder floor and is honestly flagged as missed.** Antibody–antigen is the hardest cofold class; this screen *prioritizes* wet-lab confirmation, it does not replace it. **Triage, not certification.** For each confirmed hit the tool names the confirmation assay to run first (see `verdict_table.md`).

---

## How this scales (productization — *not* built in this demo)

```
  antibody Fv ──▶ cofold vs curated self-protein reference set (few hundred)
                        │
                        ├─▶ calibrated panel (PD-1 ceiling / lysozyme floor)
                        ▼
              ranked verdict table ──▶ named confirmation assay per hit
                        │
   (scale-out, per plan.md — representative-set clustering, discovery +
    panel screens, saturation/Chao1 coverage) ──▶ 1000-candidate panel
```

The demo is the calibrated single-antibody core. The representative-set scaling path to a 1000-candidate panel is described in [`plan.md`](../plan.md) and is deliberately **out of scope here** — this demo proves the core signal is real, robust, and specific.

---

## Provenance & reproducibility

- **Data:** 24 committed cofold runs (Boltz-2.1, 5 samples each) — CIFs + PAE matrices under [`data/results/structures/`](../data/results/structures/), metrics in [`cofold_metrics.csv`](../data/results/cofold_metrics.csv), input JSONs in [`data/results/inputs/`](../data/results/inputs/).
- **Self-consistency:** `tests/test_demo.py` recomputes `cofold-fzd5`'s `PAE_IF`/`epitope_reprod` from the committed CIF+PAE and asserts it reproduces the CSV (±0.05 / ±0.03) — the structures and the extraction path agree.
- **Determinism:** two builds produce byte-identical artifacts; no committed number contradicts `cofold_metrics.csv`.
- **Verdict rule (frozen):** *confirmed* iff `PAE_IF < 9.8 Å` **and** `epitope_reprod ≥ 0.55` (mid-band of the PD-1/lysozyme calibrators); calibrators keep their ceiling/floor roles; a known off-target that fails the bar reads *missed*. See [`src/crossflag/demo/panel.py`](../src/crossflag/demo/panel.py).
