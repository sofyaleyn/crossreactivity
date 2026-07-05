# CrossFlag — demo

**From an antibody sequence alone, flag the self-proteins its CDRs will bind by mistake — before it reaches patients.** CrossFlag cofolds a therapeutic antibody's Fv (VH+VL) against a curated human surfaceome with Boltz-2 and scores each interface to prioritize likely off-targets *ahead of* a wet-lab specificity screen. This demo replays a **real-scale validation screen** on a drug that actually harmed people — SHR-1210 / camrelizumab (anti-PD-1) — and shows the frozen pipeline recovering its known off-targets out of ~1,200 self-proteins, specifically, not as generic stickiness. **Triage, not certification.**

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

Every number below traces to [`data/results/screen/validation_summary.md`](../data/results/screen/validation_summary.md) (the screen), [`data/results/screen/screen_metrics.csv`](../data/results/screen/screen_metrics.csv) (per-protein scores), [`data/results/cofold_metrics.csv`](../data/results/cofold_metrics.csv) (the calibration panel), or a numbered experiment in [`findings.md`](../findings.md). Open `dashboard.html` for the visual version; its inline `#crossflag-data` block carries the exact rows behind each mark.

---

## The ~5-minute pitch

### 1. Problem (30s)
Antibody CDR-mediated off-target binding is a leading cause of biologic drug failure. **~33% of antibody leads bind ≥1 unintended self-protein** (Norden et al., *mAbs* 2024). It is driven by surface shape, not sequence — so it is invisible to sequence inspection. **SHR-1210 / camrelizumab**, an approved anti-PD-1, caused **reactive capillary hemangioma** in patients through a CDR off-target on **VEGFR2** (plus FZD5, ULBP2). The question: could we have flagged its off-targets *before* the clinic, from sequence alone?

### 2. The result at scale — lead here (120s) → **screen figures**
We ran the frozen pipeline as a real screen, not a toy: **SHR-1210's Fv cofolded against 1,198 curated self-protein ectodomains** (1,198/1,200; 2 dropped on an API rate-limit). Counting the second anchor (401) and scramble controls, **1,603 cofolds total**, 5 samples each, for **$342.70 of a $500 budget** (~$0.20/cofold). Nothing was tuned to these results — the hit rule was frozen on the calibration panel *before* the screen ran.

**What it's screened against — the denominator matters.** The reference set is not a hand-picked decoy list. It's the human cell-surface proteome — the **SURFY surfaceome** (Bausch-Fluck et al., *PNAS* 2018): **2,886 proteins**, every self-protein an antibody's CDRs could physically reach *in vivo* (the compartment where off-target binding actually happens), **948** of them mass-spec-confirmed (CSPA). We add curated autoimmune/mimicry self-proteins → a **2,896-protein** set (99.8% with UniProt sequences), then cofold the ~1,200-protein slice above as each protein's extracellular ectodomain. Recovering the known off-targets near the top of *that* haystack is the result — against five cherry-picked decoys it would prove nothing. ([`reference-set.md`](../docs/reference-set.md))

**Specificity holds where it matters.** In the valid regime (antigen ectodomain ≥150 aa) the false-positive rate is **2.1% (14/678)**, and FPR falls monotonically with antigen size (`figures/screen_fpr_by_length.png`):

| ectodomain length | FPR |
|---|--:|
| 0–80 aa | 20.6% |
| 80–150 aa | 19.1% |
| 150–300 aa | 4.5% |
| 300–600 aa | 1.2% |
| 600–1000 aa | 0.0% |
| **≥150 aa (valid regime)** | **2.1% (14/678)** |

**Both recoverable known off-targets are enriched out of ~1,200 proteins** (`figures/screen_enrichment.png`): among same-size decoys, **FZD5** (PAE_IF 5.69 / reprod 0.660) ranks above 84% on reproducibility and 99% on PAE; **ULBP2** (5.74 / 0.893) is more reproducible than 242/243 proteins in its size class. Both flagged. (VEGFR2 is missed — owned in §4.)

**It generalizes to a second, blind antibody — zero retuning.** **ABT-736** (anti-Aβ), whose known off-target **PF4** drives a HIT-like toxicity (Steinmetz et al., *mAbs* 2021, PMID 33596779; VH/VL transcribed from patent US2009/0175847A1; PF4 = UniProt P02776), was run through the same frozen thresholds. **PF4 is flagged** — top ~9% by reproducibility, **#1 by PAE among 76 same-size decoys**. Valid-regime FPR is antibody-specific, not generic stickiness: **SHR-1210 2.1%, ABT-736 0.4%.**

**The paratope is doing the binding — negative control** (`figures/screen_scramble.png`): scrambling all six CDRs (framework preserved) collapses both off-targets to non-hits.

| antigen | wild-type | all-6-CDR scramble |
|---|---|---|
| FZD5 | hit (PAE 5.69 / reprod 0.66) | **non-hit** (9.94 / 0.26) |
| ULBP2 | hit (5.74 / 0.89) | **non-hit** (11.19 / 0.82) |

FZD5 required a *full* 6-CDR scramble to break — a partial H3/L3-only scramble did **not** — so its cross-reactivity draws on more of the paratope than the two dominant loops alone.

### 3. The calibrated core the screen rests on (60s) → **Charts A + B + C**
The screen is only trustworthy because the metric and threshold were pinned on a single antibody first ([`findings.md`](../findings.md) Exp 4). Two metrics per pair (semantics in [`data/results/README.md`](../data/results/README.md)):
- **`PAE_IF`** — interface tightness (Å; lower = tighter).
- **`epitope_reprod`** — does it hit the *same* epitope every time (Jaccard across 5 samples; higher = reproducible).
- **`ipTM` is dropped** — it over-docks (lysozyme 0.87 ≈ everything), so it is *not* in the hit rule.

Calibrated against a **PD-1 cognate ceiling** and a **lysozyme non-binder floor** (`figures/chart_a_panel.png`):

| antigen | PAE_IF (Å) | epitope_reprod | read |
|---|---|---|---|
| FZD5 | 5.69 | 0.660 | **confirmed** off-target |
| ULBP2 | 5.74 | 0.893 | **confirmed** off-target |
| PD-1 | 7.24 | 0.936 | ceiling (cognate) |
| VEGFR2 | 11.55 | 0.430 | **missed** (weakest) |
| lysozyme | 12.38 | 0.446 | floor (non-binder) |

- **Within-fold discrimination** (`figures/chart_b_family.png`, Exp 6): against **11 same-fold Frizzled decoys**, FZD5 ranks **2/12** by interface PAE — **within-family AUROC 0.909**. The screen tells the real off-target apart from its own structural family, not just "binds Frizzled-shaped things."
- **Antibody-side control** (`figures/chart_c_control.png`, Exp 8): **pembrolizumab**, a different anti-PD-1 without these off-targets, binds PD-1 but **collapses** on FZD5 (reprod 0.66→0.22) and ULBP2 (0.89→0.34) — the confirmations are SHR-1210-CDR-specific.
- **What the interface looks like:** the PAE heatmap (`figures/pae_heatmap_fzd5.png`) is a **tight, uniformly confident** block for SHR-1210×FZD5 (5.69 Å) vs a **diffuse** one for pembrolizumab×FZD5 (18.55 Å) — this *is* `PAE_IF`. The epitope map (`figures/epitope_map_fzd5.png`) shows SHR-1210 reproducing **one sharp epitope patch** while pembrolizumab scatters — this *is* `epitope_reprod`.

### 4. The honest frontier (60s)
- **We miss the off-target that caused the harm.** VEGFR2 (PAE_IF 11.55 / reprod 0.430) sits at the non-binder floor and **fails the bar — missed.** It is SHR-1210's weakest, lowest-affinity off-target ("aberrant… low affinity" per Finlay), multi-domain, with an unmapped epitope; templating it to its real 3V2A structure made the pose *worse*, not better. It is also the off-target that actually **caused the clinical capillary hemangioma** — the one case we most wanted to catch is the one we don't. We do not bury this: the screen prioritizes wet-lab confirmation, it does not stand in for it.
- **A real, bounded flaw: small antigens over-dock.** For ectodomains <150 aa the FPR is ~20%, with docking confidence *sometimes higher than the true target* (physically implausible; of 171 decoys docking tighter than a real off-target, median length 41 aa, only 3 are ≥150 aa). This is a genuine method flaw, not a test artifact — the short ectodomains are correctly extracted real domains — and it touches **~43% of the surfaceome** (1,257/2,896 proteins have <150 aa ectodomains). **Found *and* handled:** operate in the valid regime (≥150 aa) and route small-ectodomain hits through a **size-aware confidence gate** as LOW-CONFIDENCE rather than trusting or discarding them.
- **The false positives are a *characterized* failure mode, not loose ends.** The screen's top valid-regime non-target hits (SMO; and the Ig-fold receptors CD19, IL22RA1, five KIRs — KIR2DL4 / KIR2DS1 / KIR3DL3 / KIR2DL3 / KIR2DL5B) concentrate **entirely** in the fold families of the confirmed off-targets: SMO shares FZD5's Frizzled cysteine-rich domain; the rest are Ig-superfamily β-sandwiches like ULBP2 and like the PD-1 target itself. Decisively, the orthogonal *experimental* screen that found SHR-1210's real off-targets (Finlay et al., Retrogenix ~4,975 receptors, *mAbs* 2019, PMID 30541416) tested every one and scored them all **negative**. So these aren't candidate discoveries — they're the diagnostic signature of an **Ig-fold / CRD affinity bias**, the method's main non-size failure mode, now named and bounded rather than hand-waved. Full read: [`data/results/screen/novel_candidates_assessment.md`](../data/results/screen/novel_candidates_assessment.md).

### 5. Cost + reproducibility (30s)
- **Cost:** ~$0.20/cofold → the full 1,198-protein screen cost **$342.70**; a few-hundred-protein set runs **~$50–150/antibody**, vs **$10–30k** for a wet-lab specificity screen. Cofolding the whole curated set directly is cheaper than building a cheap pre-filter that works ([`findings.md`](../findings.md) §cost reframe, Exp 9).
- **Reproducibility:** the entire pitch above regenerates offline from committed data via `make demo` — see the provenance section below.

**Bottom line (from [`validation_summary.md`](../data/results/screen/validation_summary.md)): conditionally validated.** At screen scale the cofold screen reliably prioritizes off-targets for medium/large antigens (≥150 aa: 2.1% FPR, both known off-targets enriched), is paratope-specific, and generalizes across two independent antibodies. It has one real, bounded flaw — unreliable scoring for small (<150 aa) ectodomains — handled by a size-aware gate. **A prioritization tool for wet-lab specificity screens, not a replacement.**

---

## Provenance & reproducibility

- **Screen data:** the flagship run and controls live under [`data/results/screen/`](../data/results/screen/) — `screen_metrics.csv` (per-protein scores), `rig.json` (frozen thresholds), `spend_ledger.json` (budget), `checks.py` (the validation analysis), job configs under `jobs/`. Ectodomains: [`data/reference/ectodomains.csv`](../data/reference/ectodomains.csv). Complete run: **1,198/1,200** flagship + **401/401** anchor-2 cofolds scored, plus 20 pilot + 4 scramble folds (1,603 total, 5 samples each).
- **Calibration data:** the 24 committed calibration cofolds (Boltz-2.1, 5 samples each) — CIFs + PAE matrices under [`data/results/structures/`](../data/results/structures/), metrics in [`cofold_metrics.csv`](../data/results/cofold_metrics.csv), input JSONs in [`data/results/inputs/`](../data/results/inputs/).
- **Self-consistency:** `tests/test_demo.py` recomputes `cofold-fzd5`'s `PAE_IF`/`epitope_reprod` from the committed CIF+PAE and asserts it reproduces the CSV (±0.05 / ±0.03) — the structures and the extraction path agree.
- **Determinism:** two builds produce byte-identical artifacts; no committed number contradicts its source CSV. Cofolds re-run free via the boltz.bio prediction IDs (idempotency keys `cofold-*`).
- **Verdict rule (frozen):** *hit* iff `PAE_IF < 9.8 Å` **and** `epitope_reprod ≥ 0.55` (mid-band of the PD-1/lysozyme calibrators); `ipTM` is **not** used (over-docks); calibrators keep their ceiling/floor roles; a known off-target that fails the bar reads *missed*; a valid-regime hit on a <150 aa ectodomain is flagged *low-confidence*. Frozen before the screen. See [`src/crossflag/demo/panel.py`](../src/crossflag/demo/panel.py) and [`data/results/screen/rig.json`](../data/results/screen/rig.json).
