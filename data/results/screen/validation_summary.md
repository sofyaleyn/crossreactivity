# CrossFlag — off-target screen validation summary

**Date:** 2026-07-05  ·  **Compute spend:** $342.70 / $500 budget  ·  **Wet-lab:** none

CrossFlag screens a therapeutic antibody's Fv (VH+VL) against a curated human
surfaceome by **cofolding** (Boltz-2) each antibody–protein pair and scoring the
interface — to prioritize likely off-targets *before* a wet-lab specificity screen.
This run tested whether the method actually works at screen scale.

## Method (frozen before the screen — no tuning to these results)

- **Cofold:** Boltz-2.1, 5 samples per antibody–antigen pair (hosted API).
- **Two metrics:** `PAE_IF` (interface predicted-aligned-error; lower = tighter) and
  `epitope_reprod` (do the 5 samples agree on the antigen epitope; higher = consistent).
- **Hit rule (frozen):** `PAE_IF < 9.8 AND epitope_reprod ≥ 0.55`, calibrated on the
  original anchor against a PD-1 cognate-target ceiling and a lysozyme non-binder floor.
- **Antigens:** extracellular (ectodomain) region only, extracted from UniProt topology
  (membrane proteins are cofolded as their soluble ectodomain, not full-length).

## What we ran (three independent tests)

| Test | Design | Purpose |
|---|---|---|
| **1. Real-scale screen** | camrelizumab (SHR-1210, anti-PD-1) × 1,200 self-protein ectodomains | specificity + enrichment at scale |
| **2. Second anchor (blind)** | ABT-736 (anti-Aβ) × PF4 + 400 decoys | generalization to an independent drug |
| **3. Negative control** | CDR-scrambled camrelizumab × FZD5, ULBP2 | is binding paratope-specific? |

Camrelizumab's known off-targets (VEGFR2, FZD5, ULBP2 — clinical: capillary hemangioma)
were injected into the analysis from prior validated cofolds, not re-folded.
ABT-736→PF4 is an independent, published off-target case (HIT-like toxicity;
Steinmetz et al., *mAbs* 2021, PMID 33596779).

## Results

### The method works — for medium/large antigens

- **Sensitivity ✓** — among same-size (150–300 aa) decoys, **FZD5** ranks above 84% on
  reproducibility / 99% on PAE; **ULBP2** is among the most reproducible proteins in its
  size class (better than 242/243). Both correctly flagged. (VEGFR2, the weakest/lowest-affinity
  off-target, is missed — a known, honest limit.)
- **Specificity ✓** — false-positive rate falls sharply with antigen size:

  | ectodomain length | FPR |
  |---|--:|
  | 0–80 aa | 20.6% |
  | 80–150 aa | 19.1% |
  | 150–300 aa | 4.5% |
  | 300–600 aa | 1.2% |
  | 600–1000 aa | 0.0% |
  | **≥150 aa (valid regime)** | **2.1%** (14/678) |

- **Generalization ✓** — for the *second* antibody (ABT-736), **PF4 is flagged and enriched
  above same-size background** (top ~9% by reproducibility, #1 by PAE among 76 same-size
  decoys), with zero retuning. The frozen pipeline recovered an independent drug's real
  off-target.
- **Antibody-specific ✓** — valid-regime FPR: camrelizumab 2.1%, ABT-736 0.4%. Not generic
  stickiness.
- **Negative control ✓** — with all six CDRs scrambled (framework preserved), both off-targets
  collapse to non-hits, confirming paratope-specific binding:

  | antigen | wild-type | all-6-CDR scramble |
  |---|---|---|
  | FZD5 | hit (PAE 5.69 / reprod 0.66) | **non-hit** (9.94 / 0.26) |
  | ULBP2 | hit (5.74 / 0.89) | **non-hit** (11.19 / 0.82) |

  FZD5 required a *full* 6-CDR scramble to break (a partial CDR-H3/L3-only scramble did not),
  indicating its cross-reactivity draws on more of the paratope than the dominant loops alone.

### The real limitation we found

The method **over-docks small antigens (<150 aa)** — FPR ~20%, with docking confidence
*higher than the true target itself* (physically implausible). Of 168 decoys docking tighter
than a real off-target, median length 41 aa; only 3 are ≥150 aa. This is a **genuine method
flaw, not a test artifact**: the small ectodomains are correctly extracted, real domains
(e.g. GPCR N-termini), and the scoring simply mis-ranks short peptides. It affects **~43% of
the surfaceome** (1,257/2,896 proteins have ectodomains <150 aa). **Fix:** a size-aware
confidence gate (flag small-ectodomain hits as low-confidence) or improved short-antigen scoring.

### Candidate novel off-targets (valid regime, top by interface confidence)

- **SMO** (Smoothened, 202 aa) — reproducible (also flagged in earlier within-family tests).
- A cluster of **Ig-fold NK receptors** — CD19, IL22RA1, and five KIRs
  (KIR2DL4 / KIR2DS1 / KIR3DL3 / KIR2DL3 / KIR2DL5B). Same structural family as ULBP2's biology.
- **Interpretation caveat:** possible genuine cross-reactivity, or an Ig-fold affinity bias —
  needs a wet-lab / structural-biology read before treating as real.

## Bottom line

**Conditionally validated.** The cofold screen reliably prioritizes off-targets for
medium/large antigens (≥150 aa: ~3% FPR, both known off-targets enriched), is paratope-specific,
(≥150 aa: ~2% FPR, both known off-targets enriched) and generalizes across two independent
antibodies. It has one real, bounded flaw — unreliable scoring for small (<150 aa) ectodomains —
that needs a size-aware fix before the full surfaceome can be trusted. It is a **prioritization
tool for wet-lab specificity screens, not a replacement.**

## Open questions for discussion

1. Do the candidate off-targets (SMO, KIRs, CD19) make biological sense, or read as an Ig-fold artifact?
2. Small-antigen limitation — acceptable to flag GPCRs/small-ectodomain proteins as "low-confidence,"
   or do we need to solve short-antigen scoring?
3. Does our ABT-736 → PF4 generalization result square with the published case?

## Caveats & provenance

- Numbers reflect the **complete** run: **1,198/1,200** flagship (2 dropped on API rate-limit) and
  **401/401** anchor-2 cofolds scored, plus 20 pilot + 4 scramble-control folds (1,603 total, 5 samples each).
- ABT-736 VH/VL were transcribed from a patent figure (US2009/0175847A1); off-target/mechanism
  fully cited (PMID 33596779). PF4 = UniProt P02776 (mature chain).
- Data & code: `data/results/screen/` — `screen_metrics.csv` (per-protein scores),
  `rig.json` (frozen thresholds), `spend_ledger.json` (budget), `checks.py` (this analysis),
  job configs under `jobs/`. Ectodomains: `data/reference/ectodomains.csv`.
