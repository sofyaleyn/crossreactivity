# Findings — empirical viability log (2026-07-04)

Durable record of the viability investigation run on 2026-07-04. Design implications live in `pivot-spec.md`; this file is the experimental log + current verdict. All work used the SHR-1210 (camrelizumab) anti-PD-1 anchor.

> **Raw per-run data:** [`data/results/`](data/results/) — `cofold_metrics.csv` holds the full metrics for all **24 cofold runs** (ipTM max/mean, protein-ipTM, pTM, complex-pLDDT, structure-confidence, binding-confidence, PAE_IF, epitope-reprod, boltz.bio prediction IDs); plus the Tier-1 filter scores and fold-matched-benchmark scores. Human-readable tables in [`data/results/README.md`](data/results/README.md).

---

## Environment & tooling stood up

- **Compute:** Apple M2 (Metal/MPS, **no CUDA**). Local base Python 3.14 lacks ML wheels → created conda env **`crossflag-spike`** (Python 3.10) with `torch 2.12 (CPU)`, `fair-esm`, `numpy`, `biopython`, `requests`.
- **Cofolding:** hosted **boltz.bio API** (Boltz-2.1) via CLI `boltz-api` (`~/.local/bin`, OAuth device-login) — since no local GPU. Endpoints: `predictions:structure-and-binding {estimate-cost,run}`. Cost ~**$0.05/cofold** (5 samples ≈ $0.20). Total spend this session ≈ **$1.6**.
- Scratchpad (scripts, inputs, run outputs): `/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/<session>/scratchpad/`.

---

## Experiment 1 — Rung-1 method soundness (analytical + literature)

**Verdict: the original rung-1 (`cosine(antibody-PLM vector, ESM-2 antigen vector)`) is broken.** (1) Undefined as written — 512/480-d antibody vectors vs 1280-d ESM antigen vectors, no alignment. (2) Biophysically backwards — binding is complementarity, not similarity; embedding *similarity* indicates two antibodies bind the **same** antigen (epitope binning), not that an antibody binds a protein. (3) No published off-target/polyreactivity method uses unsupervised cross-PLM cosine. Full reasoning in `pivot-spec.md §Why rung-1 changed`.

## Experiment 2 — Embedding discrimination spike (16 antigens)

Same-model ESM (charitable, well-defined) paratope↔antigen cosine, WT SHR-1210 vs 3 off-targets + 13 decoys (incl. 5 hard receptor decoys):

| Formulation | VEGFR2 | FZD5 | ULBP2 | Off-targets beat hard receptor decoys? |
|---|---|---|---|---|
| whole-Fv | #1/16 | #8 | #9 | **No** |
| paratope (H3+L3) | #1/16 | #10 | #4 | **No** |

VEGFR2 ranks #1 by **fold/size confound** (largest Ig-fold receptor), not binding. FZD5/ULBP2 lose to non-target receptors.

## Experiment 3 — Anchor data verification

- WT SHR-1210 VH/VL: **real** (matches Finlay et al. mAbs 2019, PMID 30541416, Table 1).
- VEGFR2 / FZD5 / ULBP2: **real** (UniProt P35968 / Q13467 / Q9BZM5).
- **Germlined-mutant panel: NOT published — must be reconstructed** (recipe only: light chain near-germlined to IGKV1-39 + JK1). Beat 2 is therefore illustrative. Detail in `pivot-spec.md §Anchor data status`.

## Experiment 4 — Calibrated cofold panel (Boltz-2.1, 5 samples each)

Anchored by a true-target **ceiling** (PD-1, the cognate) and a non-binder **floor** (lysozyme). Metrics: ipTM (found uninformative — over-docks), **PAE-mean** (paratope↔antigen), **epitope reproducibility** (mean pairwise Jaccard of contacted antigen residues across the 5 samples). Chain identified by antigen length; heavy-atom contact cutoff 5 Å.

| complex | ipTM | PAE-mean | epitope reprod | read |
|---|---|---|---|---|
| WT × PD-1 (cognate, **ceiling**) | 0.910 | 6.41 | **0.936** | binder |
| WT × **FZD5**-CRD (off-target) | 0.952 | 4.91 | 0.660 | **confirmed** |
| WT × **ULBP2** (off-target) | 0.940 | 5.07 | 0.893 | **confirmed** |
| WT × **VEGFR2** D2-3 (off-target) | 0.883 | 10.45 | 0.430 | **floor / unconfirmed** |
| MUT × VEGFR2 D2-3 (germlined) | 0.852 | 12.88 | 0.512 | floor |
| WT × VEGFR2 D2-3 **+3V2A template** | 0.766 | 17.05 | 0.302 | worse when antigen pinned |
| MUT × VEGFR2 +template | 0.832 | 17.07 | 0.068 | floor |
| WT × lysozyme (non-binder, **floor**) | 0.869 | 11.56 | 0.446 | non-binder |

**Key results:**
- **2 of 3 known off-targets (FZD5, ULBP2) cofold as confident, reproducible binders** at/near the PD-1 ceiling and clearly above the lysozyme floor. The confirmation rung works for most off-targets.
- **VEGFR2 is the outlier at the floor** — the hardest case: lowest-affinity ("aberrant… low affinity" per Finlay), largest/multi-domain, and epitope/domain unmapped (D2-3 assumed). **Templating VEGFR2 to its real 3V2A structure made it *worse*** — untemplated, Boltz molds the deformable antigen into a spurious pose; pinned to the real fold, the antibody can't dock. WT-vs-mutant contrast was inside the noise band.
- **lysozyme (small antigen) stays at the floor**, so FZD5/ULBP2 confidence is not a small-antigen artifact.
- **ipTM is not usable** (lysozyme 0.87 ≈ everything). Use PAE-mean + epitope reproducibility vs the calibrated panel.

## Experiment 5 — Fold-matched decoy ranking benchmark (the confound-controlled ranking test)

Built a benchmark where, within each structural family, **one member is a known off-target and its close relatives were in Finlay's proteome screen and NOT flagged** (legitimately negative, same fold). Families: RTK/Ig (VEGFR2 + 14 incl. VEGFR1/3, EGFR, HER2, INSR, IGF1R…); Frizzled/CRD (FZD5 + 11 incl. all FZDs, SMO); MHC-I-like (ULBP2 + 12 incl. ULBP1/3, MICA/B, HLA). Data: `benchmark_antigens.json`.

Embedding method (same-model ESM paratope cosine, best-case), ranking each positive **within its fold-matched family**:

| family | positive | global rank | **within-family AUROC** |
|---|---|---|---|
| RTK/Ig | VEGFR2 | 7/40 (top 18%) | 0.79 |
| Frizzled | FZD5 | 34/40 | 0.55 (chance) |
| MHC-I-like | ULBP2 | 18/40 | **0.42 (below chance)** |
| **mean** | | | **0.58 ≈ chance** |

**The cheap embedding fails the confound-controlled test.** High *global* rank for VEGFR2 is the fold confound; forced to pick the true off-target out of its *own family*, embedding is at chance (0.58; ULBP2 below chance). Embedding cannot be the Stage-1 shortlisting filter.

---

## Experiment 6 — Cofold within-fold discrimination (FZD5 family)

Cofolded FZD5's 11 fold-matched Frizzled relatives (all screened-negative). Full per-member metrics in `data/results/README.md`; ranked by PAE_IF (whole-interface, lower=tighter):

| Frizzled CRD | ipTM_max | PAE_IF | epitope_reprod |
|---|---|---|---|
| SMO | 0.968 | **4.98** | 0.652 |
| **FZD5 (true off-target)** | 0.952 | **5.69** | 0.660 |
| FZD1 | 0.948 | 7.66 | 0.885 |
| FZD10 | 0.902 | 8.24 | 0.756 |
| FZD8 | 0.898 | 10.32 | 0.707 |
| FZD2 | 0.894 | 10.36 | 0.607 |
| SFRP1 / FZD7 / FZD9 / FZD4 | 0.84–0.95 | 11.0–12.1 | 0.27–0.51 |
| FZD3 / FZD6 | 0.83–0.84 | 15.4 / 16.1 | 0.23 / 0.29 |

FZD5 ranks **2/12 by PAE_IF, within-family AUROC 0.909** (epitope-reprod weaker: rank 4/12, AUROC 0.727). One within-fold false positive: **SMO** (the most divergent member, out-cofolds FZD5 on PAE). → **Cofold is a real within-fold discriminator**, not just a "cofoldable-antigen" detector.

## Experiment 7 — Surface method on the fold-matched benchmark

SASA-restricted ESM + physicochemical surface-complementarity (AlphaFold-DB structures). Best (surface-ESM) mean within-family AUROC ~0.62 (0.61–0.70 across thresholds), only marginally above the 0.58 embedding baseline, driven entirely by Frizzled; MHC-I-like stuck at 0.42 (below chance); physicochemical complementarity 0.554. → **No cheap surface proxy clears ~0.6 or robustly beats embedding. Cheap Stage-1 shortlisting fails.**

## Experiment 8 — Antibody-side specificity control (pembrolizumab)

Cofolded **pembrolizumab** (anti-PD-1, lacks SHR-1210's off-targets; Fv from PDB 5GGS) against the same 4 antigens. Whole-interface metrics (ipTM · PAE_IF · epitope-reprod):

| antigen | SHR-1210 | pembrolizumab |
|---|---|---|
| PD-1 (shared target) | 0.910 · 7.2 · 0.936 | 0.963 · 3.6 · 0.968 |
| VEGFR2 | 0.883 · 11.6 · 0.430 | 0.762 · 16.5 · 0.642 |
| FZD5 | 0.952 · 5.7 · **0.660** | 0.755 · 18.6 · **0.217** |
| ULBP2 | 0.940 · 5.7 · **0.893** | 0.901 · 13.9 · **0.341** |

Both bind the shared target PD-1; for the two off-targets SHR-1210 confirmed, **pembrolizumab collapses** (FZD5 reprod 0.66→0.22, ULBP2 0.89→0.34). → **The confirmations are SHR-1210-CDR-specific, not generic dockability. Specificity control PASSES.**

## Experiment 9 — Cheap-filter test on the real curated set (`data/reference/self_proteins.csv`, 2,896 proteins)

"Tier 1": can any cheap (non-cofold) score rank the 3 known SHR-1210 off-targets (all confirmed present) near the top of the real 2,896-protein surfaceome, so a recall-safe cutoff retains them while cutting the list? Metric = **retain-all-3** (fraction you must keep to hold all 3; lower=better).

| filter | VEGFR2 | FZD5 | ULBP2 | retain-all-3 |
|---|---|---|---|---|
| annotation "is-it-a-receptor" (FREE baseline) | 3% | 38% | 9% | **0.377** |
| embed_wholeFv (ESM) | 14% | 30% | 42% | 0.415 |
| embed_paratope (ESM) | 9% | 64% | 20% | 0.645 |
| biophys (hydrophobicity / charge / length) | — | — | — | 0.77–0.92 |
| naive CDR sequence identity | 17% | 55% | 100% | 0.996 |

**No learned filter beats the free annotation baseline (0.377); embeddings add nothing over "keep the receptors."** Even the best requires keeping **38%** of the set to retain all 3 (bottlenecked by FZD5, which every learned filter buries — embed_paratope ranks it 1867/2896). Decision rule (promising only if all 3 in top ~10–20% AND beats annotation) → **FAILED. Tier 2 ($580) not run.** Confirms the within-family benchmark (Exp 5) on real diverse data: the only cheap signal is "off-targets are receptors," captured for free by annotation. **No cheap route to a small per-candidate shortlist exists.**

## Current verdict (2026-07-04, revised)

**The cofold confirmer is validated; the cheap-triage idea is not; and the cost reality reframes the pipeline.**

| Capability | Status |
|---|---|
| Cofold **confirmation** (Stage 2) | ✅ Confirms 2/3 off-targets; discriminates within-fold (AUROC 0.909); antibody-specific (pembrolizumab control passes) — **three independent validations** |
| Cheap **embedding / surface** shortlisting (Stage 1) | ❌ Fails within-fold enrichment (~0.58–0.62 ≈ chance) |
| Antibody-only polyreactivity ("variant is sticky") | ✅ Viable, but doesn't *name* the off-target |

**Reframe:** hosted cofold ≈ **$0.20 each** → a curated few-hundred-protein reference set ≈ **$50–150/antibody**, trivial vs a wet-lab specificity screen (~$10–30k). The pivot spec's premise that cofold is "too expensive to run against every self-protein" (hence a cheap pre-filter is required) is **wrong at these prices.** Drop the broken Stage-1; cofold the whole curated set directly.

**Revised pipeline:** Stage 0 curate (bounds recall) → **cofold ALL** → rank by PAE_IF / epitope-reproducibility against the calibrated panel (PD-1 ceiling / lysozyme floor). **Limits:** VEGFR2-class (weakest) off-targets missed → imperfect recall; within-fold false positives (SMO); recall bounded by curation.

**Bottom line: viable as a validated cofold-based off-target screen for a curated reference set** — not the original embedding-triage design, but empirically defensible.

## Next steps → the build plan lives in [`plan.md`](plan.md)

The active, detailed build plan — the **representative-set off-target screen** (cluster candidates → cofold representatives to discover the shortlist → screen all candidates against it), with the module list, cost model, quantitative K-sizing, and anchor-validation steps — is in **`plan.md`**. This file (`findings.md`) is the evidence/verdict record that plan builds on. Carry-over items folded into the plan: freeze the calibrated scoring rig; resolve the VEGFR2 wrong-epitope caveat; second anchor (ABT-736→PF4) for generalization; recall/false-positive handling. Separate follow-up (out of scope): a learned cheap interaction filter distilled from cofold labels ("Option C").

## Non-goals / dropped

- Cross-PLM or same-space **embedding ranking** (broken + fails fold-matched test).
- **Surface-complementarity pre-filter** (fails fold-matched test at these prices — cofold-all is cheaper than making a filter work).
- Claiming **in-silico confirmation of the weakest off-targets** (VEGFR2-class) — honestly out of reach; wet lab confirms.

## Reproducibility

Scripts in scratchpad: `spike.py` (Exp 2), `analyze_all.py` (Exp 4 metrics), `build_benchmark.py` + `benchmark_antigens.json` (Exp 5), `analyze_pembro.py` (Exp 8), `surface_score.py` (Exp 7), `prep_*.py` (input builders). Cofold runs under `scratchpad/boltz-runs/cofold-*`; idempotency keys `cofold-*` (re-runnable without new charges). Tooling: conda env `crossflag-spike`, `boltz-api` CLI (boltz.bio, OAuth).
