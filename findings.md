# Findings — empirical viability log (2026-07-04)

Durable record of the viability investigation run on 2026-07-04. Design implications live in `pivot-spec.md`; this file is the experimental log + current verdict. All work used the SHR-1210 (camrelizumab) anti-PD-1 anchor.

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

## Current verdict (2026-07-04)

**CONDITIONAL — the project's viability now rests on one narrow, still-open question.**

| Capability | Status |
|---|---|
| Antibody-only polyreactivity ("variant is sticky") | ✅ Viable (physicochemical descriptors) — but doesn't *name* the off-target |
| Cofold **confirmation** of an off-target (Stage 2) | ✅ **Works for 2/3** (FZD5, ULBP2); VEGFR2 (low-affinity/wrong-epitope) at floor |
| Cheap **embedding shortlisting** (Stage-1 candidate) | ❌ Fails fold-matched test (AUROC 0.58 ≈ chance) |
| A cheap Stage-1 that enriches within-fold to feed the cofold | ❓ **Untested — the make-or-break question** |

Interesting inversion: embedding ranks VEGFR2 okay but cofold can't confirm it; cofold nails ULBP2/FZD5 but embedding ranks them at/below chance. The expensive cofold is the stronger, more reliable discriminator — but it can't be run proteome-wide, so a working cheap Stage-1 is required and not yet found.

## Open checks (in progress)

1. **Does cofold *discriminate* within fold?** Cofold the fold-matched decoys of a confirmed family (FZD5 vs other Frizzled CRDs) and compute cofold-based within-family AUROC. Tests whether FZD5's confident interface is *special* vs its relatives, or whether all Frizzleds cofold confidently (i.e., cofold is a real discriminator, not just a "cofoldable-antigen" detector).
2. **Does a surface method pass the benchmark where embedding failed?** Score the fold-matched benchmark with a surface-complementarity / SASA-restricted method (AlphaFold-DB structures) and compare within-family AUROC to the 0.58 embedding baseline.

## Reproducibility

Scripts in scratchpad: `spike.py` (Exp 2), `analyze_all.py` (Exp 4 metrics), `build_benchmark.py` + `benchmark_antigens.json` (Exp 5), `prep_*.py` (input builders). Cofold runs under `scratchpad/boltz-runs/cofold-*`; idempotency keys `cofold-shr1210-*` (re-runnable without new charges).
