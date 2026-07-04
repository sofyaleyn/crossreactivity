# Implementation plan — representative-set off-target screen

**Status:** active build plan (2026-07-04). Evidence basis in [`findings.md`](findings.md); supersedes the original specs in `docs/`.

## Goal & scope

De-risk a wet-lab specificity plate for a panel of **~100–1000 antibody candidates** (variants of one lead) against a **curated ~2,896-protein surfaceome** (`data/reference/self_proteins.csv`), using **Boltz-2 cofold** as the validated off-target confirmer — **without** running the full candidate × target grid.

**In scope:** cluster candidates → discover the shared off-target shortlist by cofolding representatives → screen all candidates against the shortlist → rank candidates + design the plate.

**Out of scope (follow-up project):** a learned cheap interaction filter ("Option C" — an antibody×antigen binding classifier distilled from cofold-generated labels). Justified: `findings.md` Exp 5/7/9 show no *cheap* unsupervised filter (embedding, surface, biophysical) beats a free annotation prior or reaches usable recall; a supervised model is a separate research effort needing a low-affinity/off-target training set + hard negatives. It can be bootstrapped later from labels this pipeline produces.

## Why this design (one line each; detail in findings.md)

- **Cofold confirmer is validated** — confirms 2/3 known off-targets, discriminates within-fold (AUROC 0.909), antibody-specific (pembrolizumab control). (Exp 4/6/8)
- **No cheap filter works** — all ≈ chance on fold-matched decoys and lose to free annotation on the real 2,896-set. (Exp 5/7/9)
- **So reduce the grid by candidate *similarity*, not a filter** — variants of a lead share ~all their off-targets; cofold a few representatives to discover the shortlist, then screen every candidate against it.

## Starting point (environment & reusable code — already set up)

A fresh context can build immediately; a validated toolchain and env already exist (details in `findings.md §Environment` + `§Reproducibility`):
- **Python env:** conda env `crossflag-spike` (Python 3.10) with `torch` (CPU), `fair-esm`, `numpy`, `biopython`, `requests`. Run: `conda run -n crossflag-spike python …`. (Local base Python 3.14 lacks ML wheels; no CUDA — Apple M2.)
- **Cofolding:** boltz.bio hosted API via the `boltz-api` CLI (`~/.local/bin`, OAuth session already authenticated). Pattern: `boltz-api predictions:structure-and-binding {estimate-cost,run} --input @json://<file> --model boltz-2.1 --idempotency-key <k> --root-dir <dir>`. ~$0.20/cofold (5 samples); `sampling_steps ≥ 50`; input = `protein_protein_binding` with VH/VL + antigen chains.
- **Reusable code & data** in [`data/results/`](data/results/) (adapt into `src/crossflag/`): `scripts/analyze_all.py` (cofold metrics: PAE_IF, epitope-reproducibility, chain-by-length), `scripts/prep_*.py` + `inputs/cofold_*.json` (Boltz input builders), `scripts/build_benchmark.py` (ESM embedding), `structures/` (CIFs + PAE), `cofold_metrics.csv`, `tier1_filter_scores.csv`. Cofolds re-run free via prediction IDs (idempotency keys `cofold-*`).
- **Data:** `data/reference/self_proteins.csv` (2,896 reference set + sequences), `data/anchor/` (SHR-1210 WT VH/VL, germline variants, off-target FASTAs).

## Pipeline

```
Stage 0  Inputs + calibrated scoring rig
Stage 1  Cluster candidates → K representatives        (biochemical paratope fingerprint)   ← build first
Stage 2  Discovery screen: reps × 2,896  → cofold      → off-target SHORTLIST
Stage 3  Coverage check: saturation / Chao1            → add reps until recall target met
Stage 4  Panel screen: ALL candidates × shortlist      → per-candidate off-target profile
Stage 5  Rank candidates, name off-targets, design the wet-lab plate
```

## Stage 0 — Inputs & scoring rig

- **Candidates:** VH/VL FASTA pairs (variants of a lead). Anchor demo set: `data/anchor/variants/` (WT + L1/L3/L1L3 germline).
- **Reference set:** `data/reference/self_proteins.csv` (2,896 surface proteins, sequences + 3-layer clinical-severity flags; assembly in [`reference-set.md`](docs/reference-set.md)). Cofold against the ectodomain/soluble region where the sequence has TM spans.
- **Scoring rig (frozen):** per antibody, a calibrated panel — **ceiling** = the cognate target (e.g. PD-1), **floor** = a non-binder (lysozyme). Ranking metrics: **epitope reproducibility** (Jaccard of contacted antigen residues across the 5 cofold samples; primary) and **PAE_IF** (whole antibody↔antigen interface PAE; secondary). **Do not use ipTM** (over-docks). A protein is a "hit" if its metrics sit toward the ceiling, well above the floor.

## Stage 1 — Candidate clustering (biochemical) — **the focus**

Cluster on the axis that controls off-targets: the **paratope biochemistry**. Two variants in the same cluster should share off-target profiles, so one can stand in for the others.

**Feature = biochemical paratope fingerprint (per candidate):**
1. Annotate the 6 CDRs (H1–H3, L1–L3). Proper tool: **ANARCI/IMGT**. For a variant panel of one lead (near-identical, often same length), map the lead's CDR spans onto each variant by alignment (positional) — the demo uses this.
2. Per CDR loop compute biochemical descriptors: **net charge** (K/R +1, D/E −1, H +0.1), **mean hydrophobicity** (Kyte–Doolittle), **aromatic fraction** (F/W/Y), **length**. (Optionally patch metrics from an IgFold model of the lead, restricted to solvent-exposed paratope residues — refinement, not required.)
3. Concatenate into a fixed-length vector (6 CDRs × 4 features ≈ 24-d) + global paratope net charge & total aromatic count. **Up-weight CDR-H3 and L3** (dominant paratope contributors).

**Distance & clustering:**
- Z-score standardize features across the candidate panel; apply loop weights.
- Distance = Euclidean; cluster with **k-medoids** (or agglomerative/Ward), or **sphere-exclusion** at radius *r*. Representatives = cluster **medoids**.
- Sphere-exclusion gives a guarantee: every candidate lies within *r* of a representative, so *r* bounds how different any variant is from its stand-in.

**Alternative for clean single-lead panels:** encode each candidate by its **mutations from the lead**, each weighted by (surface exposure × physicochemical Δ). Buried/framework mutations get ~0 weight. More precise than a global fingerprint when candidates are point-variants of one lead.

## Stage 3 sizing — how many representatives K (quantitative)

K is **measured, not guessed**:
1. **Calibrate radius *r*** on the anchor: SHR-1210 WT + germline variants + their known cofold profiles → how much paratope change flips an off-target → largest *r* under which the profile is stable.
2. **Off-target accumulation (rarefaction) curve:** cumulative unique off-targets vs #representatives cofolded. Because variants share off-targets, it saturates fast. **K\*** = the knee (marginal gain < threshold, e.g. <1 new off-target per 5 reps).
3. **Recall guarantee:** apply a **Chao1 / Good–Turing** unseen-estimator to the curve → estimated missed-off-target fraction; pick K\* so it's below tolerance (e.g. <5%).
4. **(Optional) sentinel bootstrap:** cofold ALL candidates × a small diverse sentinel panel (~30 proteins, ~$600) to get an *empirical* off-target fingerprint; cluster + build the saturation curve on measured behavior rather than the biochemical proxy.

## Stages 2, 4, 5

- **Stage 2 (discovery):** cofold each of the K representatives × 2,896 → union of hits = **shortlist** (~tens of proteins).
- **Stage 4 (panel):** cofold every candidate × shortlist → per-candidate off-target profile (which off-targets each variant retains/loses).
- **Stage 5 (report):** rank candidates by off-target burden (count + severity vs the calibrated panel); recommend advance/drop; name each flagged off-target + its confirmation assay; output the prioritized wet-lab plate.

## Cost model (per ~100-candidate panel)

| Step | Cofolds | ~Cost |
|---|---|---|
| Sentinel bootstrap (optional) | 100 × 30 = 3,000 | ~$600 |
| Discovery: K≈10 reps × 2,896 | ~29,000 | ~$5,800 |
| Panel: 100 × shortlist(~20) | 2,000 | ~$400 |
| **Total** | | **~$4–7k** |

vs full grid 100 × 2,896 = 289,600 cofolds ≈ **$58k**, and vs a wet-lab MPA screen ~$10–30k. Scales with K (saturation-bounded), not with panel size × set size.

## Anchor validation (before real panels)

1. **Beat 1 (real):** cofold WT SHR-1210 × all 2,896 → confirm FZD5/ULBP2 surface near the top of the calibrated ranking (~$580).
2. **Clustering sanity:** the germline variants should cluster away from WT along the mutated CDRs (built below).
3. **Radius calibration:** relate germline-variant cofold-profile change to their biochemical distance from WT.

## Modules to build

- `src/crossflag/candidate_clustering.py` — biochemical fingerprint + distance + clustering + representative selection. **← building now.**
- `src/crossflag/scoring.py` — cofold metric extraction (epitope-reproducibility, PAE_IF) + calibrated-panel call. (adapt `scratchpad/analyze_all.py`)
- `src/crossflag/cofold_client.py` — boltz.bio submit/poll/collect wrapper (adapt existing CLI usage).
- `src/crossflag/discover.py` — Stage 2 discovery + Stage 3 saturation/Chao1.
- `src/crossflag/panel.py` — Stage 4 panel screen + Stage 5 ranking/report.
- `src/crossflag/reference.py` — load/prepare the 2,896 set (ectodomain extraction).

## Risks & mitigations

- **A variant invents a new off-target** not in the reps' shortlist → cluster on the mechanism-aligned axis, calibrate *r* conservatively, and full-screen a diverse representative subset (Chao1 covers residual risk).
- **Clustering proxy error** (biochemistry ≠ exact off-target profile) → validate/replace with the empirical sentinel fingerprint.
- **Weakest off-targets missed** (VEGFR2-class, low-affinity) → known cofold limitation; flag as reduced recall for very weak binders; wet-lab is the backstop.
- **Recall bounded by the 2,896 curation** → document; expand the reference set to raise the ceiling.
