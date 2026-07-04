# Roadmap — ranked steps and tools

> **⚠️ SUPERSEDED (2026-07-04).** The Phase-1 embedding rung was found broken and the ladder reordered (cofold is the validated primary signal, not a top-of-funnel extra). See [`findings.md`](findings.md) and [`pivot-spec.md`](pivot-spec.md). Kept for history.

Read before the specs. Ranks every pipeline step and tool, marks the MVP / extension boundary, gives a 15h budget. Abbreviations in `glossary.md`.

---

## Phase overview

| Phase | Name | Goal | Ships? |
|---|---|---|---|
| **1** | MVP — embedding rank | Rank the panel by antibody-PLM ↔ antigen-PLM embedding similarity, validated on SHR-1210 | **Must ship** |
| **2** | Confirmation ladder | Add 3D surface fingerprint (rung 2) and Boltz-2 cofolding (rung 3) on the top suspects | Add if MVP done by hour 9 |
| **3** | Stretch | Second anchor, heatmap viz, broader reference set | Only if Phase 2 lands |

**Phase 1 meets the success criterion on its own** (embedding-level). Cofolding makes the demo far more convincing but is not required to have *a* demo.

---

## Steps, ranked

`[MVP]` = Phase 1. `[EXT]` = Phase 2–3.

| # | Step | Phase | Notes |
|---|---|---|---|
| 1 | **Data assembly** — SHR-1210 VH/VL, build CDR-germlining variant panel (Finlay et al.), collect VEGFR2/FZD5/ULBP2 + PDB IDs, assemble curated self-protein reference set | `[MVP]` | Do first. Verified anchor + reference set = everything. |
| 2 | **CDR extraction** — ANARCI/ANARCII → CDR annotation per variant | `[MVP]` | Localizes the paratope; used to *weight/pool* embeddings, not to string-compare. |
| 3 | **Antibody embedding** — AntiBERTy or AbLang2 → per-residue + CDR-pooled embedding per variant | `[MVP]` | The paratope representation. Replaces the old BLOSUM floor. |
| 4 | **Antigen embedding** — ESM-2 → embedding for each self-protein in the reference set (CDR-facing surface region where structure known) | `[MVP]` | The self-protein representation. |
| 5 | **Embedding-similarity ranking + calibration** — score each variant against the reference set in embedding space; calibrate vs. a benign-protein background | `[MVP]` | Rung 1. Enrichment, not raw distance. |
| 6 | **Panel ranking + assay map** — rank variants; name flagged self-proteins; recommend assay | `[MVP]` | The payoff. WT ranks risky; germlined mutants rank cleaner. |
| 7 | **LLM agent narration** — write the recommendation | `[MVP]` | Retrieval + reasoning; no training. |
| 8 | **IgFold structure + surface fingerprint** — fold top-flagged variants, fingerprint CDR surface (dMaSIF, or APBS fallback), compare vs. flagged self-protein surface | `[EXT]` | Rung 2: corroborate the shape. GPU (dMaSIF) or CPU (APBS fallback). |
| 9 | **Boltz-2 cofolding gate** — cofold top ~5 variant×self-protein pairs; score interface confidence (ipTM/PAE at CDR–epitope contacts) + affinity | `[EXT]` | Rung 3: does it actually dock? The money-shot for the demo. |
| 10 | **Second anchor / viz / broader set** | `[EXT/stretch]` | Polish. |

**Critical path:** 1 → 2 → 3 → 4 → 5 → 6 → 7. Rungs 2–3 (steps 8–9) branch off after step 7 works.

---

## Tools, ranked by tier

Full install/API/license in `tools.md`.

### Tier 0 — MVP (must have)

| Tool | Full name | Role | Access |
|---|---|---|---|
| **Python + Biopython** | — | parsing, glue | local, pip |
| **ANARCI / ANARCII** | Antibody Numbering and Antigen Receptor ClassIfication | CDR annotation | local, bioconda/pip |
| **AntiBERTy** *or* **AbLang2** | antibody-specific PLM (protein language model) | embed variant CDRs/paratope | local, pip (`antiberty` / `ablang2`) |
| **ESM-2** | Evolutionary Scale Modeling v2 (general PLM) | embed self-protein antigens | local, pip (`fair-esm` / HF) |
| **LLM** | large language model | agent narration | API |

### Tier 1 — MVP polish

| Tool | Role |
|---|---|
| CDR-pooling / attention-weighting util | turn per-residue embeddings into a paratope vector |
| background calibration | benign-protein reference distribution |

### Tier 2 — Confirmation ladder (extension)

| Tool | Full name | Role | Access | Cost / license |
|---|---|---|---|---|
| **IgFold** | antibody folding model | 3D antibody structure per variant | local, pip | CPU-ok; BSD |
| **dMaSIF** | differentiable Molecular Surface Interaction Fingerprinting | CDR surface fingerprint (rung 2) | local, PyTorch, pretrained weights | **GPU; CC-BY-NC-ND** |
| *fallback:* **APBS + MSMS** | Poisson-Boltzmann solver + surface mesher | hand-rolled surface descriptor, no GPU | local | open |
| **Boltz-2** | co-folding + affinity foundation model | rung 3: cofold top suspects, interface + affinity score | local, pip (`boltz`) | **GPU; MIT license** (commercial-ok) |

### Tier 3 — Stretch

| Tool | Role |
|---|---|
| **AbodyBuilder2 / ABodyBuilder3** | alt antibody structure cross-check |
| **SAbDab** | look up any deposited camrelizumab crystal structure |

---

## What we do NOT use (and why)

| Excluded | Reason |
|---|---|
| **BLOSUM on CDRs** | Sequence-substitution scoring is blind to surface-driven cross-reactivity (7% identity cases). Kept only as a naive baseline the tool visibly beats. |
| **IEDB linear B-cell epitopes as the reference** | Linear peptides can't represent conformational CDR off-targets. Reference set is folded self-protein surfaces instead. |
| **NetMHCpan / APE-Gen / PepSim** | HLA/pMHC tools — TCR pipeline, irrelevant to antibodies; PepSim also web-only. |
| **AlphaFold3 (for cofolding)** | Weights gated; Boltz-2 is MIT-licensed, pip-installable, gives an affinity score, and you hold credits. (AF3 is marginally better on antibody–antigen accuracy — noted, not blocking.) |
| **Whole-proteome scan** | Floods with false positives. |

---

## 15h budget

| Hours | Focus |
|---|---|
| 0–1 | Verify: ANARCI → CDR-H3 `QLYYFDYW`; AntiBERTy/AbLang2 embeds SHR-1210; ESM-2 embeds VEGFR2; Boltz-2 installs + folds a toy dimer. Fix failures first. |
| 1–6 | Steps 1–5 (data, CDR, both embeddings, ranking). Bio: reference set + assay map. MLE: embedding + ranking. |
| 6–9 | Steps 6–7 (panel rank + agent). End-to-end anchor run. **Freeze MVP.** |
| 9–13 | Phase 2: IgFold+surface on top flags, then Boltz-2 cofold WT+VEGFR2 vs. germlined mutant+VEGFR2. |
| 13–15 | Demo polish; the cofold interface-collapse visual is the closer. |
