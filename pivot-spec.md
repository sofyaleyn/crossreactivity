# Pivot spec — replacing rung-1 (embedding similarity → surface complementarity funnel)

**Status:** proposed redesign, supersedes the rung-1 mechanism in `mvp-spec.md` §4 and reorders the ladder in `extensions-spec.md`. Written after a viability assessment (2026-07-04) found the original rung-1 broken.

---

## Why rung-1 changed

The original rung-1 ranked variants by `cosine(paratope_v, antigen_i)` between an antibody-PLM embedding (AntiBERTy 512-d / AbLang2 480-d) and an ESM-2 antigen embedding (1280-d). This is **broken on three independent grounds**:

1. **Not computable as written** — cosine is undefined between vectors of different dimensionality in separately-trained spaces; no projection/alignment was specified. Even after projection, two independently-trained latent spaces share no coordinate frame, so the cosine has no defined meaning without supervised alignment on binding labels (which we don't have).
2. **Biophysically backwards** — binding is *complementarity* (shape/charge anti-correlation between paratope and epitope), not *similarity*. Embedding similarity, where it means anything, indicates two molecules bind the **same** target (epitope binning), not each other.
3. **Empirically non-specific** — the charitable single-space version (embed both sides with the *same* ESM-2 model, cosine well-defined) was tested on the SHR-1210 anchor against 3 known off-targets + 13 decoys (incl. 5 hard membrane-receptor decoys):

   | Formulation | VEGFR2 | FZD5 | ULBP2 | All 3 beat the hard receptor decoys? |
   |---|---|---|---|---|
   | whole-Fv | #1 / 16 | #8 | #9 | **No** |
   | paratope (CDR-H3+L3) | #1 / 16 | #10 | #4 | **No** |

   VEGFR2 ranks #1 only because it is the largest Ig-fold receptor in the set — a fold/size confound. EGFR, HER2, INSR, IGF1R, PDGFRA (proteins the antibody does *not* bind) rank **above** the real off-targets FZD5/ULBP2. In a realistic 200-protein set, FZD5 (~50th percentile) would fail Beat-1's top-20 requirement.

**Conclusion:** the confirmation-ladder concept and the SHR-1210 anchor are sound; the rung-1 *mechanism* must be replaced. No tweak fixes a modeling-target mismatch.

---

## The core constraint

Boltz-2 cofolding (the only signal that captures complementarity directly) is too expensive to run against every self-protein. We need a **cheap, high-recall pre-filter** to shrink the candidate set before cofolding — *without* dropping the low-homology off-targets the project exists to catch (documented cross-reactivity at as low as **7% sequence identity**; README line 17).

**Any similarity-based filter loses recall on exactly those cases.** So the pre-filter must be on the same axis as the confirmer — surface/interface complementarity — just cheaper.

### Filters considered and rejected as the primary gate

| Idea | Verdict | Reason |
|---|---|---|
| Cross-PLM / same-model embedding cosine | ❌ | Measures fold/similarity, not binding; loses low-homology hits (shown empirically). |
| Restrict to self-proteins similar to the on-target (PD-1) | ❌ primary, ✅ supplementary | Molecular-mimicry cases only; the 3 real SHR-1210 off-targets are not PD-1-like and would be filtered out. Keep as a high-precision side channel. |
| ESM antigen embedding as a binding ranker | ❌ | Same as above. **Legit role:** cluster/de-duplicate the reference set (coverage, not ranking). |

---

## The redesigned funnel

```
Stage 0  Curate + cluster the self-protein reference set        RECALL via prior
         - bias to classes where CDR off-targets occur:
           membrane receptors, surface, secreted (MPA / Norden 2024)
         - ESM-embed → cluster → keep diverse representatives
         - self-protein structures come FREE from AlphaFold DB
              │   ~hundreds of proteins, non-redundant
              ▼
Stage 1  Surface-complementarity pre-filter (was rung-2)        CHEAP, RECALL-tuned
         - fold the antibody ONCE (IgFold) → paratope surface patch
           (shape + APBS electrostatics + Kyte-Doolittle hydrophobicity)
         - match paratope patch vs each self-protein surface patch
           via MaSIF-search / dMaSIF (complementary fingerprints)
         - CPU / commercial-safe fallback: APBS+MSMS physicochemical
           patch dot-product
              │   keep top ~20–40% (cast wide; this is a recall gate)
              ▼
Stage 2  Boltz-2 cofold + affinity (rung 3, unchanged)          EXPENSIVE, PRECISION
         - cofold survivors; CDR-restricted ipTM/PAE + affinity
         - only ~10–20 pairs
```

**Panel triage, in parallel (antibody-only, no antigen named):** rank variants by validated physicochemical CDR polyreactivity descriptors — CDR-H3/L3 net charge, hydrophobicity, hydrophobic/charge patch area, aromatic (Trp/Tyr/Arg) content (TAP / CamSol style). This is the well-defined, unsupervised, hours-of-work replacement for "which variant is stickiest," separate from the antigen-specific "does WT hit VEGFR2" question, which Stages 1–2 answer.

### Why this preserves recall where similarity did not

- Self-protein structures are **pre-computed (AlphaFold DB, free)**, so Stage 1 costs only surface extraction + matching — no folding. Affordable across the *entire* curated set, so nothing is dropped by cost before the complementarity check.
- Surface complementarity **is** the binding mechanism, so a globally-unrelated protein with one matching surface patch survives — the 7%-identity case that sequence/embedding filters lose by construction.
- MaSIF-search is purpose-built to scan for interaction partners by complementary surface fingerprints; dMaSIF is the fast variant.

### Honest recall ceiling

Recall is ultimately bounded by **Stage 0** — you can only find off-targets that are in the curated set. Prior-based curation (the right protein classes) matters more than filter cleverness, and it is the project's real limitation to state on stage.

---

## Compute reality (this environment)

- **No local CUDA GPU** (Apple M2, Metal only). Boltz-2 (Stage 2) and dMaSIF (Stage 1 primary) need a GPU → use a cloud GPU (NVIDIA Boltz-2 NIM, or a rented A100) for Stage 2; use the **APBS/MSMS CPU fallback** for Stage 1 to stay local + commercial-safe.
- IgFold (fold the antibody once), APBS/MSMS, ESM clustering, and the physicochemical descriptors all run on CPU.

---

## Empirical de-risking — RUN 2026-07-04 (boltz.bio API, Boltz-2.1)

Full experimental log in **`findings.md`**. Headline results that shape this design:

**Cofold confirmation (Stage 2) works for most off-targets.** Calibrated panel (5 samples each) anchored by a true-target ceiling (PD-1) and non-binder floor (lysozyme), read on **PAE-mean + epitope reproducibility** (ipTM is uninformative — Boltz over-docks, lysozyme scores 0.87):

| complex | PAE-mean | epitope reprod | read |
|---|---|---|---|
| WT × PD-1 (cognate, ceiling) | 6.4 | 0.94 | binder |
| WT × FZD5 (off-target) | 4.9 | 0.66 | **confirmed** |
| WT × ULBP2 (off-target) | 5.1 | 0.89 | **confirmed** |
| WT × VEGFR2 D2-3 (off-target) | 10.5 | 0.43 | floor / unconfirmed |
| WT × lysozyme (floor) | 11.6 | 0.45 | non-binder |

**2 of 3 known off-targets confirm** at/near the ceiling; VEGFR2 is the outlier (lowest-affinity, multi-domain, epitope/domain unmapped — templating it to its real 3V2A structure made it *worse*, 0.43→0.30). So the Stage-2 confirmer is real, not the dead end an earlier VEGFR2-only reading suggested.

**Cheap embedding shortlisting (Stage-1 candidate) fails the confound-controlled test.** On the fold-matched decoy benchmark (`findings.md §Exp 5`), same-model ESM ranking gives **within-family AUROC 0.58 ≈ chance** (VEGFR2 0.79, FZD5 0.55, ULBP2 0.42). Global rank looks good only via the fold confound. **Embedding cannot be Stage 1.**

**Consequence for this design:** Stage 2 (cofold) is validated as the discriminator, but it can't run proteome-wide, so the funnel lives or dies on finding a **cheap Stage-1 that enriches within-fold**. The surface-complementarity proposal below is now the load-bearing untested piece. Two checks in flight: (1) does cofold *discriminate* within-fold (FZD5 vs other Frizzled CRDs); (2) does a surface method beat AUROC 0.58 on the benchmark.

---

## Anchor data status (resolved 2026-07-04)

| Item | Status |
|---|---|
| WT SHR-1210 VH/VL | ✅ Real, confirmed (matches Finlay et al. mAbs 2019, Table 1) |
| VEGFR2 / FZD5 / ULBP2 antigens | ✅ Real (UniProt P35968 / Q13467 / Q9BZM5) |
| Germlined-mutant panel (gates Beat 2) | ❌ **Not published — must be reconstructed** |

Finlay et al. (PMID 30541416) name ~15 engineered variants but publish **no full VH/VL sequences, no residue-by-residue mutation table, no accession, and no patent** for them. Only a directional recipe is given: the light chain was near-germlined to human **IGKV1-39 + JK1** (replacing JK4 to place a Trp in CDR-L3), with "almost the entire CDR-L2" revertible to germline. The lone patent hit (WO2021069670A1) is a deliberate PD1×VEGFR2 *dual-binder* — not the specificity-enhanced variants.

Implications:
- **Beat 2 is an illustration, not a result** — label all reconstructed variants clearly.
- **In the original embedding design, reconstruction is near-circular**: germline-reverting the light-chain CDRs guarantees a drop in any CDR-keyed similarity score → tests plumbing, not biology.
- **In the pivoted design this is largely absorbed**: the meaningful contrast becomes WT-vs-germlined in the **Boltz-2 cofold** against VEGFR2, which models the real interface rather than a CDR-keyed similarity — a legitimate, non-circular test. The IGKV1-39/JK1 recipe is specific enough to reconstruct reproducibly.
