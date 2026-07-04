# MVP spec (Phase 1) — embedding-based panel ranking

The minimum viable build. Rank the antibody variant panel by learned-embedding similarity between the antibody paratope and self-protein antigens. No 3D. This phase meets the success criterion at the embedding level; the confirmation ladder (`extensions-spec.md`) makes it convincing.

Abbreviations: `glossary.md`.

---

## Why embeddings, not sequence matching

CDR cross-reactivity is a surface/chemistry phenomenon — documented off-targets share as little as 7% sequence identity with the on-target. Raw sequence scoring (BLOSUM) is blind to this. Antibody-specific PLMs (protein language models) — AntiBERTy, AbLang2 — produce residue embeddings that align with paratope regions and encode binding-relevant features; general PLMs (ESM-2) do the same for arbitrary proteins. Comparing in embedding space captures similarity that sequence identity misses.

**Honest scope limit:** embedding *similarity* is not a trained *binding predictor*. Supervised antibody–antigen binding models (EPP, MambaAAI) exist but need training data we don't have in 15h. So the MVP uses unsupervised embedding similarity as a *ranking/triage* signal — and the confirmation ladder (surface + cofolding) is what turns a ranked suspicion into evidence. We rank; we don't certify.

---

## Pipeline

```
  Antibody variant panel (VH + VL per variant)
        │
        ▼
  [1] CDR extraction (ANARCI/ANARCII)                 → CDR spans per variant
        │
        ▼
  [2] Antibody embedding (AntiBERTy / AbLang2)        → CDR-pooled paratope vector
        │
  Self-protein reference set
        │
        ▼
  [3] Antigen embedding (ESM-2)                       → per-protein vector(s)
        │
        ▼
  [4] Embedding-similarity ranking + calibration      → risk per variant
        │
        ▼
  [5] Panel ranking + assay map
        │
        ▼
  [6] LLM agent narration
```

---

## [1] CDR extraction — `src/crossflag/extract/cdrs.py`

Run **ANARCI** or **ANARCII** (IMGT scheme) on each variant's VH and VL. Output all six CDR spans. Purpose here is **localization**: the CDR residue indices tell step [2] which residues to pool/weight into the paratope vector. We do not string-compare CDRs.

Verify on the anchor: CDR-H3 = `QLYYFDYW`, CDR-L3 = `QQVYSIPWT`.

**Output:** `{variant_id, vh_seq, vl_seq, cdr_spans}`.

## [2] Antibody embedding — `src/crossflag/embed/antibody.py`

Embed each variant with an **antibody-specific PLM**:
- **AntiBERTy** (`from antiberty import AntiBERTyRunner; antiberty.embed(seqs)`) → per-residue `[(L+2) × 512]`.
- or **AbLang2** (paired VH|VL, `rescoding` per-residue + `seqcoding` per-sequence).

Pool the **CDR residues** (from step [1]) into a single paratope vector — mean or attention-weighted over CDR-H3 + CDR-L3 (weight H3 highest). Optionally concatenate a whole-Fv `seqcoding` for context.

**Output:** `{variant_id, paratope_vector}`.

## [3] Antigen embedding — `src/crossflag/embed/antigen.py`

Embed each self-protein in the reference set with **ESM-2** (`esm2_t33_650M` is a good size/quality balance; smaller if compute-limited).
- If a structure/PDB is known for the self-protein, restrict to the **surface-exposed region** likely to face a CDR (solvent-accessible residues); otherwise embed the full chain and note lower confidence.
- Store per-protein vector(s). For large proteins, sliding-window sub-vectors so a local surface patch can match.

**Output:** `{protein_id, antigen_vector(s), uniprot, pdb_ref}`.

## [4] Embedding-similarity ranking + calibration — `src/crossflag/rank/embedding_rank.py`

For each variant, score similarity (cosine) between its paratope vector and every reference antigen vector. Calibrate against a **background** of benign self-proteins so the score is enrichment, not raw distance.

```
risk_v = aggregate_over_reference( cosine(paratope_v, antigen_i) )
         normalized against the benign-background distribution
```

Aggregate = max, or top-k mean. On the anchor, WT SHR-1210 should score high against VEGFR2; germlined mutants should drop toward background.

**Output:** `{variant_id, risk_score, top_flagged_protein, neighbor_list}`.

## [5] Panel ranking + assay map — `src/crossflag/rank/panel.py` + `src/crossflag/assay/map.py`

Rank variants by `risk_score`; name the cleanest variant and flagged self-protein(s).

Assay map (self-protein / biology → assay):

| Flagged self-protein / biology | Confirmation assay (priority) |
|---|---|
| VEGFR2 / angiogenic receptor | receptor proteome array (MPA / Retrogenix) → VEGFR2 SPR → HUVEC agonism assay |
| FZD5 / Wnt receptor | receptor array → FZD5-Wnt reporter |
| Any membrane receptor | cell-based protein array (MPA) → target-specific binding + function |
| Generic surface protein | SPR on recombinant protein |

MPA = Membrane Proteome Array; SPR = surface plasmon resonance; HUVEC = human umbilical vein endothelial cell.

**Output:** ranked panel + per-variant evidence (flagged protein, nearest-neighbor rationale, top-contributing CDR from embedding attribution).

## [6] LLM agent narration — `src/crossflag/agent/report.py`

LLM reads the ranked evidence, writes: which variant to advance, what each flagged variant likely cross-reacts with, which assay confirms, and the honest caveat (embedding similarity is a triage signal; the confirmation ladder is the evidence). Retrieval + reasoning; no training.

---

## Data (MVP)

### Anchor — `data/anchor/`

| Item | Value | Source |
|---|---|---|
| SHR-1210 VH / VL | see README | Thera-SAbDab / US11208484B2 |
| Off-targets | VEGFR2, FZD5, ULBP2 | Finlay et al. 2019 |
| Adverse event | capillary hemangioma; VEGFR2 agonism confirmed by antagonist rescue | same |
| Variant panel | WT + light-chain CDR-germlining mutants that ablated VEGFR2 binding | same |
| Off-target structures | VEGFR2 e.g. PDB 4ASD / 1VR2 | RCSB |

### Curated self-protein reference set — `data/curated/self_proteins.csv`

Surface-accessible human self-proteins, biased toward membrane receptors and secreted proteins where CDR off-targets are plausible (from MPA literature / Norden et al. 2024). VEGFR2, FZD5, ULBP2 explicitly included. Store sequence + PDB/AlphaFold reference.

```
protein_id | name | uniprot | sequence | pdb_or_af_ref | surface_region | is_anchor_offtarget
```

Scope: hundreds of proteins. Small + curated is the feature.

### Background — `data/background/benign_proteins.csv`

Surface self-proteins with no documented antibody off-target (e.g. abundant housekeeping/serum proteins). For calibration.

---

## Owner split

- **Bio / curation (you):** anchor data, curated reference set + surface-region annotation, background set, variant panel FASTAs, assay map, narrative.
- **MLE (partner):** ANARCI wrapper, both embedding modules, ranking + calibration, panel ranker, LLM agent, demo.
- **Converge:** hour 9, end-to-end anchor run — VEGFR2 top for WT; a germlined mutant cleaner.

## Success criterion

From the SHR-1210 panel (embeddings only): **VEGFR2 near the top for wild-type**; a CDR-germlining mutant lower; assay named.
