# Demo run — MVP acceptance test

> **⚠️ SUPERSEDED (2026-07-04).** These beats assume the embedding-ranking MVP. The validated demo is now a Boltz-2 cofold screen ranked on a calibrated panel — see [`findings.md`](findings.md) (§Next steps) and [`pivot-spec.md`](pivot-spec.md). Kept for history.

The exact sequence of runs that proves the MVP works, using the SHR-1210 anchor. Three beats, each producing a concrete artifact you show on stage. If Beat 1 fails, the MVP is broken and nothing else matters.

Prereqs: MVP built through `pipeline.py`; curated self-protein reference set assembled (see notes on size below); variant panel FASTAs prepared; ESM-2 embeddings of the reference set cached.

---

## Beat 1 — Retrospective validation (the anchor)

**The question this answers:** "Would the tool have flagged VEGFR2, FZD5, and ULBP2 as SHR-1210 off-targets, before the clinic?"

### Run

```
Input:  data/anchor/shr1210_vh.fasta + shr1210_vl.fasta   (wild-type only)
Ref:    data/reference/self_proteins.csv                    (see size requirement below)
Bg:     data/reference/background/benign_proteins.csv
```

Pipeline: CDR extraction → AntiBERTy/AbLang2 paratope embedding → cosine similarity vs. every self-protein in the reference set → background-calibrated risk → ranked flag list.

### What Beat 1 is actually measuring (read before you present it)

The paratope vector (AntiBERTy/AbLang2) and the antigen vector (ESM-2) come from two separately-trained models, and binding is *complementarity*, not resemblance — so this cosine is **not** a binding prediction, and we never present it as one. What makes it a usable signal is narrower and honest:

Both models are trained on protein sequences, so both encode the same underlying biophysical vocabulary — hydrophobicity, charge, aromatic content, local structural propensity. Two very different objects can still be described in that shared chemistry. Calibrated against the benign background (see [mvp-spec.md](mvp-spec.md) §[4]), the cross-model cosine is therefore a **fuzzy proxy** for a single claim:

> *"This paratope's chemistry is unusually enriched for the kind of chemistry present on this antigen's surface, relative to proteins we know are safe."*

That is a **smoke detector, not a binding model.** It is cheap, it runs on the entire panel in one pass, and it is *allowed to be wrong* — because every suspect it surfaces is re-checked by the surface and cofolding rungs downstream. Beat 1 earns a ranked shortlist; Beats 2–3 earn the verdict.

**The version that would make sense as a true binding signal** — a possible upgrade, not the MVP — is a nearest-neighbor search over *known antibody–self-protein binding pairs* rather than over bare self-protein sequences: *similar paratopes bind similar things.* That comparison is same-model (no cross-space mismatch) and it respects complementarity by transfer — the known binder already demonstrated the fit, so we're only flagging lookalikes. It needs a curated set of known binders we don't have in the MVP window (see [extensions-spec.md](extensions-spec.md) Phase 3).

### Pass criteria (all must hold)

1. **All three known off-targets appear in the top-N flags:**
   - VEGFR2 in top-10 (target: top-5).
   - FZD5 in top-20.
   - ULBP2 in top-20.
2. **All three score above the background 95th percentile** (they must look meaningfully enriched, not just be in the list).
3. **At least one of VEGFR2/FZD5/ULBP2 is in the top-3.**

### Reference set size requirement

For this beat to be *meaningful*, the reference set must be large enough that surfacing three specific proteins near the top is not trivial:
- **Minimum: 200 self-proteins** in the curated set, biased toward surface-exposed / membrane / secreted proteins where CDR off-targets are plausible.
- Target: 500+ if time permits.
- Anchor off-targets flagged in the CSV as `is_anchor_offtarget=True` — but that flag is used only for evaluation (compute rank), never fed into the scorer.

If the reference set has only ~20 proteins, surfacing three of them means nothing. Do not skimp here.

### Artifact for the demo

A ranked table, top ~15 rows, with VEGFR2/FZD5/ULBP2 highlighted:

```
rank  protein       risk_score   bg_percentile
1     VEGFR2   ★    0.87         99.8
2     KDR-related   0.71         98.1
3     FZD5     ★    0.68         97.4
...
9     ULBP2    ★    0.54         95.6
```

Story line: "SHR-1210 killed patients via VEGFR2. From sequence alone, the tool flags VEGFR2, FZD5, and ULBP2 — all three known off-targets — near the top of a set of 200+ self-proteins."

### If it fails

Do NOT proceed to Beat 2. Diagnose in this order:
- Are VEGFR2/FZD5/ULBP2 embedded and in the index? (Print the ref index; check.)
- Is CDR pooling correct? (Print CDR-H3, confirm `QLYYFDYW`; check pooling weights.)
- Is the background calibration inverted? (Check that higher score = more suspect.)
- Is the reference set too small? (Add proteins, retest.)
- Try the other antibody PLM (AntiBERTy ↔ AbLang2).

---

## Beat 2 — Panel optimisation (the edge)

**The question this answers:** "Does the tool correctly rank the CDR-germlined mutants as cleaner than wild-type?"

### Run

```
Input:  data/anchor/shr1210_vh.fasta + vl.fasta (WT)
        data/anchor/variants/*.fasta            (CDR-germlining mutants)
Ref:    same as Beat 1
```

Pipeline: run each variant through the same pipeline; produce a panel table with each variant's risk score against VEGFR2 (and overall).

### Pass criteria

1. **WT ranks the riskiest** in the panel by overall risk score.
2. **At least one germlined mutant has a substantially lower VEGFR2-specific score** than WT (target: drop of ≥ 30% of WT's above-background enrichment).
3. The tool names WT as the variant to drop and identifies at least one mutant as "advance."

### Artifact for the demo

Panel ranking table:

```
variant           risk_score   top_flag        VEGFR2_score    verdict
SHR-1210 WT       0.87         VEGFR2          0.87            drop
mutant_L1_v1      0.42         (background)    0.31            advance
mutant_L2_v1      0.55         FZD5            0.48            review
mutant_L3_v1      0.38         (background)    0.28            advance
```

Story line: "The tool ranks the panel: wild-type is dropped, germlined mutants advance. The published fix (Finlay et al. CDR germlining) is what the tool would have recommended."

### Caveat to state on stage

The germlined mutant sequences are (partly) reconstructed if the paper's exact sequences aren't retrievable — label reconstructed variants clearly in the demo table. The wild-type result (Beat 1) is real; the mutant result is illustrative of the workflow.

---

## Beat 3 — Cofolding confirmation (extension, optional)

**The question this answers:** "Not just similar, but actually docks."

### Run

Two Boltz-2 cofolds:
1. **WT SHR-1210 Fv + VEGFR2** (or the VEGFR2 extracellular ligand-binding domain).
2. **A germlined mutant Fv + VEGFR2.**

Each: `boltz predict <yaml> --use_msa_server`, with affinity request.

### Pass criteria

1. WT cofold shows a **confident CDR-restricted interface** (low PAE at CDR-H3 ↔ VEGFR2 contacts; reasonable ipTM for the interface).
2. Germlined mutant cofold shows the **interface collapses** (PAE at CDR contacts substantially worse; or the model docks the antibody elsewhere; or affinity score drops).

### Artifact for the demo

Two side-by-side structure views (PyMOL / Molstar screenshots), with:
- CDR-H3 highlighted.
- The VEGFR2 contact patch highlighted.
- Interface PAE color-mapped.
- Boltz-2 affinity score printed under each.

Story line: "Not just similar in embedding space — Boltz-2 predicts wild-type actually docks against VEGFR2 with high interface confidence. The germlined mutant loses the interface. That's the mechanism, from sequence alone."

### If it fails

- Try cofolding the VH+VL against the **VEGFR2 Ig-domain 2/3** (the ligand-binding region) rather than full-length — full-length is harder.
- Note that antibody–antigen is the weakest cofolding category; low confidence is a known limitation, not a bug. Report the affinity-score contrast even if interface confidence is modest.
- If time-boxed out, drop Beat 3 entirely — Beats 1+2 are sufficient for the MVP.

---

## The full demo script (~5 minutes)

1. **Problem (30s):** Antibody CDR polyspecificity killed patients on SHR-1210 (hemangioma via VEGFR2 off-target). 33% of antibody leads have this. Screening is expensive.
2. **Beat 1 (90s):** Load WT SHR-1210. Show the ranked flag list from 200+ self-proteins. VEGFR2, FZD5, ULBP2 all near the top. "We would have flagged this before the clinic."
3. **Beat 2 (90s):** Load the variant panel. Show the ranking: WT dropped, germlined mutants advanced. "Panel-level triage — we tell you which variant to advance."
4. **Beat 3 (60s, if it works):** Boltz-2 cofold contrast: WT docks VEGFR2, mutant doesn't. "Not similarity — confirmed docking."
5. **The honest frontier (30s):** This is triage, not certification. The output prioritizes wet-lab specificity screens; it does not replace them. Reference set coverage is the current limit.

## Absolute minimum for a passable demo

**Beat 1 alone.** If Beat 1 passes on the anchor with a real 200+ reference set, the MVP works and the demo is credible. Beats 2 and 3 raise the ceiling — but a working Beat 1 is the floor and the success criterion.
