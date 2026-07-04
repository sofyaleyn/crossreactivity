# Extensions spec (Phase 2–3) — the confirmation ladder

> **⚠️ SUPERSEDED (2026-07-04).** The Boltz-2 cofold "rung 3" described here is now the **primary, validated** off-target signal — not an extension. The dMaSIF/surface rung and the cheap pre-filter framing are dropped (`findings.md` shows no cheap filter works). **Build from [`plan.md`](plan.md)**; evidence in [`findings.md`](findings.md). Kept for the cofold/interface-metric rationale.

Add only after the MVP (`mvp-spec.md`) is frozen. The MVP ranks; this phase *confirms* the top suspects with two independent, increasingly expensive signals — 3D surface comparison, then cofolding. Neither runs on the whole set: rung 2 on the top ~10 flags, rung 3 on the top ~5 pairs.

Abbreviations: `glossary.md`.

---

## The ladder

```
  [MVP] embedding rank        → all variants        (cheap)
  [8]   surface fingerprint   → top ~10 flags        (medium, GPU or CPU-fallback)
  [9]   Boltz-2 cofolding     → top ~5 pairs         (expensive, GPU)
```

Each rung answers a sharper question and independently corroborates the last:
*embedding: looks similar → surface: presents a matching patch → cofold: actually docks.*

---

## [8] Rung 2 — IgFold + surface fingerprint — `src/crossflag/structure/fold.py`, `surface.py`

For the top-flagged variants:

1. **IgFold** predicts the 3D antibody structure from VH + VL (`from igfold import IgFoldRunner`). CPU-adequate for a handful; trained on antibody structures so CDR-H3 geometry is better than generic folders. Check SAbDab for any deposited camrelizumab crystal — prefer it over prediction.
2. **Surface extraction:** compute the molecular surface (MSMS), take the patch within ~4 Å of CDR-H3/L3 atoms, annotate with electrostatics (APBS), hydrophobicity (Kyte-Doolittle), H-bond geometry.
3. **Fingerprint + compare (dMaSIF):** fingerprint the CDR surface patch and the flagged self-protein's surface patch; cosine-compare.

**Fallback (no GPU):** skip dMaSIF; compare the APBS/MSMS physicochemical patches directly (dot-product over surface channels). Weaker but CPU-only and fully programmatic.

**Output:** structural corroboration score per top-flagged pair; a variant that flags by embedding *and* surface is a stronger suspect.

## [9] Rung 3 — Boltz-2 cofolding gate — `src/crossflag/structure/cofold.py`

The confirmation step, and the demo money-shot. For the top ~5 variant×self-protein pairs, cofold the antibody Fv against the flagged self-protein and read out whether they dock.

- **Tool:** Boltz-2 (`pip install boltz`; `boltz predict input.yaml --use_msa_server`). MIT-licensed (commercial-ok), open weights, jointly predicts structure **and** an affinity score.
- **Inputs:** antibody VH+VL + self-protein sequence (or the flagged domain), as a Boltz YAML with an affinity request.
- **Read-outs:**
  - **Interface confidence** — ipTM / PAE restricted to CDR–antigen contacts (a confident interface *at the CDRs* is the signal, not global pTM).
  - **Affinity score** — `affinity_pred_value` / `affinity_probability_binary` from Boltz-2's affinity head.
- **Decision:** a top suspect that cofolds with a confident CDR interface and favorable affinity is a confirmed flag; one that doesn't survives as "flagged by similarity, not confirmed by docking" — triage, not verdict.

**The anchor demo:**
- Cofold **WT SHR-1210 + VEGFR2** → expect a confident interface at CDR-H3.
- Cofold a **CDR-germlined mutant + VEGFR2** → expect the interface to collapse.
That contrast, from sequence alone, is the entire thesis in one figure.

---

## Difficulties (read before committing to Phase 2)

| Difficulty | Detail | Mitigation |
|---|---|---|
| **Antibody–antigen is the hardest cofolding case** | All co-folders (Boltz-2, AF3) are weaker here than on generic complexes; AF3 slightly edges Boltz-2 on Ab–Ag accuracy. | Treat cofold as confirmation, not ground truth. Use CDR-restricted interface confidence, not global scores. Report it as evidence, never certainty. |
| **Boltz-2 needs a GPU** | Minutes per complex on GPU; far slower on CPU. | This is *why* it's a top-k gate, not a screen. Cap at ~5 pairs. |
| **MSA step** | `--use_msa_server` calls a remote MSA server; offline needs local MSAs. | For a handful of pairs the MSA server is fine; note it's the one non-local call in the rung. |
| **dMaSIF: GPU + CC-BY-NC-ND** | Non-commercial, no-derivatives; GPU-bound. | Use the APBS/MSMS fallback if no GPU or if commercial cleanliness matters. Boltz-2 (MIT) carries the commercial-safe confirmation regardless. |
| **IgFold CDR-H3 uncertainty** | Long H3 loops are the least reliable region. | Use per-residue confidence; low-confidence H3 → lean harder on the cofold rung. |

---

## Phase 3 / stretch

| Extension | Tool | Note |
|---|---|---|
| **ABT-736 second anchor** | — | Anti-β-amyloid antibody discontinued in NHP for off-target binding to PF4 (platelet factor 4). Named antibody + self-protein; sequences may be in patents. |
| **AbodyBuilder3 cross-check** | pip | Alternative antibody structure; compare CDR-H3 conformation vs. IgFold. |
| **Panel × self-protein heatmap** | matplotlib | Rows = variants, cols = self-proteins, cells = ladder-combined risk. Clean closing visual. |
| **Broader reference set** | SAbDab + PDB | Add experimentally resolved antibody–self-protein complexes as positive controls. |

---

## Gating

| Condition | Action |
|---|---|
| MVP not frozen by hour 9 | Don't start Phase 2. Polish the ranked-panel demo. |
| MVP frozen, GPU available | Rung 2 on top flags, then rung 3 (Boltz-2) on WT+VEGFR2 and one germlined mutant. |
| MVP frozen, no GPU | Skip dMaSIF (APBS fallback); still attempt Boltz-2 if any GPU is reachable — the cofold contrast is worth more than rung 2. |
| Phase 2 done by hour 13 | One stretch item (heatmap or ABT-736). |
| Uncertain | A polished ranked panel + one clean cofold contrast beats a broad half-working ladder. |
