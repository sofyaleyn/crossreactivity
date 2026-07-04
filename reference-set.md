# Reference set — assembly instructions

This document is an executable brief for the agent that assembles CrossFlag's self-protein reference set. The output is `data/curated/self_proteins.csv` — a three-layer set consumed by the embedding pipeline (`src/crossflag/reference/build_set.py`).

Read `docs/mvp-spec.md` §Data for how this set is *used*; this file is only about how it's *built*.

---

## Why three layers

A flat surfaceome tells you a candidate *has* an off-target. A layered set tells you *how much to worry*. The layers encode clinical severity:

- **Layer A (base breadth)** — every human cell-surface protein an antibody could physically reach in vivo. Establishes the denominator and gives the demo its statistical weight.
- **Layer B (autoimmune focus)** — known human autoantigens. A hit here is a specificity concern *and* a plausible autoimmune-adverse-event story.
- **Layer C (mimicry-primed)** — self-proteins with documented pathogen mimicry (EBV, HHV-6, coxsackievirus, Campylobacter, etc.). A hit here is highest-value: chronic infections in the patient population mean pre-existing cross-reactive antibody memory may already exist, so the candidate lands into pre-primed territory.

Layers overlap — a protein can be in several. Each protein carries **layer membership flags**, not a hard filter; scoring uses these as weights (see §5).

This layering is CrossFlag's genuine wedge over ARDitox / Tope-seq: they tell you an off-target exists; CrossFlag tells you it's a known autoantigen for Sjögren's *and* an EBV-mimicry target — the patient population's chronic viral background will make it worse.

---

## Target output schema

`data/curated/self_proteins.csv`:

```
protein_id,uniprot_id,gene_symbol,name,sequence,
pdb_or_af_ref,surface_region,
layer_A,layer_B,layer_C,
autoimmune_conditions,mimicry_pathogens,mimicry_epitope,
is_anchor_offtarget,source
```

- `layer_A/B/C` — booleans.
- `autoimmune_conditions` — semicolon-separated (e.g. `SLE;Sjögren`).
- `mimicry_pathogens` — semicolon-separated (e.g. `EBV;HHV-6`).
- `mimicry_epitope` — the shared motif/pentapeptide when known (e.g. `PPPGRRP` for EBNA-1 ↔ Sm B/B').
- `is_anchor_offtarget` — `True` for VEGFR2, FZD5, ULBP2 (SHR-1210 known off-targets); used **only for evaluation**, never fed to the scorer.

---

## Layer A — Base breadth (SURFY + CSPA)

**Target size:** ~2,700 proteins.

### Sources

- **SURFY** — in silico human surfaceome, machine-learning predictor trained on CSPA. Bausch-Fluck et al., PNAS 2018. Contains **2,756 α-helical transmembrane + 130 GPI-anchored** cell-surface proteins at 5% FPR. This is the primary Layer A source.
- **CSPA** — Cell Surface Protein Atlas, mass-spec–validated. Bausch-Fluck et al., PLOS ONE 2015. **1,492 human surfaceome proteins.** Use as a validation subset — flag members with `cspa_validated=True` in the schema (extra column, optional).

### Access

- SURFY predictions: PNAS 2018 supplementary tables (Dataset S1) — direct download.
- CSPA: `https://wlab.ethz.ch/cspa/` → `data/S2_File.xlsx` (a copy is `pone.0121314.s003.xlsx` from the PLOS paper).
- Convenience wrapper (both datasets pre-packaged): `github.com/steveneschrich/surfaceome` (R package; the underlying data files are directly usable from Python).

### Assembly steps

1. Download SURFY predictions at 5% FPR → list of UniProt IDs.
2. For each UniProt ID: fetch the canonical protein sequence via UniProt REST (`https://rest.uniprot.org/uniprotkb/{id}.fasta`).
3. Add CSPA-validated subset flag if present.
4. Set `layer_A=True` for every entry.
5. Optionally, if an AlphaFold structure is available (AF2 covers ~all human proteins at `https://alphafold.ebi.ac.uk/`), record the URL/ID in `pdb_or_af_ref` and annotate `surface_region` (residue indices of extracellular / solvent-exposed regions from UniProt topology annotations). This is polish — not required for the MVP but improves antigen embeddings.

### Sanity checks

- VEGFR2 (`KDR`, UniProt `P35968`), FZD5 (`Q13467`), and ULBP2 (`Q9BZM5`) must all be present. If any is missing, add manually and log.
- Total row count: 2,500–3,000. If far outside, the SURFY threshold or the parse is wrong.

---

## Layer B — Autoimmune focus (AAgAtlas)

**Target size:** ~1,000–1,500 (v1); optionally expand toward the recent ~8,000-entry release if time permits.

### Source

- **AAgAtlas 1.0** — human autoantigen database, text-mined + manually curated. Wang et al., NAR 2017.
- v1 contains **1,126 autoantigens across 1,071 human diseases.**
- A recent expansion referenced in the 2026 mimicry literature contains **~8,045 human autoantigens** (accessed via the updated portal).

### Access

- Portal: `http://biokb.ncpsb.org/aagatlas/` (v1) and `http://biokb.ncpsb.org.cn/aagatlas_portal/` (recent expansion).
- Browse & Download tab → downloadable table (gene symbol, associated diseases, evidence).

### Assembly steps

1. Download the AAgAtlas gene/disease table.
2. Resolve each gene symbol → UniProt ID → canonical sequence (UniProt REST).
3. For each entry: set `layer_B=True`; populate `autoimmune_conditions` (semicolon-separated disease list from AAgAtlas).
4. **Merge with Layer A:** if the same UniProt ID is already in Layer A, add the `layer_B=True` flag and populate the conditions field on the existing row rather than creating a duplicate. Not every autoantigen is a cell-surface protein — intracellular autoantigens are perfectly valid entries; keep them, just with `layer_A=False`.
5. Prefer the v1 release (~1,126) for the MVP; upgrade to the recent expansion only if the v1 assembly completes before hour 4 and you have parse-time budget.

### Sanity checks

- Common lupus autoantigens like Ro/SSA (`TROVE2`), La/SSB (`SSB`), Sm B/B' (`SNRPB`) present.
- Common MS/CNS autoantigens like MBP (`MBP`), MOG (`MOG`) present.
- Row count 1,000+; if under 500, the parse dropped rows silently.

---

## Layer C — Mimicry-primed (hand-curated from mimicry literature)

**Target size:** ~50–200 self-proteins across 4–6 pathogen contexts. Small, hand-curated, high-signal.

### Why hand-curated (not a database)

No single clean database exists for "self-proteins with documented pathogen mimicry." Attempting to auto-build one from IEDB + BLAST would be a project in itself. For a 15h build, a curated set from a handful of well-cited review papers is the honest scope: high-signal, small, and defensible in the demo.

### Source papers (start here — expand if time)

1. **Suliman et al., Immunity Inflammation & Disease 2024** — clinical review of molecular-mimicry–induced autoimmunity across MS, T1D, RA, SLE, GBS, PBC, autoimmune myocarditis. Gives the top mimicry pairs across major autoimmune diseases.
2. **Almulla et al., bioRxiv 2025 / Immunity Inflammation & Disease 2026** — in silico EBV + HHV-6 vs. CNS proteins; identifies **91 mimicry pentapeptides** with 10 CNS proteins (synapsin-1 = 13 pentapeptides; MAG = 12; MBP = 9; MOG = 5). This is the cleanest single source for the EBV/HHV-6 → CNS layer.
3. **Poole et al., Frontiers Immunology 2020** — EBNA-1 ↔ SLE autoantigens (Sm B/B', Ro/SSA, dsDNA). The `PPPGRRP` ↔ `PPPGMRPP` motif is the canonical structural mimicry example.
4. **NCBI Autoimmunity chapter (NBK459460)** — molecular mimicry in autoimmunity + vaccinations; covers the EBV ↔ Ro169–180 case.
5. **Recent EBV → Sjögren's papers** on EBNA-1 ↔ La/SSB (X2AX6PG motif).

### Pathogen contexts to cover (priority order)

1. **EBV (Epstein-Barr virus)** — best-documented. Self-proteins: Sm B/B', Sm D1/D2/D3, Ro/SSA (60 kDa), La/SSB, RNP A, MBP, MOG, MAG, synapsin-1, GlialCAM, α-B-crystallin, dsDNA-binding proteins. Diseases: SLE, Sjögren's, MS.
2. **HHV-6 (human herpesvirus 6)** — CNS proteins (synapsin, myelin proteins). Diseases: MS, CFS/ME, Long COVID.
3. **Coxsackievirus B** — GAD65, cardiac myosin. Diseases: T1D, autoimmune myocarditis.
4. **Campylobacter jejuni** — gangliosides GM1/GD1a (glycan mimicry, not protein — record but flag `mimicry_epitope_type=glycan`; a protein-embedding pipeline can't score glycan mimicry, but keeping the entry documents the honest coverage gap). Disease: GBS (Guillain-Barré syndrome).
5. **Mycobacterium tuberculosis, H. influenzae, A. baumannii** — CNS antigens; secondary priority.
6. **Group A Streptococcus** — cardiac myosin, N-acetyl-β-D-glucosamine. Disease: rheumatic fever / rheumatic heart disease.

### Assembly steps

1. From each source paper, extract every named self-protein with a documented mimicry pair.
2. For each entry, record:
   - Gene symbol → UniProt ID → sequence.
   - `mimicry_pathogens` (semicolon-separated).
   - `mimicry_epitope` — the shared motif when the paper names it (e.g. `PPPGRRP`).
   - `autoimmune_conditions` — condition(s) driven by that mimicry.
3. Set `layer_C=True`. Merge with existing rows (Layer A / Layer B) rather than duplicating.
4. Record the source paper in the `source` field per entry — this is a hand-curated set and provenance matters in the demo.
5. Cross-check against AAgAtlas: most Layer C entries should also be `layer_B=True` (they're by definition autoantigens). Discrepancies are worth a note.

### Sanity checks

- EBNA-1 ↔ Sm B/B' present with motif `PPPGRRP`/`PPPGMRPP`.
- MBP, MOG, MAG, synapsin-1 all present with EBV/HHV-6 mimicry.
- GAD65 (`GAD2`) present with coxsackievirus mimicry.
- Cardiac myosin (`MYH6`/`MYH7`) present with coxsackievirus and Group A Strep mimicry.

### Honest coverage note (put this in the demo caveats)

Glycan mimicry (Campylobacter → GM1) is *not* representable in a protein-embedding pipeline. Log it as a known gap. Everything else in Layer C is protein-protein and fully representable.

---

## Background set — `data/background/benign_proteins.csv`

For score calibration. Sample from SURFY entries that are:
- Layer A only (not in B or C).
- Not implicated in any Norden et al. / MPA polyspecificity report.
- Abundant housekeeping / transport / adhesion proteins (Na/K-ATPase subunits, GLUT transporters, integrins, tetraspanins, etc.).

Target size: ~200–500. Same schema (all layer flags `False` except `layer_A`). Used to compute the enrichment normalization in `rank/embedding_rank.py`.

---

## Anchor injection (do not skip)

The anchor off-targets **must** be verifiable in the final set:

```
VEGFR2 / KDR   UniProt P35968   is_anchor_offtarget=True
FZD5           UniProt Q13467   is_anchor_offtarget=True
ULBP2          UniProt Q9BZM5   is_anchor_offtarget=True
```

If SURFY doesn't include one of these, add it manually and log. The `is_anchor_offtarget` flag is **for evaluation only** — never expose it to the scorer. It's the ground truth `test_panel_rank.py` checks against.

---

## Scoring impact — layer priority weights

The layers feed into `rank/embedding_rank.py` via a priority multiplier:

```
priority(p)  =  1.0
              + w_B if p.layer_B
              + w_C if p.layer_C

final_score(v, p)  =  embedding_similarity(v, p)  ×  priority(p)
```

Starting values: `w_B = 0.5–1.0`, `w_C = 1.0–2.0`. Tune on the anchor so VEGFR2/FZD5/ULBP2 still land in the top-N of Beat 1 (they're all Layer A; the weights should nudge, not overpower, the base similarity). The LLM agent reports the layer membership of each flagged protein.

---

## Order of operations (agent build order)

1. Layer A (SURFY → UniProt sequences) — the base.
2. Anchor injection — verify VEGFR2/FZD5/ULBP2 present; add if missing.
3. Background set assembled from Layer A.
4. Layer B (AAgAtlas → merge into Layer A rows).
5. Layer C (hand-curated from source papers → merge into existing rows).
6. Sanity checks (each layer's checklist above).
7. Freeze the CSV; embed everything (ESM-2) into the cached index.

If time is tight, ship Layers A + B + anchor + background. Layer C is the differentiator but not required for the MVP acceptance test (Beat 1). Adding it upgrades Beat 1.5 in the demo.

---

## Time budget

| Layer | Est. time | Notes |
|---|---|---|
| A | ~1.5h | download + UniProt sequence fetch (batchable) |
| Anchor injection + background | ~0.5h | small; check-driven |
| B | ~1h | AAgAtlas download + merge |
| C | ~2–3h | hand extraction from 3–5 review papers |
| Embedding index (ESM-2) | ~1h | batched inference on ~3k sequences |
| **Total** | **~6–7h** | fits in the hour 1–7 window of the 15h plan |
