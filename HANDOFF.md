# HANDOFF — build instructions for Claude Code

> **⚠️ SUPERSEDED (2026-07-04).** This build order centers the embedding rung-1, which was found broken. **Build from [`plan.md`](plan.md)** (the active representative-set cofold-screen plan), with rationale in [`pivot-spec.md`](pivot-spec.md) and evidence in [`findings.md`](findings.md). Kept for history.

This file is the executable brief. It turns the specs into an ordered build with explicit module contracts and acceptance tests. Read [docs/mvp-spec.md](docs/mvp-spec.md) and [docs/extensions-spec.md](docs/extensions-spec.md) for rationale; build from here.

**Goal:** an antibody CDR cross-reactivity triage tool, validated on the SHR-1210 → VEGFR2 anchor. Rank a variant panel by antibody-PLM ↔ antigen-PLM embedding similarity (MVP), then confirm top suspects with IgFold surface comparison and Boltz-2 cofolding (extensions).

**Non-negotiable constraints:**
- Every tool must be callable programmatically (no web-form tools). Wrap each as a skill in `skills/`.
- MVP (Phase 1) must run end-to-end and pass its acceptance test before any Phase 2 work.
- Expensive steps are top-k gated: surface on top ~10, cofolding on top ~5. Never run cofolding across the set.
- Framing everywhere is triage, not certification.

---

## Environment

- Python 3.10+. Fresh venv/conda env.
- `pip install -e ".[embed,structure,dev]"` (deps declared in `pyproject.toml`; `embed`/`structure` are the MVP and confirmation-ladder extras).
- `conda install -c bioconda anarci` (HMMER dependency; UNIX/macOS — not pip-installable, so it's not in `pyproject.toml`).
- GPU needed for Boltz-2 (rung 3) and dMaSIF (rung 2 primary). Rung-2 has a CPU fallback (APBS/MSMS). Confirm GPU early.
- Pin the exact resolved env into `requirements.lock.txt` (`pip freeze > requirements.lock.txt`) as you install.

---

## Build order

### Step 0 — scaffold + hour-zero verification
Create the package skeleton (`src/crossflag/...` per README tree), `pyproject.toml`, `tests/`. Then run the verification harness before building logic:

- ANARCI on SHR-1210 VH → CDR-H3 extracts as `QLYYFDYW`.
- AntiBERTy (or AbLang2) returns an embedding for SHR-1210 VH/VL.
- ESM-2 returns an embedding for the VEGFR2 sequence.
- `boltz predict` runs on a toy two-protein YAML and returns confidence + affinity fields.

If any fail, stop and fix — everything downstream depends on these.

### Step 1 — data (`data/`)
- `anchor/shr1210_vh.fasta`, `shr1210_vl.fasta` (sequences in README).
- `anchor/variants/` — WT + CDR-germlining mutants from Finlay et al. (encode each as a VH/VL FASTA pair; if exact mutant sequences are unavailable, generate light-chain CDR→germline revertants and label them clearly as reconstructed).
- `anchor/offtargets/` — VEGFR2, FZD5, ULBP2 sequences + PDB IDs (VEGFR2: 4ASD/1VR2).
- `reference/self_proteins.csv` — schema in mvp-spec §Data; VEGFR2/FZD5/ULBP2 flagged `is_anchor_offtarget=True`. Built by the `crossflag.reference` pipeline (see Step 5 and [docs/reference-set.md](docs/reference-set.md)); raw downloads land in `reference/raw/` and intermediates in `reference/build/` (both gitignored).
- `reference/background/benign_proteins.csv` — housekeeping/serum surface proteins.
- **Known gap:** Layer B (AAgAtlas autoantigens) is not yet built — `self_proteins.csv` is currently Layer A + the small Layer C set. See `crossflag/reference/layer_b.py`.

### Step 2 — `extract/cdrs.py`
Contract:
```
extract_cdrs(vh: str, vl: str) -> dict
  returns {"cdr_h1":str,"cdr_h2":str,"cdr_h3":str,
           "cdr_l1":str,"cdr_l2":str,"cdr_l3":str,
           "vh_cdr_indices":list[int],"vl_cdr_indices":list[int]}
```
IMGT scheme. Indices are residue positions (0-based into vh/vl) for embedding pooling.

### Step 3 — `embed/antibody.py`
Contract:
```
embed_antibody(vh: str, vl: str, cdr_indices: dict) -> np.ndarray  # paratope vector
```
Run AntiBERTy/AbLang2 → per-residue embeddings; pool CDR residues (mean or attention-weighted, weight CDR-H3 highest) into one vector. Deterministic; cache by sequence hash.

### Step 4 — `embed/antigen.py`
Contract:
```
embed_antigen(sequence: str, surface_region: list[int] | None) -> np.ndarray | list[np.ndarray]
```
ESM-2. If `surface_region` given, pool those residues; else sliding-window sub-vectors over the chain. Cache.

### Step 5 — `reference/build_set.py`
The `crossflag.reference` package assembles the layered self-protein set (see [docs/reference-set.md](docs/reference-set.md)). `build_set.py` orchestrates the layer modules in order:
```
layer_a → background → (layer_b, TODO) → layer_c → merge → embed
```
- `layer_a.py` — SURFY + CSPA base + anchor injection → `reference/build/self_proteins_layer_a.csv`.
- `background.py` — sample benign families from Layer A → `reference/background/benign_proteins.csv`.
- `layer_c.py` — fill sequences into the hand-curated mimicry seed → `reference/build/layer_c_mimicry_seed_filled.csv`.
- `merge.py` — merge layers by UniProt ID → `reference/self_proteins.csv`.
- `paths.py` — single source of truth for all `data/reference/` locations (repo-root-anchored).
Then embed every protein (step 4) and persist an embedding index (`.npz` or FAISS) into `reference/index/`. Same for the background set. Run: `python -m crossflag.reference.build_set`.

### Step 6 — `rank/embedding_rank.py`
Contract:
```
score_variant(paratope_vec, ref_index, background_index) -> dict
  returns {"risk_score":float,"neighbors":[(protein_id, cosine), ...]}
```
Cosine to every reference vector; aggregate (max or top-k mean); calibrate against background distribution (z-score or percentile). Higher = more suspect.

### Step 7 — `rank/panel.py` + `assay/map.py`
Rank all variants by `risk_score`; attach top flagged protein + assay (map table in mvp-spec §5). Emit a panel table + per-variant evidence.

### Step 8 — `agent/report.py`
LLM turns the ranked table + evidence into a recommendation (advance which variant, what it cross-reacts with, which assay, honest caveat). Retrieval + reasoning only.

### Step 9 — `pipeline.py` (MVP end-to-end)
Wire steps 2→8. Input: a directory of variant FASTA pairs. Output: ranked panel JSON + agent report.
**→ Run the MVP acceptance test (below). Freeze when green.**

### Step 10 — rung 2: `structure/fold.py` + `structure/surface.py` (EXT)
IgFold on top-N flagged variants → PDB. MSMS surface → CDR patch (within 4 Å of CDR-H3/L3) → APBS/hydrophobicity annotation → dMaSIF fingerprint + cosine vs. flagged self-protein surface. CPU fallback: APBS/MSMS patch comparison, no dMaSIF.

### Step 11 — rung 3: `structure/cofold.py` (EXT)
For top ~5 variant×self-protein pairs, write Boltz YAML (VH+VL + antigen, affinity request), run `boltz predict --use_msa_server`, parse: CDR-restricted ipTM/PAE + `affinity_pred_value`. Return a confirmation verdict per pair.
**→ Run the cofold acceptance test.**

### Step 12 — skills + notebook
Wrap each callable as a `skills/*.md` skill (input schema, invocation, output parse). Build `notebooks/anchor_validation.ipynb` demoing the full ladder on SHR-1210.

---

## Acceptance tests (`tests/`)

These are the code-level `pytest` contracts. The narrated demo and the canonical numeric
pass thresholds (the three "Beats") live in [docs/demo-run.md](docs/demo-run.md) — that is the
single source of truth for "what passing looks like."

`test_cdr_extract.py`
```
extract_cdrs(SHR1210_VH, SHR1210_VL)["cdr_h3"] == "QLYYFDYW"
```

`test_panel_rank.py` (MVP gate — must pass to freeze Phase 1)
```
run pipeline on {WT + germlined mutants}
assert VEGFR2 in top-3 flagged proteins for WT
assert risk_score(WT) > risk_score(at least one germlined mutant)
```

`test_cofold_gate.py` (Phase 2 gate)
```
cofold(WT, VEGFR2) interface_confidence > cofold(germlined_mutant, VEGFR2) interface_confidence
```

If `test_panel_rank.py` cannot pass on embeddings alone, that is the signal to tune the pooling/aggregation — not to move on.

---

## Definition of done

- **MVP:** `pipeline.py` runs from a variant-FASTA directory to a ranked panel + report; `test_panel_rank.py` green; every MVP tool wrapped as a skill.
- **Demo:** notebook shows WT SHR-1210 flags VEGFR2 (embedding), corroborated by surface (rung 2), confirmed by a confident WT–VEGFR2 cofold interface that collapses for a germlined mutant (rung 3).
- **Honesty:** every output frames results as triage/prioritization; cofold reported as confirmation evidence, not certainty (antibody–antigen is the hardest cofolding class).

## Known risks (carry into build)

- Exact Finlay mutant sequences may need reconstruction — label reconstructed variants clearly.
- Antibody–antigen cofolding is the weakest category for all models — use CDR-restricted interface confidence, never global pTM alone.
- dMaSIF is CC-BY-NC-ND (non-commercial); Boltz-2/ESM-2/IgFold are commercial-safe. Keep the commercial-safe path (APBS fallback + Boltz-2) intact.
- Embedding similarity ranks; it does not prove binding. The ladder is what earns confidence.
