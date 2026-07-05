# Tools reference

> **⚠️ TIERING SUPERSEDED (2026-07-05).** The MVP-vs-extension split below reflects the original embedding-ladder design, which was found broken. In the **validated** tool the **core engine is Boltz-2 cofold via the hosted boltz.bio API** (see the "Extension tools" section — that is now the *primary* path, not an add-on), scoring `PAE_IF` + `epitope_reprod`. The antibody/antigen embedding PLMs (AntiBERTy / AbLang2 / ESM-2) and the dMaSIF/surface rung are **no longer the primary path** — the embedding ranker failed a fold-matched benchmark and no cheap filter beats a free annotation baseline (`../findings.md`). The per-tool facts (installs, licenses, APIs) below are still accurate; only the MVP/extension *framing* is stale. Current design: [`../plan.md`](../plan.md); validated results: [`../demo/README.md`](../demo/README.md).

Every tool is programmatic (local or callable API) and skill-wrappable. Web-only tools excluded. **Skill?** = drivable without runtime human interaction.

---

## MVP tools (Phase 1) — *superseded as the primary path; see banner*

### Python + Biopython
- **Role:** parsing, FASTA/PDB handling, glue.
- **Install:** `pip install biopython`
- **Skill?** Yes. **License:** BSD-3.

### ANARCI / ANARCII
- **Full name:** Antibody Numbering and Antigen Receptor ClassIfication (ANARCII = 2025 transformer successor).
- **Role:** CDR annotation of VH/VL (and TCR) with IMGT numbering; localizes the paratope for embedding pooling.
- **Install:** `conda install -c bioconda anarci` (UNIX/macOS, needs HMMER) or `pip install anarcii`.
- **Skill?** Yes (Python API). **License:** BSD-3.

### AntiBERTy
- **Full name:** antibody-specific masked language model (BERT trained on 558M antibody sequences).
- **Role:** per-residue embeddings for antibody variants; CDR residues pooled into a paratope vector.
- **Install:** `pip install antiberty`. API: `from antiberty import AntiBERTyRunner; AntiBERTyRunner().embed(seqs)` → `[(L+2) × 512]` tensors; `return_attention=True` for attention maps.
- **Skill?** Yes. **License:** open (MIT-style). CPU-fine.

### AbLang2 (alternative to AntiBERTy)
- **Full name:** antibody language model (paired VH/VL), germline-bias-corrected, ESM-2-derived (12 layers, 480-d).
- **Role:** paired-chain embeddings (`rescoding` per-residue, `seqcoding` per-sequence).
- **Install:** `pip install ablang2`. API: `ablang2.pretrained(model_to_use='ablang2-paired')`; input `"{VH}|{VL}"`. CPU-fine (ms per pair).
- **Skill?** Yes. **License:** open. **Note:** `align=True` needs ANARCI in the env. Pick AntiBERTy *or* AbLang2 — either suffices; AbLang2 handles paired chains natively.

### ESM-2
- **Full name:** Evolutionary Scale Modeling v2 — general protein language model.
- **Role:** embed self-protein antigens in the reference set.
- **Install:** `pip install fair-esm` or via HuggingFace (`facebook/esm2_t33_650M_UR50D`). Use `t33_650M` for balance; smaller (`t12_35M`) if compute-limited.
- **Skill?** Yes. **License:** MIT. GPU helps for large proteins; CPU workable for a small set.

### LLM (large language model)
- **Role:** agent narration + recommendation.
- **Access:** API. **Skill?** Yes.

---

## Extension tools (Phase 2–3) — *Boltz-2 (below) is now the CORE engine, not an extension*

### IgFold
- **Role:** predict 3D antibody structure (VH+VL) for top-flagged variants; input to surface fingerprinting.
- **Install:** `pip install igfold`. API: `from igfold import IgFoldRunner`. CPU-adequate for a few; GPU faster.
- **Skill?** Yes. **License:** BSD-3. Trained on antibody structures (SAbDab) → better CDR loops than generic folders.

### dMaSIF
- **Full name:** differentiable Molecular Surface Interaction Fingerprinting.
- **Role:** rung-2 surface fingerprint of the CDR patch vs. self-protein surface.
- **Install:** `git clone github.com/FreyrS/dMaSIF`; pretrained weights in repo; `main_inference.py`. **GPU required.**
- **Skill?** Yes. **License:** **CC-BY-NC-ND 4.0** — non-commercial, no derivatives. Use as-is for demo; flag before commercial use.

### APBS + MSMS (rung-2 fallback)
- **Full names:** Adaptive Poisson-Boltzmann Solver; molecular surface mesher.
- **Role:** hand-rolled surface descriptor (electrostatics + hydrophobicity + H-bond) if no GPU / commercial-clean needed.
- **Install:** standalone binaries; `pip install pdb2pqr apbs`. **Skill?** Yes. **License:** open.

### Boltz-2  ← the CORE cofold engine (validated)
- **Full name:** open co-folding + binding-affinity foundation model (MIT Jameel Clinic + Recursion).
- **Role (validated):** the *primary* off-target signal — cofold the antibody Fv against **every** curated self-protein ectodomain (not a top-5 gate), 5 samples each, and score the interface with **`PAE_IF` (mean interface predicted-aligned-error, lower = tighter) + `epitope_reprod` (Jaccard of the contacted epitope across the 5 samples)**. Frozen hit rule `PAE_IF < 9.8 ∧ epitope_reprod ≥ 0.55`. **`ipTM` and the affinity head are NOT used — they over-dock** (lysozyme ipTM 0.87 ≈ everything; `../findings.md`). *(The old "rung 3 / top-5 confirmation, ranked on ipTM/affinity" framing is superseded.)*
- **Access (validated path):** the **hosted boltz.bio API** (Boltz-2.1) at ~$0.20/cofold — this is how the real-scale screen was run (`../demo/README.md`). Local `pip install boltz` (add `[cuda]` for GPU; `boltz predict input.yaml --use_msa_server`) and NVIDIA's Boltz-2 NIM (`pip install boltz2-python-client`) remain alternatives.
- **Skill?** Yes (CLI + YAML, or client API).
- **License:** **MIT — academic and commercial use.** This is why Boltz-2, not AlphaFold3, carries the confirmation rung.
- **Caveat:** antibody–antigen is the hardest cofolding category; use CDR-restricted interface confidence, treat as confirmation not certainty. `--use_msa_server` is one remote call.

### SAbDab / AbodyBuilder3 (stretch)
- **SAbDab:** Structural Antibody Database — check for a deposited camrelizumab structure (REST/sequence search, free academic).
- **AbodyBuilder3:** alt antibody structure predictor for cross-checking CDR-H3.

---

## Excluded — web-only or irrelevant

| Tool | Why |
|---|---|
| **BLOSUM (as the method)** | sequence-substitution scoring; blind to surface-driven cross-reactivity. Naive baseline only. |
| **IEDB linear epitopes (as the reference)** | linear peptides can't represent conformational off-targets; reference is folded surfaces. |
| **NetMHCpan / APE-Gen / PepSim / Expitope 2.0** | HLA/pMHC (TCR) tools; PepSim/Expitope also web-only. |
| **AlphaFold3 (for cofolding)** | gated weights; Boltz-2 is MIT + affinity + you hold credits. |

---

## Hour-zero verification

| Check | Tool |
|---|---|
| CDR-H3 == `QLYYFDYW` from SHR-1210 VH | ANARCI/ANARCII |
| Embeds SHR-1210 VH/VL, returns vectors | AntiBERTy / AbLang2 |
| Embeds VEGFR2 sequence | ESM-2 |
| VEGFR2 PDB (e.g. 4ASD) parses | Biopython/RCSB |
| LLM API key works | LLM |
| *(Phase 2)* IgFold folds SHR-1210; GPU present | IgFold |
| *(Phase 2)* `boltz predict` runs a toy protein–protein cofold + returns confidence/affinity | Boltz-2 |
