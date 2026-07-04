# Tools reference

Every tool is programmatic (local or callable API) and skill-wrappable. Web-only tools excluded. **Skill?** = drivable without runtime human interaction.

---

## MVP tools (Phase 1)

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

## Extension tools (Phase 2–3)

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

### Boltz-2  ← the confirmation gate
- **Full name:** open co-folding + binding-affinity foundation model (MIT Jameel Clinic + Recursion).
- **Role:** rung 3 — cofold top ~5 variant×self-protein pairs; read interface confidence (ipTM/PAE at CDR contacts) + affinity score.
- **Install:** `pip install boltz` (add `[cuda]` for GPU). Run: `boltz predict input.yaml --use_msa_server`. Affinity via YAML affinity request → `affinity_pred_value`, `affinity_probability_binary`.
- **Also:** NVIDIA hosts a Boltz-2 NIM with a Python client (`pip install boltz2-python-client`) if you prefer a remote endpoint over local GPU.
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
