# CrossFlag

**Panel-level CDR cross-reactivity triage for therapeutic antibodies.**

Given a panel of antibody candidates (or CDR variants of one lead), CrossFlag ranks which variants are most likely to cross-react with a self-protein, corroborates the top suspects with 3D surface comparison, confirms the highest-risk pairs by cofolding the antibody against the flagged self-protein, and names the assay to run first.

> *"Feed us your antibody panel. We embed each variant's binding region with an antibody language model, rank the panel against a self-protein reference set, confirm the top suspects by cofolding, and tell you which variant to advance and which assay confirms it — before you spend a wet-lab specificity screen."*

CDR = complementarity-determining region (the antibody loops that contact antigen). PLM = protein language model.

---

## The problem

Antibody CDR-mediated off-target binding (polyspecificity) is the largest single cause of biologic drug failure in preclinical and early clinical development. About 33% of antibody lead candidates bind at least one unintended self-protein (Norden et al., mAbs 2024). For cytotoxic modalities — ADCs (antibody-drug conjugates), bispecifics, CAR-T-targeting antibodies — a CDR that cross-reacts with a self-protein can be fatal.

The critical fact for our design: **CDR cross-reactivity is driven by surface shape and chemistry, not sequence identity.** Documented cases cross-react with as low as 7% sequence identity between on- and off-target. Any method that scores raw sequence similarity is close to blind to the real signal. CrossFlag is built around learned representations and 3D structure instead.

## How it differs from existing tools

Whole-proteome scanners (ARDitox, Dec 2025; Tope-seq, Jan 2026) take one candidate and scan the full proteome. CrossFlag differs on two axes:

1. **Curated self-protein reference set** (surfaces, not the whole proteome) — precision over recall.
2. **Panel-level agentic triage with a confirmation ladder** — rank the panel, then spend expensive compute (cofolding) only on the top suspects.

## The evidence ladder

CrossFlag stacks three independent signals, cheap to expensive, so nothing costly runs on the whole set:

```
  1. PLM-embedding similarity   →  rank the whole panel        (cheap, all variants)
  2. 3D surface fingerprint     →  corroborate the top flags   (medium, top ~10)
  3. Boltz-2 cofolding          →  confirm the top suspects    (expensive, top ~5)
```

Each rung answers a sharper question: *looks similar → presents a similar surface → actually docks.*

### What rung 1 is really measuring (and what it is not)

Rung 1 is a triage signal, not a binding predictor — and it's worth being precise about why. The paratope vector (an antibody model, AntiBERTy/AbLang2) and the antigen vector (a general protein model, ESM-2) come from two separately-trained models, and binding is *complementarity* — opposite shapes and charges fitting together — not resemblance. So a cosine between those two vectors is not a physical binding score and we never treat it as one.

What it *can* pick up is narrower and honest. Both models are trained on protein sequences, so both encode the same underlying biophysical vocabulary — hydrophobicity, charge, aromatic content, local structural propensity. Two very different objects can still be characterized by the same underlying chemistry. Calibrated against a benign background, the cross-model cosine becomes a fuzzy proxy for one claim: *"this paratope's chemistry is unusually enriched for the kind of chemistry present on this antigen's surface, relative to proteins we know are safe."*

That is a **smoke detector, not a binding model.** It's cheap, it runs on the whole panel in a single pass, and it's allowed to be wrong — because every suspect it surfaces is re-checked by the surface and cofolding rungs. Rung 1 buys a ranked shortlist; the ladder above it earns the verdict.

The version that would be a true binding signal — an upgrade beyond the MVP — is a nearest-neighbor search over *known antibody–self-protein binding pairs* rather than bare self-protein sequences: **similar paratopes bind similar things.** That comparison stays within one model (no cross-space mismatch) and respects complementarity by transfer — a known binder has already demonstrated the fit, so we're only flagging lookalikes. It requires a curated set of known binders we don't yet have (Phase 3).

## The anchor (validation) case

**SHR-1210 (camrelizumab, anti-PD-1) → VEGFR2 off-target → capillary hemangioma in patients.**

SHR-1210 is a humanized IgG4κ anti-PD-1 antibody that caused an unusual capillary hemangioma in trials — not seen with other anti-PD-1s. Receptor proteome screening traced it to CDR-mediated binding to **VEGFR2** (plus FZD5, ULBP2). VEGFR2 agonism drives the hemangioma; a VEGFR2 antagonist rescued patients, confirming the mechanism. CDR germlining ablated VEGFR2 binding while preserving PD-1 affinity — a ready-made variant panel.

This gives us: a named drug, named self-protein, CDR-confirmed mechanism, human adverse event, and a real WT-plus-mutants panel.

### Anchor sequences (public — Thera-SAbDab / US11208484B2)

```
Heavy Chain Fv (VH):
EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKG
RFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS

Light Chain Fv (VL):
DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGS
GTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK
```

Expected CDR-H3 `QLYYFDYW`, CDR-L3 `QQVYSIPWT` (ANARCI confirms). No HLA, no pMHC, no presentation gate — antibodies bind surfaces directly.

---

## Repo structure

```
crossflag/
├── README.md                       ← this file (overview + repo map)
├── HANDOFF.md                      ← build order + module contracts for Claude Code (START HERE to build)
├── pyproject.toml                  ← packaging (src layout), deps, pytest/ruff config
├── .gitignore
├── docs/
│   ├── roadmap.md                  ← ranked steps + tools; MVP → extensions  (READ FIRST)
│   ├── mvp-spec.md                 ← Phase 1: PLM-embedding panel ranking
│   ├── extensions-spec.md          ← Phase 2: surface fingerprint + Boltz-2 cofolding
│   ├── reference-set.md            ← how the curated self-protein set is BUILT (layers A/B/C)
│   ├── demo-run.md                 ← MVP acceptance script (Beats 1–3) — canonical pass thresholds
│   ├── glossary.md                 ← abbreviations expanded
│   └── tools.md                    ← every tool: role, install, API, license, skill status
├── data/                           ← DATA ONLY (no executable code lives here)
│   ├── anchor/
│   │   ├── shr1210_vh.fasta
│   │   ├── shr1210_vl.fasta
│   │   ├── variants/               ← CDR-germlining mutants (Finlay et al.)
│   │   └── offtargets/             ← VEGFR2, FZD5, ULBP2 sequences + PDB IDs
│   └── reference/                  ← curated self-protein reference set (inputs + outputs)
│       ├── raw/                    ← downloaded sources (gitignored, regenerable)
│       │   ├── surfaceome/         ← SURFY + CSPA xlsx, seq_cache.tsv
│       │   └── fasta/              ← per-accession UniProt sequence fetches
│       ├── seeds/                  ← hand-curated inputs, committed (provenance)
│       │   └── layer_c_mimicry_seed.csv
│       ├── build/                  ← intermediates (gitignored, deterministic)
│       ├── self_proteins.csv       ← curated reference set — final output (committed)
│       ├── background/
│       │   └── benign_proteins.csv ← non-cross-reactive self-proteins (calibration, committed)
│       └── index/                  ← ESM-2 embedding cache (.npz / FAISS; gitignored)
├── src/crossflag/
│   ├── extract/cdrs.py             ← ANARCI/ANARCII → CDR annotation
│   ├── embed/
│   │   ├── antibody.py             ← AntiBERTy / AbLang2 embeddings (variant side)
│   │   └── antigen.py              ← ESM-2 embeddings (self-protein side)
│   ├── reference/                  ← reference-set build pipeline (see docs/reference-set.md)
│   │   ├── paths.py                ← single source of truth for data/ locations (repo-root-anchored)
│   │   ├── build_set.py            ← orchestrator: A → anchor → background → B → C → merge → embed
│   │   ├── layer_a.py              ← SURFY + CSPA cell-surface base
│   │   ├── layer_b.py              ← AAgAtlas autoantigens (TODO: not yet implemented)
│   │   ├── layer_c.py              ← hand-curated pathogen-mimicry seed fill
│   │   ├── background.py           ← benign background sampling (calibration)
│   │   └── merge.py                ← merge layers → self_proteins.csv
│   ├── rank/
│   │   ├── embedding_rank.py       ← rung 1: embedding-space similarity + calibration
│   │   └── panel.py                ← panel assembly + ranking
│   ├── structure/
│   │   ├── fold.py                 ← IgFold antibody structure (rung 2 input)
│   │   ├── surface.py              ← dMaSIF / APBS surface fingerprint (rung 2)
│   │   └── cofold.py               ← Boltz-2 cofolding + affinity (rung 3)
│   ├── assay/map.py                ← self-protein → confirmation assay map
│   ├── agent/report.py             ← LLM narration + recommendation
│   └── pipeline.py                 ← end-to-end orchestration
├── skills/                         ← agent-skill wrappers (all programmatic)
│   ├── anarci_extract.md
│   ├── antibody_embed.md
│   ├── antigen_embed.md
│   ├── embedding_rank.md
│   ├── igfold_surface.md           (extension)
│   └── boltz_cofold.md             (extension)
├── notebooks/anchor_validation.ipynb
└── tests/
    ├── test_cdr_extract.py         ← CDR-H3 == QLYYFDYW
    ├── test_panel_rank.py          ← VEGFR2 ranks near top for WT
    └── test_cofold_gate.py         ← WT+VEGFR2 interface confidence > germlined mutant
```

## Where to start

- **To build:** read [HANDOFF.md](HANDOFF.md) — it gives Claude Code the build order and module contracts.
- **To understand:** [docs/roadmap.md](docs/roadmap.md) → [docs/mvp-spec.md](docs/mvp-spec.md) → [docs/extensions-spec.md](docs/extensions-spec.md).
- **Reference set:** [docs/reference-set.md](docs/reference-set.md) (how it's built). **Demo / acceptance:** [docs/demo-run.md](docs/demo-run.md).
- **Terminology:** [docs/glossary.md](docs/glossary.md). **Tools/install:** [docs/tools.md](docs/tools.md).

## Design constraint: everything is a programmatic skill

Every tool runs locally or via a callable API and is wrapped as an agent skill. Web-form tools (PepSim, Expitope 2.0) are excluded by design.

## Success criterion (single, binary)

On the SHR-1210 panel: CrossFlag ranks **VEGFR2 near the top of the flag list for wild-type SHR-1210**, ranks a CDR-germlining mutant as cleaner, and — for the top suspect — the Boltz-2 cofold shows a confident WT–VEGFR2 interface that collapses for the germlined mutant. From antibody sequence alone.
