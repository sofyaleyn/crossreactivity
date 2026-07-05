# CrossFlag

**Cofold-based CDR off-target triage for therapeutic antibodies.**

> **STATUS (2026-07-05): validated at screen scale.** CrossFlag is a **Boltz-2 cofold screen**. An empirical viability study ([`findings.md`](findings.md)) found the project's original embedding-similarity design broken, and validated this replacement instead: cofold an antibody Fv against a curated human self-protein surfaceome and score each interface against a calibrated ceiling/floor panel. A real-scale run ([`data/results/screen/validation_summary.md`](data/results/screen/validation_summary.md)) screened SHR-1210 against **1,198 self-protein ectodomains** and recovered its known off-targets. Verdict: **conditionally validated — a prioritization tool for wet-lab specificity screens, not a replacement.**
>
> The original 3-rung "evidence ladder" (embedding rung 1 + surface rung 2 + cofold rung 3) is **superseded and archived** at the bottom of this file; do not read it as a current claim.
>
> **Runnable demo:** **[`demo/`](demo/README.md)** — one offline command (`make demo`, no network/GPU) regenerates a dashboard + interface figures + a verdict table from committed data. Start here to *see* the result.
>
> **Read in this order:**
> 1. **[`demo/README.md`](demo/README.md)** — the ~5-minute pitch + one-command offline demo of the validated screen.
> 2. **[`plan.md`](plan.md)** — the active build plan (representative-set off-target screen). **Build from this.**
> 3. **[`findings.md`](findings.md)** — the empirical viability log: what was tested, the evidence, and why the cofold screen is the surviving design.
> 4. **[`docs/glossary.md`](docs/glossary.md)** — terminology.
>
> `mvp-spec.md`, `demo-run.md`, `HANDOFF.md`, `roadmap.md`, `tools.md`, and the "Original design" section at the bottom of this README describe the **superseded** embedding-ladder design, kept for history.

Given a therapeutic antibody (or a panel of CDR variants of one lead), CrossFlag cofolds its Fv (VH+VL) against a curated human surfaceome with Boltz-2, scores every interface, and ranks the self-proteins its CDRs are most likely to bind by mistake — so the wet-lab specificity screen tests the right suspects first.

> *"Feed us an antibody's VH/VL. We cofold it against a curated human self-protein reference set, score each interface against a calibrated true-target ceiling and non-binder floor, and hand you a ranked, paratope-specific shortlist of likely off-targets — before you spend a wet-lab specificity screen."*

CDR = complementarity-determining region (the antibody loops that contact antigen). Fv = the variable-domain fragment (VH+VL) that carries the CDRs.

---

## The problem

Antibody CDR-mediated off-target binding (polyspecificity) is the largest single cause of biologic drug failure in preclinical and early clinical development. About 33% of antibody lead candidates bind at least one unintended self-protein (Norden et al., mAbs 2024). For cytotoxic modalities — ADCs (antibody-drug conjugates), bispecifics, CAR-T-targeting antibodies — a CDR that cross-reacts with a self-protein can be fatal.

The critical fact for the method: **CDR cross-reactivity is driven by surface shape and chemistry, not sequence identity.** Documented cases cross-react with as low as 7% sequence identity between on- and off-target. Any method that scores raw sequence similarity is close to blind to the real signal — and, as the viability study found, so is cross-model embedding similarity ([`findings.md`](findings.md) Exp 1–2, 5). CrossFlag reads the actual 3D interface instead, by cofolding the antibody against each candidate self-protein and scoring how tightly and reproducibly they dock.

## How it works now

CrossFlag is a single validated stage: **cofold, then score against a calibrated panel.**

1. **Curate** a human self-protein reference surfaceome (`data/reference/self_proteins.csv`) — the **2,886-protein SURFY surfaceome** (Bausch-Fluck et al., PNAS 2018; 948 mass-spec-confirmed via CSPA) plus curated autoimmune/mimicry additions = **2,896 proteins** (99.8% with UniProt sequences; see [`docs/reference-set.md`](docs/reference-set.md)) — cofolding each against the antibody as its extracellular **ectodomain** where it spans the membrane.
2. **Cofold** the antibody Fv against every reference protein with Boltz-2 (5 samples per pair, hosted API, ~$0.20/cofold).
3. **Score** each pair on two metrics: **`PAE_IF`** (interface predicted-aligned-error — lower = tighter) and **`epitope_reprod`** (do the 5 samples agree on the same antigen epitope — higher = reproducible). **`ipTM` is dropped — it over-docks** (lysozyme non-binder scores ~0.87, like everything).
4. **Rank** against a per-antibody calibrated panel: a cognate true-target **ceiling** (e.g. PD-1) and a non-binder **floor** (lysozyme). Frozen hit rule: **`PAE_IF < 9.8 Å AND epitope_reprod ≥ 0.55`**. A valid-regime hit on a small (<150 aa) ectodomain is flagged **low-confidence** (see limits below).

The original premise that cofolding was "too expensive to run against every self-protein" (hence the need for a cheap embedding pre-filter) turned out to be **wrong at hosted prices**: a few-hundred-protein set costs ~$50–150/antibody versus ~$10–30k for a wet-lab specificity screen. So the cheap pre-filter was dropped and CrossFlag cofolds the curated set directly. For a *panel* of many candidate variants, [`plan.md`](plan.md) scales this via a representative-set approach (cluster candidates, cofold representatives to discover the shortlist, then screen all candidates against it) rather than the full candidate × target grid.

## The validated result (in brief)

A real-scale screen ([`data/results/screen/validation_summary.md`](data/results/screen/validation_summary.md)) ran the **frozen** pipeline — thresholds fixed before the screen, no tuning to these results:

- **Scale:** SHR-1210 Fv × **1,198 self-protein ectodomains** (2 dropped on an API rate-limit), plus a blind second antibody and scramble controls — **1,603 cofolds** total, 5 samples each, **$342.70** of a $500 budget.
- **Specificity:** in the valid regime (ectodomain ≥150 aa) the false-positive rate is **2.1% (14/678)**, falling monotonically with antigen size.
- **Sensitivity:** **both recoverable known off-targets recovered** — FZD5 (PAE_IF 5.69 / reprod 0.660) and ULBP2 (5.74 / 0.893) enriched out of ~1,200 proteins. **VEGFR2 is missed** — it is SHR-1210's weakest, lowest-affinity off-target, and it is also the clinically-causal one (the honest wound).
- **Generalization:** a blind second antibody, **ABT-736 (anti-Aβ) → PF4**, recovered with **zero retuning** (PF4 flagged, #1 by PAE among 76 same-size decoys).
- **Paratope-specific:** scrambling all six CDRs (framework preserved) collapses both off-targets to non-hits; within its own structural family FZD5 discriminates at **AUROC 0.909** ([`findings.md`](findings.md) Exp 6); a different anti-PD-1 (pembrolizumab) binds the shared target but collapses on the off-targets.
- **Known bounded flaw:** the method **over-docks small (<150 aa) ectodomains** (~20% FPR, ~43% of the surfaceome), handled by a **size-aware confidence gate** rather than trusting or discarding those hits.

**Bottom line: conditionally validated.** Triage / prioritization, not certification. Every quantitative claim above traces to [`findings.md`](findings.md) or [`data/results/screen/validation_summary.md`](data/results/screen/validation_summary.md).

## The anchor (validation) case

**SHR-1210 (camrelizumab, anti-PD-1) → VEGFR2 off-target → capillary hemangioma in patients.**

SHR-1210 is a humanized IgG4κ anti-PD-1 antibody that caused an unusual reactive capillary hemangioma in trials — not seen with other anti-PD-1s. Receptor proteome screening traced it to CDR-mediated binding to **VEGFR2** (plus FZD5, ULBP2). VEGFR2 agonism drives the hemangioma; a VEGFR2 antagonist rescued patients, confirming the mechanism. This gives us a named drug, named self-proteins, a CDR-confirmed mechanism, and a human adverse event — a real, published anchor to validate against.

### Anchor sequences (public — Thera-SAbDab / US11208484B2)

```
Heavy Chain Fv (VH):
EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKG
RFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS

Light Chain Fv (VL):
DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGS
GTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK
```

CDR-H3 `QLYYFDYW`, CDR-L3 `QQVYSIPWT` (ANARCI confirms). No HLA, no pMHC, no presentation gate — antibodies bind surfaces directly. The second, blind anchor (ABT-736 → PF4) is described in [`data/results/screen/validation_summary.md`](data/results/screen/validation_summary.md).

---

## Repo structure

```
crossflag/
├── README.md                       ← this file (overview + repo map)
├── plan.md                         ← ACTIVE build plan: representative-set off-target screen  (BUILD FROM THIS)
├── findings.md                     ← empirical viability log + current verdict (the plan rests on this)
├── HANDOFF.md                      ← original build order/contracts (superseded — embedding-ladder design)
├── Makefile                        ← `make setup` / `make demo` / `make test-demo`
├── pyproject.toml                  ← packaging (src layout), deps, pytest/ruff config
├── demo/                           ← runnable offline demo of the validated screen
│   ├── README.md                   ← ~5-minute pitch + one-command demo  (START HERE)
│   ├── dashboard.html              ← regenerated dashboard (embeds the data rows)
│   ├── figures/                    ← screen + calibration figures (regenerated)
│   └── verdict_table.{md,json}     ← regenerated verdict table
├── docs/                           ← ORIGINAL design docs (superseded except glossary/tools)
│   ├── glossary.md                 ← abbreviations expanded (still current)
│   ├── tools.md                    ← tool roles/install/licenses
│   ├── roadmap.md · mvp-spec.md · extensions-spec.md · demo-run.md   ← superseded embedding-ladder specs
│   └── reference-set.md            ← how the curated self-protein set is built (layers A/B/C)
├── data/                           ← DATA ONLY (no executable code lives here)
│   ├── anchor/                     ← SHR-1210 WT VH/VL, CDR-germline variants, off-target FASTAs
│   ├── anchor2/                    ← second anchor (ABT-736 → PF4)
│   ├── reference/                  ← curated self-protein reference set
│   │   ├── self_proteins.csv       ← curated reference set (~2,896 proteins, committed)
│   │   ├── ectodomains.csv         ← extracted extracellular regions used for cofolding
│   │   ├── ectodomains_report.md · background/ · seeds/ · raw/ (regenerable)
│   └── results/                    ← all experimental artifacts
│       ├── cofold_metrics.csv      ← 24 calibration cofolds (PAE_IF, epitope_reprod, ipTM, IDs)
│       ├── benchmark_*.json · tier1_filter_scores.csv   ← fold-matched + cheap-filter benchmarks
│       ├── structures/ · inputs/ · scripts/             ← CIFs+PAE, Boltz inputs, analysis code
│       ├── README.md               ← human-readable metric tables + semantics
│       └── screen/                 ← the real-scale screen
│           ├── validation_summary.md   ← screen results write-up
│           ├── screen_metrics.csv      ← per-protein scores
│           ├── rig.json                ← frozen thresholds
│           ├── spend_ledger.json       ← budget
│           ├── checks.py               ← the validation analysis
│           └── jobs/ · *_inputs/ · *_runs/   ← job configs + raw runs
├── src/crossflag/
│   ├── screen/ectodomain.py        ← extract extracellular region for cofolding
│   ├── structure/cofold.py         ← Boltz-2 cofold + interface scoring (the core stage)
│   ├── structure/fold.py · surface.py   ← (superseded rung-2 surface path)
│   ├── demo/                        ← offline demo builder
│   │   ├── build.py · run.py        ← `python -m crossflag.demo.build` orchestration
│   │   ├── panel.py · scoring.py    ← frozen hit rule + metric extraction
│   │   ├── screen_view.py · figures.py · dashboard.py · diagram.py   ← views/figures
│   │   └── paths.py                 ← repo-root-anchored data locations
│   ├── reference/                  ← reference-set build pipeline (layers A/B/C; see docs/reference-set.md)
│   ├── extract/cdrs.py             ← ANARCI/ANARCII → CDR annotation
│   ├── embed/ · rank/              ← (superseded rung-1 embedding path — kept, not used)
│   ├── assay/map.py                ← self-protein → confirmation assay map
│   ├── agent/report.py             ← LLM narration + recommendation
│   └── pipeline.py                 ← end-to-end orchestration
├── skills/                         ← agent-skill wrappers (all programmatic)
├── notebooks/anchor_validation.ipynb
└── tests/
    ├── test_demo.py                ← demo acceptance gate (recomputes PAE_IF/epitope_reprod from CIF+PAE)
    ├── test_cdr_extract.py         ← CDR-H3 == QLYYFDYW
    ├── test_panel_rank.py · test_cofold_gate.py
```

## Where to start

- **To see the result:** [`demo/README.md`](demo/README.md) → `make demo` → open `demo/dashboard.html`.
- **To build:** [`plan.md`](plan.md) (active) — the representative-set screen, module list, cost model, K-sizing, anchor-validation steps.
- **To understand the evidence:** [`findings.md`](findings.md) (9 experiments + current verdict).
- **Terminology:** [`docs/glossary.md`](docs/glossary.md). **Tools/install:** [`docs/tools.md`](docs/tools.md).

## Design constraint: everything is a programmatic skill

Every tool runs locally or via a callable API and is wrapped as an agent skill. Web-form tools (PepSim, Expitope 2.0) are excluded by design.

---

## Original design (superseded — kept for history; see [`findings.md`](findings.md) for why)

> **Everything below describes the project's ORIGINAL design and is no longer a live claim.** The viability study ([`findings.md`](findings.md)) found the embedding-similarity rung-1 method broken (fails a fold-matched benchmark at ~chance, Exp 1–2 & 5) and found that **no cheap pre-filter — embedding, surface, or biophysical — beats a free annotation baseline** (Exp 7, 9). The "cheap-to-expensive evidence ladder" narrative is dead; the cofold confirmer is the surviving, validated tool (described above). This section is retained only so the design history and its refutation are legible.

### How it *was* pitched to differ from existing tools

Whole-proteome scanners (ARDitox, Dec 2025; Tope-seq, Jan 2026) take one candidate and scan the full proteome. CrossFlag was to differ on two axes:

1. **Curated self-protein reference set** (surfaces, not the whole proteome) — precision over recall. *(This axis survives.)*
2. **Panel-level agentic triage with a confirmation ladder** — rank the panel cheaply, then spend expensive compute (cofolding) only on the top suspects. *(This axis is dead — no cheap ranking works, and cofolding the whole curated set is affordable.)*

### The evidence ladder (abandoned)

The original design stacked three independent signals, cheap to expensive:

```
  1. PLM-embedding similarity   →  rank the whole panel        (cheap, all variants)   [BROKEN — Exp 1-2,5]
  2. 3D surface fingerprint     →  corroborate the top flags   (medium, top ~10)       [NO LIFT — Exp 7]
  3. Boltz-2 cofolding          →  confirm the top suspects    (expensive, top ~5)     [VALIDATED — became the whole tool]
```

Only rung 3 survived — and, at hosted cofold prices, it absorbed the whole job (cofold the curated set directly) rather than sitting atop a cheap pre-filter.

### What rung 1 was supposed to measure (and why it failed)

Rung 1 was pitched as a triage smoke-detector, not a binding predictor: a calibrated cross-model cosine between an antibody-model paratope vector (AntiBERTy/AbLang2) and an ESM-2 antigen vector, read as *"this paratope's chemistry is unusually enriched for the chemistry on this antigen's surface, relative to safe proteins."* The intended upgrade was a nearest-neighbor search over *known* antibody–self-protein binding pairs (similar paratopes bind similar things), staying within one model to respect complementarity.

**Why it was dropped:** the cross-model cosine is undefined as written (mismatched dimensions) and biophysically backwards (binding is complementarity, not similarity). Even the charitable same-model formulation ranks the true off-targets at ~chance within their own structural family (mean within-family AUROC 0.58; ULBP2 below chance — [`findings.md`](findings.md) Exp 5), and on the real 2,896-protein set no learned filter beats the free "is-it-a-receptor" annotation baseline (Exp 9). There is no cheap route to a small per-candidate shortlist.

### Original success criterion (retired)

*On the SHR-1210 panel: rank VEGFR2 near the top for wild-type, a germlined mutant as cleaner, and confirm the top suspect by cofold.* Retired because the germlined-mutant panel is not published (only a reconstruction recipe exists, [`findings.md`](findings.md) Exp 3) and because embedding-based ranking of VEGFR2 was shown to be a fold/size confound, not signal. The current success criterion is the screen-scale validation in [`data/results/screen/validation_summary.md`](data/results/screen/validation_summary.md).
