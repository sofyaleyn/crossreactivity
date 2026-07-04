# Demo build plan — autonomous brief

**Status:** active demo-build brief (2026-07-04). Companion to [`findings.md`](findings.md) (evidence) and [`plan.md`](plan.md) (the productization pipeline — *not* built here). This file is the self-contained brief for building the **hackathon demo** from already-validated data.

---

## GOAL (north star)

> Build a self-contained, offline-reproducible **hackathon demo** of CrossFlag, following the phased build, guardrails, and per-phase acceptance gates in [`demo-plan.md`](demo-plan.md). It proves one claim: *from antibody sequence alone, a Boltz-2 cofold screen against a curated self-protein set correctly flags known off-targets of a drug that harmed patients (SHR-1210 / camrelizumab), and does so specifically — not as generic stickiness.* The demo consists of (1) a results dashboard, (2) interface-structure figures, (3) a thin replayable "tool" that reproduces the verdict table from committed data, and (4) an honest pitch narrative. It must run with **no live API calls and no GPU**, using only data already committed in `data/results/` and `data/`. Every number shown must trace to `data/results/cofold_metrics.csv` or `findings.md`. Frame everything as **triage, not certification.**

**Done when ALL of these are objectively true (verifiable):**
1. On a clean checkout with the `crossflag-spike` conda env (+ `matplotlib`), a **single command** — `python -m crossflag.demo.build` (wired to `make demo`) — regenerates `demo/dashboard.html`, `demo/figures/*.png`, and `demo/verdict_table.{md,json}` from committed data alone, with **no network and no GPU**, and exits 0.
2. `pytest tests/test_demo.py` is **green**, and it asserts: (a) the Phase-0 recompute of `cofold-fzd5` from committed CIF/PAE matches `cofold_metrics.csv` within ±0.05 PAE_IF / ±0.03 epitope_reprod; (b) `demo/verdict_table.json` marks **FZD5 and ULBP2 = confirmed, VEGFR2 = missed, lysozyme = floor**; (c) `demo/dashboard.html` exists and embeds its data inline (opens with no external requests).
3. `demo/dashboard.html` renders **three charts** — calibrated panel (FZD5/ULBP2 near the PD-1 ceiling, VEGFR2 at the floor), within-fold FZD5 rank 2/12 (AUROC 0.909), and the SHR-1210-vs-pembrolizumab specificity control — each point traceable to a `cofold_metrics.csv` row.
4. Re-running the build twice produces **byte-identical** artifacts (determinism), and **no committed number contradicts** `cofold_metrics.csv` / `findings.md`.
5. `demo/README.md` narrates the ~5-minute pitch and every claim links to a figure or a `findings.md` experiment number.

A run that cannot make all five true must **stop and report which gate failed** rather than proceed or fabricate.

---

## What is already true (do not re-derive — cite it)

The viability study is **done and validated** (`findings.md`). The demo *presents* this evidence; it does not re-run the science.

- **24 cofold runs**, full metrics in `data/results/cofold_metrics.csv`; 25 structure dirs (`data/results/structures/cofold-*/sample_{0..4}.cif` + `sample_{0..4}_pae.npz`), all committed.
- **Headline result (the demo spine):**
  - PD-1 **ceiling** (cognate target): PAE_IF 7.24, epitope-reprod 0.936.
  - **FZD5** (known off-target): PAE_IF **5.69**, reprod **0.660** — confirmed, docks as tight as the cognate.
  - **ULBP2** (known off-target): PAE_IF **5.74**, reprod **0.893** — confirmed.
  - Lysozyme **floor** (non-binder): PAE_IF 12.38, reprod 0.446.
  - **VEGFR2** (weakest off-target): PAE_IF 11.55, reprod 0.430 — sits at the floor, **honestly missed** (show this, do not hide it).
- **Within-fold discrimination:** FZD5 ranks 2/12 among Frizzled CRDs by PAE_IF; within-family AUROC **0.909** (Exp 6).
- **Antibody-specificity control:** pembrolizumab (anti-PD-1, no such off-targets) collapses on the same off-targets (FZD5 reprod 0.66→0.22, ULBP2 0.89→0.34) while both antibodies still bind PD-1 (Exp 8).
- **Cost reframe:** ~$0.20/cofold → ~$50–150/antibody for a few-hundred-protein set, vs $10–30k for a wet-lab specificity screen.
- **Problem framing:** ~33% of antibody leads bind ≥1 unintended self-protein (Norden et al., mAbs 2024); SHR-1210 caused capillary hemangioma via a VEGFR2 CDR off-target.

## Metric semantics (from `data/results/README.md` — get these right)

- **PAE_IF_mean** — mean predicted aligned error over all antibody(H+L)↔antigen(V) token pairs, Å, averaged over 5 samples. **Lower = tighter interface.**
- **epitope_reprod** — mean pairwise Jaccard of contacted antigen residues (5 Å heavy-atom) across the 5 samples. **Higher = same epitope every time.** *Primary discriminating metric.*
- **ipTM** — recorded but **not usable** (over-docks; lysozyme scores 0.87). Never rank on ipTM.

---

## Hard guardrails (these keep an autonomous run honest)

1. **Do NOT build the `plan.md` pipeline** (candidate clustering, discovery screen, saturation/Chao1, panel screen). That is the productization path — represent it as *one diagram/slide*, not code.
2. **Do NOT lead with or depend on the germlined-mutant "drop WT / advance mutant" contrast.** Those mutant sequences are **reconstructed/unpublished** and the WT-vs-mutant VEGFR2 signal is inside the noise band (`findings.md` Exp 3/4). It may appear as a clearly-labelled "illustrative, reconstructed" aside only.
3. **The headline is the trio that is real and robust:** FZD5 + ULBP2 confirmation → within-fold AUROC 0.909 → pembrolizumab collapse.
4. **Show VEGFR2 as a miss.** Honesty about the weakest off-target is a feature; it is the "triage not certification" proof point.
5. **No invented numbers.** Every figure/stat reads from `cofold_metrics.csv` or is quoted from `findings.md` with the experiment number.
6. **Offline & deterministic by default.** No live `boltz.bio` call, no network, no GPU in the default path. Any live call is behind an explicit `--live` flag and must fail gracefully to the committed-data path.
7. **Framing = triage/prioritization** everywhere; cofold = confirmation *evidence*, not certainty (antibody–antigen is the hardest cofold class).
8. **$500 is a HARD budget ceiling for all live cofolding (Phase −1) — never exceed it.** Do not plan on a fixed per-cofold price (the repo quotes both ~$0.05 and ~$0.20). Enforce the cap *mechanically*: call `estimate-cost` before every batch, keep a persistent spend ledger, and **stop submitting at $480** (a $20 hard margin under $500). Re-price from a real pilot batch, not from any number in these docs.

---

## Environment

- Conda env **`crossflag-spike`** (Python 3.10) already exists. Run things as `conda run -n crossflag-spike python …`.
- Deps present: `torch` (CPU), `fair-esm`, `numpy`, `biopython`, `requests`, `Jinja2`. **`matplotlib` is NOT in the lock** — install it into the env (`conda run -n crossflag-spike pip install matplotlib`) and append to `data/results/scripts/requirements.lock.txt`. Optional 3D render: `biotite` (nice-to-have; skip if it fails to install cleanly).
- Reusable code to adapt (do not depend on their hardcoded scratchpad paths — repoint at committed `data/results/`):
  - `data/results/scripts/analyze_all.py` — epitope/PAE extraction logic (chain-by-length, 5 Å contacts, Jaccard). **Reuse the functions; replace `BASE` and the `boltz-runs/.../archive.tar.gz` reads with the committed `data/results/structures/<run>/sample_N.cif` + `sample_N_pae.npz`.**
  - `data/results/cofold_metrics.csv` — already-extracted metrics (dashboard reads this directly; no recomputation needed).

---

## Phase −1 — Live data generation: the real-scale screen (OPTIONAL, spends money)

**Do this only if you want the stronger demo.** It converts the current *n=1 retrospective anecdote* into a *prospective screen with a measured ROC and false-positive rate*, by running the actual product at scale. It is the one part of this plan that spends `boltz.bio` credits and requires live network + auth. Everything it produces lands as **committed data** that the offline Phases 0–5 then read, so the final demo stays offline/deterministic.

### Objective
With a **frozen** scoring rig, cofold SHR-1210 against a large, diverse slice of the curated 2,896-protein surfaceome and show that the known off-targets (FZD5, ULBP2) enrich to the top of *thousands* of proteins, with an honestly-reported precision / false-positive rate — plus, budget permitting, a **second independent anchor** run blind through the same rig (n=1 → n=2 generalization).

### The budget guard (build this FIRST, before any cofold)
This is the mechanism that makes "$500 hard ceiling" true regardless of the real price.
- `src/crossflag/demo/budget.py` + a persistent ledger `data/results/screen/spend_ledger.json` holding `{committed_estimate, confirmed_actual, per_job_price_samples[]}`.
- Constants: `PLANNING_TARGET = 450`, `HARD_STOP = 480`, `CEILING = 500`. The submitter **refuses to send any batch** whose `estimate-cost` would push `committed_estimate` past `HARD_STOP`; it shrinks the batch to fit or halts.
- **Estimate → submit → reconcile** loop per batch: (1) `boltz-api … estimate-cost` for the batch; (2) check against the ledger; (3) submit only what fits, with a stable `idempotency-key` per protein (so retries/failures never double-charge); (4) record the API's actual cost back into the ledger.
- On startup, reconcile the ledger against the already-existing 24 runs (those are free to re-pull via their prediction IDs — never re-pay for them).
- **Gate for this sub-step:** a dry-run that estimates a 2,000-job screen and prints projected spend *without submitting*, and a unit test proving the guard blocks a batch that would cross `HARD_STOP`.

### Step A — Price-discovery pilot (~$5 cap)
Submit a **20-protein pilot** (random draw from `self_proteins.csv`, ectodomains) to measure the **true per-cofold price** and the **success/failure rate** on this account. Everything downstream is sized from this measured price, not from any doc estimate. Record `per_job_price` in the ledger.

### Step B — Ectodomain prep (the real engineering cost)
Most of the 2,896 set are TM proteins; you cofold the **soluble ectodomain**, not full-length (`plan.md`). Build `src/crossflag/demo/ectodomain.py`: for each target, derive the extracellular region from the `surface_region` column / UniProt topology where available, else a documented heuristic (strip predicted TM spans, keep the largest extracellular segment). Cache to `data/results/screen/ectodomains/`. Trimming also *lowers* per-job cost if price scales with antigen length — a free squeeze.

### Step C — Freeze the hit-caller (pre-register BEFORE the screen)
Commit the exact decision rule to `src/crossflag/demo/hitcaller.py` and a frozen `data/results/screen/rig.json` **before** running the screen, so the screen cannot retro-fit it. Rule (from the existing anchor calibration — do not re-tune to the screen): per protein, over 5 samples compute `PAE_IF_mean` and `epitope_reprod`; **rank** by a composite calibrated against the PD-1 ceiling (7.24 / 0.936) and lysozyme floor (12.38 / 0.446); a protein is a **"hit"** if `epitope_reprod ≥ r*` AND `PAE_IF ≤ p*`, with `r*`, `p*` set from the floor + a fixed margin and written into `rig.json` now. The threshold does **not** move after seeing results.

### Step D — Allocation (priority-ordered, adaptive to the measured price)
Spend to the cap by priority; size each tier from the pilot price so total ≤ `PLANNING_TARGET` ($450), guard at `HARD_STOP` ($480). Reuse the 24 existing runs (free) — do not re-cofold FZD5/ULBP2/VEGFR2/PD-1/lysozyme/Frizzlens.

| Prio | Item | Target size | Why this size |
|---|---|---|---|
| 1 | **SHR-1210 × curated set** (flagship) | as many as budget allows, **min 1,200**, up to full 2,896 | ROC/enrichment stabilizes by ~1,200–1,500; a huge denominator is what makes EF and FPR meaningful |
| 2 | **Second anchor × diverse decoys** (blind) | ~400–600 (incl. its known off-target) | n=1 → n=2 generalization is worth more than proteins 1,500→2,896 of anchor #1 |
| 3 | **Scrambled-CDR control** on ~5 top hits | ~5–10 | proves it's the paratope, not the scaffold |

**Squeeze rules:** if the measured price is low, expand tier 1 toward the full 2,896 **and** fund tier 2; if high, protect the **two-anchor structure** (two anchors beat one full set) by capping tier 1 at ~1,200 and keeping tier 2. Stratify the tier-1 sample across protein size/fold so enrichment isn't a size artifact (`findings.md` Exp 5 confound). 5 samples per cofold is mandatory (epitope_reprod needs it) — never cut samples to save money.

### Step E — Pre-registered pass bars (write these down before results)
- **Flagship:** FZD5 **and** ULBP2 both land in the **top 5%** of the screened pool by the frozen rig. Report exact ranks, enrichment factor, precision@K, and the **false-positive rate at the frozen threshold** — whatever it is. VEGFR2 is *expected to miss*; report its rank honestly.
- **Second anchor:** its published off-target lands in the top 5% under the **same frozen rig, no retuning**.
- **Control:** scrambled-CDR loses the off-target signal (hit → non-hit).
- **Honesty:** if the FPR is high, that is a finding, not a failure — report it; do not move the threshold.

### Step F — Outputs (feed the offline build)
Commit to `data/results/screen/`: `screen_metrics.csv` (same schema as `cofold_metrics.csv`, one row per protein + prediction_id), the `rig.json`, `spend_ledger.json` (final reconciled spend, must show ≤ $500), CIFs/PAE for the top ~15 hits (for figures), and a `screen_manifest.json` (what was screened, sampling seed, ectodomain method). The offline Phases 1–3 read these instead of / in addition to the 24-run CSV.

### Second anchor — sourcing + fallback
`findings.md` names **ABT-736 → PF4**. Sourcing the VH/VL + the validated off-target sequence is the real bottleneck (curation, not credits). **If a validated pair cannot be sourced with public sequences, skip tier 2 and reallocate its budget to expanding the flagship toward the full 2,896-set** — do not fabricate an anchor.

### Phase −1 gate
`spend_ledger.json` shows reconciled total **≤ $500** (with the guard having enforced `HARD_STOP`); `screen_metrics.csv` covers ≥ 1,200 proteins with the frozen rig; the pre-registered pass bars are evaluated and reported (pass *or* fail, honestly).

---

## Build phases

> **If Phase −1 ran,** Phase 1 Chart A upgrades from the 8-point calibrated panel to the **full-screen enrichment view** (ROC / ranked list of 1,200+ proteins with FZD5/ULBP2 at the top, VEGFR2 shown missed, FPR annotated), and the verdict table (Phase 3) is computed over the screened set. If Phase −1 was skipped, Phases 0–5 run on the committed 24-run data exactly as written below.

Each phase has a **deliverable**, **files**, and a **self-verifiable acceptance gate**. Do phases in order; do not start a phase until the prior gate is green.

### Phase 0 — Setup & self-verification (~30 min)
**Deliverable:** confirmed env + a sanity check that committed data matches `findings.md`.
- Ensure `matplotlib` in `crossflag-spike`; create `demo/` and `src/crossflag/demo/` packages.
- Load `cofold_metrics.csv`; assert the eight anchor rows match the headline numbers above (PD-1 7.24/0.936, FZD5 5.69/0.660, ULBP2 5.74/0.893, lysozyme 12.38/0.446, VEGFR2 11.55/0.430) within rounding.
- For **one** run (e.g. `cofold-fzd5`), recompute PAE_IF + epitope_reprod from the committed `sample_*.cif` + `sample_*_pae.npz` using the adapted `analyze_all.py` logic, and assert it reproduces the CSV value (±0.05 PAE, ±0.03 reprod). *This proves the committed structures + the extraction path are self-consistent — the foundation for Phases 2–3.*
- **Gate:** both assertions pass. If the recompute disagrees with the CSV, STOP and report — do not paper over it.

### Phase 1 — Results dashboard (the money chart) (~2–3 h)
**Deliverable:** `demo/dashboard.html` — a single self-contained page (data embedded inline, no external requests) that is also publishable via the Artifact tool.
- **Chart A — Calibrated panel (headline):** scatter of the SHR-1210 runs on **epitope-reprod (x) × PAE_IF (y, inverted so up = tighter)**. Mark the PD-1 ceiling and lysozyme floor as reference bands. Highlight FZD5, ULBP2 (confirmed, green) and VEGFR2 (missed, amber, annotated "weakest off-target — honestly out of reach in silico").
- **Chart B — Within-fold discrimination:** FZD5 vs its 11 Frizzled-family decoys, ranked by PAE_IF; FZD5 highlighted at rank 2/12; caption "within-family AUROC 0.909 — discriminates the true off-target from same-fold negatives."
- **Chart C — Antibody-specificity control:** grouped bars, SHR-1210 vs pembrolizumab on {PD-1, FZD5, ULBP2, VEGFR2} for epitope-reprod; show both bind PD-1, pembro collapses on the off-targets.
- A header stating the problem + the cost reframe; a footer with the honesty caveat.
- Theme-aware (light/dark), responsive, wide charts scroll in their own container. Emit as an HTML file suitable for the **Artifact** tool (favicon 🧬, title "CrossFlag — off-target cofold screen").
- **Gate:** page opens standalone in a browser with all three charts rendered from embedded data; every plotted point traces to a `cofold_metrics.csv` row (include a hidden `data-source` attribute or a small "data provenance" line).

### Phase 2 — Interface-structure figures (the "wow") (~2–3 h)
**Deliverable:** PNGs in `demo/figures/` embedded into the dashboard.
- **Required (robust, 2D):** for the key contrast **SHR-1210×FZD5 vs pembrolizumab×FZD5** (and, if time, ×ULBP2):
  - **Interface PAE heatmap** — antibody(H+L) rows × antigen(V) cols block of the PAE matrix, best sample, shared color scale. SHR-1210 shows a tight low-PAE interface block; pembro is diffuse. Directly visualizes PAE_IF.
  - **Epitope-reproducibility map** — 1D strip over antigen residues, colored by *how many of the 5 samples* contact each residue. SHR-1210 shows a sharp reproduced patch; pembro is scattered. Directly visualizes epitope_reprod.
- **Optional (nice-to-have, 3D):** cartoon render of the SHR-1210×FZD5 complex (sample 0) with CDR-H3 and the contacted antigen patch highlighted, via `biotite` or an offline PyMOL/ChimeraX script if available. If it doesn't render cleanly in autonomous mode, skip — the 2D figures carry the point.
- **Gate:** the SHR-1210 vs pembro contrast is visually obvious in the PAE heatmap and the epitope map; figures are generated by a script (`src/crossflag/demo/figures.py`), not hand-placed.

### Phase 3 — Thin replayable "tool" (~2–3 h)
**Deliverable:** a small `src/crossflag/` slice that makes the demo a *tool*, not just slides.
- `src/crossflag/demo/scoring.py` — read the committed structures for a given run and return `{PAE_IF, epitope_reprod, verdict}` where verdict is computed against the frozen calibrated panel (ceiling PD-1 / floor lysozyme). Adapt `analyze_all.py`; keep the calibration thresholds explicit and documented.
- `src/crossflag/demo/run.py` — a CLI: input = SHR-1210 VH/VL (from `data/anchor/`) + a list of reference runs; output = the ranked **verdict table** (protein, PAE_IF, epitope_reprod, read: confirmed / floor / missed, + the named confirmation assay from a small map) printed to stdout **and** written to `demo/verdict_table.{md,json}`. Default path replays committed data → fully offline & deterministic.
- **Optional `--live`:** a `cofold_client.py` that submits via the `boltz-api` CLI using the prediction IDs in `cofold_metrics.csv` (re-runs are free). Must detect missing auth/network and fall back to the committed path with a clear message — never hard-fail the demo.
- **Gate:** `python -m crossflag.demo.run` prints the verdict table offline; FZD5/ULBP2 read "confirmed", VEGFR2 reads "missed (weakest)", lysozyme "floor". Output matches the dashboard numbers.

### Phase 4 — Pitch narrative (~1–2 h)
**Deliverable:** `demo/README.md` — the ~5-minute script, plus a one-page problem/solution/cost/honesty writeup.
- Structure: Problem (30s: SHR-1210 hemangioma via VEGFR2; 33% of leads) → Result (90s: dashboard Chart A — flags FZD5/ULBP2 from sequence alone) → Specificity (60s: Charts B+C — within-fold + pembro control, "not generic stickiness") → Structures (60s: Phase-2 figures) → Cost + Honest frontier (30s: $50–150/antibody vs $10–30k; VEGFR2 missed; triage not certification).
- One diagram of the `plan.md` representative-set scaling path labelled "how this scales to a 1000-candidate panel (productization, not in this demo)."
- **Gate:** README reads end-to-end as a coherent pitch; every claim links to a figure or a `findings.md` experiment number.

### Phase 5 — Assembly & self-check (~1 h)
**Deliverable:** one entry point + a passing self-check.
- `make demo` (or `python -m crossflag.demo.build`) regenerates dashboard + figures + verdict table into `demo/` from scratch, deterministically.
- Add `tests/test_demo.py`: dashboard exists and embeds data; verdict table has FZD5/ULBP2=confirmed, VEGFR2=missed; Phase-0 recompute assertion holds. Update `README.md`'s status banner to point at `demo/` as the runnable demo.
- **Gate:** clean checkout + env → one command → all artifacts regenerate → `pytest tests/test_demo.py` green.

---

## Suggested order & rough budget
Phase 0 → 1 → 2 → 3 → 4 → 5. If time-boxed: **Phases 0–1 alone are a passable demo** (the dashboard is the pitch). Phase 2 is the biggest credibility multiplier per hour. Phase 3 makes it a "tool." Total ~10–13 h of autonomous work.

## Success criterion (single, binary)
On a clean checkout, one command regenerates a dashboard + interface figures + a verdict table that show, from committed data alone: **SHR-1210 flags FZD5 and ULBP2 as confident reproducible off-targets near the PD-1 ceiling, discriminates FZD5 within its own fold (AUROC 0.909), and does so specifically (pembrolizumab collapses)** — with VEGFR2 honestly shown as missed. Triage, not certification.
