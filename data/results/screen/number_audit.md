# CrossFlag headline-number audit

**Auditor:** independent adversarial recompute · **Date:** 2026-07-05
**Method:** every number below was recomputed from the raw CSVs with a fresh stdlib
script (scratchpad, not committed), then cross-checked against `checks.py`. The `hit`
column in `screen_metrics.csv` was verified to match the frozen rule
(`PAE_IF < 9.8 AND epitope_reprod ≥ 0.55`) on **all 1,603 rows with 0 mismatches**.

## AUDIT VERDICT: numbers reproduce — with one prose discrepancy and bookkeeping/framing caveats

The scientific/quantitative headline claims (FPR bands, valid-regime FPR, off-target
enrichment, generalization, negative control, top candidates) **all reproduce**. One
stated count is wrong (168 → actually **171**), and the summary text contains a garbled
"Bottom line" line and a `1,603 total` head-count that does not add up as written. None
of these change the conclusions, but the "168" and the "~3% FPR" phrasing should be fixed
before the pitch.

## Verdict table

| # | Claim | Recomputed | ✓/✗/⚠ | Note |
|---|---|---|---|---|
| 1 | Valid-regime (≥150aa) FPR = 2.1% (14/678) flagship | 14/678 = **2.06%** | ✓ | Rounds to 2.1%. |
| 2a | 0–80aa FPR 20.6% | 80/389 = **20.6%** | ✓ | |
| 2b | 80–150aa FPR 19.1% | 25/131 = **19.1%** | ✓ | |
| 2c | 150–300aa FPR 4.5% | 11/243 = **4.5%** | ✓ | |
| 2d | 300–600aa FPR 1.2% | 3/255 = **1.2%** | ✓ | |
| 2e | 600–1000aa FPR 0.0% | 0/180 = **0.0%** | ✓ | Bands cover all 1,198 (max len 1000, none >1000). |
| 3 | FZD5 (5.69/0.660) > ~84% reprod, ~99% PAE, is a hit | reprod > **205/243 = 84.4%**; PAE > **240/243 = 98.8%**; hit=True | ✓ | See honesty note — 84% is modest, not strong. |
| 4 | ULBP2 (5.74/0.893) more reprod than 242/243, is a hit | reprod > **242/243 = 99.6%**; hit=True | ✓ | |
| 5 | VEGFR2 is a non-hit (missed) | cofold-wt **11.55/0.43 → non-hit** | ✓ | Injected from `cofold_metrics.csv` (cofold-wt), **not present in the flagship screen CSV**; summary discloses this. |
| 6 | Anchor-2: PF4 hit; top ~9% reprod, #1 PAE of 76 same-size; FPR 0.4% | PF4 2.19/0.725 hit=True; reprod > **69/76 = 90.8%** (top 9.2%); PAE **#1/76 (100%)**; anchor2 ≥150aa FPR **1/227 = 0.44%** | ✓ | PF4 = 70aa (P02776 mature). Fully reproduces. |
| 7 | Scramble6: FZD5→non-hit (9.94/0.26), ULBP2→non-hit (11.19/0.82) | FZD5 **9.939/0.26 → non-hit**; ULBP2 **11.185/0.82 → non-hit** | ✓ | Partial-scramble rows (`scrambled`) confirm the "6-CDR needed for FZD5" note: FZD5 partial = 7.02/0.895 → **still a hit**; ULBP2 partial = 10.367/0.506 → non-hit. |
| 8 | Over-docking PAE<5.69: **168** total, median ~41aa, only 3 ≥150aa | **171** total, median **41aa**, **3** ≥150aa | ✗ | **Count is 171, not 168** — `checks.py` itself prints 171. Median and the ">=150aa = 3" reproduce. The physically-implausible finding stands; the stated count is just wrong. |
| 9 | Counts: 1,198/1,200 flagship, 401/401 anchor2, total 1,603 | flagship **1,198**, anchor2 **401**, +2 scramble6 +2 scrambled = **1,603** | ✓ / ⚠ | Scored-row counts all correct. But see honesty note: the summary's caveat says "1,198 + 401 + **20 pilot** + 4 scramble (1,603 total)" which sums to **1,623**, not 1,603 — the 1,603 excludes pilots. Also `cofold_metrics.csv` has **24** rows, not 20. Bookkeeping only. |
| 10 | Top valid-regime candidates incl. SMO + Ig-fold CD19/IL22RA1/KIRs | SMO (202aa, #1 by PAE), IL22RA1 (213), CD19 (272), **all 5 named KIRs** (2DL4/2DS1/3DL3/2DL3/2DL5B) present | ✓ | Exact gene list matches. |

**Cutoff sensitivity (flagship valid-regime FPR):** ≥130aa → 2.55% (18/707); ≥150aa →
2.06% (14/678); ≥170aa → 2.14% (14/655). **The 2.1% headline is robust to the cutoff** —
it stays ~2–2.5% across nearby thresholds and does not hinge on drawing the line at exactly 150.

## Scientific honesty section

**1. "FZD5 above 84% on reproducibility" is MODEST, not strong.** FZD5 is a *known* clinical
off-target, yet **38 of 243 same-size random decoys (15.6%) out-reproduce it** on the
`epitope_reprod` axis. The sensitivity story for FZD5 rests almost entirely on **PAE** (only
3/243 decoys dock tighter, 98.8%), not on reproducibility. ULBP2 is the genuinely strong case
(reprod beats 242/243). Pitching "84%" as a strong sensitivity result overstates the reprod
axis for FZD5 — state it as "top by PAE, mid-upper by reproducibility."

**2. Candidate off-targets are dominated by one fold class (Ig-fold) — plausible artifact.**
Of the 14 valid-regime hits, **5 are KIRs (36%)** and **≥7 are Ig-superfamily immune receptors
(50%: 5 KIRs + CD19 + IL22RA1; HLA-DRB3 is also an MHC/Ig-fold)**. A "novel off-target list"
that is half one structural family, matching the antibody's own Ig-fold, is exactly the
signature of a **fold/affinity bias rather than diverse biological discovery**. The summary
does flag this caveat honestly, but the concentration is high enough that these should be
treated as fold-bias suspects until a wet-lab/structural read, not as findings. Several hits
also sit near the decision boundary (HLA-DRB3 PAE 9.78, ADAM11 9.57, MRGPRD reprod 0.597,
KIR2DL3 0.615), so the count of 14 is threshold-fragile.

**3. Small-antigen flaw is real and correctly characterized** — 171 (not 168) decoys dock
tighter than the true off-target FZD5, median 41aa, essentially all short peptides (only 3 ≥150aa).
This is the honest self-identified limitation; the arithmetic error doesn't undermine it.

**4. Internal contradictions in `validation_summary.md`:**
- **Garbled "Bottom line" (lines 90–91):** a duplicated/contradictory parenthetical — reads
  "...medium/large antigens (≥150 aa: **~3% FPR**, both known off-targets enriched), is
  paratope-specific, (≥150 aa: **~2% FPR**, both known off-targets enriched) and generalizes...".
  Two different FPRs (~3% vs ~2%) for the same regime, one clause clearly a leftover paste. The
  correct value is **2.06%**; the "~3%" is loose/wrong and should be deleted.
- **Head-count that doesn't sum (line 105):** "1,198 flagship + 401 anchor-2 + 20 pilot + 4
  scramble (1,603 total)" — 1,198+401+20+4 = **1,623**. The real 1,603 = flagship + anchor2 +
  4 scramble; the 20 (really **24**) pilot cofolds in `cofold_metrics.csv` are **not** in that
  total. Fix the phrasing so the pilots aren't double-implied.

**5. VEGFR2 provenance (claim 5):** correctly a non-hit, but note it was never in the 1,198-protein
flagship screen — it is imported from `cofold_metrics.csv` (cofold-wt). The summary discloses this,
so it's honest, but the "screen missed VEGFR2" framing should not imply it went through the screen.

## Bottom line for the pitch
- **Safe to amplify:** ≥150aa FPR ~2% (cutoff-robust), the size-banded FPR table, PF4
  generalization (#1 by PAE, top ~9% reprod, 0.4% FPR), the negative-control collapse, the
  small-antigen flaw, and the SMO/CD19/IL22RA1/KIR candidate list.
- **Fix before pitching:** change "168" → **171**; drop the "~3% FPR" fragment in the Bottom
  line; repair the 1,603 vs 1,623 head-count wording.
- **Do not overclaim:** FZD5 sensitivity is PAE-driven (reprod only 84th pct — 38 decoys beat
  a known off-target); the candidate list is 50% Ig-fold and may be a fold/affinity artifact.
