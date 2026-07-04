# Anchor data — SHR-1210 (camrelizumab) validation case

MVP pipeline input. `pipeline.py` consumes the `variants/` directory (one VH/VL
FASTA pair per variant); `offtargets/` provides the ground-truth off-targets the
ranking is evaluated against.

## Files

| Path | Contents |
|---|---|
| `shr1210_vh.fasta` / `shr1210_vl.fasta` | Wild-type SHR-1210 Fv (Thera-SAbDab / US11208484B2). CDR-H3 `QLYYFDYW`, CDR-L3 `QQVYSIPWT`. |
| `variants/*.fasta` | Variant panel — WT + light-chain germline-reverted mutants. Each file has two records (VH, then VL). |
| `offtargets/` | VEGFR2, FZD5, ULBP2 sequences (UniProt) + `offtargets.csv` manifest with PDB refs. |

## Variant panel

Finlay et al. (2019) ablated the VEGFR2 off-target by **germlining the light
chain** while preserving PD-1 affinity — a ready-made WT-plus-mutants panel. The
exact Finlay mutant sequences are not public, so the germline revertants below
are **RECONSTRUCTED** from the IGKV1-39 / IGKJ germline and labelled as such (per
HANDOFF.md §Step 1 and Known risks). The VH is identical across all variants;
only the light-chain CDRs change.

| Variant | Light-chain change | Expected VEGFR2 risk |
|---|---|---|
| `shr1210_WT` | none (native) | **high** (VEGFR2 near top) |
| `shr1210_L1germ` | CDR-L1 `LASQTIGTWLT` → `RASQSISSYLN` | reduced |
| `shr1210_L3germ` | CDR-L3 `QQVYSIPWT` → `QQSYSTPWT` | reduced |
| `shr1210_L1L3germ` | both CDR-L1 + CDR-L3 reverted | **lowest** (VEGFR2 ablated) |

## Success criterion (test_panel_rank.py)

Run the pipeline on this panel: VEGFR2 must land in the top-3 flagged proteins
for `shr1210_WT`, and `risk_score(WT) > risk_score(shr1210_L1L3germ)`.
