# Experimental results — raw data (2026-07-04)

Durable record of the viability-study results (scratchpad is ephemeral). Narrative + interpretation in [`../../findings.md`](../../findings.md). All cofolds: **boltz.bio hosted API, Boltz-2.1, 5 samples each**, `protein_protein_binding` (antibody = binder).

## Files

| file | contents |
|---|---|
| `cofold_metrics.csv` | All 24 cofold runs × full metrics (ipTM max/mean, protein-ipTM, pTM, complex-pLDDT, structure-confidence, binding-confidence, PAE_IF, epitope-reprod, prediction_id). |
| `tier1_filter_scores.csv` | 2,896 curated proteins × 7 cheap-filter scores + `is_positive` (Exp 9). |
| `benchmark_antigens.json` | Fold-matched decoy benchmark: 40 antigens × {family, name, is_pos, sequence} (Exp 5/7). |
| `benchmark_embed_scores.json` | Embedding paratope-cosine scores on the benchmark (Exp 5). |
| `benchmark_surface_scores.json` | Surface-method scores on the benchmark (Exp 7). |
| `structures/` | Predicted-structure CIFs for all 24 cofold runs (5 samples each) — see `structures/README.md`. |
| `inputs/` | Exact Boltz-2.1 input JSONs for every run (sequences, chains, templates, sampling params). |
| `scripts/` | Raw session analysis scripts (metric extraction, cofold-input builders, benchmark/filters) — to be refactored into `src/crossflag/` per `plan.md`. |

## Metric definitions

- **ipTM_max / ipTM_mean** — Boltz interface pTM, best-of-5 and mean over the 5 samples. *Found uninformative for calling binders* (over-docks; lysozyme non-binder scores 0.87). Recorded for completeness.
- **PAE_IF_mean** — mean predicted aligned error over all antibody(H+L)↔antigen(V) token pairs, averaged across samples (Å; lower = tighter interface).
- **epitope_reprod** — mean pairwise Jaccard of the contacted antigen residues (5 Å heavy-atom) across the 5 samples (the **primary discriminating metric**; higher = same epitope every sample).
- **binding_confidence** — Boltz binding head (≈0 for all protein–protein here; not usable).

## Cofold metrics — all runs

### Calibrated panel (Exp 4) — SHR-1210 WT, PD-1 ceiling / lysozyme floor
| run | antigen | ipTM_max | ipTM_mean | PAE_IF | epitope_reprod | read |
|---|---|---|---|---|---|---|
| cofold-pd1 | PD-1 (cognate, ceiling) | 0.910 | 0.887 | 7.24 | **0.936** | binder |
| cofold-fzd5 | FZD5 CRD | 0.952 | 0.915 | 5.69 | **0.660** | confirmed |
| cofold-ulbp2 | ULBP2 | 0.940 | 0.916 | 5.74 | **0.893** | confirmed |
| cofold-wt | VEGFR2 D2-3 | 0.883 | 0.778 | 11.55 | 0.430 | floor |
| cofold-mut | VEGFR2 (germlined) | 0.852 | 0.753 | 13.41 | 0.512 | floor |
| cofold-wt-tmpl | VEGFR2 +3V2A template | 0.766 | 0.647 | 17.62 | 0.302 | worse when pinned |
| cofold-mut-tmpl | VEGFR2 germ +template | 0.832 | 0.650 | 17.87 | 0.068 | floor |
| cofold-lyz | lysozyme (non-binder, floor) | 0.869 | 0.781 | 12.38 | 0.446 | non-binder |

### FZD5-family within-fold discrimination (Exp 6) — SHR-1210 WT × Frizzled CRDs
| run | ipTM_max | ipTM_mean | PAE_IF | epitope_reprod |
|---|---|---|---|---|
| cofold-smo (SMO) | 0.968 | 0.928 | 4.98 | 0.652 |
| cofold-fzd5 (**FZD5, true off-target**) | 0.952 | 0.915 | 5.69 | 0.660 |
| cofold-fzd1 | 0.948 | 0.872 | 7.66 | 0.885 |
| cofold-fzd10 | 0.902 | 0.868 | 8.24 | 0.756 |
| cofold-fzd8 | 0.898 | 0.820 | 10.32 | 0.707 |
| cofold-fzd2 | 0.894 | 0.832 | 10.36 | 0.607 |
| cofold-sfrp1 | 0.894 | 0.792 | 12.05 | 0.266 |
| cofold-fzd7 | 0.889 | 0.807 | 11.31 | 0.299 |
| cofold-fzd9 | 0.862 | 0.813 | 11.03 | 0.513 |
| cofold-fzd4 | 0.949 | 0.802 | 11.93 | 0.509 |
| cofold-fzd3 | 0.841 | 0.740 | 15.41 | 0.232 |
| cofold-fzd6 | 0.825 | 0.738 | 16.06 | 0.289 |

FZD5 PAE_IF rank 2/12 (SMO the one within-fold false positive); PAE-based within-family AUROC = **0.909**.

### Pembrolizumab antibody-side specificity control (Exp 8)
| run | antigen | ipTM_max | ipTM_mean | PAE_IF | epitope_reprod |
|---|---|---|---|---|---|
| cofold-pembro-pd1 | PD-1 (shared target) | 0.963 | 0.958 | 3.61 | 0.968 |
| cofold-pembro-vegfr2 | VEGFR2 D2-3 | 0.762 | 0.692 | 16.46 | 0.642 |
| cofold-pembro-fzd5 | FZD5 CRD | 0.755 | 0.678 | 18.55 | 0.217 |
| cofold-pembro-ulbp2 | ULBP2 | 0.901 | 0.741 | 13.87 | 0.341 |

pembrolizumab binds PD-1 but collapses on FZD5/ULBP2 (vs SHR-1210's 0.660/0.893) → off-target confirmations are SHR-1210-CDR-specific.

### Pipeline smoke test
| run | antigen | ipTM_max |
|---|---|---|
| smoke-trpcage-aspirin | aspirin ligand | 0.457 |

Full numeric detail (protein-ipTM, pTM, pLDDT, structure-confidence, binding-confidence, prediction IDs) in `cofold_metrics.csv`.
