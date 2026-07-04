# Predicted structures — Boltz-2.1 cofolds (2026-07-04)

Predicted-structure CIFs for every cofold run, `<run>/sample_N.cif` (5 samples/run; smoke=1). Antibody = chains H+L; antigen = chain V (ligand for the smoke test). Per-run confidence metrics (ipTM, PAE_IF, epitope-reprod, prediction IDs) in [`../cofold_metrics.csv`](../cofold_metrics.csv); exact Boltz inputs in [`../inputs/`](../inputs/); analysis scripts in [`../scripts/`](../scripts/).

| run | antibody | antigen | role | samples |
|---|---|---|---|---|
| `cofold-fzd1` | SHR-1210 WT | FZD1 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd10` | SHR-1210 WT | FZD10 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd2` | SHR-1210 WT | FZD2 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd3` | SHR-1210 WT | FZD3 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd4` | SHR-1210 WT | FZD4 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd5` | SHR-1210 WT | FZD5 CRD | off-target (confirmed) | 5 |
| `cofold-fzd6` | SHR-1210 WT | FZD6 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd7` | SHR-1210 WT | FZD7 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd8` | SHR-1210 WT | FZD8 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-fzd9` | SHR-1210 WT | FZD9 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-lyz` | SHR-1210 WT | lysozyme (HEL) | non-binder / floor | 5 |
| `cofold-mut-tmpl` | SHR-1210 germlined | VEGFR2 D2-3 +3V2A template | germlined, templated | 5 |
| `cofold-mut` | SHR-1210 germlined | VEGFR2 D2-3 | germlined variant | 5 |
| `cofold-pd1` | SHR-1210 WT | PD-1 (IgV) | cognate target / ceiling | 5 |
| `cofold-pembro-fzd5` | pembrolizumab | FZD5 CRD | Ab control | 5 |
| `cofold-pembro-pd1` | pembrolizumab | PD-1 (IgV) | Ab control - shared target | 5 |
| `cofold-pembro-ulbp2` | pembrolizumab | ULBP2 ectodomain | Ab control | 5 |
| `cofold-pembro-vegfr2` | pembrolizumab | VEGFR2 D2-3 | Ab control | 5 |
| `cofold-sfrp1` | SHR-1210 WT | SFRP1 CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-smo` | SHR-1210 WT | SMO CRD | within-fold decoy (FZD5 family) | 5 |
| `cofold-ulbp2` | SHR-1210 WT | ULBP2 ectodomain | off-target (confirmed) | 5 |
| `cofold-wt-tmpl` | SHR-1210 WT | VEGFR2 D2-3 +3V2A template | off-target, templated | 5 |
| `cofold-wt` | SHR-1210 WT | VEGFR2 D2-3 | off-target | 5 |
| `smoke-trpcage-aspirin` | (Trp-cage) | aspirin (ligand) | pipeline smoke test | 1 |
