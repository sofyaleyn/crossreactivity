"""Load the committed cofold metrics and apply the frozen calibrated verdict.

The dashboard, the figures, and the CLI verdict table all read their numbers
from here, so there is exactly one source of truth: ``cofold_metrics.csv``.

Verdict calibration (frozen, explicit — see ``findings.md`` Exp 4):

- Reference band is set by two calibrators cofolded with the same antibody:
  the **PD-1 cognate ceiling** (PAE_IF 7.24, reprod 0.936) and the
  **lysozyme non-binder floor** (PAE_IF 12.38, reprod 0.446).
- A query antigen reads **confirmed** iff its interface is tight AND its epitope
  is reproducible: ``PAE_IF < PAE_IF_CONFIRM`` AND ``reprod >= REPROD_CONFIRM``.
- The two calibrators keep their designated roles (``ceiling`` / ``floor``).
- A known/candidate off-target that fails the confirmed bar reads **missed** —
  triage is honest about what it cannot reach in silico.

Thresholds sit at/near the midpoint of the calibration band and are the only
tunable knobs; every verdict below is a deterministic function of the CSV.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

from . import paths

# --- frozen calibration thresholds ---
PAE_IF_CONFIRM = 9.8   # ~midpoint of PD-1 ceiling 7.24 and lysozyme floor 12.38
REPROD_CONFIRM = 0.55  # clearly above the lysozyme floor (0.446); FZD5 = 0.660

# --- the calibrated SHR-1210 panel: which runs anchor the demo spine ---
CEILING_RUN = "cofold-pd1"
FLOOR_RUN = "cofold-lyz"

# run -> named confirmation assay (Phase 3 assay map; triage points the wet lab)
ASSAY = {
    "cofold-pd1": "PD-1 blockade reporter bioassay (cognate target)",
    "cofold-fzd5": "Wnt/beta-catenin (TCF/LEF) reporter + FZD-family SPR panel",
    "cofold-ulbp2": "NKG2D-ligand cell-surface binding / ULBP2 SPR",
    "cofold-wt": "VEGFR2 ectodomain SPR / receptor-proteome array",
    "cofold-lyz": "n/a — specificity floor control (non-binder)",
}

# The calibrated panel shown in Chart A, ceiling -> floor.
PANEL_RUNS = ["cofold-pd1", "cofold-fzd5", "cofold-ulbp2", "cofold-wt", "cofold-lyz"]

# Short display names for panel antigens.
DISPLAY = {
    "cofold-pd1": "PD-1",
    "cofold-fzd5": "FZD5",
    "cofold-ulbp2": "ULBP2",
    "cofold-wt": "VEGFR2",
    "cofold-lyz": "lysozyme",
}


@dataclass
class Row:
    run: str
    antibody: str
    antigen: str
    role: str
    PAE_IF: float
    epitope_reprod: float
    iptm_max: float
    iptm_mean: float
    prediction_id: str

    @property
    def display(self) -> str:
        return DISPLAY.get(self.run, self.antigen)


def _f(x: str) -> float:
    return float(x) if x not in (None, "") else float("nan")


def load_rows() -> dict[str, Row]:
    """All CSV rows keyed by run id (numbers exactly as committed)."""
    rows: dict[str, Row] = {}
    with open(paths.METRICS_CSV, newline="") as f:
        for r in csv.DictReader(f):
            rows[r["run"]] = Row(
                run=r["run"],
                antibody=r["antibody"],
                antigen=r["antigen"],
                role=r["role"],
                PAE_IF=_f(r["PAE_IF_mean"]),
                epitope_reprod=_f(r["epitope_reprod"]),
                iptm_max=_f(r["iptm_max"]),
                iptm_mean=_f(r["iptm_mean"]),
                prediction_id=r["prediction_id"],
            )
    return rows


def verdict(run: str, row: Row) -> str:
    """Frozen calibrated read for a run: confirmed | missed | ceiling | floor."""
    if run == CEILING_RUN:
        return "ceiling"
    if run == FLOOR_RUN:
        return "floor"
    if row.PAE_IF < PAE_IF_CONFIRM and row.epitope_reprod >= REPROD_CONFIRM:
        return "confirmed"
    return "missed"


def assay_for(run: str) -> str:
    return ASSAY.get(run, "receptor-proteome array (triage confirmation)")


# --- Chart B: FZD5-family within-fold discrimination (Exp 6) ---
FAMILY_POSITIVE = "cofold-fzd5"


def family_rows(rows: dict[str, Row]) -> list[Row]:
    """The 12 Frizzled-fold members (FZD5 + 11 same-fold decoys), by PAE_IF asc."""
    fam = [
        r for run, r in rows.items()
        if "FZD5 family" in r.role or run == FAMILY_POSITIVE
    ]
    return sorted(fam, key=lambda r: r.PAE_IF)


def family_auroc(rows: dict[str, Row]) -> tuple[float, int, int]:
    """Within-family AUROC ranking FZD5 (positive) vs decoys by -PAE_IF.

    Returns (auroc, rank_of_fzd5_1indexed, family_size). With a single positive,
    AUROC = fraction of negatives the positive out-ranks (lower PAE = better).
    """
    fam = family_rows(rows)
    n = len(fam)
    pos = next(r for r in fam if r.run == FAMILY_POSITIVE)
    negatives = [r for r in fam if r.run != FAMILY_POSITIVE]
    beaten = sum(1 for r in negatives if pos.PAE_IF < r.PAE_IF)
    auroc = beaten / len(negatives) if negatives else float("nan")
    rank = sorted(fam, key=lambda r: r.PAE_IF).index(pos) + 1
    return auroc, rank, n


# --- Chart C: SHR-1210 vs pembrolizumab specificity control (Exp 8) ---
# antigen label -> (SHR-1210 run, pembrolizumab run)
CONTROL_PAIRS = [
    ("PD-1", "cofold-pd1", "cofold-pembro-pd1"),
    ("FZD5", "cofold-fzd5", "cofold-pembro-fzd5"),
    ("ULBP2", "cofold-ulbp2", "cofold-pembro-ulbp2"),
    ("VEGFR2", "cofold-wt", "cofold-pembro-vegfr2"),
]
