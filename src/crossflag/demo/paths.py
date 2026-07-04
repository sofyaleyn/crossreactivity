"""Repo-root and committed-data path resolution for the CrossFlag demo.

Everything the demo reads lives under the committed ``data/`` tree; everything it
writes lives under ``demo/``. No scratchpad, no network, no GPU.
"""
from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    """Walk up from this file until we find the committed data tree."""
    here = Path(__file__).resolve()
    for parent in [here, *here.parents]:
        if (parent / "data" / "results" / "cofold_metrics.csv").exists():
            return parent
    raise FileNotFoundError(
        "Could not locate repo root (no data/results/cofold_metrics.csv found "
        f"walking up from {here})."
    )


ROOT = repo_root()

# --- committed inputs (read-only) ---
DATA = ROOT / "data"
RESULTS = DATA / "results"
METRICS_CSV = RESULTS / "cofold_metrics.csv"
STRUCTURES = RESULTS / "structures"
INPUTS = RESULTS / "inputs"
ANCHOR = DATA / "anchor"

# --- demo outputs (written) ---
DEMO = ROOT / "demo"
FIGURES = DEMO / "figures"
DASHBOARD = DEMO / "dashboard.html"
VERDICT_JSON = DEMO / "verdict_table.json"
VERDICT_MD = DEMO / "verdict_table.md"


def input_json_for(run: str) -> Path:
    """``cofold-fzd5`` -> ``data/results/inputs/cofold_fzd5.json``."""
    return INPUTS / (run.replace("-", "_") + ".json")


def structure_dir_for(run: str) -> Path:
    return STRUCTURES / run
