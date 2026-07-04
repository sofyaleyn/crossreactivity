"""CrossFlag verdict-table CLI (Phase 3) — the demo as a replayable tool.

Input: the SHR-1210 VH/VL anchor sequences (``data/anchor/``) + a list of
reference cofold runs. Output: the ranked verdict table (antigen, PAE_IF,
epitope_reprod, read, confirmation assay) printed to stdout AND written to
``demo/verdict_table.{md,json}``.

The default path replays committed data only — fully offline and deterministic.
A ``--live`` flag is accepted for parity with the productization story but, with
no network/auth in the demo environment, it prints a clear message and falls back
to the committed path rather than hard-failing.
"""
from __future__ import annotations

import argparse
import json
import sys

from . import panel, paths

# Reference runs shown in the verdict table, ceiling -> floor.
DEFAULT_RUNS = panel.PANEL_RUNS

READ = {
    "ceiling": "ceiling (cognate)",
    "confirmed": "confirmed",
    "missed": "missed (weakest)",
    "floor": "floor (non-binder)",
}


def build_table(runs: list[str] | None = None) -> list[dict]:
    rows = panel.load_rows()
    runs = runs or DEFAULT_RUNS
    table = []
    for run in runs:
        r = rows[run]
        v = panel.verdict(run, r)
        table.append({
            "run": run,
            "antigen": r.display,
            "PAE_IF": round(r.PAE_IF, 2),
            "epitope_reprod": round(r.epitope_reprod, 3),
            "verdict": v,
            "read": READ[v],
            "confirmation_assay": panel.assay_for(run),
            "prediction_id": r.prediction_id,
        })
    # rank by PAE_IF (tightest first) for the "triage" ordering
    return sorted(table, key=lambda t: t["PAE_IF"])


def render_md(table: list[dict]) -> str:
    header = (
        "# CrossFlag verdict table — SHR-1210 (camrelizumab)\n\n"
        "Antibody: SHR-1210 VH/VL (`data/anchor/`). Numbers replayed from "
        "committed `data/results/cofold_metrics.csv`. Triage, not certification.\n\n"
        "| antigen | PAE_IF (Å) | epitope_reprod | read | confirmation assay |\n"
        "|---|---|---|---|---|\n"
    )
    lines = [
        f"| {t['antigen']} | {t['PAE_IF']:.2f} | {t['epitope_reprod']:.3f} "
        f"| {t['read']} | {t['confirmation_assay']} |"
        for t in table
    ]
    return header + "\n".join(lines) + "\n"


def render_stdout(table: list[dict]) -> str:
    w = f"{'antigen':<10}{'PAE_IF':>8}{'reprod':>9}  {'read':<20}{'assay'}"
    out = [w, "-" * 96]
    for t in table:
        out.append(
            f"{t['antigen']:<10}{t['PAE_IF']:>8.2f}{t['epitope_reprod']:>9.3f}  "
            f"{t['read']:<20}{t['confirmation_assay']}"
        )
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="CrossFlag off-target verdict table (offline replay).")
    p.add_argument("--live", action="store_true",
                   help="attempt a live boltz.bio re-run (falls back to committed data).")
    p.add_argument("--no-write", action="store_true", help="print only; do not write demo/ files.")
    args = p.parse_args(argv)

    if args.live:
        print("[--live] no network/auth in the demo env → replaying committed data.\n",
              file=sys.stderr)

    table = build_table()
    print(render_stdout(table))

    if not args.no_write:
        paths.DEMO.mkdir(parents=True, exist_ok=True)
        paths.VERDICT_JSON.write_text(json.dumps(table, indent=2) + "\n")
        paths.VERDICT_MD.write_text(render_md(table))
        print(f"\nwrote {paths.VERDICT_JSON.relative_to(paths.ROOT)} and "
              f"{paths.VERDICT_MD.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
