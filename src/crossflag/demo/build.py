"""One-command demo build (Phase 5).

``python -m crossflag.demo.build`` regenerates, from committed data alone, with
no network and no GPU:

  demo/figures/*.png        (Phase 1 charts + Phase 2 interface figures)
  demo/verdict_table.{md,json}   (Phase 3 replayable tool output)
  demo/dashboard.html       (Phase 1 self-contained dashboard)

Deterministic: two runs produce byte-identical artifacts. Exits 0 on success.
"""
from __future__ import annotations

import sys

from . import dashboard, figures, paths, run


def main() -> int:
    print("CrossFlag demo build — offline, committed data only\n")

    print("[1/3] figures ...")
    for name in figures.build_all():
        print(f"      demo/figures/{name}")

    print("[2/3] verdict table ...")
    table = run.build_table()
    paths.DEMO.mkdir(parents=True, exist_ok=True)
    import json
    paths.VERDICT_JSON.write_text(json.dumps(table, indent=2) + "\n")
    paths.VERDICT_MD.write_text(run.render_md(table))
    print(f"      {paths.VERDICT_JSON.relative_to(paths.ROOT)}")
    print(f"      {paths.VERDICT_MD.relative_to(paths.ROOT)}")

    print("[3/3] dashboard ...")
    rel = dashboard.build()
    print(f"      {rel}")

    print("\nOK — open demo/dashboard.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
