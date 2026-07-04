"""Loop: download completed cofolds + score through the frozen rig, until all done.

Resumable and incremental (skips already-downloaded / already-scored runs).
Runs across all three jobs. Prints a progress line each pass. Exits when every
submitted run is scored, or after MAX_PASSES.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("screen", HERE / "screen.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

JOBS = ["flagship", "anchor2", "scrambled"]
MAX_PASSES = 200
SLEEP_S = 90


def n_submitted() -> int:
    t = 0
    for j in JOBS:
        f = S.SCREEN / f"submitted_{j}.json"
        if f.exists():
            t += len(json.loads(f.read_text()))
    return t


def n_scored() -> int:
    if not S.METRICS.exists():
        return 0
    return sum(1 for _ in csv.DictReader(S.METRICS.open()))


def main() -> None:
    target = n_submitted()
    for p in range(1, MAX_PASSES + 1):
        for j in JOBS:
            try:
                S.download(j)
            except Exception as e:  # noqa: BLE001
                print(f"  download {j} err: {e}")
        for j in JOBS:
            try:
                S.collect(j)
            except Exception as e:  # noqa: BLE001
                print(f"  collect {j} err: {e}")
        done = n_scored()
        target = max(target, n_submitted())
        print(f"=== pass {p}: scored {done}/{target} ===", flush=True)
        if done >= target and target > 0:
            print("ALL SCORED."); return
        time.sleep(SLEEP_S)
    print("MAX_PASSES reached.")


if __name__ == "__main__":
    main()
