"""Scaled, budget-guarded cofold screen controller (demo-plan.md Phase -1, Steps D-F).

Config-driven. A "job" is one antibody (VH/VL) x a list of antigen sequences.
Reuses the FROZEN rig (crossflag.demo.scoring + panel) — no re-tuning.

Central budget guard: a single cumulative ledger (includes the $3.90 pilot).
No batch is submitted if cumulative estimate would cross HARD_STOP ($480).
Everything is resume-safe (idempotency keys; skip already-downloaded/scored).

Usage:
  python screen.py build    <job>        # write input JSONs from jobs/<job>.json
  python screen.py estimate <job>        # FREE estimate-cost; update ledger; gate
  python screen.py submit   <job>        # start() all; record prediction ids
  python screen.py download <job>        # wait + download all completed archives
  python screen.py collect  <job>        # normalize + score through frozen rig -> screen_metrics.csv
"""
from __future__ import annotations

import json
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCREEN = ROOT / "data" / "results" / "screen"
JOBS = SCREEN / "jobs"
INPUTS = SCREEN / "screen_inputs"
RUNS = SCREEN / "screen_runs"
STRUCTURES = ROOT / "data" / "results" / "structures"
DEMO_INPUTS = ROOT / "data" / "results" / "inputs"
LEDGER = SCREEN / "spend_ledger.json"
METRICS = SCREEN / "screen_metrics.csv"

MODEL = "boltz-2.1"
HARD_STOP_USD = 480.0   # $20 hard margin under the $500 ceiling
CEILING_USD = 500.0


# ---------- ledger (single source of truth for cumulative spend) ----------
def _ledger() -> dict:
    if LEDGER.exists():
        d = json.loads(LEDGER.read_text())
    else:
        d = {}
    d.setdefault("charges", [])          # [{job, n, estimate_usd}]
    # pilot spend is recorded as its own fields; fold it in as a charge once
    if not any(c.get("job") == "pilot" for c in d["charges"]) and d.get("spend_usd"):
        d["charges"].insert(0, {"job": "pilot", "n": d.get("n_succeeded", 20),
                                "estimate_usd": d["spend_usd"]})
    return d


def _committed(d: dict) -> float:
    return round(sum(c["estimate_usd"] for c in d["charges"]), 4)


def _save(d: dict) -> None:
    d["committed_estimate_usd"] = _committed(d)
    d["hard_stop_usd"] = HARD_STOP_USD
    d["ceiling_usd"] = CEILING_USD
    d["remaining_to_hardstop_usd"] = round(HARD_STOP_USD - _committed(d), 2)
    LEDGER.write_text(json.dumps(d, indent=2))


# ---------- job config ----------
def _job(job: str) -> dict:
    return json.loads((JOBS / f"{job}.json").read_text())


def build(job: str) -> None:
    cfg = _job(job)
    H, L = cfg["antibody"]["H"], cfg["antibody"]["L"]
    INPUTS.mkdir(parents=True, exist_ok=True)
    for a in cfg["antigens"]:
        inp = {
            "entities": [
                {"type": "protein", "chain_ids": ["H"], "value": H},
                {"type": "protein", "chain_ids": ["L"], "value": L},
                {"type": "protein", "chain_ids": ["V"], "value": a["seq"]},
            ],
            "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["H", "L"]},
            "num_samples": 5,
            "model_options": {"recycling_steps": 3, "sampling_steps": 200},
        }
        (INPUTS / f"{a['run']}.json").write_text(json.dumps(inp, indent=2))
    print(f"[{job}] built {len(cfg['antigens'])} inputs")


def _estimate_one(run: str) -> float:
    r = subprocess.run(
        ["boltz-api", "predictions:structure-and-binding", "estimate-cost",
         "--input", f"@json://{INPUTS / (run + '.json')}", "--model", MODEL, "--format", "json"],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"estimate-cost {run}: {r.stderr[:200]}")
    return float(json.loads(r.stdout)["estimated_cost_usd"])


def estimate(job: str) -> None:
    from concurrent.futures import ThreadPoolExecutor
    cfg = _job(job)
    d = _ledger()
    runs = [a["run"] for a in cfg["antigens"]]
    with ThreadPoolExecutor(max_workers=16) as ex:
        costs = list(ex.map(_estimate_one, runs))
    total = round(sum(costs), 4)
    committed = _committed(d)
    print(f"[{job}] batch estimate ${total:.2f}  | already committed ${committed:.2f}"
          f"  | projected ${committed + total:.2f}  | hard-stop ${HARD_STOP_USD:.0f}")
    if committed + total > HARD_STOP_USD:
        print(f"[{job}] GATE FAIL: would cross hard-stop. NOT recording; shrink the job.")
        sys.exit(2)
    d["charges"] = [c for c in d["charges"] if c.get("job") != job]
    d["charges"].append({"job": job, "n": len(cfg["antigens"]), "estimate_usd": total})
    _save(d)
    print(f"[{job}] GATE PASS: recorded. remaining to hard-stop "
          f"${HARD_STOP_USD - _committed(d):.2f}")


def _start_one(run: str):
    import random
    import time
    for attempt in range(6):
        r = subprocess.run(
            ["boltz-api", "predictions:structure-and-binding", "start",
             "--input", f"@json://{INPUTS / (run + '.json')}", "--model", MODEL,
             "--idempotency-key", run, "--format", "json"],
            capture_output=True, text=True)
        if r.returncode == 0:
            return run, json.loads(r.stdout)["id"], None
        if "429" in r.stderr or "Too Many Requests" in r.stderr:
            time.sleep(min(2 ** attempt, 20) + random.random())
            continue
        return run, None, r.stderr[:160]
    return run, None, "429 after retries"


def submit(job: str) -> None:
    from concurrent.futures import ThreadPoolExecutor
    cfg = _job(job)
    sub_file = SCREEN / f"submitted_{job}.json"
    submitted = json.loads(sub_file.read_text()) if sub_file.exists() else {}
    todo = [a["run"] for a in cfg["antigens"] if a["run"] not in submitted]
    fails = 0
    with ThreadPoolExecutor(max_workers=5) as ex:
        for i, (run, pid, err) in enumerate(ex.map(_start_one, todo), 1):
            if pid is None:
                fails += 1; print(f"  SUBMIT FAIL {run}: {err}")
            else:
                submitted[run] = pid
            if i % 100 == 0:
                sub_file.write_text(json.dumps(submitted, indent=2))
    sub_file.write_text(json.dumps(submitted, indent=2))
    print(f"[{job}] submitted {len(submitted)}/{len(cfg['antigens'])} (fails {fails})")


def _status(pid: str) -> str:
    r = subprocess.run(
        ["boltz-api", "predictions:structure-and-binding", "retrieve", "--id", pid,
         "--format", "json"], capture_output=True, text=True)
    if r.returncode != 0:
        return "unknown"
    return json.loads(r.stdout).get("status", "unknown")


def _dl_one(run_pid):
    run, pid = run_pid
    outdir = RUNS / run
    if (outdir / "outputs" / "archive.tar.gz").exists() or (STRUCTURES / run / "sample_0.cif").exists():
        return run, "have"
    st = _status(pid)
    if st not in ("succeeded", "completed", "success"):
        return run, st  # not ready (pending/running/failed) — skip this pass
    r = subprocess.run(
        ["boltz-api", "download-results", "--id", pid, "--run-dir", str(outdir),
         "--poll-interval-seconds", "15"], capture_output=True, text=True)
    return run, ("ok" if r.returncode == 0 else f"dlfail:{r.stderr[:60]}")


def download(job: str) -> None:
    from collections import Counter
    from concurrent.futures import ThreadPoolExecutor
    submitted = json.loads((SCREEN / f"submitted_{job}.json").read_text())
    RUNS.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=10) as ex:
        res = list(ex.map(_dl_one, submitted.items()))
    c = Counter(s for _, s in res)
    have = c.get("have", 0) + c.get("ok", 0)
    print(f"[{job}] ready+downloaded {have}/{len(submitted)}  states={dict(c)}")
    return have, len(submitted)


def _normalize(run: str) -> bool:
    sdir = STRUCTURES / run
    dst_in = DEMO_INPUTS / f"{run.replace('-', '_')}.json"
    # already normalized?
    if (sdir / "sample_0.cif").exists() and (sdir / "sample_4_pae.npz").exists():
        if not dst_in.exists():
            dst_in.write_bytes((INPUTS / f"{run}.json").read_bytes())
        return True
    arc = RUNS / run / "outputs" / "archive.tar.gz"
    if not arc.exists():
        return False
    sdir.mkdir(parents=True, exist_ok=True)
    ex = RUNS / run / "outputs" / "_ex"
    with tarfile.open(arc) as t:
        t.extractall(ex)
    pred = ex / "prediction"
    for i in range(5):
        c = pred / f"sample_{i}_predicted_structure.cif"
        p = pred / f"sample_{i}_pae.npz"
        if c.exists():
            (sdir / f"sample_{i}.cif").write_bytes(c.read_bytes())
        if p.exists():
            (sdir / f"sample_{i}_pae.npz").write_bytes(p.read_bytes())
    (DEMO_INPUTS / f"{run.replace('-', '_')}.json").write_bytes(
        (INPUTS / f"{run}.json").read_bytes())
    ok = any(sdir.glob("sample_*.cif"))
    if ok:  # normalized into structures/ — drop the redundant archive+extraction (disk safety)
        import shutil
        shutil.rmtree(RUNS / run, ignore_errors=True)
    return ok


def collect(job: str) -> None:
    import csv
    from crossflag.demo import panel, scoring

    cfg = _job(job)
    meta = {a["run"]: a for a in cfg["antigens"]}
    sub_file = SCREEN / f"submitted_{job}.json"
    submitted = json.loads(sub_file.read_text())
    already = set()
    if METRICS.exists():
        already = {r["run"] for r in csv.DictReader(METRICS.open())}
    P, R = panel.PAE_IF_CONFIRM, panel.REPROD_CONFIRM
    rows = []
    for run, pid in submitted.items():
        if run in already:
            continue
        if not _normalize(run):
            continue
        try:
            s = scoring.score_run(run)
        except Exception as e:  # noqa: BLE001
            print(f"  SCORE-ERR {run}: {e}"); continue
        m = meta[run]
        hit = (s.PAE_IF < P) and (s.epitope_reprod >= R)
        rows.append({"job": job, "run": run, "gene": m.get("gene", ""),
                     "uniprot_id": m.get("uniprot_id", ""), "antigen_len": len(m["seq"]),
                     "PAE_IF": round(s.PAE_IF, 3), "epitope_reprod": round(s.epitope_reprod, 3),
                     "hit": int(hit), "prediction_id": pid})
    rows.sort(key=lambda r: r["PAE_IF"])
    # append/update METRICS csv (keyed by run)
    existing = {}
    if METRICS.exists():
        for r in csv.DictReader(METRICS.open()):
            existing[r["run"]] = r
    for r in rows:
        existing[r["run"]] = r
    cols = ["job", "run", "gene", "uniprot_id", "antigen_len", "PAE_IF",
            "epitope_reprod", "hit", "prediction_id"]
    with METRICS.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(existing.values(), key=lambda r: float(r["PAE_IF"])):
            w.writerow(r)
    hits = [r for r in rows if r["hit"]]
    print(f"[{job}] scored {len(rows)}; hits {len(hits)}"
          + (f" ({100*len(hits)/len(rows):.1f}%)" if rows else ""))
    for r in hits:
        print(f"    HIT  {r['gene']:12s} PAE_IF {r['PAE_IF']:.2f}  reprod {r['epitope_reprod']:.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: screen.py <build|estimate|submit|download|collect> <job>")
    {"build": build, "estimate": estimate, "submit": submit,
     "download": download, "collect": collect}[sys.argv[1]](sys.argv[2])
