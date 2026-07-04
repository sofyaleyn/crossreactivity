"""Step A price-discovery pilot for the off-target screen (demo-plan.md).

Subcommands:
  select    -> pick 20 ectodomain-sized proteins (length-stratified), build Boltz inputs
  estimate  -> FREE estimate-cost on all inputs; gate: total must be <= $5
  submit    -> start (non-blocking) all not-yet-submitted inputs; record prediction ids
  collect   -> poll + download + normalize + score with the FROZEN rig

No re-tuning: scoring reuses crossflag.demo.scoring + panel (frozen thresholds).
Antibody = SHR-1210 WT VH/VL (H/L taken verbatim from the validated cofold_fzd5.json).
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCREEN = ROOT / "data" / "results" / "screen"
CSV = ROOT / "data" / "reference" / "self_proteins.csv"
FZD5_INPUT = ROOT / "data" / "results" / "inputs" / "cofold_fzd5.json"
INPUTS_DIR = SCREEN / "pilot_inputs"
LEDGER = SCREEN / "spend_ledger.json"
MANIFEST = SCREEN / "pilot_manifest.json"

N = 20
LEN_MIN, LEN_MAX = 80, 400  # ectodomain-sized regime (sidesteps unbuilt Step B)
COST_CAP_USD = 5.0
MODEL = "boltz-2.1"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def antibody_chains() -> tuple[str, str]:
    d = json.loads(FZD5_INPUT.read_text())
    ents = {e["chain_ids"][0]: e["value"] for e in d["entities"]}
    return ents["H"], ents["L"]


def eligible() -> list[dict]:
    out = []
    for r in csv.DictReader(CSV.open()):
        seq = (r.get("sequence") or "").strip().upper()
        if not seq or not (LEN_MIN <= len(seq) <= LEN_MAX):
            continue
        if r.get("is_anchor_offtarget", "").strip().lower() == "true":
            continue
        if any(c not in VALID_AA for c in seq):
            continue
        out.append({"protein_id": r["protein_id"], "uniprot_id": r["uniprot_id"],
                    "gene": r["gene_symbol"], "name": r["name"], "seq": seq, "len": len(seq)})
    return out


def select() -> list[dict]:
    elig = sorted(eligible(), key=lambda d: (d["len"], d["protein_id"]))
    # length-stratified: 20 evenly-spaced picks across the size range -> price-vs-size curve
    idx = [round(i * (len(elig) - 1) / (N - 1)) for i in range(N)]
    picks = [elig[i] for i in sorted(set(idx))]
    H, L = antibody_chains()
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    manifest = []
    for p in picks:
        run = f"screen-pilot-{p['uniprot_id']}"
        inp = {
            "entities": [
                {"type": "protein", "chain_ids": ["H"], "value": H},
                {"type": "protein", "chain_ids": ["L"], "value": L},
                {"type": "protein", "chain_ids": ["V"], "value": p["seq"]},
            ],
            "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["H", "L"]},
            "num_samples": 5,
            "model_options": {"recycling_steps": 3, "sampling_steps": 200},
        }
        (INPUTS_DIR / f"{run}.json").write_text(json.dumps(inp, indent=2))
        manifest.append({"run": run, "uniprot_id": p["uniprot_id"], "gene": p["gene"],
                         "name": p["name"], "antigen_len": p["len"],
                         "idempotency_key": run})
    MANIFEST.write_text(json.dumps({"n_eligible": len(elig), "picks": manifest}, indent=2))
    print(f"eligible (len {LEN_MIN}-{LEN_MAX}aa): {len(elig)}")
    print(f"selected {len(manifest)} (length-stratified):")
    for m in manifest:
        print(f"  {m['antigen_len']:4d}aa  {m['gene']:12s} {m['uniprot_id']:8s} {m['name'][:44]}")
    return manifest


def _estimate_one(run: str) -> float:
    inp = INPUTS_DIR / f"{run}.json"
    r = subprocess.run(
        ["boltz-api", "predictions:structure-and-binding", "estimate-cost",
         "--input", f"@json://{inp}", "--model", MODEL, "--format", "json"],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"estimate-cost failed for {run}: {r.stderr[:300]}")
    return float(json.loads(r.stdout)["estimated_cost_usd"])


def estimate() -> None:
    man = json.loads(MANIFEST.read_text())["picks"]
    rows, total = [], 0.0
    for m in man:
        c = _estimate_one(m["run"])
        total += c
        rows.append({**m, "est_usd": c})
        print(f"  ${c:5.3f}  {m['antigen_len']:4d}aa  {m['gene']:12s} {m['uniprot_id']}")
    print(f"\nTOTAL estimated: ${total:.3f}  (cap ${COST_CAP_USD:.2f})")
    LEDGER.write_text(json.dumps(
        {"stage": "estimate", "cap_usd": COST_CAP_USD, "total_estimate_usd": round(total, 4),
         "per_run": rows}, indent=2))
    if total > COST_CAP_USD:
        print(f"ABORT: estimate ${total:.3f} exceeds cap ${COST_CAP_USD:.2f}")
        sys.exit(2)
    print("GATE PASSED: within cap.")


RUNS_DIR = SCREEN / "pilot_runs"
SUBMITTED = SCREEN / "pilot_submitted.json"


def submit() -> None:
    man = json.loads(MANIFEST.read_text())["picks"]
    submitted = json.loads(SUBMITTED.read_text()) if SUBMITTED.exists() else {}
    for m in man:
        run = m["run"]
        if run in submitted:
            print(f"  skip (already submitted): {run} -> {submitted[run]}")
            continue
        r = subprocess.run(
            ["boltz-api", "predictions:structure-and-binding", "start",
             "--input", f"@json://{INPUTS_DIR / (run + '.json')}",
             "--model", MODEL, "--idempotency-key", m["idempotency_key"], "--format", "json"],
            capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  SUBMIT FAILED {run}: {r.stderr[:200]}"); continue
        pid = json.loads(r.stdout)["id"]
        submitted[run] = pid
        print(f"  submitted {run} -> {pid}")
        SUBMITTED.write_text(json.dumps(submitted, indent=2))
    print(f"\n{len(submitted)}/{len(man)} submitted. ids in {SUBMITTED.name}")


def _status(pid: str) -> str:
    r = subprocess.run(
        ["boltz-api", "predictions:structure-and-binding", "retrieve", "--id", pid, "--format", "json"],
        capture_output=True, text=True)
    if r.returncode != 0:
        return f"error:{r.stderr[:80]}"
    return json.loads(r.stdout).get("status", "unknown")


def poll() -> None:
    submitted = json.loads(SUBMITTED.read_text())
    for run, pid in submitted.items():
        print(f"  {_status(pid):12s}  {run}  {pid}")


def download() -> None:
    """Download every completed prediction's raw archive to pilot_runs/<run>/."""
    submitted = json.loads(SUBMITTED.read_text())
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    for run, pid in submitted.items():
        outdir = RUNS_DIR / run
        if outdir.exists() and any(outdir.rglob("*.cif")):
            print(f"  have {run}"); continue
        r = subprocess.run(
            ["boltz-api", "download-results", "--id", pid, "--run-dir", str(outdir),
             "--poll-interval-seconds", "10"],
            capture_output=True, text=True)
        print(f"  {'ok ' if r.returncode == 0 else 'ERR'} {run}  {pid}"
              + ("" if r.returncode == 0 else f"  {r.stderr[:160]}"))


import tarfile

STRUCTURES = ROOT / "data" / "results" / "structures"
DEMO_INPUTS = ROOT / "data" / "results" / "inputs"
REPORT = SCREEN / "pilot_report.json"


def _normalize(run: str) -> bool:
    """Extract archive -> committed layout (sample_N.cif + sample_N_pae.npz + input json)."""
    arc = RUNS_DIR / run / "outputs" / "archive.tar.gz"
    if not arc.exists():
        return False
    sdir = STRUCTURES / run
    sdir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(arc) as t:
        t.extractall(RUNS_DIR / run / "outputs" / "_ex")
    pred = RUNS_DIR / run / "outputs" / "_ex" / "prediction"
    for i in range(5):
        cif = pred / f"sample_{i}_predicted_structure.cif"
        pae = pred / f"sample_{i}_pae.npz"
        if cif.exists():
            (sdir / f"sample_{i}.cif").write_bytes(cif.read_bytes())
        if pae.exists():
            (sdir / f"sample_{i}_pae.npz").write_bytes(pae.read_bytes())
    # input json in the layout paths.input_json_for() expects
    (DEMO_INPUTS / f"{run.replace('-', '_')}.json").write_bytes(
        (INPUTS_DIR / f"{run}.json").read_bytes())
    return any(sdir.glob("sample_*.cif"))


def collect() -> None:
    from crossflag.demo import panel, scoring  # frozen rig; no re-tuning

    man = {m["run"]: m for m in json.loads(MANIFEST.read_text())["picks"]}
    submitted = json.loads(SUBMITTED.read_text())
    P, R = panel.PAE_IF_CONFIRM, panel.REPROD_CONFIRM
    results = []
    for run, pid in submitted.items():
        if not _normalize(run):
            print(f"  PENDING/NO-DATA  {run}")
            results.append({**man[run], "pid": pid, "status": "no_data"})
            continue
        try:
            s = scoring.score_run(run)
        except Exception as e:  # noqa: BLE001
            print(f"  SCORE-ERR  {run}: {e}")
            results.append({**man[run], "pid": pid, "status": f"score_error:{e}"})
            continue
        hit = (s.PAE_IF < P) and (s.epitope_reprod >= R)
        results.append({**man[run], "pid": pid, "status": "ok",
                        "PAE_IF": round(s.PAE_IF, 3), "epitope_reprod": round(s.epitope_reprod, 3),
                        "hit": bool(hit)})

    ok = [r for r in results if r["status"] == "ok"]
    ok.sort(key=lambda r: r["PAE_IF"])
    print(f"\n=== PILOT RESULTS ({len(ok)}/{len(results)} scored) — frozen rig: "
          f"hit = PAE_IF<{P} AND reprod>={R} ===")
    print(f"{'gene':12s} {'len':>4s} {'PAE_IF':>7s} {'reprod':>7s}  {'verdict'}")
    for r in ok:
        v = "HIT" if r["hit"] else "-"
        print(f"{r['gene']:12s} {r['antigen_len']:4d} {r['PAE_IF']:7.2f} {r['epitope_reprod']:7.3f}  {v}")
    hits = [r for r in ok if r["hit"]]
    print(f"\nhits: {len(hits)}/{len(ok)}  ({100*len(hits)/len(ok):.0f}% flagged)"
          if ok else "\nno results scored yet")
    REPORT.write_text(json.dumps({"n": len(results), "scored": len(ok),
                                  "n_hits": len(hits), "results": results}, indent=2))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "select"
    {"select": select, "estimate": estimate, "submit": submit, "poll": poll,
     "download": download, "collect": collect}.get(cmd, lambda: sys.exit(f"unknown cmd {cmd!r}"))()
