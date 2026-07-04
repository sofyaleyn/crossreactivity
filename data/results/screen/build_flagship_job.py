"""Build jobs/flagship.json for item 1 (real-scale screen) from ectodomains.csv.

- Antibody = SHR-1210 WT VH/VL (from the validated cofold_fzd5.json).
- Antigens = stratified sample of the reference set's ectodomains.
- Excludes the 3 anchor off-targets (VEGFR2/FZD5/ULBP2) — already cofolded in
  committed data; they are injected into the ROC as positives at analysis time,
  not re-folded here.
- Stratified by ectodomain length so enrichment is not a size artifact.

Usage: python build_flagship_job.py [TARGET_N]   (default 1500)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
ECTO = ROOT / "data" / "reference" / "ectodomains.csv"
FZD5 = ROOT / "data" / "results" / "inputs" / "cofold_fzd5.json"
JOBS = ROOT / "data" / "results" / "screen" / "jobs"

ANCHOR_OFFTARGETS = {"P35968", "Q13467", "Q9BZM5"}  # VEGFR2, FZD5, ULBP2
LEN_MIN, LEN_MAX = 30, 1000
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
N_STRATA = 12


def antibody():
    d = json.loads(FZD5.read_text())
    e = {x["chain_ids"][0]: x["value"] for x in d["entities"]}
    return e["H"], e["L"]


def main(target_n: int = 1500) -> None:
    rows = []
    for r in csv.DictReader(ECTO.open()):
        seq = (r.get("ectodomain_seq") or "").strip().upper()
        if not seq or not (LEN_MIN <= len(seq) <= LEN_MAX):
            continue
        if r.get("uniprot_id") in ANCHOR_OFFTARGETS:
            continue
        if any(c not in VALID_AA for c in seq):
            continue
        rows.append({"uniprot_id": r["uniprot_id"],
                     "gene": r.get("gene_symbol", ""), "seq": seq, "len": len(seq)})
    # de-dupe by uniprot
    seen, uniq = set(), []
    for r in rows:
        if r["uniprot_id"] in seen:
            continue
        seen.add(r["uniprot_id"]); uniq.append(r)
    uniq.sort(key=lambda r: (r["len"], r["uniprot_id"]))
    n = min(target_n, len(uniq))
    # length-stratified even sampling across the sorted-by-length list
    idx = sorted({round(i * (len(uniq) - 1) / (n - 1)) for i in range(n)})
    picks = [uniq[i] for i in idx]
    antigens = [{"run": f"screen-flagship-{p['uniprot_id']}", "uniprot_id": p["uniprot_id"],
                 "gene": p["gene"], "seq": p["seq"]} for p in picks]
    H, L = antibody()
    JOBS.mkdir(parents=True, exist_ok=True)
    (JOBS / "flagship.json").write_text(json.dumps(
        {"name": "flagship", "antibody": {"H": H, "L": L}, "antigens": antigens}, indent=2))
    lens = [p["len"] for p in picks]
    print(f"eligible ectodomains {LEN_MIN}-{LEN_MAX}aa: {len(uniq)}")
    print(f"selected {len(antigens)} (length-stratified)")
    print(f"  ecto_len min/median/max: {min(lens)}/{sorted(lens)[len(lens)//2]}/{max(lens)}")
    print(f"  <=400aa (cheap): {sum(l <= 400 for l in lens)}   400-1000aa: {sum(l > 400 for l in lens)}")
    print("wrote jobs/flagship.json")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1500)
