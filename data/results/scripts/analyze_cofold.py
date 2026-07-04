import json, glob, os, tarfile, numpy as np

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
RUNS = {"WT": (f"{BASE}/boltz-runs/cofold-wt",  f"{BASE}/cofold_wt.json"),
        "MUT(germlined)": (f"{BASE}/boltz-runs/cofold-mut", f"{BASE}/cofold_mut.json")}

CDR_H3 = "QLYYFDYW"
LCDRS = {  # per-tag light-chain CDR substrings (mutant reverts these to germline)
    "WT": ["LASQTIGTWLT", "TATSLAD", "QQVYSIPWT"],
    "MUT(germlined)": ["RASQSISSYLN", "AASSLQS", "QQSYSTPWT"],
}

def locate(seq, sub, base):
    i = seq.find(sub)
    return list(range(base + i, base + i + len(sub))) if i >= 0 else []

def block(pae, a, b):
    a, b = np.array(a), np.array(b)
    return float((pae[np.ix_(a, b)].mean() + pae[np.ix_(b, a)].mean()) / 2)

def analyze(tag, run_dir, input_json):
    d = json.load(open(input_json))
    ents = {e["chain_ids"][0]: e["value"] for e in d["entities"]}
    H, L, V = ents["H"], ents["L"], ents["V"]
    offH, offL, offV = 0, len(H), len(H) + len(L)
    Vtok = list(range(offV, offV + len(V)))
    h3 = locate(H, CDR_H3, offH)
    lcdr = sum((locate(L, s, offL) for s in LCDRS[tag]), [])
    paratope = h3 + lcdr

    out_dir = os.path.join(run_dir, "outputs")
    arch = os.path.join(out_dir, "archive.tar.gz")
    if os.path.exists(arch):
        with tarfile.open(arch) as t: t.extractall(out_dir)
    metrics = json.load(open(glob.glob(f"{out_dir}/**/metrics.json", recursive=True)[0]))
    samples = metrics.get("all_sample_results", [])
    iptm = [s["metrics"].get("iptm") for s in samples]
    piptm = [s["metrics"].get("protein_iptm") for s in samples]
    plddt = [s["metrics"].get("complex_plddt") for s in samples]
    bconf = metrics.get("binding_metrics", {}).get("binding_confidence")

    pae_files = sorted(glob.glob(f"{out_dir}/**/sample_*_pae.npz", recursive=True))
    # de-dup (archive extracts to both ./ and ./files)
    seen, uniq = set(), []
    for f in pae_files:
        key = os.path.basename(f)
        if key not in seen: seen.add(key); uniq.append(f)
    rows = []
    for f in uniq:
        pae = np.load(f)["pae"]
        rows.append((block(pae, h3, Vtok), block(pae, lcdr, Vtok),
                     block(pae, paratope, Vtok), block(pae, h3 + lcdr + list(range(offH, offV)), Vtok)))
    rows = np.array(rows)  # [nsamp, 4] : h3, lcdr, paratope, wholeIF
    return {
        "tag": tag, "nsamp": len(uniq),
        "iptm_max": max([x for x in iptm if x is not None], default=None),
        "iptm_mean": float(np.mean([x for x in iptm if x is not None])) if any(iptm) else None,
        "protein_iptm_max": max([x for x in piptm if x is not None], default=None),
        "complex_plddt_max": max([x for x in plddt if x is not None], default=None),
        "binding_confidence": bconf,
        "paeH3_min": float(rows[:,0].min()), "paeH3_mean": float(rows[:,0].mean()),
        "paeLcdr_min": float(rows[:,1].min()),
        "paeParatope_min": float(rows[:,2].min()), "paeParatope_mean": float(rows[:,2].mean()),
        "paeWholeIF_min": float(rows[:,3].min()),
    }

res = {}
for tag, (rd, ij) in RUNS.items():
    if not os.path.exists(os.path.join(rd, "outputs", "archive.tar.gz")):
        print(f"[{tag}] NOT READY (no archive yet at {rd})"); continue
    res[tag] = analyze(tag, rd, ij)

if len(res) == 2:
    a, b = res["WT"], res["MUT(germlined)"]
    print("\n================  WT  vs  germlined-mutant  (SHR-1210 Fv × VEGFR2 D2-3)  ================")
    hdr = f"{'metric':<26}{'WT':>12}{'MUT':>12}   interpretation (thesis: WT binds, MUT collapses)"
    print(hdr); print("-"*len(hdr))
    def row(name, ka, better_low, note):
        va, vb = a[ka], b[ka]
        if va is None or vb is None: print(f"{name:<26}{'n/a':>12}{'n/a':>12}"); return
        good = (va < vb) if better_low else (va > vb)
        mark = "✓ supports" if good else "✗ against"
        print(f"{name:<26}{va:>12.3f}{vb:>12.3f}   {mark} — {note}")
    row("ipTM (interface, max)",       "iptm_max",          False, "higher = more confident interface")
    row("protein_ipTM (max)",          "protein_iptm_max",  False, "higher = better protein-protein interface")
    row("PAE CDR-H3↔VEGFR2 (min Å)",   "paeH3_min",         True,  "lower = tighter H3 interface")
    row("PAE L-CDRs↔VEGFR2 (min Å)",   "paeLcdr_min",       True,  "lower = tighter light-CDR interface (mutated region)")
    row("PAE paratope↔VEGFR2 (min)",   "paeParatope_min",   True,  "lower = tighter paratope interface")
    row("binding_confidence",          "binding_confidence",False, "Boltz binding head")
    print("\nfull:", json.dumps(res, indent=2))
else:
    print("\n(analysis will complete once both runs are downloaded)")
