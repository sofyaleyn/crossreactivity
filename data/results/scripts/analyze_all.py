import json, glob, os, tarfile, itertools, warnings, numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
CDR_H3 = "QLYYFDYW"
LCDRS = {"WT": ["LASQTIGTWLT", "TATSLAD", "QQVYSIPWT"],
         "MUT": ["RASQSISSYLN", "AASSLQS", "QQSYSTPWT"]}
# label -> (run_dir, input_json, lcdr_key). Ordered ceiling->floor for readability.
RUNS = [
    ("WT × PD-1        [+ positive ctrl]", "cofold-pd1", "cofold_pd1.json", "WT"),
    ("WT × VEGFR2      [off-target 1]",    "cofold-wt",  "cofold_wt.json",  "WT"),
    ("WT × FZD5-CRD    [off-target 2]",    "cofold-fzd5","cofold_fzd5.json","WT"),
    ("WT × ULBP2       [off-target 3]",    "cofold-ulbp2","cofold_ulbp2.json","WT"),
    ("MUT × VEGFR2     [germlined]",       "cofold-mut", "cofold_mut.json", "MUT"),
    ("WT × VEGFR2+tmpl [3V2A template]",   "cofold-wt-tmpl",  "cofold_wt_tmpl.json",  "WT"),
    ("MUT × VEGFR2+tmpl[3V2A template]",   "cofold-mut-tmpl", "cofold_mut_tmpl.json", "MUT"),
    ("WT × lysozyme    [- negative ctrl]", "cofold-lyz", "cofold_lyz.json", "WT"),
]
CONTACT = 5.0

def locate(seq, sub, base):
    i = seq.find(sub);  return list(range(base+i, base+i+len(sub))) if i >= 0 else []

def block(pae, a, b):
    a, b = np.array(a), np.array(b)
    return float((pae[np.ix_(a, b)].mean() + pae[np.ix_(b, a)].mean())/2)

def chain_roles(model, antigen_len):
    counts = {c.id: len([r for r in c if r.id[0]==' ']) for c in model}
    ag = [cid for cid,n in counts.items() if n == antigen_len] or [min(counts, key=lambda c: abs(counts[c]-antigen_len))]
    ab = [cid for cid in counts if cid not in ag]
    return ab, ag

def epitope(cif, antigen_len):
    m = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    ab, ag = chain_roles(m, antigen_len)
    ab_atoms = [a for c in m if c.id in ab for a in c.get_atoms() if a.element != "H"]
    ns = NeighborSearch(ab_atoms); ep = set()
    for c in m:
        if c.id not in ag: continue
        for r in c:
            if r.id[0] != ' ': continue
            if any(ns.search(a.coord, CONTACT) for a in r if a.element != "H"): ep.add((c.id, r.id[1]))
    return ep

def jac(a, b): return len(a & b)/len(a | b) if (a | b) else 0.0

def analyze(run_dir, input_json, key):
    rd = f"{BASE}/boltz-runs/{run_dir}"; out = f"{rd}/outputs"
    if not os.path.exists(f"{out}/archive.tar.gz"): return None
    d = json.load(open(f"{BASE}/{input_json}"))
    ents = {e["chain_ids"][0]: e["value"] for e in d["entities"]}
    H, L, V = ents["H"], ents["L"], ents["V"]
    offL, offV = len(H), len(H)+len(L)
    Vtok = list(range(offV, offV+len(V)))
    paratope = locate(H, CDR_H3, 0) + sum((locate(L, s, offL) for s in LCDRS[key]), [])
    with tarfile.open(f"{out}/archive.tar.gz") as t: t.extractall(out)
    metrics = json.load(open(glob.glob(f"{out}/**/metrics.json", recursive=True)[0]))
    iptm = [s["metrics"].get("iptm") for s in metrics.get("all_sample_results", [])]
    iptm = [x for x in iptm if x is not None]
    paes = []
    seen = set()
    for f in sorted(glob.glob(f"{out}/**/sample_*_pae.npz", recursive=True)):
        b = os.path.basename(f)
        if b in seen: continue
        seen.add(b); paes.append(block(np.load(f)["pae"], paratope, Vtok))
    cifs = sorted(glob.glob(f"{out}/prediction/sample_*_predicted_structure.cif"))
    eps = [epitope(c, len(V)) for c in cifs]
    jvals = [jac(a, b) for a, b in itertools.combinations(eps, 2)]
    return dict(iptm_max=max(iptm), iptm_mean=float(np.mean(iptm)),
                pae_min=float(np.min(paes)), pae_mean=float(np.mean(paes)),
                epi_size=float(np.mean([len(e) for e in eps])),
                epi_reprod=float(np.mean(jvals)) if jvals else float('nan'))

print(f"{'complex':<38}{'ipTM':>7}{'ipTM_mn':>8}{'PAE_par↔ag':>12}{'PAE_mean':>10}{'epi_sz':>8}{'epi_reprod':>11}")
print("-"*94)
for label, rd, ij, key in RUNS:
    r = analyze(rd, ij, key)
    if r is None: print(f"{label:<38}{'… not ready':>20}"); continue
    print(f"{label:<38}{r['iptm_max']:>7.3f}{r['iptm_mean']:>8.3f}{r['pae_min']:>12.2f}{r['pae_mean']:>10.2f}{r['epi_size']:>8.1f}{r['epi_reprod']:>11.3f}")
print("\nReading: PD-1 = confident-binder ceiling; lysozyme = non-binder floor.")
print("Lower PAE_par↔ag = tighter paratope interface; higher epi_reprod = same epitope every sample.")
print("Thesis holds iff WT×VEGFR2 sits near the PD-1 ceiling AND clearly above lysozyme, and MUT slides toward the floor.")
