import json, glob, os, tarfile, itertools, warnings, numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
CDR_H3 = "QLYYFDYW"
LCDRS = ["LASQTIGTWLT", "TATSLAD", "QQVYSIPWT"]  # WT SHR-1210 light-chain CDRs
CONTACT = 5.0

# (label, run_dir, input_json). FZD5 is the known off-target; rest are fold-matched decoys.
MEMBERS = [
    ("FZD5",  "cofold-fzd5",  "cofold_fzd5.json"),
    ("FZD1",  "cofold-fzd1",  "cofold_fzd1.json"),
    ("FZD2",  "cofold-fzd2",  "cofold_fzd2.json"),
    ("FZD3",  "cofold-fzd3",  "cofold_fzd3.json"),
    ("FZD4",  "cofold-fzd4",  "cofold_fzd4.json"),
    ("FZD6",  "cofold-fzd6",  "cofold_fzd6.json"),
    ("FZD7",  "cofold-fzd7",  "cofold_fzd7.json"),
    ("FZD8",  "cofold-fzd8",  "cofold_fzd8.json"),
    ("FZD9",  "cofold-fzd9",  "cofold_fzd9.json"),
    ("FZD10", "cofold-fzd10", "cofold_fzd10.json"),
    ("SMO",   "cofold-smo",   "cofold_smo.json"),
    ("SFRP1", "cofold-sfrp1", "cofold_sfrp1.json"),
]

def locate(seq, sub, base):
    i = seq.find(sub); return list(range(base+i, base+i+len(sub))) if i >= 0 else []

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

def analyze(run_dir, input_json):
    rd = f"{BASE}/boltz-runs/{run_dir}"; out = f"{rd}/outputs"
    if not os.path.exists(f"{out}/archive.tar.gz"): return None
    d = json.load(open(f"{BASE}/{input_json}"))
    ents = {e["chain_ids"][0]: e["value"] for e in d["entities"]}
    H, L, V = ents["H"], ents["L"], ents["V"]
    offL, offV = len(H), len(H)+len(L)
    Vtok = list(range(offV, offV+len(V)))
    paratope = locate(H, CDR_H3, 0) + sum((locate(L, s, offL) for s in LCDRS), [])
    with tarfile.open(f"{out}/archive.tar.gz") as t: t.extractall(out)
    paes = []; seen = set()
    for f in sorted(glob.glob(f"{out}/**/sample_*_pae.npz", recursive=True)):
        b = os.path.basename(f)
        if b in seen: continue
        seen.add(b); paes.append(block(np.load(f)["pae"], paratope, Vtok))
    cifs = sorted(glob.glob(f"{out}/prediction/sample_*_predicted_structure.cif"))
    eps = [epitope(c, len(V)) for c in cifs]
    jvals = [jac(a, b) for a, b in itertools.combinations(eps, 2)]
    return dict(alen=len(V), n_samples=len(paes),
                pae_mean=float(np.mean(paes)), pae_min=float(np.min(paes)),
                epi_size=float(np.mean([len(e) for e in eps])),
                epi_reprod=float(np.mean(jvals)) if jvals else float('nan'))

rows = []
for label, rd, ij in MEMBERS:
    r = analyze(rd, ij)
    if r is None:
        print(f"SKIPPED {label}: no archive"); continue
    r["label"] = label; rows.append(r)

# ---- Table sorted by PAE_mean (best interface first) ----
print(f"\n{'member':<8}{'AGlen':>6}{'PAE_mean':>10}{'PAE_min':>9}{'epi_sz':>8}{'epi_reprod':>12}")
print("-"*53)
for r in sorted(rows, key=lambda x: x["pae_mean"]):
    star = " *FZD5*" if r["label"]=="FZD5" else ""
    print(f"{r['label']:<8}{r['alen']:>6}{r['pae_mean']:>10.2f}{r['pae_min']:>9.2f}{r['epi_size']:>8.1f}{r['epi_reprod']:>12.3f}{star}")

# ---- within-family rank + AUROC of FZD5 vs 11 decoys ----
fzd5 = next(r for r in rows if r["label"]=="FZD5")
decoys = [r for r in rows if r["label"]!="FZD5"]
n = len(decoys)

# PAE_mean: lower is better
pae_beats = sum(1 for d in decoys if fzd5["pae_mean"] < d["pae_mean"])
pae_ties  = sum(1 for d in decoys if fzd5["pae_mean"] == d["pae_mean"])
auroc_pae = (pae_beats + 0.5*pae_ties)/n
rank_pae  = 1 + sum(1 for d in decoys if d["pae_mean"] < fzd5["pae_mean"])

# epi_reprod: higher is better
epi_beats = sum(1 for d in decoys if fzd5["epi_reprod"] > d["epi_reprod"])
epi_ties  = sum(1 for d in decoys if fzd5["epi_reprod"] == d["epi_reprod"])
auroc_epi = (epi_beats + 0.5*epi_ties)/n
rank_epi  = 1 + sum(1 for d in decoys if d["epi_reprod"] > fzd5["epi_reprod"])

print(f"\nWithin-family (FZD5 vs {n} decoys):")
print(f"  PAE_mean   : FZD5={fzd5['pae_mean']:.2f}  rank {rank_pae}/{len(rows)}  AUROC={auroc_pae:.3f}  (beats {pae_beats}/{n})")
print(f"  epi_reprod : FZD5={fzd5['epi_reprod']:.3f}  rank {rank_epi}/{len(rows)}  AUROC={auroc_epi:.3f}  (beats {epi_beats}/{n})")

json.dump(dict(rows=rows, auroc_pae=auroc_pae, auroc_epi=auroc_epi,
               rank_pae=rank_pae, rank_epi=rank_epi, n_decoys=n),
          open(f"{BASE}/frizzled_family_results.json","w"), indent=2)
