import glob, os, warnings, itertools, numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
RUNS = {"WT×VEGFR2": "cofold-wt", "MUT×VEGFR2": "cofold-mut",
        "WT×PD-1 (pos ctrl)": "cofold-pd1", "WT×lysozyme (neg ctrl)": "cofold-lyz"}
CONTACT = 5.0  # Å heavy-atom cutoff

def chain_by_len(model):
    """Map chains to roles by residue count: antibody = the two ~107-116 chains, antigen = the third."""
    chains = [(c.id, len([r for r in c if r.id[0] == ' '])) for c in model]
    chains.sort(key=lambda x: -x[1])
    # antigen is whichever isn't the two antibody chains (VH~116, VL~107). Antibody = the pair closest to 107-116.
    ab = [cid for cid, n in chains if 100 <= n <= 125]
    ag = [cid for cid, n in chains if cid not in ab]
    return ab, ag

def epitope(cif):
    m = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    ab, ag = chain_by_len(m)
    ab_atoms = [a for c in m if c.id in ab for a in c.get_atoms() if a.element != "H"]
    ns = NeighborSearch(ab_atoms)
    ep = set()
    for c in m:
        if c.id not in ag: continue
        for r in c:
            if r.id[0] != ' ': continue
            if any(ns.search(a.coord, CONTACT) for a in r if a.element != "H"):
                ep.add((c.id, r.id[1]))
    return ep

def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0

print(f"{'complex':<26}{'nsamp':>6}{'epitope_size':>14}{'epitope_reprod':>16}   (Jaccard: 1=same patch every time)")
print("-" * 90)
results = {}
for label, d in RUNS.items():
    cifs = sorted(glob.glob(f"{BASE}/boltz-runs/{d}/outputs/prediction/sample_*_predicted_structure.cif"))
    if not cifs:
        print(f"{label:<26}{'—':>6}   (not ready)"); continue
    eps = [epitope(c) for c in cifs]
    sizes = [len(e) for e in eps]
    jac = [jaccard(a, b) for a, b in itertools.combinations(eps, 2)]
    mj = float(np.mean(jac)) if jac else float("nan")
    results[label] = (len(cifs), float(np.mean(sizes)), mj)
    print(f"{label:<26}{len(cifs):>6}{np.mean(sizes):>14.1f}{mj:>16.3f}")
print("\ninterpretation: a true binder reproduces the SAME epitope across samples (high Jaccard);")
print("a non-binder docks scattered patches (low Jaccard). Compare WT×VEGFR2 to the pos/neg controls.")
