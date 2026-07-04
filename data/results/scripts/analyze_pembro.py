import json, glob, os, tarfile, itertools, warnings, numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")
BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
CONTACT = 5.0

# antigen -> (SHR1210 run_dir & input, pembro run_dir & input)
ANTIGENS = {
 "PD-1 (target)":  (("cofold-pd1","cofold_pd1.json"),   ("cofold-pembro-pd1","cofold_pembro_pd1.json")),
 "VEGFR2 (offtgt)":(("cofold-wt","cofold_wt.json"),     ("cofold-pembro-vegfr2","cofold_pembro_vegfr2.json")),
 "FZD5 (offtgt)":  (("cofold-fzd5","cofold_fzd5.json"), ("cofold-pembro-fzd5","cofold_pembro_fzd5.json")),
 "ULBP2 (offtgt)": (("cofold-ulbp2","cofold_ulbp2.json"),("cofold-pembro-ulbp2","cofold_pembro_ulbp2.json")),
}

def block(pae,a,b):
    a,b=np.array(a),np.array(b); return float((pae[np.ix_(a,b)].mean()+pae[np.ix_(b,a)].mean())/2)
def roles(model, aglen):
    cnt={c.id:len([r for r in c if r.id[0]==' ']) for c in model}
    ag=[c for c,n in cnt.items() if n==aglen] or [min(cnt,key=lambda c:abs(cnt[c]-aglen))]
    return [c for c in cnt if c not in ag], ag
def epitope(cif,aglen):
    m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]; ab,ag=roles(m,aglen)
    at=[a for c in m if c.id in ab for a in c.get_atoms() if a.element!="H"]; ns=NeighborSearch(at); ep=set()
    for c in m:
        if c.id not in ag: continue
        for r in c:
            if r.id[0]!=' ': continue
            if any(ns.search(a.coord,CONTACT) for a in r if a.element!="H"): ep.add((c.id,r.id[1]))
    return ep
def jac(a,b): return len(a&b)/len(a|b) if (a|b) else 0.0

def analyze(run_dir, input_json):
    rd=f"{BASE}/boltz-runs/{run_dir}"; out=f"{rd}/outputs"
    d=json.load(open(f"{BASE}/{input_json}")); ents={e["chain_ids"][0]:e["value"] for e in d["entities"]}
    H,L,V=ents["H"],ents["L"],ents["V"]; offV=len(H)+len(L)
    ab_tok=list(range(0,offV)); ag_tok=list(range(offV,offV+len(V)))
    with tarfile.open(f"{out}/archive.tar.gz") as t: t.extractall(out)
    met=json.load(open(glob.glob(f"{out}/**/metrics.json",recursive=True)[0]))
    iptm=[s["metrics"].get("iptm") for s in met.get("all_sample_results",[])]; iptm=[x for x in iptm if x is not None]
    seen=set(); paes=[]
    for f in sorted(glob.glob(f"{out}/**/sample_*_pae.npz",recursive=True)):
        if os.path.basename(f) in seen: continue
        seen.add(os.path.basename(f)); paes.append(block(np.load(f)["pae"],ab_tok,ag_tok))
    cifs=sorted(glob.glob(f"{out}/prediction/sample_*_predicted_structure.cif"))
    eps=[epitope(c,len(V)) for c in cifs]; jv=[jac(a,b) for a,b in itertools.combinations(eps,2)]
    return max(iptm), float(np.mean(paes)), float(np.mean(jv)) if jv else float('nan')

print(f"{'antigen':<18}{'SHR-1210 (has off-tgts)':>26}{'pembrolizumab (no off-tgts)':>30}")
print(f"{'':<18}{'ipTM  PAE_IF  epi_rep':>26}{'ipTM  PAE_IF  epi_rep':>30}")
print("-"*76)
for name,(shr,pem) in ANTIGENS.items():
    s=analyze(*shr); p=analyze(*pem)
    print(f"{name:<18}{s[0]:>7.3f}{s[1]:>8.2f}{s[2]:>8.3f}   {p[0]:>10.3f}{p[1]:>8.2f}{p[2]:>8.3f}")
print("\nlower PAE_IF (whole antibody↔antigen interface) = tighter; higher epi_rep = reproducible epitope.")
print("SPECIFICITY holds iff both bind PD-1, but only SHR-1210 (not pembrolizumab) cofolds the off-targets FZD5/ULBP2.")
