import json, glob, os, tarfile, itertools, csv, warnings, numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")
BASE="/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
RUNS=f"{BASE}/boltz-runs"; CONTACT=5.0

# run dir -> (antibody, antigen, role)
LBL={
 "cofold-pd1":("SHR-1210 WT","PD-1 (IgV)","cognate target / ceiling"),
 "cofold-wt":("SHR-1210 WT","VEGFR2 D2-3","off-target"),
 "cofold-fzd5":("SHR-1210 WT","FZD5 CRD","off-target"),
 "cofold-ulbp2":("SHR-1210 WT","ULBP2 ectodomain","off-target"),
 "cofold-mut":("SHR-1210 germlined","VEGFR2 D2-3","germlined variant"),
 "cofold-wt-tmpl":("SHR-1210 WT","VEGFR2 D2-3 +3V2A template","off-target, templated"),
 "cofold-mut-tmpl":("SHR-1210 germlined","VEGFR2 D2-3 +3V2A template","germlined, templated"),
 "cofold-lyz":("SHR-1210 WT","lysozyme (HEL)","non-binder / floor"),
 "cofold-pembro-pd1":("pembrolizumab","PD-1 (IgV)","Ab control - shared target"),
 "cofold-pembro-vegfr2":("pembrolizumab","VEGFR2 D2-3","Ab control"),
 "cofold-pembro-fzd5":("pembrolizumab","FZD5 CRD","Ab control"),
 "cofold-pembro-ulbp2":("pembrolizumab","ULBP2 ectodomain","Ab control"),
 "smoke-trpcage-aspirin":("(Trp-cage)","aspirin (ligand)","pipeline smoke test"),
}
for f in ["fzd1","fzd2","fzd3","fzd4","fzd6","fzd7","fzd8","fzd9","fzd10","smo","sfrp1"]:
    LBL[f"cofold-{f}"]=("SHR-1210 WT",f"{f.upper()} CRD","within-fold decoy (FZD5 family)")

def block(pae,a,b):
    a,b=np.array(a),np.array(b); return float((pae[np.ix_(a,b)].mean()+pae[np.ix_(b,a)].mean())/2)
def roles(model,aglen):
    cnt={c.id:len([r for r in c if r.id[0]==' ']) for c in model}
    ag=[c for c,n in cnt.items() if n==aglen] or [min(cnt,key=lambda c:abs(cnt[c]-aglen))]
    return [c for c in cnt if c not in ag],ag
def epitope(cif,aglen):
    m=MMCIFParser(QUIET=True).get_structure("x",cif)[0]; ab,ag=roles(m,aglen)
    at=[a for c in m if c.id in ab for a in c.get_atoms() if a.element!="H"]; ns=NeighborSearch(at); ep=set()
    for c in m:
        if c.id not in ag: continue
        for r in c:
            if r.id[0]==' ' and any(ns.search(a.coord,CONTACT) for a in r if a.element!="H"): ep.add((c.id,r.id[1]))
    return ep
def jac(a,b): return len(a&b)/len(a|b) if (a|b) else 0.0

def find_input(tag):
    for cand in [f"{BASE}/cofold_{tag.replace('cofold-','').replace('-','_')}.json",
                 f"{BASE}/smoke_input.json" if tag.startswith("smoke") else None]:
        if cand and os.path.exists(cand): return cand
    return None

rows=[]
for rd in sorted(glob.glob(f"{RUNS}/cofold-*")+glob.glob(f"{RUNS}/smoke-*")):
    tag=os.path.basename(rd); out=f"{rd}/outputs"
    arch=f"{out}/archive.tar.gz"
    if not os.path.exists(arch): continue
    with tarfile.open(arch) as t: t.extractall(out)
    mj=glob.glob(f"{out}/**/metrics.json",recursive=True)
    if not mj: continue
    met=json.load(open(mj[0])); samples=met.get("all_sample_results",[])
    def col(k):
        v=[s["metrics"].get(k) for s in samples]; return [x for x in v if x is not None]
    iptm=col("iptm"); ptm=col("ptm"); piptm=col("protein_iptm"); plddt=col("complex_plddt"); sc=col("structure_confidence")
    bm=met.get("binding_metrics",{})
    bconf=bm.get("binding_confidence")
    # prediction id
    pid=""
    try: pid=json.load(open(f"{rd}/run.json")).get("id","")
    except Exception: pass
    ab,ag,role=LBL.get(tag,("?","?","?"))
    # PAE_IF + epitope reproducibility (needs input json for antigen length)
    pae_if=epi=""
    ij=find_input(tag)
    if ij:
        try:
            d=json.load(open(ij)); ents={e["chain_ids"][0]:e.get("value","") for e in d["entities"]}
            if "V" in ents and "H" in ents and "L" in ents:
                offV=len(ents["H"])+len(ents["L"]); Vt=list(range(offV,offV+len(ents["V"]))); abt=list(range(0,offV))
                seen=set(); pv=[]
                for f in sorted(glob.glob(f"{out}/**/sample_*_pae.npz",recursive=True)):
                    if os.path.basename(f) in seen: continue
                    seen.add(os.path.basename(f)); pv.append(block(np.load(f)["pae"],abt,Vt))
                pae_if=round(float(np.mean(pv)),2) if pv else ""
                cifs=sorted(glob.glob(f"{out}/prediction/sample_*_predicted_structure.cif"))
                eps=[epitope(c,len(ents["V"])) for c in cifs]
                jv=[jac(a,b) for a,b in itertools.combinations(eps,2)]
                epi=round(float(np.mean(jv)),3) if jv else ""
        except Exception as e: pass
    rows.append(dict(run=tag,antibody=ab,antigen=ag,role=role,num_samples=len(samples),
        iptm_max=round(max(iptm),3) if iptm else "", iptm_mean=round(float(np.mean(iptm)),3) if iptm else "",
        protein_iptm_max=round(max(piptm),3) if piptm else "", ptm_max=round(max(ptm),3) if ptm else "",
        complex_plddt_max=round(max(plddt),3) if plddt else "", structure_confidence_max=round(max(sc),3) if sc else "",
        binding_confidence=(round(bconf,4) if isinstance(bconf,(int,float)) else ""),
        PAE_IF_mean=pae_if, epitope_reprod=epi, prediction_id=pid))

cols=["run","antibody","antigen","role","num_samples","iptm_max","iptm_mean","protein_iptm_max",
      "ptm_max","complex_plddt_max","structure_confidence_max","binding_confidence","PAE_IF_mean","epitope_reprod","prediction_id"]
outcsv=f"{BASE}/cofold_metrics.csv"
with open(outcsv,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for r in rows: w.writerow(r)
print(f"wrote {outcsv} ({len(rows)} runs)\n")
# markdown preview
print(f"{'run':22s}{'antibody':20s}{'antigen':30s}{'ipTM':>6}{'ipTMmn':>7}{'pIF':>7}{'epi':>7}")
for r in rows:
    print(f"{r['run']:22s}{r['antibody']:20s}{r['antigen']:30s}{str(r['iptm_max']):>6}{str(r['iptm_mean']):>7}{str(r['PAE_IF_mean']):>7}{str(r['epitope_reprod']):>7}")
