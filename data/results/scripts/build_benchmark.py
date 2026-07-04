import json, urllib.request, sys, numpy as np, torch, esm

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
PARATOPE_SUBS = ["QLYYFDYW", "LASQTIGTWLT", "TATSLAD", "QQVYSIPWT"]  # H3 + L1/L2/L3 (WT)

# family -> {positive_acc: name, decoys: {acc:name}}
FAMILIES = {
 "RTK/Ig-receptor": {"pos": ("P35968","VEGFR2"), "decoys": {
    "P17948":"VEGFR1","P35916":"VEGFR3","P00533":"EGFR","P04626":"HER2","P21860":"HER3",
    "P06213":"INSR","P08069":"IGF1R","P16234":"PDGFRA","P09619":"PDGFRB","P10721":"KIT",
    "P36888":"FLT3","P08581":"MET","P30530":"AXL","P11362":"FGFR1"}},
 "Frizzled/CRD": {"pos": ("Q13467","FZD5"), "decoys": {
    "Q9UP38":"FZD1","Q14332":"FZD2","Q9NPG1":"FZD3","Q9ULV1":"FZD4","O60353":"FZD6",
    "O75084":"FZD7","Q9H461":"FZD8","O00144":"FZD9","Q9ULW2":"FZD10","Q99835":"SMO","Q8N474":"SFRP1"}},
 "MHC-I-like": {"pos": ("Q9BZM5","ULBP2"), "decoys": {
    "Q9BZM6":"ULBP1","Q9BZM4":"ULBP3","Q8TD07":"ULBP4","Q6H3X3":"ULBP5","Q29983":"MICA",
    "Q29980":"MICB","P04439":"HLA-A","P13747":"HLA-E","P15813":"CD1D","Q30201":"HFE",
    "P55899":"FcRn","Q95460":"MR1"}},
}

def fetch(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30) as r:
        return "".join(l for l in r.read().decode().splitlines() if not l.startswith(">"))

print("Loading ESM-2 t12_35M ...", flush=True)
model, alph = esm.pretrained.esm2_t12_35M_UR50D(); model.eval()
bc = alph.get_batch_converter(); LAYER = 12

def per_res(seq):
    chunks = [(0,seq)] if len(seq)<=1000 else [(i,seq[i:i+1000]) for i in range(0,len(seq),500)]
    outs=[]
    for off,sub in chunks:
        _,_,tok = bc([("x",sub)])
        with torch.no_grad():
            rep = model(tok, repr_layers=[LAYER])["representations"][LAYER][0]
        outs.append(rep[1:len(sub)+1].numpy())
    return outs

def paratope_vec():
    fv=VH+VL; mat=per_res(fv)[0]; idx=[]
    for s in PARATOPE_SUBS:
        p=fv.find(s); idx+=list(range(p,p+len(s))) if p>=0 else []
    v=mat[idx].mean(0); return v/(np.linalg.norm(v)+1e-9)

def antigen_windows(seq):
    vs=[]
    for m in per_res(seq):
        v=m.mean(0); vs.append(v/(np.linalg.norm(v)+1e-9))
    return vs

q = paratope_vec()
antigens, scores, labels = {}, {}, {}   # acc -> seq / score / (family,name,is_pos)
for fam,info in FAMILIES.items():
    members = [(info["pos"][0], info["pos"][1], True)] + [(a,n,False) for a,n in info["decoys"].items()]
    for acc,name,is_pos in members:
        try:
            seq=fetch(acc); antigens[acc]=seq
            scores[acc]=max(float(np.dot(q,w)) for w in antigen_windows(seq))
            labels[acc]=(fam,name,is_pos)
            print(f"  {fam:16s} {name:8s} {acc}  len {len(seq):5d}  score {scores[acc]:.3f}", flush=True)
        except Exception as ex:
            print(f"  FAIL {acc} {name}: {ex}", file=sys.stderr)

# save dataset + scores for reuse by other methods (cofold / surface)
json.dump({acc:{"seq":antigens[acc], **dict(zip(["family","name","is_pos"],labels[acc]))} for acc in antigens},
          open(f"{BASE}/benchmark_antigens.json","w"))
json.dump({"method":"esm_paratope_cosine","scores":scores,
           "labels":{a:list(labels[a]) for a in labels}}, open(f"{BASE}/benchmark_embed_scores.json","w"), indent=2)

print("\n================  EMBEDDING RANKING (same-model ESM, best-case)  ================")
allscore = sorted(scores.items(), key=lambda kv:-kv[1])
rank_of = {a:i+1 for i,(a,_) in enumerate(allscore)}
N=len(allscore)
def fam_members(fam): return [a for a in labels if labels[a][0]==fam]
withins=[]
for fam,info in FAMILIES.items():
    pos=info["pos"][0]; mem=fam_members(fam)
    fam_sorted=sorted(mem,key=lambda a:-scores[a])
    wrank=fam_sorted.index(pos)+1
    decoys=[a for a in mem if not labels[a][2]]
    beats=sum(scores[pos]>scores[a] for a in decoys); auroc=beats/len(decoys)
    withins.append(auroc)
    print(f"{fam:16s} positive {info['pos'][1]:8s}: GLOBAL rank {rank_of[pos]}/{N} (top {100*rank_of[pos]/N:.0f}%) | "
          f"WITHIN-FAMILY rank {wrank}/{len(mem)} | beats {beats}/{len(decoys)} fold-matched decoys (AUROC {auroc:.2f})")
print(f"\nMean within-family AUROC = {np.mean(withins):.2f}   (0.5 = no binding-specific signal / pure fold confound; 1.0 = perfect)")
print("Interpretation: high GLOBAL rank + ~0.5 within-family AUROC == the 'signal' is fold similarity, not off-target binding.")
