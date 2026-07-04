import csv, sys, json, time
import numpy as np
import torch, esm

SCRATCH="/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
CSV="/Users/cheparukhin/crossreactivity/data/curated/self_proteins.csv"
POS={"P35968","Q13467","Q9BZM5"}

VH="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL="DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
H3="QLYYFDYW"; L1="LASQTIGTWLT"; L2="TATSLAD"; L3="QQVYSIPWT"

WIN=1000; STEP=500
VALID=set("ACDEFGHIKLMNPQRSTVWY")

torch.set_num_threads(8)
model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()
LAYER=12

def per_residue(seq):
    """return (L,480) per-residue embeddings for a single (<=1022) sequence."""
    _,_,toks = bc([("x",seq)])
    with torch.no_grad():
        out=model(toks,repr_layers=[LAYER],return_contacts=False)
    rep=out["representations"][LAYER][0]  # (L+2,480) incl BOS/EOS
    return rep[1:len(seq)+1].numpy()

def clean(s):
    return "".join(c for c in (s or "").upper() if c in VALID)

# ---- paratope + whole-Fv query vectors ----
vh_rep=per_residue(VH); vl_rep=per_residue(VL)
def cdr_vec(rep,seq,sub):
    i=seq.find(sub); assert i>=0, sub
    return rep[i:i+len(sub)]
loops=[cdr_vec(vh_rep,VH,H3),cdr_vec(vl_rep,VL,L1),cdr_vec(vl_rep,VL,L2),cdr_vec(vl_rep,VL,L3)]
paratope=np.concatenate(loops,axis=0).mean(axis=0)
wholefv=np.concatenate([vh_rep,vl_rep],axis=0).mean(axis=0)
def norm(v): return v/ (np.linalg.norm(v)+1e-9)
paratope_n=norm(paratope); wholefv_n=norm(wholefv)
print("query vectors built", file=sys.stderr)

# ---- load antigens ----
rows=[]
with open(CSV) as f:
    for row in csv.DictReader(f):
        rows.append((row["uniprot_id"], clean(row["sequence"])))
print("antigens",len(rows), file=sys.stderr)

t0=time.time()
results=[]
for n,(uid,seq) in enumerate(rows):
    L=len(seq)
    if L==0:
        results.append((uid,0.0,0.0)); continue
    # windows
    starts=list(range(0,max(1,L-WIN+1),STEP))
    if not starts or starts[-1]<L-WIN:
        pass
    if L<=WIN: starts=[0]
    else:
        starts=list(range(0,L-WIN+1,STEP))
        if starts[-1]!=L-WIN: starts.append(L-WIN)
    best_p=-1.0; best_w=-1.0
    for s in starts:
        sub=seq[s:s+WIN]
        # ESM max token length 1022
        if len(sub)>1022: sub=sub[:1022]
        rep=per_residue(sub)          # (w,480)
        wemb=norm(rep.mean(axis=0))
        cp=float(np.dot(paratope_n,wemb)); cw=float(np.dot(wholefv_n,wemb))
        if cp>best_p: best_p=cp
        if cw>best_w: best_w=cw
    results.append((uid,best_p,best_w))
    if (n+1)%200==0:
        el=time.time()-t0
        print(f"{n+1}/{len(rows)}  {el:.0f}s  ~{el/(n+1)*len(rows):.0f}s total", file=sys.stderr, flush=True)

with open(f"{SCRATCH}/embed_scores.json","w") as f:
    json.dump({uid:[p,w] for uid,p,w in results}, f)

def report(name, idx):
    scored=[(uid,vals[idx]) for uid,*vals in [(r[0],r[1],r[2]) for r in results]]
    s=sorted(results,key=lambda r:r[1+idx],reverse=True)
    N=len(s); rankmap={r[0]:i+1 for i,r in enumerate(s)}
    print(f"\n=== {name} ===")
    pcts=[]
    for uid in POS:
        rk=rankmap[uid]; pct=rk/N; pcts.append(pct)
        val=next(r[1+idx] for r in results if r[0]==uid)
        print(f"  {uid} rank {rk}/{N} pct {pct:.3f} cos={val:.4f}")
    worst=max(pcts)
    tk=lambda p:sum(1 for x in pcts if x<=p)
    print(f"  top5%:{tk(.05)}/3 top10%:{tk(.10)}/3 top20%:{tk(.20)}/3  RETAIN-ALL-3: {worst:.3f}")
    return worst

r1=report("embed_paratope",0)
r2=report("embed_wholeFv",1)
print(f"\nembed_paratope retain-all-3: {r1:.3f}")
print(f"embed_wholeFv  retain-all-3: {r2:.3f}")
print("DONE embed", file=sys.stderr)
