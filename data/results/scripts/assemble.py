import json, csv
SCRATCH="/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
POS={"P35968","Q13467","Q9BZM5"}
GENE={"P35968":"VEGFR2/KDR","Q13467":"FZD5","Q9BZM5":"ULBP2"}

fast=json.load(open(f"{SCRATCH}/fast_scores.json"))
emb=json.load(open(f"{SCRATCH}/embed_scores.json"))  # uid -> [paratope,wholefv]

for d in fast:
    p,w=emb.get(d["uniprot_id"],[0.0,0.0])
    d["embed_paratope"]=p; d["embed_wholeFv"]=w

# scores CSV
cols=["uniprot_id","gene","embed_paratope","embed_wholeFv","naive_cdr_identity",
"biophys_hydrophobicity","biophys_charge","biophys_length","annotation_receptor_prior","is_positive"]
with open(f"{SCRATCH}/tier1_filter_scores.csv","w",newline="") as f:
    wr=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore"); wr.writeheader()
    for d in fast: wr.writerow(d)

N=len(fast)
def analyze(key,higher=True):
    s=sorted(fast,key=lambda d:d[key],reverse=higher)
    rk={d["uniprot_id"]:i+1 for i,d in enumerate(s)}
    info={}
    pcts=[]
    for uid in POS:
        r=rk[uid]; pct=r/N; pcts.append(pct)
        info[uid]=(r,pct,next(d[key] for d in fast if d["uniprot_id"]==uid))
    worst=max(pcts)
    tk=lambda p:sum(1 for x in pcts if x<=p)
    return info,worst,(tk(.05),tk(.10),tk(.20))

filters=[("embed_paratope","embed_paratope",True),
("embed_wholeFv","embed_wholeFv",True),
("naive_cdr_identity","naive_cdr_identity",True),
("biophys_hydrophobicity(high)","biophys_hydrophobicity",True),
("biophys_charge(high)","biophys_charge",True),
("biophys_length(long)","biophys_length",True),
("annotation_receptor_prior","annotation_receptor_prior",True)]

print(f"N={N}\n")
print(f"{'FILTER':30s} {'VEGFR2':>14s} {'FZD5':>14s} {'ULBP2':>14s}  top5/10/20  RETAIN")
rows_summary=[]
for label,key,hi in filters:
    info,worst,tks=analyze(key,hi)
    def cell(uid):
        r,pct,_=info[uid]; return f"{r}({pct*100:.0f}%)"
    print(f"{label:30s} {cell('P35968'):>14s} {cell('Q13467'):>14s} {cell('Q9BZM5'):>14s}   {tks[0]}/{tks[1]}/{tks[2]}    {worst:.3f}")
    rows_summary.append((label,worst,tks))

print("\nRanked by RETAIN-ALL-3 (lower=better):")
for label,worst,tks in sorted(rows_summary,key=lambda x:x[1]):
    print(f"  {worst:.3f}  top10={tks[1]}/3  {label}")

base=[w for l,w,t in rows_summary if l=="annotation_receptor_prior"][0]
print(f"\nannotation baseline retain-all-3 = {base:.3f}")
best=min(rows_summary,key=lambda x:x[1])
print(f"best filter = {best[0]}  retain={best[1]:.3f}")
