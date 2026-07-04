import csv, math, sys

CSV="/Users/cheparukhin/crossreactivity/data/curated/self_proteins.csv"
POS={"P35968","Q13467","Q9BZM5"}

VH="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL="DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
H3="QLYYFDYW"; L1="LASQTIGTWLT"; L2="TATSLAD"; L3="QQVYSIPWT"

# Kyte-Doolittle
KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
CHG={'D':-1,'E':-1,'K':1,'R':1,'H':0.1}

def best_identity(short, seq):
    """slide short string along seq, best fraction of matching positions."""
    ls=len(short); best=0.0
    if len(seq)<ls:
        # compare against whole seq window
        n=min(len(seq),ls);
        m=sum(1 for i in range(n) if short[i]==seq[i]);
        return m/ls
    for i in range(0,len(seq)-ls+1):
        m=0
        for j in range(ls):
            if short[j]==seq[i+j]: m+=1
        f=m/ls
        if f>best: best=f
    return best

RECEPTOR_KW=["receptor","kinase","frizzled","integrin","cadherin","mhc","hla",
"ephrin","eph ","cd","ligand","channel","transporter","gpcr","semaphorin",
"notch","toll-like","interleukin","growth factor","tyrosine-protein kinase",
"g-protein","adhesion","claudin","tetraspanin","immunoglobulin","natural killer",
"nkg2","plexin","neuropilin","glypican","syndecan","selectin"]

def receptor_score(name, gene):
    t=(name or "").lower()+" "+(gene or "").lower()
    s=0
    for kw in RECEPTOR_KW:
        if kw in t: s+=1
    return s

rows=[]
with open(CSV) as f:
    r=csv.DictReader(f)
    for row in r:
        rows.append(row)
print("rows",len(rows), file=sys.stderr)

out=[]
for row in rows:
    seq=(row["sequence"] or "").upper()
    seq="".join(c for c in seq if c in KD)
    uid=row["uniprot_id"]
    name=row.get("name",""); gene=row.get("gene_symbol","")
    # filter 3: naive_cdr_identity — best of H3 and L3, and combined concat
    id_h3=best_identity(H3,seq)
    id_l3=best_identity(L3,seq)
    id_concat=best_identity(H3+L3,seq)
    cdr_score=max(id_h3,id_l3,id_concat)
    # filter 4: biophys
    L=len(seq)
    hyd=sum(KD[c] for c in seq)/L if L else 0
    charge=sum(CHG.get(c,0) for c in seq)
    # filter 5: receptor prior
    rec=receptor_score(name,gene)
    out.append(dict(uniprot_id=uid,gene=gene,name=name,length=L,
        naive_cdr_identity=cdr_score, id_h3=id_h3, id_l3=id_l3,
        biophys_hydrophobicity=hyd, biophys_charge=charge, biophys_length=L,
        annotation_receptor_prior=rec,
        is_positive=1 if uid in POS else 0))

def report(name, key, higher=True):
    s=sorted(out,key=lambda d:d[key],reverse=higher)
    N=len(s)
    rankmap={d["uniprot_id"]:(i+1) for i,d in enumerate(s)}
    print(f"\n=== {name} (higher={'more suspect' if higher else 'reversed'}) ===")
    pcts=[]
    for uid in POS:
        rk=rankmap[uid]; pct=rk/N
        pcts.append(pct)
        d=next(x for x in out if x["uniprot_id"]==uid)
        print(f"  {uid} {d['gene']:8s} rank {rk}/{N}  pct {pct:.3f}  score={d[key]:.4f}")
    worst=max(pcts)
    def topk(p): return sum(1 for x in pcts if x<=p)
    print(f"  in top5%: {topk(0.05)}/3  top10%: {topk(0.10)}/3  top20%: {topk(0.20)}/3")
    print(f"  RETAIN-ALL-3 fraction: {worst:.3f}")
    return worst

results={}
results["naive_cdr_identity"]=report("naive_cdr_identity","naive_cdr_identity")
results["biophys_hydrophobicity_high"]=report("naive_biophys hydrophobicity (high)","biophys_hydrophobicity",True)
results["biophys_hydrophobicity_low"]=report("naive_biophys hydrophobicity (low)","biophys_hydrophobicity",False)
results["biophys_charge_high"]=report("naive_biophys charge (high/positive)","biophys_charge",True)
results["biophys_charge_low"]=report("naive_biophys charge (low/negative)","biophys_charge",False)
results["biophys_length_high"]=report("naive_biophys length (long)","biophys_length",True)
results["biophys_length_low"]=report("naive_biophys length (short)","biophys_length",False)
results["annotation_receptor_prior"]=report("annotation_receptor_prior","annotation_receptor_prior")

# save partial CSV (fast filters)
import json
with open("/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad/fast_scores.json","w") as f:
    json.dump(out,f)
print("\nSUMMARY retain-all-3 fractions (lower=better):")
for k,v in sorted(results.items(),key=lambda x:x[1]):
    print(f"  {v:.3f}  {k}")
