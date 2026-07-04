import json, warnings, numpy as np, torch, esm
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
BASE="/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
VH="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL="DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
SUBS=["QLYYFDYW","LASQTIGTWLT","TATSLAD","QQVYSIPWT"]
MAXASA={'A':129.,'R':274.,'N':195.,'D':193.,'C':167.,'Q':225.,'E':223.,'G':104.,'H':224.,'I':197.,'L':201.,'K':236.,'M':224.,'F':240.,'P':159.,'S':155.,'T':172.,'W':285.,'Y':263.,'V':174.}
T3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
model,alph=esm.pretrained.esm2_t12_35M_UR50D(); model.eval(); bc=alph.get_batch_converter(); L=12
def perres(seq):
    o=[]
    for i in range(0,len(seq),1000):
        sub=seq[i:i+1000]; _,_,tok=bc([("x",sub)])
        with torch.no_grad(): rep=model(tok,repr_layers=[L])["representations"][L][0]
        o.append(rep[1:len(sub)+1].numpy())
    return np.concatenate(o,0)
fv=VH+VL; m=perres(fv); idx=[]
for s in SUBS:
    p=fv.find(s); idx+=list(range(p,p+len(s)))
q=m[idx].mean(0); q/=np.linalg.norm(q)+1e-9
ag=json.load(open(f"{BASE}/benchmark_antigens.json"))
labels={a:(v['family'],v['name'],v['is_pos']) for a,v in ag.items()}
P=PDBParser(QUIET=True); SR=ShrakeRupley()
cache={}
for acc,v in ag.items():
    seq=v['seq']; s=P.get_structure('x',f"{BASE}/af_pdbs/AF-{acc}.pdb"); SR.compute(s,level="R")
    rsa=np.full(len(seq),np.nan)
    for res in s[0]['A']:
        if not is_aa(res,standard=True): continue
        aa=T3.get(res.resname); rn=res.id[1]
        if aa and 0<=rn-1<len(seq): rsa[rn-1]=res.sasa/MAXASA[aa]
    cache[acc]=(perres(seq),rsa)
def auroc(sc):
    fams={}
    for a,l in labels.items(): fams.setdefault(l[0],[]).append(a)
    out={}
    for f,mem in fams.items():
        pos=[a for a in mem if labels[a][2]][0]; dec=[a for a in mem if not labels[a][2]]
        out[f]=sum(sc[pos]>sc[a] for a in dec)/len(dec)
    return out
fo=["RTK/Ig-receptor","Frizzled/CRD","MHC-I-like"]
print(f"{'thresh':>8s} "+" ".join(f"{x[:9]:>10s}" for x in fo)+f"{'MEAN':>8s}")
for th in [0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50]:
    sc={}
    for acc,(mat,rsa) in cache.items():
        mask=np.nan_to_num(rsa,nan=-1)>th
        if mask.sum()<5: mask=np.nan_to_num(rsa,nan=-1)>0
        vv=mat[:len(rsa)][mask[:len(rsa)]].mean(0); vv/=np.linalg.norm(vv)+1e-9
        sc[acc]=float(np.dot(q,vv))
    au=auroc(sc); mean=np.mean([au[f] for f in fo])
    print(f"{th:8.2f} "+" ".join(f"{au[f]:10.3f}" for f in fo)+f"{mean:8.3f}")
