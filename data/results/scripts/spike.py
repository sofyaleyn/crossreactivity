"""
Rung-1 discrimination spike for CrossFlag viability.
Best-case test: embed BOTH antibody Fv and antigens with the SAME ESM-2 model
(so cosine is well-defined), then check whether the 3 known SHR-1210 off-targets
rank above unrelated decoys -- especially 5 HARD membrane-receptor decoys.
If even this charitable single-space version fails, the embedding-similarity
premise is empirically unsupported.
"""
import sys, urllib.request, numpy as np, torch, esm

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
CDR_H3 = "QLYYFDYW"
CDR_L3 = "QQVYSIPWT"

POSITIVES = {"P35968": "VEGFR2", "Q13467": "FZD5", "Q9BZM5": "ULBP2"}
SOFT_DECOYS = {"P02768":"Albumin","P02787":"Transferrin","P00738":"Haptoglobin",
               "P02766":"Transthyretin","P01009":"A1AT","P69905":"HBA","P68871":"HBB","P02749":"ApoH"}
HARD_DECOYS = {"P00533":"EGFR","P06213":"INSR","P08069":"IGF1R","P04626":"HER2","P16234":"PDGFRA"}

def fetch(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        lines = r.read().decode().splitlines()
    return "".join(l for l in lines if not l.startswith(">"))

print("Loading ESM-2 t12_35M ...", flush=True)
model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()
LAYER = 12

def per_residue(seq):
    """Return per-residue embedding matrix [L, D] (windowed for long seqs)."""
    if len(seq) <= 1000:
        chunks = [(0, seq)]
    else:
        chunks = [(i, seq[i:i+1000]) for i in range(0, len(seq), 500)]
    outs = []
    for off, sub in chunks:
        _, _, toks = bc([("x", sub)])
        with torch.no_grad():
            rep = model(toks, repr_layers=[LAYER])["representations"][LAYER][0]
        outs.append((off, rep[1:len(sub)+1].numpy()))  # strip BOS/EOS
    return outs  # list of (offset, [len,D])

def mean_vec(seq):
    outs = per_residue(seq)
    allres = np.concatenate([m for _, m in outs], axis=0)
    v = allres.mean(0); return v/ (np.linalg.norm(v)+1e-9)

def window_vecs(seq):
    """One mean-pooled vector per window (for max-cosine antigen scoring)."""
    vs = []
    for _, m in per_residue(seq):
        v = m.mean(0); vs.append(v/(np.linalg.norm(v)+1e-9))
    return vs

def cdr_paratope_vec():
    fv = VH + VL
    outs = per_residue(fv)  # <=1000, single chunk
    mat = outs[0][1]
    idx = []
    for cdr, chain_seq, base in [(CDR_H3, VH, 0), (CDR_L3, VL, len(VH))]:
        p = fv.find(cdr)
        if p < 0:
            print(f"WARN: {cdr} not found", file=sys.stderr); continue
        idx += list(range(p, p+len(cdr)))
    v = mat[idx].mean(0); return v/(np.linalg.norm(v)+1e-9)

# antibody queries
q_whole = mean_vec(VH+VL)
q_para = cdr_paratope_vec()

# antigens
antigens = {**POSITIVES, **SOFT_DECOYS, **HARD_DECOYS}
labels = {**{a:("POS",n) for a,n in POSITIVES.items()},
          **{a:("soft",n) for a,n in SOFT_DECOYS.items()},
          **{a:("HARD",n) for a,n in HARD_DECOYS.items()}}

print(f"Fetching + embedding {len(antigens)} antigens ...", flush=True)
awins = {}
for acc in antigens:
    try:
        seq = fetch(acc)
        awins[acc] = window_vecs(seq)
        print(f"  {labels[acc][1]:12s} {acc}  len={len(seq)}", flush=True)
    except Exception as e:
        print(f"  FAIL {acc}: {e}", file=sys.stderr)

def score(q):
    return {acc: max(float(np.dot(q, w)) for w in wins) for acc, wins in awins.items()}

def report(name, q):
    s = score(q)
    ranking = sorted(s.items(), key=lambda kv: -kv[1])
    print(f"\n=== {name} (rank 1 = highest cosine = most 'suspect') ===")
    print(f"{'rank':>4} {'cos':>6}  {'class':<5} name")
    pos_ranks = {}
    for r,(acc,c) in enumerate(ranking,1):
        cl,nm = labels[acc]
        star = "  <-- known off-target" if cl=="POS" else ("  (hard receptor decoy)" if cl=="HARD" else "")
        print(f"{r:>4} {c:>6.3f}  {cl:<5} {nm}{star}")
        if cl=="POS": pos_ranks[nm]=r
    N=len(ranking)
    top6=[nm for nm in POSITIVES.values() if pos_ranks.get(nm,99)<=6]
    hard_ranks=[r for r,(acc,_) in enumerate(ranking,1) if labels[acc][0]=="HARD"]
    best_hard=min(hard_ranks)
    beats_hard = all(pos_ranks.get(nm,99) < best_hard for nm in POSITIVES.values())
    print(f"  -> VEGFR2 rank {pos_ranks.get('VEGFR2')}/{N}, FZD5 {pos_ranks.get('FZD5')}/{N}, ULBP2 {pos_ranks.get('ULBP2')}/{N}")
    print(f"  -> off-targets in top6: {len(top6)}/3 ; best hard-decoy rank: {best_hard} ; all 3 beat every hard decoy: {beats_hard}")
    return pos_ranks

report("F1 whole-Fv vs antigen (single ESM-2 space)", q_whole)
report("F2 CDR paratope (H3+L3) vs antigen", q_para)
print("\nDONE")
