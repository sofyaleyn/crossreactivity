import json, os, sys, urllib.request, urllib.error, warnings, numpy as np, torch, esm
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
PDBDIR = f"{BASE}/af_pdbs"; os.makedirs(PDBDIR, exist_ok=True)

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
PARATOPE_SUBS = ["QLYYFDYW", "LASQTIGTWLT", "TATSLAD", "QQVYSIPWT"]

# Kyte-Doolittle hydropathy
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
CHARGE = {'K':1.0,'R':1.0,'D':-1.0,'E':-1.0,'H':0.1}  # H partially protonated
# Tien et al. 2013 theoretical max ASA (A^2) for relative SASA
MAXASA = {'A':129.0,'R':274.0,'N':195.0,'D':193.0,'C':167.0,'Q':225.0,'E':223.0,
          'G':104.0,'H':224.0,'I':197.0,'L':201.0,'K':236.0,'M':224.0,'F':240.0,
          'P':159.0,'S':155.0,'T':172.0,'W':285.0,'Y':263.0,'V':174.0}
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
             'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
             'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
RSA_THRESH = 0.20  # relative SASA cutoff for "surface"

print("Loading ESM-2 t12_35M ...", flush=True)
model, alph = esm.pretrained.esm2_t12_35M_UR50D(); model.eval()
bc = alph.get_batch_converter(); LAYER = 12

def per_res_full(seq):
    """Per-residue ESM reps aligned to seq positions 0..N-1 (non-overlapping chunks)."""
    outs = []
    for i in range(0, len(seq), 1000):
        sub = seq[i:i+1000]
        _,_,tok = bc([("x", sub)])
        with torch.no_grad():
            rep = model(tok, repr_layers=[LAYER])["representations"][LAYER][0]
        outs.append(rep[1:len(sub)+1].numpy())
    return np.concatenate(outs, 0)  # (N, D)

def paratope_vec():
    fv = VH+VL; mat = per_res_full(fv); idx = []
    for s in PARATOPE_SUBS:
        p = fv.find(s); idx += list(range(p, p+len(s))) if p >= 0 else []
    v = mat[idx].mean(0); return v/(np.linalg.norm(v)+1e-9)

def paratope_physchem():
    res = "".join(PARATOPE_SUBS)
    q = sum(CHARGE.get(a,0.0) for a in res)
    h = np.mean([KD[a] for a in res])
    return q, h

def download(acc):
    fn = f"{PDBDIR}/AF-{acc}.pdb"
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        return fn
    # try v4 first (per task); AF-DB has since bumped to v6, so fall back to the
    # API-resolved pdbUrl for the canonical isoform.
    urls = [f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"]
    try:
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=30) as r:
            preds = json.load(r)
        for p in preds:
            if p.get("uniprotAccession") == acc and p.get("pdbUrl"):
                urls.append(p["pdbUrl"]); break
        if len(urls) == 1 and preds and preds[0].get("pdbUrl"):
            urls.append(preds[0]["pdbUrl"])
    except Exception:
        pass
    for url in urls:
        try:
            urllib.request.urlretrieve(url, fn)
            if os.path.getsize(fn) > 500:
                return fn
        except Exception:
            if os.path.exists(fn): os.remove(fn)
    if os.path.exists(fn): os.remove(fn)
    return None

def surface_mask(pdb_fn, seq):
    """Return (seq_positions_on_surface set, list of (one_letter, rsa) for surface residues).
    Maps PDB residue number i -> seq index i-1 (AF numbering == canonical seq)."""
    p = PDBParser(QUIET=True)
    struct = p.get_structure("x", pdb_fn)
    ShrakeRupley().compute(struct, level="R")
    model0 = struct[0]
    surf_pos = set(); surf_res = []
    matched = 0; total = 0
    for chain in model0:
        for res in chain:
            if not is_aa(res, standard=True): continue
            aa = THREE2ONE.get(res.resname)
            if aa is None: continue
            resnum = res.id[1]; total += 1
            idx = resnum - 1
            if 0 <= idx < len(seq) and seq[idx] == aa:
                matched += 1
            sasa = res.sasa
            rsa = sasa / MAXASA[aa]
            if rsa > RSA_THRESH:
                surf_pos.add(idx)
                surf_res.append((aa, rsa))
    return surf_pos, surf_res, matched, total

def within_family_aurocs(scores, labels):
    fams = {}
    for a,(fam,name,isp) in labels.items():
        fams.setdefault(fam, []).append(a)
    out = {}
    for fam, mem in fams.items():
        mem = [a for a in mem if a in scores]
        pos = [a for a in mem if labels[a][2]]
        if not pos: out[fam] = None; continue
        pos = pos[0]
        dec = [a for a in mem if not labels[a][2]]
        out[fam] = sum(scores[pos] > scores[a] for a in dec)/len(dec)
    return out

# ---------------- main ----------------
antigens = json.load(open(f"{BASE}/benchmark_antigens.json"))
labels = {a:(v['family'], v['name'], v['is_pos']) for a,v in antigens.items()}

q_vec = paratope_vec()
q_charge, q_hyd = paratope_physchem()
print(f"Paratope net charge={q_charge:+.1f}  mean KD hydropathy={q_hyd:+.2f}", flush=True)

score_surf_esm = {}      # (a) surface-restricted ESM cosine
score_whole_esm = {}     # sanity: whole-protein ESM cosine (single mean, no windows)
score_charge_comp = {}   # (b1) charge complementarity only
score_phys_comp = {}     # (b) charge+hydrophobicity complementarity
score_surf_frac = {}     # extra descriptor
skipped = []

for acc, meta in antigens.items():
    seq = meta['seq']
    fn = download(acc)
    if fn is None:
        skipped.append((acc, meta['name'], "no AF model"))
        print(f"  SKIP {acc} {meta['name']}: no AF model", flush=True)
        continue
    try:
        surf_pos, surf_res, matched, total = surface_mask(fn, seq)
    except Exception as e:
        skipped.append((acc, meta['name'], f"parse err {e}"))
        print(f"  SKIP {acc} {meta['name']}: {e}", flush=True)
        continue
    if len(surf_res) < 5:
        skipped.append((acc, meta['name'], "too few surface res"))
        continue

    # ESM per-residue reps aligned to seq
    mat = per_res_full(seq)  # (N,D)
    N = mat.shape[0]
    sp = sorted(i for i in surf_pos if i < N)
    # (a) surface-restricted embedding
    v_surf = mat[sp].mean(0); v_surf = v_surf/(np.linalg.norm(v_surf)+1e-9)
    score_surf_esm[acc] = float(np.dot(q_vec, v_surf))
    # whole-protein single-mean embedding (comparison)
    v_all = mat.mean(0); v_all = v_all/(np.linalg.norm(v_all)+1e-9)
    score_whole_esm[acc] = float(np.dot(q_vec, v_all))

    # physicochemical of surface
    aas = [a for a,_ in surf_res]
    net_charge = sum(CHARGE.get(a,0.0) for a in aas)
    net_charge_per = net_charge/len(aas)          # charge density
    mean_hyd = np.mean([KD[a] for a in aas])
    # complementarity: opposite charge favorable -> -(q_para*q_anti)
    charge_comp = -(q_charge/abs(q_charge)) * net_charge_per   # paratope neg -> favor +antigen
    # hydrophobic patches pair like-with-like -> product of (mean-centered) hydropathy
    hyd_comp = (q_hyd) * mean_hyd  # both hydrophobic (positive KD) favorable
    # normalize hyd term scale roughly
    score_charge_comp[acc] = float(charge_comp)
    score_phys_comp[acc] = float(charge_comp + 0.15*hyd_comp)
    score_surf_frac[acc] = float(len(surf_res)/max(matched,1))

    tag = "POS" if meta['is_pos'] else "   "
    print(f"  {tag} {meta['family']:16s} {meta['name']:8s} {acc}  Nsurf={len(surf_res):4d} "
          f"qdens={net_charge_per:+.3f} hyd={mean_hyd:+.2f} "
          f"surfESM={score_surf_esm[acc]:.3f} chgComp={charge_comp:+.3f}", flush=True)

# ---------------- report ----------------
scorers = {
    "esm_paratope_cosine (BASELINE, windowed)": json.load(open(f"{BASE}/benchmark_embed_scores.json"))["scores"],
    "surface_ESM_cosine": score_surf_esm,
    "wholeprotein_ESM_singlemean": score_whole_esm,
    "surface_charge_complementarity": score_charge_comp,
    "surface_physchem_complementarity": score_phys_comp,
}
famorder = ["RTK/Ig-receptor","Frizzled/CRD","MHC-I-like"]
print("\n" + "="*90)
print(f"{'SCORER':40s} " + " ".join(f"{f[:10]:>11s}" for f in famorder) + f"{'MEAN':>8s}")
print("-"*90)
results = {}
for name, sc in scorers.items():
    au = within_family_aurocs(sc, labels)
    row = [au.get(f) for f in famorder]
    valid = [x for x in row if x is not None]
    mean = np.mean(valid) if valid else float('nan')
    results[name] = {"per_family": au, "mean": mean}
    print(f"{name:40s} " + " ".join(f"{(x if x is not None else float('nan')):>11.3f}" for x in row) + f"{mean:>8.3f}")
print("="*90)
print(f"Skipped (missing AF model / parse): {len(skipped)}")
for s in skipped: print("   ", s)

json.dump({
    "rsa_threshold": RSA_THRESH,
    "paratope_charge": q_charge, "paratope_hydropathy": q_hyd,
    "scores": {k: v for k, v in {
        "surface_ESM_cosine": score_surf_esm,
        "wholeprotein_ESM_singlemean": score_whole_esm,
        "surface_charge_complementarity": score_charge_comp,
        "surface_physchem_complementarity": score_phys_comp,
    }.items()},
    "within_family_auroc": {k: {"per_family": v["per_family"], "mean": v["mean"]} for k,v in results.items()},
    "skipped": skipped,
    "labels": {a: list(labels[a]) for a in labels},
}, open(f"{BASE}/benchmark_surface_scores.json","w"), indent=2)
print("saved -> benchmark_surface_scores.json")
