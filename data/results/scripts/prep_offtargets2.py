import json, urllib.request

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL_WT = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"

def uni(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=40) as r:
        return json.load(r)

def feats(j, types):
    out = []
    for f in j.get("features", []):
        if f.get("type") in types:
            out.append((f["type"], f.get("description",""), f["location"]["start"]["value"], f["location"]["end"]["value"]))
    return out

# ---- FZD5 (Q13467): extracellular Frizzled cysteine-rich domain (CRD) ----
fj = uni("Q13467"); fseq = fj["sequence"]["value"]
print("FZD5 len", len(fseq))
crd = None
for t,desc,s,e in feats(fj, {"Domain","Region"}):
    print("  FZD5", t, desc, s, e)
    if ("FZ" in desc or "Frizzled" in desc or "cysteine-rich" in desc.lower()) and crd is None:
        crd = (s,e)
if crd is None:  # fallback: known FZD5 CRD ~ 27-157
    crd = (27,157)
FZD5 = fseq[crd[0]-1:crd[1]]
print(f"  -> FZD5 CRD {crd[0]}-{crd[1]} len {len(FZD5)}")

# ---- ULBP2 (Q9BZM5): mature MHC-class-I-like ectodomain (strip signal + C-term GPI region) ----
uj = uni("Q9BZM5"); useq = uj["sequence"]["value"]
print("ULBP2 len", len(useq))
sig = None; prop_start = None
for t,desc,s,e in feats(uj, {"Signal","Propeptide","Chain","Domain","Region"}):
    print("  ULBP2", t, desc, s, e)
    if t=="Signal": sig = e
    if t=="Propeptide": prop_start = s
mstart = (sig+1) if sig else 29
mend = (prop_start-1) if prop_start else 213   # drop C-terminal GPI-anchor signal
ULBP2 = useq[mstart-1:mend]
print(f"  -> ULBP2 ectodomain {mstart}-{mend} len {len(ULBP2)}")

def make_input(ag):
    return {"entities":[{"type":"protein","chain_ids":["H"],"value":VH},
                        {"type":"protein","chain_ids":["L"],"value":VL_WT},
                        {"type":"protein","chain_ids":["V"],"value":ag}],
            "binding":{"type":"protein_protein_binding","binder_chain_ids":["H","L"]},
            "num_samples":5,"model_options":{"recycling_steps":3,"sampling_steps":200}}

for tag,ag in [("fzd5",FZD5),("ulbp2",ULBP2)]:
    json.dump(make_input(ag), open(f"{BASE}/cofold_{tag}.json","w"), indent=2)
    print("wrote", f"cofold_{tag}.json | complex residues:", len(VH)+len(VL_WT)+len(ag))
