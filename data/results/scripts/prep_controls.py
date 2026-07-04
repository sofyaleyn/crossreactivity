import json, urllib.request

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL_WT = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"

def uniprot(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=40) as r:
        return json.load(r)

def domain_region(j, needle):
    for f in j.get("features", []):
        if f.get("type") == "Domain" and needle in f.get("description", ""):
            return f["location"]["start"]["value"], f["location"]["end"]["value"]
    return None

# PD-1 (PDCD1, Q15116) — camrelizumab's real target; Ig-like V-type ectodomain = confident-binder control
pj = uniprot("Q15116"); pseq = pj["sequence"]["value"]
reg = domain_region(pj, "Ig-like V-type")
if reg: s, e = reg
else:   s, e = 35, 145
# pad a little around the IgV domain to include the full folded ectodomain the mAb engages
s = max(1, s - 8); e = min(len(pseq), e + 8)
PD1 = pseq[s-1:e]
print(f"PD-1 IgV region {s}-{e} len {len(PD1)}")

# Hen egg-white lysozyme (P00698) — classic irrelevant antigen = non-binder noise floor
lj = uniprot("P00698"); lseq = lj["sequence"]["value"]
# strip signal peptide (1-18); mature 19-147
sig = None
for f in lj.get("features", []):
    if f.get("type") == "Signal":
        sig = f["location"]["end"]["value"]
mstart = (sig + 1) if sig else 19
LYZ = lseq[mstart-1:]
print(f"Lysozyme mature {mstart}-{len(lseq)} len {len(LYZ)}")

def make_input(antigen, num_samples=5):
    return {
        "entities": [
            {"type": "protein", "chain_ids": ["H"], "value": VH},
            {"type": "protein", "chain_ids": ["L"], "value": VL_WT},
            {"type": "protein", "chain_ids": ["V"], "value": antigen},
        ],
        "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["H", "L"]},
        "num_samples": num_samples,
        "model_options": {"recycling_steps": 3, "sampling_steps": 200},
    }

for tag, ag in [("pd1", PD1), ("lyz", LYZ)]:
    json.dump(make_input(ag), open(f"{BASE}/cofold_{tag}.json", "w"), indent=2)
    print("wrote", f"{BASE}/cofold_{tag}.json", "| complex residues:", len(VH)+len(VL_WT)+len(ag))
