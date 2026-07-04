import json, urllib.request

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL_WT = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"

DECOYS = [
    ("Q9UP38","FZD1"),("Q14332","FZD2"),("Q9NPG1","FZD3"),("Q9ULV1","FZD4"),
    ("O60353","FZD6"),("O75084","FZD7"),("Q9H461","FZD8"),("O00144","FZD9"),
    ("Q9ULW2","FZD10"),("Q99835","SMO"),("Q8N474","SFRP1"),
]

def uni(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60) as r:
        return json.load(r)

def feats(j, types):
    out = []
    for f in j.get("features", []):
        if f.get("type") in types:
            out.append((f["type"], f.get("description",""),
                        f["location"]["start"]["value"], f["location"]["end"]["value"]))
    return out

def make_input(ag):
    return {"entities":[{"type":"protein","chain_ids":["H"],"value":VH},
                        {"type":"protein","chain_ids":["L"],"value":VL_WT},
                        {"type":"protein","chain_ids":["V"],"value":ag}],
            "binding":{"type":"protein_protein_binding","binder_chain_ids":["H","L"]},
            "num_samples":5,"model_options":{"recycling_steps":3,"sampling_steps":200}}

results = {}
for acc, name in DECOYS:
    try:
        j = uni(acc); seq = j["sequence"]["value"]
    except Exception as e:
        print(f"{name} ({acc}) FETCH FAILED: {e}"); results[name]=None; continue
    print(f"{name} ({acc}) len {len(seq)}")
    crd = None
    # Prefer a Domain feature whose description mentions FZ/Frizzled/cysteine-rich
    for t,desc,s,e in feats(j, {"Domain"}):
        print(f"  DOMAIN '{desc}' {s}-{e}")
        if crd is None and ("FZ" in desc or "Frizzled" in desc or "cysteine-rich" in desc.lower()):
            crd = (s,e,"Domain:"+desc)
    # Fallback to Region
    if crd is None:
        for t,desc,s,e in feats(j, {"Region"}):
            if ("FZ" in desc or "Frizzled" in desc or "cysteine-rich" in desc.lower()):
                print(f"  REGION fallback '{desc}' {s}-{e}")
                crd = (s,e,"Region:"+desc); break
    if crd is None:
        print(f"  !! NO FZ DOMAIN FOUND for {name} -- skipping")
        results[name]=None; continue
    s,e,src = crd
    ag = seq[s-1:e]
    print(f"  -> {name} CRD {s}-{e} len {len(ag)} [{src}]")
    fn = f"{BASE}/cofold_{name.lower()}.json"
    json.dump(make_input(ag), open(fn,"w"), indent=2)
    results[name] = dict(acc=acc, start=s, end=e, length=len(ag), src=src,
                         complex_res=len(VH)+len(VL_WT)+len(ag), json=fn)

json.dump(results, open(f"{BASE}/decoy_domains.json","w"), indent=2)
print("\nSUMMARY")
for name,r in results.items():
    if r: print(f"  {name}: {r['start']}-{r['end']} len {r['length']}")
    else: print(f"  {name}: SKIPPED (no domain)")
