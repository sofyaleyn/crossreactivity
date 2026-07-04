import json
BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"

# Pembrolizumab Fv (from PDB 5GGS pembrolizumab-PD-1 complex; CDR-H3 RDYRFDMGFDY, CDR-L3 QHSRDLPLT)
PEMBRO_VH = "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
PEMBRO_VL = "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK"
assert "RDYRFDMGFDY" in PEMBRO_VH and "QHSRDLPLT" in PEMBRO_VL

# reuse the SAME antigen constructs SHR-1210 was tested against
SRC = {"pd1":"cofold_pd1.json", "vegfr2":"cofold_wt.json", "fzd5":"cofold_fzd5.json", "ulbp2":"cofold_ulbp2.json"}
for tag, f in SRC.items():
    d = json.load(open(f"{BASE}/{f}"))
    V = next(e for e in d["entities"] if e["chain_ids"] == ["V"])["value"]
    out = {"entities":[{"type":"protein","chain_ids":["H"],"value":PEMBRO_VH},
                       {"type":"protein","chain_ids":["L"],"value":PEMBRO_VL},
                       {"type":"protein","chain_ids":["V"],"value":V}],
           "binding":{"type":"protein_protein_binding","binder_chain_ids":["H","L"]},
           "num_samples":5,"model_options":{"recycling_steps":3,"sampling_steps":200}}
    json.dump(out, open(f"{BASE}/cofold_pembro_{tag}.json","w"), indent=2)
    print(f"cofold_pembro_{tag}.json  antigen len {len(V)}  complex {len(PEMBRO_VH)+len(PEMBRO_VL)+len(V)}")
