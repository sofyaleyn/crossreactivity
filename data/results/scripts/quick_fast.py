import json
SCRATCH="/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
POS={"P35968":"VEGFR2","Q13467":"FZD5","Q9BZM5":"ULBP2"}
fast=json.load(open(f"{SCRATCH}/fast_scores.json")); N=len(fast)
# confirm positives present
present=[u for u in POS if any(d["uniprot_id"]==u for d in fast)]
print(f"N={N}; positives present: {[POS[u] for u in present]} (missing: {[POS[u] for u in POS if u not in present]})\n")
def rep(label,key):
    s=sorted(fast,key=lambda d:d[key],reverse=True)
    rk={d['uniprot_id']:i+1 for i,d in enumerate(s)}
    pcts=[]; cells=[]
    for u in POS:
        r=rk[u]; pcts.append(r/N); cells.append(f"{POS[u]} {r}({r/N*100:.0f}%)")
    worst=max(pcts); tk=lambda p:sum(1 for x in pcts if x<=p)
    print(f"{label:26s} | {' | '.join(cells):50s} | top10%:{tk(.10)}/3 | RETAIN-ALL-3 {worst:.3f}")
for label,key in [("naive_cdr_identity","naive_cdr_identity"),
                  ("biophys_hydrophobicity","biophys_hydrophobicity"),
                  ("biophys_charge","biophys_charge"),
                  ("biophys_length","biophys_length"),
                  ("annotation_receptor(FREE)","annotation_receptor_prior")]:
    rep(label,key)
print("\n(RETAIN-ALL-3 = fraction of the 2896 you must keep to hold all 3 known off-targets; lower=better. Embedding filters pending.)")
