import json, urllib.request

# ---- Antibody: SHR-1210 (camrelizumab) WT Fv ----
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
VL_WT = "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"

# WT VL decomposed (verified vs Finlay Table 1 CDRs: L1 LASQTIGTWLT, L2 TATSLAD, L3 QQVYSIPWT)
FR1 = "DIQMTQSPSSLSASVGDRVTITC"
L1_WT, FR2 = "LASQTIGTWLT", "WYQQKPGKAPKLLIY"
L2_WT, FR3 = "TATSLAD", "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC"
L3_WT, FR4 = "QQVYSIPWT", "FGGGTKVEIK"
assert FR1+L1_WT+FR2+L2_WT+FR3+L3_WT+FR4 == VL_WT, "WT VL region decomposition mismatch"

# ---- Reconstructed germlined mutant: revert L-CDRs to IGKV1-39 germline ----
# (Finlay: light chain near-germlined to IGKV1-39; L-CDR changes drove VEGFR2 ablation.)
L1_GL = "RASQSISSYLN"   # IGKV1-39 CDR-L1
L2_GL = "AASSLQS"       # IGKV1-39 CDR-L2
L3_GL = "QQSYSTPWT"     # IGKV1-39 CDR-L3 (QQSYSTP) + Trp (JK1-derived); WT FR4 retained
VL_GL = FR1+L1_GL+FR2+L2_GL+FR3+L3_GL+FR4
print("VH        len", len(VH))
print("VL_WT     len", len(VL_WT))
print("VL_GL     len", len(VL_GL), "->", VL_GL)
print("  L-CDR reversions: L1", L1_WT, "->", L1_GL, "| L2", L2_WT, "->", L2_GL, "| L3", L3_WT, "->", L3_GL)

# ---- Antigen: VEGFR2 (P35968) Ig-like domains 2-3 (VEGF-binding region) ----
seq = None; ig_domains = []
try:
    with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/P35968.json", timeout=40) as r:
        j = json.load(r)
    seq = j["sequence"]["value"]
    for f in j.get("features", []):
        if f.get("type") == "Domain" and "Ig-like" in f.get("description", ""):
            b = f["location"]["start"]["value"]; e = f["location"]["end"]["value"]
            ig_domains.append((f["description"], b, e))
    print("\nVEGFR2 full length", len(seq))
    for d in ig_domains: print("  domain", d)
except Exception as ex:
    print("UniProt fetch failed:", ex)

# Pick Ig-like domain 2 start .. domain 3 end (VEGF-A binding site = D2-D3)
d2 = next((d for d in ig_domains if d[0].rstrip().endswith(" 2")), None)
d3 = next((d for d in ig_domains if d[0].rstrip().endswith(" 3")), None)
if seq and d2 and d3:
    start, end = d2[1], d3[2]
else:  # fallback to canonical D2-3 construct boundaries
    start, end = 121, 327
    if not seq:
        raise SystemExit("no VEGFR2 sequence; cannot proceed")
VEGFR2_D23 = seq[start-1:end]
print(f"\nVEGFR2 D2-D3 residues {start}-{end}  len {len(VEGFR2_D23)}")
print(VEGFR2_D23)

# ---- Build Boltz inputs (antibody = protein-protein binder vs antigen) ----
def make_input(vl, tag):
    return {
        "entities": [
            {"type": "protein", "chain_ids": ["H"], "value": VH},
            {"type": "protein", "chain_ids": ["L"], "value": vl},
            {"type": "protein", "chain_ids": ["V"], "value": VEGFR2_D23},
        ],
        "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["H", "L"]},
        "num_samples": 1,
        "model_options": {"recycling_steps": 3, "sampling_steps": 200},
        "_meta": {"tag": tag},
    }

base = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"
for tag, vl in [("wt", VL_WT), ("mut", VL_GL)]:
    d = make_input(vl, tag)
    d.pop("_meta")
    json.dump(d, open(f"{base}/cofold_{tag}.json", "w"), indent=2)
    print("wrote", f"{base}/cofold_{tag}.json")
print("\nComplex size (residues):", len(VH)+len(VL_WT)+len(VEGFR2_D23))
