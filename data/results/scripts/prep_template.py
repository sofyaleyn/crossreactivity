import json, urllib.request, warnings
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
warnings.filterwarnings("ignore")

BASE = "/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad"

# VEGFR2 D2-3 (must match what we cofolded: P35968 residues 141-320)
with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/P35968.json", timeout=40) as r:
    vseq = json.load(r)["sequence"]["value"]
VEGFR2_D23 = vseq[140:320]
probe = VEGFR2_D23[40:75]  # a 35-mer to identify the template chain

# Download 3V2A.pdb and find which chain is VEGFR2
pdb_path = f"{BASE}/3V2A.pdb"
urllib.request.urlretrieve("https://files.rcsb.org/download/3V2A.pdb", pdb_path)
model = PDBParser(QUIET=True).get_structure("3v2a", pdb_path)[0]
tmpl_chain = None
for c in model:
    resseq = "".join(seq1(r.resname, undef_code="") for r in c if r.id[0] == " ")
    hit = probe in resseq
    print(f"chain {c.id}: {len([r for r in c if r.id[0]==' '])} res  VEGFR2-probe={'YES' if hit else 'no'}")
    if hit and tmpl_chain is None:
        tmpl_chain = c.id
print("=> VEGFR2 template chain:", tmpl_chain)
assert tmpl_chain, "could not locate VEGFR2 chain in 3V2A"

tmpl = {
    "template_chains": [{"input_chain_id": "V", "template_chain_id": tmpl_chain}],
    "template_structure": {"type": "url", "url": "https://files.rcsb.org/download/3V2A.pdb"},
}
for tag in ["wt", "mut"]:
    d = json.load(open(f"{BASE}/cofold_{tag}.json"))
    d["templates"] = [tmpl]
    json.dump(d, open(f"{BASE}/cofold_{tag}_tmpl.json", "w"), indent=2)
    print("wrote", f"{BASE}/cofold_{tag}_tmpl.json")
