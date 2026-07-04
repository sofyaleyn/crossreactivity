"""Build the background calibration set (reference-set.md Background section).

Samples Layer-A surfaceome proteins that are abundant housekeeping / transport /
adhesion families (Na/K-ATPase, GLUT/SLC transporters, integrins, tetraspanins,
etc.) and NOT anchor off-targets. Used to compute enrichment normalization in
rank/embedding_rank.py. Output: ../background/benign_proteins.csv
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
LAYER_A = HERE / "self_proteins_layer_a.csv"
OUT_DIR = HERE.parent / "background"
OUT_CSV = OUT_DIR / "benign_proteins.csv"

TARGET_SIZE = 400

# Gene-symbol prefixes for abundant, well-characterized housekeeping surface
# families unlikely to be documented CDR off-targets.
BENIGN_PREFIXES = (
    "ATP1A", "ATP1B",   # Na/K-ATPase subunits
    "SLC2",             # GLUT transporters
    "SLC",              # broader solute carriers
    "ITGA", "ITGB",     # integrins
    "TSPAN", "CD9", "CD63", "CD81", "CD151",  # tetraspanins
    "CDH",              # cadherins
    "TFRC",             # transferrin receptor
    "BSG",              # basigin
    "CD44", "CD47", "CD59",
    "AQP",              # aquaporins
)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with LAYER_A.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    def is_benign(r: dict) -> bool:
        if r["is_anchor_offtarget"] == "True":
            return False
        if not r["sequence"]:
            return False
        gene = (r["gene_symbol"] or "").upper()
        return any(gene.startswith(p) for p in BENIGN_PREFIXES)

    picked = [r for r in rows if is_benign(r)]
    # de-dup by gene, keep breadth
    seen, deduped = set(), []
    for r in picked:
        g = r["gene_symbol"]
        if g not in seen:
            seen.add(g)
            deduped.append(r)
    benign = deduped[:TARGET_SIZE]

    for r in benign:
        r["source"] = ("Background calibration set: SURFY Layer-A housekeeping/"
                       "transport/adhesion family, no documented CDR off-target")

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(benign)

    print(f"Wrote {len(benign)} background proteins to {OUT_CSV}")
    fams: dict[str, int] = {}
    for r in benign:
        g = (r["gene_symbol"] or "?")[:4]
        fams[g] = fams.get(g, 0) + 1
    top = sorted(fams.items(), key=lambda x: -x[1])[:12]
    print("Top family prefixes:", top)


if __name__ == "__main__":
    main()
