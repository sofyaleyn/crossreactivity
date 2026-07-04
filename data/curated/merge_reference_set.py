"""Merge Layer A + Layer C into the final curated self_proteins.csv.

Per reference-set.md order of operations: Layer A is the base; Layer C rows are
merged into existing Layer-A rows by UniProt ID (setting layer_C=True and the
mimicry fields) rather than duplicated. Layer C proteins not in Layer A (e.g.
intracellular autoantigens like GAD65, Sm B/B') are appended as new rows with
layer_A=False. Anchor flags from Layer A are preserved.
"""
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
LAYER_A = HERE / "self_proteins_layer_a.csv"
LAYER_C = HERE / "layer_c_mimicry_seed.csv"
OUT = HERE / "self_proteins.csv"

SCHEMA = [
    "protein_id", "uniprot_id", "gene_symbol", "name", "sequence",
    "pdb_or_af_ref", "surface_region",
    "layer_A", "layer_B", "layer_C",
    "autoimmune_conditions", "mimicry_pathogens", "mimicry_epitope",
    "mimicry_epitope_type", "is_anchor_offtarget", "source",
    "cspa_validated", "surfy_fpr_class",
]


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def norm(row: dict) -> dict:
    return {k: row.get(k, "") for k in SCHEMA}


def merge_semicolon(a: str, b: str) -> str:
    parts = [p for p in (a.split(";") + b.split(";")) if p.strip()]
    seen, out = set(), []
    for p in parts:
        p = p.strip()
        if p and p not in seen:
            seen.add(p)
            out.append(p)
    return ";".join(out)


def main() -> None:
    layer_a = [norm(r) for r in read_rows(LAYER_A)]
    layer_c = [norm(r) for r in read_rows(LAYER_C)]

    by_uid = {r["uniprot_id"]: r for r in layer_a}

    # Collapse multiple Layer C rows per protein (e.g. MBP has 7 mimicry rows)
    # into one merged record per UniProt ID.
    c_by_uid: dict[str, dict] = {}
    for r in layer_c:
        uid = r["uniprot_id"]
        if uid not in c_by_uid:
            c_by_uid[uid] = dict(r)
        else:
            m = c_by_uid[uid]
            m["autoimmune_conditions"] = merge_semicolon(m["autoimmune_conditions"], r["autoimmune_conditions"])
            m["mimicry_pathogens"] = merge_semicolon(m["mimicry_pathogens"], r["mimicry_pathogens"])
            m["mimicry_epitope"] = (m["mimicry_epitope"] + " | " + r["mimicry_epitope"]).strip(" |")
            if "glycan" in (r["mimicry_epitope_type"], m["mimicry_epitope_type"]):
                m["mimicry_epitope_type"] = "glycan" if m["mimicry_epitope_type"] == "glycan" else m["mimicry_epitope_type"]
            m["source"] = merge_semicolon(m["source"], r["source"])

    merged_into_a, appended = 0, 0
    for uid, c in c_by_uid.items():
        if uid in by_uid:
            a = by_uid[uid]
            a["layer_C"] = "True"
            a["layer_B"] = "True" if c["layer_B"] == "True" else a["layer_B"]
            a["autoimmune_conditions"] = merge_semicolon(a["autoimmune_conditions"], c["autoimmune_conditions"])
            a["mimicry_pathogens"] = merge_semicolon(a["mimicry_pathogens"], c["mimicry_pathogens"])
            a["mimicry_epitope"] = c["mimicry_epitope"]
            a["mimicry_epitope_type"] = c["mimicry_epitope_type"]
            a["source"] = merge_semicolon(a["source"], c["source"])
            merged_into_a += 1
        else:
            # not in surfaceome (intracellular autoantigen) -> append
            by_uid[uid] = c
            appended += 1

    rows = list(by_uid.values())

    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SCHEMA)
        writer.writeheader()
        writer.writerows(rows)

    n_a = sum(1 for r in rows if r["layer_A"] == "True")
    n_b = sum(1 for r in rows if r["layer_B"] == "True")
    n_c = sum(1 for r in rows if r["layer_C"] == "True")
    n_anchor = sum(1 for r in rows if r["is_anchor_offtarget"] == "True")
    n_seq = sum(1 for r in rows if r["sequence"])
    print(f"Merged reference set: {len(rows)} unique proteins")
    print(f"  Layer A: {n_a}   Layer B: {n_b}   Layer C: {n_c}")
    print(f"  Layer C merged into A rows: {merged_into_a}, appended (non-surfaceome): {appended}")
    print(f"  anchor off-targets flagged: {n_anchor}")
    print(f"  with sequences: {n_seq}/{len(rows)}")
    for acc, nm in [("P35968", "VEGFR2"), ("Q13467", "FZD5"), ("Q9BZM5", "ULBP2")]:
        r = by_uid.get(acc)
        print(f"  {nm}: {'present' if r else 'MISSING'}, anchor={r['is_anchor_offtarget'] if r else '-'}")


if __name__ == "__main__":
    main()
