"""Merge UniProt sequences into the Layer C hand-curated mimicry seed table.

Reads layer_c_mimicry_seed.csv (rows transcribed from Suliman 2024 Table 1
and Almulla et al. 2025), fetches/reads cached FASTA files in raw_fasta/,
and writes layer_c_mimicry_seed_filled.csv with the sequence column populated.
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
RAW_FASTA = HERE / "raw_fasta"
IN_CSV = HERE / "layer_c_mimicry_seed.csv"
OUT_CSV = HERE / "layer_c_mimicry_seed_filled.csv"


def read_fasta_seq(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def load_sequences() -> dict[str, str]:
    seqs = {}
    for fasta_path in RAW_FASTA.glob("*.fasta"):
        accession = fasta_path.stem
        seqs[accession] = read_fasta_seq(fasta_path)
    return seqs


def main() -> None:
    sequences = load_sequences()

    with IN_CSV.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    filled, missing = 0, []
    for row in rows:
        uniprot_id = row["uniprot_id"]
        if uniprot_id in sequences:
            row["sequence"] = sequences[uniprot_id]
            filled += 1
        elif uniprot_id.startswith("NA"):
            row["sequence"] = ""  # glycan entry, no protein sequence by design
        else:
            missing.append((row["protein_id"], uniprot_id))

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Filled {filled}/{len(rows)} rows with sequences.")
    if missing:
        print("Missing sequences for:")
        for protein_id, uid in missing:
            print(f"  {protein_id} ({uid})")


if __name__ == "__main__":
    main()
