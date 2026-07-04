"""Build Layer A (base surfaceome breadth) of the CrossFlag reference set.

Source: SURFY in silico human surfaceome (Bausch-Fluck et al., PNAS 2018),
'surface'-labeled proteins from table S3, plus CSPA mass-spec validation flag
(Bausch-Fluck et al., PLOS ONE 2015).

Steps (per reference-set.md Layer A section):
  1. Parse SURFY table S3 -> 'surface'-labeled UniProt IDs.
  2. Flag CSPA-validated subset (cspa_validated).
  3. Batch-fetch canonical sequences from UniProt REST.
  4. Inject anchors (VEGFR2/FZD5/ULBP2) if missing; set is_anchor_offtarget.
  5. Write self_proteins_layer_a.csv in the target schema.
"""
import csv
import io
import time
import warnings

import pandas as pd
import requests

from . import paths

warnings.filterwarnings("ignore")

SURFY_XLSX = paths.SURFY_XLSX
CSPA_XLSX = paths.CSPA_XLSX
SEQ_CACHE = paths.SEQ_CACHE
OUT_CSV = paths.LAYER_A_CSV

FPR_COL = "MachineLearning FPR class (1=1%, 2=5%, 3=15%)"

ANCHORS = {
    "P35968": ("KDR", "VEGFR2"),
    "Q13467": ("FZD5", "FZD5"),
    "Q9BZM5": ("ULBP2", "ULBP2"),
}

SCHEMA = [
    "protein_id", "uniprot_id", "gene_symbol", "name", "sequence",
    "pdb_or_af_ref", "surface_region",
    "layer_A", "layer_B", "layer_C",
    "autoimmune_conditions", "mimicry_pathogens", "mimicry_epitope",
    "mimicry_epitope_type", "is_anchor_offtarget", "source",
    "cspa_validated", "surfy_fpr_class",
]


def load_surfy() -> pd.DataFrame:
    df = pd.read_excel(SURFY_XLSX, sheet_name="SurfaceomeMasterTable", header=1)
    surf = df[df["Surfaceome Label"] == "surface"].copy()
    surf = surf[["UniProt accession", "UniProt gene", "UniProt description", FPR_COL]]
    surf.columns = ["uniprot_id", "gene_symbol", "name", "surfy_fpr_class"]
    surf = surf.dropna(subset=["uniprot_id"]).drop_duplicates("uniprot_id")
    return surf


def load_cspa_accessions() -> set[str]:
    """CSPA mass-spec-validated surfaceome UniProt IDs (for cspa_validated flag)."""
    try:
        xls = pd.ExcelFile(CSPA_XLSX)
        accs: set[str] = set()
        for sheet in xls.sheet_names:
            df = pd.read_excel(CSPA_XLSX, sheet_name=sheet, header=None, nrows=200)
            # scan for any column that looks like UniProt accessions
            for col in df.columns:
                vals = df[col].astype(str)
                hits = vals.str.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
                if hits.sum() > 20:
                    full = pd.read_excel(CSPA_XLSX, sheet_name=sheet, header=None)
                    accs |= set(full[col].astype(str)[full[col].astype(str).str.fullmatch(
                        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")])
        return accs
    except Exception as e:  # noqa: BLE001
        print(f"  (CSPA parse skipped: {e})")
        return set()


def load_seq_cache() -> dict[str, str]:
    if not SEQ_CACHE.exists():
        return {}
    cache = {}
    with SEQ_CACHE.open() as f:
        for line in f:
            acc, seq = line.rstrip("\n").split("\t")
            cache[acc] = seq
    return cache


def save_seq_cache(cache: dict[str, str]) -> None:
    with SEQ_CACHE.open("w") as f:
        for acc, seq in cache.items():
            f.write(f"{acc}\t{seq}\n")


def fetch_sequences(accessions: list[str], cache: dict[str, str], batch=100) -> dict[str, str]:
    """Batch-fetch canonical sequences via UniProt REST stream endpoint."""
    todo = [a for a in accessions if a not in cache]
    print(f"  {len(cache)} cached, {len(todo)} to fetch")
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {"query": query, "format": "fasta"}
        for attempt in range(4):
            try:
                r = requests.get(url, params=params, timeout=120)
                r.raise_for_status()
                break
            except Exception as e:  # noqa: BLE001
                print(f"    retry {attempt+1} for batch {i//batch}: {e}")
                time.sleep(2 * (attempt + 1))
        else:
            print(f"    FAILED batch {i//batch}, skipping")
            continue
        # parse FASTA
        acc, seq_parts = None, []
        for line in io.StringIO(r.text):
            if line.startswith(">"):
                if acc:
                    cache[acc] = "".join(seq_parts)
                # header: >sp|P12345|NAME ...  -> accession is field 1
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) > 1 else parts[0].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if acc:
            cache[acc] = "".join(seq_parts)
        if (i // batch) % 5 == 0:
            save_seq_cache(cache)
            print(f"    fetched through batch {i//batch} ({len(cache)} total cached)")
        time.sleep(0.4)
    save_seq_cache(cache)
    return cache


def main() -> None:
    paths.ensure_dirs()
    print("Loading SURFY surface proteins...")
    surf = load_surfy()
    print(f"  {len(surf)} surface-labeled proteins")

    print("Loading CSPA validation set...")
    cspa = load_cspa_accessions()
    print(f"  {len(cspa)} CSPA accessions")

    accessions = list(surf["uniprot_id"])
    for a in ANCHORS:
        if a not in accessions:
            accessions.append(a)

    print("Fetching sequences from UniProt (batched, cached)...")
    cache = load_seq_cache()
    cache = fetch_sequences(accessions, cache)

    print("Writing Layer A CSV...")
    surf_lookup = surf.set_index("uniprot_id").to_dict("index")
    rows = []
    for acc in accessions:
        meta = surf_lookup.get(acc, {})
        gene = meta.get("gene_symbol", "")
        name = meta.get("name", "")
        fpr = meta.get("surfy_fpr_class", "")
        is_anchor = acc in ANCHORS
        if is_anchor and acc not in surf_lookup:
            gene, name = ANCHORS[acc][0], f"{ANCHORS[acc][1]} (anchor, injected)"
        rows.append({
            "protein_id": f"LAYERA_{gene or acc}_{acc}",
            "uniprot_id": acc,
            "gene_symbol": gene,
            "name": name,
            "sequence": cache.get(acc, ""),
            "pdb_or_af_ref": "",
            "surface_region": "",
            "layer_A": "True",
            "layer_B": "False",
            "layer_C": "False",
            "autoimmune_conditions": "",
            "mimicry_pathogens": "",
            "mimicry_epitope": "",
            "mimicry_epitope_type": "",
            "is_anchor_offtarget": "True" if is_anchor else "False",
            "source": "SURFY (Bausch-Fluck et al. PNAS 2018, table S3, 'surface' label)",
            "cspa_validated": "True" if acc in cspa else "False",
            "surfy_fpr_class": str(fpr),
        })

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SCHEMA)
        writer.writeheader()
        writer.writerows(rows)

    n_seq = sum(1 for r in rows if r["sequence"])
    n_cspa = sum(1 for r in rows if r["cspa_validated"] == "True")
    print(f"Wrote {len(rows)} rows, {n_seq} with sequences, {n_cspa} CSPA-validated.")
    for acc, (gene, nm) in ANCHORS.items():
        r = next((x for x in rows if x["uniprot_id"] == acc), None)
        status = "seq OK" if r and r["sequence"] else "NO SEQ"
        print(f"  anchor {nm} ({acc}): present, {status}")


if __name__ == "__main__":
    main()
