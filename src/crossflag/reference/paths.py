"""Single source of truth for reference-set data locations.

All paths are anchored to the repository root (found by walking up from this
file until a directory containing ``pyproject.toml`` is seen), NOT to the
location of any individual build script. This lets the reference build modules
live in ``src/crossflag/reference/`` while reading/writing ``data/reference/``,
and keeps every location declared in one place.

Layout (see docs/reference-set.md):
  data/reference/
  ├── raw/          downloaded sources (gitignored, regenerable)
  │   ├── surfaceome/   SURFY + CSPA xlsx, seq_cache.tsv
  │   └── fasta/        per-accession UniProt fetches
  ├── seeds/        hand-curated inputs, committed (provenance)
  ├── build/        intermediates (gitignored, deterministic)
  ├── background/   benign_proteins.csv (committed final output)
  ├── index/        ESM-2 embedding cache (gitignored)
  └── self_proteins.csv   final merged reference set (committed)
"""
from pathlib import Path


def _find_repo_root(start: Path) -> Path:
    """Walk up from ``start`` until a directory with pyproject.toml is found."""
    for parent in [start, *start.parents]:
        if (parent / "pyproject.toml").exists():
            return parent
    # Fallback: repo layout is src/crossflag/reference/paths.py -> 3 levels up.
    return start.parents[3]


REPO_ROOT = _find_repo_root(Path(__file__).resolve())

DATA = REPO_ROOT / "data"
ANCHOR = DATA / "anchor"

REFERENCE = DATA / "reference"
RAW = REFERENCE / "raw"
RAW_SURFACEOME = RAW / "surfaceome"
RAW_FASTA = RAW / "fasta"
SEEDS = REFERENCE / "seeds"
BUILD = REFERENCE / "build"
BACKGROUND = REFERENCE / "background"
INDEX = REFERENCE / "index"

# Raw sources (downloaded)
SURFY_XLSX = RAW_SURFACEOME / "surfy_table_S3.xlsx"
CSPA_XLSX = RAW_SURFACEOME / "cspa_S2_File.xlsx"
SEQ_CACHE = RAW_SURFACEOME / "seq_cache.tsv"

# Hand-curated seeds (committed)
LAYER_C_SEED = SEEDS / "layer_c_mimicry_seed.csv"

# Intermediates (gitignored)
LAYER_A_CSV = BUILD / "self_proteins_layer_a.csv"
LAYER_C_FILLED = BUILD / "layer_c_mimicry_seed_filled.csv"
MISSING_ACCS = BUILD / "missing_accs.txt"

# Final outputs (committed)
SELF_PROTEINS = REFERENCE / "self_proteins.csv"
BENIGN_PROTEINS = BACKGROUND / "benign_proteins.csv"


def ensure_dirs() -> None:
    """Create the gitignored working directories if absent."""
    for d in (RAW_SURFACEOME, RAW_FASTA, BUILD, BACKGROUND, INDEX):
        d.mkdir(parents=True, exist_ok=True)
