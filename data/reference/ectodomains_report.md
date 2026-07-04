# Ectodomain extraction report (Step B)

Input proteins: **2896**  |  rows emitted: **2896**  |  non-empty ectodomain: **2889** (99.8%)

## Method breakdown

| method | count | % |
|---|---:|---:|
| `largest_extracellular_topodom` | 2292 | 79.1% |
| `pre_first_TM_heuristic` | 317 | 10.9% |
| `mature_chain_full` | 139 | 4.8% |
| `full_fallback` | 148 | 5.1% |

## Topology / flags

- Proteins with >=1 transmembrane span: **2749**
- Sequence mismatches (CSV != UniProt canonical, forced full fallback): **0**
- Full fallbacks (all causes: fetch fail, mismatch, too-short ecto, no topology): **148**
- Fetch failures (no usable UniProt JSON): **1**
- Truncated at 1000 aa cap: **154**

## Size distribution (residues; empties excluded)

| | min | median | max |
|---|---:|---:|---:|
| orig_len | 61 | 440 | 14507 |
| ecto_len | 30 | 206 | 1000 |

- Proteins originally >1000 aa: **370**
- Ectodomains <=400 aa (cheap to cofold): **2070**
- Ectodomains 400-1000 aa: **819**
- Empty ectodomain (missing input sequence): **7**

## Derivation rules

1. **>=1 `Extracellular` topological domain** -> largest such domain (`largest_extracellular_topodom`). Ties break by earliest start.
2. **TM spans, no `Extracellular` topo domain** -> largest contiguous non-TM segment from the mature N-terminus (`pre_first_TM_heuristic`); for type-I this is the pre-first-TM region.
3. **No TM span** -> mature chain (largest `Chain` feature, else signal-stripped full length) (`mature_chain_full`).
4. **Fallbacks** (`full_fallback`): fetch failure, CSV/canonical sequence mismatch, no topology, or derived ectodomain < 30 aa -> full CSV sequence.
5. All methods cap output at 1000 aa; `truncated=True` when the cap dropped a larger extracellular domain or hard-clipped a chain.
