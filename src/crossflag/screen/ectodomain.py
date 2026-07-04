"""Step B of the CrossFlag off-target screen: ectodomain extraction.

For each self-protein we cofold an antibody against the *extracellular*
(ectodomain) region rather than the full-length chain. Multi-pass /
single-pass membrane proteins have large cytoplasmic and TM regions that are
irrelevant to a surface-binding antibody and inflate cofold cost (some chains
are up to ~14.5k aa). This module fetches UniProt topology annotations and
derives, per protein, the extracellular sequence to fold.

Deterministic: given the same input CSV and the same cached UniProt JSON, the
output CSV is byte-identical. Reruns are cheap because raw JSON is cached under
``data/reference/raw/uniprot_topology/``.

Regenerate the CSV with::

    python -m crossflag.screen.ectodomain

Everything here is free: UniProt REST + local compute only.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import requests

from ..reference import paths

# ---------------------------------------------------------------------------
# Paths (read-only import of the shared path module; new dirs derived locally
# so we do not have to edit paths.py).
# ---------------------------------------------------------------------------
RAW_TOPOLOGY = paths.RAW / "uniprot_topology"
INPUT_CSV = paths.SELF_PROTEINS
OUTPUT_CSV = paths.REFERENCE / "ectodomains.csv"
REPORT_MD = paths.REFERENCE / "ectodomains_report.md"

OUTPUT_COLUMNS = [
    "protein_id",
    "uniprot_id",
    "gene_symbol",
    "orig_len",
    "ecto_len",
    "n_tm",
    "n_extracellular_domains",
    "has_signal",
    "method",
    "truncated",
    "seq_mismatch",
    "ectodomain_seq",
]

# Extraction parameters.
MAX_ECTO_LEN = 1000  # cap: cofold cost scales with antigen length
MIN_ECTO_LEN = 30    # below this the ectodomain is not usable -> full fallback

# Fetch tuning.
N_WORKERS = 16
MAX_RETRIES = 5
TIMEOUT = 30
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"

METHOD_TOPODOM = "largest_extracellular_topodom"
METHOD_PRE_TM = "pre_first_TM_heuristic"
METHOD_MATURE = "mature_chain_full"
METHOD_FALLBACK = "full_fallback"


# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------
def _cache_path(acc: str) -> Path:
    return RAW_TOPOLOGY / f"{acc}.json"


def fetch_topology(acc: str, session: requests.Session | None = None) -> dict | None:
    """Fetch and cache the UniProt entry JSON for ``acc``.

    Returns the parsed JSON dict, or ``None`` on unrecoverable failure (e.g.
    obsolete/demerged accession, or the network being down). Cached responses
    (including a cached empty object for known-missing accessions) short-circuit
    the network. Retries 429/5xx with exponential backoff.
    """
    cache = _cache_path(acc)
    if cache.exists():
        try:
            data = json.loads(cache.read_text())
        except json.JSONDecodeError:
            data = None
        if data:
            return data
        return None  # cached miss

    sess = session or requests
    url = UNIPROT_URL.format(acc=acc)
    for attempt in range(MAX_RETRIES):
        try:
            r = sess.get(url, timeout=TIMEOUT)
        except requests.RequestException:
            time.sleep(min(2 ** attempt, 30))
            continue
        if r.status_code == 200:
            data = r.json()
            cache.write_text(json.dumps(data))
            return data
        if r.status_code in (404, 400, 410):
            # Permanently missing / obsolete accession: cache an empty miss.
            cache.write_text("{}")
            return None
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(min(2 ** attempt, 30))
            continue
        # Other 4xx: treat as a miss, don't hammer.
        cache.write_text("{}")
        return None
    return None  # exhausted retries; do NOT cache -> retried on next run


def fetch_all(accessions: list[str]) -> dict[str, dict | None]:
    """Parallel-fetch topology JSON for all accessions (cached, retrying)."""
    from concurrent.futures import as_completed

    results: dict[str, dict | None] = {}
    with requests.Session() as session:
        with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
            futures = {ex.submit(fetch_topology, acc, session): acc for acc in accessions}
            done = 0
            for fut in as_completed(futures):
                acc = futures[fut]
                try:
                    results[acc] = fut.result()
                except Exception:  # noqa: BLE001 - defensive; treat as miss
                    results[acc] = None
                done += 1
                if done % 200 == 0:
                    print(f"  fetched {done}/{len(accessions)}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Feature parsing
# ---------------------------------------------------------------------------
@dataclass
class Topology:
    canonical_seq: str = ""
    tm_spans: list[tuple[int, int]] = field(default_factory=list)
    extracellular: list[tuple[int, int]] = field(default_factory=list)
    cytoplasmic: list[tuple[int, int]] = field(default_factory=list)
    signal_end: int | None = None
    chains: list[tuple[int, int]] = field(default_factory=list)


def _loc(feature: dict) -> tuple[int, int] | None:
    loc = feature.get("location", {})
    s = loc.get("start", {}).get("value")
    e = loc.get("end", {}).get("value")
    if s is None or e is None:
        return None
    return int(s), int(e)


def parse_topology(data: dict | None) -> Topology:
    """Extract the topology features relevant to ectodomain derivation."""
    topo = Topology()
    if not data:
        return topo
    topo.canonical_seq = data.get("sequence", {}).get("value", "") or ""
    for feat in data.get("features", []):
        ftype = feat.get("type")
        span = _loc(feat)
        if span is None:
            continue
        desc = (feat.get("description") or "").lower()
        if ftype == "Transmembrane":
            topo.tm_spans.append(span)
        elif ftype == "Topological domain":
            if "extracellular" in desc:
                topo.extracellular.append(span)
            elif "cytoplasmic" in desc:
                topo.cytoplasmic.append(span)
        elif ftype == "Signal":
            topo.signal_end = span[1]
        elif ftype == "Chain":
            topo.chains.append(span)
    return topo


# ---------------------------------------------------------------------------
# Ectodomain derivation
# ---------------------------------------------------------------------------
@dataclass
class Ecto:
    ecto_seq: str
    n_tm: int
    n_extracellular_domains: int
    has_signal: bool
    method: str
    truncated: bool
    seq_mismatch: bool


def _span_len(span: tuple[int, int]) -> int:
    return span[1] - span[0] + 1


def _slice(seq: str, span: tuple[int, int]) -> str:
    # UniProt positions are 1-based inclusive.
    return seq[span[0] - 1: span[1]]


def _pick_capped_domain(seq: str, domains: list[tuple[int, int]]) -> tuple[str, bool]:
    """Choose the extracellular domain to fold, honouring the length cap.

    Prefer the largest domain that is <= MAX_ECTO_LEN. If every domain exceeds
    the cap, take the single largest and truncate it to MAX_ECTO_LEN. Returns
    (sequence, truncated_flag). Ties in span length break deterministically by
    earliest start position.
    """
    # Deterministic ordering: by descending length, then ascending start.
    ordered = sorted(domains, key=lambda s: (-_span_len(s), s[0]))
    within = [d for d in ordered if _span_len(d) <= MAX_ECTO_LEN]
    if within:
        best = within[0]
        seq_out = _slice(seq, best)
        # truncated only if a *larger* domain existed that we could not use.
        truncated = _span_len(ordered[0]) > MAX_ECTO_LEN
        return seq_out, truncated
    # All domains too big: take largest and hard-cap.
    biggest = ordered[0]
    return _slice(seq, biggest)[:MAX_ECTO_LEN], True


def derive_ectodomain(csv_seq: str, data: dict | None) -> Ecto:
    """Derive the ectodomain sequence for one protein.

    ``csv_seq`` (the sequence already stored in self_proteins.csv) is the source
    of truth for residue indexing. If it disagrees with the UniProt canonical
    sequence, feature positions are not trustworthy, so we flag ``seq_mismatch``
    and fall back to the full sequence.
    """
    topo = parse_topology(data)
    n_tm = len(topo.tm_spans)
    n_ec = len(topo.extracellular)
    has_signal = topo.signal_end is not None
    csv_seq = csv_seq or ""

    def fallback(seq_mismatch: bool) -> Ecto:
        seq = csv_seq
        truncated = False
        if len(seq) > MAX_ECTO_LEN:
            seq = seq[:MAX_ECTO_LEN]
            truncated = True
        return Ecto(seq, n_tm, n_ec, has_signal, METHOD_FALLBACK, truncated, seq_mismatch)

    # No sequence to work with, or fetch failed / obsolete accession.
    if not csv_seq or not data or not topo.canonical_seq:
        return fallback(seq_mismatch=False)

    # Verify feature indexing frame.
    if topo.canonical_seq != csv_seq:
        return fallback(seq_mismatch=True)

    seq = csv_seq
    method: str
    ecto_seq: str
    truncated = False

    if topo.extracellular:
        ecto_seq, truncated = _pick_capped_domain(seq, topo.extracellular)
        method = METHOD_TOPODOM
    elif topo.tm_spans:
        # TM spans but no explicit extracellular topo domain.
        # Rule (documented): the ectodomain is the largest contiguous segment
        # that is not inside any TM span, measured from the mature N-terminus
        # (after the signal peptide). For a type-I single-pass protein this is
        # the N-terminal region before the first TM; the rule generalises to
        # multi-pass / type-II by simply taking the largest non-TM segment.
        mature_start = (topo.signal_end + 1) if has_signal else 1
        boundaries = sorted(topo.tm_spans)
        segments: list[tuple[int, int]] = []
        cursor = mature_start
        for tm_s, tm_e in boundaries:
            if tm_s - 1 >= cursor:
                segments.append((cursor, tm_s - 1))
            cursor = max(cursor, tm_e + 1)
        if cursor <= len(seq):
            segments.append((cursor, len(seq)))
        if segments:
            ecto_seq, truncated = _pick_capped_domain(seq, segments)
        else:
            ecto_seq = ""
        method = METHOD_PRE_TM
    else:
        # No TM: secreted / GPI-anchored / soluble. Use the mature chain,
        # which excludes the signal peptide and any C-terminal propeptide
        # (e.g. GPI-anchor signal). Fall back to signal-stripped full length.
        if topo.chains:
            main_chain = sorted(topo.chains, key=lambda s: (-_span_len(s), s[0]))[0]
            ecto_seq = _slice(seq, main_chain)
        elif has_signal:
            ecto_seq = seq[topo.signal_end:]
        else:
            ecto_seq = seq
        if len(ecto_seq) > MAX_ECTO_LEN:
            ecto_seq = ecto_seq[:MAX_ECTO_LEN]
            truncated = True
        method = METHOD_MATURE

    # Minimum-length guard: an implausibly short ectodomain -> full fallback.
    if len(ecto_seq) < MIN_ECTO_LEN:
        return fallback(seq_mismatch=False)

    return Ecto(ecto_seq, n_tm, n_ec, has_signal, method, truncated, seq_mismatch=False)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def _read_input() -> list[dict]:
    rows = []
    with INPUT_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def _write_report(out_rows: list[dict], n_input: int, n_fetch_fail: int) -> None:
    def col(name):
        return [r[name] for r in out_rows]

    methods = {}
    for r in out_rows:
        methods[r["method"]] = methods.get(r["method"], 0) + 1

    n_tm_rows = sum(1 for r in out_rows if int(r["n_tm"]) > 0)
    n_mismatch = sum(1 for r in out_rows if r["seq_mismatch"] == "True")
    n_truncated = sum(1 for r in out_rows if r["truncated"] == "True")
    n_fallback = methods.get(METHOD_FALLBACK, 0)

    orig = [int(r["orig_len"]) for r in out_rows if int(r["orig_len"]) > 0]
    ecto = [int(r["ecto_len"]) for r in out_rows if int(r["ecto_len"]) > 0]
    n_orig_gt1000 = sum(1 for v in orig if v > 1000)

    n_le400 = sum(1 for v in ecto if v <= 400)
    n_400_1000 = sum(1 for v in ecto if 400 < v <= 1000)
    n_empty = sum(1 for r in out_rows if int(r["ecto_len"]) == 0)

    def stats(vals):
        if not vals:
            return "n/a", "n/a", "n/a"
        return min(vals), int(median(vals)), max(vals)

    o_min, o_med, o_max = stats(orig)
    e_min, e_med, e_max = stats(ecto)

    lines = [
        "# Ectodomain extraction report (Step B)",
        "",
        f"Input proteins: **{n_input}**  |  rows emitted: **{len(out_rows)}**  |  "
        f"non-empty ectodomain: **{len(out_rows) - n_empty}** "
        f"({100 * (len(out_rows) - n_empty) / len(out_rows):.1f}%)",
        "",
        "## Method breakdown",
        "",
        "| method | count | % |",
        "|---|---:|---:|",
    ]
    for m in [METHOD_TOPODOM, METHOD_PRE_TM, METHOD_MATURE, METHOD_FALLBACK]:
        c = methods.get(m, 0)
        lines.append(f"| `{m}` | {c} | {100 * c / len(out_rows):.1f}% |")
    lines += [
        "",
        "## Topology / flags",
        "",
        f"- Proteins with >=1 transmembrane span: **{n_tm_rows}**",
        f"- Sequence mismatches (CSV != UniProt canonical, forced full fallback): **{n_mismatch}**",
        f"- Full fallbacks (all causes: fetch fail, mismatch, too-short ecto, no topology): **{n_fallback}**",
        f"- Fetch failures (no usable UniProt JSON): **{n_fetch_fail}**",
        f"- Truncated at {MAX_ECTO_LEN} aa cap: **{n_truncated}**",
        "",
        "## Size distribution (residues; empties excluded)",
        "",
        "| | min | median | max |",
        "|---|---:|---:|---:|",
        f"| orig_len | {o_min} | {o_med} | {o_max} |",
        f"| ecto_len | {e_min} | {e_med} | {e_max} |",
        "",
        f"- Proteins originally >1000 aa: **{n_orig_gt1000}**",
        f"- Ectodomains <=400 aa (cheap to cofold): **{n_le400}**",
        f"- Ectodomains 400-1000 aa: **{n_400_1000}**",
        f"- Empty ectodomain (missing input sequence): **{n_empty}**",
        "",
        "## Derivation rules",
        "",
        "1. **>=1 `Extracellular` topological domain** -> largest such domain "
        "(`largest_extracellular_topodom`). Ties break by earliest start.",
        "2. **TM spans, no `Extracellular` topo domain** -> largest contiguous "
        "non-TM segment from the mature N-terminus (`pre_first_TM_heuristic`); "
        "for type-I this is the pre-first-TM region.",
        "3. **No TM span** -> mature chain (largest `Chain` feature, else "
        "signal-stripped full length) (`mature_chain_full`).",
        "4. **Fallbacks** (`full_fallback`): fetch failure, CSV/canonical "
        "sequence mismatch, no topology, or derived ectodomain < "
        f"{MIN_ECTO_LEN} aa -> full CSV sequence.",
        f"5. All methods cap output at {MAX_ECTO_LEN} aa; `truncated=True` when "
        "the cap dropped a larger extracellular domain or hard-clipped a chain.",
        "",
    ]
    REPORT_MD.write_text("\n".join(lines))


def main() -> None:
    RAW_TOPOLOGY.mkdir(parents=True, exist_ok=True)
    rows = _read_input()
    print(f"Loaded {len(rows)} proteins from {INPUT_CSV}", file=sys.stderr)

    accessions = [r["uniprot_id"] for r in rows if r.get("uniprot_id")]
    # De-dup preserving order for the fetch (CSV has unique uniprot_ids, but be safe).
    seen: set[str] = set()
    to_fetch = [a for a in accessions if not (a in seen or seen.add(a))]
    print(f"Fetching topology for {len(to_fetch)} accessions ({N_WORKERS} workers)...",
          file=sys.stderr)
    topo_json = fetch_all(to_fetch)

    n_fetch_fail = sum(1 for a in to_fetch if not topo_json.get(a))
    print(f"Fetch complete. {n_fetch_fail} accessions without usable JSON.", file=sys.stderr)

    out_rows: list[dict] = []
    for r in rows:
        acc = r.get("uniprot_id", "")
        csv_seq = (r.get("sequence") or "").strip()
        data = topo_json.get(acc)
        res = derive_ectodomain(csv_seq, data)
        out_rows.append({
            "protein_id": r.get("protein_id", ""),
            "uniprot_id": acc,
            "gene_symbol": r.get("gene_symbol", ""),
            "orig_len": len(csv_seq),
            "ecto_len": len(res.ecto_seq),
            "n_tm": res.n_tm,
            "n_extracellular_domains": res.n_extracellular_domains,
            "has_signal": str(res.has_signal),
            "method": res.method,
            "truncated": str(res.truncated),
            "seq_mismatch": str(res.seq_mismatch),
            "ectodomain_seq": res.ecto_seq,
        })

    with OUTPUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {OUTPUT_CSV} ({len(out_rows)} rows)", file=sys.stderr)

    _write_report(out_rows, n_input=len(rows), n_fetch_fail=n_fetch_fail)
    print(f"Wrote {REPORT_MD}", file=sys.stderr)


if __name__ == "__main__":
    main()
