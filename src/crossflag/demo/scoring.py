"""Recompute cofold interface metrics from committed structures.

Adapted from ``data/results/scripts/extract_metrics.py`` / ``analyze_all.py``,
repointed at the committed ``data/results/structures/<run>/sample_N.cif`` +
``sample_N_pae.npz`` (no scratchpad tarballs, no network).

Two metrics, matching ``data/results/README.md``:

- ``PAE_IF`` — mean predicted aligned error over ALL antibody(H+L) x antigen(V)
  token pairs, symmetrized, averaged across the 5 samples. Lower = tighter.
- ``epitope_reprod`` — mean pairwise Jaccard of contacted antigen residues
  (5 A heavy-atom) across the 5 samples. Higher = same epitope every time.

These are the numbers that feed ``cofold_metrics.csv``; Phase 0 asserts this
recompute reproduces the committed CSV within tolerance.
"""
from __future__ import annotations

import glob
import itertools
import json
import warnings
from dataclasses import dataclass

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch

from . import paths

warnings.filterwarnings("ignore")

CONTACT = 5.0  # heavy-atom contact radius, Angstrom


def _block(pae: np.ndarray, a: list[int], b: list[int]) -> float:
    """Symmetrized mean PAE over the a x b token block."""
    a_i, b_i = np.array(a), np.array(b)
    return float((pae[np.ix_(a_i, b_i)].mean() + pae[np.ix_(b_i, a_i)].mean()) / 2)


def _chain_roles(model, antigen_len: int):
    """Split chains into (antibody, antigen) by residue count == antigen length."""
    counts = {c.id: len([r for r in c if r.id[0] == " "]) for c in model}
    antigen = [c for c, n in counts.items() if n == antigen_len]
    if not antigen:
        antigen = [min(counts, key=lambda c: abs(counts[c] - antigen_len))]
    antibody = [c for c in counts if c not in antigen]
    return antibody, antigen


def _epitope(cif: str, antigen_len: int) -> set:
    """Antigen residues contacted (<=5 A heavy-atom) by any antibody atom."""
    model = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    ab, ag = _chain_roles(model, antigen_len)
    ab_atoms = [a for c in model if c.id in ab for a in c.get_atoms() if a.element != "H"]
    ns = NeighborSearch(ab_atoms)
    ep = set()
    for c in model:
        if c.id not in ag:
            continue
        for r in c:
            if r.id[0] != " ":
                continue
            if any(ns.search(a.coord, CONTACT) for a in r if a.element != "H"):
                ep.add((c.id, r.id[1]))
    return ep


def _jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else 0.0


@dataclass
class Score:
    run: str
    PAE_IF: float
    epitope_reprod: float


def sequences_for(run: str) -> dict[str, str]:
    """H/L/V chain sequences from the committed input JSON for a run."""
    d = json.loads(paths.input_json_for(run).read_text())
    return {e["chain_ids"][0]: e.get("value", "") for e in d["entities"]}


def interface_pae_block(run: str) -> tuple[np.ndarray, int, int, int]:
    """Antibody(H+L) x antigen(V) PAE submatrix for the tightest sample.

    Returns (block, n_H, n_L, n_V). Rows = antibody tokens (H then L),
    cols = antigen tokens. The 'tightest' sample is the one with the lowest
    mean interface PAE (this is what PAE_IF summarizes).
    """
    seqs = sequences_for(run)
    H, L, V = seqs["H"], seqs["L"], seqs["V"]
    off_v = len(H) + len(L)
    ab_tok = list(range(0, off_v))
    v_tok = list(range(off_v, off_v + len(V)))

    sdir = paths.structure_dir_for(run)
    pae_files = sorted(glob.glob(str(sdir / "sample_*_pae.npz")))
    best_block, best_mean = None, float("inf")
    for f in pae_files:
        pae = np.load(f)["pae"]
        blk = pae[np.ix_(ab_tok, v_tok)]
        m = float(blk.mean())
        if m < best_mean:
            best_mean, best_block = m, blk
    return best_block, len(H), len(L), len(V)


def epitope_profile(run: str) -> tuple[list[int], np.ndarray]:
    """Per-antigen-residue contact count across the 5 samples.

    Returns (residue_numbers, counts) where counts[i] in 0..5 is how many of the
    samples contact antigen residue residue_numbers[i]. This is exactly what
    epitope_reprod summarizes: a sharp reproduced patch vs a scattered one.
    """
    seqs = sequences_for(run)
    V = seqs["V"]
    sdir = paths.structure_dir_for(run)
    cif_files = sorted(glob.glob(str(sdir / "sample_*.cif")))
    epis = [_epitope(c, len(V)) for c in cif_files]
    all_pairs = sorted({p for ep in epis for p in ep}, key=lambda p: p[1])
    if not all_pairs:
        return [], np.zeros(0, dtype=int)
    resnums = [rn for (_chain, rn) in all_pairs]
    counts = np.array([sum(p in ep for ep in epis) for p in all_pairs], dtype=int)
    return resnums, counts


def score_run(run: str) -> Score:
    """Recompute PAE_IF + epitope_reprod for one run from committed data."""
    seqs = sequences_for(run)
    H, L, V = seqs["H"], seqs["L"], seqs["V"]
    off_v = len(H) + len(L)
    ab_tok = list(range(0, off_v))
    v_tok = list(range(off_v, off_v + len(V)))

    sdir = paths.structure_dir_for(run)
    pae_files = sorted(glob.glob(str(sdir / "sample_*_pae.npz")))
    cif_files = sorted(glob.glob(str(sdir / "sample_*.cif")))
    if not pae_files or not cif_files:
        raise FileNotFoundError(f"No committed structures for run {run!r} in {sdir}")

    pae_vals = [_block(np.load(f)["pae"], ab_tok, v_tok) for f in pae_files]
    epis = [_epitope(c, len(V)) for c in cif_files]
    jvals = [_jaccard(a, b) for a, b in itertools.combinations(epis, 2)]

    return Score(
        run=run,
        PAE_IF=float(np.mean(pae_vals)),
        epitope_reprod=float(np.mean(jvals)) if jvals else float("nan"),
    )
