"""Antibody paratope embedding (HANDOFF Step 3). STUB — not yet implemented.

Contract:
    embed_antibody(vh: str, vl: str, cdr_indices: dict) -> np.ndarray  # paratope vector

Run AntiBERTy/AbLang2 -> per-residue embeddings; pool CDR residues (mean or
attention-weighted, weight CDR-H3 highest) into one vector. Deterministic;
cache by sequence hash. See docs/mvp-spec.md §[2].
"""


def embed_antibody(vh: str, vl: str, cdr_indices: dict):
    raise NotImplementedError("embed_antibody: see HANDOFF Step 3 / docs/mvp-spec.md §[2].")
