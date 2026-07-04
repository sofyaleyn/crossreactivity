"""Antigen embedding (HANDOFF Step 4). STUB — not yet implemented.

Contract:
    embed_antigen(sequence: str, surface_region: list[int] | None) -> np.ndarray | list[np.ndarray]

ESM-2. If `surface_region` given, pool those residues; else sliding-window
sub-vectors over the chain. Cache. See docs/mvp-spec.md §[3].
"""


def embed_antigen(sequence: str, surface_region=None):
    raise NotImplementedError("embed_antigen: see HANDOFF Step 4 / docs/mvp-spec.md §[3].")
