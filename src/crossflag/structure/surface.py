"""Surface fingerprint (HANDOFF Step 10, rung 2). STUB — EXT.

Molecular surface (MSMS) -> CDR patch within ~4 A of CDR-H3/L3 -> APBS /
hydrophobicity annotation -> dMaSIF fingerprint + cosine vs. flagged
self-protein surface. CPU fallback: APBS/MSMS patch comparison, no dMaSIF.
See docs/extensions-spec.md §[8].
"""


def surface_fingerprint(structure_path, cdr_indices):
    raise NotImplementedError("surface_fingerprint: see HANDOFF Step 10 / docs/extensions-spec.md §[8].")
