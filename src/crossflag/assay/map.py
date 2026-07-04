"""Self-protein -> confirmation assay map (HANDOFF Step 7). STUB.

Map a flagged self-protein / biology to the confirmation assay to run first
(receptor proteome array / SPR / functional assay). Table in docs/mvp-spec.md §[5].
"""


def assay_for(protein_id: str, biology: str | None = None) -> str:
    raise NotImplementedError("assay_for: see HANDOFF Step 7 / docs/mvp-spec.md §[5] assay table.")
