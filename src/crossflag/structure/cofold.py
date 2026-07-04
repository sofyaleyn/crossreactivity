"""Boltz-2 cofolding gate (HANDOFF Step 11, rung 3). STUB — EXT.

For the top ~5 variant x self-protein pairs, write a Boltz YAML (VH+VL +
antigen, affinity request), run `boltz predict --use_msa_server`, parse
CDR-restricted ipTM/PAE + `affinity_pred_value`. Return a confirmation verdict
per pair. See docs/extensions-spec.md §[9].
"""


def cofold(vh: str, vl: str, antigen_seq: str) -> dict:
    raise NotImplementedError("cofold: see HANDOFF Step 11 / docs/extensions-spec.md §[9].")
