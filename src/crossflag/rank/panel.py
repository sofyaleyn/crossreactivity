"""Panel ranking (HANDOFF Step 7). STUB — not yet implemented.

Rank all variants by `risk_score`; attach top flagged protein + assay
(assay/map.py). Emit a panel table + per-variant evidence. See docs/mvp-spec.md §[5].
"""


def rank_panel(variant_scores: list[dict]) -> list[dict]:
    raise NotImplementedError("rank_panel: see HANDOFF Step 7 / docs/mvp-spec.md §[5].")
