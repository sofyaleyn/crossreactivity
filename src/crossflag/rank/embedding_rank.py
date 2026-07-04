"""Embedding-similarity ranking + calibration (HANDOFF Step 6, rung 1). STUB.

Contract:
    score_variant(paratope_vec, ref_index, background_index) -> dict
      returns {"risk_score": float, "neighbors": [(protein_id, cosine), ...]}

Cosine to every reference vector; aggregate (max or top-k mean); calibrate
against the benign background distribution. Higher = more suspect. Layer
priority weights per docs/reference-set.md §Scoring impact. See docs/mvp-spec.md §[4].
"""


def score_variant(paratope_vec, ref_index, background_index) -> dict:
    raise NotImplementedError("score_variant: see HANDOFF Step 6 / docs/mvp-spec.md §[4].")
