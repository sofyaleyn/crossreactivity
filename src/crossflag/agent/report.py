"""LLM agent narration (HANDOFF Step 8). STUB — not yet implemented.

The LLM reads the ranked evidence and writes: which variant to advance, what
each flagged variant likely cross-reacts with, which assay confirms, and the
honest caveat (embedding similarity is triage; the ladder is the evidence).
Retrieval + reasoning; no training. See docs/mvp-spec.md §[6].
"""


def write_report(ranked_panel: list[dict]) -> str:
    raise NotImplementedError("write_report: see HANDOFF Step 8 / docs/mvp-spec.md §[6].")
