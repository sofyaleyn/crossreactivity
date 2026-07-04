"""MVP gate (HANDOFF acceptance tests; thresholds in docs/demo-run.md Beat 1-2).

    run pipeline on {WT + germlined mutants}
    assert VEGFR2 in top-3 flagged proteins for WT
    assert risk_score(WT) > risk_score(at least one germlined mutant)

Skipped until crossflag.pipeline is implemented. This is the test that must
pass to freeze Phase 1.
"""
import pytest


@pytest.mark.skip(reason="pipeline not implemented yet (HANDOFF Step 9)")
def test_vegfr2_ranks_for_wt():
    from crossflag.pipeline import run_pipeline  # noqa: F401
    raise AssertionError("implement crossflag.pipeline, then assert VEGFR2 in top-3 for WT")
