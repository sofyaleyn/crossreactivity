"""Phase 2 gate (HANDOFF acceptance tests; thresholds in docs/demo-run.md Beat 3).

    cofold(WT, VEGFR2) interface_confidence > cofold(germlined_mutant, VEGFR2) interface_confidence

Skipped until crossflag.structure.cofold is implemented.
"""
import pytest


@pytest.mark.skip(reason="cofold not implemented yet (HANDOFF Step 11)")
def test_wt_interface_beats_germlined():
    from crossflag.structure.cofold import cofold  # noqa: F401
    raise AssertionError("implement crossflag.structure.cofold, then assert WT interface > mutant")
