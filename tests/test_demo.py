"""Acceptance tests for the CrossFlag demo (demo-plan.md Phase 5 gate).

Asserts, from committed data alone:
  (a) the Phase-0 recompute of cofold-fzd5 from committed CIF/PAE reproduces the
      CSV within +/-0.05 PAE_IF and +/-0.03 epitope_reprod;
  (b) the anchor rows match the headline numbers in findings.md;
  (c) the verdict table reads FZD5/ULBP2=confirmed, VEGFR2=missed, lysozyme=floor;
  (d) within-fold FZD5 AUROC is 0.909 (rank 2/12);
  (e) the built dashboard exists and embeds its data inline (no external requests).
"""
from __future__ import annotations

import json
import re

import pytest

from crossflag.demo import build, panel, paths, scoring


@pytest.fixture(scope="module")
def built():
    """Regenerate all artifacts once, then assert the whole thing exits 0."""
    assert build.main() == 0
    return True


# --- (a) Phase-0 recompute self-consistency ---
def test_recompute_fzd5_matches_csv():
    rows = panel.load_rows()
    csv_row = rows["cofold-fzd5"]
    s = scoring.score_run("cofold-fzd5")
    assert abs(s.PAE_IF - csv_row.PAE_IF) <= 0.05, (s.PAE_IF, csv_row.PAE_IF)
    assert abs(s.epitope_reprod - csv_row.epitope_reprod) <= 0.03, (
        s.epitope_reprod, csv_row.epitope_reprod)


# --- (b) anchor rows match findings.md headline numbers ---
@pytest.mark.parametrize("run,pae,rep", [
    ("cofold-pd1", 7.24, 0.936),
    ("cofold-fzd5", 5.69, 0.660),
    ("cofold-ulbp2", 5.74, 0.893),
    ("cofold-lyz", 12.38, 0.446),
    ("cofold-wt", 11.55, 0.430),
])
def test_anchor_rows(run, pae, rep):
    r = panel.load_rows()[run]
    assert r.PAE_IF == pytest.approx(pae, abs=0.01)
    assert r.epitope_reprod == pytest.approx(rep, abs=0.01)


# --- (c) verdict table reads ---
def test_verdict_table_json(built):
    table = json.loads(paths.VERDICT_JSON.read_text())
    by_antigen = {t["antigen"]: t["verdict"] for t in table}
    assert by_antigen["FZD5"] == "confirmed"
    assert by_antigen["ULBP2"] == "confirmed"
    assert by_antigen["VEGFR2"] == "missed"
    assert by_antigen["lysozyme"] == "floor"
    assert by_antigen["PD-1"] == "ceiling"


def test_verdict_logic_direct():
    rows = panel.load_rows()
    assert panel.verdict("cofold-fzd5", rows["cofold-fzd5"]) == "confirmed"
    assert panel.verdict("cofold-ulbp2", rows["cofold-ulbp2"]) == "confirmed"
    assert panel.verdict("cofold-wt", rows["cofold-wt"]) == "missed"
    assert panel.verdict("cofold-lyz", rows["cofold-lyz"]) == "floor"


# --- (d) within-fold discrimination ---
def test_family_auroc():
    auroc, rank, n = panel.family_auroc(panel.load_rows())
    assert auroc == pytest.approx(0.909, abs=0.001)
    assert rank == 2 and n == 12


# --- (e) dashboard exists and embeds data inline ---
def test_dashboard_self_contained(built):
    html = paths.DASHBOARD.read_text()
    assert paths.DASHBOARD.exists() and len(html) > 10_000
    # data embedded inline
    assert 'id="crossflag-data"' in html
    m = re.search(r'id="crossflag-data">(.*?)</script>', html, re.S)
    data = json.loads(m.group(1))
    assert data["family"]["auroc"] == 0.909
    assert len(data["verdict_table"]) == 5
    # all images inline; no external resource requests
    # 5 anchor figures + 3 real-scale-screen figures (when screen_metrics.csv is present)
    assert html.count("data:image/png;base64") in (5, 8)
    external = [
        u for u in re.findall(r'(?:src|href)="([^"]+)"', html)
        if not u.startswith("data:")
    ]
    assert external == [], f"external resource requests found: {external}"


def test_dashboard_numbers_match_csv(built):
    """No number in the dashboard contradicts the CSV verdict table."""
    html = paths.DASHBOARD.read_text()
    for t in json.loads(paths.VERDICT_JSON.read_text()):
        assert f'{t["PAE_IF"]:.2f}' in html
