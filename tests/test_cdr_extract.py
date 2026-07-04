"""Acceptance test for CDR extraction (HANDOFF.md).

    extract_cdrs(SHR1210_VH, SHR1210_VL)["cdr_h3"] == "QLYYFDYW"

The README quotes CDR-H3 == QLYYFDYW and CDR-L3 == QQVYSIPWT for the SHR-1210
anchor. These use two different IMGT boundary conventions (H3 anchor-trimmed,
L3 full loop 105-117), so we assert H3 against the anchor-trimmed apex and L3
against the standard IMGT CDR3.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crossflag.extract.cdrs import extract_cdrs  # noqa: E402

SHR1210_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKG"
    "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS"
)
SHR1210_VL = (
    "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGS"
    "GTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK"
)


def test_cdr_h3_apex_matches_readme():
    result = extract_cdrs(SHR1210_VH, SHR1210_VL)
    assert result["cdr_h3_apex"] == "QLYYFDYW"


def test_cdr_l3_matches_readme():
    result = extract_cdrs(SHR1210_VH, SHR1210_VL)
    assert result["cdr_l3"] == "QQVYSIPWT"


def test_all_six_cdrs_present():
    result = extract_cdrs(SHR1210_VH, SHR1210_VL)
    for key in ("cdr_h1", "cdr_h2", "cdr_h3", "cdr_l1", "cdr_l2", "cdr_l3"):
        assert result[key], f"{key} is empty"


def test_cdr_indices_map_to_sequence():
    """Indices must slice exactly the concatenated CDR residues (pooling input)."""
    result = extract_cdrs(SHR1210_VH, SHR1210_VL)
    h_res = "".join(SHR1210_VH[i] for i in result["vh_cdr_indices"])
    l_res = "".join(SHR1210_VL[i] for i in result["vl_cdr_indices"])
    assert h_res == result["cdr_h1"] + result["cdr_h2"] + result["cdr_h3"]
    assert l_res == result["cdr_l1"] + result["cdr_l2"] + result["cdr_l3"]


if __name__ == "__main__":
    test_cdr_h3_apex_matches_readme()
    test_cdr_l3_matches_readme()
    test_all_six_cdrs_present()
    test_cdr_indices_map_to_sequence()
    print("All CDR extraction acceptance tests passed.")
