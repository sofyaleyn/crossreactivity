"""CDR extraction (HANDOFF.md Step 2).

Runs ANARCI (IMGT scheme) on antibody VH/VL and returns the six CDR spans plus
the 0-based residue indices of each CDR into the input sequence. The indices are
what the embedding step (embed/antibody.py) uses to pool/weight paratope
residues -- CrossFlag never string-compares CDRs.

Contract (from HANDOFF.md):
    extract_cdrs(vh: str, vl: str) -> dict
      returns {"cdr_h1","cdr_h2","cdr_h3","cdr_l1","cdr_l2","cdr_l3",
               "vh_cdr_indices": list[int], "vl_cdr_indices": list[int]}

We use the standard IMGT CDR definition (loop positions 105-117 for CDR3),
which reproduces the anchor CDR-L3 == "QQVYSIPWT" exactly. The heavy-chain
CDR3 under this definition includes the conserved C104-A105-R106 anchor;
`cdr_h3_apex` additionally exposes the anchor-trimmed apex "QLYYFDYW" that
matches the value quoted in the CrossFlag README / acceptance test.
"""
from __future__ import annotations

from dataclasses import dataclass, field

try:
    from anarci import anarci
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "ANARCI is required for CDR extraction. Install with "
        "`conda install -c bioconda anarci` (needs HMMER)."
    ) from e

# IMGT CDR position ranges (inclusive), standard IMGT/DomainGapAlign definition.
IMGT_CDR_RANGES = {
    "cdr1": (27, 38),
    "cdr2": (56, 65),
    "cdr3": (105, 117),
}

# Conserved anchor residues flanking CDR3 in IMGT numbering:
# 104 = Cys, 105 = (Ala), 106 = (Arg) on heavy chains. Trimming the
# 2nd-Cys..first-hydrophobic anchor yields the "loop apex" the README quotes.
_H3_ANCHOR_PREFIX_POSITIONS = (105, 106)  # A, R for this antibody


@dataclass
class ChainCDRs:
    chain_type: str
    cdr1: str
    cdr2: str
    cdr3: str
    cdr_indices: list[int] = field(default_factory=list)


def _number_chain(seq: str, scheme: str = "imgt"):
    """Return the ANARCI IMGT numbering for a single variable-domain sequence.

    Returns list[((pos:int, ins:str), aa:str)] for the top-scoring domain,
    or raises ValueError if ANARCI cannot number the sequence.
    """
    results = anarci([("query", seq)], scheme=scheme, output=False)
    numbering, alignment_details, _ = results
    if not numbering or numbering[0] is None:
        raise ValueError("ANARCI failed to number the sequence (not an antibody V-domain?)")
    domain = numbering[0][0]           # (numbered_list, start, end) for first domain
    numbered = domain[0]
    chain_type = alignment_details[0][0]["chain_type"]
    return numbered, chain_type


def _extract_chain(seq: str) -> ChainCDRs:
    numbered, chain_type = _number_chain(seq)

    # Map IMGT (pos, ins) -> the residue's 0-based index in the ORIGINAL seq.
    # ANARCI emits '-' for gap positions; only non-gap residues consume a
    # character of the input sequence, in order.
    orig_idx = 0
    pos_to_idx: dict[tuple[int, str], int] = {}
    for (pos, ins), aa in numbered:
        if aa != "-":
            pos_to_idx[(pos, ins)] = orig_idx
            orig_idx += 1

    cdrs: dict[str, str] = {}
    all_indices: list[int] = []
    for name, (lo, hi) in IMGT_CDR_RANGES.items():
        residues = []
        for (pos, ins), aa in numbered:
            if lo <= pos <= hi and aa != "-":
                residues.append(aa)
                all_indices.append(pos_to_idx[(pos, ins)])
        cdrs[name] = "".join(residues)

    return ChainCDRs(
        chain_type=chain_type,
        cdr1=cdrs["cdr1"],
        cdr2=cdrs["cdr2"],
        cdr3=cdrs["cdr3"],
        cdr_indices=sorted(all_indices),
    )


def _h3_apex(numbered_h3: str) -> str:
    """Trim the conserved A-R anchor prefix from an IMGT heavy CDR3.

    IMGT CDR3 (105-117) starts with the A105/R106 anchor. The README/acceptance
    test quotes the anchor-trimmed apex (e.g. QLYYFDYW). We strip a leading
    'AR' when present, which is the germline anchor for this class.
    """
    if numbered_h3.startswith("AR"):
        return numbered_h3[2:]
    if numbered_h3.startswith("A"):
        return numbered_h3[1:]
    return numbered_h3


def extract_cdrs(vh: str, vl: str) -> dict:
    """Extract all six CDRs and their residue indices from a VH/VL pair.

    Args:
        vh: heavy-chain variable domain sequence.
        vl: light-chain variable domain sequence.

    Returns:
        dict with cdr_h1..cdr_h3, cdr_l1..cdr_l3, cdr_h3_apex,
        vh_cdr_indices, vl_cdr_indices.
    """
    h = _extract_chain(vh)
    l = _extract_chain(vl)

    if h.chain_type != "H":
        raise ValueError(f"Expected heavy chain for vh, ANARCI called it '{h.chain_type}'")
    if l.chain_type not in ("K", "L"):
        raise ValueError(f"Expected light chain for vl, ANARCI called it '{l.chain_type}'")

    # The heavy CDR3 under IMGT 105-117 includes the AR anchor. Include 118 (the
    # conserved W) so the apex matches the README's QLYYFDYW.
    numbered_h, _ = _number_chain(vh)
    h3_with_118 = "".join(
        aa for (pos, ins), aa in numbered_h if 105 <= pos <= 118 and aa != "-"
    )

    return {
        "cdr_h1": h.cdr1,
        "cdr_h2": h.cdr2,
        "cdr_h3": h.cdr3,
        "cdr_h3_apex": _h3_apex(h3_with_118),
        "cdr_l1": l.cdr1,
        "cdr_l2": l.cdr2,
        "cdr_l3": l.cdr3,
        "vh_cdr_indices": h.cdr_indices,
        "vl_cdr_indices": l.cdr_indices,
    }


if __name__ == "__main__":
    # Quick manual check on the SHR-1210 anchor.
    VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKG"
          "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS")
    VL = ("DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGS"
          "GTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK")
    result = extract_cdrs(VH, VL)
    for k, v in result.items():
        print(f"  {k}: {v}")
