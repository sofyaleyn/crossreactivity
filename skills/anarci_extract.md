---
name: anarci_extract
description: Extract the six CDR spans (IMGT) + CDR residue indices from antibody VH/VL. STUB.
---

# anarci_extract (STUB)

Wraps `crossflag.extract.cdrs.extract_cdrs`. See HANDOFF Step 2 / docs/mvp-spec.md §[1].

- **Input:** `{vh: str, vl: str}`
- **Invocation:** `python -c "from crossflag.extract.cdrs import extract_cdrs; ..."`
- **Output:** `{cdr_h1..cdr_l3: str, vh_cdr_indices: list[int], vl_cdr_indices: list[int]}`
- **Verify:** CDR-H3 apex == `QLYYFDYW`, CDR-L3 == `QQVYSIPWT` on the SHR-1210 anchor.

_Fill in exact input/output schema and invocation when the skill wiring lands (HANDOFF Step 12)._
