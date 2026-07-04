---
name: boltz_cofold
description: Cofold a variant Fv against a flagged self-protein (Boltz-2) and read out interface + affinity. STUB (extension).
---

# boltz_cofold (STUB — extension)

Wraps `crossflag.structure.cofold.cofold`. See HANDOFF Step 11 / docs/extensions-spec.md §[9].

- **Input:** `{vh: str, vl: str, antigen_seq: str}`
- **Invocation:** `boltz predict <yaml> --use_msa_server` (MIT; GPU).
- **Output:** `{interface_confidence (CDR-restricted ipTM/PAE), affinity_pred_value, verdict}`.

_Fill in exact schema + invocation when implemented (HANDOFF Step 12)._
