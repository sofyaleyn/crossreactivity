---
name: igfold_surface
description: Fold a variant (IgFold) and fingerprint its CDR surface patch vs. a flagged self-protein. STUB (extension).
---

# igfold_surface (STUB — extension)

Wraps `crossflag.structure.fold` + `crossflag.structure.surface`. See HANDOFF Step 10 / docs/extensions-spec.md §[8].

- **Input:** `{vh: str, vl: str, self_protein_ref}`
- **Output:** structural corroboration score per top-flagged pair.
- **Note:** dMaSIF is GPU + CC-BY-NC-ND; CPU fallback is APBS/MSMS patch comparison.

_Fill in exact schema + invocation when implemented (HANDOFF Step 12)._
