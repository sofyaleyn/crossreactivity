---
name: embedding_rank
description: Rank a variant against the reference set in embedding space, calibrated vs. benign background. STUB.
---

# embedding_rank (STUB)

Wraps `crossflag.rank.embedding_rank.score_variant`. See HANDOFF Step 6 / docs/mvp-spec.md §[4].

- **Input:** `{paratope_vec, ref_index, background_index}`
- **Output:** `{risk_score: float, neighbors: [(protein_id, cosine), ...]}`.

_Fill in exact schema + invocation when implemented (HANDOFF Step 12)._
