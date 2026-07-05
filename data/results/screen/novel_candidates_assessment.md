# CrossFlag Novel-Candidate Biological Assessment — SHR-1210 / Camrelizumab

**Scope:** Honest, literature-grounded read on the "valid-regime" novel off-target candidates surfaced by the CrossFlag Boltz-2 cofold screen of camrelizumab (SHR-1210, humanized IgG4κ anti-PD-1) against ~1,198 human cell-surface ectodomains. Candidates assessed: **SMO**, and the Ig-fold immune-receptor cluster **CD19, IL22RA1, KIR2DL4, KIR2DS1, KIR3DL3, KIR2DL3, KIR2DL5B**. Confirmed published off-targets used as anchors: **VEGFR2, FZD5, ULBP2** (Finlay et al., *mAbs* 2019; PMID 30541416; [DOI](https://doi.org/10.1080/19420862.2018.1550321)).

---

## BOTTOM LINE FOR THE DEMO

**None of the novel candidates has independent literature support, and every one of them falls into a fold class that already contains a confirmed off-target — which is exactly the signature of a structural-family (Ig-fold / Frizzled-CRD) affinity artifact, not new biology.** The single most important control is that the confirmed off-targets were found by an *experimental* human-receptor-proteome screen (Retrogenix cell microarray, ~4,975 receptors) that returned **only** VEGFR2, FZD5, and ULBP2 — SMO, CD19, IL22RA1 and all five KIRs were present in that near-complete surface-receptor panel and were **empirical non-hits** (PMID 30541416; [DOI](https://doi.org/10.1080/19420862.2018.1550321)). Against that ground truth, the CrossFlag novel calls are most parsimoniously read as false positives arising because camrelizumab's target PD-1 is itself an Ig-fold (IgV) antigen and all three real off-targets are Ig-superfamily/CRD folds — so the method concentrates confident interfaces in those same folds and cannot discriminate a genuine cross-reaction from fold-level β-sandwich compatibility. **This is an honest negative result: CrossFlag cleanly recovers the real off-targets but its "novel" hits are best explained as an Ig-fold/CRD affinity bias.** No candidate should be promoted as a forward prediction without wet-lab (e.g. Retrogenix/SPR) confirmation.

---

## Structural context: the anchors are already fold-biased

| Anchor (role) | Fold of ectodomain | Note |
|---|---|---|
| PD-1 (the *intended* target) | single **IgV** Ig-fold | The paratope was raised against an Ig-fold. |
| VEGFR2 (confirmed off-target) | **7 Ig-like domains** (Ig superfamily) | Ig-superfamily. |
| ULBP2 (confirmed off-target) | **MHC class-I-like** (α1/α2 platform + Ig-like) | β-sheet / Ig-adjacent NKG2D ligand. |
| FZD5 (confirmed off-target) | **Frizzled cysteine-rich domain (FZ-CRD)** | The one non-Ig anchor; a disulfide-stapled α-helical CRD. |

Two of three real off-targets are Ig-superfamily; the target itself is an Ig-fold. This matters for interpretation: **surfacing more Ig-fold antigens is what you expect from both a true signal and an artifact, so an Ig-fold candidate carries very little discriminating information.** The burden of proof is therefore high, and the null hypothesis (Ig-fold affinity bias) is the default.

---

## Candidate 1 — SMO (Smoothened), the single most confident novel hit

- **CrossFlag metrics:** PAE_IF **5.24 Å**, epitope_reprod **0.835**, antigen_len 202 aa — the *top-ranked* novel candidate, more confident than any KIR and comparable to the anchors.
- **Protein/fold:** Class-F GPCR; the screened 202-aa ectodomain is its **cysteine-rich domain (CRD)**, which is **directly homologous to the Frizzled FZ-CRD** — same 10-cysteine disulfide pattern and fold (Rana et al., *Nat Commun* 2013, PMID 24351982; Xu/Nusse family CRD homology, PMID 9096311; [web](https://www.nature.com/articles/ncomms3965)). SMO and Frizzled are sister Wnt/Hedgehog-pathway receptors.
- **Structural/epitope plausibility:** High *fold* overlap with the confirmed FZD5 hit — but that is precisely the problem. A paratope surface that docks the FZD5-CRD will geometrically fit the homologous SMO-CRD; SMO is **not an independent new epitope**, it is a within-family bleed-through from the FZD5 hit. This candidate was already flagged as a within-Frizzled-family false positive relative to FZD5 in earlier CrossFlag testing, which is consistent with this reading.
- **Independent evidence of anti-PD-1 / camrelizumab binding SMO:** **None found.** SMO was in the Retrogenix panel and was an empirical non-hit (PMID 30541416; [DOI](https://doi.org/10.1080/19420862.2018.1550321)). No polyspecificity report links any anti-PD-1 to Smoothened. Biologically, SMO surface exposure is low (Hedgehog signaling, primary-cilium-localized), an unlikely productive off-target in vivo.
- **Verdict:** **LIKELY CRD structural artifact (within-Frizzled-family cross-reaction with the FZD5 anchor).** Not an independent novel target. *Key cite:* PMID 24351982 (SMO≈FZD CRD homology); PMID 30541416 (empirical non-hit).

---

## Candidate cluster 2 — The Ig-fold NK/immune-receptor cluster: five KIRs + CD19 + IL22RA1

### The KIRs (KIR2DL4, KIR2DS1, KIR3DL3, KIR2DL3, KIR2DL5B)

- **CrossFlag metrics:** PAE_IF 7.2–8.7 Å, epitope_reprod 0.62–0.77 — a *coherent block of one gene family* clearing the hit rule together.
- **Protein/fold:** Killer-cell Ig-like receptors; ectodomains are **tandem C2-type Ig domains** (D0/D1/D2, V-shaped β-sandwich) (KIR2DL4 structure, PMID via [PMC4400354](https://pmc.ncbi.nlm.nih.gov/articles/PMC4400354/)). They are NK/T-cell receptors that bind **HLA class I** (MHC-I).
- **Structural/epitope plausibility:** This is the crux of the Ig-fold-bias case. KIRs recognize MHC-I; ULBP2 (a *confirmed* off-target) is **MHC-I-like** — so the cluster sits in the same β-sheet/MHC-I structural neighborhood as a real hit. But **the recurrence of an entire gene family (five paralogous KIRs at once) is the hallmark of a fold-level artifact, not target-specific CDR recognition.** Genuine paratope cross-reactivity tracks a *shared conformational epitope*, not a whole β-sandwich fold class. That five near-identical Ig folds all clear the threshold together argues the interface confidence is driven by generic Ig–Ig β-sheet packing (abundantly represented in the PDB that Boltz-2 was trained on), not by a specific SHR-1210 epitope.
- **Robustness flag (internal, supports artifact call):** Under the `anchor2` control condition, KIR2DL4, KIR2DL5B and KIR3DL3 collapse to PAE_IF 17–25 Å with epitope_reprod 0.17–0.44 — i.e. the KIR "hits" are **fragile / anchor-dependent** and do not reproduce when the docking anchor is perturbed. Real interfaces should be more robust.
- **Independent evidence:** **None.** All five KIRs were in the Retrogenix receptor proteome and were empirical non-hits (PMID 30541416). No literature reports any anti-PD-1 antibody binding a KIR (note: therapeutic anti-KIR antibodies such as lirilumab exist, but there is no reported PD-1-axis cross-reactivity).
- **Verdict:** **LIKELY Ig-fold artifact (whole-family recurrence + anchor-fragile + experimental non-hit).** *Key cite:* PMID 30541416 (empirical non-hits); KIR C2-Ig fold [PMC4400354].

### CD19

- **CrossFlag metrics:** PAE_IF 7.58 Å, epitope_reprod 0.76, 272 aa.
- **Protein/fold:** B-cell co-receptor; extracellular region is a **double Ig fold** — two C2-type Ig-like domains that form an elongated β-sandwich by C-terminal-half domain swapping ([PMC8470474](https://pmc.ncbi.nlm.nih.gov/articles/PMC8470474/)). Expression is **B-lineage-restricted**.
- **Structural/epitope plausibility:** Ig-fold, same β-sandwich class as PD-1 and VEGFR2 — again fold-level compatibility, no shared specific epitope with any anchor.
- **Independent evidence:** **None.** CD19 was in the Retrogenix panel and a non-hit (PMID 30541416). Camrelizumab binding CD19 would be conspicuous clinically (B-cell effects) and is not reported.
- **Verdict:** **LIKELY Ig-fold artifact.** *Key cite:* PMID 30541416; CD19 double-Ig fold [PMC8470474].

### IL22RA1

- **CrossFlag metrics:** PAE_IF **5.43 Å**, epitope_reprod 0.689, 213 aa — second-most-confident novel hit after SMO.
- **Protein/fold:** Class-II cytokine receptor α-chain; ectodomain is **two fibronectin type-III (FnIII) domains** at ~right angles ([R&D/atlas](https://atlasgeneticsoncology.org/gene/44568/il22ra1-(interleukin-22-receptor-alpha-1))). FnIII is a β-sandwich fold in the **immunoglobulin superfamily** — structurally kin to VEGFR2's Ig-like/FnIII-type extracellular architecture.
- **Structural/epitope plausibility:** FnIII β-sandwich sits in the broad Ig-superfamily basin the method favors; no specific shared epitope with an anchor. Its high confidence despite being an FnIII (not classic IgV/C2) reinforces that the bias is **β-sandwich-general**, not narrowly IgV.
- **Independent evidence:** **None.** IL22RA1 was in the Retrogenix panel and a non-hit (PMID 30541416). No anti-PD-1 polyspecificity report.
- **Verdict:** **LIKELY Ig-superfamily (FnIII β-sandwich) artifact.** *Key cite:* PMID 30541416; IL22RA1 FnIII fold (atlas/R&D).

---

## Summary verdict table

| Candidate | Ectodomain fold | Shares fold with anchor? | In Retrogenix panel? | Independent lit. evidence | Verdict |
|---|---|---|---|---|---|
| **SMO** | Frizzled CRD (FZ-CRD) | Yes — homologous to **FZD5** CRD | Yes → non-hit | None | **LIKELY CRD artifact** (within-Frizzled-family bleed from FZD5) |
| **KIR2DL4** | tandem C2 Ig | Ig / MHC-I neighborhood (ULBP2) | Yes → non-hit | None | **LIKELY Ig-fold artifact** (family recurrence; anchor-fragile) |
| **KIR2DS1** | tandem C2 Ig | " | Yes → non-hit | None | **LIKELY Ig-fold artifact** |
| **KIR3DL3** | tandem C2 Ig (D0-D1-D2) | " | Yes → non-hit | None | **LIKELY Ig-fold artifact** (anchor-fragile) |
| **KIR2DL3** | tandem C2 Ig | " | Yes → non-hit | None | **LIKELY Ig-fold artifact** |
| **KIR2DL5B** | tandem C2 Ig | " | Yes → non-hit | None | **LIKELY Ig-fold artifact** (anchor-fragile) |
| **CD19** | double (swapped) C2 Ig | Ig (PD-1/VEGFR2) | Yes → non-hit | None | **LIKELY Ig-fold artifact** |
| **IL22RA1** | 2× FnIII (Ig superfamily) | Ig-superfamily (VEGFR2) | Yes → non-hit | None | **LIKELY Ig-superfamily artifact** |

---

## Why this is still a useful demo result

1. **The recovery-vs-novel split is the story.** CrossFlag reproduces exactly the three experimentally confirmed off-targets and then, among "novel" calls, returns *only* proteins from the same fold families those anchors belong to (Frizzled-CRD for SMO; Ig/β-sandwich for the KIR–CD19–IL22RA1 cluster). That co-localization is the diagnostic signature of a **fold-affinity bias**, and being able to name it honestly is a real methodological finding.
2. **There is a hard external control.** The Finlay/Retrogenix experimental receptor-proteome screen is the ground truth, it covered these same proteins, and it did not flag any of them (PMID 30541416; [DOI](https://doi.org/10.1080/19420862.2018.1550321)). A computational "novel off-target" that an orthogonal wet screen scored as negative should not be reported as a discovery.
3. **Actionable mitigations to show the reviewer:** (a) add a fold-composition null / enrichment correction so Ig-fold and CRD antigens are scored against a fold-matched background rather than the whole proteome; (b) require **anchor-robustness** (candidates must survive the `anchor2` perturbation — the KIRs do not); (c) treat within-family homologs of a confirmed hit (SMO↔FZD5) as one event, not an independent prediction.

*No candidate is elevated to "plausible genuine cross-reactivity." If one had to be watch-listed, it would be **SMO** — solely because its anchor (FZD5) is a real off-target with a homologous CRD surface, making the physical cross-reaction mechanistically coherent even though it is not an independent new target.*

**PubMed attribution:** Key literature retrieved via PubMed. Primary source: Finlay WJJ et al., *mAbs* 2019;11(1):26–44, PMID 30541416, [DOI](https://doi.org/10.1080/19420862.2018.1550321). Supporting: Rana R et al. (SMO CRD) PMID 24351982, [DOI](https://doi.org/10.1038/ncomms3965).
