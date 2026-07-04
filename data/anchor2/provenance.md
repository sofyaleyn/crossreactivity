# Anchor 2 provenance — ABT-736 (humanized 7C6) → PF4

Second, independent antibody→self-protein off-target anchor for the CrossFlag
cofold screen. All sequences are transcribed from cited public primary sources;
none are inferred or guessed.

## The anchor at a glance

| | |
|---|---|
| **Antibody** | ABT-736 = humanized 7C6 ("7C6hum7"), a humanized IgG1 anti-amyloid-beta (Aβ) oligomer/globulomer mAb (Abbott/AbbVie) developed for Alzheimer's disease |
| **Parent murine clone** | m7C6, hybridoma ATCC deposit **PTA-7240** (raised against Aβ(20-42) "globulomer") |
| **Off-target (self-protein)** | Platelet factor 4 (PF4 / CXCL4), **UniProt P02776** |
| **Off-target form implicated** | Mature secreted chemokine, PF4 tetramer (residues 32-101). PF4 is soluble/secreted, not membrane; no membrane domain applies. |
| **Mechanism** | CDR/paratope-mediated. High-affinity off-target binding to human & cynomolgus PF4 (not mouse/rat PF4). ABT-736·PF4 immune complexes crosslink platelet FcγR → platelet activation/aggregation, thrombocytopenia and thrombosis — a HIT-like pathology — causing acute infusion reactions and chronic toxicity in cynomolgus monkeys; program terminated. Structural modelling attributes PF4 binding to a prolonged, polar heavy-chain loop (H-CDR3, e.g. Ser/Tyr residues), the same loop feature seen in the anti-PF4/HIT Fabs KKO (PDB 4R9Y) and RTO (PDB 4RAU). A backup clone (h4D10) lacking this loop did not bind PF4 and had no toxicity. |
| **Confidence** | High for the off-target claim + mechanism (peer-reviewed, dedicated paper). High for sequence identity (patent figure, cross-checked against the patent's own CDR position claims). |

## Evidence for the off-target claim (mechanism + citation)

Primary paper (PubMed / PMC, open access via mAbs):

- Steinmetz KL, et al. **"Off-target binding of an anti-amyloid beta monoclonal
  antibody to platelet factor 4 causes acute and chronic toxicity in cynomolgus
  monkeys."** *mAbs* 2021;13(1):1887628.
  - PMID: **33596779** · PMCID: **PMC7894423**
  - DOI: https://doi.org/10.1080/19420862.2021.1887628
  - Key claims used: ABT-736 = humanized version of mouse m7C6; ABT-736
    immunoprecipitates an ~8 kDa plasma protein identified by mass-spec + anti-PF4
    western as PF-4; binds human & cyno PF4 but not mouse/rat PF4; toxicity is
    off-target (not Aβ-mediated); backup h4D10 lacks PF4 binding and is non-toxic;
    structural modelling implicates a prolonged polar H-chain loop, analogous to
    the loops in anti-PF4 Fabs KKO (PDB 4R9Y) and RTO (PDB 4RAU).

Supporting structural analogues cited in that paper (public PF4·Fab complexes):
- **PDB 4R9Y** — KKO Fab · PF4 (HIT antibody)
- **PDB 4RAU** — RTO Fab · PF4
- **PDB 1F9Q** — PF4 used as the docking target in the paper's model

## Antibody VH / VL sequences (files: vh.fasta, vl.fasta)

Source: **US Patent Application US2009/0175847 A1**, "Humanized antibodies to Aβ
(20-42) globulomer and uses thereof" (Abbott; inventors incl. Barghorn, Hillen,
Hinton, Juan; pub. 2009-07-09; family member **WO2008/150949 A1**).
- Full-text: https://patents.google.com/patent/US20090175847A1/en
- Sequences transcribed from **FIGURE 2** (figure sheet US20090175847A1-20090709-D00002):
  - **FIG. 2(B) = SEQ ID NO:3** = "7C6 VH (hum7)" (variable heavy, 118 aa) → vh.fasta
  - **FIG. 2(D) = SEQ ID NO:4** = "7C6 VL (hum7)" (variable light kappa, 113 aa) → vl.fasta
- The figure underlines the CDRs; transcription was cross-validated against the
  CDR residue ranges asserted in the granted family patent **US9540432B2**
  ("Anti-Aβ globulomer 7C6 antibodies"; https://patents.google.com/patent/US9540432B2/en),
  which defines 7C6 (VH = SEQ ID NO:11, VL = SEQ ID NO:12 in that patent's numbering)
  with: VH CDR1 = res 31-35, VH CDR2 = 50-65, VH CDR3 = 98-107; VL CDR1 = 24-39,
  VL CDR2 = 55-61, VL CDR3 = 94-102. Applying those ranges to SEQ ID NO:3/4 above
  reproduces exactly the CDRs underlined in FIG. 2 (H1 SYAMS, H2 SIHNRGTIFYLDSVKG,
  H3 GRSNSYAMDY; L1 RSTQTLVHRNGDTYLE, L2 KVSNRFS, L3 FQGSHVPYT), confirming a
  consistent read.

Notes:
- ABT-736 corresponds to the wild-type humanized graft "7C6hum7" (the patent also
  describes a "7C6hum7mut" variant; the wt VH/VL SEQ ID NO:3/4 are used here).
- Humanization was CDR-grafting, so the PF4-driving CDRs are identical to those of
  the murine parent m7C6 (ATCC PTA-7240) modelled against PF4 in the mAbs paper.

## Off-target sequence (file: offtarget_PF4_P02776.fasta)

Source: **UniProtKB P02776** (PLF4_HUMAN, gene PF4), https://rest.uniprot.org/uniprotkb/P02776.fasta
- Precursor is 101 aa; UniProt annotates SIGNAL 1..31 and CHAIN 32..101
  ("Platelet factor 4").
- offtarget_PF4_P02776.fasta contains the **mature chemokine, residues 32-101 (70 aa)** —
  the physiologically relevant, secreted, tetramer-forming form that antibodies engage.
  The signal peptide (1-31) is deliberately excluded for cofolding.
- Disulfides (UniProt): 41-67, 43-83. Heparin-binding basic residues ~92-98.

## Suggested cofold jobs
- ABT-736 Fv (vh.fasta + vl.fasta) vs PF4 (offtarget_PF4_P02776.fasta) — expect a
  confident H-CDR3-loop-mediated interface if the method generalizes.
- Optional negative control: pair the same PF4 against the backup h4D10 Fv (reported
  non-binder). h4D10 VH/VL are in the same patent family (SEQ ID NOs for 4D10) if a
  matched negative is wanted later.
