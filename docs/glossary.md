# Glossary

Abbreviations and terms across CrossFlag docs.

---

## Antibody biology

| Term | Meaning |
|---|---|
| **mAb** | Monoclonal antibody. |
| **IgG / IgG4** | Immunoglobulin G / its IgG4 subclass. SHR-1210 is IgG4κ. |
| **Fv** | Variable fragment — the VH + VL pair that binds antigen. |
| **VH / VL** | Variable heavy / light chain. |
| **CDR** | Complementarity-determining region — hypervariable antibody loops that contact antigen. Three per chain: H1/H2/H3, L1/L2/L3. |
| **CDR-H3** | Third heavy-chain CDR; longest, most variable, dominates contact. |
| **paratope** | The antibody surface (CDRs) that contacts antigen. |
| **epitope** | The antigen surface recognized by the antibody. |
| **polyspecificity** | One CDR binding multiple distinct proteins — the CDR-mediated off-target problem CrossFlag targets. |
| **polyreactivity** | Non-specific binding from surface charge/hydrophobicity (distinct from polyspecificity). |
| **germlining** | Reverting CDR/framework residues to germline — reduces polyspecificity; the Finlay fix for SHR-1210. |
| **humanization** | Engineering a non-human antibody to reduce human immunogenicity. |
| **ADC** | Antibody-drug conjugate — antibody + cytotoxic payload; CDR cross-reactivity is especially dangerous. |
| **bispecific / CAR-T** | Two-antigen antibody / chimeric-antigen-receptor T-cell. |

## The anchor case

| Term | Meaning |
|---|---|
| **SHR-1210 / camrelizumab** | Anti-PD-1 IgG4κ antibody causing capillary hemangioma via CDR cross-reactivity with VEGFR2. |
| **PD-1** | Programmed cell death protein 1 — SHR-1210's intended target. |
| **VEGFR2** | Vascular endothelial growth factor receptor 2 — primary off-target; its agonism drives hemangioma. |
| **FZD5 / ULBP2** | Secondary SHR-1210 off-targets (frizzled-5; UL16 binding protein 2). |
| **Mab005** | Murine progenitor of SHR-1210; also VEGFR2-reactive → confirms CDRs are the cause. |
| **Finlay et al. 2019** | Paper identifying the off-targets and the CDR-germlining fix; source of the variant panel. |
| **ABT-736** | Anti-β-amyloid antibody discontinued in NHP for PF4 off-target; candidate second anchor. |
| **PF4** | Platelet factor 4 (ABT-736 off-target). |

## Language models & representation

| Term | Meaning |
|---|---|
| **PLM** | Protein language model — transformer trained on protein sequences; produces embeddings. |
| **embedding** | A learned numeric vector representing a residue or sequence; similar biology → nearby vectors. |
| **AntiBERTy** | Antibody-specific BERT PLM (558M antibody sequences); per-residue embeddings + attention. |
| **AbLang / AbLang2** | Antibody-specific PLMs (Oxford OPIG); AbLang2 is paired-chain, germline-bias-corrected. |
| **ESM-2** | Evolutionary Scale Modeling v2 — general-protein PLM; used for antigen embeddings. |
| **IgBERT / BALM / IgT5** | Other antibody PLMs (alternatives). |
| **rescoding / seqcoding** | Per-residue / per-sequence embedding outputs (AbLang terminology). |
| **cosine similarity** | Similarity measure between two embedding vectors. |

## Structure & cofolding

| Term | Meaning |
|---|---|
| **IgFold** | Antibody-specific structure predictor (VH+VL → 3D). |
| **AbodyBuilder2 / 3** | Alternative antibody structure predictors. |
| **cofolding / co-folding** | Predicting the joint 3D structure of two molecules together (e.g. antibody + antigen) to assess whether they bind. |
| **Boltz-2** | Open (MIT) co-folding + binding-affinity foundation model; the confirmation rung. |
| **AlphaFold3 (AF3)** | DeepMind's co-folding model; gated weights; slightly better on Ab–Ag accuracy. |
| **ipTM / pTM** | Interface / global predicted TM-score — model confidence in the (interface) structure. |
| **PAE** | Predicted aligned error — per-residue-pair confidence; low PAE at CDR–antigen contacts = confident interface. |
| **affinity score** | Boltz-2's predicted binding strength (`affinity_pred_value`, `affinity_probability_binary`). |
| **dMaSIF** | Differentiable Molecular Surface Interaction Fingerprinting — surface fingerprint for similarity. |
| **APBS / MSMS** | Poisson-Boltzmann electrostatics solver / molecular surface mesher (rung-2 fallback). |
| **Kyte-Doolittle** | Amino-acid hydrophobicity scale. |
| **pLDDT** | Per-residue structure-confidence score (IgFold/AF). |
| **PDB / SAbDab** | Protein Data Bank / Structural Antibody Database. |
| **MSA** | Multiple sequence alignment — input Boltz-2 can use (`--use_msa_server`). |

## Numbering & screening context

| Term | Meaning |
|---|---|
| **ANARCI / ANARCII** | Antibody/TCR numbering tool (ANARCII = 2025 LM version). |
| **IMGT / Kabat** | Antibody numbering schemes; CrossFlag uses IMGT. |
| **MPA** | Membrane Proteome Array — gold-standard cell-based off-target screen. |
| **SPR** | Surface plasmon resonance — binding-affinity biosensor assay. |
| **HUVEC** | Human umbilical vein endothelial cell — VEGFR2 functional assay. |
| **IEDB** | Immune Epitope Database (used as a naive linear-epitope baseline only). |

## Project

| Term | Meaning |
|---|---|
| **MVP / EXT** | Phase 1 (embedding rank) / Phase 2–3 (confirmation ladder). |
| **evidence ladder** | The three rungs: embedding rank → surface fingerprint → cofolding. |
| **agent skill** | A programmatically invokable tool wrapper; excludes web-only tools. |
| **CC-BY-NC-ND / MIT** | dMaSIF license (non-commercial) / Boltz-2, ESM-2 license (commercial-ok). |
| **ARDitox / Tope-seq** | Recent whole-proteome one-candidate tools CrossFlag differentiates from. |
