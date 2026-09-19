# EBV-MS Publication Figures & Visual Assets

This directory organizes all high-resolution figures, vector sources, supplementary plots, and exploratory graphical assets across the EBV–MS project.

---

## 1. Main Manuscript Figures (`manuscript/`)

| File | Description | Supporting Computational Source |
|---|---|---|
| `figure1.pdf` | **Study Workflow & Claim Boundary**: Overview of peptide curation, 6,400-pair screen, multi-step robustness filters, and claim boundaries. | `code/pipeline/make_publication_claim_ladder.R` |
| `figure2.pdf` | **Multi-Metric Sequence Sensitivity**: Lead-rank stability across TCR-facing identity, full-core BLOSUM62, and Grantham distance across all four alleles. | `code/pipeline/biochemical_similarity.py` |
| `figure3.pdf` | **Lead Robustness Dashboard**: Ablation grids, matched-decoy evaluation, whole-proteome null rarity (11.25M 9-mers), and two-predictor consensus. | `code/pipeline/lead_focused_robustness.py` |
| `figure4.pdf` | **Structural Ensemble Sensitivity**: AlphaFold 3 parent-versus-nested peptide C$\alpha$ RMSD distributions demonstrating register stability. | `code/pipeline/same_register_af3_analysis.py` |

---

## 2. Supplementary Figures (`supplement/`)

| File | Description |
|---|---|
| `figureS1.pdf` | Within-allele normalized BLOSUM62 score distributions for all 6,400 pairs across the four alleles. |
| `figureS2.pdf` | Multi-predictor binding percentile concordance for the initial eight audited candidate arms. |
| `figureS3.pdf` | Binding-threshold sensitivity across percentile cutoffs (1%, 2%, 5%, 10%, 20%). |

---

## 3. High-Resolution Panel Assets (`panels_v2/`)

Contains individual component plots, raw vector SVGs, and high-DPI rendered panels used in manuscript assembly.

---

## 4. Exploratory & Diagnostic Plots (`exploratory/`)

Contains preliminary diagnostic figures, score distribution tests, and early screening exploratory visualizations.
