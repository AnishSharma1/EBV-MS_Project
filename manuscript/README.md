# EBV-MS Manuscript & Supplementary Deliverables

This directory contains the submission-ready manuscript, supplement, LaTeX sources, display items, and audit records for the study:

**"Computational assessment of EBV-self peptide similarities across HLA class II contexts"**  
Author: Anish Sharma (Liberal Arts and Science Academy, Austin, TX)  
Target Journal: *BMC Bioinformatics*  
Public Code/Data: [https://github.com/AnishSharma1/EBV-MS_Project](https://github.com/AnishSharma1/EBV-MS_Project)

---

## Primary Deliverables

| File | Description | Status |
|---|---|---|
| [`paper_v9_working_main.pdf`](paper_v9_working_main.pdf) | Compiled main manuscript PDF | Submission ready |
| [`paper_v9_working_supplement.pdf`](paper_v9_working_supplement.pdf) | Compiled independent supplementary material PDF | Submission ready |
| [`main.tex`](main.tex) | Main manuscript LaTeX source | V9 Working |
| [`supplement.tex`](supplement.tex) | Independent supplement LaTeX source | V9 Working |
| [`references.tex`](references.tex) | Shared bibliography file | Audited |
| [`BMC_SUBMISSION_FIELDS.md`](BMC_SUBMISSION_FIELDS.md) | Formatted journal submission fields | Complete |

---

## Display Items

- **Main Figures** (`figures/`):
  - `figure1.pdf`: Overall study workflow and claim boundary ladder
  - `figure2.pdf`: Multi-metric sequence sensitivity across HLA class II contexts
  - `figure3.pdf`: Lead robustness dashboard (ablation, decoy matching, proteome-null, binding filters)
  - `figure4.pdf`: Corrected structural ensemble sensitivity (parent vs. nested peptide C$\alpha$ RMSD)
- **Supplementary Figures** (`supplement_figures/`):
  - `figureS1.pdf`: 6,400-pair within-allele score distributions
  - `figureS2.pdf`: Initial eight-arm binding-predictor concordance
  - `figureS3.pdf`: Binding-threshold sensitivity across percentile cutoffs
- **Supplementary Tables** (`supplement_tables/`):
  - Tables S1–S12 providing complete tabular documentation (candidate provence, 49-pair cohort, ablation grids, register scenarios, matched decoys, and benchmark status).

---

## Provenance & Audit Records

- `content_migration_ledger.csv`: Row-by-row mapping of sections, tables, and figures from V8.2 into the clean V9 two-document split.
- `source_snapshot_sha256.json`: Checksum snapshot verifying input immutability.
- `TEXT_PROVENANCE_AUDIT.md`: Direct sentence-by-sentence verification against author-drafted prose.
- `CITATION_AUDIT.md`: Verification of all reference keys against PubMed/DOI records.
- `GRAMMARLY_REPAIR_LOG.md`: Stylistic and grammatical refinement audit.

---

## Compilation

To compile the documents from source using XeLaTeX:

```bash
# Compile main manuscript
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex

# Compile supplement
xelatex -interaction=nonstopmode supplement.tex
xelatex -interaction=nonstopmode supplement.tex
```
