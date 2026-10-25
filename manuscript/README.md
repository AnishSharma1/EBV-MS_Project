# EBV-MS Manuscript & Supplementary Deliverables

This directory contains the submission-ready manuscript, supplement, LaTeX sources, display items, and audit records for the study:

**"Register-aware computation prioritizes EBV--self peptide pairs for experimental testing across HLA class II contexts"**
Author: Anish Sharma (Liberal Arts and Science Academy, Austin, TX)  
Target Journal: *BMC Bioinformatics*  
Public Code/Data: [https://github.com/AnishSharma1/EBV-MS_Project](https://github.com/AnishSharma1/EBV-MS_Project)

---

## Primary Deliverables

| File | Description | Status |
|---|---|---|
| [`main.tex`](main.tex) | Main manuscript LaTeX source | V20 (V19 plus post-audit corrections; see `V20_CHANGES.md`) |
| [`supplement.tex`](supplement.tex) | Independent supplement LaTeX source | V20 (V19 plus post-audit corrections; see `V20_CHANGES.md`) |
| [`references.tex`](references.tex) | Shared bibliography file | Audited |
| [`BMC_SUBMISSION_FIELDS.md`](BMC_SUBMISSION_FIELDS.md) | Formatted journal submission fields | Complete |

---

## Display Items

- **Main-text figures** (`figures/`): `main.tex` includes three files, which render as Figures 1-3.
  - `figure1.pdf` (Figure 1): workflow for computational prioritization; the 6,400-pair screen, the eight-pair evidence audit, and the 49-pair evaluation panel as separate tracks.
  - `figure4.pdf` (Figure 2): four-metric rank sensitivity of the two BALF5-TALDO1 pairs.
  - `figure5.pdf` (Figure 3): sensitivity of predicted peptide coordinates to input length (parent 15-mer versus nested 11-mer AlphaFold 3 models).
- **Files in `figures/` that `main.tex` does not include:** `figure2.pdf` (lead-pair TCR-facing BLOSUM62 ranks per allele), `figure3.pdf` (eight-pair evidence-audit grid), `figure6.pdf` (NetMHCIIpan-4.3 EL versus MixMHC2pred percentile concordance). They are retained for provenance; earlier versions of this README described a different four-figure set.
- **Supplementary figures** (`supplement_figures/`): `figureS1.pdf` (6,400-pair within-allele score distributions) and `figureS2.pdf` (initial eight-arm binding-predictor concordance). `supplement.tex` includes S1 and S2 only; no figure S3 exists or is referenced. The binding-threshold grid is in `output/revision_round_4/threshold_sensitivity_2026-09-20/`; a figure can be regenerated from it via `code/revision_analyses/run_computational_tightening.py` if one is added.
- **Supplementary tables:** there is no `supplement_tables/` directory. `supplement.tex` renders Tables S1-S8 inline, plus one unnumbered method-gating summary table:
  S1 top sequence-ranked pairs per allele; S2 initial eight-candidate evidence audit; S3 experimental-routing classes (49-pair cohort); S4 constrained-null results; S5 lead-pair sequences, frozen cores, and output percentiles; S6 lead rank and frozen-core sensitivity audit; S7 declared P1-P9 cores for the two BALF5-TALDO1 comparisons (LaTeX label `tab:s2a`, which renders as Table S7); S8 frozen-data ablations of ranking and routing.
- **Duplicate figure files:** `../figures/manuscript/` mirrors `figures/`. Figures 2-6 are byte-identical; `figure1.pdf` was synchronized to the compiled version on 2026-10-17 (the older copy used earlier box labels, "49 pair expanded cohort" and "provisional tiers"; it remains in Git history).

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

---

## Build Status

`paper_v20_main.pdf` / `paper_v20_supplement.pdf` (and their `.txt` extractions) are snapshots of the current sources; `paper_v20_SOURCE_SHA256.txt` records the SHA-256 of the sources, figures, and snapshots, and `code/tests/test_repository_integrity.py` fails if the sources drift from it (make the next change a V21 snapshot). The changes from V19 are listed in `V20_CHANGES.md`; `paper_v19_*` and earlier snapshots are unchanged.

XeLaTeX is not installed in the environment that produced V20, so the snapshots were built offline with `tectonic` 0.17.0 (`--only-cached`; XeTeX engine, `xdvipdfmx` driver). The main manuscript is 16 pages (matching V19) and the supplement 11 (V19: 10, because the 25-entry reference list now spills its last two entries onto an 11th page), with no unresolved citations or references. Fonts and layout were not compared page by page with V19. Rebuild with XeLaTeX (commands above, twice each) if the journal workflow requires it.

The hash files `PROVENANCE_SHA256.txt` and `FINAL_PDF_VS_LATEX_SHA256.txt` describe an earlier (v7) source and are historical.
