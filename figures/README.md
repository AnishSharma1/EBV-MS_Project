# EBV-MS Publication Figures

Vector PDFs for the manuscript and supplement. `manuscript/figures/` and `manuscript/supplement_figures/` are what `main.tex` and `supplement.tex` compile against; the files here mirror them (verified byte-identical on 2026-10-17).

---

## 1. Main-text figures (`manuscript/`)

`main.tex` includes `figure1.pdf`, `figure4.pdf`, and `figure5.pdf`, which render as Figures 1-3. The other three files are not included in the current manuscript.

| File | Rendered as | Content |
|---|---|---|
| `figure1.pdf` | Figure 1 | Workflow for computational prioritization: the 6,400-pair screen, the eight-pair evidence audit, and the 49-pair evaluation panel as separate tracks. |
| `figure4.pdf` | Figure 2 | Four-metric rank sensitivity of the two BALF5-TALDO1 pairs (TCR-facing BLOSUM62, full-core BLOSUM62, TCR-facing identity, Grantham distance). |
| `figure5.pdf` | Figure 3 | Sensitivity of predicted peptide coordinates to input length (AlphaFold 3 parent 15-mer versus nested 11-mer; 5x5 model comparisons per arm). |
| `figure2.pdf` | not included | Lead-pair TCR-facing BLOSUM62 ranks (13 and 14 of 1,600) across the four alleles. |
| `figure3.pdf` | not included | Eight-pair evidence-audit grid: binding consensus, register agreement, evidence arms, stage-1 status. |
| `figure6.pdf` | not included | NetMHCIIpan-4.3 EL versus MixMHC2pred percentile concordance for the audited arms. |

## 2. Supplementary figures (`supplement/`)

| File | Content |
|---|---|
| `figureS1.pdf` | Within-allele normalized BLOSUM62 score distributions for all 6,400 pairs. |
| `figureS2.pdf` | Binding-predictor percentile concordance for the initial eight audited candidate arms. |

There is no `figureS3.pdf`. The binding-threshold sensitivity grid exists as data in `output/revision_round_4/threshold_sensitivity_2026-09-20/`; no figure of it is included in the supplement.

## 3. Duplicate-copy note

On 2026-10-17 `figures/manuscript/figure1.pdf` was replaced with the compiled copy at `manuscript/figures/figure1.pdf`. The earlier file differed only in box wording ("49 pair expanded cohort", "provisional tiers") and remains in Git history.

## 4. Regeneration status

Figure-rendering scripts are in `code/figures/` (`build_workflow_and_rank_sensitivity.R/.py`, `build_structural_and_binding.R`, `restyle_figures_ggplot.R`, `build_panels.py`, `restyle_panels.py`, `reproduce_figures.R`). The scripts do not record which output file each figure corresponds to, and they cannot currently be run end to end from this checkout:

- several write to, or read from, absolute paths on the original author workstation;
- `reproduce_figures.R` expects `outputs/manuscript_figures_only_2026-09-06/source_data`, which is not present here;
- some input tables live in the separate `ebv_ms_publication/` workspace.

Treat the PDFs in this directory as the frozen published artifacts, not as outputs a reader can regenerate from the repository alone.
