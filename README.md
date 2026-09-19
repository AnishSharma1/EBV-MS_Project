# Computational Assessment of EBV–Self Peptide Similarities Across HLA Class II Contexts

[![Verify Computational Package](https://github.com/AnishSharma1/EBV-MS_Project/actions/workflows/verify_package.yml/badge.svg)](https://github.com/AnishSharma1/EBV-MS_Project/actions/workflows/verify_package.yml)
[![License: MIT and CC BY 4.0](https://img.shields.io/badge/License-MIT%20%26%20CC--BY--4.0-blue.svg)](LICENSE)

**Author:** Anish Sharma  
**Affiliation:** Liberal Arts and Science Academy, Austin, TX, USA  
**Correspondence:** `anishsharma0212@gmail.com`  
**Target Journal:** *BMC Bioinformatics*  
**Canonical Package:** [`computational_package/`](computational_package/)  
**Manuscript & Supplement:** [`manuscript/`](manuscript/)

---

## 1. Study Overview

Epstein–Barr virus (EBV) infection and the HLA class II risk haplotype (predominantly `HLA-DRB1*15:01`) are convergent risk factors in Multiple Sclerosis (MS) pathogenesis. Molecular mimicry between viral and myelin or central nervous system (CNS) self-peptides has long been proposed as a candidate trigger, but computational mimicry screens are often confounded by arbitrary scoring thresholds, register misassignment, or uncalibrated docking predictions.

This study implements an auditable, control-calibrated computational prioritization workflow evaluating candidate molecular mimicry across four primary HLA class II risk contexts:
- `HLA-DRB1*15:01`
- `HLA-DRB1*13:03`
- `HLA-DRB1*08:01`
- `HLA-DRB1*03:01`

We perform a systematic within-allele discovery screen of **6,400 allele-specific pairs** (40 EBV and 40 human CNS peptides per allele), prioritizing candidates by normalized BLOSUM62 similarity at T-cell receptor (TCR)-facing positions (P2, P3, P5, P7, P8). Top candidates are evaluated through multi-metric sequence sensitivity, 10,000-replicate constrained-null simulations, matched empirical decoys, proteome-wide sequence rarity filters (11.25 million 9-mers), multi-predictor binding consensus (NetMHCIIpan 4.3 and MixMHC2pred 2.1), and AlphaFold 3 structural ensemble stability across parent and nested peptide registers.

---

## 2. Study Evolution: Two-Leg Architecture

This project progressed through two distinct methodological phases:

### Leg 1: Exploratory Prototyping Phase
- **Focus**: Early heuristic screening, single-chain AlphaFold 2 and ColabFold modeling, exploratory machine learning feature extraction, and preliminary TCR docking.
- **Artifacts**: Preserved for historical provenance and transparency in [`notebooks/`](notebooks/) (`AlphaFold2.ipynb`, `Colabfold.ipynb`, `Cross Reactivity.ipynb`, etc.) and [`archive/`](archive/).
- **Critical Methodological Lessons**:
  1. Uncalibrated docking scores and generic sequence matching generated excessive false-positive candidate lists.
  2. Lack of register resolution obscured whether predicted matching residues faced the TCR or were buried in MHC anchor pockets.
  3. The classic Hy.2E11 BALF5–MBP positive-control system was clarified to be an allele-switching system (`BALF5` presented by `DRB5*01:01` [DR2a]; `MBP` presented by `DRB1*15:01` [DR2b]) rather than a single-allele DR15 mimicry event.

### Leg 2: Current Publication Architecture (Rebuilt & Calibrated)
- **Focus**: The definitive, reproducible computational pipeline supporting the submitted manuscript.
- **Artifacts**: Implemented in [`computational_package/`](computational_package/) and [`manuscript/`](manuscript/).
- **Upgrades**:
  1. Expansion to a four-allele class II universe with strict within-allele ranking (6,400 pairs).
  2. Score normalization restricted strictly to solvent-exposed, TCR-facing registers (P2, P3, P5, P7, P8).
  3. Multi-predictor binding consensus (NetMHCIIpan 4.3 + MixMHC2pred 2.1) and proteome-wide null background screens (11.25M 9-mers).
  4. AlphaFold 3 parent-versus-nested structural ensemble sensitivity verifying core register invariance.
  5. Negative control gates: APBS electrostatics and uncalibrated TCRModel2 docking were benchmarked, failed validation gates, and were strictly excluded from candidate scoring.

---

## 3. Key Prioritized Candidates

| Candidate ID | Allele | Viral Peptide | Human Self Peptide | Screen Rank | Routing Tier | Evidence & Status |
|---|---|---|---|---|---|---|
| **HP47** | `HLA-DRB1*15:01` | BALF5 627–641 (`TGGVYHFVKKHVHES`) | TALDO1 216–230 (`SVTKIYNYYKKFSYK`) | 14 / 1,600 | `strong_medium_provisional` | Multi-metric sequence robustness; high multi-predictor binding; invariant AF3 parent/nested register. Verified viral DR15 binding in IEDB (Assay 1774674). Recommended for direct recombinant binding & register mapping. |
| **HP36** | `HLA-DRB1*13:03` | BALF5 627–641 (`TGGVYHFVKKHVHES`) | TALDO1 108–122 (`DARLSFDKDAMVARA`) | 13 / 1,600 | `weak_medium_provisional` | Multi-metric sequence robustness; computational binding consensus passed. Direct binding measurement pending. |

> [!IMPORTANT]
> **Positive Control Calibration**:  
> The well-known Hy.2E11 BALF5–MBP cross-reactive system operates across distinct DR15 haplotype molecules (`BALF5` on `DRB5*01:01` [PDB 1H15]; `MBP` on `DRB1*15:01` [PDB 1BX2]) rather than a single allele. Our study preserves this biological distinction and uses it as a benchmark anchor.

---

## 4. Strict Scientific Claim Boundaries

To prevent over-interpretation of computational findings, this study enforces explicit boundaries:
- **Supported**: Auditable candidate ranking, allele-specific register evaluation, structural model quality control, and experimental routing prioritization.
- **NOT Claimed**:
  - No proof of natural *in vivo* antigen processing or presentation.
  - No measured binding affinity (computational percentiles are not $K_D$ values).
  - No proof of TCR binding, shared recognition, or T-cell activation.
  - No claim of causal molecular mimicry in MS patients.
  - APBS continuum electrostatics and uncalibrated TCRModel2 ternary docking failed control gates and were strictly excluded from candidate scoring.

---

## 5. Standardized 5-Directory Repository Architecture

```
EBV-MS_Project/
├── plans/                            # Study designs, protocols & statistical specifications
│   ├── study_scope.md                # Research question & presentation boundaries
│   ├── statistical_plan.md           # Discovery screen & constrained-null design
│   ├── publication_strategy.md       # Target journal & thesis outline
│   └── README.md                     # Directory guide
│
├── code/                             # Analysis code, test suites & staging tools
│   ├── pipeline/                     # Modular analysis scripts
│   ├── tests/                        # Repository integrity tests
│   ├── staging/                      # Staging & packaging utilities
│   └── README.md                     # Execution guide
│
├── outputs/                          # Derived tables, candidate rankings & benchmark manifests
│   ├── ebv_ms_benchmark_manifest/    # Benchmark reference manifests
│   ├── ebv_ms_model_package/         # Processed model metric outputs
│   └── README.md                     # Output catalog & evidence ledger
│
├── figures/                          # Publication figures & visual assets
│   ├── manuscript/                   # Figures 1–4 (vector PDF & SVG)
│   ├── supplement/                   # Supplementary figures S1–S3
│   └── README.md                     # Figure captions & computational sources
│
├── logs/                             # Verification logs, audit trails & QA reports
│   ├── provenance_audits/            # Sentence-level text provenance logs
│   ├── citation_audits/              # Reference verification audits
│   ├── quality_assurance/            # Language & formatting review logs
│   ├── package_verification/         # Cryptographic SHA-256 manifests
│   └── README.md                     # Audit index
│
├── manuscript/                       # Submission-ready publication deliverables
│   ├── paper_v9_working_main.pdf     # Main manuscript PDF
│   ├── paper_v9_working_supplement.pdf # Independent supplement PDF
│   ├── main.tex & supplement.tex     # Complete LaTeX sources
│   └── BMC_SUBMISSION_FIELDS.md      # Journal submission metadata
│
├── computational_package/            # Canonical reproducible package (Leg 2)
│   ├── MANIFEST.sha256               # Checksum manifest
│   ├── scripts/                      # Verified pipeline scripts
│   ├── data/                         # Curated input libraries
│   └── results/                      # Frozen score tables & decision gates
│
└── notebooks/                        # Exploratory phase Jupyter notebooks (Leg 1)
    ├── README.md                     # Leg 1 context & provenance notes
    └── *.ipynb                       # Phase 1 interactive notebooks
```

---

## 6. Reproducibility & Validation

### Automated Package Verification
The public computational package contains a deterministic integrity checker verifying all required files and SHA-256 checksums:

```bash
python3 computational_package/scripts/verify_package.py
```

Expected output:
```text
PACKAGE VALIDATION PASSED
```

### Reproducing the Screen
To execute candidate screening or inspect frozen result tables:
- **Sequence similarity & screening:** `computational_package/scripts/rank_full_pmhc_screen.py`
- **Multi-metric robustness:** `computational_package/scripts/biochemical_similarity.py`
- **Predictor consensus:** `computational_package/scripts/add_mixmhc2pred_validation.py`
- **Results summaries:** `computational_package/results/stage1_evidence_review_2026-09-03/`

---

## 7. Citation

If you use this computational pipeline, candidate scoring framework, or curated datasets, please cite:

```bibtex
@article{sharma2026ebvms,
  author = {Sharma, Anish},
  title = {Computational assessment of EBV-self peptide similarities across HLA class II contexts},
  journal = {BMC Bioinformatics},
  year = {2026},
  note = {In submission. Code available at https://github.com/AnishSharma1/EBV-MS_Project}
}
```
Or use the metadata in [`CITATION.cff`](CITATION.cff).

---

## 8. License

- **Code & Scripts**: Licensed under the [MIT License](LICENSE).
- **Manuscript Text, Figures, & Data Tables**: Licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
