# EBV-MS Codebase & Analysis Pipelines

This directory consolidates all computational code, analysis pipelines, test suites, and staging tools for the EBV–MS project.

---

## 1. Directory Structure

```
code/
├── pipeline/             # Active computational analysis and processing scripts (mirrors src/)
│   ├── prepare_datasets.py               # Rebuilds normalized tables from raw/ inputs
│   ├── rank_full_pmhc_screen.py          # Executes 6,400-pair within-allele discovery screen
│   ├── biochemical_similarity.py         # Multi-metric sequence scoring and sensitivity
│   ├── add_mixmhc2pred_validation.py     # Multi-predictor consensus (MixMHC2pred integration)
│   ├── run_null_model.py                 # 10,000-replicate constrained-null simulations
│   ├── run_focused_taldo1_register_analysis.py # TALDO1 register stability assessment
│   ├── same_register_af3_analysis.py     # AlphaFold 3 parent-vs-nested structural evaluation
│   └── ...                               # Domain-specific pipeline modules
│
├── tests/                # Automated validation suite (mirrors tests/)
│   ├── test_repository_integrity.py      # Assertions on project structure and file sizes
│   ├── test_verify_package.py            # Checksum and package manifest assertions
│   └── test_*.py                         # Numerical and unit tests
│
├── staging/              # Release packaging and staging tools
│   └── stage_github_release.py           # Validates repository size and builds GitHub releases
│
└── archive/              # Preserved early exploratory batch scripts
    ├── analyze_af_server_results.mjs
    ├── build_af_server_batch.mjs
    └── ...
```

---

## 2. Standard Execution Workflows

### Dataset Preparation:
To reprocess raw IEDB and UniProt downloads into clean, normalized row-level tables:
```bash
python3 code/pipeline/prepare_datasets.py
```

### Candidate Prioritization Screen:
To execute the multi-allele discovery screen across the 6,400 peptide pairs:
```bash
python3 code/pipeline/rank_full_pmhc_screen.py
```

### Running Validation Tests:
To run the complete repository integrity and regression test suite:
```bash
python3 -m pytest code/tests/
```

### Package Checksum Verification:
```bash
python3 computational_package/scripts/verify_package.py
```
