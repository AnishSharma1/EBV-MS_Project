# EBV-MS Complete Python Codebase & Analysis Pipelines

This directory houses the comprehensive collection of **170+ Python analysis scripts, pipelines, test suites, and staging utilities** powering the EBV–MS project.

---

## 1. Directory Structure

```
code/
├── pipeline/             # 100+ Primary analysis and data processing scripts
│   ├── prepare_datasets.py               # Rebuilds normalized tables from raw/ inputs
│   ├── rank_full_pmhc_screen.py          # Executes 6,400-pair within-allele discovery screen
│   ├── biochemical_similarity.py         # Multi-metric sequence scoring and sensitivity
│   ├── add_mixmhc2pred_validation.py     # Multi-predictor consensus (MixMHC2pred integration)
│   ├── run_null_model.py                 # 10,000-replicate constrained-null simulations
│   ├── run_focused_taldo1_register_analysis.py # TALDO1 register stability assessment
│   ├── same_register_af3_analysis.py     # AlphaFold 3 parent-vs-nested structural evaluation
│   ├── lead_focused_robustness.py        # Lead-candidate robustness dashboard generation
│   ├── multiallele_manuscript_analysis.py# Multi-allele comparative cross-context evaluation
│   └── ...                               # Domain-specific modules and figure generators
│
├── revision_analyses/    # Latest computational tightening & reviewer-requested robustness analyses
│   ├── run_proteome_sequence_null.py     # Reviewed human proteome sequence null (11.25M 9-mers)
│   ├── run_proteome_binding_filter.py    # Two-predictor top-1,000 binding prefilter
│   ├── run_iedb_proteome_binding_filter.py# IEDB binding consensus filter
│   ├── run_structural_ensemble_reanalysis.py # Corrected AF3 C-alpha RMSD & medoid calculations
│   ├── build_register_uncertainty_summary.py # 18-scenario register shift analysis
│   ├── run_metric_rank_robustness.py     # Four-metric sensitivity and rank stability
│   ├── run_threshold_sensitivity.py      # Threshold grid sensitivity (1% to 20%)
│   ├── build_structure_benchmark_inventory.py# Solved-structure benchmark inventory (144 PDBs)
│   ├── build_evidence_and_control_audit.py   # HLA Ligand Atlas & PRIDE PXD068488 audit
│   └── run_computational_tightening.py   # Integrated computational tightening runner
│
├── batch_builders/       # Structural model batch generators & AlphaFold 3 orchestration
│   ├── build_multiallele_af3_batches.py  # Builds AlphaFold 3 server JSON batches across 4 alleles
│   ├── build_missing_v2_af3_batches.py   # Reconciles and formats missing AF3 candidate runs
│   ├── validate_multiallele_af3_package.py# Validates AF3 JSON syntax and chain composition
│   ├── build_stuve_meeting_prep.py       # Data formatting for mentor technical reviews
│   ├── build_stuve_project_flow.py       # Project flow architecture generator
│   └── build_olivia_master_brief.py      # Master briefing synthesis utility
│
├── tests/                # 47 Automated validation & regression test suites
│   ├── test_repository_integrity.py      # Assertions on project structure and file sizes
│   ├── test_verify_package.py            # Checksum and package manifest assertions
│   └── test_*.py                         # Numerical regression tests and assertions
│
├── staging/              # Release packaging & GitHub repository synchronization
│   └── stage_github_release.py           # Validates repository size, checksums, and staging
│
└── archive/              # Preserved early exploratory batch scripts & prototypes
    ├── molecular_mimicry_pipeline_v3_fixed.py
    ├── analyze_af_server_results.mjs
    ├── build_af_server_batch.mjs
    ├── build_benchmark_manifest.mjs
    ├── build_tcrmodel2_fasta.mjs
    ├── analyze_hy_interfaces.py
    └── analyze_tcrmodel2_hy.py
```

---

## 2. Standard Execution Workflows

### Dataset Preparation:
To reprocess raw IEDB and UniProt downloads into clean, normalized row-level tables:
```bash
python3 code/pipeline/prepare_datasets.py
```

### Candidate Prioritization Screen:
To execute the multi-allele discovery screen across all 6,400 peptide pairs:
```bash
python3 code/pipeline/rank_full_pmhc_screen.py
```

### Reviewer Robustness & Computational Tightening:
```bash
python3 code/revision_analyses/run_computational_tightening.py
```

### Running Validation Tests:
To run the complete repository integrity and regression test suite:
```bash
python3 code/tests/test_repository_integrity.py
```

### Package Checksum Verification:
```bash
python3 computational_package/scripts/verify_package.py
```
