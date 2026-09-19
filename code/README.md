# EBV-MS Complete Codebase & Analysis Pipelines

This directory houses the comprehensive collection of **210+ analysis scripts, figure generators, test suites, and batch utilities** powering the EBV–MS project.

---

## 1. Directory Structure

```
code/
├── pipeline/             # 107 Primary analysis and data processing scripts (101 Python, 6 R)
│   ├── prepare_datasets.py               # Rebuilds normalized tables from raw/ inputs
│   ├── rank_full_pmhc_screen.py          # Executes 6,400-pair within-allele discovery screen
│   ├── biochemical_similarity.py         # Multi-metric sequence scoring and sensitivity
│   ├── add_mixmhc2pred_validation.py     # Multi-predictor consensus (MixMHC2pred integration)
│   ├── run_null_model.py                 # 10,000-replicate constrained-null simulations
│   ├── run_focused_taldo1_register_analysis.py # TALDO1 register stability assessment
│   ├── same_register_af3_analysis.py     # AlphaFold 3 parent-vs-nested structural evaluation
│   ├── lead_focused_robustness.py        # Lead-candidate robustness dashboard generation
│   ├── multiallele_manuscript_analysis.py# Multi-allele comparative cross-context evaluation
│   └── ...                               # Domain-specific modules and scoring scripts
│
├── figures/              # Manuscript & Supplement figure rendering scripts (Python & R)
│   ├── build_structural_and_binding.R    # Figure 2 / Figure 3 structural & binding panels
│   ├── build_workflow_and_rank_sensitivity.R # Workflow & ranking sensitivity ggplot figures
│   ├── build_workflow_and_rank_sensitivity.py # Python counterpart for workflow & sensitivity
│   ├── restyle_figures_ggplot.R          # High-resolution styling conforming to BMC guidelines
│   ├── reproduce_figures.R               # Standalone reproduction script for all figures
│   ├── build_panels.py                   # Multi-panel composite figure assembler
│   └── restyle_panels.py                 # Matplotlib/Seaborn panel aesthetic restyler
│
├── revision_analyses/    # 16 Computational tightening & reviewer-requested robustness analyses
│   ├── run_proteome_sequence_null.py     # Human proteome sequence null (11.25M 9-mers)
│   ├── run_proteome_binding_filter.py    # Two-predictor top-1,000 binding prefilter
│   ├── run_iedb_proteome_binding_filter.py# IEDB binding consensus filter
│   ├── run_structural_ensemble_reanalysis.py # Corrected AF3 C-alpha RMSD & medoid calculations
│   ├── build_register_uncertainty_summary.py # 18-scenario register shift analysis
│   ├── run_metric_rank_robustness.py     # Four-metric sensitivity and rank stability
│   ├── run_threshold_sensitivity.py      # Threshold grid sensitivity (1% to 20%)
│   ├── build_structure_benchmark_inventory.py# Solved-structure benchmark inventory (144 PDBs)
│   ├── build_evidence_and_control_audit.py   # HLA Ligand Atlas & PRIDE PXD068488 audit
│   ├── audit_control_manifest.py         # RCSB GraphQL audit of post-cutoff TCR-pMHC controls
│   ├── decide_predictor_gate.py          # Predeclared 9-of-12 structural recovery gate (DockQ >= 0.23)
│   ├── check_manifest_integrity.py       # SHA-256 integrity verification for control panels
│   ├── audit_tcrmodel2_components.py     # Sequence-aligned TCRModel2 diagnostic audit
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
├── tests/                # 49 Automated validation & regression test suites
│   ├── test_repository_integrity.py      # Assertions on project structure and file sizes
│   ├── test_verify_package.py            # Checksum and package manifest assertions
│   ├── test_build_missing_v2_af3_batches.py # Batch generator regression testing
│   ├── test_multiallele_validator.py     # JSON validation test for multi-allele sets
│   └── test_*.py                         # Comprehensive numerical & biological regression suites
│
├── staging/              # Release packaging & GitHub repository synchronization
│   └── stage_github_release.py           # Validates repository size, checksums, and staging
│
└── archive/              # Preserved early exploratory batch scripts, prototypes & notebooks
    ├── build_manuscript.py               # Early ReportLab PDF compilation script
    ├── colabfold_batch_cell.py           # ColabFold execution notebook snippet
    ├── molecular_mimicry_pipeline_v3_fixed.py
    ├── analyze_af_server_results.mjs
    ├── build_af_server_batch.mjs
    ├── build_benchmark_manifest.mjs
    ├── build_tcrmodel2_fasta.mjs
    ├── analyze_hy_interfaces.py
    └── analyze_tcrmodel2_hy.py
```

---

## 2. Note on Auto-Generated Files vs. Source Code

During automated homology modeling runs (e.g. via PANDORA/Modeller), execution engines automatically emit transient loop-refinement scripts (such as `MyLoop.py`, `cmd_modeller.py`, and `cmd_modeller_ini.py`) inside specific output run directories under `outputs/`. These ~300 files are generated by external packages at runtime and are preserved in local outputs for reproducibility, but are excluded from `code/` to keep author source code clean, modular, and maintainable.

Similarly, unpacked third-party dependencies (such as local `numpy` and `BioSQL` site-packages in benchmark runtimes) are kept isolated and excluded from version control.

---

## 3. Standard Execution Workflows

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

### Figure Generation:
To generate all publication figures (Figures 1–4 and Supplementary Figures):
```bash
Rscript code/figures/reproduce_figures.R
python3 code/figures/build_panels.py
```

### Reviewer Robustness & Computational Tightening:
```bash
python3 code/revision_analyses/run_computational_tightening.py
```

### Structural Recovery Gating:
```bash
python3 code/revision_analyses/decide_predictor_gate.py <results_table.csv>
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
