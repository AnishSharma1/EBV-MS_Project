# Archived Scripts and Early Plans

This directory contains standalone exploratory scripts, early batch generators, and initial draft project planning documents created during development:

| File | Purpose / Historical Context |
|---|---|
| `MS_EBV_Molecular_Mimicry_Project_Plan.docx` | Initial conceptual draft of the molecular mimicry study plan. |
| `molecular_mimicry_pipeline_v3_fixed.py` | Pipeline v3 implementation (subsequently refactored into modular components in `computational_package/scripts/`). |
| `analyze_af_server_results.mjs` | Node.js utility for parsing AlphaFold Server zip downloads. |
| `build_af_server_batch.mjs` | Node.js utility for formatting batch JSON inputs for AlphaFold Server. |
| `build_benchmark_manifest.mjs` | Manifest creation script for initial positive-control benchmarks. |
| `build_tcrmodel2_fasta.mjs` | FASTA preparation helper for TCRModel2 runs. |
| `build_project_plan.py` | Python utility for structuring project task tables. |
| `analyze_hy_interfaces.py` | Interface contact analysis script for the Hy.2E11 BALF5-MBP control system. |
| `analyze_tcrmodel2_hy.py` | Structural contact scoring script for TCRModel2 Hy.2E11 models. |
| `compare_af3_calibrators.py` | Benchmark comparison script for AlphaFold 3 control structures. |
| `compare_af3_pmhc_controls.py` | pMHC control RMSD comparison script. |
| `compare_tcrmodel2_calibrators.py` | TCRModel2 calibrator DockQ score calculation. |

---

> [!NOTE]
> These files are retained for provenance and commit history integrity. For all reproducible analyses supporting the publication, refer to [`computational_package/`](../computational_package/).
