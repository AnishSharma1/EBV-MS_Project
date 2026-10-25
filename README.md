# Register-Aware Computation Prioritizes EBV–Self Peptide Pairs for Experimental Testing Across HLA Class II Contexts

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

This study implements an auditable computational prioritization workflow, with explicit control gates, that evaluates candidate EBV–self peptide similarity across four HLA class II contexts:
- `HLA-DRB1*15:01`
- `HLA-DRB1*13:03`
- `HLA-DRB1*08:01`
- `HLA-DRB1*03:01`

We perform a within-allele screen of **6,400 allele-specific pairs** (40 EBV and 40 human self peptides per allele), ranking candidates by normalized BLOSUM62 similarity at the nominally T-cell receptor (TCR)-facing positions (P2, P3, P5, P7, P8) of declared binding cores. Top candidates are examined with four-metric rank sensitivity, a 10,000-replicate constrained-null simulation per allele, predicted-binding threshold sensitivity (NetMHCIIpan 4.3 EL, MixMHC2pred 2.1, and the IEDB-recommended method), AlphaFold 3 parent-versus-nested peptide sensitivity checks, and a PANDORA/MODELLER register stress test.

> [!NOTE]
> The reported analysis did **not** include proteome-wide background sampling or biologically matched HLA-presentable decoys; the manuscript lists both as limitations. Scripts for a proteome-wide null exist under `code/revision_analyses/` (`run_proteome_sequence_null.py`, `run_proteome_binding_filter.py`), but their outputs are not part of this checkout and are not reported in the current manuscript.

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
- **Focus**: The frozen computational package supporting the manuscript. It supports replay from frozen intermediate outputs, not one-command regeneration (see [`computational_package/docs/SCOPE_AND_LIMITS.md`](computational_package/docs/SCOPE_AND_LIMITS.md)).
- **Artifacts**: Implemented in [`computational_package/`](computational_package/) and [`manuscript/`](manuscript/).
- **Upgrades**:
  1. Expansion to a four-allele class II universe with strict within-allele ranking (6,400 pairs).
  2. Score normalization restricted to the nominally TCR-facing positions (P2, P3, P5, P7, P8) of declared cores. The score is register-conditional: alternate nine-mer windows give wide score ranges for both routed pairs.
  3. Predicted-binding threshold sensitivity (NetMHCIIpan 4.3 EL, MixMHC2pred 2.1, IEDB-recommended method). The three outputs are correlated inputs used for experimental routing, not independent measurements.
  4. AlphaFold 3 parent-versus-nested peptide sensitivity analysis. It is descriptive and does not confirm binding registers; the modeled geometry of the DRB1*13:03 TALDO1 arm remains unresolved.
  5. Negative control gates: APBS electrostatics (15 of 30 required control outputs failed, 12 missing) and uncalibrated TCRModel2 docking were benchmarked, failed validation gates, and were excluded from candidate scoring.

---

## 3. Key Prioritized Candidates

| Candidate ID | Allele | Viral Peptide | Human Self Peptide | Screen Rank | Routing Tier | Evidence & Status |
|---|---|---|---|---|---|---|
| **HP47** | `HLA-DRB1*15:01` | BALF5 627–641 (`TGGVYHFVKKHVHES`) | TALDO1 216–230 (`SVTKIYNYYKKFSYK`) | 14 / 1,600 | Experimental-routing class 1 | All three selected outputs met the permissive 20% routing rule (not stricter thresholds). Positive qualitative IEDB purified-MHC binding record for the viral arm (Assay 1774674; Hansen et al. 2007); no self-arm elution record. Routed first for experimental binding and register mapping. |
| **HP36** | `HLA-DRB1*13:03` | BALF5 627–641 (`TGGVYHFVKKHVHES`) | TALDO1 108–122 (`DARLSFDKDAMVARA`) | 13 / 1,600 | Experimental-routing class 2 | All three selected outputs met the permissive 20% routing rule only. No verified allele-specific public binding record; modeled geometry unresolved. Direct binding measurement pending. |

Routing classes are administrative labels for organizing experimental work, not biological-confidence grades (see `manuscript/supplement.tex`, SM2). The cited rank is the frozen-core, within-allele rank; the manuscript reports both pairs as non-significant after multiple-testing adjustment (maxT *p* = 1.0).

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

## 5. Repository Layout

```
EBV-MS_Project/
├── computational_package/            # Canonical frozen package (inputs, results, scripts)
│   ├── MANIFEST.sha256               # Checksum manifest (checked by verify_package.py)
│   ├── scripts/                      # Package scripts, incl. verify_package.py
│   ├── data/                         # Curated input libraries and reference sequences
│   ├── results/                      # Frozen score tables and decision gates
│   └── docs/                         # Scope, limits, and artifact policy
│
├── manuscript/                       # LaTeX sources, version snapshots, audit records
│   ├── main.tex, supplement.tex      # Current sources (V20; see manuscript/V20_CHANGES.md)
│   ├── references.tex                # Shared bibliography
│   ├── figures/                      # figure1-6.pdf; main.tex includes figure1, figure4, figure5
│   └── paper_v11 ... paper_v20 *     # Per-version PDF and text snapshots
│
├── code/                             # See code/README.md for the file inventory
│   ├── pipeline/                     # Analysis scripts (102 Python, 6 R)
│   ├── revision_analyses/            # Revision-round analyses (21 Python)
│   ├── batch_builders/               # AlphaFold 3 batch generators
│   ├── figures/                      # Figure-rendering scripts
│   ├── tests/                        # 51 test modules
│   ├── staging/                      # Release packaging
│   └── archive/                      # Early exploratory scripts
│
├── output/                           # revision_round_4/, revision_round_18/, revision_v11/ (results the V19/V20 text depends on)
├── outputs/                          # ebv_ms_benchmark_manifest/, ebv_ms_model_package/ (no README)
├── figures/                          # manuscript/ (figure1-6), supplement/ (figureS1-S2); see figures/README.md
├── plans/                            # Study designs, protocols, statistical specifications
├── logs/                             # Provenance, citation, QA, and package-verification records
├── doi_release/                      # Staged v1.1.0 archive (built around V18; not the current V20 text)
├── notebooks/                        # Exploratory-phase notebooks (Leg 1)
└── archive/, work/                   # Historical and working material; inspect provenance before reuse
```

> [!NOTE]
> Many scripts and analysis-lock files contain absolute paths from the original author workstation, and some inputs (for example the AlphaFold 3 parent-versus-nested table behind the manuscript's Figure 3) live in the separate `ebv_ms_publication/` workspace. The package supports checking frozen outputs; it does not regenerate every figure from this checkout alone.

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

### Test Suite
```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python code/tests/run_suite.py
```

51 test modules. Verified 2026-10-24 on CPython 3.14.7: **48 passed, 0 failed, 3 skipped** with the external `processed/` tier reachable, or 39 passed with 12 skipped from this repository alone. Skips are printed individually with their reason. `code/tests/README.md` records what each gate covers and which modules are blocked.

Two of these are the gates named in `CLAUDE.md`; all three run in CI alongside package verification:

```bash
python3 code/tests/test_repository_integrity.py             # 9 checks
python3 code/tests/test_peptide_coordinate_integrity.py     # 8 checks
```

### External `processed/` data tier
51 scripts under `code/`, nine tests, and every figure script read an intermediate tier that is **not in this repository**: `ebv_ms_publication/processed/` (8.0 GB, 57 directories), sitting beside this checkout. It is **not version-controlled and has no checksum manifest**, so nothing read from it is covered by `MANIFEST.sha256`. Manuscript results were computed from the frozen inputs in `computational_package/`; the tier matters for re-running scripts, not for auditing reported numbers.

`code/tests/run_suite.py` finds the tier and links it to `code/processed`, the path the scripts expect (gitignored; override with `EBV_MS_PROCESSED_ROOT`). Two copies of `ebv_ms_publication` exist on the authoring machine and have diverged in at least one file, so the tier's provenance is not yet established — see `CODEBASE_AUDIT_2026-10-24.md`.

### Inspecting the Screen
The package supports inspecting frozen result tables and replaying from frozen intermediate outputs. It does not promise one-command regeneration of every result. Relevant entry points:
- **Sequence similarity & screening:** `computational_package/scripts/rank_full_pmhc_screen.py`
- **Multi-metric robustness:** `computational_package/scripts/biochemical_similarity.py`
- **Predictor consensus:** `computational_package/scripts/add_mixmhc2pred_validation.py`
- **Results summaries:** `computational_package/results/stage1_evidence_review_2026-09-03/`

Inspect a script's arguments, inputs, and output paths before running it; several reference locations outside this checkout.

---

## 7. Citation

If you use this computational pipeline, candidate scoring framework, or curated datasets, please cite:

```bibtex
@misc{sharma2026ebvms,
  author = {Sharma, Anish},
  title = {Register-aware computation prioritizes {EBV}--self peptide pairs for experimental testing across {HLA} class {II} contexts},
  year = {2026},
  howpublished = {Manuscript and computational package},
  url = {https://github.com/AnishSharma1/EBV-MS_Project}
}
```
[`CITATION.cff`](CITATION.cff) describes the archived release v1.1.0 (built around the V18 manuscript, under its V18 title), not the current V20 text.

---

## 8. License

- **Code & Scripts**: Licensed under the [MIT License](LICENSE).
- **Manuscript Text, Figures, & Data Tables**: Licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
