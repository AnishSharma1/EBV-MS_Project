# Project Leg 1: Exploratory Phase Jupyter Notebooks

This directory preserves the interactive Jupyter notebooks from **Leg 1** (the initial exploratory and prototyping phase) of the EBV–MS project.

---

## 1. Context: Leg 1 vs. Leg 2

The EBV-MS project was developed across two distinct methodological phases:

### Leg 1 (Initial Exploratory & Prototyping Phase)
- **Timeframe & Medium**: Early development; conducted predominantly in interactive Jupyter notebooks (preserved here in `notebooks/`).
- **Scope**: Preliminary peptide extraction from EBV and CNS antigens, heuristic sequence cross-reactivity scoring, early machine learning feature exploration, single-chain AlphaFold 2 and ColabFold structural modeling, and exploratory TCR docking.
- **Findings & Methodological Evolution**:
  - Heuristic sequence matching without register control produced ambiguous candidate rankings.
  - Exploratory ternary docking without experimental receptor calibration proved uninformative for reliable cross-reactivity claims.
  - Biological review revealed that the classic Hy.2E11 BALF5–MBP positive control system operates across two co-expressed but distinct HLA-DR15 haplotype molecules (`BALF5` on `DRB5*01:01` [DR2a; PDB 1H15]; `MBP` on `DRB1*15:01` [DR2b; PDB 1BX2]) rather than a single presenting allele.
  - These lessons established the requirement for a completely rebuilt, control-calibrated computational framework (**Leg 2**).

### Leg 2 (Current Publication Architecture)
- **Timeframe & Medium**: Current production pipeline; implemented in modular, reproducible scripts in [`computational_package/`](../computational_package/) and documented in [`manuscript/`](../manuscript/).
- **Scope**:
  - Systematic evaluation across four HLA class II risk alleles (`DRB1*15:01`, `*13:03`, `*08:01`, `*03:01`) comprising 6,400 within-allele pairs.
  - Strict mathematical focus on predicted TCR-facing residue positions (P2, P3, P5, P7, P8) using normalized BLOSUM62 similarity.
  - Multi-predictor consensus filtering (NetMHCIIpan 4.3 and MixMHC2pred 2.1) and proteome-wide sequence rarity testing (11.25M 9-mers).
  - AlphaFold 3 parent-versus-nested structural ensemble sensitivity and register invariance.
  - Rigorous control-first negative gates (retiring APBS electrostatics and rejecting uncalibrated TCR docking).

---

## 2. Notebook Directory Inventory (Leg 1 Artifacts)

| Notebook | Focus in Leg 1 | Evolution into Leg 2 |
|---|---|---|
| `AlphaFold2.ipynb` | Initial AlphaFold 2 single-chain and preliminary pMHC predictions. | Superseded by AlphaFold 3 multi-seed pMHC ensembles and register sensitivity analysis. |
| `Colabfold.ipynb` | Legacy ColabFold structural predictions for preliminary comparison. | Quarantined to provenance; final structural claims use AF3 server runs with PDB controls. |
| `Cross Reactivity.ipynb` | Early exploratory heuristics for sequence cross-reactivity. | Formalized into the normalized TCR-facing BLOSUM62 equation and 10,000-replicate constrained-null simulations. |
| `Machine Learning Pipeline.ipynb` | Exploratory feature extraction, classifier testing, and embedding analysis. | Refactored into statistical sensitivity suites and multi-feature concordance rankings. |
| `Peptide Generation.ipynb` | Initial generation of peptide scanning windows from EBV and myelin antigens. | Replaced by deterministic IEDB database curation with verified UniProt provenance fields. |
| `Score Creation.ipynb` | Early prototypes for scoring similarity between peptide sequences. | Formalized into `biochemical_similarity.py` and `rank_full_pmhc_screen.py`. |
| `TCR Docking.ipynb` | Preliminary exploratory TCR docking experiments. | Evaluated through a rigorous TCRModel2 benchmark; uncalibrated docking was formally blocked from candidate scoring after control gates failed. |

---

> [!NOTE]
> **For Publication Results and Pipeline Reproduction**:  
> Please use the canonical [`computational_package/`](../computational_package/) which contains the verified, reproducible code and frozen data supporting the manuscript. These notebooks are preserved strictly for historical transparency and provenance.
