# Exploratory Phase Jupyter Notebooks

This directory preserves early exploratory workflows, initial hypothesis testing, and exploratory parameter explorations conducted during the initial phases of the EBV–MS project:

| Notebook | Focus / Exploration |
|---|---|
| `AlphaFold2.ipynb` | Initial AlphaFold 2 structural prediction experiments for single-chain and pMHC targets. |
| `Colabfold.ipynb` | Legacy ColabFold runs used for preliminary structure comparisons (superseded by AF3 server ensembles). |
| `Cross Reactivity.ipynb` | Early heuristic exploratory scripts for sequence cross-reactivity screening. |
| `Machine Learning Pipeline.ipynb` | Exploratory feature extraction, classifier testing, and exploratory embedding analysis. |
| `Peptide Generation.ipynb` | Initial generation and parsing of peptide scanning windows from EBV and CNS antigens. |
| `Score Creation.ipynb` | Early scoring logic prototypes preceding the normalized BLOSUM62 equation. |
| `TCR Docking.ipynb` | Preliminary exploratory TCR docking experiments (superseded by calibrated TCRModel2 benchmark). |

---

> [!NOTE]
> **Canonical Reproducible Pipeline**:  
> For the publication results, formal candidate prioritization, frozen tables, and verified pipeline execution, please use [`computational_package/`](../computational_package/). These notebooks are preserved strictly for historical transparency and provenance.
