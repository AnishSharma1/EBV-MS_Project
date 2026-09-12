# Scope and evidence limits

## What this package contains

- A frozen EBV/human self-peptide library and associated HLA-II prediction records.
- Predeclared control-validation and register-aware ranking records.
- A compact, control-only pMHC electrostatics provenance export.
- The current Stage-1 evidence-review recommendation table, kept distinct from the earlier V3 structural shortlist.
- TCRmodel2 calibration/recovery reports and input-protocol records.

## What the computations mean

The calculations prioritize peptide–HLA hypotheses and evaluate whether methods recover specified structural controls. Confidence scores, register hypotheses, geometrical similarity, electrostatic comparisons, or docking outputs do not demonstrate natural antigen presentation, TCR binding, recognition, activation, specificity, cross-reactivity, molecular mimicry, or disease causality.

## Current gate status

The control-first electrostatics analysis completed its technical checks, but `control_gate.json` records a failed scientific gate. Its candidate-evaluation unlock is therefore unavailable. The package keeps that result visible rather than filtering it out.

## Reproducibility limits

Several source programs depend on externally downloaded IEDB/RCSB records, AlphaFold Server results, PANDORA, APBS, or TCRmodel2. Those systems, versions, and downloads are not bundled here. The included manifests, checksums, source references, and compact outputs support auditing the reported computation; they do not promise a one-command full regeneration.
