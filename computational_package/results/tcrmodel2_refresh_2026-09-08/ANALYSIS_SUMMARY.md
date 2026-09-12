# TCRModel2 refresh: completed analysis

## Intake verification

All four server downloads were found in iCloud Downloads, copied without modifying the originals, and staged under `results/`. Each archive contains `ranked_0.pdb` through `ranked_4.pdb`, `statistics.json`, and `tcr_seqs.json`: five models per job, 20 models total.

| Job | Intended excluded template(s) | Returned pMHC templates | Exclusion check |
|---|---|---|---|
| CAL_1YMM_CORRECTED | 1YMM | 2WBJ, 6R0E, 1FYT, 4Y19 | Pass: 1YMM absent |
| CAL_2WBJ_CANONICAL | 2WBJ | 6R0E, 1FYT, 4Y19, 6CQL | Pass: 2WBJ absent |
| HY_MBP_DRB1_EXPLORATORY | 1YMM, 1ZGL, 2WBJ | 6R0E, 1FYT, 4Y19, 6CQL | Pass: all absent |
| HY_BALF5_DRB5_EXPLORATORY | 1YMM, 1ZGL, 2WBJ | 6R0E, 1FYT, 2IAM, 6CQL | Pass: all absent |

## Calibration gate

DockQ 2.1.3 was rerun with the established reference-chain normalization and `RL:RL` grouping. Its contact metrics support failure, but its automatic receptor choice reverses the requested pMHC-receptor convention for 2WBJ because the deposited native TCR group is longer. The directional component values below therefore come from the sequence-aligned pMHC-superposition audit.

| Control | Best DockQ | Native-contact recovery (`fnat`) | Gate result |
|---|---:|---:|---|
| Corrected 1YMM | 0.1358 | 0.0986 | Fail: below 0.23 |
| Canonical 2WBJ | 0.0375 | 0.0000 | Fail: below 0.23; no native contacts |

Both controls are CAPRI Incorrect across all five returned models. The sequence-aligned rank-0 pMHC C-alpha RMSDs are 0.971 A (1YMM) and 0.649 A (2WBJ); corresponding TCR-fold RMSDs are 2.242 A and 1.202 A. Yet the pMHC-aligned TCR-placement RMSDs are 15.347 A and 17.867 A. Thus the planned control gate fails because docking placement—not plainly bad component folds—is the dominant observed defect.

## Exploratory Hy.2E11 geometry

Across all five models in each condition, MBP/DRB1*15:01 had three stable peptide-contacting TCR residues and BALF5/DRB5*01:01 had five. They shared only D30S and D31I (Jaccard 0.333). This is a model-internal contact description, not evidence of TCR binding, cross-reactivity, activation, molecular mimicry, or an MS mechanism. The conditions also use different HLA beta chains, so this is not a same-allele comparison.

## Files

- `analysis/dockq_calibration/tcrmodel2_grouped_tcr_pmhc_dockq.tsv`: all ten control DockQ scores.
- `analysis/audit_tcrmodel2_components.py`: rerunnable sequence-aligned pMHC, TCR-fold, and pMHC-aligned placement audit. This supersedes `analysis/tcrmodel2_calibrator_structural_metrics.tsv`, whose residue-position matching is not valid for these trimmed model/native chains.
- `analysis/CONTROL_PEPTIDE_WINDOW_AUDIT.md`: raw-deposition-to-submitted-11-mer audit.
- `analysis/tcrmodel2_hy_contact_consensus.json`: original top-model contact summary; the all-five-model result is recorded above because it is the more stringent description.
