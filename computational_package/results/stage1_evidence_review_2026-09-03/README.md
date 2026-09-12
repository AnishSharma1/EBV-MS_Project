# EBV-MS ranking and candidate tables

Corrected on 2026-09-03. This package separates two candidate-selection tracks that answer different questions. Original research files were not changed.

## Read these first

1. `CANDIDATE_TRACKS_EXPLAINED.md` — explains why the V3 shortlist and the BALF5–TALDO1 recommendations differ.
2. `categorized_tables/01_current_stage1_evidence_review/01_stage1_assay_recommendations.csv` — the later eight-candidate evidence review: two BALF5–TALDO1 medium-priority Stage-1 candidates and six holds.
3. `categorized_tables/01_current_stage1_evidence_review/02_candidate_evidence_matrix.csv` — the evidence behind those eight decisions.
4. `categorized_tables/02_v3_structural_shortlist_and_rankings/01_v3_internal_structural_shortlist_not_final_stage1.csv` — the separate four-pair V3 structural shortlist.
5. `categorized_tables/03_original_v2_same_register/02_top_25_by_hla.csv` — the original V2 ranks, including the two BALF5–TALDO1 pairs.

## How the ordering works

- Category numbers indicate reading order.
- Lower file numbers indicate greater practical importance within a category.
- Letter suffixes group equivalent allele-specific tables.
- Category 09 contains superseded or discarded files retained for auditability.

## Categories

### 01 Current Stage-1 evidence review

An additive audit of eight preselected sequence-supported candidates. It does not rerank the complete discovery universe. The two BALF5–TALDO1 pairs are medium priority for exact-HLA binding and nested-peptide register mapping. The other six are on hold.

### 02 V3 structural shortlist and rankings

The V3 ranking universe and its separate four-pair internal structural shortlist. The file formerly called `candidate_selection_table.csv` belongs here. It should not be described as the later Stage-1 evidence recommendation table.

### 03 Original V2 same-register

The primary original V2 sequence rankings. Use these to discuss where BALF5–TALDO1 ranked before the later evidence checks.

### 04 Earlier DR15 merged development

The earlier combined DRB1*15:01 same-register analysis.

### 05 Register sensitivity and controls

Control-panel ranks and register-shift sensitivity. These explain why a high raw rank may still be held.

### 06 Alternative scoring analyses

Multifeature, RMSD, and early model-based rerankings. These are supporting sensitivity analyses.

### 07 Controls and external validation

Positive-control recovery, external benchmarks, leave-one-out robustness, and OB1A12 ternary evaluations.

### 08 Electrostatics results

Pilot, expanded-candidate, and control-sensitivity ranks. The workflow ran, but the control gate failed, so electrostatics did not support candidate prioritization.

### 09 Superseded or discarded

Earlier control-calibrated rankings and a discarded initial electrostatics analysis. Keep for provenance; do not present as current results.

## Interpretation boundary

No table proves natural presentation, TCR binding, T-cell activation, cross-reactivity, molecular mimicry, or an MS mechanism.
