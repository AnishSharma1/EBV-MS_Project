# Threshold-sensitivity QA

Status: PASS

- Frozen candidates: 49 unique pair IDs and 49 unique candidate IDs.
- Peptide arms: 98 unique arm IDs.
- Predictor records: 490 total; exactly 98 records for each of five stored predictor variants.
- Primary grid inputs: exactly three records per arm from the locked predictor set.
- Grid rows: 980 = 49 candidates x 5 cutoffs x 2 binding rules x 2 register rules.
- Exact-register agreement reproduced: 44/49 pairs.
- +/-1-register agreement reproduced: 47/49 pairs.
- Locked strict 20% exact-register gate reproduced exactly HP36 and HP47.
- Every declared and predicted 9-mer core occurred exactly once in its submitted peptide sequence.
- Percentile ranks were validated within the inclusive 0-100 range.
- Original provisional tiers were carried only as annotations and were not recomputed.

Claim boundary: Sensitivity of a computational binding/register gate only; not evidence of measured affinity, endogenous presentation, TCR recognition, cross-reactivity, molecular mimicry, or an MS mechanism.
