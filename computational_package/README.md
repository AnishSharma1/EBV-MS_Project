# Organized computational package

This is a compact public release assembled on 2026-09-11 from the project’s versioned research artifacts. It is designed for inspection, verification, and handoff; it is not a claim that all upstream external services can be replayed offline.

## Layout

- `scripts/` — source programs used to build libraries, construct pMHC inputs, evaluate controls, and analyze structural outputs.
- `data/` — frozen, compact input and evidence tables.
- `results/` — selected current result summaries, gate records, and TCRmodel2 calibration reports.
- `docs/` — scope, claims, provenance, and deliberate omissions.
- `MANIFEST.sha256` — SHA-256 hashes for all release files except itself.

## Suggested reading order

1. `docs/SCOPE_AND_LIMITS.md`
2. `data/tcell_library_v2_2026-08-22/README.md`
3. `data/high_yield_control_validation_2026-08-28/RESULTS_SUMMARY.md`
4. `results/electrostatics_controls_2026-08-30/control_gate.json`
5. `results/stage1_evidence_review_2026-09-03/01_stage1_assay_recommendations.csv`
6. `results/tcrmodel2_refresh_2026-09-08/ANALYSIS_SUMMARY.md`

## Integrity check

Run:

```bash
python3 computational_package/scripts/verify_package.py
```

It checks required files, excludes nonportable artifacts, verifies the recorded electrostatics gate state, and recomputes the manifest hashes.
