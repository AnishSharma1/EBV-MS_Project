# EBV-MS computational package

This repository preserves the computational workstream for an EBV–multiple-sclerosis molecular-mimicry project. The organized public entry point is [`computational_package/`](computational_package/README.md).

## Start here

1. Read the [scope and evidence limits](computational_package/docs/SCOPE_AND_LIMITS.md).
2. Review the [data and artifact policy](computational_package/docs/DATA_AND_ARTIFACT_POLICY.md).
3. Inspect the frozen input libraries and result summaries in `computational_package/data/` and `computational_package/results/`.
4. Run `python3 computational_package/scripts/verify_package.py` to check the public-package manifest.

The repository retains historical notebooks and earlier scripts for provenance. They are not the recommended entry point for interpreting the current package.

## Scientific boundary

The included sequence, pMHC-geometry, docking-calibration, and electrostatics analyses are computational prioritization and methodological evidence only. They do **not** establish natural presentation, TCR binding or recognition, activation, cross-reactivity, molecular mimicry, or a multiple-sclerosis mechanism.

The electrostatics control gate is preserved as a failed scientific gate. That completed computational result prevents its use for candidate ranking.
