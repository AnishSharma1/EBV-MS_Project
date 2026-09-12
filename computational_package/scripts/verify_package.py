#!/usr/bin/env python3
"""Validate the compact public EBV-MS computational package."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
REQUIRED = (
    "README.md",
    "MANIFEST.sha256",
    "docs/SCOPE_AND_LIMITS.md",
    "docs/DATA_AND_ARTIFACT_POLICY.md",
    "data/tcell_library_v2_2026-08-22/README.md",
    "data/high_yield_control_validation_2026-08-28/RESULTS_SUMMARY.md",
    "results/electrostatics_controls_2026-08-30/control_gate.json",
    "results/stage1_evidence_review_2026-09-03/01_stage1_assay_recommendations.csv",
    "results/tcrmodel2_refresh_2026-09-08/ANALYSIS_SUMMARY.md",
)
FORBIDDEN_PARTS = {
    "runtime",
    "raw_calculations",
    "aligned_models",
    "model_calculation_records",
    "sampled_surface_vectors",
    "io.mc",
    "__pycache__",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def failures() -> list[str]:
    errors = [f"missing required file: {relative}" for relative in REQUIRED if not (PACKAGE / relative).is_file()]
    for path in PACKAGE.rglob("*"):
        if set(path.parts) & FORBIDDEN_PARTS:
            errors.append(f"forbidden artifact included: {path.relative_to(PACKAGE)}")

    gate_path = PACKAGE / "results/electrostatics_controls_2026-08-30/control_gate.json"
    if gate_path.is_file():
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        expected = {
            "status": "fail",
            "candidate_evaluation_allowed": False,
            "discovery_unlock_allowed": False,
            "electrostatics_retired_from_candidate_ranking": True,
        }
        for field, value in expected.items():
            if gate.get(field) != value:
                errors.append(f"unexpected electrostatics gate {field!r}: {gate.get(field)!r}")

    manifest = PACKAGE / "MANIFEST.sha256"
    if manifest.is_file():
        for line in manifest.read_text(encoding="utf-8").splitlines():
            expected_hash, relative = line.split("  ", 1)
            target = PACKAGE / relative
            if not target.is_file():
                errors.append(f"manifest target missing: {relative}")
            elif sha256(target) != expected_hash:
                errors.append(f"checksum mismatch: {relative}")
    return errors


if __name__ == "__main__":
    errors = failures()
    if errors:
        print("PACKAGE VALIDATION FAILED")
        print("\n".join(f"- {error}" for error in errors))
        sys.exit(1)
    print("PACKAGE VALIDATION PASSED")
