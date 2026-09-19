#!/usr/bin/env python3
"""Validation tests for EBV-MS project organization and repository integrity."""

import subprocess
import sys
from pathlib import Path

# Find project root robustly whether called directly or via symlink
_curr = Path(__file__).resolve().parent
while _curr != _curr.parent:
    if (_curr / "manuscript").is_dir() and (_curr / "ANTIGRAVITY_PROJECT_BRIEF.md").is_file():
        ROOT = _curr
        break
    _curr = _curr.parent
else:
    ROOT = Path(__file__).resolve().parents[2]


def test_manuscript_deliverables_exist():
    manuscript_dir = ROOT / "manuscript"
    assert manuscript_dir.is_dir(), f"manuscript/ directory missing at {manuscript_dir}"
    
    required_files = [
        "main.tex",
        "supplement.tex",
        "paper_v9_working_main.pdf",
        "paper_v9_working_supplement.pdf",
        "BMC_SUBMISSION_FIELDS.md",
        "README.md",
        "references.tex",
    ]
    for filename in required_files:
        path = manuscript_dir / filename
        assert path.is_file(), f"Missing manuscript deliverable: {filename}"
        assert path.stat().st_size > 0, f"Empty file: {filename}"


def test_documentation_and_inventory_exist():
    assert (ROOT / "README.md").is_file(), "Root README.md missing"
    assert (ROOT / "PROJECT_DIRECTORY_MAP.md").is_file(), "PROJECT_DIRECTORY_MAP.md missing"
    assert (ROOT / "ANTIGRAVITY_PROJECT_BRIEF.md").is_file(), "ANTIGRAVITY_PROJECT_BRIEF.md missing"
    assert (ROOT / "Model_Library/STRUCTURAL_MODELS_INVENTORY.md").is_file(), "STRUCTURAL_MODELS_INVENTORY.md missing"
    assert (ROOT / "plans/README.md").is_file(), "plans/README.md missing"
    assert (ROOT / "code/README.md").is_file(), "code/README.md missing"
    assert (ROOT / "outputs/README.md").is_file(), "outputs/README.md missing"
    assert (ROOT / "figures/README.md").is_file(), "figures/README.md missing"
    assert (ROOT / "logs/README.md").is_file(), "logs/README.md missing"


def test_computational_package_passes():
    verifier = ROOT / "computational_package/scripts/verify_package.py"
    assert verifier.is_file(), "computational_package verifier script missing"
    
    proc = subprocess.run(
        [sys.executable, str(verifier)],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert proc.returncode == 0, f"computational package verification failed:\n{proc.stdout}\n{proc.stderr}"
    assert "PACKAGE VALIDATION PASSED" in proc.stdout


def test_no_forbidden_root_clutter():
    assert not (ROOT / "io.mc").exists(), "io.mc should not be at root"
    assert not (ROOT / "Rplots.pdf").exists(), "Rplots.pdf should not be at root"
    assert not (ROOT / "EBV_MS_manuscript_draft.pdf").exists(), "Old draft should not be at root"


if __name__ == "__main__":
    test_manuscript_deliverables_exist()
    test_documentation_and_inventory_exist()
    test_computational_package_passes()
    test_no_forbidden_root_clutter()
    print("ALL REPOSITORY INTEGRITY TESTS PASSED")
