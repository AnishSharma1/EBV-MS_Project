#!/usr/bin/env python3
"""Stage, validate, and prepare GitHub release for EBV-MS project."""

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

STAGING_DIR = ROOT / "git_staging_EBV-MS_Project"


def check_staging_environment():
    if not STAGING_DIR.is_dir():
        print(f"ERROR: Staging directory not found at {STAGING_DIR}")
        sys.exit(1)
    
    print("[1/4] Checking file sizes in staging (GitHub 100 MB hard limit / 50 MB warning threshold)...")
    oversized = []
    warnings = []
    for file in STAGING_DIR.rglob("*"):
        if file.is_file() and not any(p.startswith(".git") for p in file.parts):
            size_mb = file.stat().st_size / (1024 * 1024)
            if size_mb > 100:
                oversized.append((file.relative_to(STAGING_DIR), size_mb))
            elif size_mb > 50:
                warnings.append((file.relative_to(STAGING_DIR), size_mb))
    
    if oversized:
        print("ERROR: Hard-limit oversized files detected (>100 MB):")
        for f, s in oversized:
            print(f"  - {f}: {s:.2f} MB")
        sys.exit(1)
    
    if warnings:
        print("  Notice: Existing files between 50 MB and 100 MB (allowed by GitHub with warning):")
        for f, s in warnings:
            print(f"    * {f}: {s:.2f} MB")
    
    print("  File size validation passed (no files > 100 MB).")

    print("[2/4] Validating computational package integrity...")
    verifier = STAGING_DIR / "computational_package/scripts/verify_package.py"
    res = subprocess.run([sys.executable, str(verifier)], cwd=STAGING_DIR, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"ERROR: Package verification failed:\n{res.stdout}\n{res.stderr}")
        sys.exit(1)
    print("  Package verification passed.")

    print("[3/4] Staged changes summary:")
    status = subprocess.run(["git", "status", "--short"], cwd=STAGING_DIR, capture_output=True, text=True)
    lines = status.stdout.strip().split("\n")
    print(f"  Total modified / staged entries: {len(lines)}")
    for line in lines[:10]:
        print(f"    {line}")
    if len(lines) > 10:
        print(f"    ... and {len(lines) - 10} more items.")

    print("[4/4] Release Readiness Summary:")
    print("  Staging repository is healthy, validated, and ready for commit.")
    print("  To commit and push the reorganized repository to GitHub:")
    print("    cd \"git_staging_EBV-MS_Project\"")
    print("    git commit -m \"Implement standardized 5-directory repository architecture (plans, code, outputs, figures, logs)\"")
    print("    git push origin main")


if __name__ == "__main__":
    check_staging_environment()
