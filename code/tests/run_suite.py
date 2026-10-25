#!/usr/bin/env python3
"""Run every test module in ``code/tests/`` and report one honest summary.

Each file in this directory is an executable ``unittest`` module, so the suite
is run by executing them as subprocesses rather than by importing them. That
keeps each module's own ``sys.path`` bootstrap in effect and matches how the
tests are invoked by hand.

Three outcomes are reported:

``PASS``  the module exited 0.
``FAIL``  the module exited non-zero and its inputs were available.
``SKIP``  the module needs data this repository does not contain, and that data
          is not reachable on this machine. Skips are listed individually with
          their reason, never folded silently into the total.

Two groups of tests read an intermediate data tier that was never committed:
``processed/``, which lives outside this repository (see
``computational_package/docs/SCOPE_AND_LIMITS.md``). When that tier is found,
this script links it to ``code/processed`` -- the path the analysis scripts
expect -- and runs those tests normally. The link is gitignored.

Usage::

    python3 code/tests/run_suite.py
    python3 code/tests/run_suite.py --no-link      # never create the link
    python3 code/tests/run_suite.py -v             # show output of failures
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "code" / "tests"
LINK = ROOT / "code" / "processed"

# Candidate locations of the external intermediate tier, most likely first.
# Override with EBV_MS_PROCESSED_ROOT.
TIER_CANDIDATES = (
    ROOT.parent / "ebv_ms_publication" / "processed",
    Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/processed",
)

# Pass only when the external processed/ tier is reachable. Measured 2026-10-24.
REQUIRES_PROCESSED_TIER = {
    "test_candidate_expansion",
    "test_high_yield_candidate_evidence",
    "test_high_yield_control_validation",
    "test_high_yield_register_resolution",
    "test_hla2_positive_control_benchmark_v2",
    "test_hla2_positive_control_results_analysis",
    "test_lead_focused_audit_outputs",
    "test_physicochemical_similarity_upgrade",
    "test_run_pmhc_surface_electrostatics_sequence_expansion",
}

# Inputs absent from every tree on this machine, so these cannot pass anywhere
# here. Each is a real gap, recorded so it stays visible rather than being
# mistaken for a passing suite. See CODEBASE_AUDIT_2026-10-24.md section C.
MISSING_UPSTREAM_DATA = {
    "test_audit_complete_model_pipeline": "no AlphaFold 3 model root present (expects one; found none)",
    "test_build_missing_v2_af3_batches": "reads tcell_library_v2_2026-08-22/ from the CloudDocs "
    "ebv_ms_publication tree, where that file is absent",
    "test_multiallele_validator": "outputs/alphafold_multiallele_5x30_2026-08-20/ "
    "(5 AF3 job batches) exists in no tree on this machine",
}


def find_tier() -> Path | None:
    """Return the external processed/ tier, or None if it is not reachable."""
    override = os.environ.get("EBV_MS_PROCESSED_ROOT")
    candidates = (Path(override),) + TIER_CANDIDATES if override else TIER_CANDIDATES
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


def ensure_link(tier: Path | None, allow_link: bool) -> bool:
    """Point code/processed at the tier. Return True if the tier is usable."""
    if LINK.is_dir():
        return True
    if tier is None or not allow_link:
        return False
    if LINK.is_symlink():  # dangling link from an earlier run
        LINK.unlink()
    LINK.symlink_to(Path(os.path.relpath(tier, LINK.parent)))
    print(f"linked code/processed -> {tier}")
    return LINK.is_dir()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-link", action="store_true", help="never create code/processed")
    parser.add_argument("-v", "--verbose", action="store_true", help="print failure output")
    args = parser.parse_args()

    tier = find_tier()
    have_tier = ensure_link(tier, allow_link=not args.no_link)
    if not have_tier:
        print("external processed/ tier not available; dependent tests will be skipped")

    passed: list[str] = []
    failed: list[tuple[str, str]] = []
    skipped: list[tuple[str, str]] = []

    for path in sorted(TESTS.glob("test_*.py")):
        name = path.stem
        if name in MISSING_UPSTREAM_DATA:
            skipped.append((name, MISSING_UPSTREAM_DATA[name]))
            continue
        if name in REQUIRES_PROCESSED_TIER and not have_tier:
            skipped.append((name, "needs the external processed/ tier"))
            continue
        result = subprocess.run(
            [sys.executable, str(path)],
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        if result.returncode == 0:
            passed.append(name)
        else:
            failed.append((name, (result.stderr or result.stdout).strip()))

    total = len(passed) + len(failed) + len(skipped)
    for name, reason in skipped:
        print(f"SKIP  {name}: {reason}")
    for name, output in failed:
        print(f"FAIL  {name}")
        if args.verbose:
            print("\n".join("        " + line for line in output.splitlines()[-20:]))

    print(f"\n{len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped, of {total}")
    if failed and not args.verbose:
        print("re-run with -v to see failure output")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
