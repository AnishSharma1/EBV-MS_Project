#!/usr/bin/env python3
"""Validation tests for EBV-MS_Project repository integrity.

Run from any directory:  python3 code/tests/test_repository_integrity.py
These checks assert structure and cross-document consistency only. They do not
reproduce any scientific result.
"""

import hashlib
import re
import subprocess
import sys
from pathlib import Path


def _find_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in (here, *here.parents):
        if (candidate / "manuscript").is_dir() and (candidate / "computational_package").is_dir():
            return candidate
    raise RuntimeError("repository root not found: expected a directory with manuscript/ and computational_package/")


ROOT = _find_root()
MANUSCRIPT = ROOT / "manuscript"
MANIFEST_DIR = MANUSCRIPT


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _latest_snapshot_pair() -> tuple[Path, Path]:
    pattern = re.compile(r"paper_v(\d+)_main\.pdf")
    mains = sorted(
        (path for path in MANUSCRIPT.iterdir() if pattern.fullmatch(path.name)),
        key=lambda path: int(pattern.fullmatch(path.name).group(1)),
    )
    assert mains, "no paper_v<N>_main.pdf snapshot found in manuscript/"
    main = mains[-1]
    return main, main.with_name(main.name.replace("_main.pdf", "_supplement.pdf"))


def test_manuscript_deliverables_exist():
    assert MANUSCRIPT.is_dir(), f"manuscript/ directory missing at {MANUSCRIPT}"
    for filename in ("main.tex", "supplement.tex", "references.tex", "README.md", "BMC_SUBMISSION_FIELDS.md"):
        path = MANUSCRIPT / filename
        assert path.is_file(), f"Missing manuscript deliverable: {filename}"
        assert path.stat().st_size > 0, f"Empty file: {filename}"
    for snapshot in _latest_snapshot_pair():
        assert snapshot.is_file(), f"Missing version snapshot: {snapshot.name}"
        assert snapshot.stat().st_size > 0, f"Empty snapshot: {snapshot.name}"


def test_documentation_exists():
    for relative in (
        "README.md",
        "plans/README.md",
        "code/README.md",
        "figures/README.md",
        "logs/README.md",
        "computational_package/README.md",
        "computational_package/docs/SCOPE_AND_LIMITS.md",
    ):
        assert (ROOT / relative).is_file(), f"{relative} missing"


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


def test_every_cited_key_has_a_bibitem():
    defined = set(re.findall(r"\\bibitem\{([^}]+)\}", _read(MANUSCRIPT / "references.tex")))
    cited = set()
    for name in ("main.tex", "supplement.tex"):
        for group in re.findall(r"\\cite[a-z]*\{([^}]*)\}", _read(MANUSCRIPT / name)):
            cited.update(key.strip() for key in group.split(","))
    missing = sorted(cited - defined)
    assert not missing, f"cited keys with no \\bibitem in references.tex: {missing}"


def test_included_graphics_exist():
    missing = []
    for name in ("main.tex", "supplement.tex"):
        for relative in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", _read(MANUSCRIPT / name)):
            if not (MANUSCRIPT / relative).is_file():
                missing.append(f"{name}: {relative}")
    assert not missing, f"graphics included by the manuscript but absent: {missing}"


def test_readme_does_not_claim_unreported_analyses():
    readme = _read(ROOT / "README.md")
    for phrase in ("11.25 million", "11.25M", "computational_stratum"):
        assert phrase not in readme, f"README.md still contains {phrase!r}, which the manuscript does not report"


def test_readme_title_matches_manuscript_title():
    def normalize(text: str) -> str:
        return re.sub(r"\s+", " ", text.replace("--", "-").replace("\u2013", "-")).strip().lower()

    title = re.search(r"\\title\{(.+?)\}\s*$", _read(MANUSCRIPT / "main.tex"), re.M)
    assert title, "\\title{} not found in manuscript/main.tex"
    heading = re.search(r"^# (.+)$", _read(ROOT / "README.md"), re.M)
    assert heading, "README.md has no H1 heading"
    assert normalize(heading.group(1)) == normalize(title.group(1)), (
        f"README title {heading.group(1)!r} differs from main.tex title {title.group(1)!r}"
    )


def test_latest_snapshot_matches_recorded_source_hashes():
    """The newest paper_v<N> snapshot must still match the sources it was built from.

    When paper_v<N>_SOURCE_SHA256.txt exists, every file it lists (sources,
    figures, snapshot PDFs and text) must hash identically. A later source edit
    therefore requires a new snapshot version rather than silently drifting.
    """
    main, _supplement = _latest_snapshot_pair()
    version = re.fullmatch(r"paper_v(\d+)_main\.pdf", main.name).group(1)
    record = MANIFEST_DIR / f"paper_v{version}_SOURCE_SHA256.txt"
    if not record.is_file():
        return
    mismatches = []
    for line in _read(record).splitlines():
        expected, relative = line.split("  ", 1)
        path = MANIFEST_DIR / relative
        if not path.is_file():
            mismatches.append(f"{relative}: missing")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            mismatches.append(f"{relative}: changed since the paper_v{version} snapshot")
    assert not mismatches, "sources drifted from the latest snapshot (build a new version):\n" + "\n".join(
        f"  - {m}" for m in mismatches
    )


def _run_all() -> int:
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_") and callable(value)]
    failed = 0
    for test in tests:
        try:
            test()
        except AssertionError as error:
            failed += 1
            print(f"FAIL {test.__name__}\n  {error}\n")
        else:
            print(f"PASS {test.__name__}")
    if failed:
        print(f"\n{failed} of {len(tests)} repository integrity checks failed")
        return 1
    print(f"\nALL {len(tests)} REPOSITORY INTEGRITY TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(_run_all())
