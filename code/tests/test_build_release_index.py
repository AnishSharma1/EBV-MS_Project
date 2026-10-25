#!/usr/bin/env python3
"""Tests for build_release_index.py: PASS must be backed by real packages and artifacts."""

import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "revision_analyses" / "build_release_index.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_release_index_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_package(round_dir: Path, name: str, extra_files: dict[str, str]):
    """Create a package directory whose SHA256SUMS.csv correctly lists its files."""
    package = round_dir / name
    package.mkdir(parents=True)
    rows = []
    for relative, text in extra_files.items():
        path = package / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        rows.append({
            "relative_path": relative,
            "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            "bytes": len(text.encode("utf-8")),
        })
    if not rows:
        placeholder = package / "placeholder.txt"
        placeholder.write_text("x", encoding="utf-8")
        rows.append({
            "relative_path": "placeholder.txt",
            "sha256": hashlib.sha256(b"x").hexdigest(),
            "bytes": 1,
        })
    with (package / "SHA256SUMS.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "sha256", "bytes"])
        writer.writeheader()
        writer.writerows(rows)


class ReleaseIndexTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.temp = tempfile.TemporaryDirectory()
        self.round_dir = Path(self.temp.name) / "revision_round_3"
        self.out_dir = self.round_dir / "RELEASE"
        self.module.ROUND = self.round_dir
        self.module.OUT = self.out_dir

    def tearDown(self):
        self.temp.cleanup()

    def run_main(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            with self.assertRaises(SystemExit) as raised:
                self.module.main()
        return raised.exception.code, stdout.getvalue()

    def build_all_packages(self, skip_artifacts=()):
        by_package = {name: {} for name in self.module.PACKAGES}
        for _addition, _state, artifact in self.module.STATUS_ROWS:
            package, relative = artifact.split("/", 1)
            if artifact not in skip_artifacts:
                by_package[package][relative] = f"content of {artifact}\n"
        for name, files in by_package.items():
            write_package(self.round_dir, name, files)

    def test_missing_packages_exit_2_and_write_nothing(self):
        code, output = self.run_main()
        self.assertEqual(code, 2)
        self.assertIn("7 of 7 component packages are not present", output)
        self.assertFalse(self.out_dir.exists(), "index directory must not be created")

    def test_partial_packages_are_reported_by_name(self):
        write_package(self.round_dir, self.module.PACKAGES[0], {})
        absent = self.module.missing_packages()
        self.assertEqual(absent, list(self.module.PACKAGES[1:]))

    def test_all_packages_and_artifacts_present_passes(self):
        self.build_all_packages()
        with contextlib.redirect_stdout(io.StringIO()):
            self.module.main()
        qa = json.loads((self.out_dir / "RELEASE_QA.json").read_text())
        self.assertEqual(qa["release_status"], "PASS")
        self.assertEqual(qa["checksum_failures"], [])
        with (self.out_dir / "requested_additions_status.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), len(self.module.STATUS_ROWS))
        self.assertTrue(all(row["primary_artifact_present"] == "True" for row in rows))

    def test_status_row_with_absent_artifact_fails_release(self):
        absent = self.module.STATUS_ROWS[0][2]
        self.build_all_packages(skip_artifacts={absent})
        code, output = self.run_main()
        self.assertEqual(code, 1)
        qa = json.loads((self.out_dir / "RELEASE_QA.json").read_text())
        self.assertEqual(qa["release_status"], "FAIL")
        self.assertIn(f"missing_artifact:{absent}", qa["checksum_failures"])

    def test_tampered_package_file_fails_checksum(self):
        self.build_all_packages()
        artifact = self.module.STATUS_ROWS[1][2]
        package, relative = artifact.split("/", 1)
        (self.round_dir / package / relative).write_text("tampered\n", encoding="utf-8")
        code, _output = self.run_main()
        self.assertEqual(code, 1)
        qa = json.loads((self.out_dir / "RELEASE_QA.json").read_text())
        self.assertEqual(qa["release_status"], "FAIL")
        self.assertTrue(any(failure.startswith(f"{package}/checksum:") or failure.startswith(f"{package}/size:")
                            for failure in qa["checksum_failures"]))


if __name__ == "__main__":
    unittest.main()
