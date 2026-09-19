import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
BUILDER = WORKSPACE / "build_missing_v2_af3_batches.py"
MODEL_MANIFEST = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/"
    "ebv_ms_publication/processed/tcell_library_v2_2026-08-22/"
    "model_inventory_320.csv"
)
REMAINING_72 = WORKSPACE / "outputs/alphafold_v2_remaining_snapshot_2026-08-23.csv"
UNLISTED_9 = (
    WORKSPACE
    / "outputs/alphafold_v2_download_audit_2026-08-23/unlisted_missing_9.csv"
)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class MissingV2BatchBuilderTests(unittest.TestCase):
    def test_builds_the_verified_81_jobs_as_30_30_21(self):
        expected_rows = read_csv(MODEL_MANIFEST)
        expected_by_name = {row["job_name"]: row for row in expected_rows}
        missing_names = {
            row["job_name"] for row in read_csv(REMAINING_72) + read_csv(UNLISTED_9)
        }
        self.assertEqual(len(missing_names), 81)

        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, str(BUILDER), "--output-dir", directory],
                cwd=WORKSPACE,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            output = Path(directory)
            batch_files = sorted(output.glob("ebvms_v2_missing_batch_*_jobs.json"))
            self.assertEqual([path.name for path in batch_files], [
                "ebvms_v2_missing_batch_01_30_jobs.json",
                "ebvms_v2_missing_batch_02_30_jobs.json",
                "ebvms_v2_missing_batch_03_21_jobs.json",
            ])
            batches = [json.loads(path.read_text(encoding="utf-8")) for path in batch_files]

        self.assertEqual([len(batch) for batch in batches], [30, 30, 21])
        jobs = [job for batch in batches for job in batch]
        names = [job["name"] for job in jobs]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(set(names), missing_names)

        for job in jobs:
            self.assertEqual(
                set(job), {"name", "modelSeeds", "sequences", "dialect", "version"}
            )
            self.assertEqual(job["modelSeeds"], [])
            self.assertEqual(job["dialect"], "alphafoldserver")
            self.assertEqual(job["version"], 1)
            self.assertEqual(len(job["sequences"]), 3)
            chains = [entity["proteinChain"] for entity in job["sequences"]]
            self.assertTrue(all(chain["count"] == 1 for chain in chains))
            self.assertEqual(
                chains[2]["sequence"], expected_by_name[job["name"]]["peptide_sequence"]
            )


if __name__ == "__main__":
    unittest.main()
