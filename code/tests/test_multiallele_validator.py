import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
VALIDATOR = WORKSPACE / "validate_multiallele_af3_package.py"


class MultiAlleleValidatorPathTests(unittest.TestCase):
    def test_validator_succeeds_outside_workspace(self):
        with tempfile.TemporaryDirectory() as other_directory:
            result = subprocess.run(
                [sys.executable, str(VALIDATOR)],
                cwd=other_directory,
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Independent validation: PASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
