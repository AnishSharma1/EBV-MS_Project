import unittest

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "code" / "pipeline"))

from build_sequence_electrostatics_summary import classify_combined_support


class CombinedSupportTests(unittest.TestCase):
    def test_top_three_requires_register_qc_for_formal_combined_support(self):
        self.assertEqual(
            classify_combined_support(3, register_robust=True, dielectric_robust=True),
            "sequence_plus_electrostatics_supported",
        )
        self.assertEqual(
            classify_combined_support(3, register_robust=False, dielectric_robust=True),
            "sequence_plus_electrostatics_rank_supported_register_unresolved",
        )

    def test_rank_above_three_or_dielectric_instability_is_not_supported(self):
        self.assertEqual(
            classify_combined_support(4, register_robust=True, dielectric_robust=True),
            "sequence_supported_electrostatics_not_supportive",
        )
        self.assertEqual(
            classify_combined_support(1, register_robust=True, dielectric_robust=False),
            "sequence_supported_electrostatics_dielectric_unstable",
        )


if __name__ == "__main__":
    unittest.main()
