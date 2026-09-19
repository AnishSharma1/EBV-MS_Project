import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from focused_taldo1_register_analysis import classify_register_evidence, enumerate_pair_windows


class FocusedTaldo1RegisterAnalysisTests(unittest.TestCase):
    def test_consensus_requires_every_independent_predictor_family(self):
        rows = [
            {"predictor_family": "IEDB recommended", "core": "ABCDEFGHI"},
            {"predictor_family": "NetMHCIIpan 4.3", "core": "ABCDEFGHI"},
            {"predictor_family": "MixMHC2pred 2.1", "core": "BCDEFGHIJ"},
        ]
        result = classify_register_evidence(rows)
        self.assertEqual(result["status"], "predictor_register_disagreement")
        self.assertEqual(result["distinct_core_count"], 2)
        self.assertEqual(result["majority_core"], "ABCDEFGHI")
        self.assertFalse(result["experimentally_resolved"])

    def test_consensus_does_not_convert_prediction_to_experimental_resolution(self):
        rows = [
            {"predictor_family": "NetMHCIIpan 4.3", "core": "ABCDEFGHI"},
            {"predictor_family": "MixMHC2pred 2.1", "core": "ABCDEFGHI"},
        ]
        result = classify_register_evidence(rows)
        self.assertEqual(result["status"], "predictor_register_consensus")
        self.assertFalse(result["experimentally_resolved"])

    def test_all_windows_are_enumerated_for_both_15mers(self):
        rows = enumerate_pair_windows("ACDEFGHIKLMNPQR", "RSTVWYACDEFGHIK")
        self.assertEqual(len(rows), 49)
        self.assertEqual(rows[0]["left_start_1_based"], 1)
        self.assertEqual(rows[-1]["right_start_1_based"], 7)


if __name__ == "__main__":
    unittest.main()
