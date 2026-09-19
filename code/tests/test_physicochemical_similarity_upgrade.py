import math
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from amino_acid_physicochemical_distance import (  # noqa: E402
    ATCHLEY_FACTORS,
    atchley_residue_distance,
    grantham_residue_distance,
    tcr_face_atchley_distance,
    tcr_face_grantham_mismatch,
)
from build_physicochemical_similarity_upgrade import (  # noqa: E402
    DEFAULT_CONTROL_ROOT,
    DEFAULT_DISCOVERY,
    DEFAULT_FROZEN_RANKS,
    competition_ranks,
    evaluate_grantham_promotion,
    read_csv,
    run,
)


class GranthamDistanceTests(unittest.TestCase):
    def test_matches_published_table_entries(self):
        self.assertEqual(grantham_residue_distance("A", "R"), 112)
        self.assertEqual(grantham_residue_distance("I", "L"), 5)
        self.assertEqual(grantham_residue_distance("C", "W"), 215)

    def test_identity_symmetry_and_tcr_facing_mean(self):
        self.assertEqual(grantham_residue_distance("A", "A"), 0)
        self.assertEqual(
            grantham_residue_distance("W", "C"),
            grantham_residue_distance("C", "W"),
        )
        self.assertEqual(tcr_face_grantham_mismatch("AAAAAAAAA", "AAAAAAAAA"), 0.0)
        self.assertEqual(
            tcr_face_grantham_mismatch("AAAAAAAAA", "ARAAAAAAA"),
            112 / 5,
        )

    def test_rejects_noncanonical_residues_and_non_nine_mer_cores(self):
        with self.assertRaisesRegex(ValueError, "unsupported amino acid"):
            grantham_residue_distance("A", "X")
        with self.assertRaisesRegex(ValueError, "exact nine-residue cores"):
            tcr_face_grantham_mismatch("AAAAAAAA", "AAAAAAAAA")


class AtchleyDistanceTests(unittest.TestCase):
    def test_uses_published_factor_scores(self):
        self.assertEqual(ATCHLEY_FACTORS["A"], (-0.591, -1.302, -0.733, 1.570, -0.146))
        expected = math.sqrt(sum(
            (left - right) ** 2
            for left, right in zip(ATCHLEY_FACTORS["A"], ATCHLEY_FACTORS["C"])
        ))
        self.assertAlmostEqual(atchley_residue_distance("A", "C"), expected, places=12)

    def test_identity_symmetry_and_tcr_facing_mean(self):
        self.assertEqual(atchley_residue_distance("Y", "Y"), 0.0)
        self.assertAlmostEqual(
            atchley_residue_distance("W", "Y"),
            atchley_residue_distance("Y", "W"),
            places=12,
        )
        self.assertEqual(tcr_face_atchley_distance("AAAAAAAAA", "AAAAAAAAA"), 0.0)


class RankAndGateTests(unittest.TestCase):
    def test_competition_ranks_are_deterministic_under_ties(self):
        self.assertEqual(competition_ranks([0.4, 0.1, 0.1, 0.8]), [3, 1, 1, 4])
        self.assertEqual(
            competition_ranks([0.4, 0.1, 0.1, 0.8], lower_is_better=False),
            [2, 3, 3, 1],
        )

    def test_promotion_requires_lexicographic_noninferiority(self):
        current = {
            "system_capture_at_3_count": 3,
            "worst_system_rank": 2,
            "system_weighted_mrr": 0.83333333,
        }
        equal = dict(current)
        worse_worst_rank = {
            "system_capture_at_3_count": 3,
            "worst_system_rank": 3,
            "system_weighted_mrr": 1.0,
        }
        better = {
            "system_capture_at_3_count": 3,
            "worst_system_rank": 1,
            "system_weighted_mrr": 1.0,
        }
        self.assertTrue(evaluate_grantham_promotion(current, equal)["promote_grantham"])
        self.assertFalse(
            evaluate_grantham_promotion(current, worse_worst_rank)["promote_grantham"]
        )
        self.assertTrue(evaluate_grantham_promotion(current, better)["promote_grantham"])


class EndToEndTests(unittest.TestCase):
    def test_frozen_controls_gate_and_full_universe(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary)
            manifest = run(
                discovery_path=DEFAULT_DISCOVERY,
                control_root=DEFAULT_CONTROL_ROOT,
                frozen_rank_path=DEFAULT_FROZEN_RANKS,
                out=out,
            )
            self.assertTrue(manifest["control_universe_verified_against_frozen_current_metric"])
            self.assertTrue(manifest["promote_grantham"])
            self.assertEqual(manifest["grantham_objective"], [3, -1, 1.0])
            rows = read_csv(out / "all_6400_pair_physicochemical_sensitivity.csv")
            self.assertEqual(len(rows), 6400)
            self.assertEqual(
                {allele: sum(row["allele"] == allele for row in rows) for allele in manifest["discovery_pairs_by_allele"]},
                manifest["discovery_pairs_by_allele"],
            )
            self.assertTrue(all(row["robustness_label"] for row in rows))
            audited = read_csv(out / "audited_balf5_taldo1_pairs.csv")
            self.assertEqual(len(audited), 2)
            self.assertEqual(
                {(row["allele"], row["tcr_face_current_mismatch_rank"]) for row in audited},
                {("HLA-DRB1*13:03", "29"), ("HLA-DRB1*15:01", "452")},
            )


if __name__ == "__main__":
    unittest.main()
