import csv
import json
from pathlib import Path
import tempfile
import unittest


from candidate_expansion import (
    resolve_natural_context,
    select_diverse_candidates,
    summarize_candidate_predictors,
)
from build_candidate_expansion import prepare_expansion_package


class NaturalContextTests(unittest.TestCase):
    def test_known_parent_sequence_produces_exact_mixmhc_context(self):
        peptide = "DARLSFDKDAMVARA"
        records = [
            {
                "accession": "P37837",
                "protein": "Transaldolase",
                "sequence": "AAADARLSFDKDAMVARAGGG",
            }
        ]

        resolved = resolve_natural_context(peptide, records)

        self.assertEqual(resolved["source_accession"], "P37837")
        self.assertEqual(resolved["source_start_1_based"], 4)
        self.assertEqual(resolved["mixmhc_context"], "AAADARARAGGG")
        self.assertEqual(resolved["mixmhc_context_status"], "exact_parent_context")

    def test_missing_parent_remains_explicitly_not_evaluable(self):
        resolved = resolve_natural_context(
            "DARLSFDKDAMVARA",
            [{"accession": "P1", "protein": "other", "sequence": "AAAAAAAAAAAA"}],
        )

        self.assertEqual(resolved["source_accession"], "")
        self.assertEqual(resolved["mixmhc_context_status"], "not_evaluable_missing_parent")
        self.assertEqual(resolved["mixmhc_context"], "")


class CandidateSelectionTests(unittest.TestCase):
    @staticmethod
    def candidate(
        allele,
        rank,
        pair_id,
        ebv_core,
        self_core,
        ebv_protein,
        self_protein,
    ):
        return {
            "allele": allele,
            "hla_rank": str(rank),
            "pair_id": pair_id,
            "ebv_candidate_id": f"E_{pair_id}",
            "self_candidate_id": f"S_{pair_id}",
            "ebv_protein": ebv_protein,
            "self_protein": self_protein,
            "ebv_sequence": f"AAA{ebv_core}AAA",
            "self_sequence": f"BBB{self_core}BBB",
            "ebv_core": ebv_core,
            "self_core": self_core,
            "ebv_binding_percentile_rank": "5",
            "self_binding_percentile_rank": "5",
            "surface_status": "complete",
            "left_model_count": "5",
            "right_model_count": "5",
            "declared_register_status": "iedb_resolved_unique_fully_contained",
        }

    def test_selection_excludes_frozen_arms_and_same_core_duplicates(self):
        rows = [
            self.candidate("A", 1, "frozen", "ABCDEFGHI", "JKLMNOPQR", "E0", "S0"),
            self.candidate("A", 2, "best", "BCDEFGHIK", "KLMNPQRST", "E1", "S1"),
            self.candidate("A", 3, "nested_duplicate", "BCDEFGHIK", "KLMNPQRST", "E1", "S1"),
            self.candidate("B", 1, "other", "CDEFGHIKL", "LMNPQRSTV", "E2", "S2"),
        ]
        frozen = [
            {
                "pair_id": "frozen",
                "ebv_candidate_id": "E_frozen",
                "self_candidate_id": "S_frozen",
            }
        ]

        result = select_diverse_candidates(
            rows,
            frozen,
            target_count=2,
            max_per_allele=1,
            protein_cap=2,
        )

        self.assertEqual([row["pair_id"] for row in result["selected"]], ["best", "other"])
        reasons = {row["pair_id"]: row["selection_status"] for row in result["provenance"]}
        self.assertEqual(reasons["frozen"], "excluded_frozen_pair_or_arm")
        self.assertEqual(reasons["nested_duplicate"], "excluded_same_hla_core_pair_duplicate")

    def test_selection_abstains_instead_of_padding_with_repeated_cores(self):
        rows = [
            self.candidate("A", 1, "one", "ABCDEFGHI", "JKLMNOPQR", "E1", "S1"),
            self.candidate("A", 2, "two", "ABCDEFGHI", "JKLMNOPQR", "E2", "S2"),
        ]

        result = select_diverse_candidates(
            rows,
            [],
            target_count=2,
            max_per_allele=2,
            protein_cap=2,
        )

        self.assertEqual(len(result["selected"]), 1)
        self.assertEqual(result["status"], "not_evaluable_insufficient_diverse_candidates")


class ExpansionPackageTests(unittest.TestCase):
    def test_prepare_builds_ten_new_targets_and_corrected_lead_recheck(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)

            result = prepare_expansion_package(output_dir=output)

            self.assertEqual(result["new_target_count"], 10)
            self.assertEqual(result["lead_recheck_count"], 2)
            with (output / "frozen_expansion_targets.csv").open(newline="", encoding="utf-8") as handle:
                new_targets = list(csv.DictReader(handle))
            self.assertEqual(len(new_targets), 10)
            self.assertEqual(len({row["ebv_core"] for row in new_targets}), 10)
            self.assertEqual(len({row["self_core"] for row in new_targets}), 10)
            with (output / "prepared_inputs/peptide_arm_registry.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                arms = list(csv.DictReader(handle))
            self.assertEqual(len(arms), 24)
            taldo_arms = [row for row in arms if row["protein"].upper() == "TALDO1"]
            self.assertEqual(len(taldo_arms), 2)
            self.assertTrue(
                all(row["mixmhc_context_status"] == "exact_parent_context" for row in taldo_arms)
            )
            self.assertTrue(all("XXX" not in row["mixmhc_context"] for row in taldo_arms))
            protocol = json.loads((output / "protocol_lock.json").read_text(encoding="utf-8"))
            self.assertFalse(protocol["discovery_unlock_allowed"])
            self.assertEqual(protocol["status"], "prepared_not_evaluated")
            self.assertTrue((output / "SHA256SUMS.csv").exists())


class PredictorScreenTests(unittest.TestCase):
    def test_candidate_requires_both_supported_arms_and_matching_registers(self):
        arms = [
            {"target_id": "T1", "side": "ebv", "binding_consensus": True,
             "binding_supported": True, "register_consensus_matches_declared": True,
             "predictor_status": "complete", "both_predictors_above_20": False},
            {"target_id": "T1", "side": "self", "binding_consensus": False,
             "binding_supported": True, "register_consensus_matches_declared": True,
             "predictor_status": "complete", "both_predictors_above_20": False},
        ]

        summary = summarize_candidate_predictors(arms)[0]

        self.assertEqual(summary["predictor_screen_status"], "advance_with_caution")

    def test_missing_or_register_disagreement_cannot_advance(self):
        arms = [
            {"target_id": "T1", "side": "ebv", "binding_consensus": True,
             "binding_supported": True, "register_consensus_matches_declared": False,
             "predictor_status": "complete", "both_predictors_above_20": False},
            {"target_id": "T1", "side": "self", "binding_consensus": True,
             "binding_supported": True, "register_consensus_matches_declared": True,
             "predictor_status": "complete", "both_predictors_above_20": False},
        ]

        summary = summarize_candidate_predictors(arms)[0]

        self.assertEqual(summary["predictor_screen_status"], "hold_register_or_binding_conflict")


if __name__ == "__main__":
    unittest.main()
