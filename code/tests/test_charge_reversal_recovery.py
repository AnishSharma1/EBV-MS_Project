from __future__ import annotations

import unittest

from charge_reversal_recovery import (
    classify_consensus,
    make_af3_ablation_jobs,
    pandora_calibration_runs,
    pandora_candidate_runs,
)


class ChargeReversalRecoveryTests(unittest.TestCase):
    def test_af3_batch_is_exact_30_job_factorial(self) -> None:
        jobs = make_af3_ablation_jobs()
        self.assertEqual(len(jobs), 30)
        conditions = {}
        for job in jobs:
            chains = [entry["proteinChain"] for entry in job["sequences"]]
            self.assertEqual(len(chains), 3)
            self.assertEqual(len(chains[2]["sequence"]), 15)
            self.assertTrue(chains[0]["useStructureTemplate"])
            self.assertTrue(chains[1]["useStructureTemplate"])
            key = (chains[2]["useStructureTemplate"], job["modelSeeds"][0])
            conditions[key] = conditions.get(key, 0) + 1
        self.assertEqual(conditions, {(True, 104729): 10, (False, 314159): 10, (False, 104729): 10})

    def test_pandora_manifest_requests_exactly_1200_models(self) -> None:
        rows = pandora_candidate_runs()
        self.assertEqual(len(rows), 60)
        self.assertEqual(sum(row["requested_models"] for row in rows), 1200)
        self.assertEqual({row["forced_template_pdb"] for row in rows}, {"1BX2", "6CQQ"})
        self.assertEqual({row["register_label"] for row in rows}, {"minus_1", "expected", "plus_1"})

    def test_calibration_is_reciprocal_leave_one_out(self) -> None:
        rows = pandora_calibration_runs()
        self.assertEqual(len(rows), 6)
        self.assertEqual(sum(row["requested_models"] for row in rows), 120)
        self.assertTrue(all(row["forced_template_pdb"] != row["excluded_template_pdb"] for row in rows))
        expected = {row["experimental_target_pdb"]: int(row["core_start_1_based"]) for row in rows if row["register_label"] == "expected"}
        self.assertEqual(expected, {"1BX2": 5, "6CQQ": 3})

    def test_consensus_is_fail_closed_and_ordered(self) -> None:
        base = {
            "complete": True,
            "calibration_pass": True,
            "sequence_qc_pass": True,
            "structural_qc_pass": True,
            "expected_register_preferred": True,
            "expected_register_core_rmsd_max_A": 1.9,
            "af3_seeds_agree": True,
            "pandora_templates_agree": True,
            "methods_agree": True,
            "af3_contact_frequency": 0.8,
            "pandora_contact_frequency": 0.9,
        }
        self.assertEqual(classify_consensus(base), "robust_structural_contact_hypothesis")
        self.assertEqual(classify_consensus({**base, "pandora_contact_frequency": 0.79}), "no_structural_charge_contact_support")
        self.assertEqual(classify_consensus({**base, "methods_agree": False}), "method_or_template_dependent")
        self.assertEqual(classify_consensus({**base, "expected_register_core_rmsd_max_A": 2.01}), "register_confounded")
        self.assertEqual(classify_consensus({**base, "calibration_pass": False}), "not_evaluable")
        self.assertEqual(classify_consensus({}), "not_evaluable")


if __name__ == "__main__":
    unittest.main()
