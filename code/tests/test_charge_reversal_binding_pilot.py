from __future__ import annotations

import unittest

from charge_reversal_binding_pilot import (
    evaluate_fitted_results,
    fitted_results_template,
    make_plate_map,
    primary_panel,
    validate_panel,
)


class ChargeReversalBindingPilotTests(unittest.TestCase):
    def test_frozen_panel_and_single_site_mutations(self) -> None:
        panel = primary_panel()
        validate_panel(panel)
        primary = [row for row in panel if row["primary_mutation_panel"]]
        self.assertEqual(len(panel), 14)
        self.assertEqual(len(primary), 10)
        self.assertEqual({row["allele"] for row in panel}, {"HLA-DRB1*15:01"})

    def test_plate_map_is_nine_complete_unique_plates(self) -> None:
        plate_map = make_plate_map(primary_panel())
        self.assertEqual(len(plate_map), 864)
        for experiment in range(1, 4):
            for plate in range(1, 4):
                rows = [row for row in plate_map if row["plate_id"] == f"E{experiment}_P{plate}"]
                self.assertEqual(len(rows), 96)
                self.assertEqual(len({row["well"] for row in rows}), 96)

    def test_fitted_template_has_three_experiments_for_each_primary_peptide(self) -> None:
        rows = fitted_results_template(primary_panel())
        self.assertEqual(len(rows), 30)
        self.assertEqual({row["experiment_id"] for row in rows}, {1, 2, 3})
        self.assertTrue(all(row["result_status"] == "not_run" for row in rows))

    def test_gate_passes_only_stronger_reversal_with_full_qc(self) -> None:
        panel = primary_panel()
        primary = [row for row in panel if row["primary_mutation_panel"]]
        values = {
            "BALF5_WT_15": 100.0,
            "BALF5_P6_KQ": 150.0,
            "BALF5_P6_KE": 400.0,
            "BALF5_P7_KQ": 120.0,
            "BALF5_P7_KE": 150.0,
            "TALDO1_WT_15": 200.0,
            "TALDO1_P6_KQ": 250.0,
            "TALDO1_P6_KE": 300.0,
            "TALDO1_P7_KQ": 300.0,
            "TALDO1_P7_KE": 900.0,
        }
        rows = []
        for experiment in range(1, 4):
            for meta in primary:
                rows.append({
                    "experiment_id": str(experiment),
                    "sample_id": meta["sample_id"],
                    "ic50_nM": str(values[meta["sample_id"]]),
                    "kd_nM_optional": "",
                    "full_curve_fit_ok": "true",
                    "measured_or_supported_core": meta["expected_core"],
                    "register_check_pass": "true",
                    "purity_percent": "96",
                    "solubility_pass": "true",
                    "analyst_blinded": "true",
                    "exclusion_reason": "",
                    "result_status": "measured",
                })
        _, gates = evaluate_fitted_results(rows, panel)
        status = {(row["arm"], int(row["core_position"])): row["gate_status"] for row in gates}
        self.assertEqual(status[("BALF5", 6)], "pass_charge_sensitive_binding")
        self.assertEqual(status[("BALF5", 7)], "does_not_pass_charge_sensitive_binding")
        self.assertEqual(status[("TALDO1", 6)], "does_not_pass_charge_sensitive_binding")
        self.assertEqual(status[("TALDO1", 7)], "pass_charge_sensitive_binding")

    def test_incomplete_results_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly one row"):
            evaluate_fitted_results([], primary_panel())


if __name__ == "__main__":
    unittest.main()
