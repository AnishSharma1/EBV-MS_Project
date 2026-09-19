"""Contract tests for the TALDO1 HLA-specificity next leg."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hla_taldo1_next_leg import (  # noqa: E402
    classify_mutation_result,
    generate_charge_controls,
    locate_unique_core,
    pair_charge_controls,
    parse_iedb_tsv,
    summarize_mutation_result,
)
from run_hla_taldo1_next_leg import PEPTIDES  # noqa: E402


class HlaTaldo1NextLegTests(unittest.TestCase):
    def test_frozen_candidate_sequences_match_stage1_source(self):
        observed = {row["peptide_id"]: row["sequence"] for row in PEPTIDES}
        self.assertEqual(observed["TALDO1_108_122"], "DARLSFDKDAMVARA")
        self.assertEqual(observed["TALDO1_216_230"], "SVTKIYNYYKKFSYK")

    def test_locate_unique_core_returns_one_based_position(self):
        self.assertEqual(locate_unique_core("TGGVYHFVKKHVHES", "YHFVKKHVH"), 5)

    def test_locate_unique_core_rejects_missing_or_repeated_core(self):
        with self.assertRaises(ValueError):
            locate_unique_core("AAAAAAAAAA", "AAAAAAAAA")
        with self.assertRaises(ValueError):
            locate_unique_core("ACDEFGHIKLM", "YYYYYYYYY")

    def test_charge_controls_create_reversal_and_neutralization_per_formal_charge(self):
        variants = generate_charge_controls(
            peptide_id="BALF5",
            sequence="TGGVYHFVKKHVHES",
            predicted_core="YHFVKKHVH",
        )
        observed = {
            (row["core_position"], row["original_residue"], row["mutant_residue"], row["perturbation"])
            for row in variants
        }
        self.assertEqual(
            observed,
            {
                (5, "K", "E", "charge_reversal"),
                (5, "K", "Q", "neutralization"),
                (6, "K", "E", "charge_reversal"),
                (6, "K", "Q", "neutralization"),
            },
        )
        self.assertTrue(all(len(row["mutant_sequence"]) == 15 for row in variants))
        self.assertEqual(
            {row["core_position"]: row["pocket_role"] for row in variants},
            {5: "non_anchor_core", 6: "anchor"},
        )

    def test_charge_controls_label_anchor_positions(self):
        variants = generate_charge_controls(
            peptide_id="TALDO1_108_122",
            sequence="RVLSFDKDAMVARAR",
            predicted_core="LSFDKDAMV",
        )
        p4 = [row for row in variants if row["core_position"] == 4]
        p6 = [row for row in variants if row["core_position"] == 6]
        self.assertTrue(p4 and p6)
        self.assertTrue(all(row["pocket_role"] == "anchor" for row in p4 + p6))

    def test_result_classification_marks_register_shift_as_confounding(self):
        result = classify_mutation_result(
            expected_mutant_core="YHFVEKHVH",
            mutant_core="HFVEKHVHE",
            pocket_role="non_anchor_core",
            perturbation="charge_reversal",
        )
        self.assertEqual(result, "confounded_register_shift")

    def test_result_classification_never_calls_prediction_repulsion(self):
        self.assertEqual(
            classify_mutation_result(
                expected_mutant_core="YHFVEKHVH",
                mutant_core="YHFVEKHVH",
                pocket_role="non_anchor_core",
                perturbation="charge_reversal",
            ),
            "non_anchor_charge_sensitivity_same_register",
        )
        self.assertEqual(
            classify_mutation_result(
                expected_mutant_core="LSFKKDAMV",
                mutant_core="LSFKKDAMV",
                pocket_role="anchor",
                perturbation="charge_reversal",
            ),
            "anchor_charge_sensitivity_same_register",
        )

    def test_iedb_parser_maps_out_of_order_rows_by_seq_num(self):
        submissions = [
            {"seq_num": 1, "peptide_id": "BALF5", "sequence": "TGGVYHFVKKHVHES"},
            {"seq_num": 2, "peptide_id": "TALDO1", "sequence": "KELIYNYYKKFSYVI"},
        ]
        raw = (
            "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\trank\n"
            "HLA-DRB5*01:01\t2\t1\t15\t15\tYNYYKKFSY\tKELIYNYYKKFSYVI\t159.07\t13.0\n"
            "HLA-DRB5*01:01\t1\t1\t15\t15\tYHFVKKHVH\tTGGVYHFVKKHVHES\t126.27\t11.0\n"
        )
        records = parse_iedb_tsv("HLA-DRB5*01:01", submissions, raw)
        self.assertEqual([row["peptide_id"] for row in records], ["BALF5", "TALDO1"])
        self.assertEqual(records[0]["predicted_core"], "YHFVKKHVH")

    def test_iedb_parser_rejects_sequence_mismatch(self):
        submissions = [{"seq_num": 1, "peptide_id": "BALF5", "sequence": "TGGVYHFVKKHVHES"}]
        raw = (
            "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\trank\n"
            "HLA-DRB5*01:01\t1\t1\t15\t15\tYHFVKKHVH\tAAAAAAAAAAAAAAA\t126.27\t11.0\n"
        )
        with self.assertRaises(ValueError):
            parse_iedb_tsv("HLA-DRB5*01:01", submissions, raw)

    def test_mutation_summary_compares_against_expected_mutated_core(self):
        variant = generate_charge_controls(
            peptide_id="BALF5",
            sequence="TGGVYHFVKKHVHES",
            predicted_core="YHFVKKHVH",
        )[0]
        row = summarize_mutation_result(
            allele="HLA-DRB5*01:01",
            variant=variant,
            baseline={"predicted_core": "YHFVKKHVH", "rank_percentile": 11.0, "ic50_nM": 126.27},
            mutant={"predicted_core": "YHFVEKHVH", "rank_percentile": 22.0, "ic50_nM": 252.54},
        )
        self.assertEqual(row["delta_rank_percentile"], 11.0)
        self.assertEqual(row["ic50_fold_change"], 2.0)
        self.assertFalse(row["register_shifted"])
        self.assertNotIn("repulsion", row["interpretation_status"])

    def test_pair_charge_controls_compares_same_position_controls(self):
        base = {
            "allele": "HLA-DRB5*01:01",
            "peptide_id": "BALF5",
            "core_position": 5,
            "original_residue": "K",
            "pocket_role": "non_anchor_core",
            "register_shifted": False,
        }
        rows = [
            {**base, "perturbation": "charge_reversal", "mutant_residue": "E", "delta_rank_percentile": 9.0, "ic50_fold_change": 2.15},
            {**base, "perturbation": "neutralization", "mutant_residue": "Q", "delta_rank_percentile": 1.0, "ic50_fold_change": 1.10},
        ]
        paired = pair_charge_controls(rows)
        self.assertEqual(len(paired), 1)
        self.assertEqual(paired[0]["reversal_minus_neutralization_delta_rank"], 8.0)
        self.assertEqual(paired[0]["comparison_status"], "same_register_pair")


if __name__ == "__main__":
    unittest.main()
