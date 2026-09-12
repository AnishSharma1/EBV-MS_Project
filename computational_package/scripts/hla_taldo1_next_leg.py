"""Pure logic for the TALDO1 HLA-specificity and charge-sensitivity analysis."""

from __future__ import annotations

import csv
import io
from typing import Any


ANCHOR_POSITIONS = {1, 4, 6, 9}
CHARGE_REVERSAL = {"K": "E", "R": "E", "D": "K", "E": "K"}
NEUTRALIZATION = {"K": "Q", "R": "Q", "D": "N", "E": "Q"}


def locate_unique_core(sequence: str, core: str) -> int:
    """Return the one-based start of one unique core within a peptide."""
    sequence = sequence.strip().upper()
    core = core.strip().upper()
    starts = [index + 1 for index in range(len(sequence)) if sequence.startswith(core, index)]
    if len(starts) != 1:
        raise ValueError(f"expected one unique core occurrence, found {len(starts)}")
    return starts[0]


def generate_charge_controls(
    *, peptide_id: str, sequence: str, predicted_core: str
) -> list[dict[str, Any]]:
    """Generate single-site formal-charge reversal and neutralization controls."""
    sequence = sequence.strip().upper()
    predicted_core = predicted_core.strip().upper()
    core_start = locate_unique_core(sequence, predicted_core)
    variants: list[dict[str, Any]] = []
    for core_index, residue in enumerate(predicted_core, start=1):
        if residue not in CHARGE_REVERSAL:
            continue
        sequence_index = core_start + core_index - 2
        for perturbation, mapping in (
            ("charge_reversal", CHARGE_REVERSAL),
            ("neutralization", NEUTRALIZATION),
        ):
            mutant_residue = mapping[residue]
            mutant_sequence = sequence[:sequence_index] + mutant_residue + sequence[sequence_index + 1 :]
            variants.append(
                {
                    "variant_id": f"{peptide_id}__P{core_index}_{residue}{mutant_residue}__{perturbation}",
                    "peptide_id": peptide_id,
                    "wild_type_sequence": sequence,
                    "baseline_core": predicted_core,
                    "core_position": core_index,
                    "peptide_position": sequence_index + 1,
                    "pocket_role": "anchor" if core_index in ANCHOR_POSITIONS else "non_anchor_core",
                    "original_residue": residue,
                    "mutant_residue": mutant_residue,
                    "perturbation": perturbation,
                    "mutant_sequence": mutant_sequence,
                }
            )
    return variants


def classify_mutation_result(
    *, expected_mutant_core: str, mutant_core: str, pocket_role: str, perturbation: str
) -> str:
    """Label what a predictor sensitivity result can and cannot isolate."""
    if expected_mutant_core.strip().upper() != mutant_core.strip().upper():
        return "confounded_register_shift"
    if perturbation == "neutralization":
        return f"{pocket_role}_neutralization_sensitivity_same_register"
    if pocket_role == "anchor":
        return "anchor_charge_sensitivity_same_register"
    return "non_anchor_charge_sensitivity_same_register"


def parse_iedb_tsv(
    allele: str, submissions: list[dict[str, Any]], raw_text: str
) -> list[dict[str, Any]]:
    """Map one IEDB response back to frozen submissions by sequence number."""
    rows = list(csv.DictReader(io.StringIO(raw_text), delimiter="\t"))
    required = {"allele", "seq_num", "core_peptide", "peptide", "ic50", "rank"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("IEDB response is missing required columns")
    by_seq: dict[int, dict[str, str]] = {}
    for row in rows:
        seq_num = int(row["seq_num"])
        if seq_num in by_seq:
            raise ValueError(f"duplicate IEDB seq_num: {seq_num}")
        by_seq[seq_num] = row
    if len(by_seq) != len(submissions):
        raise ValueError("IEDB response does not cover every submission exactly once")
    records: list[dict[str, Any]] = []
    for submission in submissions:
        seq_num = int(submission["seq_num"])
        row = by_seq.get(seq_num)
        if row is None:
            raise ValueError(f"missing IEDB seq_num: {seq_num}")
        if row["allele"] != allele:
            raise ValueError(f"IEDB allele mismatch for seq_num {seq_num}")
        if row["peptide"] != submission["sequence"]:
            raise ValueError(f"IEDB peptide mismatch for seq_num {seq_num}")
        records.append(
            {
                "allele": allele,
                "seq_num": seq_num,
                "peptide_id": submission["peptide_id"],
                "sequence": submission["sequence"],
                "predicted_core": row["core_peptide"],
                "ic50_nM": float(row["ic50"]),
                "rank_percentile": float(row["rank"]),
            }
        )
    return records


def summarize_mutation_result(
    *,
    allele: str,
    variant: dict[str, Any],
    baseline: dict[str, Any],
    mutant: dict[str, Any],
) -> dict[str, Any]:
    """Create a descriptive mutation-versus-wild-type result record."""
    baseline_core = str(baseline["predicted_core"])
    core_index = int(variant["core_position"]) - 1
    expected_mutant_core = (
        baseline_core[:core_index]
        + str(variant["mutant_residue"])
        + baseline_core[core_index + 1 :]
    )
    mutant_core = str(mutant["predicted_core"])
    baseline_rank = float(baseline["rank_percentile"])
    mutant_rank = float(mutant["rank_percentile"])
    baseline_ic50 = float(baseline["ic50_nM"])
    mutant_ic50 = float(mutant["ic50_nM"])
    interpretation = classify_mutation_result(
        expected_mutant_core=expected_mutant_core,
        mutant_core=mutant_core,
        pocket_role=str(variant["pocket_role"]),
        perturbation=str(variant["perturbation"]),
    )
    return {
        "allele": allele,
        **variant,
        "expected_mutant_core": expected_mutant_core,
        "mutant_predicted_core": mutant_core,
        "baseline_rank_percentile": baseline_rank,
        "mutant_rank_percentile": mutant_rank,
        "delta_rank_percentile": round(mutant_rank - baseline_rank, 12),
        "baseline_ic50_nM": baseline_ic50,
        "mutant_ic50_nM": mutant_ic50,
        "ic50_fold_change": round(mutant_ic50 / baseline_ic50, 12),
        "register_shifted": mutant_core != expected_mutant_core,
        "interpretation_status": interpretation,
        "claim_boundary": (
            "Predictor sensitivity only; not evidence of physical repulsion, natural "
            "presentation, T-cell recognition, cross-reactivity, or MS mechanism."
        ),
    }


def pair_charge_controls(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Pair reversal and neutralization results at the same allele and core position."""
    grouped: dict[tuple[str, str, int], dict[str, dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["allele"]), str(row["peptide_id"]), int(row["core_position"]))
        grouped.setdefault(key, {})[str(row["perturbation"])] = row
    paired: list[dict[str, Any]] = []
    for (allele, peptide_id, core_position), controls in sorted(grouped.items()):
        reversal = controls.get("charge_reversal")
        neutral = controls.get("neutralization")
        if reversal is None or neutral is None:
            continue
        same_register = not bool(reversal["register_shifted"]) and not bool(neutral["register_shifted"])
        paired.append(
            {
                "allele": allele,
                "peptide_id": peptide_id,
                "core_position": core_position,
                "pocket_role": reversal["pocket_role"],
                "original_residue": reversal["original_residue"],
                "reversal_residue": reversal["mutant_residue"],
                "neutral_residue": neutral["mutant_residue"],
                "reversal_delta_rank": float(reversal["delta_rank_percentile"]),
                "neutralization_delta_rank": float(neutral["delta_rank_percentile"]),
                "reversal_minus_neutralization_delta_rank": round(
                    float(reversal["delta_rank_percentile"])
                    - float(neutral["delta_rank_percentile"]),
                    12,
                ),
                "reversal_ic50_fold": float(reversal["ic50_fold_change"]),
                "neutralization_ic50_fold": float(neutral["ic50_fold_change"]),
                "comparison_status": "same_register_pair" if same_register else "confounded_register_shift",
                "claim_boundary": (
                    "Descriptive predictor contrast only; not evidence of physical repulsion, "
                    "energetic mechanism, presentation, or T-cell recognition."
                ),
            }
        )
    return paired
