"""Core logic for the HY15 DRB1*15:01 charge-reversal binding pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ALLELE = "HLA-DRB1*15:01"
TOOL_ALLELE = "DRB1_15_01"
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")
CLAIM_BOUNDARY = (
    "Charge-reversal peptide-HLA binding pilot only; not evidence of natural "
    "presentation, TCR recognition, cross-reactivity, molecular mimicry, or MS mechanism."
)


def primary_panel() -> list[dict[str, Any]]:
    """Return the frozen peptide and assay-control manifest."""
    rows = [
        ("BALF5_WT_15", "BALF5", "wild_type", "TGGVYHFVKKHVHES", "VYHFVKKHV", "", "", "", True, True),
        ("BALF5_P6_KQ", "BALF5", "neutralization", "TGGVYHFVQKHVHES", "VYHFVQKHV", 6, 9, "K>Q", True, True),
        ("BALF5_P6_KE", "BALF5", "charge_reversal", "TGGVYHFVEKHVHES", "VYHFVEKHV", 6, 9, "K>E", True, True),
        ("BALF5_P7_KQ", "BALF5", "neutralization", "TGGVYHFVKQHVHES", "VYHFVKQHV", 7, 10, "K>Q", True, True),
        ("BALF5_P7_KE", "BALF5", "charge_reversal", "TGGVYHFVKEHVHES", "VYHFVKEHV", 7, 10, "K>E", True, True),
        ("TALDO1_WT_15", "TALDO1", "wild_type", "SVTKIYNYYKKFSYK", "IYNYYKKFS", "", "", "", True, True),
        ("TALDO1_P6_KQ", "TALDO1", "neutralization", "SVTKIYNYYQKFSYK", "IYNYYQKFS", 6, 10, "K>Q", True, True),
        ("TALDO1_P6_KE", "TALDO1", "charge_reversal", "SVTKIYNYYEKFSYK", "IYNYYEKFS", 6, 10, "K>E", True, True),
        ("TALDO1_P7_KQ", "TALDO1", "neutralization", "SVTKIYNYYKQFSYK", "IYNYYKQFS", 7, 11, "K>Q", True, True),
        ("TALDO1_P7_KE", "TALDO1", "charge_reversal", "SVTKIYNYYKEFSYK", "IYNYYKEFS", 7, 11, "K>E", True, True),
        ("BALF5_REGISTER_11", "BALF5", "nested_register_control", "GVYHFVKKHVH", "VYHFVKKHV", "", "", "", False, True),
        ("TALDO1_REGISTER_11", "TALDO1", "nested_register_control", "KIYNYYKKFSY", "IYNYYKKFS", "", "", "", False, True),
        ("MBP_DR15_POSITIVE", "control", "known_positive_control", "ENPVVHFFKNIVTPR", "VHFFKNIVT", "", "", "", False, True),
        ("HUMAN_BG_120112_WEAK", "control", "predicted_weak_control", "NQELRADGTVNQIEG", "RADGTVNQI", "", "", "", False, True),
    ]
    fields = (
        "sample_id", "arm", "reagent_role", "sequence", "expected_core",
        "core_position", "peptide_position_1_based", "substitution",
        "primary_mutation_panel", "include_in_binding_assay",
    )
    return [
        {
            **dict(zip(fields, row)),
            "allele": ALLELE,
            "length": len(row[3]),
            "proposed_not_ordered": True,
            "claim_boundary": CLAIM_BOUNDARY,
        }
        for row in rows
    ]


def validate_panel(rows: list[dict[str, Any]]) -> None:
    if len(rows) != 14 or len({row["sample_id"] for row in rows}) != 14:
        raise ValueError("panel must contain 14 uniquely named reagents")
    for row in rows:
        sequence = str(row["sequence"])
        core = str(row["expected_core"])
        if not sequence or set(sequence) - ALLOWED_AA:
            raise ValueError(f"invalid peptide sequence for {row['sample_id']}")
        if core not in sequence:
            raise ValueError(f"expected core is absent for {row['sample_id']}")
    primary = [row for row in rows if row["primary_mutation_panel"]]
    if len(primary) != 10:
        raise ValueError("primary mutation panel must contain exactly 10 peptides")
    by_arm = defaultdict(list)
    for row in primary:
        by_arm[row["arm"]].append(row)
    if {key: len(value) for key, value in by_arm.items()} != {"BALF5": 5, "TALDO1": 5}:
        raise ValueError("each biological arm must contain WT plus four variants")
    for arm, arm_rows in by_arm.items():
        wt = next(row for row in arm_rows if row["reagent_role"] == "wild_type")
        for row in arm_rows:
            if row is wt:
                continue
            diffs = [index for index, pair in enumerate(zip(wt["sequence"], row["sequence"]), 1) if pair[0] != pair[1]]
            if diffs != [int(row["peptide_position_1_based"])]:
                raise ValueError(f"{row['sample_id']} is not the declared single-site mutant")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_plate_map(rows: list[dict[str, Any]], seed: int = 15012026) -> list[dict[str, Any]]:
    """Create three full 96-well plates for each of three independent experiments."""
    assay_rows = [row for row in rows if row["include_in_binding_assay"]]
    concentrations = [100000.0 / (3 ** index) for index in range(10)]
    wells = [f"{letter}{column}" for letter in "ABCDEFGH" for column in range(1, 13)]
    output: list[dict[str, Any]] = []
    for experiment in range(1, 4):
        pairs = [(row, concentration) for row in assay_rows for concentration in concentrations]
        random.Random(seed + experiment).shuffle(pairs)
        by_plate: dict[int, list[dict[str, Any]]] = {1: [], 2: [], 3: []}
        for index, (row, concentration) in enumerate(pairs):
            plate_a = index % 3 + 1
            plate_b = plate_a % 3 + 1
            for technical_replicate, plate in ((1, plate_a), (2, plate_b)):
                by_plate[plate].append({
                    "experiment_id": experiment,
                    "plate_id": f"E{experiment}_P{plate}",
                    "sample_id": row["sample_id"],
                    "assay_role": row["reagent_role"],
                    "competitor_concentration_nM": round(concentration, 6),
                    "technical_replicate": technical_replicate,
                    "observed_signal": "",
                    "run_status": "not_run",
                })
        control_slots = {
            1: ["no_competitor_max_signal", "no_competitor_max_signal", "no_hla_background"],
            2: ["no_competitor_max_signal", "no_hla_background"],
            3: ["no_competitor_max_signal", "no_hla_background", "no_hla_background"],
        }
        for plate, controls in control_slots.items():
            for control in controls:
                by_plate[plate].append({
                    "experiment_id": experiment,
                    "plate_id": f"E{experiment}_P{plate}",
                    "sample_id": control,
                    "assay_role": "assay_control",
                    "competitor_concentration_nM": "",
                    "technical_replicate": "",
                    "observed_signal": "",
                    "run_status": "not_run",
                })
        for plate, plate_rows in by_plate.items():
            if len(plate_rows) != 96:
                raise ValueError(f"plate {experiment}/{plate} has {len(plate_rows)} wells")
            random.Random(seed + experiment * 10 + plate).shuffle(plate_rows)
            for well, row in zip(wells, plate_rows):
                output.append({**row, "well": well})
    return output


def fitted_results_template(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for experiment in range(1, 4):
        for row in rows:
            if not row["primary_mutation_panel"]:
                continue
            output.append({
                "experiment_id": experiment,
                "sample_id": row["sample_id"],
                "ic50_nM": "",
                "kd_nM_optional": "",
                "full_curve_fit_ok": "",
                "measured_or_supported_core": "",
                "register_check_pass": "",
                "purity_percent": "",
                "solubility_pass": "",
                "analyst_blinded": "",
                "exclusion_reason": "",
                "result_status": "not_run",
            })
    return output


def _truth(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "pass"}


def _geomean(values: Iterable[float]) -> float:
    values = list(values)
    return math.exp(sum(math.log(value) for value in values) / len(values))


def evaluate_fitted_results(
    result_rows: list[dict[str, str]], manifest_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Apply the frozen three-experiment decision rules to completed IC50 results."""
    primary = {row["sample_id"]: row for row in manifest_rows if row["primary_mutation_panel"]}
    expected = {(experiment, sample) for experiment in range(1, 4) for sample in primary}
    observed = {(int(row["experiment_id"]), row["sample_id"]) for row in result_rows}
    if observed != expected or len(result_rows) != len(expected):
        raise ValueError("results must contain exactly one row per primary peptide and experiment")
    parsed: dict[tuple[int, str], dict[str, Any]] = {}
    for row in result_rows:
        key = (int(row["experiment_id"]), row["sample_id"])
        try:
            ic50 = float(row["ic50_nM"])
        except (TypeError, ValueError):
            raise ValueError(f"missing or invalid IC50 for {key}")
        if not math.isfinite(ic50) or ic50 <= 0:
            raise ValueError(f"IC50 must be positive and finite for {key}")
        parsed[key] = {
            **row,
            "ic50_nM": ic50,
            "qc_pass": all(
                (
                    _truth(row["full_curve_fit_ok"]),
                    _truth(row["register_check_pass"]),
                    _truth(row["solubility_pass"]),
                    _truth(row["analyst_blinded"]),
                    float(row["purity_percent"]) >= 95.0,
                    row["measured_or_supported_core"] == str(primary[row["sample_id"]]["expected_core"]),
                    not row["exclusion_reason"].strip(),
                )
            ),
        }
    wt_for_arm = {"BALF5": "BALF5_WT_15", "TALDO1": "TALDO1_WT_15"}
    summaries: list[dict[str, Any]] = []
    for sample_id, meta in primary.items():
        ratios = []
        qc = []
        for experiment in range(1, 4):
            current = parsed[(experiment, sample_id)]
            wt = parsed[(experiment, wt_for_arm[meta["arm"]])]
            ratios.append(current["ic50_nM"] / wt["ic50_nM"])
            qc.extend((current["qc_pass"], wt["qc_pass"]))
        summaries.append({
            "sample_id": sample_id,
            "arm": meta["arm"],
            "reagent_role": meta["reagent_role"],
            "core_position": meta["core_position"],
            "geomean_ic50_fold_vs_wt": round(_geomean(ratios), 6),
            "direction_weaker_in_all_three": all(value > 1.0 for value in ratios),
            "all_qc_and_register_checks_pass": all(qc),
            "folds_by_experiment": ";".join(f"{value:.6g}" for value in ratios),
        })
    summary_lookup = {row["sample_id"]: row for row in summaries}
    gates: list[dict[str, Any]] = []
    for sample_id, meta in primary.items():
        if meta["reagent_role"] != "charge_reversal":
            continue
        neutral_id = sample_id[:-2] + "KQ"
        reversal = summary_lookup[sample_id]
        neutral = summary_lookup[neutral_id]
        rev_over_neutral = reversal["geomean_ic50_fold_vs_wt"] / neutral["geomean_ic50_fold_vs_wt"]
        passed = all((
            reversal["geomean_ic50_fold_vs_wt"] >= 3.0,
            rev_over_neutral >= 2.0,
            reversal["direction_weaker_in_all_three"],
            reversal["all_qc_and_register_checks_pass"],
            neutral["all_qc_and_register_checks_pass"],
        ))
        gates.append({
            "arm": meta["arm"],
            "core_position": meta["core_position"],
            "reversal_sample_id": sample_id,
            "neutralization_sample_id": neutral_id,
            "reversal_fold_vs_wt": reversal["geomean_ic50_fold_vs_wt"],
            "neutralization_fold_vs_wt": neutral["geomean_ic50_fold_vs_wt"],
            "reversal_over_neutralization_fold": round(rev_over_neutral, 6),
            "gate_status": "pass_charge_sensitive_binding" if passed else "does_not_pass_charge_sensitive_binding",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return summaries, gates


def checksum_rows(root: Path, checksum_name: str = "SHA256SUMS.csv") -> list[dict[str, Any]]:
    return [
        {
            "relative_path": str(path.relative_to(root)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != checksum_name
    ]
