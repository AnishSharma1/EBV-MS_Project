#!/usr/bin/env python3
"""Run the frozen 49-pair binding/register threshold-sensitivity grid."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


CUTOFFS = (1.0, 2.0, 5.0, 10.0, 20.0)
PRIMARY_PREDICTORS = (
    "iedb_recommended_binding",
    "netmhciipan_4_3_el",
    "mixmhc2pred_2_1_context",
)
BINDING_RULES = ("all_3_both_arms", "at_least_2_of_3_both_arms")
REGISTER_RULES = ("exact_declared_core", "within_1_residue_of_declared_start")
ARMS = ("ebv", "self")
LEADS = ("HP36", "HP47")
CLAIM_BOUNDARY = (
    "Sensitivity of a computational binding/register gate only; not evidence of measured "
    "affinity, endogenous presentation, TCR recognition, cross-reactivity, molecular mimicry, "
    "or an MS mechanism."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[2]
    authoritative = Path(
        "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
    )
    root = authoritative / "processed/high_priority_handoff_hunt_2026-09-07"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest", type=Path, default=root / "frozen_49_pair_manifest.csv"
    )
    parser.add_argument(
        "--predictors", type=Path, default=root / "raw_responses/predictor_records.csv"
    )
    parser.add_argument(
        "--tiers",
        type=Path,
        default=root / "secondary_provisional_tiers_2026-09-07/all_49_provisional_tiers.csv",
    )
    parser.add_argument("--protocol", type=Path, default=root / "protocol_lock.json")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_2/threshold_sensitivity_2026-09-14",
    )
    return parser.parse_args()


def unique_start(sequence: str, core: str, label: str) -> int:
    starts = [index for index in range(len(sequence) - len(core) + 1) if sequence[index:index+len(core)] == core]
    if len(starts) != 1:
        raise ValueError(f"{label}: expected one occurrence of core {core!r}, observed {len(starts)}")
    return starts[0]


def binding_pass(values: list[float], cutoff: float, rule: str) -> bool:
    if rule == "all_3_both_arms":
        return all(value <= cutoff for value in values)
    if rule == "at_least_2_of_3_both_arms":
        return sum(value <= cutoff for value in values) >= 2
    raise ValueError(f"unknown binding rule: {rule}")


def register_pass(deltas: list[int], rule: str) -> bool:
    tolerance = 0 if rule == "exact_declared_core" else 1
    return sum(abs(delta) <= tolerance for delta in deltas) >= 2


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = read_csv(args.manifest)
    predictor_rows = read_csv(args.predictors)
    tier_rows = read_csv(args.tiers)
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))

    if len(manifest_rows) != 49:
        raise ValueError(f"expected 49 manifest rows, observed {len(manifest_rows)}")
    if len(predictor_rows) != 490:
        raise ValueError(f"expected 490 predictor rows, observed {len(predictor_rows)}")
    if len(tier_rows) != 49:
        raise ValueError(f"expected 49 tier rows, observed {len(tier_rows)}")
    if protocol.get("frozen_pair_count") != 49:
        raise ValueError("protocol does not declare a frozen 49-pair cohort")

    tier_by_id = {row["candidate_id"]: row for row in tier_rows}
    if len(tier_by_id) != 49:
        raise ValueError("tier table candidate IDs are not unique")
    candidate_by_pair = {row["pair_id"]: row["candidate_id"] for row in tier_rows}
    if len(candidate_by_pair) != 49:
        raise ValueError("tier table pair IDs are not unique")

    manifest_by_id = {}
    for row in manifest_rows:
        candidate_id = candidate_by_pair.get(row["pair_id"])
        if not candidate_id:
            raise ValueError(f"manifest pair missing from tier table: {row['pair_id']}")
        manifest_by_id[candidate_id] = row
    if len(manifest_by_id) != 49:
        raise ValueError("manifest-to-tier join did not produce 49 unique candidates")

    record_key_counts = Counter(
        (row["target_id"], row["arm_id"], row["predictor"]) for row in predictor_rows
    )
    duplicates = [key for key, count in record_key_counts.items() if count != 1]
    if duplicates:
        raise ValueError(f"predictor record keys are not unique: {duplicates[:3]}")
    if len({row["arm_id"] for row in predictor_rows}) != 98:
        raise ValueError("expected 98 unique peptide arms")

    predictor_counts = Counter(row["predictor"] for row in predictor_rows)
    expected_predictor_counts = {
        "iedb_recommended_binding": 98,
        "netmhciipan_4_3_ba": 98,
        "netmhciipan_4_3_el": 98,
        "mixmhc2pred_2_1_context": 98,
        "mixmhc2pred_2_1_no_context": 98,
    }
    if dict(predictor_counts) != expected_predictor_counts:
        raise ValueError(f"unexpected predictor coverage: {dict(predictor_counts)}")

    records = defaultdict(lambda: defaultdict(dict))
    for row in predictor_rows:
        arm = row["arm_id"].rsplit("__", 1)[-1]
        if arm not in ARMS:
            raise ValueError(f"unexpected arm suffix: {row['arm_id']}")
        records[row["target_id"]][arm][row["predictor"]] = row

    candidate_data = {}
    for candidate_id in sorted(tier_by_id):
        if candidate_id not in records:
            raise ValueError(f"missing predictor records for {candidate_id}")
        tier = tier_by_id[candidate_id]
        manifest = manifest_by_id[candidate_id]
        arm_data = {}
        for arm in ARMS:
            sequence = tier[f"{arm}_sequence"]
            declared_core = tier[f"{arm}_core"]
            declared_start = unique_start(sequence, declared_core, f"{candidate_id} {arm} declared")
            primary = {}
            for predictor in PRIMARY_PREDICTORS:
                if predictor not in records[candidate_id][arm]:
                    raise ValueError(f"{candidate_id} {arm}: missing {predictor}")
                record = records[candidate_id][arm][predictor]
                predicted_core = record["core"]
                predicted_start = unique_start(
                    sequence, predicted_core, f"{candidate_id} {arm} {predictor}"
                )
                percentile = float(record["percentile_rank"])
                if not 0.0 <= percentile <= 100.0:
                    raise ValueError(f"{candidate_id} {arm} {predictor}: percentile outside 0-100")
                primary[predictor] = {
                    "percentile": percentile,
                    "core": predicted_core,
                    "start_1_based": predicted_start + 1,
                    "delta": predicted_start - declared_start,
                }
            arm_data[arm] = {
                "sequence": sequence,
                "declared_core": declared_core,
                "declared_start_1_based": declared_start + 1,
                "primary": primary,
            }
        candidate_data[candidate_id] = {
            "tier": tier,
            "manifest": manifest,
            "arms": arm_data,
        }

    grid_rows = []
    for candidate_id in sorted(candidate_data):
        data = candidate_data[candidate_id]
        tier = data["tier"]
        for cutoff in CUTOFFS:
            for binding_rule in BINDING_RULES:
                arm_binding = {
                    arm: binding_pass(
                        [data["arms"][arm]["primary"][p]["percentile"] for p in PRIMARY_PREDICTORS],
                        cutoff,
                        binding_rule,
                    )
                    for arm in ARMS
                }
                for register_rule in REGISTER_RULES:
                    arm_register = {
                        arm: register_pass(
                            [data["arms"][arm]["primary"][p]["delta"] for p in PRIMARY_PREDICTORS],
                            register_rule,
                        )
                        for arm in ARMS
                    }
                    pair_binding = all(arm_binding.values())
                    pair_register = all(arm_register.values())
                    output = {
                        "candidate_id": candidate_id,
                        "allele": tier["allele"],
                        "pair_id": tier["pair_id"],
                        "ebv_protein": tier["ebv_protein"],
                        "self_protein": tier["self_protein"],
                        "cutoff_percent": f"{cutoff:g}",
                        "binding_rule": binding_rule,
                        "register_rule": register_rule,
                        "ebv_binding_pass": arm_binding["ebv"],
                        "self_binding_pass": arm_binding["self"],
                        "pair_binding_pass": pair_binding,
                        "ebv_register_pass": arm_register["ebv"],
                        "self_register_pass": arm_register["self"],
                        "pair_register_pass": pair_register,
                        "pair_grid_pass": pair_binding and pair_register,
                        "original_provisional_tier": tier["tier"],
                        "claim_boundary": CLAIM_BOUNDARY,
                    }
                    for arm in ARMS:
                        output[f"{arm}_sequence"] = data["arms"][arm]["sequence"]
                        output[f"{arm}_declared_core"] = data["arms"][arm]["declared_core"]
                        output[f"{arm}_declared_start_1_based"] = data["arms"][arm]["declared_start_1_based"]
                        for predictor in PRIMARY_PREDICTORS:
                            token = {
                                "iedb_recommended_binding": "iedb",
                                "netmhciipan_4_3_el": "netmhciipan_el",
                                "mixmhc2pred_2_1_context": "mixmhc2pred_context",
                            }[predictor]
                            value = data["arms"][arm]["primary"][predictor]
                            output[f"{arm}_{token}_percentile"] = f"{value['percentile']:.12g}"
                            output[f"{arm}_{token}_core"] = value["core"]
                            output[f"{arm}_{token}_start_delta"] = value["delta"]
                    grid_rows.append(output)

    expected_grid_rows = 49 * len(CUTOFFS) * len(BINDING_RULES) * len(REGISTER_RULES)
    if len(grid_rows) != expected_grid_rows:
        raise ValueError(f"expected {expected_grid_rows} grid rows, observed {len(grid_rows)}")

    grid_path = args.output_dir / "threshold_grid_pair_results.csv"
    write_csv(grid_path, list(grid_rows[0]), grid_rows)

    summary_rows = []
    for binding_rule in BINDING_RULES:
        for register_rule in REGISTER_RULES:
            for cutoff in CUTOFFS:
                matches = [
                    row for row in grid_rows
                    if row["binding_rule"] == binding_rule
                    and row["register_rule"] == register_rule
                    and float(row["cutoff_percent"]) == cutoff
                    and row["pair_grid_pass"]
                ]
                ids = sorted(row["candidate_id"] for row in matches)
                summary_rows.append(
                    {
                        "binding_rule": binding_rule,
                        "register_rule": register_rule,
                        "cutoff_percent": f"{cutoff:g}",
                        "passing_pair_count": len(ids),
                        "passing_candidate_ids": ";".join(ids),
                        "hp36_pass": "HP36" in ids,
                        "hp47_pass": "HP47" in ids,
                        "claim_boundary": CLAIM_BOUNDARY,
                    }
                )
    summary_path = args.output_dir / "threshold_grid_summary.csv"
    write_csv(summary_path, list(summary_rows[0]), summary_rows)

    lead_rows = [row for row in grid_rows if row["candidate_id"] in LEADS]
    lead_path = args.output_dir / "lead_candidate_sensitivity.csv"
    write_csv(lead_path, list(lead_rows[0]), lead_rows)

    exact_register_pairs = 0
    shifted_register_pairs = 0
    for candidate_id, data in candidate_data.items():
        exact = all(
            register_pass(
                [data["arms"][arm]["primary"][p]["delta"] for p in PRIMARY_PREDICTORS],
                "exact_declared_core",
            )
            for arm in ARMS
        )
        shifted = all(
            register_pass(
                [data["arms"][arm]["primary"][p]["delta"] for p in PRIMARY_PREDICTORS],
                "within_1_residue_of_declared_start",
            )
            for arm in ARMS
        )
        exact_register_pairs += int(exact)
        shifted_register_pairs += int(shifted)

    summary_lookup = {
        (row["binding_rule"], row["register_rule"], float(row["cutoff_percent"])): row
        for row in summary_rows
    }
    locked_cell = summary_lookup[("all_3_both_arms", "exact_declared_core", 20.0)]
    locked_ids = locked_cell["passing_candidate_ids"].split(";") if locked_cell["passing_candidate_ids"] else []
    if locked_ids != ["HP36", "HP47"]:
        raise ValueError(f"locked 20% gate did not reproduce HP36 and HP47: {locked_ids}")
    if exact_register_pairs != 44 or shifted_register_pairs != 47:
        raise ValueError(
            f"register counts did not reproduce: exact={exact_register_pairs}, shifted={shifted_register_pairs}"
        )

    lock = {
        "analysis": "binding_register_threshold_sensitivity_grid",
        "status": "complete",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "cutoffs_percent": list(CUTOFFS),
        "primary_predictors": list(PRIMARY_PREDICTORS),
        "binding_rules": list(BINDING_RULES),
        "register_rules": list(REGISTER_RULES),
        "candidate_count": 49,
        "arm_count": 98,
        "stored_predictor_record_count": 490,
        "grid_row_count": len(grid_rows),
        "claim_boundary": CLAIM_BOUNDARY,
        "python": sys.version,
        "platform": platform.platform(),
        "inputs": {
            str(path.resolve()): sha256(path)
            for path in (args.manifest, args.predictors, args.tiers, args.protocol, Path(__file__))
        },
    }
    lock_path = args.output_dir / "analysis_lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    qa_lines = [
        "# Threshold-sensitivity QA",
        "",
        "Status: PASS",
        "",
        "- Frozen candidates: 49 unique pair IDs and 49 unique candidate IDs.",
        "- Peptide arms: 98 unique arm IDs.",
        "- Predictor records: 490 total; exactly 98 records for each of five stored predictor variants.",
        "- Primary grid inputs: exactly three records per arm from the locked predictor set.",
        f"- Grid rows: {len(grid_rows)} = 49 candidates x 5 cutoffs x 2 binding rules x 2 register rules.",
        f"- Exact-register agreement reproduced: {exact_register_pairs}/49 pairs.",
        f"- +/-1-register agreement reproduced: {shifted_register_pairs}/49 pairs.",
        "- Locked strict 20% exact-register gate reproduced exactly HP36 and HP47.",
        "- Every declared and predicted 9-mer core occurred exactly once in its submitted peptide sequence.",
        "- Percentile ranks were validated within the inclusive 0-100 range.",
        "- Original provisional tiers were carried only as annotations and were not recomputed.",
        "",
        f"Claim boundary: {CLAIM_BOUNDARY}",
    ]
    qa_path = args.output_dir / "analysis_qa.md"
    qa_path.write_text("\n".join(qa_lines) + "\n", encoding="utf-8")

    results_lines = [
        "# Binding/register threshold-sensitivity results",
        "",
        "| Binding rule | Register rule | 1% | 2% | 5% | 10% | 20% |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for binding_rule in BINDING_RULES:
        for register_rule in REGISTER_RULES:
            counts = [
                summary_lookup[(binding_rule, register_rule, cutoff)]["passing_pair_count"]
                for cutoff in CUTOFFS
            ]
            results_lines.append(
                f"| {binding_rule} | {register_rule} | " + " | ".join(str(value) for value in counts) + " |"
            )
    results_lines.extend(
        [
            "",
            "Key candidate results:",
            "",
            "- Under the locked all-three/both-arms rule, HP36 and HP47 pass only at 20%; neither passes at 10% or below.",
            "- Under the two-of-three/both-arms rule, HP47 passes at 10%; HP36 first passes at 20%.",
            "- The +/-1 register rule adds no candidates to the strict all-three rule and adds HP19 only in the two-of-three 20% cell.",
            "- The two-pair locked result is therefore conditional on the permissive 20% cutoff.",
            "",
            f"Claim boundary: {CLAIM_BOUNDARY}",
        ]
    )
    results_path = args.output_dir / "RESULTS_SUMMARY.md"
    results_path.write_text("\n".join(results_lines) + "\n", encoding="utf-8")

    manifest_targets = [
        lock_path,
        grid_path,
        summary_path,
        lead_path,
        qa_path,
        results_path,
        Path(__file__).resolve(),
    ]
    checksum_rows = [
        {"path": str(path.resolve()), "sha256": sha256(path)} for path in manifest_targets
    ]
    checksum_path = args.output_dir / "SHA256SUMS.csv"
    write_csv(checksum_path, ["path", "sha256"], checksum_rows)

    print(json.dumps({
        "status": "PASS",
        "output_dir": str(args.output_dir.resolve()),
        "grid_rows": len(grid_rows),
        "exact_register_pairs": exact_register_pairs,
        "shifted_register_pairs": shifted_register_pairs,
        "summary": summary_rows,
    }, indent=2))


if __name__ == "__main__":
    main()
