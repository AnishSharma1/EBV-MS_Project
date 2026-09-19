#!/usr/bin/env python3
"""Build a claim-bounded computational-tightening package for manuscript V8.

This script does not edit manuscript prose. It reuses the frozen 6,400-pair
universe and frozen 49-pair predictor panel to produce auditable score
definitions, ablations, register-shift sensitivity, matched empirical decoys,
corrected structural-ensemble summaries, figures, checksums, and a prospective
pre-experiment lock.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import platform
import shlex
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


TCR_POSITIONS = (1, 2, 4, 6, 7)
ANCHOR_POSITIONS = (0, 3, 5, 8)
LEADS = ("HY13_SEQ_02", "HY15_SEQ_02")
CLAIM_BOUNDARY = (
    "Computational robustness and sensitivity analysis only; not evidence of measured "
    "binding, endogenous presentation, TCR recognition, cross-reactivity, molecular "
    "mimicry, or an MS mechanism."
)


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[2]
    authoritative = Path(
        "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
    )
    categorized = workspace / (
        "outputs/EBV_MS_ranking_tables_corrected_2026-09-03/categorized_tables"
    )
    focused = authoritative / "processed/taldo1_focused_register_2026-09-05"
    expanded = authoritative / "processed/high_priority_handoff_hunt_2026-09-07"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pairs",
        type=Path,
        default=categorized / "02_v3_structural_shortlist_and_rankings/04_all_6400_v3_pairs.csv",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=categorized / "01_current_stage1_evidence_review/03_presentation_conditioned_candidate_ranks.csv",
    )
    parser.add_argument(
        "--score-source", type=Path, default=authoritative / "src/hla2_positive_control_benchmark_v2.py"
    )
    parser.add_argument(
        "--predictors", type=Path, default=expanded / "raw_responses/predictor_records.csv"
    )
    parser.add_argument("--manifest", type=Path, default=expanded / "frozen_49_pair_manifest.csv")
    parser.add_argument("--protocol", type=Path, default=expanded / "protocol_lock.json")
    parser.add_argument(
        "--tiers",
        type=Path,
        default=expanded / "secondary_provisional_tiers_2026-09-07/all_49_provisional_tiers.csv",
    )
    parser.add_argument(
        "--threshold-summary",
        type=Path,
        default=workspace / "output/revision_round_2/threshold_sensitivity_2026-09-14/threshold_grid_summary.csv",
    )
    parser.add_argument(
        "--metric-leads",
        type=Path,
        default=workspace / "output/revision_round_2/metric_rank_robustness_2026-09-14/lead_metric_rank_robustness.csv",
    )
    parser.add_argument(
        "--null-results",
        type=Path,
        default=workspace / "output/revision_round_2/null_simulation_2026-09-14/null_pair_results.csv",
    )
    parser.add_argument(
        "--parent-nested", type=Path, default=focused / "af3_analysis/parent_nested_ensemble_comparisons.csv"
    )
    parser.add_argument(
        "--job-summary", type=Path, default=focused / "af3_analysis/job_summary_13.csv"
    )
    parser.add_argument(
        "--sample-metrics", type=Path, default=focused / "af3_analysis/sample_metrics_65.csv"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_3/computational_tightening_2026-09-15",
    )
    parser.add_argument("--matched-decoys", type=int, default=100)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_sha256(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def load_score_module(path: Path):
    spec = importlib.util.spec_from_file_location("frozen_score_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load score module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def competition_ranks(values: dict[str, float], higher_is_better: bool = True) -> dict[str, int]:
    ordered = sorted(values, key=lambda key: ((-values[key]) if higher_is_better else values[key], key))
    output: dict[str, int] = {}
    position = 0
    while position < len(ordered):
        end = position + 1
        value = values[ordered[position]]
        while end < len(ordered) and math.isclose(values[ordered[end]], value, abs_tol=1e-12):
            end += 1
        for key in ordered[position:end]:
            output[key] = position + 1
        position = end
    return output


def rank_against(values: list[float], target: float) -> int:
    return 1 + sum(value > target and not math.isclose(value, target, abs_tol=1e-12) for value in values)


def sequence_counts(sequence: str) -> np.ndarray:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    counts = Counter(sequence)
    return np.asarray([counts[aa] / max(1, len(sequence)) for aa in alphabet], dtype=float)


def empirical_quantile(values: list[float], target: float) -> float:
    return (1 + sum(value >= target for value in values)) / (len(values) + 1)


def svg_escape(value: str) -> str:
    return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def build_threshold_svg(rows: list[dict[str, str]], path: Path) -> None:
    cutoffs = ["1.0", "2.0", "5.0", "10.0", "20.0"]
    order = [
        ("all_3_both_arms", "exact_declared_core"),
        ("all_3_both_arms", "within_1_residue_of_declared_start"),
        ("at_least_2_of_3_both_arms", "exact_declared_core"),
        ("at_least_2_of_3_both_arms", "within_1_residue_of_declared_start"),
    ]
    lookup = {(r["binding_rule"], r["register_rule"], r["cutoff_percent"]): int(r["passing_pair_count"]) for r in rows}
    width, height = 980, 430
    x0, y0, cw, ch = 420, 105, 94, 58
    max_value = max(lookup.values()) if lookup else 1
    labels = {
        order[0]: "3/3 predictors, exact register",
        order[1]: "3/3 predictors, register within ±1",
        order[2]: "≥2/3 predictors, exact register",
        order[3]: "≥2/3 predictors, register within ±1",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="40" y="42" font-family="Arial" font-size="21" font-weight="bold">Binding/register threshold sensitivity</text>',
        '<text x="40" y="68" font-family="Arial" font-size="13" fill="#444">Cells show passing pairs in the frozen 49-pair panel</text>',
    ]
    for col, cutoff in enumerate(cutoffs):
        parts.append(f'<text x="{x0 + col*cw + cw/2}" y="94" text-anchor="middle" font-family="Arial" font-size="13">{cutoff.replace(".0", "")}%</text>')
    for row_index, key in enumerate(order):
        y = y0 + row_index * ch
        parts.append(f'<text x="400" y="{y + 34}" text-anchor="end" font-family="Arial" font-size="13">{svg_escape(labels[key])}</text>')
        for col, cutoff in enumerate(cutoffs):
            value = lookup.get((key[0], key[1], cutoff), 0)
            intensity = value / max_value if max_value else 0
            blue = int(245 - 145 * intensity)
            fill = f'rgb({blue},{blue+5},{255})'
            x = x0 + col*cw
            parts.append(f'<rect x="{x}" y="{y}" width="{cw-4}" height="{ch-4}" fill="{fill}" stroke="#777"/>')
            parts.append(f'<text x="{x+(cw-4)/2}" y="{y+34}" text-anchor="middle" font-family="Arial" font-size="18" font-weight="bold">{value}</text>')
    parts.extend([
        '<text x="40" y="372" font-family="Arial" font-size="12" fill="#444">Result: both leading BALF5–TALDO1 pairs pass the strict rule only at the 20% cutoff.</text>',
        '<text x="40" y="397" font-family="Arial" font-size="11" fill="#666">Computational gate sensitivity only; not measured binding or biological validation.</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def build_dashboard_svg(rows: list[dict], path: Path) -> None:
    width, height = 1050, 480
    columns = [
        ("declared_score_rank", "Declared\nrank"),
        ("worst_shift_rank", "Worst ±1\nrank"),
        ("matched_decoy_upper_tail", "Matched-decoy\nupper tail"),
        ("leave_one_out_worst_rank", "Worst leave-one-\nout rank"),
    ]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="40" y="42" font-family="Arial" font-size="21" font-weight="bold">Lead-candidate computational robustness dashboard</text>',
        '<text x="40" y="68" font-family="Arial" font-size="13" fill="#444">Lower ranks are stronger; matched-decoy upper tail is descriptive and post hoc</text>',
    ]
    x_positions = [330, 510, 710, 910]
    for (field, label), x in zip(columns, x_positions):
        line1, line2 = label.split("\n")
        parts.append(f'<text x="{x}" y="105" text-anchor="middle" font-family="Arial" font-size="13">{line1}</text>')
        parts.append(f'<text x="{x}" y="122" text-anchor="middle" font-family="Arial" font-size="13">{line2}</text>')
    for idx, row in enumerate(rows):
        y = 175 + idx * 125
        label = f'{row["target_id"]} · {row["allele"].replace("HLA-", "")}'
        parts.append(f'<text x="40" y="{y+23}" font-family="Arial" font-size="15" font-weight="bold">{svg_escape(label)}</text>')
        for (field, _), x in zip(columns, x_positions):
            raw = row[field]
            if "upper_tail" in field:
                display = f'{100*float(raw):.1f}%'
            else:
                display = str(raw)
            parts.append(f'<circle cx="{x}" cy="{y+16}" r="29" fill="#dbeafe" stroke="#3568a8"/>')
            parts.append(f'<text x="{x}" y="{y+22}" text-anchor="middle" font-family="Arial" font-size="15" font-weight="bold">{display}</text>')
        parts.append(f'<line x1="40" y1="{y+72}" x2="1010" y2="{y+72}" stroke="#ddd"/>')
    parts.extend([
        '<text x="40" y="445" font-family="Arial" font-size="11" fill="#666">Robustness analyses reuse related sequence information and do not constitute independent validation.</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    score_module = load_score_module(args.score_source)
    pairs = read_csv(args.pairs)
    audit = read_csv(args.audit)
    if len(pairs) != 6400:
        raise ValueError(f"expected 6400 source pairs, observed {len(pairs)}")
    if len({row["pair_id"] for row in pairs}) != 6400:
        raise ValueError("pair IDs are not unique")
    target_to_pair = {row["target_id"]: row["pair_id"] for row in audit}
    missing_targets = set(LEADS) - set(target_to_pair)
    if missing_targets:
        raise ValueError(f"missing lead target IDs: {sorted(missing_targets)}")

    # Exact score specification and worked lead examples.
    formula = {
        "name": "normalized_BLOSUM62_similarity",
        "equation": "sum_i BLOSUM62(left_i,right_i) / sum_i max(BLOSUM62(left_i,left_i),BLOSUM62(right_i,right_i))",
        "tcr_facing_positions_1_based": [2, 3, 5, 7, 8],
        "anchor_positions_1_based": [1, 4, 6, 9],
        "negative_substitution_handling": "retained without clipping",
        "tie_rule": "competition ranks; exact numeric equality within 1e-12 in this audit",
        "claim_boundary": CLAIM_BOUNDARY,
    }
    worked = []
    pair_by_id = {row["pair_id"]: row for row in pairs}
    for target in LEADS:
        row = pair_by_id[target_to_pair[target]]
        left, right = row["ebv_predicted_core"], row["self_predicted_core"]
        substitutions = [score_module._BLOSUM62[(left[i], right[i])] for i in TCR_POSITIONS]
        denominators = [
            max(score_module._BLOSUM62[(left[i], left[i])], score_module._BLOSUM62[(right[i], right[i])])
            for i in TCR_POSITIONS
        ]
        worked.append({
            "target_id": target,
            "pair_id": row["pair_id"],
            "left_core": left,
            "right_core": right,
            "position_scores": ";".join(map(str, substitutions)),
            "position_self_maxima": ";".join(map(str, denominators)),
            "numerator": sum(substitutions),
            "denominator": sum(denominators),
            "normalized_score": score_module.blosum62_similarity(left, right, positions=TCR_POSITIONS),
        })
    formula["worked_examples"] = worked
    (args.output_dir / "scoring_formula.json").write_text(json.dumps(formula, indent=2) + "\n", encoding="utf-8")

    # Confirm the frozen source columns were generated by the documented formula.
    score_differences = []
    for row in pairs:
        recalculated_tcr = score_module.blosum62_similarity(
            row["ebv_predicted_core"], row["self_predicted_core"], positions=TCR_POSITIONS
        )
        recalculated_full = score_module.blosum62_similarity(
            row["ebv_predicted_core"], row["self_predicted_core"]
        )
        score_differences.append((
            abs(recalculated_tcr - float(row["tcr_facing_blosum62_similarity"])),
            abs(recalculated_full - float(row["full_core_blosum62_similarity"])),
        ))
    max_tcr_score_difference = max(value[0] for value in score_differences)
    max_full_score_difference = max(value[1] for value in score_differences)

    # Declared-register provenance audit.
    provenance_counts = Counter((row["allele"], row["declared_register_status"]) for row in pairs)
    provenance_rows = [
        {"allele": allele, "declared_register_status": status, "pair_rows": count}
        for (allele, status), count in sorted(provenance_counts.items())
    ]
    write_csv(args.output_dir / "declared_register_provenance_counts.csv", provenance_rows)
    if any(row["declared_register_registry_unique"].lower() != "true" for row in pairs):
        raise ValueError("at least one pair lacks a unique declared register")

    # Frozen 49-pair selection provenance, including source membership and replacements.
    tier_by_pair = {row["pair_id"]: row for row in read_csv(args.tiers)}
    selection_rows = []
    for row in read_csv(args.manifest):
        tier = tier_by_pair.get(row["pair_id"], {})
        selection_rows.append({
            "candidate_id": tier.get("candidate_id", ""),
            "pair_id": row["pair_id"],
            "allele": row["allele"],
            "universe_source": row.get("universe_source", ""),
            "retention_rank": row.get("retention_rank", ""),
            "original_hla_rank": row.get("hla_rank", ""),
            "declared_ebv_core": row.get("ebv_core", row.get("ebv_predicted_core", "")),
            "declared_self_core": row.get("self_core", row.get("self_predicted_core", "")),
            "provisional_tier": tier.get("tier", ""),
            "selection_status": "frozen_before_current_tightening_analysis",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    if len(selection_rows) != 49:
        raise ValueError("frozen 49-pair manifest did not contain 49 rows")
    write_csv(args.output_dir / "candidate_selection_provenance_49.csv", selection_rows)
    (args.output_dir / "frozen_49_protocol_copy.json").write_text(
        args.protocol.read_text(encoding="utf-8"), encoding="utf-8"
    )

    # Score ablations across the full frozen universe.
    ablation_definitions = {
        "tcr_face": TCR_POSITIONS,
        "full_core": tuple(range(9)),
        "anchors_only": ANCHOR_POSITIONS,
    }
    for omitted in TCR_POSITIONS:
        ablation_definitions[f"tcr_face_without_P{omitted+1}"] = tuple(i for i in TCR_POSITIONS if i != omitted)
    ablation_rows = []
    lead_ablation_rows = []
    for allele in sorted({row["allele"] for row in pairs}):
        subset = [row for row in pairs if row["allele"] == allele]
        metric_values: dict[str, dict[str, float]] = {}
        metric_ranks: dict[str, dict[str, int]] = {}
        for name, positions in ablation_definitions.items():
            values = {
                row["pair_id"]: score_module.blosum62_similarity(
                    row["ebv_predicted_core"], row["self_predicted_core"], positions=positions
                )
                for row in subset
            }
            metric_values[name] = values
            metric_ranks[name] = competition_ranks(values)
        for row in subset:
            output = {
                "allele": allele,
                "pair_id": row["pair_id"],
                "ebv_protein": row["ebv_protein"],
                "self_protein": row["self_protein"],
                "ebv_core": row["ebv_predicted_core"],
                "self_core": row["self_predicted_core"],
                "claim_boundary": CLAIM_BOUNDARY,
            }
            for name in ablation_definitions:
                output[f"{name}_score"] = f'{metric_values[name][row["pair_id"]]:.12g}'
                output[f"{name}_rank"] = metric_ranks[name][row["pair_id"]]
            ablation_rows.append(output)
            targets = [target for target, pair_id in target_to_pair.items() if pair_id == row["pair_id"]]
            for target in targets:
                if target in LEADS:
                    lead = dict(output)
                    lead["target_id"] = target
                    leave_out = [metric_ranks[f"tcr_face_without_P{p+1}"][row["pair_id"]] for p in TCR_POSITIONS]
                    lead["leave_one_out_best_rank"] = min(leave_out)
                    lead["leave_one_out_worst_rank"] = max(leave_out)
                    lead_ablation_rows.append(lead)
    write_csv(args.output_dir / "score_ablation_all_6400.csv", ablation_rows)
    write_csv(args.output_dir / "score_ablation_leads.csv", lead_ablation_rows)

    # +/-1 register envelope for all pairs, ranked against each allele's frozen declared-score distribution.
    register_detail = []
    register_summary = []
    for allele in sorted({row["allele"] for row in pairs}):
        subset = [row for row in pairs if row["allele"] == allele]
        reference_scores = [float(row["tcr_facing_blosum62_similarity"]) for row in subset]
        for row in subset:
            left_start = int(row["ebv_declared_core_start_1_based"]) - 1
            right_start = int(row["self_declared_core_start_1_based"]) - 1
            left_sequence, right_sequence = row["ebv_sequence"], row["self_sequence"]
            variations = []
            for left_delta in (-1, 0, 1):
                for right_delta in (-1, 0, 1):
                    ls, rs = left_start + left_delta, right_start + right_delta
                    if ls < 0 or rs < 0 or ls + 9 > len(left_sequence) or rs + 9 > len(right_sequence):
                        continue
                    left_core = left_sequence[ls:ls+9]
                    right_core = right_sequence[rs:rs+9]
                    score = score_module.blosum62_similarity(left_core, right_core, positions=TCR_POSITIONS)
                    rank = rank_against(reference_scores, score)
                    detail = {
                        "allele": allele,
                        "pair_id": row["pair_id"],
                        "left_register_delta": left_delta,
                        "right_register_delta": right_delta,
                        "left_core": left_core,
                        "right_core": right_core,
                        "score": f"{score:.12g}",
                        "rank_against_frozen_declared_universe": rank,
                        "is_declared_pair": left_delta == 0 and right_delta == 0,
                        "claim_boundary": CLAIM_BOUNDARY,
                    }
                    register_detail.append(detail)
                    variations.append(detail)
            declared = next(item for item in variations if item["is_declared_pair"])
            summary = {
                "allele": allele,
                "pair_id": row["pair_id"],
                "declared_score": declared["score"],
                "declared_score_rank": declared["rank_against_frozen_declared_universe"],
                "best_shift_score": max(float(item["score"]) for item in variations),
                "worst_shift_score": min(float(item["score"]) for item in variations),
                "best_shift_rank": min(int(item["rank_against_frozen_declared_universe"]) for item in variations),
                "worst_shift_rank": max(int(item["rank_against_frozen_declared_universe"]) for item in variations),
                "evaluated_register_combinations": len(variations),
                "declared_is_best_score": float(declared["score"]) == max(float(item["score"]) for item in variations),
                "claim_boundary": CLAIM_BOUNDARY,
            }
            targets = [target for target, pair_id in target_to_pair.items() if pair_id == row["pair_id"]]
            if targets:
                summary["target_ids"] = ";".join(sorted(targets))
            register_summary.append(summary)
    write_csv(args.output_dir / "register_shift_detail_all_6400.csv", register_detail)
    write_csv(args.output_dir / "register_shift_summary_all_6400.csv", register_summary)
    lead_register = [row for row in register_summary if any(target in row.get("target_ids", "").split(";") for target in LEADS)]
    for row in lead_register:
        row["target_id"] = next(target for target in LEADS if target in row["target_ids"].split(";"))
    write_csv(args.output_dir / "register_shift_summary_leads.csv", lead_register)

    # Deterministic post-hoc matched empirical decoys.
    matched_detail = []
    matched_summary = []
    for target in LEADS:
        observed = pair_by_id[target_to_pair[target]]
        allele_rows = [row for row in pairs if row["allele"] == observed["allele"]]
        observed_left_comp = sequence_counts("".join(observed["ebv_predicted_core"][i] for i in TCR_POSITIONS))
        observed_right_comp = sequence_counts("".join(observed["self_predicted_core"][i] for i in TCR_POSITIONS))
        candidates = []
        for row in allele_rows:
            if row["pair_id"] == observed["pair_id"]:
                continue
            if row["ebv_candidate_id"] == observed["ebv_candidate_id"] or row["self_candidate_id"] == observed["self_candidate_id"]:
                continue
            left_binding = abs(math.log1p(float(row["ebv_binding_percentile_rank"])) - math.log1p(float(observed["ebv_binding_percentile_rank"])))
            right_binding = abs(math.log1p(float(row["self_binding_percentile_rank"])) - math.log1p(float(observed["self_binding_percentile_rank"])))
            length_penalty = 0.5 * abs(len(row["ebv_sequence"]) - len(observed["ebv_sequence"])) + 0.5 * abs(len(row["self_sequence"]) - len(observed["self_sequence"]))
            left_comp = sequence_counts("".join(row["ebv_predicted_core"][i] for i in TCR_POSITIONS))
            right_comp = sequence_counts("".join(row["self_predicted_core"][i] for i in TCR_POSITIONS))
            composition = float(np.abs(left_comp - observed_left_comp).sum() + np.abs(right_comp - observed_right_comp).sum())
            distance = left_binding + right_binding + length_penalty + composition
            candidates.append((distance, row["pair_id"], row))
        candidates.sort(key=lambda item: (item[0], item[1]))
        selected = candidates[: args.matched_decoys]
        scores = [float(row["tcr_facing_blosum62_similarity"]) for _, _, row in selected]
        observed_score = float(observed["tcr_facing_blosum62_similarity"])
        for rank, (distance, _, row) in enumerate(selected, 1):
            matched_detail.append({
                "target_id": target,
                "match_rank": rank,
                "match_distance": f"{distance:.12g}",
                "decoy_pair_id": row["pair_id"],
                "decoy_ebv_protein": row["ebv_protein"],
                "decoy_self_protein": row["self_protein"],
                "decoy_score": row["tcr_facing_blosum62_similarity"],
                "decoy_ebv_binding_percentile": row["ebv_binding_percentile_rank"],
                "decoy_self_binding_percentile": row["self_binding_percentile_rank"],
                "shares_neither_peptide_arm_with_target": True,
                "claim_boundary": CLAIM_BOUNDARY,
            })
        matched_summary.append({
            "target_id": target,
            "allele": observed["allele"],
            "observed_pair_id": observed["pair_id"],
            "observed_score": observed_score,
            "matched_decoy_count": len(selected),
            "matched_decoy_score_median": statistics.median(scores),
            "matched_decoy_score_max": max(scores),
            "matched_decoy_upper_tail": empirical_quantile(scores, observed_score),
            "analysis_status": "post_hoc_descriptive_not_confirmatory",
            "matching_variables": "allele;arm-specific binding percentiles;peptide lengths;arm-specific TCR-face composition",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(args.output_dir / "matched_empirical_decoys.csv", matched_detail)
    write_csv(args.output_dir / "matched_empirical_decoy_summary.csv", matched_summary)

    # Correct interpretation of structural ensemble summaries: 25 cross-comparisons are dependent.
    job_by_name = {row["job"]: row for row in read_csv(args.job_summary)}
    structural_rows = []
    for row in read_csv(args.parent_nested):
        parent = job_by_name[row["parent_job"]]
        nested = job_by_name[row["nested_job"]]
        structural_rows.append({
            "candidate": row["candidate"],
            "arm": row["arm"],
            "parent_job": row["parent_job"],
            "nested_job": row["nested_job"],
            "parent_models": parent["models"],
            "nested_models": nested["models"],
            "independent_model_count": int(parent["models"]) + int(nested["models"]),
            "dependent_cross_comparisons": row["model_pair_comparisons"],
            "within_parent_pairwise_median_A": parent["within_job_pairwise_peptide_rmsd_median_A"],
            "within_parent_pairwise_max_A": parent["within_job_pairwise_peptide_rmsd_max_A"],
            "within_nested_pairwise_median_A": nested["within_job_pairwise_peptide_rmsd_median_A"],
            "within_nested_pairwise_max_A": nested["within_job_pairwise_peptide_rmsd_max_A"],
            "between_input_cross_median_A": row["median_A"],
            "between_input_cross_min_A": row["min_A"],
            "between_input_cross_max_A": row["max_A"],
            "statistical_status": "descriptive_dependent_cross_comparisons_not_n_equals_25",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(args.output_dir / "structural_ensemble_corrected_summary.csv", structural_rows)

    # Bring model-level peptide confidence into one directly usable supplement.
    model_rows = read_csv(args.sample_metrics)
    for row in model_rows:
        row["model_is_biological_replicate"] = "False"
        row["claim_boundary"] = CLAIM_BOUNDARY
    write_csv(args.output_dir / "af3_model_level_confidence_65.csv", model_rows)

    # Consolidated lead dashboard.
    decoy_by_target = {row["target_id"]: row for row in matched_summary}
    ablation_by_target = {row["target_id"]: row for row in lead_ablation_rows}
    register_by_target = {row["target_id"]: row for row in lead_register}
    dashboard = []
    for target in LEADS:
        dashboard.append({
            "target_id": target,
            "allele": register_by_target[target]["allele"],
            "declared_score_rank": register_by_target[target]["declared_score_rank"],
            "best_shift_rank": register_by_target[target]["best_shift_rank"],
            "worst_shift_rank": register_by_target[target]["worst_shift_rank"],
            "declared_is_best_score": register_by_target[target]["declared_is_best_score"],
            "leave_one_out_best_rank": ablation_by_target[target]["leave_one_out_best_rank"],
            "leave_one_out_worst_rank": ablation_by_target[target]["leave_one_out_worst_rank"],
            "anchor_only_rank": ablation_by_target[target]["anchors_only_rank"],
            "matched_decoy_upper_tail": decoy_by_target[target]["matched_decoy_upper_tail"],
            "matched_decoy_count": decoy_by_target[target]["matched_decoy_count"],
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(args.output_dir / "lead_robustness_dashboard.csv", dashboard)

    threshold_rows = read_csv(args.threshold_summary)
    build_threshold_svg(threshold_rows, args.output_dir / "figure_threshold_sensitivity.svg")
    build_dashboard_svg(dashboard, args.output_dir / "figure_lead_robustness_dashboard.svg")

    # Freeze computational hypotheses before experimental outcomes are available.
    predictor_rows = read_csv(args.predictors)
    tier_rows = read_csv(args.tiers)
    frozen_tiers = [row for row in tier_rows if row.get("candidate_id") in {"HP36", "HP47"}]
    frozen_predictors = [row for row in predictor_rows if row.get("target_id") in {"HP36", "HP47"}]
    prospective = {
        "status": "prepared_not_submitted",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "freeze pre-experiment computational hypotheses and thresholds",
        "experimental_results_seen_or_used": False,
        "lead_candidate_ids": ["HP36", "HP47"],
        "lead_tier_rows": frozen_tiers,
        "lead_predictor_rows": frozen_predictors,
        "binding_threshold_grid": [1, 2, 5, 10, 20],
        "primary_locked_gate": "all three predictors <=20 percent for both arms; at least two predictors recover declared core for both arms",
        "structural_interpretation": "register unresolved; structural outputs are sensitivity evidence, not confirmation",
        "score_formula_sha256": sha256(args.output_dir / "scoring_formula.json"),
        "claim_boundary": CLAIM_BOUNDARY,
    }
    prospective["protocol_sha256"] = stable_sha256(prospective)
    (args.output_dir / "prospective_preexperiment_lock.json").write_text(
        json.dumps(prospective, indent=2) + "\n", encoding="utf-8"
    )

    input_paths = [
        args.pairs, args.audit, args.score_source, args.predictors, args.manifest,
        args.protocol, args.tiers,
        args.threshold_summary, args.metric_leads, args.null_results,
        args.parent_nested, args.job_summary, args.sample_metrics,
    ]
    lock = {
        "analysis": "manuscript_v8_computational_tightening",
        "status": "complete_for_local_frozen_analyses",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "source_pair_count": len(pairs),
        "matched_decoys_per_lead": args.matched_decoys,
        "inputs": {str(path): sha256(path) for path in input_paths},
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (args.output_dir / "analysis_lock.json").write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    readme = f"""# Computational tightening package

Status: complete for the frozen local analyses described below. Manuscript prose was not edited.

## Added analyses

- exact normalized BLOSUM62 equation and reconstructable worked examples;
- declared-register provenance counts for the 6,400-pair universe;
- frozen 49-pair source-membership and deterministic-selection provenance;
- full-universe anchor/full-core/TCR-face and leave-one-position-out ablations;
- all-pair +/-1 register-shift sensitivity envelopes;
- deterministic post-hoc empirical decoys matched by allele, binding ranks, length, and exposed composition;
- corrected structural-ensemble reporting that distinguishes ten underlying models from 25 dependent cross-comparisons;
- consolidated 65-model peptide-local confidence table;
- threshold and lead-robustness figures;
- prepared-not-submitted prospective pre-experiment protocol lock.

## Deliberately not claimed

- The matched-decoy analysis is post hoc and descriptive, not confirmatory inference.
- A register shift is a sequence-level sensitivity test, not a structural or experimental register assignment.
- The saved AF3 parent/nested analysis can quantify coordinate sensitivity, but a universally calibrated residue-to-pocket structural-register caller for every allele is not established here.
- No experimental binding, presentation, TCR recognition, cross-reactivity, molecular mimicry, or MS mechanism is established.

## Source counts

- Frozen pair universe: {len(pairs)} rows.
- Score-ablation rows: {len(ablation_rows)}.
- Register-shift detail rows: {len(register_detail)}.
- Matched empirical decoys: {len(matched_detail)}.
- AF3 model-level confidence rows: {len(model_rows)}.
"""
    (args.output_dir / "README.md").write_text(readme, encoding="utf-8")

    qa_lines = [
        "# QA",
        "",
        f"- 6,400 unique pair IDs: {'PASS' if len(pairs) == len(set(row['pair_id'] for row in pairs)) == 6400 else 'FAIL'}",
        f"- Four 1,600-pair allele cohorts: {'PASS' if sorted(Counter(row['allele'] for row in pairs).values()) == [1600]*4 else 'FAIL'}",
        f"- Unique declared registers: {'PASS' if all(row['declared_register_registry_unique'].lower() == 'true' for row in pairs) else 'FAIL'}",
        f"- TCR-face score formula reproduces all frozen values (max difference {max_tcr_score_difference:.3g}): {'PASS' if max_tcr_score_difference < 1e-10 else 'FAIL'}",
        f"- Full-core score formula reproduces all frozen values (max difference {max_full_score_difference:.3g}): {'PASS' if max_full_score_difference < 1e-10 else 'FAIL'}",
        f"- Frozen selection provenance rows: {'PASS' if len(selection_rows) == 49 else 'FAIL'}",
        f"- Two lead ablation rows: {'PASS' if len(lead_ablation_rows) == 2 else 'FAIL'}",
        f"- Two lead register envelopes: {'PASS' if len(lead_register) == 2 else 'FAIL'}",
        f"- {args.matched_decoys} decoys per lead: {'PASS' if all(int(row['matched_decoy_count']) == args.matched_decoys for row in matched_summary) else 'FAIL'}",
        f"- 65 AF3 model rows: {'PASS' if len(model_rows) == 65 else 'FAIL'}",
        f"- Prospective lock marked prepared, not submitted: PASS",
    ]
    (args.output_dir / "QA.md").write_text("\n".join(qa_lines) + "\n", encoding="utf-8")

    supplement_sources = [
        ("frozen_6400_pair_universe", args.pairs),
        ("frozen_49_pair_manifest", args.manifest),
        ("frozen_49_predictor_records", args.predictors),
        ("frozen_49_provisional_tiers", args.tiers),
        ("frozen_49_protocol", args.protocol),
    ]
    supplement_rows = [
        {
            "artifact": label,
            "path": str(path.resolve()),
            "sha256": sha256(path),
            "status": "authoritative_source_not_modified",
        }
        for label, path in supplement_sources
    ]
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file() and path.name not in {"SHA256SUMS.csv", "SUPPLEMENT_INDEX.csv"}:
            supplement_rows.append({
                "artifact": path.stem,
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "status": "derived_current_package",
            })
    write_csv(args.output_dir / "SUPPLEMENT_INDEX.csv", supplement_rows)

    output_files = sorted(path for path in args.output_dir.iterdir() if path.is_file() and path.name != "SHA256SUMS.csv")
    checksum_rows = [{"path": path.name, "sha256": sha256(path)} for path in output_files]
    write_csv(args.output_dir / "SHA256SUMS.csv", checksum_rows)
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "files": len(output_files) + 1,
        "dashboard": dashboard,
    }, indent=2))


if __name__ == "__main__":
    main()
