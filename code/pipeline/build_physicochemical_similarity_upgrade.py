"""Run the frozen Grantham/Atchley physicochemical sensitivity analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from amino_acid_physicochemical_distance import (
    TCR_FACING_INDICES,
    tcr_face_atchley_distance,
    tcr_face_grantham_mismatch,
)
from hla2_positive_control_benchmark import physicochemical_mismatch
from hla2_positive_control_benchmark_v2 import blosum62_similarity


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DISCOVERY = (
    ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25/discovery/pair_summary_6400.csv"
)
DEFAULT_CONTROL_ROOT = ROOT / "processed/hla2_positive_control_benchmark_v2_pilot_2026-08-26"
DEFAULT_FROZEN_RANKS = (
    ROOT
    / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26/benchmark/method_rank_long.csv"
)
DEFAULT_OUT = ROOT / "processed/physicochemical_similarity_upgrade_2026-09-14"
ALLELES = (
    "HLA-DRB1*15:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
)
AUDITED_BALF5_TALDO1_PAIR_IDS = (
    "HLA-DRB1*13:03|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0108_0122",
    "HLA-DRB1*15:01|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0216_0230",
)
METHOD_FIELDS = {
    "current_five_property": "tcr_face_physicochemical_mismatch",
    "grantham": "tcr_face_grantham_mismatch",
    "atchley": "tcr_face_atchley_distance",
}
CLAIM_BOUNDARY = (
    "Sequence-property sensitivity analysis only; not evidence of presentation, "
    "TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism."
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] | None = None,
) -> None:
    fieldnames = list(fields or (list(rows[0]) if rows else []))
    if not fieldnames:
        raise ValueError(f"cannot write empty CSV without fields: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def competition_ranks(
    values: Sequence[float], *, lower_is_better: bool = True
) -> list[int]:
    numeric = [float(value) for value in values]
    return [
        1 + sum(
            candidate < value if lower_is_better else candidate > value
            for candidate in numeric
        )
        for value in numeric
    ]


def evaluate_grantham_promotion(
    current: Mapping[str, Any], grantham: Mapping[str, Any]
) -> dict[str, Any]:
    def objective(row: Mapping[str, Any]) -> tuple[int, int, float]:
        return (
            int(row["system_capture_at_3_count"]),
            -int(row["worst_system_rank"]),
            float(row["system_weighted_mrr"]),
        )

    current_objective = objective(current)
    grantham_objective = objective(grantham)
    promoted = grantham_objective >= current_objective
    return {
        "promote_grantham": promoted,
        "decision": (
            "promote_grantham_as_primary_physicochemical_metric"
            if promoted
            else "do_not_promote_no_primary_robustness_label"
        ),
        "predeclared_comparison_order": [
            "maximize_system_capture_at_3_count",
            "minimize_worst_system_rank",
            "maximize_system_weighted_mrr",
        ],
        "current_objective": list(current_objective),
        "grantham_objective": list(grantham_objective),
        "candidate_results_used_in_decision": False,
    }


def _control_core_lookup(control_root: Path) -> dict[str, str]:
    lookup: dict[str, str] = {}
    sources = (
        (control_root / "registry/control_ligand_registry.csv", "ligand_id", "core"),
        (control_root / "controls/control_decoy_registry.csv", "candidate_id", "predicted_core"),
    )
    for path, id_field, core_field in sources:
        for row in read_csv(path):
            identifier = str(row[id_field])
            core = str(row[core_field])
            if identifier in lookup and lookup[identifier] != core:
                raise ValueError(f"conflicting cores for control ligand {identifier}")
            if len(core) != 9:
                raise ValueError(f"invalid control core for {identifier}: {core}")
            lookup[identifier] = core
    return lookup


def build_control_metric_rows(control_root: Path) -> list[dict[str, Any]]:
    lookup = _control_core_lookup(control_root)
    universe = read_csv(control_root / "controls/comparison_universe.csv")
    if len(universe) != 208:
        raise ValueError(f"expected 208 frozen control comparisons, observed {len(universe)}")
    output = []
    for source in universe:
        left_id = str(source["left_id"])
        right_id = str(source["right_id"])
        if left_id not in lookup or right_id not in lookup:
            raise ValueError(f"missing frozen core for control comparison {source['pair_id']}")
        left_core = lookup[left_id]
        right_core = lookup[right_id]
        output.append({
            **dict(source),
            "left_core": left_core,
            "right_core": right_core,
            "tcr_face_physicochemical_mismatch": round(
                physicochemical_mismatch(left_core, right_core), 12
            ),
            "tcr_face_grantham_mismatch": round(
                tcr_face_grantham_mismatch(left_core, right_core), 12
            ),
            "tcr_face_atchley_distance": round(
                tcr_face_atchley_distance(left_core, right_core), 12
            ),
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return output


def rank_control_panels(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    panels: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for source in rows:
        panels[(str(source["positive_pair_id"]), int(source["panel_seed"]))].append(dict(source))
    if len(panels) != 8:
        raise ValueError(f"expected eight frozen control panels, observed {len(panels)}")

    annotated = []
    method_ranks = []
    for (positive_pair_id, panel_seed), group in sorted(panels.items()):
        if len(group) != 26:
            raise ValueError(f"control panel {positive_pair_id}/{panel_seed} has {len(group)} rows")
        positives = [row for row in group if str(row["comparison_role"]) == "positive"]
        if len(positives) != 1:
            raise ValueError(f"control panel {positive_pair_id}/{panel_seed} requires one positive")
        positive_pair_row_id = str(positives[0]["pair_id"])
        rank_maps: dict[str, dict[str, int]] = {}
        for method, field in METHOD_FIELDS.items():
            ranks = competition_ranks([float(row[field]) for row in group])
            rank_maps[method] = {
                str(row["pair_id"]): rank for row, rank in zip(group, ranks)
            }
            positive_rank = rank_maps[method][positive_pair_row_id]
            method_ranks.append({
                "system_id": positives[0]["system_id"],
                "positive_pair_id": positive_pair_id,
                "panel_seed": panel_seed,
                "method": method,
                "positive_rank": positive_rank,
                "capture_at_3": positive_rank <= 3,
                "comparison_count": len(group),
            })
        for row in group:
            annotated.append({
                **row,
                **{
                    f"{method}_rank": rank_maps[method][str(row["pair_id"])]
                    for method in METHOD_FIELDS
                },
            })
    return annotated, method_ranks


def summarize_control_methods(
    method_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    output = []
    for method in METHOD_FIELDS:
        rows = [row for row in method_rows if str(row["method"]) == method]
        by_system: dict[str, list[int]] = defaultdict(list)
        for row in rows:
            by_system[str(row["system_id"])].append(int(row["positive_rank"]))
        if len(rows) != 8 or len(by_system) != 3:
            raise ValueError(f"incomplete frozen controls for {method}")
        system_worst = {system: max(ranks) for system, ranks in sorted(by_system.items())}
        ranks = list(system_worst.values())
        output.append({
            "method": method,
            "control_panel_count": len(rows),
            "independent_system_count": len(ranks),
            "system_capture_at_3_count": sum(rank <= 3 for rank in ranks),
            "worst_system_rank": max(ranks),
            "system_weighted_mrr": round(sum(1.0 / rank for rank in ranks) / len(ranks), 8),
            "system_worst_ranks": ";".join(
                f"{system}:{rank}" for system, rank in system_worst.items()
            ),
            "all_panels_capture_at_3": all(int(row["positive_rank"]) <= 3 for row in rows),
        })
    return output


def verify_current_control_ranks(
    recalculated: Sequence[Mapping[str, Any]], frozen_rank_path: Path
) -> None:
    frozen = {
        (str(row["positive_pair_id"]), int(row["panel_seed"])): int(row["positive_rank"])
        for row in read_csv(frozen_rank_path)
        if str(row["method"]) == "physicochemical_only"
    }
    observed = {
        (str(row["positive_pair_id"]), int(row["panel_seed"])): int(row["positive_rank"])
        for row in recalculated
        if str(row["method"]) == "current_five_property"
    }
    if observed != frozen:
        raise ValueError("recalculated current metric ranks do not match the frozen benchmark")


def _identity(left: str, right: str) -> float:
    return sum(left[index] == right[index] for index in TCR_FACING_INDICES) / len(
        TCR_FACING_INDICES
    )


def build_discovery_rows(
    source_rows: Sequence[Mapping[str, Any]], *, promote_grantham: bool
) -> list[dict[str, Any]]:
    if len(source_rows) != 6400:
        raise ValueError(f"expected 6,400 discovery pairs, observed {len(source_rows)}")
    output = []
    for source in source_rows:
        left = str(source["ebv_predicted_core"])
        right = str(source["self_predicted_core"])
        output.append({
            **dict(source),
            "tcr_facing_blosum62_similarity": round(
                blosum62_similarity(left, right, positions=TCR_FACING_INDICES), 12
            ),
            "full_core_blosum62_similarity": round(blosum62_similarity(left, right), 12),
            "tcr_facing_sequence_identity": round(_identity(left, right), 12),
            "tcr_face_physicochemical_mismatch": round(
                physicochemical_mismatch(left, right), 12
            ),
            "tcr_face_grantham_mismatch": round(
                tcr_face_grantham_mismatch(left, right), 12
            ),
            "tcr_face_atchley_distance": round(tcr_face_atchley_distance(left, right), 12),
            "claim_boundary_physicochemical_upgrade": CLAIM_BOUNDARY,
        })

    ranking_specs = (
        ("tcr_facing_blosum62_similarity", "tcr_facing_blosum62_rank", False),
        ("full_core_blosum62_similarity", "full_core_blosum62_rank", False),
        ("tcr_facing_sequence_identity", "tcr_facing_sequence_identity_rank", False),
        ("tcr_face_physicochemical_mismatch", "tcr_face_current_mismatch_rank", True),
        ("tcr_face_grantham_mismatch", "tcr_face_grantham_rank", True),
        ("tcr_face_atchley_distance", "tcr_face_atchley_rank", True),
    )
    by_allele: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in output:
        by_allele[str(row["allele"])].append(row)
    if set(by_allele) != set(ALLELES) or any(len(by_allele[a]) != 1600 for a in ALLELES):
        raise ValueError("discovery universe must contain four allele-specific 1,600-pair cohorts")

    for allele in ALLELES:
        group = by_allele[allele]
        for field, rank_field, lower_is_better in ranking_specs:
            ranks = competition_ranks(
                [float(row[field]) for row in group], lower_is_better=lower_is_better
            )
            for row, rank in zip(group, ranks):
                row[rank_field] = rank
        top_one = math.ceil(len(group) * 0.01)
        top_five = math.ceil(len(group) * 0.05)
        for row in group:
            primary_ranks = (
                int(row["tcr_facing_blosum62_rank"]),
                int(row["full_core_blosum62_rank"]),
                int(row["tcr_facing_sequence_identity_rank"]),
                int(row["tcr_face_grantham_rank"]),
            )
            row["primary_physicochemical_metric"] = (
                "tcr_face_grantham_mismatch" if promote_grantham else ""
            )
            if not promote_grantham:
                row["robustness_metric_count_top1"] = ""
                row["robustness_metric_count_top5"] = ""
                row["robustness_label"] = "not_assigned_grantham_failed_control_gate"
                continue
            count_one = sum(rank <= top_one for rank in primary_ranks)
            count_five = sum(rank <= top_five for rank in primary_ranks)
            row["robustness_metric_count_top1"] = count_one
            row["robustness_metric_count_top5"] = count_five
            row["robustness_label"] = (
                "robust"
                if count_one >= 3
                else "partially_robust"
                if count_five >= 3
                else "not_robust"
            )
    return output


def run(
    *,
    discovery_path: Path = DEFAULT_DISCOVERY,
    control_root: Path = DEFAULT_CONTROL_ROOT,
    frozen_rank_path: Path = DEFAULT_FROZEN_RANKS,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    control_rows = build_control_metric_rows(control_root)
    ranked_controls, method_rows = rank_control_panels(control_rows)
    verify_current_control_ranks(method_rows, frozen_rank_path)
    summaries = summarize_control_methods(method_rows)
    summary_by_method = {str(row["method"]): row for row in summaries}
    gate = evaluate_grantham_promotion(
        summary_by_method["current_five_property"], summary_by_method["grantham"]
    )
    gate.update({
        "analysis_version": "EBV_MS_PHYSICOCHEMICAL_SIMILARITY_UPGRADE_2026-09-14",
        "control_universe_verified_against_frozen_current_metric": True,
        "independent_control_system_count": 3,
        "control_panel_count": 8,
        "control_comparison_count": 208,
        "grantham_source_doi": "10.1126/science.185.4154.862",
        "atchley_source_doi": "10.1073/pnas.0408677102",
        "claim_boundary": CLAIM_BOUNDARY,
    })

    discovery_rows = build_discovery_rows(
        read_csv(discovery_path), promote_grantham=bool(gate["promote_grantham"])
    )
    candidates = [
        row for row in discovery_rows
        if str(row.get("ebv_protein", "")).upper() == "BALF5"
        and str(row.get("self_protein", "")).upper() == "TALDO1"
    ]
    if not candidates:
        raise ValueError("no BALF5-TALDO1 discovery pairs were found")
    candidates_by_id = {str(row["pair_id"]): row for row in candidates}
    missing_audited = set(AUDITED_BALF5_TALDO1_PAIR_IDS) - set(candidates_by_id)
    if missing_audited:
        raise ValueError(f"missing audited BALF5-TALDO1 pairs: {sorted(missing_audited)}")
    audited_candidates = [candidates_by_id[pair_id] for pair_id in AUDITED_BALF5_TALDO1_PAIR_IDS]

    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "control_pair_metrics.csv", ranked_controls)
    write_csv(out / "control_method_positive_ranks.csv", method_rows)
    write_csv(out / "control_method_summary.csv", summaries)
    write_json(out / "grantham_promotion_gate.json", gate)
    write_csv(out / "all_6400_pair_physicochemical_sensitivity.csv", discovery_rows)
    write_csv(out / "balf5_taldo1_physicochemical_sensitivity.csv", candidates)
    write_csv(out / "audited_balf5_taldo1_pairs.csv", audited_candidates)

    manifest = {
        **gate,
        "discovery_pair_count": len(discovery_rows),
        "discovery_pairs_by_allele": {
            allele: sum(str(row["allele"]) == allele for row in discovery_rows)
            for allele in ALLELES
        },
        "balf5_taldo1_row_count": len(candidates),
        "audited_balf5_taldo1_pair_count": len(audited_candidates),
        "audited_balf5_taldo1_pair_ids": list(AUDITED_BALF5_TALDO1_PAIR_IDS),
        "ranking_scope": "within_allele_only",
        "tcr_facing_positions": ["P2", "P3", "P5", "P7", "P8"],
        "rank_tie_policy": "competition_rank",
        "top_1_percent_rank_cutoff": 16,
        "top_5_percent_rank_cutoff": 80,
        "atchley_role": "sensitivity_only",
        "inputs": {
            "discovery_pair_summary": str(discovery_path),
            "control_root": str(control_root),
            "frozen_control_method_ranks": str(frozen_rank_path),
        },
        "input_sha256": {
            "discovery_pair_summary": sha256_file(discovery_path),
            "control_comparison_universe": sha256_file(
                control_root / "controls/comparison_universe.csv"
            ),
            "control_ligand_registry": sha256_file(
                control_root / "registry/control_ligand_registry.csv"
            ),
            "control_decoy_registry": sha256_file(
                control_root / "controls/control_decoy_registry.csv"
            ),
            "frozen_control_method_ranks": sha256_file(frozen_rank_path),
        },
    }
    write_json(out / "analysis_manifest.json", manifest)
    checksums = [
        {
            "relative_path": str(path.relative_to(out)),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in sorted(out.rglob("*"), key=str)
        if path.is_file() and path.name != "SHA256SUMS.csv"
    ]
    write_csv(out / "SHA256SUMS.csv", checksums)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery", type=Path, default=DEFAULT_DISCOVERY)
    parser.add_argument("--control-root", type=Path, default=DEFAULT_CONTROL_ROOT)
    parser.add_argument("--frozen-ranks", type=Path, default=DEFAULT_FROZEN_RANKS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    manifest = run(
        discovery_path=args.discovery,
        control_root=args.control_root,
        frozen_rank_path=args.frozen_ranks,
        out=args.out,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
