#!/usr/bin/env python3
"""Reproduce frozen-data V18 rank and routing ablations.

The analyses are sensitivity checks on the archived computational workflow;
they neither estimate binding nor validate T-cell cross-reactivity.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WORKSPACE / "code" / "pipeline"))
from hla2_positive_control_benchmark_v2 import _BLOSUM62  # noqa: E402

SOURCE = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/processed 2"
)
V3 = SOURCE / "literature_grounded_hla2_rankings_v3_2026-08-27"
PANEL = SOURCE / "high_priority_handoff_hunt_2026-09-07"
OUTPUT = WORKSPACE / "output" / "revision_round_18" / "frozen_data_ablations"
POSITIONS = (1, 2, 4, 6, 7)  # zero-based P2/P3/P5/P7/P8
PREDICTORS = ("iedb_recommended_binding", "netmhciipan_4_3_el", "mixmhc2pred_2_1_context")
LEADS = {
    "HLA-DRB1*13:03|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0108_0122",
    "HLA-DRB1*15:01|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0216_0230",
}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)


def pearson(left: list[float], right: list[float]) -> float:
    mean_left, mean_right = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((x - mean_left) * (y - mean_right) for x, y in zip(left, right))
    denominator = math.sqrt(sum((x - mean_left) ** 2 for x in left) * sum((y - mean_right) ** 2 for y in right))
    return numerator / denominator


def ranked(rows: list[dict[str, str]], key) -> list[dict[str, str]]:
    return sorted(rows, key=key)


def main() -> None:
    all_rows = read(V3 / "v3_all_hla_ranked_pairs.csv")
    panel = read(PANEL / "frozen_49_pair_manifest.csv")
    tiers = read(PANEL / "secondary_provisional_tiers_2026-09-07/all_49_provisional_tiers.csv")
    records = read(PANEL / "raw_responses/predictor_records.csv")
    if len(all_rows) != 6400 or len(panel) != 49 or len(tiers) != 49 or len(records) != 490:
        raise ValueError("unexpected frozen input size")

    by_allele: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        core_a, core_b = row["ebv_core_p1_p9"], row["self_core_p1_p9"]
        row["_normalized"] = float(row["tcr_facing_blosum62_similarity"])
        row["_raw_sum"] = sum(_BLOSUM62[(core_a[i], core_b[i])] for i in POSITIONS)
        by_allele[row["allele"]].append(row)

    norm_rows, surface_rows, panel_rows, tie_rows = [], [], [], []
    surface_free_rank: dict[str, int] = {}
    for allele, rows in sorted(by_allele.items()):
        normalized = ranked(rows, lambda x: (-x["_normalized"], x["pair_id"]))
        raw = ranked(rows, lambda x: (-x["_raw_sum"], x["pair_id"]))
        norm_rank = {row["pair_id"]: index for index, row in enumerate(normalized, 1)}
        raw_rank = {row["pair_id"]: index for index, row in enumerate(raw, 1)}
        norm_rows.append({
            "allele": allele,
            "rank_correlation_normalized_vs_raw_sum": f"{pearson(list(norm_rank.values()), [raw_rank[key] for key in norm_rank]):.6f}",
            "top_16_overlap": len({x["pair_id"] for x in normalized[:16]} & {x["pair_id"] for x in raw[:16]}),
            "top_80_overlap": len({x["pair_id"] for x in normalized[:80]} & {x["pair_id"] for x in raw[:80]}),
        })
        for lead in LEADS:
            if lead in norm_rank:
                row = next(x for x in rows if x["pair_id"] == lead)
                norm_rows.append({
                    "allele": allele,
                    "rank_correlation_normalized_vs_raw_sum": "lead",
                    "top_16_overlap": f"{lead}: normalized rank {norm_rank[lead]}, raw-sum rank {raw_rank[lead]}, normalized score {row['_normalized']:.5f}, raw sum {row['_raw_sum']}",
                    "top_80_overlap": "raw mean has the same rank as raw sum",
                })
        archival = ranked(rows, lambda x: (-x["_normalized"], float(x["tcr_face_physicochemical_mismatch"]), -float(x["tcr_facing_sequence_identity"]), float(x["local_surface_percentile"]) if x["local_surface_percentile"] else math.inf, x["pair_id"]))
        no_surface = ranked(rows, lambda x: (-x["_normalized"], float(x["tcr_face_physicochemical_mismatch"]), -float(x["tcr_facing_sequence_identity"]), x["pair_id"]))
        for index, row in enumerate(no_surface, 1):
            surface_free_rank[row["pair_id"]] = index
        surface_rows.append({
            "allele": allele,
            "ordinal_positions_changed": sum(a["pair_id"] != b["pair_id"] for a, b in zip(archival, no_surface)),
            "top_16_overlap": len({x["pair_id"] for x in archival[:16]} & {x["pair_id"] for x in no_surface[:16]}),
            "top_40_overlap": len({x["pair_id"] for x in archival[:40]} & {x["pair_id"] for x in no_surface[:40]}),
            "top_80_overlap": len({x["pair_id"] for x in archival[:80]} & {x["pair_id"] for x in no_surface[:80]}),
        })
        blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            blocks[row["tcr_facing_blosum62_similarity"]].append(row)
        active = [block for block in blocks.values() if len(block) > 1]
        starting = sum(map(len, active))
        counts = {"physicochemical": 0, "identity": 0, "local_surface": 0, "pair_id": 0}
        for field, label in (("tcr_face_physicochemical_mismatch", "physicochemical"), ("tcr_facing_sequence_identity", "identity"), ("local_surface_percentile", "local_surface"), ("pair_id", "pair_id")):
            next_active = []
            for block in active:
                groups: dict[str, list[dict[str, str]]] = defaultdict(list)
                for row in block:
                    groups[row[field]].append(row)
                for subgroup in groups.values():
                    if len(subgroup) == 1:
                        counts[label] += 1
                    else:
                        next_active.append(subgroup)
            active = next_active
        tie_rows.append({"allele": allele, "primary_tied_rows": starting, **counts, "unresolved_after_pair_id": sum(map(len, active))})

    old_rank = {row["pair_id"]: int(row["primary_rank"]) for row in all_rows}
    for row in panel:
        prior, ablated = old_rank[row["pair_id"]], surface_free_rank[row["pair_id"]]
        if prior != ablated:
            panel_rows.append({"pair_id": row["pair_id"], "archival_rank": prior, "surface_free_rank": ablated, "absolute_shift": abs(prior - ablated)})

    record = defaultdict(dict)
    for row in records:
        record[(row["target_id"], row["arm_id"].rsplit("__", 1)[-1])][row["predictor"]] = row
    loo_rows = []
    for dropped in ("none", *PREDICTORS):
        keep = [name for name in PREDICTORS if name != dropped]
        passed = []
        for tier in tiers:
            pair_passes = True
            for arm in ("ebv", "self"):
                sequence, core = tier[f"{arm}_sequence"], tier[f"{arm}_core"]
                declared_start = sequence.index(core)
                for predictor in keep:
                    item = record[(tier["candidate_id"], arm)][predictor]
                    pair_passes &= float(item["percentile_rank"]) <= 20 and sequence.index(item["core"]) == declared_start
            if pair_passes:
                passed.append(tier["candidate_id"])
        loo_rows.append({"dropped_output": dropped, "retained_output_count": len(keep), "strict_20_percent_exact_core_pass_ids": ";".join(passed), "pass_count": len(passed)})

    write(OUTPUT / "normalization_ablation.csv", norm_rows)
    write(OUTPUT / "surface_key_ablation.csv", surface_rows)
    write(OUTPUT / "surface_key_panel_rank_changes.csv", panel_rows)
    write(OUTPUT / "sequential_tie_ablation.csv", tie_rows)
    write(OUTPUT / "leave_one_predictor_out.csv", loo_rows)
    print(OUTPUT)


if __name__ == "__main__":
    main()
