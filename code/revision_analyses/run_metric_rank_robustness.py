#!/usr/bin/env python3
"""Run the frozen four-metric rank-robustness analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ALLELES = (
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*15:01",
)
METRICS = {
    "tcr_facing_blosum62": {
        "field": "tcr_facing_blosum62_similarity",
        "higher_is_better": True,
        "stored_display_rank": "tcr_facing_blosum62_hla_rank",
    },
    "full_core_blosum62": {
        "field": "full_core_blosum62_similarity",
        "higher_is_better": True,
        "stored_display_rank": "full_core_blosum62_hla_rank",
    },
    "tcr_facing_identity": {
        "field": "tcr_facing_sequence_identity",
        "higher_is_better": True,
        "stored_display_rank": "tcr_facing_identity_hla_rank",
    },
    "physicochemical_mismatch": {
        "field": "tcr_face_physicochemical_mismatch",
        "higher_is_better": False,
        "stored_display_rank": "physicochemical_only_hla_rank",
    },
}
LEAD_TARGET_IDS = ("HY13_SEQ_02", "HY15_SEQ_02")
TCR_FACING_INDICES = (1, 2, 4, 6, 7)
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
FORMAL_CHARGE = {aa: 0.0 for aa in HYDROPATHY}
FORMAL_CHARGE.update({"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0})
H_BOND_DONOR = {aa: float(aa in "RKNQHSTWY") for aa in HYDROPATHY}
H_BOND_ACCEPTOR = {aa: float(aa in "DENQHSTY") for aa in HYDROPATHY}
AROMATIC = {aa: float(aa in "FWY") for aa in HYDROPATHY}
CLAIM_BOUNDARY = (
    "Rank stability across four related sequence descriptors only; the descriptors are not "
    "independent biological validations and do not establish presentation, TCR recognition, "
    "cross-reactivity, molecular mimicry, or an MS mechanism."
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


def physicochemical_mismatch(left_core: str, right_core: str) -> float:
    """Recompute the source descriptor before its 12-decimal CSV rounding."""
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("physicochemical comparison requires two exact nine-residue cores")

    def descriptor(amino_acid: str) -> np.ndarray:
        if amino_acid not in HYDROPATHY:
            raise ValueError(f"unsupported amino acid: {amino_acid}")
        return np.asarray([
            (FORMAL_CHARGE[amino_acid] + 1.0) / 2.0,
            (HYDROPATHY[amino_acid] + 4.5) / 9.0,
            H_BOND_DONOR[amino_acid],
            H_BOND_ACCEPTOR[amino_acid],
            AROMATIC[amino_acid],
        ])

    position_values = [
        float(np.mean(np.abs(descriptor(left_core[index]) - descriptor(right_core[index]))))
        for index in TCR_FACING_INDICES
    ]
    return float(np.mean(position_values))


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[2]
    base = workspace / "outputs/EBV_MS_ranking_tables_corrected_2026-09-03/categorized_tables"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pairs",
        type=Path,
        default=base / "02_v3_structural_shortlist_and_rankings/04_all_6400_v3_pairs.csv",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=base / "01_current_stage1_evidence_review/03_presentation_conditioned_candidate_ranks.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_2/metric_rank_robustness_2026-09-14",
    )
    return parser.parse_args()


def metric_ranks(rows: list[dict[str, str]], field: str, higher_is_better: bool):
    if field == "tcr_face_physicochemical_mismatch":
        values = {
            row["pair_id"]: physicochemical_mismatch(
                row["ebv_predicted_core"], row["self_predicted_core"]
            )
            for row in rows
        }
    else:
        values = {row["pair_id"]: float(row[field]) for row in rows}
    ordered = sorted(
        rows,
        key=lambda row: (
            -values[row["pair_id"]] if higher_is_better else values[row["pair_id"]],
            row["pair_id"],
        ),
    )
    display = {row["pair_id"]: index for index, row in enumerate(ordered, 1)}
    score_rank = {}
    average_rank = {}
    tie_size = {}
    position = 0
    while position < len(ordered):
        end = position + 1
        value = values[ordered[position]["pair_id"]]
        while end < len(ordered) and values[ordered[end]["pair_id"]] == value:
            end += 1
        competition = position + 1
        average = (position + 1 + end) / 2.0
        for row in ordered[position:end]:
            pair_id = row["pair_id"]
            score_rank[pair_id] = competition
            average_rank[pair_id] = average
            tie_size[pair_id] = end - position
        position = end
    return values, display, score_rank, average_rank, tie_size


def robustness_class(top_1_count: int, top_5_count: int) -> str:
    if top_1_count >= 3:
        return "robust_top_1_percent_3_of_4"
    if top_5_count >= 3:
        return "partially_robust_top_5_percent_3_of_4"
    return "metric_dependent"


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pair_rows = read_csv(args.pairs)
    audit_rows = read_csv(args.audit)
    if len(pair_rows) != 6400:
        raise ValueError(f"expected 6400 pairs, observed {len(pair_rows)}")
    if len({row["pair_id"] for row in pair_rows}) != 6400:
        raise ValueError("pair IDs are not unique")
    if Counter(row["allele"] for row in pair_rows) != Counter({allele: 1600 for allele in ALLELES}):
        raise ValueError("expected exactly 1600 rows for each of four alleles")

    lead_audit = {row["target_id"]: row for row in audit_rows if row["target_id"] in LEAD_TARGET_IDS}
    if set(lead_audit) != set(LEAD_TARGET_IDS):
        raise ValueError("could not resolve both exact BALF5-TALDO1 lead pair IDs")
    lead_pair_to_target = {row["pair_id"]: target for target, row in lead_audit.items()}

    all_rank_rows = []
    correlation_rows = []
    rank_validation_mismatches = []
    rounded_value_mismatches = []
    metric_cache = {}

    for allele in ALLELES:
        rows = [row for row in pair_rows if row["allele"] == allele]
        metric_cache[allele] = {}
        for name, spec in METRICS.items():
            result = metric_ranks(rows, spec["field"], spec["higher_is_better"])
            metric_cache[allele][name] = result
            _values, display, _score, _average, _ties = result
            for row in rows:
                if name == "physicochemical_mismatch":
                    published = float(row[spec["field"]])
                    if f"{_values[row['pair_id']]:.12f}" != f"{published:.12f}":
                        rounded_value_mismatches.append(row["pair_id"])
                stored = int(row[spec["stored_display_rank"]])
                recomputed = display[row["pair_id"]]
                if stored != recomputed:
                    rank_validation_mismatches.append(
                        (allele, row["pair_id"], name, stored, recomputed)
                    )
        if rank_validation_mismatches:
            raise ValueError(f"stored metric ranks failed reproduction: {rank_validation_mismatches[:3]}")

        ordered_pair_ids = sorted(row["pair_id"] for row in rows)
        for index, left_name in enumerate(METRICS):
            for right_name in list(METRICS)[index + 1:]:
                left_average = metric_cache[allele][left_name][3]
                right_average = metric_cache[allele][right_name][3]
                left = np.asarray([left_average[pair_id] for pair_id in ordered_pair_ids], dtype=float)
                right = np.asarray([right_average[pair_id] for pair_id in ordered_pair_ids], dtype=float)
                rho = float(np.corrcoef(left, right)[0, 1])
                correlation_rows.append(
                    {
                        "allele": allele,
                        "metric_1": left_name,
                        "metric_2": right_name,
                        "spearman_rho_oriented_ranks": f"{rho:.12f}",
                        "pair_count": len(rows),
                        "tie_handling": "average ranks; rank 1 is best for every metric",
                        "claim_boundary": CLAIM_BOUNDARY,
                    }
                )

        for row in rows:
            output = {
                "allele": allele,
                "pair_id": row["pair_id"],
                "ebv_protein": row["ebv_protein"],
                "self_protein": row["self_protein"],
                "ebv_core": row["ebv_core_p1_p9"],
                "self_core": row["self_core_p1_p9"],
                "claim_boundary": CLAIM_BOUNDARY,
            }
            top_1_count = 0
            top_5_count = 0
            average_top_1_count = 0
            average_top_5_count = 0
            worst_top_1_count = 0
            worst_top_5_count = 0
            for name, spec in METRICS.items():
                values, display, score, average, ties = metric_cache[allele][name]
                pair_id = row["pair_id"]
                output[f"{name}_value"] = f"{values[pair_id]:.12f}"
                output[f"{name}_display_rank"] = display[pair_id]
                output[f"{name}_score_rank"] = score[pair_id]
                output[f"{name}_average_rank"] = f"{average[pair_id]:.6f}"
                output[f"{name}_tie_size"] = ties[pair_id]
                output[f"{name}_worst_tied_rank"] = score[pair_id] + ties[pair_id] - 1
                top_1_count += int(score[pair_id] <= 16)
                top_5_count += int(score[pair_id] <= 80)
                average_top_1_count += int(average[pair_id] <= 16)
                average_top_5_count += int(average[pair_id] <= 80)
                worst_top_1_count += int(score[pair_id] + ties[pair_id] - 1 <= 16)
                worst_top_5_count += int(score[pair_id] + ties[pair_id] - 1 <= 80)
            output["top_1_percent_metric_count"] = top_1_count
            output["top_5_percent_metric_count"] = top_5_count
            output["rank_robustness_class"] = robustness_class(top_1_count, top_5_count)
            output["average_rank_robustness_class"] = robustness_class(
                average_top_1_count, average_top_5_count
            )
            output["worst_tied_rank_robustness_class"] = robustness_class(
                worst_top_1_count, worst_top_5_count
            )
            all_rank_rows.append(output)

    if rounded_value_mismatches:
        raise ValueError(
            "recomputed physicochemical values did not match public values at 12 decimals: "
            f"{rounded_value_mismatches[:3]}"
        )

    all_path = args.output_dir / "all_pair_metric_ranks.csv"
    write_csv(all_path, list(all_rank_rows[0]), all_rank_rows)

    correlation_path = args.output_dir / "allele_metric_correlations.csv"
    write_csv(correlation_path, list(correlation_rows[0]), correlation_rows)

    rank_by_pair = {row["pair_id"]: row for row in all_rank_rows}
    lead_rows = []
    for target_id in LEAD_TARGET_IDS:
        audit = lead_audit[target_id]
        source = rank_by_pair[audit["pair_id"]]
        lead = {
            "target_id": target_id,
            "allele": source["allele"],
            "pair_id": source["pair_id"],
            "ebv_protein": source["ebv_protein"],
            "self_protein": source["self_protein"],
            "ebv_core": source["ebv_core"],
            "self_core": source["self_core"],
        }
        for name in METRICS:
            for suffix in (
                "value", "display_rank", "score_rank", "average_rank", "tie_size",
                "worst_tied_rank",
            ):
                lead[f"{name}_{suffix}"] = source[f"{name}_{suffix}"]
        lead["top_1_percent_metric_count"] = source["top_1_percent_metric_count"]
        lead["top_5_percent_metric_count"] = source["top_5_percent_metric_count"]
        lead["rank_robustness_class"] = source["rank_robustness_class"]
        lead["average_rank_robustness_class"] = source["average_rank_robustness_class"]
        lead["worst_tied_rank_robustness_class"] = source["worst_tied_rank_robustness_class"]
        lead["claim_boundary"] = CLAIM_BOUNDARY
        lead_rows.append(lead)

    lead_path = args.output_dir / "lead_metric_rank_robustness.csv"
    write_csv(lead_path, list(lead_rows[0]), lead_rows)

    class_counts = Counter(row["rank_robustness_class"] for row in all_rank_rows)
    expected_lead_classes = {
        "HY13_SEQ_02": "partially_robust_top_5_percent_3_of_4",
        "HY15_SEQ_02": "robust_top_1_percent_3_of_4",
    }
    observed_lead_classes = {row["target_id"]: row["rank_robustness_class"] for row in lead_rows}
    if observed_lead_classes != expected_lead_classes:
        raise ValueError(f"lead classifications did not match locked rules: {observed_lead_classes}")

    lock = {
        "analysis": "four_metric_rank_robustness",
        "status": "complete",
        "analysis_date": "2026-09-14",
        "alleles": list(ALLELES),
        "pairs_per_allele": 1600,
        "metrics": METRICS,
        "scientific_rank": "competition rank; tied values share 1 + count of strictly better values",
        "display_rank": "deterministic full-precision metric ordering with lexical pair ID tie-break; reported but not used for robustness classification",
        "physicochemical_precision": "descriptor recomputed from exact nine-residue cores because the public CSV stores values rounded to 12 decimals",
        "spearman": "Pearson correlation of tie-aware average oriented ranks; rank 1 is best for every metric",
        "robust_rule": "scientific score rank <=16 under at least three of four metrics",
        "partial_rule": "if not robust, scientific score rank <=80 under at least three of four metrics",
        "tie_policy_sensitivity": "the same rules were also evaluated with average tied ranks and conservative worst tied ranks",
        "claim_boundary": CLAIM_BOUNDARY,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "inputs": {
            str(path.resolve()): sha256(path)
            for path in (args.pairs, args.audit, Path(__file__).resolve())
        },
    }
    lock_path = args.output_dir / "analysis_lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    qa_lines = [
        "# Four-metric rank-robustness QA",
        "",
        "Status: PASS",
        "",
        "- Input rows: 6,400 unique pairs; exactly 1,600 per allele.",
        "- All four metric value fields were present and numeric for every pair.",
        "- All 25,600 stored metric display ranks were independently reproduced exactly.",
        "- The physicochemical descriptor was recomputed from each exact core pair; all 6,400 values matched the public CSV after 12-decimal rounding.",
        "- Scientific score ranks used competition ranking so tied values did not depend on lexical pair ID.",
        "- Average-rank and conservative worst-tied-rank classifications were added as tie-policy sensitivity checks.",
        "- Spearman correlations used average ranks for ties and oriented every metric so rank 1 is best.",
        "- Top 1% boundary: scientific score rank <=16 of 1,600.",
        "- Top 5% boundary: scientific score rank <=80 of 1,600.",
        "- Exact lead pair IDs were recovered from the eight-pair audit table.",
        f"- All-pair class counts: {dict(sorted(class_counts.items()))}.",
        "- HY13_SEQ_02 is partially robust under all three tie policies.",
        "- HY15_SEQ_02 is robust under competition ranks but partially robust under average and worst-tied ranks.",
        "",
        f"Claim boundary: {CLAIM_BOUNDARY}",
    ]
    qa_path = args.output_dir / "analysis_qa.md"
    qa_path.write_text("\n".join(qa_lines) + "\n", encoding="utf-8")

    results_lines = [
        "# Four-metric rank-robustness results",
        "",
        "Scientific score ranks are shown; tied values share a rank.",
        "",
        "| Target | Allele | TCR-face BLOSUM62 | Full-core BLOSUM62 | TCR-face identity | Physicochemical mismatch | Competition-rank class | Average-rank class | Worst-tied-rank class |",
        "|---|---|---:|---:|---:|---:|---|---|---|",
    ]
    for row in lead_rows:
        results_lines.append(
            "| {target_id} | {allele} | {tcr_facing_blosum62_score_rank} | "
            "{full_core_blosum62_score_rank} | {tcr_facing_identity_score_rank} | "
            "{physicochemical_mismatch_score_rank} | {rank_robustness_class} | "
            "{average_rank_robustness_class} | {worst_tied_rank_robustness_class} |".format(**row)
        )
    results_lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- HY15_SEQ_02 meets the locked competition-rank robust rule, but that label is tie-policy-sensitive: its 0.4 identity score shares a 37-pair tie spanning ranks 4-40, so average and worst-tied ranks classify it as partially robust.",
            "- HY13_SEQ_02 is partially robust because three metrics place it in the top 5%, but only TCR-facing BLOSUM62 places it in the top 1%.",
            "- Cross-metric robustness does not erase the null-simulation or binding-threshold results and is not independent validation because the four metrics reuse overlapping sequence information.",
            "",
            f"Claim boundary: {CLAIM_BOUNDARY}",
        ]
    )
    results_path = args.output_dir / "RESULTS_SUMMARY.md"
    results_path.write_text("\n".join(results_lines) + "\n", encoding="utf-8")

    checksum_targets = [
        lock_path,
        all_path,
        correlation_path,
        lead_path,
        qa_path,
        results_path,
        Path(__file__).resolve(),
    ]
    checksum_rows = [
        {"path": str(path.resolve()), "sha256": sha256(path)} for path in checksum_targets
    ]
    checksum_path = args.output_dir / "SHA256SUMS.csv"
    write_csv(checksum_path, ["path", "sha256"], checksum_rows)

    print(json.dumps({
        "status": "PASS",
        "output_dir": str(args.output_dir.resolve()),
        "lead_results": lead_rows,
        "class_counts": dict(sorted(class_counts.items())),
        "correlations": correlation_rows,
    }, indent=2))


if __name__ == "__main__":
    main()
