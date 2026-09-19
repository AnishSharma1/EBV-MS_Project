#!/usr/bin/env python3
"""Run the frozen anchor-preserving TCR-face permutation null."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


TCR_INDICES = (1, 2, 4, 6, 7)
EXPECTED_ALLELES = (
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*15:01",
)
CLAIM_BOUNDARY = (
    "Sequence-score calibration only; not evidence of endogenous presentation, "
    "pMHC stability, TCR recognition, cross-reactivity, molecular mimicry, or an MS mechanism."
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


def bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [math.nan] * len(values)
    running = 1.0
    count = len(values)
    for reverse_index in range(count - 1, -1, -1):
        original_index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[2]
    authoritative = Path(
        "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pairs",
        type=Path,
        default=workspace
        / "outputs/EBV_MS_ranking_tables_corrected_2026-09-03/categorized_tables/"
        "02_v3_structural_shortlist_and_rankings/04_all_6400_v3_pairs.csv",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=workspace
        / "outputs/EBV_MS_ranking_tables_corrected_2026-09-03/categorized_tables/"
        "01_current_stage1_evidence_review/03_presentation_conditioned_candidate_ranks.csv",
    )
    parser.add_argument(
        "--score-source",
        type=Path,
        default=authoritative / "src/hla2_positive_control_benchmark_v2.py",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_2/null_simulation_2026-09-14",
    )
    parser.add_argument("--replicates", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20_260_914)
    parser.add_argument("--batch-size", type=int, default=100)
    return parser.parse_args()


def load_score_module(source: Path):
    sys.path.insert(0, str(source.parent))
    import hla2_positive_control_benchmark_v2 as module  # type: ignore

    return module


def core_maps(rows: list[dict[str, str]], allele: str):
    subset = [row for row in rows if row["allele"] == allele]
    if len(subset) != 1600:
        raise ValueError(f"{allele}: expected 1600 rows, observed {len(subset)}")

    ebv = {}
    self_peptides = {}
    pairs = {}
    for row in subset:
        ebv_id = row["ebv_candidate_id"]
        self_id = row["self_candidate_id"]
        ebv_core = row["ebv_core_p1_p9"]
        self_core = row["self_core_p1_p9"]
        if ebv_id in ebv and ebv[ebv_id] != ebv_core:
            raise ValueError(f"{allele}: EBV ID maps to multiple cores: {ebv_id}")
        if self_id in self_peptides and self_peptides[self_id] != self_core:
            raise ValueError(f"{allele}: self ID maps to multiple cores: {self_id}")
        ebv[ebv_id] = ebv_core
        self_peptides[self_id] = self_core
        pairs[row["pair_id"]] = row

    if len(ebv) != 40 or len(self_peptides) != 40 or len(pairs) != 1600:
        raise ValueError(
            f"{allele}: expected 40 EBV, 40 self, 1600 unique pairs; "
            f"observed {len(ebv)}, {len(self_peptides)}, {len(pairs)}"
        )
    return subset, ebv, self_peptides, pairs


def main() -> None:
    args = parse_args()
    if args.replicates < 1 or args.batch_size < 1:
        raise ValueError("replicates and batch size must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    score_module = load_score_module(args.score_source)
    alphabet = score_module._BLOSUM62_ALPHABET
    aa_to_index = {aa: index for index, aa in enumerate(alphabet)}
    matrix = np.asarray(
        [[score_module._BLOSUM62[(left, right)] for right in alphabet] for left in alphabet],
        dtype=np.float64,
    )
    diagonal = np.diag(matrix)

    pair_rows = read_csv(args.pairs)
    audit_rows = read_csv(args.audit)
    if len(pair_rows) != 6400:
        raise ValueError(f"expected 6400 pair rows, observed {len(pair_rows)}")
    if len(audit_rows) != 8 or len({row['target_id'] for row in audit_rows}) != 8:
        raise ValueError("audit input must contain exactly eight unique target IDs")

    audits_by_allele = defaultdict(list)
    for row in audit_rows:
        audits_by_allele[row["allele"]].append(row)

    rng = np.random.default_rng(args.seed)
    fixed_null_scores: dict[str, np.ndarray] = {
        row["target_id"]: np.empty(args.replicates, dtype=np.float64) for row in audit_rows
    }
    allele_max_scores: dict[str, np.ndarray] = {
        allele: np.empty(args.replicates, dtype=np.float64) for allele in EXPECTED_ALLELES
    }
    null_summary_path = args.output_dir / "null_replicate_summary.csv"
    summary_fields = [
        "allele",
        "replicate",
        "target_id",
        "pair_id",
        "null_fixed_pair_score",
        "null_allele_max_score",
        "null_pairs_ge_observed",
        "observed_score",
    ]

    maximum_score_difference = 0.0
    with null_summary_path.open("w", newline="", encoding="utf-8") as summary_handle:
        summary_writer = csv.DictWriter(summary_handle, fieldnames=summary_fields)
        summary_writer.writeheader()

        for allele in EXPECTED_ALLELES:
            subset, ebv_map, self_map, pair_map = core_maps(pair_rows, allele)
            ebv_ids = sorted(ebv_map)
            self_ids = sorted(self_map)
            ebv_index = {value: index for index, value in enumerate(ebv_ids)}
            self_index = {value: index for index, value in enumerate(self_ids)}

            ebv_cores = [ebv_map[value] for value in ebv_ids]
            self_cores = [self_map[value] for value in self_ids]
            if any(len(core) != 9 for core in ebv_cores + self_cores):
                raise ValueError(f"{allele}: every declared core must contain exactly nine residues")

            ebv_exposed = np.asarray(
                [[aa_to_index[core[index]] for index in TCR_INDICES] for core in ebv_cores],
                dtype=np.int16,
            )
            self_exposed = np.asarray(
                [[aa_to_index[core[index]] for index in TCR_INDICES] for core in self_cores],
                dtype=np.int16,
            )

            audit_specs = []
            for audit in audits_by_allele[allele]:
                pair_id = audit["pair_id"]
                if pair_id not in pair_map:
                    raise ValueError(f"{allele}: audited pair absent from 1600-pair universe: {pair_id}")
                observed = pair_map[pair_id]
                recomputed = score_module.blosum62_similarity(
                    observed["ebv_core_p1_p9"],
                    observed["self_core_p1_p9"],
                    positions=TCR_INDICES,
                )
                stored = float(observed["tcr_facing_blosum62_similarity"])
                maximum_score_difference = max(maximum_score_difference, abs(recomputed - stored))
                if abs(recomputed - stored) > 1e-12:
                    raise ValueError(
                        f"{pair_id}: stored score {stored} does not reproduce; got {recomputed}"
                    )
                audit_specs.append(
                    {
                        "target_id": audit["target_id"],
                        "pair_id": pair_id,
                        "ebv_i": ebv_index[observed["ebv_candidate_id"]],
                        "self_i": self_index[observed["self_candidate_id"]],
                        "observed": stored,
                        "row": observed,
                    }
                )

            completed = 0
            while completed < args.replicates:
                batch = min(args.batch_size, args.replicates - completed)
                ebv_permutation = np.argsort(rng.random((batch, 40, 5)), axis=2)
                self_permutation = np.argsort(rng.random((batch, 40, 5)), axis=2)
                shuffled_ebv = np.take_along_axis(
                    np.broadcast_to(ebv_exposed, (batch, 40, 5)), ebv_permutation, axis=2
                )
                shuffled_self = np.take_along_axis(
                    np.broadcast_to(self_exposed, (batch, 40, 5)), self_permutation, axis=2
                )

                numerator = np.zeros((batch, 40, 40), dtype=np.float64)
                denominator = np.zeros((batch, 40, 40), dtype=np.float64)
                for position in range(5):
                    left = shuffled_ebv[:, :, position]
                    right = shuffled_self[:, :, position]
                    numerator += matrix[left[:, :, None], right[:, None, :]]
                    denominator += np.maximum(
                        diagonal[left][:, :, None], diagonal[right][:, None, :]
                    )
                scores = numerator / denominator
                maxima = scores.max(axis=(1, 2))
                allele_max_scores[allele][completed : completed + batch] = maxima

                for spec in audit_specs:
                    fixed = scores[:, spec["ebv_i"], spec["self_i"]]
                    fixed_null_scores[spec["target_id"]][completed : completed + batch] = fixed
                    exceedance_counts = (scores >= spec["observed"] - 1e-12).sum(axis=(1, 2))
                    for offset in range(batch):
                        summary_writer.writerow(
                            {
                                "allele": allele,
                                "replicate": completed + offset + 1,
                                "target_id": spec["target_id"],
                                "pair_id": spec["pair_id"],
                                "null_fixed_pair_score": f"{fixed[offset]:.12f}",
                                "null_allele_max_score": f"{maxima[offset]:.12f}",
                                "null_pairs_ge_observed": int(exceedance_counts[offset]),
                                "observed_score": f"{spec['observed']:.12f}",
                            }
                        )
                completed += batch

    pair_results = []
    raw_p_values = []
    for audit in audit_rows:
        allele = audit["allele"]
        _, _, _, pair_map = core_maps(pair_rows, allele)
        observed = pair_map[audit["pair_id"]]
        observed_score = float(observed["tcr_facing_blosum62_similarity"])
        fixed = fixed_null_scores[audit["target_id"]]
        maxima = allele_max_scores[allele]
        unadjusted = (1 + int((fixed >= observed_score - 1e-12).sum())) / (
            args.replicates + 1
        )
        max_t = (1 + int((maxima >= observed_score - 1e-12).sum())) / (
            args.replicates + 1
        )
        raw_p_values.append(unadjusted)
        pair_results.append(
            {
                "target_id": audit["target_id"],
                "allele": allele,
                "pair_id": audit["pair_id"],
                "ebv_protein": observed["ebv_protein"],
                "self_protein": observed["self_protein"],
                "ebv_core": observed["ebv_core_p1_p9"],
                "self_core": observed["self_core_p1_p9"],
                "observed_score": f"{observed_score:.12f}",
                "observed_hla_rank": observed["hla_rank"],
                "observed_score_rank": observed["hla_score_rank"],
                "observed_score_tie_size": observed["hla_score_tie_size"],
                "null_median": f"{np.quantile(fixed, 0.5):.12f}",
                "null_q025": f"{np.quantile(fixed, 0.025):.12f}",
                "null_q975": f"{np.quantile(fixed, 0.975):.12f}",
                "null_percentile_le_observed": f"{100.0 * np.mean(fixed <= observed_score + 1e-12):.6f}",
                "unadjusted_permutation_p": f"{unadjusted:.8f}",
                "maxT_fwer_p": f"{max_t:.8f}",
                "monte_carlo_se_unadjusted": f"{math.sqrt(unadjusted * (1-unadjusted)/(args.replicates+1)):.8f}",
                "replicates": args.replicates,
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )

    adjusted = bh_adjust(raw_p_values)
    for row, value in zip(pair_results, adjusted):
        row["bh_q_eight_audited_pairs"] = f"{value:.8f}"
        unadjusted = float(row["unadjusted_permutation_p"])
        max_t = float(row["maxT_fwer_p"])
        if max_t < 0.05:
            interpretation = "unusual_after_allele_wide_selection_adjustment"
        elif unadjusted < 0.05:
            interpretation = "pair_specific_only_not_allele_wide"
        else:
            interpretation = "not_unusual_under_constrained_null"
        row["prespecified_interpretation"] = interpretation

    pair_result_path = args.output_dir / "null_pair_results.csv"
    result_fields = list(pair_results[0])
    write_csv(pair_result_path, result_fields, pair_results)

    input_paths = [args.pairs, args.audit, args.score_source, Path(__file__).resolve()]
    lock = {
        "analysis": "anchor_preserving_tcr_face_permutation_null",
        "status": "complete",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "replicates_per_allele": args.replicates,
        "seed": args.seed,
        "batch_size": args.batch_size,
        "alleles": list(EXPECTED_ALLELES),
        "pairs_per_allele": 1600,
        "ebv_cores_per_allele": 40,
        "self_cores_per_allele": 40,
        "fixed_anchor_positions": ["P1", "P4", "P6", "P9"],
        "permuted_positions": ["P2", "P3", "P5", "P7", "P8"],
        "score": "existing normalized BLOSUM62 function",
        "comparison": "higher score is more similar; ties count in upper-tail exceedances",
        "claim_boundary": CLAIM_BOUNDARY,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "inputs": {str(path): sha256(path) for path in input_paths},
    }
    lock_path = args.output_dir / "analysis_lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    qa_lines = [
        "# Null simulation QA",
        "",
        "Status: PASS",
        "",
        f"- Observed rows: {len(pair_rows)} total; 1,600 per allele.",
        "- Unique discovery components: 40 EBV cores and 40 self cores per allele.",
        f"- Audited targets: {len(audit_rows)} unique pairs.",
        f"- Replicates: {args.replicates} completed for each of four alleles.",
        f"- Replicate-summary rows: {args.replicates * len(audit_rows)} expected and written.",
        f"- Maximum absolute observed-score reproduction difference: {maximum_score_difference:.3e}.",
        "- Score reproduction tolerance: 1e-12.",
        "- Anchor positions remained fixed; only the five designated TCR-facing positions were permuted.",
        "- Each unique peptide was shuffled once per replicate and reused across its 40 pairings.",
        "- Upper-tail comparisons included ties.",
        "",
        f"Claim boundary: {CLAIM_BOUNDARY}",
    ]
    qa_path = args.output_dir / "analysis_qa.md"
    qa_path.write_text("\n".join(qa_lines) + "\n", encoding="utf-8")

    summary_lines = [
        "# Constrained-null simulation results",
        "",
        f"Replicates: {args.replicates} per allele. Seed: {args.seed}.",
        "",
        "| Target | Allele | Observed rank | Observed score | Pair-specific p | BH q | MaxT p | Interpretation |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in pair_results:
        summary_lines.append(
            "| {target_id} | {allele} | {observed_hla_rank} | {observed_score} | "
            "{unadjusted_permutation_p} | {bh_q_eight_audited_pairs} | {maxT_fwer_p} | "
            "{prespecified_interpretation} |".format(**row)
        )
    summary_lines.extend(["", f"Claim boundary: {CLAIM_BOUNDARY}"])
    results_summary_path = args.output_dir / "RESULTS_SUMMARY.md"
    results_summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest_targets = [
        lock_path,
        pair_result_path,
        null_summary_path,
        qa_path,
        results_summary_path,
        Path(__file__).resolve(),
    ]
    manifest_rows = [
        {"path": str(path.resolve()), "sha256": sha256(path)} for path in manifest_targets
    ]
    write_csv(args.output_dir / "SHA256SUMS.csv", ["path", "sha256"], manifest_rows)

    print(json.dumps({
        "status": "PASS",
        "output_dir": str(args.output_dir),
        "replicates_per_allele": args.replicates,
        "maximum_score_difference": maximum_score_difference,
        "pair_results": pair_results,
    }, indent=2))


if __name__ == "__main__":
    main()
