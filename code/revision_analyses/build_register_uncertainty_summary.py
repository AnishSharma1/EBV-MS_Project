#!/usr/bin/env python3
"""Summarize the nine independent +/-1 register scenarios for each lead."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/revision_round_3/computational_tightening_2026-09-15/register_shift_detail_all_6400.csv"
OUT = ROOT / "output/revision_round_3/register_uncertainty_2026-09-15"
PAIRS = {
    "HY13_SEQ_02": "HLA-DRB1*13:03|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0108_0122",
    "HY15_SEQ_02": "HLA-DRB1*15:01|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0216_0230",
}
BOUNDARY = (
    "Equal-weight register-scenario sensitivity summary, not a calibrated posterior probability. "
    "It does not establish binding, presentation, TCR recognition, or cross-reactivity."
)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    detail = []
    summary = []
    for target_id, pair_id in PAIRS.items():
        subset = [row for row in rows if row["pair_id"] == pair_id]
        if len(subset) != 9:
            raise ValueError(f"{target_id}: expected 9 register scenarios, found {len(subset)}")
        ranks = [int(row["rank_against_frozen_declared_universe"]) for row in subset]
        scores = [float(row["score"]) for row in subset]
        declared = next(row for row in subset if row["left_register_delta"] == "0" and row["right_register_delta"] == "0")
        for row in subset:
            detail.append({"target_id": target_id, **row, "scenario_weight": "0.111111111111", "claim_boundary": BOUNDARY})
        summary.append({
            "target_id": target_id,
            "allele": subset[0]["allele"],
            "scenario_count": 9,
            "declared_rank": declared["rank_against_frozen_declared_universe"],
            "best_rank": min(ranks),
            "median_rank": statistics.median(ranks),
            "worst_rank": max(ranks),
            "rank_q1": statistics.quantiles(ranks, n=4, method="inclusive")[0],
            "rank_q3": statistics.quantiles(ranks, n=4, method="inclusive")[2],
            "uniform_scenario_fraction_top_1": f"{sum(rank <= 1 for rank in ranks) / 9:.12g}",
            "uniform_scenario_fraction_top_5": f"{sum(rank <= 5 for rank in ranks) / 9:.12g}",
            "uniform_scenario_fraction_top_20": f"{sum(rank <= 20 for rank in ranks) / 9:.12g}",
            "uniform_scenario_fraction_top_100": f"{sum(rank <= 100 for rank in ranks) / 9:.12g}",
            "score_min": f"{min(scores):.12g}",
            "score_median": f"{statistics.median(scores):.12g}",
            "score_max": f"{max(scores):.12g}",
            "declared_is_best": int(declared["rank_against_frozen_declared_universe"]) == min(ranks),
            "claim_boundary": BOUNDARY,
        })
    write_csv(OUT / "lead_register_scenarios.csv", detail)
    write_csv(OUT / "lead_register_uncertainty_summary.csv", summary)
    lock = {
        "status": "complete_nine_scenario_sensitivity",
        "scenario_definition": "viral and self parent peptides shifted independently by -1, 0, or +1",
        "weighting": "uniform 1/9 for descriptive scenario fractions",
        "calibrated_probabilities": False,
        "source": str(SOURCE),
        "source_sha256": digest(SOURCE),
        "claim_boundary": BOUNDARY,
    }
    (OUT / "analysis_lock.json").write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    (OUT / "README.md").write_text(
        "# Register-uncertainty sensitivity\n\n"
        "Each lead is evaluated under all nine independent combinations of viral and self "
        "register shifts (-1, 0, +1). The reported top-k fractions use equal scenario weights "
        "for transparency. They are sensitivity summaries, not Bayesian posterior probabilities.\n",
        encoding="utf-8",
    )
    checksum_path = OUT / "SHA256SUMS.csv"
    checksums = []
    for path in sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_path):
        checksums.append({"relative_path": str(path.relative_to(OUT)), "sha256": digest(path), "bytes": path.stat().st_size})
    write_csv(checksum_path, checksums)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
