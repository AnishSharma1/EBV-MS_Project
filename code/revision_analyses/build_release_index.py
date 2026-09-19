#!/usr/bin/env python3
"""Build and verify the round-3 computational-additions release index."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ROUND = ROOT / "output/revision_round_3"
OUT = ROUND / "COMPUTATIONAL_ADDITIONS_RELEASE_2026-09-15"
PACKAGES = [
    "computational_tightening_2026-09-15",
    "proteome_sequence_null_2026-09-15",
    "proteome_binding_filter_2026-09-15",
    "proteome_binding_consensus_2026-09-15",
    "register_uncertainty_2026-09-15",
    "solved_structure_benchmark_inventory_2026-09-15",
    "evidence_and_control_audit_2026-09-15",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify_manifest(package):
    manifest = package / "SHA256SUMS.csv"
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    failures = []
    for row in rows:
        relative = row.get("relative_path") or row.get("path")
        path = package / relative
        if not path.exists():
            failures.append(f"missing:{relative}")
        elif digest(path) != row["sha256"]:
            failures.append(f"checksum:{relative}")
        elif row.get("bytes") and path.stat().st_size != int(row["bytes"]):
            failures.append(f"size:{relative}")
    return len(rows), failures


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    package_rows = []
    all_failures = []
    for name in PACKAGES:
        package = ROUND / name
        count, failures = verify_manifest(package)
        all_failures.extend(f"{name}/{failure}" for failure in failures)
        package_rows.append({
            "package": name,
            "absolute_path": str(package),
            "manifest_entry_count": count,
            "checksum_status": "PASS" if not failures else "FAIL",
            "package_manifest_sha256": digest(package / "SHA256SUMS.csv"),
        })
    write_csv(OUT / "package_index.csv", package_rows)

    status_rows = [
        ("Exact scoring equation, positions, normalization, negative-score handling, worked examples", "complete", "computational_tightening_2026-09-15/scoring_formula.json"),
        ("Complete 6,400-pair score/ablation table", "complete", "computational_tightening_2026-09-15/score_ablation_all_6400.csv"),
        ("Threshold-sensitivity analysis", "complete", "computational_tightening_2026-09-15/figure_threshold_sensitivity.svg"),
        ("Matched empirical decoys", "complete_post_hoc_descriptive", "computational_tightening_2026-09-15/matched_empirical_decoys.csv"),
        ("Independent +/-1 viral/self register scenarios", "complete_sensitivity_not_posterior", "register_uncertainty_2026-09-15/lead_register_uncertainty_summary.csv"),
        ("Corrected structural ensemble statistics and full dependent 5x5 matrices", "complete", "computational_tightening_2026-09-15/structural_ensemble_corrected_with_medoids.csv"),
        ("Model-level AF3 confidence table", "complete", "computational_tightening_2026-09-15/af3_model_level_confidence_65.csv"),
        ("General solved HLA-II structure discovery inventory", "inventory_complete_validation_pending_manual_register_ground_truth", "solved_structure_benchmark_inventory_2026-09-15/benchmark_status.json"),
        ("Reviewed-human-proteome sequence null, 11.25 million 9-mers", "complete", "proteome_sequence_null_2026-09-15/proteome_sequence_null_summary.csv"),
        ("Top-1,000-per-lead two-predictor HLA-II binding filter", "complete_sensitivity_not_full_proteome_binding_scan", "proteome_binding_consensus_2026-09-15/two_predictor_consensus_summary.csv"),
        ("HLA Ligand Atlas and submitted PXD068488 table audit", "complete_processed_source_audit", "evidence_and_control_audit_2026-09-15/lead_immunopeptidomics_summary.csv"),
        ("PXD068488 raw-spectrum reanalysis", "prepared_not_run", "evidence_and_control_audit_2026-09-15/raw_immunopeptidomics_reanalysis_lock.json"),
        ("TCR-facing score control validation", "audited_provisional_three_systems_definitive_gate_blocked", "evidence_and_control_audit_2026-09-15/tcr_score_validation_status.json"),
        ("Prospective experiment lock", "complete_preexperiment_prepared_not_submitted", "computational_tightening_2026-09-15/prospective_preexperiment_lock.json"),
        ("Direct residue-to-pocket ground-truth caller across solved structures", "pending_manual_ground_truth_and_blind_benchmark", "solved_structure_benchmark_inventory_2026-09-15/benchmark_status.json"),
        ("Experimental binding, register mapping, and cross-reactivity validation", "pending_experiment", "computational_tightening_2026-09-15/prospective_preexperiment_lock.json"),
    ]
    status = [{"requested_addition": a, "status": b, "primary_artifact": c} for a, b, c in status_rows]
    write_csv(OUT / "requested_additions_status.csv", status)

    qa = {
        "release_status": "PASS" if not all_failures else "FAIL",
        "verified_package_count": len(PACKAGES),
        "checksum_failures": all_failures,
        "completed_or_audited_items": sum(row["status"].startswith(("complete", "audited", "inventory")) for row in status),
        "prepared_not_run_items": sum("prepared_not_run" in row["status"] for row in status),
        "pending_items": sum(row["status"].startswith("pending") for row in status),
        "interpretation": (
            "The defensible pre-experiment computational additions are complete. Tasks that require manual "
            "coordinate ground truth, raw-spectrum reprocessing, or wet-lab data remain explicitly separate."
        ),
    }
    (OUT / "RELEASE_QA.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    (OUT / "README.md").write_text(
        "# V8 computational additions release\n\n"
        "This directory is the entry point for the reviewer-requested computational work. "
        "`requested_additions_status.csv` distinguishes completed analyses from inventories, "
        "prepared external work, and experiments that cannot honestly be replaced by computation.\n\n"
        "All seven component packages have independent checksum manifests. `RELEASE_QA.json` "
        "records the cross-package verification result. No manuscript prose was modified.\n",
        encoding="utf-8",
    )
    print(json.dumps(qa, indent=2))
    if all_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
