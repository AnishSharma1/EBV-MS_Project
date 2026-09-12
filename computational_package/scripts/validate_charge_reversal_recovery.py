"""Fail-closed acceptance checks for the charge-reversal recovery package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from charge_reversal_recovery import make_af3_ablation_jobs, primary_ten, validate_af3_ablation_jobs


ROOT = Path(__file__).resolve().parents[1]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate(package: Path, require_pandora_candidates: bool = False, require_af3_recovery: bool = False) -> str:
    audit = rows(package / "baseline_diagnosis" / "baseline_af3_peptide_template_audit_10.csv")
    if len(audit) != 10 or any(not row["peptide_template_pdb_ids"] for row in audit):
        raise ValueError("baseline peptide-template audit is incomplete")

    jobs = json.loads((package / "alphafold_ablation" / "af3_recovery_batch_30_jobs.json").read_text(encoding="utf-8"))
    validate_af3_ablation_jobs(jobs)
    if jobs != make_af3_ablation_jobs():
        raise ValueError("AlphaFold recovery JSON differs from deterministic generator")
    af3_status = rows(package / "alphafold_ablation" / "AF3_RECOVERY_STATUS_30.csv")
    if len(af3_status) != 30:
        raise ValueError("AlphaFold recovery status must cover 30 jobs")
    af3_summary_path = package / "alphafold_ablation" / "analysis" / "AF3_RECOVERY_SUMMARY.json"
    if af3_summary_path.exists():
        af3_summary = json.loads(af3_summary_path.read_text(encoding="utf-8"))
        if af3_summary["jobs"] != 30 or af3_summary["models"] != 150 or not af3_summary["all_exact_sequence_qc_pass"]:
            raise ValueError("AlphaFold recovery completeness or sequence QC failure")
        if af3_summary["predeclared_seed_matches"] != 29 or af3_summary["predeclared_seed_mismatches"] != 1:
            raise ValueError("the recorded AlphaFold seed discrepancy changed")
        if {row["status"] for row in af3_status} != {"returned_complete_qc_pass", "returned_complete_seed_mismatch"}:
            raise ValueError("AlphaFold returned-job status is inconsistent")
        retry_jobs = json.loads((package / "alphafold_ablation" / "af3_retry_seed_mismatch_jobs.json").read_text(encoding="utf-8"))
        if len(retry_jobs) != 1 or retry_jobs[0]["name"] != af3_summary["seed_mismatch_job_names"][0]:
            raise ValueError("AlphaFold one-job retry file does not match the failed seed gate")
    elif require_af3_recovery:
        raise ValueError("completed AlphaFold recovery analysis is required")
    elif {row["status"] for row in af3_status} != {"prepared_not_submitted"}:
        raise ValueError("external AlphaFold status must remain explicit")

    candidates = rows(package / "pandora" / "candidate_runs_60_requesting_1200_models.csv")
    calibrations = rows(package / "pandora" / "calibration_runs_6_requesting_120_models.csv")
    if len(candidates) != 60 or sum(int(row["requested_models"]) for row in candidates) != 1200:
        raise ValueError("PANDORA candidate manifest cardinality mismatch")
    if len(calibrations) != 6 or sum(int(row["requested_models"]) for row in calibrations) != 120:
        raise ValueError("PANDORA calibration manifest cardinality mismatch")
    if any(row["forced_template_pdb"] == row["experimental_target_pdb"] for row in calibrations):
        raise ValueError("leave-one-template-out leakage in calibration manifest")

    template_manifest = rows(package / "pandora" / "template_manifest.csv")
    if {row["pdb_id"] for row in template_manifest} != {"1BX2", "6CQQ"}:
        raise ValueError("template identities mismatch")
    for row in template_manifest:
        for prefix in ("raw", "prepared"):
            path = package / row[f"{prefix}_relative_path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != row[f"{prefix}_sha256"]:
                raise ValueError(f"template checksum mismatch: {path}")

    gate_path = package / "pandora" / "results" / "CALIBRATION_GATE.json"
    if gate_path.exists():
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        if gate["models_scored"] != 120 or gate["calibration_status"] != "pass":
            raise ValueError("PANDORA calibration gate did not pass")
    elif require_pandora_candidates:
        raise ValueError("PANDORA calibration results are required")

    candidate_summary_path = package / "pandora" / "results" / "PANDORA_CANDIDATE_SUMMARY.json"
    if candidate_summary_path.exists():
        summary = json.loads(candidate_summary_path.read_text(encoding="utf-8"))
        if summary["runs"] != 60 or summary["models"] != 1200 or not summary["all_exact_sequence_qc_pass"]:
            raise ValueError("PANDORA candidate ensemble acceptance failure")
    elif require_pandora_candidates:
        raise ValueError("completed PANDORA candidate analysis is required")

    consensus = rows(package / "consensus" / "consensus_classification_8.csv")
    if len(consensus) != 8:
        raise ValueError("consensus table must cover eight mutations")
    if af3_summary_path.exists():
        counts = {value: sum(row["classification"] == value for row in consensus) for value in {row["classification"] for row in consensus}}
        if counts != {"register_confounded": 7, "not_evaluable": 1}:
            raise ValueError(f"unexpected fail-closed consensus counts: {counts}")
    elif {row["classification"] for row in consensus} != {"not_evaluable"}:
        raise ValueError("consensus cannot advance before external AlphaFold recovery results")

    software = rows(package / "software_manifest.csv")
    for row in software:
        path = ROOT / row["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"software provenance mismatch: {row['path']}")

    checksums = rows(package / "SHA256SUMS.csv")
    for row in checksums:
        path = package / row["relative_path"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"package checksum mismatch: {row['relative_path']}")
    return (
        "Package-integrity QA passed: frozen 10-job template audit; 30-job/150-model AF3 return; "
        "six reciprocal PANDORA calibration runs; 60 candidate runs/1,200 requested models; "
        "eight fail-closed consensus rows; all registered checksums matched. "
        "Scientific acceptance remains failed for one predeclared AlphaFold seed mismatch."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--require-pandora-candidates", action="store_true")
    parser.add_argument("--require-af3-recovery", action="store_true")
    args = parser.parse_args()
    print(validate(args.package.resolve(), args.require_pandora_candidates, args.require_af3_recovery))


if __name__ == "__main__":
    main()
