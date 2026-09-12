"""Validate identity, status, layout, and checksums for the charge-reversal package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from charge_reversal_binding_pilot import primary_panel, validate_panel


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    args = parser.parse_args()
    package = args.package.resolve()

    expected_panel = primary_panel()
    validate_panel(expected_panel)
    expected_by_id = {row["sample_id"]: row for row in expected_panel}
    manifest = read_csv(package / "peptide_manifest.csv")
    if len(manifest) != 14 or {row["sample_id"] for row in manifest} != set(expected_by_id):
        raise ValueError("manifest identity/count mismatch")
    for row in manifest:
        expected = expected_by_id[row["sample_id"]]
        if row["sequence"] != expected["sequence"] or row["expected_core"] != expected["expected_core"]:
            raise ValueError(f"manifest sequence/core mismatch: {row['sample_id']}")

    for filename in ("iedb_prediction_results.csv", "mixmhc2pred_prediction_results.csv"):
        predictions = read_csv(package / filename)
        if len(predictions) != 14 or {row["sample_id"] for row in predictions} != set(expected_by_id):
            raise ValueError(f"prediction coverage mismatch: {filename}")
        for row in predictions:
            if row["sequence"] != expected_by_id[row["sample_id"]]["sequence"]:
                raise ValueError(f"prediction sequence mismatch: {filename}/{row['sample_id']}")

    protocol = json.loads((package / "protocol_lock.json").read_text(encoding="utf-8"))
    if protocol["hla"]["beta"] != "HLA-DRB1*15:01" or protocol["apbs_candidate_ranking"] != "prohibited_prior_control_gate_failed":
        raise ValueError("protocol identity or APBS lock mismatch")

    jobs = json.loads((package / "modeling/af3_primary_mutation_panel_10_jobs.json").read_text(encoding="utf-8"))
    model_status = read_csv(package / "modeling/MODEL_STATUS.csv")
    if len(jobs) != 10 or len(model_status) != 10:
        raise ValueError("AF3 job/status count mismatch")
    allowed_model_statuses = {"complete_exact_qc", "complete_structurally_confounded"}
    if not {row["status"] for row in model_status} <= allowed_model_statuses:
        raise ValueError("model status contains an unsupported state")
    for job in jobs:
        if job["modelSeeds"] != [314159] or len(job["sequences"]) != 3:
            raise ValueError(f"invalid AF3 job: {job['name']}")

    af3_samples = read_csv(package / "modeling/analysis/af3_sample_metrics_50.csv")
    af3_jobs = read_csv(package / "modeling/analysis/af3_job_summary_10.csv")
    comparisons = read_csv(package / "modeling/analysis/wt_mutant_register_comparisons_8.csv")
    if len(af3_samples) != 50 or len(af3_jobs) != 10 or len(comparisons) != 8:
        raise ValueError("AF3 analysis row-count mismatch")
    if any(row["sequence_layout_status"] != "pass_exact_three_chain_sequence_match" for row in af3_samples):
        raise ValueError("AF3 exact chain/sequence QC failure")
    if any(row["has_clash"] != "False" for row in af3_samples):
        raise ValueError("AF3 clash-flagged model present")
    register_counts = Counter(row["register_model_status"] for row in comparisons)
    if register_counts != Counter({"supports_expected_register": 7, "structurally_confounded": 1}):
        raise ValueError("AF3 register-status mismatch")
    provenance = json.loads((package / "modeling/analysis/archive_provenance.json").read_text(encoding="utf-8"))
    archive = Path(provenance["source_archive"])
    if not archive.is_file() or hashlib.sha256(archive.read_bytes()).hexdigest() != provenance["source_archive_sha256"]:
        raise ValueError("AF3 source archive provenance mismatch")

    software = read_csv(package / "software_manifest.csv")
    for row in software:
        if not row["path"].startswith("src/"):
            continue
        path = ROOT / row["path"]
        if (
            not path.is_file()
            or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]
            or path.stat().st_size != int(row["bytes"])
        ):
            raise ValueError(f"software provenance mismatch: {row['path']}")

    plate_map = read_csv(package / "lab_handoff/randomized_plate_map_9x96.csv")
    if len(plate_map) != 864:
        raise ValueError("plate map must contain 864 wells")
    counts = Counter(row["plate_id"] for row in plate_map)
    if len(counts) != 9 or set(counts.values()) != {96}:
        raise ValueError("plate map must contain nine complete 96-well plates")
    for plate in counts:
        wells = [row["well"] for row in plate_map if row["plate_id"] == plate]
        if len(set(wells)) != 96:
            raise ValueError(f"duplicate well on {plate}")

    fitted = read_csv(package / "lab_handoff/fitted_binding_results_template.csv")
    if len(fitted) != 30 or any(row["ic50_nM"] or row["result_status"] != "not_run" for row in fitted):
        raise ValueError("fitted-results template must be blank and marked not_run")

    checksums = read_csv(package / "SHA256SUMS.csv")
    for row in checksums:
        path = package / row["relative_path"]
        if not path.is_file():
            raise ValueError(f"missing checksummed file: {row['relative_path']}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"checksum mismatch: {row['relative_path']}")

    print(
        "QA passed: 14 exact peptides; 28 predictor records; 10 AF3 jobs/50 exact models; "
        "7/8 structurally register-supported mutants; 9 complete plates; "
        "30 blank fitted-result rows; all checksums matched"
    )


if __name__ == "__main__":
    main()
