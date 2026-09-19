#!/usr/bin/env python3
"""Build AlphaFold Server upload batches for the verified missing V2 jobs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parent
PROJECT_ROOT = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/"
    "ebv_ms_publication"
)
V2_DIR = PROJECT_ROOT / "processed/tcell_library_v2_2026-08-22"
MODEL_MANIFEST = V2_DIR / "model_inventory_320.csv"
REFERENCE_MANIFEST = V2_DIR / "reference_sequence_manifest.csv"
REMAINING_72 = WORKSPACE / "outputs/alphafold_v2_remaining_snapshot_2026-08-23.csv"
UNLISTED_9 = (
    WORKSPACE
    / "outputs/alphafold_v2_download_audit_2026-08-23/unlisted_missing_9.csv"
)
DEFAULT_OUTPUT_DIR = WORKSPACE / "outputs/alphafold_v2_missing_81_jobs_2026-08-23"

sys.path.insert(0, str(PROJECT_ROOT / "src"))
from build_tcell_library_v2 import DRA_SEQUENCE  # noqa: E402


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_jobs() -> tuple[list[dict], list[dict[str, str]]]:
    inventory = read_csv(MODEL_MANIFEST)
    if len(inventory) != 320 or len({row["job_name"] for row in inventory}) != 320:
        raise ValueError("frozen model manifest must contain 320 unique jobs")

    missing_names = {
        row["job_name"] for row in read_csv(REMAINING_72) + read_csv(UNLISTED_9)
    }
    if len(missing_names) != 81:
        raise ValueError("verified remainder must contain 81 unique jobs")

    drb_sequences = {
        row["entity"]: row["sequence"]
        for row in read_csv(REFERENCE_MANIFEST)
        if row["entity"].startswith("HLA-DRB1*")
    }
    expected_alleles = {
        "HLA-DRB1*15:01",
        "HLA-DRB1*13:03",
        "HLA-DRB1*03:01",
        "HLA-DRB1*08:01",
    }
    if set(drb_sequences) != expected_alleles:
        raise ValueError("reference manifest does not contain the four frozen DRB1 chains")

    selected_rows = [row for row in inventory if row["job_name"] in missing_names]
    if len(selected_rows) != 81:
        found = {row["job_name"] for row in selected_rows}
        raise ValueError(f"missing names absent from frozen manifest: {sorted(missing_names - found)}")

    jobs = []
    for row in selected_rows:
        chains = (DRA_SEQUENCE, drb_sequences[row["allele"]], row["peptide_sequence"])
        jobs.append({
            "name": row["job_name"],
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": sequence, "count": 1}}
                for sequence in chains
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        })
    return jobs, selected_rows


def write_package(output_dir: Path) -> None:
    jobs, rows = build_jobs()
    batches = [jobs[index:index + 30] for index in range(0, len(jobs), 30)]
    if [len(batch) for batch in batches] != [30, 30, 21]:
        raise AssertionError("expected a 30, 30, 21 partition")

    output_dir.mkdir(parents=True, exist_ok=True)
    json_paths = []
    for index, batch in enumerate(batches, start=1):
        path = output_dir / f"ebvms_v2_missing_batch_{index:02d}_{len(batch)}_jobs.json"
        path.write_text(json.dumps(batch, indent=2) + "\n", encoding="utf-8")
        json_paths.append(path)

    manifest_path = output_dir / "missing_81_job_manifest.csv"
    fields = [
        "batch", "batch_job_number", "job_name", "allele", "candidate_id",
        "peptide_sequence", "peptide_length",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        offset = 0
        for batch_index, batch in enumerate(batches, start=1):
            for batch_job_number, job in enumerate(batch, start=1):
                row = rows[offset]
                writer.writerow({
                    "batch": f"batch_{batch_index:02d}",
                    "batch_job_number": batch_job_number,
                    "job_name": job["name"],
                    "allele": row["allele"],
                    "candidate_id": row["candidate_id"],
                    "peptide_sequence": row["peptide_sequence"],
                    "peptide_length": len(row["peptide_sequence"]),
                })
                offset += 1

    checksum_paths = json_paths + [manifest_path]
    checksum_text = "".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
        for path in checksum_paths
    )
    (output_dir / "SHA256SUMS.txt").write_text(checksum_text, encoding="utf-8")
    (output_dir / "README.md").write_text(
        "# EBV-MS V2 missing-job retry package\n\n"
        "This package contains the 81 discovery jobs absent from the verified local "
        "Downloads inventory, split into AlphaFold Server uploads of 30, 30, and 21 jobs.\n\n"
        "Jobs retain frozen manifest order, exact mature HLA-DRA/DRB1 chains, full peptide "
        "sequences, `modelSeeds: []`, `dialect: alphafoldserver`, and `version: 1`.\n\n"
        "Prepared upload files are not evidence of submission or completed prediction.\n",
        encoding="utf-8",
    )
    print(f"Built 81 unique jobs in batches: 30, 30, 21")
    print(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    write_package(args.output_dir.resolve())


if __name__ == "__main__":
    main()
