#!/usr/bin/env python3
"""Consolidate immunopeptidomics and TCR-score control evidence for the two leads."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[2]
AUTHORITATIVE = Path("/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication")
EVIDENCE = AUTHORITATIVE / "processed/high_yield_candidate_evidence_2026-08-28"
PHYS = AUTHORITATIVE / "processed/physicochemical_similarity_upgrade_2026-09-14"
BENCH = AUTHORITATIVE / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26"
EXPANSION = WORKSPACE / "output/revision_round_3/expansion_sources_2026-09-15"
OUT = WORKSPACE / "output/revision_round_3/evidence_and_control_audit_2026-09-15"
TARGETS = {"HY13_SEQ_02", "HY15_SEQ_02"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    matrix_path = EVIDENCE / "candidate_evidence_matrix.csv"
    hits_path = EVIDENCE / "immunopeptidome_hits.csv"
    pxd_status_path = EVIDENCE / "raw_responses/pxd068488_status.json"
    control_summary_path = PHYS / "control_method_summary.csv"
    grantham_gate_path = PHYS / "grantham_promotion_gate.json"
    ranking_gate_path = BENCH / "benchmark/definitive_ranking_gate.json"
    specificity_gate_path = BENCH / "benchmark/specificity_gate.json"
    system_comparisons_path = BENCH / "benchmark/system_level_comparisons.csv"
    pride_project_path = EXPANSION / "pride_PXD068488_project.json"

    matrix = [row for row in read_csv(matrix_path) if row["target_id"] in TARGETS]
    hits = [row for row in read_csv(hits_path) if row["target_id"] in TARGETS]
    summary_rows = []
    for row in matrix:
        target_hits = [hit for hit in hits if hit["target_id"] == row["target_id"]]
        summary_rows.append({
            "target_id": row["target_id"],
            "allele": row["allele"],
            "ebv_peptide": row["ebv_sequence"],
            "self_peptide": row["self_sequence"],
            "hla_ligand_atlas_overlap_rows": len(target_hits),
            "atlas_exact_hla_compatible_rows": sum(hit["exact_hla_compatible"] == "True" for hit in target_hits),
            "atlas_monoallelic_rows": sum(hit["monoallelic"] == "True" for hit in target_hits),
            "iedb_exact_hla_positive_arm_count": row["iedb_exact_hla_positive_arm_count"],
            "pxd068488_processed_table_hit_count": 0,
            "raw_spectrum_reanalysis_status": "not_performed",
            "interpretation": (
                "self-core overlap observed in multiallelic HLA-II tissue data; exact lead-allele attribution absent"
                if target_hits else "no overlap found in the audited processed sources; absence is not a negative"
            ),
            "claim_boundary": "Presentation support is source- and allele-specific; no exact lead-allele presentation or cross-reactivity claim is made.",
        })
    write_csv(OUT / "lead_immunopeptidomics_summary.csv", summary_rows)
    write_csv(OUT / "lead_hla_ligand_atlas_hits.csv", hits)

    pxd_status = json.loads(pxd_status_path.read_text())
    pride = json.loads(pride_project_path.read_text())
    raw_plan = {
        "dataset": "PXD068488",
        "live_pride_project_title": pride.get("title"),
        "live_pride_project_sha256": digest(pride_project_path),
        "submitted_processed_tables_searched": 8,
        "submitted_processed_table_hit_count": pxd_status["processed_hit_count"],
        "raw_spectrum_reanalysis_status": "prepared_not_run",
        "why_not_complete_in_this_revision": (
            "A valid raw reanalysis requires file-level acquisition, a frozen search database and contaminant set, "
            "enzyme-unspecific HLA-II search settings, FDR control, and sample-to-allele provenance. Running a partial "
            "or uncalibrated search would be less defensible than reporting the processed-table audit honestly."
        ),
        "frozen_required_steps": [
            "inventory every raw and result file with size and checksum",
            "map sample identifiers to DR2a, DR2b, donor, and biological source",
            "freeze protein database, decoys, contaminants, modifications, and no-enzyme search settings",
            "control peptide-spectrum-match and peptide-level FDR",
            "search exact peptides, declared cores, and nested ligands",
            "report detections and nondetections by evaluable sample without treating absence as a negative",
        ],
        "claim_boundary": "Prepared protocol only; no new raw-spectrum result is claimed.",
    }
    (OUT / "raw_immunopeptidomics_reanalysis_lock.json").write_text(json.dumps(raw_plan, indent=2, sort_keys=True) + "\n")

    control_summary = read_csv(control_summary_path)
    system_rows = read_csv(system_comparisons_path)
    write_csv(OUT / "physicochemical_control_method_summary.csv", control_summary)
    write_csv(OUT / "tcr_system_level_comparisons.csv", system_rows)
    control_status = {
        "grantham_gate": json.loads(grantham_gate_path.read_text()),
        "definitive_ranking_gate": json.loads(ranking_gate_path.read_text()),
        "specificity_gate": json.loads(specificity_gate_path.read_text()),
        "audited_control_panels": 8,
        "strict_independent_systems": 3,
        "interpretation": (
            "Grantham is the best-performing tested physicochemical sensitivity metric on the existing controls, "
            "but the general TCR-facing scoring claim remains provisional because the registry has only three "
            "independent systems and no verified specificity negatives."
        ),
        "claim_boundary": "Control support is provisional and cannot establish TCR recognition, specificity, or cross-reactivity for the leads.",
    }
    (OUT / "tcr_score_validation_status.json").write_text(json.dumps(control_status, indent=2, sort_keys=True) + "\n")

    sources = [
        matrix_path, hits_path, pxd_status_path, control_summary_path, grantham_gate_path,
        ranking_gate_path, specificity_gate_path, system_comparisons_path, pride_project_path,
    ]
    manifest = [{"source_path": str(path), "sha256": digest(path), "bytes": path.stat().st_size} for path in sources]
    write_csv(OUT / "source_manifest.csv", manifest)
    (OUT / "README.md").write_text(
        "# Immunopeptidomics and control-validation audit\n\n"
        "This package extracts the two lead-specific presentation evidence and the current "
        "TCR-facing metric calibration status from checksum-linked project artifacts.\n\n"
        "The strongest presentation statement supported today is that nested TALDO1 ligands "
        "overlapping each lead core occur in multiallelic HLA-II tissue datasets. The exact "
        "lead HLA restrictions are not established. PXD068488 submitter tables were searched "
        "with zero lead hits, but raw spectra were not reprocessed and nondetection is not a negative.\n\n"
        "The existing score calibration contains eight panels across three independent systems. "
        "It supports Grantham as the best tested physicochemical sensitivity metric, but the "
        "predeclared benchmark itself blocks a definitive ranking claim until at least six "
        "independent systems and verified negatives are available.\n",
        encoding="utf-8",
    )
    checksum_path = OUT / "SHA256SUMS.csv"
    checksum_rows = []
    for path in sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_path):
        checksum_rows.append({"relative_path": str(path.relative_to(OUT)), "sha256": digest(path), "bytes": path.stat().st_size})
    write_csv(checksum_path, checksum_rows)
    print(json.dumps({"lead_rows": len(summary_rows), "lead_hit_rows": len(hits), "control_systems": 3, "status": "complete_audit"}, indent=2))


if __name__ == "__main__":
    main()
