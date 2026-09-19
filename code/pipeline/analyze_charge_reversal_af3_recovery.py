"""Validate and analyze the 30-job controlled AlphaFold charge-reversal return.

The outputs are structural sensitivity/QC evidence, not affinity measurements.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from analyze_charge_reversal_af3 import (
    CLAIM_BOUNDARY,
    local_contacts,
    model_index,
    normalized_job,
    peptide_slice_rmsd,
)
from analyze_af3_pmhc_downloads import parse_mmcif, peptide_hla_metrics, sequence
from charge_reversal_recovery import classify_consensus, primary_ten


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_recovery_jobs(package: Path, source: Path) -> tuple[dict[str, Any], list[dict[str, object]]]:
    submitted_path = package / "alphafold_ablation" / "af3_recovery_batch_30_jobs.json"
    submitted_raw = json.loads(submitted_path.read_text(encoding="utf-8"))
    submitted = {job["name"]: normalized_job(job) for job in submitted_raw}
    status = {row["job_name"]: row for row in read_csv(package / "alphafold_ablation" / "AF3_RECOVERY_STATUS_30.csv")}
    for row in status.values():
        state = "on" if row["peptide_templates_enabled"] == "True" else "off"
        row["condition"] = f"peptide_templates_{state}_seed_{row['seed']}"
    peptides = {row["sample_id"]: row for row in primary_ten()}
    directories = sorted(path for path in source.iterdir() if path.is_dir())
    if len(submitted) != 30 or set(submitted) != set(status) or {path.name for path in directories} != set(submitted):
        raise ValueError("returned directories, submitted batch, and 30-job status table do not match")

    jobs: dict[str, Any] = {}
    model_rows: list[dict[str, object]] = []
    for directory in directories:
        name = directory.name
        request_paths = list(directory.glob("*_job_request.json"))
        if len(request_paths) != 1:
            raise ValueError(f"expected exactly one returned request for {name}")
        returned_payload = json.loads(request_paths[0].read_text(encoding="utf-8"))
        if not isinstance(returned_payload, list) or len(returned_payload) != 1:
            raise ValueError(f"unexpected returned request shape for {name}")
        returned_job = normalized_job(returned_payload[0])
        submitted_job = submitted[name]
        returned_seed = returned_job["modelSeeds"][0]
        requested_seed = submitted_job["modelSeeds"][0]
        seed_match = returned_seed == requested_seed
        returned_except_seed = {key: value for key, value in returned_job.items() if key != "modelSeeds"}
        submitted_except_seed = {key: value for key, value in submitted_job.items() if key != "modelSeeds"}
        if returned_except_seed != submitted_except_seed:
            raise ValueError(f"returned request differs from submitted sequence/template design: {name}")
        status[name]["returned_seed"] = returned_seed
        status[name]["seed_match"] = str(seed_match)

        meta = peptides[status[name]["sample_id"]]
        requested = [entry["sequence"] for entry in submitted[name]["sequences"]]
        models: dict[int, dict[str, Any]] = {}
        for cif_path in sorted(directory.glob("*_model_*.cif"), key=model_index):
            index = model_index(cif_path)
            summary_paths = list(directory.glob(f"*_summary_confidences_{index}.json"))
            full_paths = list(directory.glob(f"*_full_data_{index}.json"))
            if len(summary_paths) != 1 or len(full_paths) != 1:
                raise ValueError(f"missing or duplicate confidence/full-data partner for {name} model {index}")
            structure = parse_mmcif(cif_path)
            observed = [sequence(structure.get(chain, [])) for chain in ("A", "B", "C")]
            if set(structure) != {"A", "B", "C"} or observed != requested:
                raise ValueError(f"exact three-chain sequence mismatch: {name} model {index}")
            summary = json.loads(summary_paths[0].read_text(encoding="utf-8"))
            metrics = peptide_hla_metrics(structure["C"], [structure["A"], structure["B"]])
            pair_iptm = [float(summary["chain_pair_iptm"][2][i]) for i in (0, 1)]
            pair_pae = [float(summary["chain_pair_pae_min"][2][i]) for i in (0, 1)]
            contacts = local_contacts(structure, int(meta["peptide_position_1_based"])) if meta["peptide_position_1_based"] else {}
            row = {
                "job_name": name,
                "sample_id": meta["sample_id"],
                "condition": status[name]["condition"],
                "peptide_templates_enabled": status[name]["peptide_templates_enabled"],
                "requested_seed": status[name]["seed"],
                "returned_seed": returned_seed,
                "seed_match": seed_match,
                "model_index": index,
                "sequence_layout_status": "pass_exact_three_chain_sequence_match",
                "peptide": requested[2],
                "ranking_score": summary["ranking_score"],
                "iptm": summary["iptm"],
                "ptm": summary["ptm"],
                "peptide_chain_iptm": summary["chain_iptm"][2],
                "peptide_hla_pair_iptm_mean": round(sum(pair_iptm) / 2, 3),
                "peptide_hla_pair_pae_min_mean_A": round(sum(pair_pae) / 2, 3),
                "has_clash": bool(summary["has_clash"]),
                **metrics,
                **contacts,
                "cif_relative_path_from_return": str(cif_path.relative_to(source)),
                "claim_boundary": CLAIM_BOUNDARY,
            }
            model_rows.append(row)
            models[index] = {"structure": structure, "summary": summary, "row": row}
        if set(models) != set(range(5)):
            raise ValueError(f"expected model indices 0-4 for {name}; found {sorted(models)}")
        jobs[name] = {"models": models, "meta": meta, "status": status[name], "requested": requested}
    if len(model_rows) != 150:
        raise ValueError(f"expected 150 models, found {len(model_rows)}")
    return jobs, model_rows


def condition_comparisons(jobs: dict[str, Any]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], str] = {}
    for name, job in jobs.items():
        grouped[(str(job["meta"]["sample_id"]), str(job["status"]["condition"]))] = name
    conditions = sorted({str(job["status"]["condition"]) for job in jobs.values()})
    rows: list[dict[str, object]] = []
    for condition in conditions:
        for arm in ("BALF5", "TALDO1"):
            wt = jobs[grouped[(f"{arm}_WT_15", condition)]]
            expected_start = str(wt["meta"]["sequence"]).index(str(wt["meta"]["expected_core"]))
            for suffix in ("P6_KQ", "P6_KE", "P7_KQ", "P7_KE"):
                mutant = jobs[grouped[(f"{arm}_{suffix}", condition)]]
                windows: dict[int, list[float]] = {start: [] for start in range(7)}
                full: list[float] = []
                for wt_model in wt["models"].values():
                    for mutant_model in mutant["models"].values():
                        full.append(peptide_slice_rmsd(wt_model["structure"], mutant_model["structure"], slice(None), slice(None)))
                        for start in range(7):
                            windows[start].append(peptide_slice_rmsd(
                                wt_model["structure"], mutant_model["structure"],
                                slice(expected_start, expected_start + 9), slice(start, start + 9),
                            ))
                window_medians = {start: median(values) for start, values in windows.items()}
                best_start = min(window_medians, key=window_medians.get)
                expected_rmsd = window_medians[expected_start]
                stable = best_start == expected_start and expected_rmsd <= 2.0
                contact_rows = [model["row"] for model in mutant["models"].values()]
                rows.append({
                    "condition": condition,
                    "peptide_templates_enabled": mutant["status"]["peptide_templates_enabled"],
                    "requested_seed": mutant["status"]["seed"],
                    "returned_seed": mutant["status"]["returned_seed"],
                    "seed_match": mutant["status"]["seed_match"],
                    "arm": arm,
                    "mutant_sample_id": mutant["meta"]["sample_id"],
                    "substitution": mutant["meta"]["substitution"],
                    "core_position": mutant["meta"]["core_position"],
                    "peptide_position_1_based": mutant["meta"]["peptide_position_1_based"],
                    "wt_mutant_model_pairs": 25,
                    "full_peptide_rmsd_median_A": round(median(full), 3),
                    "expected_core_start_1_based": expected_start + 1,
                    "best_matching_mutant_window_start_1_based": best_start + 1,
                    "expected_core_window_rmsd_median_A": round(expected_rmsd, 3),
                    "best_window_rmsd_median_A": round(window_medians[best_start], 3),
                    "register_model_status": "supports_expected_register" if stable else "structurally_confounded",
                    "sidechain_hla_contact_frequency": round(sum(int(row.get("sidechain_hla_partner_count", 0)) > 0 for row in contact_rows) / 5, 3),
                    "opposite_charge_contact_frequency": round(sum(bool(row.get("opposite_charge_atom_pair_within_4A", False)) for row in contact_rows) / 5, 3),
                    "any_clash": any(bool(row["has_clash"]) for row in contact_rows),
                    "claim_boundary": CLAIM_BOUNDARY,
                })
    if len(rows) != 24:
        raise AssertionError(f"expected 24 recovery comparisons, found {len(rows)}")
    return rows


def job_summaries(jobs: dict[str, Any]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for name, job in sorted(jobs.items()):
        models = list(job["models"].values())
        model_rows = [model["row"] for model in models]
        pairwise = [
            peptide_slice_rmsd(models[left]["structure"], models[right]["structure"], slice(None), slice(None))
            for left in range(5) for right in range(left + 1, 5)
        ]
        selected = max(model_rows, key=lambda row: float(row["ranking_score"]))
        rows.append({
            "job_name": name,
            "sample_id": job["meta"]["sample_id"],
            "condition": job["status"]["condition"],
            "peptide_templates_enabled": job["status"]["peptide_templates_enabled"],
            "requested_seed": job["status"]["seed"],
            "returned_seed": job["status"]["returned_seed"],
            "seed_match": job["status"]["seed_match"],
            "models": 5,
            "all_exact_sequence_layout_checks_pass": True,
            "any_clash": any(bool(row["has_clash"]) for row in model_rows),
            "selected_model_index": selected["model_index"],
            "selected_ranking_score": selected["ranking_score"],
            "median_peptide_chain_iptm": round(median(float(row["peptide_chain_iptm"]) for row in model_rows), 3),
            "median_peptide_mean_plddt": round(median(float(row["peptide_mean_plddt"]) for row in model_rows), 2),
            "within_job_pairwise_peptide_rmsd_median_A": round(median(pairwise), 3),
            "within_job_pairwise_peptide_rmsd_max_A": round(max(pairwise), 3),
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return rows


def mutation_consensus(package: Path, jobs: dict[str, Any], comparisons: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_mutation: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in comparisons:
        by_mutation[str(row["mutant_sample_id"])][str(row["condition"])] = row
    pandora_rows = {row["sample_id"]: row for row in read_csv(package / "pandora" / "results" / "pandora_candidate_sample_summary_10.csv")}
    calibration = json.loads((package / "pandora" / "results" / "CALIBRATION_GATE.json").read_text(encoding="utf-8"))
    status_by_sample_condition = {
        (str(job["meta"]["sample_id"]), str(job["status"]["condition"])): job for job in jobs.values()
    }
    detail_rows: list[dict[str, object]] = []
    consensus_rows: list[dict[str, object]] = []
    off_conditions = ("peptide_templates_off_seed_314159", "peptide_templates_off_seed_104729")
    for peptide in primary_ten():
        if peptide["reagent_role"] == "wild_type":
            continue
        sample_id = str(peptide["sample_id"])
        condition_rows = [by_mutation[sample_id][condition] for condition in off_conditions]
        starts = [int(row["best_matching_mutant_window_start_1_based"]) for row in condition_rows]
        af3_register_stable = all(row["register_model_status"] == "supports_expected_register" for row in condition_rows)
        af3_seeds_agree = len(set(starts)) == 1
        af3_contact = min(float(row["opposite_charge_contact_frequency"]) for row in condition_rows)
        af3_any_clash = any(bool(row["any_clash"]) for row in condition_rows)
        relevant_jobs = [status_by_sample_condition[(sample_id, condition)] for condition in off_conditions]
        af3_seed_qc_pass = all(job["status"]["seed_match"] == "True" for job in relevant_jobs)

        job_a = status_by_sample_condition[(sample_id, off_conditions[0])]
        job_b = status_by_sample_condition[(sample_id, off_conditions[1])]
        selected_a = max(job_a["models"].values(), key=lambda model: float(model["row"]["ranking_score"]))
        selected_b = max(job_b["models"].values(), key=lambda model: float(model["row"]["ranking_score"]))
        expected_start = str(peptide["sequence"]).index(str(peptide["expected_core"]))
        af3_cross_seed_rmsd = peptide_slice_rmsd(
            selected_a["structure"], selected_b["structure"],
            slice(expected_start, expected_start + 9), slice(expected_start, expected_start + 9),
        )

        pandora = pandora_rows[sample_id]
        pandora_expected = pandora["expected_register_preferred_both_templates"] == "True"
        pandora_templates_agree = pandora["preferred_register_1BX2"] == pandora["preferred_register_6CQQ"]
        pandora_rmsd = float(pandora["expected_register_top_model_cross_template_core_rmsd_A"])
        p_contacts = [
            float(pandora["expected_register_opposite_charge_contact_frequency_1BX2"] or 0),
            float(pandora["expected_register_opposite_charge_contact_frequency_6CQQ"] or 0),
        ]
        pandora_contact = min(p_contacts)
        evidence = {
            "complete": af3_seed_qc_pass,
            "calibration_pass": calibration["calibration_status"] == "pass",
            "sequence_qc_pass": True,
            "structural_qc_pass": not af3_any_clash,
            "expected_register_preferred": af3_register_stable and pandora_expected,
            "expected_register_core_rmsd_max_A": max(af3_cross_seed_rmsd, pandora_rmsd),
            "af3_seeds_agree": af3_seeds_agree,
            "pandora_templates_agree": pandora_templates_agree,
            "methods_agree": af3_register_stable == pandora_expected,
            "af3_contact_frequency": af3_contact,
            "pandora_contact_frequency": pandora_contact,
        }
        classification = classify_consensus(evidence)
        if not af3_seed_qc_pass:
            reason = "One required template-free AlphaFold job returned a different seed than the predeclared seed."
        elif classification == "register_confounded":
            reason = "The expected register did not survive both template-free AlphaFold seeds and both PANDORA templates."
        elif classification == "method_or_template_dependent":
            reason = "Register or contact evidence disagreed across seeds, templates, or methods."
        elif classification == "robust_structural_contact_hypothesis":
            reason = "The predeclared cross-method register, RMSD, and 80% contact criteria passed."
        elif classification == "no_structural_charge_contact_support":
            reason = "The register was stable, but the predeclared 80% cross-method charge-contact criterion failed."
        else:
            reason = "Completeness, calibration, sequence, or structural QC failed."
        detail_rows.append({
            "mutation_sample_id": sample_id,
            "af3_template_free_seed_314159_register_status": condition_rows[0]["register_model_status"],
            "af3_template_free_seed_104729_register_status": condition_rows[1]["register_model_status"],
            "af3_template_free_preferred_starts_1_based": ";".join(map(str, starts)),
            "af3_seeds_agree": af3_seeds_agree,
            "af3_predeclared_seed_qc_pass": af3_seed_qc_pass,
            "af3_selected_core_cross_seed_rmsd_A": round(af3_cross_seed_rmsd, 3),
            "af3_min_opposite_charge_contact_frequency": af3_contact,
            "pandora_expected_register_preferred_both_templates": pandora_expected,
            "pandora_templates_agree": pandora_templates_agree,
            "pandora_expected_core_cross_template_rmsd_A": pandora_rmsd,
            "pandora_min_opposite_charge_contact_frequency": pandora_contact,
            "classification": classification,
            "reason": reason,
            "claim_boundary": CLAIM_BOUNDARY,
        })
        consensus_rows.append({
            "mutation_sample_id": sample_id,
            "classification": classification,
            "reason": reason,
            "af3_models_observed": 15,
            "pandora_models_observed": 120,
            "models_are_biological_replicates": False,
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return detail_rows, consensus_rows


def analyze(package: Path, source: Path, archive: Path) -> None:
    output = package / "alphafold_ablation" / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    jobs, model_rows = load_recovery_jobs(package, source)
    summaries = job_summaries(jobs)
    comparisons = condition_comparisons(jobs)
    details, consensus = mutation_consensus(package, jobs, comparisons)
    write_csv(output / "af3_recovery_model_metrics_150.csv", model_rows)
    write_csv(output / "af3_recovery_job_summary_30.csv", summaries)
    write_csv(output / "af3_recovery_register_comparisons_24.csv", comparisons)
    write_csv(output / "af3_pandora_mutation_evidence_8.csv", details)
    write_csv(package / "consensus" / "consensus_classification_8.csv", consensus)

    status_path = package / "alphafold_ablation" / "AF3_RECOVERY_STATUS_30.csv"
    status_rows = read_csv(status_path)
    job_status = {name: job["status"] for name, job in jobs.items()}
    for row in status_rows:
        returned_seed = job_status[row["job_name"]]["returned_seed"]
        seed_match = job_status[row["job_name"]]["seed_match"]
        row["returned_seed"] = returned_seed
        row["seed_match"] = seed_match
        row["status"] = "returned_complete_qc_pass" if seed_match == "True" else "returned_complete_seed_mismatch"
        row["returned_archive_path"] = "alphafold_ablation/raw/folds_2026_09_08_02_30.zip"
    write_csv(status_path, status_rows)
    mismatched_names = [row["job_name"] for row in status_rows if row["seed_match"] != "True"]
    submitted_jobs = json.loads((package / "alphafold_ablation" / "af3_recovery_batch_30_jobs.json").read_text(encoding="utf-8"))
    retry_jobs = [job for job in submitted_jobs if job["name"] in mismatched_names]
    if retry_jobs:
        (package / "alphafold_ablation" / "af3_retry_seed_mismatch_jobs.json").write_text(
            json.dumps(retry_jobs, indent=2) + "\n", encoding="utf-8"
        )

    raw_dir = package / "alphafold_ablation" / "raw"
    raw_dir.mkdir(exist_ok=True)
    preserved_archive = raw_dir / "folds_2026_09_08_02_30.zip"
    shutil.copy2(archive, preserved_archive)
    provenance = {
        "source_archive": str(archive),
        "source_archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "preserved_archive_relative_path": str(preserved_archive.relative_to(package)),
        "jobs": 30,
        "models": 150,
        "confidence_files": 150,
        "full_data_files": 150,
        "exact_three_chain_sequence_and_template_checks_passed": True,
        "predeclared_seed_matches": sum(job["status"]["seed_match"] == "True" for job in jobs.values()),
        "predeclared_seed_mismatches": sum(job["status"]["seed_match"] != "True" for job in jobs.values()),
        "seed_mismatch_job_names": mismatched_names,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (output / "archive_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    classes: dict[str, int] = defaultdict(int)
    for row in consensus:
        classes[str(row["classification"])] += 1
    summary = {
        "jobs": 30,
        "models": 150,
        "all_exact_sequence_qc_pass": True,
        "predeclared_seed_matches": sum(job["status"]["seed_match"] == "True" for job in jobs.values()),
        "predeclared_seed_mismatches": sum(job["status"]["seed_match"] != "True" for job in jobs.values()),
        "seed_mismatch_job_names": mismatched_names,
        "acceptance_status": "pass" if all(job["status"]["seed_match"] == "True" for job in jobs.values()) else "fail_predeclared_seed_mismatch",
        "models_with_clashes": sum(bool(row["has_clash"]) for row in model_rows),
        "consensus_class_counts": dict(classes),
        "models_are_biological_replicates": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (output / "AF3_RECOVERY_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.package.resolve(), args.source.resolve(), args.archive.resolve())


if __name__ == "__main__":
    main()
