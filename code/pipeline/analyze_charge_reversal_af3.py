"""Validate and descriptively analyze the charge-reversal AlphaFold Server return.

The outputs are structural quality-control and model-sensitivity evidence only.
They do not measure peptide--HLA affinity or establish an electrostatic contact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np

from analyze_af3_pmhc_downloads import (
    ca_coordinates,
    kabsch,
    parse_mmcif,
    peptide_hla_metrics,
    sequence,
)


CLAIM_BOUNDARY = (
    "Descriptive AlphaFold pMHC model QC only; not measured binding, direct "
    "electrostatic-contact evidence, presentation, TCR recognition, molecular mimicry, or MS mechanism."
)
BACKBONE = {"N", "CA", "C", "O", "OXT"}
ACIDIC_OXYGENS = {"D": {"OD1", "OD2"}, "E": {"OE1", "OE2"}}
BASIC_NITROGENS = {"K": {"NZ"}, "R": {"NE", "NH1", "NH2"}}


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


def model_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.cif$", path.name)
    if not match:
        raise ValueError(f"cannot parse model index: {path.name}")
    return int(match.group(1))


def normalized_job(job: dict[str, Any]) -> dict[str, Any]:
    chains = [entry["proteinChain"] for entry in job["sequences"]]
    return {
        "name": job["name"],
        "modelSeeds": [str(value) for value in job["modelSeeds"]],
        "sequences": [
            {
                "sequence": chain["sequence"],
                "count": int(chain["count"]),
                "useStructureTemplate": bool(chain["useStructureTemplate"]),
            }
            for chain in chains
        ],
        "dialect": job["dialect"],
        "version": int(job["version"]),
    }


def transform_for_hla_fit(reference: dict[str, Any], moving: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    reference_hla = np.vstack([ca_coordinates(reference[chain][:85]) for chain in ("A", "B")])
    moving_hla = np.vstack([ca_coordinates(moving[chain][:85]) for chain in ("A", "B")])
    return kabsch(moving_hla, reference_hla)


def peptide_slice_rmsd(
    reference: dict[str, Any],
    moving: dict[str, Any],
    reference_slice: slice,
    moving_slice: slice,
) -> float:
    rotation, translation = transform_for_hla_fit(reference, moving)
    target = ca_coordinates(reference["C"])[reference_slice]
    fitted = ca_coordinates(moving["C"])[moving_slice] @ rotation + translation
    if len(target) != len(fitted):
        raise ValueError("peptide comparison lengths differ")
    return float(np.sqrt(np.mean(np.sum((target - fitted) ** 2, axis=1))))


def local_contacts(structure: dict[str, Any], peptide_position: int) -> dict[str, object]:
    residue = structure["C"][peptide_position - 1]
    sidechain = [atom for atom in residue["atoms"] if atom["element"] != "H" and atom["name"] not in BACKBONE]
    partners: list[tuple[float, str]] = []
    for chain in ("A", "B"):
        for index, hla_residue in enumerate(structure[chain], 1):
            distances = [
                float(np.linalg.norm(np.asarray(atom["xyz"]) - np.asarray(other["xyz"])))
                for atom in sidechain
                for other in hla_residue["atoms"]
                if other["element"] != "H"
            ]
            if distances and min(distances) <= 4.0:
                partners.append((min(distances), f"{chain}:{hla_residue['aa']}{index}"))

    charged_distance = ""
    aa = str(residue["aa"])
    charge_atoms = ACIDIC_OXYGENS.get(aa, BASIC_NITROGENS.get(aa, set()))
    opposite = BASIC_NITROGENS if aa in ACIDIC_OXYGENS else ACIDIC_OXYGENS if aa in BASIC_NITROGENS else {}
    distances = []
    for atom in residue["atoms"]:
        if atom["name"] not in charge_atoms:
            continue
        for chain in ("A", "B"):
            for hla_residue in structure[chain]:
                allowed = opposite.get(str(hla_residue["aa"]), set())
                for other in hla_residue["atoms"]:
                    if other["name"] in allowed:
                        distances.append(float(np.linalg.norm(np.asarray(atom["xyz"]) - np.asarray(other["xyz"]))))
    if distances:
        charged_distance = round(min(distances), 3)
    return {
        "mutated_residue_observed": aa,
        "sidechain_hla_partners_within_4A": ";".join(label for _, label in sorted(partners)),
        "sidechain_hla_partner_count": len(partners),
        "nearest_opposite_charge_atom_distance_A": charged_distance,
        "opposite_charge_atom_pair_within_4A": bool(distances and min(distances) <= 4.0),
    }


def load_jobs(source: Path, submitted_path: Path, manifest_path: Path) -> tuple[dict[str, Any], list[dict[str, object]]]:
    submitted = {job["name"]: normalized_job(job) for job in json.loads(submitted_path.read_text(encoding="utf-8"))}
    manifest = {f"cr_hy15_drb1_1501_{row['sample_id'].lower()}": row for row in read_csv(manifest_path) if row["primary_mutation_panel"] == "True"}
    if set(submitted) != set(manifest) or len(submitted) != 10:
        raise ValueError("submitted jobs and primary manifest do not match exactly")

    jobs: dict[str, Any] = {}
    sample_rows: list[dict[str, object]] = []
    directories = sorted(path for path in source.iterdir() if path.is_dir())
    if {path.name for path in directories} != set(submitted):
        raise ValueError("returned job directory names do not match the ten submitted jobs")
    for directory in directories:
        request_paths = list(directory.glob("*_job_request.json"))
        if len(request_paths) != 1:
            raise ValueError(f"expected one request file in {directory}")
        returned_job = json.loads(request_paths[0].read_text(encoding="utf-8"))[0]
        if normalized_job(returned_job) != submitted[directory.name]:
            raise ValueError(f"returned request differs from submitted request: {directory.name}")
        requested = [entry["sequence"] for entry in submitted[directory.name]["sequences"]]
        models: dict[int, dict[str, Any]] = {}
        for cif_path in sorted(directory.glob("*_model_*.cif"), key=model_index):
            index = model_index(cif_path)
            summary_matches = list(directory.glob(f"*_summary_confidences_{index}.json"))
            full_matches = list(directory.glob(f"*_full_data_{index}.json"))
            if len(summary_matches) != 1 or len(full_matches) != 1:
                raise ValueError(f"missing or duplicate JSON partner for {directory.name} model {index}")
            structure = parse_mmcif(cif_path)
            observed = [sequence(structure.get(chain, [])) for chain in ("A", "B", "C")]
            if set(structure) != {"A", "B", "C"} or observed != requested:
                raise ValueError(f"exact three-chain sequence mismatch: {directory.name} model {index}")
            summary = json.loads(summary_matches[0].read_text(encoding="utf-8"))
            metrics = peptide_hla_metrics(structure["C"], [structure["A"], structure["B"]])
            pair_iptm = [float(summary["chain_pair_iptm"][2][i]) for i in (0, 1)]
            pair_pae = [float(summary["chain_pair_pae_min"][2][i]) for i in (0, 1)]
            meta = manifest[directory.name]
            contacts = local_contacts(structure, int(meta["peptide_position_1_based"])) if meta["peptide_position_1_based"] else {}
            row = {
                "job_name": directory.name,
                "sample_id": meta["sample_id"],
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
                "cif_relative_path": str(cif_path.relative_to(source.parent.parent)),
                "claim_boundary": CLAIM_BOUNDARY,
            }
            sample_rows.append(row)
            models[index] = {"structure": structure, "summary": summary, "row": row}
        if set(models) != set(range(5)):
            raise ValueError(f"expected model indices 0-4 for {directory.name}; found {sorted(models)}")
        jobs[directory.name] = {"models": models, "meta": manifest[directory.name], "requested": requested}
    return jobs, sample_rows


def analyze(package: Path, source: Path, archive: Path) -> None:
    output = package / "modeling" / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    jobs, sample_rows = load_jobs(
        source,
        package / "modeling" / "af3_primary_mutation_panel_10_jobs.json",
        package / "peptide_manifest.csv",
    )

    job_rows: list[dict[str, object]] = []
    for name, job in sorted(jobs.items()):
        rows = [model["row"] for model in job["models"].values()]
        pairwise = [
            peptide_slice_rmsd(job["models"][left]["structure"], job["models"][right]["structure"], slice(None), slice(None))
            for left in range(5) for right in range(left + 1, 5)
        ]
        selected = max(rows, key=lambda row: float(row["ranking_score"]))
        job_rows.append({
            "job_name": name,
            "sample_id": job["meta"]["sample_id"],
            "peptide": job["requested"][2],
            "models": 5,
            "all_exact_sequence_layout_checks_pass": True,
            "any_clash": any(bool(row["has_clash"]) for row in rows),
            "selected_model_index": selected["model_index"],
            "selected_ranking_score": selected["ranking_score"],
            "median_peptide_chain_iptm": round(median(float(row["peptide_chain_iptm"]) for row in rows), 3),
            "median_peptide_hla_pair_iptm_mean": round(median(float(row["peptide_hla_pair_iptm_mean"]) for row in rows), 3),
            "median_peptide_mean_plddt": round(median(float(row["peptide_mean_plddt"]) for row in rows), 2),
            "within_job_pairwise_peptide_rmsd_median_A": round(median(pairwise), 3),
            "within_job_pairwise_peptide_rmsd_max_A": round(max(pairwise), 3),
            "claim_boundary": CLAIM_BOUNDARY,
        })

    comparison_rows: list[dict[str, object]] = []
    for arm in ("balf5", "taldo1"):
        wt_name = f"cr_hy15_drb1_1501_{arm}_wt_15"
        wt = jobs[wt_name]
        expected_start = str(wt["meta"]["sequence"]).index(str(wt["meta"]["expected_core"]))
        for suffix in ("p6_kq", "p6_ke", "p7_kq", "p7_ke"):
            name = f"cr_hy15_drb1_1501_{arm}_{suffix}"
            mutant = jobs[name]
            window_values: dict[int, list[float]] = {start: [] for start in range(7)}
            full_values: list[float] = []
            for wt_model in wt["models"].values():
                for mutant_model in mutant["models"].values():
                    full_values.append(peptide_slice_rmsd(wt_model["structure"], mutant_model["structure"], slice(None), slice(None)))
                    for start in range(7):
                        window_values[start].append(peptide_slice_rmsd(
                            wt_model["structure"], mutant_model["structure"],
                            slice(start, start + 9), slice(expected_start, expected_start + 9),
                        ))
            window_medians = {start: median(values) for start, values in window_values.items()}
            best_start = min(window_medians, key=window_medians.get)
            expected_median = window_medians[expected_start]
            stable = best_start == expected_start and expected_median <= 2.0
            contacts = [model["row"] for model in mutant["models"].values()]
            comparison_rows.append({
                "arm": arm.upper(),
                "mutant_sample_id": mutant["meta"]["sample_id"],
                "substitution": mutant["meta"]["substitution"],
                "core_position": mutant["meta"]["core_position"],
                "peptide_position_1_based": mutant["meta"]["peptide_position_1_based"],
                "wt_mutant_model_pairs": 25,
                "full_peptide_rmsd_median_A": round(median(full_values), 3),
                "full_peptide_rmsd_min_A": round(min(full_values), 3),
                "full_peptide_rmsd_max_A": round(max(full_values), 3),
                "expected_core_start_1_based": expected_start + 1,
                "best_matching_wt_window_start_1_based": best_start + 1,
                "expected_core_window_rmsd_median_A": round(expected_median, 3),
                "best_window_rmsd_median_A": round(window_medians[best_start], 3),
                "register_model_status": "supports_expected_register" if stable else "structurally_confounded",
                "register_rule": "expected WT 9-mer window must be best of seven and median RMSD <=2.0 A",
                "models_with_sidechain_hla_partner_within_4A": sum(int(row.get("sidechain_hla_partner_count", 0)) > 0 for row in contacts),
                "models_with_opposite_charge_atom_pair_within_4A": sum(bool(row.get("opposite_charge_atom_pair_within_4A", False)) for row in contacts),
                "claim_boundary": CLAIM_BOUNDARY,
            })

    write_csv(output / "af3_sample_metrics_50.csv", sample_rows)
    write_csv(output / "af3_job_summary_10.csv", job_rows)
    write_csv(output / "wt_mutant_register_comparisons_8.csv", comparison_rows)
    provenance = {
        "source_archive": str(archive.resolve()),
        "source_archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "source_archive_bytes": archive.stat().st_size,
        "extracted_result_directory": str(source.resolve()),
        "jobs": len(jobs),
        "models": len(sample_rows),
        "all_exact_job_and_chain_sequence_checks_passed": True,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (output / "archive_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    stable_count = sum(row["register_model_status"] == "supports_expected_register" for row in comparison_rows)
    clash_count = sum(bool(row["has_clash"]) for row in sample_rows)
    (output / "README.md").write_text(
        "# Charge-reversal AlphaFold structural QC\n\n"
        f"The returned archive contains all 10 exact submitted jobs and all 50 requested models. "
        f"All 50 models passed exact DRA, DRB1*15:01, and peptide sequence/layout checks; {clash_count} models were flagged for clashes.\n\n"
        f"Using the stated ensemble-calibrated geometry rule, {stable_count}/8 mutant ensembles support the expected WT P1-P9 window. "
        "Any `structurally_confounded` row must not be interpreted as an electrostatic result. "
        "Local residue contacts are descriptive and are not proof of a salt bridge or of binding affinity.\n\n"
        "Files: `af3_sample_metrics_50.csv` (all models), `af3_job_summary_10.csv` (job summaries), "
        "`wt_mutant_register_comparisons_8.csv` (WT-mutant register checks), and `archive_provenance.json` (raw linkage).\n\n"
        f"Claim boundary: {CLAIM_BOUNDARY}\n",
        encoding="utf-8",
    )
    print(f"Analyzed {len(jobs)} jobs/{len(sample_rows)} models; register support {stable_count}/8; clashes {clash_count}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.package.resolve(), args.source.resolve(), args.archive.resolve())


if __name__ == "__main__":
    main()
