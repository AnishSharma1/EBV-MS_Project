"""Validate and summarize the 1,200-model PANDORA candidate ensemble."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from analyze_charge_reversal_pandora_calibration import alignment_pairs, kabsch, residues, sequence
from build_tcell_library_v2 import DRA_SEQUENCE, DRB1_1501_SEQUENCE
from charge_reversal_binding_pilot import primary_panel


BACKBONE = {"N", "CA", "C", "O", "OXT"}
ACIDIC_OXYGENS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
BASIC_NITROGENS = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, values: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def hla_transform(reference, moving):
    pairs = alignment_pairs(residues(reference, "M"), residues(moving, "M"), 90)
    pairs += alignment_pairs(residues(reference, "N"), residues(moving, "N"), 95)
    target = np.asarray([pair[0]["CA"].coord for pair in pairs])
    source = np.asarray([pair[1]["CA"].coord for pair in pairs])
    return kabsch(source, target)


def peptide_core_rmsd(reference, moving, start: int) -> float:
    rotation, translation = hla_transform(reference, moving)
    target = []
    source = []
    for left, right in zip(residues(reference, "P")[start - 1:start + 8], residues(moving, "P")[start - 1:start + 8]):
        for atom in ("N", "CA", "C", "O"):
            if atom in left and atom in right:
                target.append(left[atom].coord)
                source.append(right[atom].coord)
    fitted = np.asarray(source) @ rotation + translation
    return float(np.sqrt(np.mean(np.sum((np.asarray(target) - fitted) ** 2, axis=1))))


def local_contact(model, peptide_position: int) -> dict[str, object]:
    residue = residues(model, "P")[peptide_position - 1]
    sidechain = [atom for atom in residue if atom.element != "H" and atom.name not in BACKBONE]
    all_hla = [atom for chain in ("M", "N") for item in residues(model, chain) for atom in item if atom.element != "H"]
    min_sidechain = min((float(np.linalg.norm(atom.coord - other.coord)) for atom in sidechain for other in all_hla), default=float("inf"))
    charge_atoms = ACIDIC_OXYGENS.get(residue.resname, BASIC_NITROGENS.get(residue.resname, set()))
    opposite = BASIC_NITROGENS if residue.resname in ACIDIC_OXYGENS else ACIDIC_OXYGENS if residue.resname in BASIC_NITROGENS else {}
    charge_distances = []
    for atom in residue:
        if atom.name not in charge_atoms:
            continue
        for chain in ("M", "N"):
            for hla_residue in residues(model, chain):
                for other in hla_residue:
                    if other.name in opposite.get(hla_residue.resname, set()):
                        charge_distances.append(float(np.linalg.norm(atom.coord - other.coord)))
    nearest_charge = min(charge_distances, default=float("inf"))
    return {
        "mutated_residue_observed": seq1(residue.resname),
        "nearest_sidechain_hla_distance_A": round(min_sidechain, 4),
        "sidechain_hla_contact_within_4A": min_sidechain <= 4.0,
        "nearest_opposite_charge_atom_distance_A": "" if not charge_distances else round(nearest_charge, 4),
        "opposite_charge_atom_pair_within_4A": nearest_charge <= 4.0,
    }


def analyze(package: Path) -> dict[str, object]:
    parser = PDBParser(QUIET=True)
    manifest = read_rows(package / "pandora" / "candidate_runs_60_requesting_1200_models.csv")
    metadata = {row["sample_id"]: row for row in primary_panel() if row["primary_mutation_panel"]}
    results = package / "pandora" / "results"
    model_rows = []
    run_rows = []
    top_models = {}
    for run in manifest:
        directory = results / run["run_id"]
        score_file = directory / "molpdf_DOPE.tsv"
        scores = []
        for line in score_file.read_text(encoding="utf-8").splitlines():
            name, molpdf, dope = line.split("\t")
            scores.append((name, float(molpdf), float(dope)))
        if len(scores) != 20:
            raise ValueError(f"expected 20 models for {run['run_id']}; found {len(scores)}")
        meta = metadata[run["sample_id"]]
        current = []
        for name, molpdf, dope in scores:
            path = directory / name
            model = parser.get_structure(name, path)
            observed = {chain: sequence(residues(model, chain)) for chain in ("M", "N", "P")}
            exact = observed == {"M": DRA_SEQUENCE, "N": DRB1_1501_SEQUENCE, "P": run["peptide"]}
            contact = local_contact(model, int(meta["peptide_position_1_based"])) if meta["peptide_position_1_based"] else {
                "mutated_residue_observed": "",
                "nearest_sidechain_hla_distance_A": "",
                "sidechain_hla_contact_within_4A": "",
                "nearest_opposite_charge_atom_distance_A": "",
                "opposite_charge_atom_pair_within_4A": "",
            }
            record = {
                "run_id": run["run_id"], "sample_id": run["sample_id"],
                "register_label": run["register_label"], "core_start_1_based": run["core_start_1_based"],
                "forced_template_pdb": run["forced_template_pdb"], "model_name": name,
                "molpdf": molpdf, "dope": dope, "exact_MNP_sequence_qc_pass": exact,
                **contact, "model_is_biological_replicate": False,
            }
            current.append(record)
            model_rows.append(record)
        ranked = sorted(current, key=lambda row: float(row["molpdf"]))
        top_models[(run["sample_id"], run["register_label"], run["forced_template_pdb"])] = parser.get_structure(
            run["run_id"] + "_top", directory / ranked[0]["model_name"]
        )
        run_rows.append({
            "run_id": run["run_id"], "sample_id": run["sample_id"], "register_label": run["register_label"],
            "core_start_1_based": run["core_start_1_based"], "forced_template_pdb": run["forced_template_pdb"],
            "models": len(current), "all_exact_MNP_sequence_qc_pass": all(row["exact_MNP_sequence_qc_pass"] for row in current),
            "top_five_molpdf_median": round(median(float(row["molpdf"]) for row in ranked[:5]), 5),
            "all_model_molpdf_median": round(median(float(row["molpdf"]) for row in ranked), 5),
            "sidechain_contact_frequency": round(sum(bool(row.get("sidechain_hla_contact_within_4A")) for row in current) / 20, 3) if meta["peptide_position_1_based"] else "",
            "opposite_charge_contact_frequency": round(sum(bool(row.get("opposite_charge_atom_pair_within_4A")) for row in current) / 20, 3) if meta["peptide_position_1_based"] else "",
        })
    write_rows(results / "pandora_candidate_model_metrics_1200.csv", model_rows)
    write_rows(results / "pandora_candidate_run_summary_60.csv", run_rows)

    sample_rows = []
    for sample_id, meta in metadata.items():
        sample_runs = [row for row in run_rows if row["sample_id"] == sample_id]
        preferred = {}
        for template in ("1BX2", "6CQQ"):
            template_runs = [row for row in sample_runs if row["forced_template_pdb"] == template]
            preferred[template] = min(template_runs, key=lambda row: float(row["top_five_molpdf_median"]))["register_label"]
        expected_start = 4 if meta["arm"] == "BALF5" else 5
        cross_template = peptide_core_rmsd(
            top_models[(sample_id, "expected", "1BX2")],
            top_models[(sample_id, "expected", "6CQQ")],
            expected_start,
        )
        expected_runs = [row for row in sample_runs if row["register_label"] == "expected"]
        sample_rows.append({
            "sample_id": sample_id, "arm": meta["arm"],
            "preferred_register_1BX2": preferred["1BX2"], "preferred_register_6CQQ": preferred["6CQQ"],
            "expected_register_preferred_both_templates": preferred == {"1BX2": "expected", "6CQQ": "expected"},
            "expected_register_top_model_cross_template_core_rmsd_A": round(cross_template, 4),
            "cross_template_core_rmsd_le_2A": cross_template <= 2.0,
            "expected_register_sidechain_contact_frequency_1BX2": next(row["sidechain_contact_frequency"] for row in expected_runs if row["forced_template_pdb"] == "1BX2"),
            "expected_register_sidechain_contact_frequency_6CQQ": next(row["sidechain_contact_frequency"] for row in expected_runs if row["forced_template_pdb"] == "6CQQ"),
            "expected_register_opposite_charge_contact_frequency_1BX2": next(row["opposite_charge_contact_frequency"] for row in expected_runs if row["forced_template_pdb"] == "1BX2"),
            "expected_register_opposite_charge_contact_frequency_6CQQ": next(row["opposite_charge_contact_frequency"] for row in expected_runs if row["forced_template_pdb"] == "6CQQ"),
            "pandora_register_status": "supports_expected_register" if preferred == {"1BX2": "expected", "6CQQ": "expected"} and cross_template <= 2.0 else "register_or_template_confounded",
            "scores_are_affinities": False,
        })
    write_rows(results / "pandora_candidate_sample_summary_10.csv", sample_rows)
    summary = {
        "runs": len(run_rows), "models": len(model_rows),
        "all_exact_sequence_qc_pass": all(row["all_exact_MNP_sequence_qc_pass"] for row in run_rows),
        "samples_supporting_expected_register": sum(row["pandora_register_status"] == "supports_expected_register" for row in sample_rows),
        "models_are_biological_replicates": False,
        "interpretation": "structural method sensitivity only; molpdf/DOPE are not binding affinities",
    }
    (results / "PANDORA_CANDIDATE_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.package.resolve()), indent=2))


if __name__ == "__main__":
    main()
