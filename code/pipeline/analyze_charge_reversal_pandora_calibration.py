"""Score reciprocal PANDORA register calibration against experimental structures."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from build_tcell_library_v2 import DRA_SEQUENCE, DRB1_1501_SEQUENCE


BACKBONE = ("N", "CA", "C", "O")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, values: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def residues(structure, chain: str):
    return [residue for residue in structure[0][chain] if residue.id[0] == " "]


def sequence(residue_list) -> str:
    return "".join(seq1(residue.resname, custom_map={"MSE": "M"}) for residue in residue_list)


def alignment_pairs(reference_residues, model_residues, cap: int) -> list[tuple[object, object]]:
    reference_residues = reference_residues[:cap]
    model_residues = model_residues[:cap]
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(sequence(reference_residues), sequence(model_residues))[0]
    pairs = []
    for (r0, r1), (m0, m1) in zip(alignment.aligned[0], alignment.aligned[1]):
        for r_index, m_index in zip(range(r0, r1), range(m0, m1)):
            if "CA" in reference_residues[r_index] and "CA" in model_residues[m_index]:
                pairs.append((reference_residues[r_index], model_residues[m_index]))
    return pairs


def kabsch(moving: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    moving_center = moving.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (moving - moving_center).T @ (target - target_center)
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    translation = target_center - moving_center @ rotation
    return rotation, translation


def core_rmsd(reference, model, start: int) -> float:
    hla_pairs = alignment_pairs(residues(reference, "M"), residues(model, "M"), 90)
    hla_pairs += alignment_pairs(residues(reference, "N"), residues(model, "N"), 95)
    target = np.asarray([pair[0]["CA"].coord for pair in hla_pairs])
    moving = np.asarray([pair[1]["CA"].coord for pair in hla_pairs])
    rotation, translation = kabsch(moving, target)
    ref_peptide = residues(reference, "P")[start - 1:start + 8]
    model_peptide = residues(model, "P")[start - 1:start + 8]
    target_atoms = []
    moving_atoms = []
    for ref_residue, model_residue in zip(ref_peptide, model_peptide):
        for atom in BACKBONE:
            if atom in ref_residue and atom in model_residue:
                target_atoms.append(ref_residue[atom].coord)
                moving_atoms.append(model_residue[atom].coord)
    target_array = np.asarray(target_atoms)
    fitted = np.asarray(moving_atoms) @ rotation + translation
    return float(np.sqrt(np.mean(np.sum((target_array - fitted) ** 2, axis=1))))


def analyze(package: Path) -> dict[str, object]:
    parser = PDBParser(QUIET=True)
    manifest = rows(package / "pandora" / "calibration_runs_6_requesting_120_models.csv")
    results = package / "pandora" / "results"
    model_rows = []
    run_rows = []
    for run in manifest:
        directory = results / run["run_id"]
        scores = {}
        for line in (directory / "molpdf_DOPE.tsv").read_text(encoding="utf-8").splitlines():
            name, molpdf, dope = line.split("\t")
            scores[name] = (float(molpdf), float(dope))
        target_id = run["experimental_target_pdb"]
        reference_path = package / "pandora" / "templates" / "pmhc_only" / f"{target_id}_MNP_pandora.pdb"
        reference = parser.get_structure(f"ref_{target_id}", reference_path)
        expected_start = 5 if target_id == "1BX2" else 3
        current = []
        for name, (molpdf, dope) in scores.items():
            model_path = directory / name
            model = parser.get_structure(name, model_path)
            rmsd = core_rmsd(reference, model, expected_start)
            exact = {
                chain: sequence(residues(model, chain)) for chain in ("M", "N", "P")
            } == {"M": DRA_SEQUENCE, "N": DRB1_1501_SEQUENCE, "P": run["peptide"]}
            record = {
                "run_id": run["run_id"], "experimental_target_pdb": target_id,
                "forced_template_pdb": run["forced_template_pdb"], "register_label": run["register_label"],
                "model_name": name, "molpdf": molpdf, "dope": dope,
                "experimental_core_backbone_rmsd_A": round(rmsd, 4),
                "exact_MNP_sequence_qc_pass": exact,
                "model_is_biological_replicate": False,
            }
            current.append(record)
            model_rows.append(record)
        ranked = sorted(current, key=lambda value: float(value["molpdf"]))
        top_five = ranked[:5]
        run_rows.append({
            "run_id": run["run_id"], "experimental_target_pdb": target_id,
            "forced_template_pdb": run["forced_template_pdb"], "register_label": run["register_label"],
            "models": len(current), "top_five_min_core_backbone_rmsd_A": min(float(value["experimental_core_backbone_rmsd_A"]) for value in top_five),
            "top_five_median_core_backbone_rmsd_A": round(median(float(value["experimental_core_backbone_rmsd_A"]) for value in top_five), 4),
            "any_top_five_below_2A": any(float(value["experimental_core_backbone_rmsd_A"]) < 2.0 for value in top_five),
            "all_exact_MNP_sequence_qc_pass": all(bool(value["exact_MNP_sequence_qc_pass"]) for value in current),
            "template_leakage_check": "pass" if run["forced_template_pdb"] != target_id and not (directory / f"{target_id}.pdb").exists() else "fail",
        })
    write(results / "calibration_model_metrics_120.csv", model_rows)
    write(results / "calibration_register_summary_6.csv", run_rows)
    target_gates = []
    for target_id in ("1BX2", "6CQQ"):
        target_runs = [row for row in run_rows if row["experimental_target_pdb"] == target_id]
        expected = next(row for row in target_runs if row["register_label"] == "expected")
        decoys = [row for row in target_runs if row["register_label"] != "expected"]
        passed = (
            bool(expected["any_top_five_below_2A"])
            and all(float(expected["top_five_median_core_backbone_rmsd_A"]) < float(row["top_five_median_core_backbone_rmsd_A"]) for row in decoys)
            and all(row["template_leakage_check"] == "pass" for row in target_runs)
            and all(bool(row["all_exact_MNP_sequence_qc_pass"]) for row in target_runs)
            and all(int(row["models"]) == 20 for row in target_runs)
        )
        target_gates.append({
            "experimental_target_pdb": target_id,
            "expected_register_top_five_min_rmsd_A": expected["top_five_min_core_backbone_rmsd_A"],
            "expected_register_top_five_median_rmsd_A": expected["top_five_median_core_backbone_rmsd_A"],
            "both_adjacent_decoy_medians_worse": all(float(expected["top_five_median_core_backbone_rmsd_A"]) < float(row["top_five_median_core_backbone_rmsd_A"]) for row in decoys),
            "at_least_one_expected_top_five_model_below_2A": expected["any_top_five_below_2A"],
            "template_leakage_checks_pass": all(row["template_leakage_check"] == "pass" for row in target_runs),
            "exact_sequence_checks_pass": all(bool(row["all_exact_MNP_sequence_qc_pass"]) for row in target_runs),
            "calibration_status": "pass" if passed else "fail",
        })
    write(results / "calibration_gate_by_target_2.csv", target_gates)
    overall = {
        "calibration_status": "pass" if all(row["calibration_status"] == "pass" for row in target_gates) else "fail",
        "candidate_modeling_unlocked": all(row["calibration_status"] == "pass" for row in target_gates),
        "models_scored": len(model_rows),
        "target_gates": target_gates,
        "rule": "each expected register needs a top-five model below 2 A, a lower top-five median than both adjacent decoys, complete 20-model runs, and no target-template leakage",
    }
    (results / "CALIBRATION_GATE.json").write_text(json.dumps(overall, indent=2) + "\n", encoding="utf-8")
    return overall


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.package.resolve()), indent=2))


if __name__ == "__main__":
    main()
