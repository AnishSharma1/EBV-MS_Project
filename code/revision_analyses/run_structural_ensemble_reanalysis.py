#!/usr/bin/env python3
"""Re-express parent/nested AF3 ensembles without treating 5x5 RMSDs as n=25."""

from __future__ import annotations

import argparse
import csv
import math
import shlex
from collections import defaultdict
from pathlib import Path

import numpy as np


# The frozen focused-register analysis used peptide C-alpha coordinates and
# described them as backbone RMSD. Reuse that exact coordinate definition so
# this reanalysis reproduces the archived 5x5 summaries rather than silently
# changing the metric to all-backbone atoms.
BACKBONE = ("CA",)
CLAIM_BOUNDARY = (
    "Descriptive AlphaFold coordinate sensitivity only; models are not biological replicates "
    "and do not establish binding, register identity, presentation, or cross-reactivity."
)


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[2]
    focused = Path(
        "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/"
        "processed/taldo1_focused_register_2026-09-05/af3_analysis"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metrics", type=Path, default=focused / "sample_metrics_65.csv")
    parser.add_argument("--parent-nested", type=Path, default=focused / "parent_nested_ensemble_comparisons.csv")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_3/computational_tightening_2026-09-15",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_mmcif_atoms(path: Path) -> dict[str, dict[int, dict[str, np.ndarray]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers: list[str] = []
    start = None
    for index, line in enumerate(lines):
        if line.strip() == "loop_":
            candidate = []
            cursor = index + 1
            while cursor < len(lines) and lines[cursor].startswith("_atom_site."):
                candidate.append(lines[cursor].strip())
                cursor += 1
            if candidate:
                headers = candidate
                start = cursor
                break
    if start is None:
        raise ValueError(f"atom_site loop not found: {path}")
    index = {name: offset for offset, name in enumerate(headers)}
    required = [
        "_atom_site.group_PDB", "_atom_site.label_atom_id", "_atom_site.label_asym_id",
        "_atom_site.label_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
    ]
    missing = [name for name in required if name not in index]
    if missing:
        raise ValueError(f"missing atom_site fields {missing}: {path}")
    output: dict[str, dict[int, dict[str, np.ndarray]]] = defaultdict(lambda: defaultdict(dict))
    for line in lines[start:]:
        if line.startswith("#") or line.startswith("loop_") or line.startswith("_"):
            break
        if not line.strip():
            continue
        values = shlex.split(line, posix=True)
        if len(values) != len(headers):
            raise ValueError(f"unexpected atom row width {len(values)} != {len(headers)}: {path}")
        if values[index["_atom_site.group_PDB"]] != "ATOM":
            continue
        chain = values[index["_atom_site.label_asym_id"]]
        residue = int(values[index["_atom_site.label_seq_id"]])
        atom = values[index["_atom_site.label_atom_id"]]
        output[chain][residue][atom] = np.asarray([
            float(values[index["_atom_site.Cartn_x"]]),
            float(values[index["_atom_site.Cartn_y"]]),
            float(values[index["_atom_site.Cartn_z"]]),
        ])
    return output


def coords(atoms, chain: str, residues: list[int], atom_names: tuple[str, ...]) -> np.ndarray:
    output = []
    for residue in residues:
        for atom in atom_names:
            if atom not in atoms[chain][residue]:
                raise ValueError(f"missing {chain}:{residue}:{atom}")
            output.append(atoms[chain][residue][atom])
    return np.asarray(output, dtype=float)


def hla_fit(reference, mobile):
    reference_hla = np.vstack([
        coords(reference, "A", sorted(reference["A"])[:85], ("CA",)),
        coords(reference, "B", sorted(reference["B"])[:85], ("CA",)),
    ])
    mobile_hla = np.vstack([
        coords(mobile, "A", sorted(mobile["A"])[:85], ("CA",)),
        coords(mobile, "B", sorted(mobile["B"])[:85], ("CA",)),
    ])
    ref_center = reference_hla.mean(axis=0)
    mob_center = mobile_hla.mean(axis=0)
    covariance = (mobile_hla - mob_center).T @ (reference_hla - ref_center)
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return rotation, ref_center - mob_center @ rotation


def peptide_rmsd(reference, mobile, reference_residues: list[int], mobile_residues: list[int]) -> float:
    rotation, translation = hla_fit(reference, mobile)
    ref_coords = coords(reference, "C", reference_residues, BACKBONE)
    mob_coords = coords(mobile, "C", mobile_residues, BACKBONE) @ rotation + translation
    return float(np.sqrt(np.mean(np.sum((ref_coords - mob_coords) ** 2, axis=1))))


def pairwise(models: list[dict], residues: list[int]) -> tuple[list[dict], int]:
    rows = []
    values: dict[tuple[int, int], float] = {}
    for left in range(len(models)):
        for right in range(left + 1, len(models)):
            value = peptide_rmsd(models[left]["atoms"], models[right]["atoms"], residues, residues)
            values[left, right] = value
            rows.append({
                "left_model_index": models[left]["model_index"],
                "right_model_index": models[right]["model_index"],
                "rmsd_A": f"{value:.6f}",
            })
    mean_by_model = {}
    for index in range(len(models)):
        model_values = [value for (left, right), value in values.items() if index in (left, right)]
        mean_by_model[index] = statistics_mean(model_values)
    medoid = min(mean_by_model, key=lambda index: (mean_by_model[index], models[index]["model_index"]))
    return rows, medoid


def statistics_mean(values: list[float]) -> float:
    return sum(values) / len(values)


def summary(values: list[float]) -> tuple[float, float, float]:
    return min(values), float(np.median(values)), max(values)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sample_rows = read_csv(args.sample_metrics)
    by_job: dict[str, list[dict]] = defaultdict(list)
    for row in sample_rows:
        path = Path(row["cif"])
        if not path.exists():
            raise FileNotFoundError(path)
        by_job[row["job"]].append({
            "model_index": int(row["model_index"]),
            "path": str(path),
            "atoms": parse_mmcif_atoms(path),
        })
    for models in by_job.values():
        models.sort(key=lambda row: row["model_index"])

    within_rows = []
    cross_rows = []
    corrected_rows = []
    for comparison in read_csv(args.parent_nested):
        parent_models = by_job[comparison["parent_job"]]
        nested_models = by_job[comparison["nested_job"]]
        start = int(comparison["nested_start_in_parent_1based"])
        parent_residues = list(range(start, start + 11))
        nested_residues = list(range(1, 12))
        parent_full_residues = list(range(1, len(comparison["parent_peptide"]) + 1))
        nested_full_residues = list(range(1, len(comparison["nested_peptide"]) + 1))
        parent_pairs, parent_medoid = pairwise(parent_models, parent_full_residues)
        nested_pairs, nested_medoid = pairwise(nested_models, nested_full_residues)
        for row in parent_pairs:
            within_rows.append({"candidate": comparison["candidate"], "arm": comparison["arm"], "ensemble": "parent", **row, "claim_boundary": CLAIM_BOUNDARY})
        for row in nested_pairs:
            within_rows.append({"candidate": comparison["candidate"], "arm": comparison["arm"], "ensemble": "nested", **row, "claim_boundary": CLAIM_BOUNDARY})
        cross_values = []
        for parent in parent_models:
            for nested in nested_models:
                value = peptide_rmsd(parent["atoms"], nested["atoms"], parent_residues, nested_residues)
                cross_values.append(value)
                cross_rows.append({
                    "candidate": comparison["candidate"],
                    "arm": comparison["arm"],
                    "parent_model_index": parent["model_index"],
                    "nested_model_index": nested["model_index"],
                    "rmsd_A": f"{value:.6f}",
                    "comparison_is_independent": False,
                    "claim_boundary": CLAIM_BOUNDARY,
                })
        parent_values = [float(row["rmsd_A"]) for row in parent_pairs]
        nested_values = [float(row["rmsd_A"]) for row in nested_pairs]
        medoid_value = peptide_rmsd(
            parent_models[parent_medoid]["atoms"], nested_models[nested_medoid]["atoms"],
            parent_residues, nested_residues,
        )
        parent_min, parent_median, parent_max = summary(parent_values)
        nested_min, nested_median, nested_max = summary(nested_values)
        cross_min, cross_median, cross_max = summary(cross_values)
        corrected_rows.append({
            "candidate": comparison["candidate"],
            "arm": comparison["arm"],
            "parent_model_count": 5,
            "nested_model_count": 5,
            "underlying_model_count": 10,
            "within_parent_unique_pairs": 10,
            "within_parent_rmsd_min_A": f"{parent_min:.6f}",
            "within_parent_rmsd_median_A": f"{parent_median:.6f}",
            "within_parent_rmsd_max_A": f"{parent_max:.6f}",
            "within_nested_unique_pairs": 10,
            "within_nested_rmsd_min_A": f"{nested_min:.6f}",
            "within_nested_rmsd_median_A": f"{nested_median:.6f}",
            "within_nested_rmsd_max_A": f"{nested_max:.6f}",
            "cross_input_dependent_pairs": 25,
            "cross_input_rmsd_min_A": f"{cross_min:.6f}",
            "cross_input_rmsd_median_A": f"{cross_median:.6f}",
            "cross_input_rmsd_max_A": f"{cross_max:.6f}",
            "parent_medoid_model_index": parent_models[parent_medoid]["model_index"],
            "nested_medoid_model_index": nested_models[nested_medoid]["model_index"],
            "medoid_to_medoid_shared_11mer_rmsd_A": f"{medoid_value:.6f}",
            "matched_seed_comparison_status": "not_evaluated_no_verified_cross_job_seed_identity",
            "peptide_coordinate_definition": "C-alpha atoms only",
            "statistical_status": "descriptive_ensembles_not_biological_replicates",
            "claim_boundary": CLAIM_BOUNDARY,
        })

    write_csv(args.output_dir / "structural_within_ensemble_pairwise_rmsd.csv", within_rows)
    write_csv(args.output_dir / "structural_parent_nested_full_5x5_matrices.csv", cross_rows)
    write_csv(args.output_dir / "structural_ensemble_corrected_with_medoids.csv", corrected_rows)
    print(f"wrote {len(corrected_rows)} corrected ensemble rows and {len(cross_rows)} cross-input comparisons")


if __name__ == "__main__":
    main()
