"""Analyze the focused TALDO1 parent/nested AlphaFold Server pMHC ensemble.

The outputs describe model confidence, peptide-HLA contacts, within-job pose
consistency, and parent-versus-nested structural sensitivity. They do not infer
presentation, TCR recognition, cross-reactivity, or disease mechanism.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from statistics import median

import numpy as np

from analyze_af3_pmhc_downloads import (
    ca_coordinates,
    kabsch,
    parse_mmcif,
    peptide_hla_metrics,
    sequence,
)


SOURCE = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads/"
    "folds_2026_09_05_16_53"
)
OUT = Path(__file__).resolve().parents[1] / "processed" / "taldo1_focused_register_2026-09-05" / "af3_analysis"


PARENT_NESTED = [
    ("HY13_SEQ_02", "ebv", "taldo1_hy13_drb1_1303_ebv_parent", "taldo1_hy13_drb1_1303_ebv_nested"),
    ("HY13_SEQ_02", "self", "taldo1_hy13_drb1_1303_self_parent", "taldo1_hy13_drb1_1303_self_nested"),
    ("HY15_SEQ_02", "ebv", "taldo1_hy15_drb1_1501_ebv_parent", "taldo1_hy15_drb1_1501_ebv_nested"),
    ("HY15_SEQ_02", "self", "taldo1_hy15_drb1_1501_self_parent", "taldo1_hy15_drb1_1501_self_nested"),
    ("HY15_SEQ_02_DRB5", "ebv", "taldo1_hy15_drb5_0101_ebv_parent", "taldo1_hy15_drb5_0101_ebv_nested"),
    ("HY15_SEQ_02_DRB5", "self_mixmhc", "taldo1_hy15_drb5_0101_self_parent", "taldo1_hy15_drb5_0101_self_nested_mixmhc"),
    ("HY15_SEQ_02_DRB5", "self_iedb", "taldo1_hy15_drb5_0101_self_parent", "taldo1_hy15_drb5_0101_self_nested_iedb"),
]


SAME_PEPTIDE_CROSS_HLA = [
    ("ebv_parent_1303_vs_1501", "taldo1_hy13_drb1_1303_ebv_parent", "taldo1_hy15_drb1_1501_ebv_parent"),
    ("ebv_parent_1501_vs_drb5", "taldo1_hy15_drb1_1501_ebv_parent", "taldo1_hy15_drb5_0101_ebv_parent"),
    ("ebv_parent_1303_vs_drb5", "taldo1_hy13_drb1_1303_ebv_parent", "taldo1_hy15_drb5_0101_ebv_parent"),
    ("self_parent_1501_vs_drb5", "taldo1_hy15_drb1_1501_self_parent", "taldo1_hy15_drb5_0101_self_parent"),
    ("ebv_nested_1303_vs_drb5", "taldo1_hy13_drb1_1303_ebv_nested", "taldo1_hy15_drb5_0101_ebv_nested"),
]


def index_from_name(path: Path) -> int:
    match = re.search(r"_(\d+)\.(?:cif|json)$", path.name)
    if not match:
        raise ValueError(f"Cannot recover model index from {path.name}")
    return int(match.group(1))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_jobs() -> dict[str, dict[str, object]]:
    jobs: dict[str, dict[str, object]] = {}
    for directory in sorted(path for path in SOURCE.iterdir() if path.is_dir()):
        request_paths = list(directory.glob("*_job_request.json"))
        if len(request_paths) != 1:
            raise ValueError(f"Expected one request in {directory}")
        request = json.loads(request_paths[0].read_text(encoding="utf-8"))[0]
        name = str(request["name"])
        requested = [entry["proteinChain"]["sequence"] for entry in request["sequences"]]
        models: dict[int, dict[str, object]] = {}
        for cif_path in sorted(directory.glob("*_model_*.cif"), key=index_from_name):
            index = index_from_name(cif_path)
            summary_path = next(directory.glob(f"*_summary_confidences_{index}.json"))
            parsed = parse_mmcif(cif_path)
            observed = [sequence(parsed.get(chain, [])) for chain in ("A", "B", "C")]
            if observed != requested:
                raise ValueError(f"Sequence mismatch in {name}, model {index}")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            models[index] = {"structure": parsed, "summary": summary, "cif": str(cif_path)}
        if len(models) != 5:
            raise ValueError(f"Expected five models for {name}; found {len(models)}")
        jobs[name] = {"request": request, "requested": requested, "models": models}
    if len(jobs) != 13:
        raise ValueError(f"Expected 13 jobs; found {len(jobs)}")
    return jobs


def peptide_rmsd(reference: dict[str, object], moving: dict[str, object], reference_slice: slice | None = None) -> float:
    ref_structure = reference["structure"]
    mov_structure = moving["structure"]
    ref_hla = np.vstack([ca_coordinates(ref_structure[chain][:85]) for chain in ("A", "B")])
    mov_hla = np.vstack([ca_coordinates(mov_structure[chain][:85]) for chain in ("A", "B")])
    rotation, translation = kabsch(mov_hla, ref_hla)
    ref_peptide = ca_coordinates(ref_structure["C"])
    if reference_slice is not None:
        ref_peptide = ref_peptide[reference_slice]
    mov_peptide = ca_coordinates(mov_structure["C"]) @ rotation + translation
    if len(ref_peptide) != len(mov_peptide):
        raise ValueError("Peptide comparison lengths differ")
    return float(np.sqrt(np.mean(np.sum((mov_peptide - ref_peptide) ** 2, axis=1))))


def summarize(values: list[float]) -> dict[str, float]:
    return {
        "median_A": round(median(values), 3),
        "min_A": round(min(values), 3),
        "max_A": round(max(values), 3),
    }


def main() -> None:
    jobs = load_jobs()
    OUT.mkdir(parents=True, exist_ok=True)

    sample_rows: list[dict[str, object]] = []
    job_rows: list[dict[str, object]] = []
    for name, job in sorted(jobs.items()):
        per_job = []
        for index, model in sorted(job["models"].items()):
            structure = model["structure"]
            summary = model["summary"]
            metrics = peptide_hla_metrics(structure["C"], [structure["A"], structure["B"]])
            peptide_pair_iptm = [float(summary["chain_pair_iptm"][2][i]) for i in (0, 1)]
            peptide_pair_pae = [float(summary["chain_pair_pae_min"][2][i]) for i in (0, 1)]
            row = {
                "job": name,
                "model_index": index,
                "peptide": job["requested"][2],
                "ranking_score": summary.get("ranking_score", ""),
                "iptm": summary.get("iptm", ""),
                "ptm": summary.get("ptm", ""),
                "peptide_chain_iptm": summary["chain_iptm"][2],
                "peptide_hla_pair_iptm_mean": round(sum(peptide_pair_iptm) / 2, 3),
                "peptide_hla_pair_pae_min_mean_A": round(sum(peptide_pair_pae) / 2, 3),
                "has_clash": summary.get("has_clash", ""),
                **metrics,
                "cif": model["cif"],
            }
            sample_rows.append(row)
            per_job.append(row)
        within = []
        for left in range(5):
            for right in range(left + 1, 5):
                within.append(peptide_rmsd(job["models"][left], job["models"][right]))
        selected = max(per_job, key=lambda row: float(row["ranking_score"]))
        peptide_selected = max(per_job, key=lambda row: float(row["peptide_chain_iptm"]))
        job_rows.append({
            "job": name,
            "peptide": job["requested"][2],
            "models": 5,
            "selected_model_index": selected["model_index"],
            "selected_ranking_score": selected["ranking_score"],
            "selected_iptm": selected["iptm"],
            "server_selected_peptide_chain_iptm": selected["peptide_chain_iptm"],
            "best_peptide_model_index": peptide_selected["model_index"],
            "best_peptide_chain_iptm": peptide_selected["peptide_chain_iptm"],
            "median_peptide_chain_iptm": round(median(float(row["peptide_chain_iptm"]) for row in per_job), 3),
            "median_peptide_hla_pair_iptm_mean": round(median(float(row["peptide_hla_pair_iptm_mean"]) for row in per_job), 3),
            "median_peptide_hla_pair_pae_min_mean_A": round(median(float(row["peptide_hla_pair_pae_min_mean_A"]) for row in per_job), 3),
            "selected_peptide_mean_plddt": selected["peptide_mean_plddt"],
            "median_peptide_mean_plddt": round(median(float(row["peptide_mean_plddt"]) for row in per_job), 2),
            "median_peptide_contact_residues": median(float(row["peptide_residues_with_any_hla_contact"]) for row in per_job),
            "within_job_pairwise_peptide_rmsd_median_A": round(median(within), 3),
            "within_job_pairwise_peptide_rmsd_max_A": round(max(within), 3),
            "any_clash": any(bool(row["has_clash"]) for row in per_job),
        })

    register_rows: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []
    parent_model_best_window_rows: list[dict[str, object]] = []
    for candidate, arm, parent_name, nested_name in PARENT_NESTED:
        parent, nested = jobs[parent_name], jobs[nested_name]
        parent_peptide, nested_peptide = parent["requested"][2], nested["requested"][2]
        start = parent_peptide.index(nested_peptide)
        rmsds = []
        for parent_index, parent_model in parent["models"].items():
            for nested_index, nested_model in nested["models"].items():
                rmsds.append(peptide_rmsd(parent_model, nested_model, slice(start, start + len(nested_peptide))))
        register_rows.append({
            "candidate": candidate,
            "arm": arm,
            "parent_job": parent_name,
            "nested_job": nested_name,
            "parent_peptide": parent_peptide,
            "nested_peptide": nested_peptide,
            "nested_start_in_parent_1based": start + 1,
            "model_pair_comparisons": len(rmsds),
            **summarize(rmsds),
            "interpretation": "Parent-versus-nested peptide pose sensitivity after equivalent HLA-groove alignment; descriptive modeling result only.",
        })
        for window_start in range(len(parent_peptide) - len(nested_peptide) + 1):
            window_rmsds = [
                peptide_rmsd(parent_model, nested_model, slice(window_start, window_start + len(nested_peptide)))
                for parent_model in parent["models"].values()
                for nested_model in nested["models"].values()
            ]
            window_rows.append({
                "candidate": candidate,
                "arm": arm,
                "nested_job": nested_name,
                "nested_peptide": nested_peptide,
                "parent_window_start_1based": window_start + 1,
                "parent_window": parent_peptide[window_start:window_start + len(nested_peptide)],
                "is_sequence_identical_hypothesis": parent_peptide[window_start:window_start + len(nested_peptide)] == nested_peptide,
                "model_pair_comparisons": len(window_rmsds),
                **summarize(window_rmsds),
            })
        for parent_index, parent_model in parent["models"].items():
            choices = []
            for window_start in range(len(parent_peptide) - len(nested_peptide) + 1):
                values = [
                    peptide_rmsd(parent_model, nested_model, slice(window_start, window_start + len(nested_peptide)))
                    for nested_model in nested["models"].values()
                ]
                choices.append((median(values), window_start, values))
            best_median, best_start, best_values = min(choices)
            parent_model_best_window_rows.append({
                "candidate": candidate,
                "arm": arm,
                "nested_job": nested_name,
                "parent_model_index": parent_index,
                "parent_model_ranking_score": parent_model["summary"].get("ranking_score", ""),
                "best_parent_window_start_1based": best_start + 1,
                "best_parent_window": parent_peptide[best_start:best_start + len(nested_peptide)],
                "best_window_median_rmsd_A": round(best_median, 3),
                "best_window_min_rmsd_A": round(min(best_values), 3),
                "best_window_max_rmsd_A": round(max(best_values), 3),
            })

    cross_hla_rows: list[dict[str, object]] = []
    for comparison, left_name, right_name in SAME_PEPTIDE_CROSS_HLA:
        left, right = jobs[left_name], jobs[right_name]
        if left["requested"][2] != right["requested"][2]:
            raise ValueError(f"Cross-HLA peptide mismatch for {comparison}")
        rmsds = [
            peptide_rmsd(left_model, right_model)
            for left_model in left["models"].values()
            for right_model in right["models"].values()
        ]
        cross_hla_rows.append({
            "comparison": comparison,
            "left_job": left_name,
            "right_job": right_name,
            "peptide": left["requested"][2],
            "model_pair_comparisons": len(rmsds),
            **summarize(rmsds),
            "interpretation": "Same-peptide pose sensitivity across remodeled HLA complexes; descriptive modeling result only.",
        })

    write_csv(OUT / "sample_metrics_65.csv", sample_rows)
    write_csv(OUT / "job_summary_13.csv", job_rows)
    write_csv(OUT / "parent_nested_ensemble_comparisons.csv", register_rows)
    write_csv(OUT / "parent_window_scan.csv", window_rows)
    write_csv(OUT / "parent_model_best_window.csv", parent_model_best_window_rows)
    write_csv(OUT / "same_peptide_cross_hla_comparisons.csv", cross_hla_rows)
    summary = {
        "source": str(SOURCE),
        "jobs": len(jobs),
        "models": len(sample_rows),
        "all_sequence_layout_checks_passed": True,
        "parent_nested_comparisons": register_rows,
        "same_peptide_cross_hla_comparisons": cross_hla_rows,
        "claim_boundary": "Computational pMHC structure sensitivity only; not evidence of presentation, TCR binding, cross-reactivity, molecular mimicry, or MS mechanism.",
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    submitted_jobs = [jobs[name]["request"] for name in sorted(jobs)]
    (OUT / "submitted_jobs_reconstructed_13.json").write_text(
        json.dumps(submitted_jobs, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "README.md").write_text(
        "# TALDO1 focused-register AF3 analysis\n\n"
        "All 13 jobs and all 65 models passed exact three-chain sequence checks. "
        "The analysis uses all five models per job. Parent-versus-nested RMSD compares the shared 11-mer "
        "after fitting equivalent HLA-groove C-alpha atoms; cross-HLA rows compare only identical peptide sequences.\n\n"
        "These are descriptive model-sensitivity outputs. They do not establish peptide presentation, TCR binding, "
        "cross-reactivity, molecular mimicry, or an MS mechanism.\n",
        encoding="utf-8",
    )
    print(f"Analyzed {len(jobs)} jobs and {len(sample_rows)} models into {OUT}")


if __name__ == "__main__":
    main()
