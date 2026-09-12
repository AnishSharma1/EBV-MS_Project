"""Evaluate the two AlphaFold 3 OB.1A12 TCR-pMHC calibration jobs.

The AF3 request chain order is TCRa, TCRb, peptide, DRA, DRB. Models are
canonicalized to DRA, DRB, peptide, TCRa, TCRb and scored against the
experimental 1YMM and 2WBJ complexes. DockQ treats pMHC as receptor and the
paired TCR as ligand.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path
from statistics import median

import gemmi

from evaluate_ob1a12_ternary_models import (
    fitted_chain_rmsd,
    fit_to_reference,
    mean_plddt,
    parse_pdb,
    sequence,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads/"
    "folds_2026_09_05_22_09"
)
OUT = ROOT / "processed" / "olivia_literature_and_tcr_audit_2026-09-05" / "af3_tcr_retry_analysis"
OLD = Path("/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package/results_analysis/dockq_calibration")

CASES = {
    "tcr_1ymm": {
        "reference_id": "1YMM",
        "native": OLD / "CAL_1YMM" / "CAL_1YMM.native.canonical.pdb",
    },
    "tcr_2wbj": {
        "reference_id": "2WBJ",
        "native": OLD / "CAL_2WBJ_EXCLUDED" / "CAL_2WBJ_EXCLUDED.native.canonical.pdb",
    },
}

MODEL_TO_CANONICAL = {"A": "D", "B": "E", "C": "C", "D": "A", "E": "B"}


def quality_class(value: float) -> str:
    if value < 0.23:
        return "Incorrect"
    if value < 0.49:
        return "Acceptable"
    if value < 0.80:
        return "Medium"
    return "High"


def atom_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.startswith("ATOM")]


def canonicalize_af3(cif: Path, target: Path, temporary: Path) -> None:
    structure = gemmi.read_structure(str(cif))
    structure.write_pdb(str(temporary))
    buckets = {chain: [] for chain in "ABCDE"}
    for line in atom_lines(temporary):
        canonical = MODEL_TO_CANONICAL.get(line[21])
        if canonical:
            buckets[canonical].append(f"{line[:21]}{canonical}{line[22:]}")
    output = []
    for chain in "ABCDE":
        output.extend(buckets[chain])
        output.append("TER")
    output.append("END")
    target.write_text("\n".join(output) + "\n", encoding="utf-8")


def write_grouped(source: Path, target: Path) -> None:
    residue_ids = {"R": {}, "L": {}}
    counts = {"R": 0, "L": 0}
    output = []
    last_group = None
    for line in atom_lines(source):
        source_chain = line[21]
        group = "R" if source_chain in {"A", "B", "C"} else "L"
        if last_group is not None and group != last_group:
            output.append("TER")
        last_group = group
        old_id = (source_chain, line[22:27])
        if old_id not in residue_ids[group]:
            counts[group] += 1
            residue_ids[group][old_id] = counts[group]
        number = residue_ids[group][old_id]
        output.append(f"{line[:21]}{group}{number:4d} {line[27:]}")
    output.extend(["TER", "END"])
    target.write_text("\n".join(output) + "\n", encoding="utf-8")


def dockq_score(executable: Path, model: Path, native: Path, output_json: Path) -> dict[str, float]:
    subprocess.run(
        [str(executable), str(model), str(native), "--mapping", "RL:RL", "--json", str(output_json)],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    result = next(iter(payload["best_result"].values()))
    return {key: float(result[key]) for key in ("DockQ", "iRMSD", "LRMSD", "fnat", "fnonnat", "clashes")}


def template_ids(directory: Path) -> list[str]:
    ids = set()
    for path in directory.glob("templates/*_template_hit_*_chains_*.cif"):
        match = re.search(r"^_entry\.id\s+(\S+)", path.read_text(encoding="utf-8", errors="ignore"), re.MULTILINE)
        if match:
            ids.add(match.group(1).upper())
    return sorted(ids)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dockq", type=Path, required=True)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    template_rows = []

    for folder_name, spec in CASES.items():
        folder = SOURCE / folder_name
        native = spec["native"]
        if not native.exists():
            raise FileNotFoundError(native)
        request_path = next(folder.glob("*_job_request.json"))
        request = json.loads(request_path.read_text(encoding="utf-8"))[0]
        requested = [entry["proteinChain"]["sequence"] for entry in request["sequences"]]
        ids = template_ids(folder)
        template_rows.append({
            "job": request["name"],
            "reference_id": spec["reference_id"],
            "template_ids": ";".join(ids),
            "direct_reference_present": spec["reference_id"] in ids,
            "interpretation": "Per-chain AF3 template inventory; absence of the native ID does not make this equivalent to a TCRmodel2 template-exclusion run.",
        })

        native_grouped = OUT / f"{folder_name}.native.grouped.pdb"
        write_grouped(native, native_grouped)
        native_parsed = parse_pdb(native)
        cifs = sorted(folder.glob("*_model_*.cif"), key=lambda path: int(re.search(r"_(\d+)\.cif$", path.name).group(1)))
        if len(cifs) != 5:
            raise ValueError(f"Expected five models in {folder}")
        for cif in cifs:
            model_index = int(re.search(r"_(\d+)\.cif$", cif.name).group(1))
            summary_path = next(folder.glob(f"*_summary_confidences_{model_index}.json"))
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            canonical = OUT / f"{folder_name}.model_{model_index}.canonical.pdb"
            temporary = OUT / f"{folder_name}.model_{model_index}.raw.pdb"
            canonicalize_af3(cif, canonical, temporary)
            temporary.unlink()
            parsed = parse_pdb(canonical)
            observed_sequences = [sequence(parsed[chain]) for chain in ("D", "E", "C", "A", "B")]
            expected_lengths = [len(requested[i]) for i in (0, 1, 2, 3, 4)]
            observed_lengths = [len(item) for item in observed_sequences]
            if observed_lengths != expected_lengths:
                raise ValueError(f"Chain length mismatch for {folder_name} model {model_index}: {observed_lengths}")
            if observed_sequences != requested:
                raise ValueError(f"Chain sequence mismatch for {folder_name} model {model_index}")

            grouped = OUT / f"{folder_name}.model_{model_index}.grouped.pdb"
            write_grouped(canonical, grouped)
            dockq_json = OUT / f"{folder_name}.model_{model_index}.dockq.json"
            dockq = dockq_score(args.dockq, grouped, native_grouped, dockq_json)

            rotation, translation, hla_rmsd = fit_to_reference(parsed, native_parsed, ["A", "B"])
            peptide_rmsd = fitted_chain_rmsd(parsed["C"], native_parsed["C"], rotation, translation)
            tcra_rmsd = fitted_chain_rmsd(parsed["D"], native_parsed["D"], rotation, translation)
            tcrb_rmsd = fitted_chain_rmsd(parsed["E"], native_parsed["E"], rotation, translation)
            pair_iptm = summary["chain_pair_iptm"]
            pair_pae = summary["chain_pair_pae_min"]
            tcr_pmhc_pairs = [(i, j) for i in (0, 1) for j in (2, 3, 4)]
            peptide_tcr_pairs = [(0, 2), (1, 2)]
            rows.append({
                "job": request["name"],
                "reference_id": spec["reference_id"],
                "model_index": model_index,
                "sequences_match_request": True,
                "ranking_score": summary["ranking_score"],
                "global_iptm": summary["iptm"],
                "tcr_pmhc_pair_iptm_mean": round(sum(pair_iptm[i][j] for i, j in tcr_pmhc_pairs) / len(tcr_pmhc_pairs), 4),
                "peptide_tcr_pair_iptm_mean": round(sum(pair_iptm[i][j] for i, j in peptide_tcr_pairs) / len(peptide_tcr_pairs), 4),
                "tcr_pmhc_pair_pae_min_mean_A": round(sum(pair_pae[i][j] for i, j in tcr_pmhc_pairs) / len(tcr_pmhc_pairs), 4),
                "tcr_alpha_mean_plddt": mean_plddt(parsed["D"]),
                "tcr_beta_mean_plddt": mean_plddt(parsed["E"]),
                "peptide_mean_plddt": mean_plddt(parsed["C"]),
                "hla_fit_rmsd_A": round(float(hla_rmsd), 3),
                "peptide_rmsd_after_hla_fit_A": round(float(peptide_rmsd), 3) if peptide_rmsd is not None else "",
                "tcr_alpha_rmsd_after_hla_fit_A": round(float(tcra_rmsd), 3) if tcra_rmsd is not None else "",
                "tcr_beta_rmsd_after_hla_fit_A": round(float(tcrb_rmsd), 3) if tcrb_rmsd is not None else "",
                **{key: round(value, 4) for key, value in dockq.items()},
                "CAPRI_class": quality_class(dockq["DockQ"]),
                "canonical_pdb": str(canonical),
            })

    write_csv(OUT / "af3_tcr_retry_model_metrics_10.csv", rows)
    write_csv(OUT / "af3_tcr_retry_template_audit.csv", template_rows)
    case_summary = []
    for job in sorted({row["job"] for row in rows}):
        selected = [row for row in rows if row["job"] == job]
        case_summary.append({
            "job": job,
            "models": len(selected),
            "dockq_median": round(median(float(row["DockQ"]) for row in selected), 4),
            "dockq_min": min(float(row["DockQ"]) for row in selected),
            "dockq_max": max(float(row["DockQ"]) for row in selected),
            "capri_classes": ";".join(sorted({str(row["CAPRI_class"]) for row in selected})),
            "ranking_score_median": round(median(float(row["ranking_score"]) for row in selected), 4),
            "tcr_pmhc_pair_iptm_median": round(median(float(row["tcr_pmhc_pair_iptm_mean"]) for row in selected), 4),
            "peptide_tcr_pair_iptm_median": round(median(float(row["peptide_tcr_pair_iptm_mean"]) for row in selected), 4),
            "hla_fit_rmsd_median_A": round(median(float(row["hla_fit_rmsd_A"]) for row in selected), 3),
            "tcr_alpha_rmsd_median_A": round(median(float(row["tcr_alpha_rmsd_after_hla_fit_A"]) for row in selected), 3),
            "tcr_beta_rmsd_median_A": round(median(float(row["tcr_beta_rmsd_after_hla_fit_A"]) for row in selected), 3),
            "peptide_rmsd_median_A": round(median(float(row["peptide_rmsd_after_hla_fit_A"]) for row in selected), 3),
        })
    write_csv(OUT / "af3_tcr_retry_case_summary.csv", case_summary)
    print(f"Scored {len(rows)} models into {OUT}")


if __name__ == "__main__":
    main()
