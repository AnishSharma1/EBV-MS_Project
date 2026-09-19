"""Create the dated, auditable computational-recovery package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import urllib.request
from itertools import combinations
from pathlib import Path
from typing import Any

from charge_reversal_binding_pilot import checksum_rows, write_csv, write_json
from charge_reversal_recovery import (
    AF3_CONDITIONS,
    CLAIM_BOUNDARY,
    audit_baseline_templates,
    empty_consensus_rows,
    make_af3_ablation_jobs,
    pandora_calibration_runs,
    pandora_candidate_runs,
    primary_ten,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / "processed" / "charge_reversal_binding_pilot_2026-09-07"
DEFAULT_OUT = ROOT / "processed" / "charge_reversal_computational_recovery_2026-09-07"
PANDORA_COMMIT = "618bbb573645bd428697666475a2e47a2e418d01"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "charge-reversal-recovery/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        destination.write_bytes(response.read())


def strip_pmhc_copy(raw: Path, output: Path) -> dict[str, Any]:
    """Keep the first pMHC copy and map A/B/C to PANDORA's M/N/P convention."""
    kept: list[str] = []
    counts = {chain: 0 for chain in "ABC"}
    residue_numbers: dict[tuple[str, str, str], int] = {}
    next_residue = {chain: 0 for chain in "ABC"}
    for line in raw.read_text(encoding="utf-8", errors="replace").splitlines():
        record = line[:6].strip()
        if record in {"HEADER", "TITLE", "COMPND", "SOURCE", "REMARK", "CRYST1"}:
            kept.append(line)
        elif record in {"ATOM", "TER"} and len(line) > 21 and line[21] in counts:
            source_chain = line[21]
            mapped = {"A": "M", "B": "N", "C": "P"}[source_chain]
            residue_key = (source_chain, line[22:26], line[26:27])
            if record != "TER" and residue_key not in residue_numbers:
                next_residue[source_chain] += 1
                residue_numbers[residue_key] = next_residue[source_chain]
            new_number = residue_numbers.get(residue_key, next_residue[source_chain])
            kept.append(line[:21] + mapped + f"{new_number:4d}" + " " + line[27:])
            if record == "ATOM":
                counts[line[21]] += 1
    kept.append("END")
    output.write_text("\n".join(kept) + "\n", encoding="utf-8")
    if any(value == 0 for value in counts.values()):
        raise ValueError(f"missing A/B/C chain atoms in {raw.name}: {counts}")
    return {"chain_atom_counts": counts, "selection": "first A/B/C pMHC copy remapped to M/N/P; TCR and duplicate pMHC copies excluded"}


def baseline_pairwise(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    by_arm = {
        arm: [row for row in rows if str(row["sample_id"]).startswith(arm)]
        for arm in ("BALF5", "TALDO1")
    }
    for arm, arm_rows in by_arm.items():
        for left, right in combinations(arm_rows, 2):
            left_set = set(str(left["peptide_template_pdb_ids"]).split(";"))
            right_set = set(str(right["peptide_template_pdb_ids"]).split(";"))
            union = left_set | right_set
            output.append({
                "arm": arm,
                "left_sample_id": left["sample_id"],
                "right_sample_id": right["sample_id"],
                "template_set_jaccard": round(len(left_set & right_set) / len(union), 3),
                "same_template_set": left_set == right_set,
                "absolute_delta_median_peptide_plddt": round(abs(float(left["median_peptide_mean_plddt"]) - float(right["median_peptide_mean_plddt"])), 3),
                "absolute_delta_within_job_peptide_rmsd_A": round(abs(float(left["within_job_pairwise_peptide_rmsd_median_A"]) - float(right["within_job_pairwise_peptide_rmsd_median_A"])), 3),
                "interpretation": "descriptive association; seed and sequence also differ",
                "claim_boundary": CLAIM_BOUNDARY,
            })
    return output


def build_package(out: Path, baseline: Path) -> None:
    if out.exists():
        raise FileExistsError(f"output already exists: {out}")
    analysis = out / "baseline_diagnosis"
    af3 = out / "alphafold_ablation"
    pandora = out / "pandora"
    templates_raw = pandora / "templates" / "raw"
    templates_pmhc = pandora / "templates" / "pmhc_only"
    environment = pandora / "environment"
    consensus = out / "consensus"
    for directory in (analysis, af3, templates_raw, templates_pmhc, environment, consensus):
        directory.mkdir(parents=True, exist_ok=True)

    baseline_return = baseline / "modeling" / "af3_return_2026-09-07"
    baseline_summary = baseline / "modeling" / "analysis" / "af3_job_summary_10.csv"
    audit = audit_baseline_templates(baseline_return, baseline_summary)
    write_csv(analysis / "baseline_af3_peptide_template_audit_10.csv", audit)
    write_csv(analysis / "baseline_template_outcome_pairwise_20.csv", baseline_pairwise(audit))
    shutil.copy2(
        baseline / "modeling" / "analysis" / "af3_sample_metrics_50.csv",
        analysis / "baseline_model_contacts_and_qc_50.csv",
    )
    shutil.copy2(
        baseline / "modeling" / "analysis" / "wt_mutant_register_comparisons_8.csv",
        analysis / "baseline_register_assessment_8.csv",
    )
    (analysis / "BASELINE_DIAGNOSIS.md").write_text(
        "# Frozen AlphaFold baseline diagnosis\n\n"
        "The 10-job/50-model return is preserved in the original pilot package and linked by checksum here. "
        "All jobs enabled structure templates for the peptide. The audit confirms that variants did not all receive "
        "the same peptide-template set or alignment coverage. BALF5 WT/P7 also had lower peptide confidence or greater "
        "pose variability than the P6 mutants. This is a plausible technical confounder, not proof that templates caused "
        "the structural differences: sequence and seed-sensitive inference remain entangled in the frozen run.\n\n"
        "The matched 2x2 design is therefore required. Global AlphaFold ranking scores are retained only as provenance "
        "and are not biological evidence.\n",
        encoding="utf-8",
    )

    jobs = make_af3_ablation_jobs()
    write_json(af3 / "af3_recovery_batch_30_jobs.json", jobs)
    job_rows = []
    peptide_by_sequence = {row["sequence"]: row["sample_id"] for row in primary_ten()}
    for job in jobs:
        chain = job["sequences"][2]["proteinChain"]
        job_rows.append({
            "job_name": job["name"],
            "sample_id": peptide_by_sequence[chain["sequence"]],
            "peptide_templates_enabled": chain["useStructureTemplate"],
            "hla_templates_enabled": True,
            "seed": job["modelSeeds"][0],
            "requested_models": 5,
            "status": "prepared_not_submitted",
            "returned_archive_path": "",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(af3 / "AF3_RECOVERY_STATUS_30.csv", job_rows)
    factorial = [{
        "condition": "existing_baseline",
        "peptide_templates_enabled": True,
        "seed": 314159,
        "jobs": 10,
        "models": 50,
        "status": "complete_exact_sequence_qc",
    }] + [{
        "condition": label,
        "peptide_templates_enabled": enabled,
        "seed": seed,
        "jobs": 10,
        "models": 50,
        "status": "prepared_not_submitted",
    } for label, enabled, seed in AF3_CONDITIONS]
    write_csv(af3 / "af3_factorial_design_40_jobs.csv", factorial)
    (af3 / "UPLOAD_AND_RETURN_QC.md").write_text(
        "# AlphaFold Server recovery batch\n\n"
        "Upload `af3_recovery_batch_30_jobs.json` as one batch. Before submission, confirm exactly 30 jobs. "
        "The two HLA chains must retain templates; peptide-template state and seed must match the status table. "
        "Each exact 15-mer must return five model CIF, summary-confidence JSON, and full-data JSON files.\n\n"
        "Do not merge results with the frozen baseline until all 150 new models pass exact three-chain sequence, "
        "completeness, clash, and register checks. Partial returns remain `not_evaluable`.\n",
        encoding="utf-8",
    )

    template_records = []
    for pdb_id in ("1BX2", "6CQQ"):
        raw = templates_raw / f"{pdb_id}.pdb"
        download(f"https://files.rcsb.org/download/{pdb_id}.pdb", raw)
        prepared = templates_pmhc / f"{pdb_id}_MNP_pandora.pdb"
        selection = strip_pmhc_copy(raw, prepared)
        template_records.append({
            "pdb_id": pdb_id,
            "raw_rcsb_url": f"https://files.rcsb.org/download/{pdb_id}.pdb",
            "raw_relative_path": str(raw.relative_to(out)),
            "raw_sha256": sha256(raw),
            "prepared_relative_path": str(prepared.relative_to(out)),
            "prepared_sha256": sha256(prepared),
            "chain_selection": selection["selection"],
            "chain_atom_counts": json.dumps(selection["chain_atom_counts"], sort_keys=True),
            "template_role": "primary" if pdb_id == "1BX2" else "sensitivity_check",
        })
    write_csv(pandora / "template_manifest.csv", template_records)
    write_csv(pandora / "candidate_runs_60_requesting_1200_models.csv", pandora_candidate_runs())
    write_csv(pandora / "calibration_runs_6_requesting_120_models.csv", pandora_calibration_runs())
    (pandora / "PANDORA_STATUS.md").write_text(
        "# PANDORA status\n\n"
        "The exact candidate and leave-one-template-out calibration manifests are prepared. The raw RCSB files and "
        "pMHC-only A/B/C copies are checksummed. Execution is gated on a working pinned PANDORA/MODELLER environment.\n\n"
        "Interpretation stays locked until both experimental-template calibration targets recover the expected register "
        "and at least one top-five model is below 2 A core-backbone RMSD. Any target-template leakage fails calibration.\n",
        encoding="utf-8",
    )
    (environment / "environment.yml").write_text(
        "name: charge-reversal-pandora\nchannels:\n  - csb-nijmegen\n  - salilab\n  - bioconda\n  - conda-forge\n"
        "dependencies:\n  - python=3.10.21\n  - modeller=10.8\n  - muscle=5.1.0\n  - blast=2.16.0\n  - biopython=1.88\n  - pdb2sql=0.5.3\n  - pip\n",
        encoding="utf-8",
    )
    write_json(environment / "software_lock.json", {
        "pandora_repository": "https://github.com/X-lab-3D/PANDORA.git",
        "pandora_commit": PANDORA_COMMIT,
        "python": "3.10.21",
        "modeller": "10.8",
        "muscle": "5.1.0 package; executable reports 5.2.osxarm64",
        "blast": "2.16.0",
        "license_key_stored_in_package": False,
        "platform_note": "Upstream recommends Python 3.11.10 and BLAST 2.10. On osx-arm64 MODELLER 10.8 resolves with Python 3.10, and the available native BLAST pin is 2.16.0. Runtime imports and MODELLER license initialization were tested.",
    })
    (environment / "INSTALL.md").write_text(
        "# Pinned environment installation\n\n"
        "Use a separate Conda environment from `environment.yml`. Set `KEY_MODELLER` only in the process environment; "
        "never place the license key in this package. Clone PANDORA, check out the commit in `software_lock.json`, install "
        "it editable, then fetch its template database. Record the resulting explicit package list before modeling.\n",
        encoding="utf-8",
    )

    write_csv(consensus / "consensus_classification_8.csv", empty_consensus_rows())
    write_csv(consensus / "per_model_metrics_template.csv", [{
        "method": "",
        "condition_or_template": "",
        "sample_id": "",
        "register_start_1_based": "",
        "model_id": "",
        "sequence_qc_pass": "",
        "clash_qc_pass": "",
        "core_backbone_rmsd_A": "",
        "contact_present_within_4A": "",
        "exclusion_reason": "",
        "note": "blank schema; models are technical ensemble members, not biological replicates",
    }])
    write_json(consensus / "decision_rules.json", {
        "classifications": [
            "robust_structural_contact_hypothesis", "no_structural_charge_contact_support",
            "method_or_template_dependent", "register_confounded", "not_evaluable",
        ],
        "contact_frequency_min_each_method": 0.8,
        "expected_register_core_rmsd_max_A": 2.0,
        "calibration_required": True,
        "all_ensembles_required": True,
        "p_values_prohibited": True,
        "model_scores_are_not_affinities": True,
        "apbs_ranking_prohibited": True,
        "claim_boundary": CLAIM_BOUNDARY,
    })

    baseline_zip = baseline.with_suffix(".zip")
    sources = []
    for label, path in (("frozen_baseline_directory_checksum_table", baseline / "SHA256SUMS.csv"), ("frozen_baseline_zip", baseline_zip)):
        sources.append({"source": label, "absolute_path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size})
    write_csv(out / "source_manifest.csv", sources)
    software_rows = []
    for path in (
        ROOT / "src" / "charge_reversal_recovery.py",
        ROOT / "src" / "prepare_charge_reversal_recovery.py",
        ROOT / "src" / "run_charge_reversal_pandora.py",
        ROOT / "src" / "analyze_charge_reversal_pandora_calibration.py",
        ROOT / "src" / "analyze_charge_reversal_pandora_candidates.py",
        ROOT / "src" / "validate_charge_reversal_recovery.py",
    ):
        software_rows.append({"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size})
    write_csv(out / "software_manifest.csv", software_rows)
    (out / "LITERATURE_RATIONALE.md").write_text(
        "# Literature rationale\n\n"
        "PANDORA v2.0 benchmarked 136 pMHC-II structures and reported median peptide backbone RMSDs of 0.42 A "
        "for the binding core and 0.88 A for the full peptide. The same study found incorrect automated anchor predictions "
        "in 33/136 structures and shifted AlphaFold cores in two of four direct examples, motivating explicit adjacent-register "
        "enumeration and method comparison (PMCID: PMC10739464; PMID: 38143769).\n\n"
        "Computed mutation structures remain secondary evidence. They cannot turn a modeled K-to-E geometry into measured "
        "binding or direct-contact proof (PMID: 34222328).\n",
        encoding="utf-8",
    )
    (out / "README.md").write_text(
        "# BALF5-TALDO1 computational charge-reversal recovery\n\n"
        "## Current result\n\n"
        "The frozen 10-job/50-model AlphaFold run is diagnosed and preserved. The exact 30-job matched AlphaFold batch, "
        "60 PANDORA candidate runs requesting 1,200 models, and six calibration runs requesting 120 models are prepared. "
        "No consensus contact classification is yet evaluable because the new ensembles have not completed.\n\n"
        "## What is ready\n\n"
        "- `baseline_diagnosis/`: actual peptide template IDs, coverage/mappings, confidence, and pose variability.\n"
        "- `alphafold_ablation/`: one exact external upload batch plus the complete 2x2 design/status.\n"
        "- `pandora/`: checksummed 1BX2/6CQQ inputs, calibration and 1,200-model manifests, pinned environment specification.\n"
        "- `consensus/`: fail-closed decision rules and current `not_evaluable` classifications.\n\n"
        "## Boundary\n\n"
        f"{CLAIM_BOUNDARY}\n",
        encoding="utf-8",
    )
    write_csv(out / "SHA256SUMS.csv", checksum_rows(out))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    args = parser.parse_args()
    build_package(args.output.resolve(), args.baseline.resolve())
    print(json.dumps({"output": str(args.output.resolve()), "status": "prepared", "new_af3_jobs": 30, "pandora_candidate_models_requested": 1200}, indent=2))


if __name__ == "__main__":
    main()
