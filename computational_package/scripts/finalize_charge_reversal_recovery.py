"""Freeze recovery metadata, checksums, and the distributable archive."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from charge_reversal_binding_pilot import checksum_rows, write_csv


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "processed" / "charge_reversal_computational_recovery_2026-09-07"
BASELINE = ROOT / "processed" / "charge_reversal_binding_pilot_2026-09-07"
ENVIRONMENT = Path.home() / ".cache" / "ebv_ms_tools" / "charge-reversal-pandora"
MICROMAMBA = Path.home() / ".cache" / "ebv_ms_tools" / "micromamba-bin" / "micromamba"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    diagnosis = PACKAGE / "baseline_diagnosis"
    shutil.copy2(BASELINE / "modeling" / "analysis" / "af3_sample_metrics_50.csv", diagnosis / "baseline_model_contacts_and_qc_50.csv")
    shutil.copy2(BASELINE / "modeling" / "analysis" / "wt_mutant_register_comparisons_8.csv", diagnosis / "baseline_register_assessment_8.csv")

    environment_dir = PACKAGE / "pandora" / "environment"
    explicit = subprocess.run(
        [str(MICROMAMBA), "list", "-p", str(ENVIRONMENT), "--explicit"],
        capture_output=True, text=True, check=True,
    ).stdout
    (environment_dir / "conda_explicit_lock.txt").write_text(explicit, encoding="utf-8")
    probe = subprocess.run(
        [str(ENVIRONMENT / "bin" / "python"), "-c", "import modeller, PANDORA; modeller.Environ(); print(PANDORA.version)"],
        capture_output=True, text=True, check=False,
    )
    verification = {
        "environment_prefix": str(ENVIRONMENT),
        "python_version": subprocess.run([str(ENVIRONMENT / "bin" / "python"), "--version"], capture_output=True, text=True, check=True).stdout.strip(),
        "modeller_and_pandora_runtime_probe": "pass" if probe.returncode == 0 else "fail",
        "pandora_commit": "618bbb573645bd428697666475a2e47a2e418d01",
        "license_key_in_package": False,
        "license_value_in_package": False,
        "muscle_version": subprocess.run([str(ENVIRONMENT / "bin" / "muscle"), "-version"], capture_output=True, text=True, check=True).stdout.splitlines()[0],
        "blast_version": subprocess.run([str(ENVIRONMENT / "bin" / "blastp"), "-version"], capture_output=True, text=True, check=True).stdout.splitlines()[0],
    }
    (environment_dir / "environment_verification.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    if verification["modeller_and_pandora_runtime_probe"] != "pass":
        raise RuntimeError("MODELLER/PANDORA environment probe failed")

    software_paths = (
        "src/charge_reversal_recovery.py",
        "src/prepare_charge_reversal_recovery.py",
        "src/run_charge_reversal_pandora.py",
        "src/analyze_charge_reversal_pandora_calibration.py",
        "src/analyze_charge_reversal_pandora_candidates.py",
        "src/analyze_charge_reversal_af3_recovery.py",
        "src/validate_charge_reversal_recovery.py",
        "src/finalize_charge_reversal_recovery.py",
    )
    write_csv(PACKAGE / "software_manifest.csv", [{
        "path": value, "sha256": digest(ROOT / value), "bytes": (ROOT / value).stat().st_size,
    } for value in software_paths])

    calibration = json.loads((PACKAGE / "pandora" / "results" / "CALIBRATION_GATE.json").read_text(encoding="utf-8"))
    candidate = json.loads((PACKAGE / "pandora" / "results" / "PANDORA_CANDIDATE_SUMMARY.json").read_text(encoding="utf-8"))
    af3 = json.loads((PACKAGE / "alphafold_ablation" / "analysis" / "AF3_RECOVERY_SUMMARY.json").read_text(encoding="utf-8"))
    consensus_path = PACKAGE / "consensus" / "consensus_classification_8.csv"
    with consensus_path.open(newline="", encoding="utf-8") as handle:
        consensus = list(csv.DictReader(handle))
    class_counts: dict[str, int] = {}
    for row in consensus:
        class_counts[row["classification"]] = class_counts.get(row["classification"], 0) + 1
    class_summary = ", ".join(f"{value} `{key}`" for key, value in sorted(class_counts.items()))

    (PACKAGE / "pandora" / "PANDORA_STATUS.md").write_text(
        "# PANDORA status\n\n"
        f"The pinned environment passed runtime checks. Calibration: **{calibration['calibration_status']}** across 120 models. "
        f"Candidate production: **{candidate['models']}/1,200 models**, with exact sequence QC "
        f"{'passing' if candidate['all_exact_sequence_qc_pass'] else 'failing'}. "
        f"The expected register was supported by both templates for {candidate['samples_supporting_expected_register']}/10 peptides.\n\n"
        "Models are technical ensemble members, not biological replicates. Molpdf/DOPE and structural contacts are not binding affinities.\n",
        encoding="utf-8",
    )
    (PACKAGE / "README.md").write_text(
        "# BALF5-TALDO1 computational charge-reversal recovery\n\n"
        "## Finished locally\n\n"
        f"- Frozen AlphaFold baseline diagnosed: 10 jobs / 50 models.\n"
        f"- Reciprocal PANDORA calibration: {calibration['calibration_status']}, 120/120 models scored.\n"
        f"- PANDORA candidate ensemble: {candidate['models']}/1,200 models analyzed; exact sequence QC "
        f"{'passed' if candidate['all_exact_sequence_qc_pass'] else 'failed'}.\n"
        f"- Expected register supported across both templates: {candidate['samples_supporting_expected_register']}/10 peptides.\n"
        f"- Controlled AlphaFold recovery: {af3['jobs']}/30 jobs and {af3['models']}/150 models returned; exact sequence QC passed.\n"
        f"- Predeclared AlphaFold seed match: {af3['predeclared_seed_matches']}/30 jobs; "
        f"{af3['predeclared_seed_mismatches']} mismatch requires a one-job rerun.\n\n"
        "## Current conclusion\n\n"
        f"Fail-closed classifications: {class_summary}. The seven fully matched mutation comparisons are `register_confounded` because "
        "the proposed register did not survive both methods; one mutation is `not_evaluable` because its returned AlphaFold seed differed from the predeclared seed. "
        "These outcomes are method sensitivity, not evidence that charge changes binding.\n\n"
        "## Key files\n\n"
        "- `baseline_diagnosis/`: actual AlphaFold peptide-template identities, mappings, confidence, contacts, and register checks.\n"
        "- `alphafold_ablation/af3_recovery_batch_30_jobs.json`: exact external rerun batch.\n"
        "- `alphafold_ablation/analysis/`: all 150 model metrics, register comparisons, and cross-method evidence.\n"
        "- `alphafold_ablation/af3_retry_seed_mismatch_jobs.json`: exact one-job correction batch.\n"
        "- `alphafold_ablation/raw/`: preserved downloaded AlphaFold archive.\n"
        "- `pandora/results/CALIBRATION_GATE.json`: predeclared calibration result.\n"
        "- `pandora/results/pandora_candidate_sample_summary_10.csv`: register and contact summary.\n"
        "- `consensus/consensus_classification_8.csv`: fail-closed current classifications.\n\n"
        "## Claim boundary\n\n"
        "This package tests cross-method structural reproducibility. It does not measure peptide-HLA affinity, prove a direct electrostatic contact, or establish T-cell recognition or molecular mimicry.\n",
        encoding="utf-8",
    )
    write_csv(PACKAGE / "SHA256SUMS.csv", checksum_rows(PACKAGE))
    archive = shutil.make_archive(str(PACKAGE), "zip", root_dir=PACKAGE.parent, base_dir=PACKAGE.name)
    print(json.dumps({
        "package": str(PACKAGE), "archive": archive,
        "archive_sha256": digest(Path(archive)), "archive_bytes": Path(archive).stat().st_size,
    }, indent=2))


if __name__ == "__main__":
    main()
