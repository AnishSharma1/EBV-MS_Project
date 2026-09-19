"""Add an independent MixMHC2pred baseline check to a TALDO1 HLA package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PACKAGE = ROOT / "processed" / "taldo1_hla_next_leg_2026-09-04"
DEFAULT_BINARY = Path.home() / ".cache/ebv_ms_tools/mixmhc2pred_v2.1.beta1.2/MixMHC2pred"
ALLELE_MAP = {
    "HLA-DRB1*13:03": "DRB1_13_03",
    "HLA-DRB1*15:01": "DRB1_15_01",
    "HLA-DRB5*01:01": "DRB5_01_01",
    "HLA-DQA1*01:02/DQB1*06:02": "DQA1_01_02__DQB1_06_02",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_mixmhc(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    args = parser.parse_args()
    package = args.package.resolve()
    baseline = read_csv(package / "baseline_hla_predictions.csv")
    baseline_lookup = {(row["allele"], row["peptide_id"]): row for row in baseline}
    peptides: dict[str, tuple[str, str]] = {}
    for row in baseline:
        peptides.setdefault(row["sequence"], (row["peptide_id"], row["protein"]))

    prepared = package / "prepared_inputs/mixmhc2pred_baseline_no_context.txt"
    prepared.write_text("\n".join(peptides) + "\n", encoding="utf-8")
    raw = package / "raw_responses/mixmhc2pred_2_1_baseline_no_context.tsv"
    command = [str(args.binary), "-i", str(prepared), "-o", str(raw), "-a", *ALLELE_MAP.values(), "--no_context"]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    (package / "raw_responses/mixmhc2pred_2_1_baseline.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (package / "raw_responses/mixmhc2pred_2_1_baseline.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"MixMHC2pred failed with exit {completed.returncode}")

    by_peptide = {row["Peptide"]: row for row in parse_mixmhc(raw)}
    output: list[dict[str, object]] = []
    for sequence, (peptide_id, protein) in peptides.items():
        row = by_peptide[sequence]
        for hla, tool_allele in ALLELE_MAP.items():
            p1 = int(row[f"CoreP1_{tool_allele}"])
            predicted_core = sequence[p1 - 1:p1 + 8]
            iedb_core = baseline_lookup[(hla, peptide_id)]["predicted_core"]
            output.append({
                "allele": hla,
                "peptide_id": peptide_id,
                "protein": protein,
                "sequence": sequence,
                "predicted_core": predicted_core,
                "iedb_predicted_core": iedb_core,
                "register_agreement_with_iedb": predicted_core == iedb_core,
                "rank_percentile": float(row[f"%Rank_{tool_allele}"]),
                "core_start_1_based": p1,
                "predictor": "MixMHC2pred v2.1-beta1",
                "context_mode": "no_context",
                "claim_boundary": "Independent binding-prediction check only; not evidence of presentation, T-cell recognition, cross-reactivity, or MS mechanism.",
            })
    write_csv(package / "mixmhc2pred_baseline_validation.csv", output)

    readme = package / "README.md"
    text = readme.read_text(encoding="utf-8")
    marker = "## How to read the files"
    addition = (
        "## Independent predictor check\n\n"
        "MixMHC2pred v2.1-beta1, run without peptide-flanking context, independently placed both arms of HY15_SEQ_02 at or below rank 20 for DRB1*15:01 (3.15 and 3.49) and DRB5*01:01 (9.26 and 14.3). It also placed both HY13_SEQ_02 arms at or below rank 20 for DRB1*13:03 (0.23 and 12.2). The native HY13/DRB1*13:03 and HY15/DRB1*15:01 cores agree with IEDB. For HY15 on DRB5*01:01, the BALF5 core agrees, but TALDO1 216-230 does not: IEDB predicts YNYYKKFSY while MixMHC2pred predicts VTKIYNYYK. Thus two predictors corroborate baseline binding priority, while the DRB5 TALDO1 register remains unresolved and must be mapped before selecting register-specific charge mutations. These predictions do not validate natural presentation or cross-reactivity.\n\n"
    )
    if "## Independent predictor check" in text:
        before, remainder = text.split("## Independent predictor check", 1)
        _, after = remainder.split(marker, 1)
        text = before + addition + marker + after
    else:
        text = text.replace(marker, addition + marker)
    readme.write_text(text, encoding="utf-8")

    checksum = package / "SHA256SUMS.csv"
    rows = []
    for path in sorted(item for item in package.rglob("*") if item.is_file() and item.name != checksum.name):
        rows.append({"relative_path": str(path.relative_to(package)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
    write_csv(checksum, rows)
    print(f"wrote {len(output)} validation rows and refreshed {len(rows)} checksums")


if __name__ == "__main__":
    main()
