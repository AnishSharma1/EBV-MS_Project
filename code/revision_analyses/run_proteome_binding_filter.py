#!/usr/bin/env python3
"""Apply local MixMHC2pred to the sequence-prefiltered proteome null.

This is a two-stage screen: every reviewed-human 9-mer was scored for sequence
similarity first, then the top 1,000 per lead were evaluated for predicted HLA-II
binding. It is not a whole-proteome binding scan and cannot establish presentation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path


ALLELE_TOOL = {
    "HLA-DRB1*13:03": "DRB1_13_03",
    "HLA-DRB1*15:01": "DRB1_15_01",
}
CONTROLS = {
    "HY13_SEQ_02": {
        "allele": "HLA-DRB1*13:03",
        "viral": "TGGVYHFVKKHVHES",
        "taldo1": "DARLSFDKDAMVARA",
    },
    "HY15_SEQ_02": {
        "allele": "HLA-DRB1*15:01",
        "viral": "TGGVYHFVKKHVHES",
        "taldo1": "SVTKIYNYYKKFSYK",
    },
}
BOUNDARY = (
    "Two-stage screen of the top 1,000 sequence-similar human windows per lead using "
    "local MixMHC2pred 2.1-beta1 without flanking context. This is not a whole-proteome "
    "binding scan or evidence of natural presentation, TCR recognition, cross-reactivity, "
    "or disease mechanism."
)


def parse_args():
    workspace = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=workspace / "output/revision_round_3/proteome_sequence_null_2026-09-15/binding_prediction_manifest_top_1000_per_lead.csv",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path.home() / ".cache/ebv_ms_tools/mixmhc2pred_v2.1.beta1.2/MixMHC2pred",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_3/proteome_binding_filter_2026-09-15",
    )
    return parser.parse_args()


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_mix(path: Path):
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def main():
    args = parse_args()
    out = args.output_dir
    raw = out / "raw_responses"
    prepared = out / "prepared_inputs"
    raw.mkdir(parents=True, exist_ok=True)
    prepared.mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.manifest)
    if len(rows) != 2000:
        raise ValueError(f"expected 2,000 manifest rows, found {len(rows)}")
    if not args.binary.exists():
        raise FileNotFoundError(args.binary)

    results = []
    summaries = []
    for target_id, control in CONTROLS.items():
        allele = control["allele"]
        tool_allele = ALLELE_TOOL[allele]
        subset = [row for row in rows if row["target_id"] == target_id]
        inputs = [(row["submission_id"], "human_prefilter", row["parent_15mer"], row) for row in subset]
        inputs += [
            (f"{target_id}_VIRAL_CONTROL", "viral_control", control["viral"], {}),
            (f"{target_id}_TALDO1_CONTROL", "taldo1_control", control["taldo1"], {}),
        ]
        input_path = prepared / f"{target_id.lower()}_no_context.txt"
        input_path.write_text("\n".join(sequence for _, _, sequence, _ in inputs) + "\n", encoding="utf-8")
        output_path = raw / f"{target_id.lower()}_mixmhc2pred_no_context.tsv"
        completed = subprocess.run(
            [str(args.binary), "-i", str(input_path), "-o", str(output_path), "-a", tool_allele, "--no_context"],
            text=True,
            capture_output=True,
            check=False,
        )
        (raw / f"{target_id.lower()}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (raw / f"{target_id.lower()}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"MixMHC2pred failed for {target_id}: {completed.stderr.strip()}")
        predicted = parse_mix(output_path)
        if len(predicted) != len(inputs):
            raise ValueError(f"{target_id}: expected {len(inputs)} predictions, got {len(predicted)}")

        target_results = []
        for (submission_id, record_type, sequence, source), pred in zip(inputs, predicted):
            if pred.get("Peptide") != sequence:
                raise ValueError(f"sequence mismatch for {submission_id}")
            raw_p1 = pred[f"CoreP1_{tool_allele}"]
            raw_rank = pred[f"%Rank_{tool_allele}"]
            if raw_p1 in {"", "NA"} or raw_rank in {"", "NA"}:
                raise ValueError(f"missing prediction for {submission_id}: P1={raw_p1}, rank={raw_rank}")
            p1 = int(float(raw_p1))
            rank = float(raw_rank)
            result = {
                "submission_id": submission_id,
                "target_id": target_id,
                "record_type": record_type,
                "allele": allele,
                "sequence_prefilter_rank": source.get("sequence_prefilter_rank", ""),
                "accession": source.get("accession", ""),
                "peptide": sequence,
                "predicted_core": sequence[p1 - 1:p1 + 8],
                "predicted_core_start_1_based": p1,
                "mixmhc2pred_rank_percentile": f"{rank:.12g}",
                "passes_rank_2": rank <= 2,
                "passes_rank_5": rank <= 5,
                "passes_rank_10": rank <= 10,
                "passes_rank_20": rank <= 20,
                "predictor": "MixMHC2pred 2.1-beta1 no_context",
                "claim_boundary": BOUNDARY,
            }
            results.append(result)
            target_results.append(result)

        human = [row for row in target_results if row["record_type"] == "human_prefilter"]
        controls = {row["record_type"]: row for row in target_results if row["record_type"] != "human_prefilter"}
        summary = {
            "target_id": target_id,
            "allele": allele,
            "human_sequence_prefilter_count": len(human),
            "human_pass_rank_2": sum(row["passes_rank_2"] for row in human),
            "human_pass_rank_5": sum(row["passes_rank_5"] for row in human),
            "human_pass_rank_10": sum(row["passes_rank_10"] for row in human),
            "human_pass_rank_20": sum(row["passes_rank_20"] for row in human),
            "viral_control_rank": controls["viral_control"]["mixmhc2pred_rank_percentile"],
            "viral_control_core": controls["viral_control"]["predicted_core"],
            "taldo1_control_rank": controls["taldo1_control"]["mixmhc2pred_rank_percentile"],
            "taldo1_control_core": controls["taldo1_control"]["predicted_core"],
            "human_with_rank_le_taldo1": sum(
                float(row["mixmhc2pred_rank_percentile"]) <= float(controls["taldo1_control"]["mixmhc2pred_rank_percentile"])
                for row in human
            ),
            "claim_boundary": BOUNDARY,
        }
        summaries.append(summary)

    write_csv(out / "mixmhc2pred_all_predictions.csv", results)
    write_csv(out / "two_stage_binding_filter_summary.csv", summaries)
    lock = {
        "status": "complete_local_two_stage_screen",
        "manifest": str(args.manifest),
        "manifest_sha256": digest(args.manifest),
        "binary": str(args.binary),
        "binary_sha256": digest(args.binary),
        "predictor": "MixMHC2pred 2.1-beta1 no_context",
        "human_rows": 2000,
        "control_rows": 4,
        "claim_boundary": BOUNDARY,
    }
    (out / "analysis_lock.json").write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    readme = (
        "# Proteome sequence-plus-binding screen\n\n"
        "This package applies local MixMHC2pred 2.1-beta1 to the 1,000 highest-scoring "
        "human sequence windows for each BALF5 lead. Viral and TALDO1 peptides are explicit "
        "controls. It is a two-stage sensitivity analysis, not a whole-proteome HLA-binding "
        "scan and not immunopeptidomics evidence.\n\n"
        "The summary reports how many prefiltered human windows pass 2%, 5%, 10%, and 20% "
        "predicted-rank thresholds and how TALDO1 compares within that deliberately enriched set.\n"
    )
    (out / "README.md").write_text(readme, encoding="utf-8")
    checksum_rows = []
    checksum_path = out / "SHA256SUMS.csv"
    for path in sorted(p for p in out.rglob("*") if p.is_file() and p != checksum_path):
        checksum_rows.append({"relative_path": str(path.relative_to(out)), "sha256": digest(path), "bytes": path.stat().st_size})
    write_csv(checksum_path, checksum_rows)
    counts = Counter(row["record_type"] for row in results)
    print(json.dumps({"rows": len(results), "counts": counts, "summaries": summaries}, indent=2, default=dict))


if __name__ == "__main__":
    main()
