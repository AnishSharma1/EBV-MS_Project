#!/usr/bin/env python3
"""Run NetMHCIIpan 4.3 EL on the sequence-prefiltered proteome screen."""

from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "output/revision_round_3/proteome_sequence_null_2026-09-15/binding_prediction_manifest_top_1000_per_lead.csv"
MIX_RESULTS = ROOT / "output/revision_round_3/proteome_binding_filter_2026-09-15/mixmhc2pred_all_predictions.csv"
OUT = ROOT / "output/revision_round_3/proteome_binding_consensus_2026-09-15"
ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
METHOD = "netmhciipan_el-4.3"
CONTROLS = {
    "HY13_SEQ_02": [("viral_control", "TGGVYHFVKKHVHES"), ("taldo1_control", "DARLSFDKDAMVARA")],
    "HY15_SEQ_02": [("viral_control", "TGGVYHFVKKHVHES"), ("taldo1_control", "SVTKIYNYYKKFSYK")],
}
BOUNDARY = (
    "Two-stage top-1,000 sequence-prefilter screen using NetMHCIIpan 4.3 EL and local "
    "MixMHC2pred 2.1-beta1. Not a whole-proteome binding scan or evidence of presentation, "
    "TCR recognition, cross-reactivity, or disease mechanism."
)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_batch(records, allele, raw_path):
    if raw_path.exists():
        data = raw_path.read_bytes()
    else:
        fasta = "".join(f">{i}\n{row['peptide']}\n" for i, row in enumerate(records, 1))
        error = None
        for attempt in range(3):
            try:
                response = requests.post(
                    ENDPOINT,
                    data={"method": METHOD, "sequence_text": fasta, "allele": allele, "length": "asis"},
                    timeout=240,
                )
                response.raise_for_status()
                data = response.content
                if not data.strip():
                    raise RuntimeError("empty response")
                raw_path.write_bytes(data)
                break
            except Exception as exc:
                error = exc
                time.sleep(2 * (attempt + 1))
        else:
            raise RuntimeError(f"IEDB batch failed after 3 attempts: {error}")
    text = data.decode("utf-8", "replace")
    table = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    by_seq = {}
    for row in table:
        seq_num = int(row["seq_num"])
        current = by_seq.get(seq_num)
        if current is None or float(row["rank"]) < float(current["rank"]):
            by_seq[seq_num] = row
    if len(by_seq) != len(records):
        raise ValueError(f"expected {len(records)} predictions, found {len(by_seq)} in {raw_path.name}")
    return [by_seq[i] for i in range(1, len(records) + 1)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw_responses"
    raw_dir.mkdir(exist_ok=True)
    manifest = read_csv(MANIFEST)
    mix = {row["submission_id"]: row for row in read_csv(MIX_RESULTS)}
    predictions = []
    for target_id in sorted(CONTROLS):
        source = [row for row in manifest if row["target_id"] == target_id]
        allele = source[0]["allele"]
        records = [
            {"submission_id": row["submission_id"], "record_type": "human_prefilter", "peptide": row["parent_15mer"], "sequence_prefilter_rank": row["sequence_prefilter_rank"], "accession": row["accession"]}
            for row in source
        ]
        records.extend(
            {"submission_id": f"{target_id}_{kind.upper()}", "record_type": kind, "peptide": peptide, "sequence_prefilter_rank": "", "accession": ""}
            for kind, peptide in CONTROLS[target_id]
        )
        for batch_index in range(0, len(records), 100):
            batch = records[batch_index:batch_index + 100]
            raw_path = raw_dir / f"{target_id.lower()}_{batch_index // 100 + 1:02d}.tsv"
            returned = fetch_batch(batch, allele, raw_path)
            for source_row, pred in zip(batch, returned):
                predictions.append({
                    "submission_id": source_row["submission_id"],
                    "target_id": target_id,
                    "record_type": source_row["record_type"],
                    "allele": allele,
                    "sequence_prefilter_rank": source_row["sequence_prefilter_rank"],
                    "accession": source_row["accession"],
                    "peptide": source_row["peptide"],
                    "netmhciipan_el_core": pred["core_peptide"],
                    "netmhciipan_el_rank": pred["rank"],
                    "raw_response": str(raw_path.relative_to(OUT)),
                    "claim_boundary": BOUNDARY,
                })

    combined = []
    for row in predictions:
        mix_row = mix[row["submission_id"]]
        net_rank = float(row["netmhciipan_el_rank"])
        mix_rank = float(mix_row["mixmhc2pred_rank_percentile"])
        combined.append({
            **row,
            "mixmhc2pred_core": mix_row["predicted_core"],
            "mixmhc2pred_rank": mix_row["mixmhc2pred_rank_percentile"],
            "core_agreement": row["netmhciipan_el_core"] == mix_row["predicted_core"],
            "both_pass_rank_2": net_rank <= 2 and mix_rank <= 2,
            "both_pass_rank_5": net_rank <= 5 and mix_rank <= 5,
            "both_pass_rank_10": net_rank <= 10 and mix_rank <= 10,
            "both_pass_rank_20": net_rank <= 20 and mix_rank <= 20,
        })
    write_csv(OUT / "two_predictor_predictions.csv", combined)

    summaries = []
    for target_id in sorted(CONTROLS):
        target = [row for row in combined if row["target_id"] == target_id]
        human = [row for row in target if row["record_type"] == "human_prefilter"]
        controls = {row["record_type"]: row for row in target if row["record_type"] != "human_prefilter"}
        summaries.append({
            "target_id": target_id,
            "allele": target[0]["allele"],
            "human_prefilter_count": len(human),
            "both_predictors_pass_rank_2": sum(row["both_pass_rank_2"] for row in human),
            "both_predictors_pass_rank_5": sum(row["both_pass_rank_5"] for row in human),
            "both_predictors_pass_rank_10": sum(row["both_pass_rank_10"] for row in human),
            "both_predictors_pass_rank_20": sum(row["both_pass_rank_20"] for row in human),
            "predictor_core_agreement_count": sum(row["core_agreement"] for row in human),
            "viral_netmhciipan_rank": controls["viral_control"]["netmhciipan_el_rank"],
            "viral_mixmhc2pred_rank": controls["viral_control"]["mixmhc2pred_rank"],
            "taldo1_netmhciipan_rank": controls["taldo1_control"]["netmhciipan_el_rank"],
            "taldo1_mixmhc2pred_rank": controls["taldo1_control"]["mixmhc2pred_rank"],
            "claim_boundary": BOUNDARY,
        })
    write_csv(OUT / "two_predictor_consensus_summary.csv", summaries)
    lock = {
        "status": "complete_two_predictor_top_1000_screen",
        "manifest_sha256": digest(MANIFEST),
        "mixmhc2pred_results_sha256": digest(MIX_RESULTS),
        "netmhciipan_method": METHOD,
        "iedb_endpoint": ENDPOINT,
        "response_file_count": len(list(raw_dir.glob("*.tsv"))),
        "claim_boundary": BOUNDARY,
    }
    (OUT / "analysis_lock.json").write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    (OUT / "README.md").write_text(
        "# Two-predictor sequence-plus-binding proteome screen\n\n"
        "NetMHCIIpan 4.3 EL and MixMHC2pred 2.1-beta1 were applied independently to the "
        "top 1,000 sequence-similar human windows for each lead. The summary uses a strict "
        "both-predictors rule at 2%, 5%, 10%, and 20% rank cutoffs. This remains a top-window "
        "sensitivity analysis rather than a complete proteome-wide binding prediction.\n",
        encoding="utf-8",
    )
    checksum_path = OUT / "SHA256SUMS.csv"
    checksum_rows = []
    for path in sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_path):
        checksum_rows.append({"relative_path": str(path.relative_to(OUT)), "sha256": digest(path), "bytes": path.stat().st_size})
    write_csv(checksum_path, checksum_rows)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
