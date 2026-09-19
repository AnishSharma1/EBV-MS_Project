"""Build a focused alternate-register package for two BALF5-TALDO1 leads and DRB5."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from focused_taldo1_register_analysis import classify_register_evidence, enumerate_pair_windows


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "processed/taldo1_focused_register_2026-09-05"
ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
MIX_BINARY = Path.home() / ".cache/ebv_ms_tools/mixmhc2pred_v2.1.beta1.2/MixMHC2pred"
CLAIM = (
    "Computational peptide-HLA register sensitivity only; not experimental register resolution, "
    "natural presentation, TCR recognition, cross-reactivity, molecular mimicry, or MS mechanism."
)
TARGETS = [
    {
        "target_id": "HY13_SEQ_02",
        "source_target_id": "HY13_SEQ_02",
        "allele": "HLA-DRB1*13:03",
        "mix_allele": "DRB1_13_03",
        "ebv_sequence": "TGGVYHFVKKHVHES",
        "self_sequence": "DARLSFDKDAMVARA",
    },
    {
        "target_id": "HY15_SEQ_02",
        "source_target_id": "HY15_SEQ_02",
        "allele": "HLA-DRB1*15:01",
        "mix_allele": "DRB1_15_01",
        "ebv_sequence": "TGGVYHFVKKHVHES",
        "self_sequence": "SVTKIYNYYKKFSYK",
    },
    {
        "target_id": "HY15_SEQ_02_DRB5",
        "source_target_id": "HY15_SEQ_02",
        "allele": "HLA-DRB5*01:01",
        "mix_allele": "DRB5_01_01",
        "ebv_sequence": "TGGVYHFVKKHVHES",
        "self_sequence": "SVTKIYNYYKKFSYK",
    },
]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def post_iedb(allele: str, method: str, sequences: list[tuple[str, str]]) -> str:
    fasta = "\n".join(f">{name}\n{sequence}" for name, sequence in sequences)
    body = urllib.parse.urlencode({"method": method, "sequence_text": fasta, "allele": allele, "length": "15"}).encode()
    request = urllib.request.Request(ENDPOINT, data=body, headers={"Content-Type": "application/x-www-form-urlencoded"}, method="POST")
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8")


def parse_iedb(target: dict[str, str], method: str, raw: str, raw_path: Path, out: Path) -> list[dict[str, Any]]:
    table = list(csv.DictReader(io.StringIO(raw), delimiter="\t"))
    by_seq = {int(row["seq_num"]): row for row in table}
    output = []
    family = "IEDB recommended" if method == "recommended_binding" else "NetMHCIIpan 4.3"
    for seq_num, arm in enumerate(("ebv", "self"), start=1):
        row = by_seq[seq_num]
        sequence = target[f"{arm}_sequence"]
        if row["peptide"] != sequence or row["allele"] != target["allele"]:
            raise ValueError(f"IEDB response mismatch: {target['target_id']} {method} {arm}")
        core = row["core_peptide"]
        output.append({
            "target_id": target["target_id"], "source_target_id": target["source_target_id"],
            "allele": target["allele"], "arm": arm, "sequence": sequence,
            "predictor_family": family, "method": method, "core": core,
            "core_start_1_based": sequence.index(core) + 1,
            "rank_percentile": row["rank"],
            "score_type": "ic50_nM" if "ic50" in row else "el_score",
            "score_value": row.get("ic50", row.get("score", "")),
            "raw_response_path": str(raw_path.relative_to(out)), "claim_boundary": CLAIM,
        })
    return output


def parse_mix(target: dict[str, str], raw_path: Path, out: Path) -> list[dict[str, Any]]:
    with raw_path.open(newline="", encoding="utf-8") as handle:
        table = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if len(table) != 2:
        raise ValueError("MixMHC2pred must return exactly two peptide rows")
    output = []
    key = target["mix_allele"]
    for arm, row in zip(("ebv", "self"), table):
        sequence = target[f"{arm}_sequence"]
        if row["Peptide"] != sequence:
            raise ValueError(f"MixMHC2pred response mismatch: {target['target_id']} {arm}")
        start = int(row[f"CoreP1_{key}"])
        output.append({
            "target_id": target["target_id"], "source_target_id": target["source_target_id"],
            "allele": target["allele"], "arm": arm, "sequence": sequence,
            "predictor_family": "MixMHC2pred 2.1", "method": "MixMHC2pred 2.1-beta1 no_context",
            "core": sequence[start - 1:start + 8], "core_start_1_based": start,
            "rank_percentile": row[f"%Rank_{key}"], "score_type": "percentile_rank",
            "score_value": row[f"%Rank_{key}"], "raw_response_path": str(raw_path.relative_to(out)),
            "claim_boundary": CLAIM,
        })
    return output


def nested_11mer(sequence: str, core_start: int) -> tuple[int, str]:
    start = max(1, min(core_start - 1, len(sequence) - 10))
    return start, sequence[start - 1:start + 10]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--reuse-raw-from", type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(out)
    (out / "raw_responses").mkdir(parents=True)
    (out / "prepared_inputs").mkdir()
    built_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    cached_manifest: dict[tuple[str, str], dict[str, str]] = {}
    if args.reuse_raw_from:
        cached_manifest_path = args.reuse_raw_from / "raw_response_manifest.csv"
        if cached_manifest_path.is_file():
            with cached_manifest_path.open(newline="", encoding="utf-8") as handle:
                cached_manifest = {(row["target_id"], row["method"]): row for row in csv.DictReader(handle)}

    predictor_rows: list[dict[str, Any]] = []
    raw_manifest: list[dict[str, Any]] = []
    for target in TARGETS:
        sequences = [("ebv", target["ebv_sequence"]), ("self", target["self_sequence"])]
        slug = target["target_id"].lower()
        for method in ("recommended_binding", "netmhciipan_el-4.3", "netmhciipan_ba-4.3"):
            raw_path = out / "raw_responses" / f"{slug}__{method.replace('.', '_')}.tsv"
            cached = args.reuse_raw_from / "raw_responses" / raw_path.name if args.reuse_raw_from else None
            if cached and cached.is_file():
                shutil.copy2(cached, raw_path)
                raw = raw_path.read_text(encoding="utf-8")
                retrieved_utc = cached_manifest.get((target["target_id"], method), {}).get("retrieved_utc", "unknown_cached_retrieval_time")
            else:
                raw = post_iedb(target["allele"], method, sequences)
                raw_path.write_text(raw, encoding="utf-8")
                retrieved_utc = built_utc
            predictor_rows.extend(parse_iedb(target, method, raw, raw_path, out))
            raw_manifest.append({"target_id": target["target_id"], "method": method, "raw_response_path": str(raw_path.relative_to(out)), "sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(), "retrieved_utc": retrieved_utc})

        mix_input = out / "prepared_inputs" / f"{slug}__mixmhc2pred.txt"
        mix_input.write_text("\n".join(sequence for _, sequence in sequences) + "\n", encoding="utf-8")
        mix_raw = out / "raw_responses" / f"{slug}__mixmhc2pred_no_context.tsv"
        command = [str(MIX_BINARY), "-i", str(mix_input), "-o", str(mix_raw), "-a", target["mix_allele"], "--no_context"]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        (out / "raw_responses" / f"{slug}__mixmhc2pred.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (out / "raw_responses" / f"{slug}__mixmhc2pred.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"MixMHC2pred failed for {target['target_id']}")
        predictor_rows.extend(parse_mix(target, mix_raw, out))
        raw_manifest.append({"target_id": target["target_id"], "method": "MixMHC2pred 2.1-beta1 no_context", "raw_response_path": str(mix_raw.relative_to(out)), "sha256": hashlib.sha256(mix_raw.read_bytes()).hexdigest(), "retrieved_utc": built_utc})

    predictor_rows.sort(key=lambda row: (row["target_id"], row["arm"], row["predictor_family"], row["method"]))
    write_csv(out / "predictor_register_records.csv", predictor_rows)
    write_csv(out / "raw_response_manifest.csv", raw_manifest)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in predictor_rows:
        grouped[(row["target_id"], row["arm"])].append(row)
    evidence_rows = []
    for target in TARGETS:
        for arm in ("ebv", "self"):
            rows = grouped[(target["target_id"], arm)]
            result = classify_register_evidence(rows)
            cores = Counter(row["core"] for row in rows if row["method"] in {"recommended_binding", "netmhciipan_el-4.3", "MixMHC2pred 2.1-beta1 no_context"})
            evidence_rows.append({
                "target_id": target["target_id"], "source_target_id": target["source_target_id"],
                "allele": target["allele"], "arm": arm, "sequence": target[f"{arm}_sequence"],
                **result, "core_family_votes": ";".join(f"{core}:{count}" for core, count in sorted(cores.items())),
                "claim_boundary": CLAIM,
            })
    write_csv(out / "register_evidence_summary.csv", evidence_rows)

    window_rows = []
    target_summaries = []
    old_summary_path = ROOT / "processed/high_yield_register_resolution_2026-08-28/target_register_summary.csv"
    with old_summary_path.open(newline="", encoding="utf-8") as handle:
        old_summary = {row["target_id"]: row for row in csv.DictReader(handle)}
    for target in TARGETS:
        arm_rows = {arm: grouped[(target["target_id"], arm)] for arm in ("ebv", "self")}
        arm_evidence = {row["arm"]: row for row in evidence_rows if row["target_id"] == target["target_id"]}
        family_cores = {}
        for arm in ("ebv", "self"):
            for row in arm_rows[arm]:
                family_cores.setdefault((arm, row["predictor_family"]), row["core"])
        rows = enumerate_pair_windows(target["ebv_sequence"], target["self_sequence"])
        for row in rows:
            families = sorted(family for family in {v[1] for v in family_cores} if family_cores.get(("ebv", family)) == row["left_core"] and family_cores.get(("self", family)) == row["right_core"])
            row.update({
                "target_id": target["target_id"], "allele": target["allele"],
                "paired_predictor_family_count": len(families),
                "paired_predictor_families": ";".join(families),
                "is_majority_core_pair": row["left_core"] == arm_evidence["ebv"]["majority_core"] and row["right_core"] == arm_evidence["self"]["majority_core"],
                "claim_boundary": CLAIM,
            })
        rows.sort(key=lambda row: (-float(row["tcr_facing_blosum62_similarity"]), float(row["tcr_face_physicochemical_mismatch"]), -float(row["tcr_facing_sequence_identity"]), int(row["left_start_1_based"]), int(row["right_start_1_based"])))
        for rank, row in enumerate(rows, start=1):
            row["within_target_sequence_similarity_rank"] = rank
        window_rows.extend(rows)
        majority_row = next(row for row in rows if row["is_majority_core_pair"])
        supported_rows = sorted(
            (row for row in rows if int(row["paired_predictor_family_count"]) > 0),
            key=lambda row: int(row["within_target_sequence_similarity_rank"]),
        )
        prior = old_summary.get(target["source_target_id"]) if target["target_id"] == target["source_target_id"] else None
        target_summaries.append({
            "target_id": target["target_id"], "source_target_id": target["source_target_id"], "allele": target["allele"],
            "ebv_register_status": arm_evidence["ebv"]["status"], "ebv_majority_core": arm_evidence["ebv"]["majority_core"],
            "self_register_status": arm_evidence["self"]["status"], "self_majority_core": arm_evidence["self"]["majority_core"],
            "pair_register_status": "predictor_register_consensus" if arm_evidence["ebv"]["status"] == arm_evidence["self"]["status"] == "predictor_register_consensus" else "predictor_register_disagreement",
            "majority_core_pair_sequence_similarity_rank_of_49": majority_row["within_target_sequence_similarity_rank"],
            "predictor_supported_core_pair_ranks": ";".join(
                f"{row['left_core']}|{row['right_core']}|rank={row['within_target_sequence_similarity_rank']}|families={row['paired_predictor_families']}"
                for row in supported_rows
            ),
            "prior_n3_register_sensitivity_status": prior["register_resolution_status"] if prior else "not_available_for_DRB5",
            "experimentally_resolved": False, "recommended_action": "nested_peptide_binding_and_register_mapping_before_charge_mutants",
            "claim_boundary": CLAIM,
        })
    write_csv(out / "all_window_sequence_sensitivity.csv", window_rows)
    write_csv(out / "target_register_summary.csv", target_summaries)

    panel_rows = []
    for target in TARGETS:
        for arm in ("ebv", "self"):
            sequence = target[f"{arm}_sequence"]
            panel_rows.append({"target_id": target["target_id"], "allele": target["allele"], "arm": arm, "reagent_role": "parent_reference", "source_core": "", "sequence": sequence, "length": len(sequence), "proposed_not_ordered": True, "claim_boundary": CLAIM})
            unique = {}
            for row in grouped[(target["target_id"], arm)]:
                unique.setdefault(row["core"], set()).add(row["predictor_family"])
            for core, families in sorted(unique.items()):
                core_start = sequence.index(core) + 1
                nested_start, nested = nested_11mer(sequence, core_start)
                panel_rows.append({"target_id": target["target_id"], "allele": target["allele"], "arm": arm, "reagent_role": "predicted_register_nested_11mer", "source_core": core, "sequence": nested, "length": len(nested), "proposed_not_ordered": True, "claim_boundary": CLAIM})
    write_csv(out / "proposed_nested_peptide_panel.csv", panel_rows)

    (out / "protocol_lock.json").write_text(json.dumps({
        "package_id": "taldo1_focused_register_2026-09-05", "targets": TARGETS,
        "predictor_families": ["IEDB recommended", "NetMHCIIpan 4.3", "MixMHC2pred 2.1"],
        "all_window_rule": "enumerate all fully contained 9mer pairs; sequence metrics only",
        "register_consensus_rule": "all predictor families must select the same 9mer core",
        "nested_peptides_ordered": False, "experimentally_resolved": False, "claim_boundary": CLAIM,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Focused BALF5-TALDO1 register analysis", "",
        f"Package built {built_utc}. Exact retrieval times for each live IEDB response and the local MixMHC2pred runs are recorded in `raw_response_manifest.csv`. The two native-HLA leads and the DRB5 restriction test were analyzed separately.", "",
        "| Target | HLA | BALF5 register | TALDO1 register | Similarity rank | Result |",
        "|---|---|---|---|---:|---|",
    ]
    for row in target_summaries:
        lines.append(f"| {row['target_id']} | {row['allele']} | `{row['ebv_majority_core']}` | `{row['self_majority_core']}` | {row['majority_core_pair_sequence_similarity_rank_of_49']}/49 | {row['pair_register_status']} |")
    lines += [
        "", "HY13_SEQ_02 and HY15_SEQ_02 have full predictor-family register agreement on both arms. Their earlier exhaustive window test nevertheless classified both as `declared_window_only`, so agreement among predictors does not prove the physical register.", "",
        "For HY15_SEQ_02, an unsupported shifted pair (`GVYHFVKKH` / `KIYNYYKKF`) has the highest raw sequence-similarity rank, while the core pair chosen by every predictor family ranks 2/49. This is why sequence similarity alone must not select an HLA register.", "",
        "For HY15_SEQ_02_DRB5, IEDB recommended binding and both NetMHCIIpan 4.3 modes select TALDO1 core `YNYYKKFSY` (sequence-similarity rank 8/49); MixMHC2pred selects `VTKIYNYYK` (rank 37/49). BALF5 remains `YHFVKKHVH` across all three families. The DRB5 TALDO1 arm is therefore predictor-discordant and should be tested with both natural nested 11-mers before charge mutations are selected.", "",
        "`all_window_sequence_sensitivity.csv` contains all 49 window pairs for each target (147 total). These ranks describe sequence resemblance within one target and are not binding predictions or a discovery rerank.", "",
        "No peptide was ordered and no experiment was submitted.", "", "## Claim boundary", "", CLAIM, "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")

    checksum_rows = []
    for path in sorted(item for item in out.rglob("*") if item.is_file() and item.name != "SHA256SUMS.csv"):
        checksum_rows.append({"relative_path": str(path.relative_to(out)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
    write_csv(out / "SHA256SUMS.csv", checksum_rows)
    print(f"wrote {len(checksum_rows) + 1} files to {out}")


if __name__ == "__main__":
    main()
