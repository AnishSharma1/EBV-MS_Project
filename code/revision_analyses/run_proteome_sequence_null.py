#!/usr/bin/env python3
"""Screen BALF5 cores against every valid 9-mer in the reviewed human proteome.

This is deliberately a sequence-first proteome null. It prepares a smaller
top-window manifest for HLA-II binding prediction but does not label the full
proteome HLA-filtered unless those predictions are actually completed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import importlib.util
import json
from collections import Counter
from pathlib import Path


TCR_POSITIONS = (1, 2, 4, 6, 7)
TARGETS = {
    "HY13_SEQ_02": {
        "allele": "HLA-DRB1*13:03",
        "viral_core": "YHFVKKHVH",
        "taldo1_core": "LSFDKDAMV",
    },
    "HY15_SEQ_02": {
        "allele": "HLA-DRB1*15:01",
        "viral_core": "VYHFVKKHV",
        "taldo1_core": "IYNYYKKFS",
    },
}
CLAIM_BOUNDARY = (
    "Proteome-wide sequence-similarity null before HLA-binding filtering; not evidence of "
    "HLA binding, presentation, TCR recognition, cross-reactivity, or disease mechanism."
)


def parse_args():
    workspace = Path(__file__).resolve().parents[2]
    authoritative = Path(
        "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--human-fasta",
        type=Path,
        default=authoritative / "processed/high_yield_candidate_evidence_2026-08-28/raw_responses/human_reviewed_canonical.fasta",
    )
    parser.add_argument(
        "--score-source", type=Path, default=authoritative / "src/hla2_positive_control_benchmark_v2.py"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "output/revision_round_3/proteome_sequence_null_2026-09-15",
    )
    parser.add_argument("--top-windows", type=int, default=1000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fasta(path: Path):
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def accession(header: str) -> str:
    parts = header.split("|")
    return parts[1] if len(parts) >= 3 else header.split()[0]


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parent_15mer(sequence: str, core_start: int) -> tuple[str, int]:
    start = min(max(core_start - 3, 0), max(len(sequence) - 15, 0))
    peptide = sequence[start:start+15]
    return peptide, core_start - start + 1


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("score", args.score_source)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    valid = set(module._BLOSUM62_ALPHABET)
    records = list(parse_fasta(args.human_fasta))
    summaries = []
    top_output = []
    histogram_output = []
    binding_manifest = []

    for target_id, target in TARGETS.items():
        viral = target["viral_core"]
        taldo1_score = module.blosum62_similarity(viral, target["taldo1_core"], positions=TCR_POSITIONS)
        total = 0
        greater = 0
        equal = 0
        skipped = 0
        heap = []
        histogram = Counter()
        taldo1_occurrences = []
        serial = 0
        for header, sequence in records:
            acc = accession(header)
            for start in range(max(0, len(sequence) - 8)):
                core = sequence[start:start+9]
                if len(core) != 9:
                    continue
                if not set(core) <= valid:
                    skipped += 1
                    continue
                score = module.blosum62_similarity(viral, core, positions=TCR_POSITIONS)
                total += 1
                if score > taldo1_score + 1e-12:
                    greater += 1
                elif abs(score - taldo1_score) <= 1e-12:
                    equal += 1
                histogram[f"{round(score, 2):.2f}"] += 1
                peptide15, core_start_in_15 = parent_15mer(sequence, start)
                item = (score, serial, acc, header, start + 1, core, peptide15, core_start_in_15)
                serial += 1
                if len(heap) < args.top_windows:
                    heapq.heappush(heap, item)
                elif (score, -serial) > (heap[0][0], -heap[0][1]):
                    heapq.heapreplace(heap, item)
                if core == target["taldo1_core"] and ("TALDO1" in header or "TAL_HUMAN" in header):
                    taldo1_occurrences.append({"accession": acc, "start_1_based": start + 1, "header": header})
        ordered = sorted(heap, key=lambda item: (-item[0], item[2], item[4], item[1]))
        score_rank = greater + 1
        percentile = greater / total if total else None
        summaries.append({
            "target_id": target_id,
            "allele": target["allele"],
            "viral_core": viral,
            "taldo1_core": target["taldo1_core"],
            "taldo1_score": f"{taldo1_score:.12g}",
            "valid_human_9mers": total,
            "skipped_noncanonical_9mers": skipped,
            "strictly_greater_count": greater,
            "equal_score_count": equal,
            "competition_rank": score_rank,
            "strictly_greater_fraction": f"{percentile:.12g}",
            "taldo1_occurrence_count": len(taldo1_occurrences),
            "taldo1_occurrences_json": json.dumps(taldo1_occurrences, sort_keys=True),
            "hla_binding_filter_status": "not_yet_applied",
            "claim_boundary": CLAIM_BOUNDARY,
        })
        for display_rank, item in enumerate(ordered, 1):
            score, _, acc, header, start_1_based, core, peptide15, core_start_in_15 = item
            row = {
                "target_id": target_id,
                "allele": target["allele"],
                "display_rank": display_rank,
                "score": f"{score:.12g}",
                "accession": acc,
                "protein_header": header,
                "core_start_1_based": start_1_based,
                "human_core": core,
                "parent_15mer": peptide15,
                "core_start_in_parent_1_based": core_start_in_15,
                "claim_boundary": CLAIM_BOUNDARY,
            }
            top_output.append(row)
            binding_manifest.append({
                "submission_id": f"{target_id}_HUMAN_{display_rank:04d}",
                "target_id": target_id,
                "allele": target["allele"],
                "accession": acc,
                "parent_15mer": peptide15,
                "core_9mer": core,
                "core_start_in_parent_1_based": core_start_in_15,
                "sequence_prefilter_rank": display_rank,
                "submission_status": "prepared_not_submitted",
                "claim_boundary": CLAIM_BOUNDARY,
            })
        for score_bin, count in sorted(histogram.items(), key=lambda item: float(item[0])):
            histogram_output.append({
                "target_id": target_id,
                "score_bin_rounded_0_01": score_bin,
                "window_count": count,
            })

    write_csv(args.output_dir / "proteome_sequence_null_summary.csv", summaries)
    write_csv(args.output_dir / "top_1000_human_windows_per_lead.csv", top_output)
    write_csv(args.output_dir / "binding_prediction_manifest_top_1000_per_lead.csv", binding_manifest)
    write_csv(args.output_dir / "proteome_score_histogram.csv", histogram_output)
    lock = {
        "analysis": "BALF5_core_vs_reviewed_human_proteome_sequence_null",
        "status": "complete_sequence_stage_binding_stage_prepared_not_submitted",
        "human_fasta": str(args.human_fasta),
        "human_fasta_sha256": sha256(args.human_fasta),
        "human_protein_records": len(records),
        "score_source": str(args.score_source),
        "score_source_sha256": sha256(args.score_source),
        "top_windows_per_target": args.top_windows,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    lock["protocol_sha256"] = hashlib.sha256(json.dumps(lock, sort_keys=True).encode()).hexdigest()
    (args.output_dir / "analysis_lock.json").write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "README.md").write_text(
        "# Proteome-scale sequence null\n\n"
        "Every valid nine-residue window in the frozen reviewed human canonical FASTA was scored "
        "against the two BALF5 declared cores. This is not yet the HLA-binding-filtered null proposed "
        "by the reviewer. The top 1,000 windows per allele were exported as a prepared-not-submitted "
        "binding-prediction manifest so that HLA filtering can be completed without mislabeling the "
        "sequence-first result.\n",
        encoding="utf-8",
    )
    files = [p for p in args.output_dir.iterdir() if p.is_file() and p.name != "SHA256SUMS.csv"]
    write_csv(args.output_dir / "SHA256SUMS.csv", [
        {"path": p.name, "sha256": sha256(p)} for p in sorted(files)
    ])
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
