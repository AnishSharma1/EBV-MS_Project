"""Build the additive, diversity-constrained high-yield candidate expansion."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from candidate_expansion import (
    CLAIM_BOUNDARY,
    normalize_candidate_rows,
    resolve_natural_context,
    select_diverse_candidates,
)


ROOT = Path(__file__).resolve().parents[1]
V3_DIR = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
EIGHT_DOSSIER = ROOT / "processed/high_yield_candidate_evidence_2026-08-28"
DEFAULT_OUT = ROOT / "processed/high_yield_candidate_expansion_2026-08-28"
HUMAN_FASTA = EIGHT_DOSSIER / "raw_responses/human_reviewed_canonical.fasta"
EBV_FASTA = EIGHT_DOSSIER / "raw_responses/ebv_uniprot_sequences.fasta"

ALLELE_TOKENS = {
    "HLA-DRB1*03:01": "03",
    "HLA-DRB1*08:01": "08",
    "HLA-DRB1*13:03": "13",
    "HLA-DRB1*15:01": "15",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] = (),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_fasta(path: Path) -> list[dict[str, str]]:
    records = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append(_fasta_record(header, "".join(sequence)))
            header = line[1:].strip()
            sequence = []
        else:
            sequence.append(line.strip())
    if header:
        records.append(_fasta_record(header, "".join(sequence)))
    return records


def _fasta_record(header: str, sequence: str) -> dict[str, str]:
    token = header.split()[0]
    accession = token.split("|")[1] if token.count("|") >= 2 else token
    return {
        "accession": accession,
        "protein": header,
        "sequence": re.sub(r"[^A-Za-z]", "", sequence).upper(),
    }


def _frozen_new_targets(selected: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    allele_counts: dict[str, int] = {}
    rows = []
    for candidate in selected:
        allele = str(candidate["allele"])
        allele_counts[allele] = allele_counts.get(allele, 0) + 1
        target_id = f"EXP{ALLELE_TOKENS[allele]}_{allele_counts[allele]:02d}"
        rows.append(
            {
                "target_id": target_id,
                "expansion_role": "new_candidate",
                "allele": allele,
                "upstream_hla_rank": candidate["hla_rank"],
                "upstream_rank_scope": candidate.get("rank_scope", ""),
                "selection_universe": candidate["selection_universe"],
                "pair_id": candidate["pair_id"],
                "ebv_candidate_id": candidate["ebv_candidate_id"],
                "ebv_protein": candidate["ebv_protein"],
                "ebv_sequence": candidate["ebv_sequence"],
                "ebv_core": candidate["ebv_core"],
                "ebv_binding_percentile_rank": candidate["ebv_binding_percentile_rank"],
                "self_candidate_id": candidate["self_candidate_id"],
                "self_protein": candidate["self_protein"],
                "self_sequence": candidate["self_sequence"],
                "self_core": candidate["self_core"],
                "self_binding_percentile_rank": candidate["self_binding_percentile_rank"],
                "tcr_facing_blosum62_similarity": candidate.get("tcr_facing_blosum62_similarity", ""),
                "tcr_face_physicochemical_mismatch": candidate.get(
                    "tcr_face_physicochemical_mismatch", ""
                ),
                "tcr_facing_sequence_identity": candidate.get("tcr_facing_sequence_identity", ""),
                "local_surface_percentile": candidate.get("local_surface_percentile", ""),
                "left_model_count": candidate["left_model_count"],
                "right_model_count": candidate["right_model_count"],
                "declared_register_status": candidate["declared_register_status"],
                "surface_status": candidate["surface_status"],
                "computational_pair_marker": "*",
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
    return rows


def _lead_recheck_targets() -> list[dict[str, Any]]:
    recommendations = read_csv(EIGHT_DOSSIER / "stage1_assay_recommendations.csv")
    registry = {
        row["target_id"]: row for row in read_csv(EIGHT_DOSSIER / "frozen_candidate_registry.csv")
    }
    rows = []
    for recommendation in recommendations:
        if recommendation["stage1_status"] != "stage1_medium_priority":
            continue
        source = registry[recommendation["target_id"]]
        rows.append(
            {
                "target_id": recommendation["target_id"],
                "expansion_role": "existing_medium_priority_lead_context_recheck",
                "allele": recommendation["allele"],
                "upstream_hla_rank": "",
                "upstream_rank_scope": "frozen_eight_candidate_dossier",
                "selection_universe": "frozen_eight_candidate_dossier",
                "pair_id": recommendation["pair_id"],
                "ebv_candidate_id": source["ebv_candidate_id"],
                "ebv_protein": recommendation["ebv_protein"],
                "ebv_sequence": recommendation["ebv_sequence"],
                "ebv_core": recommendation["ebv_core"],
                "ebv_binding_percentile_rank": source["ebv_binding_percentile_rank"],
                "self_candidate_id": source["self_candidate_id"],
                "self_protein": recommendation["self_protein"],
                "self_sequence": recommendation["self_sequence"],
                "self_core": recommendation["self_core"],
                "self_binding_percentile_rank": source["self_binding_percentile_rank"],
                "left_model_count": source["left_model_count"],
                "right_model_count": source["right_model_count"],
                "declared_register_status": "frozen_declared_core",
                "surface_status": source["surface_status"],
                "computational_pair_marker": "*",
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
    rows.sort(key=lambda row: row["target_id"])
    if len(rows) != 2:
        raise ValueError(f"expected exactly two medium-priority lead rechecks, found {len(rows)}")
    return rows


def _arm_registry(
    targets: Sequence[Mapping[str, Any]],
    human_records: Sequence[Mapping[str, Any]],
    ebv_records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows = []
    for target in targets:
        for side, records in (("ebv", ebv_records), ("self", human_records)):
            sequence = str(target[f"{side}_sequence"]).upper()
            core = str(target[f"{side}_core"]).upper()
            context = resolve_natural_context(sequence, records)
            rows.append(
                {
                    "arm_id": f"{target['target_id']}__{side}",
                    "target_id": target["target_id"],
                    "expansion_role": target["expansion_role"],
                    "pair_id": target["pair_id"],
                    "side": side,
                    "allele": target["allele"],
                    "candidate_id": target[f"{side}_candidate_id"],
                    "protein": target[f"{side}_protein"],
                    "sequence": sequence,
                    "core": core,
                    "declared_core_start_1_based": sequence.index(core) + 1,
                    **context,
                    "claim_boundary": CLAIM_BOUNDARY,
                }
            )
    rows.sort(key=lambda row: row["arm_id"])
    return rows


def _predictor_manifest(arms: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    groups: dict[tuple[str, int], list[Mapping[str, Any]]] = {}
    for arm in arms:
        groups.setdefault((str(arm["allele"]), len(str(arm["sequence"]))), []).append(arm)
    for (allele, length), group in sorted(groups.items()):
        batch_id = (
            allele.replace("HLA-", "").replace("*", "_").replace(":", "_")
            + f"__len{length}"
        )
        for seq_num, arm in enumerate(sorted(group, key=lambda row: str(row["arm_id"])), start=1):
            rows.append(
                {
                    "batch_id": batch_id,
                    "seq_num": seq_num,
                    "arm_id": arm["arm_id"],
                    "target_id": arm["target_id"],
                    "allele": allele,
                    "peptide_length": length,
                    "exact_sequence": arm["sequence"],
                    "declared_core": arm["core"],
                    "mixmhc_context": arm["mixmhc_context"],
                    "mixmhc_context_status": arm["mixmhc_context_status"],
                    "netmhciipan_methods": "4.3_EL;4.3_BA",
                    "mixmhc2pred_methods": "2.1_context;2.1_no_context",
                    "analysis_status": "prepared_not_run",
                }
            )
    return rows


def _write_checksums(output_dir: Path) -> None:
    files = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.csv"
    )
    write_csv(
        output_dir / "SHA256SUMS.csv",
        [
            {"relative_path": str(path.relative_to(output_dir)), "sha256": sha256_file(path)}
            for path in files
        ],
        ("relative_path", "sha256"),
    )


def prepare_expansion_package(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    required = [
        V3_DIR / "v3_all_hla_ranked_pairs.csv",
        V3_DIR / "combined_drb1501_v3_ranked_pairs.csv",
        EIGHT_DOSSIER / "frozen_candidate_registry.csv",
        EIGHT_DOSSIER / "stage1_assay_recommendations.csv",
        HUMAN_FASTA,
        EBV_FASTA,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing required expansion inputs: " + "; ".join(missing))
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = normalize_candidate_rows(
        read_csv(V3_DIR / "v3_all_hla_ranked_pairs.csv"),
        read_csv(V3_DIR / "combined_drb1501_v3_ranked_pairs.csv"),
    )
    frozen = read_csv(EIGHT_DOSSIER / "frozen_candidate_registry.csv")
    selection = select_diverse_candidates(rows, frozen)
    if selection["status"] != "complete":
        raise ValueError(
            "the locked rules did not yield ten diverse candidates: "
            f"{selection['selected_target_count']}"
        )
    new_targets = _frozen_new_targets(selection["selected"])
    lead_targets = _lead_recheck_targets()
    human_records = parse_fasta(HUMAN_FASTA)
    ebv_records = parse_fasta(EBV_FASTA)
    arms = _arm_registry(new_targets + lead_targets, human_records, ebv_records)

    write_csv(output_dir / "frozen_expansion_targets.csv", new_targets)
    write_csv(output_dir / "existing_lead_recheck_targets.csv", lead_targets)
    write_csv(output_dir / "selection_provenance.csv", selection["provenance"])
    write_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv", arms)
    write_csv(output_dir / "prepared_inputs/predictor_query_manifest.csv", _predictor_manifest(arms))
    with (output_dir / "prepared_inputs/mixmhc2pred_context.tsv").open("w", encoding="utf-8") as handle:
        for arm in arms:
            handle.write(f"{arm['sequence']}\t{arm['mixmhc_context']}\n")
    with (output_dir / "prepared_inputs/mixmhc2pred_no_context.tsv").open("w", encoding="utf-8") as handle:
        for arm in arms:
            handle.write(f"{arm['sequence']}\n")

    old_arms = {
        row["arm_id"]: row
        for row in read_csv(EIGHT_DOSSIER / "prepared_inputs/peptide_arm_registry.csv")
    }
    correction_rows = []
    for arm in arms:
        if arm["expansion_role"] != "existing_medium_priority_lead_context_recheck":
            continue
        old = old_arms[arm["arm_id"]]
        correction_rows.append(
            {
                "arm_id": arm["arm_id"],
                "protein": arm["protein"],
                "sequence": arm["sequence"],
                "old_context": old["mixmhc_context"],
                "old_context_status": old["mixmhc_context_status"],
                "corrected_context": arm["mixmhc_context"],
                "corrected_context_status": arm["mixmhc_context_status"],
                "source_accession": arm["source_accession"],
                "correction_status": (
                    "corrected_in_additive_recheck"
                    if old["mixmhc_context_status"] != arm["mixmhc_context_status"]
                    else "unchanged_exact_context"
                ),
            }
        )
    write_csv(output_dir / "lead_context_correction.csv", correction_rows)

    source_checksums = {str(path.relative_to(ROOT)): sha256_file(path) for path in required}
    allele_counts: dict[str, int] = {}
    for row in new_targets:
        allele_counts[row["allele"]] = allele_counts.get(row["allele"], 0) + 1
    context_missing = [arm["arm_id"] for arm in arms if arm["mixmhc_context_status"] != "exact_parent_context"]
    protocol = {
        "protocol_id": "high_yield_candidate_expansion_2026-08-28",
        "status": "prepared_not_evaluated",
        "selection_method": "frozen_V3_rank_order_with_prospective_diversity_constraints",
        "new_target_count": len(new_targets),
        "new_target_count_by_hla": allele_counts,
        "lead_recheck_count": len(lead_targets),
        "peptide_arm_count": len(arms),
        "selection_rules": {
            "original_binding_percentile_max_each_arm": 20,
            "complete_model_count_each_arm": 5,
            "declared_register": "IEDB-resolved unique fully contained nonamer",
            "original_eight_pair_and_arm_reuse": "forbidden",
            "same_hla_core_pair_duplicate": "keep best frozen V3 rank only",
            "global_exact_ebv_core_reuse": "forbidden",
            "global_exact_self_core_reuse": "forbidden",
            "global_ebv_self_protein_family_reuse": "forbidden",
            "maximum_use_per_source_protein": 2,
            "maximum_candidates_per_hla": 3,
            "geometry_inspected_during_selection": False,
            "new_composite_created": False,
        },
        "diversity_result": {
            "eligible_before_same_hla_core_deduplication": selection[
                "eligible_count_before_core_deduplication"
            ],
            "eligible_after_same_hla_core_deduplication": selection[
                "eligible_count_after_core_deduplication"
            ],
            "honest_diverse_shortlist_size": len(new_targets),
            "near_duplicates_added_to_reach_a_larger_count": False,
        },
        "lead_context_recheck_reason": (
            "The two medium-priority TALDO1 arms were resolved in the cached canonical "
            "human proteome after the original prepare stage and require corrected natural-flank runs."
        ),
        "missing_natural_context_arm_ids": context_missing,
        "predictor_status": "prepared_not_run",
        "existing_packages_modified": False,
        "discovery_unlock_allowed": False,
        "weights_frozen": False,
        "specificity_claim_allowed": False,
        "cross_reactivity_claim_allowed": False,
        "molecular_mimicry_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
        "source_checksums": source_checksums,
    }
    write_json(output_dir / "protocol_lock.json", protocol)
    write_json(
        output_dir / "selection_summary.json",
        {
            "selection_status": selection["status"],
            "new_target_count": len(new_targets),
            "new_target_count_by_hla": allele_counts,
            "lead_recheck_count": len(lead_targets),
            "all_natural_contexts_resolved": not context_missing,
            "claim_boundary": CLAIM_BOUNDARY,
        },
    )
    (output_dir / "README.md").write_text(
        """# High-Yield Candidate Expansion

This additive package freezes ten fresh HLA-specific hypotheses from the existing V3 universes. It does not rerank the discovery library. Selection used the frozen V3 order only after requiring both original binding percentiles at or below 20, complete five-model ensembles, and a uniquely contained declared register.

An initially larger list would have required repeated cores, nested versions of the same peptide, or repeated EBV-human protein families. Those were not added. The ten retained candidates are therefore a diversity-constrained screening set, not ten proven independent biological systems.

The two medium-priority BALF5-TALDO1 leads are included separately for corrected MixMHC2pred natural-context runs. The original eight-candidate dossier remains unchanged.

## State

- New candidates: 10.
- Existing lead rechecks: 2.
- Peptide arms prepared: 24.
- Independent predictor and provenance results: not yet collected.
- Discovery unlock, specificity, TCR cross-reactivity, molecular mimicry, and MS-mechanism claims: not allowed.

## Files

- `frozen_expansion_targets.csv`: the ten new hypotheses.
- `existing_lead_recheck_targets.csv`: the two prior leads, kept separate.
- `selection_provenance.csv`: every eligibility and exclusion decision.
- `lead_context_correction.csv`: old versus corrected lead contexts.
- `prepared_inputs/`: exact arm registry and predictor inputs.
- `protocol_lock.json`: frozen rules, source checksums, and claim boundary.
""",
        encoding="utf-8",
    )
    _write_checksums(output_dir)
    return {
        "new_target_count": len(new_targets),
        "lead_recheck_count": len(lead_targets),
        "peptide_arm_count": len(arms),
        "all_natural_contexts_resolved": not context_missing,
        "output_dir": str(output_dir),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(prepare_expansion_package(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
