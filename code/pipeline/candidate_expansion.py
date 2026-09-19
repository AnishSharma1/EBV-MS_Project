"""Pure selection and provenance logic for the candidate expansion screen."""

from __future__ import annotations

from collections import Counter
from itertools import permutations
from typing import Any, Iterable, Mapping, Optional, Sequence


CLAIM_BOUNDARY = (
    "This expansion identifies peptide-HLA hypotheses for independent predictor and "
    "provenance review. It does not establish natural presentation, TCR recognition, "
    "specificity, cross-reactivity, molecular mimicry, probability, false-discovery "
    "rate, or an MS mechanism."
)


def _number(value: Any) -> Optional[float]:
    try:
        if value in (None, "", "NA", "N/A"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _integer(value: Any) -> Optional[int]:
    number = _number(value)
    return int(number) if number is not None and number.is_integer() else None


def normalize_candidate_row(row: Mapping[str, Any], universe: str) -> dict[str, Any]:
    """Map V3 and combined-DR15 rows into one immutable selection schema."""
    combined = "ebv_core_p1_p9" in row
    normalized = dict(row)
    normalized.update(
        {
            "hla_rank": str(row.get("primary_rank") if combined else row.get("hla_rank", "")),
            "ebv_core": str(
                row.get("ebv_core_p1_p9") if combined else row.get("ebv_predicted_core", "")
            ).upper(),
            "self_core": str(
                row.get("self_core_p1_p9") if combined else row.get("self_predicted_core", "")
            ).upper(),
            "selection_universe": universe,
        }
    )
    return normalized


def normalize_candidate_rows(
    v3_rows: Iterable[Mapping[str, Any]],
    combined_dr15_rows: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows = [
        normalize_candidate_row(row, "v3_1600_pair_hla_universe")
        for row in v3_rows
        if row.get("allele") != "HLA-DRB1*15:01"
    ]
    rows.extend(
        normalize_candidate_row(row, "combined_drb1501_2043_pair_universe")
        for row in combined_dr15_rows
    )
    rows.sort(
        key=lambda row: (
            str(row.get("allele", "")),
            _integer(row.get("hla_rank")) or 10**9,
            str(row.get("pair_id", "")),
        )
    )
    return rows


def _eligibility_reason(
    row: Mapping[str, Any],
    frozen_pair_ids: set[str],
    frozen_ebv_ids: set[str],
    frozen_self_ids: set[str],
) -> str:
    if (
        str(row.get("pair_id", "")) in frozen_pair_ids
        or str(row.get("ebv_candidate_id", "")) in frozen_ebv_ids
        or str(row.get("self_candidate_id", "")) in frozen_self_ids
    ):
        return "excluded_frozen_pair_or_arm"
    ebv_binding = _number(row.get("ebv_binding_percentile_rank"))
    self_binding = _number(row.get("self_binding_percentile_rank"))
    if (
        ebv_binding is None
        or self_binding is None
        or ebv_binding > 20.0
        or self_binding > 20.0
    ):
        return "excluded_original_binding_percentile_above_20_or_missing"
    if (
        row.get("surface_status") != "complete"
        or _integer(row.get("left_model_count")) != 5
        or _integer(row.get("right_model_count")) != 5
    ):
        return "excluded_incomplete_five_model_ensemble"
    if row.get("declared_register_status") != "iedb_resolved_unique_fully_contained":
        return "excluded_declared_register_not_unique_or_contained"
    ebv_core = str(row.get("ebv_core", "")).upper()
    self_core = str(row.get("self_core", "")).upper()
    ebv_sequence = str(row.get("ebv_sequence", "")).upper()
    self_sequence = str(row.get("self_sequence", "")).upper()
    if (
        len(ebv_core) != 9
        or len(self_core) != 9
        or ebv_core not in ebv_sequence
        or self_core not in self_sequence
    ):
        return "excluded_invalid_or_uncontained_core"
    if _integer(row.get("hla_rank")) is None:
        return "excluded_missing_hla_rank"
    return "eligible"


def _select_for_order(
    rows_by_allele: Mapping[str, Sequence[Mapping[str, Any]]],
    allele_order: Sequence[str],
    *,
    target_count: int,
    max_per_allele: int,
    protein_cap: int,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    allele_counts: Counter[str] = Counter()
    ebv_protein_counts: Counter[str] = Counter()
    self_protein_counts: Counter[str] = Counter()
    ebv_cores: set[str] = set()
    self_cores: set[str] = set()
    protein_families: set[tuple[str, str]] = set()
    selected_ids: set[str] = set()

    changed = True
    while changed and len(selected) < target_count:
        changed = False
        for allele in allele_order:
            if len(selected) >= target_count:
                break
            if allele_counts[allele] >= max_per_allele:
                continue
            for candidate in rows_by_allele[allele]:
                pair_id = str(candidate["pair_id"])
                ebv_protein = str(candidate["ebv_protein"]).upper()
                self_protein = str(candidate["self_protein"]).upper()
                family = (ebv_protein, self_protein)
                if pair_id in selected_ids:
                    continue
                if candidate["ebv_core"] in ebv_cores or candidate["self_core"] in self_cores:
                    continue
                if family in protein_families:
                    continue
                if ebv_protein_counts[ebv_protein] >= protein_cap:
                    continue
                if self_protein_counts[self_protein] >= protein_cap:
                    continue
                selected.append(dict(candidate))
                selected_ids.add(pair_id)
                allele_counts[allele] += 1
                ebv_protein_counts[ebv_protein] += 1
                self_protein_counts[self_protein] += 1
                ebv_cores.add(str(candidate["ebv_core"]))
                self_cores.add(str(candidate["self_core"]))
                protein_families.add(family)
                changed = True
                break
    return selected


def select_diverse_candidates(
    rows: Sequence[Mapping[str, Any]],
    frozen_targets: Sequence[Mapping[str, Any]],
    *,
    target_count: int = 10,
    max_per_allele: int = 3,
    protein_cap: int = 2,
) -> dict[str, Any]:
    """Select a deterministic, diversity-constrained expansion without rescoring."""
    frozen_pair_ids = {str(row.get("pair_id", "")) for row in frozen_targets}
    frozen_ebv_ids = {str(row.get("ebv_candidate_id", "")) for row in frozen_targets}
    frozen_self_ids = {str(row.get("self_candidate_id", "")) for row in frozen_targets}
    provenance: list[dict[str, Any]] = []
    eligible: list[dict[str, Any]] = []
    for input_index, raw in enumerate(rows, start=1):
        row = dict(raw)
        reason = _eligibility_reason(row, frozen_pair_ids, frozen_ebv_ids, frozen_self_ids)
        provenance.append(
            {
                "input_index": input_index,
                "allele": row.get("allele", ""),
                "pair_id": row.get("pair_id", ""),
                "hla_rank": row.get("hla_rank", ""),
                "ebv_candidate_id": row.get("ebv_candidate_id", ""),
                "self_candidate_id": row.get("self_candidate_id", ""),
                "ebv_core": row.get("ebv_core", ""),
                "self_core": row.get("self_core", ""),
                "selection_status": reason,
            }
        )
        if reason == "eligible":
            eligible.append(row)

    eligible.sort(
        key=lambda row: (
            str(row["allele"]),
            _integer(row["hla_rank"]) or 10**9,
            str(row["pair_id"]),
        )
    )
    deduplicated: list[dict[str, Any]] = []
    retained_core_pairs: set[tuple[str, str, str]] = set()
    duplicate_ids: set[str] = set()
    for row in eligible:
        key = (str(row["allele"]), str(row["ebv_core"]), str(row["self_core"]))
        if key in retained_core_pairs:
            duplicate_ids.add(str(row["pair_id"]))
            continue
        retained_core_pairs.add(key)
        deduplicated.append(row)
    for row in provenance:
        if row["pair_id"] in duplicate_ids and row["selection_status"] == "eligible":
            row["selection_status"] = "excluded_same_hla_core_pair_duplicate"

    rows_by_allele: dict[str, list[dict[str, Any]]] = {}
    for row in deduplicated:
        rows_by_allele.setdefault(str(row["allele"]), []).append(row)
    alleles = sorted(rows_by_allele)
    alternatives = []
    for allele_order in permutations(alleles):
        selected = _select_for_order(
            rows_by_allele,
            allele_order,
            target_count=target_count,
            max_per_allele=max_per_allele,
            protein_cap=protein_cap,
        )
        counts = Counter(str(row["allele"]) for row in selected)
        alternatives.append(
            (
                -len(selected),
                -min((counts[allele] for allele in alleles), default=0),
                sum(_integer(row["hla_rank"]) or 10**9 for row in selected),
                tuple(sorted(str(row["pair_id"]) for row in selected)),
                allele_order,
                selected,
            )
        )
    best = min(alternatives) if alternatives else (0, 0, 0, (), (), [])
    selected = list(best[-1])
    selected.sort(
        key=lambda row: (
            str(row["allele"]),
            _integer(row["hla_rank"]) or 10**9,
            str(row["pair_id"]),
        )
    )
    selected_ids = {str(row["pair_id"]) for row in selected}
    for row in provenance:
        if row["pair_id"] in selected_ids:
            row["selection_status"] = "selected"
        elif row["selection_status"] == "eligible":
            row["selection_status"] = "eligible_not_selected_diversity_or_quota"
    provenance.sort(key=lambda row: int(row["input_index"]))
    return {
        "status": (
            "complete" if len(selected) == target_count
            else "not_evaluable_insufficient_diverse_candidates"
        ),
        "selected": selected,
        "provenance": provenance,
        "eligible_count_before_core_deduplication": len(eligible),
        "eligible_count_after_core_deduplication": len(deduplicated),
        "requested_target_count": target_count,
        "selected_target_count": len(selected),
    }


def resolve_natural_context(
    peptide: str,
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Resolve one peptide to a parent sequence and construct MixMHC2pred context."""
    peptide = peptide.strip().upper()
    matches = []
    for record in records:
        sequence = str(record.get("sequence", "")).upper()
        start = sequence.find(peptide)
        if start >= 0:
            matches.append(
                (
                    str(record.get("accession", "")),
                    str(record.get("protein", "")),
                    sequence,
                    start,
                )
            )
    if not matches:
        return {
            "source_accession": "",
            "source_protein_record": "",
            "source_start_1_based": "",
            "mixmhc_context": "",
            "mixmhc_context_status": "not_evaluable_missing_parent",
        }
    accession, protein, sequence, start = sorted(matches, key=lambda item: (item[0], item[3]))[0]
    end = start + len(peptide)
    upstream = sequence[max(0, start - 3) : start].rjust(3, "-")
    downstream = sequence[end : end + 3].ljust(3, "-")
    return {
        "source_accession": accession,
        "source_protein_record": protein,
        "source_start_1_based": start + 1,
        "mixmhc_context": upstream + peptide[:3] + peptide[-3:] + downstream,
        "mixmhc_context_status": "exact_parent_context",
    }


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def summarize_candidate_predictors(
    arm_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Aggregate independent predictor support without producing a composite score."""
    target_ids = sorted({str(row.get("target_id", "")) for row in arm_rows})
    summaries = []
    for target_id in target_ids:
        rows = [row for row in arm_rows if str(row.get("target_id", "")) == target_id]
        by_side = {str(row.get("side", "")): row for row in rows}
        ebv = by_side.get("ebv")
        self_arm = by_side.get("self")
        if not ebv or not self_arm:
            status = "not_evaluable_missing_arm"
        else:
            complete = all(str(row.get("predictor_status", "")) == "complete" for row in (ebv, self_arm))
            register_matches = all(
                _truth(row.get("register_consensus_matches_declared")) for row in (ebv, self_arm)
            )
            supported = all(_truth(row.get("binding_supported")) for row in (ebv, self_arm))
            strict = all(_truth(row.get("binding_consensus")) for row in (ebv, self_arm))
            if not complete:
                status = "not_evaluable_missing_predictor"
            elif strict and register_matches:
                status = "advance_strict_predictor_support"
            elif supported and register_matches:
                status = "advance_with_caution"
            else:
                status = "hold_register_or_binding_conflict"
        summaries.append(
            {
                "target_id": target_id,
                "pair_id": (ebv or self_arm or {}).get("pair_id", ""),
                "allele": (ebv or self_arm or {}).get("allele", ""),
                "expansion_role": (ebv or self_arm or {}).get("expansion_role", ""),
                "ebv_binding_consensus": bool(ebv and _truth(ebv.get("binding_consensus"))),
                "self_binding_consensus": bool(
                    self_arm and _truth(self_arm.get("binding_consensus"))
                ),
                "ebv_binding_supported": bool(ebv and _truth(ebv.get("binding_supported"))),
                "self_binding_supported": bool(
                    self_arm and _truth(self_arm.get("binding_supported"))
                ),
                "ebv_register_matches_declared": bool(
                    ebv and _truth(ebv.get("register_consensus_matches_declared"))
                ),
                "self_register_matches_declared": bool(
                    self_arm and _truth(self_arm.get("register_consensus_matches_declared"))
                ),
                "predictor_screen_status": status,
                "single_composite_score": "not_created",
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
    return summaries
