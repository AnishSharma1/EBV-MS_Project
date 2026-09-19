"""Pure logic for the focused BALF5-TALDO1 alternate-register analysis."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Sequence

from build_same_register_hla_rankings_v2 import sequence_metrics


def classify_register_evidence(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Summarize predictor-family core agreement without calling it experimental proof."""
    by_family: dict[str, str] = {}
    for row in rows:
        family = str(row["predictor_family"])
        core = str(row["core"]).strip().upper()
        if len(core) != 9:
            raise ValueError(f"predictor core must contain nine residues: {core}")
        previous = by_family.setdefault(family, core)
        if previous != core:
            raise ValueError(f"one predictor family returned multiple cores: {family}")
    if not by_family:
        raise ValueError("at least one predictor family is required")
    counts = Counter(by_family.values())
    majority_core, majority_count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    return {
        "status": "predictor_register_consensus" if len(counts) == 1 else "predictor_register_disagreement",
        "predictor_family_count": len(by_family),
        "distinct_core_count": len(counts),
        "majority_core": majority_core,
        "majority_family_count": majority_count,
        "experimentally_resolved": False,
    }


def enumerate_pair_windows(left_sequence: str, right_sequence: str) -> list[dict[str, Any]]:
    """Enumerate every fully contained nine-residue window pair and sequence metric."""
    left_sequence = left_sequence.strip().upper()
    right_sequence = right_sequence.strip().upper()
    if len(left_sequence) < 9 or len(right_sequence) < 9:
        raise ValueError("both sequences must contain at least nine residues")
    rows: list[dict[str, Any]] = []
    for left_start in range(len(left_sequence) - 8):
        left_core = left_sequence[left_start:left_start + 9]
        for right_start in range(len(right_sequence) - 8):
            right_core = right_sequence[right_start:right_start + 9]
            rows.append({
                "left_start_1_based": left_start + 1,
                "right_start_1_based": right_start + 1,
                "left_core": left_core,
                "right_core": right_core,
                **sequence_metrics(left_core, right_core),
            })
    return rows
