"""Run and summarize the independent predictor screen for the expansion package."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

from build_candidate_expansion import DEFAULT_OUT, read_csv, write_csv, write_json, _write_checksums
from build_high_yield_candidate_evidence import _fetch_netmhcii, _run_mixmhc
from candidate_expansion import CLAIM_BOUNDARY, summarize_candidate_predictors
from high_yield_candidate_evidence import summarize_predictor_evidence


DEFAULT_MIX_BINARY = Path(
    "/Users/anishsharma/.cache/ebv_ms_tools/"
    "mixmhc2pred_v2.1.beta1.2/MixMHC2pred"
)


def run_predictor_screen(
    output_dir: Path = DEFAULT_OUT,
    mixmhc_binary: Path = DEFAULT_MIX_BINARY,
) -> dict[str, Any]:
    protocol_path = output_dir / "protocol_lock.json"
    if not protocol_path.exists():
        raise FileNotFoundError("prepare the candidate expansion package before running predictors")
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    net_rows, net_failures = _fetch_netmhcii(output_dir)
    mix_rows, mix_failures = _run_mixmhc(output_dir, mixmhc_binary)
    predictor_rows = sorted(net_rows + mix_rows, key=lambda row: (row["arm_id"], row["predictor"]))
    predictor_fields = [
        "arm_id", "target_id", "allele", "predictor", "percentile_rank", "core",
        "core_reliability", "orientation", "score", "context_status", "raw_response",
    ]
    write_csv(output_dir / "raw_responses/predictor_records.csv", predictor_rows, predictor_fields)

    arm_screen = []
    for arm in arms:
        records = [row for row in predictor_rows if row["arm_id"] == arm["arm_id"]]
        arm_screen.append(
            {
                **arm,
                **summarize_predictor_evidence(records, arm["core"]),
                "predictor_record_count": len(records),
                "single_composite_score": "not_created",
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
    write_csv(output_dir / "peptide_arm_predictor_screen.csv", arm_screen)

    candidate_rows = summarize_candidate_predictors(arm_screen)
    target_details = {
        row["target_id"]: row
        for row in (
            read_csv(output_dir / "frozen_expansion_targets.csv")
            + read_csv(output_dir / "existing_lead_recheck_targets.csv")
        )
    }
    candidate_screen = []
    for row in candidate_rows:
        detail = target_details[row["target_id"]]
        candidate_screen.append(
            {
                **row,
                "ebv_protein": detail["ebv_protein"],
                "ebv_sequence": detail["ebv_sequence"],
                "ebv_core": detail["ebv_core"],
                "self_protein": detail["self_protein"],
                "self_sequence": detail["self_sequence"],
                "self_core": detail["self_core"],
                "upstream_hla_rank": detail["upstream_hla_rank"],
            }
        )
    candidate_screen.sort(key=lambda row: row["target_id"])
    write_csv(output_dir / "candidate_predictor_screen.csv", candidate_screen)

    failures = net_failures + mix_failures
    arm_status_counts = Counter(row["predictor_status"] for row in arm_screen)
    candidate_status_counts = Counter(row["predictor_screen_status"] for row in candidate_screen)
    new_status_counts = Counter(
        row["predictor_screen_status"]
        for row in candidate_screen
        if row["expansion_role"] == "new_candidate"
    )
    lead_status_counts = Counter(
        row["predictor_screen_status"]
        for row in candidate_screen
        if row["expansion_role"] != "new_candidate"
    )
    status = (
        "complete"
        if not failures and len(predictor_rows) == 4 * len(arms)
        and arm_status_counts.get("complete", 0) == len(arms)
        else "not_evaluable_incomplete_predictor_results"
    )
    gate = {
        "status": status,
        "predictor_record_count": len(predictor_rows),
        "expected_predictor_record_count": 4 * len(arms),
        "arm_status_counts": dict(sorted(arm_status_counts.items())),
        "candidate_status_counts": dict(sorted(candidate_status_counts.items())),
        "new_candidate_status_counts": dict(sorted(new_status_counts.items())),
        "lead_recheck_status_counts": dict(sorted(lead_status_counts.items())),
        "failures": failures,
        "decision_role": "screen_for_full_evidence_dossier_only",
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "cross_reactivity_claim_allowed": False,
        "molecular_mimicry_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "predictor_screen_gate.json", gate)

    lines = [
        "# Independent Predictor Screen",
        "",
        f"Status: `{status}`.",
        "",
        "This screen uses NetMHCIIpan 4.3 EL and MixMHC2pred 2.1 with corrected natural context. It does not change V3 ranks or create a composite score.",
        "",
        "## Counts",
        "",
    ]
    for key, value in sorted(new_status_counts.items()):
        lines.append(f"- New candidates `{key}`: {value}.")
    for key, value in sorted(lead_status_counts.items()):
        lines.append(f"- Existing lead rechecks `{key}`: {value}.")
    lines.extend(
        [
            "",
            "Advancement means only that both peptide arms retain independent binding support and the two predictors agree with the declared register. It is not evidence of natural presentation or TCR cross-reactivity.",
        ]
    )
    (output_dir / "PREDICTOR_SCREEN_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    _write_checksums(output_dir)
    return gate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--mixmhc-binary", type=Path, default=DEFAULT_MIX_BINARY)
    args = parser.parse_args()
    print(json.dumps(run_predictor_screen(args.output_dir, args.mixmhc_binary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
