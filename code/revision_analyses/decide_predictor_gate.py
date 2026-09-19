#!/usr/bin/env python3
"""Apply the predeclared, predictor-specific 9-of-12 structural recovery gate."""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKQ_ACCEPTABLE = 0.23
REQUIRED_CONTROLS = 12
REQUIRED_RECOVERIES = 9


def as_float(row: dict, key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing/invalid {key} for {row.get('control_id')}/{row.get('predictor')}") from exc


def main(path: Path) -> int:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    grouped: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row.get("status") != "scored":
            continue
        grouped[row["predictor"]][row["control_id"]].append(row)

    report = {"rule": "at least 9 of exactly 12 controls with DockQ >= 0.23 and native_contact_recovery > 0", "predictors": {}}
    for predictor, by_control in grouped.items():
        if len(by_control) != REQUIRED_CONTROLS:
            report["predictors"][predictor] = {"decision": "INCOMPLETE", "scored_controls": len(by_control)}
            continue
        recovered = []
        for control, models in sorted(by_control.items()):
            best = max(models, key=lambda r: as_float(r, "dockq"))
            passes = any(as_float(m, "dockq") >= DOCKQ_ACCEPTABLE and as_float(m, "native_contact_recovery") > 0 for m in models)
            if passes:
                recovered.append(control)
            report.setdefault("per_control", {}).setdefault(predictor, {})[control] = {
                "best_dockq": as_float(best, "dockq"), "any_rank_passes": passes,
            }
        report["predictors"][predictor] = {
            "recovered_controls": recovered,
            "n_recovered": len(recovered),
            "decision": "PASS" if len(recovered) >= REQUIRED_RECOVERIES else "FAIL",
        }
    report["candidate_pose_generation_allowed"] = any(x.get("decision") == "PASS" for x in report["predictors"].values())
    target = ROOT / "manifests" / "predictor_gate_decision.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: decide_predictor_gate.py pose_metrics.tsv")
    raise SystemExit(main(Path(sys.argv[1])))
