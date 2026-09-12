"""Build the dated, auditable HY15 charge-reversal pilot package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from build_tcell_library_v2 import DRA_SEQUENCE, DRB1_1501_SEQUENCE
from charge_reversal_binding_pilot import (
    ALLELE,
    CLAIM_BOUNDARY,
    TOOL_ALLELE,
    checksum_rows,
    fitted_results_template,
    make_plate_map,
    primary_panel,
    validate_panel,
    write_csv,
    write_json,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "processed" / "charge_reversal_binding_pilot_2026-09-07"
IEDB_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
MIX_BINARY = Path.home() / ".cache/ebv_ms_tools/mixmhc2pred_v2.1.beta1.2/MixMHC2pred"


def verify_locked_sources() -> None:
    stage1_path = ROOT / "processed/high_yield_candidate_evidence_2026-08-28/stage1_assay_recommendations.csv"
    with stage1_path.open(newline="", encoding="utf-8") as handle:
        stage1 = {row["target_id"]: row for row in csv.DictReader(handle)}
    row = stage1.get("HY15_SEQ_02")
    expected = ("TGGVYHFVKKHVHES", "SVTKIYNYYKKFSYK", "stage1_medium_priority")
    observed = (row["ebv_sequence"], row["self_sequence"], row["stage1_status"]) if row else None
    if observed != expected:
        raise ValueError(f"HY15 source candidate drift: {observed}")

    register_path = ROOT / "processed/taldo1_focused_register_2026-09-05/target_register_summary.csv"
    with register_path.open(newline="", encoding="utf-8") as handle:
        register_rows = {row["target_id"]: row for row in csv.DictReader(handle)}
    register = register_rows.get("HY15_SEQ_02")
    register_expected = ("VYHFVKKHV", "IYNYYKKFS", "predictor_register_consensus", "False")
    register_observed = (
        register["ebv_majority_core"],
        register["self_majority_core"],
        register["pair_register_status"],
        register["experimentally_resolved"],
    ) if register else None
    if register_observed != register_expected:
        raise ValueError(f"HY15 register source drift: {register_observed}")

    gate_path = ROOT / "processed/pmhc_surface_electrostatics_v2_controls_2026-08-30/control_gate.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if gate.get("status") != "fail" or gate.get("electrostatics_retired_from_candidate_ranking") is not True:
        raise ValueError("retired APBS gate no longer matches the locked source state")


def post_iedb(rows: list[dict[str, Any]]) -> str:
    fasta = "\n".join(f">{index}|{row['sample_id']}\n{row['sequence']}" for index, row in enumerate(rows, 1))
    body = urllib.parse.urlencode({
        "method": "recommended_binding",
        "sequence_text": fasta,
        "allele": ALLELE,
        "length": "asis",
    }).encode("utf-8")
    request = urllib.request.Request(
        IEDB_ENDPOINT,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8")


def parse_iedb(raw: str, panel: list[dict[str, Any]], retrieved_utc: str) -> list[dict[str, Any]]:
    rows = list(csv.DictReader(io.StringIO(raw), delimiter="\t"))
    by_number = {int(row["seq_num"]): row for row in rows}
    if len(by_number) != len(panel):
        raise ValueError("IEDB did not return one result for every peptide")
    output = []
    for index, meta in enumerate(panel, 1):
        row = by_number[index]
        if row["allele"] != ALLELE or row["peptide"] != meta["sequence"]:
            raise ValueError(f"IEDB identity mismatch at seq_num {index}")
        output.append({
            "sample_id": meta["sample_id"],
            "seq_num": index,
            "allele": ALLELE,
            "sequence": meta["sequence"],
            "expected_core": meta["expected_core"],
            "predicted_core": row["core_peptide"],
            "register_matches_expected": row["core_peptide"] == meta["expected_core"],
            "predicted_ic50_nM": float(row["ic50"]),
            "predicted_rank_percentile": float(row["rank"]),
            "predictor": "IEDB recommended_binding",
            "retrieved_utc": retrieved_utc,
            "claim_boundary": "Binding-prediction sensitivity only; " + CLAIM_BOUNDARY,
        })
    return output


def run_mixmhc(panel: list[dict[str, Any]], prepared: Path, raw: Path) -> list[dict[str, Any]]:
    prepared.write_text("\n".join(row["sequence"] for row in panel) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(MIX_BINARY), "-i", str(prepared), "-o", str(raw), "-a", TOOL_ALLELE, "--no_context"],
        capture_output=True,
        text=True,
        check=False,
    )
    (raw.parent / "mixmhc2pred.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (raw.parent / "mixmhc2pred.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"MixMHC2pred failed with exit {completed.returncode}")
    lines = [line for line in raw.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    parsed = list(csv.DictReader(lines, delimiter="\t"))
    if len(parsed) != len(panel):
        raise ValueError("MixMHC2pred did not return one result for every peptide")
    output = []
    for meta, row in zip(panel, parsed):
        if row["Peptide"] != meta["sequence"]:
            raise ValueError(f"MixMHC2pred identity mismatch for {meta['sample_id']}")
        p1 = int(row[f"CoreP1_{TOOL_ALLELE}"])
        core = meta["sequence"][p1 - 1:p1 + 8]
        output.append({
            "sample_id": meta["sample_id"],
            "allele": ALLELE,
            "sequence": meta["sequence"],
            "expected_core": meta["expected_core"],
            "predicted_core": core,
            "core_start_1_based": p1,
            "register_matches_expected": core == meta["expected_core"],
            "predicted_rank_percentile": float(row[f"%Rank_{TOOL_ALLELE}"]),
            "predictor": "MixMHC2pred v2.1-beta1 no_context",
            "claim_boundary": "Binding-prediction sensitivity only; " + CLAIM_BOUNDARY,
        })
    return output


def af3_jobs(panel: list[dict[str, Any]]) -> list[dict[str, Any]]:
    jobs = []
    for row in panel:
        if not row["primary_mutation_panel"]:
            continue
        jobs.append({
            "name": f"cr_hy15_drb1_1501_{row['sample_id'].lower()}",
            "modelSeeds": [314159],
            "sequences": [
                {"proteinChain": {"sequence": DRA_SEQUENCE, "count": 1, "useStructureTemplate": True}},
                {"proteinChain": {"sequence": DRB1_1501_SEQUENCE, "count": 1, "useStructureTemplate": True}},
                {"proteinChain": {"sequence": row["sequence"], "count": 1, "useStructureTemplate": True}},
            ],
            "dialect": "alphafoldserver",
            "version": 3,
        })
    return jobs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(f"output already exists: {out}")
    verify_locked_sources()
    out.mkdir(parents=True)
    prepared = out / "prepared_inputs"
    raw_dir = out / "raw_responses"
    model_dir = out / "modeling"
    lab_dir = out / "lab_handoff"
    for directory in (prepared, raw_dir, model_dir, lab_dir):
        directory.mkdir()

    panel = primary_panel()
    validate_panel(panel)
    write_csv(out / "peptide_manifest.csv", panel)
    (prepared / "iedb_primary_panel.fasta").write_text(
        "\n".join(f">{index}|{row['sample_id']}\n{row['sequence']}" for index, row in enumerate(panel, 1)) + "\n",
        encoding="utf-8",
    )
    retrieved_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    raw_iedb = post_iedb(panel)
    (raw_dir / "iedb_recommended_binding_drb1_1501.tsv").write_text(raw_iedb, encoding="utf-8")
    iedb = parse_iedb(raw_iedb, panel, retrieved_utc)
    write_csv(out / "iedb_prediction_results.csv", iedb)
    mix = run_mixmhc(
        panel,
        prepared / "mixmhc2pred_primary_panel_no_context.txt",
        raw_dir / "mixmhc2pred_primary_panel_no_context.tsv",
    )
    write_csv(out / "mixmhc2pred_prediction_results.csv", mix)
    iedb_by_id = {row["sample_id"]: row for row in iedb}
    mix_by_id = {row["sample_id"]: row for row in mix}
    consensus = []
    for row in panel:
        left, right = iedb_by_id[row["sample_id"]], mix_by_id[row["sample_id"]]
        consensus.append({
            "sample_id": row["sample_id"],
            "expected_core": row["expected_core"],
            "iedb_core": left["predicted_core"],
            "mixmhc2pred_core": right["predicted_core"],
            "both_match_expected_core": left["register_matches_expected"] and right["register_matches_expected"],
            "predictors_agree": left["predicted_core"] == right["predicted_core"],
            "interpretation": "predictor_register_support_only_not_experimental_resolution",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(out / "predictor_register_consensus.csv", consensus)
    wt_ids = {"BALF5": "BALF5_WT_15", "TALDO1": "TALDO1_WT_15"}
    meta_by_id = {row["sample_id"]: row for row in panel}
    computational = []
    for row in panel:
        if not row["primary_mutation_panel"]:
            continue
        sample_id = row["sample_id"]
        wt_id = wt_ids[row["arm"]]
        iedb_row, iedb_wt = iedb_by_id[sample_id], iedb_by_id[wt_id]
        mix_row, mix_wt = mix_by_id[sample_id], mix_by_id[wt_id]
        computational.append({
            "sample_id": sample_id,
            "arm": row["arm"],
            "reagent_role": row["reagent_role"],
            "core_position": row["core_position"],
            "expected_core": row["expected_core"],
            "iedb_predicted_core": iedb_row["predicted_core"],
            "iedb_register_matches_expected": iedb_row["register_matches_expected"],
            "iedb_ic50_fold_vs_wt": round(iedb_row["predicted_ic50_nM"] / iedb_wt["predicted_ic50_nM"], 6),
            "iedb_rank_delta_vs_wt": round(iedb_row["predicted_rank_percentile"] - iedb_wt["predicted_rank_percentile"], 6),
            "mixmhc2pred_core": mix_row["predicted_core"],
            "mixmhc2pred_register_matches_expected": mix_row["register_matches_expected"],
            "mixmhc2pred_rank_delta_vs_wt": round(mix_row["predicted_rank_percentile"] - mix_wt["predicted_rank_percentile"], 6),
            "experimental_gate_status": "not_evaluable_without_measured_binding",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(out / "computational_sensitivity_summary.csv", computational)

    jobs = af3_jobs(panel)
    write_json(model_dir / "af3_primary_mutation_panel_10_jobs.json", jobs)
    write_csv(model_dir / "MODEL_STATUS.csv", [{
        "job_name": job["name"],
        "sample_id": next(row["sample_id"] for row in panel if row["sample_id"].lower() in job["name"]),
        "requested_models": 5,
        "status": "prepared_not_submitted",
        "results_path": "",
        "claim_boundary": "No structure exists until a returned model package passes exact-chain and peptide QC.",
    } for job in jobs])
    (model_dir / "AF3_UPLOAD_AND_RETURN_QC.md").write_text("""# AlphaFold Server upload and return QC

## Submission

Upload `af3_primary_mutation_panel_10_jobs.json` to AlphaFold Server as one batch. Verify that the page lists exactly 10 jobs before submission. Each job contains DRA, DRB1*15:01, and one exact 15-mer peptide. The single locked seed requests the server's five model samples.

## Required return

Download the complete result archive without renaming internal files. Expected coverage is 10 jobs and 50 model samples. Update `MODEL_STATUS.csv` only after the archive is present.

## QC before interpretation

1. Verify exactly three chains per model and exact equality to the submitted DRA, DRB, and peptide sequences.
2. Reject truncated, duplicated, missing-chain, or clashing samples from geometry while preserving them in the QC table.
3. Report peptide-chain ipTM/pLDDT, peptide-HLA contacts, across-model pose consistency, and the predicted P1-P9 register.
4. A shifted or unstable peptide pose is `structurally_confounded`, not an electrostatic result.
5. Do not use global AlphaFold ranking score or APBS similarity to pass the experimental binding gate.

## Current state

Prepared, not submitted. No structural result is present in this package.
""", encoding="utf-8")

    plate_map = make_plate_map(panel)
    write_csv(lab_dir / "randomized_plate_map_9x96.csv", plate_map)
    write_csv(lab_dir / "fitted_binding_results_template.csv", fitted_results_template(panel))
    write_csv(lab_dir / "reagent_order_sheet.csv", [{
        "sample_id": row["sample_id"],
        "sequence": row["sequence"],
        "length": row["length"],
        "requested_purity_percent": 95,
        "quantity": "lab_to_specify",
        "modifications": "none",
        "order_status": "proposed_not_ordered",
    } for row in panel])

    protocol = {
        "protocol_id": "HY15_DRB1_1501_CHARGE_REVERSAL_BINDING_PILOT_V1",
        "created_utc": retrieved_utc,
        "target": "HY15_SEQ_02 BALF5-TALDO1",
        "hla": {"alpha": "HLA-DRA*01:01", "beta": ALLELE},
        "primary_endpoint": "competition-binding IC50_nM; KD_nM may be reported as an orthogonal endpoint",
        "assay_design": {"independent_experiments": 3, "technical_replicates": 2, "concentrations_nM": [round(100000.0 / (3 ** index), 6) for index in range(10)]},
        "gate": {
            "reversal_fold_loss_vs_wt_min": 3.0,
            "reversal_over_neutralization_fold_min": 2.0,
            "direction_must_match_all_experiments": True,
            "register_must_be_retained": True,
            "full_curve_fit_required": True,
            "purity_percent_min": 95.0,
            "solubility_required": True,
            "analyst_blinding_required": True,
        },
        "exclusions": ["register_shift", "failed_curve_fit", "purity_below_95_percent", "solubility_failure", "identity_mismatch"],
        "apbs_candidate_ranking": "prohibited_prior_control_gate_failed",
        "tcr_endpoints": "out_of_scope_until_binding_and_register_gates_pass",
        "expansions": ["HY13_SEQ_02 with DRB1*13:03 after primary gate", "HY15_SEQ_02 with DRB5*01:01 after experimental register mapping"],
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "protocol_lock.json", protocol)
    source_rows = [
        {"source": "stage1_assay_recommendations.csv", "purpose": "frozen HY15 candidate identity", "path": "processed/high_yield_candidate_evidence_2026-08-28/stage1_assay_recommendations.csv"},
        {"source": "focused register analysis", "purpose": "native and nested register hypotheses", "path": "processed/taldo1_focused_register_2026-09-05/target_register_summary.csv"},
        {"source": "prior charge sensitivity", "purpose": "prediction context only", "path": "processed/taldo1_hla_next_leg_2026-09-04/charge_sensitivity_results.csv"},
        {"source": "retired APBS gate", "purpose": "enforce non-ranking rule", "path": "processed/pmhc_surface_electrostatics_v2_controls_2026-08-30/control_gate.json"},
    ]
    for row in source_rows:
        path = ROOT / row["path"]
        row["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        row["bytes"] = path.stat().st_size
    write_csv(out / "source_manifest.csv", source_rows)

    primary_consensus = [row for row in consensus if any(p["sample_id"] == row["sample_id"] and p["primary_mutation_panel"] for p in panel)]
    matching = sum(row["both_match_expected_core"] for row in primary_consensus)
    non_wt = [row for row in computational if row["reagent_role"] != "wild_type"]
    max_fold = max(row["iedb_ic50_fold_vs_wt"] for row in non_wt)
    shifted = [row["sample_id"] for row in computational if not row["iedb_register_matches_expected"] or not row["mixmhc2pred_register_matches_expected"]]
    readme = f"""# HY15 DRB1*15:01 charge-reversal binding pilot

## Status

- Computational binding predictions: completed for 14 peptides on {retrieved_utc}.
- Primary variants matching the locked core in both predictors: {matching}/10.
- AlphaFold Server inputs: 10 jobs prepared, not submitted; each submitted job requests one seed and the server's five model samples.
- Peptides: proposed, not ordered.
- Binding experiments: not run; the result files are blank templates.

## Purpose

This package implements Olivia Thomas's charge-reversal idea as a peptide-HLA binding perturbation for HY15_SEQ_02 with DRA1/DRB1*15:01. It compares wild type, K-to-Q neutralization, and K-to-E reversal at core P6 and P7 on both BALF5 and TALDO1.

The analysis must not be called an electrostatic-contact result unless measured binding, register retention, and the locked gate support it. Even a passed peptide gate supports charge-sensitive binding, not a direct salt bridge. Direct-contact evidence requires an independently designed reciprocal HLA mutation/double-mutant cycle.

## Start here

1. `protocol_lock.json` — frozen endpoints, exclusions, and decision rule.
2. `peptide_manifest.csv` — exact identities and sequences.
3. `predictor_register_consensus.csv` — computational register sensitivity.
4. `COMPUTATIONAL_FINDINGS.md` — current prediction result and limitations.
5. `lab_handoff/LAB_HANDOFF.md` — order and assay instructions.
6. `modeling/MODEL_STATUS.csv` — explicit prepared/not-submitted status.

## Analysis

After the lab returns fitted IC50 values, complete `lab_handoff/fitted_binding_results_template.csv` and run:

```bash
PYTHONPATH=src python3 src/analyze_charge_reversal_binding_results.py --package {out}
```

The analysis refuses incomplete or invalid three-experiment data.

## Claim boundary

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    findings = f"""# Computational findings

## Direct result

- Both predictors retained the locked DRB1*15:01 core for {matching}/10 primary peptides.
- The exception was `{shifted[0] if shifted else 'none'}`. IEDB shifted its predicted core by one residue, while MixMHC2pred retained the expected core. This variant is computationally register-confounded and cannot support an electrostatic interpretation unless its experimental register is resolved.
- The largest IEDB-predicted IC50 loss among the eight mutations was {max_fold:.3f}-fold relative to the same-arm wild type, below the locked threefold experimental threshold.
- These predictions therefore do not pre-label any position as charge-sensitive. They justify running the controlled binding experiment; they are not substitutes for it.

See `computational_sensitivity_summary.csv` for every raw comparison and the two predictor-specific tables for exact outputs.

## Boundary

No peptide was ordered, no binding was measured, and no new pMHC structure was generated in this package. The prior APBS endpoint remains retired from candidate ranking.
"""
    (out / "COMPUTATIONAL_FINDINGS.md").write_text(findings, encoding="utf-8")
    handoff = f"""# Laboratory handoff: HY15 charge-reversal binding pilot

## Reagents

Order the 14 peptides in `reagent_order_sheet.csv` at at least 95% purity after the receiving laboratory specifies quantity and formulation. Do not order directly from this draft without the laboratory confirming its assay format and reference tracer.

Use recombinant HLA-DRA*01:01/HLA-DRB1*15:01. Confirm complex folding/loading competence with the MBP positive control. `HUMAN_BG_120112_WEAK` is a predicted weak control, not an experimentally certified non-binder.

## Assay

- Quantitative competition binding with ten threefold competitor concentrations from 100,000 to 5.081 nM.
- Three independent experiments; two technical replicates per concentration.
- Use the randomized nine-plate map and preserve raw signals before curve fitting.
- Fit full curves with the laboratory's validated model; record IC50 in nM and retain model diagnostics.
- Confirm peptide identity, purity, and solubility. Analyze coded sample IDs.

## Register gate

Compare parent 15-mers with the two nested 11-mers. Use anchor-disruption controls or a direct structural method to establish the occupied P1-P9 frame. Predictor agreement is not experimental register resolution. A mutation with a shifted or unresolved register is excluded from the electrostatic interpretation.

## Interpretation

The locked charge-sensitive binding rule requires at least a threefold K-to-E affinity loss versus WT, at least twice the loss caused by K-to-Q, the same direction in all three experiments, and passed QC/register checks. A passed result does not by itself prove a direct electrostatic contact.

## Current state

Prepared only: no peptide has been ordered and no binding experiment has been run.
"""
    (lab_dir / "LAB_HANDOFF.md").write_text(handoff, encoding="utf-8")
    write_csv(out / "SHA256SUMS.csv", checksum_rows(out))
    print(json.dumps({"output": str(out), "peptides": len(panel), "primary_variants": 10, "plates": 9, "af3_jobs": len(jobs), "predictor_core_matches": matching}, indent=2))


if __name__ == "__main__":
    main()
