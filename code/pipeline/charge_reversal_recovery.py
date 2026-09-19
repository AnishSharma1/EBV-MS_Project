"""Deterministic preparation and fail-closed logic for charge-reversal recovery."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

from build_tcell_library_v2 import DRA_SEQUENCE, DRB1_1501_SEQUENCE
from charge_reversal_binding_pilot import primary_panel


CLAIM_BOUNDARY = (
    "Method-consensus structural sensitivity only; not measured affinity, direct "
    "electrostatic-contact proof, presentation, TCR recognition, molecular mimicry, or MS mechanism."
)
AF3_CONDITIONS = (
    ("peptide_templates_on_seed_104729", True, 104729),
    ("peptide_templates_off_seed_314159", False, 314159),
    ("peptide_templates_off_seed_104729", False, 104729),
)
PANDORA_TEMPLATES = ("1BX2", "6CQQ")
REGISTER_STARTS = {"BALF5": (3, 4, 5), "TALDO1": (4, 5, 6)}
EXPECTED_START = {"BALF5": 4, "TALDO1": 5}
CONSENSUS_CLASSES = {
    "robust_structural_contact_hypothesis",
    "no_structural_charge_contact_support",
    "method_or_template_dependent",
    "register_confounded",
    "not_evaluable",
}


def primary_ten() -> list[dict[str, Any]]:
    return [row for row in primary_panel() if row["primary_mutation_panel"]]


def make_af3_ablation_jobs() -> list[dict[str, Any]]:
    """Build the exact 30-job external AlphaFold Server upload batch."""
    jobs: list[dict[str, Any]] = []
    for label, peptide_templates, seed in AF3_CONDITIONS:
        suffix = f"pt_{'on' if peptide_templates else 'off'}_s{seed}"
        for row in primary_ten():
            jobs.append({
                "name": f"crr_{row['sample_id'].lower()}_{suffix}",
                "modelSeeds": [seed],
                "sequences": [
                    {"proteinChain": {"sequence": DRA_SEQUENCE, "count": 1, "useStructureTemplate": True}},
                    {"proteinChain": {"sequence": DRB1_1501_SEQUENCE, "count": 1, "useStructureTemplate": True}},
                    {"proteinChain": {"sequence": row["sequence"], "count": 1, "useStructureTemplate": peptide_templates}},
                ],
                "dialect": "alphafoldserver",
                "version": 3,
                "recoveryCondition": label,
            })
    validate_af3_ablation_jobs(jobs)
    return jobs


def validate_af3_ablation_jobs(jobs: list[dict[str, Any]]) -> None:
    if len(jobs) != 30 or len({job["name"] for job in jobs}) != 30:
        raise ValueError("AlphaFold recovery batch must have 30 unique jobs")
    expected_sequences = {row["sequence"] for row in primary_ten()}
    counts: dict[tuple[bool, int], int] = {}
    for job in jobs:
        chains = [entry["proteinChain"] for entry in job["sequences"]]
        if len(chains) != 3 or chains[0]["sequence"] != DRA_SEQUENCE or chains[1]["sequence"] != DRB1_1501_SEQUENCE:
            raise ValueError(f"chain order or HLA sequence mismatch: {job['name']}")
        if chains[2]["sequence"] not in expected_sequences or len(chains[2]["sequence"]) != 15:
            raise ValueError(f"peptide mismatch or truncation: {job['name']}")
        if not chains[0]["useStructureTemplate"] or not chains[1]["useStructureTemplate"]:
            raise ValueError(f"HLA templates must remain enabled: {job['name']}")
        key = (bool(chains[2]["useStructureTemplate"]), int(job["modelSeeds"][0]))
        counts[key] = counts.get(key, 0) + 1
    if counts != {(True, 104729): 10, (False, 314159): 10, (False, 104729): 10}:
        raise ValueError(f"incorrect template/seed factorial: {counts}")


def pandora_candidate_runs() -> list[dict[str, Any]]:
    """Return 60 runs requesting 1,200 candidate models in total."""
    rows: list[dict[str, Any]] = []
    for peptide in primary_ten():
        arm = str(peptide["arm"])
        for start in REGISTER_STARTS[arm]:
            register = "expected" if start == EXPECTED_START[arm] else ("minus_1" if start < EXPECTED_START[arm] else "plus_1")
            anchors = [start, start + 3, start + 5, start + 8]
            for template in PANDORA_TEMPLATES:
                rows.append({
                    "run_id": f"{peptide['sample_id']}__start{start}__{template}",
                    "sample_id": peptide["sample_id"],
                    "arm": arm,
                    "peptide": peptide["sequence"],
                    "mhc_class": "II",
                    "hla_alpha": "HLA-DRA*01:01",
                    "hla_beta": "HLA-DRB1*15:01",
                    "core_start_1_based": start,
                    "core_9mer": peptide["sequence"][start - 1:start + 8],
                    "register_label": register,
                    "anchors_1_based": ",".join(map(str, anchors)),
                    "forced_template_pdb": template,
                    "requested_models": 20,
                    "status": "prepared_not_run",
                    "claim_boundary": CLAIM_BOUNDARY,
                })
    if len(rows) != 60 or sum(row["requested_models"] for row in rows) != 1200:
        raise AssertionError("PANDORA candidate manifest cardinality failure")
    return rows


def pandora_calibration_runs() -> list[dict[str, Any]]:
    """Prepare reciprocal leave-one-template-out expected and adjacent-register controls."""
    targets = (
        ("1BX2", "6CQQ", "ENPVVHFFKNIVTP", 5, "1BX2 experimental complex rebuilt from 6CQQ only"),
        ("6CQQ", "1BX2", "RFYKTLRAEQASQ", 3, "6CQQ experimental complex rebuilt from 1BX2 only"),
    )
    rows: list[dict[str, Any]] = []
    for target, allowed_template, peptide, expected_start, note in targets:
        for delta, label in ((-1, "minus_1"), (0, "expected"), (1, "plus_1")):
            start = expected_start + delta
            rows.append({
                "run_id": f"calibration_{target}__start{start}__template_{allowed_template}",
                "experimental_target_pdb": target,
                "forced_template_pdb": allowed_template,
                "excluded_template_pdb": target,
                "peptide": peptide,
                "core_start_1_based": start,
                "register_label": label,
                "anchors_1_based": ",".join(map(str, (start, start + 3, start + 5, start + 8))),
                "requested_models": 20,
                "top_five_core_backbone_rmsd_threshold_A": 2.0,
                "status": "prepared_not_run",
                "leakage_rule": "target PDB must be absent from template and alignment inputs",
                "note": note,
                "claim_boundary": CLAIM_BOUNDARY,
            })
    return rows


def _template_pdb_id(cif_path: Path) -> str:
    head = cif_path.read_text(encoding="utf-8", errors="replace")[:4096]
    match = re.search(r"^_entry\.id\s+(\S+)", head, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"missing template PDB ID: {cif_path}")
    return match.group(1).strip("'\"").upper()


def audit_baseline_templates(return_dir: Path, job_summary_csv: Path) -> list[dict[str, Any]]:
    """Extract the actual peptide-chain template identities and mappings per frozen job."""
    with job_summary_csv.open(newline="", encoding="utf-8") as handle:
        summaries = {row["job_name"]: row for row in csv.DictReader(handle)}
    rows: list[dict[str, Any]] = []
    for job_dir in sorted(path for path in return_dir.iterdir() if path.is_dir()):
        maps = list((job_dir / "templates").glob("*chains_c_query_to_hit.json"))
        if len(maps) != 1:
            raise ValueError(f"expected one peptide template map for {job_dir.name}")
        mappings = json.loads(maps[0].read_text(encoding="utf-8"))
        hit_rows = []
        for hit in mappings:
            cif_path = job_dir / "templates" / hit["name"]
            query_indices = sorted(int(value) for value in hit["queryIndices"])
            hit_rows.append({
                "pdb": _template_pdb_id(cif_path),
                "coverage": len(set(query_indices)) / 15.0,
                "query": ",".join(str(value + 1) for value in query_indices),
                "template": ",".join(str(int(value) + 1) for value in hit["templateIndices"]),
            })
        summary = summaries[job_dir.name]
        rows.append({
            "job_name": job_dir.name,
            "sample_id": summary["sample_id"],
            "peptide_template_count": len(hit_rows),
            "peptide_template_pdb_ids": ";".join(hit["pdb"] for hit in hit_rows),
            "peptide_template_query_coverages": ";".join(f"{hit['coverage']:.3f}" for hit in hit_rows),
            "peptide_template_query_positions_1_based": " | ".join(hit["query"] for hit in hit_rows),
            "peptide_template_positions_1_based": " | ".join(hit["template"] for hit in hit_rows),
            "median_peptide_chain_iptm": summary["median_peptide_chain_iptm"],
            "median_peptide_mean_plddt": summary["median_peptide_mean_plddt"],
            "within_job_pairwise_peptide_rmsd_median_A": summary["within_job_pairwise_peptide_rmsd_median_A"],
            "diagnostic_scope": "template-outcome association only; matched ablation is required for causal attribution",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    if len(rows) != 10:
        raise ValueError("baseline audit must cover exactly ten jobs")
    return rows


def classify_consensus(evidence: dict[str, Any]) -> str:
    """Apply the predeclared consensus rules with fail-closed precedence."""
    required = (
        "complete", "calibration_pass", "sequence_qc_pass", "structural_qc_pass",
        "expected_register_preferred", "expected_register_core_rmsd_max_A",
        "af3_seeds_agree", "pandora_templates_agree", "methods_agree",
        "af3_contact_frequency", "pandora_contact_frequency",
    )
    if any(key not in evidence for key in required):
        return "not_evaluable"
    if not all(bool(evidence[key]) for key in ("complete", "calibration_pass", "sequence_qc_pass", "structural_qc_pass")):
        return "not_evaluable"
    if not bool(evidence["expected_register_preferred"]) or float(evidence["expected_register_core_rmsd_max_A"]) > 2.0:
        return "register_confounded"
    if not all(bool(evidence[key]) for key in ("af3_seeds_agree", "pandora_templates_agree", "methods_agree")):
        return "method_or_template_dependent"
    if float(evidence["af3_contact_frequency"]) >= 0.8 and float(evidence["pandora_contact_frequency"]) >= 0.8:
        return "robust_structural_contact_hypothesis"
    return "no_structural_charge_contact_support"


def empty_consensus_rows() -> list[dict[str, Any]]:
    return [{
        "mutation_sample_id": row["sample_id"],
        "classification": "not_evaluable",
        "reason": "new AlphaFold and PANDORA ensembles are not yet complete",
        "af3_models_expected": 15,
        "pandora_models_expected": 120,
        "models_are_biological_replicates": False,
        "claim_boundary": CLAIM_BOUNDARY,
    } for row in primary_ten() if row["reagent_role"] != "wild_type"]
