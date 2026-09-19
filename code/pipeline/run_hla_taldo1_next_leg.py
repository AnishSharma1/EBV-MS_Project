"""Build the dated TALDO1 HLA-specificity and charge-sensitivity package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

from hla_taldo1_next_leg import (
    generate_charge_controls,
    pair_charge_controls,
    parse_iedb_tsv,
    summarize_mutation_result,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "processed" / "taldo1_hla_next_leg_2026-09-04"
IEDB_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_METHOD = "recommended_binding"
CLAIM_BOUNDARY = (
    "Exploratory peptide-HLA binding-prediction sensitivity only; not evidence of "
    "physical repulsion, natural presentation, T-cell recognition, cross-reactivity, "
    "molecular mimicry, protection from MS, or an MS mechanism."
)

PEPTIDES = [
    {
        "peptide_id": "BALF5_627_641",
        "protein": "EBV BALF5",
        "coordinates_1_based": "627-641",
        "sequence": "TGGVYHFVKKHVHES",
    },
    {
        "peptide_id": "TALDO1_108_122",
        "protein": "human TALDO1",
        "coordinates_1_based": "108-122",
        "sequence": "DARLSFDKDAMVARA",
    },
    {
        "peptide_id": "TALDO1_216_230",
        "protein": "human TALDO1",
        "coordinates_1_based": "216-230",
        "sequence": "SVTKIYNYYKKFSYK",
    },
]

ALLELES = [
    "HLA-DRB1*13:03",
    "HLA-DRB1*15:01",
    "HLA-DRB5*01:01",
    "HLA-DQA1*01:02/DQB1*06:02",
]

PAIRS = [
    {
        "target_id": "HY13_SEQ_02",
        "native_hla": "HLA-DRB1*13:03",
        "ebv_peptide_id": "BALF5_627_641",
        "self_peptide_id": "TALDO1_108_122",
    },
    {
        "target_id": "HY15_SEQ_02",
        "native_hla": "HLA-DRB1*15:01",
        "ebv_peptide_id": "BALF5_627_641",
        "self_peptide_id": "TALDO1_216_230",
    },
]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def slug(value: str) -> str:
    return (
        value.lower()
        .replace("hla-", "")
        .replace("*", "_")
        .replace(":", "_")
        .replace("/", "__")
    )


def post_iedb(allele: str, submissions: list[dict[str, Any]]) -> str:
    fasta = "\n".join(f">{row['peptide_id']}\n{row['sequence']}" for row in submissions)
    body = urllib.parse.urlencode(
        {
            "method": IEDB_METHOD,
            "sequence_text": fasta,
            "allele": allele,
            "length": "asis",
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        IEDB_ENDPOINT,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8")


def verify_source_candidates() -> None:
    source = ROOT / "processed" / "high_yield_candidate_evidence_2026-08-28" / "stage1_assay_recommendations.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = {row["target_id"]: row for row in csv.DictReader(handle)}
    expected = {
        "HY13_SEQ_02": ("TGGVYHFVKKHVHES", "DARLSFDKDAMVARA", "stage1_medium_priority"),
        "HY15_SEQ_02": ("TGGVYHFVKKHVHES", "SVTKIYNYYKKFSYK", "stage1_medium_priority"),
    }
    for target_id, values in expected.items():
        row = rows.get(target_id)
        observed = (row["ebv_sequence"], row["self_sequence"], row["stage1_status"]) if row else None
        if observed != values:
            raise ValueError(f"source candidate drift for {target_id}: {observed}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(f"output already exists: {out}")
    out.mkdir(parents=True)
    raw_out = out / "raw_responses"
    raw_out.mkdir()
    prepared_out = out / "prepared_inputs"
    prepared_out.mkdir()

    verify_source_candidates()
    retrieved_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    baseline_records: list[dict[str, Any]] = []
    variants_by_allele: dict[str, list[dict[str, Any]]] = {}
    raw_manifest: list[dict[str, Any]] = []

    for allele in ALLELES:
        submissions = [
            {"seq_num": index, "peptide_id": row["peptide_id"], "sequence": row["sequence"]}
            for index, row in enumerate(PEPTIDES, start=1)
        ]
        raw = post_iedb(allele, submissions)
        raw_path = raw_out / f"baseline_{slug(allele)}.tsv"
        raw_path.write_text(raw, encoding="utf-8")
        records = parse_iedb_tsv(allele, submissions, raw)
        peptide_meta = {row["peptide_id"]: row for row in PEPTIDES}
        for record in records:
            baseline_records.append(
                {
                    **record,
                    "protein": peptide_meta[record["peptide_id"]]["protein"],
                    "coordinates_1_based": peptide_meta[record["peptide_id"]]["coordinates_1_based"],
                    "prediction_method_requested": IEDB_METHOD,
                    "prediction_endpoint": IEDB_ENDPOINT,
                    "retrieved_utc": retrieved_utc,
                    "raw_response_path": str(raw_path.relative_to(ROOT)),
                    "claim_boundary": CLAIM_BOUNDARY,
                }
            )
        variants: list[dict[str, Any]] = []
        for record in records:
            variants.extend(
                generate_charge_controls(
                    peptide_id=record["peptide_id"],
                    sequence=record["sequence"],
                    predicted_core=record["predicted_core"],
                )
            )
        for index, variant in enumerate(variants, start=1):
            variant["seq_num"] = index
        variants_by_allele[allele] = variants
        write_csv(prepared_out / f"charge_controls_{slug(allele)}.csv", variants)
        raw_manifest.append(
            {
                "allele": allele,
                "request_type": "baseline",
                "record_count": len(submissions),
                "raw_response_path": str(raw_path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                "retrieved_utc": retrieved_utc,
            }
        )

    baseline_lookup = {(row["allele"], row["peptide_id"]): row for row in baseline_records}
    mutation_records: list[dict[str, Any]] = []
    for allele in ALLELES:
        variants = variants_by_allele[allele]
        submissions = [
            {
                "seq_num": int(row["seq_num"]),
                "peptide_id": row["variant_id"],
                "sequence": row["mutant_sequence"],
            }
            for row in variants
        ]
        raw = post_iedb(allele, submissions)
        raw_path = raw_out / f"charge_controls_{slug(allele)}.tsv"
        raw_path.write_text(raw, encoding="utf-8")
        mutant_rows = parse_iedb_tsv(allele, submissions, raw)
        mutant_lookup = {row["peptide_id"]: row for row in mutant_rows}
        for variant in variants:
            record = summarize_mutation_result(
                allele=allele,
                variant=variant,
                baseline=baseline_lookup[(allele, variant["peptide_id"])],
                mutant=mutant_lookup[variant["variant_id"]],
            )
            record.update(
                {
                    "prediction_method_requested": IEDB_METHOD,
                    "prediction_endpoint": IEDB_ENDPOINT,
                    "retrieved_utc": retrieved_utc,
                    "raw_response_path": str(raw_path.relative_to(ROOT)),
                }
            )
            mutation_records.append(record)
        raw_manifest.append(
            {
                "allele": allele,
                "request_type": "charge_controls",
                "record_count": len(submissions),
                "raw_response_path": str(raw_path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                "retrieved_utc": retrieved_utc,
            }
        )

    pair_rows: list[dict[str, Any]] = []
    for pair in PAIRS:
        for allele in ALLELES:
            ebv = baseline_lookup[(allele, pair["ebv_peptide_id"])]
            self_row = baseline_lookup[(allele, pair["self_peptide_id"])]
            supported_count = sum(float(row["rank_percentile"]) <= 20.0 for row in (ebv, self_row))
            pair_rows.append(
                {
                    "target_id": pair["target_id"],
                    "native_hla": pair["native_hla"],
                    "tested_hla": allele,
                    "ebv_peptide_id": pair["ebv_peptide_id"],
                    "ebv_core": ebv["predicted_core"],
                    "ebv_rank_percentile": ebv["rank_percentile"],
                    "self_peptide_id": pair["self_peptide_id"],
                    "self_core": self_row["predicted_core"],
                    "self_rank_percentile": self_row["rank_percentile"],
                    "descriptive_rank_status": (
                        "both_arms_rank_le_20" if supported_count == 2
                        else "one_arm_rank_le_20" if supported_count == 1
                        else "neither_arm_rank_le_20"
                    ),
                    "threshold_source": "existing dossier binding-supported descriptive threshold",
                    "claim_boundary": CLAIM_BOUNDARY,
                }
            )

    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in mutation_records:
        grouped[(row["allele"], row["peptide_id"], row["perturbation"])].append(row)
    mutation_summary: list[dict[str, Any]] = []
    for (allele, peptide_id, perturbation), rows in sorted(grouped.items()):
        same_register = [row for row in rows if not row["register_shifted"]]
        mutation_summary.append(
            {
                "allele": allele,
                "peptide_id": peptide_id,
                "perturbation": perturbation,
                "variant_count": len(rows),
                "same_register_count": len(same_register),
                "register_shift_count": len(rows) - len(same_register),
                "median_delta_rank_same_register": (
                    statistics.median(row["delta_rank_percentile"] for row in same_register)
                    if same_register else ""
                ),
                "median_ic50_fold_same_register": (
                    statistics.median(row["ic50_fold_change"] for row in same_register)
                    if same_register else ""
                ),
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )

    write_csv(out / "baseline_hla_predictions.csv", baseline_records)
    write_csv(out / "candidate_pair_hla_comparison.csv", pair_rows)
    write_csv(out / "charge_sensitivity_results.csv", mutation_records)
    write_csv(out / "charge_sensitivity_summary.csv", mutation_summary)
    paired_charge_rows = pair_charge_controls(mutation_records)
    write_csv(out / "paired_charge_control_comparison.csv", paired_charge_rows)
    write_csv(out / "raw_response_manifest.csv", raw_manifest)

    class_i_feasibility = [
        {
            "allele": allele,
            "allele_class": "HLA class I",
            "ms_association": "reported protective association",
            "current_15mer_directly_comparable": False,
            "required_next_analysis": (
                "separate 8-11mer processing/presentation screen over EBV BALF5 and TALDO1; "
                "use class-I-specific controls and CD8 claim boundaries"
            ),
            "decision": "defer_from_current_HLA-II_binding_panel",
            "reason": "protective association does not make this allele a negative binding control",
        }
        for allele in ("HLA-B*38:01", "HLA-B*44:02")
    ]
    write_csv(out / "class_i_protective_hla_feasibility.csv", class_i_feasibility)

    protocol = {
        "analysis": "TALDO1 HLA specificity and formal-charge sensitivity",
        "alleles": ALLELES,
        "anchor_positions": [1, 4, 6, 9],
        "formal_charge_reversal": {"K": "E", "R": "E", "D": "K", "E": "K"},
        "neutralization_controls": {"K": "Q", "R": "Q", "D": "N", "E": "Q"},
        "histidine_rule": "not mutated because its charge is pH-dependent near physiological conditions",
        "mutation_rule": "one formally charged residue per variant within each allele-specific predicted 9-mer core",
        "register_rule": "if the returned core differs from the expected mutated core, classify as confounded_register_shift",
        "scoring_rule": "report raw rank, IC50, deltas, and fold changes without a new composite score",
        "prediction_method": IEDB_METHOD,
        "prediction_endpoint": IEDB_ENDPOINT,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "protocol_lock.json", protocol)

    sources = [
        {
            "topic": "DRB1*15:01 and DRB5*01:01 co-expression",
            "citation": "PMID:16111772",
            "url": "https://pubmed.ncbi.nlm.nih.gov/16111772/",
        },
        {
            "topic": "Hy.2E11 cross-reactivity across DRB5*01:01 and DRB1*15:01",
            "citation": "PMID:12244309; DOI:10.1038/ni835",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12244309/",
        },
        {
            "topic": "DR15 molecules jointly shape autoreactive repertoire",
            "citation": "PMCID:PMC7707104",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7707104/",
        },
        {
            "topic": "TALDO1 autoantigenicity and oligodendrocyte expression",
            "citation": "PMCID:PMC2191732; PMID:7964452",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2191732/",
        },
        {
            "topic": "TALDO1 antibody epitopes and viral peptide cross-reactivity",
            "citation": "PMID:10491006",
            "url": "https://pubmed.ncbi.nlm.nih.gov/10491006/",
        },
        {
            "topic": "HLA-A*02:01-restricted TALDO1 168-176 CD8 response",
            "citation": "PMID:16339578",
            "url": "https://pubmed.ncbi.nlm.nih.gov/16339578/",
        },
        {
            "topic": "HLA-B*38:01 and B*44:02 protective association",
            "citation": "PMID:26343388; DOI:10.1038/ng.3395",
            "url": "https://pubmed.ncbi.nlm.nih.gov/26343388/",
        },
    ]
    write_csv(out / "literature_source_manifest.csv", sources)

    brief = """# Experimental discussion brief

## What the computational perturbation can test

Each variant changes one formally charged residue in the allele-specific predicted HLA-II core. Every charge reversal has a same-position neutralization control. A larger binding change for reversal than neutralization, with the predicted register retained, would be consistent with a charge-sensitive contribution. It would not demonstrate physical repulsion. Similar effects from both variants would be consistent with a general side-chain or packing effect. A shifted predicted core makes the comparison inconclusive.

## Recommended experimental order

1. Measure binding of each wild-type BALF5 and TALDO1 peptide to the exact HLA molecule.
2. Use nested peptides or a direct register-mapping method to establish the occupied P1-P9 frame.
3. Test a compact set of paired charge-reversal and neutralization variants chosen from same-register computational results.
4. Only after binding and register are established, consider pMHC tetramers or functional T-cell assays for recognition and cross-reactivity.

DRB1*15:01 and DRB5*01:01 are both expressed on the DR15 haplotype, so a comparison tests restriction and presentation behavior but cannot genetically separate the effects of their co-inheritance. DQA1*01:02/DQB1*06:02 is a linked, secondary HLA-II context. HLA-B*38:01 and B*44:02 require a separate class-I/CD8 workflow.
"""
    (out / "EXPERIMENTAL_DISCUSSION_BRIEF.md").write_text(brief, encoding="utf-8")

    checklist = """# Manuscript claim and evidence checklist

| Proposed statement | Evidence status | Allowed wording |
|---|---|---|
| TALDO1 is biologically relevant to MS research | Published protein-level evidence | Prior studies reported TALDO1 expression in oligodendrocytes and immune reactivity in subsets of people with MS. |
| Olivia Thomas independently mentioned TALDO1 | Oral, unpublished context | Describe as mentor-reported protein-level convergence unless a citable dataset or exact peptide is supplied. |
| HY13 TALDO1 108-122 matches prior work | Partial regional connection | It overlaps residues 108-115 of a reported 101-115 antibody-reactive region; this is not exact epitope replication. |
| HY15 TALDO1 216-230 matches prior work | Regional proximity only | It lies immediately before a reported 231-245 antibody-reactive region; do not call it an overlap. |
| Either candidate matches the published TALDO1 168-176 CD8 epitope | Unsupported | State that the published HLA-A*02:01-restricted 168-176 epitope is distinct from both candidates. |
| BALF5-TALDO1 is a confirmed molecular mimic | Unsupported | Call both rows Stage-1 computational experimental leads pending exact-HLA binding and register mapping. |
| Charge reversal proves peptide-HLA repulsion | Unsupported | Call it predictor sensitivity; physical interaction requires experimental or validated energetic evidence. |
| DRB5 analysis removes dual-inheritance confounding | Unsupported | Say it tests HLA-specific presentation behavior within the co-inherited DR15 context. |
| B*38:01/B*44:02 are negative controls | Unsupported | Describe them as protective-association alleles requiring a separate HLA-I hypothesis. |
"""
    (out / "MANUSCRIPT_CLAIM_CHECKLIST.md").write_text(checklist, encoding="utf-8")

    pair_lookup = {(row["target_id"], row["tested_hla"]): row for row in pair_rows}
    readme = f"""# TALDO1 HLA next leg

This package extends the two Stage-1 BALF5-TALDO1 leads across DRB1*13:03, DRB1*15:01, DRB5*01:01, and DQA1*01:02/DQB1*06:02 using one live IEDB `recommended_binding` retrieval dated {retrieved_utc}. It adds single-site charge-reversal and neutralization sensitivity controls without changing the earlier candidate rankings.

## Baseline result

| Candidate | Tested HLA | BALF5 rank | TALDO1 rank | Descriptive status |
|---|---|---:|---:|---|
"""
    for pair in PAIRS:
        for allele in ALLELES:
            row = pair_lookup[(pair["target_id"], allele)]
            readme += (
                f"| {pair['target_id']} | {allele} | {float(row['ebv_rank_percentile']):.1f} | "
                f"{float(row['self_rank_percentile']):.1f} | {row['descriptive_rank_status']} |\n"
            )
    readme += f"""

The shared BALF5 peptide and TALDO1 216-230 both fall at or below the existing descriptive rank-20 threshold on DRB5*01:01. Neither DQA1*01:02/DQB1*06:02 candidate pair has both arms at or below that threshold in this run. These are single-predictor hypotheses and do not establish presentation.

## How to read the files

- `candidate_pair_hla_comparison.csv`: compact HLA comparison for the two candidate pairs.
- `charge_sensitivity_results.csv`: every mutation, raw prediction, register check, and claim boundary.
- `charge_sensitivity_summary.csv`: descriptive summaries only; no new composite score.
- `paired_charge_control_comparison.csv`: same-position reversal versus neutralization contrasts.
- `EXPERIMENTAL_DISCUSSION_BRIEF.md`: proposed experimental sequence and interpretation logic.
- `MANUSCRIPT_CLAIM_CHECKLIST.md`: language that is supported, suggestive, or unsupported.
- `class_i_protective_hla_feasibility.csv`: why B*38:01/B*44:02 remain a separate class-I branch.
- `raw_responses/`, `raw_response_manifest.csv`, and `protocol_lock.json`: exact provenance and locked rules.

## Claim boundary

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")

    checksum_rows = []
    for path in sorted(item for item in out.rglob("*") if item.is_file()):
        if path.name == "SHA256SUMS.csv":
            continue
        checksum_rows.append(
            {
                "relative_path": str(path.relative_to(out)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "bytes": path.stat().st_size,
            }
        )
    write_csv(out / "SHA256SUMS.csv", checksum_rows)
    print(f"wrote {len(checksum_rows) + 1} files to {out}")


if __name__ == "__main__":
    main()
