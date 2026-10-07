#!/usr/bin/env python3
"""Retrieve the frozen candidate arms from public experimental-evidence sources.

The IEDB IQ-API is queried literally by peptide and requested HLA allele.
Its response is saved in full so the compact census can be regenerated or
audited.  Other sources are represented as explicit, source-appropriate
review tasks: immunopeptidomics supports natural-presentation context but
does not replace the locked exact-allele MHC/ligand gate.

This script intentionally does not change candidate ranks.
"""

from __future__ import annotations

import csv
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/revision_round_4/threshold_sensitivity_2026-09-20/threshold_grid_pair_results.csv"
OUT = ROOT / "output/revision_round_4/cross_source_public_evidence_2026-09-20"
API = "https://query-api.iedb.org/mhc_search"
CHECKED_AT = "2026-09-20"
SELECT = (
    "elution_id,structure_id,linear_sequence,mhc_allele_name,qualitative_measure,"
    "quantitative_measure,assay_description,pubmed_id,reference_authors,"
    "reference_titles,reference_dates,journal_name,mhc_allele_resolution"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def query_iedb(peptide: str, allele: str) -> tuple[list[dict[str, Any]], str]:
    params = {
        "linear_sequence": f"eq.{peptide}",
        "mhc_allele_name": f"eq.{allele}",
        "select": SELECT,
        "limit": "100",
    }
    url = f"{API}?{urlencode(params)}"
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "EBV-MS-public-evidence-audit/1.0"})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=45) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if not isinstance(payload, list):
                raise ValueError("IEDB response was not a JSON list")
            return payload, url
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def is_positive(record: dict[str, Any]) -> bool:
    return str(record.get("qualitative_measure", "")).strip().lower() == "positive"


def source_tasks(peptide: str, arm: str) -> list[dict[str, str]]:
    # The HLA Ligand Atlas is a normal-tissue resource, so EBV peptides are
    # deliberately not queried there.  SysteMHC and PeptideAtlas remain useful
    # for checking naturally presented peptides, not purified binding.
    tasks = [
        {
            "source": "SysteMHC Atlas",
            "source_url": "https://systemhc.sjtu.edu.cn/explore",
            "evidence_domain": "natural_mhc_presentation",
            "review_scope": "exact_peptide; retain sample HLA typing and allele-assignment ambiguity",
        },
        {
            "source": "Human HLA PeptideAtlas",
            "source_url": "https://peptideatlas.org/builds/human/hla/",
            "evidence_domain": "natural_mhc_presentation",
            "review_scope": "exact_peptide; retain dataset and HLA typing when available",
        },
    ]
    if arm == "self":
        tasks.insert(0, {
            "source": "HLA Ligand Atlas",
            "source_url": "https://hla-ligand-atlas.org/welcome",
            "evidence_domain": "normal_tissue_mhc_presentation",
            "review_scope": "exact_peptide; retain tissue, donor HLA typing, and allele-assignment ambiguity",
        })
    return [dict(task, peptide=peptide, arm=arm) for task in tasks]


def literature_task(peptide: str, allele: str, arm: str) -> dict[str, str]:
    """Build a reviewable PubMed query without treating search hits as assays."""
    allele_short = allele.replace("HLA-", "")
    aliases = f'"{allele}" OR "{allele_short}" OR "HLA-DR" OR MHC'
    term = f'"{peptide}"[Title/Abstract] AND ({aliases})'
    return {
        "source": "PubMed primary-literature recovery",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/",
        "peptide": peptide,
        "requested_hla": allele,
        "arm": arm,
        "pubmed_query": term,
        "evidence_domain": "primary_experimental_assay_review",
        "review_scope": "retain only primary experimental peptide-HLA assays; record PMID, assay type, exact allele, outcome, and value/unit",
        "review_status": "pending_targeted_review",
        "gate_policy": "qualifies_only_if_positive_exact_sequence_exact_hla_mhc_ligand_assay",
    }


def main() -> None:
    candidates: dict[str, dict[str, str]] = {}
    for row in read_rows(SOURCE):
        candidates.setdefault(row["candidate_id"], row)
    if len(candidates) != 49:
        raise ValueError(f"Expected 49 candidates, found {len(candidates)}")

    OUT.mkdir(parents=True, exist_ok=True)
    raw_path = OUT / "iedb_iq_api_raw_records.jsonl"
    arms: list[dict[str, str]] = []
    source_review: list[dict[str, str]] = []
    literature_review: list[dict[str, str]] = []
    with raw_path.open("w", encoding="utf-8") as raw:
        for candidate_id in sorted(candidates, key=lambda value: int(value[2:])):
            candidate = candidates[candidate_id]
            for arm in ("ebv", "self"):
                peptide = candidate[f"{arm}_sequence"]
                allele = candidate["allele"]
                records, query_url = query_iedb(peptide, allele)
                raw.write(json.dumps({
                    "candidate_id": candidate_id,
                    "arm": arm,
                    "peptide": peptide,
                    "requested_hla": allele,
                    "query_url": query_url,
                    "records": records,
                }, sort_keys=True) + "\n")
                positives = [record for record in records if is_positive(record)]
                assay_ids = ";".join(str(record.get("elution_id", "")) for record in positives)
                pmids = ";".join(str(record.get("pubmed_id", "")) for record in positives if record.get("pubmed_id"))
                arms.append({
                    "candidate_id": candidate_id,
                    "arm": arm,
                    "peptide": peptide,
                    "requested_hla": allele,
                    "retrieval_source": "IEDB IQ-API",
                    "retrieval_status": "evaluated",
                    "evidence_domain": "mhc_ligand_or_binding",
                    "allele_resolution": "exact_literal_query",
                    "returned_record_count": str(len(records)),
                    "positive_exact_hla_record_count": str(len(positives)),
                    "qualifying_iedb_mhc_assay_ids": assay_ids,
                    "qualifying_pubmed_ids": pmids,
                    "qualifies_locked_mhc_gate": str(bool(positives)),
                    "decision_effect": "qualifies_public_evidence_gate" if positives else "reviewed_no_qualifying_record_no_penalty",
                    "query_url": query_url,
                    "checked_at": CHECKED_AT,
                })
                source_review.extend(source_tasks(peptide, arm))
                literature_review.append(literature_task(peptide, allele, arm))
                # Be polite to a public database while retaining deterministic output.
                time.sleep(0.08)

    arm_fields = list(arms[0])
    write_rows(OUT / "iedb_iq_api_arm_census.csv", arms, arm_fields)
    review_fields = ["source", "source_url", "evidence_domain", "review_scope", "peptide", "arm"]
    deduplicated = {(r["source"], r["peptide"], r["arm"]): r for r in source_review}
    review_rows = sorted(deduplicated.values(), key=lambda r: (r["source"], r["peptide"], r["arm"]))
    for row in review_rows:
        row["review_status"] = "pending_targeted_review"
        row["gate_policy"] = "context_only_never_substitutes_for_exact_hla_mhc_gate"
    review_fields += ["review_status", "gate_policy"]
    write_rows(OUT / "presentation_source_review_manifest.csv", review_rows, review_fields)
    literature_fields = [
        "source", "source_url", "peptide", "requested_hla", "arm", "pubmed_query",
        "evidence_domain", "review_scope", "review_status", "gate_policy",
    ]
    literature_deduplicated = {
        (r["peptide"], r["requested_hla"], r["arm"]): r for r in literature_review
    }
    write_rows(
        OUT / "literature_recovery_manifest.csv",
        sorted(literature_deduplicated.values(), key=lambda r: (r["peptide"], r["requested_hla"], r["arm"])),
        literature_fields,
    )

    counts = Counter(row["decision_effect"] for row in arms)
    (OUT / "status.json").write_text(json.dumps({
        "status": "complete_iedb_iq_api_retrieval; presentation_review_manifest_created",
        "candidate_count": 49,
        "arm_count": len(arms),
        "iedb_iq_api_evaluated_arm_count": len(arms),
        "iedb_positive_exact_hla_mhc_arm_count": sum(r["qualifies_locked_mhc_gate"] == "True" for r in arms),
        "decision_counts": dict(sorted(counts.items())),
        "rank_changes": "none",
        "gate_rule": "Only positive exact-sequence/exact-allele MHC/ligand evidence qualifies.",
        "absence_policy": "A reviewed no-record result is not evidence of absence and carries no ranking penalty.",
        "presentation_policy": "Atlas evidence is recorded as presentation context and cannot substitute for the locked binding gate.",
        "literature_recovery_unique_peptide_allele_arm_count": len(literature_deduplicated),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# Cross-source public-evidence recovery\n\n"
        "`iedb_iq_api_arm_census.csv` is a literal exact-peptide/exact-allele retrieval against IEDB's public Query API. "
        "The full returned records are preserved in `iedb_iq_api_raw_records.jsonl`.\n\n"
        "`literature_recovery_manifest.csv` provides one exact peptide-and-HLA primary-literature query per unique frozen peptide/allele/arm combination. Search hits must be reviewed down to an experimental assay before they can affect the gate.\n\n"
        "`presentation_source_review_manifest.csv` defines the next targeted review across HLA Ligand Atlas (self peptides only), "
        "SysteMHC Atlas, and Human HLA PeptideAtlas. These sources address natural presentation and are deliberately kept separate "
        "from the locked exact-allele binding/ligand gate. No candidate rank is changed by this recovery.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
