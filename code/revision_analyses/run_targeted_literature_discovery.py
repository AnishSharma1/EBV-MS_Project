#!/usr/bin/env python3
"""Run the predeclared peptide-HLA PubMed discovery queries.

This is a discovery layer only.  A PubMed hit is never treated as an assay
result until its primary paper is inspected for peptide, exact HLA, assay
type, and outcome.
"""

from __future__ import annotations

import csv
import json
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "output/revision_round_4/cross_source_public_evidence_2026-09-20/literature_recovery_manifest.csv"
OUT = ROOT / "output/revision_round_4/cross_source_public_evidence_2026-09-20"
ENTREZ = Path("/Users/anishsharma/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/ncbi-entrez-skill/scripts/ncbi_entrez.py")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    rows = read_rows(INPUT)
    result_rows: list[dict[str, str]] = []
    raw_path = OUT / "pubmed_targeted_discovery_raw.jsonl"
    with raw_path.open("w", encoding="utf-8") as raw:
        for row in rows:
            # Older epitope papers commonly omit HLA names from abstracts.  Use
            # exact sequence alone for discovery, while retaining the original
            # peptide-HLA query in the output as provenance for the stricter pass.
            sequence_only_query = f'"{row["peptide"]}"[Title/Abstract]'
            request = {
                "endpoint": "esearch",
                "params": {"db": "pubmed", "term": sequence_only_query, "retmode": "json", "retmax": 20},
                "max_items": 20,
                "max_depth": 3,
                "timeout_sec": 30,
            }
            completed = subprocess.run(
                ["python3", str(ENTREZ)], input=json.dumps(request), text=True,
                capture_output=True, check=False,
            )
            try:
                response = json.loads(completed.stdout)
            except json.JSONDecodeError:
                response = {"ok": False, "error": {"message": completed.stderr.strip() or "invalid helper output"}}
            ids = []
            if response.get("ok"):
                summary = response.get("summary", {})
                if isinstance(summary, dict):
                    ids = summary.get("idlist", [])
            raw.write(json.dumps({"input": row, "response": response}, sort_keys=True) + "\n")
            result_rows.append({
                "peptide": row["peptide"],
                "requested_hla": row["requested_hla"],
                "arm": row["arm"],
                "pubmed_query": row["pubmed_query"],
                "pubmed_sequence_only_query": sequence_only_query,
                "discovery_status": "candidate_citations_require_primary_assay_review" if ids else "no_title_abstract_candidate_citation_no_penalty",
                "pubmed_candidate_count": str(len(ids)),
                "pubmed_candidate_ids": ";".join(str(value) for value in ids),
            })
            time.sleep(0.34)
    fields = list(result_rows[0])
    with (OUT / "pubmed_targeted_discovery.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(result_rows)


if __name__ == "__main__":
    main()
