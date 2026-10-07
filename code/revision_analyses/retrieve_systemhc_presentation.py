#!/usr/bin/env python3
"""Exact-peptide presentation-context census in the public SysteMHC Atlas.

SysteMHC records MHC-associated peptides found by mass spectrometry.  A hit is
not converted to binding evidence; allele assignment is retained as reported
because many samples are polyallelic.
"""

from __future__ import annotations

import csv
import html
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/revision_round_4/cross_source_public_evidence_2026-09-20/presentation_source_review_manifest.csv"
OUT = ROOT / "output/revision_round_4/cross_source_public_evidence_2026-09-20"
URL = "https://systemhc.sjtu.edu.cn/explore"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def query(peptide: str) -> tuple[str, str, str]:
    request_url = f"{URL}?{urlencode({'search_hit': peptide})}"
    request = Request(request_url, headers={"User-Agent": "EBV-MS-public-evidence-audit/1.0"})
    try:
        with urlopen(request, timeout=60) as response:
            return peptide, request_url, response.read().decode("utf-8", errors="replace")
    except Exception as exc:
        return peptide, request_url, json.dumps({"error": str(exc)})


def rows_from_html(peptide: str, page: str) -> list[dict[str, str]]:
    if page.startswith('{"error":'):
        return []
    body = re.search(r'<tbody>(.*?)</tbody>', page, flags=re.S | re.I)
    if not body:
        return []
    parsed: list[dict[str, str]] = []
    for row_html in re.findall(r'<tr[^>]*>(.*?)</tr>', body.group(1), flags=re.S | re.I):
        cells = [html.unescape(re.sub(r'<[^>]+>', '', cell)).strip() for cell in re.findall(r'<td[^>]*>(.*?)</td>', row_html, flags=re.S | re.I)]
        # The public results table orders these fields as SysteMHC ID, SampleID,
        # protein, peptide, modified peptide, PTM, position, top allele, class,
        # spectra count, probability, annotation score.
        if len(cells) >= 12 and cells[3].upper() == peptide.upper():
            parsed.append({
                "systemhc_id": cells[0], "sample_id": cells[1], "protein_id": cells[2],
                "peptide": cells[3], "top_allele": cells[7], "mhc_class": cells[8],
                "spectra_count": cells[9], "probability": cells[10], "annotation_score": cells[11],
            })
    return parsed


def main() -> None:
    manifest = [row for row in read_rows(SOURCE) if row["source"] == "SysteMHC Atlas"]
    peptides = sorted({row["peptide"] for row in manifest})
    raw: dict[str, dict[str, str]] = {}
    # Two concurrent read-only requests keeps pressure low on a public atlas.
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(query, peptide): peptide for peptide in peptides}
        for future in as_completed(futures):
            peptide, request_url, page = future.result()
            raw[peptide] = {"query_url": request_url, "page": page}
            time.sleep(0.15)
    (OUT / "systemhc_exact_peptide_raw.json").write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    evidence = []
    for row in manifest:
        results = rows_from_html(row["peptide"], raw[row["peptide"]]["page"])
        if raw[row["peptide"]]["page"].startswith('{"error":'):
            status = "source_access_failed_no_inference"
        elif results:
            status = "exact_peptide_presentation_records_found_context_only"
        else:
            status = "no_exact_peptide_record_in_systemhc_no_penalty"
        if not results:
            results = [{"systemhc_id": "", "sample_id": "", "protein_id": "", "peptide": row["peptide"], "top_allele": "", "mhc_class": "", "spectra_count": "", "probability": "", "annotation_score": ""}]
        for result in results:
            evidence.append({
                "source": "SysteMHC Atlas", "arm": row["arm"], "peptide": row["peptide"],
                "retrieval_status": status, "evidence_domain": "natural_mhc_presentation_context_only",
                "allele_assignment": result["top_allele"] or "not_applicable", "mhc_class": result["mhc_class"],
                "systemhc_id": result["systemhc_id"], "sample_id": result["sample_id"], "protein_id": result["protein_id"],
                "spectra_count": result["spectra_count"], "probability": result["probability"],
                "annotation_score": result["annotation_score"], "qualifies_locked_mhc_gate": "False",
                "decision_effect": "context_only_no_ranking_change", "query_url": raw[row["peptide"]]["query_url"],
            })
    fields = list(evidence[0])
    with (OUT / "systemhc_exact_peptide_presentation_census.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(evidence)


if __name__ == "__main__":
    main()
