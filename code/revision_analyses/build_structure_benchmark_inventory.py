#!/usr/bin/env python3
"""Build a transparent HLA-II solved-structure benchmark inventory from RCSB metadata."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/revision_round_3/expansion_sources_2026-09-15"
OUT = ROOT / "output/revision_round_3/solved_structure_benchmark_inventory_2026-09-15"
BOUNDARY = (
    "RCSB metadata inventory only. Candidate peptide chains are length-filtered, not yet "
    "manually verified or assigned to P1-P9 pockets; no predictor accuracy is claimed."
)


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_sequence(value):
    return "".join((value or "").split()).upper()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    page1 = json.loads((SOURCE / "rcsb_hla_dra_search.json").read_text())
    page2 = json.loads((SOURCE / "rcsb_hla_dra_search_page2.json").read_text())
    ids = [row["identifier"] for row in page1["result_set"] + page2["result_set"]]
    if len(ids) != 144 or len(set(ids)) != 144:
        raise ValueError(f"expected 144 unique entries, found {len(ids)} / {len(set(ids))}")
    metadata = json.loads((SOURCE / "rcsb_hla_dra_entry_metadata.json").read_text())["data"]["entries"]
    by_id = {entry["rcsb_id"]: entry for entry in metadata}
    if set(ids) != set(by_id):
        raise ValueError("search and metadata entry sets differ")

    entries = []
    peptides = []
    for pdb_id in ids:
        entry = by_id[pdb_id]
        title = (entry.get("struct") or {}).get("title", "")
        methods = ";".join(row.get("method", "") for row in (entry.get("exptl") or []))
        resolution_values = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        resolution = min(resolution_values) if resolution_values else ""
        descriptions = []
        candidate_count = 0
        for entity in entry.get("polymer_entities") or []:
            entity_poly = entity.get("entity_poly") or {}
            seq = clean_sequence(entity_poly.get("pdbx_seq_one_letter_code_can"))
            description = (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description", "") or ""
            descriptions.append(description)
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            refs = identifiers.get("reference_sequence_identifiers") or []
            uniprot = sorted({ref.get("database_accession", "") for ref in refs if ref.get("database_name") == "UniProt"})
            if 9 <= len(seq) <= 30 and "class ii" not in description.lower() and "hla" not in description.lower():
                candidate_count += 1
                peptides.append({
                    "pdb_id": pdb_id,
                    "candidate_peptide_sequence": seq,
                    "length": len(seq),
                    "description": description,
                    "asym_ids": ";".join(identifiers.get("asym_ids") or []),
                    "auth_asym_ids": ";".join(identifiers.get("auth_asym_ids") or []),
                    "uniprot_accessions": ";".join(uniprot),
                    "manual_peptide_identity_status": "pending_manual_verification",
                    "experimental_register_status": "pending_residue_to_pocket_annotation",
                    "claim_boundary": BOUNDARY,
                })
        joined = " | ".join(descriptions).lower()
        has_tcr = "t-cell receptor" in joined or "t cell receptor" in joined or "tcr" in joined
        resolution_ok = bool(resolution != "" and float(resolution) <= 4.0)
        entries.append({
            "pdb_id": pdb_id,
            "title": title,
            "experimental_method": methods,
            "resolution_angstrom": resolution,
            "candidate_peptide_chain_count": candidate_count,
            "contains_tcr_description": has_tcr,
            "initial_metadata_screen": "candidate" if candidate_count == 1 and resolution_ok else "manual_review",
            "coordinate_download_status": "prepared_not_downloaded",
            "manual_register_annotation_status": "pending",
            "claim_boundary": BOUNDARY,
        })

    write_csv(OUT / "rcsb_hla_dra_entry_inventory.csv", entries)
    write_csv(OUT / "candidate_peptide_chains.csv", peptides)
    summary = {
        "status": "inventory_complete_benchmark_not_run",
        "query_definition": "experimental RCSB entries with a polymer entity cross-referenced to UniProt P01903 (HLA-DRA)",
        "entry_count": len(entries),
        "candidate_peptide_chain_rows": len(peptides),
        "entries_with_exactly_one_candidate_peptide_chain": sum(row["candidate_peptide_chain_count"] == 1 for row in entries),
        "entries_passing_initial_metadata_screen": sum(row["initial_metadata_screen"] == "candidate" for row in entries),
        "entries_with_tcr_description": sum(row["contains_tcr_description"] for row in entries),
        "required_before_accuracy_estimation": [
            "download and parse coordinates",
            "manually verify peptide chains and HLA allele assignments",
            "annotate experimental P1-P9 register from peptide-pocket geometry",
            "deduplicate near-identical structures and freeze train/test exclusions",
            "run each register predictor blind to the annotations",
        ],
        "claim_boundary": BOUNDARY,
        "source_checksums": {
            path.name: digest(path)
            for path in [
                SOURCE / "rcsb_hla_dra_search.json",
                SOURCE / "rcsb_hla_dra_search_page2.json",
                SOURCE / "rcsb_hla_dra_entry_metadata.json",
            ]
        },
    }
    (OUT / "benchmark_status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# Solved HLA-II structure benchmark inventory\n\n"
        "This is the completed discovery and metadata-audit stage for the reviewer-requested "
        "general register benchmark. It contains every experimental RCSB entry returned for "
        "HLA-DRA UniProt P01903 and flags likely short peptide chains.\n\n"
        "It is not yet an accuracy benchmark. Length alone cannot establish which polymer is "
        "the bound peptide, and a reliable ground truth requires coordinate-level assignment "
        "of peptide residues to the P1-P9 pockets. The pending steps are frozen in "
        "`benchmark_status.json`; this prevents an inventory from being mislabeled as validation.\n",
        encoding="utf-8",
    )
    checksums = []
    checksum_path = OUT / "SHA256SUMS.csv"
    for path in sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_path):
        checksums.append({"relative_path": str(path.relative_to(OUT)), "sha256": digest(path), "bytes": path.stat().st_size})
    write_csv(checksum_path, checksums)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
