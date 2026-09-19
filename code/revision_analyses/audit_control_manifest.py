#!/usr/bin/env python3
"""Audit a proposed post-cutoff class-II TCR-pMHC control panel against RCSB.

The script intentionally refuses to write a frozen panel unless every row has
an unambiguous paired TCR, its peptide is deposited as a polymer sequence, the
release date is after the cutoff, and no paired alpha/beta sequence repeats.
It uses only the public RCSB Data API and the standard library.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

CUTOFF = date(2021, 9, 30)
ROOT = Path(__file__).resolve().parents[1]


def api(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "control-panel-audit/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def digest(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def polymer_entities(pdb_id: str) -> list[dict]:
    entry = api(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    ids = entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    entities = []
    for entity_id in ids:
        record = api(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        description = record.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        identifiers = record["rcsb_polymer_entity_container_identifiers"]
        entities.append(
            {
                "entity_id": str(entity_id),
                "description": description,
                "sequence": record["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", ""),
                "auth_asym_ids": identifiers.get("auth_asym_ids", []),
            }
        )
    return entry, entities


def main(manifest: Path) -> int:
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    results, failed, pair_digests = [], [], set()
    for row in rows:
        pdb = row["pdb_id"].upper()
        try:
            entry, entities = polymer_entities(pdb)
            release = date.fromisoformat(entry["rcsb_accession_info"]["initial_release_date"][:10])
            tcrs = [
                x for x in entities
                if "t cell receptor" in x["description"].lower().replace("-", " ")
                or "tcr" in x["description"].lower()
            ]
            # Description is deposited annotation; alpha/beta naming is audited rather than inferred by chain order.
            alpha = [x for x in tcrs if "alpha" in x["description"].lower()]
            beta = [x for x in tcrs if "beta" in x["description"].lower()]
            peptide = [x for x in entities if x["sequence"] == row["peptide_sequence"]]
            if release <= CUTOFF:
                raise ValueError(f"release date {release} is not after {CUTOFF}")
            if len(alpha) != 1 or len(beta) != 1:
                raise ValueError(f"expected one alpha and one beta entity, found {len(alpha)} alpha / {len(beta)} beta")
            if len(peptide) != 1:
                raise ValueError(f"expected one exact peptide entity, found {len(peptide)}")
            pair_digest = digest(alpha[0]["sequence"] + "|" + beta[0]["sequence"])
            if pair_digest in pair_digests:
                raise ValueError("paired TCR sequence duplicates an earlier proposed control")
            pair_digests.add(pair_digest)
            citation = entry.get("rcsb_primary_citation", {})
            resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined", [""])
            row.update(
                {
                    "release_date": str(release),
                    "resolution_a": resolution[0] if resolution else "",
                    "tcr_alpha_chain": ",".join(alpha[0]["auth_asym_ids"]),
                    "tcr_beta_chain": ",".join(beta[0]["auth_asym_ids"]),
                    "peptide_chain": ",".join(peptide[0]["auth_asym_ids"]),
                    "tcr_alpha_sha256": digest(alpha[0]["sequence"]),
                    "tcr_beta_sha256": digest(beta[0]["sequence"]),
                    "pair_sha256": pair_digest,
                    "primary_pubmed_id": str(citation.get("pdbx_database_id_PubMed", "")),
                    "primary_doi": citation.get("pdbx_database_id_DOI", ""),
                    "chain_map_status": "RCSB_entity_audited",
                    "tcr_pair_status": "unique_deposited_pair_audited",
                    "selection_status": "frozen_pending_register_citation",
                    "audit_status": "PASS",
                }
            )
            # Kept in the local audit JSON and emitted to a separate FASTA on
            # successful freeze; checksums alone are not sufficient to rerun a pose generator.
            row["_tcr_alpha_sequence"] = alpha[0]["sequence"]
            row["_tcr_beta_sequence"] = beta[0]["sequence"]
        except Exception as exc:  # keep full panel evidence; no partial freeze
            row["audit_status"] = f"FAIL: {exc}"
            failed.append(f"{row['control_id']} ({pdb}): {exc}")
        results.append(row)

    audit_path = ROOT / "manifests" / "control_audit.json"
    audit_path.write_text(json.dumps(results, indent=2) + "\n")
    if failed:
        print("PANEL NOT FROZEN:", *failed, sep="\n- ", file=sys.stderr)
        print(f"Detailed audit written to {audit_path}", file=sys.stderr)
        return 2

    frozen = ROOT / "manifests" / "control_panel_frozen.tsv"
    fields = [key for key in results[0] if not key.startswith("_")]
    with frozen.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows([{key: row.get(key, "") for key in fields} for row in results])
    fasta = ROOT / "manifests" / "control_tcr_sequences.fasta"
    with fasta.open("w") as handle:
        for row in results:
            for chain, field in (("alpha", "_tcr_alpha_sequence"), ("beta", "_tcr_beta_sequence")):
                handle.write(f">{row['control_id']}|{row['pdb_id']}|TCR_{chain}|chains={row[f'tcr_{chain}_chain']}\n")
                handle.write(row[field] + "\n")
    print(f"FROZEN: {len(results)} controls -> {frozen}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: audit_control_manifest.py CONTROL_PANEL.tsv")
    raise SystemExit(main(Path(sys.argv[1])))
