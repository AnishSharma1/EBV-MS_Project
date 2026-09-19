#!/usr/bin/env python3
"""Independently validate the generated five-batch AlphaFold Server package."""

from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parent
PACKAGE = WORKSPACE / "outputs" / "alphafold_multiallele_5x30_2026-08-20"
ZIP_PATH = PACKAGE.parent / "alphafold_multiallele_5x30_2026-08-20.zip"
PROJECT = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
)
SOURCE_MANIFEST = PROJECT / "processed" / "pmhc_candidate_manifest.csv"
PRIOR_AF3_JSON = (
    PROJECT
    / "af3_migration_2026-08-10"
    / "recalibrated_30_job_batches"
    / "pmhc_90_one_seed_alphafoldserver.json"
)

EXPECTED_HLA = {
    "DRB1*13:03:01:01": ("IPD-IMGT/HLA HLA00799", 189, "2393660fc0d90e20dc0dc1f6e51805c3abe048e79e1868ab43d9589be9f43584"),
    "DRB1*03:01:01:01": ("IPD-IMGT/HLA HLA00671", 189, "cd8a608f9b400f7209a626850ed4fac571268d3e0c0d13d060464501efc4c7fd"),
    "DRB1*08:01:01:01": ("IPD-IMGT/HLA HLA00723", 189, "ed34bdce188109eac93be802e36da8c5f004205a77ccc2666d66e8deb48d3263"),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def validate() -> list[str]:
    json_files = sorted(PACKAGE.glob("ebvms_multiallele_batch_*_30_jobs.json"))
    assert len(json_files) == 5

    all_jobs: list[dict] = []
    for path in json_files:
        jobs = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(jobs, list) and len(jobs) == 30
        assert Counter(job["name"].split("_")[1] for job in jobs) == {
            "drb1303": 10,
            "drb0301": 10,
            "drb0801": 10,
        }
        for job in jobs:
            assert set(job) == {"name", "modelSeeds", "sequences", "dialect", "version"}
            assert job["dialect"] == "alphafoldserver"
            assert job["version"] == 1
            assert job["modelSeeds"] == []
            assert len(job["sequences"]) == 3
            for entity in job["sequences"]:
                assert set(entity) == {"proteinChain"}
                assert set(entity["proteinChain"]) == {"sequence", "count"}
                assert entity["proteinChain"]["count"] == 1
        all_jobs.extend(jobs)

    assert len(all_jobs) == 150
    assert len({job["name"] for job in all_jobs}) == 150

    peptide_alleles: dict[tuple[str, str], set[str]] = defaultdict(set)
    for job in all_jobs:
        _, allele, candidate = job["name"].split("_", 2)
        peptide = job["sequences"][2]["proteinChain"]["sequence"]
        peptide_alleles[(candidate, peptide)].add(allele)
    assert len(peptide_alleles) == 50
    assert all(
        alleles == {"drb1303", "drb0301", "drb0801"}
        for alleles in peptide_alleles.values()
    )

    job_rows = read_csv(PACKAGE / "job_manifest_150.csv")
    assert len(job_rows) == 150
    assert [row["job_name"] for row in job_rows] == [job["name"] for job in all_jobs]

    panel_rows = read_csv(PACKAGE / "peptide_panel_manifest.csv")
    assert len(panel_rows) == 50
    assert Counter(row["arm_group"] for row in panel_rows) == {"EBV": 25, "CNS/self": 25}
    source = {row["candidate_id"]: row for row in read_csv(SOURCE_MANIFEST)}
    for row in panel_rows:
        if row["candidate_id"] == "GLIALCAM_370_389_UNMODIFIED":
            assert row["peptide_sequence"] == "ATGRTHSSPPRAPSSPGRSR"
            continue
        original = source[row["candidate_id"]]
        assert row["peptide_sequence"] == original["peptide"]
        assert int(row["peptide_length"]) == int(original["peptide_length"])

    hla_rows = read_csv(PACKAGE / "hla_sequence_manifest.csv")
    assert len(hla_rows) == 4
    prior_jobs = json.loads(PRIOR_AF3_JSON.read_text(encoding="utf-8"))
    prior_dra = prior_jobs[0]["sequences"][0]["proteinChain"]["sequence"]
    assert hla_rows[0]["sequence"] == prior_dra
    assert hla_rows[0]["sha256"] == hashlib.sha256(prior_dra.encode()).hexdigest()
    for row in hla_rows[1:]:
        accession, length, digest = EXPECTED_HLA[row["allele_or_name"]]
        assert row["source_accession"] == accession
        assert int(row["sequence_length"]) == length
        assert row["sha256"] == digest
        assert hashlib.sha256(row["sequence"].encode()).hexdigest() == digest

    checksum_lines = (PACKAGE / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    for line in checksum_lines:
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((PACKAGE / name).read_bytes()).hexdigest() == digest

    with zipfile.ZipFile(ZIP_PATH) as archive:
        prefix = PACKAGE.name + "/"
        zipped_files = {
            name[len(prefix) :]
            for name in archive.namelist()
            if name.startswith(prefix) and not name.endswith("/")
        }
    disk_files = {path.name for path in PACKAGE.iterdir() if path.is_file()}
    assert zipped_files == disk_files

    return [
        "Independent validation: PASS",
        "5 JSON files; 30 jobs per file; 150 unique jobs total",
        "50 unique peptides crossed once with each of 3 alleles",
        "25 EBV and 25 CNS/self peptides",
        "49 peptides match the frozen project manifest; GlialCAM 370-389 verified separately",
        "DRB sequence identities, IPD accessions, lengths, and hashes match the frozen sequence manifest",
        "ZIP contents exactly match the package directory; all listed SHA-256 checksums pass",
    ]


if __name__ == "__main__":
    for message in validate():
        print(message)
