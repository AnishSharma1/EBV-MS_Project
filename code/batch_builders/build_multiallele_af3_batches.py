#!/usr/bin/env python3
"""Build five balanced 30-job AlphaFold Server pMHC discovery batches."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication"
)
SOURCE_MANIFEST = PROJECT_ROOT / "processed" / "pmhc_candidate_manifest.csv"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs" / "alphafold_multiallele_5x30_2026-08-20"
ZIP_PATH = OUTPUT_DIR.parent / "alphafold_multiallele_5x30_2026-08-20.zip"

DRA_SEQUENCE = (
    "EHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKAN"
    "LEIMTKRSNYTPITNVPPEVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETV"
    "FLPREDHLFRKFHYLPFLPSTEDVYDCRVEHWGLDEPLLKHWEFD"
)

ALLELES = [
    {
        "allele": "HLA-DRB1*13:03",
        "allele_code": "drb1303",
        "four_field_name": "DRB1*13:03:01:01",
        "ipd_accession": "HLA00799",
        "sequence": (
            "DTRPRFLEYSTSECHFFNGTERVRFLDRYFYNQEEYVRFDSDVGEYRAVTELGRPSAEYWNSQKDI"
            "LEDKRAAVDTYCRHNYGVGESFTVQRRVHPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFR"
            "NGQEEKTGVVSTGLIHNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRA"
        ),
    },
    {
        "allele": "HLA-DRB1*03:01",
        "allele_code": "drb0301",
        "four_field_name": "DRB1*03:01:01:01",
        "ipd_accession": "HLA00671",
        "sequence": (
            "DTRPRFLEYSTSECHFFNGTERVRYLDRYFHNQEENVRFDSDVGEFRAVTELGRPDAEYWNSQKDL"
            "LEQKRGRVDNYCRHNYGVVESFTVQRRVHPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFR"
            "NGQEEKTGVVSTGLIHNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRA"
        ),
    },
    {
        "allele": "HLA-DRB1*08:01",
        "allele_code": "drb0801",
        "four_field_name": "DRB1*08:01:01:01",
        "ipd_accession": "HLA00723",
        "sequence": (
            "DTRPRFLEYSTGECYFFNGTERVRFLDRYFYNQEEYVRFDSDVGEYRAVTELGRPSAEYWNSQKDF"
            "LEDRRALVDTYCRHNYGVGESFTVQRRVHPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFR"
            "NGQEEKTGVVSTGLIHNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWSA"
        ),
    },
]

# Ten peptides per batch: five EBV and five CNS/self peptides. Every peptide is
# modeled against all three alleles, producing 30 jobs in each upload file.
BATCH_PANELS = [
    [
        "EBV_TCELL_950",
        "EBV_TCELL_63843",
        "EBV_TCELL_2268683",
        "EBV_TCELL_2268933",
        "EBV_TCELL_119155",
        "HUMAN_MYELIN_112214",
        "HUMAN_MYELIN_114806",
        "HUMAN_MYELIN_5516",
        "HUMAN_MYELIN_112226",
        "GLIALCAM_370_389_UNMODIFIED",
    ],
    [
        "EBV_TCELL_2268741",
        "EBV_TCELL_2268720",
        "EBV_TCELL_2268934",
        "EBV_TCELL_1862913",
        "EBV_TCELL_1393562",
        "HUMAN_MYELIN_112191",
        "HUMAN_MYELIN_112319",
        "HUMAN_MYELIN_35803",
        "HUMAN_MYELIN_40353",
        "HUMAN_MYELIN_118650",
    ],
    [
        "EBV_TCELL_149795",
        "EBV_TCELL_1392779",
        "EBV_MHC_15168",
        "EBV_MHC_38633",
        "EBV_MHC_45380",
        "HUMAN_MYELIN_112782",
        "HUMAN_MYELIN_114566",
        "HUMAN_MYELIN_40762",
        "HUMAN_MYELIN_51191",
        "HUMAN_MYELIN_63975",
    ],
    [
        "EBV_MHC_53359",
        "EBV_MHC_48976",
        "EBV_MHC_12325",
        "EBV_MHC_26480",
        "EBV_MHC_3005",
        "HUMAN_MYELIN_115061",
        "HUMAN_MYELIN_115622",
        "HUMAN_MYELIN_115641",
        "HUMAN_MYELIN_112354",
        "HUMAN_MYELIN_73159",
    ],
    [
        "EBV_MHC_1997",
        "EBV_MHC_23636",
        "EBV_MHC_20020",
        "EBV_MHC_41628",
        "EBV_MHC_11422",
        "HUMAN_MYELIN_115648",
        "HUMAN_MYELIN_116924",
        "HUMAN_MYELIN_119927",
        "HUMAN_MYELIN_112603",
        "HUMAN_MYELIN_112756",
    ],
]

GLIALCAM_RECORD = {
    "candidate_id": "GLIALCAM_370_389_UNMODIFIED",
    "arm": "Human CNS calibration",
    "evidence_tier": "Published EBNA1-GlialCAM calibration peptide",
    "peptide": "ATGRTHSSPPRAPSSPGRSR",
    "peptide_length": "20",
    "source_antigen": "Glial cell adhesion molecule (GlialCAM)",
    "source_accession": "UniProt:Q14CZ8; residues 370-389",
    "iedb_assay_id": "",
    "iedb_epitope_id": "",
    "pubmed_id": "35025605",
    "hla": "HLA-DRB1*15:01/DRB5*01:01 published context",
    "mhc_class": "II",
    "modeling_status": "cross-allele calibration only",
}

AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def load_source_records() -> dict[str, dict[str, str]]:
    with SOURCE_MANIFEST.open(newline="", encoding="utf-8-sig") as handle:
        records = {row["candidate_id"]: row for row in csv.DictReader(handle)}
    records[GLIALCAM_RECORD["candidate_id"]] = GLIALCAM_RECORD
    return records


def make_job(candidate: dict[str, str], allele: dict[str, str]) -> dict:
    name = f"ebvms_{allele['allele_code']}_{candidate['candidate_id'].lower()}"
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": DRA_SEQUENCE, "count": 1}},
            {"proteinChain": {"sequence": allele["sequence"], "count": 1}},
            {"proteinChain": {"sequence": candidate["peptide"], "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_jobs(all_jobs: list[dict], job_rows: list[dict], panel_rows: list[dict]) -> list[str]:
    checks: list[str] = []
    assert len(BATCH_PANELS) == 5
    assert all(len(panel) == 10 for panel in BATCH_PANELS)
    assert len({item for panel in BATCH_PANELS for item in panel}) == 50
    checks.append("PASS: five panels contain 10 unique peptides each (50 unique peptides total)")

    assert len(DRA_SEQUENCE) == 178 and set(DRA_SEQUENCE) <= AA
    assert all(len(a["sequence"]) == 189 and set(a["sequence"]) <= AA for a in ALLELES)
    checks.append("PASS: HLA-DRA is 178 aa; each allele-specific DRB extracellular chain is 189 aa")

    assert len(all_jobs) == 150
    assert len({job["name"] for job in all_jobs}) == 150
    checks.append("PASS: 150 jobs and 150 unique AlphaFold Server job names")

    for job in all_jobs:
        assert job["dialect"] == "alphafoldserver"
        assert job["version"] == 1
        assert job["modelSeeds"] == []
        assert len(job["sequences"]) == 3
        for entity in job["sequences"]:
            chain = entity["proteinChain"]
            assert chain["count"] == 1
            assert chain["sequence"] and set(chain["sequence"]) <= AA
    checks.append("PASS: every job uses the validated three-protein-chain AlphaFold Server schema")

    batch_counts = Counter(row["batch"] for row in job_rows)
    assert batch_counts == {f"batch_{i:02d}": 30 for i in range(1, 6)}
    checks.append("PASS: each upload batch contains exactly 30 jobs")

    for batch in sorted(batch_counts):
        rows = [row for row in job_rows if row["batch"] == batch]
        assert Counter(row["allele"] for row in rows) == {a["allele"]: 10 for a in ALLELES}
        assert Counter(row["arm_group"] for row in rows) == {"EBV": 15, "CNS/self": 15}
    checks.append("PASS: every batch has 10 jobs per allele and a 15 EBV / 15 CNS-self balance")

    assert Counter(row["allele"] for row in job_rows) == {a["allele"]: 50 for a in ALLELES}
    assert Counter(row["arm_group"] for row in panel_rows) == {"EBV": 25, "CNS/self": 25}
    checks.append("PASS: full package has 50 jobs per allele and a 25 EBV / 25 CNS-self peptide panel")
    return checks


def main() -> None:
    source_records = load_source_records()
    selected_ids = [item for panel in BATCH_PANELS for item in panel]
    missing = sorted(set(selected_ids) - set(source_records))
    if missing:
        raise SystemExit(f"Missing source-manifest candidates: {missing}")

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    panel_rows: list[dict] = []
    job_rows: list[dict] = []
    all_jobs: list[dict] = []

    for batch_number, candidate_ids in enumerate(BATCH_PANELS, start=1):
        batch_name = f"batch_{batch_number:02d}"
        jobs: list[dict] = []
        for peptide_slot, candidate_id in enumerate(candidate_ids, start=1):
            candidate = source_records[candidate_id]
            peptide = candidate["peptide"].strip().upper()
            if not peptide or set(peptide) - AA:
                raise SystemExit(f"Invalid peptide sequence for {candidate_id}: {peptide}")
            if int(candidate["peptide_length"]) != len(peptide):
                raise SystemExit(f"Length mismatch for {candidate_id}")
            candidate = dict(candidate, peptide=peptide)
            arm_group = "EBV" if candidate["arm"] == "EBV" else "CNS/self"
            panel_rows.append(
                {
                    "batch": batch_name,
                    "peptide_slot": peptide_slot,
                    "candidate_id": candidate_id,
                    "arm_group": arm_group,
                    "original_arm": candidate["arm"],
                    "evidence_tier": candidate["evidence_tier"],
                    "source_antigen": candidate["source_antigen"],
                    "source_accession": candidate["source_accession"],
                    "peptide_sequence": peptide,
                    "peptide_length": len(peptide),
                    "pubmed_id": candidate.get("pubmed_id", ""),
                    "original_hla_context": candidate.get("hla", ""),
                    "selection_basis": (
                        "all Tier-1 EBV T-cell candidates retained"
                        if candidate_id.startswith("EBV_TCELL")
                        else "source-coverage transfer panel; no new-allele AF3 score used"
                    ),
                }
            )
            for allele in ALLELES:
                job = make_job(candidate, allele)
                jobs.append(job)
                all_jobs.append(job)
                job_rows.append(
                    {
                        "global_job_number": len(all_jobs),
                        "batch": batch_name,
                        "batch_job_number": len(jobs),
                        "job_name": job["name"],
                        "candidate_id": candidate_id,
                        "arm_group": arm_group,
                        "allele": allele["allele"],
                        "allele_four_field_name": allele["four_field_name"],
                        "ipd_accession": allele["ipd_accession"],
                        "peptide_sequence": peptide,
                        "peptide_length": len(peptide),
                        "model_seeds": "auto (one server seed)",
                    }
                )
        output_file = OUTPUT_DIR / f"ebvms_multiallele_{batch_name}_30_jobs.json"
        output_file.write_text(json.dumps(jobs, indent=2) + "\n", encoding="utf-8")

    panel_fields = [
        "batch", "peptide_slot", "candidate_id", "arm_group", "original_arm",
        "evidence_tier", "source_antigen", "source_accession", "peptide_sequence",
        "peptide_length", "pubmed_id", "original_hla_context", "selection_basis",
    ]
    job_fields = [
        "global_job_number", "batch", "batch_job_number", "job_name", "candidate_id",
        "arm_group", "allele", "allele_four_field_name", "ipd_accession",
        "peptide_sequence", "peptide_length", "model_seeds",
    ]
    write_csv(OUTPUT_DIR / "peptide_panel_manifest.csv", panel_rows, panel_fields)
    write_csv(OUTPUT_DIR / "job_manifest_150.csv", job_rows, job_fields)

    hla_rows = [
        {
            "chain": "HLA-DRA",
            "allele_or_name": "mature HLA-DRA extracellular chain (project reference)",
            "source_accession": "existing validated DRB1*15:01 project input",
            "sequence_length": len(DRA_SEQUENCE),
            "sha256": sha256_text(DRA_SEQUENCE),
            "sequence": DRA_SEQUENCE,
        }
    ]
    for allele in ALLELES:
        hla_rows.append(
            {
                "chain": "HLA-DRB",
                "allele_or_name": allele["four_field_name"],
                "source_accession": f"IPD-IMGT/HLA {allele['ipd_accession']}",
                "sequence_length": len(allele["sequence"]),
                "sha256": sha256_text(allele["sequence"]),
                "sequence": allele["sequence"],
            }
        )
    write_csv(
        OUTPUT_DIR / "hla_sequence_manifest.csv",
        hla_rows,
        ["chain", "allele_or_name", "source_accession", "sequence_length", "sha256", "sequence"],
    )

    checks = validate_jobs(all_jobs, job_rows, panel_rows)
    (OUTPUT_DIR / "VALIDATION_REPORT.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")

    readme = f"""# EBV-MS multi-allele AlphaFold Server package

## What this package is

- Five upload-ready JSON files, each containing exactly 30 three-chain pMHC jobs.
- 150 unique complexes total: 50 frozen peptides crossed once with each of HLA-DRB1*13:03, HLA-DRB1*03:01, and HLA-DRB1*08:01.
- Every batch contains 10 peptides x 3 alleles, with 15 EBV jobs and 15 CNS/self jobs.
- `modelSeeds: []` requests one automatic AlphaFold Server seed per unique complex.

## Scientific boundary

This is a targeted cross-allele **discovery/transfer screen** of the project's existing source-verified peptide panel. It is not the completed allele-specific IEDB universe, it does not establish allele-specific presentation, and it is not a robustness replicate set. The three alleles must be analyzed independently; do not pool raw structural scores across alleles.

## Sequence convention

- Chain A: the same 178-aa mature extracellular HLA-DRA sequence used in the existing project.
- Chain B: a 189-aa mature extracellular DRB chain, using IPD-IMGT/HLA first-listed four-field alleles and the same N-terminal `DTR...`/extracellular truncation convention as the existing project.
- Chain C: the exact full peptide sequence from the frozen source manifest; GlialCAM 370-389 is unmodified in this discovery screen.
- IPD accessions: HLA00799 (DRB1*13:03:01:01), HLA00671 (DRB1*03:01:01:01), and HLA00723 (DRB1*08:01:01:01).

## Upload order

Upload `batch_01` through `batch_05` separately. Each is already at the requested 30-job limit. Preserve the downloaded job names so results can be joined directly to `job_manifest_150.csv`.

## Files

- `ebvms_multiallele_batch_01_30_jobs.json` through `batch_05`: AlphaFold Server uploads
- `job_manifest_150.csv`: one row per complex
- `peptide_panel_manifest.csv`: one row per frozen peptide
- `hla_sequence_manifest.csv`: exact HLA sequences, accessions, lengths, and checksums
- `VALIDATION_REPORT.txt`: machine-run integrity checks

## Required downstream gates

Before interpretation: verify complete downloads, apply pMHC/register quality checks, fit the HLA groove, compare only predeclared exposed positions, and retain the project's claim boundary. These outputs can prioritize hypotheses; they cannot by themselves establish natural presentation, shared-TCR binding, T-cell activation, cross-reactivity, molecular mimicry, or an EBV-driven MS mechanism.
"""
    (OUTPUT_DIR / "README.md").write_text(readme, encoding="utf-8")

    checksum_lines = []
    for path in sorted(OUTPUT_DIR.iterdir()):
        if path.name == "SHA256SUMS.txt" or not path.is_file():
            continue
        checksum_lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (OUTPUT_DIR / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")

    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUTPUT_DIR.iterdir()):
            archive.write(path, arcname=f"{OUTPUT_DIR.name}/{path.name}")

    print(f"Created {OUTPUT_DIR}")
    print(f"Created {ZIP_PATH}")
    for check in checks:
        print(check)


if __name__ == "__main__":
    main()
