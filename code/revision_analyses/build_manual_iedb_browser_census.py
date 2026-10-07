#!/usr/bin/env python3
"""Materialize the 98-arm IEDB public-browser census from an audited UI pass.

The IEDB Query API was unavailable for the original retrieval.  This script
does not call an API: its small, explicit observation ledger is transcribed
from exact-peptide searches in IEDB's public browser results interface on
2026-09-20.  Each distinct peptide was searched once; all requested alleles
for that peptide were then compared literally against the displayed rows.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/revision_round_4/threshold_sensitivity_2026-09-20/threshold_grid_pair_results.csv"
OUT = ROOT / "output/revision_round_4/iedb_manual_browser_census_2026-09-20"
CHECKED_AT = "2026-09-20"
IEDB_URL = "https://www.iedb.org/"

# Exact-sequence public result existed for these peptides.  Every peptide not
# listed here visibly returned zero epitopes/antigens/assays/references.
EXACT_SEQUENCE_HITS = {
    "EFIWMCMTVRHRCQA", "DARLSFDKDAMVARA", "TVFYNIPPMPL",
    "LGERVTLTSYWRRVS", "MDIAIHHPWIRRPFF", "RLYELLSDPRSALGL",
    "FFRDHSYQEEA", "ENPVVHFFKNIVTPRTP", "MANYIVRQSRGDRGL",
    "RTELLKDAIGEGK", "ILCFVMAARQRLQDI", "AEGLRALLARSHVER",
    "QEPRRVRRRVLVQQE", "LTAYHVSTAPTGSWF", "VTLTSYWRRVSLGPE",
    "QHYREVAAAKSSE", "PYYVVDLSVRGM", "QQFLRLLGATYVLRV",
    "VTTFVHYANFHNYFVPN", "MEVGWYRPPFSRVVHLYRNGK",
    "PVVHFFKNIVTPRTP", "TGGVYHFVKKHVHES", "VTNILIYNGWYADS",
    "IMNILRIYYSPSIM", "SRDELLHTRAASLLY", "MLGQDDFIKFKSPLV",
}

# Only the following row meets the locked exact-sequence/exact-allele positive
# MHC/ligand criterion.  A literal allele match is required: broad labels such
# as HLA class II, HLA-DR2, or HLA-DR13 do not qualify.
MHC_POSITIVE = {
    ("TGGVYHFVKKHVHES", "HLA-DRB1*15:01"): {
        "iedb_epitope_id": "63843",
        "iedb_assay_id": "1774674",
        "assay_type": "purified MHC/competitive/radioactivity qualitative binding",
        "qualitative_outcome": "Positive",
        "reference": "B E Hansen; Tissue Antigens; 2007",
    },
}

# Kept as non-qualifying context so this evidence is visible rather than
# silently discarded.  The gate requires an MHC/ligand assay, not a T-cell
# response experiment.
TCELL_POSITIVE = {
    ("ENPVVHFFKNIVTPRTP", "HLA-DRB1*15:01"): "1666849;1686870",
    ("AEGLRALLARSHVER", "HLA-DRB1*15:01"): "12579318",
    ("TGGVYHFVKKHVHES", "HLA-DRB1*15:01"): "1686220;1774672",
    ("VTNILIYNGWYADS", "HLA-DRB1*15:01"): "22930779",
    ("IMNILRIYYSPSIM", "HLA-DRB1*15:01"): "22930778",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    by_candidate: dict[str, dict[str, str]] = {}
    for row in read_rows(SOURCE):
        by_candidate.setdefault(row["candidate_id"], row)
    if len(by_candidate) != 49:
        raise ValueError(f"Expected 49 candidates, found {len(by_candidate)}")

    rows: list[dict[str, str]] = []
    for candidate_id in sorted(by_candidate, key=lambda value: int(value[2:])):
        source = by_candidate[candidate_id]
        for arm in ("ebv", "self"):
            peptide = source[f"{arm}_sequence"]
            allele = source["allele"]
            key = (peptide, allele)
            exact_sequence = peptide in EXACT_SEQUENCE_HITS
            mhc = MHC_POSITIVE.get(key)
            tcell_ids = TCELL_POSITIVE.get(key, "")
            if mhc:
                classification = "exact_sequence_exact_hla_positive_mhc_ligand"
                decision = "qualifies_public_evidence_gate"
            elif tcell_ids:
                classification = "exact_sequence_exact_hla_positive_tcell_only"
                decision = "does_not_qualify_mhc_ligand_gate"
            elif exact_sequence:
                classification = "exact_sequence_no_positive_exact_hla_mhc_ligand"
                decision = "does_not_qualify_mhc_ligand_gate"
            else:
                classification = "no_exact_sequence_public_record"
                decision = "not_evidence_of_absence_or_penalty"
            rows.append({
                "candidate_id": candidate_id,
                "arm": arm,
                "peptide": peptide,
                "requested_hla": allele,
                "manual_retrieval_status": "evaluated_public_browser",
                "exact_sequence_public_record": str(exact_sequence),
                "exact_hla_positive_mhc_ligand": str(bool(mhc)),
                "iedb_epitope_id": mhc["iedb_epitope_id"] if mhc else "",
                "iedb_mhc_assay_ids": mhc["iedb_assay_id"] if mhc else "",
                "iedb_exact_hla_tcell_assay_ids": tcell_ids,
                "classification": classification,
                "decision_effect": decision,
                "browser_check_date": CHECKED_AT,
                "source_url": IEDB_URL,
            })
    if len(rows) != 98:
        raise ValueError(f"Expected 98 arms, found {len(rows)}")
    write_rows(OUT / "iedb_manual_browser_census.csv", rows)
    candidate_rows = []
    for candidate_id in sorted(by_candidate, key=lambda value: int(value[2:])):
        arms = [row for row in rows if row["candidate_id"] == candidate_id]
        qualifying = [row for row in arms if row["exact_hla_positive_mhc_ligand"] == "True"]
        candidate_rows.append({
            "candidate_id": candidate_id,
            "allele": by_candidate[candidate_id]["allele"],
            "ebv_peptide": next(row["peptide"] for row in arms if row["arm"] == "ebv"),
            "self_peptide": next(row["peptide"] for row in arms if row["arm"] == "self"),
            "public_evidence_gate": "pass" if qualifying else "not_passed",
            "qualifying_arm_ids": ";".join(row["arm"] for row in qualifying),
            "qualifying_iedb_mhc_assay_ids": ";".join(row["iedb_mhc_assay_ids"] for row in qualifying),
            "non_penalizing_no_record_arm_count": str(sum(row["classification"] == "no_exact_sequence_public_record" for row in arms)),
        })
    write_rows(OUT / "candidate_public_evidence_gate.csv", candidate_rows)
    counts = Counter(row["classification"] for row in rows)
    report = {
        "status": "complete_manual_public_browser_census",
        "candidate_count": 49,
        "arm_count": 98,
        "unique_exact_peptide_searches": len({row["peptide"] for row in rows}),
        "classification_counts": dict(sorted(counts.items())),
        "qualifying_mhc_ligand_arm_count": sum(row["exact_hla_positive_mhc_ligand"] == "True" for row in rows),
        "rule": "Positive exact-sequence/exact-allele MHC/ligand evidence qualifies; T-cell-only evidence is retained but does not qualify.",
        "absence_policy": "No exact-sequence record is recorded as no public record, not as negative evidence or a ranking penalty.",
        "source": "Manual exact-peptide queries in the IEDB public browser interface; no Query API call was used.",
    }
    (OUT / "census_status.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# IEDB manual public-browser census\n\n"
        "All 98 frozen candidate arms were evaluated through 35 exact-peptide searches in IEDB's public browser interface. "
        "The browser interface was available even though the prior Query API retrieval failed.\n\n"
        "One arm, HP47 EBV (BALF5 627-641, `TGGVYHFVKKHVHES`) at HLA-DRB1*15:01, has a positive exact-sequence/exact-allele MHC binding record (IEDB assay 1774674). "
        "Exact-HLA T-cell-only records are retained separately and do not satisfy the locked MHC/ligand criterion. "
        "No-hit searches remain non-penalizing: they are not evidence of absence.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
