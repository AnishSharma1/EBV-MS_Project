#!/usr/bin/env python3
"""Verify frozen panel peptides against their stored UniProt reference sequences.

Every peptide in the frozen 80-peptide panel declares an accession and 1-based
source coordinates. This test slices the stored reference sequence at those
coordinates and requires an exact match, so a coordinate or sequence edit cannot
pass unnoticed.

It also checks the declared P1-P9 binding cores of the two routed BALF5-TALDO1
leads, since the manuscript's primary score is computed from those cores.

Offline by design: it reads the reference FASTA committed to the package, never
the live UniProt API. Refresh that FASTA with
``code/pipeline/fetch_reference_sequences.py`` when the curated panel changes.
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/pipeline"))

from fetch_reference_sequences import KNOWN_COORDINATE_EXCEPTIONS  # noqa: E402

PACKAGE = ROOT / "computational_package"
REFERENCE_DIR = PACKAGE / "data/reference_sequences_2026-10-16"
FASTA = REFERENCE_DIR / "uniprot_reference_proteins.fasta"
PANEL = PACKAGE / "data/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv"
LEGACY_MANIFEST = PACKAGE / "data/tcell_library_v2_2026-08-22/reference_sequence_manifest.csv"

# Declared V3 cores used for the primary TCR-facing BLOSUM62 score, with the
# 1-based offset of the core inside its 15-mer parent peptide.
DECLARED_LEAD_CORES = (
    # candidate, accession, peptide start, peptide, core offset, core
    ("HP36 BALF5", "P03198", 627, "TGGVYHFVKKHVHES", 5, "YHFVKKHVH"),
    ("HP36 TALDO1", "P37837", 108, "DARLSFDKDAMVARA", 4, "LSFDKDAMV"),
    ("HP47 BALF5", "P03198", 627, "TGGVYHFVKKHVHES", 4, "VYHFVKKHV"),
    ("HP47 TALDO1", "P37837", 216, "SVTKIYNYYKKFSYK", 5, "IYNYYKKFS"),
)


# Peptide labels written as \\pep{PROTEIN}{start-end} in the manuscript that are
# intentionally NOT canonical UniProt coordinates.
MANUSCRIPT_LABEL_EXCEPTIONS = {
    # Hy.2E11 control: literature 18.5 kDa MBP numbering for ENPVVHFFKNIVTPR, the
    # peptide chain of PDB 1BX2. Canonical P02686 217-231.
    ("MBP", 85, 99): "literature numbering; canonical P02686 217-231",
}


# Each exception label must be glossed in the manuscript with its canonical
# coordinates, so a reader checking it against UniProt is not misled.
MANUSCRIPT_LABEL_GLOSSES = {
    ("MBP", 85, 99): "217--231",
}


def read_fasta(path: Path) -> dict[str, str]:
    """Parse a multi-record FASTA keyed by UniProt accession."""
    records: dict[str, str] = {}
    accession: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if accession is not None:
                records[accession] = "".join(chunks)
            # sp|P37837|TALDO_HUMAN ... -> P37837
            parts = line[1:].split("|")
            accession = parts[1] if len(parts) >= 2 else parts[0].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if accession is not None:
        records[accession] = "".join(chunks)
    return records


def test_reference_fasta_present():
    assert FASTA.is_file(), f"reference FASTA missing: {FASTA}"
    records = read_fasta(FASTA)
    assert records, "reference FASTA parsed to zero records"


def test_every_panel_accession_has_a_reference():
    records = read_fasta(FASTA)
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8")))
    missing = sorted({r["accession"] for r in panel} - set(records))
    assert not missing, f"panel accessions with no stored reference sequence: {missing}"


def test_panel_peptides_match_reference_coordinates():
    records = read_fasta(FASTA)
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8")))
    assert panel, "frozen panel is empty"

    failures = []
    unused_exceptions = set(KNOWN_COORDINATE_EXCEPTIONS)
    for row in panel:
        accession = row["accession"]
        reference = records.get(accession)
        if reference is None:
            continue  # reported by test_every_panel_accession_has_a_reference
        start = int(row["source_start_1_based"])
        end = int(row["source_end_1_based"])
        peptide = row["sequence"]

        if end - start + 1 != len(peptide):
            failures.append(
                f"{accession} {start}-{end}: span {end - start + 1} != peptide length {len(peptide)}"
            )
            continue
        if end > len(reference):
            failures.append(
                f"{accession} {start}-{end}: end past reference length {len(reference)}"
            )
            continue
        observed = reference[start - 1 : end]
        if observed == peptide:
            continue

        key = (accession, start, peptide)
        if key in KNOWN_COORDINATE_EXCEPTIONS:
            unused_exceptions.discard(key)
            continue
        failures.append(
            f"{accession} {start}-{end}: declared {peptide} but reference has {observed}"
        )

    assert not failures, "peptide coordinate mismatches:\n" + "\n".join(f"  - {f}" for f in failures)
    assert not unused_exceptions, (
        "these entries are listed as known exceptions but now match the reference; "
        f"remove them from KNOWN_COORDINATE_EXCEPTIONS: {sorted(unused_exceptions)}"
    )


def test_known_exceptions_are_real_sequences():
    """Each exempted peptide must still exist somewhere in its own reference."""
    records = read_fasta(FASTA)
    failures = []
    for (accession, start, peptide), reason in KNOWN_COORDINATE_EXCEPTIONS.items():
        reference = records.get(accession, "")
        if "variant" in reason:
            continue  # deliberately differs from the reference by design
        if peptide not in reference:
            failures.append(f"{accession} {start} {peptide}: absent from reference ({reason})")
    assert not failures, "exempted peptides missing from their reference:\n" + "\n".join(
        f"  - {f}" for f in failures
    )


def test_declared_lead_cores_are_consistent():
    """The two routed leads' scoring cores must sit where the records say."""
    records = read_fasta(FASTA)
    failures = []
    for label, accession, start, peptide, offset, core in DECLARED_LEAD_CORES:
        reference = records.get(accession)
        if reference is None:
            failures.append(f"{label}: no reference sequence for {accession}")
            continue
        if reference[start - 1 : start - 1 + len(peptide)] != peptide:
            failures.append(f"{label}: parent peptide not at {accession}:{start}")
            continue
        observed = peptide[offset - 1 : offset - 1 + len(core)]
        if observed != core:
            failures.append(
                f"{label}: declared core {core} at offset {offset} but peptide gives {observed}"
            )
    assert not failures, "declared core mismatches:\n" + "\n".join(f"  - {f}" for f in failures)


def test_legacy_manifest_agrees_with_reference_fasta():
    """The pre-existing human-protein manifest must not contradict the FASTA."""
    records = read_fasta(FASTA)
    failures = []
    for row in csv.DictReader(LEGACY_MANIFEST.open(encoding="utf-8")):
        entity = row["entity"]
        stored = row["sequence"]
        recomputed = hashlib.sha256(stored.encode("utf-8")).hexdigest()
        if recomputed != row["sha256"]:
            failures.append(f"{entity}: stored sha256 does not match its own stored sequence")
        if int(row["sequence_length"]) != len(stored):
            failures.append(f"{entity}: sequence_length != len(sequence)")
        reference = records.get(entity)
        if reference is not None and reference != stored:
            failures.append(f"{entity}: manifest sequence differs from reference FASTA")
    assert not failures, "reference manifest inconsistencies:\n" + "\n".join(
        f"  - {f}" for f in failures
    )


def test_symbol_corrections_are_consistent():
    """Recorded protein-symbol errata must agree with the panel and the UniProt record."""
    corrections_path = REFERENCE_DIR / "symbol_corrections.csv"
    assert corrections_path.is_file(), "symbol_corrections.csv missing"
    corrections = list(csv.DictReader(corrections_path.open(encoding="utf-8")))
    assert corrections, "symbol_corrections.csv has no rows"
    panel = {row["candidate_id"]: row for row in csv.DictReader(PANEL.open(encoding="utf-8"))}
    headers = {}
    for line in FASTA.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            headers[line.split("|")[1]] = line
    failures = []
    for row in corrections:
        panel_row = panel.get(row["candidate_id"])
        if panel_row is None:
            failures.append(f"{row['candidate_id']}: not in the frozen panel")
            continue
        if panel_row["accession"] != row["accession"]:
            failures.append(f"{row['candidate_id']}: accession differs from the panel")
        recorded = row["panel_protein_symbol"].split("_")[0]
        accepted = {recorded, row["uniprot_gene_name"], row["uniprot_orf_name"]}
        if panel_row["protein_symbol"].split("_")[0] not in accepted:
            failures.append(f"{row['candidate_id']}: panel symbol {panel_row['protein_symbol']!r} is none of {sorted(accepted)}")
        header = headers.get(row["accession"], "")
        if f"GN={row['uniprot_gene_name']} " not in header + " ":
            failures.append(f"{row['candidate_id']}: UniProt header for {row['accession']} has no GN={row['uniprot_gene_name']}")
        if row["panel_protein_symbol"] == row["uniprot_gene_name"]:
            failures.append(f"{row['candidate_id']}: correction row does not correct anything")
    assert not failures, "symbol correction inconsistencies:\n" + "\n".join(f"  - {f}" for f in failures)


def test_manuscript_pep_labels_match_panel_peptides():
    """Every \\pep{PROTEIN}{start-end} label must name the peptide at those UniProt coordinates.

    The label must resolve to a peptide that is in the frozen panel, so a wrong
    coordinate in the manuscript text (as "EBNA1 87-101" once was, for the
    peptide at P03211 482-496) cannot pass.
    """
    records = read_fasta(FASTA)
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8")))
    accessions_by_symbol: dict[str, set[str]] = {}
    for row in panel:
        accessions_by_symbol.setdefault(row["protein_symbol"].split("_")[0], set()).add(row["accession"])
    panel_sequences = {row["sequence"] for row in panel}
    pattern = re.compile(r"\\pep\{([^}]+)\}\{(\d+)-+(\d+)\}")
    failures = []
    checked = 0
    for name in ("main.tex", "supplement.tex"):
        text = (ROOT / "manuscript" / name).read_text(encoding="utf-8")
        for symbol, start, end in sorted(set(pattern.findall(text))):
            start, end = int(start), int(end)
            if (symbol, start, end) in MANUSCRIPT_LABEL_EXCEPTIONS:
                continue
            accessions = accessions_by_symbol.get(symbol)
            if not accessions:
                failures.append(f"{name}: {symbol} {start}-{end}: no panel protein with this symbol")
                continue
            checked += 1
            slices = {acc: records[acc][start - 1 : end] for acc in accessions}
            if not any(sl in panel_sequences for sl in slices.values()):
                failures.append(
                    f"{name}: {symbol} {start}-{end}: UniProt slice {slices} is not a panel peptide"
                )
    for (symbol, start, end), canonical in MANUSCRIPT_LABEL_GLOSSES.items():
        label = f"\\pep{{{symbol}}}{{{start}-{end}}}"
        for name in ("main.tex", "supplement.tex"):
            text = (ROOT / "manuscript" / name).read_text(encoding="utf-8")
            for match in re.finditer(re.escape(label), text):
                if canonical not in text[match.end() : match.end() + 220]:
                    failures.append(f"{name}: {symbol} {start}-{end} is not glossed with canonical {canonical}")
    assert checked > 20, f"only {checked} manuscript peptide labels were checked; the pattern may be stale"
    assert not failures, "manuscript peptide labels that do not match UniProt coordinates:\n" + "\n".join(
        f"  - {f}" for f in failures
    )


if __name__ == "__main__":
    tests = [
        test_reference_fasta_present,
        test_every_panel_accession_has_a_reference,
        test_panel_peptides_match_reference_coordinates,
        test_known_exceptions_are_real_sequences,
        test_declared_lead_cores_are_consistent,
        test_legacy_manifest_agrees_with_reference_fasta,
        test_manuscript_pep_labels_match_panel_peptides,
        test_symbol_corrections_are_consistent,
    ]
    failed = 0
    for test in tests:
        try:
            test()
        except AssertionError as error:
            failed += 1
            print(f"FAIL {test.__name__}\n{error}\n")
        else:
            print(f"PASS {test.__name__}")
    if failed:
        print(f"\n{failed} of {len(tests)} checks failed")
        sys.exit(1)
    print(f"\nALL {len(tests)} PEPTIDE COORDINATE CHECKS PASSED")
