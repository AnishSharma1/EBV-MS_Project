#!/usr/bin/env python3
"""Fetch the UniProt reference sequences behind the frozen 80-peptide panel.

Reads every distinct accession from the frozen panel, retrieves those records
from the UniProt REST API in a single streamed query, and writes a FASTA plus a
provenance record into the computational package. The FASTA is what
``code/tests/test_peptide_coordinate_integrity.py`` checks peptide coordinates
against, so the test stays offline and deterministic.

Usage::

    python3 code/pipeline/fetch_reference_sequences.py            # write files
    python3 code/pipeline/fetch_reference_sequences.py --check    # verify only

Re-run this only when the curated panel changes. Refreshing it rewrites
``SHA256SUMS.csv`` and requires a corresponding, documented update to
``computational_package/MANIFEST.sha256``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "computational_package"
PANEL = PACKAGE / "data/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv"
OUT_DIR = PACKAGE / "data/reference_sequences_2026-10-16"
FASTA = OUT_DIR / "uniprot_reference_proteins.fasta"
COVERAGE = OUT_DIR / "coordinate_verification.csv"
PROVENANCE = OUT_DIR / "retrieval_provenance.json"
SUMS = OUT_DIR / "SHA256SUMS.csv"

API = "https://rest.uniprot.org/uniprotkb/stream"
TIMEOUT = 180

# Peptides whose declared coordinates do not index the canonical UniProt record.
# Each was checked individually on 2026-10-16: the peptide is real, and the
# discrepancy is a numbering convention or a documented strain variant. Keyed by
# (accession, start, declared sequence) so any NEW mismatch is still reported as
# a failure. None involves the two routed BALF5-TALDO1 leads, which match the
# canonical records exactly. See the directory README for the full write-up.
KNOWN_COORDINATE_EXCEPTIONS = {
    # MBP: declared coordinates follow 18.5 kDa "classic" MBP numbering; the
    # canonical P02686 entry is the 304-aa golli isoform. Constant offset +133.
    ("P02686", 84, "ENPVVHFFKNIVTPR"): "18.5 kDa MBP numbering; canonical P02686 217-231",
    ("P02686", 84, "ENPVVHFFKNIVTP"): "18.5 kDa MBP numbering; canonical P02686 217-230",
    ("P02686", 84, "ENPVVHFFKNIVTPRTP"): "18.5 kDa MBP numbering; canonical P02686 217-233",
    ("P02686", 86, "PVVHFFKNIVTPRTP"): "18.5 kDa MBP numbering; canonical P02686 219-233",
    # MOG: mature-protein numbering (signal peptide 1-29 removed), offset +30.
    ("Q16653", 98, "FFRDHSYQEEA"): "mature MOG numbering; canonical Q16653 128-138",
    # EBNA2: present verbatim 3 residues earlier in B95-8; repeat length varies
    # between strains.
    ("P12978", 139, "MANYIVRQSRGDRGL"): "strain offset; present at B95-8 P12978 136-150",
    # LMP2: study-variant residues retained over the B95-8 reference, as recorded
    # in the panel's reference_strain_or_isoform field.
    ("P13285", 72, "DYQPLGTQDPSLYLG"): "study variant; B95-8 P13285 72-86 reads SDYQPLGTQDQSLYL",
    ("P13285", 73, "DYQPLGNQDPSLYLG"): "study variant; B95-8 P13285 73-87 reads DYQPLGTQDQSLYLG",
}


def panel_rows() -> list[dict[str, str]]:
    with PANEL.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def accessions(rows) -> list[str]:
    return sorted({row["accession"] for row in rows})


def fetch(accs: list[str], attempts: int = 5) -> str:
    """GET the records, retrying transient transport failures with backoff.

    The request is a read-only GET, so retrying is safe. HTTP 4xx other than 429
    is a request error and is not retried.
    """
    query = " OR ".join(f"accession:{a}" for a in accs)
    url = f"{API}?{urllib.parse.urlencode({'query': query, 'format': 'fasta'})}"
    request = urllib.request.Request(url, headers={"User-Agent": "EBV-MS-package/1.0"})

    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:  # noqa: S310
                return response.read().decode("utf-8")
        except urllib.error.HTTPError as error:
            retryable = error.code == 429 or error.code >= 500
            if not retryable or attempt == attempts:
                raise
            delay = min(float(error.headers.get("Retry-After") or 0), 60.0) or _backoff(attempt)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            if attempt == attempts:
                raise
            delay = _backoff(attempt)
            print(f"  attempt {attempt}/{attempts} failed ({error}); retrying in {delay:.1f}s")
        time.sleep(delay)

    raise RuntimeError("unreachable")


def _backoff(attempt: int) -> float:
    """Exponential backoff with jitter, capped at 30 seconds."""
    return min(2.0**attempt, 30.0) * (0.5 + random.random() / 2)


def parse_fasta(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    accession: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if accession is not None:
                records[accession] = "".join(chunks)
            parts = line[1:].split("|")
            accession = parts[1] if len(parts) >= 2 else parts[0].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if accession is not None:
        records[accession] = "".join(chunks)
    return records


def verify(rows, records) -> list[dict[str, object]]:
    results = []
    for row in rows:
        acc = row["accession"]
        start = int(row["source_start_1_based"])
        end = int(row["source_end_1_based"])
        reference = records.get(acc, "")
        observed = reference[start - 1 : end] if reference else ""
        matched = observed == row["sequence"]
        reason = KNOWN_COORDINATE_EXCEPTIONS.get((acc, start, row["sequence"]), "")
        results.append(
            {
                "candidate_id": row["candidate_id"],
                "accession": acc,
                "protein_symbol": row["protein_symbol"],
                "kingdom": row["kingdom"],
                "start_1_based": start,
                "end_1_based": end,
                "declared_sequence": row["sequence"],
                "reference_slice": observed,
                "reference_length": len(reference),
                "match": matched,
                "known_exception": "" if matched else ("yes" if reason else "no"),
                "exception_reason": "" if matched else reason,
            }
        )
    return results


def write_outputs(accs, text, records, results) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FASTA.write_text(text, encoding="utf-8")

    with COVERAGE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    PROVENANCE.write_text(
        json.dumps(
            {
                "source": "UniProtKB REST API",
                "endpoint": API,
                "query": "accession OR-list of every accession in the frozen 80-peptide panel",
                "format": "fasta",
                "requested_accessions": accs,
                "returned_accessions": sorted(records),
                "retrieved_utc_date": date.today().isoformat(),
                "panel_file": str(PANEL.relative_to(PACKAGE)),
                "peptides_checked": len(results),
                "peptides_matching": sum(r["match"] for r in results),
                "note": (
                    "Canonical UniProt sequences retrieved to verify the declared 1-based "
                    "peptide coordinates in the frozen panel. Sequence versions (SV) are "
                    "recorded in the FASTA headers; a UniProt sequence revision can shift "
                    "coordinates, so re-run this check before reusing the coordinates."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    with SUMS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "sha256", "bytes"])
        for path in sorted(p for p in OUT_DIR.iterdir() if p.name != SUMS.name):
            data = path.read_bytes()
            writer.writerow([path.name, hashlib.sha256(data).hexdigest(), len(data)])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify against the stored FASTA without contacting UniProt",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="regenerate the derived records from the stored FASTA, without re-fetching",
    )
    args = parser.parse_args()

    rows = panel_rows()
    accs = accessions(rows)

    if args.check or args.rebuild:
        if not FASTA.is_file():
            print(f"stored FASTA missing: {FASTA}")
            return 1
        text = FASTA.read_text(encoding="utf-8")
        records = parse_fasta(text)
    else:
        print(f"requesting {len(accs)} accessions from UniProt ...")
        text = fetch(accs)
        records = parse_fasta(text)
        missing = sorted(set(accs) - set(records))
        if missing:
            print(f"UniProt returned no record for: {missing}")
            return 1

    if not args.check:
        write_outputs(accs, text, records, verify(rows, records))
        print(f"wrote {OUT_DIR.relative_to(ROOT)} ({len(records)} records)")

    results = verify(rows, records)
    mismatched = [r for r in results if not r["match"]]
    unexplained = [r for r in mismatched if r["known_exception"] != "yes"]

    for r in mismatched:
        tag = "KNOWN   " if r["known_exception"] == "yes" else "MISMATCH"
        detail = f" [{r['exception_reason']}]" if r["exception_reason"] else ""
        print(
            f"{tag} {r['accession']} {r['start_1_based']}-{r['end_1_based']}: "
            f"declared {r['declared_sequence']} / reference {r['reference_slice']}{detail}"
        )

    matched = len(results) - len(mismatched)
    print(
        f"{matched}/{len(results)} panel peptides match their coordinates; "
        f"{len(mismatched) - len(unexplained)} documented exceptions, "
        f"{len(unexplained)} unexplained"
    )
    return 1 if unexplained else 0


if __name__ == "__main__":
    sys.exit(main())
