"""Analyze completed experimental IC50 rows using the frozen pilot gate."""

from __future__ import annotations

import argparse
from pathlib import Path

from charge_reversal_binding_pilot import (
    checksum_rows,
    evaluate_fitted_results,
    read_csv,
    write_csv,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    args = parser.parse_args()
    package = args.package.resolve()
    manifest = read_csv(package / "peptide_manifest.csv")
    for row in manifest:
        row["primary_mutation_panel"] = row["primary_mutation_panel"].lower() == "true"
    results = read_csv(package / "lab_handoff/fitted_binding_results_template.csv")
    try:
        summaries, gates = evaluate_fitted_results(results, manifest)
    except ValueError as exc:
        raise SystemExit(f"analysis blocked: {exc}") from None
    write_csv(package / "experimental_binding_summary.csv", summaries)
    write_csv(package / "charge_sensitive_binding_gates.csv", gates)
    write_csv(package / "SHA256SUMS.csv", checksum_rows(package))
    print(f"analyzed {len(results)} fitted rows; wrote {len(gates)} locked gate decisions")


if __name__ == "__main__":
    main()
