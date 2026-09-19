"""Execute prepared PANDORA runs with explicit templates and anchors."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

from build_tcell_library_v2 import DRA_SEQUENCE, DRB1_1501_SEQUENCE


TEMPLATE_META = {
    "1BX2": {"peptide": "ENPVVHFFKNIVTP", "anchors": [5, 8, 10, 13]},
    "6CQQ": {"peptide": "RFYKTLRAEQASQ", "anchors": [3, 6, 8, 11]},
}
ALLELES = ["HLA-DRA*01:01", "HLA-DRB1*15:01"]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def configure_runtime(package: Path):
    import PANDORA
    from PANDORA import Align, Database, Template

    environment_bin = str(Path(sys.executable).resolve().parent)
    os.environ["PATH"] = environment_bin + os.pathsep + os.environ.get("PATH", "")

    original_align = Align.Align

    class RetainTargetTerminiAlign(original_align):
        """Prevent PANDORA from silently deleting target residues absent in a template."""

        def __init__(self, target, template, clip_C_domain=False, remove_terms=True):
            super().__init__(target, template, clip_C_domain=clip_C_domain, remove_terms=False)

    Align.Align = RetainTargetTerminiAlign

    runtime = Path.home() / ".cache" / "ebv_ms_tools" / "charge-reversal-pandora-runtime" / "data"
    pdb_dir = runtime / "PDBs" / "pMHCII"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for pdb_id in TEMPLATE_META:
        source = package / "pandora" / "templates" / "pmhc_only" / f"{pdb_id}_MNP_pandora.pdb"
        shutil.copy2(source, pdb_dir / f"{pdb_id}.pdb")
    PANDORA.PANDORA_data = str(runtime)
    database = Database.Database()
    templates = {
        pdb_id: Template(
            id=pdb_id,
            peptide=meta["peptide"],
            allele_type=ALLELES,
            MHC_class="II",
            anchors=meta["anchors"],
        )
        for pdb_id, meta in TEMPLATE_META.items()
    }
    return database, templates


def run_manifest(
    package: Path,
    manifest: Path,
    models_override: int | None = None,
    limit: int | None = None,
    jobs_per_run: int = 4,
    offset: int = 0,
    resume: bool = False,
    template_only: str | None = None,
) -> Path:
    from PANDORA import Pandora, Target

    database, templates = configure_runtime(package)
    rows = read_rows(manifest)
    if template_only:
        rows = [row for row in rows if row["forced_template_pdb"] == template_only]
    rows = rows[offset:]
    if limit is not None:
        rows = rows[:limit]
    results_dir = package / "pandora" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    scratch_results = Path.home() / ".cache" / "ebv_ms_tools" / "charge-reversal-pandora-runtime" / "runs"
    scratch_results.mkdir(parents=True, exist_ok=True)
    status_rows: list[dict[str, object]] = []
    execution_label = manifest.stem + (f"_offset{offset}_limit{limit}" if offset or limit is not None else "")
    for row in rows:
        run_id = row["run_id"]
        anchors = [int(value) for value in row["anchors_1_based"].split(",")]
        requested = models_override or int(row["requested_models"])
        started = time.time()
        status = "complete"
        error = ""
        produced = 0
        packaged_target = results_dir / run_id
        existing = len(list(packaged_target.glob("*BL*.pdb"))) if packaged_target.exists() else 0
        if resume and existing == requested:
            status_rows.append({
                "run_id": run_id, "forced_template_pdb": row["forced_template_pdb"],
                "register_label": row["register_label"], "anchors_1_based": row["anchors_1_based"],
                "requested_models": requested, "produced_models": existing,
                "status": "complete_existing", "elapsed_seconds": 0.0, "error": "",
            })
            write_rows(results_dir / f"{execution_label}_execution_status.csv", status_rows)
            continue
        try:
            scratch_target = scratch_results / run_id
            if scratch_target.exists():
                shutil.rmtree(scratch_target)
            target = Target(
                id=run_id,
                peptide=row["peptide"],
                allele_type=ALLELES,
                MHC_class="II",
                M_chain_seq=DRA_SEQUENCE,
                N_chain_seq=DRB1_1501_SEQUENCE,
                anchors=anchors,
                output_dir=str(scratch_results),
            )
            case = Pandora.Pandora(target, database=database, template=templates[row["forced_template_pdb"]])
            case.model(
                n_loop_models=requested,
                n_jobs=min(jobs_per_run, requested),
                loop_refinement="slow",
                pickle_out=False,
                benchmark=False,
                verbose=False,
            )
            scratch_target = scratch_results / target.id
            produced = len(list(scratch_target.glob("*BL*.pdb")))
            if produced != requested:
                status = "incomplete_model_count"
                error = f"expected {requested}, found {produced} BL PDB models"
            packaged_target = results_dir / target.id
            if packaged_target.exists():
                shutil.rmtree(packaged_target)
            shutil.copytree(scratch_target, packaged_target)
        except Exception as exc:  # preserve per-run failure and continue batch
            status = "failed"
            error = f"{type(exc).__name__}: {exc}"
        status_rows.append({
            "run_id": run_id,
            "forced_template_pdb": row["forced_template_pdb"],
            "register_label": row["register_label"],
            "anchors_1_based": row["anchors_1_based"],
            "requested_models": requested,
            "produced_models": produced,
            "status": status,
            "elapsed_seconds": round(time.time() - started, 3),
            "error": error,
        })
        write_rows(results_dir / f"{execution_label}_execution_status.csv", status_rows)
    status_path = results_dir / f"{execution_label}_execution_status.csv"
    (results_dir / f"{execution_label}_execution_summary.json").write_text(json.dumps({
        "manifest": str(manifest.resolve()),
        "runs_attempted": len(status_rows),
        "runs_complete": sum(row["status"] == "complete" for row in status_rows),
        "models_produced": sum(int(row["produced_models"]) for row in status_rows),
        "models_are_biological_replicates": False,
    }, indent=2) + "\n", encoding="utf-8")
    return status_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--models-override", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--jobs-per-run", type=int, default=4)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--template-only", choices=("1BX2", "6CQQ"))
    args = parser.parse_args()
    path = run_manifest(
        args.package.resolve(), args.manifest.resolve(), args.models_override,
        args.limit, args.jobs_per_run, args.offset, args.resume, args.template_only,
    )
    print(path)


if __name__ == "__main__":
    main()
