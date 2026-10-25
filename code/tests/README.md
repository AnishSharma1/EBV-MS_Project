# Test suite

51 test modules. Each one is an executable `unittest` module, so it can be run
on its own or through the runner.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python code/tests/run_suite.py
```

Verified on 2026-10-24 (CPython 3.14.7, macOS arm64, versions as pinned in
`requirements.txt`):

| | count | |
|---|---|---|
| pass | 48 | with the external `processed/` tier reachable |
| pass | 39 | from this repository alone |
| skip | 3 | inputs absent from every tree on this machine |

The runner exits non-zero only on a real failure. Skips are printed
individually with their reason; they are never folded into the pass count.

## Two gates are the documented ones

These are the checks named in `CLAUDE.md`, and both run in CI:

```sh
python3 computational_package/scripts/verify_package.py     # 219 checksums
python3 code/tests/test_repository_integrity.py             # 9 checks
python3 code/tests/test_peptide_coordinate_integrity.py     # 8 checks
```

The third gate verifies every frozen-panel peptide against its stored UniProt
reference sequence, so a coordinate or sequence edit cannot pass unnoticed.

## The external `processed/` tier

Nine tests read an intermediate data tier that was never committed to this
repository: `ebv_ms_publication/processed/` (8.0 GB, 57 directories), which sits
beside this checkout rather than inside it. It is not under version control and
carries no checksum manifest — treat anything read from it as unversioned.

`run_suite.py` looks for it, and when it is found links it to `code/processed`,
the path the analysis scripts expect. The link is gitignored. Point the runner
elsewhere with `EBV_MS_PROCESSED_ROOT`, or pass `--no-link` to skip those tests
instead.

Without the tier the suite still reports 39 passed, 0 failed, 12 skipped.

## Why each module bootstraps its own path

The modules under test live in `code/pipeline/`. Tests import them by bare name,
so each test file inserts that directory on `sys.path`:

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "code" / "pipeline"))
```

Before 2026-10-24 these lines pointed at `code/src`, the layout used before the
code was reorganised into `code/` subdirectories, and 23 files had no bootstrap
at all. That is why the suite reported 3 of 51 passing: the failures were import
errors, not analysis defects.

## Known-blocked modules

Recorded in `MISSING_UPSTREAM_DATA` in `run_suite.py` so they stay visible:

- `test_audit_complete_model_pipeline` — expects one AlphaFold 3 model root; none present.
- `test_build_missing_v2_af3_batches` — reads `tcell_library_v2_2026-08-22/` from the
  CloudDocs `ebv_ms_publication` tree, where that file is absent. The file does exist
  in the tree beside this checkout, so the two trees have diverged; which one is
  authoritative is a provenance decision and has not been made here.
- `test_multiallele_validator` — needs `outputs/alphafold_multiallele_5x30_2026-08-20/`
  (5 AlphaFold 3 job batches), which exists in no tree on this machine.
