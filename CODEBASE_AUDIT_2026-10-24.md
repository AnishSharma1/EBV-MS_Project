# EBV–MS codebase audit — 24 Oct 2026

Follow-up to `CODEBASE_AUDIT_2026-10-16.md`, eight days on. Scope: the canonical
`EBV-MS_Project/` checkout, with the three sibling workspaces inspected where
they bear on it.

Verdict in one line: **the science is in good shape and the headline validation
now passes; the gap was reproducibility infrastructure, and most of it is now
closed.** The test suite went from 3 of 51 passing to 48, and the cause of the
failures was never analysis defects — it was import paths left behind by the
`code/` reorganisation.

---

## A. What the 16 Oct audit flagged and is now fixed

Confirmed by running the checks, not by reading the diff.

| Item | Status |
|---|---|
| TALDO1 108–122 / 216–230 unverified against UniProt P37837 (the one open scientific risk) | **closed.** `computational_package/data/reference_sequences_2026-10-16/` stores the FASTA, retrieval provenance, and a coordinate verification table; `code/tests/test_peptide_coordinate_integrity.py` re-checks every panel peptide against it (8 checks, passing) |
| README claimed an 11.25M-9-mer proteome null the manuscript disclaims | **removed** (0 occurrences in `README.md`) |
| `computational_stratum_A/B` invented vocabulary | **removed** (0 occurrences) |
| `figures/README.md` claimed a non-existent `figureS3.pdf` | **corrected**, now states explicitly that no such figure exists |
| `manuscript/figures/figure1.pdf` diverged from `figures/manuscript/figure1.pdf` | **resolved**, byte-identical (`d561c526…`) |
| `test_repository_integrity.py` asserted paths from another workspace | **rewritten** to this repo; 9/9 passing |
| `Rplots.pdf` in the repository root | **removed** (one copy remains at `manuscript/Rplots.pdf`, gitignored) |
| MixMHC2pred v2 used but only the 2019 v1 paper cited | **fixed** in V20 (`mix2023`) |
| Manifest did not cover the new reference sequences | **covered**, 213 → 219 checksums, `verify_package.py` passing |

V20 itself is well handled: `manuscript/V20_CHANGES.md` records every text change
with its evidence, states plainly that no result, figure, table, or dataset
changed, and documents that the snapshots were built with `tectonic` 0.17.0
rather than the XeLaTeX named in the README.

---

## B. Two findings in the 16 Oct audit that were wrong

Recorded because both overstated the damage.

**B1. The test suite was described as unrunnable for want of `pytest`.** It was
runnable; it was misconfigured. 44 of 51 modules import nothing beyond the
standard library. The real cause was that tests import their targets by bare
name (`import build_rmsd_rerankings`) while the targets live in
`code/pipeline/` — moved there by `2f5e26f` and `fa58f8f`. 23 modules carried a
`sys.path` line pointing at `code/src`, a directory this layout never had, and
22 more had no bootstrap at all. Setting `PYTHONPATH` alone took the suite from
3 passing to 19.

**B2. "The published figures are currently unregenerable by anyone, you
included."** The three input directories that audit reported as missing
everywhere are present in `ebv_ms_publication/processed/`. Verified file by
file, including `v3_all_hla_ranked_pairs.csv` (13 MB) and
`frozen_candidate_registry.csv`. The figures are regenerable; their inputs are
merely outside the repository and unversioned, which is a different and smaller
problem — see section C.

---

## C. The central remaining gap: the external `processed/` tier

51 scripts under `code/`, nine test modules, and every figure script read
`ebv_ms_publication/processed/` — **8.0 GB across 57 directories, sitting beside
this checkout, under no version control and with no checksum manifest.** It is
the single largest unprotected dependency in the project.

Nothing read from it is covered by `MANIFEST.sha256` or `verify_package.py`.
Manuscript results were computed from the frozen inputs inside
`computational_package/`, so this does not affect auditing the reported numbers;
it affects anyone re-running the analysis or figure scripts.

Two further facts about it:

1. **Two copies of `ebv_ms_publication` exist on the authoring machine** — one
   beside this checkout (15 GB) and one directly under CloudDocs (18 GB) — and
   they have diverged in at least one file:
   `processed/tcell_library_v2_2026-08-22/model_inventory_320.csv` is present in
   the first and absent from the second. Several
   `code/revision_analyses/` scripts hardcode the second path. **Which tree is
   authoritative has not been established, and this audit did not decide it.**
2. `computational_package/docs/SCOPE_AND_LIMITS.md` honestly declines to promise
   regeneration and names the external *systems* (IEDB, RCSB, AlphaFold Server,
   PANDORA, APBS, TCRmodel2), but it does not name this directory. A draft
   paragraph documenting it was written and then reverted, because that file is
   inside the checksummed package and its hash is also recorded in
   `logs/package_verification/MANIFEST.sha256` and in the frozen `doi_release/`
   v1.1.0 staging. Re-hashing it belongs in a deliberate package update, not in
   a test-infrastructure commit. **Recommended: fold it into the same update
   that refreshes `CITATION.cff` and `doi_release/` to V20.** The dependency is
   documented meanwhile in `README.md` §6 and `code/tests/README.md`.

---

## D. What this commit changed

Test and infrastructure files only. No analysis logic, no data, no manuscript
text, no frozen package contents.

- **23 test modules**: `parents[1] / "src"` → `parents[2] / "code" / "pipeline"`.
- **22 test modules**: added the same bootstrap, which they lacked entirely.
- `code/tests/test_multiallele_validator.py`: validator path corrected to
  `code/batch_builders/`, where the script now lives.
- `code/pipeline/build_hla2_positive_control_benchmark_v2.py`: the two
  `geometry_module_sha256` / `v2_contract_module_sha256` provenance fields hashed
  `ROOT / "src/…"` and so could not be computed at all; repointed to
  `ROOT / "pipeline/…"`. This is the only edit outside `code/tests/`, and it
  restores a provenance field rather than changing a result.
- **New** `requirements.txt`: numpy 2.5.4, pandas 3.1.0, scipy 1.18.1,
  matplotlib 3.11.2, scikit-learn 1.9.1, biopython 1.88, requests 2.34.2,
  pytest 9.1.1. Pinned exactly, and resolved by installing them rather than
  recalled from memory; verified on CPython 3.14.7.
- **New** `code/tests/run_suite.py`: runs all 51 modules, reports PASS/FAIL/SKIP,
  exits non-zero only on a real failure. Locates the external tier and links it
  to `code/processed`; declares the tier-dependent and known-blocked modules
  explicitly so skips stay visible instead of being absorbed into a pass count.
- **New** `code/tests/README.md`; `README.md` §6 gains a test-suite and
  external-tier section; `.gitignore` ignores `code/processed`.
- `.github/workflows/verify_package.yml`: previously ran only
  `verify_package.py`, the weakest gate. Now runs repository integrity and
  peptide coordinate integrity too, then installs `requirements.txt` and runs the
  suite — so the coordinate guarantee closed on 16 Oct can no longer regress
  silently.

### Suite status, measured

| | result |
|---|---|
| before | 3 passed, 48 failed |
| from this repository alone | 39 passed, 0 failed, 12 skipped |
| with the external tier reachable | **48 passed, 0 failed, 3 skipped** |

The three skips are inputs absent from every tree on this machine, recorded in
`MISSING_UPSTREAM_DATA` in `run_suite.py`: no AlphaFold 3 model root for
`test_audit_complete_model_pipeline`; the diverged-tree file above for
`test_build_missing_v2_af3_batches`; and
`outputs/alphafold_multiallele_5x30_2026-08-20/` (5 AF3 job batches) for
`test_multiallele_validator`, which exists nowhere on this machine.

---

## E. Still open

Not addressed here, in rough priority order.

1. **`code/revision_analyses/build_release_index.py`** still declares seven
   packages under `output/revision_round_3/`, a directory that does not exist,
   and still writes `"release_status": "PASS"`. Its hardcoded
   `requested_additions_status.csv` marks the proteome null, matched decoys and
   AF3 confidence tables "complete" with no artifacts behind them. V19/V20 cite
   none of it, so no manuscript claim depends on it — but it is the file in this
   repository most likely to be mistaken for evidence. Delete it or restore the
   artifacts.
2. **Target journal undecided.** `README.md` and `manuscript/README.md` say BMC
   Bioinformatics; `manuscript/V20_CHANGES.md` and
   `scientific_reports_cover_letter_v19.md` say Scientific Reports. V20 states
   the decision was deliberately deferred. Nothing records which is current.
3. **`CITATION.cff` and `doi_release/` are two manuscript versions stale** —
   both v1.1.0 / 2026-09-20 / V18-era against a V20 manuscript. No DOI is minted.
   Bundle with the `SCOPE_AND_LIMITS.md` update from section C2.
4. **`EBV-MS_Project_repo/` regression trap is still live**: 10 modified files,
   `main.tex` title five versions stale
   ("Computational assessment of EBV-self peptide similarities…" against the
   current "Register-aware computation prioritizes…"). Committing that tree would
   silently revert the manuscript. Archive or discard it.
5. **Figure scripts still write outside the repository** to absolute
   `/Users/anishsharma/…/ebv_ms_publication/manuscript_figures_v2` paths, and
   `code/figures/build_panels.py` reads from `~/Documents/New project/…`. The
   inputs exist (section B2), but no script can run on another machine. The R
   scripts have no `renv.lock`.
6. **Orphan figures**: `figure2.pdf`, `figure3.pdf`, `figure6.pdf` are present
   and unused — `main.tex` includes only figure1, figure4, figure5.
7. **Minor**: supplement label `tab:s2a` renders as Table S7; both `output/` and
   `outputs/` exist, with only `outputs/` in the README's directory map while
   `output/` holds the rounds V19/V20 depend on.

Three scientific items from the 16 Oct audit's own "not verified" list remain
unverified and were not attempted here: the 49-pair panel arithmetic
(49 source rows → 47 unique → 2 replacements), the AlphaFold 3 RMSD values
behind Figure 3, and the references in `references.tex` beyond the five checked.

---

## F. What I did not verify

Stated plainly so this is not mistaken for a clean bill of health.

- **No PDF was rebuilt.** The three gates and the suite were run; XeLaTeX is not
  installed in this environment, and V20's own note that fonts and layout were
  never compared page-by-page with V19 still stands.
- **The 27 dependency-blocked modules now pass, but only on CPython 3.14.7 with
  the pinned versions.** pandas 3.1.0 is a major version ahead of what much of
  this code was presumably written against; the suite passing is evidence it is
  compatible, not proof across other versions. CI will now exercise this.
- **The external tier's contents were not checksummed or compared between the
  two `ebv_ms_publication` trees** beyond the one diverging file found. A full
  comparison of 8 GB against 18 GB was not attempted.
- **No analysis was re-run and no result recomputed.** This audit verified
  infrastructure. The scientific recomputations in the 16 Oct audit's section A
  were not repeated.
