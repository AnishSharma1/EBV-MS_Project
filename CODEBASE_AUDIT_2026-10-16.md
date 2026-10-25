# EBV–MS codebase audit — 16 Oct 2026

Scope: all four workspaces under `EBV-MS_Molecular_Mimicry_Project/`, with
`EBV-MS_Project/` treated as the canonical checkout (HEAD `3d1dfbe`, two commits
ahead of `EBV-MS_Project_repo` at `fa58f8f`).

Verdict in one line: **the science is sound and independently reproducible; the
documentation around it is not.** Every headline number I could recompute from
frozen data came out exactly right. The failures are all in the README/manifest
layer — documentation describing analyses that were dropped, artifacts that are
absent, and vocabulary that exists nowhere but the README.

---

## A. What I verified as correct

These were checked against the frozen data, not taken on trust.

| Claim | Check | Result |
|---|---|---|
| 6,400 within-allele pairs, 1,600 per allele × 4 alleles | counted `within_allele_pair_universe_6400.csv` | **exact** |
| 40 EBV + 40 human self peptides | counted `frozen_v2_80_peptide_panel.csv` `kingdom` | **exact** (40/40) |
| DR13 lead score 0.28125, rank 13/1,600 | recomputed BLOSUM62 at P2/P3/P5/P7/P8 from declared cores `YHFVKKHVH` vs `LSFDKDAMV` → 9/32 | **exact** |
| DR15 lead score 0.31429, rank 14/1,600 | recomputed from `VYHFVKKHV` vs `IYNYYKKFS` → 11/35 | **exact** |
| "Both pairs had 40.0% identity at the five TCR-exposed positions" | recomputed | **exact** (2/5 both) |
| "full-core identity ~22.2% in \*13:03 and 33.3% in \*15:01" | recomputed (2/9, 3/9) | **exact** |
| Table S4 Benjamini–Hochberg q-values | recomputed step-up BH over the 8 p-values | **exact** to 5 d.p. |
| "No audited pair had BH q < 0.05"; "both leads maxT p = 1.0" | read Table S4 | **correct** |
| IEDB Assay 1774674 (BALF5 627–641 / DRB1\*15:01) | traced to raw API record with PMID 17559580, Hansen et al., *Tissue Antigens* 2007 | **fully provenanced** |
| Thomas et al., ANO2 T cells, *Cell* 2026 | verified against Karolinska/KI publication record | **real; author initials correct** |
| PANDORA v2.0 reference | verified: *Front Immunol* 2023;14:1285899 | **exact** |
| `verify_package.py` | ran it | **PASSES** (213 manifest checksums) |
| TALDO1 coordinates, UniProt P37837, 337 aa | cross-checked panel entries at 1–15 / 108–122 / 216–230 / 323–337 | **internally consistent** |

`computational_package/docs/SCOPE_AND_LIMITS.md` is honest and well-calibrated —
it explicitly declines to promise one-command regeneration. The manuscript's own
claim boundaries (`main.tex` §limitations) are disciplined.

---

## B. Critical — fix before anyone else reads the repo

### B1. The root README describes analyses the manuscript explicitly disclaims

`README.md:25` and `README.md:47` state the screen used:

> "matched empirical decoys, proteome-wide sequence rarity filters (**11.25 million 9-mers**)"
> "proteome-wide null background screens (11.25M 9-mers)"

`manuscript/main.tex:168` (the V19 limitations paragraph) says the opposite:

> "HLA-DQ and HLA-DP loci, broader haplotype effects, and **proteome-wide background sampling were omitted**."

The 11.25M claim appears in the V14–V17 supplements only; it was **removed from
V18 and V19**. `grep "11.25" manuscript/supplement.tex` returns zero. The README
is still advertising the V14-era method set.

This is the single most damaging item: a reviewer comparing README to paper sees
the repository claiming a proteome-wide null the paper says it did not do.

### B2. The seven "round-3 computational additions" packages do not exist

`code/revision_analyses/build_release_index.py` declares seven packages under
`output/revision_round_3/` and writes `"release_status": "PASS"` plus prose
asserting "All seven component packages have independent checksum manifests."

Searched the entire workspace, all four trees:

```
computational_tightening_2026-09-15        -> 0 directories found
proteome_sequence_null_2026-09-15          -> 0
proteome_binding_consensus_2026-09-15      -> 0
register_uncertainty_2026-09-15            -> 0
solved_structure_benchmark_inventory_...   -> 0
evidence_and_control_audit_2026-09-15      -> 0
proteome_binding_filter_2026-09-15         -> 0
```

`output/revision_round_3/` itself does not exist. The script would crash at
`verify_manifest()` on the first missing `SHA256SUMS.csv`. Its hardcoded
`requested_additions_status.csv` — which marks the proteome null, matched decoys,
corrected structural ensembles and AF3 confidence tables all "complete" — is
therefore a **table of assertions with no artifacts behind it**. This is the
provenance record most likely to be mistaken for evidence.

Note this does **not** compromise the manuscript: V19 does not cite these
results. It compromises the repository's self-description.

### B3. Figures cannot be regenerated anywhere, including on your own machine

Every figure script in `code/figures/` reads or writes outside the repo:

- `build_workflow_and_rank_sensitivity.R:1`, `build_structural_and_binding.R:1`,
  `restyle_figures_ggplot.R:2`, `build_workflow_and_rank_sensitivity.py:5` →
  write to `/Users/anishsharma/.../CloudDocs/ebv_ms_publication/manuscript_figures_v2`
- `build_panels.py:25-26` → read `/Users/anishsharma/Documents/New project/...`
- `restyle_figures_ggplot.R:6`, `build_structural_and_binding.R:23` → read
  relative `processed/...` paths that exist in neither checkout

I checked the external workspace. The three input directories those scripts need
are **also missing there**:

```
MISSING  high_yield_candidate_evidence_2026-08-28
MISSING  literature_grounded_hla2_rankings_v3_2026-08-27
MISSING  pmhc_surface_electrostatics_v2_controls_2026-08-30
```

So the published figures are currently unregenerable by anyone, you included.

### B4. `EBV-MS_Project_repo/` is a regression trap

`git status` there shows eight modified files including `manuscript/main.tex`.
Those modifications look like newer work. They are **older**:

| | `EBV-MS_Project` (canonical) | `EBV-MS_Project_repo` (working tree) |
|---|---|---|
| Date | September 21, 2026 | September 14, 2026 |
| Title | "Register-aware computation prioritizes EBV–self peptide pairs…" | "Computational assessment of EBV–self peptide similarities…" |
| Abstract | single-paragraph, V19 | BMC Background/Results/Conclusions, pre-V18 |

Committing that working tree would silently revert the manuscript by five
versions. Recommend discarding it or stashing it with a loud name.

Related: the root `README.md` BibTeX block and `CITATION.cff` both still carry
the **old title** from this stale copy. Three titles are in circulation —
README/BibTeX, CITATION.cff/manuscript README, and the actual V19 `\title{}`.

---

## C. Documentation that asserts things that do not exist

| Claim | Location | Reality |
|---|---|---|
| `figureS3.pdf` — "Binding-threshold sensitivity" | `figures/README.md:27` | Does not exist. Only S1, S2. (`manuscript/README.md` was corrected for this in `3d1dfbe`; `figures/README.md` was not.) |
| `figures/panels_v2/`, `figures/exploratory/` | `figures/README.md:31,37` | Neither directory exists. |
| `manuscript/supplement_tables/` with Tables S1–S7 | `manuscript/README.md` | Directory does not exist; the supplement actually renders **S1–S8** (8 tables). |
| Main figures are "figure1–figure4" with given captions | `manuscript/README.md`, `figures/README.md` | `main.tex` includes only **figure1, figure4, figure5** (rendering as Figures 1, 2, 3). figure2/3/6 are orphans. The README captions describe a figure set the paper no longer uses. |
| `outputs/README.md` — "Output catalog & evidence ledger" | root `README.md` §5 | Missing. |
| Routing tiers `computational_stratum_A` / `computational_stratum_B` | root `README.md:57-58` | **Invented vocabulary.** The string appears in no script, no data file, and no manuscript version. The data use `weak_medium_provisional` / `strong_medium_provisional`; the manuscript uses "experimental-routing class 1/2/3". Both leads carry identical `stage1_medium_priority` status, so the A/B distinction is unsupported as well as unsourced. |
| "210+ analysis scripts" | `code/README.md:3` | 200 (190 `.py` + 10 `.R`). |
| "16 Computational tightening… analyses" | `code/README.md:31` | 21 files in `code/revision_analyses/`. |
| "deposited to Zenodo" | `.gitignore:` doi_release comment | No DOI exists. `BMC_SUBMISSION_FIELDS.md` says one is "recommended"; `ZENODO_RELEASE_CHECKLIST.md` is still a to-do list. |
| "IEDB assay 1774674 original Hansen 2007 publication remains unresolved" | `manuscript/CITATION_NOTES.md:3` | Stale — it *was* resolved. `iedb_iq_api_raw_records.jsonl` carries PMID 17559580 with full author list. |

---

## D. Broken validation

`CLAUDE.md` instructs running two validation commands. One of them fails.

```
python3 code/tests/test_repository_integrity.py
  AssertionError: PROJECT_DIRECTORY_MAP.md missing
```

The test was copied from the `ebv_ms_publication` workspace at the CloudDocs
root without adaptation. It asserts four paths that exist **there** but not here
(`PROJECT_DIRECTORY_MAP.md`, `ANTIGRAVITY_PROJECT_BRIEF.md`,
`Model_Library/STRUCTURAL_MODELS_INVENTORY.md`, `outputs/README.md`), and its
`test_no_forbidden_root_clutter()` asserts `Rplots.pdf` is absent from the root —
it is present at both `./Rplots.pdf` and `./manuscript/Rplots.pdf`.

Its ROOT-discovery loop also searches for `ANTIGRAVITY_PROJECT_BRIEF.md`, never
finds it, and silently falls through to `parents[2]`. That happens to be correct
here, but only by accident.

Compounding this: **no test runner and no dependency manifest.** No
`requirements.txt`, `pyproject.toml`, `environment.yml`, or `renv.lock` anywhere
in the tree. Neither available interpreter (`python3.14`, `/usr/bin/python3.9`)
has `pytest`, `numpy`, or `pandas`. The 49 test files cannot be run as a suite.
For the R scripts, `R-ENV-001` applies — no `renv.lock`, and the figure scripts
hardcode absolute output paths rather than using `here()`.

---

## E. Provenance hazards (lower severity, worth fixing)

1. **Duplicate, divergent figure1.** `manuscript/figures/figure1.pdf` and
   `figures/manuscript/figure1.pdf` have **different SHA-256**. figures 2–6 are
   byte-identical between the two locations. The paper compiles against
   `manuscript/figures/`, so `figures/manuscript/figure1.pdf` — the copy the
   "publication figures" README points at — is the stale one.

2. **Two core definitions for the same peptide arm.**
   `frozen_v2_80_peptide_panel.csv` gives TALDO1 216–230 `proposed_core =
   KIYNYYKKF`; the V3 declared core used for ranking is `IYNYYKKFS` (one residue
   left). Scoring under the panel's value gives 0.11429, not 0.31429. The
   manuscript uses the correct (declared) core and documents core freezing — but
   a reader auditing from the frozen panel file will reproduce the wrong number.
   Worth an explicit note in the panel README distinguishing the two columns.

3. **`output/` and `outputs/` both exist.** Only `outputs/` appears in the root
   README's directory map, yet `output/` holds revision rounds 4, 11 and 18 —
   the results the V19 manuscript actually depends on.

4. **Absolute paths everywhere.** 3,803 references to
   `/Users/anishsharma/Library/...` outside `.git`, plus `~/.cache/ebv_ms_tools/`
   solver paths and 86 references to
   `/Users/anishsharma/Documents/Codex/2026-08-01/i-al/work/...`, which **no
   longer exists**. `output/revision_round_4/threshold_sensitivity_2026-09-20/analysis_lock.json`
   pins an input by absolute path into `ebv_ms_publication/processed 2/` — a
   Finder-duplicate directory outside every checkout.

5. **Supplement label `tab:s2a` renders as Table S7.** Numbering is correct in
   the PDF (S1–S8), but the label name implies a position it does not occupy.
   Rename to `tab:s7` to avoid a future miscitation.

6. **Target journal is ambiguous.** Root README, `manuscript/README.md`,
   `BMC_SUBMISSION_FIELDS.md` and `ZENODO_RELEASE_CHECKLIST.md` all say *BMC
   Bioinformatics*. The newest cover letter,
   `scientific_reports_cover_letter_v19.md` (Sept 21), is addressed to
   *Scientific Reports*. `supplement.tex` carries the Scientific-Reports-style
   title. Nothing records which is current.

7. **DOI release is one manuscript version stale.** `doi_release/` is
   `v1.1.0` built around `paper_v18`; `CITATION.cff` (v1.1.0, 2026-09-20)
   describes "V18 sensitivity analyses". The manuscript is V19.

8. **Citation version gap.** MixMHC2pred **2.1** is used throughout, but
   `references.tex` cites only Racle et al. 2019 (*Nat Biotechnol*), which is
   v1. The v2 methods paper is not cited.

---

## F. What I did not verify

Stated plainly so this is not mistaken for a clean bill of health:

- **TALDO1 108–122 and 216–230 sequences were not checked against UniProt
  P37837 directly.** No proteome FASTA is stored in the repo and
  `reference_sequence_manifest.csv` records only lengths and SHA-256, not the
  residues. Internal consistency is good (length 337 matches P37837, the 1–15
  and 323–337 entries look like genuine termini), but the two lead self-peptide
  coordinates remain unconfirmed against the primary record. **This is the one
  remaining factual risk in the lead candidates and should be closed first** —
  if either coordinate is off, both leads change.
- The 49-pair panel arithmetic (49 source rows → 47 unique → 2 replacements) was
  not recomputed from the panel file.
- AlphaFold 3 RMSD values behind Figure 3 were not recomputed from coordinates.
- PDF rendering of main/supplement was not re-run (no XeLaTeX check performed);
  the existing `.log` files show no undefined references or citations.
- Remaining references in `references.tex` beyond the four I checked were not
  individually verified against primary sources.

---

## G. Suggested order of work

1. Verify TALDO1 108–122 / 216–230 against UniProt P37837 and store the FASTA in
   the package. (Closes the only open scientific risk.)
2. Rewrite root `README.md` §1 and §2 against V19 — delete the proteome-null and
   matched-decoy claims, delete `computational_stratum_A/B`, fix the title and
   BibTeX to match `main.tex`.
3. Either restore `output/revision_round_3/` or delete
   `build_release_index.py` and its "PASS" assertions.
4. Fix `code/tests/test_repository_integrity.py` to assert against this repo;
   add `requirements.txt` (or `pyproject.toml`) and `renv.lock`.
5. Correct `figures/README.md` and `manuscript/README.md` figure/table
   inventories; remove or regenerate the orphan figure2/3/6 and the divergent
   `figures/manuscript/figure1.pdf`.
6. Decide the target journal and make every document agree.
7. Discard or clearly archive the `EBV-MS_Project_repo` working tree.
