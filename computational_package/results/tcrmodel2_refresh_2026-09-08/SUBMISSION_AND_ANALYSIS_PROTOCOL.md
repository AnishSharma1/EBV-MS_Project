# Submission and analysis protocol

## Submission

Use the official TCRModel2 **TCR-pMHC complex / MHC Class II** workflow at [tcrmodel.ibbr.umd.edu](https://tcrmodel.ibbr.umd.edu/).

For each row in `SUBMISSION_MANIFEST.tsv`:

1. Upload the named five-record FASTA file.
2. Enter the exact comma-separated template-exclusion list from `pMHC_templates_to_exclude`.
3. Leave Amber minimization off so all four jobs use the same server setting.
4. Save the job URL in `server_job_url` and, when complete, copy the server-returned pMHC template list into `returned_template_list`.
5. Download the archive into `results/<submission_id>/` without renaming the server files. Record that directory in `result_archive`.

## Intake checks

- Match each uploaded FASTA to `SHA256SUMS.txt` before upload.
- Confirm each peptide is 11 residues and the intended excluded template does not appear in the returned pMHC template list.
- Confirm the server archive contains five ranked structures and its input/template report.

## Control scoring

1. Stage the downloaded controls under the existing `tcrmodel2_results` directory as `CAL_1YMM` and `CAL_2WBJ_EXCLUDED`, retaining the server's five `ranked_*.pdb` filenames.
2. Run the existing `score_tcrmodel_docking_quality.py` script unchanged, using DockQ 2.1.3. It normalizes the reference structures and scores pMHC as receptor and TCR as ligand with mapping `RL:RL`.
3. Re-run `compare_tcrmodel2_calibrators.py` for pMHC-aligned TCR-placement RMSD. Do not use whole-complex RMSD as the pass/fail criterion.
4. A repaired control passes only if at least one model has DockQ >= 0.23 (CAPRI Acceptable) and `fnat > 0`. Both controls must pass before any Hy.2E11 geometry is discussed beyond an exploratory structural hypothesis.

## Exploratory Hy.2E11 analysis

- Use the existing contact-recurrence script on one MBP/DRB1*15:01 and one BALF5/DRB5*01:01 job.
- Report contact recurrence separately from calibration. The two conditions differ in their HLA beta chain, so contact overlap is not a same-allele cross-reactivity test.
- Do not infer TCR binding, cross-reactivity, activation, molecular mimicry, or MS mechanism.

## Current execution state

The uploads are prepared but not submitted: the local browser's enforced security policy blocked access to the TCRModel2 server. No job IDs, server templates, downloaded models, or scoring outputs exist yet.
