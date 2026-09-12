# TCRModel2 control-failure technical debugging report

**Purpose.** Explain what the two repaired class-II controls do—and do not—show before the Tom Resnik call. This is a structural-validation report, not evidence that either candidate peptide binds an HLA, binds a TCR, cross-reacts, activates T cells, or contributes to MS.

## Bottom line

The repaired controls still fail the predeclared gate: none of the ten ranked models is CAPRI Acceptable (DockQ >= 0.23) with nonzero native-contact recovery. After an independent, sequence-aligned audit, the main explanation is **not** wholesale pMHC or TCR-fold failure. Both rank-0 pMHCs and TCR folds are close to their deposited structures; the folded TCR is placed in the wrong location/orientation on the pMHC.

There are two important evaluation caveats, neither of which rescues the controls:

1. The existing DockQ wrapper correctly preserves the intended R/L chain mapping, but DockQ automatically assigns the longer native group as “receptor.” For 2WBJ, that reverses the requested pMHC-receptor/TCR-ligand convention and inflates the reported placement RMSD (about 55 Å). The contact metrics, iRMSD, and “Incorrect” classification nevertheless remain failed.
2. In 1YMM, the deposited peptide is a 15-mer (`ENPVVHFFKNIVTPR`) but its terminal Arg lacks resolved coordinates. The 11-mer submitted here (`VHFFKNIVTPR`) is the intended C-terminal 11-residue window; sequence alignment compares the ten coordinate-observed overlapping peptide residues. It does not compare an unobserved terminal Arg. This is a limitation to label, not a reason to claim the result is invalid.

## Diagnosis hierarchy

| Rank | Diagnosis | Evidence | What remains uncertain |
|---|---|---|---|
| 1 | **TCR docking-placement failure** | Sequence-aligned pMHC-to-native alignment leaves the rank-0 TCR 15.35 Å (1YMM) or 17.87 Å (2WBJ) from its native placement. 2WBJ recovers zero native contacts in every rank. | Whether another sampling/modeling method can recover these particular poses. |
| 2 | **Evaluation implementation caveat for 2WBJ placement RMSD** | DockQ calls the longer 2WBJ native TCR group the receptor despite `RL:RL`; therefore its LRMSD is not the requested pMHC-aligned TCR-placement statistic. | A reimplementation of DockQ with forced receptor identity would be needed for a replacement CAPRI-style LRMSD, but would not change zero contact recovery. |
| 3 | **1YMM peptide-window observability caveat** | The raw 1YMM coordinates resolve 14 of 15 deposited peptide residues and omit terminal Arg. | The exact influence of the missing terminal flank on docking cannot be measured from this deposition. |
| 4 | **No simple direct-template-copy explanation** | Excluded pMHC PDBs are absent from the returned pMHC lists; returned TCR-chain templates are not 1YMM or 2WBJ. | The public “exclude” option is not a whole-pipeline training-data holdout. |
| 5 | **Not primarily pMHC reconstruction or TCR-fold failure** | Rank-0 pMHC C-alpha RMSD: 0.97 Å (1YMM), 0.65 Å (2WBJ). TCR-fold RMSD after separate TCR alignment: 2.24 Å, 1.20 Å. | These are global matched-residue measures, not proof that every peptide side chain/register detail is right. |

## Local forensic audit

### Inputs, residues, and chain correspondence

The raw deposited structures were normalized to pMHC = `R` and TCR = `L`, then the model and native sequences were globally aligned before coordinates were compared. This preserves the existing normalizer's chain choices, including the special 2WBJ handling that separates its peptide from the fused peptide/TCR-beta source chain.

| Control | Native/model matched pMHC residues | Identity among matched | Native/model matched TCR residues | Identity among matched | Audit finding |
|---|---:|---:|---:|---:|---|
| 1YMM | 173 / 185 model residues | 100.0% | 214 / 230 | 100.0% | Model sequence maps unambiguously to native variable/domain subset. Ten submitted peptide residues have coordinates in common; deposited terminal Arg is unresolved. |
| 2WBJ | 172 / 185 model residues | 100.0% | 223 / 230 | 99.1% | Model sequence maps unambiguously to the native subset after peptide/TCR-beta separation. All 11 submitted peptide residues are represented. |

This confirms that the fresh component analysis uses the same biological chains and sequence-overlap subset in model and reference. It also identifies the limits of that comparison instead of silently treating full-length deposited chains and server-trimmed domains as identical.

### DockQ rerun: existing scorer, unchanged normalization

The existing scorer was rerun unchanged: it groups pMHC as `R`, TCR as `L`, calls DockQ with `--mapping RL:RL`, and reports `best_result`. The raw results are retained in [`analysis/dockq_calibration/tcrmodel2_grouped_tcr_pmhc_dockq.tsv`](analysis/dockq_calibration/tcrmodel2_grouped_tcr_pmhc_dockq.tsv).

| Control | Best DockQ | Native-contact recovery (best fnat) | All five CAPRI calls | Gate |
|---|---:|---:|---|---|
| CAL_1YMM_CORRECTED | 0.1358 | 0.0986 | Incorrect | Fail |
| CAL_2WBJ_CANONICAL | 0.0375 | 0.0000 | Incorrect | Fail |

For 1YMM, the DockQ placement term (15.31–18.41 Å) agrees with the independent pMHC-aligned placement calculation. For 2WBJ, the raw DockQ placement term is 54.82–55.38 Å because DockQ selects the longer native TCR group as receptor; it is not the planned directional placement statistic. Crucially, 2WBJ still has `fnat = 0` in all five ranks, and its independent pMHC-aligned TCR displacement is 14.91–17.87 Å. The control therefore fails under either treatment.

### Component decomposition

The previous component-RMSD helper used incompatible residue-position matching and produced implausibly large pMHC RMSDs. It should **not** be used. The replacement audit uses only sequence-aligned, identical-residue C-alpha pairs and reports these values:

| Control | Rank | pMHC matched C-alpha RMSD (Å) | TCR-fold RMSD after TCR alignment (Å) | TCR placement RMSD after pMHC alignment (Å) |
|---|---|---:|---:|---:|
| 1YMM | 0 | 0.971 | 2.242 | 15.347 |
| 1YMM | 1 | 0.991 | 2.259 | 16.117 |
| 1YMM | 2 | 0.662 | 1.950 | 18.444 |
| 1YMM | 3 | 0.992 | 2.120 | 15.499 |
| 1YMM | 4 | 0.972 | 1.964 | 18.341 |
| 2WBJ | 0 | 0.649 | 1.202 | 17.867 |
| 2WBJ | 1 | 0.690 | 1.898 | 14.911 |
| 2WBJ | 2 | 0.749 | 1.527 | 16.756 |
| 2WBJ | 3 | 0.696 | 1.680 | 16.206 |
| 2WBJ | 4 | 0.708 | 1.195 | 17.827 |

Interpretation: these controls make sensible pMHC-like scaffolds and roughly correct individual TCR structures, but not the correct *relative interface*. That is exactly the failure mode a single global confidence score can conceal.

### Template and leakage audit

| Job | Requested pMHC exclusion | Returned pMHC templates | Returned TCR alpha/beta templates | Result |
|---|---|---|---|---|
| CAL_1YMM_CORRECTED | 1YMM | 2WBJ, 6R0E, 1FYT, 4Y19 | 6OVN_A, 6DFX_G, 4EN3_A, 4MJI_D / 4DZB_B, 3O4L_E, 3MFG_B, 1KTK_E | 1YMM absent; no direct 1YMM TCR-chain template. |
| CAL_2WBJ_CANONICAL | 2WBJ | 6R0E, 1FYT, 4Y19, 6CQL | same returned TCR list | 2WBJ absent; no direct 2WBJ TCR-chain template. |

The source code defines `ignore_pdbs_string` as “Do not use these pdbs as pmhc templates”; it does not claim to exclude separate TCR templates or the model's learned parameters. It also defaults `max_template_date` to 2100-01-01. Therefore, these runs rule out a **direct returned pMHC-template copy** of the withheld control but are not a date-held-out, training-held-out benchmark. The 1YMM run used 2WBJ as a structurally related pMHC template: that plausibly helped the good pMHC scaffold, but it plainly did not solve TCR placement. [TCRmodel2 source](https://raw.githubusercontent.com/piercelab/tcrmodel2/main/run_tcrmodel2.py)

### Why confidence and DockQ differ

Rank-0 global pLDDT is high (89.69 for 1YMM; 90.32 for 2WBJ), while whole-model confidence is 0.8397 and 0.7960. Both are below TCRmodel2's suggested 0.85 high-specificity threshold for a likely accurate model. The primary paper reports benchmark-level association between confidence and quality—not a guarantee that a particular high-local-confidence structure has the right interface. Here, a locally well-folded pMHC and TCR are confidently assembled into a wrong complex. [Yin et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10320165/)

## What the external evidence says

TCRmodel2 itself describes the challenge as predicting TCR recognition with flexible, diverse CDR loops and broad docking orientations; its class-II input convention is an 11-mer made of a 9-residue core plus one flank on each end. [TCRmodel2 documentation](https://github.com/piercelab/tcrmodel2)

This local pattern agrees with recent independent evidence: a held-out benchmark found that pMHC quality can remain strong while full TCR-pMHC quality varies widely, identifying TCR docking orientation as a principal remaining bottleneck. It also cautions that post-hoc docking refinement is not reliably corrective. [Held-out benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC13240595/) A recent multi-tool comparison likewise identifies CDR3 loops, TCR docking orientation, TCR-peptide interfaces, and class-II MHC-peptide interfaces as continuing hard cases. [Comparative analysis](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00298)

For class II, peptide register is an additional upstream uncertainty: “11 residues” is a server input convention, not independent proof that the biological core/register is correct. CDR3 flexibility and the broad class-II docking landscape compound that uncertainty. The present control failure is therefore compatible with known limitations; it does **not** establish that TCRmodel2 is universally unusable.

## Recommended next tests

1. **Repair the evaluator before interpreting any threshold.** Keep the current normalizer and raw DockQ output as an audit record, but add a sequence-aligned directional placement metric with pMHC explicitly used as the superposition scaffold. Treat 2WBJ's raw DockQ LRMSD as non-directional/unsuitable for placement claims. Retain fnat and iRMSD as failed diagnostics.
2. **Use a strict, date-held-out class-II panel.** Preselect deposited complexes released after the server/model template cutoff, then exclude each target PDB, use unrelated TCRs and pMHCs, lock peptide register and sequences before runs, download all five ranks, and score blindly with the repaired fixed pipeline. Include a positive benchmark, a deliberately difficult but biologically valid benchmark, and negative/decoy inputs only for clearly defined questions.
3. **Use a second method only as sensitivity analysis.** A second predictor can test whether an interface is stable to modeling assumptions; agreement is not evidence that a TCR binds, and disagreement is an uncertainty result—not a vote.
4. **Keep experimental ordering explicit.** First test exact-HLA peptide binding and map the register/core. Only when an appropriate paired-TCR system exists should TCR binding and activation be tested. Structural models may help choose experiments after validation, but cannot establish cross-reactivity or mechanism.

## Questions to read to Tom

1. “Does this pattern—good pMHC and TCR folds but 15–18 Å error in pMHC-aligned TCR placement—look like the normal docking-orientation failure mode to you?”
2. “What minimum, date-held-out class-II control panel and recovery criterion would you require before a TCR-pMHC model could prioritize an EBV/self experiment?”
3. “For the 1YMM peptide window and the 2WBJ fused-chain normalization, would you change the benchmark or scoring design?”
4. “Do the different DRB1*15:01 and DRB5*01:01 alleles make the current Hy.2E11 comparison structurally useful, or should it wait for matched-HLA evidence?”
5. “If the next control panel is still inconclusive, what is the most informative experimental next step: peptide-HLA binding/register mapping, or a paired-TCR assay?”

## Audit trail

- Archives and five ranked PDBs per job: [`results/`](results/)
- Job URLs, input peptides, exclusions, and returned templates: [`SUBMISSION_MANIFEST.tsv`](SUBMISSION_MANIFEST.tsv)
- Unchanged-normalizer DockQ outputs: [`analysis/dockq_calibration/`](analysis/dockq_calibration/)
- Rerunnable sequence-aligned component audit: [`analysis/audit_tcrmodel2_components.py`](analysis/audit_tcrmodel2_components.py)
- Peptide-window/register-boundary audit: [`analysis/CONTROL_PEPTIDE_WINDOW_AUDIT.md`](analysis/CONTROL_PEPTIDE_WINDOW_AUDIT.md)
- Original submission/analysis rules: [`SUBMISSION_AND_ANALYSIS_PROTOCOL.md`](SUBMISSION_AND_ANALYSIS_PROTOCOL.md)
- Archive checksums: [`ARCHIVE_SHA256SUMS.txt`](ARCHIVE_SHA256SUMS.txt)
