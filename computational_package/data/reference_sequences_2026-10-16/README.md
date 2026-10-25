# UniProt reference sequences (retrieved 2026-10-16)

Canonical UniProt records for every accession cited by
`data/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv`, stored so that
the panel's declared 1-based peptide coordinates can be checked offline.

## Contents

| File | Description |
|---|---|
| `uniprot_reference_proteins.fasta` | 31 canonical UniProtKB records (20 EBV, 11 human). Headers retain the UniProt sequence version (`SV=`). |
| `coordinate_verification.csv` | Per-peptide result: declared sequence, reference slice at the declared coordinates, and match flag. |
| `retrieval_provenance.json` | Endpoint, query, requested/returned accessions, retrieval date, and match counts. |
| `symbol_corrections.csv` | Panel protein symbols that disagree with UniProt for the same accession (currently BFRF3 for P03187), recorded without editing the frozen panel. |
| `SHA256SUMS.csv` | Checksums for the files above. |

## Provenance

Retrieved from the UniProtKB REST API (`https://rest.uniprot.org/uniprotkb/stream`)
in a single streamed FASTA query on 2026-10-16 by
`code/pipeline/fetch_reference_sequences.py`. Regenerate with:

```sh
python3 code/pipeline/fetch_reference_sequences.py          # re-fetch and rewrite
python3 code/pipeline/fetch_reference_sequences.py --check  # verify, no network
```

Refreshing these files requires a corresponding documented update to
`computational_package/MANIFEST.sha256`.

## Verification result

**72 of 80 panel peptides match their declared coordinates exactly** against the
canonical record.

Both routed BALF5–TALDO1 leads verify exactly, including their declared P1–P9
scoring cores:

| Arm | Accession | Declared | Canonical slice | Core | Core offset |
|---|---|---|---|---|---|
| BALF5 627–641 | P03198 (1,015 aa, SV=1) | `TGGVYHFVKKHVHES` | identical | `YHFVKKHVH` / `VYHFVKKHV` | 5 / 4 |
| TALDO1 108–122 | P37837 (337 aa, SV=2) | `DARLSFDKDAMVARA` | identical | `LSFDKDAMV` | 4 |
| TALDO1 216–230 | P37837 (337 aa, SV=2) | `SVTKIYNYYKKFSYK` | identical | `IYNYYKKFS` | 5 |

TALDO1 1–15 and 323–337 also match, confirming the full 337-residue span.

## The eight non-matching peptides

None is a fabricated sequence. Six are present verbatim in their own reference at
a different offset (a numbering convention); two are deliberately retained strain
variants. None belongs to a routed lead. They are listed as explicit exceptions in
`code/tests/test_peptide_coordinate_integrity.py`, so any *new* mismatch fails.

| Protein | Declared | Peptide | Canonical position | Cause |
|---|---|---|---|---|
| MBP (P02686) | 84–98 | `ENPVVHFFKNIVTPR` | 217–231 | 18.5 kDa MBP numbering (offset +133) |
| MBP (P02686) | 84–97 | `ENPVVHFFKNIVTP` | 217–230 | same |
| MBP (P02686) | 84–100 | `ENPVVHFFKNIVTPRTP` | 217–233 | same |
| MBP (P02686) | 86–100 | `PVVHFFKNIVTPRTP` | 219–233 | same |
| MOG (Q16653) | 98–108 | `FFRDHSYQEEA` | 128–138 | mature-protein numbering (offset +30) |
| EBNA2 (P12978) | 139–153 | `MANYIVRQSRGDRGL` | 136–150 | strain offset (−3) in B95-8 |
| LMP2 (P13285) | 72–86 | `DYQPLGTQDPSLYLG` | absent | study variant; B95-8 reads `SDYQPLGTQDQSLYL` |
| LMP2 (P13285) | 73–87 | `DYQPLGNQDPSLYLG` | absent | study variant; B95-8 reads `DYQPLGTQDQSLYLG` |

### Known inconsistencies these expose

1. **`reference_strain_or_isoform` is wrong for MBP and MOG.** Those rows read
   "canonical reviewed UniProt entry", but their coordinates are *not* indexed
   against the canonical entry. The field should name the isoform and numbering
   convention actually used.

2. **MBP is numbered three ways across the project**, with no cross-walk:
   - `manuscript/main.tex` and `supplement.tex` §SM8: `MBP 85-99`
     (18.5 kDa literature numbering for the Hy.2E11 control)
   - `manuscript/supplement.tex` line 138: `MBP 217-233`
     (canonical P02686 numbering — correct for `ENPVVHFFKNIVTPRTP`)
   - this panel: `84-98` / `84-100` (18.5 kDa numbering, one lower than the
     literature's 85–99)

   All three denote the same region. Pick one convention, state it, and record
   the mapping.

3. **The EBV arm had no stored reference sequence before this directory
   existed.** `data/tcell_library_v2_2026-08-22/reference_sequence_manifest.csv`
   covers 11 human proteins and 5 HLA chains — none of the 20 EBV accessions.
   The viral coordinates were therefore previously unverifiable offline.

4. **`protein_symbol` is wrong for P03187.** The panel row `EBV_IEDB_1bcd6d24071a` (`IMNILRIYYSPSIM`, 289-302) is labelled `BFRF3`, but UniProt P03187 is the 364-residue triplex capsid protein 1 (gene names TRX1/BORF1). UniProt's BFRF3 is P14348, a different 176-residue protein that does not contain this peptide. The accession, coordinates, and sequence are correct, and the panel's `protein_name` column correctly reads "Triplex capsid protein 1"; only the symbol is wrong. It originates in a hardcoded accession-to-symbol map in `code/pipeline/build_tcell_library_v2.py` (line 105). The manuscript text never uses "BFRF3", so no manuscript claim is affected. The label was **not** changed in the panel, builder, or derived tables: the panel's SHA-256 is pinned in `../high_yield_control_validation_2026-08-28/protocol_lock.json` (`frozen_before_control_geometry_read: true`) and in package manifests, so relabeling would invalidate a pre-analysis lock. The correction is recorded in `symbol_corrections.csv` (use BORF1/TRX1 for P03187) and checked by `code/tests/test_peptide_coordinate_integrity.py`.

5. **A manuscript coordinate label was wrong and has been corrected.** `supplement.tex` (Tables S1, S3) named the viral arm of the rank-1 DRB1*08:01 and DRB1*15:01 pairs "EBNA1 87-101". The peptide in the data is `AEGLRALLARSHVER`, which is residues **482-496** of UniProt P03211 (declared 482-496 in the frozen panel and in the stage-1 evidence table); P03211 87-101 is a Gly-Ala repeat (`GGTGAGAGAGGAGAG`) that is not a panel peptide. All four occurrences refer to this same peptide, including the class-2 EBNA1-MBP 217-233 pair (`output/revision_round_4/threshold_sensitivity_2026-09-20/threshold_grid_pair_results.csv`). The label was corrected to 482-496 in `supplement.tex` and `manuscript/main_bmc_bioinformatics.tex` on 2026-10-17. It persists in `main_google_docs_aligned.tex` and the `paper_v11`-`paper_v19` snapshots, which are provenance copies and were not edited. `code/tests/test_peptide_coordinate_integrity.py` now requires every `\pep{PROTEIN}{start-end}` label in `main.tex` and `supplement.tex` to resolve to a panel peptide.

## Numbering crosswalk (checked 2026-10-17)

The same sequence carries different coordinates in the frozen panel, the manuscript, and the canonical UniProt record. All coordinates below were checked against `uniprot_reference_proteins.fasta`; the frozen panel itself was not edited, because it is part of a checksummed frozen release.

| Sequence | Panel label | Canonical UniProt | Manuscript label | Comment |
|---|---|---|---|---|
| MBP `ENPVVHFFKNIVTPR` | 84-98 | P02686 217-231 | "MBP 85-99" (`main.tex` methods; `supplement.tex` SM8) | Same 15-mer. It is also the peptide chain of PDB 1BX2 (checked against the RCSB sequence record 2026-10-17). The panel label is one lower than the literature label. Offset to canonical: +133 from the panel label, +132 from the literature label. |
| MBP `ENPVVHFFKNIVTPRTP` | 84-100 | P02686 217-233 | "MBP 217-233" (`supplement.tex`, Table S3) | The manuscript uses canonical numbering here. |
| MOG `FFRDHSYQEEA` | 98-108 | Q16653 128-138 | not labelled in the manuscript | Mature-protein numbering in the panel (offset +30). |
| MOG `RTELLKDAIGEGK` | 97-109 | Q16653 97-109 | "MOG 97-109" (`supplement.tex`) | No offset; panel and manuscript agree. |
| EBNA2 `MANYIVRQSRGDRGL` | 139-153 | P12978 136-150 | "EBNA2 136-150" (`supplement.tex`) | The manuscript uses the B95-8 coordinates; the panel label is +3. |

The panel field `reference_strain_or_isoform` reads "canonical reviewed UniProt entry" for the MBP and MOG rows, which is inaccurate for the four MBP rows and the MOG 98-108 row above: their coordinates are not indexed against the canonical entry. The "MBP 85-99" label is the numbering used in the DR2-MBP structural literature: Gauthier et al., *PNAS* 1998;95:11828-11833 (PMID 9751750), describe the immunodominant human MBP peptide in the DR2 complex as "residues 85-99"; the peptide chain of PDB 1BX2 is `ENPVVHFFKNIVTPR`, residues 217-231 of P02686. The manuscript now carries a gloss giving the canonical coordinates at the first mention of the label.

## Change log

- 2026-10-16: directory created (UniProt retrieval and coordinate verification).
- 2026-10-17: added the numbering crosswalk and inconsistencies 4 and 5 above, and the Gauthier et al. source for the MBP 85-99 label. Added `symbol_corrections.csv` (BFRF3 erratum). `SHA256SUMS.csv` and `../../MANIFEST.sha256` were updated for these README edits only; no other file in this directory changed.

## Scope

These sequences support coordinate auditing only. They do not establish binding,
natural presentation, TCR recognition, cross-reactivity, molecular mimicry, or
any MS mechanism.
