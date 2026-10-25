# V20: changes from V19

V20 (2026-10-17) is V19 plus corrections found in a post-V19 integrity audit. **No result, figure, table, dataset, or analysis changed.** `paper_v19_*` are untouched. Snapshots: `paper_v20_main.pdf/.txt`, `paper_v20_supplement.pdf/.txt`; `paper_v20_SOURCE_SHA256.txt` records the hashes of the sources, figures, and snapshots they were built from.

## Text changes

| Where | V19 | V20 | Evidence |
|---|---|---|---|
| `main.tex`, `supplement.tex` `\date` | September 21, 2026 | October 17, 2026 | |
| `supplement.tex` Tables S1, S3 (4 places) | viral arm "EBNA1 87-101" | "EBNA1 482-496" | Peptide `AEGLRALLARSHVER` is P03211 482-496 (live UniProt; frozen panel; stage-1 table). P03211 87-101 is a Gly-Ala repeat not in the panel. |
| `main.tex`, `supplement.tex`, first mention of "MBP 85-99" | bare label | gloss: "residue numbering as in the DR2--MBP structural literature [Gauthier et al.]; UniProt P02686 residues 217--231" | Gauthier et al. 1998 abstract gives "residues 85-99"; PDB 1BX2 peptide chain `ENPVVHFFKNIVTPR` = P02686 217-231. |
| `main.tex` Results, evidence review | no mention of IEDB's TALDO1 record | one sentence: IEDB lists a cellular MHC mass-spectrometry ligand record for TALDO1 108-122 (epitope 2008648), MHC recorded only as generic HLA class II, so it supports no specific allele | Live IEDB query; project census `output/revision_round_4/iedb_manual_browser_validation_2026-09-20/README.md` |
| `main.tex` limitations | "absent self-peptide records" | "absent allele-specific self-peptide records" | same |

## Reference changes (23 to 25 entries)

PRIDE pages D543-D551 to D543-D553; Science Advances article number `adj6367` to `eadj6367` (both checked against Crossref). Added `mix2023` (Racle et al., *Immunity* 2023, MixMHC2pred v2) and `gauthier` (Gauthier et al., *PNAS* 1998). Reference numbers after [7] shift accordingly.

## How V20 was checked against V19

Extracted text of each V20 PDF was diffed against the V19 snapshot (citation numbers ignored): similarity about 0.99 for both documents (0.989 main, 0.987 supplement), and every difference is one of the changes above or a line-break hyphenation artifact of text extraction.

## Build

Built offline with `tectonic` 0.17.0 (`--only-cached`; XeTeX engine, `xdvipdfmx` driver, the producer of the V19 PDFs), not the XeLaTeX named in `README.md`; fonts and layout were not compared page by page with V19. Main: 16 pages (V19: 16). Supplement: 11 pages (V19: 10): the 25-entry reference list now spills its last two entries onto an 11th page. No unresolved citations or references.

## Not part of V20

- The panel symbol `BFRF3` for P03187 (UniProt: TRX1/BORF1) is a data-label error, not a manuscript claim. It is recorded in `computational_package/data/reference_sequences_2026-10-16/symbol_corrections.csv` and left unedited because the panel hash is pinned in a protocol lock.
- The old "EBNA1 87-101" label persists in `main_google_docs_aligned.tex` and the `paper_v11`-`paper_v19` snapshots, which are provenance copies. `main_bmc_bioinformatics.tex` carries the same correction.
- Target journal (BMC Bioinformatics vs Scientific Reports) is undecided; submission documents were not changed.
