# Citation audit

Audit date: 2026-09-12. Scope: all externally sourced scientific, database, software, benchmark, and direct-assay claims in `main.tex`. Project-generated counts, rankings, model outputs, and classifications are identified as results of this study; they do not require external literature citations.

| Claim family | Manuscript location | Supporting reference(s) | Audit status |
|---|---|---|---|
| EBV–MS epidemiology and post-seroconversion neurofilament signal | Abstract; Introduction; Discussion | Bjornevik et al. 2022 | Cited |
| MS clinical description | Introduction | NINDS | Cited |
| EBNA1–GlialCAM antibody cross-reactivity | Abstract; Introduction; Discussion | Lanz et al. 2022 | Cited |
| HLA-II groove, core register, anchors, and TCR-facing positions | Introduction; Methods | Brown et al. 1993; Sewell 2012 | Cited |
| IEDB and UniProt data sources | Methods | Vita et al. 2025; UniProt Consortium 2025 | Cited |
| Sequence, binding, structural, and docking tools | Methods | Henikoff & Henikoff; Nilsson et al.; Racle et al.; Abramson et al.; Jurrus et al.; Basu & Wallner; Collins et al.; Yin et al.; Parizi et al.; Sali & Blundell | Cited |
| Hy.2E11 BALF5–MBP control | Methods | Lang et al. 2002 | Cited |
| PRIDE as an immunopeptidome evidence source | Results | Perez-Riverol et al. 2025 | Cited |
| Exact BALF5 627–641 / HLA-DRB1*15:01 direct-binding claim | Abstract; Results 3.5; Discussion | IEDB Assay 1774674; Hansen et al. 2007 | Cited; scope explicitly limited to the viral arm |
| EBNA1–ANO2 prior literature connection | Discussion | Thomas, Olsson et al. 2026 | Cited |
| HLA-DM peptide editing limitation | Discussion | Denzin & Cresswell 1995 | Cited |

## Direct-assay provenance note

- **Record:** IEDB MHC Ligand Assay 1774674.
- **Analyte / allele:** BALF5 627–641 (`TGGVYHFVKKHVHES`) in HLA-DRB1*15:01.
- **Outcome:** positive purified-MHC competitive radioactivity binding assay, qualitative.
- **Linked source:** Hansen BE, Rasmussen AH, Jakobsen BK, Ryder LP, Svejgaard A. *Tissue Antigens*. 2007;70:42–52. DOI: 10.1111/j.1399-0039.2007.00849.x.
- **Claim boundary:** this record supports binding of the viral BALF5 arm only. It does not demonstrate binding or presentation of TALDO1 216–230, shared TCR recognition, functional activation, or molecular mimicry.


## Addendum 2026-10-17

All 22 DOI-bearing entries in `references.tex` were re-checked against Crossref (title, journal, volume, pages, first author). Two bibliographic errors were found and fixed (PRIDE pages; Science Advances article number) and the MixMHC2pred v2 methods paper was added as `mix2023`. The `ninds` web reference could not be fetched (HTTP 403) and remains unverified. PDB `1H15` (HLA-DRA*01:01/DRB5*01:01 with an EBV DNA polymerase peptide, Lang et al. 2002) and `1BX2` (HLA-DRA*01:01/DRB1*15:01 with an MBP peptide; *J Exp Med* 1998;188:1511-1520) were confirmed against RCSB entry records. Later the same day `gauthier` (Gauthier et al., *PNAS* 1998) was added and verified, and the claim it supports (the "MBP 85-99" label) was checked against its abstract. See `CITATION_NOTES.md` for details.
