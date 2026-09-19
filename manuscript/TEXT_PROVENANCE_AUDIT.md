# Text provenance audit

## Scope

This audit compares the current LaTeX prose with the supplied Google Docs text. It is a provenance check, not an AI-detector score. “Not directly traceable” means the paragraph does not closely match the supplied snapshot; it does **not** establish who wrote it, because the snapshot may omit later user revisions.

- Canonical source: `/Users/anishsharma/.codex/attachments/e1b5ba5c-a523-43df-b035-2493ea7773d9/pasted-text.txt`
- Source SHA-256: `b419e8d32cf9c3f3a507eafcf9b2e069f1d18b6e271d889a798300eada4e9219`
- Current source: `main.tex`
- Current source SHA-256: `70d21a8fc62c971ca312e2f5bc8d30119c5025cf78b6fd6fb062cef9012a9990`
- User-source paragraphs analysed: 66
- Current prose paragraphs analysed: 39

## Summary

| Classification | Paragraphs |
|---|---:|
| Direct match | 0 |
| Partial match | 8 |
| Not directly traceable | 31 |

## Paragraph-level results

| Current paragraph | Classification | Best source paragraph | Similarity | Opening text |
|---:|---|---:|---:|---|
| 1 | Not directly traceable | 22 | 0.32 | center 1Liberal Arts and Science Academy, 1012 Arthur Stiles Road, Austin, TX 78721, USA\\ *Correspondence: anishsharma0212@gmail.com center |
| 2 | Not directly traceable | 1 | 0.30 | abstract Background Over the years, researchers have noticed that a prior Epstein--Barr virus (EBV) infection is present in almost all multiple sclerosis pat... |
| 3 | Partial match | 1 | 0.85 | Results In this study, a multi-stage computational pipeline incorporating NetMHCIIpan, MixMHC2pred, TCR-facing BLOSUM62 similarity scoring, and 3D structural... |
| 4 | Not directly traceable | 52 | 0.33 | Conclusions The study results provide a prioritized shortlist for experimental binding, register-mapping, and T-cell cross-reactivity testing. abstract |
| 5 | Not directly traceable | 3 | 0.60 | Multiple sclerosis (MS) is an autoimmune disease that can cause degeneration of motor, sensory, visual, and cognitive abilities through immune-mediated damag... |
| 6 | Partial match | 4 | 0.85 | Previous studies such as Lanz et al. identified cross-reactive antibodies ; comprehensive mapping of T-cell cross-reactivity remains incompletely defined. T-... |
| 7 | Partial match | 5 | 0.75 | The computational pipeline developed here integrates peptide--HLA binding predictions, register-aware sequence alignments, and predicted peptide--HLA structu... |
| 8 | Partial match | 6 | 0.92 | This study describes a computational prioritization pipeline designed to identify high-interest EBV--human self-peptide candidate pairs for subsequent experi... |
| 9 | Partial match | 7 | 0.76 | An EBV repertoire was curated from the Immune Epitope Database (IEDB) by filtering EBV source peptides with positive T-cell assay outcomes from MHC class II ... |
| 10 | Partial match | 8 | 0.74 | The initial sequence screen was limited to HLA-DRB1*03:01, HLA-DRB1*08:01, HLA-DRB1*13:03, and HLA-DRB1*15:01. Each HLA allele was evaluated independently; p... |
| 11 | Partial match | 9 | 0.73 | Sequence similarity was calculated with BLOSUM62 . Because all non-peptide components are identical within an allele-bound comparison, scoring was restricted... |
| 12 | Partial match | 11 | 0.77 | Eight candidate pairs selected from the 6,400-comparison matrix underwent a standardized manual evidence audit. Predicted binding support, predicted register... |
| 13 | Not directly traceable | 12 | 0.69 | For the initial eight-pair review, peptide--HLA class II percentile ranks were predicted with NetMHCIIpan-4.3 EL and MixMHC2pred context predictions . Full b... |
| 14 | Not directly traceable | 13 | 0.58 | Binding-register analysis defined the register as the predicted P1--P9 core within the open-ended HLA class II groove . NetMHCIIpan-4.3 EL and MixMHC2pred co... |
| 15 | Not directly traceable | 10 | 0.68 | Heterotrimeric pMHC complexes (HLA-DRA, relevant HLA-DR beta chain, and target peptide) were modeled using the AlphaFold 3 web server, with five returned mod... |
| 16 | Not directly traceable | 14 | 0.65 | The literature-derived Hy.2E11 cross-reactive benchmark, BALF5 627--641--MBP 85--99, served as a structural-comparison reference . The benchmark contains HLA... |
| 17 | Not directly traceable | 16 | 0.54 | Poisson--Boltzmann/APBS pMHC surface-potential comparison was explored as an exploratory surface-resemblance descriptor . Its score was excluded from primary... |
| 18 | Not directly traceable | 17 | 0.63 | Separately, TCRModel2 generated five-component TCR--pMHC models comprising HLA-alpha, HLA-beta, peptide, TCR-alpha, and TCR-beta chains . The calibration ref... |
| 19 | Not directly traceable | 18 | 0.65 | Reproducibility materials in the project GitHub repository retain software environments, random seeds, input logs, and SHA-256 dataset checksums. Ranking was... |
| 20 | Not directly traceable | 20 | 0.65 | A total of 6,400 allele-specific peptide-pair comparisons were evaluated across four target alleles, with 1,600 pairs per allele. Sequence-similarity metrics... |
| 21 | Not directly traceable | 29 | 0.36 | BALF5 627--641--TALDO1 108--122 scored 0.281 in HLA-DRB1*13:03 (rank 13/1,600), while BALF5 627--641--TALDO1 216--230 scored 0.314 in HLA-DRB1*15:01 (rank 14... |
| 22 | Not directly traceable | 32 | 0.41 | Eight candidates selected from the 6,400-pair matrix underwent deep evidence review, yielding two medium-priority candidates, six holds, and zero high-priori... |
| 23 | Not directly traceable | 33 | 0.35 | Both BALF5--TALDO1 pairs showed predicted core-register agreement, although the registers remain unconfirmed experimentally. Neither pair satisfied full bind... |
| 24 | Not directly traceable | 38 | 0.44 | The two medium-priority pairs were structurally assessed alongside a DRB5*01:01 restriction test. Thirteen AlphaFold 3 jobs, each returning five models, gene... |
| 25 | Not directly traceable | 41 | 0.38 | Models of shorter 11-mers frequently converged to consistent placements, but these placements were rarely retained when parent 15-mers were modeled, producin... |
| 26 | Not directly traceable | 42 | 0.62 | PANDORA coupled with MODELLER was used for cross-method structural reproducibility . The reciprocal PANDORA controls 1BX2 and 6CQQ were separate from the Hy.... |
| 27 | Not directly traceable | 43 | 0.62 | The subsequent ensemble contained 1,200 models spanning BALF5, TALDO1, and eight P6/P7 charge-substitution variants across two structural templates and three... |
| 28 | Not directly traceable | 46 | 0.45 | A separate expanded cohort of 49 pairs was recovered from the 6,400-pair ranking. NetMHCIIpan, MixMHC2pred, and IEDB-recommended binding prediction profiles ... |
| 29 | Not directly traceable | 48 | 0.63 | For public-evidence assessment, 98 automated IEDB API queries returned zero retrievable records because the endpoint was unavailable. The resulting frozen st... |
| 30 | Not directly traceable | 49 | 0.55 | A provisional shortlist contained seven candidates: one strong-medium, two weak-medium, and four strong-low. The remaining 42 candidates were untiered becaus... |
| 31 | Not directly traceable | 59 | 0.66 | This multi-stage computational screening and evidence-review workflow prioritized EBV--human peptide pairs for later experimental validation. Its central fin... |
| 32 | Not directly traceable | 61 | 0.54 | BALF5--TALDO1 did not have the highest raw sequence-similarity score; EBNA1--ANO2 and EBNA3B--MOG scored higher. Raw sequence alignment was only the first ev... |
| 33 | Not directly traceable | 62 | 0.38 | Agreement among sequence-based register predictors did not guarantee stable structural coordinates. Shorter 11-mer inputs often converged in the models analy... |
| 34 | Not directly traceable | 63 | 0.55 | Bjornevik et al. established the epidemiological link between EBV infection and MS, and Lanz et al. reported antibody cross-reactivity between EBNA1 and Glia... |
| 35 | Not directly traceable | 64 | 0.55 | Electrostatic surface-potential comparison was excluded from primary ranking because its control gate did not discriminate the calibration control from backg... |
| 36 | Not directly traceable | 65 | 0.56 | The search space was limited to pre-curated viral epitopes, selected myelin/CNS proteins, and four HLA-DRB1 alleles. HLA-DQ and HLA-DP loci were omitted. Ran... |
| 37 | Not directly traceable | 66 | 0.69 | Peptide--HLA binding can be quantified for viral and self peptides against recombinant HLA-DRB1 molecules by fluorescence polarization or surface plasmon res... |
| 38 | Not directly traceable | 60 | 0.63 | The first experimental priority is HLA-DRB1*15:01 BALF5 627--641--TALDO1 216--230, because of its strong-medium classification, its passing of expanded bindi... |
| 39 | Not directly traceable | 51 | 0.28 | AF3: AlphaFold 3; APBS: Adaptive Poisson--Boltzmann Solver; EBV: Epstein--Barr virus; HLA: human leukocyte antigen; IEDB: Immune Epitope Database; MS: multip... |

## Action rule

To create a text-only version traceable to the supplied Google Docs draft, retain only “Direct match” paragraphs, or restore the corresponding source paragraphs for every other row. Formatting elements, section headings, references, tables, figures, and declarations should be reviewed separately because they are not prose paragraphs in the supplied source.
