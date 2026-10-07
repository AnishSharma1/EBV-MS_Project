# V18 target-journal submission guide

Prepared 21 September 2026 for the V18 manuscript, *Register-aware computational prioritization of EBV–self peptide pairs across HLA class II contexts*. This is a pre-submission checklist, not a claim that a journal will send the work for review. Publisher instructions and fee schedules change; recheck the linked official page immediately before uploading.

## Current V18 baseline

The current main manuscript is a 16-page LaTeX PDF with a structured abstract, four displayed figures, a supplementary manuscript, public GitHub link, and a **prepared but not yet published** DOI archive. The work is a deliberately cautious computational prioritization study: it does not claim experimental peptide presentation, stable register occupancy, T-cell cross-reactivity, molecular mimicry, or MS causality.

That framing is scientifically appropriate. It also matters for journal selection: journals that expect independent validation of public-data computational analyses are a poor match until an external cohort or experimental validation is added.

## One cross-journal package to finish first

1. **Publish the frozen release before submission.** Deposit the exact source code, inputs where redistributable, derived tables, environments/locks, manifests, failure artifacts, and instructions in a versioned archive. Replace “a public DOI archive has been prepared” with the minted DOI, version, and access date. A prepared ZIP is not a citable DOI archive.
2. **Make the reproducibility path executable.** Retain the GitHub repository, but give the archival release a `README`, a top-level reproduction command, expected outputs/checksums, software versions, licenses for third-party inputs, and a clear statement of anything that cannot be redistributed.
3. **Prepare transparent declarations.** Confirm author affiliations, corresponding-author details, funding, competing interests, author contributions, acknowledgments, ethics status, and data/code availability. Do not add claims that no data exist when derived data and code are available.
4. **Prepare an AI-use disclosure where required.** Bioinformatics Advances explicitly asks for disclosure of AI use in content, code, or data; state accurately what was used and how author(s) verified the final work. Do not list an AI system as an author.
5. **Keep the conclusion discipline.** Across all targets, retain the distinction between computational routing and biological validation. The two highest-impact targets will likely require stronger independent validation, not merely a reformat.

## At-a-glance fit and work needed

| Journal | Plausible article route | V18 fit | Work needed before upload |
|---|---|---|---|
| Scientific Reports | Article | Good potential fit; concise-rewrite needed | Use an unstructured 200-word abstract, 20-word title, Nature reference style, data statement, AI-use documentation, and fit the recommended 4,500-word/11-typeset-page presentation. |
| BMC Bioinformatics | Research Article | Good methodological fit | Reformat for BMC, add line/page numbers and double spacing, use Vancouver references, and foreground reusable computational contribution and reproducibility. |
| Bioinformatics Advances | Original Article | Good scope, but severe length mismatch | Reduce to no more than 8 pages/6 figures/6 tables/50 references; unstructured 200-word abstract; author–date references; mandatory AI disclosure. |
| Viruses | Article | Reasonable viral-immunology fit | Lead with EBV biology and experimental testability; use the MDPI front/back matter and numbered reference style; publish archive and complete declarations. |
| HLA | Original Article | Topically close but evidence-sensitive | Center HLA-II biology, use correct WHO HLA nomenclature, provide data statement, and make the computational-only status explicit. A concise pre-submission synopsis is available. |
| Frontiers in Bioinformatics | Original Research | Format fit, but current evidence-bar mismatch | Add appropriate independent validation before submission; make materials/data FAIR and repository-addressable; prepare contribution-to-field statement. |

## 1. Scientific Reports (Nature Portfolio)

**Official guidance:** <https://www.nature.com/srep/author-instructions/submission-guidelines>

### What to prepare

- Submit as an **Article**. Scientific Reports does not generally impose strict word/page limits, but recommends no more than **11 typeset pages** and **4,500 main-text words** (excluding Abstract, Methods, references, and figure legends).
- Use a title of no more than **20 words** that states the main message in one scientifically accurate sentence. Supply an **unstructured abstract of no more than 200 words**, without citations or figure references, and up to **six keywords/key phrases**. Graphical abstracts are not supported.
- A common main-text organization is Introduction, Results, Discussion, and Methods. Follow with Nature-style numbered references, Acknowledgements, Author Contributions, the mandatory Data Availability Statement, Additional Information/Competing Interests, figure legends, and editable tables. Footnotes are not allowed.
- Keep display items to **eight figures and/or tables** combined. Figure legends must be no more than 350 words each. Submit high-resolution individual figure files and combine the supplementary material as a single separate file where possible.
- Include a cover letter with the corresponding author’s contact information, why the work suits *Scientific Reports*, suitable/excluded reviewers, and any prior editorial discussion. Initial submissions may be Word, LaTeX, or PDF; for LaTeX, use the Springer Nature template and ensure the complete source compiles without warnings. Revised manuscripts must provide editable Word or LaTeX text rather than PDF.
- Add a transparent Methods disclosure of any LLM use. The journal states that LLMs cannot be authors and that their use must be documented in Methods or a suitable alternative location.

### V18-specific assessment

- This is a materially more realistic target than *Science Advances*: the editorial question is whether the computational study is methodologically sound, clearly reported, and supported by accessible underlying materials—not whether it demonstrates a field-defining advance. The current careful limitations language helps.
- V18’s four main figures fit the eight-display cap. The 16-page PDF should be condensed toward the journal’s 11-typeset-page/4,500-main-text-word recommendation, while allowing Methods to remain detailed. Move extended sensitivity grids, audit tables, and supporting implementation detail into the existing supplement and archival release.
- Shorten the present title to 20 words or fewer, convert the structured abstract to one unstructured paragraph, and convert the reference list to Nature’s numbered style. Confirm the final reference count is at or below the journal’s 60-reference guideline.
- Put the archived code/data DOI in the mandatory Data Availability Statement before References. The current declaration that an archive has merely been prepared is not sufficient for a submission that promises accessible materials.

## 2. BMC Bioinformatics (Springer Nature)

**Official guidance:** <https://link.springer.com/journal/12859/submission-guidelines>  
**Research Article page:** <https://link.springer.com/journal/12859/submission-guidelines/research-article>

### What to prepare

- BMC Bioinformatics welcomes novel computational algorithms, software, models, and tools for biological data, including sequence and structural analysis. Its assessment emphasizes a defined and sound question, appropriate methods/analysis, and relevant reporting standards rather than perceived impact alone.
- LaTeX is accepted. Prepare a double-spaced manuscript with line and page numbering, SI units, and no manual page breaks.
- Use the journal’s Research Article structure and **Vancouver** numbered reference style. Put URLs in the numbered reference list rather than only in the narrative text.
- Confirm the current open-access fee/licence route at submission; the displayed fee can change.

### V18-specific assessment

- This is the cleanest direct format match. Recast the title, abstract, and opening paragraphs around the **register-aware reproducible prioritization workflow**, not a claim that a disease mechanism has been found.
- Add line numbers and page numbers to the submission build; do not alter the scientific content merely to satisfy a house style.
- Preserve the limitations and add an explicit code/data release DOI once published. The existing GitHub link is useful but should not be the only reproducibility endpoint.
- A BMC version may retain substantially more methods detail in the main text than the 8-page Bioinformatics Advances route.

## 3. Bioinformatics Advances (Oxford University Press)

**Official guidance:** <https://academic.oup.com/bioinformaticsadvances/pages/author-guidelines>

### What to prepare

- The manuscript belongs most naturally in **Immunoinformatics** or **Structural Bioinformatics**, subject to the journal’s current section choices.
- For an Original Article, the official guide specifies: **up to 8 pages**, **200-word abstract**, **up to 5 keywords**, **up to 50 references**, and **no more than 6 figures and 6 tables**. Required core sections are Abstract, Introduction, Methods, Discussion, and References.
- Use the OUP Modern Large template and **author–date** citations. Upload a PDF for initial review; supply source files after acceptance as requested.
- Include full author/affiliation information, corresponding-author mailing and email address, a running head of no more than 50 characters, figure alt text, and ORCID(s) required by the journal.
- Disclose any AI use affecting content, code, data, images, or similar research outputs in the cover letter and in Methods or Acknowledgments as applicable; document LLM use in supplementary material when required by the policy.

### V18-specific assessment

- The V18 main PDF is 16 pages, so this requires a real **short-form rewrite**, not a cosmetic format conversion. Move audit detail, calibration, complete threshold grids, structural sensitivity detail, and extended limitations to the supplement while leaving the central method, main result, validation limits, and one or two decisive figures in the article.
- The existing four figures fit; check that all tables plus figures remain within the separate six-and-six caps. Reduce references to 50 or fewer and convert the current numbered style.
- Preserve the critical caveat in the shortened Discussion: the output is candidate prioritization, not evidence of T-cell cross-reactivity.

## 4. Viruses (MDPI)

**Official guidance:** <https://www.mdpi.com/journal/viruses/instructions>  
**Publisher layout guide:** <https://www.mdpi.com/authors/layout>

### What to prepare

- Use the journal’s current article template and its front matter, main-text, and back-matter order. For a research paper, prepare title/author/affiliation information, abstract, keywords, main text, supplementary-material declaration, author contributions, funding, institutional-review-board and consent statements where applicable, data-availability statement, acknowledgments, conflicts statement, and references.
- The abstract should stand alone: avoid citations, table/figure references, and unsupported abbreviations. Verify the current article-type word/figure limits in the live portal because the publisher page is periodically updated.
- Use MDPI’s numbered in-text citation/reference convention and ensure all software, databases, and deposited materials are citable.
- Make the data/methods sufficiently detailed for replication and be ready to provide underlying data to editors/referees if requested.

### V18-specific assessment

- The biological entry point is strong: EBV, HLA-II, and MS relevance should lead the title/abstract while the claims remain computational and hypothesis-generating.
- Bring the current Declarations into the MDPI back-matter order, and expand the data statement after the archive DOI exists. Use CRediT-compatible wording for the single-author contribution statement.
- Do not present the unvalidated peptide pairs as viral immune escape, natural presentation, or a causal MS mechanism; frame them as experimentally testable candidates.

## 5. HLA (Wiley)

**Official guidance:** <https://onlinelibrary.wiley.com/page/journal/20592310/homepage/forauthors.html>

### What to prepare

- Choose **Original Article** unless the paper is radically condensed. The journal also explicitly welcomes Brief Communications on MHC peptide motifs, but that format is not a substitute for a full computational-methods paper.
- The guidance specifies a 250-word abstract for Original Articles. It requires ORCID and a data-availability statement for accepted articles; published-result data should be archived in an appropriate public repository.
- Use official WHO HLA nomenclature consistently. If reporting novel nucleotide or amino-acid sequences, deposit them in GenBank/ENA before submission and provide accession numbers. V18 analyzes existing peptide records, so this sequence-deposition condition should not apply unless new sequences are added.
- Submit figures, supplemental files, and appendices as separate files as directed by the Research Exchange workflow. A one-page pre-submission synopsis may be sent to the Editor-in-Chief to assess suitability.

### V18-specific assessment

- The HLA-II register focus makes this thematically attractive, particularly if the paper emphasizes peptide motif/register handling and HLA-restricted interpretation.
- It is not yet a secure fit: the key biological claims remain unvalidated and one lead has limited direct evidence only for the viral arm. Submit only with a plainly computational title/abstract and a strong reproducibility release; consider the optional synopsis before investing in a journal-specific rebuild.
- Audit every allele name and presentation statement against WHO nomenclature and retain the distinction between HLA-DRB1*15:01 and the DR15 haplotype-level benchmark.

## 6. Frontiers in Bioinformatics

**Official journal guidance:** <https://www.frontiersin.org/journals/bioinformatics/for-authors/author-guidelines>  
**Publisher policies:** <https://www.frontiersin.org/guidelines/policies-and-publication-ethics>

### What to prepare

- Submit as **Original Research** in the best-fitting section, likely Protein Bioinformatics or Integrative Bioinformatics, after checking the active section scope.
- Use a single-paragraph abstract with no citations, figures, or tables. Include the manuscript word count plus number of figures/tables on the first page; follow the article-type word limit in the current portal.
- Frontiers accepts Word or LaTeX. For LaTeX, upload the `.tex`, PDF, bibliography, and all required source files. Supply concise, meaningful alt text for figures.
- Prepare a 200-word plain-language “contribution to the field” statement for the submission form; Frontiers does not use a conventional cover letter.
- A data-availability statement is required for every article. Materials necessary to reproduce the work should be identified in Methods and preferably placed in a persistent online repository. Cite datasets in the reference list where appropriate.
- The publisher states that manuscripts based solely on bioinformatics analyses or public data must include **appropriate validation** (e.g., an independent clinical/patient cohort, biological validation, or specified multi-database validation). It may regard a manuscript without it as unsuitable.

### V18-specific assessment

- **This is the decisive barrier.** V18 is a computational analysis of public data and reports that its proposed pairs require later binding, register mapping, and T-cell testing. Its internal calibration/sensitivity analyses improve transparency but are not clearly one of the policy’s accepted independent-validation examples.
- The right next move is experimental peptide–HLA binding/register validation or a genuinely independent validation dataset, then a frozen archive DOI. Reformatting alone should wait until that evidence is available.
- If validation is added, retain the current cautious language and make the repository path FAIR, fully versioned, and visible in both Methods and the data-availability statement.

## Recommended submission order today

1. **Scientific Reports** and **BMC Bioinformatics** are the strongest present-tense fits for a carefully bounded, reproducible computational workflow.
2. **Bioinformatics Advances** is a strong scope match if a concise 8-page article can retain the central method and its limitations.
3. **Viruses** is viable when the viral-immunology motivation is foregrounded and all data/declarations are finalized.
4. **HLA** is conceptually close but should follow a positive response to a one-page synopsis or additional HLA-focused evidence.
5. **Frontiers in Bioinformatics** should wait for independent validation; its stated hurdle is scientific evidence rather than copyediting.

## Final pre-upload checklist

- [ ] Minted archival DOI/version replaces “prepared DOI archive” everywhere.
- [ ] GitHub URL, archive DOI, licenses, file manifest, and data/code availability statement agree exactly.
- [ ] All main-text values, supplement values, and source tables have been checked against frozen outputs.
- [ ] Target-specific abstract, title page, declarations, reference style, figure/table cap, and supplement packaging are applied.
- [ ] Author/affiliation/correspondence/ORCID information is complete and approved.
- [ ] Cover letter, reviewer suggestions, contribution-to-field statement, or HLA synopsis are prepared only for the intended target.
- [ ] Any AI-use disclosure is accurate and complies with the selected journal’s live policy.
- [ ] The target journal’s live portal is rechecked for article-type limits, fees, template revision, and required forms on the day of submission.
