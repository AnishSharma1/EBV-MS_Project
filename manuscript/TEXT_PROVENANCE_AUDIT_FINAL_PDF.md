# Final-PDF provenance comparison

## Source

- Authoritative PDF extraction: `FINAL_GOOGLE_DOCS_PDF_EXTRACTED.txt`
- Current manuscript: `main.tex`
- Current narrative tokens analysed: 3299

## Consecutive phrase overlap

| Minimum matching phrase length | Current-text tokens in a matching phrase | Coverage |
|---:|---:|---:|
| 5 words | 1120 | 33.9% |
| 8 words | 663 | 20.1% |
| 12 words | 459 | 13.9% |

## Interpretation

This is a literal phrase-provenance metric, not an AI-detector score or a determination of who wrote any passage. The 8-word line is the most useful conservative estimate: it counts current manuscript words that participate in at least one consecutive eight-word phrase found in the final Google Docs PDF. It does not count paraphrases, user edits made after the PDF, headings, tables, figure captions, or citations as matching prose.
