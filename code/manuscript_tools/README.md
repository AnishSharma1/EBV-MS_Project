# Manuscript .docx tools

Builds editable Word copies of the manuscript and supplement from the LaTeX
sources, for collaborators who comment in Word rather than in LaTeX.

```sh
python3 code/manuscript_tools/build_docx.py
```

This writes `manuscript/paper_v<N>_main.docx` and `..._supplement.docx`, where
`<N>` is the latest `paper_v<N>_main.pdf` snapshot. It reads the LaTeX sources
and never modifies them, so `test_repository_integrity.py` still passes.

Requires `pandoc` and `pdftoppm` (Poppler). The .docx files are a derived
convenience copy: **`main.tex` and `supplement.tex` remain the source of
truth**, and edits made in Word have to be carried back into the LaTeX by
hand.

## Why a preprocessing step is needed

Pandoc cannot convert these sources directly. It exits 0 with no warning while
losing content four separate ways, so a direct conversion looks fine and is
not. `prepare_for_docx.py` rewrites the LaTeX to avoid each one:

| Failure | Effect on the output |
|---|---|
| `\cite` inside a `thebibliography` document | **Every citation vanishes.** "optic nerves [2]" becomes "optic nerves ." |
| `thebibliography` | Collapses to one run-on paragraph, leaking its `{99}` widest-label argument |
| `\resizebox{..}{..}{<tabular>}` | **The whole table is discarded** — this silently removed Tables S4, S5 and S6 |
| `\ref` | Numbers only `\label`ed floats, so `\ref{fig:ranksens}` rendered "Figure 1" where the PDF says "Figure 2" |

Two more differences are handled because Word has no equivalent of the LaTeX
mechanism:

- Figures are PDFs, which Word cannot display, so they are rasterised to PNG.
- `\renewcommand{\thefigure}{S\arabic{figure}}` is ignored, so "Figure S1." and
  "Table S1." prefixes are injected into captions in source order. `\caption*`
  is unnumbered in LaTeX and is skipped, matching the compiled PDF.

`shrink_wide_tables.py` then does what `\scriptsize` + `\resizebox` did: without
it a nine-column table gets nine equal columns and wraps mid-word.

## Verifying a build

`compare_to_pdf.py` diffs the generated .docx against `paper_v<N>_*.txt`, which
is `pdftotext` of the tracked PDF snapshot and therefore the reference for what
the version actually says:

```sh
python3 code/manuscript_tools/compare_to_pdf.py \
    manuscript/paper_v20_main.docx manuscript/paper_v20_main.txt
```

It reports prose words present in the PDF but missing from the .docx. Expect a
small residue of differences that are not content: the PDF text carries page
numbers and LaTeX section numbers, and hyphenates across line breaks.

When checking table values by hand, read the `w:tc` cells rather than searching
flattened document text — cell text is split across runs, and the minus sign in
the tables is an ASCII `-`, not U+2212.

## Known limitations

- Headings carry Word heading styles but not LaTeX's "2.1" numbering. This is
  what lets a reviewer restructure sections cleanly; apply a multilevel list in
  Word if the numbers are wanted.
- Equations convert to native OMML and are editable, but round-tripping them
  back to LaTeX is not automated.
- The `.docx` files are not part of the frozen package and are not recorded in
  `paper_v<N>_SOURCE_SHA256.txt`; rebuild them after any source change.

Unit tests: `code/tests/test_manuscript_docx_tools.py`.
