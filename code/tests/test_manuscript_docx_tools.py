#!/usr/bin/env python3
"""Tests for the LaTeX-to-.docx conversion helpers.

Run from any directory:  python3 code/tests/test_manuscript_docx_tools.py

These cover the rewriting rules in ``code/manuscript_tools/`` that work around
pandoc's silent content loss. Each rule is tested on a small fixture and, where
it matters, against the real manuscript sources. No .docx is built here and no
scientific result is reproduced.
"""

import re
import sys
from pathlib import Path


def _find_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in (here, *here.parents):
        if (candidate / "manuscript").is_dir() and (candidate / "code").is_dir():
            return candidate
    raise SystemExit("could not locate the repository root")


ROOT = _find_root()
MANUSCRIPT = ROOT / "manuscript"
TOOLS = ROOT / "code" / "manuscript_tools"
sys.path.insert(0, str(TOOLS))

import prepare_for_docx as prep           # noqa: E402
import shrink_wide_tables as shrink       # noqa: E402

REFERENCES = r"""
\begin{thebibliography}{99}
\bibitem{alpha} Alpha A. First work. \textit{Journal}. 2020;1:1-2.
\bibitem{beta} Beta B. Second work. \textit{Journal}. 2021;2:3-4.
\bibitem{gamma} Gamma C. Third work. \textit{Journal}. 2022;3:5-6.
\end{thebibliography}
"""


# --- citation numbering ------------------------------------------------

def test_citation_numbers_follow_bibitem_order():
    numbers = prep.citation_numbers(REFERENCES)
    assert numbers == {"alpha": 1, "beta": 2, "gamma": 3}, numbers


def test_single_and_grouped_cites_resolve():
    text, missing = prep.resolve_citations(
        r"one \cite{alpha} and two \cite{beta,gamma}.",
        prep.citation_numbers(REFERENCES))
    assert not missing, missing
    assert text == "one [1] and two [2, 3].", text


def test_unknown_cite_key_is_reported():
    _, missing = prep.resolve_citations(
        r"\cite{nope}", prep.citation_numbers(REFERENCES))
    assert missing == ["nope"], missing


def test_real_cite_keys_match_the_compiled_pdf_numbering():
    """The .docx must use the same numbers the tracked PDF shows."""
    numbers = prep.citation_numbers(_read(MANUSCRIPT / "references.tex"))
    # Spot values read from paper_v20_main.txt: "optic nerves [2]" is \cite{ninds}
    # and "may contact the TCR [23, 24]" is \cite{brown,sewell}.
    assert numbers["ninds"] == 2, numbers.get("ninds")
    assert numbers["brown"] == 23, numbers.get("brown")
    assert numbers["sewell"] == 24, numbers.get("sewell")


def test_every_cite_in_the_sources_resolves():
    numbers = prep.citation_numbers(_read(MANUSCRIPT / "references.tex"))
    for name in ("main.tex", "supplement.tex"):
        _, missing = prep.resolve_citations(_read(MANUSCRIPT / name), numbers)
        assert not missing, f"{name}: unresolved {sorted(set(missing))}"


# --- bibliography ------------------------------------------------------

def test_bibliography_becomes_one_paragraph_per_reference():
    block = prep.bibliography_paragraphs(REFERENCES)
    assert r"\section*{References}" in block
    assert "[1] Alpha A." in block
    assert "[3] Gamma C." in block
    # the {99} widest-label argument must not survive into the text
    assert "99" not in block.replace("thebibliography", ""), block


def test_bibliography_keeps_every_real_reference():
    text = _read(MANUSCRIPT / "references.tex")
    block = prep.bibliography_paragraphs(text)
    expected = len(re.findall(r"\\bibitem\{", text))
    found = len(re.findall(r"^\[\d+\] ", block, re.M))
    assert found == expected, f"{found} paragraphs for {expected} bibitems"
    assert f"[{expected}] " in block


# --- box macros that make pandoc drop tables ---------------------------

def test_resizebox_unwraps_to_its_payload():
    """The tabular must survive; pandoc discards the table if the macro does not.

    The payload keeps the ``%`` that ended the macro's opening line, which is
    just a LaTeX comment and is ignored downstream.
    """
    text, count = prep.unwrap_macro(
        r"\resizebox{\textwidth}{!}{%" "\n" r"\begin{tabular}{ll}a&b\\\end{tabular}%"
        "\n}", "resizebox", 3, keep=2)
    assert count == 1, count
    assert r"\resizebox" not in text, text
    assert r"\begin{tabular}{ll}a&b" in text, text
    assert r"\end{tabular}" in text, text
    before = text.split(r"\begin{tabular}")[0].strip()
    assert before in ("", "%"), f"unexpected leading content: {before!r}"


def test_shortstack_optional_argument_is_skipped():
    """\\shortstack[l]{..} has a positional arg before its braced one."""
    text, count = prep.unwrap_macro(
        r"x & \shortstack[l]{Candidate\\pair} & y", "shortstack", 1, keep=0,
        transform=lambda s: " ".join(s.replace(r"\\", " ").split()))
    assert count == 1
    assert text == "x & Candidate pair & y", text


def test_shortstack_linebreak_would_not_end_the_row():
    text, _ = prep.unwrap_macro(
        r"\shortstack{Monte Carlo\\$p$}", "shortstack", 1, keep=0,
        transform=lambda s: " ".join(s.replace(r"\\", " ").split()))
    assert r"\\" not in text, text


def test_balanced_arg_handles_nesting_and_escapes():
    content, end = prep.balanced_arg(r"{a{b}c\}d}tail", 0)
    assert content == r"a{b}c\}d", content
    assert end == len(r"{a{b}c\}d}")


def test_unbalanced_braces_raise_rather_than_corrupt():
    try:
        prep.balanced_arg("{no end", 0)
    except ValueError:
        return
    raise AssertionError("unbalanced braces should raise")


# --- float numbering and cross references ------------------------------

SAMPLE = r"""
\begin{figure}\includegraphics{a.pdf}\caption{First figure.}\label{fig:one}\end{figure}
\begin{table}\caption*{Unnumbered summary.}\begin{tabular}{l}x\\\end{tabular}\end{table}
\begin{table}\caption{First table.}\label{tab:one}\end{table}
\begin{figure}\includegraphics{b.pdf}\caption{Second figure.}\label{fig:two}\end{figure}
See Figure~\ref{fig:two} and Table~\ref{tab:one}.
"""


def test_captions_are_numbered_in_source_order():
    text, figures, tables, _ = prep.number_captions(SAMPLE, "")
    assert (figures, tables) == (2, 1), (figures, tables)
    assert "Figure 1. First figure." in text
    assert "Figure 2. Second figure." in text
    assert "Table 1. First table." in text


def test_unnumbered_caption_is_left_alone():
    text, _, tables, _ = prep.number_captions(SAMPLE, "")
    assert tables == 1, "caption* must not advance the table counter"
    assert r"\caption*{Unnumbered summary.}" in text


def test_supplement_prefix_renumbers_both_counters():
    text, _, _, labels = prep.number_captions(SAMPLE, "S")
    assert "Figure S1. First figure." in text
    assert "Table S1. First table." in text
    assert labels["fig:two"] == "S2", labels


def test_refs_resolve_to_the_float_number():
    text, _, _, labels = prep.number_captions(SAMPLE, "")
    assert labels == {"fig:one": "1", "tab:one": "1", "fig:two": "2"}, labels
    text, unknown = prep.resolve_refs(text, labels)
    assert not unknown, unknown
    assert "See Figure 2 and Table 1." in text


def test_unknown_ref_is_reported():
    _, unknown = prep.resolve_refs(r"see \ref{nope}", {})
    assert unknown == ["nope"], unknown


def test_real_cross_references_match_the_compiled_pdf():
    """paper_v20_main.txt says "(Figure 2)"; pandoc alone produced "Figure 1"."""
    text = _read(MANUSCRIPT / "main.tex")
    text, _, _, labels = prep.number_captions(text, "")
    resolved, unknown = prep.resolve_refs(text, labels)
    assert not unknown, unknown
    assert labels["fig:ranksens"] == "2", labels
    assert labels["fig:length"] == "3", labels
    assert "all 6,400 pairs (Figure 2)" in resolved


# --- wide table shrinking ----------------------------------------------

TABLE_XML = (
    "<w:tbl><w:tblPr><w:tblW w:type=\"auto\" w:w=\"0\" /></w:tblPr>"
    "<w:tblGrid><w:gridCol w:w=\"900\"/><w:gridCol w:w=\"900\"/></w:tblGrid>"
    "<w:tr><w:tc><w:p><w:r><w:rPr><w:b/><w:sz w:val=\"24\"/></w:rPr>"
    "<w:t>Header</w:t></w:r></w:p></w:tc>"
    "<w:tc><w:p><w:r><w:t>Plain</w:t></w:r></w:p></w:tc></w:tr></w:tbl>"
)


def test_shrink_sets_one_size_on_every_run():
    out = shrink.shrink_table(TABLE_XML, 16)
    assert out.count('<w:sz w:val="16"/>') == 2, out
    assert 'w:val="24"' not in out, "the old size must be replaced"


def test_shrink_preserves_text_and_other_formatting():
    out = shrink.shrink_table(TABLE_XML, 16)
    assert "<w:t>Header</w:t>" in out and "<w:t>Plain</w:t>" in out
    assert "<w:b/>" in out, "bold must survive"


def test_shrink_enables_autofit_and_clears_fixed_widths():
    out = shrink.shrink_table(TABLE_XML, 16)
    assert '<w:tblLayout w:type="autofit" />' in out
    assert 'w:gridCol w:w=' not in out, "fixed widths block autofit"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _run_all() -> int:
    tests = [value for name, value in sorted(globals().items())
             if name.startswith("test_") and callable(value)]
    failed = 0
    for test in tests:
        try:
            test()
        except AssertionError as error:
            failed += 1
            print(f"FAIL {test.__name__}\n  {error}\n")
        else:
            print(f"PASS {test.__name__}")
    if failed:
        print(f"\n{failed} of {len(tests)} docx tool checks failed")
        return 1
    print(f"\nALL {len(tests)} DOCX TOOL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(_run_all())
