#!/usr/bin/env python3
"""Check the generated .docx carries the same prose as the compiled V20 PDF.

Compares against manuscript/paper_v20_*.txt, which is pdftotext of the
tracked PDF snapshot, so the PDF is the reference for what V20 says.
Reports sentences present in one and not the other.
"""

import argparse
import collections
import re
import sys
import zipfile
from pathlib import Path


def docx_text(path: Path) -> str:
    """Visible text in document order, including OMML equation characters.

    Equations live in <m:t>, not <w:t>; reading only <w:t> makes every
    formula look like it went missing.
    """
    d = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    # A table caption is carried in a w:tblPr attribute; drop attributes so
    # their text is not mistaken for body content.
    d = re.sub(r"<w:tblPr>.*?</w:tblPr>", "", d, flags=re.S)
    d = re.sub(r"<w:tab[^>]*/>", " ", d)
    d = d.replace("</w:p>", "</w:p>\n").replace("</w:tc>", " </w:tc>")

    out = []
    for para in d.split("\n"):
        out.append("".join(
            m.group(1) for m in
            re.finditer(r"<(?:w|m):t(?:\s[^>]*)?>(.*?)</(?:w|m):t>",
                        para, re.S)))
    text = "\n".join(out)
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
                 ("&quot;", '"'), ("&apos;", "'")):
        text = text.replace(a, b)
    return text


def normalise(text: str) -> str:
    """Flatten the differences that are purely typesetting."""
    text = text.replace("­", "")
    # pdftotext hyphenates across line breaks
    text = re.sub(r"-\n(?=[a-z])", "", text)
    text = text.replace("\n", " ")
    # unify the dashes, quotes and spaces the two renderers choose differently
    for ch in "‐‑‒–—−":
        text = text.replace(ch, "-")
    for ch in "‘’ʼ":
        text = text.replace(ch, "'")
    for ch in "“”":
        text = text.replace(ch, '"')
    text = text.replace(" ", " ").replace(" ", " ")
    text = text.replace("Å", "A")           # angstrom sign vs letter
    text = re.sub(r"\s+", " ", text)
    return text


def sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[.:;])\s+", text)
    return [p.strip() for p in parts if len(p.strip()) >= 40]


def fingerprint(sentence: str) -> str:
    """Compare on letters and digits only, so markup and spacing drop out.

    Also drops the section and page numbers LaTeX typesets into the PDF but
    Word carries as heading styles, which are not prose differences.
    """
    return re.sub(r"[^a-z0-9]", "", sentence.lower())


def strip_pdf_furniture(text: str) -> str:
    """Remove LaTeX section numbers and stray page numbers from the PDF text."""
    # "2.3 TCR-facing ..." / "1 Background" at a heading
    text = re.sub(r"(?<![\w.])\d+(?:\.\d+)*\s+(?=[A-Z][a-z])", "", text)
    return text


def words(text: str) -> collections.Counter:
    """Multiset of alphabetic words, which excludes page and section numbers."""
    return collections.Counter(
        w for w in re.findall(r"[A-Za-z][A-Za-z'-]{3,}", text))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx", type=Path)
    ap.add_argument("pdf_text", type=Path)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    doc = normalise(docx_text(args.docx))
    ref = normalise(args.pdf_text.read_text(encoding="utf-8", errors="replace"))

    # Primary signal: prose words in the PDF that the docx does not carry.
    dw, rw = words(doc), words(ref)
    dropped = rw - dw
    extra = dw - rw

    # Secondary: sentence-level, with section numbering removed from both so
    # LaTeX's auto-numbered headings do not register as content differences.
    doc_flat = fingerprint(strip_pdf_furniture(doc))
    ref_flat = fingerprint(strip_pdf_furniture(ref))
    missing = [s for s in sentences(strip_pdf_furniture(ref))
               if fingerprint(s) not in doc_flat]

    print(f"=== {args.label or args.docx.name} ===")
    print(f"  reference words      : {sum(rw.values())}")
    print(f"  docx words           : {sum(dw.values())}")
    print(f"  words dropped        : {sum(dropped.values())}"
          f"  {dict(dropped.most_common(12)) if dropped else ''}")
    print(f"  words added          : {sum(extra.values())}"
          f"  {dict(extra.most_common(12)) if extra else ''}")
    print(f"  sentences not found  : {len(missing)}"
          f" / {len(sentences(ref))}")
    for s in missing[:10]:
        print(f"      - {s[:160]}")
    return 1 if sum(dropped.values()) or missing else 0


if __name__ == "__main__":
    sys.exit(main())
