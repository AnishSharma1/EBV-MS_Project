#!/usr/bin/env python3
"""Make wide tables legible in the .docx.

The LaTeX wrapped the widest tables in \\scriptsize and \\resizebox so they
fit the text block. Pandoc carries over neither, so a 9-column table gets
nine equal columns and every cell wraps mid-word. This shrinks the type in
tables past a column threshold and lets Word fit columns to their content,
which is what \\resizebox was doing by other means.

Edits only table runs; body text, captions and narrow tables are untouched.
"""

import argparse
import re
import shutil
import zipfile
from pathlib import Path

TBL_RE = re.compile(r"<w:tbl>.*?</w:tbl>", re.S)
GRIDCOL_RE = re.compile(r"<w:gridCol[^>]*/>")
# a run's properties block, where font size belongs
RPR_RE = re.compile(r"<w:rPr>(.*?)</w:rPr>", re.S)
SZ_RE = re.compile(r"<w:sz[^>]*/>|<w:szCs[^>]*/>")


def shrink_table(tbl: str, half_points: int) -> str:
    """Set every run in the table to the given size and enable autofit."""
    size = f'<w:sz w:val="{half_points}"/><w:szCs w:val="{half_points}"/>'

    def fix_rpr(m: re.Match[str]) -> str:
        inner = SZ_RE.sub("", m.group(1))
        return f"<w:rPr>{inner}{size}</w:rPr>"

    tbl = RPR_RE.sub(fix_rpr, tbl)
    # runs with no rPr at all need one
    tbl = re.sub(r"<w:r>(?!<w:rPr>)", f"<w:r><w:rPr>{size}</w:rPr>", tbl)

    # Let Word size columns to content rather than keeping them equal.
    tbl = tbl.replace(
        '<w:tblW w:type="auto" w:w="0" />',
        '<w:tblW w:type="auto" w:w="0" /><w:tblLayout w:type="autofit" />')
    tbl = GRIDCOL_RE.sub("<w:gridCol/>", tbl)
    return tbl


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx", type=Path)
    ap.add_argument("--min-cols", type=int, default=6,
                    help="only tables with at least this many columns")
    ap.add_argument("--half-points", type=int, default=16, help="16 = 8pt")
    args = ap.parse_args()

    src = zipfile.ZipFile(args.docx)
    doc = src.read("word/document.xml").decode("utf-8")

    changed = []

    def maybe(m: re.Match[str]) -> str:
        tbl = m.group(0)
        cols = len(GRIDCOL_RE.findall(tbl))
        if cols < args.min_cols:
            return tbl
        changed.append(cols)
        return shrink_table(tbl, args.half_points)

    new = TBL_RE.sub(maybe, doc)

    tmp = args.docx.with_suffix(".tmp.docx")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
        for item in src.infolist():
            data = src.read(item.filename)
            if item.filename == "word/document.xml":
                data = new.encode("utf-8")
            out.writestr(item, data)
    src.close()
    shutil.move(tmp, args.docx)

    print(f"{args.docx.name}: shrank {len(changed)} table(s) "
          f"with {changed} columns to {args.half_points/2:g}pt + autofit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
