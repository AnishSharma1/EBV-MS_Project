#!/usr/bin/env python3
"""Build editable .docx copies of the manuscript and supplement.

Run from any directory::

    python3 code/manuscript_tools/build_docx.py

Writes ``manuscript/paper_v<N>_main.docx`` and ``..._supplement.docx``.
The LaTeX sources are read only; nothing tracked is modified.

Pandoc cannot convert these sources directly without losing content -- see
``README.md`` in this directory for the four silent failures this pipeline
works around. The steps are:

1. rasterise the figure PDFs, which Word cannot display, to PNG;
2. rewrite the LaTeX (citations, bibliography, float numbering, cross
   references) with ``prepare_for_docx.py``;
3. run pandoc against a reference document that carries the manuscript font;
4. shrink the widest tables with ``shrink_wide_tables.py``.

Requires ``pandoc`` and ``pdftoppm`` (Poppler) on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Figures are rasterised at this resolution; 300 dpi keeps the vector line art
# crisp when Word scales it into the text block.
FIGURE_DPI = 300
# Tables at or beyond this many columns get the smaller type that \resizebox
# was providing in LaTeX.
WIDE_TABLE_COLUMNS = 6
TABLE_HALF_POINTS = 16  # 8pt


def find_root() -> Path:
    for candidate in (HERE, *HERE.parents):
        if (candidate / "manuscript").is_dir() and (candidate / "code").is_dir():
            return candidate
    raise SystemExit("could not locate the repository root")


def require(tool: str) -> None:
    if shutil.which(tool) is None:
        raise SystemExit(f"{tool} is not on PATH; install it and re-run")


def current_version(manuscript: Path) -> str:
    """Latest paper_v<N>_main.pdf snapshot in the manuscript directory."""
    versions = []
    for pdf in manuscript.glob("paper_v*_main.pdf"):
        tag = pdf.name[len("paper_v"):].split("_", 1)[0]
        if tag.isdigit():
            versions.append(int(tag))
    if not versions:
        raise SystemExit("no paper_v<N>_main.pdf snapshot found")
    return str(max(versions))


def rasterise_figures(manuscript: Path, build: Path) -> int:
    """Convert every figure the two sources include into a PNG beside it."""
    import re

    count = 0
    pattern = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")
    for name in ("main.tex", "supplement.tex"):
        text = (manuscript / name).read_text(encoding="utf-8")
        for relative in pattern.findall(text):
            source = manuscript / relative
            if source.suffix.lower() != ".pdf":
                continue
            target = build / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                ["pdftoppm", "-png", "-r", str(FIGURE_DPI), "-singlefile",
                 str(source), str(target.with_suffix(""))],
                check=True)
            count += 1
    return count


def reference_docx(build: Path, font: str) -> Path:
    """Pandoc's default reference document, restyled to the manuscript font.

    Word takes heading fonts from the theme and body fonts from the styles,
    so both have to be set or the headings come out in pandoc's sans default.
    """
    import re

    base = build / "reference_base.docx"
    base.write_bytes(subprocess.run(
        ["pandoc", "--print-default-data-file", "reference.docx"],
        check=True, capture_output=True).stdout)

    unpacked = build / "reference_unpacked"
    with zipfile.ZipFile(base) as z:
        z.extractall(unpacked)

    styles = unpacked / "word" / "styles.xml"
    text = styles.read_text(encoding="utf-8")
    for attr in ("ascii", "hAnsi", "cs"):
        text = re.sub(rf'w:{attr}="(?!{font})[^"]*"', f'w:{attr}="{font}"', text)
    styles.write_text(text, encoding="utf-8")

    theme = unpacked / "word" / "theme" / "theme1.xml"
    if theme.exists():
        text = theme.read_text(encoding="utf-8")
        text = re.sub(
            r'(<a:(?:majorFont|minorFont)>\s*<a:latin[^>]*typeface=")[^"]*(")',
            rf'\1{font}\2', text)
        theme.write_text(text, encoding="utf-8")

    out = build / "reference.docx"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(unpacked.rglob("*")):
            if path.is_file():
                z.write(path, path.relative_to(unpacked).as_posix())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--font", default="Times New Roman",
                    help="body and heading font for the .docx")
    ap.add_argument("--keep-build", action="store_true",
                    help="keep the intermediate build directory")
    args = ap.parse_args()

    require("pandoc")
    require("pdftoppm")

    root = find_root()
    manuscript = root / "manuscript"
    version = current_version(manuscript)

    build_dir = Path(tempfile.mkdtemp(prefix="docx_build_"))
    try:
        figures = rasterise_figures(manuscript, build_dir)
        print(f"rasterised {figures} figure(s) at {FIGURE_DPI} dpi")

        reference = reference_docx(build_dir, args.font)
        print(f"reference document styled as {args.font}")

        written = []
        for name, prefix in (("main", ""), ("supplement", "S")):
            source = manuscript / f"{name}.tex"
            cmd = [sys.executable, "-I", str(HERE / "prepare_for_docx.py"),
                   str(source),
                   "--references", str(manuscript / "references.tex"),
                   "--outdir", str(build_dir)]
            if prefix:
                cmd += ["--caption-prefix", prefix]
            subprocess.run(cmd, check=True)

            target = manuscript / f"paper_v{version}_{name}.docx"
            subprocess.run(
                ["pandoc", f"{name}.tex", "-f", "latex", "-t", "docx",
                 "--reference-doc", str(reference),
                 "--resource-path", f".:{build_dir}",
                 "-o", str(target)],
                check=True, cwd=build_dir)

            subprocess.run(
                [sys.executable, "-I", str(HERE / "shrink_wide_tables.py"),
                 str(target),
                 "--min-cols", str(WIDE_TABLE_COLUMNS),
                 "--half-points", str(TABLE_HALF_POINTS)],
                check=True)
            written.append(target)

        print("\nwrote:")
        for path in written:
            print(f"  {path.relative_to(root)}  ({path.stat().st_size:,} bytes)")
        print("\nVerify against the compiled PDF with:")
        for path in written:
            name = path.stem.split("_", 2)[-1]
            print(f"  python3 code/manuscript_tools/compare_to_pdf.py "
                  f"{path.relative_to(root)} "
                  f"manuscript/paper_v{version}_{name}.txt")
    finally:
        if args.keep_build:
            print(f"\nbuild directory kept at {build_dir}")
        else:
            shutil.rmtree(build_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
