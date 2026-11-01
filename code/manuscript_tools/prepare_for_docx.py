#!/usr/bin/env python3
"""Preprocess the V20 LaTeX sources so pandoc produces a faithful .docx.

Pandoc's LaTeX reader silently drops \\cite when the bibliography is a
thebibliography environment, flattens that environment into one paragraph,
and leaves \\includegraphics pointing at PDFs that Word cannot display.
It also cannot honour \\renewcommand{\\thefigure}{S\\arabic{figure}}.

This rewrites the source to resolve those four things before conversion.
Nothing here edits the tracked sources; output goes to a build directory.
"""

import argparse
import re
import shutil
from pathlib import Path

CITE_RE = re.compile(r"\\cite\{([^}]*)\}")
BIBITEM_RE = re.compile(r"\\bibitem\{([^}]*)\}")


def balanced_arg(tex: str, open_brace: int) -> tuple[str, int]:
    """Return the content of the braced group at open_brace and its end index."""
    assert tex[open_brace] == "{", tex[open_brace : open_brace + 20]
    depth = 0
    i = open_brace
    while i < len(tex):
        c = tex[i]
        if c == "\\" and i + 1 < len(tex):
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return tex[open_brace + 1 : i], i + 1
        i += 1
    raise ValueError(f"unbalanced braces from offset {open_brace}")


def unwrap_macro(tex: str, macro: str, nargs: int, keep: int,
                 transform=None) -> tuple[str, int]:
    """Replace \\macro{..}..{..} with its ``keep``-th argument (0-based).

    Pandoc's LaTeX reader drops the whole enclosing environment when it meets
    a box macro it does not model, which silently deletes tables. Unwrapping
    to the payload keeps the content and loses only the scaling.
    """
    count = 0
    while True:
        at = tex.find("\\" + macro)
        if at < 0:
            return tex, count
        i = at + len(macro) + 1
        args = []
        for _ in range(nargs):
            while i < len(tex):
                # skip whitespace, a trailing comment, and optional [..] args
                if tex[i] in " \t\n":
                    i += 1
                elif tex[i] == "%":
                    i = tex.index("\n", i) + 1
                elif tex[i] == "[":
                    i = tex.index("]", i) + 1
                else:
                    break
            arg, i = balanced_arg(tex, i)
            args.append(arg)
        payload = args[keep]
        if transform is not None:
            payload = transform(payload)
        tex = tex[:at] + payload + tex[i:]
        count += 1


def citation_numbers(references_tex: str) -> dict[str, int]:
    """Map each \\bibitem key to its 1-based position, as LaTeX numbers them."""
    return {key: n for n, key in enumerate(BIBITEM_RE.findall(references_tex), 1)}


def bibliography_paragraphs(references_tex: str) -> str:
    """Rewrite thebibliography as one numbered paragraph per reference.

    Pandoc merges the environment's items into a single run-on paragraph, so
    the entries are emitted as separate paragraphs carrying explicit [n]
    labels instead.
    """
    body = references_tex.split(r"\begin{thebibliography}", 1)[1]
    body = body.split(r"\end{thebibliography}", 1)[0]
    body = body.split("}", 1)[1]  # drop the {99} widest-label argument

    entries = [chunk.strip() for chunk in BIBITEM_RE.split(body)]
    entries = entries[1:]  # leading text before the first \bibitem is empty

    out = ["\\section*{References}", ""]
    for n, (_key, text) in enumerate(zip(entries[::2], entries[1::2]), 1):
        collapsed = " ".join(text.split())
        out.append(f"[{n}] {collapsed}")
        out.append("")
    return "\n".join(out)


def resolve_citations(tex: str, numbers: dict[str, int]) -> tuple[str, list[str]]:
    """Replace \\cite{a,b} with the [1, 2] form the compiled PDF renders."""
    missing: list[str] = []

    def repl(match: re.Match[str]) -> str:
        keys = [k.strip() for k in match.group(1).split(",") if k.strip()]
        nums = []
        for key in keys:
            if key not in numbers:
                missing.append(key)
                continue
            nums.append(numbers[key])
        return "[" + ", ".join(str(n) for n in nums) + "]"

    return CITE_RE.sub(repl, tex), missing


def number_captions(tex: str, prefix: str) -> tuple[str, int, int, dict[str, str]]:
    """Prepend 'Figure N.' / 'Table N.' to each numbered caption.

    Counters advance in source order, exactly as LaTeX does. \\caption*{} is
    unnumbered in LaTeX and is left alone. ``prefix`` is "" for the main text
    and "S" for the supplement, which renumbers both counters as S1, S2, ...

    Also returns the \\label -> number map, because pandoc numbers only the
    floats that carry a \\label and so renders \\ref with the wrong value.
    """
    figures = tables = 0
    labels: dict[str, str] = {}
    out: list[str] = []
    # Split on the environment openers, captions and labels, keeping order.
    tokens = re.split(
        r"(\\begin\{figure\}|\\begin\{table\}|\\caption\*?\{|\\label\{[^}]*\})",
        tex)
    current = None      # which float we are inside
    number = None       # the number that float was given
    for token in tokens:
        if token == r"\begin{figure}":
            current, number = "figure", None
        elif token == r"\begin{table}":
            current, number = "table", None
        elif token == r"\caption{":
            if current == "figure":
                figures += 1
                number = f"{prefix}{figures}"
                token = token + f"Figure {number}. "
            elif current == "table":
                tables += 1
                number = f"{prefix}{tables}"
                token = token + f"Table {number}. "
        elif token.startswith(r"\label{"):
            key = token[len(r"\label{"):-1]
            if number is not None:
                labels[key] = number
        out.append(token)
    return "".join(out), figures, tables, labels


def resolve_refs(tex: str, labels: dict[str, str]) -> tuple[str, list[str]]:
    """Replace \\ref{key} with the float's number."""
    unknown: list[str] = []

    def repl(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in labels:
            unknown.append(key)
            return match.group(0)
        return labels[key]

    # ~ before \ref is a non-breaking space in the source.
    tex = re.sub(r"~(?=\\ref\{)", " ", tex)
    return re.sub(r"\\ref\{([^}]*)\}", repl, tex), unknown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="main.tex or supplement.tex")
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--caption-prefix", default="", help='"S" for the supplement')
    args = parser.parse_args()

    tex = args.source.read_text(encoding="utf-8")
    references_tex = args.references.read_text(encoding="utf-8")
    numbers = citation_numbers(references_tex)

    tex, missing = resolve_citations(tex, numbers)
    if missing:
        raise SystemExit(f"unknown citation keys: {sorted(set(missing))}")

    cites_left = CITE_RE.search(tex)
    if cites_left:
        raise SystemExit(f"unresolved \\cite remains: {cites_left.group(0)}")

    # \resizebox makes pandoc discard the enclosing table entirely; \shortstack
    # cells come through empty. Unwrap both to their text payload.
    tex, n_resize = unwrap_macro(tex, "resizebox", 3, keep=2)
    # \shortstack uses \\ as its internal line break; left in place it would
    # end the table row, so it becomes a space.
    tex, n_stack = unwrap_macro(
        tex, "shortstack", 1, keep=0,
        transform=lambda s: " ".join(s.replace(r"\\", " ").split()))

    # Word cannot render PDF images; the figures were rasterised to PNG.
    tex = tex.replace(".pdf}", ".png}")

    tex, figures, tables, labels = number_captions(tex, args.caption_prefix)

    tex, unknown = resolve_refs(tex, labels)
    if unknown:
        raise SystemExit(f"\\ref to unknown labels: {sorted(set(unknown))}")

    tex = tex.replace(
        r"\input{references.tex}", bibliography_paragraphs(references_tex))
    if r"\input{" in tex:
        raise SystemExit("an \\input other than references.tex is unhandled")

    args.outdir.mkdir(parents=True, exist_ok=True)
    target = args.outdir / args.source.name
    target.write_text(tex, encoding="utf-8")

    print(f"{args.source.name}: {len(numbers)} references, "
          f"{figures} figures, {tables} tables, "
          f"{n_resize} resizebox and {n_stack} shortstack unwrapped, "
          f"refs resolved {labels} -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
