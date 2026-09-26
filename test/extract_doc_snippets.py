#!/usr/bin/env python3
"""Extract C++ code blocks from the Sphinx sources into compilable units.

Every ``code-block:: cpp`` (or ``c++``) block in ``docs/**/*.rst`` becomes
part of one generated translation unit per document.  The blocks of a
document are placed, in order, in nested scopes inside a single ``bool``
function, so a later block can use names declared by an earlier one (as a
tutorial does) and can also redeclare them (as a self-contained reference
example does).  ``#include`` lines are hoisted to file scope.  A document
whose blocks show no ``#include`` lines is compiled against the umbrella
header ``<silex/silex.hpp>`` and ``<cassert>``; a document whose blocks do
show includes must compile with exactly those includes.

``#line`` directives map compiler diagnostics back to the ``.rst`` lines.

The build compiles the generated files as an object library; see the
``silex-doc-snippets`` test in the top-level ``CMakeLists.txt``.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


DIRECTIVE = re.compile(
    r"^(?P<indent>\s*)\.\.\s+(code-block|code|sourcecode)::\s*(cpp|c\+\+)\s*$")
OPTION = re.compile(r"^\s*:[\w-]+:")
INCLUDE = re.compile(r"^\s*#\s*include\b")


class Block:
    def __init__(self, first_line: int, lines: list[str]) -> None:
        self.first_line = first_line
        self.lines = lines


def indentation(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def extract_blocks(path: Path) -> list[Block]:
    lines = path.read_text(encoding="utf-8").splitlines()
    blocks: list[Block] = []
    i = 0
    while i < len(lines):
        match = DIRECTIVE.match(lines[i])
        if match is None:
            i += 1
            continue
        directive_indent = len(match.group("indent"))
        i += 1
        while i < len(lines) and OPTION.match(lines[i]):
            i += 1
        body: list[tuple[int, str]] = []
        while i < len(lines):
            line = lines[i]
            if line.strip() and indentation(line) <= directive_indent:
                break
            body.append((i + 1, line))
            i += 1
        while body and not body[0][1].strip():
            body.pop(0)
        while body and not body[-1][1].strip():
            body.pop()
        if not body:
            raise ValueError(f"{path}: empty C++ code block")
        content_indent = min(indentation(text) for _, text in body
                             if text.strip())
        blocks.append(Block(body[0][0],
                            [text[content_indent:] for _, text in body]))
    return blocks


def unit_name(rel: Path) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(rel.with_suffix("")))


def render(rel: Path, blocks: list[Block]) -> str:
    includes: list[str] = []
    for block in blocks:
        for line in block.lines:
            if INCLUDE.match(line) and line.strip() not in includes:
                includes.append(line.strip())
    if not includes:
        includes = ["#include <cassert>", "#include <silex/silex.hpp>"]

    out = [
        f"// Generated from {rel.as_posix()} by test/extract_doc_snippets.py.",
        "// Do not edit; change the documentation instead.",
        "",
        *includes,
        "",
        "namespace {",
        "",
        f"[[maybe_unused]] bool doc_snippet_{unit_name(rel)}() {{",
    ]
    for block in blocks:
        # Hoisted includes become blank lines, so one #line per block keeps
        # every later line in step with the .rst source.
        out.append("{")
        out.append(f'#line {block.first_line} "{rel.as_posix()}"')
        for line in block.lines:
            out.append("" if INCLUDE.match(line) else line)
    out.append("}" * len(blocks))
    out.append("    return true;")
    out.append("}")
    out.append("")
    out.append("}  // namespace")
    out.append("")
    return "\n".join(out)


def write_if_changed(path: Path, text: str) -> None:
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True,
                        help="source root containing docs/")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--list-file", type=Path, required=True,
                        help="file receiving the generated source paths")
    args = parser.parse_args()

    root = args.root.resolve()
    output_dir = args.output_dir.resolve()
    generated: list[Path] = []
    block_count = 0
    try:
        for doc in sorted((root / "docs").rglob("*.rst")):
            blocks = extract_blocks(doc)
            if not blocks:
                continue
            rel = doc.relative_to(root)
            target = output_dir / (unit_name(rel) + ".cpp")
            write_if_changed(target, render(rel, blocks))
            generated.append(target)
            block_count += len(blocks)
    except (OSError, ValueError) as error:
        print(f"extract_doc_snippets: {error}", file=sys.stderr)
        return 1

    for stale in output_dir.glob("*.cpp") if output_dir.exists() else []:
        if stale not in generated:
            stale.unlink()

    write_if_changed(args.list_file,
                     "".join(f"{path.as_posix()}\n" for path in generated))
    print(f"extract_doc_snippets: {block_count} C++ blocks from "
          f"{len(generated)} documents")
    return 0


if __name__ == "__main__":
    sys.exit(main())
