#!/usr/bin/env python3
"""Lightweight restricted C++ profile checks for Silex."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


TOKEN_CHECKS = (
    ("exceptions", re.compile(r"\b(throw|try|catch)\b")),
    ("rtti", re.compile(r"\b(dynamic_cast|typeid)\b")),
    ("casts requiring review", re.compile(r"\b(reinterpret_cast|const_cast)\s*<")),
)

PUBLIC_STL_PATTERNS = (
    re.compile(r"\bstd::(string|vector|map|unordered_map|set|unordered_set|span)\b"),
    re.compile(r"^\s*#\s*include\s*<(string|vector|map|unordered_map|set|unordered_set)>",
               re.MULTILINE),
)

PUBLIC_STL_ALLOWLIST = {
    "include/silex/class_group.hpp",
    "include/silex/detail/class_relation_module_context.hpp",
    "include/silex/factor_base.hpp",
    "include/silex/factored_element.hpp",
    "include/silex/factored.hpp",
    "include/silex/fmpz_smat.hpp",
    "include/silex/ideal_factorization.hpp",
    "include/silex/relation.hpp",
}

CLASS_GROUP_PUBLIC_INTERNAL_PATTERNS = (
    re.compile(r"\bPartialRelationEntry\b"),
    re.compile(r"\bRelationAdmissionCache\b"),
    re.compile(r"\bSubfactorBaseSchedule\b"),
    re.compile(r"\bRelationCompletionState\b"),
    re.compile(r"\bFactorBaseHonestyScanAudit\b"),
    re.compile(r"\bRelationSearchControlState\b"),
    re.compile(r"\bRelationFinishState\b"),
    re.compile(r"\bRandomIdealSearchState\b"),
    re.compile(r"\bClassUnitRunAudit\b"),
    re.compile(r"\bClassUnitFallback(Stage|Kind)\b"),
    re.compile(r"\btry_certify_class_unit_with_bf_audit\b"),
    re.compile(
        r"\b(ideal_lattice_short_element|weighted_ideal_lattice_short_element|"
        r"reduce_ideal_(lattice|product|signed_power))\b"
    ),
    re.compile(r"\bbuild_weighted_ideal_lattice_reduction\b"),
    re.compile(
        r"\b(factor_base_honesty_|initialize_relation_admission_cache|"
        r"next_relation_random_exponent|"
        r"rebuild_relation_factor_base_and_replay)\w*\b"
    ),
    re.compile(r"\brecord_class_unit_fallback_\b"),
)


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//.*", "", text)
    return text


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def check_tokens(path: Path, rel: str, failures: list[str]) -> None:
    raw = path.read_text(encoding="utf-8")
    text = strip_comments(raw)
    for label, pattern in TOKEN_CHECKS:
        for match in pattern.finditer(text):
            failures.append(f"{rel}:{line_number(text, match.start())}: forbidden {label}: "
                            f"{match.group(0)}")


def check_public_header(path: Path, rel: str, failures: list[str]) -> None:
    if rel in PUBLIC_STL_ALLOWLIST:
        return
    raw = path.read_text(encoding="utf-8")
    text = strip_comments(raw)
    for pattern in PUBLIC_STL_PATTERNS:
        for match in pattern.finditer(text):
            failures.append(f"{rel}:{line_number(text, match.start())}: public header exposes "
                            f"restricted STL ABI surface: {match.group(0)}")


def check_class_group_public_boundary(
        path: Path, rel: str, failures: list[str]) -> None:
    raw = path.read_text(encoding="utf-8")
    text = strip_comments(raw)
    for pattern in CLASS_GROUP_PUBLIC_INTERNAL_PATTERNS:
        for match in pattern.finditer(text):
            failures.append(
                f"{rel}:{line_number(text, match.start())}: class-group public header "
                f"exposes private implementation contract: {match.group(0)}"
            )


def iter_existing(root: Path, patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(root.glob(pattern))
    return sorted(path for path in paths if path.is_file())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()

    root = args.root.resolve()
    checked_paths = iter_existing(root, [
        "src/**/*.cpp",
        "include/silex/**/*.hpp",
        "include/silex/**/*.h",
    ])
    public_headers = iter_existing(root, [
        "include/silex/**/*.hpp",
        "include/silex/**/*.h",
    ])

    failures: list[str] = []
    for path in checked_paths:
        check_tokens(path, str(path.relative_to(root)), failures)
    for path in public_headers:
        rel = str(path.relative_to(root))
        check_public_header(path, rel, failures)
        if rel == "include/silex/class_group.hpp":
            check_class_group_public_boundary(path, rel, failures)

    if failures:
        print("Restricted C++ profile check failed:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    print(f"Restricted C++ profile check passed ({len(checked_paths)} files scanned).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
