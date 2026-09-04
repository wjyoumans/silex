#!/usr/bin/env python3
"""Compare two single-field class/unit JSON result files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


DEFAULT_PATHS = (
    "success",
    "class_group.order",
    "class_group.certification",
    "unit_group.free_rank",
    "unit_group.certification",
)


def get_path(data: dict[str, Any], path: str) -> Any:
    current: Any = data
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--actual", required=True, type=Path)
    parser.add_argument("--path", action="append", dest="paths")
    args = parser.parse_args()

    expected = json.loads(args.expected.read_text())
    actual = json.loads(args.actual.read_text())
    paths = args.paths or list(DEFAULT_PATHS)

    failures: list[str] = []
    for path in paths:
        expected_value = get_path(expected, path)
        actual_value = get_path(actual, path)
        if expected_value != actual_value:
            failures.append(
                f"{path}: expected {expected_value!r}, actual {actual_value!r}"
            )

    if failures:
        print("class/unit instance comparison failed")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print("class/unit instance comparison passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
