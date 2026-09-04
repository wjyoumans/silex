#!/usr/bin/env python3
"""Run manifest class/unit rows that should pass inside the fast gate."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

FIELD_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def validated_field_id(value: object) -> str:
    if not isinstance(value, str) or FIELD_ID_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "field id must match [A-Za-z0-9][A-Za-z0-9_.-]*"
        )
    return value


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_fields(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text())
    fields = data.get("fields", [])
    if not isinstance(fields, list):
        raise ValueError("manifest field 'fields' must be a list")
    for row in fields:
        if not isinstance(row, dict):
            raise ValueError("manifest fields must contain objects")
        validated_field_id(row.get("id"))
    return fields


def selected_fields(
    fields: list[dict[str, Any]],
    field_ids: list[str],
    statuses: set[str],
) -> list[dict[str, Any]]:
    if field_ids:
        by_id = {validated_field_id(row.get("id")): row for row in fields}
        missing = [field_id for field_id in field_ids if field_id not in by_id]
        if missing:
            raise ValueError("unknown field id(s): " + ", ".join(missing))
        return [by_id[field_id] for field_id in field_ids]
    return [row for row in fields if str(row.get("status")) in statuses]


def main() -> int:
    root = repo_root()
    root_resolved = root.resolve()
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=root / "test/data/class_unit_fields.json")
    parser.add_argument("--build-dir", type=Path, default=root / "build")
    parser.add_argument("--exe")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--status", action="append")
    parser.add_argument("--field-id", action="append", default=[])
    parser.add_argument("--out-dir", type=Path, default=root / "build/class-unit-priority-gate")
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()

    try:
        fields = load_fields(args.manifest)
        statuses = set(args.status or ["must_pass_fast"])
        rows = selected_fields(fields, args.field_id, statuses)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"class/unit priority gate setup failed: {exc}", file=sys.stderr)
        return 2

    if not rows:
        print("class/unit priority gate: no selected rows")
        return 0

    script = root / "tools/bench/run-class-unit-instance.py"
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    for row in rows:
        field_id = validated_field_id(row.get("id"))
        timeout = float(row.get("timeout_seconds", args.timeout))
        out = args.out_dir / f"{field_id}.json"
        cmd = [
            sys.executable,
            str(script),
            "--manifest",
            str(args.manifest),
            "--field-id",
            field_id,
            "--build-dir",
            str(args.build_dir),
            "--timeout",
            str(timeout),
            "--expect-success",
            "1",
            "--out",
            str(out),
        ]
        if args.exe:
            cmd.extend(["--exe", args.exe])
        if args.log:
            cmd.append("--log")
        if args.verbose:
            cmd.append("--verbose")
        if args.profile:
            cmd.append("--profile")

        completed = subprocess.run(cmd, cwd=root, check=False, text=True)
        try:
            result = json.loads(out.read_text())
        except (OSError, json.JSONDecodeError):
            result = {"success": False, "error": "missing or invalid result JSON"}

        status = "pass" if completed.returncode == 0 else "fail"
        elapsed = result.get("elapsed_ms")
        elapsed_text = f" elapsed_ms={elapsed:.3f}" if isinstance(elapsed, float) else ""
        try:
            result_path = out.relative_to(root_resolved)
        except ValueError:
            result_path = out
        print(f"{status}: {field_id}{elapsed_text} result={result_path}")
        if completed.returncode != 0:
            failures.append(field_id)

    if failures:
        print("class/unit priority gate failed: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
