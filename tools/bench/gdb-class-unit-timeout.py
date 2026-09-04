#!/usr/bin/env python3
"""Capture a gdb backtrace for one class/unit instance after a timeout."""

from __future__ import annotations

import argparse
import json
import re
import shutil
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


def load_field(path: Path, field_id: str) -> dict[str, Any]:
    data = json.loads(path.read_text())
    for row in data.get("fields", []):
        if row.get("id") == field_id:
            return row
    raise ValueError(f"unknown field id: {field_id}")


def timeout_arg(seconds: float) -> str:
    if seconds == int(seconds):
        return f"{int(seconds)}s"
    return f"{seconds}s"


def main() -> int:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=root / "test/data/class_unit_fields.json")
    parser.add_argument("--field-id")
    parser.add_argument("--coeffs")
    parser.add_argument(
        "--mode", choices=["proven", "grh", "candidate", "unknown"]
    )
    parser.add_argument("--build-dir", type=Path, default=root / "build")
    parser.add_argument("--exe")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--precision", type=int)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-relations", type=int)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--out")
    args = parser.parse_args()

    field_label = "field"
    if args.field_id:
        try:
            field_label = validated_field_id(args.field_id)
            row = load_field(args.manifest, field_label)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"gdb class/unit setup failed: {exc}", file=sys.stderr)
            return 2
        if args.coeffs is None:
            args.coeffs = ",".join(str(value) for value in row["coefficients_low_to_high"])
        if args.mode is None:
            args.mode = str(row.get("mode", "proven"))
        if "timeout_seconds" in row and args.timeout == 60.0:
            args.timeout = float(row["timeout_seconds"])

    if args.coeffs is None:
        print("missing --coeffs or --field-id", file=sys.stderr)
        return 2
    if args.mode is None:
        args.mode = "proven"

    exe = Path(args.exe) if args.exe else args.build_dir / "silex-class-unit-instance"
    if not exe.exists():
        print(f"missing executable: {exe}", file=sys.stderr)
        return 2
    if shutil.which("timeout") is None:
        print("missing timeout executable", file=sys.stderr)
        return 2
    if shutil.which("gdb") is None:
        print("missing gdb executable", file=sys.stderr)
        return 2

    program_args = [str(exe), "--coeffs", args.coeffs, "--mode", args.mode]
    if args.precision is not None:
        program_args.extend(["--precision", str(args.precision)])
    if args.max_candidates is not None:
        program_args.extend(["--max-candidates", str(args.max_candidates)])
    if args.max_relations is not None:
        program_args.extend(["--max-relations", str(args.max_relations)])
    if args.profile:
        program_args.append("--profile")

    cmd = [
        "timeout",
        "-s",
        "INT",
        timeout_arg(args.timeout),
        "gdb",
        "-batch",
        "-ex",
        "set debuginfod enabled off",
        "-ex",
        "set pagination off",
        "-ex",
        "run",
        "-ex",
        "thread apply all bt",
        "--args",
        *program_args,
    ]

    out = Path(args.out) if args.out else root / "build/class-unit-backtraces" / f"{field_label}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(cmd, cwd=root, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out.write_text(
        "$ " + " ".join(cmd) + "\n\n"
        + "returncode=" + str(completed.returncode) + "\n\n"
        + "stdout:\n" + completed.stdout + "\n"
        + "stderr:\n" + completed.stderr + "\n"
    )
    print(f"wrote gdb timeout backtrace output to {out.relative_to(root)}")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
