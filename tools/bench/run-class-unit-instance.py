#!/usr/bin/env python3
"""Run one class/unit instance with a process timeout and JSON output."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def find_field(manifest: dict[str, Any], field_id: str) -> dict[str, Any] | None:
    for collection in ("fields", "benchmark_warmups"):
        for row in manifest.get(collection, []):
            if row.get("id") == field_id:
                return row
    return None


def coeffs_from_field(row: dict[str, Any]) -> str:
    coeffs = row.get("coefficients_low_to_high")
    if not isinstance(coeffs, list) or not coeffs:
        raise ValueError("field row is missing coefficients_low_to_high")
    return ",".join(str(value) for value in coeffs)


def write_result(result: dict[str, Any], out: Path | None) -> None:
    text = json.dumps(result, indent=2) + "\n"
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    else:
        sys.stdout.write(text)


def timeout_stream_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def add_optional_arg(cmd: list[str], name: str, value: object | None) -> None:
    if value is not None:
        cmd.extend([name, str(value)])


def main() -> int:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe")
    parser.add_argument("--build-dir", type=Path, default=root / "build")
    parser.add_argument("--manifest", type=Path, default=root / "test/data/class_unit_fields.json")
    parser.add_argument("--field-id", help="field id from the manifest")
    parser.add_argument("--coeffs", help="low-to-high coefficients")
    warmup = parser.add_mutually_exclusive_group()
    warmup.add_argument("--warmup-field-id", help="unmeasured warmup field id")
    warmup.add_argument("--warmup-coeffs", help="unmeasured low-to-high coefficients")
    parser.add_argument("--mode", choices=["proven", "grh"])
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--precision", type=int)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-relations", type=int)
    parser.add_argument("--bf-cutoff", type=int)
    parser.add_argument("--factor-base-bound", type=int)
    parser.add_argument("--expect-class-order", type=int)
    parser.add_argument("--expect-unit-rank", type=int)
    parser.add_argument("--expect-success", choices=["0", "1", "false", "true"])
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--trace", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    manifest: dict[str, Any] | None = None
    field_row: dict[str, Any] | None = None
    if args.field_id:
        try:
            manifest = load_manifest(args.manifest)
        except OSError as exc:
            write_result(
                {"success": False, "error": f"could not read manifest: {exc}"},
                args.out,
            )
            return 2
        except json.JSONDecodeError as exc:
            write_result(
                {"success": False, "error": f"invalid manifest JSON: {exc}"},
                args.out,
            )
            return 2
        field_row = find_field(manifest, args.field_id)
        if field_row is None:
            write_result(
                {"success": False, "error": f"unknown field id: {args.field_id}"},
                args.out,
            )
            return 2
        if args.coeffs is None:
            try:
                args.coeffs = coeffs_from_field(field_row)
            except ValueError as exc:
                write_result({"success": False, "error": str(exc)}, args.out)
                return 2
        if args.mode is None and field_row.get("mode") is not None:
            args.mode = str(field_row["mode"])
        if args.timeout is None and field_row.get("timeout_seconds") is not None:
            args.timeout = float(field_row["timeout_seconds"])
        if args.expect_success is None and field_row.get("expected_success") is not None:
            args.expect_success = "1" if bool(field_row["expected_success"]) else "0"
        expecting_success = args.expect_success not in {"0", "false"}
        if expecting_success and args.expect_class_order is None and field_row.get("expected_class_order") is not None:
            args.expect_class_order = int(field_row["expected_class_order"])
        if expecting_success and args.expect_unit_rank is None and field_row.get("expected_unit_rank") is not None:
            args.expect_unit_rank = int(field_row["expected_unit_rank"])

    if args.warmup_field_id:
        try:
            if manifest is None:
                manifest = load_manifest(args.manifest)
            warmup_row = find_field(manifest, args.warmup_field_id)
            if warmup_row is None:
                raise ValueError(f"unknown warmup field id: {args.warmup_field_id}")
            args.warmup_coeffs = coeffs_from_field(warmup_row)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            write_result({"success": False, "error": str(exc)}, args.out)
            return 2

    if args.coeffs is None:
        write_result(
            {"success": False, "error": "missing --coeffs or --field-id"},
            args.out,
        )
        return 2
    if args.mode is None:
        args.mode = "proven"
    if args.timeout is None:
        args.timeout = 60.0

    exe = Path(args.exe) if args.exe else args.build_dir / "silex-class-unit-instance"
    if not exe.exists():
        write_result(
            {
                "success": False,
                "error": f"missing executable: {exe}",
                "hint": "build the silex-class-unit-instance target first",
                "field_id": args.field_id,
            },
            args.out,
        )
        return 2

    cmd = [str(exe), "--coeffs", args.coeffs, "--mode", args.mode]
    add_optional_arg(cmd, "--warmup-coeffs", args.warmup_coeffs)
    add_optional_arg(cmd, "--precision", args.precision)
    add_optional_arg(cmd, "--max-candidates", args.max_candidates)
    add_optional_arg(cmd, "--max-relations", args.max_relations)
    add_optional_arg(cmd, "--bf-cutoff", args.bf_cutoff)
    add_optional_arg(cmd, "--factor-base-bound", args.factor_base_bound)
    add_optional_arg(cmd, "--expect-class-order", args.expect_class_order)
    add_optional_arg(cmd, "--expect-unit-rank", args.expect_unit_rank)
    add_optional_arg(cmd, "--expect-success", args.expect_success)
    if args.log:
        cmd.append("--log")
    if args.trace:
        cmd.append("--trace")
    if args.verbose:
        cmd.append("--verbose")
    if args.profile:
        cmd.append("--profile")

    try:
        completed = subprocess.run(
            cmd,
            cwd=root,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=args.timeout,
        )
    except subprocess.TimeoutExpired as exc:
        result = {
            "success": False,
            "timeout": True,
            "timeout_seconds": args.timeout,
            "field_id": args.field_id,
            "cmd": cmd,
            "stdout": timeout_stream_text(exc.stdout),
            "stderr": timeout_stream_text(exc.stderr),
        }
        write_result(result, args.out)
        return 124

    text = completed.stdout
    if not text.strip():
        result = {
            "success": False,
            "returncode": completed.returncode,
            "field_id": args.field_id,
            "cmd": cmd,
            "stderr": completed.stderr,
        }
        text = json.dumps(result, indent=2) + "\n"

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    else:
        sys.stdout.write(text)
    if completed.stderr:
        sys.stderr.write(completed.stderr)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
