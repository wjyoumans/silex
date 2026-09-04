#!/usr/bin/env python3
"""Exercise the public-API backend operation instance."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )

SCOPES = {
    "maximal_order": "field_construction+equation_order+maximal_order",
    "ideal_multiply": "ideal_multiply_only",
    "element_square_root": "element_is_square_with_root_only",
}
READY_MARKER = "__SILEX_BENCH_SILEX_READY__"
TARGET_DONE_MARKER = "__SILEX_BENCH_SILEX_TARGET_DONE__"
TARGET_NONCE = "0123456789abcdef0123456789abcdef"


def run_instance(
    executable: Path,
    coefficients: str,
    operation: str,
    *,
    warmup_coefficients: str | None = None,
    expected_returncode: int = 0,
) -> dict[str, Any]:
    command = [
        str(executable),
        "--coeffs",
        coefficients,
        "--operation",
        operation,
    ]
    if warmup_coefficients is not None:
        command.extend(["--warmup-coeffs", warmup_coefficients])
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=20.0,
    )
    assert completed.returncode == expected_returncode, (
        completed.returncode,
        completed.stdout,
        completed.stderr,
    )
    assert completed.stderr == ""
    return json.loads(completed.stdout)


def run_marked_instance(
    executable: Path,
    coefficients: str,
    operation: str,
    *,
    warmup_coefficients: str | None = None,
) -> dict[str, Any]:
    command = [
        str(executable),
        "--coeffs",
        coefficients,
        "--operation",
        operation,
        "--marked-protocol",
    ]
    if warmup_coefficients is not None:
        command.extend(["--warmup-coeffs", warmup_coefficients])
    process = subprocess.Popen(
        command,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdin is not None
    assert process.stdout is not None
    process.stdin.write("run-warmup\n")
    process.stdin.flush()
    assert process.stdout.readline().strip() == READY_MARKER
    process.stdin.write(TARGET_NONCE + "\n")
    process.stdin.flush()
    assert process.stdout.readline().strip() == f"{TARGET_DONE_MARKER}:{TARGET_NONCE}"
    stdout, stderr = process.communicate("publish-json\n", timeout=20.0)
    assert process.returncode == 0, (process.returncode, stdout, stderr)
    assert stderr == ""
    return json.loads(stdout)


def assert_invalid_marked_nonce(executable: Path) -> None:
    completed = subprocess.run(
        [
            str(executable),
            "--coeffs=-5,0,1",
            "--operation=maximal_order",
            "--marked-protocol",
        ],
        input=f"run-warmup\n{'A' * 32}\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=20.0,
    )
    assert completed.returncode == 5
    lines = completed.stdout.splitlines()
    assert lines[0] == READY_MARKER
    assert TARGET_DONE_MARKER not in completed.stdout
    payload = json.loads("\n".join(lines[1:]))
    assert payload["error"] == "marked protocol target nonce is invalid"


def check_common(payload: dict[str, Any], operation: str) -> None:
    assert payload["engine"] == "silex"
    assert payload["operation"] == operation
    assert payload["success"] is True
    assert payload["engine_thread_count"] == 1
    assert payload["error"] is None
    assert payload["source"] == "silex_public_api"
    assert payload["timing_scope"] == SCOPES[operation]
    assert payload["timing_clock"] == {
        "cpu": "std_clock_process_cpu",
        "wall": "steady_clock",
    }
    assert payload["target_cpu_ms"] is not None
    assert payload["target_cpu_ms"] >= 0.0
    assert payload["target_wall_ms"] is not None
    assert payload["target_wall_ms"] >= 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    args = parser.parse_args()

    fields = [
        ("-5,0,1", ["-5", "0", "1"], "5"),
        ("47,0,1", ["47", "0", "1"], "-47"),
    ]
    for coefficients, normalized_coefficients, discriminant in fields:
        maximal = run_instance(args.exe, coefficients, "maximal_order")
        check_common(maximal, "maximal_order")
        assert maximal["coefficients_low_to_high"] == normalized_coefficients
        assert maximal["maximal_order_discriminant"] == discriminant
        assert maximal["ideal_norm"] is None
        assert maximal["root_found"] is None
        assert maximal["root_verified"] is None
        assert maximal["warmup"] == {"used": False, "degree": None}

        ideal = run_instance(args.exe, coefficients, "ideal_multiply")
        check_common(ideal, "ideal_multiply")
        assert ideal["maximal_order_discriminant"] is None
        assert ideal["ideal_norm"] == "36"
        assert ideal["root_found"] is None
        assert ideal["root_verified"] is None

        square_root = run_instance(
            args.exe, coefficients, "element_square_root"
        )
        check_common(square_root, "element_square_root")
        assert square_root["maximal_order_discriminant"] is None
        assert square_root["ideal_norm"] is None
        assert square_root["root_found"] is True
        assert square_root["root_verified"] is True

    warmed = run_instance(
        args.exe,
        "-5,0,1",
        "ideal_multiply",
        warmup_coefficients="47,0,1",
    )
    check_common(warmed, "ideal_multiply")
    assert warmed["warmup"] == {"used": True, "degree": 2}

    marked = run_marked_instance(
        args.exe,
        "-5,0,1",
        "ideal_multiply",
        warmup_coefficients="47,0,1",
    )
    check_common(marked, "ideal_multiply")
    assert marked["ideal_norm"] == "36"
    assert marked["warmup"] == {"used": True, "degree": 2}
    assert_invalid_marked_nonce(args.exe)

    repeated_warmup = run_instance(
        args.exe,
        "-5,0,1",
        "maximal_order",
        warmup_coefficients="-5,0,1",
        expected_returncode=2,
    )
    assert repeated_warmup["success"] is False
    assert repeated_warmup["error"] == (
        "--warmup-coeffs must define a distinct field"
    )
    assert repeated_warmup["target_cpu_ms"] is None
    assert repeated_warmup["target_wall_ms"] is None

    wrong_degree = run_instance(
        args.exe,
        "-5,0,1",
        "maximal_order",
        warmup_coefficients="-2,0,0,1",
        expected_returncode=2,
    )
    assert wrong_degree["success"] is False
    assert wrong_degree["error"] == (
        "--warmup-coeffs must have the same degree as --coeffs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
