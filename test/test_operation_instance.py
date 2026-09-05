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
    "maximal_order": "maximal_order_only",
    "ideal_multiply": "ideal_multiplication_only",
    "element_square_root": "number_field_element_is_square_only",
}
READY_MARKER = "__SILEX_BENCH_SILEX_READY__"
TARGET_DONE_MARKER = "__SILEX_BENCH_SILEX_TARGET_DONE__"
TARGET_NONCE = "0123456789abcdef0123456789abcdef"


def run_instance(
    executable: Path,
    coefficients: str,
    operation: str,
    *,
    expected_returncode: int = 0,
) -> dict[str, Any]:
    command = [
        str(executable),
        "--coeffs",
        coefficients,
        "--operation",
        operation,
    ]
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
) -> dict[str, Any]:
    command = [
        str(executable),
        "--coeffs",
        coefficients,
        "--operation",
        operation,
        "--marked-protocol",
    ]
    process = subprocess.Popen(
        command,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdin is not None
    assert process.stdout is not None
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
        input=f"{'A' * 32}\n",
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


def assert_marked_phase_failures(executable: Path) -> None:
    for operation in SCOPES:
        for supplied, done in (("", False), ("a" * 33 + "\n", False),
                               (TARGET_NONCE + "\n", True)):
            completed = subprocess.run(
                [str(executable), "--coeffs=-5,0,1", f"--operation={operation}",
                 "--marked-protocol"], input=supplied, text=True,
                capture_output=True, timeout=20.0,
            )
            assert completed.returncode == 5
            lines = completed.stdout.splitlines()
            assert lines.pop(0) == READY_MARKER
            if done:
                assert lines.pop(0) == f"{TARGET_DONE_MARKER}:{TARGET_NONCE}"
            else:
                assert TARGET_DONE_MARKER not in completed.stdout
            payload = json.loads("\n".join(lines))
            assert payload["success"] is False
            assert payload["error"] == (
                "marked protocol missing final phase input" if done
                else "marked protocol target nonce is invalid")
            assert (payload["target_cpu_ms"] is not None) == done


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    args = parser.parse_args()
    assert_marked_phase_failures(args.exe)

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
        assert "warmup" not in maximal

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

    # The maximal-order marked protocol reaches READY only after field and
    # equation-order preparation, then brackets the maximal-order call itself.
    marked_maximal = run_marked_instance(
        args.exe,
        "-5,0,1",
        "maximal_order",
    )
    check_common(marked_maximal, "maximal_order")
    assert marked_maximal["coefficients_low_to_high"] == ["-5", "0", "1"]
    assert marked_maximal["maximal_order_discriminant"] == "5"
    assert marked_maximal["ideal_norm"] is None
    assert marked_maximal["root_found"] is None
    assert marked_maximal["root_verified"] is None
    assert "warmup" not in marked_maximal

    marked = run_marked_instance(
        args.exe,
        "-5,0,1",
        "ideal_multiply",
    )
    check_common(marked, "ideal_multiply")
    assert marked["ideal_norm"] == "36"
    assert "warmup" not in marked

    marked_square_root = run_marked_instance(
        args.exe,
        "-5,0,1",
        "element_square_root",
    )
    check_common(marked_square_root, "element_square_root")
    assert marked_square_root["root_found"] is True
    assert marked_square_root["root_verified"] is True
    assert_invalid_marked_nonce(args.exe)

    removed_warmup = subprocess.run(
        [
            str(args.exe),
            "--coeffs=-5,0,1",
            "--operation=maximal_order",
            "--warmup-coeffs=47,0,1",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=20.0,
    )
    assert removed_warmup.returncode == 2
    removed_payload = json.loads(removed_warmup.stdout)
    assert removed_payload["success"] is False
    assert removed_payload["error"] == "unknown argument: --warmup-coeffs=47,0,1"
    assert removed_payload["target_cpu_ms"] is None
    assert removed_payload["target_wall_ms"] is None
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
