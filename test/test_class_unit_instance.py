#!/usr/bin/env python3
"""Replay the fixed class/unit inventory and marked adapter protocol."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )

READY_MARKER = "__SILEX_BENCH_SILEX_READY__"
TARGET_DONE_MARKER = "__SILEX_BENCH_SILEX_TARGET_DONE__"
TARGET_NONCE = "0123456789abcdef0123456789abcdef"


def run_json(cmd: list[str], root: Path) -> dict[str, object]:
    completed = subprocess.run(
        cmd,
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=30.0,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed ({completed.returncode}): {completed.stderr}\n"
            f"{completed.stdout}"
        )
    return json.loads(completed.stdout)


def run_marked_json(cmd: list[str], root: Path) -> dict[str, object]:
    process = subprocess.Popen(
        [*cmd, "--marked-protocol"],
        cwd=root,
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
    stdout, stderr = process.communicate("publish-json\n", timeout=30.0)
    if process.returncode != 0:
        raise AssertionError(
            f"marked command failed ({process.returncode}): {stderr}\n{stdout}"
        )
    return json.loads(stdout)


def assert_invalid_marked_nonce(cmd: list[str], root: Path) -> None:
    completed = subprocess.run(
        [*cmd, "--marked-protocol"],
        cwd=root,
        input=f"{'A' * 32}\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=30.0,
    )
    assert completed.returncode == 5
    lines = completed.stdout.splitlines()
    assert lines[0] == READY_MARKER
    assert TARGET_DONE_MARKER not in completed.stdout
    payload = json.loads("\n".join(lines[1:]))
    assert payload["error"] == "marked protocol target nonce is invalid"


def assert_rejected(cmd: list[str], root: Path) -> None:
    completed = subprocess.run(
        cmd,
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=30.0,
    )
    assert completed.returncode != 0


def assert_marked_phase_failures(exe: Path, root: Path) -> None:
    for supplied, done in (("", False), ("a" * 33 + "\n", False),
                           (TARGET_NONCE + "\n", True)):
        completed = subprocess.run(
            [str(exe), "--coeffs=-5,0,1", "--mode=proven", "--marked-protocol"],
            input=supplied, text=True, capture_output=True, cwd=root, timeout=30.0,
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    assert_marked_phase_failures(args.exe, root)
    instance_script = root / "tools/bench/run-class-unit-instance.py"
    manifest = json.loads(args.manifest.read_text())

    proven_rows = [
        row
        for row in manifest["fields"]
        if row.get("status") == "must_pass_fast"
    ]
    assert len(manifest["fields"]) == 22
    assert len(proven_rows) == 17
    exact_rows = {
        "degree_one_proven",
        "imaginary_quadratic_disc3_proven",
        "imaginary_quadratic_disc3_alternate_proven",
        "imaginary_quadratic_disc4_proven",
        "imaginary_quadratic_47_proven",
    }
    proven_instances: dict[str, dict[str, object]] = {}
    for row in proven_rows:
        assert row["timeout_seconds"] == 20
        proven_instance = run_json(
            [
                sys.executable,
                str(instance_script),
                "--exe",
                str(args.exe),
                "--manifest",
                str(args.manifest),
                "--field-id",
                row["id"],
            ],
            root,
        )
        assert proven_instance["success"] is True
        assert proven_instance["final_result_published"] is True
        assert proven_instance["field_defined"] is True
        assert proven_instance["maximal_order_defined"] is True
        assert proven_instance["class_group"]["has_presentation"] is True
        assert proven_instance["unit_group"]["is_set"] is True
        assert proven_instance["certification_status"] == "proven"
        assert proven_instance["class_group_proof_status"] == "proven"
        assert proven_instance["unit_group_proof_status"] == "proven"
        assert proven_instance["regulator_proof_status"] == "verified"
        assert proven_instance["class_group"][
            "factor_base_generation_status"
        ] == "verified"
        assert proven_instance["class_group"][
            "relation_saturation_status"
        ] == ("not_checked" if row["id"] == "real_quadratic_210_proven"
              else "verified")
        assert proven_instance["class_group"]["unit_proof_status"] == (
            "verified"
        )
        assert proven_instance["class_group"]["regulator_proof_status"] == (
            "verified"
        )
        assert proven_instance["class_group"]["order"] == str(
            row["expected_class_order"]
        )
        assert proven_instance["unit_group"]["free_rank"] == row[
            "expected_unit_rank"
        ]
        assert proven_instance["class_group"]["certification"] == "proven"
        assert proven_instance["unit_group"]["certification"] == "proven"
        r1, r2 = proven_instance["signature"]
        assert r1 + 2 * r2 == len(row["coefficients_low_to_high"]) - 1
        assert r1 + r2 - 1 == row["expected_unit_rank"]
        if "expected_signature" in row:
            assert proven_instance["signature"] == row["expected_signature"]
            assert proven_instance["maximal_order_discriminant"] == row[
                "maximal_order_discriminant"
            ]
        if row["expected_class_order"] == 1:
            assert proven_instance["class_group"]["invariants"] == []
        if "expected_class_invariants" in row:
            assert proven_instance["class_group"]["invariants"] == [
                str(value) for value in row["expected_class_invariants"]
            ]
        if row["id"] in exact_rows:
            # Exact completion needs no analytic or BF receipt.
            for component in (
                "analytic_class_regulator_status", "zeta_bf_proof_status"
            ):
                assert proven_instance["class_group"][component] == "not_checked"
        if row["id"] == "real_quadratic_210_proven":
            # The canonical Dirichlet index-one gate certifies this pair
            # without a relation-saturation or BF proof attempt.
            assert proven_instance["class_group"][
                "analytic_class_regulator_status"
            ] == "verified"
            assert proven_instance["class_group"][
                "zeta_bf_proof_status"
            ] == "not_checked"
        proven_instances[row["id"]] = proven_instance

    default_instance = proven_instances["real_quadratic_5_proven"]
    assert "algorithm" not in default_instance
    assert "fallback_used" not in default_instance
    assert "policy" not in default_instance

    grh_rows = [
        row
        for row in manifest["fields"]
        if row.get("status") == "grh_certification"
    ]
    assert len(grh_rows) == 5
    for row in grh_rows:
        grh_instance = run_json(
            [
                sys.executable,
                str(instance_script),
                "--exe",
                str(args.exe),
                "--manifest",
                str(args.manifest),
                "--field-id",
                row["id"],
                "--timeout",
                "20",
            ],
            root,
        )
        assert grh_instance["success"] is True
        assert grh_instance["final_result_published"] is True
        assert grh_instance["certification_status"] == "grh"
        assert grh_instance["class_group_proof_status"] == "grh"
        assert grh_instance["unit_group_proof_status"] == "grh"
        assert grh_instance["regulator_proof_status"] == "not_checked"
        assert grh_instance["class_group"]["order"] == str(
            row["expected_class_order"]
        )
        assert grh_instance["class_group"]["invariants"] == [
            str(value) for value in row["expected_class_invariants"]
        ]
        assert grh_instance["unit_group"]["free_rank"] == 0

    for removed_option in (
        "--coordinate-radius=2",
        "--ideal-radius=1",
        "--target-kernel-units=1",
        "--post-finite-budget=1",
    ):
        assert_rejected(
            [
                str(args.exe),
                "--coeffs=-5,0,1",
                "--mode=proven",
                removed_option,
            ],
            root,
        )

    marked_instance = run_marked_json(
        [
            str(args.exe),
            "--coeffs=-5,0,1",
            "--mode=proven",
        ],
        root,
    )
    assert marked_instance["success"] is True
    assert marked_instance["engine_thread_count"] == 1
    assert "warmup" not in marked_instance
    assert marked_instance["measurement_timing"] == {
        "algorithm_clock": "std_clock_process_cpu",
        "algorithm_scope": "class_unit_transaction_only",
        "component_clock": "steady_clock",
        "preparation_scope": (
            "field_construction+equation_order+maximal_order+compute_options"
        ),
        "preparation_excluded": True,
        "finalization_excluded": True,
        "target_cpu_ms": marked_instance["measurement_timing"][
            "target_cpu_ms"
        ],
        "target_wall_ms": marked_instance["measurement_timing"][
            "target_wall_ms"
        ],
    }
    assert_invalid_marked_nonce(
        [
            str(args.exe),
            "--coeffs=-5,0,1",
            "--mode=proven",
        ],
        root,
    )

    assert_rejected(
        [
            str(args.exe),
            "--coeffs=-5,0,1",
            "--mode=proven",
            "--warmup-coeffs=2,2,1",
        ],
        root,
    )

    instance = run_json(
        [
            sys.executable,
            str(instance_script),
            "--exe",
            str(args.exe),
            "--manifest",
            str(args.manifest),
            "--field-id",
            "real_quadratic_5_proven",
            "--timeout",
            "20",
        ],
        root,
    )
    assert instance["success"] is True
    assert "warmup" not in instance
    assert instance["component_timing_ms"]["total"] > 0.0
    assert instance["measurement_timing"]["algorithm_clock"] == (
        "std_clock_process_cpu"
    )
    assert instance["measurement_timing"]["target_cpu_ms"] > 0.0
    assert instance["measurement_timing"]["target_wall_ms"] > 0.0
    assert instance["signature"] == [2, 0]
    assert instance["maximal_order_discriminant"] == "5"

    assert_rejected(
        [
            sys.executable,
            str(instance_script),
            "--exe",
            str(args.exe),
            "--manifest",
            str(args.manifest),
            "--field-id",
            "real_quadratic_5_proven",
            "--warmup-coeffs=2,2,1",
        ],
        root,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
