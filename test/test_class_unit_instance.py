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
    assert len(manifest["fields"]) == 23
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
        # The instance JSON always reports the GRH dependence of the
        # analytic hR check (ClassGroup::analytic_class_regulator_
        # certification(), decision 2026-09-26 "T-012 public gates and
        # GRH accessor"), whether or not that check was exercised.
        assert "analytic_class_regulator_certification" in (
            proven_instance["class_group"]
        )
        if proven_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "not_checked":
            assert proven_instance["class_group"][
                "analytic_class_regulator_certification"
            ] == "unknown"
        if row["id"] in exact_rows:
            # Exact completion needs no analytic or BF receipt.
            for component in (
                "analytic_class_regulator_status", "zeta_bf_proof_status"
            ):
                assert proven_instance["class_group"][component] == "not_checked"
        if row["id"] == "real_quadratic_210_proven":
            # The canonical Dirichlet index-one gate certifies this pair
            # without a relation-saturation or BF proof attempt. The gate
            # is unconditional, so the GRH dependence field reads
            # `proven`, not the result's overall `proven` certification.
            assert proven_instance["class_group"][
                "analytic_class_regulator_status"
            ] == "verified"
            assert proven_instance["class_group"][
                "analytic_class_regulator_certification"
            ] == "proven"
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
        # A `grh` request accepts the GRH-conditional index-one bound in
        # `try_validate_refine_loop` without ever calling the certify helper
        # that records `analytic_class_regulator_status`/`_certification`
        # (src/order_unit/validation.cpp); this is a known library gap (see
        # notes/reviews/T-044-reviewer.md, "GRH-mode gap"), not something
        # this task changes. Every `grh_certification` row therefore leaves
        # the analytic check itself unrecorded here; its `grh` label comes
        # from the GRH-dependent factor-base bound instead, not from
        # `analytic_class_regulator_status`/`_certification`. Pinning this
        # as `not_checked`/`unknown` means a future fix to that gap shows up
        # as a test change here.
        assert grh_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "not_checked"
        assert grh_instance["class_group"][
            "analytic_class_regulator_certification"
        ] == "unknown"

    # An actually GRH-conditional (`"grh"`-valued) analytic hR check is
    # exercised by a dedicated fixture row instead: `--zeta-bf-audit`
    # (tools/class_unit_instance.cpp) runs ClassGroupContext::
    # try_certify_class_unit_with_zeta_bf(...) post-hoc, mirroring
    # test/t-order-unit.cpp's
    # test_belabas_friedman_class_regulator_is_grh_conditional (a `proven`
    # request settles a nontrivial-degree field by saturation, then a
    # Belabas-Friedman audit recorded afterward is labelled `grh` without
    # replacing the saturation proof or changing the overall `proven`
    # result).
    analytic_grh_rows = [
        row
        for row in manifest["fields"]
        if row.get("status") == "analytic_grh_certification"
    ]
    assert len(analytic_grh_rows) == 1
    for row in analytic_grh_rows:
        assert row["mode"] == "proven"
        coeffs = ",".join(
            str(value) for value in row["coefficients_low_to_high"]
        )
        analytic_instance = run_json(
            [
                str(args.exe),
                "--coeffs",
                coeffs,
                "--mode",
                row["mode"],
                "--zeta-bf-audit",
            ],
            root,
        )
        assert analytic_instance["success"] is True
        assert analytic_instance["final_result_published"] is True
        # The overall result is unaffected: still `proven`, via saturation.
        assert analytic_instance["certification_status"] == "proven"
        assert analytic_instance["class_group_proof_status"] == "proven"
        assert analytic_instance["unit_group_proof_status"] == "proven"
        assert analytic_instance["class_group"]["order"] == str(
            row["expected_class_order"]
        )
        assert analytic_instance["unit_group"]["free_rank"] == row[
            "expected_unit_rank"
        ]
        assert analytic_instance["class_group"][
            "relation_saturation_status"
        ] == "verified"
        # The post-hoc audit is what newly reports the analytic check's own
        # GRH-conditional certification.
        assert analytic_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "verified"
        assert analytic_instance["class_group"][
            "analytic_class_regulator_certification"
        ] == "grh"
        assert analytic_instance["class_group"][
            "zeta_bf_proof_status"
        ] == "verified"
        # The audit's outcome is visible on its own: it ran and succeeded,
        # distinct from "not requested" or "requested but skipped/failed".
        assert analytic_instance["zeta_bf_audit"] == {
            "requested": True,
            "ran": True,
            "succeeded": True,
            "skip_reason": None,
            "wall_ms": analytic_instance["zeta_bf_audit"]["wall_ms"],
        }
        assert analytic_instance["zeta_bf_audit"]["wall_ms"] >= 0.0

        # Without `--zeta-bf-audit` the same field settles by saturation
        # alone and never exercises the analytic check (regression check
        # for the flag's default-off behavior).
        unaudited_instance = run_json(
            [
                str(args.exe),
                "--coeffs",
                coeffs,
                "--mode",
                row["mode"],
            ],
            root,
        )
        assert unaudited_instance["success"] is True
        assert unaudited_instance["certification_status"] == "proven"
        assert unaudited_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "not_checked"
        assert unaudited_instance["class_group"][
            "analytic_class_regulator_certification"
        ] == "unknown"
        assert unaudited_instance["zeta_bf_audit"] == {
            "requested": False,
            "ran": False,
            "succeeded": None,
            "skip_reason": None,
            "wall_ms": None,
        }

        # `--zeta-bf-audit` must not change published certification labels
        # when the transaction itself did not publish `proven` for both the
        # class group and the units. `try_certify_class_unit_with_zeta_bf`'s
        # own hR is unconditional in degree one (and on the quadratic
        # L(1, chi) route), so it would otherwise publish `proven` outright
        # even for a `grh`-requested run, mismatching the top-level
        # `certification_status` (already fixed by the transaction, before
        # any audit runs). Checked here on this row's own higher-degree
        # field under `--mode grh` (regression for T-044 round-2 blocker
        # B1); the degree-one case is checked separately below.
        grh_audited_instance = run_json(
            [
                str(args.exe),
                "--coeffs",
                coeffs,
                "--mode",
                "grh",
                "--zeta-bf-audit",
            ],
            root,
        )
        assert grh_audited_instance["success"] is True
        assert grh_audited_instance["certification_status"] == "grh"
        assert grh_audited_instance["class_group"]["certification"] == "grh"
        assert grh_audited_instance["unit_group"]["certification"] == "grh"
        assert grh_audited_instance["zeta_bf_audit"] == {
            "requested": True,
            "ran": False,
            "succeeded": None,
            "skip_reason": "transaction_certification_not_proven",
            "wall_ms": None,
        }

    degree_one_rows = [
        row
        for row in manifest["fields"]
        if row["id"] == "degree_one_proven"
    ]
    assert len(degree_one_rows) == 1
    for row in degree_one_rows:
        coeffs = ",".join(
            str(value) for value in row["coefficients_low_to_high"]
        )
        # Degree one is exactly the case where the Belabas-Friedman hR is
        # unconditional, so `try_certify_class_unit_with_zeta_bf` publishes
        # `proven` outright rather than only recording a check. Gating the
        # audit on the transaction having already published `proven` for
        # both objects (B1) means a `grh`-requested degree-one run keeps its
        # `grh` labels even with `--zeta-bf-audit`.
        degree_one_grh_audited = run_json(
            [
                str(args.exe),
                "--coeffs",
                coeffs,
                "--mode",
                "grh",
                "--zeta-bf-audit",
            ],
            root,
        )
        assert degree_one_grh_audited["success"] is True
        assert degree_one_grh_audited["degree"] == 1
        assert degree_one_grh_audited["certification_status"] == "grh"
        assert degree_one_grh_audited["class_group"]["certification"] == (
            "grh"
        )
        assert degree_one_grh_audited["unit_group"]["certification"] == (
            "grh"
        )
        assert degree_one_grh_audited["zeta_bf_audit"] == {
            "requested": True,
            "ran": False,
            "succeeded": None,
            "skip_reason": "transaction_certification_not_proven",
            "wall_ms": None,
        }

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
