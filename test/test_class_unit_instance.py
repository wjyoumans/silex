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


def assert_unit_regulator_enclosure(instance, expected_proof_status):
    """The unit_group object reports the regulator enclosure and its label."""
    unit_group = instance["unit_group"]
    midpoint = unit_group["regulator_midpoint"]
    radius = unit_group["regulator_radius"]
    assert isinstance(midpoint, (int, float))
    assert midpoint > 0.0
    assert 0.0 <= radius < 1e-20 * max(1.0, midpoint)
    assert abs(float(unit_group["regulator_decimal"]) - midpoint) <= (
        1e-12 * midpoint
    )
    assert unit_group["regulator_proof_status"] == expected_proof_status
    assert unit_group["regulator_proof_status"] == instance[
        "regulator_proof_status"
    ]


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


def assert_fail_closed_inputs(exe: Path, root: Path) -> None:
    # A reducible defining polynomial is rejected at field construction
    # (NumberField::by_polynomial requires irreducibility over Q), so no
    # maximal order or class/unit transaction exists in either mode.
    # x^3 + x + 30 = (x + 3)(x^2 - 3x + 10).
    for mode in ("proven", "grh"):
        instance = run_json(
            [str(exe), "--coeffs", "30,1,0,1", "--mode", mode], root
        )
        assert instance["success"] is False
        assert instance["final_result_published"] is False
        assert instance["field_defined"] is False
        assert instance["maximal_order_defined"] is False
        assert instance["failure_reason"] == "input_or_options_unavailable"


def assert_grh_generation_below_minkowski(exe: Path, root: Path) -> None:
    # A grh request of positive unit rank may take factor-base generation
    # from GRH: its base contains every prime ideal of norm at most the GRH
    # bound min(BDF, Bach), so the analytic index-one test may accept below
    # the proven generation bound (factor_base_class_group_bound). These fields
    # used to fail closed for want of Minkowski-type coverage; they now
    # publish exactly the grh label. GRH coverage is never reported as
    # unconditional evidence: generation stays `unavailable` and relation
    # saturation is not `verified`. x^3 + x + 200: |D| = 1080004, bound
    # 100 < 295, class group Z/2. x^2 - 100003: D = 400012, bound
    # 100 < 316, h = 1.
    for coeffs, used, requested, order, invariants in (
        ("200,1,0,1", "100", "295", "2", ["2"]),
        ("-100003,0,1", "100", "316", "1", []),
    ):
        instance = run_json(
            [str(exe), "--coeffs", coeffs, "--mode", "grh"], root
        )
        assert instance["success"] is True
        assert instance["final_result_published"] is True
        assert instance["certification_status"] == "grh"
        assert instance["class_group_proof_status"] == "grh"
        assert instance["unit_group_proof_status"] == "grh"
        assert instance["factor_base_bound"] == used
        assert instance["requested_factor_base_bound"] == requested
        class_group = instance["class_group"]
        assert class_group["order"] == order
        assert class_group["invariants"] == invariants
        assert class_group["certification"] == "grh"
        assert instance["unit_group"]["certification"] == "grh"
        assert class_group["factor_base_generation_status"] == "unavailable"
        assert class_group["relation_saturation_status"] != "verified"
        # The real-quadratic hR comes from an unconditional L(1, chi)
        # evaluation, the cubic one from a Belabas-Friedman bound; the
        # published label is grh either way, since generation rests on GRH.
        assert class_group["analytic_class_regulator_status"] == "verified"
        assert class_group["analytic_class_regulator_certification"] == (
            "proven" if instance["degree"] == 2 else "grh"
        )
        # Provenance: generation rests on the BDF bound (50 for both
        # fields, below Bach's), completeness on the hR that accepted.
        assert class_group["factor_base_generation_certification"] == "grh"
        assert class_group["factor_base_generation_basis"] == "bdf"
        assert class_group[
            "factor_base_generation_certification_bound"
        ] == "50"
        if instance["degree"] == 2:
            assert class_group[
                "class_unit_completeness_certification"
            ] == "proven"
            assert class_group["class_unit_completeness_basis"] == (
                "unconditional_analytic"
            )
        else:
            assert class_group[
                "class_unit_completeness_certification"
            ] == "grh"
            assert class_group["class_unit_completeness_basis"] == (
                "belabas_friedman"
            )


def assert_proven_and_degree_one_provenance(exe: Path, root: Path) -> None:
    # A proven result: generation to the Minkowski-type bound and
    # unconditional completeness. The degree-one grh route is the exact
    # route relabelled grh: both components proven, label grh.
    for coeffs, mode, bound in (
        ("200,1,0,1", "proven", "295"),
        ("0,1", "grh", "1"),
    ):
        instance = run_json(
            [str(exe), "--coeffs", coeffs, "--mode", mode], root
        )
        assert instance["success"] is True
        assert instance["certification_status"] == mode
        class_group = instance["class_group"]
        assert class_group["factor_base_generation_certification"] == (
            "proven"
        )
        assert class_group["factor_base_generation_basis"] == (
            "minkowski_type"
        )
        assert class_group[
            "factor_base_generation_certification_bound"
        ] == bound
        assert class_group["class_unit_completeness_certification"] == (
            "proven"
        )
        assert class_group["class_unit_completeness_basis"] == (
            "unconditional_certification"
        )


def assert_grh_minkowski_equality_boundary(exe: Path, root: Path) -> None:
    # The equality boundary of unconditional generation in grh mode:
    # record_factor_base_generation_ (src/class_group/class_group.cpp)
    # verifies with `>=`, so a GRH-sized factor base that reaches the
    # Minkowski-type bound exactly (not just strictly above it) records
    # unconditional (`verified`) generation, while the published label
    # stays grh. x^2 - 40001: 40001 = 13 * 17 * 181 is
    # squarefree and 1 mod 4, so D = 40001 and the real-quadratic branch of
    # factor_base_class_group_bound gives floor(sqrt(40001) / 2) = 100,
    # exactly the grh policy's selected bound here. h = 32, class group
    # Z/16 x Z/2 (GP 2.17 quadclassunit(40001)).
    instance = run_json(
        [str(exe), "--coeffs", "-40001,0,1", "--mode", "grh"], root
    )
    assert instance["success"] is True
    assert instance["final_result_published"] is True
    assert instance["certification_status"] == "grh"
    assert instance["factor_base_bound"] == "100"
    assert instance["requested_factor_base_bound"] == "100"
    assert instance["class_group"]["factor_base_generation_status"] == (
        "verified"
    )
    assert instance["class_group"]["order"] == "32"
    assert instance["class_group"]["invariants"] == ["2", "16"]
    # Generation verified to the Minkowski-type bound is reported proven,
    # which tells this result apart from one whose generation rests on GRH
    # although both carry grh labels.
    class_group = instance["class_group"]
    assert class_group["factor_base_generation_certification"] == "proven"
    assert class_group["factor_base_generation_basis"] == "minkowski_type"
    assert class_group["factor_base_generation_certification_bound"] == "100"
    assert class_group["class_unit_completeness_certification"] == "proven"
    assert class_group["class_unit_completeness_basis"] == (
        "unconditional_analytic"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    assert_marked_phase_failures(args.exe, root)
    assert_fail_closed_inputs(args.exe, root)
    assert_grh_generation_below_minkowski(args.exe, root)
    assert_grh_minkowski_equality_boundary(args.exe, root)
    assert_proven_and_degree_one_provenance(args.exe, root)
    instance_script = root / "tools/bench/run-class-unit-instance.py"
    manifest = json.loads(args.manifest.read_text())

    proven_rows = [
        row
        for row in manifest["fields"]
        if row.get("status") == "must_pass_fast"
    ]
    assert len(manifest["fields"]) == 29
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
        ] == "verified"
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
        assert_unit_regulator_enclosure(proven_instance, "verified")
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
        # certification()), whether or not that check was exercised.
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
            # without an ell-local saturation or BF proof attempt; index
            # one proves saturation at every prime, so the saturation
            # status above is `verified`. The gate is unconditional, so
            # the GRH dependence field reads
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
    assert len(grh_rows) == 11
    grh_bounds: dict[str, tuple[int, int]] = {}
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
        assert grh_instance["unit_group"]["free_rank"] == row[
            "expected_unit_rank"
        ]
        assert grh_instance["class_group"]["certification"] == "grh"
        assert grh_instance["unit_group"]["certification"] == "grh"
        assert_unit_regulator_enclosure(
            grh_instance, grh_instance["regulator_proof_status"]
        )
        assert grh_instance["unit_group"]["regulator_proof_status"] == (
            "not_checked"
        )
        # `factor_base_bound` is the bound the transaction used (the GRH
        # policy's selected bound), `requested_factor_base_bound` the
        # tool-side request (the Minkowski-type bound).
        used = int(grh_instance["factor_base_bound"])
        minkowski = int(grh_instance["requested_factor_base_bound"])
        assert used >= 2
        assert minkowski >= 2
        grh_bounds[row["id"]] = (used, minkowski)
        provenance = grh_instance["class_group"]
        generation_bound = int(
            provenance["factor_base_generation_certification_bound"]
        )
        if used >= minkowski:
            assert provenance["factor_base_generation_certification"] == (
                "proven"
            )
            assert provenance["factor_base_generation_basis"] == (
                "minkowski_type"
            )
            assert generation_bound == minkowski
        else:
            assert provenance["factor_base_generation_certification"] == (
                "grh"
            )
            assert provenance["factor_base_generation_basis"] in (
                "bdf",
                "bach",
            )
            assert generation_bound <= used
        if row["expected_unit_rank"] == 0:
            # Rank-zero grh rows are imaginary quadratic and take the exact
            # imaginary-quadratic `grh` route: the index comes from the
            # exact class number, not from an analytic hR, and its GRH
            # dependence is in factor-base generation (checked only up to
            # the GRH bound). That route records no analytic check, so
            # these rows stay `not_checked`/`unknown`.
            assert grh_instance["class_group"][
                "analytic_class_regulator_status"
            ] == "not_checked"
            assert grh_instance["class_group"][
                "analytic_class_regulator_certification"
            ] == "unknown"
            assert provenance["class_unit_completeness_certification"] == (
                "proven"
            )
            assert provenance["class_unit_completeness_basis"] == (
                "exact_class_number"
            )
            continue
        # Positive unit rank: the analytic index-one test accepts, with a
        # Belabas-Friedman hR (GRH-conditional) above degree two and an
        # unconditional L(1, chi) hR for real quadratic fields. Generation is
        # unconditional (`verified`) only when the GRH-sized base reaches
        # the Minkowski-type bound; below it, generation rests on GRH and
        # is never reported `verified`. Relation saturation is not proved
        # by a GRH-conditional acceptance.
        assert grh_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "verified"
        assert grh_instance["class_group"][
            "analytic_class_regulator_certification"
        ] == ("proven" if grh_instance["degree"] == 2 else "grh")
        assert grh_instance["class_group"][
            "factor_base_generation_status"
        ] == ("verified" if used >= minkowski else "unavailable")
        assert grh_instance["class_group"][
            "relation_saturation_status"
        ] != "verified"
        if grh_instance["degree"] == 2:
            assert provenance["class_unit_completeness_certification"] == (
                "proven"
            )
            assert provenance["class_unit_completeness_basis"] == (
                "unconditional_analytic"
            )
        else:
            assert provenance["class_unit_completeness_certification"] == (
                "grh"
            )
            assert provenance["class_unit_completeness_basis"] == (
                "belabas_friedman"
            )
    # Each GRH-generation row is below the Minkowski-type bound; the
    # x^2 - 10007 row sits exactly on it.
    for field_id in (
        "cubic_x3_x_200_grh",
        "real_quadratic_100003_grh",
        "quartic_disc1412343_grh",
        "quintic_disc401370255_grh",
        "imaginary_quadratic_100003_grh",
    ):
        used, minkowski = grh_bounds[field_id]
        assert used < minkowski
    used, minkowski = grh_bounds["real_quadratic_10007_grh"]
    assert used == minkowski

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
        # The row asks for the audit itself, so the bench wrapper passes
        # `--zeta-bf-audit` when the row is run by id.
        assert row["zeta_bf_audit"] is True
        coeffs = ",".join(
            str(value) for value in row["coefficients_low_to_high"]
        )
        analytic_instance = run_json(
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
        # own Belabas-Friedman hR is unconditional in degree one, so it
        # would otherwise publish `proven` outright even for a
        # `grh`-requested run, mismatching the top-level
        # `certification_status` (already fixed by the transaction, before
        # any audit runs). Checked here on this row's own higher-degree
        # field under `--mode grh`; the degree-one case is checked
        # separately below.
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
        # The `grh` transaction itself records the analytic index-one check
        # that accepted the pair: a Belabas-Friedman hR, so `grh`, together
        # with that evaluation's BF audit data (the record's contents are
        # checked in t-class-unit-matrix.cpp). The record is informational
        # and leaves the labels and the unit/regulator proof states alone.
        assert grh_audited_instance["class_group"][
            "analytic_class_regulator_status"
        ] == "verified"
        assert grh_audited_instance["class_group"][
            "analytic_class_regulator_certification"
        ] == "grh"
        assert grh_audited_instance["class_group"][
            "zeta_bf_proof_status"
        ] == "verified"
        assert grh_audited_instance["class_group_proof_status"] == "grh"
        assert grh_audited_instance["unit_group_proof_status"] == "grh"
        assert grh_audited_instance["regulator_proof_status"] == (
            "not_checked"
        )
        assert grh_audited_instance["class_group"]["unit_proof_status"] == (
            "not_checked"
        )
        # `factor_base_bound` reports the bound the transaction used, not the
        # request: the grh policy raises the requested bound on this row.
        assert int(grh_audited_instance["factor_base_bound"]) > int(
            grh_audited_instance["requested_factor_base_bound"]
        )

    # A real-quadratic `grh` transaction is accepted by the analytic
    # index-one test against the unconditional L(1, chi_D) hR, so the
    # recorded check reads `proven` while both labels stay `grh`.
    real_quadratic_grh = run_json(
        [str(args.exe), "--coeffs", "-5,0,1", "--mode", "grh"],
        root,
    )
    assert real_quadratic_grh["success"] is True
    assert real_quadratic_grh["certification_status"] == "grh"
    assert real_quadratic_grh["class_group"]["certification"] == "grh"
    assert real_quadratic_grh["unit_group"]["certification"] == "grh"
    assert real_quadratic_grh["class_group"]["order"] == "1"
    assert real_quadratic_grh["unit_group"]["free_rank"] == 1
    assert real_quadratic_grh["class_group"][
        "analytic_class_regulator_status"
    ] == "verified"
    assert real_quadratic_grh["class_group"][
        "analytic_class_regulator_certification"
    ] == "proven"
    assert real_quadratic_grh["class_group"][
        "zeta_bf_proof_status"
    ] == "not_checked"

    # The same with a nontrivial class group: Q(sqrt(10)), h = 2, class
    # group Z/2, unit rank 1 (GP 2.17 quadclassunit(40)). The regulator is
    # checked against log(3 + sqrt(10)) in t-class-unit-matrix.cpp.
    real_quadratic_h2_grh = run_json(
        [str(args.exe), "--coeffs", "-10,0,1", "--mode", "grh"],
        root,
    )
    assert real_quadratic_h2_grh["success"] is True
    assert real_quadratic_h2_grh["certification_status"] == "grh"
    assert real_quadratic_h2_grh["class_group"]["certification"] == "grh"
    assert real_quadratic_h2_grh["unit_group"]["certification"] == "grh"
    assert real_quadratic_h2_grh["class_group"]["order"] == "2"
    assert real_quadratic_h2_grh["class_group"]["invariants"] == ["2"]
    assert real_quadratic_h2_grh["unit_group"]["free_rank"] == 1
    assert real_quadratic_h2_grh["class_group"][
        "analytic_class_regulator_status"
    ] == "verified"
    assert real_quadratic_h2_grh["class_group"][
        "analytic_class_regulator_certification"
    ] == "proven"
    assert real_quadratic_h2_grh["class_group"][
        "zeta_bf_proof_status"
    ] == "not_checked"

    # A requested audit on x^2 + 5 (h = 2) in `proven` mode runs and
    # succeeds. The exact imaginary-quadratic route stores a verified
    # relation-saturation record at ell = 2, and the GRH-conditional BF
    # audit keeps `proven` by promoting from that stored record.
    iq_audit = run_json(
        [str(args.exe), "--coeffs", "5,0,1", "--mode", "proven",
         "--zeta-bf-audit"],
        root,
    )
    assert iq_audit["success"] is True
    assert iq_audit["certification_status"] == "proven"
    assert iq_audit["class_group"]["certification"] == "proven"
    assert iq_audit["class_group"]["order"] == "2"
    assert iq_audit["zeta_bf_audit"]["requested"] is True
    assert iq_audit["zeta_bf_audit"]["ran"] is True
    assert iq_audit["zeta_bf_audit"]["succeeded"] is True
    assert iq_audit["zeta_bf_audit"]["skip_reason"] is None

    # A requested audit that runs but does not succeed is reported as
    # ran/failed, distinct from not requested or skipped. With a BF cutoff
    # of 1 the evaluation cannot meet the error target, so the audit gate
    # fails; the audit is rolled back and the transaction stays `proven`.
    failed_audit = run_json(
        [str(args.exe), "--coeffs", "5,0,1", "--mode", "proven",
         "--zeta-bf-audit", "--bf-cutoff", "1"],
        root,
    )
    assert failed_audit["success"] is True
    assert failed_audit["certification_status"] == "proven"
    assert failed_audit["zeta_bf_audit"]["requested"] is True
    assert failed_audit["zeta_bf_audit"]["ran"] is True
    assert failed_audit["zeta_bf_audit"]["succeeded"] is False
    assert failed_audit["zeta_bf_audit"]["skip_reason"] is None

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
        # Degree one takes the exact degree-one route even for a `grh`
        # request, so no analytic index-one check is used or recorded.
        assert degree_one_grh_audited["class_group"][
            "analytic_class_regulator_status"
        ] == "not_checked"
        assert degree_one_grh_audited["class_group"][
            "analytic_class_regulator_certification"
        ] == "unknown"

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
