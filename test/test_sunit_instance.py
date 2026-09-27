#!/usr/bin/env python3
"""Smoke-test the native S-unit fixture and instance protocol."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )

PRIME_INDEX_CONVENTION = (
    "zero-based authoritative manifest order of exact two-generator ideals "
    "(p, beta_power_basis) within each rational-prime decomposition"
)


def run_json(command: list[str], root: Path) -> dict[str, Any]:
    completed = subprocess.run(
        command,
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


def check_timeout_protocol(root: Path, manifest: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="silex-sunit-timeout-") as temporary:
        temporary_path = Path(temporary)
        fake_executable = temporary_path / "sleeping-instance"
        fake_executable.write_text(
            f"#!{sys.executable}\n"
            "import sys\n"
            "import time\n"
            "sys.stdout.write('partial stdout\\n')\n"
            "sys.stdout.flush()\n"
            "sys.stderr.write('partial stderr\\n')\n"
            "sys.stderr.flush()\n"
            "time.sleep(60)\n",
            encoding="utf-8",
        )
        fake_executable.chmod(0o755)
        output = temporary_path / "timeout.json"
        completed = subprocess.run(
            [
                sys.executable,
                str(root / "tools/bench/run-sunit-instance.py"),
                "--exe",
                str(fake_executable),
                "--manifest",
                str(manifest),
                "--field-id",
                "cubic_x3_minus_2_empty_s",
                "--timeout",
                "0.1",
                "--out",
                str(output),
            ],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=5.0,
        )
        assert completed.returncode == 124, (
            completed.returncode,
            completed.stdout,
            completed.stderr,
        )
        payload = json.loads(output.read_text(encoding="utf-8"))
        assert payload["success"] is False
        assert payload["timeout"] is True
        assert payload["stdout"] == "partial stdout\n"
        assert payload["stderr"] == "partial stderr\n"


def check_repeated_rational_prime_selectors(
    root: Path, executable: Path, manifest: dict[str, Any]
) -> None:
    with tempfile.TemporaryDirectory(prefix="silex-sunit-selectors-") as temporary:
        manifest_copy = json.loads(json.dumps(manifest))
        row = next(
            entry
            for entry in manifest_copy["fields"]
            if entry["id"] == "real_quadratic_5_split_11"
        )
        row["selected_primes"] = [
            {"p": 11, "index": 1},
            {"p": 11, "index": 0},
        ]
        manifest_path = Path(temporary) / "sunit-fields.json"
        manifest_path.write_text(json.dumps(manifest_copy), encoding="utf-8")

        payload = run_json(
            [
                sys.executable,
                str(root / "tools/bench/run-sunit-instance.py"),
                "--exe",
                str(executable),
                "--manifest",
                str(manifest_path),
                "--field-id",
                row["id"],
                "--timeout",
                "20",
            ],
            root,
        )
        selected = payload["sunit"]["selected_primes"]
        complete = payload["sunit"]["canonical_prime_decompositions"]
        assert [entry["canonical_index"] for entry in selected] == [1, 0]
        assert [entry["canonical_index"] for entry in complete] == [0, 1]
        assert selected == list(reversed(complete))


def check_outside_prime_search_extends_beyond_fixed_prefix(
    root: Path, executable: Path, manifest: dict[str, Any]
) -> None:
    with tempfile.TemporaryDirectory(prefix="silex-sunit-outside-") as temporary:
        manifest_copy = json.loads(json.dumps(manifest))
        row = next(
            entry
            for entry in manifest_copy["fields"]
            if entry["id"] == "real_quadratic_5_split_11"
        )
        rational_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
        row["selected_primes"] = [
            {"p": rational_prime, "index": "all"}
            for rational_prime in rational_primes
        ]
        witnesses: list[dict[str, Any]] = []
        for rational_prime in rational_primes:
            roots = [
                residue
                for residue in range(rational_prime)
                if (residue * residue - residue - 1) % rational_prime == 0
            ]
            if not roots:
                witnesses.append(
                    {
                        "p": rational_prime,
                        "canonical_index": 0,
                        "e": 1,
                        "f": 2,
                        "beta_power_basis": [str(rational_prime), "0"],
                    }
                )
                continue
            ramification_index = 2 if len(roots) == 1 else 1
            for canonical_index, residue in enumerate(roots):
                witnesses.append(
                    {
                        "p": rational_prime,
                        "canonical_index": canonical_index,
                        "e": ramification_index,
                        "f": 1,
                        "beta_power_basis": [str(-residue), "1"],
                    }
                )
        row["prime_ideal_witnesses"] = witnesses
        row["expected"] = {"nonunit_rank": len(witnesses)}
        manifest_path = Path(temporary) / "sunit-fields.json"
        manifest_path.write_text(json.dumps(manifest_copy), encoding="utf-8")

        completed = subprocess.run(
            [
                sys.executable,
                str(root / "tools/bench/run-sunit-instance.py"),
                "--exe",
                str(executable),
                "--manifest",
                str(manifest_path),
                "--field-id",
                row["id"],
                "--timeout",
                "30",
            ],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=35.0,
        )
        assert completed.returncode == 0, (
            completed.returncode,
            completed.stdout,
            completed.stderr,
        )
        payload = json.loads(completed.stdout)
        assert payload["success"] is True
        assert payload["sunit"]["membership"]["outside_outcome"] == "not_sunit"


def assert_sunit_instance(instance: dict[str, Any], row: dict[str, Any]) -> None:
    assert instance["field_id"] == row["id"]
    assert instance["success"] is True
    assert instance["timeout"] is False
    assert instance["manifest_expectations_match"] is True
    assert instance["engine"] == "silex"
    assert instance["maximal_order_discriminant"] == row[
        "maximal_order_discriminant"
    ]
    assert instance["final_result_published"] is True
    for key in ("certification_status", "class_group_proof_status",
                "unit_group_proof_status"):
        assert instance[key] == "proven"
    assert instance["regulator_proof_status"] == "verified"
    assert instance["class_group"]["certification"] == "proven"
    assert instance["unit_group"]["certification"] == "proven"
    for key in ("factor_base_generation_status", "unit_proof_status",
                "regulator_proof_status"):
        assert instance["class_group"][key] == "verified"
    dirichlet_route = row["id"] == "real_quadratic_210_first_over_2"
    assert instance["class_group"]["relation_saturation_status"] == (
        "not_checked" if dirichlet_route else "verified"
    )
    if dirichlet_route:
        # The source pair uses the canonical Dirichlet index-one gate;
        # S-unit publication preserves its unused component receipts.
        assert instance["class_group"][
            "analytic_class_regulator_status"
        ] == "verified"
        assert instance["class_group"][
            "analytic_class_regulator_certification"
        ] == "proven"
        assert instance["class_group"]["zeta_bf_proof_status"] == "not_checked"
    assert "algorithm" not in instance
    assert "fallback_used" not in instance
    assert "--class-unit-route" not in instance["cmd"]

    sunit = instance["sunit"]
    assert sunit["success"] is True
    assert sunit["final_result_published"] is True
    assert sunit["failure_stage"] is None
    assert "source_fallback_used" not in sunit
    assert "computation_fallback_used" not in sunit
    s_class = sunit["s_class_group"]
    s_units = sunit["s_unit_group"]
    expected = row["expected"]
    assert s_class["order"] == expected["s_class_order"]
    assert s_class["invariants"] == expected["s_class_invariants"]
    for key in ("torsion_order", "ordinary_free_rank", "nonunit_rank",
                "free_rank", "valuation_lattice_index"):
        assert s_units[key] == expected[key]
    for group in (s_class, s_units):
        assert group["certification_status"] == "proven"
        assert group["proof_status"] == "verified"
    assert s_units["regulator_proof_status"] == "verified"
    assert math.isfinite(s_units["regulator_midpoint"])
    assert s_units["regulator_midpoint"] > 0
    assert len(s_units["valuation_matrix"]) == expected["nonunit_rank"]
    assert all(len(values) == expected["nonunit_rank"]
               for values in s_units["valuation_matrix"])

    complete = sunit["canonical_prime_decompositions"]
    witnesses = row["prime_ideal_witnesses"]
    assert len(complete) == len(witnesses)
    for descriptor, witness in zip(complete, witnesses):
        assert descriptor["p"] == str(witness["p"])
        for key in ("canonical_index", "e", "f", "beta_power_basis"):
            assert descriptor[key] == witness[key]
    selected = []
    for selector in row["selected_primes"]:
        decomposition = [entry for entry in complete
                         if entry["p"] == str(selector["p"])]
        assert [entry["canonical_index"] for entry in decomposition] == list(
            range(len(decomposition))
        )
        assert sum(entry["e"] * entry["f"] for entry in decomposition) == (
            len(row["coefficients_low_to_high"]) - 1
        )
        selected.extend(decomposition if selector["index"] == "all" else
                        [decomposition[selector["index"]]])
    assert sunit["selected_primes"] == selected
    assert len(selected) == expected["nonunit_rank"]
    assert sunit["membership"] == {
        "status": "verified",
        "mixed_round_trip_verified": True,
        "verified_round_trip_count": 1,
        "mixed_outcome": "verified",
        "outside_support_rejected": True,
        "outside_outcome": "not_sunit",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(args.manifest.read_text())
    fields = {row["id"]: row for row in manifest["fields"]}
    assert manifest["schema_version"] == 2
    assert manifest["prime_index_convention"] == PRIME_INDEX_CONVENTION
    assert len(manifest["fields"]) == len(fields) == 6
    assert set(fields) == {
        "cubic_x3_minus_2_empty_s",
        "cubic_x3_minus_2_regression",
        "imaginary_quadratic_23_first_over_2",
        "real_quadratic_210_first_over_2",
        "real_quadratic_5_split_11",
        "real_quadratic_5_ramified_5",
    }
    assert fields["cubic_x3_minus_2_empty_s"]["selected_primes"] == []
    assert fields["cubic_x3_minus_2_regression"]["expected"][
        "nonunit_rank"
    ] == 5
    assert fields["imaginary_quadratic_23_first_over_2"]["expected"][
        "valuation_lattice_index"
    ] == "3"
    assert fields["real_quadratic_210_first_over_2"]["expected"][
        "s_class_invariants"
    ] == ["2"]

    check_timeout_protocol(root, args.manifest)
    check_repeated_rational_prime_selectors(root, args.exe, manifest)
    check_outside_prime_search_extends_beyond_fixed_prefix(root, args.exe, manifest)

    instances: dict[str, dict[str, Any]] = {}
    for field_id, row in fields.items():
        instance = run_json(
            [
                sys.executable,
                str(root / "tools/bench/run-sunit-instance.py"),
                "--exe", str(args.exe),
                "--manifest", str(args.manifest),
                "--field-id", field_id,
                "--timeout", "20",
            ],
            root,
        )
        assert_sunit_instance(instance, row)
        instances[field_id] = instance

    real_quadratic = instances["real_quadratic_210_first_over_2"]
    assert real_quadratic["class_group"]["order"] == "4"
    assert real_quadratic["class_group"]["invariants"] == ["2", "2"]
    assert real_quadratic["unit_group"]["free_rank"] == 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
