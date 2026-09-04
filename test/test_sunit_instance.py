#!/usr/bin/env python3
"""Smoke-test the native S-unit fixture and instance protocol."""

from __future__ import annotations

import argparse
import json
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
    assert len(fields) == 6
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

    instance = run_json(
        [
            sys.executable,
            str(root / "tools/bench/run-sunit-instance.py"),
            "--exe",
            str(args.exe),
            "--manifest",
            str(args.manifest),
            "--field-id",
            "cubic_x3_minus_2_empty_s",
            "--timeout",
            "20",
        ],
        root,
    )
    assert instance["success"] is True
    assert instance["timeout"] is False
    assert instance["manifest_expectations_match"] is True
    assert instance["engine"] == "silex"
    assert "algorithm" not in instance
    assert "fallback_used" not in instance
    assert instance["sunit"]["selected_primes"] == []
    assert instance["sunit"]["canonical_prime_decompositions"] == []
    assert instance["sunit"]["membership"]["status"] == "verified"
    assert "source_fallback_used" not in instance["sunit"]
    assert "computation_fallback_used" not in instance["sunit"]
    assert "--class-unit-route" not in instance["cmd"]

    real_quadratic = run_json(
        [
            sys.executable,
            str(root / "tools/bench/run-sunit-instance.py"),
            "--exe",
            str(args.exe),
            "--manifest",
            str(args.manifest),
            "--field-id",
            "real_quadratic_210_first_over_2",
            "--timeout",
            "20",
        ],
        root,
    )
    assert real_quadratic["success"] is True
    assert real_quadratic["timeout"] is False
    assert real_quadratic["manifest_expectations_match"] is True
    assert "algorithm" not in real_quadratic
    assert "--class-unit-route" not in real_quadratic["cmd"]
    assert "fallback_used" not in real_quadratic
    assert real_quadratic["certification_status"] == "proven"
    assert real_quadratic["class_group"]["order"] == "4"
    assert real_quadratic["class_group"]["invariants"] == ["2", "2"]
    assert real_quadratic["unit_group"]["free_rank"] == 1
    assert real_quadratic["unit_group"]["certification"] == "proven"
    real_quadratic_sunit = real_quadratic["sunit"]
    assert real_quadratic_sunit["final_result_published"] is True
    assert "source_fallback_used" not in real_quadratic_sunit
    assert "computation_fallback_used" not in real_quadratic_sunit
    assert real_quadratic_sunit["s_class_group"]["order"] == "2"
    assert real_quadratic_sunit["s_class_group"]["invariants"] == ["2"]
    assert (
        real_quadratic_sunit["s_class_group"]["proof_status"] == "verified"
    )
    assert real_quadratic_sunit["s_unit_group"]["ordinary_free_rank"] == 1
    assert real_quadratic_sunit["s_unit_group"]["nonunit_rank"] == 1
    assert real_quadratic_sunit["s_unit_group"]["free_rank"] == 2
    assert (
        real_quadratic_sunit["s_unit_group"]["valuation_lattice_index"] == "2"
    )
    assert [
        descriptor["canonical_index"]
        for descriptor in real_quadratic_sunit["selected_primes"]
    ] == [0]
    real_quadratic_decomposition = real_quadratic_sunit[
        "canonical_prime_decompositions"
    ]
    assert [
        descriptor["canonical_index"]
        for descriptor in real_quadratic_decomposition
    ] == list(range(len(real_quadratic_decomposition)))
    assert real_quadratic_sunit["selected_primes"][0] == (
        real_quadratic_decomposition[0]
    )
    assert sum(
        descriptor["e"] * descriptor["f"]
        for descriptor in real_quadratic_decomposition
    ) == 2
    assert (
        real_quadratic_sunit["s_unit_group"]["proof_status"] == "verified"
    )
    assert real_quadratic_sunit["membership"]["status"] == "verified"

    split_prime = run_json(
        [
            sys.executable,
            str(root / "tools/bench/run-sunit-instance.py"),
            "--exe",
            str(args.exe),
            "--manifest",
            str(args.manifest),
            "--field-id",
            "real_quadratic_5_split_11",
            "--timeout",
            "20",
        ],
        root,
    )
    assert split_prime["success"] is True
    assert [
        descriptor["canonical_index"]
        for descriptor in split_prime["sunit"]["selected_primes"]
    ] == [0, 1]
    assert split_prime["sunit"]["selected_primes"] == split_prime["sunit"][
        "canonical_prime_decompositions"
    ]

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
