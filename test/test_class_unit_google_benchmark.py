#!/usr/bin/env python3
"""Require semantic success from the tagged-release class/unit benchmarks."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

PROVEN_CERTIFICATION = 3.0
VERIFIED_PROOF = 2.0

EXACT_RESULT_COUNTERS = (
    "class_order_exact",
    "class_invariant_count_exact",
    "class_invariants_exact",
    "unit_rank_exact",
    "expected_unit_rank_matches_signature",
    "torsion_order_exact",
    "torsion_generator_available",
    "torsion_generator_in_order",
    "torsion_generator_unit_norm",
    "torsion_generator_exact_order",
    "free_generators_available",
    "free_generators_in_order",
    "free_generators_unit_norm",
    "stored_regulator_available",
    "recomputed_regulator_available",
    "regulators_overlap",
    "stored_regulator_in_expected_interval",
    "recomputed_regulator_in_expected_interval",
    "unit_output_exact",
)

RELEASE_BENCHMARK_NAMES = {
    "BM_class_unit_0_1_0_degree_one_proven",
    "BM_class_unit_0_1_0_real_quadratic_proven",
    "BM_class_unit_0_1_0_imag_quadratic_proven",
    "BM_class_unit_0_1_0_cubic_trivial_proven",
    "BM_class_unit_0_1_0_cubic_nontrivial_proven",
    "BM_class_unit_0_1_0_quartic_cyclotomic_proven",
    "BM_class_unit_0_1_0_quartic_noncyclotomic_proven",
    "BM_class_unit_0_1_0_quartic_disc70640_proven",
    "BM_class_unit_0_1_0_quartic_disc223479_proven",
    "BM_class_unit_0_1_0_quartic_disc35019_proven",
    "BM_class_unit_0_1_0_quartic_disc1412343_proven",
    "BM_class_unit_0_1_0_quartic_x4_minus_x_minus_1_proven",
    "BM_class_unit_0_1_0_quintic_proven",
    "BM_class_unit_0_1_0_quintic_disc11119_proven",
    "BM_class_unit_0_1_0_quintic_disc401370255_proven",
    "BM_class_unit_0_1_0_quintic_disc57895_proven",
    "BM_class_unit_0_1_0_sextic_proven",
    "BM_class_unit_random_matrix/3/4/0/1/0/iterations:1",
}

RELEASE_FILTER = (
    "^(BM_class_unit_0_1_0_.*_proven|"
    "BM_class_unit_random_matrix/3/4/0/1/0.*)$"
)


def require_counter(
    failures: list[str],
    row: dict[str, object],
    name: str,
    expected: float,
) -> None:
    value = row.get(name)
    if not isinstance(value, (int, float)) or value != expected:
        failures.append(
            f"{row.get('name', '<unnamed>')}: {name}={value!r}, "
            f"expected {expected:g}"
        )


def validate(payload: object) -> list[str]:
    if not isinstance(payload, dict):
        return ["Google Benchmark output is not a JSON object"]
    raw_rows = payload.get("benchmarks")
    if not isinstance(raw_rows, list):
        return ["Google Benchmark JSON has no benchmarks array"]
    if not all(isinstance(row, dict) for row in raw_rows):
        return ["Google Benchmark JSON contains a non-object benchmark row"]

    rows = raw_rows
    actual_names = {
        str(row.get("name", "<unnamed>")) for row in rows
    }
    failures: list[str] = []
    missing = sorted(RELEASE_BENCHMARK_NAMES - actual_names)
    unexpected = sorted(actual_names - RELEASE_BENCHMARK_NAMES)
    if missing:
        failures.append("missing release benchmark rows: " + ", ".join(missing))
    if unexpected:
        failures.append(
            "unexpected release benchmark rows: "
            + ", ".join(str(name) for name in unexpected)
        )
    if len(rows) != len(actual_names):
        failures.append("release benchmark output contains duplicate row names")

    for row in rows:
        name = str(row.get("name", "<unnamed>"))
        if row.get("error_occurred") is True:
            failures.append(
                f"{name}: Google Benchmark error: "
                f"{row.get('error_message', '<no message>')}"
            )
            continue
        require_counter(failures, row, "success", 1.0)
        require_counter(failures, row, "failure_reason", 0.0)
        require_counter(
            failures, row, "requested_cert", PROVEN_CERTIFICATION
        )
        require_counter(failures, row, "class_cert", PROVEN_CERTIFICATION)
        require_counter(failures, row, "unit_cert", PROVEN_CERTIFICATION)
        require_counter(failures, row, "fb_checked", VERIFIED_PROOF)
        require_counter(failures, row, "relation_saturation", VERIFIED_PROOF)
        require_counter(failures, row, "unit_proof", VERIFIED_PROOF)
        require_counter(failures, row, "regulator_proof", VERIFIED_PROOF)
        for counter in EXACT_RESULT_COUNTERS:
            require_counter(failures, row, counter, 1.0)
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=180.0)
    args = parser.parse_args()

    command = [
        str(args.exe),
        f"--benchmark_filter={RELEASE_FILTER}",
        "--benchmark_min_time=1x",
        "--benchmark_repetitions=1",
        "--benchmark_report_aggregates_only=false",
        "--benchmark_format=json",
    ]
    try:
        completed = subprocess.run(
            command,
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=args.timeout,
        )
    except subprocess.TimeoutExpired as exc:
        print(
            f"class/unit Google Benchmark gate timed out after "
            f"{args.timeout:g} seconds\n{exc.stderr or ''}",
            file=sys.stderr,
        )
        return 1

    if completed.returncode != 0:
        print(
            f"Google Benchmark exited with {completed.returncode}:\n"
            f"{completed.stderr}\n{completed.stdout}",
            file=sys.stderr,
        )
        return 1
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        print(
            f"invalid Google Benchmark JSON: {exc}\n"
            f"stderr:\n{completed.stderr}\nstdout:\n{completed.stdout}",
            file=sys.stderr,
        )
        return 1

    failures = validate(payload)
    if failures:
        print("class/unit Google Benchmark release gate failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        if completed.stderr:
            print(f"benchmark stderr:\n{completed.stderr}", file=sys.stderr)
        return 1

    print(f"validated {len(RELEASE_BENCHMARK_NAMES)} release benchmark rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
