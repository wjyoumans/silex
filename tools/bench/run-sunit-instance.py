#!/usr/bin/env python3
"""Run one Silex S-class/S-unit fixture with bounded process time."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


PRIME_INDEX_CONVENTION = (
    "zero-based authoritative manifest order of exact two-generator ideals "
    "(p, beta_power_basis) within each rational-prime decomposition"
)
MAX_RATIONAL_PRIME = (1 << 63) - 1
MAX_WITNESS_DEGREE = 16
MAX_RATIONAL_TEXT_CHARS = 8192


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_field(path: Path, field_id: str) -> dict[str, Any]:
    manifest = json.loads(path.read_text())
    if (
        type(manifest) is not dict
        or manifest.get("schema_version") != 2
        or manifest.get("prime_index_convention") != PRIME_INDEX_CONVENTION
        or type(manifest.get("fields")) is not list
    ):
        raise ValueError("invalid S-unit manifest schema")
    seen_ids: set[str] = set()
    selected_row: dict[str, Any] | None = None
    for row in manifest["fields"]:
        if type(row) is not dict or type(row.get("id")) is not str or not row["id"]:
            raise ValueError("invalid S-unit manifest field")
        if row["id"] in seen_ids:
            raise ValueError(f"duplicate S-unit field id: {row['id']}")
        seen_ids.add(row["id"])
        coefficient_arg(row)
        witness_args(row)
        if row.get("id") == field_id:
            selected_row = row
    if selected_row is not None:
        return selected_row
    raise ValueError(f"unknown S-unit field id: {field_id}")


def coefficient_arg(row: dict[str, Any]) -> str:
    values = row.get("coefficients_low_to_high")
    if (
        type(values) is not list
        or not 2 <= len(values) <= MAX_WITNESS_DEGREE + 1
        or any(type(value) is not int for value in values)
        or values[-1] != 1
    ):
        raise ValueError("field row is missing coefficients_low_to_high")
    return ",".join(str(value) for value in values)


def rational_prime_is_valid(value: Any) -> bool:
    if type(value) is not int or value < 2 or value > MAX_RATIONAL_PRIME:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if value in small_primes:
        return True
    if any(value % prime == 0 for prime in small_primes):
        return False
    odd_part = value - 1
    power_of_two = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        power_of_two += 1
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(power_of_two - 1):
            witness = pow(witness, 2, value)
            if witness == value - 1:
                break
        else:
            return False
    return True


def canonical_rational_text(value: Any) -> bool:
    if type(value) is not str or not value or len(value) > MAX_RATIONAL_TEXT_CHARS:
        return False
    if "/" not in value:
        return (
            re.fullmatch(r"(?:0|-?[1-9][0-9]*)", value) is not None
            and value != "-0"
        )
    if value.count("/") != 1:
        return False
    numerator, denominator = value.split("/", 1)
    if (
        re.fullmatch(r"-?[1-9][0-9]*", numerator) is None
        or re.fullmatch(r"[1-9][0-9]*", denominator) is None
    ):
        return False
    parsed = Fraction(int(numerator), int(denominator))
    return parsed.denominator > 1 and math.gcd(abs(parsed.numerator), parsed.denominator) == 1 and f"{parsed.numerator}/{parsed.denominator}" == value


def witness_args(row: dict[str, Any]) -> list[str]:
    coefficients = row.get("coefficients_low_to_high")
    selectors = row.get("selected_primes")
    witnesses = row.get("prime_ideal_witnesses")
    discriminant = row.get("maximal_order_discriminant")
    if (
        type(coefficients) is not list
        or type(selectors) is not list
        or type(witnesses) is not list
        or type(discriminant) is not str
        or re.fullmatch(r"-?[1-9][0-9]*", discriminant) is None
    ):
        raise ValueError("field row has invalid prime witness metadata")
    degree = len(coefficients) - 1
    selectors_by_prime: dict[int, list[str | int]] = {}
    for selector in selectors:
        if type(selector) is not dict or set(selector) != {"p", "index"}:
            raise ValueError("field row has invalid selected_primes")
        p = selector.get("p")
        index = selector.get("index")
        if (
            type(p) is not int
            or not rational_prime_is_valid(p)
            or not (index == "all" or (type(index) is int and index >= 0))
            or index in selectors_by_prime.setdefault(p, [])
        ):
            raise ValueError("field row has invalid selected_primes")
        selectors_by_prime[p].append(index)
    if any("all" in values and values != ["all"] for values in selectors_by_prime.values()):
        raise ValueError("field row has contradictory selected_primes")

    witnesses_by_prime: dict[int, list[dict[str, Any]]] = {}
    for witness in witnesses:
        if type(witness) is not dict or set(witness) != {
            "p",
            "canonical_index",
            "e",
            "f",
            "beta_power_basis",
        }:
            raise ValueError("field row has invalid prime_ideal_witnesses")
        p = witness.get("p")
        index = witness.get("canonical_index")
        e = witness.get("e")
        f = witness.get("f")
        beta = witness.get("beta_power_basis")
        if (
            type(p) is not int
            or not rational_prime_is_valid(p)
            or type(index) is not int
            or index < 0
            or type(e) is not int
            or e < 1
            or type(f) is not int
            or f < 1
            or type(beta) is not list
            or len(beta) != degree
            or any(not canonical_rational_text(value) for value in beta)
        ):
            raise ValueError("field row has invalid prime_ideal_witnesses")
        witnesses_by_prime.setdefault(p, []).append(witness)
    if list(witnesses_by_prime) != list(selectors_by_prime):
        raise ValueError("prime witnesses do not cover selected rational primes")

    selection_indices: dict[tuple[int, int], int] = {}
    selected_count = 0
    for p, selected_indices in selectors_by_prime.items():
        complete = witnesses_by_prime[p]
        if (
            [witness["canonical_index"] for witness in complete]
            != list(range(len(complete)))
            or sum(witness["e"] * witness["f"] for witness in complete) != degree
        ):
            raise ValueError("prime witnesses are not a complete decomposition")
        expanded_indices = (
            range(len(complete))
            if selected_indices == ["all"]
            else selected_indices
        )
        for index in expanded_indices:
            if type(index) is not int or index >= len(complete):
                raise ValueError("field row has invalid selected_primes")
            selection_indices[(p, index)] = selected_count
            selected_count += 1

    arguments: list[str] = []
    for p in selectors_by_prime:
        complete = witnesses_by_prime[p]
        for witness in complete:
            selection_index = selection_indices.get(
                (p, witness["canonical_index"]), -1
            )
            arguments.append(
                f"{p}:{witness['canonical_index']}:{selection_index}:"
                + ",".join(witness["beta_power_basis"])
            )
    expected = row.get("expected")
    if type(expected) is not dict or expected.get("nonunit_rank") != selected_count:
        raise ValueError("prime witness selection rank does not match expectations")
    return arguments


def expected_matches(payload: dict[str, Any], row: dict[str, Any]) -> bool:
    expected = row.get("expected", {})
    sunit = payload.get("sunit", {})
    s_class = sunit.get("s_class_group", {})
    s_units = sunit.get("s_unit_group", {})
    actual = {
        "s_class_order": s_class.get("order"),
        "s_class_invariants": s_class.get("invariants"),
        "torsion_order": s_units.get("torsion_order"),
        "ordinary_free_rank": s_units.get("ordinary_free_rank"),
        "nonunit_rank": s_units.get("nonunit_rank"),
        "free_rank": s_units.get("free_rank"),
        "valuation_lattice_index": s_units.get("valuation_lattice_index"),
    }
    return all(actual.get(key) == value for key, value in expected.items())


def write_result(payload: dict[str, Any], out: Path | None) -> None:
    text = json.dumps(payload, indent=2) + "\n"
    if out is None:
        sys.stdout.write(text)
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)


def timeout_stream_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def main() -> int:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--build-dir", type=Path, default=root / "build/default")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "test/data/sunit_fields.json",
    )
    parser.add_argument("--field-id", required=True)
    parser.add_argument("--precision", type=int, default=128)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    try:
        row = load_field(args.manifest, args.field_id)
        coeffs = coefficient_arg(row)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        write_result({"success": False, "error": str(exc)}, args.out)
        return 2

    exe = args.exe or args.build_dir / "silex-class-unit-instance"
    if not exe.exists():
        write_result(
            {"success": False, "error": f"missing executable: {exe}"},
            args.out,
        )
        return 2

    cmd = [
        str(exe),
        "--coeffs",
        coeffs,
        "--mode",
        "proven",
        "--precision",
        str(args.precision),
        "--compute-sunit",
    ]
    for witness in witness_args(row):
        cmd.extend(["--s-prime-witness", witness])
    option_names = {
        "factor_base_bound": "--factor-base-bound",
        "max_candidates": "--max-candidates",
        "max_relations": "--max-relations",
        "bf_cutoff": "--bf-cutoff",
    }
    for key, option in option_names.items():
        value = row.get("class_unit_options", {}).get(key)
        if value is not None:
            cmd.extend([option, str(value)])
    if args.profile:
        cmd.append("--profile")

    try:
        completed = subprocess.run(
            cmd,
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=args.timeout,
        )
    except subprocess.TimeoutExpired as exc:
        write_result(
            {
                "success": False,
                "timeout": True,
                "timeout_seconds": args.timeout,
                "field_id": args.field_id,
                "cmd": cmd,
                "stdout": timeout_stream_text(exc.stdout),
                "stderr": timeout_stream_text(exc.stderr),
            },
            args.out,
        )
        return 124

    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        write_result(
            {
                "success": False,
                "error": f"invalid Silex JSON output: {exc}",
                "returncode": completed.returncode,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
            },
            args.out,
        )
        return 2

    matches = expected_matches(payload, row)
    payload["field_id"] = args.field_id
    payload["manifest_expectations_match"] = matches
    payload["cmd"] = cmd
    payload["stderr"] = completed.stderr
    payload["timeout"] = False
    sunit = payload.get("sunit")
    payload["success"] = (
        completed.returncode == 0
        and payload.get("success") is True
        and type(sunit) is dict
        and sunit.get("success") is True
        and matches
    )
    write_result(payload, args.out)
    return 0 if payload["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
