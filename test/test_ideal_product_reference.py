#!/usr/bin/env python3
"""Compare native product lattices with fixed exact power-basis expectations."""

import argparse
import json
import subprocess
from fractions import Fraction
from pathlib import Path


def inverse(matrix):
    n = len(matrix)
    augmented = [list(row) + [Fraction(i == j) for j in range(n)]
                 for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if augmented[i][j])
        augmented[j], augmented[pivot] = augmented[pivot], augmented[j]
        scale = augmented[j][j]
        augmented[j] = [x / scale for x in augmented[j]]
        for i in range(n):
            if i != j:
                scale = augmented[i][j]
                augmented[i] = [x - scale * y for x, y in zip(augmented[i], augmented[j])]
    return [row[n:] for row in augmented]


def integral_coordinates(left, right_inverse):
    n = len(left)
    return all(sum(left[i][k] * right_inverse[k][j] for k in range(n)).denominator == 1
               for i in range(n) for j in range(n))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.fixtures.read_text())
    assert data["schema_version"] == 1
    assert {(case["degree"], case["kind"]) for case in data["cases"]} == {
        (degree, 3) for degree in (2, 3, 4, 6, 8, 12)}
    completed = subprocess.run([args.exe, "--dump-product-fixtures"],
                               check=True, capture_output=True, text=True, timeout=60)
    actual = {}
    for line in completed.stdout.splitlines():
        degree, kind, text = line.split(" ", 2)
        key = int(degree), int(kind)
        assert key not in actual
        actual[key] = [[Fraction(entry) for entry in row.split(",")]
                       for row in text.removeprefix("[").removesuffix("]").split(";")]
    assert len(actual) == 30
    for case in data["cases"]:
        key = case["degree"], case["kind"]
        expected = [[Fraction(entry) for entry in row] for row in case["power_basis_rows"]]
        result = actual[key]
        assert len(result) == key[0] and all(len(row) == key[0] for row in result)
        # Mutual integral containment is exact lattice equality, independent of
        # the chosen row/column HNF convention or a unimodular basis change.
        assert integral_coordinates(result, inverse(expected)), key
        assert integral_coordinates(expected, inverse(result)), key
    print(f"verified {len(data['cases'])} exact reference lattices")


if __name__ == "__main__":
    main()
