#!/usr/bin/env python3
"""Verify the compact source-neutral class/unit replay manifest."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )


def load_tool(path: Path, module_name: str) -> object:
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_field_id_validation(root: Path) -> None:
    tools = (
        load_tool(
            root / "tools/bench/run-class-unit-priority-gate.py",
            "silex_priority_gate",
        ),
        load_tool(
            root / "tools/bench/gdb-class-unit-timeout.py",
            "silex_gdb_timeout",
        ),
    )
    for tool in tools:
        validator = getattr(tool, "validated_field_id")
        for unsafe in (
            "../escape",
            "../../CMakeLists.txt",
            "/tmp/absolute",
            "",
            None,
            7,
        ):
            try:
                validator(unsafe)
            except ValueError:
                pass
            else:
                raise AssertionError(f"unsafe field id accepted: {unsafe!r}")
        for safe in ("field", "field-1", "field_2", "field.3"):
            assert validator(safe) == safe


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    assert manifest["schema_version"] == 1
    assert "source-neutral" in manifest["fixture_contract"]
    fields = manifest["fields"]
    assert isinstance(fields, list)
    assert len(fields) == 17
    assert len({row["id"] for row in fields}) == len(fields)
    assert all(row["expected_success"] is True for row in fields)
    assert all(row["mode"] in {"proven", "grh"} for row in fields)

    proven = [row for row in fields if row["status"] == "must_pass_fast"]
    grh = [row for row in fields if row["status"] == "grh_certification"]
    assert len(proven) == 12
    assert len(grh) == 5
    assert all(row["mode"] == "proven" for row in proven)
    assert all(row["mode"] == "grh" for row in grh)
    assert {row["id"] for row in proven} == {
        "real_quadratic_5_proven",
        "imaginary_quadratic_47_proven",
        "cubic_disc81_proven",
        "quartic_x4_plus_1_proven",
        "cubic_disc643_nontrivial_proven",
        "cubic_seeded_h4_proven",
        "quartic_disc892_proven",
        "quartic_disc70640_proven",
        "quartic_disc223479_proven",
        "quartic_disc35019_proven",
        "quartic_disc1412343_proven",
        "quintic_disc401370255_proven",
    }
    assert {tuple(row["coefficients_low_to_high"]) for row in grh} == {
        (34, 0, 1),
        (46, 0, 1),
        (66, 0, 1),
        (185, 0, 1),
        (47, 0, 1),
    }

    campaign_only_keys = {
        "benchmark_role",
        "optimization_external_engines",
        "pa" + "ri_elapsed_ms",
    }
    assert all(not campaign_only_keys.intersection(row) for row in fields)
    serialized = json.dumps(manifest).lower()
    external_names = ("pa" + "ri", "hec" + "ke", "magma")
    assert all(name not in serialized for name in external_names)
    check_field_id_validation(Path(__file__).resolve().parents[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
