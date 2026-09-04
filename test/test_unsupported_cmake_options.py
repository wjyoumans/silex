#!/usr/bin/env python3
"""Exercise fail-closed handling for unsupported CMake options."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

from cmake_test_driver import nested_cmake_configure_command

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )


OPTIONS = {
    "spdlog": "SILEX_WITH_SPDLOG",
    "tracy": "SILEX_WITH_TRACY",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--unsupported-backend", choices=tuple(OPTIONS), required=True)
    args = parser.parse_args()

    option = OPTIONS[args.unsupported_backend]
    cmake = str(args.cmake.resolve())
    source_dir = args.source_dir.resolve()
    outer_build_dir = args.build_dir.resolve()
    with tempfile.TemporaryDirectory(
        prefix=f"unsupported-{args.unsupported_backend}-", dir=outer_build_dir
    ) as scratch:
        configure_command = nested_cmake_configure_command(
            cmake,
            outer_build_dir,
            source_dir,
            Path(scratch),
        )
        configure_command.extend(
            [
                "-DSILEX_BUILD_TESTS=OFF",
                "-DSILEX_BUILD_EXAMPLES=OFF",
                "-DSILEX_BUILD_BENCHMARKS=OFF",
                f"-D{option}=ON",
            ]
        )
        completed = subprocess.run(
            configure_command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=60.0,
        )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, (option, completed.stdout, completed.stderr)
    assert f"{option} is not implemented" in output, output
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
