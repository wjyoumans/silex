#!/usr/bin/env python3
"""Build and run a clean add_subdirectory consumer."""

from __future__ import annotations

import argparse
import shutil
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


def run(command: list[str], *, timeout: float = 180.0) -> None:
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
    )
    assert completed.returncode == 0, (
        command,
        completed.returncode,
        completed.stdout,
        completed.stderr,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--config", default="")
    args = parser.parse_args()

    cmake = str(args.cmake.resolve())
    source_dir = args.source_dir.resolve()
    outer_build_dir = args.build_dir.resolve()
    fixture = source_dir / "test/cmake/add-subdirectory-consumer"
    assert (fixture / "CMakeLists.txt").is_file()

    with tempfile.TemporaryDirectory(
        prefix="add-subdirectory-consumer-", dir=outer_build_dir
    ) as scratch:
        scratch_dir = Path(scratch)
        consumer_source = scratch_dir / "source"
        consumer_build = scratch_dir / "build"
        shutil.copytree(fixture, consumer_source)
        configure_command = nested_cmake_configure_command(
            cmake,
            outer_build_dir,
            consumer_source,
            consumer_build,
        )
        configure_command.append(f"-DSILEX_SOURCE_DIR={source_dir}")
        if args.config:
            configure_command.append(f"-DCMAKE_BUILD_TYPE={args.config}")
        run(configure_command)
        build_command = [
            cmake,
            "--build",
            str(consumer_build),
            "--parallel",
            "2",
        ]
        if args.config:
            build_command.extend(("--config", args.config))
        run(build_command)

        executable = consumer_build / "silex-add-subdirectory-consumer"
        if args.config:
            configured = consumer_build / args.config / executable.name
            if configured.exists():
                executable = configured
        run([str(executable)])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
