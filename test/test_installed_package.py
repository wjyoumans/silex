#!/usr/bin/env python3
"""Validate default and development-backend installed Silex packages."""

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


def build_command(cmake: str, build_dir: Path, config: str) -> list[str]:
    command = [cmake, "--build", str(build_dir), "--parallel", "2"]
    if config:
        command.extend(("--config", config))
    return command


def installed_executable(build_dir: Path, config: str) -> Path:
    executable = build_dir / "silex-install-consumer"
    if config:
        configured = build_dir / config / executable.name
        if configured.exists():
            executable = configured
    return executable


def assert_package_is_relocatable(
    prefix: Path, source_dir: Path, producer_build: Path
) -> None:
    cmake_files = sorted(prefix.rglob("*.cmake"))
    assert cmake_files, prefix
    forbidden = (str(source_dir), str(producer_build))
    for path in cmake_files:
        contents = path.read_text(encoding="utf-8")
        assert all(value not in contents for value in forbidden), (
            path,
            forbidden,
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--config", default="")
    parser.add_argument(
        "--backend", choices=("default", "fplll", "flatter"), default="default"
    )
    parser.add_argument(
        "--layout", choices=("default", "nondefault"), default="default"
    )
    parser.add_argument("--fplll-root", default="")
    parser.add_argument("--flatter-root", default="")
    args = parser.parse_args()

    cmake = str(args.cmake.resolve())
    source_dir = args.source_dir.resolve()
    outer_build_dir = args.build_dir.resolve()
    assert (source_dir / "CMakeLists.txt").is_file()
    assert (source_dir / "test/cmake/install-consumer/CMakeLists.txt").is_file()
    if args.layout == "nondefault":
        assert args.backend == "default"

    with tempfile.TemporaryDirectory(
        prefix=f"installed-{args.backend}-{args.layout}-consumer-",
        dir=outer_build_dir,
    ) as scratch:
        scratch_dir = Path(scratch)
        prefix = scratch_dir / "prefix"
        producer_build = outer_build_dir
        package_dir: Path | None = None
        expected_include = prefix / "include/silex/silex.hpp"
        expected_notice_dir = prefix / "share/doc/Silex"

        if args.layout == "nondefault":
            producer_build = scratch_dir / "producer"
            producer_configure = nested_cmake_configure_command(
                cmake,
                outer_build_dir,
                source_dir,
                producer_build,
            )
            producer_configure.extend(
                [
                    "-DSILEX_BUILD_TESTS=OFF",
                    "-DSILEX_BUILD_RELEASE_TESTS=OFF",
                    "-DSILEX_BUILD_EXAMPLES=OFF",
                    "-DCMAKE_INSTALL_INCLUDEDIR=silex-include",
                    "-DCMAKE_INSTALL_LIBDIR=silex-lib",
                    "-DCMAKE_INSTALL_DOCDIR=silex-doc",
                ]
            )
            if args.config:
                producer_configure.append(
                    f"-DCMAKE_BUILD_TYPE={args.config}"
                )
            run(producer_configure)
            run(build_command(cmake, producer_build, args.config))
            package_dir = prefix / "silex-lib/cmake/Silex"
            expected_include = prefix / "silex-include/silex/silex.hpp"
            expected_notice_dir = prefix / "silex-doc"

        install_command = [cmake, "--install", str(producer_build)]
        if args.config:
            install_command.extend(("--config", args.config))
        run([*install_command, "--prefix", str(prefix)])

        assert expected_include.is_file(), expected_include
        for notice_name in (
            "LICENSE",
            "NOTICE.md",
            "CITATION.cff",
            "THIRD_PARTY_NOTICES.md",
        ):
            assert (expected_notice_dir / notice_name).is_file(), (
                expected_notice_dir / notice_name
            )
        assert_package_is_relocatable(prefix, source_dir, producer_build)

        consumer_source = scratch_dir / "consumer-source"
        consumer_build = scratch_dir / "consumer-build"
        shutil.copytree(
            source_dir / "test/cmake/install-consumer", consumer_source
        )
        configure_command = nested_cmake_configure_command(
            cmake,
            outer_build_dir,
            consumer_source,
            consumer_build,
        )
        if package_dir is None:
            configure_command.append(f"-DCMAKE_PREFIX_PATH={prefix}")
        else:
            configure_command.append(f"-DSilex_DIR={package_dir}")
        for backend, root in (
            ("FPLLL", args.fplll_root),
            ("FLATTER", args.flatter_root),
        ):
            if root:
                configure_command.append(
                    f"-DSilex_{backend}_ROOT={Path(root).resolve()}"
                )
        if args.config:
            configure_command.append(f"-DCMAKE_BUILD_TYPE={args.config}")
        run(configure_command)
        run(build_command(cmake, consumer_build, args.config))
        run([str(installed_executable(consumer_build, args.config))])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
