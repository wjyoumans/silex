"""Shared CMake invocation support for nested test-driver builds."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from pathlib import Path


def cmake_cache_value(build_dir: Path, name: str) -> str | None:
    cache = build_dir / "CMakeCache.txt"
    assert cache.is_file(), f"missing enclosing CMake cache: {cache}"
    prefix = f"{name}:"
    for line in cache.read_text(encoding="utf-8").splitlines():
        if not line.startswith(prefix):
            continue
        declaration, separator, value = line.partition("=")
        assert separator and declaration.partition(":")[0] == name
        return value
    return None


def cmake_cache_bool(build_dir: Path, name: str) -> bool:
    value = cmake_cache_value(build_dir, name)
    assert value is not None, f"missing {name} in enclosing CMake cache"
    normalized = value.strip().upper()
    if normalized in {"1", "ON", "YES", "TRUE", "Y"}:
        return True
    if normalized in {
        "",
        "0",
        "OFF",
        "NO",
        "FALSE",
        "N",
        "IGNORE",
        "NOTFOUND",
    } or normalized.endswith("-NOTFOUND"):
        return False
    try:
        numeric = Decimal(normalized)
    except InvalidOperation:
        numeric = None
    if numeric is not None and numeric.is_finite():
        return numeric != 0
    raise AssertionError(
        f"unexpected {name} value in enclosing CMake cache: {value!r}"
    )


def nested_cmake_configure_command(
    cmake: str,
    outer_build_dir: Path,
    source_dir: Path,
    build_dir: Path,
) -> list[str]:
    """Start a nested configure with the enclosing generator and toolchain."""
    generator = cmake_cache_value(outer_build_dir, "CMAKE_GENERATOR")
    assert generator, "missing CMAKE_GENERATOR in enclosing CMake cache"
    command = [cmake, "-G", generator]

    generator_platform = cmake_cache_value(
        outer_build_dir, "CMAKE_GENERATOR_PLATFORM"
    )
    if generator_platform:
        command.extend(("-A", generator_platform))
    generator_toolset = cmake_cache_value(
        outer_build_dir, "CMAKE_GENERATOR_TOOLSET"
    )
    if generator_toolset:
        command.extend(("-T", generator_toolset))

    command.extend(("-S", str(source_dir), "-B", str(build_dir)))

    generator_instance = cmake_cache_value(
        outer_build_dir, "CMAKE_GENERATOR_INSTANCE"
    )
    if generator_instance:
        command.append(f"-DCMAKE_GENERATOR_INSTANCE={generator_instance}")
    make_program = cmake_cache_value(outer_build_dir, "CMAKE_MAKE_PROGRAM")
    if make_program:
        command.append(f"-DCMAKE_MAKE_PROGRAM={make_program}")

    toolchain_file = cmake_cache_value(outer_build_dir, "CMAKE_TOOLCHAIN_FILE")
    if toolchain_file:
        command.append(f"-DCMAKE_TOOLCHAIN_FILE={toolchain_file}")
    else:
        cxx_compiler = cmake_cache_value(
            outer_build_dir, "CMAKE_CXX_COMPILER"
        )
        assert cxx_compiler, "missing CMAKE_CXX_COMPILER in enclosing CMake cache"
        command.append(f"-DCMAKE_CXX_COMPILER={cxx_compiler}")
    return command


def enclosing_ctest_command(outer_build_dir: Path) -> str:
    ctest = cmake_cache_value(outer_build_dir, "CMAKE_CTEST_COMMAND")
    assert ctest, "missing CMAKE_CTEST_COMMAND in enclosing CMake cache"
    return ctest
