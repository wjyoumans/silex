#!/usr/bin/env python3
"""Build and test a clean archive of the tracked Silex source tree."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

from cmake_test_driver import (
    cmake_cache_bool,
    cmake_cache_value,
    enclosing_ctest_command,
    nested_cmake_configure_command,
)

if not __debug__:
    raise RuntimeError(
        "this assertion-based test must not run with Python optimization"
    )

REQUIRED_SOURCE_ARCHIVE_MEMBERS = {
    Path("README.md"),
    Path("CONTRIBUTING.md"),
    Path("RELEASE_NOTES.md"),
    Path("LICENSE"),
    Path("NOTICE.md"),
    Path("CITATION.cff"),
    Path("SECURITY.md"),
    Path("SUPPORT.md"),
    Path("CODE_OF_CONDUCT.md"),
    Path("THIRD_PARTY_NOTICES.md"),
    Path(".github/ISSUE_TEMPLATE/bug.yml"),
    Path(".github/ISSUE_TEMPLATE/config.yml"),
    Path(".github/ISSUE_TEMPLATE/documentation.yml"),
    Path(".github/ISSUE_TEMPLATE/feature.yml"),
    Path(".github/pull_request_template.md"),
    Path(".github/workflows/ci.yml"),
    Path(".github/workflows/pages.yml"),
    Path("docs/index.rst"),
    Path("docs/conf.py"),
    Path("docs/requirements.txt"),
    Path("docs/tutorial/index.rst"),
    Path("docs/reference/index.rst"),
    Path("docs/development/index.rst"),
    Path("docs/releases/index.rst"),
    Path("bench/benchmark_contract.hpp"),
    Path("test/cmake_test_driver.py"),
    Path("test/test_native_benchmark.py"),
    Path("tools/bench/clean-benchmark-results.py"),
    Path("tools/bench/native-benchmark.py"),
}


def required_source_archive_members(project_version: str) -> set[Path]:
    return REQUIRED_SOURCE_ARCHIVE_MEMBERS | {
        Path(f"docs/releases/{project_version}.rst")
    }


def run(
    command: list[str], *, cwd: Path | None = None, timeout: float = 280.0
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
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
    return completed


def staged_index(source_dir: Path) -> str:
    completed = run(
        ["git", "ls-files", "--stage", "-z"], cwd=source_dir, timeout=30.0
    )
    assert completed.stdout
    return completed.stdout


def tracked_paths(index: str, project_version: str) -> list[Path]:
    paths: list[Path] = []
    for entry in (value for value in index.split("\0") if value):
        metadata, separator, name = entry.partition("\t")
        fields = metadata.split()
        assert separator and len(fields) == 3, entry
        _mode, _object_id, stage = fields
        assert stage == "0", f"source archive cannot contain unmerged path: {name!r}"
        paths.append(Path(name))
    assert paths
    assert len(paths) == len(set(paths)), "Git index contains duplicate paths"
    assert all(not path.is_absolute() and ".." not in path.parts for path in paths)
    missing = required_source_archive_members(project_version) - set(paths)
    assert not missing, f"required source-archive members are untracked: {missing}"
    legacy_docs = [path for path in paths if path.parts[0] == "doc"]
    assert not legacy_docs, (
        "legacy doc/ paths remain tracked after the docs/ migration: "
        + ", ".join(str(path) for path in legacy_docs)
    )
    return paths


def require_index_worktree_match(source_dir: Path) -> None:
    completed = run(
        [
            "git",
            "diff",
            "--name-only",
            "-z",
            "--no-ext-diff",
            "--ignore-submodules=none",
            "--",
        ],
        cwd=source_dir,
        timeout=30.0,
    )
    changed = [value for value in completed.stdout.split("\0") if value]
    assert not changed, (
        "source-archive validation requires tracked working-tree bytes to "
        "match the staged Git index; unstaged paths: "
        + ", ".join(repr(value) for value in changed)
    )


def materialize_index(
    source_dir: Path, destination: Path, project_version: str
) -> tuple[list[Path], str]:
    index = staged_index(source_dir)
    paths = tracked_paths(index, project_version)
    require_index_worktree_match(source_dir)
    assert staged_index(source_dir) == index, "Git index changed before checkout"
    destination.mkdir()
    run(
        [
            "git",
            "checkout-index",
            "--all",
            f"--prefix={destination.resolve()}{os.sep}",
        ],
        cwd=source_dir,
        timeout=30.0,
    )
    assert staged_index(source_dir) == index, "Git index changed during checkout"
    for relative in paths:
        path = destination / relative
        assert path.exists() or path.is_symlink(), relative
    return paths, index


def verify_version_metadata(source_dir: Path, project_version: str) -> None:
    cmake_text = (source_dir / "CMakeLists.txt").read_text(encoding="utf-8")
    project_match = re.search(
        r"project\s*\(\s*Silex\b(?P<body>.*?)\)",
        cmake_text,
        flags=re.DOTALL,
    )
    assert project_match, "CMake project(Silex ...) declaration is missing"
    cmake_version_match = re.search(
        r"\bVERSION\s+([0-9]+(?:\.[0-9]+){2})\b",
        project_match.group("body"),
    )
    assert cmake_version_match, "CMake project version is missing or malformed"
    assert cmake_version_match.group(1) == project_version, (
        "CMake project version does not match --project-version: "
        f"{cmake_version_match.group(1)!r} != {project_version!r}"
    )

    citation_text = (source_dir / "CITATION.cff").read_text(encoding="utf-8")
    citation_version_match = re.search(
        r"^version:\s*['\"]?([^'\"\s#]+)['\"]?\s*(?:#.*)?$",
        citation_text,
        flags=re.MULTILINE,
    )
    assert citation_version_match, "CITATION.cff has no top-level version"
    assert citation_version_match.group(1) == project_version, (
        "CITATION.cff version does not match --project-version: "
        f"{citation_version_match.group(1)!r} != {project_version!r}"
    )

    release_notes = (source_dir / "RELEASE_NOTES.md").read_text(
        encoding="utf-8"
    )
    assert re.search(
        rf"^#\s+Silex\s+{re.escape(project_version)}\s+release notes\s*$",
        release_notes,
        flags=re.IGNORECASE | re.MULTILINE,
    ), f"RELEASE_NOTES.md is not headed by the {project_version} release"

    release_page = source_dir / "docs" / "releases" / f"{project_version}.rst"
    release_page_text = release_page.read_text(encoding="utf-8")
    assert re.search(
        rf"\bSilex\s+{re.escape(project_version)}\b", release_page_text
    ), f"{release_page.relative_to(source_dir)} does not identify its release"


def create_archive(
    index_source: Path,
    paths: list[Path],
    archive: Path,
    root_name: str,
) -> None:
    with tarfile.open(archive, "w:gz") as output:
        for relative in paths:
            path = index_source / relative
            output.add(path, arcname=Path(root_name) / relative, recursive=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--project-version", required=True)
    parser.add_argument("--config", default="")
    args = parser.parse_args()

    assert re.fullmatch(r"[0-9]+(?:\.[0-9]+){2}", args.project_version), (
        "--project-version must be a three-component numeric version"
    )

    cmake = str(args.cmake.resolve())
    source_dir = args.source_dir.resolve()
    outer_build_dir = args.build_dir.resolve()
    assert (source_dir / "CMakeLists.txt").is_file()
    build_benchmarks = cmake_cache_bool(
        outer_build_dir, "SILEX_BUILD_BENCHMARKS"
    )
    ctest = enclosing_ctest_command(outer_build_dir)

    with tempfile.TemporaryDirectory(
        prefix="source-archive-", dir=outer_build_dir
    ) as scratch:
        scratch_dir = Path(scratch)
        root_name = f"silex-{args.project_version}"
        archive = scratch_dir / f"{root_name}.tar.gz"
        index_source = scratch_dir / "index"
        index_paths, index = materialize_index(
            source_dir, index_source, args.project_version
        )
        verify_version_metadata(index_source, args.project_version)
        create_archive(index_source, index_paths, archive, root_name)
        require_index_worktree_match(source_dir)
        assert staged_index(source_dir) == index, (
            "Git index changed while creating the source archive"
        )

        extract_dir = scratch_dir / "extract"
        extract_dir.mkdir()
        with tarfile.open(archive, "r:gz") as source:
            members = source.getmembers()
            assert members
            assert all(
                Path(member.name).parts
                and Path(member.name).parts[0] == root_name
                and ".git" not in Path(member.name).parts
                and "build" not in Path(member.name).parts[1:]
                for member in members
            )
            source.extractall(extract_dir)

        archive_source = extract_dir / root_name
        archive_build = scratch_dir / "build"
        for relative in required_source_archive_members(args.project_version):
            assert (archive_source / relative).is_file(), relative
        assert not (archive_source / ".git").exists()
        verify_version_metadata(archive_source, args.project_version)

        configure_command = nested_cmake_configure_command(
            cmake,
            outer_build_dir,
            archive_source,
            archive_build,
        )
        configure_command.extend(
            [
                "-DSILEX_BUILD_TESTS=ON",
                "-DSILEX_BUILD_RELEASE_TESTS=OFF",
                "-DSILEX_BUILD_EXAMPLES=ON",
                "-DSILEX_BUILD_DOCS=ON",
                f"-DSILEX_DOCS_RELEASE={args.project_version}",
                "-DSILEX_DOCS_CHANNEL=stable",
                f"-DSILEX_BUILD_BENCHMARKS={'ON' if build_benchmarks else 'OFF'}",
            ]
        )
        if build_benchmarks:
            benchmark_dir = cmake_cache_value(outer_build_dir, "benchmark_DIR")
            assert benchmark_dir and not benchmark_dir.endswith("-NOTFOUND"), (
                "benchmark-enabled source-archive validation requires the "
                "resolved benchmark_DIR from the enclosing build"
            )
            configure_command.append(f"-Dbenchmark_DIR={benchmark_dir}")
        if args.config:
            configure_command.append(f"-DCMAKE_BUILD_TYPE={args.config}")
        run(configure_command)
        build_command = [
            cmake,
            "--build",
            str(archive_build),
            "--parallel",
            "2",
        ]
        if args.config:
            build_command.extend(("--config", args.config))
        run(build_command)

        docs_build_command = [
            cmake,
            "--build",
            str(archive_build),
            "--parallel",
            "2",
            "--target",
            "docs-html",
        ]
        if args.config:
            docs_build_command.extend(("--config", args.config))
        run(docs_build_command)

        test_command = [
            ctest,
            "--test-dir",
            str(archive_build),
            "--output-on-failure",
        ]
        if args.config:
            test_command.extend(("-C", args.config))
        if build_benchmarks:
            # The enclosing Release suite executes every semantic row from the
            # same index-matched sources.  The archive replay still builds all
            # managed benchmark executables and retains the lightweight
            # manifest and registered-row inventory audits without recursively
            # duplicating the expensive semantic gates.
            test_command.extend(
                ("--exclude-regex", "(google-benchmark-success|semantic-gate)$")
            )
        run(test_command)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
