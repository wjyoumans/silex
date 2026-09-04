#!/usr/bin/env python3
"""Remove complete ignored native-benchmark scratch artifact families."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


ARTIFACT_SUFFIXES = (
    ".metadata.json",
    ".google.json",
    ".stdout.txt",
    ".stderr.txt",
    ".json",
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def is_result_artifact(path: Path) -> bool:
    return any(path.name.endswith(suffix) for suffix in ARTIFACT_SUFFIXES)


def kept_result_family(raw_path: Path) -> set[Path]:
    """Return the raw result and companions produced for that exact path."""
    stem = raw_path.name[:-5] if raw_path.name.endswith(".json") else raw_path.name
    return {
        raw_path,
        raw_path.with_name(stem + ".metadata.json"),
        raw_path.with_name(stem + ".google.json"),
        raw_path.with_name(stem + ".stdout.txt"),
        raw_path.with_name(stem + ".stderr.txt"),
    }


def ignored_result_artifacts(root: Path) -> list[Path]:
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "ls-files",
            "--others",
            "--ignored",
            "--exclude-standard",
            "-z",
            "--",
            "build/benchmark-results",
        ],
        check=True,
        stdout=subprocess.PIPE,
    )
    paths = []
    for encoded_path in completed.stdout.split(b"\0"):
        if not encoded_path:
            continue
        path = root / os.fsdecode(encoded_path)
        if is_result_artifact(path) and path.is_file():
            paths.append(path)
    return sorted(paths)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=repo_root())
    parser.add_argument("--apply", action="store_true")
    parser.add_argument(
        "--keep",
        action="append",
        default=[],
        metavar="RESULT",
        help="preserve one exact result path and all native-helper companions",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    keep = set()
    for item in args.keep:
        keep.update(kept_result_family((root / item).resolve()))
    files = [
        path
        for path in ignored_result_artifacts(root)
        if path.resolve() not in keep
    ]

    action = "removing" if args.apply else "would remove"
    print(f"{action} {len(files)} ignored benchmark scratch artifact(s)")
    for path in files:
        print(path.relative_to(root))
        if args.apply:
            path.unlink()
    if not args.apply:
        print("dry run only; pass --apply to delete these files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
