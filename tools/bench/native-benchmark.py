#!/usr/bin/env python3
"""Run, validate, and compare Silex-native Google Benchmark rows."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import signal
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any


SCHEMA_VERSION = 1
MANIFEST_KIND = "silex-native-benchmark-executables"
RUN_METADATA_KIND = "silex-native-benchmark-run"
COUNTER_NAME_PATTERN = re.compile(r"[A-Za-z_][A-Za-z0-9_.-]*")
DIGEST_PATTERN = re.compile(r"[0-9A-Fa-f]{64}")
TIME_UNITS = {"ns", "us", "ms", "s"}
# JSON inputs are evidence artifacts, not an unbounded data transport.  The
# limit is deliberately generous for the complete native benchmark inventory.
MAX_JSON_BYTES = 64 * 1024 * 1024
# Each subprocess stream is retained in full only up to this byte count.  The
# readers continue draining while termination begins, so a noisy child cannot
# deadlock the workflow or consume unbounded memory/disk.
MAX_PROCESS_LOG_BYTES = 16 * 1024 * 1024
MAX_CONSOLE_LOG_BYTES = 64 * 1024
PROCESS_IO_CHUNK_BYTES = 64 * 1024
PROCESS_TERMINATION_GRACE_SECONDS = 0.25
REQUIRED_BOOLEAN_OPTIONS = {
    "logging": False,
    "debug_checks": False,
    "profiling": False,
    "sanitizers": False,
    "frame_pointers": False,
    "exceptions": False,
    "rtti": False,
    "shared_libraries": True,
    "position_independent_code": True,
    "unity_build": True,
    "fplll": False,
    "flatter": False,
}
KNOWN_OPTIONAL_DEPENDENCIES = {"fplll", "flatter"}
REQUIRED_STATUS_COUNTERS = (
    ("success", 1),
    ("failure_reason", 0),
)
OWNED_BENCHMARK_FLAGS = (
    "--benchmark_filter",
    "--benchmark_min_time",
    "--benchmark_out",
    "--benchmark_out_format",
    "--benchmark_repetitions",
    "--benchmark_report_aggregates_only",
)
COMPARE_TOOL_CANDIDATES = (
    Path("/usr/local/share/googlebenchmark/tools/compare.py"),
    Path("/usr/share/googlebenchmark/tools/compare.py"),
    Path("/usr/local/share/benchmark/tools/compare.py"),
    Path("/usr/share/benchmark/tools/compare.py"),
)


class ToolError(RuntimeError):
    """A user-facing workflow error."""


def _file_identity(stat_result: os.stat_result) -> tuple[int, int, int, int]:
    return (
        stat_result.st_dev,
        stat_result.st_ino,
        stat_result.st_size,
        stat_result.st_mtime_ns,
    )


def open_regular_evidence(path: Path) -> tuple[Any, os.stat_result]:
    """Open one evidence path without blocking on a FIFO or other special file."""
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NONBLOCK", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise ToolError(f"could not open evidence file {path}: {exc}") from exc
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ToolError(f"evidence path is not a regular file: {path}")
        return os.fdopen(descriptor, "rb"), before
    except BaseException:
        os.close(descriptor)
        raise


def read_bytes_snapshot(path: Path, *, limit: int) -> tuple[bytes, os.stat_result]:
    """Read one stable regular-file snapshot without exceeding ``limit``."""
    try:
        stream, before = open_regular_evidence(path)
        with stream:
            if before.st_size > limit:
                raise ToolError(
                    f"file exceeds the {limit}-byte evidence limit: {path}"
                )
            contents = stream.read(limit + 1)
            after = os.fstat(stream.fileno())
    except ToolError:
        raise
    except OSError as exc:
        raise ToolError(f"could not read file {path}: {exc}") from exc
    if len(contents) > limit:
        raise ToolError(f"file exceeds the {limit}-byte evidence limit: {path}")
    if (
        _file_identity(before) != _file_identity(after)
        or len(contents) != after.st_size
    ):
        raise ToolError(f"file changed while it was read: {path}")
    return contents, after


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_json(path: Path) -> object:
    payload, _contents, _snapshot = read_json_snapshot(path)
    return payload


def read_json_snapshot(
    path: Path,
) -> tuple[object, bytes, dict[str, object]]:
    """Read and identify one stable byte snapshot of a JSON file."""
    try:
        contents, after = read_bytes_snapshot(path, limit=MAX_JSON_BYTES)
        payload = json.loads(contents)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ToolError(f"invalid JSON in {path}: {exc}") from exc
    return (
        payload,
        contents,
        {
            "path": str(path),
            "sha256": hashlib.sha256(contents).hexdigest(),
            "size": after.st_size,
            "mtime_ns": after.st_mtime_ns,
        },
    )


def publish_bytes(path: Path, contents: bytes, *, overwrite: bool) -> None:
    """Publish bytes exclusively, or atomically replace after explicit consent."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not overwrite:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        flags |= getattr(os, "O_CLOEXEC", 0)
        try:
            descriptor = os.open(path, flags, 0o666)
        except FileExistsError as exc:
            raise ToolError(
                f"refusing to overwrite existing artifact: {path}; "
                "pass --overwrite to replace it"
            ) from exc
        created_identity: tuple[int, int] | None = None
        try:
            opened = os.fstat(descriptor)
            created_identity = (opened.st_dev, opened.st_ino)
            with os.fdopen(descriptor, "wb") as stream:
                descriptor = -1
                stream.write(contents)
                stream.flush()
                os.fsync(stream.fileno())
        except BaseException:
            if descriptor >= 0:
                os.close(descriptor)
            try:
                current = path.lstat()
                if created_identity == (current.st_dev, current.st_ino):
                    path.unlink()
            except OSError:
                pass
            raise
        return

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(contents)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def write_json(path: Path, payload: object, *, overwrite: bool) -> None:
    contents = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    publish_bytes(path, contents, overwrite=overwrite)


def hashed_file_snapshot(path: Path) -> dict[str, object]:
    """Hash one regular file through the same descriptor used for identity."""
    digest = hashlib.sha256()
    try:
        stream, before = open_regular_evidence(path)
        with stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
            after = os.fstat(stream.fileno())
    except ToolError:
        raise
    except OSError as exc:
        raise ToolError(f"could not hash evidence file {path}: {exc}") from exc
    if _file_identity(before) != _file_identity(after):
        raise ToolError(f"file changed while its identity was captured: {path}")
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
    }


def file_sha256(path: Path) -> str:
    return str(hashed_file_snapshot(path)["sha256"])


def file_snapshot(path: Path) -> dict[str, object]:
    return hashed_file_snapshot(path)


def companion_path(path: Path, suffix: str) -> Path:
    stem = path.name[:-5] if path.name.endswith(".json") else path.name
    return path.with_name(stem + suffix)


def helper_invocation() -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        *sys.argv[1:],
    ]


def is_number(value: object) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError):
        return False


def is_nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def values_equal(left: object, right: object) -> bool:
    if is_number(left) and is_number(right):
        return float(left) == float(right)
    return type(left) is type(right) and left == right


def parse_counter_requirement(text: str) -> tuple[str, object]:
    name, separator, encoded = text.partition("=")
    if separator == "" or COUNTER_NAME_PATTERN.fullmatch(name) is None:
        raise argparse.ArgumentTypeError(
            "counter requirements must use NAME=JSON_VALUE"
        )
    value = parse_scalar(encoded)
    return name, value


def parse_scalar(encoded: str) -> object:
    try:
        value = json.loads(encoded)
    except json.JSONDecodeError:
        value = encoded
    if isinstance(value, (dict, list)):
        raise argparse.ArgumentTypeError(
            "counter requirements must use a scalar JSON value"
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise argparse.ArgumentTypeError("counter values must be finite")
    return value


def parse_row_counter_requirement(text: str) -> tuple[str, str, object]:
    row_and_counter, separator, encoded = text.partition("=")
    row, row_separator, counter = row_and_counter.rpartition(":")
    if (
        separator == ""
        or row_separator == ""
        or not row
        or COUNTER_NAME_PATTERN.fullmatch(counter) is None
    ):
        raise argparse.ArgumentTypeError(
            "row counter requirements must use ROW:NAME=JSON_VALUE"
        )
    return row, counter, parse_scalar(encoded)


def parse_counter_name(text: str) -> str:
    if COUNTER_NAME_PATTERN.fullmatch(text) is None:
        raise argparse.ArgumentTypeError(
            "counter names must match [A-Za-z_][A-Za-z0-9_.-]*"
        )
    return text


def normalized_counter_requirements(
    requirements: list[tuple[str, object]],
    *,
    mandatory: tuple[tuple[str, object], ...] = (),
) -> list[tuple[str, object]]:
    """Merge counter requirements and reject contradictory duplicates."""
    merged: list[tuple[str, object]] = []
    positions: dict[str, int] = {}
    for name, expected in [*mandatory, *requirements]:
        position = positions.get(name)
        if position is None:
            positions[name] = len(merged)
            merged.append((name, expected))
            continue
        previous = merged[position][1]
        if not values_equal(previous, expected):
            raise ToolError(
                f"counter {name} has conflicting requirements "
                f"{previous!r} and {expected!r}"
            )
    return merged


def canonical_row_name(row: dict[str, object]) -> str | None:
    run_name = row.get("run_name")
    if isinstance(run_name, str) and run_name:
        return run_name
    name = row.get("name")
    if isinstance(name, str) and name:
        return name
    return None


def validate_payload(
    payload: object,
    *,
    expected_rows: list[str] | None = None,
    expected_repetitions: int | None = None,
    required_counters: list[tuple[str, object]] | None = None,
    required_row_counters: list[tuple[str, str, object]] | None = None,
) -> dict[str, object]:
    """Validate raw Google Benchmark JSON and Silex semantic conventions."""

    errors: list[str] = []
    if not isinstance(payload, dict):
        return {
            "ok": False,
            "errors": ["Google Benchmark output is not a JSON object"],
            "benchmark_names": [],
            "sample_count": 0,
        }

    context = payload.get("context")
    if not isinstance(context, dict):
        errors.append("Google Benchmark JSON has no context object")
    else:
        schema_version = context.get("json_schema_version")
        if type(schema_version) is not int or schema_version != 1:
            errors.append(
                "unsupported or missing Google Benchmark "
                f"json_schema_version: {schema_version!r}"
            )
    rows = payload.get("benchmarks")
    if not isinstance(rows, list):
        return {
            "ok": False,
            "errors": errors
            + ["Google Benchmark JSON has no benchmarks array"],
            "benchmark_names": [],
            "sample_count": 0,
        }
    if not rows:
        errors.append("Google Benchmark JSON contains no benchmark rows")

    requirements = required_counters or []
    row_requirements: dict[str, list[tuple[str, object]]] = {}
    for row_name, counter, expected in required_row_counters or []:
        row_requirements.setdefault(row_name, []).append((counter, expected))
    samples_by_name: dict[str, list[dict[str, object]]] = {}
    repetition_indexes: dict[str, set[int]] = {}
    for index, value in enumerate(rows):
        if not isinstance(value, dict):
            errors.append(f"benchmark row {index} is not an object")
            continue
        row = value
        name = canonical_row_name(row)
        display_name = name or f"row {index}"
        if name is None:
            errors.append(f"benchmark row {index} has no nonempty name")

        if "error_occurred" in row:
            error_occurred = row["error_occurred"]
            if error_occurred is True:
                message = row.get("error_message", "<no message>")
                errors.append(f"{display_name}: benchmark error: {message}")
            elif error_occurred is not False:
                errors.append(
                    f"{display_name}: error_occurred must be a JSON boolean"
                )

        run_type = row.get("run_type", "iteration")
        if run_type == "aggregate":
            continue
        if run_type != "iteration":
            errors.append(f"{display_name}: unsupported run_type {run_type!r}")
            continue
        if name is None:
            continue
        raw_name = row.get("name")
        run_name = row.get("run_name")
        if (
            isinstance(raw_name, str)
            and isinstance(run_name, str)
            and raw_name != run_name
        ):
            errors.append(
                f"{display_name}: iteration name and run_name differ"
            )

        samples_by_name.setdefault(name, []).append(row)
        for field in ("real_time", "cpu_time"):
            measured = row.get(field)
            if not is_number(measured) or float(measured) < 0.0:
                errors.append(
                    f"{display_name}: {field} must be a finite nonnegative number"
                )
        iterations = row.get("iterations")
        if (
            not is_number(iterations)
            or float(iterations) <= 0.0
            or float(iterations) != int(float(iterations))
        ):
            errors.append(
                f"{display_name}: iterations must be a positive integer"
            )
        if row.get("time_unit") not in TIME_UNITS:
            errors.append(
                f"{display_name}: unsupported time_unit {row.get('time_unit')!r}"
            )

        repetition_index = row.get("repetition_index")
        if repetition_index is not None:
            if (
                not is_number(repetition_index)
                or float(repetition_index) < 0.0
                or float(repetition_index) != int(float(repetition_index))
            ):
                errors.append(
                    f"{display_name}: repetition_index must be a nonnegative integer"
                )
            else:
                numeric_index = int(float(repetition_index))
                seen = repetition_indexes.setdefault(name, set())
                if numeric_index in seen:
                    errors.append(
                        f"{display_name}: duplicate repetition_index {numeric_index}"
                    )
                seen.add(numeric_index)

        # Silex benchmark rows use these names as semantic status fields.
        # Generic validation rejects failed values whenever they are present;
        # run, gate, and compare additionally require both fields.
        if "success" in row and not values_equal(row["success"], 1):
            errors.append(
                f"{display_name}: success={row['success']!r}, expected 1"
            )
        if "failure_reason" in row and not values_equal(
            row["failure_reason"], 0
        ):
            errors.append(
                f"{display_name}: failure_reason={row['failure_reason']!r}, expected 0"
            )
        for counter, expected in requirements:
            if counter not in row:
                errors.append(f"{display_name}: missing counter {counter}")
            elif not values_equal(row[counter], expected):
                errors.append(
                    f"{display_name}: {counter}={row[counter]!r}, "
                    f"expected {expected!r}"
                )
        for counter, expected in row_requirements.get(name, []):
            if counter not in row:
                errors.append(f"{display_name}: missing counter {counter}")
            elif not values_equal(row[counter], expected):
                errors.append(
                    f"{display_name}: {counter}={row[counter]!r}, "
                    f"expected {expected!r}"
                )

    actual_names = sorted(samples_by_name)
    if rows and not samples_by_name:
        errors.append(
            "Google Benchmark JSON has no raw iteration rows; "
            "do not use aggregates-only output"
        )
    for name, samples in samples_by_name.items():
        has_repetition_indexes = [
            "repetition_index" in row for row in samples
        ]
        if len(samples) > 1 and not all(has_repetition_indexes):
            errors.append(
                f"{name}: every repeated raw sample must have a "
                "repetition_index"
            )
        if any(has_repetition_indexes):
            actual_indexes = repetition_indexes.get(name, set())
            expected_indexes = set(range(len(samples)))
            if actual_indexes != expected_indexes:
                errors.append(
                    f"{name}: repetition indexes {sorted(actual_indexes)!r}, "
                    f"expected {sorted(expected_indexes)!r}"
                )
        has_repetitions = ["repetitions" in row for row in samples]
        if any(has_repetitions):
            if not all(has_repetitions) or any(
                not values_equal(row.get("repetitions"), len(samples))
                for row in samples
            ):
                errors.append(
                    f"{name}: repetitions metadata does not match "
                    f"the {len(samples)} raw samples"
                )
    if expected_rows is not None:
        expected_names = set(expected_rows)
        actual_name_set = set(actual_names)
        missing = sorted(expected_names - actual_name_set)
        unexpected = sorted(actual_name_set - expected_names)
        if missing:
            errors.append("missing expected benchmark rows: " + ", ".join(missing))
        if unexpected:
            errors.append("unexpected benchmark rows: " + ", ".join(unexpected))
        if len(expected_names) != len(expected_rows):
            errors.append("expected benchmark row list contains duplicates")
    missing_requirement_rows = sorted(set(row_requirements) - set(actual_names))
    if missing_requirement_rows:
        errors.append(
            "rows required by counter rules are missing: "
            + ", ".join(missing_requirement_rows)
        )

    if expected_repetitions is not None:
        for name, samples in samples_by_name.items():
            if len(samples) != expected_repetitions:
                errors.append(
                    f"{name}: found {len(samples)} raw samples, "
                    f"expected {expected_repetitions}"
                )

    return {
        "ok": not errors,
        "errors": errors,
        "benchmark_names": actual_names,
        "sample_count": sum(len(samples) for samples in samples_by_name.values()),
        "samples_per_benchmark": {
            name: len(samples) for name, samples in sorted(samples_by_name.items())
        },
    }


def load_manifest(path: Path) -> dict[str, object]:
    payload = read_json(path)
    if not isinstance(payload, dict):
        raise ToolError(f"benchmark manifest {path} is not a JSON object")
    schema_version = payload.get("schema_version")
    if type(schema_version) is not int or schema_version != SCHEMA_VERSION:
        raise ToolError(
            f"benchmark manifest {path} has unsupported schema_version "
            f"{schema_version!r}"
        )
    if payload.get("kind") != MANIFEST_KIND:
        raise ToolError(
            f"benchmark manifest {path} has invalid kind "
            f"{payload.get('kind')!r}"
        )
    targets = payload.get("targets")
    if not isinstance(targets, list) or not all(
        isinstance(row, dict) for row in targets
    ):
        raise ToolError(f"benchmark manifest {path} has no targets array")
    return payload


def resolve_executable(
    exe_argument: Path | None,
    manifest_path: Path | None,
    target: str | None,
) -> tuple[Path, str | None, dict[str, object] | None]:
    if exe_argument is not None and target is not None:
        raise ToolError("use either --exe or --target, not both")
    if exe_argument is not None and manifest_path is not None:
        raise ToolError("--manifest may only be used with --target")
    if exe_argument is None and target is None:
        raise ToolError("one of --exe or --target is required")

    manifest: dict[str, object] | None = None
    resolved_target = target
    if target is not None:
        if manifest_path is None:
            raise ToolError("--target requires --manifest")
        manifest = load_manifest(manifest_path)
        matches = [
            row
            for row in manifest["targets"]  # type: ignore[index]
            if isinstance(row, dict) and row.get("target") == target
        ]
        if not matches:
            raise ToolError(f"unknown benchmark target in manifest: {target}")
        if len(matches) != 1 or not isinstance(matches[0].get("path"), str):
            raise ToolError(f"invalid manifest entry for benchmark target: {target}")
        executable = Path(str(matches[0]["path"]))
    else:
        executable = Path(exe_argument)  # type: ignore[arg-type]

    executable = executable.expanduser().resolve()
    if not executable.is_file():
        raise ToolError(f"benchmark executable does not exist: {executable}")
    if os.name != "nt" and not os.access(executable, os.X_OK):
        raise ToolError(f"benchmark executable is not executable: {executable}")
    return executable, resolved_target, manifest


def git_metadata(root: Path) -> dict[str, object]:
    def git(*arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", *arguments],
            cwd=root,
            check=False,
            text=True,
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    try:
        revision = git("rev-parse", "HEAD")
        status = git("status", "--short", "--untracked-files=all")
    except OSError:
        return {"available": False}
    if revision.returncode != 0 or status.returncode != 0:
        return {"available": False}
    status_lines = status.stdout.splitlines()
    return {
        "available": True,
        "commit": revision.stdout.strip(),
        "dirty": bool(status_lines),
        "status": status_lines,
    }


def environment_metadata() -> dict[str, object]:
    affinity: list[int] | None = None
    if hasattr(os, "sched_getaffinity"):
        try:
            affinity = sorted(os.sched_getaffinity(0))  # type: ignore[attr-defined]
        except OSError:
            affinity = None
    processor = platform.processor().strip()
    cpu_model = processor
    if not cpu_model:
        try:
            for line in Path("/proc/cpuinfo").read_text(errors="replace").splitlines():
                name, separator, value = line.partition(":")
                if separator and name.strip() in {"model name", "Hardware"}:
                    cpu_model = value.strip()
                    if cpu_model:
                        break
        except OSError:
            pass
    if not cpu_model:
        cpu_model = platform.machine()
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": processor,
        "cpu_model": cpu_model,
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "cpu_affinity": affinity,
    }


def source_root_from_manifest(
    manifest: dict[str, object] | None,
) -> Path:
    if manifest is None:
        return repo_root()
    value = manifest.get("source_dir")
    if not isinstance(value, str) or not value:
        raise ToolError("benchmark manifest has no source_dir")
    candidate = Path(value).expanduser().resolve()
    if not candidate.is_dir():
        raise ToolError(
            f"benchmark manifest source_dir does not exist: {candidate}"
        )
    return candidate


def benchmark_arguments(args: argparse.Namespace, output: Path) -> list[str]:
    extras = list(args.benchmark_arg or [])
    for argument in extras:
        if any(
            argument == flag or argument.startswith(flag + "=")
            for flag in OWNED_BENCHMARK_FLAGS
        ):
            raise ToolError(
                f"{argument} is managed by the native benchmark helper"
            )
    return [
        f"--benchmark_filter={args.filter}",
        f"--benchmark_min_time={args.min_time}",
        f"--benchmark_repetitions={args.repetitions}",
        "--benchmark_report_aggregates_only=false",
        f"--benchmark_out={output}",
        "--benchmark_out_format=json",
        *extras,
    ]


def prepare_outputs(paths: list[Path], overwrite: bool) -> None:
    if len(set(paths)) != len(paths):
        raise ToolError("artifact output paths must be distinct")
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        raise ToolError(
            "refusing to overwrite existing artifact(s): "
            + ", ".join(str(path) for path in existing)
            + "; pass --overwrite to replace them"
        )
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def prepare_stream_outputs(paths: list[Path], overwrite: bool) -> None:
    """Make stream paths available for the exclusive opens in the supervisor."""
    if not overwrite:
        return
    for path in paths:
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            raise ToolError(f"could not replace process log {path}: {exc}") from exc


def discard_stale_output(path: Path, *, overwrite: bool) -> None:
    """Remove an old artifact when an explicit overwrite produced no replacement."""
    if not overwrite:
        return
    try:
        path.unlink()
    except FileNotFoundError:
        pass
    except OSError as exc:
        raise ToolError(f"could not discard stale artifact {path}: {exc}") from exc


def reject_output_input_overlap(
    outputs: list[Path], inputs: list[Path]
) -> None:
    overlap = set(outputs) & {path.resolve() for path in inputs}
    if overlap:
        raise ToolError(
            "artifact output path would overwrite an input: "
            + ", ".join(str(path) for path in sorted(overlap))
        )


def read_bounded_log_text(path: Path) -> str:
    """Read a small head/tail view without loading an entire retained log."""
    try:
        retained, _after = read_bytes_snapshot(
            path, limit=MAX_PROCESS_LOG_BYTES
        )
        size = len(retained)
        if size <= MAX_CONSOLE_LOG_BYTES:
            contents = retained
        else:
            half = MAX_CONSOLE_LOG_BYTES // 2
            contents = (
                retained[:half]
                + b"\n... bounded console view; see full capped artifact ...\n"
                + retained[-half:]
            )
    except ToolError as exc:
        raise ToolError(f"could not read process log {path}: {exc}") from exc
    return contents.decode(errors="replace")


def terminate_process_tree(
    process: subprocess.Popen[bytes], *, posix_group: bool
) -> None:
    """Terminate a supervised process group, with a direct-process fallback."""
    if posix_group:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except (OSError, ProcessLookupError):
            pass
    elif process.poll() is None:
        try:
            process.terminate()
        except OSError:
            pass

    try:
        process.wait(timeout=PROCESS_TERMINATION_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        pass

    # Send the hard signal to the group even when the leader exited during the
    # grace period: descendants may still hold a pipe open.
    if posix_group:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            pass
    if process.poll() is None:
        try:
            process.kill()
        except OSError:
            pass
    try:
        process.wait(timeout=PROCESS_TERMINATION_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        pass


def run_bounded_process(
    command: list[str],
    *,
    cwd: Path,
    timeout: float,
    stdout_path: Path,
    stderr_path: Path,
    log_limit: int = MAX_PROCESS_LOG_BYTES,
    result_path: Path | None = None,
    result_limit: int = MAX_JSON_BYTES,
) -> dict[str, object]:
    """Run a child with bounded logs, result capture, and shutdown."""
    if log_limit <= 0:
        raise ToolError("process log limit must be positive")
    if result_path is not None and result_limit <= 0:
        raise ToolError("process result limit must be positive")
    states = {
        "stdout": {"bytes_seen": 0, "truncated": False, "error": None},
        "stderr": {"bytes_seen": 0, "truncated": False, "error": None},
    }
    limit_event = threading.Event()
    error_event = threading.Event()
    result_state: dict[str, object] = {
        "bytes_observed": 0,
        "limit_exceeded": False,
        "error": None,
    }
    done_events = {"stdout": threading.Event(), "stderr": threading.Event()}
    log_handles: dict[str, Any] = {}
    process: subprocess.Popen[bytes] | None = None
    posix_group = os.name == "posix"
    started = time.monotonic()

    def inspect_result() -> None:
        if result_path is None or result_state["error"] is not None:
            return
        try:
            observed = result_path.lstat()
        except FileNotFoundError:
            return
        except OSError as exc:
            result_state["error"] = str(exc)
            return
        if not stat.S_ISREG(observed.st_mode):
            result_state["error"] = (
                f"result path is not a regular file: {result_path}"
            )
            return
        result_state["bytes_observed"] = max(
            int(result_state["bytes_observed"]), observed.st_size
        )
        if observed.st_size > result_limit:
            result_state["limit_exceeded"] = True

    def drain(name: str, stream: Any) -> None:
        state = states[name]
        destination = log_handles[name]
        try:
            while True:
                chunk = stream.read(PROCESS_IO_CHUNK_BYTES)
                if not chunk:
                    break
                state["bytes_seen"] = int(state["bytes_seen"]) + len(chunk)
                retained = destination.tell()
                remaining = max(0, log_limit - retained)
                if remaining:
                    destination.write(chunk[:remaining])
                    destination.flush()
                if len(chunk) > remaining:
                    state["truncated"] = True
                    limit_event.set()
                # Continue draining/discarding while the supervisor terminates
                # the child, preventing a full pipe from deadlocking shutdown.
        except (OSError, ValueError) as exc:
            state["error"] = str(exc)
            error_event.set()
        finally:
            try:
                stream.close()
            except (OSError, ValueError):
                pass
            done_events[name].set()

    try:
        log_handles["stdout"] = stdout_path.open("xb")
        log_handles["stderr"] = stderr_path.open("xb")
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=posix_group,
        )
        if process.stdout is None or process.stderr is None:
            raise ToolError("could not create subprocess output pipes")
        threads = [
            threading.Thread(
                target=drain,
                args=("stdout", process.stdout),
                daemon=True,
                name="native-benchmark-stdout",
            ),
            threading.Thread(
                target=drain,
                args=("stderr", process.stderr),
                daemon=True,
                name="native-benchmark-stderr",
            ),
        ]
        for thread in threads:
            thread.start()

        deadline = started + timeout
        termination_reason: str | None = None
        while True:
            inspect_result()
            if bool(result_state["limit_exceeded"]):
                termination_reason = "result_limit"
                break
            if result_state["error"] is not None:
                termination_reason = "result_error"
                break
            readers_done = all(event.is_set() for event in done_events.values())
            if process.poll() is not None and readers_done:
                break
            if limit_event.is_set():
                termination_reason = "output_limit"
                break
            if error_event.is_set():
                termination_reason = "log_error"
                break
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                termination_reason = "timeout"
                break
            limit_event.wait(timeout=min(remaining, 0.025))

        # Catch a final write made immediately before a clean child exit.
        inspect_result()
        if termination_reason is None and bool(result_state["limit_exceeded"]):
            termination_reason = "result_limit"
        if termination_reason is None and result_state["error"] is not None:
            termination_reason = "result_error"

        if termination_reason is not None:
            terminate_process_tree(process, posix_group=posix_group)
        for event in done_events.values():
            event.wait(PROCESS_TERMINATION_GRACE_SECONDS)
        if not all(event.is_set() for event in done_events.values()):
            # A detached descendant may retain a pipe on platforms without
            # process-tree control. Closing our read ends keeps this function
            # bounded; the incomplete drain invalidates the result.
            for stream in (process.stdout, process.stderr):
                try:
                    stream.close()
                except (OSError, ValueError):
                    pass
            error_event.set()
            for name, event in done_events.items():
                if not event.is_set() and states[name]["error"] is None:
                    states[name]["error"] = "process log drain did not finish"
        for thread in threads:
            thread.join(timeout=PROCESS_TERMINATION_GRACE_SECONDS)

        elapsed = time.monotonic() - started
        forced = termination_reason is not None
        result_returncode = None if forced else process.returncode
        return {
            "returncode": result_returncode,
            "observed_returncode": process.returncode,
            "elapsed_seconds": elapsed,
            "timed_out": termination_reason == "timeout",
            "output_limit_exceeded": bool(limit_event.is_set()),
            "result_limit_exceeded": bool(result_state["limit_exceeded"]),
            "result_capture_error": result_state["error"],
            "termination_reason": termination_reason,
            "log_capture_error": (
                "; ".join(
                    str(state["error"])
                    for state in states.values()
                    if state["error"] is not None
                )
                or None
            ),
            "stdout_bytes_seen": states["stdout"]["bytes_seen"],
            "stderr_bytes_seen": states["stderr"]["bytes_seen"],
            "stdout_truncated": states["stdout"]["truncated"],
            "stderr_truncated": states["stderr"]["truncated"],
            "log_byte_limit_per_stream": log_limit,
            "result_bytes_observed": result_state["bytes_observed"],
            "result_byte_limit": result_limit if result_path is not None else None,
        }
    finally:
        if process is not None and process.poll() is None:
            terminate_process_tree(process, posix_group=posix_group)
        for handle in log_handles.values():
            try:
                handle.close()
            except OSError:
                pass


def execution_metadata(
    *,
    args: argparse.Namespace,
    target: str | None,
    manifest: dict[str, object] | None,
    manifest_snapshot: dict[str, object] | None,
    executable_snapshot: dict[str, object],
    helper_snapshot: dict[str, object],
    command: list[str],
    raw_output: Path,
    raw_output_sha256: str | None,
    raw_output_published: bool,
    stdout_path: Path,
    stderr_path: Path,
    process_result: dict[str, object],
    dependency_verification: dict[str, object] | None,
    validation: dict[str, object],
    started_at_utc: str,
    repository: dict[str, object],
    environment: dict[str, object],
) -> dict[str, object]:
    source_root = source_root_from_manifest(manifest)
    build = manifest.get("build", {}) if manifest is not None else {}
    return {
        "schema_version": SCHEMA_VERSION,
        "kind": RUN_METADATA_KIND,
        "mode": args.command,
        "label": args.label,
        "notes": list(args.note or []),
        "started_at_utc": started_at_utc,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository": repository,
        "environment": environment,
        "build": build,
        "dependency_verification": dependency_verification,
        "helper": helper_snapshot,
        "manifest": manifest_snapshot,
        "executable": {"target": target, **executable_snapshot},
        "selection": {
            "filter": args.filter,
            "min_time": args.min_time,
            "repetitions": args.repetitions,
            "benchmark_args": list(args.benchmark_arg or []),
            "expected_rows": list(args.expect_row or []),
            "required_counters": {
                name: expected for name, expected in args.require_counter
            },
            "required_row_counters": [
                {"row": row, "counter": counter, "value": expected}
                for row, counter, expected in args.require_row_counter
            ],
        },
        "helper_invocation": helper_invocation(),
        "command": command,
        "working_directory": str(source_root),
        "artifacts": {
            "raw_json": str(raw_output),
            "raw_json_sha256": raw_output_sha256,
            "raw_json_published": raw_output_published,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
        },
        "process": process_result,
        "validation": validation,
    }


def run_benchmark(args: argparse.Namespace) -> int:
    args.require_counter = normalized_counter_requirements(
        list(args.require_counter),
        mandatory=REQUIRED_STATUS_COUNTERS,
    )
    output = args.output.expanduser().resolve()
    metadata_path = companion_path(output, ".metadata.json")
    stdout_path = companion_path(output, ".stdout.txt")
    stderr_path = companion_path(output, ".stderr.txt")
    if output == metadata_path:
        raise ToolError("raw JSON and metadata paths must be distinct")
    executable, target, manifest = resolve_executable(
        args.exe, args.manifest, args.target
    )
    source_root = source_root_from_manifest(manifest)
    started_at_utc = datetime.now(timezone.utc).isoformat()
    repository = git_metadata(source_root)
    environment = environment_metadata()
    output_paths = [output, metadata_path, stdout_path, stderr_path]
    input_paths = [executable, Path(__file__).resolve()]
    if args.manifest is not None:
        input_paths.append(args.manifest.expanduser().resolve())
    reject_output_input_overlap(output_paths, input_paths)
    executable_snapshot = file_snapshot(executable)
    helper_snapshot = file_snapshot(Path(__file__).resolve())
    manifest_snapshot = (
        file_snapshot(args.manifest.expanduser().resolve())
        if args.manifest is not None
        else None
    )
    dependency_verification = (
        verify_manifest_dependency_files(manifest)
        if manifest is not None
        else None
    )
    prepare_outputs(output_paths, args.overwrite)
    prepare_stream_outputs([stdout_path, stderr_path], args.overwrite)
    raw_bytes: bytes | None = None
    with tempfile.TemporaryDirectory(prefix="silex-native-run-") as directory:
        private_output = Path(directory) / "benchmark.json"
        command = [
            str(executable),
            *benchmark_arguments(args, private_output),
        ]
        process_result = run_bounded_process(
            command,
            cwd=source_root,
            timeout=args.timeout,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            result_path=private_output,
            result_limit=MAX_JSON_BYTES,
        )

        validation: dict[str, object]
        if process_result["result_limit_exceeded"]:
            validation = {
                "ok": False,
                "errors": [
                    "benchmark result exceeded the JSON evidence byte limit; "
                    "the private result was not published"
                ],
                "benchmark_names": [],
                "sample_count": 0,
            }
        elif process_result["result_capture_error"] is not None:
            validation = {
                "ok": False,
                "errors": [
                    "benchmark result capture failed: "
                    + str(process_result["result_capture_error"])
                ],
                "benchmark_names": [],
                "sample_count": 0,
            }
        else:
            try:
                raw_bytes, _raw_after = read_bytes_snapshot(
                    private_output, limit=MAX_JSON_BYTES
                )
            except ToolError as exc:
                validation = {
                    "ok": False,
                    "errors": [str(exc)],
                    "benchmark_names": [],
                    "sample_count": 0,
                }
            else:
                try:
                    payload = json.loads(raw_bytes)
                except (json.JSONDecodeError, UnicodeDecodeError) as exc:
                    validation = {
                        "ok": False,
                        "errors": [f"invalid benchmark JSON: {exc}"],
                        "benchmark_names": [],
                        "sample_count": 0,
                    }
                else:
                    validation = validate_payload(
                        payload,
                        expected_rows=(
                            list(args.expect_row)
                            if args.expect_row is not None
                            else None
                        ),
                        expected_repetitions=args.repetitions,
                        required_counters=list(args.require_counter),
                        required_row_counters=list(args.require_row_counter),
                    )

    if process_result["timed_out"]:
        validation["errors"] = list(validation["errors"]) + [
            f"benchmark timed out after {args.timeout:g} seconds"
        ]
        validation["ok"] = False
    if process_result["output_limit_exceeded"]:
        validation["errors"] = list(validation["errors"]) + [
            "benchmark exceeded the per-stream process log byte limit"
        ]
        validation["ok"] = False
    if process_result["result_limit_exceeded"]:
        validation["ok"] = False
    if process_result["result_capture_error"] is not None:
        validation["ok"] = False
    if process_result["log_capture_error"] is not None:
        validation["errors"] = list(validation["errors"]) + [
            "benchmark log capture failed: "
            + str(process_result["log_capture_error"])
        ]
        validation["ok"] = False
    returncode = process_result["returncode"]
    if returncode != 0:
        validation["errors"] = list(validation["errors"]) + [
            (
                f"benchmark executable exited with {returncode}"
                if returncode is not None
                else "benchmark executable did not complete successfully"
            )
        ]
        validation["ok"] = False

    if raw_bytes is not None and validation["ok"]:
        publish_bytes(output, raw_bytes, overwrite=args.overwrite)
        raw_output_sha256 = hashlib.sha256(raw_bytes).hexdigest()
        raw_output_published = True
    else:
        discard_stale_output(output, overwrite=args.overwrite)
        raw_output_sha256 = None
        raw_output_published = False

    metadata = execution_metadata(
        args=args,
        target=target,
        manifest=manifest,
        manifest_snapshot=manifest_snapshot,
        executable_snapshot=executable_snapshot,
        helper_snapshot=helper_snapshot,
        command=command,
        raw_output=output,
        raw_output_sha256=raw_output_sha256,
        raw_output_published=raw_output_published,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        process_result=process_result,
        dependency_verification=dependency_verification,
        validation=validation,
        started_at_utc=started_at_utc,
        repository=repository,
        environment=environment,
    )
    write_json(metadata_path, metadata, overwrite=args.overwrite)

    if not validation["ok"]:
        print(f"native benchmark {args.command} failed:", file=sys.stderr)
        for error in validation["errors"]:
            print(f"- {error}", file=sys.stderr)
        print(f"metadata: {metadata_path}", file=sys.stderr)
        return 1
    count = validation["sample_count"]
    names = len(validation["benchmark_names"])
    print(
        f"validated {count} raw sample(s) across {names} benchmark row(s); "
        f"metadata: {metadata_path}"
    )
    return 0


def samples_by_name(payload: object) -> dict[str, list[dict[str, object]]]:
    if not isinstance(payload, dict):
        raise ToolError("validated benchmark payload is not an object")
    rows = payload.get("benchmarks")
    if not isinstance(rows, list):
        raise ToolError("validated benchmark payload has no benchmarks array")
    result: dict[str, list[dict[str, object]]] = {}
    for value in rows:
        if not isinstance(value, dict):
            raise ToolError("validated benchmark payload has a non-object row")
        if value.get("run_type", "iteration") != "iteration":
            continue
        name = canonical_row_name(value)
        if name is None:
            raise ToolError("validated benchmark row has no name")
        result.setdefault(name, []).append(value)
    return result


def compare_semantics(
    baseline: object,
    candidate: object,
    counter_names: list[str],
) -> list[str]:
    failures: list[str] = []
    baseline_rows = samples_by_name(baseline)
    candidate_rows = samples_by_name(candidate)
    baseline_names = set(baseline_rows)
    candidate_names = set(candidate_rows)
    if baseline_names != candidate_names:
        missing = sorted(baseline_names - candidate_names)
        added = sorted(candidate_names - baseline_names)
        if missing:
            failures.append("candidate is missing rows: " + ", ".join(missing))
        if added:
            failures.append("candidate has extra rows: " + ", ".join(added))

    for name in sorted(baseline_names & candidate_names):
        left_rows = baseline_rows[name]
        right_rows = candidate_rows[name]
        if len(left_rows) != len(right_rows):
            failures.append(
                f"{name}: sample counts differ: "
                f"baseline={len(left_rows)}, candidate={len(right_rows)}"
            )
        left_units = {row.get("time_unit") for row in left_rows}
        right_units = {row.get("time_unit") for row in right_rows}
        if left_units != right_units:
            failures.append(
                f"{name}: time units differ: baseline={sorted(left_units, key=str)!r}, "
                f"candidate={sorted(right_units, key=str)!r}"
            )
        for counter in counter_names:
            left_values = [row.get(counter, None) for row in left_rows]
            right_values = [row.get(counter, None) for row in right_rows]
            if any(counter not in row for row in left_rows):
                failures.append(f"{name}: baseline is missing counter {counter}")
                continue
            if any(counter not in row for row in right_rows):
                failures.append(f"{name}: candidate is missing counter {counter}")
                continue
            if any(
                not values_equal(left_values[0], value)
                for value in left_values[1:]
            ):
                failures.append(
                    f"{name}: baseline counter {counter} changes across repetitions"
                )
                continue
            if any(
                not values_equal(right_values[0], value)
                for value in right_values[1:]
            ):
                failures.append(
                    f"{name}: candidate counter {counter} changes across repetitions"
                )
                continue
            if not values_equal(left_values[0], right_values[0]):
                failures.append(
                    f"{name}: counter {counter} differs: "
                    f"baseline={left_values[0]!r}, candidate={right_values[0]!r}"
                )
    return failures


def compare_tool_candidates(argument: Path | None) -> list[Path]:
    if argument is not None:
        return [argument]
    env_value = os.environ.get("SILEX_GOOGLE_BENCHMARK_COMPARE")
    candidates = ([Path(env_value)] if env_value else []) + list(
        COMPARE_TOOL_CANDIDATES
    )
    executable = shutil.which("compare.py")
    if executable:
        candidates.append(Path(executable))
    return candidates


def discover_compare_tool(argument: Path | None) -> Path:
    candidates = compare_tool_candidates(argument)
    for candidate in candidates:
        path = candidate.expanduser().resolve()
        if path.is_file():
            return path
    raise ToolError(
        "could not find Google Benchmark compare.py; pass --compare-tool "
        "or set SILEX_GOOGLE_BENCHMARK_COMPARE"
    )


def comparison_artifact_ref(
    raw_snapshot: dict[str, object],
    metadata_snapshot: dict[str, object] | None,
) -> dict[str, object]:
    return {
        **raw_snapshot,
        "metadata": (
            metadata_snapshot.get("path")
            if metadata_snapshot is not None
            else None
        ),
        "metadata_sha256": (
            metadata_snapshot.get("sha256")
            if metadata_snapshot is not None
            else None
        ),
    }


def load_run_metadata(
    raw_path: Path,
    raw_sha256: str,
) -> tuple[
    dict[str, object] | None,
    list[str],
    dict[str, object] | None,
]:
    metadata_path = companion_path(raw_path, ".metadata.json")
    if not metadata_path.is_file():
        return None, [f"missing run metadata: {metadata_path}"], None
    try:
        payload, _contents, metadata_snapshot = read_json_snapshot(metadata_path)
    except ToolError as exc:
        return None, [str(exc)], None
    if not isinstance(payload, dict):
        return (
            None,
            [f"run metadata is not an object: {metadata_path}"],
            metadata_snapshot,
        )

    errors: list[str] = []
    schema_version = payload.get("schema_version")
    if type(schema_version) is not int or schema_version != SCHEMA_VERSION:
        errors.append(
            f"unsupported run metadata schema_version: "
            f"{schema_version!r}"
        )
    if payload.get("kind") != RUN_METADATA_KIND:
        errors.append(f"invalid run metadata kind: {payload.get('kind')!r}")
    if payload.get("mode") != "run":
        errors.append(
            f"timing comparison requires run metadata, got mode "
            f"{payload.get('mode')!r}"
        )
    if not isinstance(payload.get("label"), str) or not str(
        payload.get("label", "")
    ).strip():
        errors.append("run metadata has no nonempty label")
    for field in ("helper", "executable", "manifest"):
        snapshot = payload.get(field)
        if not isinstance(snapshot, dict) or not is_nonempty_string(
            snapshot.get("sha256")
        ):
            errors.append(f"run metadata has no {field} file digest")

    validation = payload.get("validation")
    if not isinstance(validation, dict) or validation.get("ok") is not True:
        errors.append("run metadata does not record successful validation")
    process = payload.get("process")
    if (
        not isinstance(process, dict)
        or process.get("returncode") != 0
        or process.get("timed_out") is not False
        or process.get("output_limit_exceeded") is not False
        or process.get("result_limit_exceeded") is not False
        or process.get("result_capture_error") is not None
        or process.get("log_capture_error") is not None
    ):
        errors.append("run metadata does not record a successful process")
    artifacts = payload.get("artifacts")
    recorded_digest = (
        artifacts.get("raw_json_sha256")
        if isinstance(artifacts, dict)
        else None
    )
    if not isinstance(artifacts, dict) or artifacts.get(
        "raw_json_published"
    ) is not True:
        errors.append("run metadata does not record a published raw result")
    if recorded_digest != raw_sha256:
        errors.append(
            "run metadata raw_json_sha256 does not match the comparison input"
        )

    dependency_verification = payload.get("dependency_verification")
    verification_files = (
        dependency_verification.get("files")
        if isinstance(dependency_verification, dict)
        else None
    )
    if not isinstance(verification_files, list) or not verification_files:
        errors.append("run metadata has no dependency verification records")
    else:
        recorded_entries: set[tuple[str, str, str]] = set()
        for entry in verification_files:
            if not isinstance(entry, dict):
                errors.append("run metadata has an invalid dependency verification")
                continue
            dependency_name = entry.get("dependency")
            path_text = entry.get("path")
            expected_digest = entry.get("expected_sha256")
            actual_digest = entry.get("sha256")
            if (
                not isinstance(dependency_name, str)
                or not isinstance(path_text, str)
                or not is_absolute_evidence_path(path_text)
                or not isinstance(expected_digest, str)
                or DIGEST_PATTERN.fullmatch(expected_digest) is None
                or actual_digest != expected_digest
            ):
                errors.append(
                    "run metadata has an invalid dependency verification record"
                )
                continue
            recorded_entries.add(
                (dependency_name, path_text, expected_digest.lower())
            )
        build = payload.get("build")
        if isinstance(build, dict):
            try:
                expected_entries = set(
                    manifest_dependency_digest_entries({"build": build})
                )
                if recorded_entries != expected_entries:
                    errors.append(
                        "run metadata dependency verification inventory differs "
                        "from its build manifest"
                    )
            except ToolError as exc:
                errors.append(str(exc))

    repository = payload.get("repository")
    if not isinstance(repository, dict) or repository.get("available") is not True:
        errors.append("run metadata has no Git provenance")
    elif (
        not is_nonempty_string(repository.get("commit"))
        or not isinstance(repository.get("dirty"), bool)
        or not isinstance(repository.get("status"), list)
    ):
        errors.append("run metadata has incomplete Git provenance")
    return payload, errors, metadata_snapshot


def metadata_build_for_comparison(metadata: dict[str, object]) -> object:
    build = metadata.get("build")
    if not isinstance(build, dict):
        return build
    return {
        key: value
        for key, value in build.items()
        if key != "google_benchmark_compare"
    }


def parse_cmake_boolean(
    value: object, *, field: str, allow_unset: bool = False
) -> tuple[bool | None, str | None]:
    """Parse a serialized CMake boolean without treating unknown text as false."""
    if not isinstance(value, str):
        return None, f"{field} is not a CMake boolean string"
    normalized = value.strip().upper()
    if allow_unset and normalized == "":
        return None, None
    if normalized in {"1", "ON", "YES", "TRUE", "Y"}:
        return True, None
    if normalized in {"0", "OFF", "NO", "FALSE", "N"}:
        return False, None
    return None, f"{field} has invalid CMake boolean {value!r}"


def is_absolute_evidence_path(value: str) -> bool:
    return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()


def parse_library_digest_list(
    value: object, *, field: str
) -> list[tuple[str, str]]:
    """Parse the manifest's semicolon-delimited absolute PATH=SHA256 list."""
    if not is_nonempty_string(value):
        raise ToolError(f"{field} is empty")
    entries: list[tuple[str, str]] = []
    seen: set[str] = set()
    for encoded in str(value).split(";"):
        path_text, separator, digest = encoded.rpartition("=")
        if (
            separator == ""
            or not path_text
            or not is_absolute_evidence_path(path_text)
            or DIGEST_PATTERN.fullmatch(digest) is None
        ):
            raise ToolError(
                f"{field} entry must use absolute PATH=64HEX: {encoded!r}"
            )
        if path_text in seen:
            raise ToolError(f"{field} repeats dependency path {path_text!r}")
        seen.add(path_text)
        entries.append((path_text, digest.lower()))
    return entries


def validate_comparable_build_metadata(
    build: dict[str, object], side: str, *, require_release: bool = True
) -> list[str]:
    """Validate the build provenance required for timing evidence."""
    errors: list[str] = []
    required_build_keys = {
        "configuration",
        "cmake_generator",
        "cmake_generator_platform",
        "cmake_generator_toolset",
        "cxx_compiler",
        "cxx_compiler_id",
        "cxx_compiler_version",
        "cxx_compiler_launcher",
        "cxx_compiler_target",
        "cxx_flags",
        "executable_linker_flags",
        "interprocedural_optimization",
        "sysroot",
        "toolchain_file",
        "flint_version",
        "dependencies",
        "optional_dependencies",
        "options",
    }
    missing = sorted(required_build_keys - set(build))
    if missing:
        errors.append(
            f"{side} build metadata is missing: {', '.join(missing)}"
        )

    if require_release and build.get("configuration") != "Release":
        errors.append(
            f"{side} timing build is not Release: "
            f"{build.get('configuration')!r}"
        )
    elif not is_nonempty_string(build.get("configuration")):
        errors.append(f"{side} build metadata has no configuration")

    for field in (
        "cmake_generator",
        "cxx_compiler",
        "cxx_compiler_id",
        "cxx_compiler_version",
        "flint_version",
    ):
        if field in build and not is_nonempty_string(build[field]):
            errors.append(
                f"{side} build metadata has no nonempty {field}"
            )

    for field in (
        "cxx_flags",
        "executable_linker_flags",
        "interprocedural_optimization",
    ):
        values = build.get(field)
        if field not in build:
            continue
        if not isinstance(values, dict):
            errors.append(f"{side} build metadata {field} is not an object")
            continue
        if "release" not in values or not isinstance(values["release"], str):
            errors.append(
                f"{side} build metadata {field} has no Release setting"
            )
        elif field == "cxx_flags" and not values["release"].strip():
            errors.append(
                f"{side} build metadata has no nonempty Release C++ flags"
            )

    options = build.get("options")
    option_values: dict[str, bool | None] = {}
    if "options" in build and (
        not isinstance(options, dict) or not options
    ):
        errors.append(f"{side} build metadata has no populated options object")
    elif isinstance(options, dict):
        missing_options = sorted(
            (set(REQUIRED_BOOLEAN_OPTIONS) | {"unity_build_batch_size"})
            - set(options)
        )
        if missing_options:
            errors.append(
                f"{side} build metadata options are missing: "
                + ", ".join(missing_options)
            )
        for option, allow_unset in REQUIRED_BOOLEAN_OPTIONS.items():
            if option not in options:
                continue
            parsed, error = parse_cmake_boolean(
                options[option],
                field=f"{side} build option {option}",
                allow_unset=allow_unset,
            )
            option_values[option] = parsed
            if error is not None:
                errors.append(error)
        batch_size = options.get("unity_build_batch_size")
        if "unity_build_batch_size" in options and (
            not isinstance(batch_size, str)
            or (
                batch_size.strip()
                and (
                    re.fullmatch(r"[1-9][0-9]*", batch_size.strip()) is None
                )
            )
        ):
            errors.append(
                f"{side} build option unity_build_batch_size is invalid"
            )

    optional_dependencies = build.get("optional_dependencies")
    if "optional_dependencies" in build and not isinstance(
        optional_dependencies, dict
    ):
        errors.append(
            f"{side} build metadata optional_dependencies is not an object"
        )
    elif isinstance(optional_dependencies, dict):
        optional_names = set(optional_dependencies)
        if optional_names != KNOWN_OPTIONAL_DEPENDENCIES:
            missing_optional = sorted(
                KNOWN_OPTIONAL_DEPENDENCIES - optional_names
            )
            unknown_optional = sorted(
                optional_names - KNOWN_OPTIONAL_DEPENDENCIES
            )
            if missing_optional:
                errors.append(
                    f"{side} build metadata optional_dependencies are missing: "
                    + ", ".join(missing_optional)
                )
            if unknown_optional:
                errors.append(
                    f"{side} build metadata has unknown optional dependencies: "
                    + ", ".join(unknown_optional)
                )
        for dependency_name, dependency in optional_dependencies.items():
            if not isinstance(dependency, dict):
                errors.append(
                    f"{side} optional dependency {dependency_name} "
                    "metadata is not an object"
                )
                continue
            enabled, enabled_error = parse_cmake_boolean(
                dependency.get("enabled"),
                field=(
                    f"{side} optional dependency {dependency_name} enabled"
                ),
            )
            if enabled_error is not None:
                errors.append(enabled_error)
                continue
            option_enabled = option_values.get(dependency_name)
            if option_enabled is not None and enabled != option_enabled:
                errors.append(
                    f"{side} optional dependency {dependency_name} enabled "
                    "state disagrees with its build option"
                )
            if not enabled:
                digest_value = dependency.get("library_digests")
                if is_nonempty_string(digest_value):
                    try:
                        parse_library_digest_list(
                            digest_value,
                            field=(
                                f"{side} optional dependency {dependency_name} "
                                "library_digests"
                            ),
                        )
                    except ToolError as exc:
                        errors.append(str(exc))
                continue
            for field in ("include_dirs", "libraries", "library_digests"):
                if not is_nonempty_string(dependency.get(field)):
                    errors.append(
                        f"{side} enabled optional dependency "
                        f"{dependency_name} has no nonempty {field}"
                    )
            if is_nonempty_string(dependency.get("library_digests")):
                try:
                    parse_library_digest_list(
                        dependency["library_digests"],
                        field=(
                            f"{side} optional dependency {dependency_name} "
                            "library_digests"
                        ),
                    )
                except ToolError as exc:
                    errors.append(str(exc))

    dependencies = build.get("dependencies")
    if "dependencies" in build and not isinstance(dependencies, dict):
        errors.append(f"{side} build metadata dependencies is not an object")
        return errors
    if not isinstance(dependencies, dict):
        return errors

    required_dependency_fields = {
        "flint": (
            "version",
            "include_dirs",
            "library_dirs",
            "libraries",
            "library_digests",
        ),
        "google_benchmark": (
            "version",
            "libraries",
            "library_digests",
        ),
    }
    for dependency_name, fields in required_dependency_fields.items():
        dependency = dependencies.get(dependency_name)
        if not isinstance(dependency, dict):
            errors.append(
                f"{side} build metadata has no {dependency_name} "
                "dependency object"
            )
            continue
        for field in fields:
            if not is_nonempty_string(dependency.get(field)):
                errors.append(
                    f"{side} {dependency_name} dependency has no nonempty "
                    f"{field}"
                )
        digest_value = dependency.get("library_digests")
        if is_nonempty_string(digest_value):
            try:
                parse_library_digest_list(
                    digest_value,
                    field=(
                        f"{side} {dependency_name} dependency library_digests"
                    ),
                )
            except ToolError as exc:
                errors.append(str(exc))
    return errors


def manifest_dependency_digest_entries(
    manifest: dict[str, object],
) -> list[tuple[str, str, str]]:
    build = manifest.get("build")
    if not isinstance(build, dict):
        raise ToolError("benchmark manifest has no build object")
    errors = validate_comparable_build_metadata(
        build, "manifest", require_release=False
    )
    if errors:
        raise ToolError("; ".join(errors))

    result: list[tuple[str, str, str]] = []
    dependencies = build.get("dependencies")
    if not isinstance(dependencies, dict):
        raise ToolError("benchmark manifest dependencies are not an object")
    for dependency_name in ("flint", "google_benchmark"):
        dependency = dependencies.get(dependency_name)
        if not isinstance(dependency, dict):
            raise ToolError(
                f"benchmark manifest has no {dependency_name} dependency"
            )
        for path_text, digest in parse_library_digest_list(
            dependency.get("library_digests"),
            field=f"manifest {dependency_name} dependency library_digests",
        ):
            result.append((dependency_name, path_text, digest))

    optional_dependencies = build.get("optional_dependencies")
    if not isinstance(optional_dependencies, dict):
        raise ToolError(
            "benchmark manifest optional_dependencies are not an object"
        )
    for dependency_name in sorted(KNOWN_OPTIONAL_DEPENDENCIES):
        dependency = optional_dependencies.get(dependency_name)
        if not isinstance(dependency, dict):
            raise ToolError(
                f"benchmark manifest has no optional dependency {dependency_name}"
            )
        enabled, error = parse_cmake_boolean(
            dependency.get("enabled"),
            field=f"manifest optional dependency {dependency_name} enabled",
        )
        if error is not None:
            raise ToolError(error)
        if not enabled:
            continue
        for path_text, digest in parse_library_digest_list(
            dependency.get("library_digests"),
            field=(
                f"manifest optional dependency {dependency_name} "
                "library_digests"
            ),
        ):
            result.append((dependency_name, path_text, digest))

    seen: dict[str, str] = {}
    for _dependency_name, path_text, digest in result:
        previous = seen.get(path_text)
        if previous is not None and previous != digest:
            raise ToolError(
                f"manifest records conflicting digests for {path_text}"
            )
        seen[path_text] = digest
    return result


def verify_manifest_dependency_files(
    manifest: dict[str, object],
) -> dict[str, object]:
    verified: list[dict[str, object]] = []
    for dependency_name, path_text, expected_digest in (
        manifest_dependency_digest_entries(manifest)
    ):
        path = Path(path_text)
        try:
            snapshot = file_snapshot(path)
        except OSError as exc:
            raise ToolError(
                f"could not verify manifest dependency {path}: {exc}"
            ) from exc
        actual_digest = str(snapshot["sha256"])
        if actual_digest.lower() != expected_digest:
            raise ToolError(
                f"manifest dependency digest changed for {path}: "
                f"expected {expected_digest}, found {actual_digest}"
            )
        verified.append(
            {
                "dependency": dependency_name,
                **snapshot,
                "expected_sha256": expected_digest,
            }
        )
    return {
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": verified,
    }


def compare_run_metadata(
    baseline: dict[str, object],
    candidate: dict[str, object],
    *,
    current_helper_digest: str | None = None,
) -> list[str]:
    errors: list[str] = []
    baseline_helper = baseline.get("helper")
    candidate_helper = candidate.get("helper")
    baseline_helper_digest = (
        baseline_helper.get("sha256")
        if isinstance(baseline_helper, dict)
        else None
    )
    candidate_helper_digest = (
        candidate_helper.get("sha256")
        if isinstance(candidate_helper, dict)
        else None
    )
    if not is_nonempty_string(baseline_helper_digest) or not is_nonempty_string(
        candidate_helper_digest
    ):
        errors.append("run metadata has no helper digest")
    elif baseline_helper_digest != candidate_helper_digest:
        errors.append("baseline and candidate validation helpers differ")
    else:
        current_digest = current_helper_digest or file_sha256(
            Path(__file__).resolve()
        )
        if baseline_helper_digest != current_digest:
            errors.append(
                "run metadata validation helper does not match the current helper"
            )

    baseline_build = metadata_build_for_comparison(baseline)
    candidate_build = metadata_build_for_comparison(candidate)
    if not isinstance(baseline_build, dict) or not isinstance(
        candidate_build, dict
    ):
        errors.append("run metadata has no build object")
    else:
        for side, build in (
            ("baseline", baseline_build),
            ("candidate", candidate_build),
        ):
            errors.extend(validate_comparable_build_metadata(build, side))
        if baseline_build != candidate_build:
            errors.append("baseline and candidate build metadata differ")

    baseline_executable = baseline.get("executable")
    candidate_executable = candidate.get("executable")
    baseline_target = (
        baseline_executable.get("target")
        if isinstance(baseline_executable, dict)
        else None
    )
    candidate_target = (
        candidate_executable.get("target")
        if isinstance(candidate_executable, dict)
        else None
    )
    if not isinstance(baseline_target, str) or not baseline_target:
        errors.append("baseline metadata has no benchmark target")
    if not isinstance(candidate_target, str) or not candidate_target:
        errors.append("candidate metadata has no benchmark target")
    if baseline_target != candidate_target:
        errors.append(
            "baseline and candidate benchmark targets differ: "
            f"{baseline_target!r} != {candidate_target!r}"
        )

    if baseline.get("selection") != candidate.get("selection"):
        errors.append("baseline and candidate benchmark selections differ")

    environment_keys = (
        "platform",
        "machine",
        "cpu_model",
        "cpu_count",
        "cpu_affinity",
    )
    for key in environment_keys:
        baseline_environment = baseline.get("environment")
        candidate_environment = candidate.get("environment")
        if (
            not isinstance(baseline_environment, dict)
            or key not in baseline_environment
            or not isinstance(candidate_environment, dict)
            or key not in candidate_environment
        ):
            errors.append(f"run metadata is missing environment field {key}")
        elif baseline_environment[key] != candidate_environment[key]:
            errors.append(
                f"baseline and candidate environment field {key} differs"
            )
    return errors


def compare_raw_contexts(
    baseline: object, candidate: object
) -> list[str]:
    if not isinstance(baseline, dict) or not isinstance(candidate, dict):
        return ["raw benchmark payload is not an object"]
    baseline_context = baseline.get("context")
    candidate_context = candidate.get("context")
    if not isinstance(baseline_context, dict) or not isinstance(
        candidate_context, dict
    ):
        return ["raw benchmark context is not an object"]
    errors: list[str] = []
    comparable_keys = (
        "host_name",
        "num_cpus",
        "cpu_scaling_enabled",
        "aslr_enabled",
        "caches",
        "library_version",
        "library_build_type",
        "json_schema_version",
    )
    for key in comparable_keys:
        if key not in baseline_context or key not in candidate_context:
            errors.append(f"raw benchmark context is missing field {key}")
        elif baseline_context[key] != candidate_context[key]:
            errors.append(
                f"baseline and candidate raw context field {key} differs"
            )
    if baseline_context.get("library_build_type") != "release":
        errors.append("baseline Google Benchmark library is not a release build")
    if candidate_context.get("library_build_type") != "release":
        errors.append("candidate Google Benchmark library is not a release build")
    return errors


def validate_comparator_statistics(
    payload: object,
    expected_samples: dict[str, int],
    expected_units: dict[str, str],
    *,
    require_utest: bool,
) -> list[str]:
    if not isinstance(payload, list) or not payload:
        return ["Google Benchmark comparison report is empty or not an array"]
    errors: list[str] = []
    names: dict[str, int] = {}
    expected_names = set(expected_samples)
    if set(expected_units) != expected_names:
        errors.append("comparison unit inventory does not match selected rows")
    aggregate_suffixes = {
        "_mean": "mean",
        "_median": "median",
        "_stddev": "stddev",
        "_cv": "cv",
    }
    for index, value in enumerate(payload):
        if not isinstance(value, dict):
            errors.append(f"comparison row {index} is not an object")
            continue
        name = value.get("name")
        if not isinstance(name, str) or not name:
            errors.append(f"comparison row {index} has no nonempty name")
            continue
        names[name] = names.get(name, 0) + 1
        if names[name] > 1:
            errors.append(f"comparison report has duplicate row {name}")

        expected_aggregate: str | None = None
        associated_base: str | None = None
        if name == "OVERALL_GEOMEAN":
            expected_aggregate = "geomean"
        elif name not in expected_names:
            for base_name in expected_names:
                for suffix, aggregate_name in aggregate_suffixes.items():
                    if name == base_name + suffix:
                        expected_aggregate = aggregate_name
                        associated_base = base_name
                        break
                if expected_aggregate is not None:
                    break
            if expected_aggregate is None:
                errors.append(f"comparison report has unexpected row {name}")

        measurements = value.get("measurements")
        if not isinstance(measurements, list) or not measurements:
            errors.append(f"comparison row {name} has no measurements")
            continue
        if name in expected_names:
            expected_count = expected_samples[name]
            if len(measurements) != expected_count:
                errors.append(
                    f"comparison row {name} has {len(measurements)} "
                    f"measurements, expected {expected_count}"
                )
            if value.get("run_type") != "iteration":
                errors.append(
                    f"comparison row {name} is not an iteration row"
                )
            if value.get("aggregate_name") != "":
                errors.append(
                    f"comparison row {name} has invalid aggregate_name"
                )
        elif expected_aggregate is not None:
            if len(measurements) != 1:
                errors.append(
                    f"comparison aggregate row {name} must have one measurement"
                )
            if value.get("run_type") != "aggregate":
                errors.append(
                    f"comparison aggregate row {name} is not aggregate output"
                )
            if value.get("aggregate_name") != expected_aggregate:
                errors.append(
                    f"comparison aggregate row {name} has invalid "
                    "aggregate_name"
                )
        time_unit = value.get("time_unit")
        if name in expected_names:
            expected_unit = expected_units.get(name)
            if time_unit != expected_unit:
                errors.append(
                    f"comparison row {name} has time_unit {time_unit!r}, "
                    f"expected {expected_unit!r}"
                )
        elif associated_base is not None:
            expected_unit = expected_units.get(associated_base)
            if time_unit != expected_unit:
                errors.append(
                    f"comparison aggregate row {name} has time_unit "
                    f"{time_unit!r}, expected {expected_unit!r}"
                )
        elif expected_aggregate == "geomean" and time_unit not in TIME_UNITS:
            errors.append(
                f"comparison aggregate row {name} has invalid time_unit "
                f"{time_unit!r}"
            )
        for measurement in measurements:
            if not isinstance(measurement, dict):
                errors.append(f"comparison row {name} has a non-object measurement")
                continue
            for field in (
                "real_time",
                "cpu_time",
                "real_time_other",
                "cpu_time_other",
                "time",
                "cpu",
            ):
                if not is_number(measurement.get(field)):
                    errors.append(
                        f"comparison row {name} has invalid measurement {field}"
                    )
                elif field in {
                    "real_time",
                    "cpu_time",
                    "real_time_other",
                    "cpu_time_other",
                } and float(measurement[field]) < 0.0:
                    errors.append(
                        f"comparison row {name} has negative measurement {field}"
                    )
        if name in expected_names and require_utest:
            utest = value.get("utest")
            expected_count = expected_samples[name]
            if not isinstance(utest, dict):
                errors.append(f"comparison row {name} has no U-test result")
            else:
                if utest.get("have_optimal_repetitions") is not True:
                    errors.append(
                        f"comparison row {name} lacks optimal U-test repetitions"
                    )
                for field in ("cpu_pvalue", "time_pvalue"):
                    pvalue = utest.get(field)
                    if (
                        not is_number(pvalue)
                        or float(pvalue) < 0.0
                        or float(pvalue) > 1.0
                    ):
                        errors.append(
                            f"comparison row {name} has invalid U-test {field}"
                        )
                for field in (
                    "nr_of_repetitions",
                    "nr_of_repetitions_other",
                ):
                    if not values_equal(utest.get(field), expected_count):
                        errors.append(
                            f"comparison row {name} has invalid U-test {field}"
                        )
    missing = sorted(expected_names - set(names))
    if missing:
        errors.append(
            "comparison report is missing benchmark rows: " + ", ".join(missing)
        )
    return errors


def comparison_options(args: argparse.Namespace) -> dict[str, object]:
    return {
        "required_counters": {
            name: value for name, value in args.require_counter
        },
        "required_row_counters": [
            {"row": row, "counter": counter, "value": value}
            for row, counter, value in args.require_row_counter
        ],
        "compared_counters": list(args.compare_counter),
        "allow_fewer_repetitions": args.allow_fewer_repetitions,
        "allow_incomparable_metadata": args.allow_incomparable_metadata,
        "mann_whitney_enabled": not args.no_utest,
        "mann_whitney_alpha": None if args.no_utest else args.alpha,
        "comparator_timeout_seconds": args.timeout,
    }


def compare_benchmarks(args: argparse.Namespace) -> int:
    args.require_counter = normalized_counter_requirements(
        list(args.require_counter), mandatory=REQUIRED_STATUS_COUNTERS
    )
    helper_snapshot = file_snapshot(Path(__file__).resolve())
    baseline_path = args.baseline.expanduser().resolve()
    candidate_path = args.candidate.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if baseline_path == candidate_path:
        raise ToolError("baseline and candidate must be distinct raw JSON files")
    if output in (baseline_path, candidate_path):
        raise ToolError("comparison output may not overwrite an input run")
    google_output = companion_path(output, ".google.json")
    stdout_path = companion_path(output, ".stdout.txt")
    stderr_path = companion_path(output, ".stderr.txt")
    baseline, baseline_bytes, baseline_snapshot = read_json_snapshot(
        baseline_path
    )
    candidate, candidate_bytes, candidate_snapshot = read_json_snapshot(
        candidate_path
    )
    output_paths = [output, google_output, stdout_path, stderr_path]
    input_paths = [baseline_path, candidate_path, Path(__file__).resolve()]
    for raw_path in (baseline_path, candidate_path):
        input_paths.extend(
            [
                companion_path(raw_path, ".metadata.json"),
                companion_path(raw_path, ".stdout.txt"),
                companion_path(raw_path, ".stderr.txt"),
            ]
        )
    input_paths.extend(compare_tool_candidates(args.compare_tool))
    reject_output_input_overlap(output_paths, input_paths)

    requirements = list(args.require_counter)
    row_requirements = list(args.require_row_counter)
    baseline_validation = validate_payload(
        baseline,
        required_counters=requirements,
        required_row_counters=row_requirements,
    )
    candidate_validation = validate_payload(
        candidate,
        required_counters=requirements,
        required_row_counters=row_requirements,
    )
    (
        baseline_metadata,
        baseline_metadata_errors,
        baseline_metadata_snapshot,
    ) = load_run_metadata(
        baseline_path,
        str(baseline_snapshot["sha256"]),
    )
    (
        candidate_metadata,
        candidate_metadata_errors,
        candidate_metadata_snapshot,
    ) = load_run_metadata(
        candidate_path,
        str(candidate_snapshot["sha256"]),
    )
    metadata_errors = [
        f"baseline: {error}" for error in baseline_metadata_errors
    ] + [f"candidate: {error}" for error in candidate_metadata_errors]
    if baseline_metadata is not None and candidate_metadata is not None:
        metadata_errors.extend(
            compare_run_metadata(
                baseline_metadata,
                candidate_metadata,
                current_helper_digest=str(helper_snapshot["sha256"]),
            )
        )
    if baseline_validation["ok"] and candidate_validation["ok"]:
        metadata_errors.extend(compare_raw_contexts(baseline, candidate))
    semantic_failures: list[str] = []
    if baseline_validation["ok"] and candidate_validation["ok"]:
        semantic_failures = compare_semantics(
            baseline, candidate, list(args.compare_counter)
        )

    sample_counts = list(
        dict(baseline_validation.get("samples_per_benchmark", {})).values()
    ) + list(dict(candidate_validation.get("samples_per_benchmark", {})).values())
    fewer_repetition_override_used = bool(
        args.allow_fewer_repetitions
        and sample_counts
        and min(sample_counts) < 9
    )
    metadata_override_used = bool(
        args.allow_incomparable_metadata and metadata_errors
    )
    evidence_class = (
        "tooling_only"
        if fewer_repetition_override_used
        or metadata_override_used
        or args.no_utest
        else "routine"
    )
    if (
        not args.allow_fewer_repetitions
        and sample_counts
        and min(sample_counts) < 9
    ):
        semantic_failures.append(
            "timing comparison requires at least 9 raw repetitions per row; "
            "use --allow-fewer-repetitions only for tooling checks or exploration"
        )

    validation_ok = bool(
        baseline_validation["ok"]
        and candidate_validation["ok"]
        and not semantic_failures
        and (not metadata_errors or args.allow_incomparable_metadata)
    )
    if not validation_ok:
        prepare_outputs(output_paths, args.overwrite)
        for stale in (google_output, stdout_path, stderr_path):
            discard_stale_output(stale, overwrite=args.overwrite)
        report = {
            "schema_version": SCHEMA_VERSION,
            "kind": "silex-native-benchmark-comparison",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "helper": helper_snapshot,
            "helper_invocation": helper_invocation(),
            "baseline": comparison_artifact_ref(
                baseline_snapshot, baseline_metadata_snapshot
            ),
            "candidate": comparison_artifact_ref(
                candidate_snapshot, candidate_metadata_snapshot
            ),
            "comparison_options": comparison_options(args),
            "overrides_used": {
                "fewer_repetitions": fewer_repetition_override_used,
                "incomparable_metadata": metadata_override_used,
            },
            "evidence_class": evidence_class,
            "validation": {
                "ok": False,
                "baseline": baseline_validation,
                "candidate": candidate_validation,
                "semantic_errors": semantic_failures,
                "metadata": {
                    "ok": not metadata_errors,
                    "override": args.allow_incomparable_metadata,
                    "errors": metadata_errors,
                },
            },
            "verdict": "invalid",
        }
        write_json(output, report, overwrite=args.overwrite)
        print("native benchmark comparison rejected:", file=sys.stderr)
        for failure in list(baseline_validation["errors"]):
            print(f"- baseline: {failure}", file=sys.stderr)
        for failure in list(candidate_validation["errors"]):
            print(f"- candidate: {failure}", file=sys.stderr)
        for failure in semantic_failures:
            print(f"- {failure}", file=sys.stderr)
        for failure in metadata_errors:
            print(f"- metadata: {failure}", file=sys.stderr)
        return 1

    compare_tool = discover_compare_tool(args.compare_tool)
    reject_output_input_overlap(output_paths, [compare_tool])
    prepare_outputs(output_paths, args.overwrite)
    prepare_stream_outputs([stdout_path, stderr_path], args.overwrite)
    compare_tool_snapshot = file_snapshot(compare_tool)
    comparator_bytes: bytes | None = None
    comparator_payload: object | None = None
    comparator_read_error: str | None = None
    with tempfile.TemporaryDirectory(prefix="silex-native-compare-") as directory:
        comparison_directory = Path(directory)
        comparator_baseline = comparison_directory / "baseline.json"
        comparator_candidate = comparison_directory / "candidate.json"
        private_google_output = comparison_directory / "comparison.json"
        comparator_baseline.write_bytes(baseline_bytes)
        comparator_candidate.write_bytes(candidate_bytes)
        command = [
            sys.executable,
            str(compare_tool),
            "--no-color",
            "--dump_to_json",
            str(private_google_output),
        ]
        if args.no_utest:
            command.append("--no-utest")
        else:
            command.extend(["--alpha", str(args.alpha)])
        command.extend(
            [
                "benchmarks",
                str(comparator_baseline),
                str(comparator_candidate),
            ]
        )
        process_result = run_bounded_process(
            command,
            cwd=repo_root(),
            timeout=args.timeout,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            result_path=private_google_output,
            result_limit=MAX_JSON_BYTES,
        )
        if (
            not process_result["result_limit_exceeded"]
            and process_result["result_capture_error"] is None
        ):
            try:
                comparator_bytes, _comparator_after = read_bytes_snapshot(
                    private_google_output, limit=MAX_JSON_BYTES
                )
            except ToolError as exc:
                comparator_read_error = str(exc)
            else:
                try:
                    comparator_payload = json.loads(comparator_bytes)
                except (json.JSONDecodeError, UnicodeDecodeError) as exc:
                    comparator_read_error = f"invalid comparator JSON: {exc}"

    stdout = read_bounded_log_text(stdout_path)
    stderr = read_bounded_log_text(stderr_path)

    statistics: object = []
    comparator_error: str | None = None
    if process_result["timed_out"]:
        comparator_error = (
            "Google Benchmark compare.py timed out after "
            f"{args.timeout:g} seconds"
        )
    elif process_result["output_limit_exceeded"]:
        comparator_error = (
            "Google Benchmark compare.py exceeded the per-stream "
            "process log byte limit"
        )
    elif process_result["result_limit_exceeded"]:
        comparator_error = (
            "Google Benchmark compare.py result exceeded the JSON evidence "
            "byte limit; the private result was not published"
        )
    elif process_result["result_capture_error"] is not None:
        comparator_error = (
            "Google Benchmark compare.py result capture failed: "
            + str(process_result["result_capture_error"])
        )
    elif process_result["log_capture_error"] is not None:
        comparator_error = (
            "Google Benchmark compare.py log capture failed: "
            + str(process_result["log_capture_error"])
        )
    elif process_result["returncode"] != 0:
        comparator_error = (
            "Google Benchmark compare.py exited with "
            f"{process_result['returncode']}"
        )
    elif comparator_read_error is not None:
        comparator_error = comparator_read_error
    elif comparator_payload is None:
        comparator_error = "Google Benchmark compare.py produced no JSON report"
    else:
        try:
            statistics = comparator_payload
            baseline_units = {
                name: str(rows[0]["time_unit"])
                for name, rows in samples_by_name(baseline).items()
            }
            statistic_errors = validate_comparator_statistics(
                statistics,
                dict(baseline_validation["samples_per_benchmark"]),
                baseline_units,
                require_utest=(
                    not args.no_utest
                    and not fewer_repetition_override_used
                ),
            )
            if statistic_errors:
                comparator_error = "; ".join(statistic_errors)
        except ToolError as exc:
            comparator_error = str(exc)

    if comparator_bytes is not None and comparator_error is None:
        publish_bytes(
            google_output, comparator_bytes, overwrite=args.overwrite
        )
        comparator_sha256 = hashlib.sha256(comparator_bytes).hexdigest()
        comparator_published = True
    else:
        discard_stale_output(google_output, overwrite=args.overwrite)
        comparator_sha256 = None
        comparator_published = False

    report = {
        "schema_version": SCHEMA_VERSION,
        "kind": "silex-native-benchmark-comparison",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "helper": helper_snapshot,
        "helper_invocation": helper_invocation(),
        "baseline": comparison_artifact_ref(
            baseline_snapshot, baseline_metadata_snapshot
        ),
        "candidate": comparison_artifact_ref(
            candidate_snapshot, candidate_metadata_snapshot
        ),
        "comparison_options": comparison_options(args),
        "overrides_used": {
            "fewer_repetitions": fewer_repetition_override_used,
            "incomparable_metadata": metadata_override_used,
        },
        "evidence_class": evidence_class,
        "validation": {
            "ok": comparator_error is None,
            "baseline": baseline_validation,
            "candidate": candidate_validation,
            "semantic_errors": [],
            "metadata": {
                "ok": not metadata_errors,
                "override": args.allow_incomparable_metadata,
                "errors": metadata_errors,
            },
        },
        "comparator": {
            **compare_tool_snapshot,
            "command": command,
            "working_directory": str(repo_root()),
            **process_result,
            "timeout_seconds": args.timeout,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "raw_report": str(google_output),
            "raw_report_sha256": comparator_sha256,
            "raw_report_published": comparator_published,
            "error": comparator_error,
        },
        "statistics": statistics,
        # Statistical output is evidence, not a repository-wide numeric gate.
        "verdict": "measurement_only" if comparator_error is None else "invalid",
    }
    write_json(output, report, overwrite=args.overwrite)
    if stdout:
        sys.stdout.write(stdout)
    if metadata_override_used:
        print(
            "warning: incomparable metadata override used; report is tooling-only",
            file=sys.stderr,
        )
    if fewer_repetition_override_used:
        print(
            "warning: fewer-repetitions override used; report is tooling-only",
            file=sys.stderr,
        )
    if args.no_utest:
        print(
            "warning: Mann-Whitney testing disabled; report is tooling-only",
            file=sys.stderr,
        )
    if comparator_error is not None:
        print(comparator_error, file=sys.stderr)
        if stderr:
            sys.stderr.write(stderr)
        return 1
    print(f"comparison report: {output}")
    return 0


def add_execution_arguments(
    parser: argparse.ArgumentParser, *, gate: bool
) -> None:
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--target")
    parser.add_argument("--filter", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--label", default="semantic-gate" if gate else None,
                        required=not gate)
    parser.add_argument("--timeout", type=float, default=60.0 if gate else 600.0)
    parser.add_argument("--expect-row", action="append")
    parser.add_argument(
        "--require-counter",
        action="append",
        default=[],
        type=parse_counter_requirement,
        metavar="NAME=VALUE",
    )
    parser.add_argument(
        "--require-row-counter",
        action="append",
        default=[],
        type=parse_row_counter_requirement,
        metavar="ROW:NAME=VALUE",
    )
    parser.add_argument("--benchmark-arg", action="append", default=[])
    parser.add_argument("--note", action="append", default=[])
    parser.add_argument("--overwrite", action="store_true")
    if gate:
        parser.set_defaults(repetitions=1, min_time="1x")
    else:
        parser.add_argument("--repetitions", type=int, default=9)
        parser.add_argument("--min-time", default="0.1s")


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Silex-native Google Benchmark workflow helper"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser(
        "run", help="collect validated raw timing samples and metadata"
    )
    add_execution_arguments(run, gate=False)
    run.set_defaults(handler=run_benchmark)

    gate = subparsers.add_parser(
        "gate", help="run a one-iteration semantic benchmark gate"
    )
    add_execution_arguments(gate, gate=True)
    gate.set_defaults(handler=run_benchmark)

    compare = subparsers.add_parser(
        "compare", help="validate and statistically compare two raw runs"
    )
    compare.add_argument("--baseline", required=True, type=Path)
    compare.add_argument("--candidate", required=True, type=Path)
    compare.add_argument("--output", required=True, type=Path)
    compare.add_argument("--compare-tool", type=Path)
    compare.add_argument(
        "--require-counter",
        action="append",
        default=[],
        type=parse_counter_requirement,
        metavar="NAME=VALUE",
    )
    compare.add_argument(
        "--require-row-counter",
        action="append",
        default=[],
        type=parse_row_counter_requirement,
        metavar="ROW:NAME=VALUE",
    )
    compare.add_argument(
        "--compare-counter",
        action="append",
        default=[],
        type=parse_counter_name,
        metavar="NAME",
    )
    compare.add_argument("--allow-fewer-repetitions", action="store_true")
    compare.add_argument(
        "--allow-incomparable-metadata",
        action="store_true",
        help="allow missing or mismatched provenance only for tooling checks",
    )
    compare.add_argument("--no-utest", action="store_true")
    compare.add_argument("--alpha", type=float, default=0.05)
    compare.add_argument("--timeout", type=float, default=600.0)
    compare.add_argument("--overwrite", action="store_true")
    compare.set_defaults(handler=compare_benchmarks)
    return parser


def main() -> int:
    parser = create_parser()
    args = parser.parse_args()
    if hasattr(args, "repetitions") and args.repetitions <= 0:
        parser.error("--repetitions must be positive")
    if hasattr(args, "timeout") and (
        not math.isfinite(args.timeout) or args.timeout <= 0
    ):
        parser.error("--timeout must be finite and positive")
    if hasattr(args, "filter") and not args.filter:
        parser.error("--filter must be nonempty")
    if hasattr(args, "label") and not str(args.label).strip():
        parser.error("--label must be nonempty")
    if hasattr(args, "alpha") and not 0.0 < args.alpha < 1.0:
        parser.error("--alpha must be between 0 and 1")
    try:
        return int(args.handler(args))
    except (OSError, ToolError) as exc:
        print(f"native benchmark workflow failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
