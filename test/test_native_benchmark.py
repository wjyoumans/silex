#!/usr/bin/env python3
"""Unit and generated-manifest checks for the native benchmark helper."""

from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import importlib.util
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.dont_write_bytecode = True


EXPECTED_TARGETS = {
    "b-silex-status",
    "b-silex-fmpz_smat_hnf_ctx",
    "b-silex-fmpz_smat_kernels",
    "b-silex-fmpz_smat_rank",
    "b-silex-lat_hnf",
    "b-silex-lat_intersection",
    "b-silex-lat_lll",
    "b-silex-nf_clgp",
    "b-silex-nf_clgp_factor_base_honesty",
    "b-silex-nf_fac_elt",
    "b-silex-nf_fac_elt_compact_reconstruction",
    "b-silex-nf_idl",
    "b-silex-nf_ord_maximal_order",
    "b-silex-nf_ord_pmaximal_overorder",
    "b-silex-nf_prime_idl",
    "b-silex-zeta_bf_linear_factor_count",
    "b-silex-zeta_bf_prime_scratch",
}


def load_tool(path: Path) -> object:
    spec = importlib.util.spec_from_file_location("silex_native_benchmark", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load Python tool: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def benchmark_payload(
    *,
    repetitions: int = 9,
    answer: int = 42,
    success: int = 1,
    names: tuple[str, ...] = ("BM_kernel", "BM_control"),
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for name_index, name in enumerate(names):
        for repetition in range(repetitions):
            rows.append(
                {
                    "name": name,
                    "run_name": name,
                    "run_type": "iteration",
                    "repetitions": repetitions,
                    "repetition_index": repetition,
                    "iterations": 1,
                    "real_time": float(100 + name_index * 10 + repetition),
                    "cpu_time": float(90 + name_index * 10 + repetition),
                    "time_unit": "ns",
                    "success": float(success),
                    "failure_reason": 0.0 if success else 7.0,
                    "answer": float(answer),
                }
            )
    return {
        "context": {
            "date": "2026-01-01T00:00:00+00:00",
            "host_name": "fixture",
            "json_schema_version": 1,
        },
        "benchmarks": rows,
    }


def run_metadata() -> dict[str, object]:
    return {
        "helper": {"sha256": "helper-digest"},
        "build": {
            "configuration": "Release",
            "cmake_generator": "Ninja",
            "cmake_generator_platform": "",
            "cmake_generator_toolset": "",
            "cxx_compiler": "/usr/bin/c++",
            "cxx_compiler_id": "GNU",
            "cxx_compiler_version": "16",
            "cxx_compiler_launcher": "",
            "cxx_compiler_target": "",
            "cxx_flags": {"all": "", "release": "-O3 -DNDEBUG"},
            "executable_linker_flags": {"all": "", "release": ""},
            "interprocedural_optimization": {"all": "", "release": ""},
            "sysroot": "",
            "toolchain_file": "",
            "flint_version": "3.6.0",
            "dependencies": {
                "flint": {
                    "version": "3.6.0",
                    "include_dirs": "/usr/include",
                    "library_dirs": "/usr/lib",
                    "libraries": "/usr/lib/libflint.so",
                    "library_digests": "/usr/lib/libflint.so=" + "a" * 64,
                },
                "google_benchmark": {
                    "version": "1.9.5",
                    "libraries": "/usr/lib/libbenchmark.a",
                    "library_digests": (
                        "/usr/lib/libbenchmark.a=" + "b" * 64
                    ),
                },
            },
            "optional_dependencies": {
                "fplll": {
                    "enabled": "OFF",
                    "include_dirs": "",
                    "libraries": "",
                    "library_digests": "",
                },
                "flatter": {
                    "enabled": "OFF",
                    "include_dirs": "",
                    "libraries": "",
                    "library_digests": "",
                },
            },
            "options": {
                "logging": "OFF",
                "debug_checks": "OFF",
                "profiling": "OFF",
                "sanitizers": "OFF",
                "frame_pointers": "OFF",
                "exceptions": "OFF",
                "rtti": "OFF",
                "shared_libraries": "",
                "position_independent_code": "",
                "unity_build": "",
                "unity_build_batch_size": "",
                "fplll": "OFF",
                "flatter": "OFF",
            },
        },
        "executable": {"target": "b-silex-fixture"},
        "selection": {"filter": "BM_kernel"},
        "environment": {
            "platform": "fixture-os",
            "machine": "fixture-machine",
            "cpu_model": "fixture-cpu",
            "cpu_count": None,
            "cpu_affinity": None,
        },
    }


def manifest_build(library: Path, *, configuration: str) -> dict[str, object]:
    build = copy.deepcopy(run_metadata()["build"])
    build["configuration"] = configuration
    resolved = library.resolve()
    digest = hashlib.sha256(resolved.read_bytes()).hexdigest()
    encoded_digest = f"{resolved}={digest}"
    for dependency_name in ("flint", "google_benchmark"):
        dependency = build["dependencies"][dependency_name]
        dependency["libraries"] = str(resolved)
        dependency["library_digests"] = encoded_digest
    build["dependencies"]["flint"]["include_dirs"] = str(resolved.parent)
    build["dependencies"]["flint"]["library_dirs"] = str(resolved.parent)
    return build


def comparison_measurement() -> dict[str, float]:
    return {
        "real_time": 100.0,
        "cpu_time": 90.0,
        "real_time_other": 101.0,
        "cpu_time_other": 91.0,
        "time": 0.01,
        "cpu": 0.02,
    }


def comparison_row(
    name: str,
    *,
    measurements: int = 9,
    aggregate_name: str | None = None,
) -> dict[str, object]:
    aggregate = aggregate_name is not None
    return {
        "name": name,
        "measurements": [
            comparison_measurement() for _ in range(measurements)
        ],
        "run_type": "aggregate" if aggregate else "iteration",
        "aggregate_name": aggregate_name or "",
        "time_unit": "s" if aggregate_name == "geomean" else "ns",
        "utest": (
            {}
            if aggregate
            else {
                "have_optimal_repetitions": True,
                "cpu_pvalue": 0.5,
                "time_pvalue": 0.5,
                "nr_of_repetitions": measurements,
                "nr_of_repetitions_other": measurements,
            }
        ),
    }


class NativeBenchmarkTests(unittest.TestCase):
    tool_path: Path
    tool: object
    cleaner_path: Path
    cleaner: object

    @classmethod
    def setUpClass(cls) -> None:
        cls.tool_path = Path(ARGS.tool).resolve()
        cls.tool = load_tool(cls.tool_path)
        cls.cleaner_path = cls.tool_path.with_name("clean-benchmark-results.py")
        cls.cleaner = load_tool(cls.cleaner_path)

    def compare_metadata(
        self, baseline: dict[str, object], candidate: dict[str, object]
    ) -> list[str]:
        return self.tool.compare_run_metadata(
            baseline,
            candidate,
            current_helper_digest="helper-digest",
        )

    def test_raw_payload_validation(self) -> None:
        result = self.tool.validate_payload(
            benchmark_payload(),
            expected_rows=["BM_kernel", "BM_control"],
            expected_repetitions=9,
            required_counters=[("answer", 42)],
            required_row_counters=[("BM_kernel", "answer", 42)],
        )
        self.assertTrue(result["ok"], result["errors"])
        self.assertEqual(result["sample_count"], 18)
        self.assertEqual(
            result["samples_per_benchmark"],
            {"BM_control": 9, "BM_kernel": 9},
        )

    def test_errors_and_failed_status_are_rejected(self) -> None:
        payload = benchmark_payload(repetitions=1, success=0)
        first = payload["benchmarks"][0]
        first["error_occurred"] = True
        first["error_message"] = "fixture skip"
        result = self.tool.validate_payload(payload, expected_repetitions=1)
        self.assertFalse(result["ok"])
        joined = "\n".join(result["errors"])
        self.assertIn("fixture skip", joined)
        self.assertIn("success=0.0", joined)
        self.assertIn("failure_reason=7.0", joined)

    def test_counter_requirements_reject_conflicting_duplicates(self) -> None:
        self.assertEqual(
            self.tool.normalized_counter_requirements(
                [("success", 1.0)],
                mandatory=self.tool.REQUIRED_STATUS_COUNTERS,
            ),
            [("success", 1), ("failure_reason", 0)],
        )
        with self.assertRaisesRegex(
            self.tool.ToolError, "success.*conflicting requirements"
        ):
            self.tool.normalized_counter_requirements(
                [("success", 0)],
                mandatory=self.tool.REQUIRED_STATUS_COUNTERS,
            )

    def test_git_metadata_without_git_is_unavailable(self) -> None:
        previous_path = os.environ.get("PATH")
        os.environ["PATH"] = ""
        try:
            with tempfile.TemporaryDirectory() as directory:
                self.assertEqual(
                    self.tool.git_metadata(Path(directory)),
                    {"available": False},
                )
        finally:
            if previous_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = previous_path

    def test_aggregate_only_output_is_rejected(self) -> None:
        payload = benchmark_payload(repetitions=1)
        row = payload["benchmarks"][0]
        row["name"] = "BM_kernel_median"
        row["run_type"] = "aggregate"
        row["aggregate_name"] = "median"
        payload["benchmarks"] = [row]
        result = self.tool.validate_payload(payload)
        self.assertFalse(result["ok"])
        self.assertTrue(
            any("aggregates-only" in error for error in result["errors"])
        )

    def test_repetition_indexes_must_match_raw_samples(self) -> None:
        payload = benchmark_payload(repetitions=2, names=("BM_kernel",))
        payload["benchmarks"][1]["repetition_index"] = 3
        result = self.tool.validate_payload(payload)
        self.assertFalse(result["ok"])
        self.assertTrue(
            any("repetition indexes" in error for error in result["errors"])
        )

    def test_repetition_indexes_are_required_for_repeated_samples(self) -> None:
        payload = benchmark_payload(repetitions=2, names=("BM_kernel",))
        for row in payload["benchmarks"]:
            del row["repetition_index"]
        result = self.tool.validate_payload(payload)
        self.assertFalse(result["ok"])
        self.assertTrue(
            any(
                "every repeated raw sample must have a repetition_index"
                in error
                for error in result["errors"]
            )
        )

    def test_unsupported_google_benchmark_schema_is_rejected(self) -> None:
        for value in (999, True, None):
            with self.subTest(value=value):
                payload = benchmark_payload(repetitions=1)
                if value is None:
                    del payload["context"]["json_schema_version"]
                else:
                    payload["context"]["json_schema_version"] = value
                result = self.tool.validate_payload(payload)
                self.assertFalse(result["ok"])
                self.assertTrue(
                    any(
                        "json_schema_version" in error
                        for error in result["errors"]
                    )
                )

    def test_oversized_numeric_value_is_rejected_without_overflow(self) -> None:
        payload = benchmark_payload(repetitions=1, names=("BM_kernel",))
        payload["benchmarks"][0]["real_time"] = 10**400
        result = self.tool.validate_payload(payload)
        self.assertFalse(result["ok"])
        self.assertTrue(
            any("real_time" in error for error in result["errors"])
        )

    def test_json_inputs_are_size_bounded(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "oversized.json"
            path.write_bytes(b" " * 33)
            original_limit = self.tool.MAX_JSON_BYTES
            self.tool.MAX_JSON_BYTES = 32
            try:
                with self.assertRaisesRegex(
                    self.tool.ToolError, "evidence limit"
                ):
                    self.tool.read_json(path)
            finally:
                self.tool.MAX_JSON_BYTES = original_limit

    def test_evidence_inputs_reject_special_files_without_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(
                self.tool.ToolError, "not a regular file"
            ):
                self.tool.file_snapshot(root)

            if not hasattr(os, "mkfifo"):
                return
            fifo = root / "baseline.fifo"
            os.mkfifo(fifo)
            candidate = root / "candidate.json"
            candidate.write_text(json.dumps(benchmark_payload()))
            completed = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "compare",
                    "--baseline",
                    str(fifo),
                    "--candidate",
                    str(candidate),
                    "--output",
                    str(root / "comparison.json"),
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=2.0,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("not a regular file", completed.stderr)

    def test_publication_is_exclusive_or_atomic_after_explicit_overwrite(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "evidence.json"
            destination.write_bytes(b"original")
            with self.assertRaisesRegex(
                self.tool.ToolError, "refusing to overwrite"
            ):
                self.tool.publish_bytes(
                    destination, b"candidate", overwrite=False
                )
            self.assertEqual(destination.read_bytes(), b"original")

            self.tool.publish_bytes(destination, b"candidate", overwrite=True)
            self.assertEqual(destination.read_bytes(), b"candidate")
            self.assertEqual(
                list(destination.parent.glob(f".{destination.name}.*.tmp")),
                [],
            )

    def test_bounded_process_logs_and_output_limit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stdout_path = root / "normal.stdout"
            stderr_path = root / "normal.stderr"
            result = self.tool.run_bounded_process(
                [
                    sys.executable,
                    "-c",
                    (
                        "import os; "
                        "os.write(1, b'o' * 2048); "
                        "os.write(2, b'e' * 2048)"
                    ),
                ],
                cwd=root,
                timeout=2.0,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                log_limit=4096,
            )
            self.assertEqual(result["returncode"], 0)
            self.assertFalse(result["output_limit_exceeded"])
            self.assertEqual(stdout_path.stat().st_size, 2048)
            self.assertEqual(stderr_path.stat().st_size, 2048)

            stdout_path = root / "limited.stdout"
            stderr_path = root / "limited.stderr"
            started = time.monotonic()
            result = self.tool.run_bounded_process(
                [
                    sys.executable,
                    "-c",
                    (
                        "import os,time; "
                        "os.write(1, b'x' * 1048576); time.sleep(5)"
                    ),
                ],
                cwd=root,
                timeout=2.0,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                log_limit=4096,
            )
            self.assertLess(time.monotonic() - started, 1.5)
            self.assertTrue(result["output_limit_exceeded"])
            self.assertEqual(result["termination_reason"], "output_limit")
            self.assertIsNone(result["returncode"])
            self.assertLessEqual(stdout_path.stat().st_size, 4096)

    def test_bounded_process_result_limit_terminates_and_records_state(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_path = root / "private-result.json"
            started = time.monotonic()
            result = self.tool.run_bounded_process(
                [
                    sys.executable,
                    "-c",
                    (
                        "from pathlib import Path; import sys,time; "
                        "Path(sys.argv[1]).write_bytes(b'x' * 1048576); "
                        "time.sleep(5)"
                    ),
                    str(result_path),
                ],
                cwd=root,
                timeout=2.0,
                stdout_path=root / "result-limit.stdout",
                stderr_path=root / "result-limit.stderr",
                result_path=result_path,
                result_limit=4096,
            )
            self.assertLess(time.monotonic() - started, 1.5)
            self.assertTrue(result["result_limit_exceeded"])
            self.assertEqual(result["termination_reason"], "result_limit")
            self.assertIsNone(result["returncode"])
            self.assertGreater(result["result_bytes_observed"], 4096)
            self.assertEqual(result["result_byte_limit"], 4096)

    def test_run_discards_private_result_that_exceeds_limit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "oversized-benchmark"
            executable.write_text(
                f"#!{sys.executable}\n"
                "import pathlib,sys,time\n"
                "argument = next(value for value in sys.argv "
                "if value.startswith('--benchmark_out='))\n"
                "pathlib.Path(argument.split('=', 1)[1]).write_bytes(b'x' * 1048576)\n"
                "time.sleep(5)\n"
            )
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            output = root / "published.json"
            args = argparse.Namespace(
                output=output,
                exe=executable,
                manifest=None,
                target=None,
                filter="BM_oversized",
                min_time="1x",
                repetitions=1,
                benchmark_arg=[],
                expect_row=None,
                require_counter=[],
                require_row_counter=[],
                command="gate",
                label="semantic-gate",
                note=[],
                overwrite=False,
                timeout=2.0,
            )
            original_limit = self.tool.MAX_JSON_BYTES
            self.tool.MAX_JSON_BYTES = 4096
            try:
                with contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(self.tool.run_benchmark(args), 1)
            finally:
                self.tool.MAX_JSON_BYTES = original_limit
            self.assertFalse(output.exists())
            metadata = json.loads(
                (root / "published.metadata.json").read_text()
            )
            self.assertTrue(metadata["process"]["result_limit_exceeded"])
            self.assertEqual(
                metadata["process"]["termination_reason"], "result_limit"
            )
            self.assertFalse(metadata["artifacts"]["raw_json_published"])
            self.assertIsNone(metadata["artifacts"]["raw_json_sha256"])

    def test_run_does_not_publish_valid_json_after_nonzero_exit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "failed-benchmark"
            payload = json.dumps(
                benchmark_payload(
                    repetitions=1,
                    names=("BM_failed",),
                )
            )
            executable.write_text(
                f"#!{sys.executable}\n"
                "import pathlib,sys\n"
                "argument = next(value for value in sys.argv "
                "if value.startswith('--benchmark_out='))\n"
                f"pathlib.Path(argument.split('=', 1)[1]).write_text({payload!r})\n"
                "raise SystemExit(7)\n"
            )
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            output = root / "published.json"
            output.write_bytes(b"stale result")
            args = argparse.Namespace(
                output=output,
                exe=executable,
                manifest=None,
                target=None,
                filter="BM_failed",
                min_time="1x",
                repetitions=1,
                benchmark_arg=[],
                expect_row=["BM_failed"],
                require_counter=[("success", 1), ("failure_reason", 0)],
                require_row_counter=[],
                command="gate",
                label="semantic-gate",
                note=[],
                overwrite=True,
                timeout=2.0,
            )
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.tool.run_benchmark(args), 1)

            self.assertFalse(output.exists())
            metadata = json.loads(
                (root / "published.metadata.json").read_text()
            )
            self.assertFalse(metadata["validation"]["ok"])
            self.assertEqual(metadata["process"]["returncode"], 7)
            self.assertFalse(metadata["artifacts"]["raw_json_published"])
            self.assertIsNone(metadata["artifacts"]["raw_json_sha256"])
            self.assertTrue((root / "published.stdout.txt").is_file())
            self.assertTrue((root / "published.stderr.txt").is_file())

    def test_run_requires_status_counters_without_explicit_flags(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "missing-status-benchmark"
            payload = benchmark_payload(
                repetitions=1,
                names=("BM_missing_status",),
            )
            del payload["benchmarks"][0]["success"]
            del payload["benchmarks"][0]["failure_reason"]
            encoded_payload = json.dumps(payload)
            executable.write_text(
                f"#!{sys.executable}\n"
                "import pathlib,sys\n"
                "argument = next(value for value in sys.argv "
                "if value.startswith('--benchmark_out='))\n"
                "pathlib.Path(argument.split('=', 1)[1]).write_text("
                f"{encoded_payload!r})\n"
            )
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            output = root / "published.json"
            args = argparse.Namespace(
                output=output,
                exe=executable,
                manifest=None,
                target=None,
                filter="BM_missing_status",
                min_time="1x",
                repetitions=1,
                benchmark_arg=[],
                expect_row=["BM_missing_status"],
                require_counter=[],
                require_row_counter=[],
                command="run",
                label="candidate",
                note=[],
                overwrite=False,
                timeout=2.0,
            )
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.tool.run_benchmark(args), 1)

            self.assertFalse(output.exists())
            metadata = json.loads(
                (root / "published.metadata.json").read_text()
            )
            self.assertEqual(
                metadata["selection"]["required_counters"],
                {"failure_reason": 0, "success": 1},
            )
            errors = "\n".join(metadata["validation"]["errors"])
            self.assertIn("missing counter success", errors)
            self.assertIn("missing counter failure_reason", errors)
            self.assertFalse(metadata["artifacts"]["raw_json_published"])

    @unittest.skipUnless(os.name == "posix", "POSIX process-group test")
    def test_bounded_process_timeout_stops_descendants(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            heartbeat = root / "heartbeat.txt"
            child_code = (
                "from pathlib import Path\n"
                "import time\n"
                f"path = Path({str(heartbeat)!r})\n"
                "while True:\n"
                "    path.write_text(str(time.monotonic()))\n"
                "    time.sleep(0.01)\n"
            )
            parent_code = (
                "import subprocess,sys,time\n"
                f"subprocess.Popen([sys.executable, '-c', {child_code!r}])\n"
                "time.sleep(5)\n"
            )
            result = self.tool.run_bounded_process(
                [sys.executable, "-c", parent_code],
                cwd=root,
                timeout=0.2,
                stdout_path=root / "timeout.stdout",
                stderr_path=root / "timeout.stderr",
                log_limit=4096,
            )
            self.assertTrue(result["timed_out"])
            self.assertTrue(heartbeat.is_file())
            stopped_value = heartbeat.read_text()
            time.sleep(0.1)
            self.assertEqual(heartbeat.read_text(), stopped_value)

    def test_semantic_counter_change_is_rejected(self) -> None:
        baseline = benchmark_payload()
        candidate = benchmark_payload(answer=43)
        failures = self.tool.compare_semantics(
            baseline, candidate, ["answer"]
        )
        self.assertEqual(len(failures), 2)
        self.assertTrue(all("counter answer differs" in row for row in failures))

    def test_metadata_accepts_equal_unsupported_environment_fields(self) -> None:
        baseline = run_metadata()
        candidate = run_metadata()
        self.assertEqual(
            self.compare_metadata(baseline, candidate),
            [],
        )
        candidate["build"]["cxx_flags"]["release"] = "-O2 -DNDEBUG"
        self.assertIn(
            "baseline and candidate build metadata differ",
            self.compare_metadata(baseline, candidate),
        )
        baseline = run_metadata()
        candidate = run_metadata()
        baseline["helper"]["sha256"] = ""
        self.assertIn(
            "run metadata has no helper digest",
            self.compare_metadata(baseline, candidate),
        )

    def test_manifest_requires_exact_kind_and_integer_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "manifest.json"
            payload = {
                "schema_version": 1,
                "kind": "silex-native-benchmark-executables",
                "targets": [],
            }
            manifest.write_text(json.dumps(payload))
            self.assertEqual(self.tool.load_manifest(manifest), payload)

            for field, value, message in (
                ("kind", "not-a-silex-manifest", "invalid kind"),
                ("schema_version", True, "unsupported schema_version"),
            ):
                with self.subTest(field=field, value=value):
                    invalid = dict(payload)
                    invalid[field] = value
                    manifest.write_text(json.dumps(invalid))
                    with self.assertRaisesRegex(self.tool.ToolError, message):
                        self.tool.load_manifest(manifest)

    def test_generated_manifest_check_is_fail_closed_under_optimization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "invalid-manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "schema_version": 999,
                        "kind": "wrong",
                        "targets": [],
                    }
                )
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    "-O",
                    str(Path(__file__).resolve()),
                    "--tool",
                    str(self.tool_path),
                    "--manifest",
                    str(manifest),
                    "--check-manifest-only",
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("unsupported schema_version", completed.stderr)

    def test_metadata_rejects_empty_critical_build_provenance(self) -> None:
        cases = (
            (("cmake_generator",), "cmake_generator"),
            (("cxx_compiler",), "cxx_compiler"),
            (("cxx_compiler_id",), "cxx_compiler_id"),
            (("cxx_compiler_version",), "cxx_compiler_version"),
            (("flint_version",), "flint_version"),
            (("cxx_flags", "release"), "Release C++ flags"),
            (("options",), "populated options object"),
            (("dependencies", "flint", "include_dirs"), "include_dirs"),
            (("dependencies", "flint", "library_dirs"), "library_dirs"),
            (("dependencies", "flint", "libraries"), "libraries"),
            (
                ("dependencies", "flint", "library_digests"),
                "library_digests",
            ),
            (
                ("dependencies", "google_benchmark", "libraries"),
                "libraries",
            ),
            (
                ("dependencies", "google_benchmark", "library_digests"),
                "library_digests",
            ),
        )
        for path, expected_error in cases:
            with self.subTest(path=path):
                baseline = run_metadata()
                candidate = run_metadata()
                value = baseline["build"]
                for key in path[:-1]:
                    value = value[key]
                value[path[-1]] = {} if path == ("options",) else ""
                errors = self.compare_metadata(baseline, candidate)
                self.assertTrue(
                    any(expected_error in error for error in errors), errors
                )

    def test_metadata_rejects_malformed_boolean_and_digest_provenance(self) -> None:
        baseline = run_metadata()
        candidate = run_metadata()
        baseline["build"]["optional_dependencies"]["fplll"][
            "enabled"
        ] = "MAYBE"
        errors = self.compare_metadata(baseline, candidate)
        self.assertTrue(
            any("invalid CMake boolean" in error for error in errors), errors
        )

        baseline = run_metadata()
        candidate = run_metadata()
        del baseline["build"]["optional_dependencies"]["flatter"]
        errors = self.compare_metadata(baseline, candidate)
        self.assertTrue(
            any(
                "optional_dependencies are missing: flatter" in error
                for error in errors
            ),
            errors,
        )

        baseline = run_metadata()
        candidate = run_metadata()
        baseline["build"]["dependencies"]["flint"][
            "library_digests"
        ] = "/usr/lib/libflint.so=x"
        errors = self.compare_metadata(baseline, candidate)
        self.assertTrue(
            any("absolute PATH=64HEX" in error for error in errors), errors
        )

        baseline = run_metadata()
        candidate = run_metadata()
        baseline["helper"]["sha256"] = "old-helper"
        candidate["helper"]["sha256"] = "old-helper"
        self.assertIn(
            "run metadata validation helper does not match the current helper",
            self.compare_metadata(baseline, candidate),
        )

    def test_manifest_dependencies_are_rehashed_for_nonrelease_gates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            dependency = Path(directory) / "libfixture.so"
            dependency.write_bytes(b"first dependency image")
            manifest = {
                "build": manifest_build(
                    dependency, configuration="RelWithDebInfo"
                )
            }
            verification = self.tool.verify_manifest_dependency_files(manifest)
            self.assertEqual(len(verification["files"]), 2)
            dependency.write_bytes(b"changed dependency image")
            with self.assertRaisesRegex(
                self.tool.ToolError, "dependency digest changed"
            ):
                self.tool.verify_manifest_dependency_files(manifest)

    def test_disabled_optional_dependency_paths_may_be_empty(self) -> None:
        baseline = run_metadata()
        candidate = run_metadata()
        optional_dependencies = {
            "fplll": {
                "enabled": "OFF",
                "include_dirs": "",
                "libraries": "",
                "library_digests": "",
            },
            "flatter": {
                "enabled": "OFF",
                "include_dirs": "",
                "libraries": "",
                "library_digests": "",
            },
        }
        baseline["build"]["optional_dependencies"] = optional_dependencies
        candidate["build"]["optional_dependencies"] = optional_dependencies
        self.assertEqual(
            self.compare_metadata(baseline, candidate), []
        )

        baseline["build"]["optional_dependencies"]["fplll"][
            "enabled"
        ] = "ON"
        errors = self.compare_metadata(baseline, candidate)
        self.assertTrue(
            any(
                "enabled optional dependency fplll has no nonempty "
                "include_dirs" in error
                for error in errors
            ),
            errors,
        )

    def test_comparator_statistics_require_exact_selected_rows(self) -> None:
        valid = [
            comparison_row("BM_kernel"),
            comparison_row(
                "BM_kernel_mean", measurements=1, aggregate_name="mean"
            ),
            comparison_row(
                "BM_kernel_median", measurements=1, aggregate_name="median"
            ),
            comparison_row(
                "BM_kernel_stddev", measurements=1, aggregate_name="stddev"
            ),
            comparison_row(
                "BM_kernel_cv", measurements=1, aggregate_name="cv"
            ),
            comparison_row(
                "OVERALL_GEOMEAN", measurements=1, aggregate_name="geomean"
            ),
        ]
        self.assertEqual(
            self.tool.validate_comparator_statistics(
                valid,
                {"BM_kernel": 9},
                {"BM_kernel": "ns"},
                require_utest=True,
            ),
            [],
        )

        missing_unit = comparison_row("BM_kernel")
        del missing_unit["time_unit"]
        negative_time = comparison_row("BM_kernel")
        negative_time["measurements"][0]["real_time"] = -1.0
        invalid_base_aggregate = comparison_row("BM_kernel")
        invalid_base_aggregate["aggregate_name"] = "mean"
        cases = (
            (
                [comparison_row("BM_kernel", measurements=8)],
                "8 measurements, expected 9",
            ),
            (
                [comparison_row("BM_kernel"), comparison_row("BM_kernel")],
                "duplicate row BM_kernel",
            ),
            (
                [comparison_row("BM_kernel"), comparison_row("BM_extra")],
                "unexpected row BM_extra",
            ),
            (
                [
                    comparison_row("BM_kernel"),
                    comparison_row(
                        "BM_kernel_mean",
                        measurements=2,
                        aggregate_name="mean",
                    ),
                ],
                "aggregate row BM_kernel_mean must have one measurement",
            ),
            ([missing_unit], "time_unit None"),
            ([negative_time], "negative measurement real_time"),
            ([invalid_base_aggregate], "invalid aggregate_name"),
        )
        for payload, expected_error in cases:
            with self.subTest(expected_error=expected_error):
                errors = self.tool.validate_comparator_statistics(
                    payload,
                    {"BM_kernel": 9},
                    {"BM_kernel": "ns"},
                    require_utest=True,
                )
                self.assertTrue(
                    any(expected_error in error for error in errors), errors
                )

    def test_gate_uses_manifest_and_writes_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "fake-benchmark"
            executable.write_text(
                """#!/usr/bin/env python3
import json
import sys
from pathlib import Path

def value(prefix):
    return next(arg.split('=', 1)[1] for arg in sys.argv if arg.startswith(prefix))

repetitions = int(value('--benchmark_repetitions='))
output = Path(value('--benchmark_out='))
rows = []
for name_index, name in enumerate(('BM_kernel', 'BM_control')):
    for repetition in range(repetitions):
        rows.append({
            'name': name,
            'run_name': name,
            'run_type': 'iteration',
            'repetitions': repetitions,
            'repetition_index': repetition,
            'iterations': 1,
            'real_time': 100.0 + name_index + repetition,
            'cpu_time': 90.0 + name_index + repetition,
            'time_unit': 'ns',
            'success': 1.0,
            'failure_reason': 0.0,
            'answer': 42.0,
        })
output.write_text(json.dumps({
    'context': {'json_schema_version': 1},
    'benchmarks': rows,
}))
print('fake benchmark complete')
"""
            )
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "kind": "silex-native-benchmark-executables",
                        "source_dir": str(root),
                        "build": manifest_build(
                            executable, configuration="RelWithDebInfo"
                        ),
                        "targets": [
                            {"target": "b-silex-fixture", "path": str(executable)}
                        ],
                    }
                )
            )
            output = root / "gate.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "gate",
                    "--manifest",
                    str(manifest),
                    "--target",
                    "b-silex-fixture",
                    "--filter=^(BM_kernel|BM_control)$",
                    "--expect-row",
                    "BM_kernel",
                    "--expect-row",
                    "BM_control",
                    "--require-counter",
                    "answer=42",
                    "--require-row-counter",
                    "BM_kernel:answer=42",
                    "--output",
                    str(output),
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            metadata_path = root / "gate.metadata.json"
            metadata = json.loads(metadata_path.read_text())
            self.assertEqual(metadata["mode"], "gate")
            self.assertEqual(metadata["executable"]["target"], "b-silex-fixture")
            self.assertEqual(metadata["working_directory"], str(root))
            self.assertTrue(metadata["validation"]["ok"])
            self.assertEqual(metadata["validation"]["sample_count"], 2)
            private_output_argument = next(
                argument
                for argument in metadata["command"]
                if argument.startswith("--benchmark_out=")
            )
            private_output = Path(private_output_argument.split("=", 1)[1])
            self.assertNotEqual(private_output, output)
            self.assertFalse(private_output.exists())
            self.assertTrue(metadata["artifacts"]["raw_json_published"])
            self.assertFalse(metadata["process"]["result_limit_exceeded"])
            self.assertTrue((root / "gate.stdout.txt").is_file())
            self.assertTrue((root / "gate.stderr.txt").is_file())

            original_output = output.read_text()
            invalid = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "gate",
                    "--manifest",
                    str(manifest),
                    "--target",
                    "b-silex-missing",
                    "--filter=BM_kernel",
                    "--output",
                    str(output),
                    "--overwrite",
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(invalid.returncode, 2)
            self.assertEqual(output.read_text(), original_output)

            invalid_timeout = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "gate",
                    "--exe",
                    str(executable),
                    "--filter=BM_kernel",
                    "--timeout",
                    "nan",
                    "--output",
                    str(root / "invalid-timeout.json"),
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(invalid_timeout.returncode, 2)
            self.assertIn("finite and positive", invalid_timeout.stderr)

    def test_compare_delegates_statistics_without_numeric_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline.json"
            candidate = root / "candidate.json"
            baseline.write_text(json.dumps(benchmark_payload()))
            candidate.write_text(json.dumps(benchmark_payload()))
            comparator = root / "compare.py"
            comparator.write_text(
                """import json
import sys
from pathlib import Path

output = Path(sys.argv[sys.argv.index('--dump_to_json') + 1])
measurement = {
    'real_time': 100.0,
    'cpu_time': 90.0,
    'real_time_other': 101.0,
    'cpu_time_other': 91.0,
    'time': 0.01,
    'cpu': 0.02,
}
utest = {
    'have_optimal_repetitions': True,
    'cpu_pvalue': 0.5,
    'time_pvalue': 0.5,
    'nr_of_repetitions': 9,
    'nr_of_repetitions_other': 9,
}
output.write_text(json.dumps([
    {
        'name': 'BM_kernel',
        'measurements': [measurement for _ in range(9)],
        'run_type': 'iteration',
        'aggregate_name': '',
        'time_unit': 'ns',
        'utest': utest,
    },
    {
        'name': 'BM_control',
        'measurements': [measurement for _ in range(9)],
        'run_type': 'iteration',
        'aggregate_name': '',
        'time_unit': 'ns',
        'utest': utest,
    },
]))
print('fixture comparison')
"""
            )
            report_path = root / "comparison.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "compare",
                    "--baseline",
                    str(baseline),
                    "--candidate",
                    str(candidate),
                    "--output",
                    str(report_path),
                    "--compare-tool",
                    str(comparator),
                    "--require-counter",
                    "success=1",
                    "--compare-counter",
                    "answer",
                    "--allow-incomparable-metadata",
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(report_path.read_text())
            self.assertEqual(report["verdict"], "measurement_only")
            self.assertEqual(report["evidence_class"], "tooling_only")
            self.assertTrue(report["validation"]["ok"])
            self.assertFalse(report["validation"]["metadata"]["ok"])
            self.assertTrue(
                report["comparison_options"]["allow_incomparable_metadata"]
            )
            self.assertEqual(
                report["comparison_options"]["required_counters"],
                {"failure_reason": 0, "success": 1},
            )
            self.assertEqual(report["statistics"][0]["name"], "BM_kernel")
            self.assertTrue((root / "comparison.google.json").is_file())
            self.assertTrue(report["comparator"]["raw_report_published"])
            dump_index = report["comparator"]["command"].index(
                "--dump_to_json"
            )
            private_report = Path(
                report["comparator"]["command"][dump_index + 1]
            )
            self.assertNotEqual(
                private_report, root / "comparison.google.json"
            )
            self.assertFalse(private_report.exists())

    def test_compare_requires_status_counters_without_explicit_flags(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = benchmark_payload(repetitions=1)
            for row in payload["benchmarks"]:
                del row["success"]
                del row["failure_reason"]
            baseline = root / "baseline.json"
            candidate = root / "candidate.json"
            baseline.write_text(json.dumps(payload))
            candidate.write_text(json.dumps(payload))
            output = root / "comparison.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(self.tool_path),
                    "compare",
                    "--baseline",
                    str(baseline),
                    "--candidate",
                    str(candidate),
                    "--output",
                    str(output),
                    "--allow-fewer-repetitions",
                    "--allow-incomparable-metadata",
                    "--no-utest",
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 1)
            report = json.loads(output.read_text())
            self.assertFalse(report["validation"]["ok"])
            errors = "\n".join(report["validation"]["baseline"]["errors"])
            self.assertIn("missing counter success", errors)
            self.assertIn("missing counter failure_reason", errors)
            self.assertEqual(
                report["comparison_options"]["required_counters"],
                {"failure_reason": 0, "success": 1},
            )

    def test_compare_timeout_invalidates_and_records_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline.json"
            candidate = root / "candidate.json"
            baseline.write_text(json.dumps(benchmark_payload()))
            candidate.write_text(json.dumps(benchmark_payload()))
            comparator = root / "compare.py"
            comparator.write_text(
                "import time\n"
                "time.sleep(2)\n"
            )
            report_path = root / "comparison.json"
            command = [
                sys.executable,
                str(self.tool_path),
                "compare",
                "--baseline",
                str(baseline),
                "--candidate",
                str(candidate),
                "--output",
                str(report_path),
                "--compare-tool",
                str(comparator),
                "--allow-incomparable-metadata",
                "--timeout",
                "0.05",
            ]
            completed = subprocess.run(
                command,
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            report = json.loads(report_path.read_text())
            self.assertEqual(report["verdict"], "invalid")
            self.assertFalse(report["validation"]["ok"])
            self.assertTrue(report["comparator"]["timed_out"])
            self.assertIsNone(report["comparator"]["returncode"])
            self.assertEqual(report["comparator"]["timeout_seconds"], 0.05)
            self.assertIn("timed out", report["comparator"]["error"])

            invalid_timeout = subprocess.run(
                [*command[:-1], "nan"],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(invalid_timeout.returncode, 2)
            self.assertIn("finite and positive", invalid_timeout.stderr)

    def test_compare_publishes_only_a_clean_validated_report(self) -> None:
        for case, returncode in (("invalid-statistics", 0), ("nonzero", 7)):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                baseline = root / "baseline.json"
                candidate = root / "candidate.json"
                baseline.write_text(json.dumps(benchmark_payload()))
                candidate.write_text(json.dumps(benchmark_payload()))
                comparator = root / "compare.py"
                comparator.write_text(
                    "import json,pathlib,sys\n"
                    "output = pathlib.Path("
                    "sys.argv[sys.argv.index('--dump_to_json') + 1])\n"
                    "output.write_text(json.dumps([]))\n"
                    f"raise SystemExit({returncode})\n"
                )
                report_path = root / "comparison.json"
                google_path = root / "comparison.google.json"
                google_path.write_bytes(b"stale report")
                completed = subprocess.run(
                    [
                        sys.executable,
                        str(self.tool_path),
                        "compare",
                        "--baseline",
                        str(baseline),
                        "--candidate",
                        str(candidate),
                        "--output",
                        str(report_path),
                        "--compare-tool",
                        str(comparator),
                        "--allow-incomparable-metadata",
                        "--overwrite",
                    ],
                    check=False,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                self.assertEqual(completed.returncode, 1, completed.stderr)
                self.assertFalse(google_path.exists())
                report = json.loads(report_path.read_text())
                self.assertEqual(report["verdict"], "invalid")
                self.assertFalse(report["validation"]["ok"])
                self.assertFalse(
                    report["comparator"]["raw_report_published"]
                )
                self.assertIsNone(report["comparator"]["raw_report_sha256"])
                self.assertTrue((root / "comparison.stdout.txt").is_file())
                self.assertTrue((root / "comparison.stderr.txt").is_file())

    def test_compare_discards_private_report_that_exceeds_limit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline.json"
            candidate = root / "candidate.json"
            baseline.write_text(json.dumps(benchmark_payload()))
            candidate.write_text(json.dumps(benchmark_payload()))
            comparator = root / "compare.py"
            comparator.write_text(
                "import pathlib,sys,time\n"
                "output = pathlib.Path(sys.argv[sys.argv.index('--dump_to_json') + 1])\n"
                "output.write_bytes(b'x' * 1048576)\n"
                "time.sleep(5)\n"
            )
            report_path = root / "comparison.json"
            args = argparse.Namespace(
                baseline=baseline,
                candidate=candidate,
                output=report_path,
                compare_tool=comparator,
                require_counter=[],
                require_row_counter=[],
                compare_counter=[],
                allow_fewer_repetitions=False,
                allow_incomparable_metadata=True,
                no_utest=False,
                alpha=0.05,
                timeout=2.0,
                overwrite=False,
                command="compare",
            )
            original_limit = self.tool.MAX_JSON_BYTES
            self.tool.MAX_JSON_BYTES = 8192
            try:
                with contextlib.redirect_stdout(io.StringIO()), \
                     contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(self.tool.compare_benchmarks(args), 1)
            finally:
                self.tool.MAX_JSON_BYTES = original_limit
            report = json.loads(report_path.read_text())
            self.assertTrue(report["comparator"]["result_limit_exceeded"])
            self.assertEqual(
                report["comparator"]["termination_reason"], "result_limit"
            )
            self.assertFalse(report["comparator"]["raw_report_published"])
            self.assertFalse((root / "comparison.google.json").exists())

    @unittest.skipIf(shutil.which("git") is None, "Git is not available")
    def test_cleaner_dry_run_apply_and_reserved_unicode_keep(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(
                ["git", "init", "--quiet", str(root)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            (root / ".gitignore").write_text("/build/benchmark-results/\n")
            subprocess.run(
                ["git", "-C", str(root), "add", ".gitignore"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            results = root / "build" / "benchmark-results"
            results.mkdir(parents=True)
            keep_raw = results / "résumé.metadata.json"
            drop_raw = results / "Δ.google.json"
            keep_family = self.cleaner.kept_result_family(keep_raw)
            drop_family = self.cleaner.kept_result_family(drop_raw)
            for path in keep_family | drop_family:
                path.write_text(f"fixture {path.name}\n")
            unrelated = results / "notes.bin"
            unrelated.write_bytes(b"preserve")

            command = [
                sys.executable,
                str(self.cleaner_path),
                "--root",
                str(root),
            ]
            dry_run = subprocess.run(
                command,
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(dry_run.returncode, 0, dry_run.stderr)
            self.assertIn("would remove 10", dry_run.stdout)
            self.assertTrue(all(path.is_file() for path in keep_family | drop_family))

            applied = subprocess.run(
                [
                    *command,
                    "--apply",
                    "--keep",
                    str(keep_raw.relative_to(root)),
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(applied.returncode, 0, applied.stderr)
            self.assertTrue(all(path.is_file() for path in keep_family))
            self.assertTrue(all(not path.exists() for path in drop_family))
            self.assertTrue(unrelated.is_file())


def check_generated_manifest(path: Path) -> None:
    payload = json.loads(path.read_text())
    def require(condition: bool, message: str) -> None:
        if not condition:
            raise RuntimeError(f"invalid generated benchmark manifest: {message}")

    require(isinstance(payload, dict), "top level is not an object")
    require(
        type(payload.get("schema_version")) is int
        and payload["schema_version"] == 1,
        "unsupported schema_version",
    )
    require(
        payload.get("kind") == "silex-native-benchmark-executables",
        "invalid kind",
    )
    require(isinstance(payload.get("source_dir"), str), "invalid source_dir")
    require(isinstance(payload.get("build_dir"), str), "invalid build_dir")
    build = payload.get("build")
    require(isinstance(build, dict), "build is not an object")
    if not isinstance(build, dict):
        return
    for field in (
        "options",
        "cxx_flags",
        "executable_linker_flags",
        "interprocedural_optimization",
        "dependencies",
        "optional_dependencies",
    ):
        require(isinstance(build.get(field), dict), f"invalid build.{field}")
    dependencies = build.get("dependencies")
    require(isinstance(dependencies, dict), "invalid dependencies")
    if not isinstance(dependencies, dict):
        return
    for dependency_name in ("flint", "google_benchmark"):
        dependency = dependencies.get(dependency_name)
        require(
            isinstance(dependency, dict),
            f"missing {dependency_name} dependency",
        )
        if isinstance(dependency, dict):
            require(
                isinstance(dependency.get("library_digests"), str)
                and bool(dependency["library_digests"]),
                f"empty {dependency_name} library_digests",
            )
    options = build.get("options")
    require(isinstance(options, dict), "invalid options")
    if isinstance(options, dict):
        require("shared_libraries" in options, "missing shared_libraries")
        require(
            "position_independent_code" in options,
            "missing position_independent_code",
        )
    rows = payload.get("targets")
    require(isinstance(rows, list), "targets is not an array")
    if not isinstance(rows, list):
        return
    require(all(isinstance(row, dict) for row in rows), "invalid target entry")
    target_names = {
        row.get("target") for row in rows if isinstance(row, dict)
    }
    require(target_names == EXPECTED_TARGETS, "target inventory differs")
    require(len(rows) == len(EXPECTED_TARGETS), "duplicate target entry")
    require(
        all(
            isinstance(row, dict)
            and isinstance(row.get("path"), str)
            and Path(row["path"]).is_file()
            for row in rows
        ),
        "target executable is missing",
    )
    for row in rows:
        require(isinstance(row, dict), "invalid target entry")
        if not isinstance(row, dict):
            continue
        release_rows = row.get("release_semantic_rows")
        excluded_rows = row.get("excluded_semantic_rows")
        require(
            isinstance(release_rows, list) and bool(release_rows),
            "target has no release semantic rows",
        )
        require(
            isinstance(release_rows, list)
            and all(isinstance(name, str) and name for name in release_rows),
            "invalid release semantic row",
        )
        if not isinstance(release_rows, list):
            continue
        require(
            len(release_rows) == len(set(release_rows)),
            "duplicate release semantic row",
        )
        require(isinstance(excluded_rows, list), "invalid excluded rows")
        if not isinstance(excluded_rows, list):
            continue
        require(
            all(
                isinstance(exclusion, dict)
                and isinstance(exclusion.get("name"), str)
                and bool(exclusion["name"])
                and isinstance(exclusion.get("reason"), str)
                and bool(exclusion["reason"])
                for exclusion in excluded_rows
            ),
            "invalid excluded semantic row",
        )
        excluded_names = [exclusion["name"] for exclusion in excluded_rows]
        require(
            len(excluded_names) == len(set(excluded_names)),
            "duplicate excluded semantic row",
        )
        require(
            set(release_rows).isdisjoint(excluded_names),
            "row is both released and excluded",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--check-manifest-only", action="store_true")
    return parser.parse_args()


ARGS = parse_args()


if __name__ == "__main__":
    if ARGS.check_manifest_only:
        if ARGS.manifest is None:
            raise SystemExit("--check-manifest-only requires --manifest")
        check_generated_manifest(ARGS.manifest)
    else:
        result = unittest.main(argv=[sys.argv[0]], exit=False)
        if not result.result.wasSuccessful():
            raise SystemExit(1)
        if ARGS.manifest is not None:
            check_generated_manifest(ARGS.manifest)
