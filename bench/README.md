# Benchmarks

Silex uses Google Benchmark for native microbenchmarks and performance
regression tracking. Keep workloads explicit, reproducible, and separate from
correctness tests and profiling diagnostics.

Every runtime change requires performance-impact triage before editing. A
localized path needs a directly affected row and an unaffected control; an
algorithm, routing, or shared-infrastructure change also needs the smallest
representative end-to-end row. Record a concrete rationale when a change has
no runtime effect. See `docs/development/benchmarking.rst` for the full policy.

Configure and build with:

```sh
cmake --preset benchmarks
cmake --build --preset benchmark-targets
```

Google Benchmark is optional for ordinary builds. When
`SILEX_BUILD_BENCHMARKS=ON`—including the `benchmarks` and `profile-release`
presets—configuration fails if Google Benchmark is unavailable, so benchmark
semantic gates cannot disappear silently.

Before creating or moving a tag, build everything needed by the Release test
suite and run it:

```sh
cmake --build --preset benchmarks
ctest --preset benchmarks
```

The generated manifest classifies every registered row as release-semantic or
as an explicit exclusion with a reason. CTest audits that inventory and runs a
one-iteration semantic gate for every release row. A missing, duplicate,
vanished, or unclassified row fails the suite; timing thresholds are never
applied in CI.

CMake writes the executable and build-metadata manifest to
`build/benchmarks/bench/silex-native-benchmarks-Release.json`. Use the native
helper for semantic gates and timing collection rather than hand-assembling
Google Benchmark output flags:

```sh
mkdir -p build/benchmark-results
python3 tools/bench/native-benchmark.py gate \
  --manifest build/benchmarks/bench/silex-native-benchmarks-Release.json \
  --target b-silex-status \
  --filter '^(BM_status_ok|BM_version_string)$' \
  --expect-row BM_status_ok \
  --expect-row BM_version_string \
  --output build/benchmark-results/status-gate.json
```

Both `run` and `gate` reject invalid or aggregate-only JSON, missing and
unexpected rows, Google Benchmark skip/error rows, incomplete or duplicate
repetition indices, and missing or failed Silex semantic status counters.
Every release-semantic row must publish `success=1` and `failure_reason=0`
after its postcondition.
Those two counters are mandatory automatically, and `compare` revalidates
them in existing evidence. Use
`--require-counter NAME=VALUE` or
`--require-row-counter ROW:NAME=VALUE` to state additional contracts.

Current targets:

- `b-silex-status`: lightweight status/version scaffold.
- `b-silex-fmpz_smat_kernels`: sparse-times-dense and transpose kernels.
- `b-silex-fmpz_smat_rank`: one-shot and persistent modular-rank kernels.
- `b-silex-fmpz_smat_hnf_ctx`: incremental HNF-context kernels.
- `b-silex-lat_hnf`: exact HNF kernels.
- `b-silex-lat_lll`: native LLL plus optional backend comparison rows.
- `b-silex-lat_intersection`: exact lattice intersection.
- `b-silex-nf_ord_maximal_order`: maximal-order construction.
- `b-silex-nf_ord_pmaximal_overorder`: p-maximal overorder construction.
- `b-silex-nf_clgp`: class-group and paired class/unit rows.
- `b-silex-nf_clgp_factor_base_honesty`: focused internal factor-base rows.
- `b-silex-nf_fac_elt`: factored and compact element operations.
- `b-silex-nf_fac_elt_compact_reconstruction`: coordinate-bound, direct, and
  bounded-CRT compact reconstruction.
- `b-silex-nf_idl`: ideal factor-over-base operations.
- `b-silex-nf_prime_idl`: prime-ideal reduction and valuation operations.
- `b-silex-zeta_bf_linear_factor_count`: finite-field linear-factor kernels.
- `b-silex-zeta_bf_prime_scratch`: per-prime zeta scratch kernels.

For timing, build or preserve the baseline before implementation and use
separate fresh Release trees with identical compiler, dependency, and option
settings. `run` preserves raw samples and writes companion metadata, stdout,
and stderr files. It defaults to the nine repetitions needed for a routine
statistical comparison:

```sh
python3 tools/bench/native-benchmark.py run \
  --manifest build/benchmarks/bench/silex-native-benchmarks-Release.json \
  --target b-silex-nf_clgp \
  --filter '^BM_class_unit_0_1_0_real_quadratic_proven$' \
  --expect-row BM_class_unit_0_1_0_real_quadratic_proven \
  --label candidate-WORKTREE \
  --output build/benchmark-results/candidate.json

python3 tools/bench/native-benchmark.py compare \
  --baseline build/benchmark-results/baseline.json \
  --candidate build/benchmark-results/candidate.json \
  --output build/benchmark-results/comparison.json
```

The comparison revalidates both inputs and companion provenance, requires
resolved FLINT and Google Benchmark library digests, rejects mismatched build
or host metadata, and delegates ratios and Mann--Whitney statistics to Google
Benchmark's `compare.py`. The comparator is executable Python: use only the
copy discovered in a known Google Benchmark installation or pass a trusted
local file with `--compare-tool`. Its execution is bounded by `--timeout`, and
its path and digest are recorded. A timeout, nonzero exit, malformed report,
or row/sample mismatch invalidates the comparison.

The helper accepts only regular evidence files, caps each process log at 16
MiB and each JSON result at 64 MiB, and writes child-produced JSON to a private
path first. Validated JSON is then published without overwriting an existing
artifact unless `--overwrite` was explicitly supplied.

The workflow never applies a universal percentage threshold. A repeatable
slowdown outside control drift must be fixed or explicitly accepted; an
inconclusive result requires a controlled rerun and is not a pass. Keep raw
output in ignored build directories or outside this repository. The comparison
tool requires Python, NumPy, and SciPy. Preview or remove complete ignored
result families, including metadata and process logs, with:

```sh
python3 tools/bench/clean-benchmark-results.py
python3 tools/bench/clean-benchmark-results.py --apply
```
