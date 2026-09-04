Benchmarking
============

Silex uses Google Benchmark for native timing checkpoints.  Benchmarks answer
timing questions, profiling attributes cost, and tests plus debug checks
establish correctness.  Benchmark support is optional for ordinary builds and
required only when ``SILEX_BUILD_BENCHMARKS=ON``.

Build and verify
----------------

Configure the Release benchmark tree, build its registered targets, and run
the semantic suite with:

.. code-block:: console

   cmake --preset benchmarks
   cmake --build --preset benchmark-targets
   ctest --preset benchmarks

Configuration fails if Google Benchmark is unavailable.  The generated
manifest classifies every registered row as release-semantic or records an
explicit exclusion.  CTest audits that inventory and runs each release row
once; it does not apply wall-clock performance thresholds.

The manifest is written to
``build/benchmarks/bench/silex-native-benchmarks-Release.json``.  The native
helper validates subprocess results, row identity, repetition structure, and
Silex semantic counters.  A small gate can be run independently:

.. code-block:: console

   mkdir -p build/benchmark-results
   python3 tools/bench/native-benchmark.py gate \
     --manifest build/benchmarks/bench/silex-native-benchmarks-Release.json \
     --target b-silex-status \
     --filter '^(BM_status_ok|BM_version_string)$' \
     --expect-row BM_status_ok \
     --expect-row BM_version_string \
     --output build/benchmark-results/status-gate.json

Every release-semantic row must publish ``success=1`` and
``failure_reason=0`` after its postcondition.  Additional contracts may be
stated with ``--require-counter NAME=VALUE`` or
``--require-row-counter ROW:NAME=VALUE``.

Collect and compare timing evidence
-----------------------------------

Select only the affected rows and controls.  The ``run`` command defaults to
nine repetitions, preserves the raw iteration rows, and writes companion
metadata and process logs:

.. code-block:: console

   python3 tools/bench/native-benchmark.py run \
     --manifest build/benchmarks/bench/silex-native-benchmarks-Release.json \
     --target b-silex-nf_clgp \
     --filter '^BM_class_unit_0_1_0_real_quadratic_proven$' \
     --expect-row BM_class_unit_0_1_0_real_quadratic_proven \
     --label candidate-WORKTREE \
     --output build/benchmark-results/candidate.json

Collect baseline and candidate results in separate fresh Release build trees
with matching compiler, dependency, build, host, and affinity settings.  Then
compare them with:

.. code-block:: console

   python3 tools/bench/native-benchmark.py compare \
     --baseline build/benchmark-results/baseline.json \
     --candidate build/benchmark-results/candidate.json \
     --output build/benchmark-results/comparison.json

The comparison revalidates both inputs and their provenance before delegating
ratios and Mann--Whitney statistics to Google Benchmark's ``compare.py``.
NumPy and SciPy are required for that comparison support.  There is no
repository-wide percentage threshold: report a repeatable slowdown outside
control drift as a regression, and treat an inconclusive result as requiring a
controlled rerun.

Native workload surface
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Target
     - Workload
   * - ``b-silex-lat_hnf``, ``b-silex-lat_lll``,
       ``b-silex-lat_intersection``
     - Exact lattice kernels and reduction quality.
   * - ``b-silex-nf_ord_maximal_order``,
       ``b-silex-nf_ord_pmaximal_overorder``
     - Direct quadratic rows and native generic smoke rows.
   * - ``b-silex-nf_clgp``
     - Public class-group and paired class/unit rows.
   * - ``b-silex-nf_fac_elt``,
       ``b-silex-nf_fac_elt_compact_reconstruction``
     - Factored-element roots and compact reconstruction kernels.
   * - ``b-silex-nf_idl``, ``b-silex-nf_prime_idl``
     - Factor-over-base, local reduction, and valuation rows.
   * - ``b-silex-zeta_bf_linear_factor_count``,
       ``b-silex-zeta_bf_prime_scratch``
     - Zeta finite-field and per-prime scratch kernels.
   * - ``b-silex-fmpz_smat_*``, ``b-silex-status``
     - Native regression and smoke coverage.

Replay and profiling adapters
-----------------------------

Noninstalled source-tree adapters provide bounded reproduction and profiling,
not benchmark history.  For example:

.. code-block:: console

   python3 tools/bench/run-class-unit-instance.py \
     --field-id cubic_disc81_proven \
     --timeout 60 \
     --build-dir build/profile-release

Raw evidence belongs in ignored build trees or outside the repository.  See
:doc:`development/benchmarking` for contributor-facing impact classification,
baseline, evidence, and benchmark-authoring requirements.
