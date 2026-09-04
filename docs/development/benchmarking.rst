Benchmarking and Performance Evidence
=====================================

Benchmarks measure performance, profiling attributes cost, and tests/debug
checks establish correctness.  Keep those roles separate.  The native targets
and helper commands are introduced in :doc:`../benchmarking`.

Impact classification
---------------------

Before editing runtime implementation, installed inline/template code, backend
routing, or performance-affecting flags, record the starting commit and dirty
state, affected symbols and callers, then classify the change:

``no runtime effect``
   Give a concrete rationale; do not time unrelated code.

``localized runtime path``
   Select one affected native microbenchmark and one unaffected control.

``algorithm, routing, or shared infrastructure``
   Add the smallest representative end-to-end row to the affected row and
   control.

``optimization or performance claim``
   Use the preceding coverage plus controlled repeated timing and profiling.

Build or preserve the baseline before implementation.  Baseline and candidate
must use separate fresh Release trees with the same compiler, dependencies,
options, and filters.  Use ``profile-release`` for attribution and
``benchmarks`` for timing.

Semantic release gate
---------------------

Before creating or moving a release tag:

.. code-block:: console

   cmake --preset benchmarks
   cmake --build --preset benchmarks
   ctest --preset benchmarks

Google Benchmark is optional for ordinary builds but mandatory for this gate.
The generated manifest classifies every row as release-semantic or explicitly
excluded.  Inventory checks reject missing, duplicate, vanished, or
unclassified rows, and semantic gates execute every release row once.  Every
release-semantic row must publish ``success=1`` and ``failure_reason=0`` after
its postcondition.

Timing comparisons
------------------

Use ``tools/bench/native-benchmark.py run`` for raw repetitions and
``compare`` for matching baseline/candidate evidence.  Routine comparisons use
at least nine repetitions.  Preserve:

* repository commit and dirty status;
* compiler, build type, flags, and dependency versions/digests;
* host, operating system, and affinity information;
* exact commands, filters, selected rows, and raw JSON;
* affected, end-to-end, and control results;
* medians, variability, ratios, and statistical output; and
* deterministic mathematical/proof counters that establish comparability.

The comparator must reject mismatched rows, repetitions, time units, host,
toolchain, options, dependencies, affinity, or helper metadata.  NumPy and
SciPy are required for Google Benchmark's comparison support.  A trusted
``compare.py`` path and its digest belong in the evidence.

There is no repository-wide percentage threshold.  Report one of ``no
regression detected``, ``regression``, or ``inconclusive--controlled rerun
required``.  A repeatable slowdown outside control drift must be fixed or
explicitly accepted; inconclusive is not a pass.  Do not collect final
CPU-sensitive evidence while other CPU-intensive work is active.

Benchmark authoring
-------------------

Use deterministic fixtures, keep setup and verification outside the timed
region, expose stable row names and semantic counters, and add bounded
representative sizes.  An optimization declares its target workload and
acceptance or stop criterion before implementation.  A component row explains
cost but does not replace end-to-end evidence.

Silex-specific replay adapters are noninstalled development tools.  Keep raw
results and companion metadata under ignored build trees or outside the
repository; only stable workloads, interfaces, and policy conclusions belong
in public documentation.
