Instrumentation
===============

Silex keeps user-facing progress, developer logging, correctness checks,
profiling, and reproducible benchmarking independent.  Enabling one channel
must not silently enable another, and instrumentation must not change
mathematical behavior.

Channels and gates
------------------

``DiagnosticsContext`` owns runtime channel configuration.  Library call sites
use the Silex-owned ``SILEX_VERBOSE``, ``SILEX_LOG``, ``SILEX_DEBUG_CHECK``,
``SILEX_PROFILE_SCOPE``, and ``SILEX_PROFILE_EVENT`` facades instead of direct
printing, third-party logging types, or ad hoc timers.

``SILEX_ENABLE_LOGGING``, ``SILEX_ENABLE_DEBUG_CHECKS``, and
``SILEX_ENABLE_PROFILING`` are independent compile-time gates.  Disabled
facades must not evaluate expensive arguments.  Verbose progress remains
runtime-controlled because it is the user-facing channel and is cheap when
disabled.  Google Benchmark is a separate build surface; see
:doc:`benchmarking`.

The 0.1.0 implementation has callback hooks rather than mandatory logging or
profiling dependencies.  spdlog and Tracy are not implemented backends, and
explicit requests for their unsupported CMake options fail configuration.
Optional-backend types must not enter installed headers.

Diagnostic example
------------------

Build the compiled example with all library-side diagnostic gates enabled:

.. code-block:: console

   cmake -S . -B build/diagnostics \
     -DSILEX_BUILD_EXAMPLES=ON \
     -DSILEX_ENABLE_LOGGING=ON \
     -DSILEX_ENABLE_DEBUG_CHECKS=ON \
     -DSILEX_ENABLE_PROFILING=ON
   cmake --build build/diagnostics --target example-diagnostics-hot-paths

Exercise the runtime channels separately:

.. code-block:: console

   SILEX_DIAGNOSTICS_VERBOSE=progress \
     build/diagnostics/examples/example-diagnostics-hot-paths
   SILEX_DIAGNOSTICS_LOG=detail \
     build/diagnostics/examples/example-diagnostics-hot-paths
   SILEX_DIAGNOSTICS_DEBUG=cheap \
     build/diagnostics/examples/example-diagnostics-hot-paths
   SILEX_DIAGNOSTICS_PROFILE=1 \
     build/diagnostics/examples/example-diagnostics-hot-paths

Profiling callback durations are inclusive and may double-count nested scopes.
Use them for attribution and quick inspection, not as performance evidence;
use Google Benchmark or an external profiler for controlled measurement.

Review requirements
-------------------

New instrumentation must use the central facade, preserve failure behavior,
remain cheap when disabled, and carry a concrete workload need.  Tests must
show that channels can be enabled independently, logging does not trigger
debug checks, and expensive checks run only when explicitly requested.  Add
finer scopes or counters only after a benchmark or investigation identifies a
specific attribution gap.
