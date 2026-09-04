Getting Started
===============

Silex is a source-built restricted native C++ library.  User code should
prefer domain objects such as ``NumberField``, ``Order``, ``Element``,
``Ideal``, ``ClassGroupContext``, and ``OrderUnitGroup`` plus Silex-owned
FLINT RAII support types.

Prerequisites
-------------

The current build expects:

* a C++20 compiler;
* CMake 3.20 or newer;
* ``pkg-config``;
* FLINT at least the version requested by ``SILEX_REQUIRED_FLINT_VERSION``;
* Python 3 when tests are enabled;
* Google Benchmark when benchmark targets are enabled;
* Python 3.12 or newer and the pinned Sphinx version in
  ``docs/requirements.txt`` when documentation targets are enabled.

The repository provides supported source-tree, ``add_subdirectory``, and
default installed-package consumption for the FLINT-backed configuration.
See :doc:`support_matrix` for the exact qualified stack and limitations.

Install And Packaging Status
----------------------------

Silex currently supports three CMake consumption modes:

* for local development, configure and link against the build tree;
* for a parent CMake project, use ``add_subdirectory(path/to/silex)`` and
  link ``Silex::silex``;
* for installed consumption, install to a prefix and use
  ``find_package(Silex CONFIG REQUIRED)`` before linking ``Silex::silex``.

The downstream CMake project under ``test/cmake/install-consumer`` configures,
builds, and runs against an installed prefix.  Automated release tests cover
the default layout, nondefault ``GNUInstallDirs`` include/library/documentation
directories, a clean ``add_subdirectory`` consumer, and a clean tracked source
archive.  The exported package depends on ``pkg-config`` FLINT discovery, as
the source build does.  ``LICENSE``, ``THIRD_PARTY_NOTICES.md``, ``NOTICE.md``,
and ``CITATION.cff`` install under ``CMAKE_INSTALL_DOCDIR``.

Configure And Build
-------------------

A normal development build enables tests and examples:

.. code-block:: console

   cmake -S . -B build \
     -DSILEX_BUILD_TESTS=ON \
     -DSILEX_BUILD_EXAMPLES=ON
   cmake --build build -j
   ctest --test-dir build --output-on-failure

The main library target is ``Silex::silex``.  In a parent CMake project that
vendors this repository, use:

.. code-block:: cmake

   add_subdirectory(path/to/silex)
   target_link_libraries(my_target PRIVATE Silex::silex)

The library exports public include paths for ``include/`` and the generated
build configuration header in the build tree.

To smoke-test installed consumption manually:

.. code-block:: console

   cmake -S . -B build/install-smoke -DSILEX_BUILD_TESTS=OFF -DSILEX_BUILD_EXAMPLES=OFF
   cmake --build build/install-smoke --target silex
   cmake --install build/install-smoke --prefix /tmp/silex-install
   cmake -S test/cmake/install-consumer -B /tmp/silex-consumer \
     -DCMAKE_PREFIX_PATH=/tmp/silex-install
   cmake --build /tmp/silex-consumer

Build Options
-------------

Common build options:

``SILEX_BUILD_TESTS``
   Build CTest tests.  Enabled by default.  The active source-name check also
   requires a Git index at the source root and is not registered in a
   generated source archive; CMake reports that omission.

``SILEX_BUILD_RELEASE_TESTS``
   Add the installed-package, nondefault-layout, ``add_subdirectory``, and
   source-archive release checks when tests are enabled.  Disabled by default;
   enable it for release-candidate validation.
   The source-archive check requires a Git index at the source root and is not
   registered when configuring a generated source archive; CMake reports that
   omission and the other release checks still run.  Disable the option for a
   nested or vendored build that should run only the native suite.

``SILEX_BUILD_EXAMPLES``
   Build compiled example programs under ``examples/``.  Enabled by default.

``SILEX_BUILD_BENCHMARKS``
   Build Google Benchmark targets under ``bench/``.  Disabled by default.

``SILEX_BUILD_DOCS``
   Add the strict ``docs-html`` and ``docs-linkcheck`` targets.  Disabled by
   default; configuration fails clearly if ``sphinx-build`` is unavailable
   when this option is enabled.

``SILEX_ENABLE_LOGGING``
   Compile logging facade support.  Enabled by default.

``SILEX_ENABLE_DEBUG_CHECKS``
   Compile debug-check support.  Disabled by default.

``SILEX_ENABLE_PROFILING``
   Compile profiling facade support.  Disabled by default.

``SILEX_ENABLE_SANITIZERS``
   Enable AddressSanitizer and UndefinedBehaviorSanitizer for supported
   compilers.  Disabled by default.

``SILEX_ENABLE_EXCEPTIONS`` and ``SILEX_ENABLE_RTTI``
   Enable exceptions or RTTI in the core target.  Both are disabled by default
   while the restricted C++ profile remains in force.

Examples
--------

Examples are built with ``SILEX_BUILD_EXAMPLES=ON``:

.. code-block:: console

   cmake --build build --target example-field-order-ideal-basics
   ./build/examples/example-field-order-ideal-basics

Examples are also registered as CTest tests when tests are enabled:

.. code-block:: console

   ctest --test-dir build -R silex-example --output-on-failure

See :doc:`examples` for the current example list.

Diagnostics Builds
------------------

Logging, debug checks, profiling, and verbose output are separate channels.
Compile the channels you need, then enable the example callbacks at runtime:

.. code-block:: console

   cmake -S . -B build/diagnostics \
     -DSILEX_BUILD_EXAMPLES=ON \
     -DSILEX_ENABLE_LOGGING=ON \
     -DSILEX_ENABLE_DEBUG_CHECKS=ON \
     -DSILEX_ENABLE_PROFILING=ON
   cmake --build build/diagnostics --target example-diagnostics-hot-paths
   SILEX_DIAGNOSTICS_PROFILE=1 \
     ./build/diagnostics/examples/example-diagnostics-hot-paths

Profiling callbacks add overhead.  Use them for attribution and quick
inspection, not as replacement for benchmark timing.

Benchmarks
----------

Build benchmarks explicitly:

.. code-block:: console

   cmake --preset benchmarks
   cmake --build --preset benchmark-targets
   mkdir -p build/benchmark-results
   python3 tools/bench/native-benchmark.py gate \
     --manifest build/benchmarks/bench/silex-native-benchmarks-Release.json \
     --target b-silex-status \
     --filter '^(BM_status_ok|BM_version_string)$' \
     --expect-row BM_status_ok \
     --expect-row BM_version_string \
     --output build/benchmark-results/status-gate.json

Before creating or moving a tag, build the full Release suite and run its
ordinary tests, benchmark inventory audit, and one-iteration semantic gates:

.. code-block:: console

   cmake --build --preset benchmarks
   ctest --preset benchmarks

Use the native helper's ``run`` and ``compare`` commands for performance
comparisons.  They preserve raw repetitions and record build, machine, command,
and commit metadata.  See :doc:`benchmarking` for the complete workflow and
controlled-rerun policy.

Documentation
-------------

Install the shared documentation dependency and build the development channel
with:

.. code-block:: console

   python3 -m pip install -r docs/requirements.txt
   cmake --preset docs
   cmake --build --preset docs
   cmake --build --preset docs-linkcheck

The presets use ``sphinx-build -n -W --keep-going`` so unresolved references
and warnings fail the HTML or link-check build.  Direct builds may set the
published version and channel explicitly:

.. code-block:: console

   SILEX_DOCS_RELEASE=0.1.1 SILEX_DOCS_CHANNEL=stable \
     sphinx-build -n -W --keep-going -b html docs build/docs/html

Use :doc:`reference/module_map` to map installed headers to the conceptual API
guide, :doc:`legal_and_provenance` for source anchors and license pointers, and
:doc:`releases/index` for channel semantics and versioned release notes.
