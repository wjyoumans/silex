Silex Documentation
===================

Silex is a native C++20 library for computational algebraic number theory over
FLINT-backed exact arithmetic.  This is the |docs-channel| documentation for
release |silex-release|.

.. only:: dev

   These pages describe the development branch.  Interfaces may differ from a
   tagged release; use the stable channel when reproducing a released result.

.. only:: stable

   These pages describe a tagged release.  The :doc:`support_matrix` defines
   its qualified platform and mathematical boundary.

The public declarations live under ``include/silex/``.  Start with
:doc:`getting_started`, then use the tutorials and the
:doc:`reference/index` to find the relevant domain types.  The reference is a
conceptual API guide, not generated symbol documentation; the headers shipped
with the selected release remain authoritative for signatures.

Build Quick Reference
---------------------

Configure, build, and test the default source tree with:

.. code-block:: console

   cmake -S . -B build -DSILEX_BUILD_TESTS=ON -DSILEX_BUILD_EXAMPLES=ON
   cmake --build build -j
   ctest --test-dir build --output-on-failure

Build public documentation and benchmark targets explicitly:

.. code-block:: console

   python3 -m pip install -r docs/requirements.txt
   cmake --preset docs
   cmake --build --preset docs
   cmake --build --preset docs-linkcheck

   cmake --preset benchmarks
   cmake --build --preset benchmark-targets

See :doc:`reference/module_map` for a header-to-module index,
:doc:`legal_and_provenance` for licenses and source anchors, and
:doc:`releases/index` for release-channel semantics and versioned notes.

.. toctree::
   :maxdepth: 2

   getting_started
   support_matrix
   tutorial/index
   reference/index
   examples
   benchmarking
   algorithms_and_sources
   legal_and_provenance
   releases/index
   development/index
