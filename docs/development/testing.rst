Testing and Review
==================

Tests establish mathematical parity, ownership and lifetime safety, failure
preservation, documentation accuracy, and package behavior.  Add the smallest
focused test first, then run broader gates proportional to the change.

Development loop
----------------

The source-archive release test intentionally requires tracked working-tree
bytes to match the Git index.  During ordinary unstaged development, configure
a separate tree with release-only package/archive checks disabled:

.. code-block:: console

   cmake -S . -B build/dev \
     -DSILEX_BUILD_TESTS=ON \
     -DSILEX_BUILD_RELEASE_TESTS=OFF \
     -DSILEX_BUILD_EXAMPLES=ON
   cmake --build build/dev
   ctest --test-dir build/dev --output-on-failure

Before review, stage the exact intended paths and run the full applicable
preset so the source-archive and installed-consumer tests inspect the same
bytes that will be committed.  Release-only gates are not substitutes for the
focused development suite.

Test expectations
-----------------

Behavior changes need edge and failure cases, permitted aliasing cases,
mathematical parent/lifetime cases, and explicit proof or certification
status.  Prefer known mathematical results and differential checks against the
source identified in :doc:`source_fidelity`.  Keep public API tests, internal
algorithm tests, FLINT-wrapper ownership tests, adapter tests, diagnostics
tests, and performance measurements conceptually separate.

Documentation changes must run:

.. code-block:: console

   python3 -m pip install -r docs/requirements.txt
   cmake --preset docs
   cmake --build --preset docs
   cmake --build --preset docs-linkcheck

Use direct Sphinx for a channel/version override:

.. code-block:: console

   SILEX_DOCS_RELEASE=0.1.0 SILEX_DOCS_CHANNEL=dev \
     sphinx-build -n -W --keep-going -b html docs build/docs/html

Compiled examples are the executable source of truth for tutorial code.  When
a tutorial needs a focused excerpt, add stable ``literalinclude`` markers to
the compiled example and keep the excerpt inside its tested path.

Review evidence
---------------

A review handoff records changed paths, exact commands and results, source or
provenance implications, performance-impact classification, and remaining
risk.  Review confirms canonical forms, cache invalidation, parent identity,
failure behavior, optional-dependency boundaries, and instrumentation
separation.  Keep chronological logs and raw benchmark artifacts outside the
tracked source tree; publish only durable conclusions.
