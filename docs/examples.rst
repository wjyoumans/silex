Examples
========

The example programs under ``examples/`` are small native C++ programs that
exercise implemented module surfaces without requiring users to manage raw
FLINT storage directly.  They are built when ``SILEX_BUILD_EXAMPLES=ON`` and
are registered as CTest tests when ``SILEX_BUILD_TESTS=ON``.

Build and run them with a configured build tree:

.. code-block:: console

   cmake --build build --target example-log-unit-lattice
   ./build/examples/example-log-unit-lattice

Current examples:

:download:`factored_compact_roots.cpp <../examples/factored_compact_roots.cpp>`
   Factored number-field elements, compact presentations, structural roots,
   and exact square/power checks.

:download:`field_order_ideal_basics.cpp <../examples/field_order_ideal_basics.cpp>`
   A C++-first walkthrough from ``NumberField`` to ``Order``, ``Ideal``, and
   ``FractionalIdeal`` using principal ideals and norm/containment queries.

:download:`element_arithmetic_basics.cpp <../examples/element_arithmetic_basics.cpp>`
   Exact ``Element`` arithmetic in a quadratic field, including trace/norm,
   multiplication, inversion, and signed integer powers through named methods.

:download:`log_unit_lattice.cpp <../examples/log_unit_lattice.cpp>`
   Embeddings, logarithmic unit matrices, regulator computation, and
   independence checks.

:download:`local_relation_basics.cpp <../examples/local_relation_basics.cpp>`
   Prime decomposition, copied factor-base primes, rational-prime block
   access, relation rows, relation matrices, and residue-field arithmetic
   through C++ domain objects.

:download:`maximal_order_zeta.cpp <../examples/maximal_order_zeta.cpp>`
   Cubic maximal orders, order index/discriminant access, and
   Belabas-Friedman zeta class-regulator audit output as a named result
   record.

:download:`class_unit_computation.cpp <../examples/class_unit_computation.cpp>`
   Paired class-group candidate and order-unit subgroup computation through
   native ``ClassGroupContext`` and ``OrderUnitGroup`` objects, including
   coarse certification and component proof-status queries.

:download:`diagnostics_hot_paths.cpp <../examples/diagnostics_hot_paths.cpp>`
   Callback-based verbose, log, debug-check, and profiling output for
   compact-root, relation, saturation, class-group, and unit-group hot paths.
   Runtime switches are documented in
   :download:`examples/README.md <../examples/README.md>` and
   :doc:`getting_started`.

:download:`rank_zero_order_units.cpp <../examples/rank_zero_order_units.cpp>`
   Proven rank-zero order-unit computation and torsion access.

:download:`real_quadratic_order_units.cpp <../examples/real_quadratic_order_units.cpp>`
   Exact continued-fraction fundamental-unit publication for a
   release-qualified maximal real-quadratic example through
   ``OrderUnitGroup``.

:download:`relation_kernel_order_units.cpp <../examples/relation_kernel_order_units.cpp>`
   Relation-kernel compact unit extraction from class-group relation data.

:download:`supplied_order_units.cpp <../examples/supplied_order_units.cpp>`
   Caller-supplied compact order-unit subgroup validation, regulator storage,
   and certification metadata.
