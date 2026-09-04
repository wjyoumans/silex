# Silex examples

Examples are small native C++ programs that exercise implemented module
surfaces without raw FLINT handles in user code.

Number fields and orders in examples use factory construction, such as
`NumberField::by_polynomial`, `NumberField::quadratic`, and
`Order::equation_order`.  Mutating `define_*` methods remain compatibility and
scratch-object helpers, but they are not the ordinary example style.

Examples prefer owned-return helpers such as `Element::to_fmpq_poly`,
`ClassGroupContext::order`, `OrderUnitGroup::torsion_order`,
`OrderUnitGroup::regulator`, and zeta BF audit result records for ordinary
queries.  Explicit output buffers still appear where the example is exercising
matrix/log computations, proof-driver APIs, or operations without an
owned-return helper.

- `factored_compact_roots.cpp`: factored number-field elements, compact
  presentations, exact structural roots, and compact root strategies.
- `field_order_ideal_basics.cpp`: C++-first construction of a number field,
  equation order, principal integral ideal from an `OrderElement`, and
  principal fractional ideal from an ambient `Element`.
- `element_arithmetic_basics.cpp`: exact `Element` arithmetic, trace/norm,
  inversion, and signed-power calls in a quadratic field.
- `log_unit_lattice.cpp`: embeddings, logarithmic unit matrices, regulator
  computation, and independence checks.
- `local_relation_basics.cpp`: prime decomposition, copied factor-base primes,
  rational-prime block access, relation rows/matrices, and residue-field
  arithmetic through C++ domain objects.
- `maximal_order_zeta.cpp`: a cubic maximal order, order index and
  discriminant access, and Belabas-Friedman zeta class-regulator audit output
  as a named result record.
- `class_unit_computation.cpp`: paired class-group and order-unit transaction
  through the native C++ API, including coarse
  certification and component proof-status queries.
- `diagnostics_hot_paths.cpp`: callback-based verbose/log/debug/profile output
  for compact-root, relation, saturation, class-group, and unit-group hot
  paths.  Enable channels independently with `SILEX_DIAGNOSTICS_VERBOSE`,
  `SILEX_DIAGNOSTICS_LOG`, `SILEX_DIAGNOSTICS_DEBUG`, and
  `SILEX_DIAGNOSTICS_PROFILE`.
- `rank_zero_order_units.cpp`: proven rank-zero order-unit computation and
  torsion access.
- `real_quadratic_order_units.cpp`: exact continued-fraction fundamental-unit
  publication for a release-qualified maximal real-quadratic example through
  `OrderUnitGroup`.
- `relation_kernel_order_units.cpp`: relation-kernel compact unit extraction
  from class-group relation data.
- `supplied_order_units.cpp`: caller-supplied compact order-unit subgroup
  validation, regulator storage, and certification metadata.

## Diagnostics example

Build a diagnostics-enabled tree when you want library-side debug checks or
profile scopes:

```sh
cmake -S . -B build/diagnostics \
  -DSILEX_BUILD_EXAMPLES=ON \
  -DSILEX_ENABLE_LOGGING=ON \
  -DSILEX_ENABLE_DEBUG_CHECKS=ON \
  -DSILEX_ENABLE_PROFILING=ON
cmake --build build/diagnostics --target example-diagnostics-hot-paths
```

Run channels independently:

```sh
SILEX_DIAGNOSTICS_VERBOSE=progress build/diagnostics/examples/example-diagnostics-hot-paths
SILEX_DIAGNOSTICS_LOG=detail build/diagnostics/examples/example-diagnostics-hot-paths
SILEX_DIAGNOSTICS_DEBUG=cheap build/diagnostics/examples/example-diagnostics-hot-paths
SILEX_DIAGNOSTICS_PROFILE=1 build/diagnostics/examples/example-diagnostics-hot-paths
```

When profiling is enabled, scope-end lines include elapsed time:

```text
[profile:relation:end] set_generator: relation.set_generator elapsed_us=128.843
profile inclusive scope time = 10.0482 ms
```

The summary is an inclusive sum of completed profile scopes, so nested scopes
are counted in both the child and parent timing. Use it for quick callback
inspection; use Google Benchmark or an external profiler for stable
performance measurements.

The example also accepts combinations, for example:

```sh
SILEX_DIAGNOSTICS_VERBOSE=progress \
SILEX_DIAGNOSTICS_LOG=detail \
SILEX_DIAGNOSTICS_DEBUG=cheap \
SILEX_DIAGNOSTICS_PROFILE=1 \
build/diagnostics/examples/example-diagnostics-hot-paths
```
