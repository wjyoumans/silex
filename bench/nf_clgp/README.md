# nf_clgp benchmarks

This directory contains native class-group and paired class/unit benchmark
rows. The main target uses the public C++ API; explicitly named internal
targets may include implementation headers to isolate a source-backed
checkpoint without exposing it as stable API.

The `b-silex-nf_clgp` target includes:

- deterministic cubic, quartic, quintic, and bounded sextic candidate rows;
- scaled factor-base variants for selected fields;
- degree-one and real/imaginary quadratic paired class/unit rows;
- deterministic random-sweep, random-matrix, and boundary workloads; and
- public counters for success, factor-base size, relation counts and rank,
  relation sources, kernel witnesses, class order, and proof metadata.

Private diagnostic counters must remain in separately named internal targets.
They must use existing internal boundaries, publish no production result, and
remain outside the stable API and default CTest suite.

The internal `b-silex-nf_clgp_factor_base_honesty` target compares exact
witness predicates on frozen quartic and quintic batches (672 and 3937
candidates). `BM_honesty_selected_*` measures the selected predicate, including
its fallback. Every timed candidate must match its frozen classification,
not merely the final match count. The `--dump-witness-fixtures` development
mode emits source-neutral JSON lines containing polynomial coefficients,
factor-base and required-prime lattices in the polynomial power basis, and
each candidate's expected classification for an independent exact oracle.

## Reproduce a field

Run an exact failing or slow field with a process timeout before promoting it
to a benchmark:

```sh
python tools/bench/run-class-unit-instance.py \
  --field-id cubic_disc81_proven \
  --timeout 60 \
  --build-dir build/profile-release
```

Use the fast manifest gate for rows marked `must_pass_fast`:

```sh
python tools/bench/run-class-unit-priority-gate.py \
  --build-dir build/profile-release
```

For a bounded backtrace of a stuck process:

```sh
python tools/bench/gdb-class-unit-timeout.py \
  --field-id cubic_disc81_proven \
  --build-dir build/profile-release \
  --timeout 60
```

Promote a field to normal tests when it must pass inside the one-minute test
ceiling. Promote it to this benchmark only when it is stable performance
coverage rather than a transient investigation target.

## Noninstalled adapter

The honesty witness target includes `BM_honesty_scalar_selected_*` and
`BM_honesty_search_selected_*` rows. They respectively time the frozen scalar
prefix and the mixed scalar/enumerated batch through the search representation
boundary. The existing `BM_honesty_selected_*` rows retain the general
order-element path as a control. Setup checks each scalar against the exact
scalar kernel and reference classifications; positive and failure cases are
covered by the native factor-base honesty test.

`BM_honesty_lattice_*` times a batch of `Lat::lll_reduce` calls on all searched
prime ideals from the same quartic/quintic fixtures. Field construction,
decomposition, input preparation, witness validation, and reference checks
are outside timing. Each reduced basis must equal the same-build reference
with a zero transformation argument, an identity-initialized transformation,
and no transformation argument. Setup also checks the unimodular row identity,
canonical HNF, and reducedness independently. The `lattices` counter records
the number of reductions per batch; these are not per-lattice timing rows.

The `silex-class-unit-instance` executable is a source-tree development tool.
It is built when tests are enabled or
`SILEX_BUILD_BENCHMARK_ADAPTERS=ON`, and is not part of the installed package.
Its JSON output records the requested and published certification labels,
proof metadata, failure stage, publication state, component timings, and
transaction-phase timings.

Only compact, source-neutral correctness and certification fixtures belong in
`test/data/class_unit_fields.json`. External-engine corpora and campaign
configuration do not belong in this repository.

Validation uses focused CTest fixtures, the priority gate, and task-local
Google Benchmark JSON. Historical run narratives and raw campaign artifacts
are not part of the public source tree.
