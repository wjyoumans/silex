# nf_ord Benchmarks

This directory contains order-kernel benchmarks and small native C++ coverage
for the generic maximalization path.

Current targets:

- `b-silex-nf_ord_maximal_order`
  - `BM_order_maximal_order_quadratic_fast` covers `Q(sqrt(47))`.
  - `BM_order_maximal_order_generic_quadratic` exercises the generic
    p-radical/multiplier-ring path on `x^2 - 5` without quadratic field
    metadata.
  - `BM_order_maximal_order_stable_cubic` exercises the discriminant-factor
    loop on `x^3 - 2`.
- `b-silex-nf_ord_pmaximal_overorder`
  - `BM_order_pmaximal_overorder_quadratic_fast` covers `Q(sqrt(47))`,
    `p = 2`.
  - `BM_order_pmaximal_overorder_generic_quadratic` exercises the generic
    p-maximal path on `x^2 - 5`, `p = 2`.

The generic-path entries are small coverage targets for the native C++
order/ideal maximalization stack. Use external comparisons only for rows with
matching workloads.

Profiling/logging/debug hooks for these algorithms should be added through a
diagnostics-aware options/context API. Avoid adding global diagnostics state or
inert profile scopes that callers cannot enable.
