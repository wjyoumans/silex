Fmpz Sparse Matrix Benchmarks
=============================

This directory contains initial ``fmpz_smat`` benchmark coverage.

These workloads use small sparse-matrix test fixtures through the native C++
API:

- ``BM_fmpz_smat_mul_fmpz_mat`` mirrors the sparse-times-dense multiplication
  matrix from ``t-kernels.c``.
- ``BM_fmpz_smat_transpose`` mirrors the transpose matrix from ``t-kernels.c``.
- ``BM_fmpz_smat_rank_mod_prime_ui`` mirrors the one-shot modular rank matrix
  from ``t-modular.c``.
- ``BM_fmpz_smat_mod_rank_add_rows`` mirrors the incremental modular-rank rows
  from ``t-modular.c``.
- ``BM_fmpz_smat_hnf_ctx_add_rows`` mirrors the dense-row HNF context workload
  from ``t-hnf_ctx.c``.

These are fixed-size smoke/regression kernels.  They are not broad performance
claims and should be extended with larger workloads before optimization
decisions.
