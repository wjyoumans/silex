#include <benchmark/benchmark.h>

#include "fmpz_smat_bench_data.hpp"
#include "benchmark_contract.hpp"

namespace {

using silex::bench::fmpz_smat::set_sparse_kernel_matrix;
using silex::bench::fmpz_smat::set_sparse_kernel_right;
using silex::fmpz_smat::SparseMat;

void BM_fmpz_smat_mul_fmpz_mat(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    SparseMat matrix(4);
    fmpz_mat_t right;
    fmpz_mat_t out;
    fmpz_mat_t dense;
    fmpz_mat_t expected;

    fmpz_mat_init(right, 4, 2);
    fmpz_mat_init(out, 3, 2);
    fmpz_mat_init(dense, 3, 4);
    fmpz_mat_init(expected, 3, 2);

    set_sparse_kernel_matrix(matrix);
    set_sparse_kernel_right(right);
    matrix.get_fmpz_mat(dense);
    fmpz_mat_mul(expected, dense, right);

    for (auto _ : state) {
        matrix.mul_fmpz_mat(out, right);
        benchmark::DoNotOptimize(fmpz_mat_entry(out, 0, 0));
    }

    benchmark::ClobberMemory();
    if (fmpz_mat_equal(out, expected) == 0) {
        silex::bench_contract::fail(
                state, "sparse-dense product differs from FLINT reference",
                silex::bench_contract::FailureReason::reference_mismatch);
        fmpz_mat_clear(expected);
        fmpz_mat_clear(dense);
        fmpz_mat_clear(out);
        fmpz_mat_clear(right);
        return;
    }
    state.counters["rows"] = static_cast<double>(fmpz_mat_nrows(out));
    state.counters["columns"] = static_cast<double>(fmpz_mat_ncols(out));
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(expected);
    fmpz_mat_clear(dense);
    fmpz_mat_clear(out);
    fmpz_mat_clear(right);
}

void BM_fmpz_smat_transpose(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    SparseMat matrix(4);
    SparseMat transposed(0);
    fmpz_mat_t dense;
    fmpz_mat_t expected;
    fmpz_mat_t observed;

    fmpz_mat_init(dense, 3, 4);
    fmpz_mat_init(expected, 4, 3);
    fmpz_mat_init(observed, 4, 3);

    set_sparse_kernel_matrix(matrix);
    matrix.get_fmpz_mat(dense);
    fmpz_mat_transpose(expected, dense);

    for (auto _ : state) {
        matrix.transpose(transposed);
        benchmark::DoNotOptimize(transposed.nrows());
    }

    benchmark::ClobberMemory();
    transposed.get_fmpz_mat(observed);
    const bool valid = transposed.nrows() == 4 &&
                       transposed.ncols() == 3 &&
                       transposed.nnz() == matrix.nnz() &&
                       fmpz_mat_equal(observed, expected) != 0;
    if (!valid) {
        silex::bench_contract::fail(
                state, "sparse transpose differs from FLINT reference",
                silex::bench_contract::FailureReason::reference_mismatch);
        fmpz_mat_clear(observed);
        fmpz_mat_clear(expected);
        fmpz_mat_clear(dense);
        return;
    }
    state.counters["rows"] = static_cast<double>(transposed.nrows());
    state.counters["columns"] = static_cast<double>(transposed.ncols());
    state.counters["nonzeros"] = static_cast<double>(transposed.nnz());
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(observed);
    fmpz_mat_clear(expected);
    fmpz_mat_clear(dense);
}

}  // namespace

BENCHMARK(BM_fmpz_smat_mul_fmpz_mat);
BENCHMARK(BM_fmpz_smat_transpose);

BENCHMARK_MAIN();
