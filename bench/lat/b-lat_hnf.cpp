#include <benchmark/benchmark.h>

#include <silex/lat.hpp>

#include "benchmark_contract.hpp"

namespace {

void set_entry_si(fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

void BM_lat_hnf(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::lat::Lat lattice(4);
    silex::lat::Lat hnf(4);
    fmpz_mat_t basis;
    fmpz_mat_t expected;

    fmpz_mat_init(basis, 4, 4);
    fmpz_mat_init(expected, 4, 4);

    set_entry_si(basis, 0, 0, 6);
    set_entry_si(basis, 1, 0, 2);
    set_entry_si(basis, 1, 1, 9);
    set_entry_si(basis, 2, 0, 5);
    set_entry_si(basis, 2, 1, 1);
    set_entry_si(basis, 2, 2, 7);
    set_entry_si(basis, 3, 0, 3);
    set_entry_si(basis, 3, 1, 4);
    set_entry_si(basis, 3, 2, 2);
    set_entry_si(basis, 3, 3, 11);
    if (!lattice.set_basis(basis)) {
        silex::bench_contract::fail(
                state, "lattice setup failed",
                silex::bench_contract::FailureReason::setup);
        fmpz_mat_clear(expected);
        fmpz_mat_clear(basis);
        return;
    }
    fmpz_mat_hnf(expected, basis);

    bool operation_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = lattice.hnf(hnf);
        operation_ok = iteration_ok && operation_ok;
        benchmark::DoNotOptimize(hnf.nrows());
    }

    benchmark::ClobberMemory();
    const bool valid = operation_ok && hnf.is_hnf() && hnf.nrows() == 4 &&
                       fmpz_mat_equal(hnf.basis_ref().raw(), expected) != 0;
    if (!valid) {
        silex::bench_contract::fail(
                state, "lattice HNF failed its exact reference check",
                operation_ok
                        ? silex::bench_contract::FailureReason::reference_mismatch
                        : silex::bench_contract::FailureReason::operation);
        fmpz_mat_clear(expected);
        fmpz_mat_clear(basis);
        return;
    }
    state.counters["rows"] = static_cast<double>(hnf.nrows());
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(expected);
    fmpz_mat_clear(basis);
}

}  // namespace

BENCHMARK(BM_lat_hnf);
BENCHMARK_MAIN();
