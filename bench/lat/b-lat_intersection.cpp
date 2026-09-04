#include <benchmark/benchmark.h>

#include <silex/lat.hpp>

#include "benchmark_contract.hpp"

namespace {

void set_entry_si(fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

void BM_lat_intersection(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::lat::Lat left(2);
    silex::lat::Lat right(2);
    silex::lat::Lat intersection(2);
    fmpz_mat_t left_basis;
    fmpz_mat_t right_basis;

    fmpz_mat_init(left_basis, 2, 2);
    fmpz_mat_init(right_basis, 2, 2);

    set_entry_si(left_basis, 0, 0, 1);
    set_entry_si(left_basis, 0, 1, 1);
    set_entry_si(left_basis, 1, 0, 1);
    set_entry_si(left_basis, 1, 1, -1);
    if (!left.set_basis(left_basis)) {
        silex::bench_contract::fail(
                state, "left lattice setup failed",
                silex::bench_contract::FailureReason::setup);
        fmpz_mat_clear(right_basis);
        fmpz_mat_clear(left_basis);
        return;
    }

    set_entry_si(right_basis, 0, 0, 2);
    set_entry_si(right_basis, 1, 1, 2);
    if (!right.set_basis(right_basis)) {
        silex::bench_contract::fail(
                state, "right lattice setup failed",
                silex::bench_contract::FailureReason::setup);
        fmpz_mat_clear(right_basis);
        fmpz_mat_clear(left_basis);
        return;
    }

    bool operation_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = left.intersection(intersection, right);
        operation_ok = iteration_ok && operation_ok;
        benchmark::DoNotOptimize(intersection.nrows());
    }

    benchmark::ClobberMemory();
    const bool valid = operation_ok && intersection.nrows() == 2 &&
                       fmpz_mat_equal(intersection.basis_ref().raw(),
                                      right.basis_ref().raw()) != 0;
    if (!valid) {
        silex::bench_contract::fail(
                state, "lattice intersection failed its exact reference check",
                operation_ok
                        ? silex::bench_contract::FailureReason::reference_mismatch
                        : silex::bench_contract::FailureReason::operation);
        fmpz_mat_clear(right_basis);
        fmpz_mat_clear(left_basis);
        return;
    }
    state.counters["rows"] = static_cast<double>(intersection.nrows());
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(right_basis);
    fmpz_mat_clear(left_basis);
}

}  // namespace

BENCHMARK(BM_lat_intersection);
BENCHMARK_MAIN();
