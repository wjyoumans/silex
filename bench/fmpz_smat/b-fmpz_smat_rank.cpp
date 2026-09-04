#include <benchmark/benchmark.h>

#include "fmpz_smat_bench_data.hpp"
#include "benchmark_contract.hpp"

#include <array>

namespace {

using silex::bench::fmpz_smat::set_mod_rank_context_rows;
using silex::bench::fmpz_smat::set_rank_matrix;
using silex::fmpz_smat::ModRankContext;
using silex::fmpz_smat::SparseMat;

void BM_fmpz_smat_rank_mod_prime_ui(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    SparseMat matrix(2);
    fmpz_mat_t dense;
    slong rank = -1;

    fmpz_mat_init(dense, 2, 2);
    set_rank_matrix(dense);

    matrix.set_fmpz_mat(dense);

    bool all_operations_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = matrix.rank_mod_prime_ui(&rank, 3);
        all_operations_ok = iteration_ok && all_operations_ok;
        benchmark::DoNotOptimize(rank);
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        silex::bench_contract::fail(
                state, "timed sparse modular-rank operation failed",
                silex::bench_contract::FailureReason::operation);
        fmpz_mat_clear(dense);
        return;
    }
    rank = -1;
    if (!matrix.rank_mod_prime_ui(&rank, 3) || rank != 2) {
        silex::bench_contract::fail(
                state, "sparse modular rank result is not two",
                silex::bench_contract::FailureReason::invariant);
        fmpz_mat_clear(dense);
        return;
    }
    state.counters["rank"] = static_cast<double>(rank);
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(dense);
}

void BM_fmpz_smat_mod_rank_add_rows(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    ModRankContext context;
    fmpz_mat_t dense;
    bool independent = false;

    fmpz_mat_init(dense, 5, 5);
    set_mod_rank_context_rows(dense);

    const bool ok = context.set_prime_ui(5, 7);

    bool all_operations_ok = true;
    if (!ok) {
        silex::bench_contract::fail(
                state, "fmpz_smat modular rank context setup failed",
                silex::bench_contract::FailureReason::setup);
    } else {
        for (auto _ : state) {
            context.reset();
            for (slong row = 0; row < 5; ++row) {
                const bool iteration_ok =
                        context.add_fmpz_mat_row(&independent, dense, row);
                all_operations_ok = iteration_ok && all_operations_ok;
            }
            benchmark::DoNotOptimize(context.rank());
        }
    }

    benchmark::ClobberMemory();
    if (ok) {
        if (!all_operations_ok) {
            silex::bench_contract::fail(
                    state, "timed modular-rank row insertion failed",
                    silex::bench_contract::FailureReason::operation);
            fmpz_mat_clear(dense);
            return;
        }
        context.reset();
        constexpr std::array<bool, 5> expected_independence = {
                true, true, true, true, true};
        bool valid = true;
        for (slong row = 0; row < 5; ++row) {
            independent = false;
            valid = context.add_fmpz_mat_row(&independent, dense, row) &&
                    independent == expected_independence[static_cast<std::size_t>(row)] &&
                    valid;
        }
        valid = valid && context.rank() == 5;
        if (!valid) {
            silex::bench_contract::fail(
                    state, "modular rank row sequence failed validation",
                    silex::bench_contract::FailureReason::invariant);
            fmpz_mat_clear(dense);
            return;
        }
        state.counters["rank"] = static_cast<double>(context.rank());
        silex::bench_contract::succeed(state);
    }
    fmpz_mat_clear(dense);
}

}  // namespace

BENCHMARK(BM_fmpz_smat_rank_mod_prime_ui);
BENCHMARK(BM_fmpz_smat_mod_rank_add_rows);

BENCHMARK_MAIN();
