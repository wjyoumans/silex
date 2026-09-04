#include <benchmark/benchmark.h>

#include "fmpz_smat_bench_data.hpp"
#include "benchmark_contract.hpp"

#include <array>

namespace {

using silex::bench::fmpz_smat::set_hnf_context_rows;
using silex::fmpz_smat::HnfContext;

void BM_fmpz_smat_hnf_ctx_add_rows(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    HnfContext context;
    fmpz_mat_t dense;
    fmpz_mat_t hnf;
    fmpz_t index;
    bool independent = false;

    fmpz_mat_init(dense, 3, 2);
    fmpz_mat_init(hnf, 2, 2);
    fmpz_init(index);
    set_hnf_context_rows(dense);

    bool all_operations_ok = true;
    for (auto _ : state) {
        const bool reset_ok = context.reset(2, 0);
        all_operations_ok = reset_ok && all_operations_ok;
        for (slong row = 0; row < 3; ++row) {
            const bool add_ok =
                    context.add_fmpz_mat_row(&independent, dense, row);
            all_operations_ok = add_ok && all_operations_ok;
        }
        const bool index_ok = context.full_rank_index(index);
        const bool hnf_ok = context.get_hnf(hnf);
        all_operations_ok = index_ok && hnf_ok && all_operations_ok;
        benchmark::DoNotOptimize(fmpz_mat_entry(hnf, 0, 0));
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        silex::bench_contract::fail(
                state, "timed incremental HNF operation failed",
                silex::bench_contract::FailureReason::operation);
        fmpz_clear(index);
        fmpz_mat_clear(hnf);
        fmpz_mat_clear(dense);
        return;
    }
    constexpr std::array<bool, 3> expected_independence = {true, true, false};
    bool valid = context.reset(2, 0);
    for (slong row = 0; row < 3; ++row) {
        independent = false;
        valid = context.add_fmpz_mat_row(&independent, dense, row) &&
                independent == expected_independence[static_cast<std::size_t>(row)] &&
                valid;
    }
    valid = valid && context.rank() == 2 && context.full_rank_index(index) &&
            fmpz_equal_si(index, 10) != 0 && context.get_hnf(hnf) &&
            fmpz_equal_si(fmpz_mat_entry(hnf, 0, 0), 2) != 0 &&
            fmpz_is_zero(fmpz_mat_entry(hnf, 0, 1)) != 0 &&
            fmpz_is_zero(fmpz_mat_entry(hnf, 1, 0)) != 0 &&
            fmpz_equal_si(fmpz_mat_entry(hnf, 1, 1), 5) != 0;
    if (!valid) {
        silex::bench_contract::fail(
                state, "incremental HNF context failed exact validation",
                silex::bench_contract::FailureReason::invariant);
        fmpz_clear(index);
        fmpz_mat_clear(hnf);
        fmpz_mat_clear(dense);
        return;
    }
    state.counters["rank"] = static_cast<double>(context.rank());
    state.counters["index"] = fmpz_get_d(index);
    silex::bench_contract::succeed(state);
    fmpz_clear(index);
    fmpz_mat_clear(hnf);
    fmpz_mat_clear(dense);
}

}  // namespace

BENCHMARK(BM_fmpz_smat_hnf_ctx_add_rows);

BENCHMARK_MAIN();
