#include <benchmark/benchmark.h>

#include "../benchmark_contract.hpp"
#include "../../test/ideal_product_fixtures.hpp"

namespace {
namespace contract = silex::bench_contract;

void BM_ideal_multiply(benchmark::State& state) {
    contract::initialize(state);
    const slong degree = state.range(0);
    const slong kind = state.range(1);
    silex::Order order;
    silex::Ideal left, right;
    if (!silex::test::ideal_product::fixture(order, left, right, degree, kind)) {
        contract::fail(state, "ideal fixture setup failed", contract::FailureReason::setup);
        return;
    }
    silex::Ideal expected(order), result(order);
    if (!silex::test::ideal_product::exhaustive_product(expected, left, right)) {
        contract::fail(state, "exact product oracle failed", contract::FailureReason::setup);
        return;
    }
    bool ok = true;
    for (auto _ : state) {
        ok = result.multiply(left, right) && ok;
        benchmark::DoNotOptimize(result);
    }
    benchmark::ClobberMemory();
    if (!ok || !result.equal(expected)) {
        contract::fail(state, "product lattice differs", contract::FailureReason::invariant);
        return;
    }
    state.counters["degree"] = static_cast<double>(degree);
    state.counters["fixture_kind"] = static_cast<double>(kind);
    state.counters["exact_product"] = 1;
    contract::succeed(state);
}

void multiplication_arguments(benchmark::Benchmark* bench) {
    for (slong degree : {2, 3, 4, 6, 8, 12}) {
        for (slong kind = 0; kind < 5; ++kind) bench->Args({degree, kind});
    }
    bench->Args({2, 5});
}

}  // namespace

BENCHMARK(BM_ideal_multiply)->Apply(multiplication_arguments);
BENCHMARK_MAIN();
