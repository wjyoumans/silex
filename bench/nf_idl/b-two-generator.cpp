#include <benchmark/benchmark.h>

#include "../benchmark_contract.hpp"
#include "../../test/ideal_product_fixtures.hpp"
#include "../../src/ideal/two_generator_internal.hpp"

namespace {
namespace contract = silex::bench_contract;

void BM_ideal_two_generator(benchmark::State& state) {
    contract::initialize(state);
    const slong degree = state.range(0), operation = state.range(1);
    silex::Order order;
    silex::Ideal left, right;
    const bool setup = operation == 2
            ? silex::test::ideal_product::random_search_fixture(order, right)
            : silex::test::ideal_product::fixture(order, left, right, degree, 3);
    if (!setup) {
        contract::fail(state, "component fixture failed", contract::FailureReason::setup);
        return;
    }
    silex::flint::Fmpz scalar;
    silex::OrderElement beta(order);
    silex::Ideal result(order), expected(order);
    silex::detail::TwoGeneratorSearchOptions options;
    options.random_trials = operation == 2 ? 0 : 32;
    if (operation == 1 &&
        (!silex::detail::ideal_two_generator(silex::flint::FmpzRef(scalar), beta, right) ||
         !silex::test::ideal_product::exhaustive_product(expected, left, right))) {
        contract::fail(state, "kernel setup failed", contract::FailureReason::setup);
        return;
    }
    bool ok = true;
    for (auto _ : state) {
        if (operation == 1) {
            ok = silex::detail::multiply_integral_ideal_by_two_generator(
                    result, left, silex::flint::FmpzConstRef(scalar), beta) && ok;
        } else {
            const bool found = silex::detail::ideal_two_generator(
                    silex::flint::FmpzRef(scalar), beta, right, options);
            ok = (found == (operation == 0)) && ok;
        }
        benchmark::DoNotOptimize(ok);
    }
    if (operation == 0) {
        ok = silex::detail::set_known_two_generator_ideal(
                result, silex::flint::FmpzConstRef(scalar), beta) && result.equal(right) && ok;
    } else if (operation == 1) {
        ok = result.equal(expected) && ok;
    }
    if (!ok) {
        contract::fail(state, "component postcondition failed", contract::FailureReason::invariant);
        return;
    }
    state.counters["operation"] = static_cast<double>(operation);
    contract::succeed(state);
}

void component_arguments(benchmark::Benchmark* bench) {
    for (slong degree : {2, 3, 4, 6, 8, 12}) {
        bench->Args({degree, 0});
        bench->Args({degree, 1});
    }
    bench->Args({3, 2});
}
}  // namespace

BENCHMARK(BM_ideal_two_generator)->Apply(component_arguments);
BENCHMARK_MAIN();
