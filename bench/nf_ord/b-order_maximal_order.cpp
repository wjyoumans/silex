#include <benchmark/benchmark.h>

#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>

#include "../benchmark_contract.hpp"

namespace {
namespace sflint = silex::flint;
namespace contract = silex::bench_contract;

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

void BM_order_maximal_order_quadratic_fast(benchmark::State& state) {
    contract::initialize(state);

    sflint::Fmpz radicand;
    sflint::fmpz_set_si(radicand, 47);

    silex::NumberField field;
    if (!field.define_quadratic(sflint::FmpzConstRef(radicand))) {
        contract::fail(state, "quadratic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order order;
    silex::Order maximal(field);
    if (!order.define_equation_order(field)) {
        contract::fail(state, "order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    bool all_operations_ok = true;
    bool all_results_maximal = true;
    for (auto _ : state) {
        const bool iteration_ok = maximal.maximal_order(order);
        const bool iteration_maximal = maximal.is_maximal();
        all_operations_ok = iteration_ok && all_operations_ok;
        all_results_maximal = iteration_maximal && all_results_maximal;
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed quadratic maximal-order call failed",
                       contract::FailureReason::operation);
        return;
    }
    if (!all_results_maximal) {
        contract::fail(state, "timed quadratic result was not maximal",
                       contract::FailureReason::invariant);
        return;
    }
    if (!maximal.maximal_order(order)) {
        contract::fail(state, "quadratic maximal-order computation failed",
                       contract::FailureReason::operation);
        return;
    }

    sflint::Fmpz discriminant;
    if (!maximal.discriminant(sflint::FmpzRef(discriminant))) {
        contract::fail(state, "quadratic maximal-order discriminant failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["discriminant"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(discriminant));
    if (!maximal.maximality_known() || !maximal.is_maximal() ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(discriminant), 188)) {
        contract::fail(state, "quadratic maximal-order invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

void BM_order_maximal_order_generic_quadratic(benchmark::State& state) {
    contract::initialize(state);

    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial))) {
        contract::fail(state, "generic quadratic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order order;
    silex::Order maximal(field);
    if (!order.define_equation_order(field)) {
        contract::fail(state, "order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    bool all_operations_ok = true;
    bool all_results_maximal = true;
    for (auto _ : state) {
        const bool iteration_ok = maximal.maximal_order(order);
        const bool iteration_maximal = maximal.is_maximal();
        all_operations_ok = iteration_ok && all_operations_ok;
        all_results_maximal = iteration_maximal && all_results_maximal;
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(
                state, "timed generic quadratic maximal-order call failed",
                contract::FailureReason::operation);
        return;
    }
    if (!all_results_maximal) {
        contract::fail(state,
                       "timed generic quadratic result was not maximal",
                       contract::FailureReason::invariant);
        return;
    }
    if (!maximal.maximal_order(order)) {
        contract::fail(state,
                       "generic quadratic maximal-order computation failed",
                       contract::FailureReason::operation);
        return;
    }

    sflint::Fmpz discriminant;
    if (!maximal.discriminant(sflint::FmpzRef(discriminant))) {
        contract::fail(state,
                       "generic quadratic maximal-order discriminant failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["discriminant"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(discriminant));
    if (!maximal.maximality_known() || !maximal.is_maximal() ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(discriminant), 5)) {
        contract::fail(state,
                       "generic quadratic maximal-order invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

void BM_order_maximal_order_stable_cubic(benchmark::State& state) {
    contract::initialize(state);

    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial))) {
        contract::fail(state, "cubic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order order;
    silex::Order maximal(field);
    if (!order.define_equation_order(field)) {
        contract::fail(state, "order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    bool all_operations_ok = true;
    bool all_results_maximal = true;
    for (auto _ : state) {
        const bool iteration_ok = maximal.maximal_order(order);
        const bool iteration_maximal = maximal.is_maximal();
        all_operations_ok = iteration_ok && all_operations_ok;
        all_results_maximal = iteration_maximal && all_results_maximal;
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed cubic maximal-order call failed",
                       contract::FailureReason::operation);
        return;
    }
    if (!all_results_maximal) {
        contract::fail(state, "timed cubic result was not maximal",
                       contract::FailureReason::invariant);
        return;
    }
    if (!maximal.maximal_order(order)) {
        contract::fail(state, "cubic maximal-order computation failed",
                       contract::FailureReason::operation);
        return;
    }

    sflint::Fmpz discriminant;
    if (!maximal.discriminant(sflint::FmpzRef(discriminant))) {
        contract::fail(state, "cubic maximal-order discriminant failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["discriminant"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(discriminant));
    if (!maximal.maximality_known() || !maximal.is_maximal() ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(discriminant), -108)) {
        contract::fail(state, "cubic maximal-order invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_order_maximal_order_quadratic_fast);
BENCHMARK(BM_order_maximal_order_generic_quadratic);
BENCHMARK(BM_order_maximal_order_stable_cubic);

BENCHMARK_MAIN();
