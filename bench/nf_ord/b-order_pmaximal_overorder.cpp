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

void BM_order_pmaximal_overorder_quadratic_fast(benchmark::State& state) {
    contract::initialize(state);

    sflint::Fmpz radicand;
    sflint::Fmpz prime;
    sflint::fmpz_set_si(radicand, 47);
    sflint::fmpz_set_ui(prime, 2);

    silex::NumberField field;
    if (!field.define_quadratic(sflint::FmpzConstRef(radicand))) {
        contract::fail(state, "quadratic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order order;
    silex::Order pmaximal(field);
    if (!order.define_equation_order(field)) {
        contract::fail(state, "order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    bool all_operations_ok = true;
    bool all_results_maximal = true;
    for (auto _ : state) {
        const bool iteration_ok = pmaximal.pmaximal_overorder(
                order, sflint::FmpzConstRef(prime));
        const bool iteration_maximal = pmaximal.is_maximal();
        all_operations_ok = iteration_ok && all_operations_ok;
        all_results_maximal = iteration_maximal && all_results_maximal;
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed quadratic p-maximal call failed",
                       contract::FailureReason::operation);
        return;
    }
    if (!all_results_maximal) {
        contract::fail(state, "timed quadratic result was not maximal",
                       contract::FailureReason::invariant);
        return;
    }
    if (!pmaximal.pmaximal_overorder(order,
                                     sflint::FmpzConstRef(prime))) {
        contract::fail(state, "quadratic p-maximal computation failed",
                       contract::FailureReason::operation);
        return;
    }

    sflint::Fmpz discriminant;
    if (!pmaximal.discriminant(sflint::FmpzRef(discriminant))) {
        contract::fail(state, "quadratic p-maximal discriminant failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["discriminant"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(discriminant));
    if (!pmaximal.maximality_known() || !pmaximal.is_maximal() ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(discriminant), 188)) {
        contract::fail(state, "quadratic p-maximal invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

void BM_order_pmaximal_overorder_generic_quadratic(benchmark::State& state) {
    contract::initialize(state);

    sflint::FmpqPoly polynomial;
    sflint::Fmpz prime;
    sflint::Fmpz discriminant;
    poly_x2_minus(polynomial, 5);
    sflint::fmpz_set_ui(prime, 2);

    silex::NumberField field;
    if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial))) {
        contract::fail(state, "generic quadratic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order order;
    silex::Order pmaximal(field);
    if (!order.define_equation_order(field)) {
        contract::fail(state, "order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    bool all_operations_ok = true;
    bool all_discriminants_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = pmaximal.pmaximal_overorder(
                order, sflint::FmpzConstRef(prime));
        const bool iteration_discriminant_ok =
                pmaximal.discriminant(sflint::FmpzRef(discriminant));
        all_operations_ok = iteration_ok && all_operations_ok;
        all_discriminants_ok =
                iteration_discriminant_ok && all_discriminants_ok;
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok || !all_discriminants_ok) {
        contract::fail(
                state,
                !all_operations_ok
                        ? "timed generic quadratic p-maximal call failed"
                        : "timed generic quadratic discriminant call failed",
                contract::FailureReason::operation);
        return;
    }
    if (!pmaximal.pmaximal_overorder(order, sflint::FmpzConstRef(prime)) ||
        !pmaximal.discriminant(sflint::FmpzRef(discriminant))) {
        contract::fail(state, "generic quadratic p-maximal computation failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["discriminant"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(discriminant));
    if (!pmaximal.maximality_known() || !pmaximal.is_maximal() ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(discriminant), 5)) {
        contract::fail(state,
                       "generic quadratic p-maximal invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_order_pmaximal_overorder_quadratic_fast);
BENCHMARK(BM_order_pmaximal_overorder_generic_quadratic);

BENCHMARK_MAIN();
