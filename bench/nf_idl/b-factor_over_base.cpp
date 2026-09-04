#include <benchmark/benchmark.h>

#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include "../benchmark_contract.hpp"

namespace {
namespace sflint = silex::flint;
namespace contract = silex::bench_contract;

bool set_quadratic_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);
    return true;
}

void BM_ideal_factor_over_base_quadratic(benchmark::State& state) {
    contract::initialize(state);

    sflint::FmpqPoly polynomial;
    if (!set_quadratic_polynomial(polynomial)) {
        contract::fail(state, "quadratic polynomial setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::NumberField field;
    if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial))) {
        contract::fail(state, "quadratic field setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Order equation_order;
    silex::Order maximal_order(field);
    if (!equation_order.define_equation_order(field) ||
        !maximal_order.maximal_order(equation_order)) {
        contract::fail(state, "maximal order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz bound;
    sflint::fmpz_set_si(bound, 7);

    silex::FactorBase base(maximal_order);
    if (!base.build(sflint::FmpzConstRef(bound))) {
        contract::fail(state, "factor base setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::OrderElement generator(maximal_order);
    silex::Ideal ideal(maximal_order);
    if (!generator.set_si(7) || !ideal.set_principal(generator)) {
        contract::fail(state, "principal ideal setup failed",
                       contract::FailureReason::setup);
        return;
    }

    sflint::FmpzMat exponents(1, base.length());
    bool all_operations_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = silex::ideal_factor_over_base(
                sflint::FmpzMatRef(exponents), ideal, base);
        all_operations_ok = iteration_ok && all_operations_ok;
        benchmark::DoNotOptimize(sflint::fmpz_mat_ncols(exponents));
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed ideal factorization over base failed",
                       contract::FailureReason::operation);
        return;
    }
    if (!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents), ideal,
                                       base)) {
        contract::fail(state, "ideal factorization over base failed",
                       contract::FailureReason::operation);
        return;
    }

    state.counters["factor_base_length"] =
            static_cast<double>(base.length());
    const bool factorization_matches =
            sflint::fmpz_mat_nrows(exponents) == 1 &&
            sflint::fmpz_mat_ncols(exponents) == 5 &&
            sflint::fmpz_equal_si(sflint::fmpz_mat_entry(exponents, 0, 0),
                                  0) &&
            sflint::fmpz_equal_si(sflint::fmpz_mat_entry(exponents, 0, 1),
                                  0) &&
            sflint::fmpz_equal_si(sflint::fmpz_mat_entry(exponents, 0, 2),
                                  0) &&
            sflint::fmpz_equal_si(sflint::fmpz_mat_entry(exponents, 0, 3),
                                  1) &&
            sflint::fmpz_equal_si(sflint::fmpz_mat_entry(exponents, 0, 4),
                                  1);
    if (base.length() != 5 || !factorization_matches) {
        contract::fail(state, "ideal factorization invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_ideal_factor_over_base_quadratic);
BENCHMARK_MAIN();
