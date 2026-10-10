#include <benchmark/benchmark.h>

#include <flint/fmpq.h>

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include "benchmark_contract.hpp"
#include "zeta/zeta_internal.hpp"

// Times the GRH generation bound min(BDF, Bach) on inputs where the
// Belabas-Diaz y Diaz-Friedman criterion cannot be computed, so that
// grh_factor_base_bound_with_diagnostics falls back to Bach's bound alone
// (Bach 1990, Theorem 4).  The fallback fixture is the non-monic
// presentation 4x^3 + x + 100 of the field of x^3 + x + 200
// (|D| = 1080004, Bach bound 2317); the control is the monic presentation,
// where the BDF search finds 19.  Each row checks the path it takes
// (bach_selected) and the bound before timing.

namespace {
namespace sflint = silex::flint;

struct Fixture {
    silex::NumberField field;
    silex::Order maximal;
};

bool monic_fixture(Fixture& out) {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 200);
    out.field = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    if (!out.field.is_defined()) {
        return false;
    }
    silex::Order equation = silex::Order::equation_order(out.field);
    out.maximal = silex::Order(out.field);
    return equation.is_defined() && out.maximal.maximal_order(equation);
}

bool nonmonic_fixture(Fixture& out) {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 100);
    out.field = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    if (!out.field.is_defined()) {
        return false;
    }
    // Z[2 beta] = Z[alpha], with basis 1, 2 beta, 4 beta^2.
    sflint::FmpqMat basis(3, 3);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 1, 1), 2, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 2, 2), 4, 1);
    silex::Order alpha = silex::Order::from_basis(
            out.field, sflint::FmpqMatConstRef(basis));
    out.maximal = silex::Order(out.field);
    return alpha.is_defined() && out.maximal.maximal_order(alpha);
}

void run(benchmark::State& state, bool monic) {
    silex::bench_contract::initialize(state);
    Fixture fixture;
    if (!(monic ? monic_fixture(fixture) : nonmonic_fixture(fixture))) {
        silex::bench_contract::fail(
                state, "fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    const slong expected = monic ? 19 : 2317;
    sflint::Fmpz bound;
    bool bach_selected = monic;
    if (!silex::detail::grh_factor_base_bound_with_diagnostics(
                sflint::FmpzRef(bound), fixture.maximal, nullptr,
                &bach_selected) ||
        !sflint::fmpz_equal_si(bound, expected) ||
        bach_selected == monic) {
        silex::bench_contract::fail(
                state, "bound or fallback path differs from the reference",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    for (auto _ : state) {
        bool selected = false;
        bool ok = silex::detail::grh_factor_base_bound_with_diagnostics(
                sflint::FmpzRef(bound), fixture.maximal, nullptr, &selected);
        benchmark::DoNotOptimize(ok);
        benchmark::DoNotOptimize(selected);
        if (!ok || selected == monic) {
            silex::bench_contract::fail(
                    state, "path changed during timing",
                    silex::bench_contract::FailureReason::invariant);
            return;
        }
    }
    state.counters["bach_selected"] = bach_selected ? 1.0 : 0.0;
    state.counters["bound"] = static_cast<double>(expected);
    silex::bench_contract::succeed(state);
}

void BM_grh_factor_base_bound_bach_fallback(benchmark::State& state) {
    run(state, false);
}

void BM_grh_factor_base_bound_bdf_control(benchmark::State& state) {
    run(state, true);
}

BENCHMARK(BM_grh_factor_base_bound_bach_fallback);
BENCHMARK(BM_grh_factor_base_bound_bdf_control);

}  // namespace

BENCHMARK_MAIN();
