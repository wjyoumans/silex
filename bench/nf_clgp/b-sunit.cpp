// In-process S-unit rows (BM_sunit_*).  Each row times only
// compute_sunit_groups; the class group, the ordinary units and the prime
// selection are built once in setup, outside the timed region.
//
// The fields are the "index": "all" rows of test/data/sunit_fields.json that
// take under a second: cubic_x3_minus_2_regression, real_quadratic_5_split_11
// and real_quadratic_5_ramified_5, plus cubic_x3_minus_2_empty_s as the
// empty-S control.  "all" rows select every prime ideal above each listed
// rational prime, so no witness matching is needed.  The two manifest rows
// that select a single prime ideal by witness (imaginary_quadratic_23 and
// real_quadratic_210) are left to tools/bench/run-sunit-instance.py.
// Polynomials and primes are copied from the manifest; the expected S-class
// order and ranks are checked once per row against its "expected" block.
#include <benchmark/benchmark.h>

#include <silex/class_group.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/sunit.hpp>

#include "benchmark_contract.hpp"
#include "../../test/test_support.hpp"

#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;
namespace bench = silex::bench_contract;

struct SUnitCase {
    std::initializer_list<slong> coefficients_low_to_high;
    std::initializer_list<slong> rational_primes;
    slong s_class_order;
    slong ordinary_free_rank;
    slong nonunit_rank;
};

struct SUnitFixture {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    std::vector<silex::PrimeIdeal> selected;
    bool ok = false;
};

bool build_fixture(SUnitFixture& out, const SUnitCase& c) {
    sflint::FmpqPoly polynomial;
    slong degree = 0;
    for (const slong coefficient : c.coefficients_low_to_high) {
        sflint::fmpq_poly_set_coeff_si(polynomial, degree++, coefficient);
    }
    out.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    if (!out.order.maximal_order(equation) || !out.order.is_maximal()) {
        return false;
    }
    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    sflint::Fmpz bound;
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(bound), out.order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    if (!out.units.compute_with_class_group(
                out.class_group, out.order, sflint::FmpzConstRef(bound),
                options, 192)) {
        return false;
    }
    for (const slong rational_prime : c.rational_primes) {
        sflint::Fmpz p;
        sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
        silex::PrimeIdealList decomposition;
        if (!silex::decompose_prime(
                    decomposition, out.order, sflint::FmpzConstRef(p))) {
            return false;
        }
        for (slong i = 0; i < decomposition.size(); ++i) {
            const silex::PrimeIdeal* source = decomposition.at(i);
            out.selected.emplace_back(out.order);
            if (source == nullptr || !out.selected.back().set(*source)) {
                return false;
            }
        }
    }
    return true;
}

bool matches_expected(const SUnitCase& c,
                      const silex::SClassGroup& s_class_group,
                      const silex::SUnitGroup& s_unit_group) {
    sflint::Fmpz order;
    return s_class_group.order(sflint::FmpzRef(order)) &&
           sflint::fmpz_equal_si(order, c.s_class_order) &&
           s_unit_group.ordinary_free_rank() == c.ordinary_free_rank &&
           s_unit_group.nonunit_rank() == c.nonunit_rank;
}

void run_sunit(benchmark::State& state, const SUnitCase& c) {
    bench::initialize(state);
    SUnitFixture fixture;
    if (!build_fixture(fixture, c)) {
        bench::fail(state, "fixture setup failed", bench::FailureReason::setup);
        return;
    }
    const silex::PrimeIdealSpan selected(
            fixture.selected.data(),
            static_cast<std::size_t>(fixture.selected.size()));
    silex::SUnitComputeOptions options;
    options.regulator_precision = 192;

    // One untimed semantic gate before the timed loop.
    {
        silex::SUnitComputeResult result;
        silex::SClassGroup s_class_group;
        silex::SUnitGroup s_unit_group;
        if (!silex::compute_sunit_groups(
                    result, s_class_group, s_unit_group,
                    fixture.class_group, fixture.units, selected, options)) {
            bench::fail(state, "compute_sunit_groups failed",
                        bench::FailureReason::operation);
            return;
        }
        if (!matches_expected(c, s_class_group, s_unit_group)) {
            bench::fail(state, "S-unit result differs from the manifest",
                        bench::FailureReason::reference_mismatch);
            return;
        }
        state.counters["selected_primes"] =
                static_cast<double>(s_unit_group.selected_prime_count());
        state.counters["free_rank"] =
                static_cast<double>(s_unit_group.free_rank());
    }

    for (auto _ : state) {
        silex::SUnitComputeResult result;
        silex::SClassGroup s_class_group;
        silex::SUnitGroup s_unit_group;
        const bool ok = silex::compute_sunit_groups(
                result, s_class_group, s_unit_group, fixture.class_group,
                fixture.units, selected, options);
        benchmark::DoNotOptimize(s_unit_group);
        if (!ok) {
            bench::fail(state, "compute_sunit_groups failed",
                        bench::FailureReason::operation);
            return;
        }
    }
    bench::succeed(state);
}

// cubic_x3_minus_2_empty_s
void BM_sunit_cubic_x3_minus_2_empty_s(benchmark::State& state) {
    run_sunit(state, {{-2, 0, 0, 1}, {}, 1, 1, 0});
}

// cubic_x3_minus_2_regression
void BM_sunit_cubic_x3_minus_2_primes_2_3_5_7(benchmark::State& state) {
    run_sunit(state, {{-2, 0, 0, 1}, {2, 3, 5, 7}, 1, 1, 5});
}

// real_quadratic_5_split_11
void BM_sunit_real_quadratic_5_split_11(benchmark::State& state) {
    run_sunit(state, {{-1, -1, 1}, {11}, 1, 1, 2});
}

// real_quadratic_5_ramified_5
void BM_sunit_real_quadratic_5_ramified_5(benchmark::State& state) {
    run_sunit(state, {{-1, -1, 1}, {5}, 1, 1, 1});
}

}  // namespace

BENCHMARK(BM_sunit_cubic_x3_minus_2_empty_s);
BENCHMARK(BM_sunit_cubic_x3_minus_2_primes_2_3_5_7);
BENCHMARK(BM_sunit_real_quadratic_5_split_11);
BENCHMARK(BM_sunit_real_quadratic_5_ramified_5);
BENCHMARK_MAIN();
