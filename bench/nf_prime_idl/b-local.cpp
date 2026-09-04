#include <benchmark/benchmark.h>

#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include "../benchmark_contract.hpp"

namespace {
namespace sflint = silex::flint;
namespace contract = silex::bench_contract;

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

bool define_maximal_quadratic(silex::NumberField& field,
                              silex::Order& maximal_order,
                              slong radicand) noexcept {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, radicand);

    silex::Order equation_order;
    return field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial)) &&
           equation_order.define_equation_order(field) &&
           maximal_order.define(field) &&
           maximal_order.maximal_order(equation_order);
}

bool define_marked_equation_quadratic(silex::NumberField& field,
                                      silex::Order& order,
                                      slong radicand) noexcept {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, radicand);
    if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial)) ||
        !order.define_equation_order(field)) {
        return false;
    }
    order.set_maximality(true);
    return true;
}

bool set_one_plus_theta_over_two(silex::Element& element) noexcept {
    sflint::Fmpq half;
    sflint::FmpqPoly polynomial;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, half);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, half);
    return element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

void BM_prime_ideal_reduce_nf_elem_quadratic(benchmark::State& state) {
    contract::initialize(state);

    silex::NumberField field;
    silex::Order maximal_order(field);
    if (!define_maximal_quadratic(field, maximal_order, 5)) {
        contract::fail(state, "maximal order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Element alpha(field);
    silex::OrderElement alpha_order(maximal_order);
    if (!set_one_plus_theta_over_two(alpha) ||
        !alpha_order.set_element(alpha)) {
        contract::fail(state, "order element setup failed",
                       contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz p;
    sflint::fmpz_set_si(p, 11);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, maximal_order,
                                sflint::FmpzConstRef(p)) ||
        primes.size() <= 0 || primes.at(0) == nullptr) {
        contract::fail(state, "prime decomposition failed",
                       contract::FailureReason::setup);
        return;
    }

    const silex::PrimeIdeal& prime = *primes.at(0);
    sflint::FmpzPoly reduced;
    bool all_operations_ok = true;
    for (auto _ : state) {
        const bool iteration_ok =
                prime.reduce(sflint::FmpzPolyRef(reduced), alpha_order);
        all_operations_ok = iteration_ok && all_operations_ok;
        benchmark::DoNotOptimize(fmpz_poly_get_coeff_si(reduced.raw(), 0));
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed prime-ideal reduction failed",
                       contract::FailureReason::operation);
        return;
    }
    if (!prime.reduce(sflint::FmpzPolyRef(reduced), alpha_order)) {
        contract::fail(state, "prime-ideal reduction failed",
                       contract::FailureReason::operation);
        return;
    }

    const slong residue = fmpz_poly_get_coeff_si(reduced.raw(), 0);
    state.counters["residue"] = static_cast<double>(residue);
    if (fmpz_poly_degree(reduced.raw()) != 0 ||
        (residue != 4 && residue != 8)) {
        contract::fail(state, "prime-ideal reduction invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

void BM_prime_ideal_valuation_quadratic(benchmark::State& state) {
    contract::initialize(state);

    silex::NumberField field;
    silex::Order maximal_order(field);
    if (!define_maximal_quadratic(field, maximal_order, 2)) {
        contract::fail(state, "maximal order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz p;
    sflint::fmpz_set_si(p, 3);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, maximal_order,
                                sflint::FmpzConstRef(p)) ||
        primes.size() <= 0 || primes.at(0) == nullptr) {
        contract::fail(state, "prime decomposition failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::OrderElement generator(maximal_order);
    silex::Ideal ideal(maximal_order);
    if (!generator.set_si(9) || !ideal.set_principal(generator)) {
        contract::fail(state, "principal ideal setup failed",
                       contract::FailureReason::setup);
        return;
    }

    const silex::PrimeIdeal& prime = *primes.at(0);
    slong value = -1;
    bool all_operations_ok = true;
    for (auto _ : state) {
        value = -1;
        const bool iteration_ok = prime.valuation(value, ideal);
        all_operations_ok = iteration_ok && all_operations_ok;
        benchmark::DoNotOptimize(value);
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed prime-ideal valuation failed",
                       contract::FailureReason::operation);
        return;
    }
    value = -1;
    if (!prime.valuation(value, ideal)) {
        contract::fail(state, "prime-ideal valuation failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["valuation"] = static_cast<double>(value);
    if (value != 2) {
        contract::fail(state, "prime-ideal valuation invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

void BM_prime_ideal_valuation_order_element_quadratic_split(
        benchmark::State& state) {
    contract::initialize(state);

    silex::NumberField field;
    silex::Order order(field);
    if (!define_marked_equation_quadratic(field, order, 2)) {
        contract::fail(state, "equation order setup failed",
                       contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz p;
    sflint::fmpz_set_si(p, 7);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)) ||
        primes.size() <= 0) {
        contract::fail(state, "prime decomposition failed",
                       contract::FailureReason::setup);
        return;
    }

    silex::Element theta(field);
    silex::Element three(field);
    silex::Element theta_minus_three(field);
    silex::Element square(field);
    silex::OrderElement generator(order);
    if (!theta.gen() || !three.set_si(3) ||
        !theta_minus_three.subtract(theta, three) ||
        !square.multiply(theta_minus_three, theta_minus_three) ||
        !generator.set_element(square)) {
        contract::fail(state, "order element setup failed",
                       contract::FailureReason::setup);
        return;
    }

    const silex::PrimeIdeal* prime = nullptr;
    sflint::FmpzPoly reduced;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* candidate = primes.at(i);
        if (candidate != nullptr &&
            candidate->reduce(sflint::FmpzPolyRef(reduced), theta) &&
            fmpz_poly_degree(reduced.raw()) == 0 &&
            fmpz_poly_get_coeff_si(reduced.raw(), 0) == 3) {
            prime = candidate;
            break;
        }
    }
    if (prime == nullptr) {
        contract::fail(state, "target split prime not found",
                       contract::FailureReason::setup);
        return;
    }

    slong value = -1;
    bool all_operations_ok = true;
    for (auto _ : state) {
        value = -1;
        const bool iteration_ok = prime->valuation(value, generator);
        all_operations_ok = iteration_ok && all_operations_ok;
        benchmark::DoNotOptimize(value);
    }

    benchmark::ClobberMemory();
    if (!all_operations_ok) {
        contract::fail(state, "timed order-element valuation failed",
                       contract::FailureReason::operation);
        return;
    }
    value = -1;
    if (!prime->valuation(value, generator)) {
        contract::fail(state, "order-element valuation failed",
                       contract::FailureReason::operation);
        return;
    }
    state.counters["valuation"] = static_cast<double>(value);
    if (value != 2) {
        contract::fail(state, "order-element valuation invariant failed",
                       contract::FailureReason::invariant);
        return;
    }
    contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_prime_ideal_reduce_nf_elem_quadratic);
BENCHMARK(BM_prime_ideal_valuation_quadratic);
BENCHMARK(BM_prime_ideal_valuation_order_element_quadratic_split);
BENCHMARK_MAIN();
