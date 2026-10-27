#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/ideal_factorization.hpp>

#include "ideal/ideal_access_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"

#include <flint/fmpz_factor.h>

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

bool set_rational_principal(silex::Ideal& ideal, slong value) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr) {
        return false;
    }

    silex::OrderElement generator(*order);
    return generator.set_si(value) && ideal.set_principal(generator);
}

bool set_fractional_principal(silex::FractionalIdeal& ideal,
                              slong numerator,
                              slong denominator) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr || order->parent() == nullptr) {
        return false;
    }

    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator,
                        static_cast<ulong>(denominator));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);

    silex::Element element(*order->parent());
    return element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)) &&
           ideal.set_principal(element);
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

slong count_prime(const silex::IdealFactorization& factorization,
                  slong rational_prime,
                  slong ramification_index,
                  slong residue_degree,
                  slong exponent) noexcept {
    const silex::Order* parent = factorization.parent();
    if (parent == nullptr) {
        return 0;
    }

    silex::PrimeIdeal prime(*parent);
    sflint::Fmpz p;
    slong count = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exp = -1;
        if (!factorization.prime(prime, i) ||
            !factorization.exponent(exp, i) ||
            !prime.rational_prime(sflint::FmpzRef(p))) {
            return -1;
        }
        if (sflint::fmpz_equal_si(p, rational_prime) &&
            prime.ramification_index() == ramification_index &&
            prime.residue_degree() == residue_degree && exp == exponent) {
            ++count;
        }
    }
    return count;
}

slong exponent_of_prime(const silex::IdealFactorization& factorization,
                        const silex::PrimeIdeal& target) noexcept {
    const silex::Order* parent = factorization.parent();
    if (!silex::same_order_parent(parent, target.parent())) {
        return -1;
    }

    silex::PrimeIdeal prime(*parent);
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = -1;
        if (!factorization.prime(prime, i) ||
            !factorization.exponent(exponent, i)) {
            return -1;
        }
        if (prime.equal(target)) {
            return exponent;
        }
    }
    return 0;
}

bool row_matches_factorization(const sflint::FmpzMat& row,
                               const silex::FactorBase& base,
                               const silex::IdealFactorization& factorization)
        noexcept {
    const silex::Order* parent = base.parent();
    if (parent == nullptr ||
        !silex::same_order_parent(factorization.parent(), parent) ||
        sflint::fmpz_mat_nrows(row) != 1 ||
        sflint::fmpz_mat_ncols(row) != base.length()) {
        return false;
    }

    silex::PrimeIdeal prime(*parent);
    for (slong i = 0; i < base.length(); ++i) {
        if (!base.prime(prime, i)) {
            return false;
        }
        const slong exponent = exponent_of_prime(factorization, prime);
        if (exponent < 0 ||
            !sflint::fmpz_equal_si(sflint::fmpz_mat_entry(row, 0, i),
                                   exponent)) {
            return false;
        }
    }
    return true;
}

// Reference for the one non-base factor a required-prime witness may carry:
// a positive power of the only prime ideal above its rational prime q that
// the base lacks, with e = 1, when the base holds another ideal above q and
// q is not the required prime's rational prime.  The relation
// qO_K = Q_1^e_1 ... Q_g^e_g then puts that class in the subgroup of the base
// ideals above q.
bool lacking_prime_is_admissible(const silex::PrimeIdeal& factor,
                                 slong exponent,
                                 const silex::FactorBase& base,
                                 const silex::PrimeIdeal& required_prime)
        noexcept {
    sflint::Fmpz q;
    sflint::Fmpz required_p;
    if (exponent <= 0 || factor.ramification_index() != 1 ||
        !factor.rational_prime(sflint::FmpzRef(q)) ||
        !required_prime.rational_prime(sflint::FmpzRef(required_p)) ||
        sflint::fmpz_equal(sflint::FmpzConstRef(q),
                           sflint::FmpzConstRef(required_p))) {
        return false;
    }
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, *factor.parent(),
                                sflint::FmpzConstRef(q))) {
        return false;
    }
    slong in_base = 0;
    slong outside = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) {
            return false;
        }
        if (base.contains(*prime)) {
            ++in_base;
        } else {
            ++outside;
        }
    }
    return in_base > 0 && outside == 1;
}

bool full_factorization_has_required_prime(
        const silex::Ideal& ideal,
        const silex::FactorBase& base,
        const silex::PrimeIdeal& required_prime) noexcept {
    const silex::Order* parent = ideal.parent();
    if (parent == nullptr ||
        !silex::same_order_parent(base.parent(), parent) ||
        !silex::same_order_parent(required_prime.parent(), parent)) {
        return false;
    }

    silex::IdealFactorization factorization(*parent);
    silex::PrimeIdeal factor(*parent);
    if (!factorization.factor(ideal)) {
        return false;
    }

    slong required_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        if (factor.equal(required_prime)) {
            required_exponent += exponent;
        } else if (!base.contains(factor) &&
                   !lacking_prime_is_admissible(factor, exponent, base,
                                                required_prime)) {
            return false;
        }
    }
    return required_exponent == 1;
}

bool first_prime_above(silex::PrimeIdeal& out,
                       const silex::Order& order,
                       slong rational_prime) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList primes;
    const silex::PrimeIdeal* first = nullptr;
    return silex::decompose_prime(
                   primes, order, sflint::FmpzConstRef(p)) &&
           primes.size() > 0 && (first = primes.at(0)) != nullptr &&
           out.set(*first);
}

bool required_prime_result_matches(
        const silex::Ideal& ideal,
        const silex::FactorBase& base,
        const silex::PrimeIdeal& required_prime,
        bool expected) noexcept {
    bool matches = !expected;
    return silex::detail::ideal_factor_over_base_with_required_prime(
                   matches, ideal, base, required_prime) &&
           matches == expected &&
           matches == full_factorization_has_required_prime(
                              ideal, base, required_prime);
}

int test_degree_one_factorization() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::Ideal ideal(order);
    silex::Ideal reconstructed(order);
    silex::IdealFactorization factorization;
    silex::IdealFactorization copy(order);
    silex::IdealFactorization other(order);

    assert(set_rational_principal(ideal, 60));
    assert(factorization.factor(ideal));
    assert(factorization.length() == 3);
    assert(count_prime(factorization, 2, 1, 1, 2) == 1);
    assert(count_prime(factorization, 3, 1, 1, 1) == 1);
    assert(count_prime(factorization, 5, 1, 1, 1) == 1);
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.equal(ideal));

    assert(copy.set(factorization));
    assert(copy.length() == 3);
    assert(count_prime(copy, 2, 1, 1, 2) == 1);
    assert(count_prime(copy, 3, 1, 1, 1) == 1);
    assert(count_prime(copy, 5, 1, 1, 1) == 1);
    swap(copy, other);
    assert(copy.length() == 0);
    assert(other.length() == 3);
    assert(count_prime(other, 2, 1, 1, 2) == 1);

    assert(set_rational_principal(ideal, 1));
    assert(factorization.factor(ideal));
    assert(factorization.length() == 0);
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.equal(ideal));

    return 0;
}

int test_quadratic_factorization() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::Ideal ideal(order);
    silex::Ideal reconstructed(order);
    silex::IdealFactorization factorization(order);

    assert(set_rational_principal(ideal, 9));
    assert(factorization.factor(ideal));
    assert(factorization.length() == 1);
    assert(count_prime(factorization, 3, 1, 2, 2) == 1);
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.equal(ideal));

    assert(set_rational_principal(ideal, 7));
    assert(factorization.factor(ideal));
    assert(factorization.length() == 2);
    assert(count_prime(factorization, 7, 1, 1, 1) == 2);
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.equal(ideal));

    assert(set_rational_principal(ideal, 2));
    assert(factorization.factor(ideal));
    assert(factorization.length() == 1);
    assert(count_prime(factorization, 2, 2, 1, 2) == 1);
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.equal(ideal));

    return 0;
}

int test_factor_failure_preserves_output() {
    silex::NumberField rational_field = degree_one_field();
    silex::Order rational_order = silex::test::equation_order(rational_field);

    silex::NumberField quadratic = quadratic_field(5);
    silex::Order nonmaximal = silex::test::equation_order(quadratic);
    assert(!nonmaximal.is_maximal());

    silex::Ideal good(rational_order);
    silex::Ideal bad(nonmaximal);
    silex::IdealFactorization factorization;
    assert(set_rational_principal(good, 5));
    assert(factorization.factor(good));
    assert(factorization.length() == 1);
    assert(count_prime(factorization, 5, 1, 1, 1) == 1);

    assert(set_rational_principal(bad, 9));
    assert(!factorization.factor(bad));
    assert(factorization.length() == 1);
    assert(count_prime(factorization, 5, 1, 1, 1) == 1);

    return 0;
}

int test_factor_over_base_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 3);

    silex::Ideal ideal(order);
    sflint::FmpzMat exponents(1, 3);
    sflint::FmpzMat wrong(2, 3);
    bool smooth = false;

    assert(set_rational_principal(ideal, 30));
    smooth = false;
    assert(silex::ideal_is_smooth(smooth, ideal, base));
    assert(smooth);
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                         ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));
    assert(mat_entry_is_si(exponents, 0, 1, 1));
    assert(mat_entry_is_si(exponents, 0, 2, 1));

    assert(set_rational_principal(ideal, 210));
    smooth = true;
    assert(silex::ideal_is_smooth(smooth, ideal, base));
    assert(!smooth);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                          ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));
    assert(mat_entry_is_si(exponents, 0, 1, 1));
    assert(mat_entry_is_si(exponents, 0, 2, 1));

    assert(set_rational_principal(ideal, 2000006));
    smooth = true;
    assert(silex::ideal_is_smooth(smooth, ideal, base));
    assert(!smooth);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                          ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));
    assert(mat_entry_is_si(exponents, 0, 1, 1));
    assert(mat_entry_is_si(exponents, 0, 2, 1));

    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(wrong),
                                          ideal, base));

    return 0;
}

int test_factor_over_base_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 5);

    silex::Ideal ideal(order);
    sflint::FmpzMat exponents(1, 5);
    bool smooth = false;

    assert(set_rational_principal(ideal, 14));
    assert(silex::ideal_is_smooth(smooth, ideal, base));
    assert(smooth);
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                         ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 2));
    assert(mat_entry_is_si(exponents, 0, 1, 0));
    assert(mat_entry_is_si(exponents, 0, 2, 0));
    assert(mat_entry_is_si(exponents, 0, 3, 1));
    assert(mat_entry_is_si(exponents, 0, 4, 1));

    silex::IdealFactorization full_factorization(order);
    assert(full_factorization.factor(ideal));
    assert(row_matches_factorization(exponents, base, full_factorization));

    assert(set_rational_principal(ideal, 22));
    smooth = true;
    assert(silex::ideal_is_smooth(smooth, ideal, base));
    assert(!smooth);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                          ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 2));
    assert(mat_entry_is_si(exponents, 0, 3, 1));
    assert(mat_entry_is_si(exponents, 0, 4, 1));

    return 0;
}

int test_factor_over_base_with_required_prime() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::PrimeIdeal required_prime(order);
    silex::PrimeIdeal other_prime(order);
    silex::PrimeIdeal base_prime(order);
    assert(first_prime_above(required_prime, order, 11));
    assert(first_prime_above(other_prime, order, 13));
    assert(base.prime(base_prime, 0));
    assert(!base.contains(required_prime));
    assert(!base.contains(other_prime));

    silex::Ideal required_ideal(order);
    silex::Ideal other_ideal(order);
    silex::Ideal base_ideal(order);
    silex::Ideal product(order);
    assert(required_prime.get_ideal(required_ideal));
    assert(other_prime.get_ideal(other_ideal));
    assert(base_prime.get_ideal(base_ideal));

    assert(required_prime_result_matches(
            required_ideal, base, required_prime, true));

    assert(product.multiply(required_ideal, base_ideal));
    assert(required_prime_result_matches(
            product, base, required_prime, true));

    assert(product.multiply(required_ideal, required_ideal));
    assert(required_prime_result_matches(
            product, base, required_prime, false));

    assert(product.multiply(required_ideal, other_ideal));
    assert(required_prime_result_matches(
            product, base, required_prime, false));

    return 0;
}

int test_order_element_factor_over_base_with_one_large_prime() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());
    assert(!order.is_equation_order());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::OrderElement element(order);
    silex::Element ambient(field);
    sflint::Fmpq norm;
    sflint::FmpzMat row(1, base.length());
    silex::PrimeIdeal large_prime(order);
    silex::detail::OneLargePrimeFactorStatus status =
            silex::detail::OneLargePrimeFactorStatus::no_candidate;

    assert(element.set_si(6));
    assert(element.get_element(ambient));
    assert(ambient.norm(sflint::FmpqRef(norm)));
    assert(silex::detail::order_element_factor_over_base_with_one_large_prime(
            status, sflint::FmpzMatRef(row), large_prime, element,
            sflint::FmpqConstRef(norm), base));
    assert(status == silex::detail::OneLargePrimeFactorStatus::found);
    assert(mat_entry_is_si(row, 0, 0, 1));
    assert(!base.contains(large_prime));

    silex::Ideal principal(order);
    silex::IdealFactorization full(order);
    assert(principal.set_principal(element));
    assert(full.factor(principal));
    assert(row_matches_factorization(row, base, full));
    assert(exponent_of_prime(full, large_prime) == 1);

    // 3 is inert, so the q-part left at q = 3 is 9, not 3: the large prime
    // has degree 2 and is built even when deferral is offered.
    silex::detail::OneLargePrimeDeferral deferral;
    silex::PrimeIdeal offered(order);
    sflint::FmpzMat offered_row(1, base.length());
    assert(silex::detail::order_element_factor_over_base_with_one_large_prime(
            status, sflint::FmpzMatRef(offered_row), offered, element,
            sflint::FmpqConstRef(norm), base, nullptr, &deferral));
    assert(status == silex::detail::OneLargePrimeFactorStatus::found);
    assert(offered.equal(large_prime));
    assert(offered.residue_degree() == 2);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(offered_row),
                                  sflint::FmpzMatConstRef(row)));

    const slong no_candidate_values[] = {2, 9, 21};
    for (slong value : no_candidate_values) {
        assert(element.set_si(value));
        assert(element.get_element(ambient));
        assert(ambient.norm(sflint::FmpqRef(norm)));
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(row), large_prime,
                               element, sflint::FmpqConstRef(norm), base));
        assert(status ==
               silex::detail::OneLargePrimeFactorStatus::no_candidate);
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(row), large_prime,
                               element, sflint::FmpqConstRef(norm), base,
                               nullptr, &deferral));
        assert(status ==
               silex::detail::OneLargePrimeFactorStatus::no_candidate);
    }

    return 0;
}

int test_factor_over_base_parent_failure_preserves_outputs() {
    silex::NumberField rational_field = degree_one_field();
    silex::Order rational_order = silex::test::equation_order(rational_field);

    silex::NumberField quadratic = quadratic_field(2);
    silex::Order quadratic_order = silex::test::equation_order(quadratic);

    silex::FactorBase base(rational_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Ideal rational_ideal(rational_order);
    silex::Ideal quadratic_ideal(quadratic_order);
    sflint::FmpzMat exponents(1, 1);
    assert(set_rational_principal(rational_ideal, 2));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                         rational_ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));

    bool smooth = true;
    assert(set_rational_principal(quadratic_ideal, 2));
    assert(!silex::ideal_is_smooth(smooth, quadratic_ideal, base));
    assert(smooth);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                          quadratic_ideal, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));

    silex::FractionalIdeal quadratic_fractional(quadratic_order);
    assert(set_fractional_principal(quadratic_fractional, 1, 2));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(exponents),
                                          quadratic_fractional, base));
    assert(mat_entry_is_si(exponents, 0, 0, 1));

    return 0;
}

int test_fractional_factor_over_base_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::FractionalIdeal ideal(order);
    sflint::FmpzMat row(1, 1);

    assert(set_fractional_principal(ideal, 1, 2));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                         ideal, base));
    assert(mat_entry_is_si(row, 0, 0, -1));

    assert(set_fractional_principal(ideal, 3, 2));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                          ideal, base));
    assert(mat_entry_is_si(row, 0, 0, -1));

    assert(set_fractional_principal(ideal, 4, 1));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                         ideal, base));
    assert(mat_entry_is_si(row, 0, 0, 2));

    return 0;
}

int test_fractional_factor_over_base_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::FractionalIdeal ideal(order);
    sflint::FmpzMat row(1, 1);

    assert(set_fractional_principal(ideal, 1, 2));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                         ideal, base));
    assert(mat_entry_is_si(row, 0, 0, -2));

    assert(set_fractional_principal(ideal, 3, 2));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                          ideal, base));
    assert(mat_entry_is_si(row, 0, 0, -2));

    return 0;
}

silex::NumberField field_from_coefficients(const slong* coefficients,
                                           slong length) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < length; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool set_prime_power(silex::FractionalIdeal& out,
                     const silex::PrimeIdeal& prime,
                     slong exponent) noexcept {
    const silex::Order* order = prime.parent();
    if (order == nullptr) {
        return false;
    }
    silex::Ideal integral(*order);
    silex::FractionalIdeal fractional(*order);
    sflint::Fmpz power;
    sflint::fmpz_set_si(sflint::FmpzRef(power), exponent);
    return prime.get_ideal(integral) && fractional.set_integral(integral) &&
           out.pow_fmpz(fractional, sflint::FmpzConstRef(power));
}

bool multiply_prime_power(silex::FractionalIdeal& inout,
                          const silex::PrimeIdeal& prime,
                          slong exponent) noexcept {
    const silex::Order* order = prime.parent();
    if (order == nullptr) {
        return false;
    }
    silex::FractionalIdeal factor(*order);
    return set_prime_power(factor, prime, exponent) &&
           inout.multiply(inout, factor);
}

// Primes above `rational_prime` with residue degree `residue_degree`, in
// decomposition order.
slong primes_above_with_degree(silex::PrimeIdeal* out,
                               slong capacity,
                               const silex::Order& order,
                               slong rational_prime,
                               slong residue_degree) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order, sflint::FmpzConstRef(p))) {
        return -1;
    }
    slong count = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) {
            return -1;
        }
        if (prime->residue_degree() != residue_degree) {
            continue;
        }
        if (count < capacity && !out[count].set(*prime)) {
            return -1;
        }
        ++count;
    }
    return count;
}

// Every base entry other than the listed (prime, exponent) pairs is zero.
bool row_is(const sflint::FmpzMat& row,
            const silex::FactorBase& base,
            const silex::PrimeIdeal* primes,
            const slong* exponents,
            slong count) noexcept {
    if (sflint::fmpz_mat_nrows(row) != 1 ||
        sflint::fmpz_mat_ncols(row) != base.length()) {
        return false;
    }
    for (slong i = 0; i < base.length(); ++i) {
        const silex::PrimeIdeal* base_prime = base.prime_at(i);
        if (base_prime == nullptr) {
            return false;
        }
        slong expected = 0;
        for (slong j = 0; j < count; ++j) {
            if (base_prime->equal(primes[j])) {
                expected = exponents[j];
            }
        }
        if (!sflint::fmpz_equal_si(sflint::fmpz_mat_entry(row, 0, i),
                                   expected)) {
            return false;
        }
    }
    return true;
}

void fill_row(sflint::FmpzMat& row, slong value) noexcept {
    for (slong i = 0; i < sflint::fmpz_mat_ncols(row); ++i) {
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, i), value);
    }
}

bool row_is_filled(const sflint::FmpzMat& row, slong value) noexcept {
    for (slong i = 0; i < sflint::fmpz_mat_ncols(row); ++i) {
        if (!sflint::fmpz_equal_si(sflint::fmpz_mat_entry(row, 0, i),
                                   value)) {
            return false;
        }
    }
    return true;
}

// The cubic x^3 - 2 (maximal order Z[theta]), base bound 5: a - theta has
// norm a^3 - 2.  At a = 5 (3 * 41) and a = 9 (727, which splits completely)
// the one large prime is the degree-one prime (p, theta - a), found among
// the primes above p that can contain the element.  At a = 7 (11 * 31) and
// for 31 (three primes of norm 31) there are two large primes.  With
// deferral offered, the found cases leave a q-part of exactly q, so they are
// deferred and the deferred build gives the same prime; the others, and
// (5 - theta)^2 (one degree-one prime with valuation 2), are rejected.
int test_order_element_factor_one_large_prime_cubic() {
    const slong coefficients[] = {-2, 0, 0, 1};
    silex::NumberField field = field_from_coefficients(coefficients, 4);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::OrderElement element(order);
    silex::Element ambient(field);
    sflint::FmpqPoly polynomial;
    sflint::Fmpq norm;
    sflint::FmpzMat row(1, base.length());
    silex::PrimeIdeal large_prime(order);
    silex::detail::OneLargePrimeFactorStatus status =
            silex::detail::OneLargePrimeFactorStatus::no_candidate;
    auto set_a_minus_theta = [&](slong a) {
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, a);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
        assert(ambient.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        assert(element.set_element(ambient));
        assert(ambient.norm(sflint::FmpqRef(norm)));
    };

    const slong found_values[] = {5, 9};
    for (slong a : found_values) {
        set_a_minus_theta(a);
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(row), large_prime,
                               element, sflint::FmpqConstRef(norm), base));
        assert(status == silex::detail::OneLargePrimeFactorStatus::found);
        assert(!base.contains(large_prime));
        assert(large_prime.residue_degree() == 1);

        silex::Ideal principal(order);
        silex::IdealFactorization full(order);
        assert(principal.set_principal(element));
        assert(full.factor(principal));
        assert(row_matches_factorization(row, base, full));
        assert(exponent_of_prime(full, large_prime) == 1);

        silex::detail::OneLargePrimeDeferral deferral;
        silex::PrimeIdeal deferred(order);
        sflint::FmpzMat deferred_row(1, base.length());
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(deferred_row),
                               deferred, element, sflint::FmpqConstRef(norm),
                               base, nullptr, &deferral));
        assert(status ==
               silex::detail::OneLargePrimeFactorStatus::found_deferred);
        assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(deferred_row),
                                      sflint::FmpzMatConstRef(row)));
        sflint::Fmpz large_norm;
        assert(large_prime.norm(sflint::FmpzRef(large_norm)));
        assert(fmpz_equal(deferral.rational_prime.raw(), large_norm.raw()));
        assert(deferral.norm_valuation == 1);
        assert(silex::detail::order_element_deferred_large_prime(
                deferred, element,
                sflint::FmpzConstRef(deferral.rational_prime),
                deferral.norm_valuation, base));
        assert(deferred.equal(large_prime));
    }

    auto rejected_both_ways = [&]() {
        silex::detail::OneLargePrimeDeferral deferral;
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(row), large_prime,
                               element, sflint::FmpqConstRef(norm), base));
        assert(status ==
               silex::detail::OneLargePrimeFactorStatus::no_candidate);
        assert(silex::detail::
                       order_element_factor_over_base_with_one_large_prime(
                               status, sflint::FmpzMatRef(row), large_prime,
                               element, sflint::FmpqConstRef(norm), base,
                               nullptr, &deferral));
        assert(status ==
               silex::detail::OneLargePrimeFactorStatus::no_candidate);
    };

    set_a_minus_theta(7);
    rejected_both_ways();

    assert(element.set_si(31));
    assert(element.get_element(ambient));
    assert(ambient.norm(sflint::FmpqRef(norm)));
    rejected_both_ways();

    set_a_minus_theta(5);
    silex::OrderElement square(order);
    assert(square.multiply(element, element));
    assert(element.set(square));
    assert(element.get_element(ambient));
    assert(ambient.norm(sflint::FmpqRef(norm)));
    rejected_both_ways();
    return 0;
}

// Q(2^(1/3)): 5 = P Q with f(P) = 1 and f(Q) = 2.  The norm-bounded base at
// 5 keeps P and drops Q, so its block above 5 is incomplete.
int test_required_prime_through_lacking_conjugate() {
    // x^3 - 2: 5 = R1 R2 with f(R1) = 1, f(R2) = 2, and 31 = P1 P2 P3.
    const slong coefficients[] = {-2, 0, 0, 1};
    silex::NumberField field = field_from_coefficients(coefficients, 4);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::PrimeIdeal r1(order);
    silex::PrimeIdeal r2(order);
    silex::PrimeIdeal p31[] = {silex::PrimeIdeal(order),
                               silex::PrimeIdeal(order),
                               silex::PrimeIdeal(order)};
    assert(primes_above_with_degree(&r1, 1, order, 5, 1) == 1);
    assert(primes_above_with_degree(&r2, 1, order, 5, 2) == 1);
    assert(primes_above_with_degree(p31, 3, order, 31, 1) == 3);

    // The norm-bounded base holds R1 and lacks R2.
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    silex::FactorBase partial(order);
    assert(partial.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(partial.contains(r1) && !partial.contains(r2));
    assert(!partial.contains(p31[0]) && !partial.contains(p31[1]));

    silex::Ideal required_ideal(order);
    silex::Ideal r2_ideal(order);
    silex::Ideal other_ideal(order);
    silex::Ideal product(order);
    assert(p31[0].get_ideal(required_ideal));
    assert(r2.get_ideal(r2_ideal));
    assert(p31[1].get_ideal(other_ideal));

    // P1 R2: R2 is the only ideal above 5 outside the base, with e = 1.
    assert(product.multiply(required_ideal, r2_ideal));
    assert(required_prime_result_matches(product, partial, p31[0], true));
    // P1 P2: P2 lies above the required prime's own rational prime.
    assert(product.multiply(required_ideal, other_ideal));
    assert(required_prime_result_matches(product, partial, p31[0], false));
    return 0;
}

int test_fractional_factor_over_incomplete_block() {
    const slong coefficients[] = {-2, 0, 0, 1};
    silex::NumberField field = field_from_coefficients(coefficients, 4);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    // The rows below compare against both primes, so they live in one array
    // from the start and the names refer into it.
    silex::PrimeIdeal primes[] = {silex::PrimeIdeal(order),
                                  silex::PrimeIdeal(order)};
    silex::PrimeIdeal& p_prime = primes[0];
    silex::PrimeIdeal& q_prime = primes[1];
    assert(primes_above_with_degree(&p_prime, 1, order, 5, 1) == 1);
    assert(primes_above_with_degree(&q_prime, 1, order, 5, 2) == 1);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    silex::FactorBase partial(order);
    silex::FactorBase full(order);
    assert(partial.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(full.build(sflint::FmpzConstRef(bound)));
    assert(partial.contains(p_prime) && !partial.contains(q_prime));
    assert(full.contains(p_prime) && full.contains(q_prime));

    silex::FractionalIdeal ideal(order);
    sflint::FmpzMat row(1, partial.length());
    sflint::FmpzMat full_row(1, full.length());

    // P^{-1} = Q / 5: the numerator Q is outside the base, but P^{-1} is
    // smooth over {P}.
    assert(set_prime_power(ideal, p_prime, -1));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                         partial));
    const slong minus_one[] = {-1};
    assert(row_is(row, partial, &p_prime, minus_one, 1));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(full_row), ideal,
                                         full));
    assert(row_is(full_row, full, &p_prime, minus_one, 1));

    // P^2 Q^{-1} is smooth over the complete base and not over {P}.
    assert(set_prime_power(ideal, p_prime, 2));
    assert(multiply_prime_power(ideal, q_prime, -1));
    const slong mixed_exponents[] = {2, -1};
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(full_row), ideal,
                                         full));
    assert(row_is(full_row, full, primes, mixed_exponents, 2));
    fill_row(row, 7);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                          partial));
    assert(row_is_filled(row, 7));

    // Q^{-1} and P / 7 are not smooth over {P}; the output is unchanged.
    assert(set_prime_power(ideal, q_prime, -1));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                          partial));
    assert(row_is_filled(row, 7));
    assert(set_fractional_principal(ideal, 1, 7));
    assert(multiply_prime_power(ideal, p_prime, 1));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                          partial));
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(full_row), ideal,
                                          full));
    assert(row_is_filled(row, 7));

    // 2 P^{-1} / 3 is smooth over both bases.
    assert(set_fractional_principal(ideal, 2, 3));
    assert(multiply_prime_power(ideal, p_prime, -1));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                         partial));
    silex::Ideal integral(order);
    silex::IdealFactorization two(order);
    silex::IdealFactorization three(order);
    assert(set_rational_principal(integral, 2) && two.factor(integral));
    assert(set_rational_principal(integral, 3) && three.factor(integral));
    for (slong i = 0; i < partial.length(); ++i) {
        const silex::PrimeIdeal* prime = partial.prime_at(i);
        assert(prime != nullptr);
        slong expected = exponent_of_prime(two, *prime) -
                         exponent_of_prime(three, *prime);
        if (prime->equal(p_prime)) {
            expected -= 1;
        }
        assert(sflint::fmpz_equal_si(sflint::fmpz_mat_entry(row, 0, i),
                                     expected));
    }

    return 0;
}

// x^5 - 2x^4 - 2x^3 - x^2 - 2 has Z[x] maximal and 3 = P Q R with
// f(P) = 1 and f(Q) = f(R) = 2.  Over the norm-bounded base at 3 the block
// above 3 is {P}.  Q R^{-1} has norm 1 and v_P = 0, so norm accounting
// alone would accept it; the non-base valuations must reject it.
int test_fractional_factor_over_base_cancelling_nonbase_primes() {
    const slong coefficients[] = {-2, 0, -1, -2, -2, 1};
    silex::NumberField field = field_from_coefficients(coefficients, 6);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::PrimeIdeal p_prime(order);
    silex::PrimeIdeal pair[2] = {silex::PrimeIdeal(order),
                                 silex::PrimeIdeal(order)};
    assert(primes_above_with_degree(&p_prime, 1, order, 3, 1) == 1);
    assert(primes_above_with_degree(pair, 2, order, 3, 2) == 2);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    silex::FactorBase partial(order);
    silex::FactorBase full(order);
    assert(partial.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(full.build(sflint::FmpzConstRef(bound)));
    assert(partial.contains(p_prime));
    assert(!partial.contains(pair[0]) && !partial.contains(pair[1]));
    assert(full.contains(pair[0]) && full.contains(pair[1]));

    silex::FractionalIdeal ideal(order);
    sflint::FmpzMat row(1, partial.length());
    sflint::FmpzMat full_row(1, full.length());
    assert(set_prime_power(ideal, pair[0], 1));
    assert(multiply_prime_power(ideal, pair[1], -1));

    fill_row(row, 7);
    assert(!silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                          partial));
    assert(row_is_filled(row, 7));

    const slong exponents[] = {1, -1};
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(full_row), ideal,
                                         full));
    assert(row_is(full_row, full, pair, exponents, 2));

    // P^{-1} Q^0 R^0 remains smooth over {P}.
    assert(set_prime_power(ideal, p_prime, -1));
    const slong minus_one[] = {-1};
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row), ideal,
                                         partial));
    assert(row_is(row, partial, &p_prime, minus_one, 1));

    return 0;
}

silex::IdealFactorization local_ideal_factorization() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 4));

    silex::IdealFactorization factorization(order);
    assert(factorization.factor(ideal));
    return factorization;
}

// The norm accounting that IdealFactorization::factor checks in place of a
// reconstruction (buch2.c:can_factor/divide_p_id): every factor has positive
// exponent and lies over some p | N(I), and sum_{P | p} f_P v_P = v_p(N(I))
// for every p | N(I).
bool norm_accounting_holds(const silex::IdealFactorization& factorization,
                           const silex::Ideal& ideal) noexcept {
    const silex::Order* parent = factorization.parent();
    sflint::Fmpz norm;
    sflint::Fmpz p;
    sflint::Fmpz factor_p;
    if (parent == nullptr || !ideal.norm(sflint::FmpzRef(norm)) ||
        fmpz_sgn(norm.raw()) <= 0) {
        return false;
    }
    sflint::FmpzFactor norm_factorization;
    fmpz_factor(norm_factorization.raw(), norm.raw());
    silex::PrimeIdeal prime(*parent);
    slong matched = 0;
    for (slong i = 0; i < norm_factorization.raw()->num; ++i) {
        fmpz_factor_get_fmpz(p.raw(), norm_factorization.raw(), i);
        slong accounted = 0;
        for (slong j = 0; j < factorization.length(); ++j) {
            slong exponent = 0;
            if (!factorization.prime(prime, j) ||
                !factorization.exponent(exponent, j) || exponent <= 0 ||
                !prime.rational_prime(sflint::FmpzRef(factor_p))) {
                return false;
            }
            if (fmpz_equal(factor_p.raw(), p.raw()) != 0) {
                accounted += prime.residue_degree() * exponent;
                ++matched;
            }
        }
        if (accounted != static_cast<slong>(norm_factorization.raw()->exp[i])) {
            return false;
        }
    }
    return matched == factorization.length();
}

int test_factor_norm_accounting() {
    // Q(sqrt 2): 2 ramifies, 3 is inert (f = 2), 7 splits.
    {
        silex::NumberField field = quadratic_field(2);
        silex::Order order = silex::test::equation_order(field);
        assert(order.is_maximal());
        silex::Ideal ideal(order);
        silex::Ideal reconstructed(order);
        silex::IdealFactorization factorization(order);
        assert(set_rational_principal(ideal, 42));
        assert(factorization.factor(ideal));
        assert(factorization.length() == 4);
        assert(count_prime(factorization, 2, 2, 1, 2) == 1);
        assert(count_prime(factorization, 3, 1, 2, 1) == 1);
        assert(count_prime(factorization, 7, 1, 1, 1) == 2);
        assert(norm_accounting_holds(factorization, ideal));
        assert(factorization.reconstruct(reconstructed));
        assert(reconstructed.equal(ideal));
    }

    // Q(cbrt 2): 2 and 3 are totally ramified, 7 is inert (f = 3), and
    // 5 = P P' and 41 = Q Q' with f(P) = f(Q) = 1, f(P') = f(Q') = 2.
    const slong coefficients[] = {-2, 0, 0, 1};
    silex::NumberField field = field_from_coefficients(coefficients, 4);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    // I = P^3 P' (14), built by ideal multiplication (no stored generator).
    {
        silex::PrimeIdeal degree_one(order);
        silex::PrimeIdeal degree_two(order);
        assert(primes_above_with_degree(&degree_one, 1, order, 5, 1) == 1);
        assert(primes_above_with_degree(&degree_two, 1, order, 5, 2) == 1);
        silex::Ideal p(order);
        silex::Ideal p_prime(order);
        silex::Ideal fourteen(order);
        silex::Ideal left(order);
        silex::Ideal right(order);
        silex::Ideal ideal(order);
        assert(degree_one.get_ideal(p));
        assert(degree_two.get_ideal(p_prime));
        assert(set_rational_principal(fourteen, 14));
        assert(left.multiply(p, p));
        assert(right.multiply(left, p));
        assert(left.multiply(right, p_prime));
        assert(ideal.multiply(left, fourteen));

        silex::IdealFactorization factorization(order);
        silex::Ideal reconstructed(order);
        assert(factorization.factor(ideal));
        assert(factorization.length() == 4);
        assert(count_prime(factorization, 5, 1, 1, 3) == 1);
        assert(count_prime(factorization, 5, 1, 2, 1) == 1);
        assert(count_prime(factorization, 2, 3, 1, 3) == 1);
        assert(count_prime(factorization, 7, 1, 3, 1) == 1);
        assert(norm_accounting_holds(factorization, ideal));
        assert(factorization.reconstruct(reconstructed));
        assert(reconstructed.equal(ideal));
    }

    // I = (14 (5 - theta)), with its generator stored; N(5 - theta) = 3 * 41.
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, 70);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -14);
        silex::Element ambient(field);
        silex::OrderElement generator(order);
        silex::Ideal ideal(order);
        assert(ambient.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        assert(generator.set_element(ambient));
        assert(ideal.set_principal(generator));

        silex::IdealFactorization factorization(order);
        silex::Ideal reconstructed(order);
        assert(factorization.factor(ideal));
        assert(factorization.length() == 4);
        assert(count_prime(factorization, 2, 3, 1, 3) == 1);
        assert(count_prime(factorization, 3, 3, 1, 1) == 1);
        assert(count_prime(factorization, 7, 1, 3, 1) == 1);
        assert(count_prime(factorization, 41, 1, 1, 1) == 1);
        assert(norm_accounting_holds(factorization, ideal));
        assert(factorization.reconstruct(reconstructed));
        assert(reconstructed.equal(ideal));
    }
    return 0;
}

int test_factor_rejects_zero_ideal() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    // No public route stores the zero ideal.
    silex::Ideal zero(order);
    silex::OrderElement zero_generator(order);
    sflint::FmpzMat zero_hnf(2, 2);
    assert(zero_generator.set_si(0));
    assert(!zero.set_principal(zero_generator));
    assert(!zero.set_hnf(sflint::FmpzMatConstRef(zero_hnf)));

    // The test-only hook stores it; factor() must reject it, leaving the
    // previous output in place, although fmpz_factor(0) has no primes.
    assert(silex::detail::IdealAccess::set_hnf_direct(
            zero, sflint::FmpzMatConstRef(zero_hnf)));
    sflint::Fmpz norm;
    assert(zero.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_is_zero(sflint::FmpzConstRef(norm)));

    silex::Ideal good(order);
    silex::IdealFactorization factorization(order);
    assert(set_rational_principal(good, 7));
    assert(factorization.factor(good));
    assert(factorization.length() == 2);
    assert(!factorization.factor(zero));
    assert(factorization.length() == 2);
    assert(count_prime(factorization, 7, 1, 1, 1) == 2);
    return 0;
}

int test_keeps_parent_order_alive() {
    silex::IdealFactorization factorization = local_ideal_factorization();
    assert(factorization.is_defined());
    assert(factorization.parent() != nullptr);
    assert(factorization.parent()->parent() != nullptr);
    assert(factorization.length() == 1);

    silex::PrimeIdeal prime(*factorization.parent());
    slong exponent = 0;
    assert(factorization.prime(prime, 0));
    assert(factorization.exponent(exponent, 0));
    assert(exponent == 2);
    assert(silex::same_order_parent(prime.parent(), factorization.parent()));

    silex::Ideal reconstructed(*factorization.parent());
    assert(factorization.reconstruct(reconstructed));
    assert(reconstructed.parent() != nullptr);
    assert(silex::same_order_parent(reconstructed.parent(),
                                    factorization.parent()));
    return 0;
}

}  // namespace

int test_factor_rejects_norm_accounting_mismatch() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::Ideal good(order);
    silex::IdealFactorization factorization(order);
    assert(set_rational_principal(good, 7));
    assert(factorization.factor(good));
    assert(factorization.length() == 2);

    // Overstate one valuation: sum f_P v_P no longer equals v_p(N(I)), so
    // factor() must reject and leave the previous output in place.  The
    // norm accounting runs before the debug-check reconstruct, which
    // therefore cannot mask this.
    silex::detail::TestHookAccess::
            set_ideal_factorization_valuation_perturbation(1);
    silex::Ideal other(order);
    assert(set_rational_principal(other, 3));
    assert(!factorization.factor(other));
    assert(factorization.length() == 2);
    assert(count_prime(factorization, 7, 1, 1, 1) == 2);

    // The perturbation was one-shot: the same call now succeeds.
    assert(factorization.factor(other));
    assert(factorization.length() == 1);
    return 0;
}

int main() {
    assert(test_degree_one_factorization() == 0);
    assert(test_quadratic_factorization() == 0);
    assert(test_factor_failure_preserves_output() == 0);
    assert(test_factor_over_base_degree_one() == 0);
    assert(test_factor_over_base_quadratic() == 0);
    assert(test_factor_over_base_with_required_prime() == 0);
    assert(test_required_prime_through_lacking_conjugate() == 0);
    assert(test_order_element_factor_over_base_with_one_large_prime() == 0);
    assert(test_order_element_factor_one_large_prime_cubic() == 0);
    assert(test_factor_over_base_parent_failure_preserves_outputs() == 0);
    assert(test_fractional_factor_over_base_degree_one() == 0);
    assert(test_fractional_factor_over_base_quadratic() == 0);
    assert(test_fractional_factor_over_incomplete_block() == 0);
    assert(test_fractional_factor_over_base_cancelling_nonbase_primes() == 0);
    assert(test_factor_norm_accounting() == 0);
    assert(test_factor_rejects_zero_ideal() == 0);
    assert(test_factor_rejects_norm_accounting_mismatch() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    return 0;
}
