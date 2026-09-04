#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal_factorization.hpp>

#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "test_support.hpp"

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
        } else if (!base.contains(factor)) {
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

silex::IdealFactorization local_ideal_factorization() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 4));

    silex::IdealFactorization factorization(order);
    assert(factorization.factor(ideal));
    return factorization;
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

int main() {
    assert(test_degree_one_factorization() == 0);
    assert(test_quadratic_factorization() == 0);
    assert(test_factor_failure_preserves_output() == 0);
    assert(test_factor_over_base_degree_one() == 0);
    assert(test_factor_over_base_quadratic() == 0);
    assert(test_factor_over_base_with_required_prime() == 0);
    assert(test_order_element_factor_over_base_with_one_large_prime() == 0);
    assert(test_factor_over_base_parent_failure_preserves_outputs() == 0);
    assert(test_fractional_factor_over_base_degree_one() == 0);
    assert(test_fractional_factor_over_base_quadratic() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    return 0;
}
