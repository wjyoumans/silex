#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::Order order_by_polynomial(silex::NumberField& field,
                                  sflint::FmpqPoly& polynomial) noexcept {
    field = field_by_polynomial(polynomial);
    return silex::test::equation_order(field);
}

void fmpz_set_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
}

void zpoly_set_coeff_si(sflint::FmpzPoly& polynomial,
                        slong index,
                        slong value) noexcept {
    fmpz_poly_set_coeff_si(polynomial.raw(), index, value);
}

bool fmpz_poly_coeff_is_si(const sflint::FmpzPoly& polynomial,
                           slong index,
                           slong value) noexcept {
    sflint::Fmpz coeff;
    fmpz_poly_get_coeff_fmpz(coeff.raw(), polynomial.raw(), index);
    return sflint::fmpz_equal_si(coeff, value);
}

bool fmpz_poly_is_const_si(const sflint::FmpzPoly& polynomial,
                           slong value) noexcept {
    return fmpz_poly_degree(polynomial.raw()) <= 0 &&
           fmpz_poly_coeff_is_si(polynomial, 0, value);
}

void element_from_coeff(silex::Element& element,
                        slong index,
                        slong num,
                        slong den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_set_si(coeff, num, den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, index, coeff);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

silex::ResidueField degree_one_residue_field(slong characteristic) noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    sflint::Fmpz p;
    fmpz_set_si(p, characteristic);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*prime));
    return residue_field;
}

silex::ResidueField quadratic_residue_field_25() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    sflint::Fmpz p;
    fmpz_set_si(p, 5);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*prime));
    assert(residue_field.degree() == 2);
    sflint::FmpzPoly modulus;
    assert(residue_field.modulus(sflint::FmpzPolyRef(modulus)));
    assert(fmpz_poly_degree(modulus.raw()) == 2);
    assert(fmpz_poly_coeff_is_si(modulus, 2, 1));
    assert(fmpz_poly_coeff_is_si(modulus, 1, 0));
    assert(fmpz_poly_coeff_is_si(modulus, 0, 3));
    return residue_field;
}

silex::ResidueFieldElement local_residue_field_element() noexcept {
    silex::ResidueField residue_field = degree_one_residue_field(5);
    silex::ResidueFieldElement element(residue_field);
    sflint::FmpzPoly input;
    zpoly_set_coeff_si(input, 0, 2);
    assert(element.set_polynomial(sflint::FmpzPolyConstRef(input)));
    return element;
}

silex::ResidueFieldQuotientLog local_quotient_log() noexcept {
    silex::ResidueField residue_field = degree_one_residue_field(5);
    silex::ResidueFieldQuotientLog cache(residue_field);
    sflint::Fmpz ell;
    fmpz_set_si(ell, 2);
    assert(cache.set_ell(sflint::FmpzConstRef(ell)));
    return cache;
}

int test_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    fmpz_set_si(p, 5);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    silex::ResidueField residue_field;
    silex::ResidueField copy;
    silex::ResidueField other;
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(residue_field.set_prime(*prime));
    assert(residue_field.degree() == 1);

    sflint::Fmpz out_integer;
    assert(residue_field.characteristic(sflint::FmpzRef(out_integer)));
    assert(sflint::fmpz_equal_si(out_integer, 5));
    auto owned_characteristic = residue_field.characteristic();
    assert(owned_characteristic.has_value());
    assert(sflint::fmpz_equal_si(*owned_characteristic, 5));
    assert(residue_field.cardinality(sflint::FmpzRef(out_integer)));
    assert(sflint::fmpz_equal_si(out_integer, 5));
    auto owned_cardinality = residue_field.cardinality();
    assert(owned_cardinality.has_value());
    assert(sflint::fmpz_equal_si(*owned_cardinality, 5));

    sflint::FmpzPoly modulus;
    assert(residue_field.modulus(sflint::FmpzPolyRef(modulus)));
    assert(fmpz_poly_degree(modulus.raw()) == 1);
    auto owned_modulus = residue_field.modulus();
    assert(owned_modulus.has_value());
    assert(fmpz_poly_degree(owned_modulus->raw()) == 1);

    assert(copy.set(residue_field));
    assert(copy.degree() == 1);
    copy.swap(other);
    assert(other.degree() == 1);

    silex::ResidueFieldElement a(residue_field);
    silex::ResidueFieldElement b(residue_field);
    silex::ResidueFieldElement c(residue_field);
    silex::ResidueFieldElement zero(residue_field);
    sflint::FmpzPoly input;
    sflint::FmpzPoly output;

    zpoly_set_coeff_si(input, 0, 7);
    zpoly_set_coeff_si(input, 1, 3);
    assert(a.set_polynomial(sflint::FmpzPolyConstRef(input)));
    assert(a.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 2));
    auto owned_output = a.polynomial();
    assert(owned_output.has_value());
    assert(fmpz_poly_is_const_si(*owned_output, 2));

    assert(b.invert(a));
    assert(c.multiply(a, b));
    assert(c.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    assert(zero.zero());
    assert(b.one());
    assert(!b.invert(zero));
    assert(b.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    assert(a.add(a, a));
    assert(a.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 4));

    assert(a.subtract(a, b));
    assert(a.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 3));

    assert(a.negate(a));
    assert(a.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 2));

    return 0;
}

int test_field_element_keeps_parent_alive() {
    silex::ResidueFieldElement element = local_residue_field_element();
    assert(element.is_defined());
    assert(element.parent() != nullptr);
    assert(element.parent()->parent_order() != nullptr);
    assert(element.parent()->parent_order()->parent() != nullptr);

    sflint::FmpzPoly output;
    assert(element.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 2));

    silex::ResidueFieldElement square(*element.parent());
    assert(square.multiply(element, element));
    assert(square.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 4));

    return 0;
}

int test_quotient_log_keeps_parent_alive() {
    silex::ResidueFieldQuotientLog cache = local_quotient_log();
    assert(cache.is_defined());
    assert(cache.is_set());
    assert(cache.parent() != nullptr);
    assert(cache.parent()->parent_order() != nullptr);

    silex::ResidueFieldElement generator;
    assert(cache.generator(generator));
    assert(generator.parent() != nullptr);
    assert(generator.parent()->equal(*cache.parent()));

    silex::ResidueFieldElement element(*cache.parent());
    sflint::FmpzPoly input;
    zpoly_set_coeff_si(input, 0, 2);
    assert(element.set_polynomial(sflint::FmpzPolyConstRef(input)));

    sflint::Fmpz log;
    assert(cache.apply(sflint::FmpzRef(log), element));
    assert(sflint::fmpz_equal_si(log, 1));

    return 0;
}

int test_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    fmpz_set_si(p, 3);
    silex::PrimeIdealList primes3;
    assert(silex::decompose_prime(primes3, order, sflint::FmpzConstRef(p)));
    assert(primes3.size() == 1);

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*primes3.at(0)));
    assert(residue_field.degree() == 2);
    sflint::Fmpz cardinality;
    assert(residue_field.cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 9));

    silex::ResidueFieldElement theta(residue_field);
    silex::ResidueFieldElement product(residue_field);
    silex::ResidueFieldElement inverse(residue_field);
    silex::Element x(field);
    sflint::FmpzPoly output;
    assert(x.gen());
    assert(theta.set_element(x));
    assert(product.multiply(theta, theta));
    assert(product.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 2));

    assert(inverse.invert(theta));
    assert(product.multiply(theta, inverse));
    assert(product.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    sflint::FmpzPoly input;
    zpoly_set_coeff_si(input, 3, 1);
    zpoly_set_coeff_si(input, 1, 4);
    zpoly_set_coeff_si(input, 0, 8);
    assert(product.set_polynomial(sflint::FmpzPolyConstRef(input)));
    assert(product.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 2));

    silex::Element theta_plus_one_over_two(field);
    {
        sflint::FmpqPoly half_poly;
        sflint::Fmpq half;
        sflint::fmpq_set_si(half, 1, 2);
        sflint::fmpq_poly_set_coeff_fmpq(half_poly, 0, half);
        sflint::fmpq_poly_set_coeff_fmpq(half_poly, 1, half);
        assert(theta_plus_one_over_two.set_fmpq_poly(
                sflint::FmpqPolyConstRef(half_poly)));
    }

    fmpz_set_si(p, 11);
    silex::PrimeIdealList primes11;
    assert(silex::decompose_prime(primes11, order, sflint::FmpzConstRef(p)));
    assert(primes11.size() == 2);
    for (slong i = 0; i < primes11.size(); ++i) {
        assert(residue_field.set_prime(*primes11.at(i)));
        assert(residue_field.degree() == 1);
        assert(residue_field.cardinality(sflint::FmpzRef(cardinality)));
        assert(sflint::fmpz_equal_si(cardinality, 11));

        silex::ResidueFieldElement split_theta(residue_field);
        silex::ResidueFieldElement split_half(residue_field);
        silex::ResidueFieldElement split_product(residue_field);
        silex::ResidueFieldElement split_inverse(residue_field);
        assert(split_theta.set_element(x));
        assert(split_theta.get_polynomial(sflint::FmpzPolyRef(output)));
        assert(fmpz_poly_degree(output.raw()) <= 0);
        const slong r = fmpz_poly_get_coeff_si(output.raw(), 0);
        assert(((r * r - 5) % 11 + 11) % 11 == 0);

        assert(split_half.set_element(theta_plus_one_over_two));
        assert(split_half.get_polynomial(sflint::FmpzPolyRef(output)));
        assert(fmpz_poly_degree(output.raw()) <= 0);
        const slong h = fmpz_poly_get_coeff_si(output.raw(), 0);
        assert(((2 * h - 1 - r) % 11 + 11) % 11 == 0);

        assert(split_product.multiply(split_theta, split_theta));
        assert(split_product.get_polynomial(sflint::FmpzPolyRef(output)));
        assert(fmpz_poly_is_const_si(output, 5));

        assert(split_inverse.invert(split_theta));
        assert(split_product.multiply(split_theta, split_inverse));
        assert(split_product.get_polynomial(sflint::FmpzPolyRef(output)));
        assert(fmpz_poly_is_const_si(output, 1));
    }

    return 0;
}

int test_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    fmpz_set_si(p, 3);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*primes.at(0)));
    sflint::Fmpz cardinality;
    assert(residue_field.cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 9));

    sflint::FmpqPoly unsupported_polynomial;
    poly_x2_minus(unsupported_polynomial, 12);
    silex::NumberField unsupported_field;
    silex::Order unsupported_order;
    unsupported_order = order_by_polynomial(unsupported_field,
                                            unsupported_polynomial);
    assert(unsupported_field.backend_kind() ==
           silex::NumberFieldBackendKind::generic);
    silex::Order maximal_order(unsupported_field);
    assert(maximal_order.maximal_order(unsupported_order));
    fmpz_set_si(p, 2);
    silex::PrimeIdealList explicit_primes;
    assert(silex::decompose_prime(explicit_primes, maximal_order,
                                  sflint::FmpzConstRef(p)));
    assert(!residue_field.set_prime(*explicit_primes.at(0)));
    assert(residue_field.cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 9));

    silex::ResidueFieldElement a(residue_field);
    assert(a.one());
    silex::Element theta_over_three(field);
    element_from_coeff(theta_over_three, 1, 1, 3);
    assert(!a.set_element(theta_over_three));
    sflint::FmpzPoly output;
    assert(a.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    return 0;
}

int test_multiplicative_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    fmpz_set_si(p, 5);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*primes.at(0)));

    silex::ResidueFieldElement two(residue_field);
    silex::ResidueFieldElement four(residue_field);
    silex::ResidueFieldElement x(residue_field);
    silex::ResidueFieldElement y(residue_field);
    silex::ResidueFieldElement g(residue_field);
    silex::ResidueFieldElement undefined_generator;
    silex::ResidueFieldElement quotient_generator(residue_field);
    silex::ResidueFieldElement zero(residue_field);
    silex::ResidueFieldQuotientLog cache(residue_field);
    silex::ResidueFieldQuotientLog cache_copy;
    silex::ResidueFieldQuotientLog cache_other(residue_field);
    silex::ResidueFieldQuotientLog cache_unset(residue_field);
    sflint::FmpzPoly input;
    sflint::FmpzPoly output;
    sflint::Fmpz e;
    sflint::Fmpz ord;
    sflint::Fmpz log;

    zpoly_set_coeff_si(input, 0, 2);
    assert(two.set_polynomial(sflint::FmpzPolyConstRef(input)));
    zpoly_set_coeff_si(input, 0, 4);
    assert(four.set_polynomial(sflint::FmpzPolyConstRef(input)));

    assert(two.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 4));

    fmpz_set_si(e, 3);
    assert(x.pow_fmpz(two, sflint::FmpzConstRef(e)));
    assert(x.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 3));

    fmpz_set_si(e, -1);
    assert(x.pow_fmpz(two, sflint::FmpzConstRef(e)));
    assert(x.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 3));

    fmpz_set_si(e, 2);
    assert(x.pow_fmpz(two, sflint::FmpzConstRef(e)));
    assert(x.pow_fmpz(x, sflint::FmpzConstRef(e)));
    assert(x.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    assert(residue_field.multiplicative_generator(g));
    assert(g.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 4));
    assert(residue_field.multiplicative_generator(undefined_generator));
    assert(undefined_generator.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 4));
    fmpz_set_si(e, 3);
    assert(y.pow_fmpz(g, sflint::FmpzConstRef(e)));
    assert(y.discrete_log(sflint::FmpzRef(log), g));
    assert(sflint::fmpz_equal_si(log, 3));

    assert(four.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 2));
    assert(!two.discrete_log(sflint::FmpzRef(log), four));
    assert(sflint::fmpz_equal_si(log, 3));
    assert(two.degree_one_scalar(sflint::FmpzRef(log)));
    assert(sflint::fmpz_equal_si(log, 2));

    fmpz_set_si(e, 2);
    assert(two.quotient_log_mod_prime(sflint::FmpzRef(log),
                                      sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 1));
    assert(four.quotient_log_mod_prime(sflint::FmpzRef(log),
                                       sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 0));

    fmpz_set_si(log, 66);
    assert(!cache_unset.apply(sflint::FmpzRef(log), two));
    assert(sflint::fmpz_equal_si(log, 66));
    assert(!cache_unset.generator(g));

    assert(cache.set_ell(sflint::FmpzConstRef(e)));
    assert(cache.is_set());
    assert(cache.ell(sflint::FmpzRef(log)));
    assert(sflint::fmpz_equal_si(log, 2));
    assert(cache.generator(g));
    assert(g.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 4));
    assert(cache.quotient_generator(quotient_generator));
    assert(quotient_generator.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 2));
    assert(cache.apply(sflint::FmpzRef(log), two));
    assert(sflint::fmpz_equal_si(log, 1));
    assert(cache.apply(sflint::FmpzRef(log), four));
    assert(sflint::fmpz_equal_si(log, 0));

    assert(cache_copy.set(cache));
    assert(cache_copy.apply(sflint::FmpzRef(log), two));
    assert(sflint::fmpz_equal_si(log, 1));

    fmpz_set_si(log, 77);
    fmpz_set_si(e, 3);
    assert(!two.quotient_log_mod_prime(sflint::FmpzRef(log),
                                       sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 77));
    assert(!cache.set_ell(sflint::FmpzConstRef(e)));
    assert(cache.is_set());
    assert(cache.apply(sflint::FmpzRef(log), two));
    assert(sflint::fmpz_equal_si(log, 1));
    fmpz_set_si(e, 4);
    assert(!two.quotient_log_mod_prime(sflint::FmpzRef(log),
                                       sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 1));
    assert(!cache_copy.set_ell(sflint::FmpzConstRef(e)));
    assert(cache_copy.is_set());

    assert(zero.zero());
    assert(x.one());
    fmpz_set_si(e, -1);
    assert(!x.pow_fmpz(zero, sflint::FmpzConstRef(e)));
    assert(x.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    fmpz_set_si(ord, 99);
    assert(!zero.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 99));

    fmpz_set_si(log, 77);
    assert(!zero.discrete_log(sflint::FmpzRef(log), g));
    assert(sflint::fmpz_equal_si(log, 77));
    fmpz_set_si(e, 2);
    assert(!zero.quotient_log_mod_prime(sflint::FmpzRef(log),
                                        sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 77));
    assert(!cache.apply(sflint::FmpzRef(log), zero));
    assert(sflint::fmpz_equal_si(log, 77));

    cache.swap(cache_other);
    assert(cache_other.apply(sflint::FmpzRef(log), two));
    assert(sflint::fmpz_equal_si(log, 1));

    silex::ResidueFieldElement rational_image(residue_field);
    silex::Element rational(field);
    element_from_coeff(rational, 0, 3, 2);
    assert(rational_image.set_element(rational));
    assert(rational_image.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 4));

    element_from_coeff(rational, 0, 1, 5);
    assert(!rational_image.set_element(rational));
    assert(rational_image.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 4));

    silex::Element two_element(field);
    silex::Element three_element(field);
    silex::Element five_element(field);
    assert(two_element.set_si(2));
    assert(three_element.set_si(3));
    assert(five_element.set_si(5));

    silex::FactoredElement factored(field);
    assert(factored.push(two_element, 3));
    assert(factored.push(three_element, -1));
    assert(rational_image.set_factored_element(factored));
    assert(rational_image.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    assert(rational_image.one());
    assert(factored.one());
    assert(factored.push(five_element, -1));
    assert(!rational_image.set_factored_element(factored));
    assert(rational_image.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_is_const_si(output, 1));

    return 0;
}

int test_multiplicative_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    fmpz_set_si(p, 3);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    silex::ResidueField residue_field;
    assert(residue_field.set_prime(*primes.at(0)));

    silex::ResidueFieldElement g(residue_field);
    silex::ResidueFieldElement h(residue_field);
    silex::ResidueFieldElement one(residue_field);
    silex::ResidueFieldQuotientLog cache(residue_field);
    sflint::Fmpz e;
    sflint::Fmpz ord;
    sflint::Fmpz log;

    assert(residue_field.multiplicative_generator(g));
    assert(g.multiplicative_order(sflint::FmpzRef(ord)));
    assert(sflint::fmpz_equal_si(ord, 8));
    assert(one.one());

    fmpz_set_si(e, 8);
    assert(h.pow_fmpz(g, sflint::FmpzConstRef(e)));
    assert(h.equal(one));

    fmpz_set_si(e, 4);
    assert(h.pow_fmpz(g, sflint::FmpzConstRef(e)));
    assert(!h.equal(one));

    fmpz_set_si(e, 5);
    assert(h.pow_fmpz(g, sflint::FmpzConstRef(e)));
    assert(h.discrete_log(sflint::FmpzRef(log), g));
    assert(sflint::fmpz_equal_si(log, 5));

    fmpz_set_si(e, 2);
    assert(h.quotient_log_mod_prime(sflint::FmpzRef(log),
                                    sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 1));
    assert(cache.set_ell(sflint::FmpzConstRef(e)));
    assert(cache.apply(sflint::FmpzRef(log), h));
    assert(sflint::fmpz_equal_si(log, 1));

    fmpz_set_si(e, 4);
    assert(h.pow_fmpz(g, sflint::FmpzConstRef(e)));
    fmpz_set_si(e, 2);
    assert(h.quotient_log_mod_prime(sflint::FmpzRef(log),
                                    sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 0));
    assert(cache.apply(sflint::FmpzRef(log), h));
    assert(sflint::fmpz_equal_si(log, 0));

    fmpz_set_si(log, 88);
    fmpz_set_si(e, 3);
    assert(!h.quotient_log_mod_prime(sflint::FmpzRef(log),
                                     sflint::FmpzConstRef(e)));
    assert(sflint::fmpz_equal_si(log, 88));
    assert(!cache.set_ell(sflint::FmpzConstRef(e)));
    assert(cache.apply(sflint::FmpzRef(log), h));
    assert(sflint::fmpz_equal_si(log, 0));

    return 0;
}

void check_quotient_logs_mod_three(const silex::ResidueField& residue_field,
                                  const silex::ResidueField& foreign_field,
                                  slong characteristic,
                                  slong cardinality,
                                  slong primary_order) noexcept {
    sflint::Fmpz value;
    assert(residue_field.characteristic(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, characteristic));
    assert(residue_field.cardinality(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, cardinality));
    assert(foreign_field.characteristic(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, characteristic));
    assert(foreign_field.cardinality(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, cardinality));
    assert(!foreign_field.equal(residue_field));

    sflint::Fmpz ell;
    fmpz_set_si(ell, 3);
    silex::ResidueFieldQuotientLog cache(residue_field);
    assert(cache.set_ell(sflint::FmpzConstRef(ell)));
    silex::ResidueFieldElement generator(residue_field);
    silex::ResidueFieldElement quotient_generator(residue_field);
    assert(cache.generator(generator));
    assert(generator.multiplicative_order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, primary_order));
    assert(cache.quotient_generator(quotient_generator));
    assert(quotient_generator.multiplicative_order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 3));

    silex::ResidueFieldElement power(residue_field);
    sflint::Fmpz exponent;
    fmpz_set_si(exponent, primary_order / 3);
    assert(power.pow_fmpz(generator, sflint::FmpzConstRef(exponent)));
    assert(power.equal(quotient_generator));

    // The stateless map uses a full multiplicative generator and therefore
    // has its own normalization, independent of the cached primary generator.
    silex::ResidueFieldElement full_generator(residue_field);
    silex::ResidueFieldElement stateless_generator(residue_field);
    assert(residue_field.multiplicative_generator(full_generator));
    assert(full_generator.multiplicative_order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, cardinality - 1));
    sflint::Fmpz projection_exponent;
    fmpz_set_si(projection_exponent, (cardinality - 1) / 3);
    assert(stateless_generator.pow_fmpz(
            full_generator, sflint::FmpzConstRef(projection_exponent)));
    assert(stateless_generator.multiplicative_order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 3));

    const auto cached_log = [&cache, &ell](
            const silex::ResidueFieldElement& element) noexcept {
        sflint::Fmpz log;
        assert(cache.apply(sflint::FmpzRef(log), element));
        assert(sflint::fmpz_sgn(sflint::FmpzConstRef(log)) >= 0);
        assert(sflint::fmpz_cmp(sflint::FmpzConstRef(log),
                               sflint::FmpzConstRef(ell)) < 0);
        return sflint::fmpz_get_si(sflint::FmpzConstRef(log));
    };
    const auto stateless_log = [&ell](
            const silex::ResidueFieldElement& element) noexcept {
        sflint::Fmpz log;
        assert(element.quotient_log_mod_prime(sflint::FmpzRef(log),
                                             sflint::FmpzConstRef(ell)));
        assert(sflint::fmpz_sgn(sflint::FmpzConstRef(log)) >= 0);
        assert(sflint::fmpz_cmp(sflint::FmpzConstRef(log),
                               sflint::FmpzConstRef(ell)) < 0);
        return sflint::fmpz_get_si(sflint::FmpzConstRef(log));
    };

    // For q - 1 = c * 3^v, the defining projection gives L(generator) = c
    // modulo 3.  The cofactors are 2, 2, and 8 for F7, F19, and F25.
    assert((cardinality - 1) / primary_order % 3 == 2);
    assert(cached_log(generator) == 2);
    assert(stateless_log(full_generator) == 1);

    // Enumerate constants in F7/F19 and a + b*theta in F25, theta^2 = 2.
    const auto set_nonzero = [characteristic, cardinality](
            silex::ResidueFieldElement& element, slong index) noexcept {
        assert(index > 0 && index < cardinality);
        sflint::FmpzPoly input;
        zpoly_set_coeff_si(input, 0, index % characteristic);
        zpoly_set_coeff_si(input, 1, index / characteristic);
        assert(element.set_polynomial(sflint::FmpzPolyConstRef(input)));
        assert(!element.is_zero());
    };
    silex::ResidueFieldElement x(residue_field);
    silex::ResidueFieldElement y(residue_field);
    silex::ResidueFieldElement product(residue_field);
    silex::ResidueFieldElement inverse(residue_field);
    silex::ResidueFieldElement image(residue_field);
    silex::ResidueFieldElement one(residue_field);
    assert(one.one());
    slong kernel_size = 0;
    for (slong i = 1; i < cardinality; ++i) {
        set_nonzero(x, i);
        const slong log_x = cached_log(x);
        const slong stateless_log_x = stateless_log(x);
        assert(image.pow_fmpz(x, sflint::FmpzConstRef(projection_exponent)));
        fmpz_set_si(exponent, log_x);
        assert(power.pow_fmpz(quotient_generator,
                              sflint::FmpzConstRef(exponent)));
        assert(power.equal(image));
        fmpz_set_si(exponent, stateless_log_x);
        assert(power.pow_fmpz(stateless_generator,
                              sflint::FmpzConstRef(exponent)));
        assert(power.equal(image));
        assert((log_x == 0) == (stateless_log_x == 0));
        assert((log_x == 0) == image.equal(one));
        if (log_x == 0) {
            ++kernel_size;
        }

        assert(inverse.invert(x));
        assert(cached_log(inverse) == (3 - log_x) % 3);
        assert(stateless_log(inverse) == (3 - stateless_log_x) % 3);
        for (slong j = 1; j < cardinality; ++j) {
            set_nonzero(y, j);
            assert(product.multiply(x, y));
            assert(cached_log(product) == (log_x + cached_log(y)) % 3);
            assert(stateless_log(product) ==
                   (stateless_log_x + stateless_log(y)) % 3);
        }
    }
    assert(kernel_size == (cardinality - 1) / 3);

    const auto check_cache_preserved = [&]() noexcept {
        assert(cache.is_set());
        assert(cache.parent() != nullptr);
        assert(cache.parent()->equal(residue_field));
        assert(cache.ell(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, 3));
        assert(cache.generator(power));
        assert(power.equal(generator));
        assert(cache.quotient_generator(power));
        assert(power.equal(quotient_generator));
        assert(cached_log(generator) == 2);
    };
    silex::ResidueFieldElement zero(residue_field);
    assert(zero.zero());
    fmpz_set_si(value, 77);
    assert(!cache.apply(sflint::FmpzRef(value), zero));
    assert(sflint::fmpz_equal_si(value, 77));
    assert(!zero.quotient_log_mod_prime(sflint::FmpzRef(value),
                                       sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_equal_si(value, 77));
    check_cache_preserved();

    silex::ResidueFieldElement foreign(foreign_field);
    assert(foreign.one());
    fmpz_set_si(value, 77);
    assert(!cache.apply(sflint::FmpzRef(value), foreign));
    assert(sflint::fmpz_equal_si(value, 77));
    check_cache_preserved();

    const slong invalid_ells[] = {9, 5};
    for (slong invalid_ell : invalid_ells) {
        sflint::Fmpz bad_ell;
        fmpz_set_si(bad_ell, invalid_ell);
        fmpz_set_si(value, 77);
        assert(!generator.quotient_log_mod_prime(sflint::FmpzRef(value),
                                                sflint::FmpzConstRef(bad_ell)));
        assert(sflint::fmpz_equal_si(value, 77));
        assert(!cache.set_ell(sflint::FmpzConstRef(bad_ell)));
        check_cache_preserved();
    }
}

int test_quotient_logs_mod_three() {
    const slong characteristics[] = {7, 19};
    for (slong characteristic : characteristics) {
        silex::ResidueField residue_field =
                degree_one_residue_field(characteristic);
        silex::ResidueField foreign_field =
                degree_one_residue_field(characteristic);
        assert(residue_field.degree() == 1);
        check_quotient_logs_mod_three(residue_field, foreign_field,
                                      characteristic, characteristic,
                                      characteristic == 19 ? 9 : 3);
    }

    silex::ResidueField residue_field = quadratic_residue_field_25();
    silex::ResidueField foreign_field = quadratic_residue_field_25();
    check_quotient_logs_mod_three(residue_field, foreign_field, 5, 25, 3);
    return 0;
}

void element_linear_over(silex::Element& out,
                         slong constant,
                         slong alpha,
                         slong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, alpha);
    fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), denominator);
    assert(out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

slong residue_constant(const silex::ResidueFieldElement& element) noexcept {
    sflint::FmpzPoly output;
    assert(element.get_polynomial(sflint::FmpzPolyRef(output)));
    assert(fmpz_poly_degree(output.raw()) <= 0);
    return fmpz_poly_get_coeff_si(output.raw(), 0);
}

int test_generic_nonpower_basis_set_element() {
    // x^2 - 12 on the generic backend: the maximal order [1, alpha/2] is not
    // a power basis, and the residue polynomial is stored in alpha.
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 12);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    assert(field.backend_kind() == silex::NumberFieldBackendKind::generic);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!maximal.is_equation_order());

    silex::Element alpha(field);
    silex::Element half_alpha(field);
    silex::Element third_alpha(field);
    silex::Element alpha_over_p(field);
    silex::Element epsilon(field);
    silex::Element epsilon_squared(field);
    assert(alpha.gen());
    element_linear_over(half_alpha, 0, 1, 2);
    element_linear_over(third_alpha, 0, 1, 3);
    element_linear_over(alpha_over_p, 0, 1, 11);
    element_linear_over(epsilon, 4, 1, 2);  // 2 + alpha/2
    assert(epsilon_squared.multiply(epsilon, epsilon));

    silex::OrderElement half_order(maximal);
    assert(half_order.set_element(half_alpha));

    sflint::Fmpz p;
    fmpz_set_si(p, 11);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    bool saw_plus_one = false;
    bool saw_minus_one = false;
    for (slong i = 0; i < primes.size(); ++i) {
        silex::ResidueField residue_field;
        assert(residue_field.set_prime(*primes.at(i)));
        assert(residue_field.degree() == 1);

        silex::ResidueFieldElement image(residue_field);
        assert(image.set_element(alpha));
        const slong root = residue_constant(image);
        assert(root == 1 || root == 10);
        const bool plus = root == 1;
        saw_plus_one = saw_plus_one || plus;
        saw_minus_one = saw_minus_one || !plus;

        assert(image.set_element(half_alpha));
        assert(residue_constant(image) == (plus ? 6 : 5));
        assert(image.set_order_element(half_order));
        assert(residue_constant(image) == (plus ? 6 : 5));

        // Denominators prime to p are inverted: 1/3 = 4, -1/3 = 7 mod 11.
        assert(image.set_element(third_alpha));
        assert(residue_constant(image) == (plus ? 4 : 7));

        // epsilon = 2 + alpha/2 -> 8 or 7; epsilon^2 = 7 + 2 alpha -> 9 or 5.
        assert(image.set_element(epsilon));
        assert(residue_constant(image) == (plus ? 8 : 7));
        assert(image.set_element(epsilon_squared));
        assert(residue_constant(image) == (plus ? 9 : 5));

        // A denominator divisible by p fails and preserves the output.
        assert(!image.set_element(alpha_over_p));
        assert(residue_constant(image) == (plus ? 9 : 5));
    }
    assert(saw_plus_one && saw_minus_one);
    return 0;
}

void poly_cubic_disc1724(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -8);
}

int test_set_element_homomorphism_on_index_maximal_order() {
    // Cubic field with [O : Z[alpha]] = 2; check that set_element is
    // multiplicative and additive, including elements with denominators
    // prime to p, at every prime above p not dividing the index.
    sflint::FmpqPoly polynomial;
    poly_cubic_disc1724(polynomial);
    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    sflint::FmpqMat basis(3, 3);
    assert(maximal.get_basis(sflint::FmpqMatRef(basis)));
    std::vector<silex::Element> elements;
    for (slong i = 0; i < 3; ++i) {
        sflint::FmpqPoly row;
        for (slong j = 0; j < 3; ++j) {
            fmpq_poly_set_coeff_fmpq(row.raw(), j,
                                     fmpq_mat_entry(basis.raw(), i, j));
        }
        elements.emplace_back(field);
        assert(elements.back().set_fmpq_poly(sflint::FmpqPolyConstRef(row)));
    }
    {
        silex::Element mixed(field);
        silex::Element scaled(field);
        assert(mixed.add(elements[1], elements[2]));
        assert(scaled.scalar_div_si(mixed, 3));
        elements.push_back(std::move(scaled));
        silex::Element shifted(field);
        silex::Element two(field);
        assert(two.set_si(2));
        assert(shifted.subtract(elements[2], two));
        elements.push_back(std::move(shifted));
    }

    const slong rational_primes[] = {3, 5, 7, 11, 13};
    slong checked = 0;
    for (const slong q : rational_primes) {
        sflint::Fmpz p;
        fmpz_set_si(p, q);
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, maximal,
                                      sflint::FmpzConstRef(p)));
        for (slong i = 0; i < primes.size(); ++i) {
            silex::ResidueField residue_field;
            assert(residue_field.set_prime(*primes.at(i)));
            silex::ResidueFieldElement left_image(residue_field);
            silex::ResidueFieldElement right_image(residue_field);
            silex::ResidueFieldElement combined_image(residue_field);
            silex::ResidueFieldElement expected(residue_field);
            for (const silex::Element& left : elements) {
                if (q == 3 && &left == &elements[3]) {
                    assert(!left_image.set_element(left));
                    continue;
                }
                assert(left_image.set_element(left));
                silex::OrderElement left_order(maximal);
                if (left_order.set_element(left)) {
                    silex::ResidueFieldElement order_image(residue_field);
                    assert(order_image.set_order_element(left_order));
                    assert(order_image.equal(left_image));
                }
                for (const silex::Element& right : elements) {
                    if (q == 3 && &right == &elements[3]) {
                        continue;
                    }
                    silex::Element product(field);
                    silex::Element sum(field);
                    assert(product.multiply(left, right));
                    assert(sum.add(left, right));
                    assert(right_image.set_element(right));

                    assert(combined_image.set_element(product));
                    assert(expected.multiply(left_image, right_image));
                    assert(combined_image.equal(expected));

                    assert(combined_image.set_element(sum));
                    assert(expected.add(left_image, right_image));
                    assert(combined_image.equal(expected));
                }
            }
            ++checked;
        }
    }
    assert(checked >= 5);
    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_field_element_keeps_parent_alive() == 0);
    assert(test_quotient_log_keeps_parent_alive() == 0);
    assert(test_quadratic() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_multiplicative_degree_one() == 0);
    assert(test_multiplicative_quadratic() == 0);
    assert(test_quotient_logs_mod_three() == 0);
    assert(test_generic_nonpower_basis_set_element() == 0);
    assert(test_set_element_homomorphism_on_index_maximal_order() == 0);
    return 0;
}
