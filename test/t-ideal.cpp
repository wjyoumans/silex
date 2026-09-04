#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

bool same_parent(const silex::Order* parent,
                 const silex::Order& order) noexcept {
    return parent != nullptr && parent->has_same_data(order);
}

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

void mat_entry_si(sflint::FmpzMat& matrix,
                  slong row,
                  slong col,
                  slong value) noexcept {
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(matrix, row, col), value);
}

void qmat_entry_si(sflint::FmpqMat& matrix,
                   slong row,
                   slong col,
                   slong value) noexcept {
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(matrix, row, col), value, 1);
}

void qmat_entry_frac_si(sflint::FmpqMat& matrix,
                        slong row,
                        slong col,
                        slong num,
                        slong den) noexcept {
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(matrix, row, col), num, den);
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool mat_is_diag2(const sflint::FmpzMat& matrix, slong diagonal) noexcept {
    return mat_entry_is_si(matrix, 0, 0, diagonal) &&
           mat_entry_is_si(matrix, 0, 1, 0) &&
           mat_entry_is_si(matrix, 1, 0, 0) &&
           mat_entry_is_si(matrix, 1, 1, diagonal);
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

void element_set_pseudorandom_integral(silex::Element& element,
                                       slong seed) noexcept {
    sflint::FmpqPoly polynomial;
    for (slong j = 0; j < 3; ++j) {
        const slong c = ((seed * (13 + 5 * j) + 7 * j * j + 3) % 11) - 5;
        sflint::fmpq_poly_set_coeff_si(polynomial, j, c);
    }
    if (fmpq_poly_is_zero(polynomial.raw()) != 0) {
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, seed + 1);
    }
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

void element_set_fractional_pseudorandom(silex::Element& element,
                                         slong seed,
                                         slong den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    for (slong j = 0; j < 3; ++j) {
        const slong c = ((seed * (17 + 3 * j) + 5 * j * j + 9) % 13) - 6;
        sflint::fmpq_set_si(coeff, c, den);
        sflint::fmpq_poly_set_coeff_fmpq(polynomial, j, coeff);
    }
    if (fmpq_poly_is_zero(polynomial.raw()) != 0) {
        sflint::fmpq_set_si(coeff, seed + 1, den);
        sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coeff);
    }
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

void element_one_plus_theta_over_two(silex::Element& element) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_set_si(coeff, 1, 2);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coeff);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coeff);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool fmpq_is_si(const sflint::Fmpq& value,
                slong num,
                slong den) noexcept {
    sflint::Fmpq expected;
    sflint::fmpq_set_si(expected, num, den);
    return sflint::fmpq_equal(value, expected);
}

silex::FractionalIdeal fractional_from_hnf(const silex::Order& order,
                                           slong hnf_entry,
                                           slong den_value) noexcept {
    silex::Ideal numerator(order);
    sflint::FmpzMat hnf(1, 1);
    mat_entry_si(hnf, 0, 0, hnf_entry);
    assert(numerator.set_hnf(sflint::FmpzMatConstRef(hnf)));

    silex::FractionalIdeal ideal(order);
    sflint::Fmpz den;
    sflint::fmpz_set_ui(den, static_cast<ulong>(den_value));
    assert(ideal.set_integral_den(numerator, sflint::FmpzConstRef(den)));
    return ideal;
}

bool set_rational_principal(silex::Ideal& ideal, slong value) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr) {
        return false;
    }

    silex::OrderElement generator(*order);
    return generator.set_si(value) && ideal.set_principal(generator);
}

int test_degree_one_integral_arithmetic() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal two(order);
    silex::Ideal three(order);
    silex::Ideal out(order);
    sflint::FmpzMat hnf(1, 1);
    sflint::Fmpz norm;

    mat_entry_si(hnf, 0, 0, 2);
    assert(two.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(two.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 2));

    silex::OrderElement value(order);
    assert(value.set_si(4));
    assert(two.contains(value));
    assert(value.set_si(3));
    assert(!two.contains(value));

    assert(value.set_si(3));
    assert(three.set_principal(value));
    assert(three.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 3));

    silex::Ideal negative_six(order);
    silex::Ideal six_from_principal(order);
    mat_entry_si(hnf, 0, 0, -6);
    assert(negative_six.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(value.set_si(6));
    assert(six_from_principal.set_principal(value));
    assert(negative_six.equal(six_from_principal));
    assert(negative_six.get_hnf(sflint::FmpzMatRef(hnf)));
    assert(mat_entry_is_si(hnf, 0, 0, 6));

    silex::Ideal before(order);
    assert(before.set(three));
    assert(value.set_si(0));
    assert(!three.set_principal(value));
    assert(three.equal(before));

    silex::Ideal four(order);
    assert(value.set_si(4));
    assert(four.set_principal(value));
    assert(two.contains(four));
    assert(!four.contains(two));
    assert(!two.contains(three));

    assert(out.add(two, three));
    assert(out.is_one());
    assert(out.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 1));

    assert(out.intersect(two, three));
    assert(out.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 6));
    assert(out.get_hnf(sflint::FmpzMatRef(hnf)));
    assert(mat_entry_is_si(hnf, 0, 0, 6));

    assert(out.multiply(two, three));
    assert(out.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 6));

    assert(two.multiply(two, three));
    assert(two.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 6));
    return 0;
}

int test_integral_arithmetic_nongalois_cubic() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_element(field);
    silex::Element theta_squared_element(field);
    assert(theta_element.gen());
    assert(theta_squared_element.multiply(theta_element, theta_element));

    silex::OrderElement theta(order);
    silex::OrderElement theta_squared(order);
    silex::OrderElement two(order);
    assert(theta.set_element(theta_element));
    assert(theta_squared.set_element(theta_squared_element));
    assert(two.set_si(2));

    silex::Ideal left(order);
    silex::Ideal right(order);
    silex::Ideal out(order);
    silex::Ideal expected(order);
    sflint::Fmpz norm_left;
    sflint::Fmpz norm_right;
    sflint::Fmpz norm_product;

    assert(left.set_principal(theta));
    assert(right.set_principal(theta_squared));

    assert(out.add(left, right));
    assert(out.equal(left));

    assert(out.intersect(left, right));
    assert(out.equal(right));

    assert(out.multiply(left, right));
    assert(expected.set_principal(two));
    assert(out.equal(expected));

    assert(left.norm(sflint::FmpzRef(norm_left)));
    assert(right.norm(sflint::FmpzRef(norm_right)));
    assert(out.norm(sflint::FmpzRef(norm_product)));
    assert(sflint::fmpz_equal_si(norm_left, 2));
    assert(sflint::fmpz_equal_si(norm_right, 4));
    assert(sflint::fmpz_equal_si(norm_product, 8));
    sflint::fmpz_mul(sflint::FmpzRef(norm_left),
                     sflint::FmpzConstRef(norm_left),
                     sflint::FmpzConstRef(norm_right));
    assert(sflint::fmpz_equal(norm_product, norm_left));

    assert(right.contains(two));
    assert(!right.contains(theta));
    return 0;
}

int test_integral_cubic_swap_and_self_set() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal left(order);
    silex::Ideal right(order);
    silex::OrderElement three(order);
    sflint::Fmpz norm;

    assert(left.one());
    assert(three.set_si(3));
    assert(right.set_principal(three));

    left.swap(right);
    assert(left.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 27));
    assert(right.is_one());

    assert(left.set(left));
    assert(left.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 27));
    return 0;
}

int test_integral_arithmetic_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_element(field);
    assert(theta_element.gen());

    silex::OrderElement theta(order);
    silex::OrderElement two(order);
    assert(theta.set_element(theta_element));
    assert(two.set_si(2));

    silex::Ideal theta_ideal(order);
    silex::Ideal two_ideal(order);
    silex::Ideal out(order);
    silex::Ideal unit(order);
    sflint::Fmpz norm_theta;
    sflint::Fmpz norm_two;
    sflint::Fmpz norm_product;

    assert(theta_ideal.set_principal(theta));
    assert(two_ideal.set_principal(two));
    assert(unit.one());

    assert(out.add(theta_ideal, two_ideal));
    assert(out.equal(theta_ideal));

    assert(out.intersect(theta_ideal, two_ideal));
    assert(out.equal(two_ideal));

    assert(out.multiply(theta_ideal, theta_ideal));
    assert(out.equal(two_ideal));

    assert(theta_ideal.norm(sflint::FmpzRef(norm_theta)));
    assert(two_ideal.norm(sflint::FmpzRef(norm_two)));
    assert(out.multiply(theta_ideal, two_ideal));
    assert(out.norm(sflint::FmpzRef(norm_product)));
    sflint::fmpz_mul(sflint::FmpzRef(norm_theta),
                     sflint::FmpzConstRef(norm_theta),
                     sflint::FmpzConstRef(norm_two));
    assert(sflint::fmpz_equal(norm_product, norm_theta));

    assert(out.multiply(out, unit));
    assert(out.norm(sflint::FmpzRef(norm_product)));
    assert(sflint::fmpz_equal(norm_product, norm_theta));
    return 0;
}

int test_integral_arithmetic_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta_element(field);
    assert(theta_element.gen());

    silex::OrderElement theta(order);
    silex::OrderElement five(order);
    assert(theta.set_element(theta_element));
    assert(five.set_si(5));

    silex::Ideal theta_ideal(order);
    silex::Ideal five_ideal(order);
    silex::Ideal out(order);
    sflint::Fmpz norm_left;
    sflint::Fmpz norm_right;
    sflint::Fmpz norm_product;

    assert(theta_ideal.set_principal(theta));
    assert(five_ideal.set_principal(five));

    assert(out.add(theta_ideal, five_ideal));
    assert(out.equal(theta_ideal));

    assert(out.intersect(theta_ideal, five_ideal));
    assert(out.equal(five_ideal));

    assert(out.multiply(theta_ideal, theta_ideal));
    assert(out.equal(five_ideal));

    assert(theta_ideal.norm(sflint::FmpzRef(norm_left)));
    assert(five_ideal.norm(sflint::FmpzRef(norm_right)));
    assert(out.multiply(theta_ideal, five_ideal));
    assert(out.norm(sflint::FmpzRef(norm_product)));
    sflint::fmpz_mul(sflint::FmpzRef(norm_left),
                     sflint::FmpzConstRef(norm_left),
                     sflint::FmpzConstRef(norm_right));
    assert(sflint::fmpz_equal(norm_product, norm_left));
    return 0;
}

int test_integral_principal_pseudorandom_identities() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element a(field);
    silex::Element b(field);
    silex::Element ab(field);
    silex::OrderElement a_order(order);
    silex::OrderElement b_order(order);
    silex::OrderElement ab_order(order);
    silex::Ideal left(order);
    silex::Ideal right(order);
    silex::Ideal product(order);
    silex::Ideal expected(order);
    sflint::Fmpz norm_left;
    sflint::Fmpz norm_right;
    sflint::Fmpz norm_product;

    for (slong t = 1; t <= 10; ++t) {
        element_set_pseudorandom_integral(a, t);
        element_set_pseudorandom_integral(b, 17 + 3 * t);
        assert(ab.multiply(a, b));

        assert(a_order.set_element(a));
        assert(b_order.set_element(b));
        assert(ab_order.set_element(ab));
        assert(left.set_principal(a_order));
        assert(right.set_principal(b_order));
        assert(expected.set_principal(ab_order));

        assert(product.multiply(left, right));
        assert(product.equal(expected));

        assert(left.norm(sflint::FmpzRef(norm_left)));
        assert(right.norm(sflint::FmpzRef(norm_right)));
        assert(product.norm(sflint::FmpzRef(norm_product)));
        sflint::fmpz_mul(sflint::FmpzRef(norm_left),
                         sflint::FmpzConstRef(norm_left),
                         sflint::FmpzConstRef(norm_right));
        assert(sflint::fmpz_equal(norm_product, norm_left));
    }

    return 0;
}

int test_principal_times_integral_direct_matches_fractional_product() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta(field);
    silex::Element theta_over_two(field);
    assert(theta.gen());
    element_from_coeff(theta_over_two, 1, 1, 2);

    silex::OrderElement two_element(order);
    silex::OrderElement theta_order(order);
    assert(two_element.set_si(2));
    assert(theta_order.set_element(theta));

    silex::Ideal two_ideal(order);
    silex::Ideal expected(order);
    silex::Ideal direct(order);
    silex::Ideal generic(order);
    assert(two_ideal.set_principal(two_element));
    assert(expected.set_principal(theta_order));

    assert(silex::detail::multiply_integral_ideal_by_element(
            direct, two_ideal, theta_over_two));
    assert(direct.equal(expected));

    silex::FractionalIdeal principal(order);
    silex::FractionalIdeal integral(order);
    silex::FractionalIdeal product(order);
    sflint::Fmpz denominator;
    assert(principal.set_principal(theta_over_two));
    assert(integral.set_integral(two_ideal));
    assert(product.multiply(principal, integral));
    assert(product.get_integral_den(generic, sflint::FmpzRef(denominator)));
    assert(sflint::fmpz_is_one(denominator));
    assert(direct.equal(generic));

    return 0;
}

int test_two_generator_product_matches_generic() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta(field);
    silex::OrderElement theta_order(order);
    silex::OrderElement one(order);
    silex::OrderElement theta_plus_one(order);
    assert(theta.gen());
    assert(theta_order.set_element(theta));
    assert(one.set_si(1));
    assert(theta_plus_one.add(theta_order, one));

    sflint::Fmpz two;
    sflint::Fmpz three;
    sflint::fmpz_set_ui(two, 2);
    sflint::fmpz_set_ui(three, 3);

    silex::Ideal left(order);
    silex::Ideal right(order);
    silex::Ideal generic(order);
    silex::Ideal direct(order);
    silex::Ideal before(order);
    assert(silex::detail::set_known_two_generator_ideal(
            left, sflint::FmpzConstRef(two), theta_plus_one));
    assert(silex::detail::set_known_two_generator_ideal(
            right, sflint::FmpzConstRef(three), theta_order));
    assert(generic.multiply(left, right));
    assert(silex::detail::multiply_integral_ideal_by_two_generator(
            direct, left, sflint::FmpzConstRef(three), theta_order));
    assert(direct.equal(generic));

    assert(before.set(direct));
    assert(silex::detail::multiply_integral_ideal_by_two_generator(
            left, left, sflint::FmpzConstRef(three), theta_order));
    assert(left.equal(generic));

    sflint::Fmpz zero;
    assert(!silex::detail::multiply_integral_ideal_by_two_generator(
            direct, generic, sflint::FmpzConstRef(zero), theta_order));
    assert(direct.equal(before));
    return 0;
}

int test_principal_times_integral_direct_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta(field);
    silex::Element theta_over_two(field);
    assert(theta.gen());
    element_from_coeff(theta_over_two, 1, 1, 2);

    silex::OrderElement two_element(order);
    silex::OrderElement theta_order(order);
    assert(two_element.set_si(2));
    assert(theta_order.set_element(theta));

    silex::Ideal two_ideal(order);
    silex::Ideal expected(order);
    silex::Ideal direct(order);
    assert(two_ideal.set_principal(two_element));
    assert(expected.set_principal(theta_order));
    assert(silex::detail::multiply_integral_ideal_by_element(
            direct, two_ideal, theta_over_two));
    assert(direct.equal(expected));

    return 0;
}

int test_principal_times_integral_direct_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_over_three(field);
    element_from_coeff(theta_over_three, 1, 1, 3);

    silex::OrderElement two_element(order);
    assert(two_element.set_si(2));

    silex::Ideal two_ideal(order);
    silex::Ideal direct(order);
    silex::Ideal before(order);
    assert(two_ideal.set_principal(two_element));
    assert(direct.one());
    assert(before.one());

    assert(!silex::detail::multiply_integral_ideal_by_element(
            direct, two_ideal, theta_over_three));
    assert(direct.equal(before));

    return 0;
}

int test_integral_coprime_add_to_one_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal three(order);
    silex::Ideal five(order);
    silex::Ideal six(order);
    silex::Ideal ten(order);
    silex::Ideal thirty_five(order);
    assert(set_rational_principal(three, 3));
    assert(set_rational_principal(five, 5));
    assert(set_rational_principal(six, 6));
    assert(set_rational_principal(ten, 10));
    assert(set_rational_principal(thirty_five, 35));

    bool coprime = false;
    assert(six.is_coprime(coprime, thirty_five));
    assert(coprime);
    assert(six.is_coprime(coprime, ten));
    assert(!coprime);

    silex::OrderElement left(order);
    silex::OrderElement right(order);
    silex::OrderElement sum(order);
    assert(three.add_to_one(left, right, five));
    assert(three.contains(left));
    assert(five.contains(right));
    assert(sum.add(left, right));
    assert(sum.equal_si(1));

    assert(left.set_si(123));
    assert(right.set_si(456));
    assert(!six.add_to_one(left, right, ten));
    assert(left.equal_si(123));
    assert(right.equal_si(456));

    return 0;
}

int test_integral_add_to_one_swapped_quadratic_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 1, 1);
    qmat_entry_si(basis, 1, 0, 1);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Ideal three(order);
    silex::Ideal five(order);
    assert(set_rational_principal(three, 3));
    assert(set_rational_principal(five, 5));

    bool coprime = false;
    assert(three.is_coprime(coprime, five));
    assert(coprime);

    silex::OrderElement left(order);
    silex::OrderElement right(order);
    silex::OrderElement sum(order);
    assert(three.add_to_one(left, right, five));
    assert(three.contains(left));
    assert(five.contains(right));
    assert(sum.add(left, right));
    assert(sum.equal_si(1));

    return 0;
}

int test_multiplier_ring_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal principal(order);
    silex::OrderElement seven(order);
    assert(seven.set_si(7));
    assert(principal.set_principal(seven));

    silex::Order ring(field);
    sflint::Fmpz index;
    assert(principal.multiplier_ring(ring));
    assert(order_index(sflint::FmpzRef(index), order, ring));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(order_index(sflint::FmpzRef(index), ring, order));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(!ring.maximality_known());
    return 0;
}

int test_p_radical_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal radical(order);
    sflint::Fmpz p;
    sflint::Fmpz norm;
    sflint::FmpzMat hnf(1, 1);
    sflint::fmpz_set_ui(p, 7);
    assert(order.p_radical(radical, sflint::FmpzConstRef(p)));
    assert(radical.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 7));
    assert(radical.get_hnf(sflint::FmpzMatRef(hnf)));
    assert(mat_entry_is_si(hnf, 0, 0, 7));
    return 0;
}

int test_multiplier_ring_quadratic_p_radical_shape() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal radical(order);
    sflint::Fmpz p;
    sflint::fmpz_set_ui(p, 2);
    assert(order.p_radical(radical, sflint::FmpzConstRef(p)));

    silex::Order ring(field);
    sflint::Fmpz index;
    sflint::Fmpz discriminant;
    assert(radical.multiplier_ring(ring));
    assert(order_index(sflint::FmpzRef(index), order, ring));
    assert(sflint::fmpz_equal_si(index, 2));
    assert(ring.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, 5));

    silex::Element alpha(field);
    element_one_plus_theta_over_two(alpha);
    assert(ring.contains(alpha));
    assert(!ring.maximality_known());
    return 0;
}

int test_p_radical_stable_multiplier_rings() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField sqrt2;
    sqrt2 = field_by_polynomial(polynomial);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);

    silex::Ideal radical(sqrt2_order);
    silex::Order ring(sqrt2);
    sflint::Fmpz p;
    sflint::Fmpz index;
    sflint::Fmpz discriminant;
    sflint::fmpz_set_ui(p, 2);
    assert(sqrt2_order.p_radical(radical, sflint::FmpzConstRef(p)));
    assert(radical.multiplier_ring(ring));
    assert(order_index(sflint::FmpzRef(index), sqrt2_order, ring));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(order_index(sflint::FmpzRef(index), ring, sqrt2_order));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(ring.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, 8));

    sflint::FmpqPoly polynomial5;
    poly_x2_minus(polynomial5, 5);
    silex::NumberField sqrt5;
    sqrt5 = field_by_polynomial(polynomial5);
    silex::Order sqrt5_order;
    sqrt5_order = silex::test::equation_order(sqrt5);
    silex::Ideal radical3(sqrt5_order);
    silex::Order ring3(sqrt5);
    sflint::fmpz_set_ui(p, 3);
    assert(sqrt5_order.p_radical(radical3, sflint::FmpzConstRef(p)));
    assert(radical3.multiplier_ring(ring3));
    assert(order_index(sflint::FmpzRef(index), sqrt5_order, ring3));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(order_index(sflint::FmpzRef(index), ring3, sqrt5_order));
    assert(sflint::fmpz_equal_si(index, 1));
    return 0;
}

int test_multiplier_ring_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta_element(field);
    assert(theta_element.gen());
    silex::OrderElement theta(order);
    assert(theta.set_element(theta_element));

    silex::Ideal principal(order);
    assert(principal.set_principal(theta));

    silex::Order ring(field);
    sflint::Fmpz index;
    sflint::Fmpz discriminant;
    assert(principal.multiplier_ring(ring));
    assert(order_index(sflint::FmpzRef(index), order, ring));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(order_index(sflint::FmpzRef(index), ring, order));
    assert(sflint::fmpz_equal_si(index, 1));
    assert(ring.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, 5));
    assert(!ring.maximality_known());
    return 0;
}

int test_p_radical_explicit_order_cases() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Ideal radical(order);
    silex::Ideal expected(order);
    silex::Element rational(field);
    silex::OrderElement generator(order);
    sflint::Fmpz p;
    sflint::Fmpz norm;

    sflint::fmpz_set_ui(p, 2);
    assert(order.p_radical(radical, sflint::FmpzConstRef(p)));
    element_from_coeff(rational, 0, 2, 1);
    assert(generator.set_element(rational));
    assert(expected.set_principal(generator));
    assert(radical.equal(expected));
    assert(radical.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 4));

    sflint::fmpz_set_ui(p, 5);
    assert(order.p_radical(radical, sflint::FmpzConstRef(p)));
    silex::Element theta(field);
    assert(theta.gen());
    assert(generator.set_element(theta));
    assert(expected.set_principal(generator));
    assert(radical.equal(expected));
    assert(radical.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 5));

    sflint::fmpz_set_ui(p, 3);
    assert(order.p_radical(radical, sflint::FmpzConstRef(p)));
    element_from_coeff(rational, 0, 3, 1);
    assert(generator.set_element(rational));
    assert(expected.set_principal(generator));
    assert(radical.equal(expected));
    assert(radical.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 9));

    return 0;
}

int test_fractional_arithmetic_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::FractionalIdeal half = fractional_from_hnf(order, 1, 2);
    silex::FractionalIdeal third = fractional_from_hnf(order, 1, 3);
    silex::FractionalIdeal out(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq qnorm;
    sflint::Fmpz exponent;

    assert(out.add(half, third));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 6));
    expected = fractional_from_hnf(order, 1, 6);
    assert(out.equal(expected));

    assert(out.intersect(half, third));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 1));
    assert(expected.one());
    assert(out.equal(expected));

    assert(out.multiply(half, third));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 6));
    expected = fractional_from_hnf(order, 1, 6);
    assert(out.equal(expected));

    assert(out.colon(half, third));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 3, 2));
    expected = fractional_from_hnf(order, 3, 2);
    assert(out.equal(expected));

    assert(out.colon(half, half));
    assert(out.one());

    assert(expected.one());
    assert(out.colon(half, expected));
    assert(out.equal(half));

    assert(out.colon(expected, half));
    expected = fractional_from_hnf(order, 2, 1);
    assert(out.equal(expected));

    silex::FractionalIdeal alias_colon = fractional_from_hnf(order, 1, 2);
    assert(expected.one());
    assert(alias_colon.colon(expected, alias_colon));
    expected = fractional_from_hnf(order, 2, 1);
    assert(alias_colon.equal(expected));

    assert(out.invert(half));
    expected = fractional_from_hnf(order, 2, 1);
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 2, 1));

    silex::FractionalIdeal three = fractional_from_hnf(order, 3, 1);
    silex::FractionalIdeal one_third = fractional_from_hnf(order, 1, 3);
    assert(out.invert(three));
    assert(out.equal(one_third));
    silex::FractionalIdeal product(order);
    assert(product.multiply(three, out));
    assert(product.equal(fractional_from_hnf(order, 1, 1)));

    assert(three.invert(three));
    assert(three.equal(one_third));

    sflint::fmpz_zero(sflint::FmpzRef(exponent));
    assert(out.pow_fmpz(expected, sflint::FmpzConstRef(exponent)));
    assert(out.one());

    sflint::fmpz_set_ui(exponent, 1);
    assert(out.pow_fmpz(expected, sflint::FmpzConstRef(exponent)));
    assert(out.equal(expected));

    sflint::fmpz_set_ui(exponent, 2);
    assert(out.pow_fmpz(expected, sflint::FmpzConstRef(exponent)));
    expected = fractional_from_hnf(order, 4, 1);
    assert(out.equal(expected));

    sflint::fmpz_set_ui(exponent, 3);
    assert(out.pow_fmpz(half, sflint::FmpzConstRef(exponent)));
    expected = fractional_from_hnf(order, 1, 8);
    assert(out.equal(expected));

    sflint::fmpz_set_si(exponent, -1);
    assert(out.pow_fmpz(half, sflint::FmpzConstRef(exponent)));
    expected = fractional_from_hnf(order, 2, 1);
    assert(out.equal(expected));

    sflint::fmpz_set_si(exponent, -2);
    assert(expected.pow_fmpz(expected, sflint::FmpzConstRef(exponent)));
    assert(expected.equal(fractional_from_hnf(order, 1, 4)));

    silex::FractionalIdeal alias_add = fractional_from_hnf(order, 1, 2);
    assert(alias_add.add(alias_add, third));
    expected = fractional_from_hnf(order, 1, 6);
    assert(alias_add.equal(expected));

    silex::FractionalIdeal alias_intersect = fractional_from_hnf(order, 1, 3);
    silex::FractionalIdeal unit(order);
    assert(unit.one());
    assert(alias_intersect.intersect(alias_intersect, unit));
    assert(alias_intersect.equal(unit));

    assert(half.multiply(half, third));
    expected = fractional_from_hnf(order, 1, 6);
    assert(half.equal(expected));
    return 0;
}

int test_fractional_arithmetic_nongalois_cubic() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_over_two(field);
    silex::Element theta_over_three(field);
    silex::Element theta_over_six(field);
    silex::Element theta_squared_over_six(field);
    silex::Element three_over_two(field);
    element_from_coeff(theta_over_two, 1, 1, 2);
    element_from_coeff(theta_over_three, 1, 1, 3);
    element_from_coeff(theta_over_six, 1, 1, 6);
    element_from_coeff(theta_squared_over_six, 2, 1, 6);
    element_from_coeff(three_over_two, 0, 3, 2);

    silex::Element theta(field);
    assert(theta.gen());

    silex::FractionalIdeal left(order);
    silex::FractionalIdeal right(order);
    silex::FractionalIdeal out(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq qnorm;

    assert(left.set_principal(theta_over_two));
    assert(right.set_principal(theta_over_three));

    assert(out.add(left, right));
    assert(expected.set_principal(theta_over_six));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 108));

    assert(out.intersect(left, right));
    assert(expected.set_principal(theta));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 2, 1));

    assert(out.multiply(left, right));
    assert(expected.set_principal(theta_squared_over_six));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 54));

    assert(out.set_principal(three_over_two));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 27, 8));

    return 0;
}

int test_fractional_arithmetic_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta(field);
    silex::Element theta_over_two(field);
    silex::Element theta_over_three(field);
    silex::Element theta_over_six(field);
    silex::Element one_over_three(field);
    assert(theta.gen());
    element_from_coeff(theta_over_two, 1, 1, 2);
    element_from_coeff(theta_over_three, 1, 1, 3);
    element_from_coeff(theta_over_six, 1, 1, 6);
    element_from_coeff(one_over_three, 0, 1, 3);

    silex::FractionalIdeal left(order);
    silex::FractionalIdeal right(order);
    silex::FractionalIdeal out(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq qnorm;

    assert(left.set_principal(theta_over_two));
    assert(right.set_principal(theta_over_three));

    assert(out.add(left, right));
    assert(expected.set_principal(theta_over_six));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 18));

    assert(out.intersect(left, right));
    assert(expected.set_principal(theta));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 2, 1));

    assert(out.multiply(left, right));
    assert(expected.set_principal(one_over_three));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 9));
    return 0;
}

int test_fractional_colon_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta(field);
    silex::Element two(field);
    silex::Element theta_over_two(field);
    assert(theta.gen());
    element_from_coeff(two, 0, 2, 1);
    element_from_coeff(theta_over_two, 1, 1, 2);

    silex::FractionalIdeal left(order);
    silex::FractionalIdeal right(order);
    silex::FractionalIdeal out(order);
    silex::FractionalIdeal product(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq qnorm;

    assert(left.set_principal(theta));
    assert(right.set_principal(two));

    assert(out.colon(left, right));
    assert(expected.set_principal(theta_over_two));
    assert(out.equal(expected));
    assert(product.multiply(out, right));
    assert(product.equal(left));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));

    assert(left.one());
    assert(out.colon(left, expected));
    assert(expected.set_principal(theta));
    assert(out.equal(expected));
    return 0;
}

int test_fractional_colon_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta_over_two(field);
    silex::Element theta_over_three(field);
    silex::Element three_over_two(field);
    element_from_coeff(theta_over_two, 1, 1, 2);
    element_from_coeff(theta_over_three, 1, 1, 3);
    element_from_coeff(three_over_two, 0, 3, 2);

    silex::FractionalIdeal left(order);
    silex::FractionalIdeal right(order);
    silex::FractionalIdeal out(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq qnorm;

    assert(left.set_principal(theta_over_two));
    assert(right.set_principal(theta_over_three));

    assert(out.colon(left, right));
    assert(expected.set_principal(three_over_two));
    assert(out.equal(expected));
    assert(out.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 9, 4));

    return 0;
}

int test_fractional_principal_pseudorandom_identities() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element a(field);
    silex::Element b(field);
    silex::Element ab(field);
    silex::FractionalIdeal left(order);
    silex::FractionalIdeal right(order);
    silex::FractionalIdeal product(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpq norm_left;
    sflint::Fmpq norm_right;
    sflint::Fmpq norm_product;

    for (slong t = 1; t <= 10; ++t) {
        element_set_fractional_pseudorandom(a, t, 2 + (t % 3));
        element_set_fractional_pseudorandom(b, 19 + 5 * t, 3 + (t % 2));
        assert(ab.multiply(a, b));

        assert(left.set_principal(a));
        assert(right.set_principal(b));
        assert(expected.set_principal(ab));

        assert(product.multiply(left, right));
        assert(product.equal(expected));

        assert(left.norm(sflint::FmpqRef(norm_left)));
        assert(right.norm(sflint::FmpqRef(norm_right)));
        assert(product.norm(sflint::FmpqRef(norm_product)));
        fmpq_mul(norm_left.raw(), norm_left.raw(), norm_right.raw());
        assert(sflint::fmpq_equal(norm_product, norm_left));
    }

    return 0;
}

int test_fractional_principal_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_over_two(field);
    element_from_coeff(theta_over_two, 1, 1, 2);

    silex::FractionalIdeal principal(order);
    silex::FractionalIdeal copy(order);
    sflint::Fmpq qnorm;

    assert(principal.set_principal(theta_over_two));
    assert(principal.contains(theta_over_two));
    assert(principal.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));

    assert(copy.set(principal));
    assert(copy.equal(principal));
    assert(copy.one());
    copy.swap(principal);
    assert(copy.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));
    assert(principal.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 1));
    assert(copy.set(copy));
    assert(copy.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));

    return 0;
}

int test_fractional_principal_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta_over_two(field);
    element_from_coeff(theta_over_two, 1, 1, 2);

    silex::FractionalIdeal principal(order);
    sflint::Fmpq qnorm;

    assert(principal.set_principal(theta_over_two));
    assert(principal.contains(theta_over_two));
    assert(principal.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 5, 4));

    return 0;
}

int test_quadratic_principal_and_validation() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderElement theta(order);
    silex::Element theta_element(field);
    assert(theta_element.gen());
    assert(theta.set_element(theta_element));

    silex::Ideal theta_ideal(order);
    sflint::Fmpz norm;
    assert(theta_ideal.set_principal(theta));
    assert(theta_ideal.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 2));
    assert(theta_ideal.contains(theta));

    silex::OrderElement one(order);
    assert(one.one());
    assert(!theta_ideal.contains(one));

    silex::Ideal out(order);
    assert(out.one());
    sflint::FmpzMat invalid(2, 2);
    mat_entry_si(invalid, 0, 0, 1);
    mat_entry_si(invalid, 1, 1, 2);
    assert(!out.set_hnf(sflint::FmpzMatConstRef(invalid)));
    assert(out.is_one());

    sflint::FmpzMat wrong_shape(1, 2);
    assert(!out.set_hnf(sflint::FmpzMatConstRef(wrong_shape)));
    assert(out.is_one());

    return 0;
}

int test_set_hnf_validation() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    sflint::FmpzMat hnf(2, 2);
    sflint::FmpzMat got(2, 2);
    sflint::Fmpz norm;

    mat_entry_si(hnf, 0, 0, 2);
    mat_entry_si(hnf, 1, 1, 3);
    assert(!ideal.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(!ideal.has_hnf());
    assert(!ideal.get_hnf(sflint::FmpzMatRef(got)));
    assert(!ideal.norm(sflint::FmpzRef(norm)));

    mat_entry_si(hnf, 1, 1, 0);
    assert(!ideal.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(!ideal.has_hnf());
    assert(!ideal.get_hnf(sflint::FmpzMatRef(got)));
    assert(!ideal.norm(sflint::FmpzRef(norm)));

    mat_entry_si(hnf, 0, 0, 2);
    mat_entry_si(hnf, 0, 1, 0);
    mat_entry_si(hnf, 1, 0, 0);
    mat_entry_si(hnf, 1, 1, 2);
    assert(ideal.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(ideal.get_hnf(sflint::FmpzMatRef(got)));
    assert(mat_is_diag2(got, 2));
    assert(ideal.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 4));

    assert(ideal.set_hnf(sflint::FmpzMatConstRef(got)));
    assert(ideal.get_hnf(sflint::FmpzMatRef(got)));
    assert(mat_is_diag2(got, 2));
    assert(ideal.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 4));

    mat_entry_si(hnf, 0, 0, 2);
    mat_entry_si(hnf, 0, 1, 0);
    mat_entry_si(hnf, 1, 0, 0);
    mat_entry_si(hnf, 1, 1, 3);
    assert(!ideal.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(ideal.get_hnf(sflint::FmpzMatRef(got)));
    assert(mat_is_diag2(got, 2));
    assert(ideal.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 4));

    return 0;
}

int test_integral_noninvertible_norm_trap() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -4);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    silex::Ideal product(order);
    sflint::FmpzMat hnf(2, 2);
    sflint::FmpzMat got(2, 2);
    sflint::Fmpz norm;

    mat_entry_si(hnf, 0, 0, 2);
    mat_entry_si(hnf, 1, 1, 1);
    assert(ideal.set_hnf(sflint::FmpzMatConstRef(hnf)));
    assert(ideal.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 2));

    assert(product.multiply(ideal, ideal));
    assert(product.get_hnf(sflint::FmpzMatRef(got)));
    assert(mat_entry_is_si(got, 0, 0, 4));
    assert(mat_entry_is_si(got, 0, 1, 0));
    assert(mat_entry_is_si(got, 1, 0, 0));
    assert(mat_entry_is_si(got, 1, 1, 2));
    assert(product.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 8));
    return 0;
}

int test_fractional_inverse_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Element theta_element(field);
    assert(theta_element.gen());
    silex::OrderElement theta(order);
    assert(theta.set_element(theta_element));

    silex::Ideal theta_num(order);
    assert(theta_num.set_principal(theta));

    silex::FractionalIdeal theta_over_two(order);
    silex::FractionalIdeal inverse(order);
    silex::FractionalIdeal product(order);
    silex::FractionalIdeal expected(order);
    sflint::Fmpz den;

    sflint::fmpz_set_ui(den, 2);
    assert(theta_over_two.set_integral_den(theta_num,
                                           sflint::FmpzConstRef(den)));
    assert(inverse.invert(theta_over_two));
    assert(expected.set_integral(theta_num));
    assert(inverse.equal(expected));

    assert(product.multiply(theta_over_two, inverse));
    assert(expected.one());
    assert(product.equal(expected));

    assert(inverse.invert(inverse));
    assert(inverse.equal(theta_over_two));
    return 0;
}

int test_fractional_inverse_nontrivial_order_basis() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::Element theta_element(field);
    assert(theta_element.gen());
    silex::Element two_theta_element(field);
    assert(two_theta_element.add(theta_element, theta_element));

    silex::OrderElement theta(order);
    silex::OrderElement two_theta(order);
    assert(theta.set_element(theta_element));
    assert(two_theta.set_element(two_theta_element));

    silex::Ideal theta_num(order);
    silex::Ideal two_theta_num(order);
    assert(theta_num.set_principal(theta));
    assert(two_theta_num.set_principal(two_theta));

    silex::FractionalIdeal theta_over_two(order);
    silex::FractionalIdeal expected(order);
    silex::FractionalIdeal inverse(order);
    silex::FractionalIdeal product(order);
    sflint::Fmpz den;

    sflint::fmpz_set_ui(den, 2);
    assert(theta_over_two.set_integral_den(theta_num,
                                           sflint::FmpzConstRef(den)));
    sflint::fmpz_set_ui(den, 5);
    assert(expected.set_integral_den(two_theta_num,
                                     sflint::FmpzConstRef(den)));

    assert(inverse.invert(theta_over_two));
    assert(inverse.equal(expected));
    assert(product.multiply(theta_over_two, inverse));
    assert(expected.one());
    assert(product.equal(expected));
    return 0;
}

int test_fractional_inverse_maximal_quadratic_matches_colon() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order equation = silex::test::equation_order(field);

    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::Element theta_element(field);
    assert(theta_element.gen());
    silex::OrderElement theta(order);
    assert(theta.set_element(theta_element));

    silex::Ideal theta_num(order);
    assert(theta_num.set_principal(theta));

    silex::FractionalIdeal theta_over_two(order);
    sflint::Fmpz den;
    sflint::fmpz_set_ui(den, 2);
    assert(theta_over_two.set_integral_den(theta_num,
                                           sflint::FmpzConstRef(den)));

    silex::FractionalIdeal inverse(order);
    silex::FractionalIdeal expected(order);
    silex::FractionalIdeal one(order);
    silex::FractionalIdeal product(order);
    assert(one.one());

    assert(inverse.invert(theta_over_two));
    assert(expected.colon(one, theta_over_two));
    assert(inverse.equal(expected));

    assert(product.multiply(theta_over_two, inverse));
    assert(product.equal(one));
    return 0;
}

int test_fractional_integral_denominator() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal unit(order);
    silex::Ideal two(order);
    sflint::FmpzMat hnf(1, 1);
    assert(unit.one());
    mat_entry_si(hnf, 0, 0, 2);
    assert(two.set_hnf(sflint::FmpzMatConstRef(hnf)));

    silex::FractionalIdeal half(order);
    sflint::Fmpz den;
    sflint::fmpz_set_ui(den, 4);
    assert(half.set_integral_den(two, sflint::FmpzConstRef(den)));

    silex::Ideal numerator(order);
    assert(half.get_integral_den(numerator, sflint::FmpzRef(den)));
    assert(numerator.is_one());
    assert(sflint::fmpz_equal_si(den, 2));

    silex::FractionalIdeal direct_half(order);
    sflint::fmpz_set_ui(den, 2);
    assert(direct_half.set_integral_den(unit, sflint::FmpzConstRef(den)));
    assert(half.equal(direct_half));

    sflint::Fmpq qnorm;
    assert(half.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));
    const silex::FractionalIdeal& const_half = half;
    assert(const_half.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 1, 2));

    silex::OrderElement six_generator(order);
    assert(six_generator.set_si(6));
    silex::Ideal six(order);
    assert(six.set_principal(six_generator));

    silex::FractionalIdeal three(order);
    sflint::fmpz_set_ui(den, 2);
    assert(three.set_integral_den(six, sflint::FmpzConstRef(den)));
    assert(three.get_integral_den(numerator, sflint::FmpzRef(den)));
    assert(numerator.get_hnf(sflint::FmpzMatRef(hnf)));
    assert(mat_entry_is_si(hnf, 0, 0, 3));
    assert(sflint::fmpz_equal_si(den, 1));
    assert(three.norm(sflint::FmpqRef(qnorm)));
    assert(fmpq_is_si(qnorm, 3, 1));

    silex::OrderElement three_generator(order);
    assert(three_generator.set_si(3));
    silex::Ideal three_integral(order);
    assert(three_integral.set_principal(three_generator));
    silex::FractionalIdeal integral_three(order);
    assert(integral_three.set_integral(three_integral));
    assert(three.equal(integral_three));

    silex::FractionalIdeal known_three(order);
    assert(silex::detail::set_integral_ideal_known_hnf(known_three,
                                                       three_integral));
    assert(known_three.equal(integral_three));
    assert(known_three.get_integral_den(numerator, sflint::FmpzRef(den)));
    assert(numerator.get_hnf(sflint::FmpzMatRef(hnf)));
    assert(mat_entry_is_si(hnf, 0, 0, 3));
    assert(sflint::fmpz_equal_si(den, 1));

    silex::FractionalIdeal known_before(order);
    silex::Ideal undefined_ideal;
    assert(known_before.set(known_three));
    assert(!silex::detail::set_integral_ideal_known_hnf(known_before,
                                                        undefined_ideal));
    assert(known_before.equal(known_three));

    silex::Element rational(field);
    element_from_coeff(rational, 0, 1, 2);
    assert(half.contains(rational));
    element_from_coeff(rational, 0, 1, 1);
    assert(half.contains(rational));
    element_from_coeff(rational, 0, 1, 3);
    assert(!half.contains(rational));

    silex::FractionalIdeal before(order);
    assert(before.set(three));
    sflint::fmpz_zero(sflint::FmpzRef(den));
    assert(!three.set_integral_den(three_integral, sflint::FmpzConstRef(den)));
    assert(three.equal(before));

    sflint::fmpz_set_si(den, -1);
    assert(!three.set_integral_den(three_integral, sflint::FmpzConstRef(den)));
    assert(three.equal(before));

    return 0;
}

int test_fractional_noninvertible_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -4);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal noninvertible_num(order);
    sflint::FmpzMat hnf(2, 2);
    mat_entry_si(hnf, 0, 0, 2);
    mat_entry_si(hnf, 1, 1, 1);
    assert(noninvertible_num.set_hnf(sflint::FmpzMatConstRef(hnf)));

    silex::FractionalIdeal noninvertible(order);
    silex::FractionalIdeal one(order);
    silex::FractionalIdeal colon(order);
    silex::FractionalIdeal product(order);
    silex::FractionalIdeal before(order);
    sflint::Fmpz exponent;
    assert(noninvertible.set_integral(noninvertible_num));
    assert(one.one());

    assert(colon.colon(one, noninvertible));
    assert(product.multiply(noninvertible, colon));
    assert(!product.equal(one));

    assert(before.set(colon));
    assert(!colon.invert(noninvertible));
    assert(colon.equal(before));

    assert(before.set(noninvertible));
    assert(!noninvertible.invert(noninvertible));
    assert(noninvertible.equal(before));

    assert(before.set(one));
    sflint::fmpz_set_si(exponent, -1);
    assert(!one.pow_fmpz(noninvertible, sflint::FmpzConstRef(exponent)));
    assert(one.equal(before));
    return 0;
}

int test_move_clear_and_redefine() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal source(order);
    assert(source.one());
    silex::Ideal moved(std::move(source));
    assert(moved.is_defined());
    assert(same_parent(moved.parent(), order));
    assert(moved.degree() == 1);
    assert(moved.is_one());
    assert(!source.is_defined());
    assert(source.parent() == nullptr);
    assert(source.degree() == 0);

    silex::Ideal assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(same_parent(assigned.parent(), order));
    assert(assigned.is_one());
    assert(!moved.is_defined());
    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(assigned.degree() == 0);
    assert(!assigned.one());
    assert(assigned.define(order));
    assert(assigned.one());

    silex::FractionalIdeal frac(order);
    assert(frac.one());
    silex::FractionalIdeal frac_moved(std::move(frac));
    assert(frac_moved.is_defined());
    assert(same_parent(frac_moved.parent(), order));
    assert(frac_moved.degree() == 1);
    assert(frac_moved.has_integral_denominator());
    assert(!frac.is_defined());
    assert(frac.parent() == nullptr);
    assert(frac.degree() == 0);

    silex::FractionalIdeal frac_assigned;
    frac_assigned = std::move(frac_moved);
    assert(frac_assigned.is_defined());
    assert(same_parent(frac_assigned.parent(), order));
    assert(frac_assigned.has_integral_denominator());
    assert(!frac_moved.is_defined());
    frac_assigned.clear();
    assert(!frac_assigned.is_defined());
    assert(frac_assigned.parent() == nullptr);
    assert(frac_assigned.degree() == 0);
    assert(!frac_assigned.one());
    assert(frac_assigned.define(order));
    assert(frac_assigned.one());

    return 0;
}

int test_failure_preserves_outputs() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::FmpqPoly other_polynomial;
    poly_x2_minus(other_polynomial, 3);
    silex::NumberField other_field;
    other_field = field_by_polynomial(other_polynomial);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);

    silex::Ideal out(order);
    silex::Ideal left(order);
    silex::Ideal other(other_order);
    assert(out.one());
    assert(left.one());
    assert(other.one());
    assert(!out.add(left, other));
    assert(out.is_one());
    assert(!out.intersect(left, other));
    assert(out.is_one());
    assert(!out.multiply(left, other));
    assert(out.is_one());

    silex::FractionalIdeal frac(order);
    assert(frac.one());
    sflint::Fmpz den;
    sflint::fmpz_zero(den);
    assert(!frac.set_integral_den(left, sflint::FmpzConstRef(den)));
    assert(frac.has_integral_denominator());
    assert(frac.get_integral_den(left, sflint::FmpzRef(den)));
    assert(sflint::fmpz_is_one(den));

    sflint::fmpz_set_si(den, -1);
    assert(!frac.set_integral_den(left, sflint::FmpzConstRef(den)));
    assert(frac.has_integral_denominator());
    assert(frac.get_integral_den(left, sflint::FmpzRef(den)));
    assert(sflint::fmpz_is_one(den));

    silex::Element zero(field);
    assert(zero.zero());
    assert(!frac.set_principal(zero));
    assert(frac.get_integral_den(left, sflint::FmpzRef(den)));
    assert(sflint::fmpz_is_one(den));

    silex::Element other_element(other_field);
    assert(other_element.one());
    assert(!frac.set_principal(other_element));
    assert(frac.get_integral_den(left, sflint::FmpzRef(den)));
    assert(sflint::fmpz_is_one(den));

    sflint::Fmpz before;
    sflint::Fmpz after;
    assert(order.discriminant(sflint::FmpzRef(before)));
    assert(!left.multiplier_ring(order));
    assert(order.discriminant(sflint::FmpzRef(after)));
    assert(sflint::fmpz_equal(before, after));

    assert(other_order.discriminant(sflint::FmpzRef(before)));
    assert(!left.multiplier_ring(other_order));
    assert(other_order.discriminant(sflint::FmpzRef(after)));
    assert(sflint::fmpz_equal(before, after));

    assert(left.is_one());
    sflint::Fmpz p;
    sflint::fmpz_set_ui(p, 2);
    assert(!order.p_radical(other, sflint::FmpzConstRef(p)));
    assert(other.is_one());

    sflint::fmpz_set_ui(p, 1);
    assert(!order.p_radical(left, sflint::FmpzConstRef(p)));
    assert(left.is_one());

    return 0;
}

silex::Ideal local_integral_ideal() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);
    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(ideal.one());
    return ideal;
}

silex::FractionalIdeal local_fractional_ideal() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);
    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::FractionalIdeal ideal(order);
    assert(ideal.one());
    return ideal;
}

int test_keeps_parent_order_alive() {
    silex::Ideal ideal = local_integral_ideal();
    assert(ideal.is_defined());
    assert(ideal.parent() != nullptr);
    assert(ideal.parent()->parent() != nullptr);
    assert(ideal.degree() == 1);
    assert(ideal.is_one());

    silex::Ideal copy(*ideal.parent());
    assert(copy.is_defined());
    assert(copy.set(ideal));
    assert(copy.is_one());

    silex::FractionalIdeal fractional = local_fractional_ideal();
    assert(fractional.is_defined());
    assert(fractional.parent() != nullptr);
    assert(fractional.parent()->parent() != nullptr);
    assert(fractional.degree() == 1);
    assert(fractional.has_integral_denominator());

    silex::FractionalIdeal fractional_copy(*fractional.parent());
    assert(fractional_copy.is_defined());
    assert(fractional_copy.set(fractional));
    assert(fractional_copy.has_integral_denominator());
    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one_integral_arithmetic() == 0);
    assert(test_integral_arithmetic_nongalois_cubic() == 0);
    assert(test_integral_cubic_swap_and_self_set() == 0);
    assert(test_integral_arithmetic_quadratic() == 0);
    assert(test_integral_arithmetic_nontrivial_order_basis() == 0);
    assert(test_integral_principal_pseudorandom_identities() == 0);
    assert(test_principal_times_integral_direct_matches_fractional_product() == 0);
    assert(test_two_generator_product_matches_generic() == 0);
    assert(test_principal_times_integral_direct_nontrivial_order_basis() == 0);
    assert(test_principal_times_integral_direct_failure_preserves_output() == 0);
    assert(test_integral_coprime_add_to_one_degree_one() == 0);
    assert(test_integral_add_to_one_swapped_quadratic_basis() == 0);
    assert(test_multiplier_ring_degree_one() == 0);
    assert(test_p_radical_degree_one() == 0);
    assert(test_multiplier_ring_quadratic_p_radical_shape() == 0);
    assert(test_p_radical_stable_multiplier_rings() == 0);
    assert(test_multiplier_ring_nontrivial_order_basis() == 0);
    assert(test_p_radical_explicit_order_cases() == 0);
    assert(test_fractional_arithmetic_degree_one() == 0);
    assert(test_fractional_arithmetic_nongalois_cubic() == 0);
    assert(test_fractional_arithmetic_quadratic() == 0);
    assert(test_fractional_colon_quadratic() == 0);
    assert(test_fractional_colon_nontrivial_order_basis() == 0);
    assert(test_fractional_principal_pseudorandom_identities() == 0);
    assert(test_fractional_principal_quadratic() == 0);
    assert(test_fractional_principal_nontrivial_order_basis() == 0);
    assert(test_quadratic_principal_and_validation() == 0);
    assert(test_set_hnf_validation() == 0);
    assert(test_integral_noninvertible_norm_trap() == 0);
    assert(test_fractional_inverse_quadratic() == 0);
    assert(test_fractional_inverse_nontrivial_order_basis() == 0);
    assert(test_fractional_inverse_maximal_quadratic_matches_colon() == 0);
    assert(test_fractional_integral_denominator() == 0);
    assert(test_fractional_noninvertible_preserves_output() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_failure_preserves_outputs() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    return 0;
}
