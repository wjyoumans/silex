#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>

#include "test_support.hpp"
#include "ideal_product_fixtures.hpp"
#include "../src/ideal/hnf_membership_internal.hpp"
#include "../src/ideal/ideal_access_internal.hpp"
#include "../src/ideal/two_generator_internal.hpp"
#include "../src/order/order_internal.hpp"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <string_view>
#include <utility>
#include <vector>

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

// Checks detail::set_principal_ideal_known_norm against Ideal::set_principal
// for pseudorandom integral generators of `order`, with the true norm (the
// determinant-modular route), and with a norm that is a proper multiple, a
// noninteger or zero (the full-HNF route).  Both must store the same HNF.
void check_principal_known_norm_matches(silex::Order& order,
                                        slong degree,
                                        slong seeds) noexcept {
    const silex::NumberField* field = order.parent();
    assert(field != nullptr);
    silex::Element alpha(*field);
    silex::OrderElement generator(order);
    silex::Ideal expected(order);
    silex::Ideal out(order);
    sflint::FmpzMat expected_hnf(degree, degree);
    sflint::FmpzMat out_hnf(degree, degree);
    sflint::Fmpq norm;
    sflint::Fmpq wrong_norm;

    for (slong seed = 1; seed <= seeds; ++seed) {
        sflint::FmpqPoly polynomial;
        for (slong j = 0; j < degree; ++j) {
            const slong c = ((seed * (7 + 3 * j) + j * j + 2) % 9) - 4;
            sflint::fmpq_poly_set_coeff_si(polynomial, j, c);
        }
        if (fmpq_poly_is_zero(polynomial.raw()) != 0) {
            sflint::fmpq_poly_set_coeff_si(polynomial, 0, seed + 1);
        }
        assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        if (!generator.set_element(alpha)) {
            continue;
        }
        assert(alpha.norm(sflint::FmpqRef(norm)));
        if (fmpq_is_zero(norm.raw()) != 0) {
            continue;
        }
        assert(expected.set_principal(generator));
        assert(expected.get_hnf(sflint::FmpzMatRef(expected_hnf)));

        const auto check = [&](const sflint::Fmpq& value) noexcept {
            silex::Ideal fresh(order);
            assert(silex::detail::set_principal_ideal_known_norm(
                    fresh, generator, sflint::FmpqConstRef(value)));
            assert(fresh.get_hnf(sflint::FmpzMatRef(out_hnf)));
            assert(fmpz_mat_equal(out_hnf.raw(), expected_hnf.raw()) != 0);
            assert(fresh.equal(expected));
            const silex::OrderElement* known =
                    fresh.known_principal_generator();
            assert(known != nullptr && known->equal(generator));
        };
        check(norm);
        fmpq_neg(wrong_norm.raw(), norm.raw());
        check(wrong_norm);
        fmpq_mul_si(wrong_norm.raw(), norm.raw(), 3);
        check(wrong_norm);
        fmpz_mul_ui(fmpq_denref(wrong_norm.raw()),
                    fmpq_denref(wrong_norm.raw()), 2);
        fmpq_canonicalise(wrong_norm.raw());
        check(wrong_norm);
        fmpq_zero(wrong_norm.raw());
        check(wrong_norm);
    }

    // A zero generator fails and leaves the output unchanged.
    assert(generator.set_si(0));
    assert(out.set_principal(generator) == false);
    assert(generator.set_si(2));
    assert(out.set_principal(generator));
    assert(out.get_hnf(sflint::FmpzMatRef(expected_hnf)));
    assert(generator.set_si(0));
    sflint::fmpq_set_si(norm, 1, 1);
    assert(!silex::detail::set_principal_ideal_known_norm(
            out, generator, sflint::FmpqConstRef(norm)));
    assert(out.get_hnf(sflint::FmpzMatRef(out_hnf)));
    assert(fmpz_mat_equal(out_hnf.raw(), expected_hnf.raw()) != 0);
}

int test_principal_known_norm_matches_set_principal() {
    {
        sflint::FmpqPoly polynomial;
        poly_x3_minus(polynomial, 2);
        silex::NumberField field;
        field = field_by_polynomial(polynomial);
        silex::Order order;
        order = silex::test::equation_order(field);
        check_principal_known_norm_matches(order, 3, 12);
    }
    {
        // Z[(1 + sqrt 5)/2]: an order basis other than the power basis.
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
        check_principal_known_norm_matches(order, 2, 12);
    }
    {
        // Z[3 sqrt 5]: a nonmaximal order.
        sflint::FmpqPoly polynomial;
        poly_x2_minus(polynomial, 5);
        silex::NumberField field;
        field = field_by_polynomial(polynomial);
        silex::Order order(field);
        sflint::FmpqMat basis(2, 2);
        qmat_entry_si(basis, 0, 0, 1);
        qmat_entry_si(basis, 1, 1, 3);
        assert(order.set_basis(sflint::FmpqMatConstRef(basis)));
        check_principal_known_norm_matches(order, 2, 12);
    }
    {
        // Degree 7: x^7 - 3x + 3, Eisenstein at 3.
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 7, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, 3);
        silex::NumberField field;
        field = field_by_polynomial(polynomial);
        silex::Order order;
        order = silex::test::equation_order(field);
        check_principal_known_norm_matches(order, 7, 12);
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

int test_two_generator_search_contracts() {
    silex::Order original;
    silex::Ideal original_left, original_right;
    assert(silex::test::ideal_product::fixture(original, original_left, original_right, 3, 3));
    sflint::FmpqMat basis(3, 3), transform(3, 3), changed(3, 3);
    assert(original.get_basis(sflint::FmpqMatRef(basis)));
    const slong entries[3][3] = {{7, 5, -21}, {-3, -2, 9}, {-4, -4, 13}};
    for (slong i = 0; i < 3; ++i)
        for (slong j = 0; j < 3; ++j)
            fmpq_set_si(fmpq_mat_entry(transform.raw(), i, j), entries[i][j], 1);
    sflint::Fmpq determinant;
    fmpq_mat_det(determinant.raw(), transform.raw());
    assert(fmpq_is_one(determinant.raw()));
    fmpq_mat_mul(changed.raw(), transform.raw(), basis.raw());
    silex::Order changed_order = silex::Order::from_basis(*original.parent(), sflint::FmpqMatConstRef(changed));
    assert(changed_order.is_defined());
    // Unimodular change of the already established maximal order; maximality
    // is recomputed, and the basis is kept.
    silex::Order order = silex::test::verified_maximal_order(changed_order);
    sflint::FmpqPoly polynomial;
    silex::test::ideal_product::polynomial(polynomial, 3);
    silex::Ideal left(order), right(order), rebuilt(order), expected(order), actual(order);
    assert(silex::test::ideal_product::evaluation_ideal(left, polynomial, 2));
    assert(silex::test::ideal_product::evaluation_ideal(right, polynomial, 3));
    silex::OrderElement beta(order);
    assert(beta.set_si(7));
    sflint::Fmpz minimum;
    fmpz_set_si(minimum.raw(), 99);
    silex::detail::TwoGeneratorSearchReport report;
    silex::detail::TwoGeneratorSearchOptions options;
    options.random_trials = 0;
    assert(!silex::detail::ideal_two_generator(sflint::FmpzRef(minimum), beta, left, options, &report));
    assert(report.basis_trials > 0 && report.random_trials == 0 && !report.verified);
    assert(beta.equal_si(7) && fmpz_equal_si(minimum.raw(), 99));
    sflint::FmpzMat preserved(3, 3), snapshot(3, 3);
    fmpz_mat_one(preserved.raw());
    fmpz_mat_set(snapshot.raw(), preserved.raw());
    assert(!silex::detail::try_ideal_product_hnf(sflint::FmpzMatRef(preserved), left, left, options));
    assert(fmpz_mat_equal(preserved.raw(), snapshot.raw()));
    options.random_trials = 32;
    assert(silex::detail::ideal_two_generator(sflint::FmpzRef(minimum), beta, left, options, &report));
    assert(report.verified && report.random_trials > 0 && report.random_trials <= 32);
    assert(fmpz_equal_si(minimum.raw(), 6));
    assert(silex::detail::set_known_two_generator_ideal(rebuilt, sflint::FmpzConstRef(minimum), beta));
    assert(rebuilt.equal(left));
    assert(silex::test::ideal_product::exhaustive_product(expected, left, right));
    assert(actual.multiply(left, right) && actual.equal(expected));
    assert(actual.known_principal_generator() == nullptr);
    assert(right.multiply(left, right) && right.equal(expected));

    // The parent check precedes search and publication.
    silex::OrderElement wrong_parent(original);
    assert(wrong_parent.set_si(9));
    assert(!silex::detail::ideal_two_generator(sflint::FmpzRef(minimum), wrong_parent, left));
    assert(wrong_parent.equal_si(9) && fmpz_equal_si(minimum.raw(), 6));
    options.random_trials = 33;
    assert(!silex::detail::ideal_two_generator(sflint::FmpzRef(minimum), beta, left, options));
    assert(rebuilt.equal(left));
    return 0;
}

int test_general_product_contents_and_fallback() {
    for (slong kind : {3, 5}) {
        silex::Order order;
        silex::Ideal left, right;
        assert(silex::test::ideal_product::fixture(order, left, right, 2, kind));
        silex::Ideal expected(order), scaled_left(order), scaled_right(order), result(order);
        assert(silex::test::ideal_product::exhaustive_product(expected, left, right));
        sflint::FmpzMat h(2, 2), k(2, 2), expected_hnf(2, 2);
        assert(left.get_hnf(sflint::FmpzMatRef(h)));
        assert(right.get_hnf(sflint::FmpzMatRef(k)));
        assert(expected.get_hnf(sflint::FmpzMatRef(expected_hnf)));
        sflint::Fmpz content;
        fmpz_one(content.raw());
        fmpz_mul_2exp(content.raw(), content.raw(), 80);
        fmpz_add_ui(content.raw(), content.raw(), 7);
        fmpz_mat_scalar_mul_fmpz(h.raw(), h.raw(), content.raw());
        fmpz_mat_scalar_mul_si(k.raw(), k.raw(), 3);
        assert(scaled_left.set_hnf(sflint::FmpzMatConstRef(h)));
        assert(scaled_right.set_hnf(sflint::FmpzMatConstRef(k)));
        fmpz_mul_ui(content.raw(), content.raw(), 3);
        fmpz_mat_scalar_mul_fmpz(expected_hnf.raw(), expected_hnf.raw(), content.raw());
        assert(expected.set_hnf(sflint::FmpzMatConstRef(expected_hnf)));
        assert(result.multiply(scaled_left, scaled_right) && result.equal(expected));
        assert(scaled_right.multiply(scaled_left, scaled_right) && scaled_right.equal(expected));

        silex::FractionalIdeal a(order), b(order), product(order), rational_expected(order);
        sflint::Fmpz denominator;
        fmpz_set_ui(denominator.raw(), 6);
        assert(a.set_integral_den(left, sflint::FmpzConstRef(denominator)));
        fmpz_set_ui(denominator.raw(), 10);
        assert(b.set_integral_den(right, sflint::FmpzConstRef(denominator)));
        assert(silex::test::ideal_product::exhaustive_product(expected, left, right));
        fmpz_set_ui(denominator.raw(), 60);
        assert(rational_expected.set_integral_den(expected, sflint::FmpzConstRef(denominator)));
        assert(product.multiply(a, b) && product.equal(rational_expected));
        assert(a.multiply(a, b) && a.equal(rational_expected));
    }
    return 0;
}

// set_known_two_generator_ideal takes the HNF of a O + beta O modulo a
// (fmpz_mat_hnf_modular_eldiv).  Compare it with the generic fmpz_mat_hnf of
// the 2n x n rows a I and the coordinates of beta w_i, over prime, composite,
// unit and multiword scalars, in orders with trivial and nontrivial index.
void check_known_two_generator_hnf(silex::Order& order,
                                   sflint::FmpzConstRef scalar,
                                   const sflint::FmpzMat& coordinates) noexcept {
    const slong n = order.degree();
    silex::OrderElement beta(order), basis_element(order), product(order);
    assert(beta.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
    sflint::FmpzMat rows(2 * n, n), reference(2 * n, n), unit(1, n), row(1, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i), scalar.raw());
        fmpz_mat_zero(unit.raw());
        fmpz_one(fmpz_mat_entry(unit.raw(), 0, i));
        assert(basis_element.set_coordinates(sflint::FmpzMatConstRef(unit)));
        assert(product.multiply(beta, basis_element));
        assert(product.get_coordinates(sflint::FmpzMatRef(row)));
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(rows.raw(), n + i, j),
                     fmpz_mat_entry(row.raw(), 0, j));
        }
    }
    fmpz_mat_hnf(reference.raw(), rows.raw());
    sflint::FmpzMat expected(n, n), actual(n, n);
    sflint::FmpzMatConstWindow top(reference, 0, 0, n, n);
    fmpz_mat_set(expected.raw(), top.raw());
    for (slong i = n; i < 2 * n; ++i) {
        for (slong j = 0; j < n; ++j) {
            assert(fmpz_is_zero(fmpz_mat_entry(reference.raw(), i, j)));
        }
    }

    silex::Ideal ideal(order);
    assert(silex::detail::set_known_two_generator_ideal(ideal, scalar, beta));
    assert(ideal.get_hnf(sflint::FmpzMatRef(actual)));
    if (!fmpz_mat_equal(actual.raw(), expected.raw())) {
        std::fprintf(stderr, "known two-generator HNF mismatch: n=%ld a=", n);
        fmpz_fprint(stderr, scalar.raw());
        std::fprintf(stderr, "\n");
        fmpz_mat_print_pretty(actual.raw());
        fmpz_mat_print_pretty(expected.raw());
        std::fflush(stdout);
    }
    assert(fmpz_mat_equal(actual.raw(), expected.raw()));
}

int test_known_two_generator_hnf_matches_generic() {
    sflint::FmpqPoly polynomial;
    std::vector<silex::Order> orders;

    // Q(sqrt 5) with basis 1, (1 + theta) / 2.
    poly_x2_minus(polynomial, 5);
    {
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order order(field);
        sflint::FmpqMat basis(2, 2);
        qmat_entry_si(basis, 0, 0, 1);
        qmat_entry_frac_si(basis, 1, 0, 1, 2);
        qmat_entry_frac_si(basis, 1, 1, 1, 2);
        assert(order.set_basis(sflint::FmpqMatConstRef(basis)));
        orders.push_back(std::move(order));
    }
    // x^3 - 2: the equation order.
    poly_x3_minus(polynomial, 2);
    orders.push_back(silex::test::equation_order(field_by_polynomial(polynomial)));
    // x^3 - 10: 10 = 1 mod 9, so the maximal order has index 3 over Z[theta].
    poly_x3_minus(polynomial, 10);
    {
        silex::Order equation =
                silex::test::equation_order(field_by_polynomial(polynomial));
        silex::Order order(*equation.parent());
        assert(order.maximal_order(equation));
        assert(!order.is_equation_order());
        orders.push_back(std::move(order));
    }
    // Phi_15, degree 8, and x^6 + x + 1 (equation orders).
    sflint::fmpq_poly_zero(polynomial);
    for (slong e : {8, 7, 5, 4, 3, 1, 0}) {
        sflint::fmpq_poly_set_coeff_si(polynomial, e, (e == 7 || e == 4 || e == 1) ? -1 : 1);
    }
    orders.push_back(silex::test::equation_order(field_by_polynomial(polynomial)));
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 6, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    orders.push_back(silex::test::equation_order(field_by_polynomial(polynomial)));

    std::vector<sflint::Fmpz> scalars;
    for (ulong a : {1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL, 11UL, 12UL, 13UL, 30UL,
                    97UL, 128UL, 4093UL * 4093UL,
                    18446744073709551557UL}) {  // 2^64 - 59, full-word nmod path
        scalars.emplace_back();
        fmpz_set_ui(scalars.back().raw(), a);
    }
    // Multiword scalars take FLINT's fmpz path: 2^64 + 13 and 6 (2^64 + 13).
    scalars.emplace_back();
    fmpz_one(scalars.back().raw());
    fmpz_mul_2exp(scalars.back().raw(), scalars.back().raw(), 64);
    fmpz_add_ui(scalars.back().raw(), scalars.back().raw(), 13);
    scalars.emplace_back();
    fmpz_mul_ui(scalars.back().raw(), scalars[scalars.size() - 2].raw(), 6);

    ulong state = 0x2545f491UL;
    for (silex::Order& order : orders) {
        const slong n = order.degree();
        for (const sflint::Fmpz& scalar : scalars) {
            for (slong trial = 0; trial < 5; ++trial) {
                sflint::FmpzMat coordinates(1, n);
                for (slong j = 0; j < n; ++j) {
                    state = state * 6364136223846793005UL + 1442695040888963407UL;
                    slong value = static_cast<slong>((state >> 33) % 41) - 20;
                    if (trial == 0) value = 0;                 // beta = 0
                    if (trial == 1) value = (j == 0) ? 1 : 0;  // beta = 1
                    if (trial == 4 && j == 0) value += 1000003;
                    fmpz_set_si(fmpz_mat_entry(coordinates.raw(), 0, j), value);
                }
                check_known_two_generator_hnf(order, sflint::FmpzConstRef(scalar),
                                              coordinates);
            }
        }
    }
    return 0;
}

// The generic branch of try_ideal_product_hnf takes the HNF of the rows
// [mk H; H M_beta] modulo mh mk (fmpz_mat_hnf_modular_eldiv).  Compare it with
// the plain fmpz_mat_hnf of all n^2 basis products, over maximal orders with
// trivial and nontrivial index, ideals with content, equal minima, squares and
// word and multiword moduli.  Outside a maximal order the branch is refused.
// The modulus needs only that h k is an ideal, not maximality: on a
// non-maximal order whose flag is forced through the test hook, every product
// the branch returns is still the plain product.
void plain_product_hnf(sflint::FmpzMat& out, const silex::Ideal& left,
                       const silex::Ideal& right) noexcept {
    const silex::Order& order = *left.parent();
    const slong n = order.degree();
    sflint::FmpzMat h(n, n), k(n, n), row(1, n), rows(n * n, n), reduced(n * n, n);
    assert(left.get_hnf(sflint::FmpzMatRef(h)));
    assert(right.get_hnf(sflint::FmpzMatRef(k)));
    silex::OrderElement a(order), b(order), ab(order);
    for (slong i = 0; i < n; ++i) {
        sflint::FmpzMatConstWindow hi(h, i, 0, i + 1, n);
        fmpz_mat_set(row.raw(), hi.raw());
        assert(a.set_coordinates(sflint::FmpzMatConstRef(row)));
        for (slong j = 0; j < n; ++j) {
            sflint::FmpzMatConstWindow kj(k, j, 0, j + 1, n);
            fmpz_mat_set(row.raw(), kj.raw());
            assert(b.set_coordinates(sflint::FmpzMatConstRef(row)));
            assert(ab.multiply(a, b) && ab.get_coordinates(sflint::FmpzMatRef(row)));
            for (slong c = 0; c < n; ++c) {
                fmpz_set(fmpz_mat_entry(rows.raw(), i * n + j, c),
                         fmpz_mat_entry(row.raw(), 0, c));
            }
        }
    }
    fmpz_mat_hnf(reduced.raw(), rows.raw());
    sflint::FmpzMatConstWindow top(reduced, 0, 0, n, n);
    fmpz_mat_set(out.raw(), top.raw());
}

int test_general_product_hnf_matches_plain_hnf() {
    sflint::FmpqPoly polynomial;
    slong generic = 0, multiword = 0;
    // x^2 - 5 (index 2), x^3 - 10 (index 3), x^3 - 2, the fixture quartic,
    // Phi_15 and x^6 + x + 1.
    for (slong which = 0; which < 6; ++which) {
        sflint::fmpq_poly_zero(polynomial);
        if (which == 0) {
            poly_x2_minus(polynomial, 5);
        } else if (which == 1) {
            poly_x3_minus(polynomial, 10);
        } else if (which == 2) {
            poly_x3_minus(polynomial, 2);
        } else if (which == 3) {
            silex::test::ideal_product::polynomial(polynomial, 4);
        } else if (which == 4) {
            for (slong e : {8, 7, 5, 4, 3, 1, 0}) {
                sflint::fmpq_poly_set_coeff_si(polynomial, e,
                        (e == 7 || e == 4 || e == 1) ? -1 : 1);
            }
        } else {
            sflint::fmpq_poly_set_coeff_si(polynomial, 6, 1);
            sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
            sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
        }
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order order(field);
        assert(order.maximal_order(equation) && order.is_maximal());
        const slong n = order.degree();

        // Evaluation ideals (f(k), theta - k); the last point makes f(k)
        // exceed one word, so mh mk takes FLINT's multiword path.
        std::vector<silex::Ideal> ideals;
        const slong large = (n == 2) ? (WORD(1) << 33) : (n == 3) ? WORD(3000000)
                          : (n == 4) ? (WORD(1) << 17) : WORD(3000);
        for (slong point : {2, 3, -2, 7, 11}) {
            ideals.emplace_back(order);
            assert(silex::test::ideal_product::evaluation_ideal(ideals.back(), polynomial, point));
        }
        ideals.emplace_back(order);
        assert(silex::test::ideal_product::evaluation_ideal(ideals.back(), polynomial, large));
        {
            sflint::Fmpz minimum;
            silex::OrderElement beta(order);
            assert(silex::detail::ideal_two_generator(sflint::FmpzRef(minimum), beta,
                                                      ideals.back()));
            assert(!fmpz_abs_fits_ui(minimum.raw()));
        }
        // A product without metadata, and an ideal with content 6.
        ideals.emplace_back(order);
        assert(ideals.back().multiply(ideals[0], ideals[1]));
        {
            sflint::FmpzMat h(n, n);
            assert(ideals[0].get_hnf(sflint::FmpzMatRef(h)));
            fmpz_mat_scalar_mul_si(h.raw(), h.raw(), 6);
            ideals.emplace_back(order);
            assert(ideals.back().set_hnf(sflint::FmpzMatConstRef(h)));
        }

        sflint::FmpzMat expected(n, n), actual(n, n);
        for (const silex::Ideal& left : ideals) {
            for (const silex::Ideal& right : ideals) {
                if (left.is_one() || right.is_one()) continue;
                plain_product_hnf(expected, left, right);
                silex::detail::TwoGeneratorSearchReport report;
                fmpz_mat_zero(actual.raw());
                if (!silex::detail::try_ideal_product_hnf(sflint::FmpzMatRef(actual), left, right,
                                                          {}, &report)) continue;
                if (report.verified) {
                    ++generic;
                    if (&left == &ideals[5] || &right == &ideals[5]) ++multiword;
                }
                if (!fmpz_mat_equal(actual.raw(), expected.raw())) {
                    std::fprintf(stderr, "general product HNF mismatch: n=%ld\n", n);
                    fmpz_mat_print_pretty(actual.raw());
                    fmpz_mat_print_pretty(expected.raw());
                    std::fflush(stdout);
                }
                assert(fmpz_mat_equal(actual.raw(), expected.raw()));
            }
        }

        // The equation order of x^3 - 10 is not maximal: the generic branch
        // is refused and the output is unchanged.
        if (which == 1) {
            assert(!equation.is_maximal());
            silex::Ideal a(equation), b(equation), product(equation);
            assert(silex::test::ideal_product::evaluation_ideal(a, polynomial, 2));
            assert(silex::test::ideal_product::evaluation_ideal(b, polynomial, 3));
            sflint::FmpzMat out(3, 3), snapshot(3, 3);
            fmpz_mat_one(out.raw());
            fmpz_mat_set(snapshot.raw(), out.raw());
            assert(!silex::detail::try_ideal_product_hnf(sflint::FmpzMatRef(out), a, b));
            assert(fmpz_mat_equal(out.raw(), snapshot.raw()));
            plain_product_hnf(expected, a, b);
            assert(product.multiply(a, b) && product.get_hnf(sflint::FmpzMatRef(actual)));
            assert(fmpz_mat_equal(actual.raw(), expected.raw()));

            // Forced flag: the branch runs on Z[10^(1/3)], index 3 in the
            // maximal order, and its products stay exact.
            std::vector<silex::Ideal> forced;
            for (slong point : {2, 3, -2, 7, 11}) {
                forced.emplace_back(equation);
                assert(silex::test::ideal_product::evaluation_ideal(forced.back(), polynomial, point));
            }
            // The evaluation ideals are principal, so invertible.  Add the
            // singular prime P = (3, theta - 1), which is not: P^2 =
            // <9, 3u, u^2> with u = theta - 1 has index 27, not N(P)^2 = 9.
            {
                silex::Element theta(field), shifted(field);
                silex::OrderElement u(equation);
                sflint::Fmpz three;
                fmpz_set_ui(three.raw(), 3);
                assert(theta.gen() && shifted.add_si(theta, -1) && u.set_element(shifted));
                forced.emplace_back(equation);
                assert(silex::detail::set_known_two_generator_ideal(
                        forced.back(), sflint::FmpzConstRef(three), u));
            }
            const silex::Ideal& singular = forced.back();
            silex::detail::OrderAccess::set_maximality_unchecked(equation, true);
            slong forced_generic = 0;
            bool singular_square = false;
            for (const silex::Ideal& left : forced) {
                for (const silex::Ideal& right : forced) {
                    if (left.is_one() || right.is_one()) continue;
                    silex::detail::TwoGeneratorSearchReport report;
                    if (!silex::detail::try_ideal_product_hnf(sflint::FmpzMatRef(actual), left,
                                                              right, {}, &report)) continue;
                    if (report.verified) ++forced_generic;
                    plain_product_hnf(expected, left, right);
                    assert(fmpz_mat_equal(actual.raw(), expected.raw()));
                    if (&left == &singular && &right == &singular && report.verified) {
                        sflint::Fmpz index;
                        fmpz_mat_det(index.raw(), expected.raw());
                        assert(fmpz_equal_ui(index.raw(), 27));
                        singular_square = true;
                    }
                }
            }
            silex::detail::OrderAccess::set_maximality_unchecked(equation, false);
            assert(forced_generic > 0 && singular_square);
        }
    }
    assert(generic >= 100 && multiword > 0);
    return 0;
}

// extract_two_generator accepts a candidate beta only when the HNF of the
// rows [a I; M_beta mod a] equals the ideal HNF (base4.c ok_elt).  Replay the
// whole search here with the plain fmpz_mat_hnf of those rows and compare the
// outcome, trial counts and returned element of ideal_two_generator, over
// maximal orders with trivial and nontrivial index, minimum 1, scalar ideals,
// ideals with content, and word and multiword minima.
void reduced_multiplication_rows(sflint::FmpzMat& out, silex::Order& order,
                                 const sflint::FmpzMat& coordinates,
                                 const sflint::Fmpz& modulus) noexcept {
    const slong n = order.degree();
    silex::OrderElement beta(order), basis_element(order), product(order);
    assert(beta.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
    sflint::FmpzMat unit(1, n), row(1, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_mat_zero(unit.raw());
        fmpz_one(fmpz_mat_entry(unit.raw(), 0, i));
        assert(basis_element.set_coordinates(sflint::FmpzMatConstRef(unit)));
        assert(product.multiply(beta, basis_element));
        assert(product.get_coordinates(sflint::FmpzMatRef(row)));
        for (slong j = 0; j < n; ++j) {
            fmpz_mod(fmpz_mat_entry(out.raw(), i, j),
                     fmpz_mat_entry(row.raw(), 0, j), modulus.raw());
        }
    }
}

bool plain_presentation_hnf_matches(const sflint::FmpzMat& h,
                                    const sflint::Fmpz& minimum,
                                    const sflint::FmpzMat& matrix) noexcept {
    const slong n = fmpz_mat_nrows(h.raw());
    sflint::FmpzMat rows(2 * n, n), reduced(2 * n, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i), minimum.raw());
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(rows.raw(), n + i, j),
                     fmpz_mat_entry(matrix.raw(), i, j));
        }
    }
    fmpz_mat_hnf(reduced.raw(), rows.raw());
    sflint::FmpzMatConstWindow top(reduced, 0, 0, n, n);
    return fmpz_mat_equal(top.raw(), h.raw()) != 0;
}

// Least positive integer in the ideal: the denominator of the coordinates of
// 1 in the basis given by the rows of h.
void plain_ideal_minimum(sflint::Fmpz& minimum, const sflint::FmpzMat& h,
                         const sflint::FmpzMat& one) noexcept {
    const slong n = fmpz_mat_nrows(h.raw());
    sflint::FmpqMat basis(n, n), inverse(n, n), target(1, n), solution(1, n);
    fmpq_mat_set_fmpz_mat(basis.raw(), h.raw());
    assert(fmpq_mat_inv(inverse.raw(), basis.raw()));
    fmpq_mat_set_fmpz_mat(target.raw(), one.raw());
    fmpq_mat_mul(solution.raw(), target.raw(), inverse.raw());
    fmpz_one(minimum.raw());
    for (slong j = 0; j < n; ++j) {
        fmpz_lcm(minimum.raw(), minimum.raw(),
                 fmpq_denref(fmpq_mat_entry(solution.raw(), 0, j)));
    }
}

struct PlainSearch {
    bool found = false;
    silex::detail::TwoGeneratorSearchReport report;
};

PlainSearch plain_two_generator_search(sflint::FmpzMat& found_row, silex::Order& order,
                                       const sflint::FmpzMat& h,
                                       const sflint::Fmpz& minimum,
                                       const silex::detail::TwoGeneratorSearchOptions& options) noexcept {
    const slong n = order.degree();
    PlainSearch search;
    sflint::FmpzMat row(1, n), matrix(n, n);
    std::vector<sflint::FmpzMat> rows, matrices;
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, j), fmpz_mat_entry(h.raw(), i, j));
        }
        reduced_multiplication_rows(matrix, order, row, minimum);
        if (fmpz_mat_is_zero(matrix.raw())) continue;
        if (options.basis_candidates) {
            ++search.report.basis_trials;
            if (plain_presentation_hnf_matches(h, minimum, matrix)) {
                fmpz_mat_set(found_row.raw(), row.raw());
                search.found = true;
                return search;
            }
        }
        rows.emplace_back(1, n);
        matrices.emplace_back(n, n);
        fmpz_mat_set(rows.back().raw(), row.raw());
        fmpz_mat_set(matrices.back().raw(), matrix.raw());
    }
    if (rows.empty()) {
        fmpz_mat_zero(matrix.raw());
        search.found = plain_presentation_hnf_matches(h, minimum, matrix);
        if (search.found) fmpz_mat_zero(found_row.raw());
        return search;
    }
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 0x53494c45UL, 0x58UL);
    sflint::Fmpz coefficient;
    for (slong trial = 0; trial < options.random_trials && !search.found; ++trial) {
        ++search.report.random_trials;
        fmpz_mat_zero(row.raw());
        for (std::size_t k = 0; k < rows.size(); ++k) {
            fmpz_randm(coefficient.raw(), state, minimum.raw());
            fmpz_mat_scalar_addmul_fmpz(row.raw(), rows[k].raw(), coefficient.raw());
        }
        fmpz_mat_scalar_mod_fmpz(row.raw(), row.raw(), minimum.raw());
        reduced_multiplication_rows(matrix, order, row, minimum);
        if (!fmpz_mat_is_zero(matrix.raw()) &&
            plain_presentation_hnf_matches(h, minimum, matrix)) {
            fmpz_mat_set(found_row.raw(), row.raw());
            search.found = true;
        }
    }
    flint_rand_clear(state);
    return search;
}

int test_two_generator_search_matches_plain_hnf() {
    sflint::FmpqPoly polynomial;
    slong searches = 0, rejected = 0, multiword = 0, unit_minimum = 0;
    for (slong which = 0; which < 6; ++which) {
        sflint::fmpq_poly_zero(polynomial);
        if (which == 0) {
            poly_x2_minus(polynomial, 5);
        } else if (which == 1) {
            poly_x3_minus(polynomial, 10);
        } else if (which == 2) {
            poly_x3_minus(polynomial, 2);
        } else if (which == 3) {
            silex::test::ideal_product::polynomial(polynomial, 4);
        } else if (which == 4) {
            for (slong e : {8, 7, 5, 4, 3, 1, 0}) {
                sflint::fmpq_poly_set_coeff_si(polynomial, e,
                        (e == 7 || e == 4 || e == 1) ? -1 : 1);
            }
        } else {
            sflint::fmpq_poly_set_coeff_si(polynomial, 6, 1);
            sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
            sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
        }
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order order(field);
        assert(order.maximal_order(equation) && order.is_maximal());
        const slong n = order.degree();

        std::vector<silex::Ideal> ideals;
        const slong large = (n == 2) ? (WORD(1) << 33) : (n == 3) ? WORD(3000000)
                          : (n == 4) ? (WORD(1) << 17) : WORD(3000);
        for (slong point : {WORD(2), WORD(3), WORD(-2), WORD(7), WORD(11), large}) {
            ideals.emplace_back(order);
            assert(silex::test::ideal_product::evaluation_ideal(ideals.back(), polynomial, point));
        }
        ideals.emplace_back(order);
        assert(ideals.back().multiply(ideals[0], ideals[1]));
        ideals.emplace_back(order);
        assert(ideals.back().multiply(ideals[0], ideals[0]));
        // O (minimum 1), 6 O, (2^64 + 13) O, and 6 times an evaluation ideal.
        sflint::Fmpz scale;
        fmpz_one(scale.raw());
        fmpz_mul_2exp(scale.raw(), scale.raw(), 64);
        fmpz_add_ui(scale.raw(), scale.raw(), 13);
        for (slong kind = 0; kind < 4; ++kind) {
            sflint::FmpzMat h(n, n);
            if (kind == 3) {
                assert(ideals[0].get_hnf(sflint::FmpzMatRef(h)));
            } else {
                fmpz_mat_one(h.raw());
            }
            if (kind == 2) {
                fmpz_mat_scalar_mul_fmpz(h.raw(), h.raw(), scale.raw());
            } else if (kind != 0) {
                fmpz_mat_scalar_mul_si(h.raw(), h.raw(), 6);
            }
            ideals.emplace_back(order);
            assert(ideals.back().set_hnf(sflint::FmpzMatConstRef(h)));
        }

        sflint::FmpzMat one_row(1, n), h(n, n), expected_row(1, n), actual_row(1, n);
        silex::OrderElement one(order);
        assert(one.one() && one.get_coordinates(sflint::FmpzMatRef(one_row)));
        for (const silex::Ideal& ideal : ideals) {
            assert(ideal.get_hnf(sflint::FmpzMatRef(h)));
            sflint::Fmpz minimum;
            plain_ideal_minimum(minimum, h, one_row);
            for (bool basis_candidates : {true, false}) {
                silex::detail::TwoGeneratorSearchOptions options;
                options.basis_candidates = basis_candidates;
                const PlainSearch expected = plain_two_generator_search(
                        expected_row, order, h, minimum, options);
                silex::detail::TwoGeneratorSearchReport report;
                sflint::Fmpz scalar;
                silex::OrderElement beta(order);
                const bool found = silex::detail::ideal_two_generator(
                        sflint::FmpzRef(scalar), beta, ideal, options, &report);
                ++searches;
                const slong trials = report.basis_trials + report.random_trials;
                rejected += trials - (found && trials > 0 ? 1 : 0);
                if (!fmpz_abs_fits_ui(minimum.raw())) ++multiword;
                if (fmpz_is_one(minimum.raw())) ++unit_minimum;
                if (found != expected.found ||
                    report.basis_trials != expected.report.basis_trials ||
                    report.random_trials != expected.report.random_trials) {
                    std::fprintf(stderr, "two-generator search mismatch: n=%ld\n", n);
                    std::fflush(stderr);
                }
                assert(found == expected.found && report.verified == found);
                assert(report.basis_trials == expected.report.basis_trials);
                assert(report.random_trials == expected.report.random_trials);
                if (!found) continue;
                assert(fmpz_equal(scalar.raw(), minimum.raw()));
                assert(beta.get_coordinates(sflint::FmpzMatRef(actual_row)));
                assert(fmpz_mat_equal(actual_row.raw(), expected_row.raw()));
            }
        }
    }
    assert(searches >= 100 && rejected >= 40 && multiword >= 20 && unit_minimum >= 12);
    return 0;
}

// Forward substitution in ideal_minimum against the fmpz_mat_solve route and
// the rational lcm of denominators; both must agree on the return value.
void check_ideal_minimum(const sflint::FmpzMat& h, const sflint::FmpzMat& one,
                         bool nonsingular) noexcept {
    sflint::Fmpz minimum, generic, plain;
    fmpz_set_si(minimum.raw(), -7);
    fmpz_set_si(generic.raw(), -7);
    const bool ok = silex::detail::ideal_minimum(
            minimum, sflint::FmpzMatConstRef(one), sflint::FmpzMatConstRef(h));
    const bool generic_ok = silex::detail::ideal_minimum_generic(
            generic, sflint::FmpzMatConstRef(one), sflint::FmpzMatConstRef(h));
    assert(ok == generic_ok && ok == nonsingular);
    if (!ok) return;
    assert(fmpz_equal(minimum.raw(), generic.raw()));
    plain_ideal_minimum(plain, h, one);
    assert(fmpz_equal(minimum.raw(), plain.raw()));
}

int test_ideal_minimum_forward_substitution_matches_solve() {
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 0x4d494e49UL, 0x4d554dUL);
    // Non-square-free and prime-power diagonals, a multiword prime power and
    // units; entry 0 marks a random diagonal below 2^70.
    const slong diagonals[] = {1, 1, 1, 2, 3, 4, 8, 9, 12, 16, 18, 25, 27, 36, 49, 0};
    const slong choices = static_cast<slong>(sizeof(diagonals) / sizeof(diagonals[0]));
    slong integral = 0, fractional = 0, square_diagonal = 0, multiword = 0;
    slong signed_cases = 0, dense = 0, singular = 0;
    for (slong trial = 0; trial < 1500; ++trial) {
        const slong n = 1 + static_cast<slong>(n_randint(state, 8));
        const slong form = trial % 3;  // 0 HNF, 1 HNF with prime-power powers, 2 signed
        sflint::FmpzMat h(n, n), one(1, n), z(1, n);
        bool square = false;
        for (slong j = 0; j < n; ++j) {
            fmpz* d = fmpz_mat_entry(h.raw(), j, j);
            const slong pick = diagonals[n_randint(state, static_cast<ulong>(choices))];
            if (pick == 0) {
                fmpz_randbits(d, state, 1 + static_cast<flint_bitcnt_t>(n_randint(state, 70)));
                fmpz_abs(d, d);
                if (fmpz_is_zero(d)) fmpz_one(d);
            } else {
                fmpz_set_si(d, pick);
            }
            if (form == 1 && n_randint(state, 3) == 0) {
                fmpz_ui_pow_ui(d, n_randint(state, 2) == 0 ? 2 : 3,
                               1 + n_randint(state, 70));
            }
            if (form == 2 && n_randint(state, 4) == 0) fmpz_neg(d, d);
            // 4 | d or 9 | d covers the listed non-square-free choices.
            square = square || fmpz_fdiv_ui(d, 4) == 0 || fmpz_fdiv_ui(d, 9) == 0;
            if (!fmpz_abs_fits_ui(d)) ++multiword;
            for (slong i = 0; i < j; ++i) {
                fmpz* entry = fmpz_mat_entry(h.raw(), i, j);
                if (form == 2) {
                    fmpz_randtest(entry, state, 40);
                } else {
                    sflint::Fmpz bound;
                    fmpz_abs(bound.raw(), d);
                    fmpz_randm(entry, state, bound.raw());
                }
            }
        }
        if (form != 2) assert(fmpz_mat_is_in_hnf(h.raw()));
        if (form == 2) ++signed_cases;
        if (square) ++square_diagonal;
        const slong target = trial % 4;
        if (target == 0) {
            fmpz_one(fmpz_mat_entry(one.raw(), 0, 0));
        } else if (target == 1) {
            // c = z H: x = z is integral, every division is exact, minimum 1.
            for (slong j = 0; j < n; ++j) fmpz_randtest(fmpz_mat_entry(z.raw(), 0, j), state, 20);
            fmpz_mat_mul(one.raw(), z.raw(), h.raw());
        } else if (target == 2) {
            fmpz_mat_randtest(one.raw(), state, 30);
        } else {
            // c = z H + e_k: integral plus one basis-vector perturbation.
            for (slong j = 0; j < n; ++j) fmpz_randtest(fmpz_mat_entry(z.raw(), 0, j), state, 10);
            fmpz_mat_mul(one.raw(), z.raw(), h.raw());
            fmpz* entry = fmpz_mat_entry(one.raw(), 0,
                    static_cast<slong>(n_randint(state, static_cast<ulong>(n))));
            fmpz_add_ui(entry, entry, 1);
        }
        check_ideal_minimum(h, one, true);
        sflint::Fmpz minimum;
        assert(silex::detail::ideal_minimum(minimum, sflint::FmpzMatConstRef(one),
                                            sflint::FmpzMatConstRef(h)));
        if (target == 1) assert(fmpz_is_one(minimum.raw()));
        ++(fmpz_is_one(minimum.raw()) ? integral : fractional);

        if (n >= 2) {
            // Not triangular: a swapped row pair takes the generic solve.
            sflint::FmpzMat swapped(n, n);
            fmpz_mat_set(swapped.raw(), h.raw());
            fmpz_mat_swap_rows(swapped.raw(), nullptr, 0, n - 1);
            check_ideal_minimum(swapped, one, true);
            ++dense;
            // Triangular but singular: both routes report failure.
            sflint::FmpzMat zeroed(n, n);
            fmpz_mat_set(zeroed.raw(), h.raw());
            const slong k = static_cast<slong>(n_randint(state, static_cast<ulong>(n)));
            for (slong j = 0; j < n; ++j) fmpz_zero(fmpz_mat_entry(zeroed.raw(), k, j));
            check_ideal_minimum(zeroed, one, false);
            ++singular;
        }
    }
    flint_rand_clear(state);
    assert(integral >= 300 && fractional >= 600 && square_diagonal >= 500 &&
           multiword >= 100 && signed_cases >= 400 && dense >= 1000 && singular >= 1000);

    // Stored HNFs of ideals in maximal orders, and their primitive parts.
    sflint::FmpqPoly polynomial;
    slong ideals_checked = 0, nonsquarefree = 0;
    for (slong which = 0; which < 3; ++which) {
        sflint::fmpq_poly_zero(polynomial);
        if (which == 0) {
            poly_x2_minus(polynomial, 5);
        } else if (which == 1) {
            poly_x3_minus(polynomial, 10);
        } else {
            silex::test::ideal_product::polynomial(polynomial, 4);
        }
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order order(field);
        assert(order.maximal_order(equation) && order.is_maximal());
        const slong n = order.degree();
        std::vector<silex::Ideal> ideals;
        for (slong point : {WORD(2), WORD(3), WORD(-2), WORD(7), WORD(12)}) {
            ideals.emplace_back(order);
            assert(silex::test::ideal_product::evaluation_ideal(ideals.back(), polynomial, point));
        }
        const std::size_t base = ideals.size();
        for (std::size_t i = 0; i < base; ++i) {
            ideals.emplace_back(order);
            assert(ideals.back().multiply(ideals[i], ideals[i]));
            ideals.emplace_back(order);
            assert(ideals.back().multiply(ideals[ideals.size() - 2], ideals[i]));
            ideals.emplace_back(order);
            assert(ideals.back().multiply(ideals[i], ideals[(i + 1) % base]));
        }
        sflint::FmpzMat one(1, n), h(n, n);
        silex::OrderElement unity(order);
        assert(unity.one() && unity.get_coordinates(sflint::FmpzMatRef(one)));
        for (const silex::Ideal& ideal : ideals) {
            assert(ideal.get_hnf(sflint::FmpzMatRef(h)));
            for (slong j = 0; j < n; ++j) {
                const fmpz* d = fmpz_mat_entry(h.raw(), j, j);
                if (fmpz_fdiv_ui(d, 4) == 0 || fmpz_fdiv_ui(d, 9) == 0) {
                    ++nonsquarefree;
                    break;
                }
            }
            check_ideal_minimum(h, one, true);
            sflint::Fmpz content;
            fmpz_mat_content(content.raw(), h.raw());
            fmpz_mat_scalar_divexact_fmpz(h.raw(), h.raw(), content.raw());
            check_ideal_minimum(h, one, true);
            ++ideals_checked;
        }
    }
    assert(ideals_checked >= 60 && nonsquarefree >= 10);
    return 0;
}

int test_nonprincipal_product_and_associativity() {
    sflint::FmpqPoly f;
    poly_x2_minus(f, -5);
    silex::NumberField field = field_by_polynomial(f);
    silex::Order equation = silex::Order::equation_order(field), order(field);
    assert(order.maximal_order(equation));
    silex::Element theta(field), shifted(field);
    silex::OrderElement beta(order), two(order);
    sflint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    assert(theta.gen() && shifted.add_si(theta, 1) && beta.set_element(shifted));
    silex::Ideal ideal(order), square(order), expected(order), left(order), right(order);
    // (2, 1+theta) in discriminant -20 has norm 2 and is nonprincipal:
    // a^2+5b^2=2 has no integral solutions. Its square is (2).
    assert(silex::detail::set_known_two_generator_ideal(ideal, sflint::FmpzConstRef(p), beta));
    sflint::Fmpz norm;
    assert(ideal.norm(sflint::FmpzRef(norm)) && fmpz_equal_ui(norm.raw(), 2));
    assert(two.set_si(2) && expected.set_principal(two));
    assert(square.multiply(ideal, ideal) && square.equal(expected));
    assert(left.multiply(square, ideal));
    assert(right.multiply(ideal, square) && right.equal(left));
    assert(square.set(ideal) && square.multiply(square, square) && square.equal(expected));
    return 0;
}

// Noninstalled differential-test output, in the common polynomial power basis.
int dump_product_fixtures() {
    for (slong degree : {2, 3, 4, 6, 8, 12}) {
        for (slong kind = 0; kind < 5; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            silex::Ideal product(order);
            assert(product.multiply(left, right));
            sflint::FmpzMat h(degree, degree);
            sflint::FmpqMat basis(degree, degree), rows(degree, degree), coordinates(degree, degree);
            assert(product.get_hnf(sflint::FmpzMatRef(h)));
            assert(order.get_basis(sflint::FmpqMatRef(basis)));
            fmpq_mat_set_fmpz_mat(rows.raw(), h.raw());
            fmpq_mat_mul(coordinates.raw(), rows.raw(), basis.raw());
            std::printf("%ld %ld [", degree, kind);
            for (slong i = 0; i < degree; ++i) {
                if (i) std::printf(";");
                for (slong j = 0; j < degree; ++j) {
                    if (j) std::printf(",");
                    fmpq_print(fmpq_mat_entry(coordinates.raw(), i, j));
                }
            }
            std::printf("]\n");
        }
    }
    return 0;
}

}  // namespace

// Deterministic pseudorandom source for the differential containment tests.
struct ContainsRng {
    std::uint64_t state;

    slong next(slong lo, slong hi) noexcept {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<slong>((state >> 33) %
                                       static_cast<std::uint64_t>(hi - lo + 1));
    }
};

// The containment rule Ideal::contains and Ideal::set_hnf used before the
// echelon back-substitution path in Lat::contains: a row lies in the lattice
// of the n x n basis iff HNF(basis) and HNF(basis; row) agree after trailing
// zero rows are dropped.
bool reference_hnf_contains(const sflint::FmpzMat& basis,
                            const sflint::FmpzMat& row) {
    const slong n = fmpz_mat_ncols(basis.raw());
    sflint::FmpzMat hnf(n, n);
    sflint::FmpzMat extended(n + 1, n);
    sflint::FmpzMat extended_hnf(n + 1, n);
    fmpz_mat_hnf(hnf.raw(), basis.raw());
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(extended.raw(), i, j),
                     fmpz_mat_entry(basis.raw(), i, j));
        }
    }
    for (slong j = 0; j < n; ++j) {
        fmpz_set(fmpz_mat_entry(extended.raw(), n, j),
                 fmpz_mat_entry(row.raw(), 0, j));
    }
    fmpz_mat_hnf(extended_hnf.raw(), extended.raw());
    if (fmpz_mat_is_zero_row(extended_hnf.raw(), n) == 0) {
        return false;
    }
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (fmpz_equal(fmpz_mat_entry(hnf.raw(), i, j),
                           fmpz_mat_entry(extended_hnf.raw(), i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

// Reference acceptance rule of Ideal::set_hnf: the input normalizes to a
// nonsingular HNF whose lattice is closed under multiplication by every order
// basis element, with membership decided by reference_hnf_contains.
bool reference_set_hnf_accepts(const silex::Order& order,
                               const sflint::FmpzMat& input,
                               sflint::FmpzMat& normalized) {
    const slong n = order.degree();
    fmpz_mat_hnf(normalized.raw(), input.raw());
    sflint::Fmpz det;
    fmpz_mat_det(det.raw(), normalized.raw());
    if (fmpz_is_zero(det.raw()) != 0) {
        return false;
    }
    silex::OrderElement gen(order);
    silex::OrderElement omega(order);
    silex::OrderElement product(order);
    sflint::FmpzMat row(1, n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                     fmpz_mat_entry(normalized.raw(), i, j));
        }
        assert(gen.set_coordinates(sflint::FmpzMatConstRef(row)));
        for (slong j = 0; j < n; ++j) {
            fmpz_mat_zero(row.raw());
            fmpz_one(fmpz_mat_entry(row.raw(), 0, j));
            assert(omega.set_coordinates(sflint::FmpzMatConstRef(row)));
            assert(product.multiply(gen, omega));
            assert(product.get_coordinates(sflint::FmpzMatRef(row)));
            if (!reference_hnf_contains(normalized, row)) {
                return false;
            }
        }
    }
    return true;
}

// Differential test of Ideal::contains(OrderElement) and Ideal::set_hnf
// against the HNF-comparison rule, on the product fixtures (maximal and
// equation orders, degrees 2-8) and their sums and products.  Elements are
// ideal combinations, perturbed combinations, arbitrary rows and halved
// combinations; set_hnf candidates are the ideal HNF, a unimodular scramble
// of it, twice it, a perturbed entry, a random triangular matrix and a
// singular matrix.
int test_contains_and_set_hnf_match_hnf_reference() {
    ContainsRng rng{0x7179c0de};
    slong members = 0;
    slong nonmembers = 0;
    slong accepted = 0;
    slong rejected = 0;
    for (slong degree : {2, 3, 4, 6, 8}) {
        for (slong kind = 0; kind < 5; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            silex::Ideal sum(order), product(order);
            assert(sum.add(left, right));
            assert(product.multiply(left, right));
            const slong n = order.degree();
            for (const silex::Ideal* ideal : {&left, &right, &sum, &product}) {
                sflint::FmpzMat hnf(n, n);
                assert(ideal->get_hnf(sflint::FmpzMatRef(hnf)));
                // The stored HNF is full-rank upper triangular, so these
                // samples take the triangular-solve route of contains.
                assert(silex::detail::is_full_rank_upper_triangular(hnf.raw()));

                sflint::FmpzMat row(1, n);
                sflint::Fmpz c;
                silex::OrderElement element(order);
                for (slong sample = 0; sample < 16; ++sample) {
                    const slong shape = sample % 4;
                    fmpz_mat_zero(row.raw());
                    if (shape != 2) {
                        for (slong i = 0; i < n; ++i) {
                            fmpz_set_si(c.raw(), rng.next(-4, 4));
                            for (slong j = 0; j < n; ++j) {
                                fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j),
                                            fmpz_mat_entry(hnf.raw(), i, j),
                                            c.raw());
                            }
                        }
                    }
                    if (shape == 1) {
                        fmpz_add_si(fmpz_mat_entry(row.raw(), 0, rng.next(0, n - 1)),
                                    fmpz_mat_entry(row.raw(), 0, rng.next(0, n - 1)),
                                    rng.next(1, 3));
                    } else if (shape == 2) {
                        for (slong j = 0; j < n; ++j) {
                            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, j),
                                        rng.next(-30, 30));
                        }
                    } else if (shape == 3) {
                        bool even = true;
                        for (slong j = 0; j < n; ++j) {
                            even = even && fmpz_is_even(fmpz_mat_entry(row.raw(), 0, j)) != 0;
                        }
                        for (slong j = 0; even && j < n; ++j) {
                            fmpz_divexact_si(fmpz_mat_entry(row.raw(), 0, j),
                                             fmpz_mat_entry(row.raw(), 0, j), 2);
                        }
                    }
                    assert(element.set_coordinates(sflint::FmpzMatConstRef(row)));
                    const bool expected = reference_hnf_contains(hnf, row);
                    if (ideal->contains(element) != expected) {
                        return 1;
                    }
                    ++(expected ? members : nonmembers);
                }

                for (slong candidate_kind = 0; candidate_kind < 6; ++candidate_kind) {
                    sflint::FmpzMat candidate(n, n);
                    fmpz_mat_set(candidate.raw(), hnf.raw());
                    if (candidate_kind == 1) {
                        for (slong step = 0; step < 3 * n; ++step) {
                            const slong target = rng.next(0, n - 1);
                            const slong source = rng.next(0, n - 1);
                            if (target == source) {
                                continue;
                            }
                            fmpz_set_si(c.raw(), rng.next(-2, 2));
                            for (slong j = 0; j < n; ++j) {
                                fmpz_addmul(fmpz_mat_entry(candidate.raw(), target, j),
                                            fmpz_mat_entry(candidate.raw(), source, j),
                                            c.raw());
                            }
                        }
                    } else if (candidate_kind == 2) {
                        fmpz_mat_scalar_mul_si(candidate.raw(), candidate.raw(), 2);
                    } else if (candidate_kind == 3) {
                        const slong i = rng.next(0, n - 2);
                        const slong j = rng.next(i + 1, n - 1);
                        fmpz_add_ui(fmpz_mat_entry(candidate.raw(), i, j),
                                    fmpz_mat_entry(candidate.raw(), i, j), 1);
                    } else if (candidate_kind == 4) {
                        fmpz_mat_zero(candidate.raw());
                        for (slong i = 0; i < n; ++i) {
                            const slong diagonal = rng.next(1, 3);
                            fmpz_set_si(fmpz_mat_entry(candidate.raw(), i, i), diagonal);
                            for (slong k = 0; k < i; ++k) {
                                fmpz_set_si(fmpz_mat_entry(candidate.raw(), k, i),
                                            rng.next(0, diagonal - 1));
                            }
                        }
                    } else if (candidate_kind == 5) {
                        fmpz_mat_zero(candidate.raw());
                        for (slong j = 0; j < n; ++j) {
                            fmpz_set(fmpz_mat_entry(candidate.raw(), 0, j),
                                     fmpz_mat_entry(hnf.raw(), 0, j));
                        }
                    }
                    sflint::FmpzMat normalized(n, n);
                    const bool expected = reference_set_hnf_accepts(order, candidate, normalized);
                    silex::Ideal out(order);
                    if (out.set_hnf(sflint::FmpzMatConstRef(candidate)) != expected) {
                        return 1;
                    }
                    if (expected) {
                        sflint::FmpzMat got(n, n);
                        if (!out.get_hnf(sflint::FmpzMatRef(got)) ||
                            fmpz_mat_equal(got.raw(), normalized.raw()) == 0) {
                            return 1;
                        }
                    }
                    ++(expected ? accepted : rejected);
                }
            }
        }
    }
    return members > 0 && nonmembers > 0 && accepted > 0 && rejected > 0 ? 0 : 1;
}

// Differential test of Ideal::contains(const Ideal&) against the
// HNF-comparison rule: J is in I iff every row of J's HNF is in the lattice of
// I's HNF.  The ideals are the product fixtures (maximal and equation orders,
// degrees 2-8), their sums, products and doubles and the unit ideal, compared
// in every ordered pair, so both inclusions and non-inclusions occur.
int test_contains_ideal_matches_hnf_reference() {
    slong included = 0;
    slong excluded = 0;
    for (slong degree : {2, 3, 4, 6, 8}) {
        for (slong kind = 0; kind < 6; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            const slong n = order.degree();
            silex::Ideal sum(order), product(order), two(order),
                    double_left(order), double_sum(order), unit(order);
            silex::OrderElement two_element(order);
            assert(sum.add(left, right));
            assert(product.multiply(left, right));
            assert(two_element.set_si(2));
            assert(two.set_principal(two_element));
            assert(double_left.multiply(two, left));
            assert(double_sum.multiply(two, sum));
            assert(unit.one());
            const silex::Ideal* ideals[] = {&left, &right, &sum, &product, &two,
                                            &double_left, &double_sum, &unit};
            for (const silex::Ideal* outer : ideals) {
                sflint::FmpzMat outer_hnf(n, n);
                assert(outer->get_hnf(sflint::FmpzMatRef(outer_hnf)));
                // Stored HNFs take the triangular-solve route of contains.
                assert(silex::detail::is_full_rank_upper_triangular(outer_hnf.raw()));
                for (const silex::Ideal* inner : ideals) {
                    sflint::FmpzMat inner_hnf(n, n);
                    assert(inner->get_hnf(sflint::FmpzMatRef(inner_hnf)));
                    bool expected = true;
                    sflint::FmpzMat row(1, n);
                    for (slong i = 0; expected && i < n; ++i) {
                        for (slong j = 0; j < n; ++j) {
                            fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                                     fmpz_mat_entry(inner_hnf.raw(), i, j));
                        }
                        expected = reference_hnf_contains(outer_hnf, row);
                    }
                    if (outer->contains(*inner) != expected) {
                        return 1;
                    }
                    ++(expected ? included : excluded);
                }
            }
        }
    }
    return included > 0 && excluded > 0 ? 0 : 1;
}

// A basis of the row lattice of the n x n HNF `hnf` (n >= 2) that is not
// upper triangular: the rows in reverse order, so row 0 has a zero first
// entry, followed by random unimodular row additions.
void reversed_scrambled_basis(sflint::FmpzMat& out,
                              const sflint::FmpzMat& hnf,
                              ContainsRng& rng) {
    const slong n = fmpz_mat_nrows(hnf.raw());
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(out.raw(), i, j),
                     fmpz_mat_entry(hnf.raw(), n - 1 - i, j));
        }
    }
    sflint::Fmpz c;
    for (slong step = 0; step < 2 * n; ++step) {
        const slong target = rng.next(1, n - 1);
        const slong source = rng.next(0, n - 1);
        if (target == source) {
            continue;
        }
        fmpz_set_si(c.raw(), rng.next(-2, 2));
        for (slong j = 0; j < n; ++j) {
            fmpz_addmul(fmpz_mat_entry(out.raw(), target, j),
                        fmpz_mat_entry(out.raw(), source, j), c.raw());
        }
    }
}

// The defensive lattice fallbacks of both Ideal::contains overloads, which no
// public construction route reaches: every route stores a full-rank upper
// triangular HNF.  The noninstalled detail::IdealAccess hook stores a
// non-triangular basis of a fixture ideal's lattice (rows reversed, so row 0
// has a zero first entry, then scrambled by unimodular row additions), and
// both overloads are compared with the HNF-comparison rule on that basis and
// with the triangular route of the fixture ideal itself.
int test_contains_lattice_fallback_matches_hnf_reference() {
    ContainsRng rng{0xfa11bac4};
    slong members = 0;
    slong nonmembers = 0;
    slong included = 0;
    slong excluded = 0;
    for (slong degree : {2, 3, 4, 6}) {
        for (slong kind = 0; kind < 5; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            const slong n = order.degree();
            silex::Ideal sum(order), product(order), unit(order);
            assert(sum.add(left, right));
            assert(product.multiply(left, right));
            assert(unit.one());
            const silex::Ideal* ideals[] = {&left, &right, &sum, &product, &unit};
            for (const silex::Ideal* outer : ideals) {
                sflint::FmpzMat hnf(n, n);
                assert(outer->get_hnf(sflint::FmpzMatRef(hnf)));
                sflint::FmpzMat basis(n, n);
                reversed_scrambled_basis(basis, hnf, rng);
                // The stored basis fails the triangular shape check, so both
                // overloads take the lattice fallback.
                assert(!silex::detail::is_full_rank_upper_triangular(basis.raw()));
                silex::Ideal stored(order);
                assert(silex::detail::IdealAccess::set_hnf_direct(
                        stored, sflint::FmpzMatConstRef(basis)));

                sflint::FmpzMat row(1, n);
                sflint::Fmpz c;
                silex::OrderElement element(order);
                for (slong sample = 0; sample < 12; ++sample) {
                    const slong shape = sample % 3;
                    fmpz_mat_zero(row.raw());
                    if (shape != 2) {
                        for (slong i = 0; i < n; ++i) {
                            fmpz_set_si(c.raw(), rng.next(-4, 4));
                            for (slong j = 0; j < n; ++j) {
                                fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j),
                                            fmpz_mat_entry(hnf.raw(), i, j),
                                            c.raw());
                            }
                        }
                    }
                    if (shape == 1) {
                        fmpz_add_si(fmpz_mat_entry(row.raw(), 0, rng.next(0, n - 1)),
                                    fmpz_mat_entry(row.raw(), 0, rng.next(0, n - 1)),
                                    rng.next(1, 3));
                    } else if (shape == 2) {
                        for (slong j = 0; j < n; ++j) {
                            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, j),
                                        rng.next(-30, 30));
                        }
                    }
                    assert(element.set_coordinates(sflint::FmpzMatConstRef(row)));
                    const bool expected = reference_hnf_contains(basis, row);
                    if (stored.contains(element) != expected ||
                        outer->contains(element) != expected) {
                        return 1;
                    }
                    ++(expected ? members : nonmembers);
                }

                for (const silex::Ideal* inner : ideals) {
                    sflint::FmpzMat inner_hnf(n, n);
                    assert(inner->get_hnf(sflint::FmpzMatRef(inner_hnf)));
                    bool expected = true;
                    for (slong i = 0; expected && i < n; ++i) {
                        for (slong j = 0; j < n; ++j) {
                            fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                                     fmpz_mat_entry(inner_hnf.raw(), i, j));
                        }
                        expected = reference_hnf_contains(basis, row);
                    }
                    if (stored.contains(*inner) != expected ||
                        outer->contains(*inner) != expected) {
                        return 1;
                    }
                    ++(expected ? included : excluded);
                }
                if (!stored.contains(stored)) {
                    return 1;
                }
            }
        }
    }
    return members > 0 && nonmembers > 0 && included > 0 && excluded > 0 ? 0 : 1;
}

// The noninstalled detail::IdealAccess::set_known_ideal_hnf hook, which
// library code uses for lattices that are ideals by construction, stores the
// same normalized HNF as the checked Ideal::set_hnf for any full-rank basis
// of an ideal's lattice, and fails like set_hnf, leaving the output
// unchanged, on a singular or misshapen matrix.  It does not test that the
// lattice is an ideal in release builds, so no non-ideal input is passed.
int test_set_known_ideal_hnf_matches_set_hnf() {
    ContainsRng rng{0x5e7c0de1};
    for (slong degree : {2, 3, 4, 6}) {
        for (slong kind = 0; kind < 5; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            const slong n = order.degree();
            silex::Ideal product(order);
            assert(product.multiply(left, right));
            const silex::Ideal* ideals[] = {&left, &right, &product};
            for (const silex::Ideal* source : ideals) {
                sflint::FmpzMat hnf(n, n);
                assert(source->get_hnf(sflint::FmpzMatRef(hnf)));
                sflint::FmpzMat basis(n, n);
                reversed_scrambled_basis(basis, hnf, rng);
                for (const sflint::FmpzMat* input : {&hnf, &basis}) {
                    silex::Ideal checked(order), known(order);
                    assert(checked.set_hnf(sflint::FmpzMatConstRef(*input)));
                    assert(silex::detail::IdealAccess::set_known_ideal_hnf(
                            known, sflint::FmpzMatConstRef(*input)));
                    sflint::FmpzMat want(n, n), got(n, n);
                    assert(checked.get_hnf(sflint::FmpzMatRef(want)));
                    assert(known.get_hnf(sflint::FmpzMatRef(got)));
                    assert(fmpz_mat_equal(want.raw(), got.raw()) != 0);
                    assert(fmpz_mat_equal(got.raw(), hnf.raw()) != 0);
                    assert(known.equal(*source));
                    sflint::Fmpz want_norm, got_norm;
                    assert(checked.norm(sflint::FmpzRef(want_norm)));
                    assert(known.norm(sflint::FmpzRef(got_norm)));
                    assert(fmpz_equal(want_norm.raw(), got_norm.raw()) != 0);
                }

                silex::Ideal kept(order);
                assert(kept.set(*source));
                sflint::FmpzMat singular(n, n);
                fmpz_mat_set(singular.raw(), hnf.raw());
                for (slong j = 0; j < n; ++j) {
                    fmpz_zero(fmpz_mat_entry(singular.raw(), n - 1, j));
                }
                assert(!silex::detail::IdealAccess::set_known_ideal_hnf(
                        kept, sflint::FmpzMatConstRef(singular)));
                sflint::FmpzMat wrong_shape(n - 1, n);
                assert(!silex::detail::IdealAccess::set_known_ideal_hnf(
                        kept, sflint::FmpzMatConstRef(wrong_shape)));
                assert(kept.equal(*source));
            }
        }
    }
    silex::Ideal undefined;
    sflint::FmpzMat identity(2, 2);
    fmpz_mat_one(identity.raw());
    assert(!silex::detail::IdealAccess::set_known_ideal_hnf(
            undefined, sflint::FmpzMatConstRef(identity)));
    return 0;
}

// A fractional ideal is built from a matrix only through the checked
// Ideal::set_hnf followed by FractionalIdeal::set_integral_den, which takes an
// Ideal and so does not repeat the ideal test.  On the product fixtures:
// set_hnf rejects every non-ideal HNF (decided by reference_set_hnf_accepts)
// and leaves the ideal unchanged, so set_integral_den still sees the earlier
// ideal; set_integral_den on gI with denominator d stores the checked
// set_hnf of gI/c over d/c, c = gcd(content(gI), d); and on maximal orders
// invert gives an inverse whose inverse is the input.
int test_fractional_from_validated_numerator() {
    ContainsRng rng{0x7a11d4e1};
    slong rejected = 0;
    slong reduced = 0;
    for (slong degree : {2, 3, 4, 6}) {
        for (slong kind = 0; kind < 6; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            const slong n = order.degree();
            silex::Ideal product(order);
            assert(product.multiply(left, right));
            const silex::Ideal* ideals[] = {&left, &right, &product};
            for (const silex::Ideal* source : ideals) {
                sflint::FmpzMat hnf(n, n);
                assert(source->get_hnf(sflint::FmpzMatRef(hnf)));

                for (slong trial = 0; trial < 8; ++trial) {
                    sflint::FmpzMat candidate(n, n);
                    fmpz_mat_set(candidate.raw(), hnf.raw());
                    const slong i = rng.next(0, n - 1);
                    const slong j = rng.next(0, n - 1);
                    fmpz_add_ui(fmpz_mat_entry(candidate.raw(), i, j),
                                fmpz_mat_entry(candidate.raw(), i, j),
                                static_cast<ulong>(rng.next(1, 3)));
                    sflint::FmpzMat normalized(n, n);
                    if (reference_set_hnf_accepts(order, candidate, normalized)) {
                        continue;
                    }
                    silex::Ideal kept(order);
                    assert(kept.set(*source));
                    assert(!kept.set_hnf(sflint::FmpzMatConstRef(candidate)));
                    assert(kept.equal(*source));
                    silex::FractionalIdeal from_kept(order), from_source(order);
                    assert(from_kept.set_integral(kept));
                    assert(from_source.set_integral(*source));
                    assert(from_kept.equal(from_source));
                    ++rejected;
                }

                for (slong g : {1, 6}) {
                    sflint::FmpzMat scaled(n, n);
                    fmpz_mat_scalar_mul_si(scaled.raw(), hnf.raw(), g);
                    silex::Ideal numerator(order);
                    assert(numerator.set_hnf(sflint::FmpzMatConstRef(scaled)));
                    for (slong d : {1, 4, 9}) {
                        sflint::Fmpz den, content;
                        fmpz_set_si(den.raw(), d);
                        fmpz_mat_content(content.raw(), scaled.raw());
                        fmpz_gcd(content.raw(), content.raw(), den.raw());
                        sflint::FmpzMat divided(n, n);
                        fmpz_mat_scalar_divexact_fmpz(divided.raw(), scaled.raw(),
                                                      content.raw());
                        sflint::Fmpz want_den;
                        fmpz_divexact(want_den.raw(), den.raw(), content.raw());
                        silex::Ideal want(order);
                        assert(want.set_hnf(sflint::FmpzMatConstRef(divided)));
                        if (!fmpz_is_one(content.raw())) {
                            ++reduced;
                        }

                        silex::FractionalIdeal frac(order);
                        assert(frac.set_integral_den(numerator,
                                                     sflint::FmpzConstRef(den)));
                        silex::Ideal got(order);
                        sflint::Fmpz got_den;
                        assert(frac.get_integral_den(got, sflint::FmpzRef(got_den)));
                        sflint::FmpzMat want_hnf(n, n), got_hnf(n, n);
                        assert(want.get_hnf(sflint::FmpzMatRef(want_hnf)));
                        assert(got.get_hnf(sflint::FmpzMatRef(got_hnf)));
                        assert(fmpz_mat_equal(want_hnf.raw(), got_hnf.raw()) != 0);
                        assert(fmpz_equal(want_den.raw(), got_den.raw()) != 0);

                        if (kind != 5) {
                            silex::FractionalIdeal inverse(order), back(order);
                            silex::FractionalIdeal unit(order), one(order);
                            assert(inverse.invert(frac));
                            assert(back.invert(inverse));
                            assert(back.equal(frac));
                            assert(unit.multiply(frac, inverse));
                            assert(one.one());
                            assert(unit.equal(one));
                        }
                    }
                }
            }
        }
    }
    return rejected > 0 && reduced > 0 ? 0 : 1;
}

// Failure returns of the public Ideal::contains overloads: an undefined ideal
// or element, an ideal without an HNF, and an element or ideal of another
// order (a different order of the same field, and an order of another field)
// give false.  The probes are 0 and the unit ideal's own data, which a
// successful test would accept, so each false is a failure return.
int test_contains_failure_returns() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);
    silex::NumberField field;
    field = field_by_polynomial(polynomial);
    silex::Order order;
    order = silex::test::equation_order(field);

    // The maximal order Z[(1 + x) / 2] of the same field.
    silex::Order maximal(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 1, 2);
    qmat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(maximal.set_basis(sflint::FmpqMatConstRef(basis)));

    sflint::FmpqPoly other_polynomial;
    poly_x2_minus(other_polynomial, 3);
    silex::NumberField other_field;
    other_field = field_by_polynomial(other_polynomial);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);

    silex::Ideal unit(order);
    assert(unit.one());
    silex::OrderElement zero(order);
    assert(zero.set_si(0));
    assert(unit.contains(zero));
    assert(unit.contains(unit));

    silex::Ideal undefined;
    assert(!undefined.is_defined());
    assert(!undefined.contains(zero));

    silex::Ideal without_hnf(order);
    assert(without_hnf.is_defined() && !without_hnf.has_hnf());
    assert(!without_hnf.contains(zero));

    silex::OrderElement undefined_element;
    assert(!undefined_element.is_defined());
    assert(!unit.contains(undefined_element));

    silex::OrderElement maximal_zero(maximal);
    assert(maximal_zero.set_si(0));
    assert(!unit.contains(maximal_zero));

    silex::OrderElement other_zero(other_order);
    assert(other_zero.set_si(0));
    assert(!unit.contains(other_zero));

    silex::Ideal maximal_unit(maximal);
    silex::Ideal other_unit(other_order);
    assert(maximal_unit.one() && other_unit.one());
    assert(!unit.contains(undefined));
    assert(!undefined.contains(unit));
    assert(!unit.contains(without_hnf));
    assert(!without_hnf.contains(unit));
    assert(!unit.contains(maximal_unit));
    assert(!maximal_unit.contains(unit));
    assert(!unit.contains(other_unit));
    assert(!other_unit.contains(unit));
    return 0;
}

// The triangular solve behind Ideal::contains, on bases no ideal stores:
// square upper triangular matrices with negative or large diagonal entries
// and unreduced entries above the diagonal, compared with the HNF-comparison
// rule.  Determinants that fit in a word (the route modulo |det B|) and that
// do not (the exact route) both occur, with single-word and multiword
// vectors.  Lower triangular, singular and non-square matrices fail the shape
// check, so Ideal::contains keeps the lattice test for them.
int test_triangular_basis_contains_matches_hnf_reference() {
    sflint::FmpzMat lower(2, 2);
    sflint::FmpzMat singular(2, 2);
    sflint::FmpzMat wide(2, 3);
    fmpz_set_si(fmpz_mat_entry(lower.raw(), 0, 0), 2);
    fmpz_set_si(fmpz_mat_entry(lower.raw(), 1, 0), 1);
    fmpz_set_si(fmpz_mat_entry(lower.raw(), 1, 1), 3);
    fmpz_set_si(fmpz_mat_entry(singular.raw(), 0, 0), 2);
    fmpz_set_si(fmpz_mat_entry(singular.raw(), 0, 1), 1);
    fmpz_set_si(fmpz_mat_entry(wide.raw(), 0, 0), 1);
    fmpz_set_si(fmpz_mat_entry(wide.raw(), 1, 1), 1);
    assert(!silex::detail::is_full_rank_upper_triangular(lower.raw()));
    assert(!silex::detail::is_full_rank_upper_triangular(singular.raw()));
    assert(!silex::detail::is_full_rank_upper_triangular(wide.raw()));

    ContainsRng rng{0x7a1a5017};
    slong members = 0;
    slong nonmembers = 0;
    for (slong n : {1, 2, 3, 5, 8}) {
        for (slong trial = 0; trial < 24; ++trial) {
            sflint::FmpzMat basis(n, n);
            for (slong i = 0; i < n; ++i) {
                slong diagonal = 0;
                while (diagonal == 0) {
                    diagonal = rng.next(-6, 6);
                }
                fmpz_set_si(fmpz_mat_entry(basis.raw(), i, i), diagonal);
                for (slong j = i + 1; j < n; ++j) {
                    fmpz_set_si(fmpz_mat_entry(basis.raw(), i, j),
                                rng.next(-40, 40));
                }
            }
            if (trial % 6 == 5) {
                fmpz_mul_2exp(fmpz_mat_entry(basis.raw(), 0, 0),
                              fmpz_mat_entry(basis.raw(), 0, 0), 80);
            }
            // |det B| near the word size: two pivots scaled by 2^31 give a
            // determinant that fits in a word for some bases and not others.
            if (trial % 6 == 4 && n >= 2) {
                fmpz_mul_2exp(fmpz_mat_entry(basis.raw(), 0, 0),
                              fmpz_mat_entry(basis.raw(), 0, 0), 31);
                fmpz_mul_2exp(fmpz_mat_entry(basis.raw(), n - 1, n - 1),
                              fmpz_mat_entry(basis.raw(), n - 1, n - 1), 31);
            }
            assert(silex::detail::is_full_rank_upper_triangular(basis.raw()));

            sflint::FmpzMat row(1, n);
            sflint::FmpzMat residual(1, n);
            sflint::Fmpz c;
            for (slong sample = 0; sample < 8; ++sample) {
                fmpz_mat_zero(row.raw());
                if (sample % 2 == 0) {
                    // An integer combination of the rows, perturbed in one
                    // entry on every fourth sample.
                    for (slong i = 0; i < n; ++i) {
                        fmpz_set_si(c.raw(), rng.next(-5, 5));
                        // Multiword coefficients exercise the reduction
                        // modulo |det B|.
                        if (sample >= 4) {
                            fmpz_mul_2exp(c.raw(), c.raw(), 90);
                            fmpz_add_si(c.raw(), c.raw(), rng.next(-5, 5));
                        }
                        for (slong j = 0; j < n; ++j) {
                            fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j),
                                        fmpz_mat_entry(basis.raw(), i, j),
                                        c.raw());
                        }
                    }
                    if (sample % 4 == 2) {
                        fmpz_add_si(fmpz_mat_entry(row.raw(), 0,
                                                   rng.next(0, n - 1)),
                                    fmpz_mat_entry(row.raw(), 0,
                                                   rng.next(0, n - 1)),
                                    1);
                    }
                } else {
                    for (slong j = 0; j < n; ++j) {
                        fmpz_set_si(fmpz_mat_entry(row.raw(), 0, j),
                                    rng.next(-60, 60));
                        if (sample >= 4) {
                            fmpz_mul_2exp(fmpz_mat_entry(row.raw(), 0, j),
                                          fmpz_mat_entry(row.raw(), 0, j),
                                          70);
                            fmpz_add_si(fmpz_mat_entry(row.raw(), 0, j),
                                        fmpz_mat_entry(row.raw(), 0, j),
                                        rng.next(-60, 60));
                        }
                    }
                }
                fmpz_mat_set(residual.raw(), row.raw());
                const bool expected = reference_hnf_contains(basis, row);
                if (silex::detail::triangular_basis_contains(
                            basis.raw(), fmpz_mat_row(residual.raw(), 0)) !=
                    expected) {
                    return 1;
                }
                ++(expected ? members : nonmembers);
            }
        }
    }
    // Fixed 2 x 2 bases with d = |det B| on each side of the word limit:
    // d = (2^32 - 1)(2^32 + 1) = 2^64 - 1 takes the word-size screen and
    // d = 2^32 2^32 = 2^64 takes the exact loop.  Each has a member (an
    // integer combination of the rows, with a multiword coefficient) and a
    // non-member (the same row with one entry shifted by 1).
    for (int above = 0; above < 2; ++above) {
        sflint::FmpzMat basis(2, 2);
        fmpz_set_ui(fmpz_mat_entry(basis.raw(), 0, 0),
                    above ? 4294967296UL : 4294967295UL);
        fmpz_set_ui(fmpz_mat_entry(basis.raw(), 1, 1),
                    above ? 4294967296UL : 4294967297UL);
        fmpz_set_si(fmpz_mat_entry(basis.raw(), 0, 1), 123456789);
        sflint::FmpzMat row(1, 2);
        sflint::Fmpz c0;
        sflint::Fmpz c1;
        fmpz_set_si(c0.raw(), 7);
        fmpz_mul_2exp(c0.raw(), c0.raw(), 90);
        fmpz_set_si(c1.raw(), -5);
        for (slong j = 0; j < 2; ++j) {
            fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j),
                        fmpz_mat_entry(basis.raw(), 0, j), c0.raw());
            fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j),
                        fmpz_mat_entry(basis.raw(), 1, j), c1.raw());
        }
        for (int member = 1; member >= 0; --member) {
            if (!member) {
                fmpz_add_ui(fmpz_mat_entry(row.raw(), 0, 1),
                            fmpz_mat_entry(row.raw(), 0, 1), 1);
            }
            assert(reference_hnf_contains(basis, row) == (member != 0));
            sflint::FmpzMat residual(1, 2);
            fmpz_mat_set(residual.raw(), row.raw());
            if (silex::detail::triangular_basis_contains(
                        basis.raw(), fmpz_mat_row(residual.raw(), 0)) !=
                (member != 0)) {
                return 1;
            }
        }
    }
    assert(members > 0 && nonmembers > 0);
    return 0;
}

// zero_skipping_mul is the H M product of the two-generator ideal product; it
// must return exactly what fmpz_mat_mul does.  The cases are random ideal-style
// HNFs (unit pivots give zero columns above them), HNFs with zero rows or
// columns, 1 x 1 and empty shapes, rectangular factors and dense signed left
// factors, with multiword entries on both sides.
void check_zero_skipping_mul(const sflint::FmpzMat& left, const sflint::FmpzMat& right,
                             flint_rand_t state) noexcept {
    const slong rows = fmpz_mat_nrows(left.raw());
    const slong columns = fmpz_mat_ncols(right.raw());
    sflint::FmpzMat expected(rows, columns), product(rows, columns);
    fmpz_mat_mul(expected.raw(), left.raw(), right.raw());
    // The output is overwritten, not accumulated into.
    fmpz_mat_randtest(product.raw(), state, 80);
    silex::detail::zero_skipping_mul(sflint::FmpzMatRef(product),
                                     sflint::FmpzMatConstRef(left),
                                     sflint::FmpzMatConstRef(right));
    assert(fmpz_mat_equal(product.raw(), expected.raw()));
    // The ideal product writes into the lower window of a 2n-row matrix.
    sflint::FmpzMat stacked(2 * rows, columns);
    fmpz_mat_randtest(stacked.raw(), state, 40);
    sflint::FmpzMat top(rows, columns);
    for (slong i = 0; i < rows; ++i)
        for (slong j = 0; j < columns; ++j)
            fmpz_set(fmpz_mat_entry(top.raw(), i, j), fmpz_mat_entry(stacked.raw(), i, j));
    {
        sflint::FmpzMatWindow lower(stacked, rows, 0, 2 * rows, columns);
        silex::detail::zero_skipping_mul(sflint::FmpzMatRef(lower.raw()),
                                         sflint::FmpzMatConstRef(left),
                                         sflint::FmpzMatConstRef(right));
    }
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < columns; ++j) {
            assert(fmpz_equal(fmpz_mat_entry(stacked.raw(), i, j),
                              fmpz_mat_entry(top.raw(), i, j)));
            assert(fmpz_equal(fmpz_mat_entry(stacked.raw(), rows + i, j),
                              fmpz_mat_entry(expected.raw(), i, j)));
        }
    }
}

int test_zero_skipping_mul_matches_fmpz_mat_mul() {
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 0x5a45524fUL, 0x534b4950UL);
    slong hnf = 0, unit_pivots = 0, zero_rows = 0, zero_columns = 0, one_by_one = 0;
    slong multiword = 0, dense = 0;
    for (slong trial = 0; trial < 600; ++trial) {
        const slong n = trial < 20 ? 1 : 1 + static_cast<slong>(n_randint(state, 24));
        const slong form = trial % 4;  // 0, 1 HNF; 2 HNF with zero rows/columns; 3 dense
        const flint_bitcnt_t bits = 1 + static_cast<flint_bitcnt_t>(n_randint(state, 130));
        sflint::FmpzMat h(n, n);
        if (form == 3) {
            fmpz_mat_randtest(h.raw(), state, bits);
            ++dense;
        } else {
            for (slong j = 0; j < n; ++j) {
                fmpz* d = fmpz_mat_entry(h.raw(), j, j);
                if (n_randint(state, 2) == 0) {
                    fmpz_one(d);
                    ++unit_pivots;
                } else {
                    fmpz_randbits(d, state, bits);
                    fmpz_abs(d, d);
                    if (fmpz_cmp_ui(d, 2) < 0) fmpz_set_ui(d, 2);
                }
                for (slong i = 0; i < j; ++i) {
                    fmpz_randm(fmpz_mat_entry(h.raw(), i, j), state, d);
                }
            }
            assert(fmpz_mat_is_in_hnf(h.raw()));
            ++hnf;
            if (form == 2) {
                const slong k = static_cast<slong>(n_randint(state, static_cast<ulong>(n)));
                if (trial % 8 == 2) {
                    for (slong j = 0; j < n; ++j) fmpz_zero(fmpz_mat_entry(h.raw(), k, j));
                    ++zero_rows;
                } else {
                    for (slong i = 0; i < n; ++i) fmpz_zero(fmpz_mat_entry(h.raw(), i, k));
                    ++zero_columns;
                }
            }
        }
        if (n == 1) ++one_by_one;
        if (fmpz_mat_max_bits(h.raw()) > FLINT_BITS - 2 ||
            fmpz_mat_max_bits(h.raw()) < -(FLINT_BITS - 2)) ++multiword;
        // Square, as in the ideal product, and one rectangular right factor.
        sflint::FmpzMat m(n, n);
        fmpz_mat_randtest(m.raw(), state, 1 + n_randint(state, 100));
        check_zero_skipping_mul(h, m, state);
        const slong columns = static_cast<slong>(n_randint(state, 5));
        sflint::FmpzMat wide(n, columns);
        fmpz_mat_randtest(wide.raw(), state, 1 + n_randint(state, 70));
        check_zero_skipping_mul(h, wide, state);
    }
    // All-zero and empty shapes.
    for (slong n : {0, 1, 3}) {
        sflint::FmpzMat zero(n, n), m(n, n);
        fmpz_mat_randtest(m.raw(), state, 20);
        check_zero_skipping_mul(zero, m, state);
        check_zero_skipping_mul(m, zero, state);
    }
    flint_rand_clear(state);
    assert(hnf >= 400 && unit_pivots >= 1000 && zero_rows >= 50 && zero_columns >= 50 &&
           one_by_one >= 20 && multiword >= 100 && dense >= 100);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--dump-product-fixtures")
        return dump_product_fixtures();
    assert(argc == 1);
    assert(test_two_generator_search_contracts() == 0);
    assert(test_general_product_contents_and_fallback() == 0);
    assert(test_nonprincipal_product_and_associativity() == 0);
    assert(test_known_two_generator_hnf_matches_generic() == 0);
    assert(test_general_product_hnf_matches_plain_hnf() == 0);
    assert(test_two_generator_search_matches_plain_hnf() == 0);
    assert(test_ideal_minimum_forward_substitution_matches_solve() == 0);
    assert(test_zero_skipping_mul_matches_fmpz_mat_mul() == 0);
    for (slong degree : {2, 3, 4, 6, 8, 12}) {
        for (slong kind = 0; kind < 5; ++kind) {
            silex::Order order;
            silex::Ideal left, right;
            assert(silex::test::ideal_product::fixture(order, left, right, degree, kind));
            silex::Ideal expected(order), product(order), reversed(order);
            assert(silex::test::ideal_product::exhaustive_product(expected, left, right));
            assert(product.multiply(left, right));
            if (!product.equal(expected)) {
                std::fprintf(stderr, "product mismatch: degree=%ld kind=%ld\n", degree, kind);
                sflint::FmpzMat actual(degree, degree), wanted(degree, degree);
                assert(product.get_hnf(sflint::FmpzMatRef(actual)));
                assert(expected.get_hnf(sflint::FmpzMatRef(wanted)));
                fmpz_mat_print_pretty(actual.raw());
                fmpz_mat_print_pretty(wanted.raw());
                std::fflush(stdout);
            }
            assert(product.equal(expected));
            assert(reversed.multiply(right, left));
            assert(reversed.equal(expected));
            if (kind == 4) {
                assert(reversed.set(left));
                assert(reversed.multiply(reversed, reversed));
                assert(reversed.equal(expected));
            }
            assert(left.multiply(left, right));
            assert(left.equal(expected));
        }
    }
    assert(test_degree_one_integral_arithmetic() == 0);
    assert(test_integral_arithmetic_nongalois_cubic() == 0);
    assert(test_integral_cubic_swap_and_self_set() == 0);
    assert(test_integral_arithmetic_quadratic() == 0);
    assert(test_integral_arithmetic_nontrivial_order_basis() == 0);
    assert(test_integral_principal_pseudorandom_identities() == 0);
    assert(test_principal_known_norm_matches_set_principal() == 0);
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
    assert(test_contains_and_set_hnf_match_hnf_reference() == 0);
    assert(test_contains_ideal_matches_hnf_reference() == 0);
    assert(test_contains_lattice_fallback_matches_hnf_reference() == 0);
    assert(test_set_known_ideal_hnf_matches_set_hnf() == 0);
    assert(test_fractional_from_validated_numerator() == 0);
    assert(test_contains_failure_returns() == 0);
    assert(test_triangular_basis_contains_matches_hnf_reference() == 0);
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
