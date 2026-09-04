#include <silex/element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>

namespace sflint = silex::flint;

namespace {

bool qcoeff_equal_si(
        const sflint::FmpqPoly& polynomial, slong i, slong expected) noexcept {
    sflint::Fmpq coeff;
    sflint::fmpq_poly_get_coeff_fmpq(
            sflint::FmpqRef(coeff), polynomial, i);
    return sflint::fmpq_equal_si(coeff, expected);
}

silex::NumberField quadratic_field() noexcept {
    return silex::test::quadratic_field(5);
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField full_residue_degree_cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField imaginary_quadratic_field() noexcept {
    return silex::test::quadratic_field(-1);
}

void set_rational(silex::Element& element, slong numerator, ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator, denominator);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

}  // namespace

int main() {
    silex::NumberField field = quadratic_field();
    silex::NumberField same_model = quadratic_field();

    silex::Element empty;
    assert(!empty.is_defined());
    assert(empty.parent() == nullptr);
    assert(empty.raw_flint_element() == nullptr);
    assert(!empty.zero());
    assert(!empty.equal_si(0));

    silex::Element x(field);
    assert(x.is_defined());
    assert(x.parent() != nullptr && x.parent()->has_same_data(field));
    assert(x.raw_flint_element() != nullptr);
    assert(x.equal_si(0));

    assert(x.one());
    assert(x.equal_si(1));
    assert(x.zero());
    assert(x.equal_si(0));
    assert(x.set_si(7));
    assert(x.equal_si(7));
    sflint::Fmpz large_integer;
    sflint::fmpz_one(sflint::FmpzRef(large_integer));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(large_integer),
                          sflint::FmpzConstRef(large_integer), 80);
    sflint::fmpz_add_ui(sflint::FmpzRef(large_integer),
                        sflint::FmpzConstRef(large_integer), 123);
    assert(!sflint::fmpz_fits_si(sflint::FmpzConstRef(large_integer)));
    assert(x.set_fmpz(sflint::FmpzConstRef(large_integer)));
    sflint::FmpqPoly large_poly;
    sflint::Fmpq large_coeff;
    assert(x.get_fmpq_poly(sflint::FmpqPolyRef(large_poly)));
    assert(sflint::fmpq_poly_degree(large_poly) == 0);
    sflint::fmpq_poly_get_coeff_fmpq(
            sflint::FmpqRef(large_coeff), large_poly, 0);
    assert(sflint::fmpz_equal(sflint::fmpq_num_ref(large_coeff),
                              sflint::FmpzConstRef(large_integer)));
    assert(sflint::fmpz_is_one(sflint::fmpq_den_ref(large_coeff)));

    silex::Element theta(field);
    assert(theta.gen());

    sflint::FmpqPoly theta_poly;
    assert(theta.get_fmpq_poly(sflint::FmpqPolyRef(theta_poly)));
    assert(sflint::fmpq_poly_degree(theta_poly) == 1);
    assert(qcoeff_equal_si(theta_poly, 0, 0));
    assert(qcoeff_equal_si(theta_poly, 1, 1));

    sflint::FmpqPoly input_poly;
    sflint::fmpq_poly_set_coeff_si(input_poly, 0, 3);
    sflint::fmpq_poly_set_coeff_si(input_poly, 1, 2);
    assert(x.set_fmpq_poly(sflint::FmpqPolyConstRef(input_poly)));

    sflint::FmpqPoly output_poly;
    assert(x.get_fmpq_poly(sflint::FmpqPolyRef(output_poly)));
    assert(sflint::fmpq_poly_degree(output_poly) == 1);
    assert(qcoeff_equal_si(output_poly, 0, 3));
    assert(qcoeff_equal_si(output_poly, 1, 2));
    auto owned_poly = x.to_fmpq_poly();
    assert(owned_poly.has_value());
    assert(sflint::fmpq_poly_degree(*owned_poly) == 1);
    assert(qcoeff_equal_si(*owned_poly, 0, 3));
    assert(qcoeff_equal_si(*owned_poly, 1, 2));

    silex::Element y(field);
    assert(y.set(x));
    assert(y.equal(x));

    silex::Element moved(std::move(y));
    assert(moved.is_defined());
    assert(!y.is_defined());
    assert(moved.equal(x));

    silex::Element assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(!moved.is_defined());
    assert(assigned.equal(x));

    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(assigned.raw_flint_element() == nullptr);
    assert(!assigned.to_fmpq_poly().has_value());
    assert(!assigned.set(x));
    assert(assigned.define(field));
    assert(assigned.set(x));
    assert(assigned.equal(x));

    silex::Element z(field);
    assert(z.set_si(-4));
    swap(assigned, z);
    assert(z.equal(x));
    assert(assigned.equal_si(-4));
    assert(assigned.parent() != nullptr &&
           assigned.parent()->has_same_data(field));

    silex::Element other_parent(same_model);
    assert(other_parent.set_si(7));
    assert(!other_parent.has_same_parent(assigned));
    assert(!other_parent.set(assigned));
    assert(other_parent.equal_si(7));
    assert(!other_parent.equal(assigned));

    sflint::Fmpq trace;
    sflint::Fmpq norm;
    assert(x.trace(sflint::FmpqRef(trace)));
    assert(sflint::fmpq_equal_si(trace, 6));
    assert(x.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, -11));

    silex::Element conjugate(field);
    assert(x.conjugate(conjugate));
    assert(conjugate.get_fmpq_poly(sflint::FmpqPolyRef(output_poly)));
    assert(qcoeff_equal_si(output_poly, 0, 3));
    assert(qcoeff_equal_si(output_poly, 1, -2));
    assert(conjugate.trace(sflint::FmpqRef(trace)));
    assert(sflint::fmpq_equal_si(trace, 6));
    assert(conjugate.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, -11));

    sflint::FmpqPoly generic_polynomial;
    sflint::fmpq_poly_set_coeff_si(generic_polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(generic_polynomial, 0, -12);

    silex::NumberField generic_field = silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(generic_polynomial));
    assert(generic_field.backend_kind() ==
            silex::NumberFieldBackendKind::generic);
    silex::Element generic_x(generic_field);
    silex::Element generic_out(generic_field);
    assert(generic_x.set_fmpq_poly(sflint::FmpqPolyConstRef(input_poly)));
    assert(generic_x.trace(sflint::FmpqRef(trace)));
    assert(sflint::fmpq_equal_si(trace, 6));
    assert(generic_x.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, -39));

    assert(generic_out.set_si(7));
    assert(!generic_x.conjugate(generic_out));
    assert(generic_out.equal_si(7));

    silex::Element half_theta(field);
    assert(half_theta.scalar_div_si(theta, 2));

    silex::Element product(field);
    assert(product.multiply(half_theta, half_theta));
    assert(product.norm(sflint::FmpqRef(norm)));
    sflint::Fmpq expected_norm;
    sflint::fmpq_set_si(expected_norm, 25, 16);
    assert(sflint::fmpq_equal(norm, expected_norm));

    silex::Element shifted(field);
    silex::Element negated(field);
    silex::Element cancelled(field);
    assert(shifted.add_si(half_theta, 3));
    assert(negated.negate(shifted));
    assert(cancelled.add(shifted, negated));
    assert(cancelled.equal_si(0));

    silex::Element half(field);
    assert(half.set_si_over_si(1, 2));
    assert(shifted.add_si(half, 3));
    assert(shifted.get_fmpq_poly(sflint::FmpqPolyRef(output_poly)));
    assert(qcoeff_equal_si(output_poly, 1, 0));
    sflint::Fmpq rational_constant;
    sflint::Fmpq expected_constant;
    sflint::fmpq_poly_get_coeff_fmpq(
            sflint::FmpqRef(rational_constant), output_poly, 0);
    sflint::fmpq_set_si(expected_constant, 7, 2);
    assert(sflint::fmpq_equal(rational_constant, expected_constant));
    assert(half.set_si(9));
    assert(!half.set_si_over_si(1, 0));
    assert(half.equal_si(9));
    assert(half_theta.set_si(11));
    assert(!half_theta.scalar_div_si(theta, 0));
    assert(half_theta.equal_si(11));
    assert(!half_theta.scalar_div_si(other_parent, 2));
    assert(half_theta.equal_si(11));
    assert(half_theta.scalar_div_si(theta, 2));

    silex::Element twice_theta(field);
    assert(twice_theta.add(theta, theta));
    assert(cancelled.subtract(twice_theta, theta));
    assert(cancelled.equal(theta));

    assert(twice_theta.set_si(42));
    assert(!twice_theta.add(theta, other_parent));
    assert(twice_theta.equal_si(42));
    assert(!twice_theta.multiply(theta, other_parent));
    assert(twice_theta.equal_si(42));

    silex::Element two(field);
    silex::Element zero(field);
    silex::Element expected(field);
    silex::Element power_out(field);
    sflint::Fmpz exponent;
    assert(two.set_si(2));
    assert(zero.zero());

    assert(power_out.invert(two));
    sflint::FmpqPoly half_poly;
    sflint::fmpq_poly_set_coeff_si(half_poly, 0, 1);
    sflint::fmpq_poly_scalar_div_ui(half_poly, half_poly, 2);
    assert(expected.set_fmpq_poly(sflint::FmpqPolyConstRef(half_poly)));
    assert(power_out.equal(expected));

    sflint::fmpz_set_si(exponent, 0);
    assert(power_out.pow_fmpz(zero, sflint::FmpzConstRef(exponent)));
    assert(power_out.equal_si(1));

    sflint::fmpz_set_si(exponent, 3);
    assert(power_out.pow_fmpz(two, sflint::FmpzConstRef(exponent)));
    assert(power_out.equal_si(8));

    sflint::fmpz_set_si(exponent, -1);
    assert(power_out.pow_fmpz(two, sflint::FmpzConstRef(exponent)));
    assert(power_out.equal(expected));

    assert(power_out.set_si(7));
    assert(!power_out.invert(zero));
    assert(power_out.equal_si(7));
    assert(!power_out.pow_fmpz(zero, sflint::FmpzConstRef(exponent)));
    assert(power_out.equal_si(7));

    sflint::fmpz_set_si(exponent, 4);
    assert(two.pow_fmpz(two, sflint::FmpzConstRef(exponent)));
    assert(two.equal_si(16));

    silex::Element one_plus_theta(field);
    assert(one_plus_theta.add_si(theta, 1));
    sflint::fmpz_set_si(exponent, 2);
    assert(power_out.pow_fmpz(one_plus_theta, sflint::FmpzConstRef(exponent)));
    assert(expected.multiply(one_plus_theta, one_plus_theta));
    assert(power_out.equal(expected));

    sflint::fmpz_set_si(exponent, -1);
    assert(power_out.pow_fmpz(one_plus_theta, sflint::FmpzConstRef(exponent)));
    assert(z.multiply(one_plus_theta, power_out));
    assert(z.equal_si(1));

    sflint::fmpz_set_si(exponent, 2);
    assert(power_out.pow_fmpz(one_plus_theta, sflint::FmpzConstRef(exponent)));
    sflint::fmpz_set_si(exponent, 3);
    assert(z.pow_fmpz(one_plus_theta, sflint::FmpzConstRef(exponent)));
    assert(expected.multiply(power_out, z));
    sflint::fmpz_set_si(exponent, 5);
    assert(power_out.pow_fmpz(one_plus_theta, sflint::FmpzConstRef(exponent)));
    assert(power_out.equal(expected));

    silex::NumberField rational_field = degree_one_field();
    silex::Element rational_x(rational_field);
    silex::Element rational_root(rational_field);
    silex::Element rational_check(rational_field);
    bool is_square = false;
    bool is_power = false;

    set_rational(rational_x, 9, 16);
    assert(rational_x.is_square(is_square, rational_root));
    assert(is_square);
    assert(rational_check.multiply(rational_root, rational_root));
    assert(rational_check.equal(rational_x));
    set_rational(rational_check, 3, 4);
    assert(rational_root.equal(rational_check));

    assert(rational_x.set_si(0));
    assert(rational_x.is_square(is_square, rational_root));
    assert(is_square);
    assert(rational_root.equal_si(0));

    assert(rational_root.set_si(7));
    is_square = true;
    assert(rational_x.set_si(2));
    assert(rational_x.is_square(is_square, rational_root));
    assert(!is_square);
    assert(rational_root.equal_si(7));

    is_square = true;
    assert(rational_x.set_si(-4));
    assert(rational_x.is_square(is_square, rational_root));
    assert(!is_square);
    assert(rational_root.equal_si(7));

    set_rational(rational_x, 27, 8);
    sflint::fmpz_set_ui(exponent, 3);
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_check.pow_fmpz(rational_root,
                                   sflint::FmpzConstRef(exponent)));
    assert(rational_check.equal(rational_x));
    set_rational(rational_check, 3, 2);
    assert(rational_root.equal(rational_check));

    assert(rational_x.set_si(-8));
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_root.equal_si(-2));

    assert(rational_root.set_si(7));
    sflint::fmpz_set_ui(exponent, 2);
    is_power = true;
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(rational_root.equal_si(7));

    assert(rational_x.set_si(16));
    sflint::fmpz_set_ui(exponent, 3);
    is_power = true;
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(rational_root.equal_si(7));

    assert(rational_x.set_si(0));
    sflint::fmpz_set_ui(exponent, 5);
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_root.equal_si(0));

    assert(rational_x.set_si(1));
    sflint::fmpz_one(exponent);
    sflint::fmpz_mul_2exp(exponent, exponent, 70);
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_root.equal_si(1));

    assert(rational_x.set_si(2));
    assert(rational_root.set_si(7));
    is_power = true;
    assert(rational_x.is_power(is_power, rational_root,
                               sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(rational_root.equal_si(7));

    assert(rational_x.set_si(11));
    sflint::fmpz_one(exponent);
    assert(rational_x.is_power(is_power, rational_x,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_x.equal_si(11));

    assert(rational_root.set_si(7));
    sflint::fmpz_zero(exponent);
    is_power = true;
    assert(!rational_x.is_power(is_power, rational_root,
                                sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(rational_root.equal_si(7));

    assert(one_plus_theta.add_si(theta, 1));
    assert(expected.multiply(one_plus_theta, one_plus_theta));
    assert(power_out.set(expected));
    assert(power_out.is_square(is_square, power_out));
    assert(is_square);
    assert(z.multiply(power_out, power_out));
    assert(z.equal(expected));

    assert(power_out.set_si(5));
    assert(power_out.is_square(is_square, power_out));
    assert(is_square);
    assert(z.multiply(power_out, power_out));
    assert(z.equal_si(5));

    assert(power_out.set_si(2));
    assert(z.set_si(7));
    is_square = true;
    assert(power_out.is_square(is_square, z));
    assert(!is_square);
    assert(z.equal_si(7));

    sflint::fmpz_set_ui(exponent, 2);
    assert(expected.multiply(one_plus_theta, one_plus_theta));
    assert(expected.is_power(is_power, power_out,
                             sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(z.pow_fmpz(power_out, sflint::FmpzConstRef(exponent)));
    assert(z.equal(expected));

    sflint::fmpz_set_ui(exponent, 3);
    assert(power_out.set_si(7));
    is_power = true;
    assert(expected.is_power(is_power, power_out,
                             sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(power_out.equal_si(7));

    silex::NumberField imaginary_field = imaginary_quadratic_field();
    silex::Element imaginary_theta(imaginary_field);
    silex::Element imaginary_x(imaginary_field);
    silex::Element imaginary_root(imaginary_field);
    silex::Element imaginary_check(imaginary_field);
    assert(imaginary_theta.gen());
    assert(imaginary_x.set_si(-1));
    assert(imaginary_x.is_square(is_square, imaginary_root));
    assert(is_square);
    assert(imaginary_check.multiply(imaginary_root, imaginary_root));
    assert(imaginary_check.equal(imaginary_x));
    assert(imaginary_root.equal(imaginary_theta));

    silex::NumberField cubic = cubic_field();
    silex::Element cubic_x(cubic);
    silex::Element cubic_root(cubic);
    silex::Element cubic_check(cubic);
    assert(cubic_x.set_si(8));
    assert(cubic_root.set_si(7));
    sflint::fmpz_set_ui(exponent, 3);
    is_power = false;
    assert(cubic_x.is_power(is_power, cubic_root,
                            sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(cubic_root.equal_si(2));
    assert(cubic_check.pow_fmpz(cubic_root, sflint::FmpzConstRef(exponent)));
    assert(cubic_check.equal(cubic_x));

    silex::Element cubic_theta(cubic);
    silex::Element cubic_shift(cubic);
    silex::Element cubic_cube(cubic);
    assert(cubic_theta.gen());
    assert(cubic_shift.add_si(cubic_theta, 1));
    assert(cubic_cube.pow_fmpz(cubic_shift, sflint::FmpzConstRef(exponent)));
    assert(cubic_root.set_si(7));
    is_power = false;
    assert(cubic_cube.is_power(is_power, cubic_root,
                               sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(cubic_check.pow_fmpz(cubic_root, sflint::FmpzConstRef(exponent)));
    assert(cubic_check.equal(cubic_cube));

    assert(cubic_x.set_si(16));
    assert(cubic_root.set_si(7));
    is_power = true;
    assert(cubic_x.is_power(is_power, cubic_root,
                            sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(cubic_check.pow_fmpz(cubic_root, sflint::FmpzConstRef(exponent)));
    assert(cubic_check.equal(cubic_x));

    assert(cubic_x.set_si(4));
    sflint::fmpz_set_ui(exponent, 2);
    assert(cubic_x.is_power(is_power, cubic_root,
                            sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(cubic_root.equal_si(2));

    assert(cubic_x.set_si(4));
    assert(cubic_root.set_si(7));
    is_square = false;
    assert(cubic_x.is_square(is_square, cubic_root));
    assert(is_square);
    assert(cubic_root.equal_si(2));

    assert(cubic_x.set_si(2));
    assert(cubic_root.set_si(7));
    is_square = true;
    assert(cubic_x.is_square(is_square, cubic_root));
    assert(!is_square);
    assert(cubic_root.equal_si(7));

    assert(cubic_x.gen());
    assert(cubic_root.set_si(7));
    is_square = true;
    assert(cubic_x.is_square(is_square, cubic_root));
    assert(!is_square);
    assert(cubic_root.equal_si(7));

    sflint::fmpz_set_ui(exponent, 2);
    is_power = true;
    assert(cubic_x.is_power(is_power, cubic_root,
                            sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(cubic_root.equal_si(7));

    silex::NumberField full_residue_cubic =
        full_residue_degree_cubic_field();
    silex::Element full_residue_theta(full_residue_cubic);
    silex::Element full_residue_shift(full_residue_cubic);
    silex::Element full_residue_square(full_residue_cubic);
    silex::Element full_residue_root(full_residue_cubic);
    silex::Element full_residue_check(full_residue_cubic);
    assert(full_residue_theta.gen());
    assert(full_residue_shift.add_si(full_residue_theta, 2));
    assert(full_residue_square.multiply(full_residue_shift,
                                        full_residue_shift));
    assert(full_residue_root.set_si(7));
    is_square = false;
    assert(full_residue_square.is_square(is_square, full_residue_root));
    assert(is_square);
    assert(full_residue_check.multiply(full_residue_root,
                                       full_residue_root));
    assert(full_residue_check.equal(full_residue_square));

    assert(full_residue_root.set_si(7));
    sflint::fmpz_set_ui(exponent, 2);
    is_power = false;
    assert(full_residue_square.is_power(
            is_power, full_residue_root, sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(full_residue_check.pow_fmpz(
            full_residue_root, sflint::FmpzConstRef(exponent)));
    assert(full_residue_check.equal(full_residue_square));

    sflint::fmpz_set_ui(exponent, 3);
    is_power = true;
    assert(!cubic_x.is_power(is_power, cubic_root,
                             sflint::FmpzConstRef(exponent)));
    assert(is_power);
    assert(cubic_root.equal_si(7));

    return 0;
}
