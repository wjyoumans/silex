#include <silex/element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"
#include "../src/element/element_internal.hpp"

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/nf_elem.h>

#include <cassert>
#include <cstdint>
#include <initializer_list>
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

silex::NumberField sqrt_two_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

// Regression: the Hensel and residue-disproof prime loops must skip primes
// p with p | n.  They previously skipped only p with n | p, so exponents
// such as 10, 15 and 20 selected p = 5 and fell back to unsupported.
void test_power_skips_primes_dividing_exponent() noexcept {
    silex::NumberField field = sqrt_two_field();
    silex::Element theta(field);
    silex::Element base(field);
    assert(theta.gen());
    assert(base.add_si(theta, 1));

    for (slong n : {2, 3, 5, 6, 7, 10, 14, 15, 20, 21}) {
        sflint::Fmpz exponent;
        sflint::fmpz_set_si(sflint::FmpzRef(exponent), n);
        silex::Element power(field);
        silex::Element root(field);
        silex::Element check(field);
        assert(power.pow_fmpz(base, sflint::FmpzConstRef(exponent)));

        bool is_power = false;
        assert(power.is_power(is_power, root,
                              sflint::FmpzConstRef(exponent)));
        assert(is_power);
        assert(check.pow_fmpz(root, sflint::FmpzConstRef(exponent)));
        assert(check.equal(power));
    }
}

// Regression: a huge exponent must not build the dense residue polynomial
// y^n - a.  For these inputs the lifted candidates have norm other than +-1,
// so the norm pre-check also rejects them without forming candidate^n in K.
// Candidates of norm +-1 (units) are bounded by the height pre-check; see
// test_power_huge_exponent_unit_candidates.  1 + 5*theta has norm -49, so it
// is not a 2^40-th power.  Unsupported is acceptable; a definite answer must
// be false.
void test_power_huge_exponent_is_bounded() noexcept {
    silex::NumberField field = sqrt_two_field();
    silex::Element value(field);
    silex::Element root(field);
    sflint::FmpqPoly one_plus_five_theta;
    sflint::fmpq_poly_set_coeff_si(one_plus_five_theta, 0, 1);
    sflint::fmpq_poly_set_coeff_si(one_plus_five_theta, 1, 5);
    assert(value.set_fmpq_poly(
            sflint::FmpqPolyConstRef(one_plus_five_theta)));
    assert(root.set_si(7));

    sflint::Fmpz exponent;
    sflint::fmpz_one(sflint::FmpzRef(exponent));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(exponent),
                          sflint::FmpzConstRef(exponent), 40);
    bool is_power = true;
    if (value.is_power(is_power, root, sflint::FmpzConstRef(exponent))) {
        assert(!is_power);
    }
    assert(root.equal_si(7));

    // The fundamental unit 1 + theta is not a 2^40-th power either.
    sflint::FmpqPoly one_plus_theta;
    sflint::fmpq_poly_set_coeff_si(one_plus_theta, 0, 1);
    sflint::fmpq_poly_set_coeff_si(one_plus_theta, 1, 1);
    assert(value.set_fmpq_poly(sflint::FmpqPolyConstRef(one_plus_theta)));
    is_power = true;
    if (value.is_power(is_power, root, sflint::FmpzConstRef(exponent))) {
        assert(!is_power);
    }
    assert(root.equal_si(7));
}

// Regression: the residue-field membership disproof is valid at primes
// p | n.  In Q(sqrt 2), -theta is not a 10th power: at p = 5 (inert),
// (-theta)^((25 - 1)/gcd(10, 24)) = theta^12 = 2^6 = 4 != 1 in F_25.  Skipping
// p = 5 would reach p = 7 first, where no disproof exists.
void test_power_disproof_at_prime_dividing_exponent() noexcept {
    silex::NumberField field = sqrt_two_field();
    sflint::Fmpz exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), 10);

    for (slong constant : {0, -3}) {
        silex::Element value(field);
        silex::Element root(field);
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
        assert(value.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        assert(root.set_si(7));
        bool is_power = true;
        assert(value.is_power(is_power, root,
                              sflint::FmpzConstRef(exponent)));
        assert(!is_power);
        assert(root.equal_si(7));
    }
}

silex::NumberField field_from_coefficients(
        std::initializer_list<slong> coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    slong i = 0;
    for (slong c : coefficients) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i++, c);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::Element element_from_coefficients(
        const silex::NumberField& field,
        std::initializer_list<slong> coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    slong i = 0;
    for (slong c : coefficients) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i++, c);
    }
    silex::Element element(field);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    return element;
}

// Any definite answer must be correct; a failed query leaves root unchanged.
void assert_power_answer_consistent(const silex::Element& value,
                                    const sflint::Fmpz& exponent,
                                    bool expected_power) noexcept {
    const silex::NumberField& field = *value.parent();
    silex::Element root(field);
    assert(root.set_si(7));
    bool is_power = !expected_power;
    if (value.is_power(is_power, root, sflint::FmpzConstRef(exponent))) {
        assert(is_power == expected_power);
        if (is_power) {
            silex::Element check(field);
            assert(check.pow_fmpz(root, sflint::FmpzConstRef(exponent)));
            assert(check.equal(value));
        } else {
            assert(root.equal_si(7));
        }
    } else {
        assert(root.equal_si(7));
    }
}

// Regression: a lifted candidate of norm +-1 passes the norm
// pre-check, and verifying it formed candidate^n for n near 2^40 or 2^60,
// which did not finish.  The height pre-check (n log M(P_c) <= log M(P_a) + 1)
// now rejects such candidates before powering.  None of these units is an
// n-th power; the calls must return promptly, unsupported or false.
void test_power_huge_exponent_unit_candidates() noexcept {
    sflint::Fmpz two_40;
    sflint::fmpz_one(sflint::FmpzRef(two_40));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(two_40),
                          sflint::FmpzConstRef(two_40), 40);
    sflint::Fmpz three_25;
    sflint::fmpz_set_si(sflint::FmpzRef(three_25), 847288609443);
    // 1 + 24 * 5^24
    sflint::Fmpz one_plus_24_five_24;
    sflint::fmpz_set_si(sflint::FmpzRef(one_plus_24_five_24),
                        1430511474609375001);

    silex::NumberField sqrt_two = sqrt_two_field();
    assert_power_answer_consistent(
            element_from_coefficients(sqrt_two, {1, 1}),
            one_plus_24_five_24, false);
    assert_power_answer_consistent(
            element_from_coefficients(sqrt_two, {-3, -2}), two_40, false);
    assert_power_answer_consistent(
            element_from_coefficients(sqrt_two, {-3, 2}), two_40, false);

    silex::NumberField sqrt_three = field_from_coefficients({-3, 0, 1});
    assert_power_answer_consistent(
            element_from_coefficients(sqrt_three, {2, 1}), two_40, false);
    assert_power_answer_consistent(
            element_from_coefficients(sqrt_three, {2, -1}), two_40, false);

    silex::NumberField golden = field_from_coefficients({-1, -1, 1});
    for (auto coefficients : {std::pair<slong, slong>{-3, 2},
                              std::pair<slong, slong>{-1, -2},
                              std::pair<slong, slong>{1, 2},
                              std::pair<slong, slong>{3, -2}}) {
        assert_power_answer_consistent(
                element_from_coefficients(
                        golden, {coefficients.first, coefficients.second}),
                three_25, false);
    }
}

// The height pre-check must keep genuine roots: units and roots of unity
// (Mahler measure 1, where binary powering stays cheap for huge n) must be
// found with a verified root, and a non-power of norm 1 must not be reported
// as a power.
void assert_is_power_with_root(const silex::Element& value,
                               const sflint::Fmpz& exponent) noexcept {
    const silex::NumberField& field = *value.parent();
    silex::Element root(field);
    bool is_power = false;
    assert(value.is_power(is_power, root, sflint::FmpzConstRef(exponent)));
    assert(is_power);
    silex::Element check(field);
    assert(check.pow_fmpz(root, sflint::FmpzConstRef(exponent)));
    assert(check.equal(value));
}

void test_power_height_check_keeps_roots() noexcept {
    silex::NumberField golden = field_from_coefficients({-1, -1, 1});
    silex::Element phi = element_from_coefficients(golden, {0, 1});
    for (slong n : {3, 27, 243, 729}) {
        sflint::Fmpz exponent;
        sflint::fmpz_set_si(sflint::FmpzRef(exponent), n);
        silex::Element power(golden);
        assert(power.pow_fmpz(phi, sflint::FmpzConstRef(exponent)));
        assert_is_power_with_root(power, exponent);
    }

    // In Q(i): i = i^n for n = 1 (mod 4), including n = 1 + 4 * 5^25.
    silex::NumberField gaussian = imaginary_quadratic_field();
    silex::Element i_unit = element_from_coefficients(gaussian, {0, 1});
    for (slong n : {slong{5}, slong{13}, slong{1192092895507812501}}) {
        sflint::Fmpz exponent;
        sflint::fmpz_set_si(sflint::FmpzRef(exponent), n);
        assert_is_power_with_root(i_unit, exponent);
    }
    // -1 is not a 2^40-th power in Q(i).
    sflint::Fmpz two_40;
    sflint::fmpz_one(sflint::FmpzRef(two_40));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(two_40),
                          sflint::FmpzConstRef(two_40), 40);
    {
        silex::Element minus_one = element_from_coefficients(gaussian, {-1});
        silex::Element root(gaussian);
        assert(root.set_si(7));
        bool is_power = true;
        assert(minus_one.is_power(is_power, root,
                                  sflint::FmpzConstRef(two_40)));
        assert(!is_power);
        assert(root.equal_si(7));
    }

    // Non-integral coverage: ((1 + i)/2)^3 is found as a cube.
    sflint::Fmpz three;
    sflint::fmpz_set_si(sflint::FmpzRef(three), 3);
    {
        silex::Element half(gaussian);
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
        ::fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), 2);
        assert(half.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        silex::Element cube(gaussian);
        assert(cube.pow_fmpz(half, sflint::FmpzConstRef(three)));
        assert_is_power_with_root(cube, three);
    }

    // (3 + 4i)/5 has norm 1 but is not integral, and it is not a
    // (2^40 + 1)-th power.  The leading coefficient 5 of its primitive
    // characteristic polynomial is not a (2^40 + 1)-th power, so the Hensel
    // path rejects it before rescaling and the call is unsupported today.
    // Only consistency is checked: it must not report a power.
    silex::Element fraction(gaussian);
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, 3);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, 4);
        ::fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), 5);
        assert(fraction.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    }
    sflint::Fmpz two_40_plus_one;
    sflint::fmpz_add_ui(sflint::FmpzRef(two_40_plus_one),
                        sflint::FmpzConstRef(two_40), 1);
    assert_power_answer_consistent(fraction, two_40_plus_one, false);
}

silex::Element element_from_fraction(
        const silex::NumberField& field,
        std::initializer_list<slong> numerator,
        slong denominator) noexcept {
    silex::Element element = element_from_coefficients(field, numerator);
    assert(element.scalar_div_si(element, denominator));
    return element;
}

// Regression: non-integral powers whose roots have norm +-1 were
// never found.  The Hensel reconstruction recovers f'(theta) c in Z[theta],
// which needs an integral root c.  The reference is_power
// (src/NumField/NfAbs/Elem.jl) solves y^n = a d^n for the denominator d of a
// and returns y/d; Silex now does the same for non-integral a.  The found
// roots have norm +-1, so they also pass the height pre-check with a nonzero
// leading-coefficient term on both sides.  For (3 + 4i)/5 that term, log 5,
// is the whole of log M(P), since both conjugates have absolute value 1.
void assert_power_of_root_found(const silex::Element& base,
                                slong n) noexcept {
    const silex::NumberField& field = *base.parent();
    sflint::Fmpz exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), n);
    silex::Element power(field);
    assert(power.pow_fmpz(base, sflint::FmpzConstRef(exponent)));
    assert_is_power_with_root(power, exponent);
}

void test_power_non_integral_unit_norm_roots() noexcept {
    silex::NumberField gaussian = imaginary_quadratic_field();
    const silex::Element three_four =
            element_from_fraction(gaussian, {3, 4}, 5);
    for (slong n : {3, 4, 5, 6, 7}) {
        assert_power_of_root_found(three_four, n);
    }
    assert_power_of_root_found(element_from_fraction(gaussian, {5, 12}, 13),
                               3);
    assert_power_of_root_found(element_from_fraction(gaussian, {7, 24}, 25),
                               3);

    silex::NumberField sqrt_two = sqrt_two_field();
    for (slong n : {3, 5}) {
        assert_power_of_root_found(
                element_from_fraction(sqrt_two, {11, 6}, 7), n);
    }

    silex::NumberField sqrt_seven = field_from_coefficients({-7, 0, 1});
    assert_power_of_root_found(element_from_fraction(sqrt_seven, {4, 1}, 3),
                               3);

    // Q(zeta_8) = Q[x]/(x^4 + 1), with zeta_8^2 = i.  Squares here go
    // through the Hensel square path, not the quadratic formula.
    silex::NumberField zeta_eight = field_from_coefficients({1, 0, 0, 0, 1});
    for (slong n : {2, 3}) {
        assert_power_of_root_found(
                element_from_fraction(zeta_eight, {3, 0, 4}, 5), n);
    }
    assert_power_of_root_found(
            element_from_fraction(zeta_eight, {5, 0, 12}, 13), 2);

    // Denominator in Z[theta] but integral: ((1 + sqrt 5)/2)^7 needs no
    // rescaling and stays found.
    silex::NumberField sqrt_five = field_from_coefficients({-5, 0, 1});
    assert_power_of_root_found(element_from_fraction(sqrt_five, {1, 1}, 2),
                               7);

    // A rational radicand whose root is not rational, so the rational
    // shortcut in is_power_rational_constant cannot decide it and the query
    // falls through to the leading-coefficient code.  -1/4 is rational, but
    // a rational 4th power is nonnegative, so is_power_rational_constant
    // finds no rational 4th root and does not settle the query; the true
    // root, (1 + i)/2, is irrational.  Its primitive characteristic
    // polynomial over Q(i) is (4x + 1)^2 = 16x^2 + 8x + 1 (content 1),
    // giving lc = 16 = 2^4 and bits = 5 > 4, so the pre-filter passes and
    // primitive_characteristic_leading_coefficient runs.
    assert_power_of_root_found(element_from_fraction(gaussian, {1, 1}, 2), 4);

    // A non-integral root of non-unit norm: N((1 + 2i)/3) = 5/9, not +-1.
    // This exercises the norm-based verification in verify_power_root after
    // rescaling, not only the height check that guards unit-norm candidates.
    assert_power_of_root_found(element_from_fraction(gaussian, {1, 2}, 3), 3);

    // The power-basis denominator can exceed the denominator over the
    // maximal order O_K.  In Q(sqrt 5), O_K = Z[(1 + sqrt 5)/2], and
    // (1 + sqrt 5)/6 = ((1 + sqrt 5)/2)/3 has O_K-denominator 3, but its
    // power-basis (Z[sqrt 5]) denominator is 6.  Silex rescales by the
    // (possibly larger) power-basis denominator, per the documented
    // deviation from the reference's maximal-order choice; the exact check
    // in assert_power_of_root_found still requires the result to be exactly
    // correct.
    assert_power_of_root_found(
            element_from_fraction(sqrt_five, {1, 1}, 6), 3);

    // Non-powers stay non-powers: (3 + 4i)/5 is not a square or a cube, and
    // ((3 + 4i)/5)^2 is not a cube.  A definite answer must be false.
    sflint::Fmpz exponent;
    for (slong n : {2, 3}) {
        sflint::fmpz_set_si(sflint::FmpzRef(exponent), n);
        assert_power_answer_consistent(three_four, exponent, false);
    }
    silex::Element square(gaussian);
    assert(square.multiply(three_four, three_four));
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), 3);
    assert_power_answer_consistent(square, exponent, false);

    // A non-power that passes the leading-coefficient pre-filter at n = 3,
    // exercising the rescaled residue disproof (pure_power_hensel_prime_
    // search's definite-false return) rather than an early unsupported
    // return from the filter, or (at n = 2, in a quadratic-backend field)
    // the quadratic formula in is_square_quadratic, which never reaches the
    // Hensel path at all.  (-9 + 13i)/125 = (1 + i) * ((2 + i)/5)^3.  Its
    // primitive characteristic polynomial is 125x^2 + 18x + 2 (content 1),
    // so lc = 125 = 5^3 and bits = 7 > 3: the filter cannot reject it before
    // rescaling.  It is not a cube in Q(i): N(1 + i) = 2 is not a cube, so
    // 1 + i is not a cube, and a cube is a cube of a cube only if its
    // factor is; ((2 + i)/5)^3 is manifestly a cube, so a itself is a cube
    // only if 1 + i is.  The answer must be a definite false.
    {
        silex::Element non_power =
                element_from_fraction(gaussian, {-9, 13}, 125);
        silex::Element root(gaussian);
        assert(root.set_si(7));
        sflint::fmpz_set_si(sflint::FmpzRef(exponent), 3);
        bool is_power = true;
        assert(non_power.is_power(is_power, root,
                                  sflint::FmpzConstRef(exponent)));
        assert(!is_power);
        assert(root.equal_si(7));
    }
}

// x^8 - 4x^7 + 4x^6 + 2x^5 + 69x^4 - 146x^3 - 524x^2 + 598x + 1201, a
// defining polynomial of Q(sqrt(-3), sqrt(-5), sqrt(-7)) (the grh benchmark
// row multiquadratic_m3_m5_m7).  Its rational squares are the products of
// -3, -5 and -7 times rational squares.
silex::NumberField multiquadratic_m3_m5_m7_field() noexcept {
    return field_from_coefficients({1201, 598, -524, -146, 69, 2, 4, -4, 1});
}

void assert_definitely_not_power(const silex::Element& value,
                                 slong exponent) noexcept {
    sflint::Fmpz n;
    sflint::fmpz_set_si(sflint::FmpzRef(n), exponent);
    silex::Element root(*value.parent());
    assert(root.set_si(7));
    bool is_power = true;
    assert(value.is_power(is_power, root, sflint::FmpzConstRef(n)));
    assert(!is_power);
    assert(root.equal_si(7));
}

// Elements whose residue is an n-th power at the first factor of the first
// good prime, so the Hensel lift runs and verifies no candidate.  Before the
// pre-lift residue disproof these queries failed (unknown); a non-power
// residue at a further prime now proves "not a power".  Each value is
// checked to be no n-th power by hand: -1, 6 are not squares in the real
// field Q(2^(1/3)) of odd degree (a rational square of K has square norm,
// N(-1) = -1, N(6) = 216 = 6^3; and 6 = b^2 forces sqrt(6) in K, impossible
// in a cubic field); 1 + theta has norm 3 and 3 + theta norm 29, neither a
// rational square or cube.  2, 3, -1, 6 are not squares in
// Q(sqrt(-3), sqrt(-5), sqrt(-7)), whose rational square classes are those
// of -3, -5, -7 and their products, and 2 (1 + theta)^2 is a square there
// only if 2 is.
void test_power_prelift_residue_disproof() noexcept {
    silex::NumberField cubic = cubic_field();
    for (std::initializer_list<slong> c :
         {std::initializer_list<slong>{-1}, {6}, {1, 1}, {3, 1}}) {
        assert_definitely_not_power(element_from_coefficients(cubic, c), 2);
    }
    for (std::initializer_list<slong> c :
         {std::initializer_list<slong>{3}, {1, 1}, {3, 1}}) {
        assert_definitely_not_power(element_from_coefficients(cubic, c), 3);
    }

    silex::NumberField field = multiquadratic_m3_m5_m7_field();
    for (slong c : {2, 3, -1, 6}) {
        silex::Element value(field);
        assert(value.set_si(c));
        bool is_square = true;
        silex::Element root(field);
        assert(value.is_square(is_square, root));
        assert(!is_square);
    }
    silex::Element one_plus_theta = element_from_coefficients(field, {1, 1});
    silex::Element structured(field);
    silex::Element two(field);
    assert(two.set_si(2));
    assert(structured.multiply(one_plus_theta, one_plus_theta));
    assert(structured.multiply(structured, two));
    assert_definitely_not_power(structured, 2);

    // Genuine squares still answer yes with a verified root.
    for (slong c : {-3, 15, -105}) {
        silex::Element value(field);
        assert(value.set_si(c));
        assert(structured.multiply(one_plus_theta, one_plus_theta));
        assert(value.multiply(value, structured));
        sflint::Fmpz two;
        sflint::fmpz_set_si(sflint::FmpzRef(two), 2);
        assert_is_power_with_root(value, two);
    }
}

// The exact answers of the lift-only path (pre-lift disproof disabled) and
// of the production path, each followed by the residue disproof that
// Element::is_power runs after a failed lift.
struct PowerAnswer {
    bool known = false;
    bool is_power = false;
};

PowerAnswer hensel_power_answer(const silex::Element& value,
                                slong exponent,
                                silex::detail::PreLiftResidueDisproof mode,
                                silex::Element& root) noexcept {
    PowerAnswer answer;
    bool flag = false;
    if (silex::detail::pure_power_hensel_root(flag, root, value, exponent,
                                              nullptr, mode)) {
        answer.known = true;
        answer.is_power = flag;
        return answer;
    }
    flag = true;
    if (silex::detail::pure_power_residue_disproves(flag, value, exponent,
                                                    nullptr)) {
        assert(!flag);
        answer.known = true;
        answer.is_power = false;
    }
    return answer;
}

ulong next_test_random(ulong& state) noexcept {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    return state >> 33;
}

// Randomised differential test of the pre-lift residue disproof against the
// lift-only answers: it may only turn an unknown into "not a power".  Every
// definite answer stays, a found root is the same verified root, and no
// unknown or "no" becomes "yes".
void test_power_prelift_disproof_differential() noexcept {
    struct Case {
        silex::NumberField field;
        slong exponent;
        int samples;
    };
    Case cases[] = {
            {cubic_field(), 2, 60},
            {cubic_field(), 3, 60},
            {field_from_coefficients({1, 0, -10, 0, 1}), 2, 60},
            {field_from_coefficients({-3, 1, 0, 0, 0, 1}), 5, 30},
            {multiquadratic_m3_m5_m7_field(), 2, 16},
    };
    ulong state = 237;
    int unknown_to_no = 0;
    for (Case& test_case : cases) {
        const silex::NumberField& field = test_case.field;
        const slong degree = field.degree();
        sflint::Fmpz n;
        sflint::fmpz_set_si(sflint::FmpzRef(n), test_case.exponent);
        for (int sample = 0; sample < test_case.samples; ++sample) {
            sflint::FmpqPoly coefficients;
            for (slong i = 0; i < degree; ++i) {
                sflint::fmpq_poly_set_coeff_si(
                        coefficients, i,
                        static_cast<slong>(next_test_random(state) % 15) - 7);
            }
            silex::Element base(field);
            assert(base.set_fmpq_poly(
                    sflint::FmpqPolyConstRef(coefficients)));
            if (base.equal_si(0)) {
                continue;
            }
            // Kinds: b, b^n, c b^n with c in {2, 3, -1, 6}, and the
            // non-integral (b/3)^n and c (b/3)^n.
            const ulong kind = next_test_random(state) % 5;
            const slong scalars[] = {2, 3, -1, 6};
            const slong scalar =
                    scalars[next_test_random(state) % 4];
            silex::Element value(field);
            if (kind == 0) {
                assert(value.set(base));
            } else {
                if (kind >= 3) {
                    assert(base.scalar_div_si(base, 3));
                }
                assert(value.pow_fmpz(base, sflint::FmpzConstRef(n)));
                if (kind == 2 || kind == 4) {
                    silex::Element factor(field);
                    assert(factor.set_si(scalar));
                    assert(value.multiply(value, factor));
                }
            }

            silex::Element old_root(field);
            silex::Element new_root(field);
            const PowerAnswer old_answer = hensel_power_answer(
                    value, test_case.exponent,
                    silex::detail::PreLiftResidueDisproof::disabled,
                    old_root);
            const PowerAnswer new_answer = hensel_power_answer(
                    value, test_case.exponent,
                    silex::detail::PreLiftResidueDisproof::enabled, new_root);
            if (old_answer.known) {
                assert(new_answer.known);
                assert(new_answer.is_power == old_answer.is_power);
                if (old_answer.is_power) {
                    assert(new_root.equal(old_root));
                }
            } else if (new_answer.known) {
                assert(!new_answer.is_power);
                ++unknown_to_no;
            }
            if (new_answer.known && new_answer.is_power) {
                silex::Element check(field);
                assert(check.pow_fmpz(new_root, sflint::FmpzConstRef(n)));
                assert(check.equal(value));
            }
            if (kind == 1 || kind == 3) {
                // A genuine power is never disproved.
                assert(!new_answer.known || new_answer.is_power);
            }
        }
    }
    assert(unknown_to_no > 0);
}

void set_rational(silex::Element& element, slong numerator, ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator, denominator);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

// Replaces the FLINT context behind an already-defined field with a context
// for `reducible`, bypassing NumberField construction checks.  The public
// constructors reject reducible polynomials, so this white-box path is the
// only way to exercise Element::invert's own zero-divisor guard.
void reinitialize_flint_field(silex::NumberField& field,
                              const sflint::FmpqPoly& reducible) noexcept {
    nf_struct* raw = field.raw_flint_field();
    assert(raw != nullptr);
    nf_clear(raw);
    nf_init(raw, reducible.raw());
}

void check_invert_rejects_zero_divisor(silex::NumberField field,
                                       const sflint::FmpqPoly& reducible,
                                       slong root) noexcept {
    reinitialize_flint_field(field, reducible);

    silex::Element zero_divisor(field);
    silex::Element theta(field);
    silex::Element out(field);
    assert(theta.gen());
    assert(zero_divisor.add_si(theta, -root));
    assert(!zero_divisor.equal_si(0));

    assert(out.set_si(7));
    assert(!out.invert(zero_divisor));
    assert(out.equal_si(7));

    sflint::Fmpz exponent;
    sflint::fmpz_set_si(exponent, -1);
    assert(!out.pow_fmpz(zero_divisor, sflint::FmpzConstRef(exponent)));
    assert(out.equal_si(7));

    // Aliased invert (out == input): Element::invert writes into a separate
    // temporary and only swaps it in on success (src/element/element.cpp,
    // Element::invert), so a rejected zero divisor must fail the same way
    // and leave the aliased element's own value untouched.
    silex::Element aliased(field);
    assert(aliased.add_si(theta, -root));
    assert(!aliased.invert(aliased));
    assert(aliased.equal(zero_divisor));
}

void check_invert_round_trip(const silex::NumberField& field) noexcept {
    silex::Element theta(field);
    silex::Element value(field);
    silex::Element inverse(field);
    silex::Element product(field);
    assert(theta.gen());
    assert(value.add_si(theta, 1));
    assert(inverse.invert(value));
    assert(product.multiply(value, inverse));
    assert(product.equal_si(1));
    assert(value.invert(value));
    assert(value.equal(inverse));
}

void test_invert_zero_divisors() noexcept {
    // Nonmonic, nonintegral irreducible cubic (numerator 21x^3 + 14x + 30 is
    // 2-Eisenstein), plus the quadratic and linear representations.
    sflint::FmpqPoly nonmonic;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, 1, 2);
    sflint::fmpq_poly_set_coeff_fmpq(nonmonic, 3, coefficient);
    sflint::fmpq_set_si(coefficient, 1, 3);
    sflint::fmpq_poly_set_coeff_fmpq(nonmonic, 1, coefficient);
    sflint::fmpq_set_si(coefficient, 5, 7);
    sflint::fmpq_poly_set_coeff_fmpq(nonmonic, 0, coefficient);
    check_invert_round_trip(silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(nonmonic)));
    check_invert_round_trip(cubic_field());
    check_invert_round_trip(quadratic_field());
    check_invert_round_trip(imaginary_quadratic_field());
    check_invert_round_trip(degree_one_field());

    // Quadratic representation: x^2 - 4, zero divisor theta - 2.
    sflint::FmpqPoly quadratic_reducible;
    sflint::fmpq_poly_set_coeff_si(quadratic_reducible, 2, 1);
    sflint::fmpq_poly_set_coeff_si(quadratic_reducible, 0, -4);
    check_invert_rejects_zero_divisor(
            quadratic_field(), quadratic_reducible, 2);

    // The same quadratic (degree two) representation, but starting from a
    // NumberField whose own backend_kind() is generic rather than the
    // Silex quadratic backend: x^2 - 12 is irreducible but its radicand 12
    // is not squarefree, so by_polynomial() never installs the quadratic
    // backend for it (T-015).  Element::invert's own dispatch instead reads
    // FLINT's own nf_struct flag (element.cpp:262, NF_QUADRATIC), which is
    // set from the reinitialized degree-two polynomial regardless of
    // backend_kind(), so this exercises the same guard independently of
    // Silex's backend classification.
    sflint::FmpqPoly generic_quadratic_seed;
    sflint::fmpq_poly_set_coeff_si(generic_quadratic_seed, 2, 1);
    sflint::fmpq_poly_set_coeff_si(generic_quadratic_seed, 0, -12);
    silex::NumberField generic_quadratic_field =
            silex::test::field_by_polynomial(
                    sflint::FmpqPolyConstRef(generic_quadratic_seed));
    assert(generic_quadratic_field.backend_kind() ==
           silex::NumberFieldBackendKind::generic);
    check_invert_rejects_zero_divisor(
            generic_quadratic_field, quadratic_reducible, 2);

    // Generic representation: (x^2 + 1)(x - 3), zero divisor theta - 3.
    sflint::FmpqPoly cubic_reducible;
    sflint::fmpq_poly_set_coeff_si(cubic_reducible, 3, 1);
    sflint::fmpq_poly_set_coeff_si(cubic_reducible, 2, -3);
    sflint::fmpq_poly_set_coeff_si(cubic_reducible, 1, 1);
    sflint::fmpq_poly_set_coeff_si(cubic_reducible, 0, -3);
    check_invert_rejects_zero_divisor(cubic_field(), cubic_reducible, 3);

    // Nonmonic generic representation: 2(x - 1)(x^2 + 2), zero divisor
    // theta - 1.
    sflint::FmpqPoly nonmonic_reducible;
    sflint::fmpq_poly_set_coeff_si(nonmonic_reducible, 3, 2);
    sflint::fmpq_poly_set_coeff_si(nonmonic_reducible, 2, -2);
    sflint::fmpq_poly_set_coeff_si(nonmonic_reducible, 1, 4);
    sflint::fmpq_poly_set_coeff_si(nonmonic_reducible, 0, -4);
    check_invert_rejects_zero_divisor(silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(nonmonic)), nonmonic_reducible, 1);
}


// Element::norm computes the norm in fields x^n + c (n even, n >= 4) by
// repeated Graeffe halving of the defining polynomial (see
// src/element/element.cpp).  The result must be the same canonical rational
// as FLINT's nf_elem_norm, which every other field uses.
std::uint64_t next_random(std::uint64_t& state) noexcept {
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31U);
}

void random_integer(fmpz_t out, std::uint64_t& state, ulong max_bits) noexcept {
    const ulong bits = static_cast<ulong>(next_random(state) % (max_bits + 1));
    fmpz_zero(out);
    for (ulong done = 0; done < bits; done += 32) {
        fmpz_mul_2exp(out, out, 32);
        fmpz_add_ui(out, out, static_cast<ulong>(next_random(state) >> 32U));
    }
    if (bits > 0) {
        fmpz_fdiv_r_2exp(out, out, bits);
    }
    if ((next_random(state) & 1U) != 0) {
        fmpz_neg(out, out);
    }
}

void check_norm_matches_flint(const silex::NumberField& field,
                              const sflint::FmpqPoly& polynomial) noexcept {
    silex::Element element(field);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    sflint::Fmpq norm;
    sflint::Fmpq expected;
    assert(element.norm(sflint::FmpqRef(norm)));
    nf_elem_norm(expected.raw(), element.raw_flint_element(),
                 field.raw_flint_field());
    assert(fmpq_equal(norm.raw(), expected.raw()) != 0);
}

silex::NumberField binomial_field(slong degree, const fmpz_t constant) noexcept {
    sflint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_si(polynomial.raw(), degree, 1);
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, constant);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

slong check_binomial_norms(slong degree,
                           const fmpz_t constant,
                           slong count,
                           std::uint64_t& state) noexcept {
    silex::NumberField field = binomial_field(degree, constant);
    sflint::FmpqPoly polynomial;
    sflint::Fmpz numerator;
    sflint::Fmpz denominator;
    sflint::Fmpq coefficient;
    for (slong i = 0; i < count; ++i) {
        fmpq_poly_zero(polynomial.raw());
        const slong length =
                static_cast<slong>(next_random(state) % (degree + 1));
        const ulong shape = static_cast<ulong>(next_random(state) % 4);
        const ulong max_bits = shape == 0 ? 4 : (shape == 1 ? 40 : 200);
        for (slong k = 0; k < length; ++k) {
            if (shape == 3 && next_random(state) % 3 != 0) {
                continue;
            }
            random_integer(numerator.raw(), state, max_bits);
            random_integer(denominator.raw(), state, 24);
            fmpz_abs(denominator.raw(), denominator.raw());
            fmpz_add_ui(denominator.raw(), denominator.raw(), 1);
            fmpq_set_fmpz_frac(coefficient.raw(), numerator.raw(),
                               denominator.raw());
            fmpq_poly_set_coeff_fmpq(polynomial.raw(), k, coefficient.raw());
        }
        check_norm_matches_flint(field, polynomial);
    }
    return count;
}

void test_binomial_field_norm_matches_flint() noexcept {
    // Known values in Q(alpha), alpha^12 = -7: N(alpha) = 7 and
    // N(1 + alpha) = f(-1) = 8.
    sflint::Fmpz constant;
    fmpz_set_si(constant.raw(), 7);
    silex::NumberField field = binomial_field(12, constant.raw());
    silex::Element alpha(field);
    silex::Element shifted(field);
    sflint::Fmpq norm;
    assert(alpha.gen());
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, 7));
    assert(shifted.add_si(alpha, 1));
    assert(shifted.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, 8));
    assert(shifted.zero());
    assert(shifted.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, 0));

    // Irreducible x^n + c (Eisenstein, or cyclotomic for c = 1), plus odd
    // and non-binomial controls that take the generic path.
    struct Case {
        slong degree;
        slong constant;
    };
    const Case cases[] = {
            {4, 2},   {4, -6},  {4, 1},   {6, 5},   {6, -3},  {8, -2},
            {8, 5},   {8, 1},   {10, 5},  {12, 7},  {12, -10}, {14, 3},
            {16, -2}, {16, 1},  {18, 6},  {24, -5}, {7, -3},  {5, 2},
    };
    std::uint64_t state = 0x5eed1e5ULL;
    slong checked = 0;
    for (const Case& entry : cases) {
        fmpz_set_si(constant.raw(), entry.constant);
        checked += check_binomial_norms(
                entry.degree, constant.raw(), 600, state);
    }

    // Constants beyond one word: 3 (2^70 + 1) and -2 (2^64 + 1), Eisenstein
    // at 3 and 2.
    fmpz_one(constant.raw());
    fmpz_mul_2exp(constant.raw(), constant.raw(), 70);
    fmpz_add_ui(constant.raw(), constant.raw(), 1);
    fmpz_mul_ui(constant.raw(), constant.raw(), 3);
    checked += check_binomial_norms(8, constant.raw(), 300, state);
    fmpz_one(constant.raw());
    fmpz_mul_2exp(constant.raw(), constant.raw(), 64);
    fmpz_add_ui(constant.raw(), constant.raw(), 1);
    fmpz_mul_si(constant.raw(), constant.raw(), -2);
    checked += check_binomial_norms(6, constant.raw(), 300, state);

    // Non-binomial even-degree control: x^4 + x + 2 (irreducible) takes the
    // generic path.
    sflint::FmpqPoly control;
    fmpq_poly_set_coeff_si(control.raw(), 4, 1);
    fmpq_poly_set_coeff_si(control.raw(), 1, 1);
    fmpq_poly_set_coeff_si(control.raw(), 0, 2);
    silex::NumberField control_field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(control));
    sflint::FmpqPoly element_poly;
    fmpq_poly_set_coeff_si(element_poly.raw(), 3, -2);
    fmpq_poly_set_coeff_si(element_poly.raw(), 0, 5);
    check_norm_matches_flint(control_field, element_poly);

    assert(checked >= 10000);
}
}  // namespace

int main() {
    test_binomial_field_norm_matches_flint();
    test_invert_zero_divisors();
    test_power_skips_primes_dividing_exponent();
    test_power_huge_exponent_is_bounded();
    test_power_disproof_at_prime_dividing_exponent();
    test_power_huge_exponent_unit_candidates();
    test_power_height_check_keeps_roots();
    test_power_non_integral_unit_norm_roots();
    test_power_prelift_residue_disproof();
    test_power_prelift_disproof_differential();
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

    // theta = 2^(1/3) is not a cube: N(theta) = 2 is not a rational cube.
    // The first good prime leaves it undecided after the lift; the pre-lift
    // residue disproof at a further prime now answers "not a power" (this
    // query used to fail), and the root stays untouched.
    sflint::fmpz_set_ui(exponent, 3);
    is_power = true;
    assert(cubic_x.is_power(is_power, cubic_root,
                            sflint::FmpzConstRef(exponent)));
    assert(!is_power);
    assert(cubic_root.equal_si(7));

    return 0;
}
