#include <silex/flint/fmpq.hpp>
#include <silex/unit.hpp>

#include "element/element_internal.hpp"
#include "test_support.hpp"
#include "unit/unit_internal.hpp"

#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>

#include <cassert>
#include <vector>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField pure_quadratic_polynomial_field(
        slong radicand) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -radicand);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField x4_plus_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_trivial_good_prime_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_zeta4_subfield_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_zeta6_subfield_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField generic_x2_plus_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

bool element_pow_si(silex::Element& out,
                    const silex::Element& input,
                    slong exponent) noexcept {
    sflint::Fmpz exp;
    fmpz_set_si(exp.raw(), exponent);
    return out.pow_fmpz(input, sflint::FmpzConstRef(exp));
}

bool set_quadratic_coeffs(silex::Element& out,
                          slong constant_num,
                          ulong constant_den,
                          slong linear_num,
                          ulong linear_den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;

    sflint::fmpq_set_si(coefficient, constant_num, constant_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    sflint::fmpq_set_si(coefficient, linear_num, linear_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool set_integral_quadratic_coeffs_str(silex::Element& out,
                                       const char* constant,
                                       const char* linear) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::Fmpz integer;

    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), constant)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), linear)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

int test_rank() {
    silex::NumberField degree_one = degree_one_field();
    silex::NumberField real_quadratic = quadratic_field(2);
    silex::NumberField imaginary_quadratic = quadratic_field(-1);

    slong rank = -1;
    assert(silex::unit_rank(rank, degree_one));
    assert(rank == 0);
    assert(silex::unit_rank(rank, real_quadratic));
    assert(rank == 1);
    assert(silex::unit_rank(rank, imaginary_quadratic));
    assert(rank == 0);

    silex::NumberField unset;
    rank = 77;
    assert(!silex::unit_rank(rank, unset));
    assert(rank == 77);
    return 0;
}

int test_quadratic_fundamental_unit() {
    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Element eps2(sqrt2);
    silex::Element expected2(sqrt2);
    assert(silex::quadratic_fundamental_unit(eps2, sqrt2));
    assert(set_quadratic_coeffs(expected2, 1, 1, 1, 1));
    assert(eps2.equal(expected2));

    silex::NumberField polynomial_sqrt2 = pure_quadratic_polynomial_field(2);
    assert(polynomial_sqrt2.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    silex::Element generic_eps2(polynomial_sqrt2);
    silex::Element generic_expected2(polynomial_sqrt2);
    assert(silex::quadratic_fundamental_unit(generic_eps2,
                                             polynomial_sqrt2));
    assert(set_quadratic_coeffs(generic_expected2, 1, 1, 1, 1));
    assert(generic_eps2.equal(generic_expected2));

    silex::NumberField sqrt5 = quadratic_field(5);
    silex::Element eps5(sqrt5);
    silex::Element expected5(sqrt5);
    assert(silex::quadratic_fundamental_unit(eps5, sqrt5));
    assert(set_quadratic_coeffs(expected5, 1, 2, 1, 2));
    assert(eps5.equal(expected5));

    silex::NumberField polynomial_sqrt5 = pure_quadratic_polynomial_field(5);
    assert(polynomial_sqrt5.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    silex::Element generic_eps5(polynomial_sqrt5);
    silex::Element generic_expected5(polynomial_sqrt5);
    assert(silex::quadratic_fundamental_unit(generic_eps5,
                                             polynomial_sqrt5));
    assert(set_quadratic_coeffs(generic_expected5, 1, 2, 1, 2));
    assert(generic_eps5.equal(generic_expected5));

    silex::NumberField sqrt13 = quadratic_field(13);
    silex::Element eps13(sqrt13);
    silex::Element expected13(sqrt13);
    assert(silex::quadratic_fundamental_unit(eps13, sqrt13));
    assert(set_quadratic_coeffs(expected13, 3, 2, 1, 2));
    assert(eps13.equal(expected13));

    silex::NumberField sqrt17 = quadratic_field(17);
    silex::Element eps17(sqrt17);
    silex::Element expected17(sqrt17);
    assert(silex::quadratic_fundamental_unit(eps17, sqrt17));
    assert(set_quadratic_coeffs(expected17, 4, 1, 1, 1));
    assert(eps17.equal(expected17));

    silex::NumberField polynomial_sqrt210 =
            pure_quadratic_polynomial_field(210);
    silex::Element eps210(polynomial_sqrt210);
    silex::Element expected210(polynomial_sqrt210);
    assert(silex::quadratic_fundamental_unit(eps210,
                                             polynomial_sqrt210));
    assert(set_quadratic_coeffs(expected210, 29, 1, 2, 1));
    assert(eps210.equal(expected210));

    silex::NumberField sqrt17345 = quadratic_field(17345);
    silex::Element eps17345(sqrt17345);
    silex::Element expected17345(sqrt17345);
    sflint::Fmpq norm17345;
    assert(silex::quadratic_fundamental_unit(eps17345, sqrt17345));
    assert(set_integral_quadratic_coeffs_str(
            expected17345, "301977958012", "2292915721"));
    assert(eps17345.equal(expected17345));
    assert(eps17345.norm(sflint::FmpqRef(norm17345)));
    assert(sflint::fmpq_equal_si(norm17345, -1));

    silex::NumberField polynomial_sqrt17345 =
            pure_quadratic_polynomial_field(17345);
    assert(polynomial_sqrt17345.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    silex::Element generic_eps17345(polynomial_sqrt17345);
    silex::Element generic_expected17345(polynomial_sqrt17345);
    assert(silex::quadratic_fundamental_unit(generic_eps17345,
                                             polynomial_sqrt17345));
    assert(set_integral_quadratic_coeffs_str(
            generic_expected17345, "301977958012", "2292915721"));
    assert(generic_eps17345.equal(generic_expected17345));
    assert(generic_eps17345.norm(sflint::FmpqRef(norm17345)));
    assert(sflint::fmpq_equal_si(norm17345, -1));

    silex::NumberField qi = quadratic_field(-1);
    silex::Element imaginary(qi);
    assert(imaginary.set_si(7));
    assert(!silex::quadratic_fundamental_unit(imaginary, qi));
    assert(imaginary.equal_si(7));

    silex::NumberField degree_one = degree_one_field();
    silex::Element generic(degree_one);
    assert(generic.set_si(7));
    assert(!silex::quadratic_fundamental_unit(generic, degree_one));
    assert(generic.equal_si(7));

    silex::NumberField nonsquarefree =
            pure_quadratic_polynomial_field(12);
    silex::Element unsupported(nonsquarefree);
    assert(unsupported.set_si(7));
    assert(!silex::quadratic_fundamental_unit(unsupported, nonsquarefree));
    assert(unsupported.equal_si(7));

    sflint::FmpqPoly shifted_polynomial;
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 0, 1);
    silex::NumberField shifted = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(shifted_polynomial));
    silex::Element shifted_output(shifted);
    assert(shifted_output.set_si(7));
    assert(!silex::quadratic_fundamental_unit(shifted_output, shifted));
    assert(shifted_output.equal_si(7));

    return 0;
}

int test_roots_of_unity() {
    silex::NumberField degree_one = degree_one_field();
    sflint::Fmpz order;
    silex::Element generator(degree_one);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), generator,
                                  degree_one));
    assert(sflint::fmpz_equal_si(order, 2));
    assert(generator.equal_si(-1));

    silex::NumberField qi = quadratic_field(-1);
    silex::Element i(qi);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), i, qi));
    assert(sflint::fmpz_equal_si(order, 4));
    silex::Element i2(qi);
    assert(element_pow_si(i2, i, 2));
    assert(i2.equal_si(-1));

    silex::NumberField qzeta6 = quadratic_field(-3);
    silex::Element zeta6(qzeta6);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), zeta6, qzeta6));
    assert(sflint::fmpz_equal_si(order, 6));
    silex::Element zeta6_cubed(qzeta6);
    assert(element_pow_si(zeta6_cubed, zeta6, 3));
    assert(zeta6_cubed.equal_si(-1));

    silex::NumberField x4p1 = x4_plus_one_field();
    silex::Element zeta8(x4p1);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), zeta8, x4p1));
    assert(sflint::fmpz_equal_si(order, 8));
    silex::Element zeta8_4(x4p1);
    assert(element_pow_si(zeta8_4, zeta8, 4));
    assert(zeta8_4.equal_si(-1));

    silex::NumberField quartic_trivial = quartic_trivial_good_prime_field();
    silex::Element trivial_generator(quartic_trivial);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), trivial_generator,
                                  quartic_trivial));
    assert(sflint::fmpz_equal_si(order, 2));
    assert(trivial_generator.equal_si(-1));
    assert(silex::root_of_unity_order(sflint::FmpzRef(order),
                                       quartic_trivial));
    assert(sflint::fmpz_equal_si(order, 2));
    assert(silex::root_of_unity_generator(trivial_generator,
                                           quartic_trivial));
    assert(trivial_generator.equal_si(-1));

    silex::NumberField quartic_zeta4 = quartic_zeta4_subfield_field();
    silex::Element zeta4(quartic_zeta4);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), zeta4,
                                  quartic_zeta4));
    assert(sflint::fmpz_equal_si(order, 4));
    silex::Element zeta4_squared(quartic_zeta4);
    assert(element_pow_si(zeta4_squared, zeta4, 2));
    assert(zeta4_squared.equal_si(-1));

    silex::NumberField quartic_zeta6 = quartic_zeta6_subfield_field();
    silex::Element zeta6_quartic(quartic_zeta6);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), zeta6_quartic,
                                  quartic_zeta6));
    assert(sflint::fmpz_equal_si(order, 6));
    silex::Element zeta6_quartic_cubed(quartic_zeta6);
    assert(element_pow_si(zeta6_quartic_cubed, zeta6_quartic, 3));
    assert(zeta6_quartic_cubed.equal_si(-1));

    silex::NumberField generic_x2_plus_one = generic_x2_plus_one_field();
    silex::Element generic_generator(generic_x2_plus_one);
    assert(generic_generator.set_si(7));
    sflint::fmpz_set_ui(sflint::FmpzRef(order), 17);
    assert(silex::roots_of_unity(sflint::FmpzRef(order),
                                  generic_generator,
                                  generic_x2_plus_one));
    assert(sflint::fmpz_equal_si(order, 4));
    silex::Element generic_generator_squared(generic_x2_plus_one);
    assert(element_pow_si(generic_generator_squared, generic_generator, 2));
    assert(generic_generator_squared.equal_si(-1));

    silex::Element wrong_parent(degree_one);
    assert(!silex::root_of_unity_generator(wrong_parent, qi));
    assert(wrong_parent.parent() != nullptr &&
           wrong_parent.parent()->has_same_data(degree_one));
    return 0;
}

silex::NumberField integer_polynomial_field(const slong* coefficients,
                                            slong length) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < length; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// True when generator has exact multiplicative order `order`: generator^order
// is one and generator^(order/l) is not one for every prime l | order.
bool has_exact_order(const silex::Element& generator, ulong order) noexcept {
    silex::Element power(*generator.parent());
    sflint::Fmpz exponent;
    fmpz_set_ui(exponent.raw(), order);
    if (!power.pow_fmpz(generator, sflint::FmpzConstRef(exponent)) ||
        !power.equal_si(1)) {
        return false;
    }
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, order, 1);
    for (int i = 0; i < factors.num; ++i) {
        fmpz_set_ui(exponent.raw(), order / factors.p[i]);
        if (!power.pow_fmpz(generator, sflint::FmpzConstRef(exponent)) ||
            power.equal_si(1)) {
            return false;
        }
    }
    return true;
}

void assert_roots_of_unity(const silex::NumberField& field,
                           ulong expected) noexcept {
    sflint::Fmpz order;
    silex::Element generator(field);
    assert(silex::roots_of_unity(sflint::FmpzRef(order), generator, field));
    assert(fmpz_equal_ui(order.raw(), expected) != 0);
    assert(has_exact_order(generator, expected));

    sflint::Fmpz order_only;
    assert(silex::root_of_unity_order(sflint::FmpzRef(order_only), field));
    assert(fmpz_equal_ui(order_only.raw(), expected) != 0);

    silex::Element generator_only(field);
    assert(silex::root_of_unity_generator(generator_only, field));
    assert(has_exact_order(generator_only, expected));
}

// Q(zeta_n) defined by Phi_n for n = 3..30.  The expected w values are from
// an external reference, GP 2.17.4, not from Silex:
//
//     for(n=3,30, print(n, " ", nfrootsof1(nfinit(polcyclo(n)))[1]))
//
// This used to return w = 6 for Q(zeta_9) (true 18) and fail for every n
// with a prime factor >= 5.
int test_roots_of_unity_cyclotomic_fields() {
    const ulong expected_w[31] = {
            0,  0,  0,  6,  4,  10, 6,  14, 8,  18, 10,
            22, 12, 26, 14, 30, 16, 34, 18, 38, 20, 42,
            22, 46, 24, 50, 26, 54, 28, 58, 30};
    for (ulong n = 3; n <= 30; ++n) {
        sflint::FmpzPoly cyclotomic;
        fmpz_poly_cyclotomic(cyclotomic.raw(), n);
        sflint::FmpqPoly polynomial;
        fmpq_poly_set_fmpz_poly(polynomial.raw(), cyclotomic.raw());
        silex::NumberField field = silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(polynomial));
        assert_roots_of_unity(field, expected_w[n]);
    }
    return 0;
}

// Cyclotomic fields given by polynomials that are not translates of a
// cyclotomic polynomial, so the search must build the prime-power roots.
// Each polynomial is charpoly(Mod(x + 2*x^2, polcyclo(m))) in GP 2.17.4, and
// the expected w is GP's nfrootsof1 of it.
int test_roots_of_unity_noncyclotomic_presentations() {
    // Q(zeta_9): w = 18 (3-part 9 found through a cube root of zeta_3).
    const slong zeta9[] = {57, -36, 36, 9, 6, 0, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta9, 7), 18);

    // Q(zeta_16): w = 16 (square roots from zeta_4 up to zeta_16).
    const slong zeta16[] = {257, 16, 80, 128, 32, 0, 0, 0, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta16, 9), 16);

    // Q(zeta_12): w = 12.
    const slong zeta12[] = {13, -8, 11, -4, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta12, 5), 12);

    // Q(zeta_5): w = 10.  This used to fail closed (no root finder
    // for zeta_5 outside cyclotomic-translate presentations); the Hensel
    // search for Phi_5 now finds it.
    const slong zeta5[] = {11, 7, 9, 3, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta5, 5), 10);

    // Q(zeta_7): w = 14, and Q(zeta_11) (degree 10): w = 22.
    const slong zeta7[] = {43, 19, 25, 27, 9, 3, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta7, 7), 14);
    const slong zeta11[] = {683, 235, 137, 251, 377, 243, 81, 27, 9, 3, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta11, 11), 22);
    return 0;
}

// Degree 10 or more, where the exact square and power roots are not
// supported, so the 2- and 3-power roots come from the Hensel search.  The
// polynomials are charpoly(Mod(x + 2*x^2, polcyclo(m))) and
// charpoly(Mod(x + x^2, polcyclo(27))), charpoly(Mod(x + x^3, polcyclo(32)))
// in GP 2.17.4, with w from nfrootsof1.
int test_roots_of_unity_high_degree_prime_powers() {
    // Q(zeta_27), degree 18: w = 54.
    const slong zeta27[] = {261633, -9180, -54756, -118512, -53136, 57024,
                            88704,  69120, 20736,  513,     18,     108,
                            240,    144,   0,      0,       0,      0,
                            1};
    assert_roots_of_unity(integer_polynomial_field(zeta27, 19), 54);
    const slong zeta27b[] = {1, 9, 108, 516, 1278, 1782, 1386, 540, 81, 2,
                             9, 27, 30, 9, 0, 0, 0, 0, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta27b, 19), 54);

    // Q(zeta_32), degree 16: w = 32.
    const slong zeta32[] = {65537, 32, 416, 2816, 10560, 21504, 21504, 8192,
                            512,   0,  0,   0,    0,     0,     0,     0,
                            1};
    assert_roots_of_unity(integer_polynomial_field(zeta32, 17), 32);
    const slong zeta32b[] = {4, 0, -32, 0, 128, 0, 192, 0, 140,
                             0, 16, 0, 0, 0, 0, 0, 1};
    assert_roots_of_unity(integer_polynomial_field(zeta32b, 17), 32);
    return 0;
}

silex::NumberField rational_polynomial_field(const slong* numerators,
                                             const slong* denominators,
                                             slong length) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < length; ++i) {
        fmpq_set_si(coeff.raw(), numerators[i],
                    static_cast<ulong>(denominators[i]));
        fmpq_poly_set_coeff_fmpq(polynomial.raw(), i, coeff.raw());
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// Defining polynomials that are not monic and integral.  Each is a monic
// polynomial above with x replaced by c x (times a constant), and the expected
// w is GP 2.17.4 nfrootsof1(nfinit(P)).
int test_roots_of_unity_nonmonic_presentations() {
    // w = 2: 3x^4 + x + 5 and 2x^4 + 3x^2 + x + 7, both totally complex.
    const slong w2a[] = {5, 1, 0, 0, 3};
    assert_roots_of_unity(integer_polynomial_field(w2a, 5), 2);
    const slong w2b[] = {7, 1, 3, 0, 2};
    assert_roots_of_unity(integer_polynomial_field(w2b, 5), 2);

    // Cyclotomic translate: Phi_5(2x), w = 10.
    const slong phi5_2x[] = {1, 2, 4, 8, 16};
    assert_roots_of_unity(integer_polynomial_field(phi5_2x, 5), 10);

    // w > 2 through the monic model: Q(zeta_5) as f(3x) and f(2x)/16 for
    // f = x^4 + 3x^3 + 9x^2 + 7x + 11 (w = 10), Q(zeta_12) as g(5x) (w = 12)
    // and Q(zeta_9) as h(2x) (w = 18).
    const slong zeta5_3x[] = {11, 21, 81, 81, 81};
    assert_roots_of_unity(integer_polynomial_field(zeta5_3x, 5), 10);
    const slong zeta5_half_num[] = {11, 7, 9, 3, 1};
    const slong zeta5_half_den[] = {16, 8, 4, 2, 1};
    assert_roots_of_unity(
            rational_polynomial_field(zeta5_half_num, zeta5_half_den, 5), 10);
    const slong zeta12_5x[] = {13, -40, 275, -500, 625};
    assert_roots_of_unity(integer_polynomial_field(zeta12_5x, 5), 12);
    const slong zeta9_2x[] = {57, -72, 144, 72, 96, 0, 64};
    assert_roots_of_unity(integer_polynomial_field(zeta9_2x, 7), 18);

    // Degree 16 and 20: P(3x) for P = charpoly(Mod(x + 2*x^2, polcyclo(m)))
    // with m = 32 (w = 32) and m = 25 (w = 50).  The monic model is P itself
    // (theta' = 3 theta), where the Hensel search succeeds; the model
    // theta' = a_n theta has coefficients of several hundred bits, for which
    // the lifting bound cannot be evaluated.
    const slong zeta32_3x[] = {65537,   96,      3744,     76032,   855360,
                               5225472, 15676416, 17915904, 3359232, 0,
                               0,       0,        0,        0,       0,
                               0,       43046721};
    assert_roots_of_unity(integer_polynomial_field(zeta32_3x, 17), 32);
    const slong zeta25_3x[] = {
            1016801,   -924360,    -4705560,   9331200,    17463600,
            -54719469, -1844370,   640441080,  1316136600, 897544800,
            772892361, 1534093020, 754646220,  637729200,  1913187600,
            473513931, 430467210,  2582803260, 0,          0,
            3486784401};
    assert_roots_of_unity(integer_polynomial_field(zeta25_3x, 21), 50);
    return 0;
}

// Q(zeta_40) as T = charpoly(Mod(x + 3*x^2 + x^3, polcyclo(40))) in GP 2.17.4
// (the "40_3" presentation), with w = 40 from nfrootsof1.  The good-prime
// gcd stays at 120 for more than n + 20 = 36 consecutive good primes, where
// the stopping rule of reference `guess_roots` alone stops, and 120 cannot be
// reduced, so this used to fail closed.  phi(120) = 32 does not divide
// n = 16, so the reset of reference `_torsion_group_order_divisor` keeps the
// search going until the gcd reaches 40.  -T(3x) is the same field with a
// non-monic defining polynomial with negative leading coefficient.
int test_roots_of_unity_phi_reset() {
    const slong zeta40[] = {4669921, -6094020, 817562, 2113680, -400340,
                            -802560, 252122,   524760, 183534,  60,
                            -8492,   0,        435,    0,       -22,
                            0,       1};
    assert_roots_of_unity(integer_polynomial_field(zeta40, 17), 40);
    const slong zeta40_m3x[] = {
            -4669921,   18282060,    -7358058,  -57069360, 32427540,
            195022080,  -183796938,  -1147650120, -1204166574, -1180980,
            501444108,  0,           -231176835, 0,           105225318,
            0,          -43046721};
    assert_roots_of_unity(integer_polynomial_field(zeta40_m3x, 17), 40);
    return 0;
}

// Asserts that all three public entry points fail and leave their outputs
// unchanged.
void assert_roots_of_unity_fail_closed(
        const silex::NumberField& field) noexcept {
    sflint::Fmpz order;
    fmpz_set_ui(order.raw(), 17);
    silex::Element generator(field);
    assert(generator.set_si(7));
    assert(!silex::roots_of_unity(sflint::FmpzRef(order), generator, field));
    assert(fmpz_equal_ui(order.raw(), 17) != 0);
    assert(generator.equal_si(7));

    assert(!silex::root_of_unity_order(sflint::FmpzRef(order), field));
    assert(fmpz_equal_ui(order.raw(), 17) != 0);

    assert(!silex::root_of_unity_generator(generator, field));
    assert(generator.equal_si(7));
}

// KNOWN LIMITATION: a field whose good-prime gcd stalls at a strict multiple
// of w that the reset does not catch fails closed through the public API.
// T = charpoly(Mod(x + 2*x^2 + x^3, Q)) in GP 2.17.4 for Q the first
// polcompositum of polcyclo(5) and x^2 + 11, so K = Q(zeta_5, sqrt(-11)) with
// w = 10 (nfrootsof1).  The gcd stays at 30 for more than n + 20 = 28
// consecutive good primes.  phi(30) = 8 divides n = 8, so the reset does not
// apply, and v_3(disc T) = 4 meets the degree and ramification bound, so the
// reduction keeps the factor 3.  K has no zeta_3, so no root of exact order
// 30 is certified.  The threshold 5n = 40 of reference
// `_torsion_group_order_divisor` would reach 10, but Silex keeps n + 20.
// T(3x) is the same field with a non-monic defining polynomial and fails
// closed the same way.  A change of stopping rule may flip these to w = 10,
// never to a wrong w.
int test_roots_of_unity_public_fail_closed() {
    const slong stalled[] = {6813116442581, 1047154163766, 94346195848,
                             4956165560,    147295001,     2623880,
                             31172,         242,           1};
    assert_roots_of_unity_fail_closed(
            integer_polynomial_field(stalled, 9));
    const slong stalled_3x[] = {6813116442581, 3141462491298, 849115762632,
                                133816470120,  11930895081,   637602840,
                                22724388,      529254,        6561};
    assert_roots_of_unity_fail_closed(
            integer_polynomial_field(stalled_3x, 9));
    return 0;
}

// The Hensel search for a root of Phi_m on its own: it finds primitive m-th
// roots of unity that exist and reports failure for those that do not.
int test_cyclotomic_root_hensel() {
    const slong zeta9[] = {57, -36, 36, 9, 6, 0, 1};
    silex::NumberField field9 = integer_polynomial_field(zeta9, 7);
    silex::Element root(field9);
    assert(silex::detail::cyclotomic_root_hensel(root, field9, 9));
    assert(has_exact_order(root, 9));
    assert(silex::detail::cyclotomic_root_hensel(root, field9, 3));
    assert(has_exact_order(root, 3));
    assert(silex::detail::cyclotomic_root_hensel(root, field9, 18));
    assert(has_exact_order(root, 18));

    // Q(zeta_9) has no primitive 4th, 5th or 27th root of unity; the root is
    // left unchanged.
    assert(root.set_si(7));
    assert(!silex::detail::cyclotomic_root_hensel(root, field9, 4));
    assert(!silex::detail::cyclotomic_root_hensel(root, field9, 5));
    assert(!silex::detail::cyclotomic_root_hensel(root, field9, 27));
    assert(root.equal_si(7));

    // A non-monic defining polynomial is outside the helper's contract.
    const slong zeta5_3x[] = {11, 21, 81, 81, 81};
    silex::NumberField nonmonic = integer_polynomial_field(zeta5_3x, 5);
    silex::Element nonmonic_root(nonmonic);
    assert(!silex::detail::cyclotomic_root_hensel(nonmonic_root, nonmonic, 5));
    return 0;
}

// The search publishes w only when it finds a root of unity whose exact order
// is the supplied good-prime bound; a bound it cannot attain fails closed.
int test_roots_of_unity_fail_closed_against_bound() {
    silex::NumberField qi = quadratic_field(-1);
    const slong zeta9[] = {57, -36, 36, 9, 6, 0, 1};
    silex::NumberField field9 = integer_polynomial_field(zeta9, 7);

    sflint::Fmpz order;
    silex::Element generator(qi);
    assert(silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator, qi, 4));
    assert(fmpz_equal_ui(order.raw(), 4) != 0);
    assert(has_exact_order(generator, 4));

    // A bound that is a strict multiple of w: Q(i) has no zeta_3 or zeta_8.
    fmpz_set_ui(order.raw(), 17);
    assert(generator.set_si(7));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator, qi, 12));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator, qi, 8));
    assert(fmpz_equal_ui(order.raw(), 17) != 0);
    assert(generator.equal_si(7));

    silex::Element generator9(field9);
    assert(silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 18));
    assert(fmpz_equal_ui(order.raw(), 18) != 0);
    assert(has_exact_order(generator9, 18));
    // The old zeta_4/zeta_6 search returned 6 here.  A bound of 54 is a multiple
    // of the true w = 18 that the search cannot attain.
    silex::Element saved9(field9);
    assert(saved9.set(generator9));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 54));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 36));
    assert(fmpz_equal_ui(order.raw(), 18) != 0);
    assert(generator9.equal(saved9));

    // A proper divisor of w violates the precondition that the bound is a
    // multiple of w.  Q(zeta_9) contains zeta_6, so without the guard the
    // exact-order check would accept 6 and publish the wrong w; the guard
    // rejects any bound that is not a multiple of the proven good-prime
    // bound.
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 6));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 2));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator9, field9, 0));
    assert(fmpz_equal_ui(order.raw(), 18) != 0);
    assert(generator9.equal(saved9));
    fmpz_set_ui(order.raw(), 17);
    assert(generator.set_si(7));
    assert(!silex::detail::roots_of_unity_for_bound(
            sflint::FmpzRef(order), generator, qi, 2));
    assert(fmpz_equal_ui(order.raw(), 17) != 0);
    assert(generator.equal_si(7));
    return 0;
}

// |value - reference| < 2^-90, with reference a decimal string from GP.
bool arb_near_decimal(const sflint::Arb& value, const char* reference) noexcept {
    sflint::Arb expected;
    sflint::Arb difference;
    sflint::Arf upper;
    if (::arb_set_str(expected.raw(), reference, 256) != 0) {
        return false;
    }
    ::arb_sub(difference.raw(), value.raw(), expected.raw(), 256);
    ::arb_abs(difference.raw(), difference.raw());
    ::arb_get_ubound_arf(upper.raw(), difference.raw(), 256);
    return ::arf_is_finite(upper.raw()) != 0 &&
           ::arf_cmp_2exp_si(upper.raw(), -90) < 0;
}

// The bound is an exact point equal to the lower endpoint of an enclosure
// of ten_thousandths / 10^4: at most that value and within 2^-50 of it.
bool arb_is_lower_point_of(const sflint::Arb& value,
                           ulong ten_thousandths) noexcept {
    sflint::Arb expected;
    sflint::Arb difference;
    ::arb_set_ui(expected.raw(), ten_thousandths);
    ::arb_div_ui(expected.raw(), expected.raw(), 10000, 256);
    ::arb_sub(difference.raw(), expected.raw(), value.raw(), 256);
    sflint::Arf lower;
    sflint::Arf upper;
    ::arb_get_lbound_arf(lower.raw(), difference.raw(), 256);
    ::arb_get_ubound_arf(upper.raw(), difference.raw(), 256);
    return ::arb_is_exact(value.raw()) != 0 &&
           ::arf_sgn(lower.raw()) >= 0 &&
           ::arf_cmp_2exp_si(upper.raw(), -50) < 0;
}

// floor(regulator / bound) for a GP regulator string and an exact bound.
slong floor_regulator_quotient(const char* regulator,
                               const sflint::Arb& bound) noexcept {
    sflint::Arb value;
    sflint::Arf upper;
    sflint::Fmpz result;
    if (::arb_set_str(value.raw(), regulator, 256) != 0) {
        return -1;
    }
    ::arb_div(value.raw(), value.raw(), bound.raw(), 256);
    ::arb_get_ubound_arf(upper.raw(), value.raw(), 256);
    ::arf_get_fmpz(result.raw(), upper.raw(), ARF_RND_FLOOR);
    return ::fmpz_get_si(result.raw());
}

int test_lower_regulator_bound() {
    silex::NumberField real_quadratic = quadratic_field(2);
    sflint::Arb bound;
    assert(silex::unit_lower_regulator_bound(sflint::ArbRef(bound),
                                             real_quadratic, 128));
    assert(sflint::arb_is_positive(bound));
    // Friedman 1989 Table 6, (r1, r2) = (2, 0): R >= 0.48.  With
    // R(Q(sqrt2)) = log(1 + sqrt2) = 0.8813... the index bound is one.
    assert(arb_is_lower_point_of(bound, 4800));

    sflint::Arb sentinel;
    sflint::arb_set_si(sflint::ArbRef(sentinel), 123);
    assert(!silex::unit_lower_regulator_bound(sflint::ArbRef(sentinel),
                                              real_quadratic, 0));
    assert(sflint::arb_contains_si(sentinel, 123));

    // Unit rank zero: still a positive bound (R = 1 by convention), which
    // no unit-index caller consumes.
    silex::NumberField rational = degree_one_field();
    assert(silex::unit_lower_regulator_bound(sflint::ArbRef(bound), rational,
                                             128));
    assert(sflint::arb_is_positive(bound));
    silex::NumberField gaussian = generic_x2_plus_one_field();
    assert(silex::unit_lower_regulator_bound(sflint::ArbRef(bound), gaussian,
                                             128));
    assert(sflint::arb_is_positive(bound));
    return 0;
}

// Each regulator lower bound term against GP 2.17 (\p 40), e.g.
//   zl(r1,r2,g) = log((1+g)*(1+2*g)/2) + (r1+r2)*lngamma(1+g)
//       + r2*lngamma(3/2+g) - (r1+r2)*log(2) - r2*log(Pi)/2
//       - (1+g)*((r1+r2)*psi((1+g)/2) + r2*psi(1+g/2) + 2/g + 1/(1+g));
//   Z(r1,r2,w) = w/50*exp((46*r1+10*r2)/100);
//   F(r1,r2,w) = w*31/10000*exp((241*(r1+2*r2)+497*r1)/1000);
//   S(r1,r2,w,g) = w*exp(zl(r1,r2,g));
int test_lower_regulator_bound_terms() {
    sflint::Fmpz two;
    sflint::Fmpz ten;
    ::fmpz_set_ui(two.raw(), 2);
    ::fmpz_set_ui(ten.raw(), 10);
    const sflint::FmpzConstRef w2(two.raw());
    const sflint::FmpzConstRef w10(ten.raw());
    sflint::Arb value;

    // Friedman 1989 Theorem B.
    assert(silex::detail::regulator_lower_bound_friedman_minimum(
            sflint::ArbRef(value), 128));
    assert(arb_near_decimal(value, "0.2052"));

    // Zimmert 1981 Korollar (i): R >= 0.02 w exp(0.46 r1 + 0.1 r2).
    assert(silex::detail::regulator_lower_bound_zimmert_corollary(
            sflint::ArbRef(value), 3, 0, w2, 128));
    assert(arb_near_decimal(value,
                            "0.1589960650997899247563667112375822200772"));
    assert(silex::detail::regulator_lower_bound_zimmert_corollary(
            sflint::ArbRef(value), 2, 1, w2, 128));
    assert(arb_near_decimal(value,
                            "0.1109277905585719166719679910858177542534"));
    assert(silex::detail::regulator_lower_bound_zimmert_corollary(
            sflint::ArbRef(value), 3, 1, w2, 128));
    assert(arb_near_decimal(value,
                            "0.1757178272367502674293515057326898228977"));
    assert(silex::detail::regulator_lower_bound_zimmert_corollary(
            sflint::ArbRef(value), 0, 2, w10, 128));
    assert(arb_near_decimal(value,
                            "0.2442805516320339667842143989279348340615"));

    // Friedman 1989 p. 620 Corollary: R/w > 0.0031 exp(0.241 n + 0.497 r1).
    assert(silex::detail::regulator_lower_bound_friedman_corollary(
            sflint::ArbRef(value), 3, 0, w2, 128));
    assert(arb_near_decimal(value,
                            "0.05674396414638643833575044923010525274285"));
    assert(silex::detail::regulator_lower_bound_friedman_corollary(
            sflint::ArbRef(value), 2, 1, w2, 128));
    assert(arb_near_decimal(value,
                            "0.04392788412136035000632128871714391496935"));
    assert(silex::detail::regulator_lower_bound_friedman_corollary(
            sflint::ArbRef(value), 8, 0, w2, 128));
    assert(arb_near_decimal(value,
                            "2.272303359991809705171606228068175887917"));
    assert(silex::detail::regulator_lower_bound_friedman_corollary(
            sflint::ArbRef(value), 0, 2, w10, 128));
    assert(arb_near_decimal(value,
                            "0.08128708960401178432231306997820516539408"));

    // Zimmert 1981 Satz 3 at gamma = 1 and gamma = 3/5.
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 3, 0, w2, 1, 1, 128));
    assert(arb_near_decimal(value,
                            "0.1613163261295697741200235837340640423589"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 2, 1, w2, 1, 1, 128));
    assert(arb_near_decimal(value,
                            "0.1124720982727918166884165012387947991104"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 0, 3, w2, 1, 1, 128));
    assert(arb_near_decimal(value,
                            "0.05467360953127803854474959348325015224870"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 0, 2, w10, 1, 1, 128));
    assert(arb_near_decimal(value,
                            "0.2471998004427213873853834918008501855207"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 8, 0, w2, 3, 5, 128));
    assert(arb_near_decimal(value,
                            "2.296361101036429348993838953587253457251"));

    // At gamma = 1, R/w >= 0.0202138 exp(0.461284 r1 + 0.100622 r2); the
    // printed Korollar (i) constants are these rounded down.
    sflint::Arb base;
    sflint::Arb step;
    sflint::Arb ratio;
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(base), 3, 0, w2, 1, 1, 128));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(step), 4, 0, w2, 1, 1, 128));
    ::arb_div(ratio.raw(), step.raw(), base.raw(), 128);
    ::arb_log(ratio.raw(), ratio.raw(), 128);
    assert(arb_near_decimal(ratio,
                            "0.4612841492431204117957920587066282940088"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(step), 3, 1, w2, 1, 1, 128));
    ::arb_div(ratio.raw(), step.raw(), base.raw(), 128);
    ::arb_log(ratio.raw(), ratio.raw(), 128);
    assert(arb_near_decimal(ratio,
                            "0.1006221288341864432385257187103119968916"));

    // Zimmert Tabelle 2 (p. 375) lists lower bounds for 2R/w, rounded down,
    // at his chosen gamma: (3,0) 0.2129 at 1.58, (2,1) 0.1306 at 1.37,
    // (3,1) 0.1809 at 1.09.  With w = 2 these are bounds for R.
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 3, 0, w2, 158, 100, 128));
    assert(arb_near_decimal(value,
                            "0.2129474323813426813568508772937398306256"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 2, 1, w2, 137, 100, 128));
    assert(arb_near_decimal(value,
                            "0.1306738533620172259551575231187187248846"));
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 3, 1, w2, 109, 100, 128));
    assert(arb_near_decimal(value,
                            "0.1809380062002522306708648373948250271538"));

    // Friedman 1989 Table 6 (p. 621).
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 3, 0, 128));
    assert(arb_near_decimal(value, "0.52"));
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 2, 1, 128));
    assert(arb_near_decimal(value, "0.36"));
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 3, 1, 128));
    assert(arb_near_decimal(value, "0.62"));
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 6, 0, 128));
    assert(arb_near_decimal(value, "3.22"));
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 0, 9, 128));
    assert(arb_near_decimal(value, "0.47"));
    // (0,3) is 0.27 except for three sextics; capped at 0.2052.
    assert(silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 0, 3, 128));
    assert(arb_near_decimal(value, "0.2052"));
    // Signatures outside the table, and rank-zero entries, are absent.
    sflint::arb_set_si(sflint::ArbRef(value), 123);
    assert(!silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 8, 0, 128));
    assert(!silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 2, 6, 128));
    assert(!silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 1, 0, 128));
    assert(!silex::detail::regulator_lower_bound_friedman_table6(
            sflint::ArbRef(value), 0, 1, 128));
    assert(sflint::arb_contains_si(value, 123));

    // Invalid inputs.
    assert(!silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(value), 3, 0, w2, 0, 1, 128));
    assert(!silex::detail::regulator_lower_bound_zimmert_corollary(
            sflint::ArbRef(value), 0, 0, w2, 128));
    assert(!silex::detail::regulator_lower_bound_friedman_corollary(
            sflint::ArbRef(value), 3, 0, w2, 0));
    assert(sflint::arb_contains_si(value, 123));
    return 0;
}

int test_lower_regulator_bound_maximum() {
    sflint::Fmpz two;
    ::fmpz_set_ui(two.raw(), 2);
    const sflint::FmpzConstRef w2(two.raw());
    sflint::Arb value;

    // Signatures in Friedman's Table 6: the table entry dominates.
    assert(silex::detail::regulator_lower_bound_from_signature(
            sflint::ArbRef(value), 3, 0, w2, 256));
    assert(arb_is_lower_point_of(value, 5200));
    assert(silex::detail::regulator_lower_bound_from_signature(
            sflint::ArbRef(value), 2, 1, w2, 256));
    assert(arb_is_lower_point_of(value, 3600));
    assert(silex::detail::regulator_lower_bound_from_signature(
            sflint::ArbRef(value), 3, 1, w2, 256));
    assert(arb_is_lower_point_of(value, 6200));
    assert(silex::detail::regulator_lower_bound_from_signature(
            sflint::ArbRef(value), 0, 3, w2, 256));
    assert(arb_is_lower_point_of(value, 2052));

    // (8,0) is not in the table; Satz 3 at gamma = 3/5 is the best of the
    // fixed set (GP: 2.29636...), above Friedman's corollary (2.27230...).
    assert(silex::detail::regulator_lower_bound_from_signature(
            sflint::ArbRef(value), 8, 0, w2, 256));
    sflint::Arb satz3;
    assert(silex::detail::regulator_lower_bound_zimmert_satz3(
            sflint::ArbRef(satz3), 8, 0, w2, 3, 5, 64));
    sflint::Arf expected;
    ::arb_get_lbound_arf(expected.raw(), satz3.raw(), 64);
    assert(::arb_is_exact(value.raw()) != 0);
    assert(::arf_equal(arb_midref(value.raw()), expected.raw()) != 0);

    // Five fields that need the unbounded auxiliary-prime search to prove
    // their unit index (GP 2.17 bnfinit(f, 1); w = 2 and bnfcertify = 1 for
    // each).  Index bound floor(R / R_lower) uses Friedman Table 6:
    // (3,0) 0.52, (2,1) 0.36, (3,1) 0.62.
    struct FriedmanIndexBoundField {
        slong coefficients[6];
        slong length;
        ulong table_entry;
        const char* regulator;
        slong index_bound;
    };
    const FriedmanIndexBoundField fields[] = {
            {{-1, -8, -4, 1}, 4, 5200,
             "12.68082022713452093603368239737775544244", 24},
            {{-3, 2, 1, 4, 1}, 5, 3600,
             "14.46098754389304356185892633589732684899", 40},
            {{3, -3, -3, -4, 1}, 5, 3600, "49.8132070223", 138},
            {{4, 1, 3, -8, 1}, 5, 3600, "28.9150416062", 80},
            {{3, -3, -3, -6, -7, 1}, 6, 6200, "734.654033002", 1184},
    };
    for (const FriedmanIndexBoundField& entry : fields) {
        silex::NumberField field =
                integer_polynomial_field(entry.coefficients, entry.length);
        assert(silex::unit_lower_regulator_bound(sflint::ArbRef(value), field,
                                                 256));
        assert(arb_is_lower_point_of(value, entry.table_entry));
        assert(floor_regulator_quotient(entry.regulator, value) ==
               entry.index_bound);
    }
    return 0;
}

int test_log_matrix_regulator_and_independence() {
    silex::NumberField real_quadratic = quadratic_field(2);
    silex::EmbeddingContext embeddings(real_quadratic);
    assert(embeddings.is_defined());

    silex::Element theta(real_quadratic);
    silex::Element unit(real_quadratic);
    assert(theta.gen());
    assert(unit.add_si(theta, 1));

    std::vector<silex::Element> units;
    units.emplace_back(real_quadratic);
    assert(units[0].set(unit));

    sflint::ArbMat logs(1, 2);
    assert(silex::unit_log_matrix(sflint::ArbMatRef(logs), embeddings,
                                  silex::ElementSpan(units.data(), units.size()),
                                  silex::LogEmbeddingMode::product,
                                  128));
    assert(!sflint::arb_contains_zero(arb_mat_entry(logs.raw(), 0, 0)));
    assert(!sflint::arb_contains_zero(arb_mat_entry(logs.raw(), 0, 1)));

    sflint::ArbVec direct(2);
    assert(silex::logarithmic_embedding(sflint::ArbVecRef(direct),
                                        embeddings, unit,
                                        silex::LogEmbeddingMode::product,
                                        128));
    assert(sflint::arb_overlaps(arb_mat_entry(logs.raw(), 0, 0),
                                sflint::ArbConstRef(direct.data() + 0)));
    assert(sflint::arb_overlaps(arb_mat_entry(logs.raw(), 0, 1),
                                sflint::ArbConstRef(direct.data() + 1)));

    sflint::Arb regulator;
    assert(silex::unit_regulator(sflint::ArbRef(regulator), embeddings,
                                 silex::ElementSpan(units.data(), units.size()),
                                 128));
    assert(sflint::arb_is_positive(regulator));

    bool independent = false;
    silex::EmbeddingContext unset_embeddings;
    assert(silex::units_independent(independent, unset_embeddings,
                                    silex::ElementSpan(), 64));
    assert(independent);

    assert(silex::units_independent(
            independent, embeddings,
            silex::ElementSpan(units.data(), units.size()), 128));
    assert(independent);

    units.emplace_back(real_quadratic);
    assert(units[1].set(unit));
    assert(silex::units_independent(
            independent, embeddings,
            silex::ElementSpan(units.data(), units.size()), 128));
    assert(!independent);

    sflint::ArbMat sentinel(1, 2);
    sflint::arb_set_si(sflint::arb_mat_entry_ref(sentinel, 0, 0), 123);
    sflint::arb_set_si(sflint::arb_mat_entry_ref(sentinel, 0, 1), 456);
    assert(!silex::unit_log_matrix(sflint::ArbMatRef(sentinel), embeddings,
                                   silex::ElementSpan(units.data(), units.size()),
                                   silex::LogEmbeddingMode::product, 0));
    assert(sflint::arb_contains_si(
            sflint::arb_mat_entry_ref(sentinel, 0, 0), 123));
    assert(sflint::arb_contains_si(
            sflint::arb_mat_entry_ref(sentinel, 0, 1), 456));

    return 0;
}

}  // namespace

int main() {
    test_rank();
    test_quadratic_fundamental_unit();
    test_roots_of_unity();
    test_roots_of_unity_cyclotomic_fields();
    test_roots_of_unity_noncyclotomic_presentations();
    test_roots_of_unity_high_degree_prime_powers();
    test_roots_of_unity_nonmonic_presentations();
    test_cyclotomic_root_hensel();
    test_roots_of_unity_fail_closed_against_bound();
    test_roots_of_unity_phi_reset();
    test_roots_of_unity_public_fail_closed();
    test_lower_regulator_bound();
    test_lower_regulator_bound_terms();
    test_lower_regulator_bound_maximum();
    test_log_matrix_regulator_and_independence();
    return 0;
}
