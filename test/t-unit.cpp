#include <silex/unit.hpp>

#include "test_support.hpp"

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

int test_lower_regulator_bound() {
    silex::NumberField real_quadratic = quadratic_field(2);
    sflint::Arb bound;
    assert(silex::unit_lower_regulator_bound(sflint::ArbRef(bound),
                                             real_quadratic, 128));
    assert(sflint::arb_is_positive(bound));

    sflint::Arb sentinel;
    sflint::arb_set_si(sflint::ArbRef(sentinel), 123);
    assert(!silex::unit_lower_regulator_bound(sflint::ArbRef(sentinel),
                                              real_quadratic, 0));
    assert(sflint::arb_contains_si(sentinel, 123));
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
    test_lower_regulator_bound();
    test_log_matrix_regulator_and_independence();
    return 0;
}
