#include "class_group/ideal_lattice_reduction_internal.hpp"
#include "class_group/ideal_lattice_lll_internal.hpp"
#include "class_group/ideal_t2_enumeration_internal.hpp"
#include "test_support.hpp"

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {
namespace sflint = silex::flint;

constexpr slong kReductionPrecision = 256;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
};

FieldSetup cubic_x3_minus_2() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

bool degree_one_prime_ideal_above_five(silex::Ideal& out,
                                       const silex::Order& order) noexcept {
    sflint::Fmpz five;
    sflint::fmpz_set_ui(sflint::FmpzRef(five), UWORD(5));
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order,
                                sflint::FmpzConstRef(five))) {
        return false;
    }

    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime != nullptr && prime->ramification_index() == 1 &&
            prime->residue_degree() == 1) {
            return prime->get_ideal(out);
        }
    }
    return false;
}

bool build_large_ideal(silex::Ideal& out,
                       const silex::Ideal& prime_ideal,
                       const silex::Order& order,
                       const silex::NumberField& field) noexcept {
    sflint::Fmpz scalar;
    sflint::fmpz_one(sflint::FmpzRef(scalar));
    ::fmpz_mul_2exp(scalar.raw(), scalar.raw(), 102);

    silex::Element scalar_element(field);
    silex::OrderElement scalar_generator(order);
    silex::Ideal scalar_ideal(order);
    return scalar_element.is_defined() && scalar_generator.is_defined() &&
           scalar_ideal.is_defined() &&
           scalar_element.set_fmpz(sflint::FmpzConstRef(scalar)) &&
           scalar_generator.set_element(scalar_element) &&
           scalar_ideal.set_principal(scalar_generator) &&
           out.multiply(scalar_ideal, prime_ideal);
}

bool reduced_ideal_identity_holds(const silex::Ideal& input,
                                  const silex::Ideal& reduced,
                                  const silex::Element& multiplier) noexcept {
    const silex::Order* order = input.parent();
    if (order == nullptr) {
        return false;
    }

    silex::FractionalIdeal principal(*order);
    silex::FractionalIdeal input_fractional(*order);
    silex::FractionalIdeal actual(*order);
    silex::FractionalIdeal expected(*order);
    return principal.is_defined() && input_fractional.is_defined() &&
           actual.is_defined() && expected.is_defined() &&
           principal.set_principal(multiplier) &&
           input_fractional.set_integral(input) &&
           actual.multiply(principal, input_fractional) &&
           expected.set_integral(reduced) && actual.equal(expected);
}

bool reduced_product_identity_holds(const silex::Ideal& left,
                                    const silex::Ideal& right,
                                    const silex::Ideal& reduced,
                                    const silex::Element& multiplier) noexcept {
    const silex::Order* order = left.parent();
    if (order == nullptr) {
        return false;
    }

    silex::FractionalIdeal principal(*order);
    silex::FractionalIdeal left_fractional(*order);
    silex::FractionalIdeal right_fractional(*order);
    silex::FractionalIdeal product(*order);
    silex::FractionalIdeal actual(*order);
    silex::FractionalIdeal expected(*order);
    return principal.is_defined() && left_fractional.is_defined() &&
           right_fractional.is_defined() && product.is_defined() &&
           actual.is_defined() && expected.is_defined() &&
           principal.set_principal(multiplier) &&
           left_fractional.set_integral(left) &&
           right_fractional.set_integral(right) &&
           product.multiply(left_fractional, right_fractional) &&
           actual.multiply(principal, product) &&
           expected.set_integral(reduced) && actual.equal(expected);
}

bool reduced_power_identity_holds(const silex::Ideal& input,
                                  slong exponent,
                                  const silex::Ideal& reduced,
                                  const silex::FactoredElement& multiplier)
        noexcept {
    const silex::Order* order = input.parent();
    const silex::NumberField* field =
            order == nullptr ? nullptr : order->parent();
    if (field == nullptr) {
        return false;
    }

    sflint::Fmpz exponent_value;
    sflint::fmpz_set_si(sflint::FmpzRef(exponent_value), exponent);
    silex::Element multiplier_value(*field);
    silex::FractionalIdeal input_fractional(*order);
    silex::FractionalIdeal expected(*order);
    silex::FractionalIdeal reduced_fractional(*order);
    silex::FractionalIdeal principal(*order);
    silex::FractionalIdeal actual(*order);
    return multiplier_value.is_defined() && input_fractional.is_defined() &&
           expected.is_defined() && reduced_fractional.is_defined() &&
           principal.is_defined() && actual.is_defined() &&
           input_fractional.set_integral(input) &&
           expected.pow_fmpz(input_fractional,
                             sflint::FmpzConstRef(exponent_value)) &&
           multiplier.evaluate(multiplier_value) &&
           reduced_fractional.set_integral(reduced) &&
           principal.set_principal(multiplier_value) &&
           actual.multiply(reduced_fractional, principal) &&
           actual.equal(expected);
}

bool integral_ideal_contains_element(const silex::Ideal& ideal,
                                     const silex::Element& element) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr) {
        return false;
    }

    silex::OrderElement order_element(*order);
    return order_element.is_defined() && order_element.set_element(element) &&
           ideal.contains(order_element);
}

bool reduction_matrices_have_expected_properties(
        const silex::Ideal& ideal,
        const sflint::FmpzMat& reduced_basis,
        const sflint::FmpzMat& transform,
        const sflint::FmpzMat& scaled_gram) noexcept {
    const slong degree = ideal.degree();
    if (degree <= 0 || sflint::fmpz_mat_nrows(reduced_basis) != degree ||
        sflint::fmpz_mat_ncols(reduced_basis) != degree ||
        sflint::fmpz_mat_nrows(transform) != degree ||
        sflint::fmpz_mat_ncols(transform) != degree ||
        sflint::fmpz_mat_nrows(scaled_gram) != degree ||
        sflint::fmpz_mat_ncols(scaled_gram) != degree) {
        return false;
    }

    sflint::FmpzMat original_basis(degree, degree);
    sflint::FmpzMat expected_reduced_basis(degree, degree);
    if (!ideal.get_hnf(sflint::FmpzMatRef(original_basis))) {
        return false;
    }
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(expected_reduced_basis),
                         sflint::FmpzMatConstRef(transform),
                         sflint::FmpzMatConstRef(original_basis));
    if (!sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(expected_reduced_basis),
                sflint::FmpzMatConstRef(reduced_basis))) {
        return false;
    }

    sflint::Fmpz determinant;
    sflint::fmpz_mat_det(sflint::FmpzRef(determinant),
                         sflint::FmpzMatConstRef(transform));
    if (!sflint::fmpz_is_pm1(sflint::FmpzConstRef(determinant))) {
        return false;
    }
    for (slong i = 0; i < degree; ++i) {
        if (sflint::fmpz_sgn(sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(scaled_gram), i, i)) <= 0) {
            return false;
        }
    }
    return true;
}

void reset_element_sentinels(silex::Ideal& reduced,
                             silex::Element& multiplier,
                             const silex::Ideal& reduced_sentinel,
                             const silex::Element& multiplier_sentinel) {
    assert(reduced.set(reduced_sentinel));
    assert(multiplier.set(multiplier_sentinel));
}

int test_ideal_lattice_result_contracts() {
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const slong degree = order.degree();

    silex::Ideal prime_ideal(order);
    assert(prime_ideal.is_defined());
    assert(degree_one_prime_ideal_above_five(prime_ideal, order));

    silex::detail::relation_search::IdealLatticeLllData general_data;
    assert(silex::detail::relation_search::build_ideal_lattice_lll_data(
            general_data, prime_ideal, kReductionPrecision));
    assert(sflint::fmpz_mat_nrows(general_data.basis) == degree);
    assert(sflint::fmpz_mat_ncols(general_data.basis) == degree);
    assert(sflint::fmpz_mat_nrows(general_data.scaled_gram) == degree);
    assert(sflint::fmpz_mat_ncols(general_data.scaled_gram) == degree);
    assert(sflint::fmpz_sgn(
                   sflint::FmpzConstRef(general_data.gram_denominator)) > 0);

    silex::detail::relation_search::FiniteIdealT2EnumerationData finite_data;
    assert(silex::detail::relation_search::
                   build_finite_ideal_t2_enumeration_data_with_retry(
                           finite_data, prime_ideal, nullptr, nullptr));
    assert(sflint::fmpz_mat_nrows(finite_data.basis) == degree);
    assert(sflint::fmpz_mat_ncols(finite_data.basis) == degree);
    assert(finite_data.quadratic_form_data.size() ==
           static_cast<std::size_t>(degree * degree));
    assert(std::isfinite(finite_data.initial_bound_value));
    assert(finite_data.initial_bound_value > 0.0);

    // These are mutation sentinels, not golden mathematical outputs.
    general_data.basis = sflint::FmpzMat(1, 2);
    general_data.scaled_gram = sflint::FmpzMat(2, 1);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(general_data.basis, 0, 0), 11);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(general_data.basis, 0, 1), 13);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(general_data.scaled_gram, 0, 0), 17);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(general_data.scaled_gram, 1, 0), 19);
    sflint::fmpz_set_si(sflint::FmpzRef(general_data.gram_denominator), 23);
    assert(!silex::detail::relation_search::build_ideal_lattice_lll_data(
            general_data, prime_ideal, 0));
    assert(sflint::fmpz_mat_nrows(general_data.basis) == 1);
    assert(sflint::fmpz_mat_ncols(general_data.basis) == 2);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(general_data.basis), 0, 0),
            11));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(general_data.basis), 0, 1),
            13));
    assert(sflint::fmpz_mat_nrows(general_data.scaled_gram) == 2);
    assert(sflint::fmpz_mat_ncols(general_data.scaled_gram) == 1);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(general_data.scaled_gram), 0, 0),
            17));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(general_data.scaled_gram), 1, 0),
            19));
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(general_data.gram_denominator), 23));

    finite_data.basis = sflint::FmpzMat(1, 2);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(finite_data.basis, 0, 0), 29);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(finite_data.basis, 0, 1), 31);
    finite_data.quadratic_form_data = {37.0, 41.0};
    finite_data.initial_bound_value = 43.0;
    silex::Ideal no_hnf_ideal(order);
    assert(no_hnf_ideal.is_defined());
    assert(!no_hnf_ideal.has_hnf());
    assert(!silex::detail::relation_search::
                    build_finite_ideal_t2_enumeration_data_with_retry(
                            finite_data, no_hnf_ideal, nullptr, nullptr));
    assert(sflint::fmpz_mat_nrows(finite_data.basis) == 1);
    assert(sflint::fmpz_mat_ncols(finite_data.basis) == 2);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(finite_data.basis), 0, 0),
            29));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(finite_data.basis), 0, 1),
            31));
    assert(finite_data.quadratic_form_data.size() == 2);
    assert(finite_data.quadratic_form_data[0] == 37.0);
    assert(finite_data.quadratic_form_data[1] == 41.0);
    assert(finite_data.initial_bound_value == 43.0);
    return 0;
}

int test_reduction_identities_and_mutation() {
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const silex::NumberField& field = setup.field;

    silex::Ideal prime_ideal(order);
    silex::Ideal large_ideal(order);
    assert(prime_ideal.is_defined());
    assert(large_ideal.is_defined());
    assert(degree_one_prime_ideal_above_five(prime_ideal, order));
    assert(build_large_ideal(large_ideal, prime_ideal, order, field));

    silex::Ideal reduced_large(order);
    silex::Element large_multiplier(field);
    assert(reduced_large.is_defined());
    assert(large_multiplier.is_defined());
    assert(silex::detail::reduce_ideal_lattice(
            reduced_large, large_multiplier, large_ideal,
            kReductionPrecision));
    // These nonidentity checks only prove that the fixture reaches the
    // reduction path; the reconstruction identity below is the contract.
    assert(!reduced_large.equal(large_ideal));
    assert(!large_multiplier.equal_si(1));
    assert(reduced_ideal_identity_holds(
            large_ideal, reduced_large, large_multiplier));

    silex::Ideal reduced_product(order);
    silex::Element product_multiplier(field);
    assert(reduced_product.is_defined());
    assert(product_multiplier.is_defined());
    assert(silex::detail::reduce_ideal_product(
            reduced_product, product_multiplier, prime_ideal, prime_ideal,
            kReductionPrecision));
    assert(reduced_product_identity_holds(
            prime_ideal, prime_ideal, reduced_product, product_multiplier));

    constexpr std::array<slong, 5> exponents{-1, 0, 1, 2, 3};
    for (slong exponent : exponents) {
        silex::Ideal reduced_power(order);
        silex::FactoredElement power_multiplier(field);
        assert(reduced_power.is_defined());
        assert(power_multiplier.is_defined());
        assert(silex::detail::reduce_ideal_signed_power(
                reduced_power, power_multiplier, prime_ideal, exponent,
                kReductionPrecision));
        assert(reduced_power_identity_holds(
                prime_ideal, exponent, reduced_power, power_multiplier));
    }
    return 0;
}

int test_short_elements_and_reduction_matrices() {
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const silex::NumberField& field = setup.field;
    const slong degree = order.degree();

    silex::Ideal prime_ideal(order);
    assert(prime_ideal.is_defined());
    assert(degree_one_prime_ideal_above_five(prime_ideal, order));

    sflint::FmpzMat zero_weights(1, degree);
    sflint::FmpzMat nonzero_weights(1, degree);
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(zero_weights));
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(nonzero_weights));
    sflint::fmpz_one(sflint::fmpz_mat_entry(nonzero_weights, 0, 0));

    silex::Element unweighted_short(field);
    silex::Element zero_weighted_short(field);
    silex::Element weighted_short(field);
    assert(unweighted_short.is_defined());
    assert(zero_weighted_short.is_defined());
    assert(weighted_short.is_defined());
    assert(silex::detail::ideal_lattice_short_element(
            unweighted_short, prime_ideal, kReductionPrecision));
    assert(integral_ideal_contains_element(prime_ideal, unweighted_short));
    assert(silex::detail::weighted_ideal_lattice_short_element(
            zero_weighted_short, prime_ideal,
            sflint::FmpzMatConstRef(zero_weights), kReductionPrecision));
    assert(zero_weighted_short.equal(unweighted_short));
    assert(silex::detail::weighted_ideal_lattice_short_element(
            weighted_short, prime_ideal,
            sflint::FmpzMatConstRef(nonzero_weights), kReductionPrecision));
    assert(integral_ideal_contains_element(prime_ideal, weighted_short));

    sflint::Fmpz two;
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    silex::FractionalIdeal fractional_ideal(order);
    assert(fractional_ideal.is_defined());
    assert(fractional_ideal.set_integral_den(
            prime_ideal, sflint::FmpzConstRef(two)));

    silex::Ideal fractional_numerator(order);
    sflint::Fmpz fractional_denominator;
    assert(fractional_numerator.is_defined());
    assert(fractional_ideal.get_integral_den(
            fractional_numerator,
            sflint::FmpzRef(fractional_denominator)));
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(fractional_denominator), 2));

    silex::Element fractional_numerator_short(field);
    silex::Element fractional_short(field);
    silex::Element expected_fractional_short(field);
    assert(fractional_numerator_short.is_defined());
    assert(fractional_short.is_defined());
    assert(expected_fractional_short.is_defined());
    assert(silex::detail::weighted_ideal_lattice_short_element(
            fractional_numerator_short, fractional_numerator,
            sflint::FmpzMatConstRef(nonzero_weights), kReductionPrecision));
    assert(silex::detail::weighted_ideal_lattice_short_element(
            fractional_short, fractional_ideal,
            sflint::FmpzMatConstRef(nonzero_weights), kReductionPrecision));
    assert(fractional_ideal.contains(fractional_short));
    assert(expected_fractional_short.scalar_div_si(
            fractional_numerator_short, 2));
    assert(fractional_short.equal(expected_fractional_short));

    sflint::FmpzMat zero_basis(degree, degree);
    sflint::FmpzMat zero_transform(degree, degree);
    sflint::FmpzMat zero_scaled_gram(degree, degree);
    sflint::Fmpz zero_gram_denominator;
    assert(silex::detail::relation_search::
                   build_weighted_ideal_lattice_reduction(
                           zero_basis, zero_transform, zero_scaled_gram,
                           zero_gram_denominator, prime_ideal,
                           sflint::FmpzMatConstRef(zero_weights),
                           kReductionPrecision, nullptr, nullptr));
    assert(reduction_matrices_have_expected_properties(
            prime_ideal, zero_basis, zero_transform, zero_scaled_gram));

    sflint::Fmpz expected_zero_denominator;
    sflint::fmpz_one(sflint::FmpzRef(expected_zero_denominator));
    sflint::fmpz_mul_2exp(
            sflint::FmpzRef(expected_zero_denominator),
            sflint::FmpzConstRef(expected_zero_denominator),
            static_cast<ulong>(kReductionPrecision / 2));
    assert(sflint::fmpz_equal(
            sflint::FmpzConstRef(zero_gram_denominator),
            sflint::FmpzConstRef(expected_zero_denominator)));

    sflint::FmpzMat weighted_basis(degree, degree);
    sflint::FmpzMat weighted_transform(degree, degree);
    sflint::FmpzMat weighted_scaled_gram(degree, degree);
    sflint::Fmpz weighted_gram_denominator;
    assert(silex::detail::relation_search::
                   build_weighted_ideal_lattice_reduction(
                           weighted_basis, weighted_transform,
                           weighted_scaled_gram, weighted_gram_denominator,
                           prime_ideal,
                           sflint::FmpzMatConstRef(nonzero_weights),
                           kReductionPrecision, nullptr, nullptr));
    assert(reduction_matrices_have_expected_properties(
            prime_ideal, weighted_basis, weighted_transform,
            weighted_scaled_gram));
    assert(sflint::fmpz_sgn(
                   sflint::FmpzConstRef(weighted_gram_denominator)) > 0);
    return 0;
}

int test_deterministic_failures_preserve_outputs() {
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const silex::NumberField& field = setup.field;

    silex::Ideal prime_ideal(order);
    assert(prime_ideal.is_defined());
    assert(degree_one_prime_ideal_above_five(prime_ideal, order));

    silex::Element seven(field);
    silex::Ideal reduced(order);
    silex::Element multiplier(field);
    assert(seven.is_defined());
    assert(reduced.is_defined());
    assert(multiplier.is_defined());
    assert(seven.set_si(7));

    reset_element_sentinels(reduced, multiplier, prime_ideal, seven);
    assert(!silex::detail::reduce_ideal_lattice(
            reduced, multiplier, prime_ideal, 0));
    assert(reduced.equal(prime_ideal));
    assert(multiplier.equal(seven));

    const silex::NumberField other_field = silex::test::quadratic_field(3);
    const silex::Order other_order = silex::test::equation_order(other_field);
    silex::Ideal other_ideal(other_order);
    assert(other_ideal.is_defined());
    assert(other_ideal.one());

    reset_element_sentinels(reduced, multiplier, prime_ideal, seven);
    assert(!silex::detail::reduce_ideal_product(
            reduced, multiplier, prime_ideal, other_ideal,
            kReductionPrecision));
    assert(reduced.equal(prime_ideal));
    assert(multiplier.equal(seven));

    assert(reduced.set(prime_ideal));
    silex::FactoredElement factored_multiplier(field);
    assert(factored_multiplier.is_defined());
    assert(factored_multiplier.set_element(seven));
    assert(!silex::detail::reduce_ideal_signed_power(
            reduced, factored_multiplier, prime_ideal,
            std::numeric_limits<slong>::min(), kReductionPrecision));
    assert(reduced.equal(prime_ideal));
    silex::Element evaluated(field);
    assert(evaluated.is_defined());
    assert(factored_multiplier.evaluate(evaluated));
    assert(evaluated.equal(seven));

    sflint::FmpzMat valid_weights(1, order.degree());
    sflint::FmpzMat invalid_weights(1, order.degree() - 1);
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(valid_weights));
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(invalid_weights));

    silex::Element short_element(field);
    assert(short_element.is_defined());
    assert(short_element.set(seven));
    assert(!silex::detail::ideal_lattice_short_element(
            short_element, prime_ideal, 0));
    assert(short_element.equal(seven));

    assert(short_element.set(seven));
    assert(!silex::detail::weighted_ideal_lattice_short_element(
            short_element, prime_ideal,
            sflint::FmpzMatConstRef(valid_weights), 0));
    assert(short_element.equal(seven));

    assert(short_element.set(seven));
    assert(!silex::detail::weighted_ideal_lattice_short_element(
            short_element, prime_ideal,
            sflint::FmpzMatConstRef(invalid_weights), kReductionPrecision));
    assert(short_element.equal(seven));

    silex::FractionalIdeal fractional_ideal(order);
    assert(fractional_ideal.is_defined());
    assert(fractional_ideal.set_integral(prime_ideal));

    assert(short_element.set(seven));
    assert(!silex::detail::weighted_ideal_lattice_short_element(
            short_element, fractional_ideal,
            sflint::FmpzMatConstRef(valid_weights), 0));
    assert(short_element.equal(seven));

    assert(short_element.set(seven));
    assert(!silex::detail::weighted_ideal_lattice_short_element(
            short_element, fractional_ideal,
            sflint::FmpzMatConstRef(invalid_weights), kReductionPrecision));
    assert(short_element.equal(seven));
    return 0;
}

}  // namespace

int main() {
    assert(test_ideal_lattice_result_contracts() == 0);
    assert(test_reduction_identities_and_mutation() == 0);
    assert(test_short_elements_and_reduction_matrices() == 0);
    assert(test_deterministic_failures_preserve_outputs() == 0);
    return 0;
}
