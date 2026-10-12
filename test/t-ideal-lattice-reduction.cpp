#include "class_group/class_group_storage_internal.hpp"
#include "class_group/ideal_lattice_reduction_internal.hpp"
#include "class_group/ideal_lattice_lll_internal.hpp"
#include "class_group/ideal_t2_enumeration_internal.hpp"
#include "test_support.hpp"

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstring>
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

FieldSetup quadratic_x2_minus_6() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -6);

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

FieldSetup cubic_x3_plus_x_plus_200() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 200);

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

// The real subfield of Q(zeta_136), degree 32, by the minimal polynomial of
// 2 cos(2 pi / 136) (benchmark row cyclotomic_real_136_deg32).  Its power
// basis embeddings span about 2^0 to 2^31.
FieldSetup real_cyclotomic_136() noexcept {
    constexpr std::array<slong, 33> coefficients = {
            1,       0, -144,    0, 3396,    0, -31248, 0, 149106,
            0, -425128, 0,  786798,  0, -995824, 0, 891311, 0,
            -575832, 0, 271214, 0, -93128, 0, 23051, 0, -4004, 0,
            463,     0, -32,    0, 1};
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }

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

// The field of the monic polynomial with low-to-high coefficients
// `coefficients` (the leading 1 included), with its maximal order.
FieldSetup monic_field(const slong* coefficients, slong count) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < count; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }

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

// The degree-one unramified prime ideal above p with the given index among
// those in decompose_prime order.
bool degree_one_prime_ideal_above(silex::Ideal& out,
                                  const silex::Order& order, ulong p,
                                  slong index = 0) noexcept {
    sflint::Fmpz prime_number;
    sflint::fmpz_set_ui(sflint::FmpzRef(prime_number), p);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order,
                                sflint::FmpzConstRef(prime_number))) {
        return false;
    }

    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime != nullptr && prime->ramification_index() == 1 &&
            prime->residue_degree() == 1) {
            if (index == 0) {
                return prime->get_ideal(out);
            }
            --index;
        }
    }
    return false;
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

bool same_lattice_as_ideal(sflint::FmpzMatConstRef basis,
                           const silex::Ideal& ideal) noexcept {
    const slong degree = sflint::fmpz_mat_nrows(basis);
    if (degree <= 0 || sflint::fmpz_mat_ncols(basis) != degree) {
        return false;
    }
    sflint::FmpzMat ideal_hnf(degree, degree);
    sflint::FmpzMat basis_hnf(degree, degree);
    sflint::FmpzMat expected_hnf(degree, degree);
    if (!ideal.get_hnf(sflint::FmpzMatRef(ideal_hnf))) {
        return false;
    }
    ::fmpz_mat_hnf(basis_hnf.raw(), basis.raw());
    ::fmpz_mat_hnf(expected_hnf.raw(), ideal_hnf.raw());
    return sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(basis_hnf),
                                  sflint::FmpzMatConstRef(expected_hnf));
}

bool rounded_embedding_data_has_contract(
        const silex::detail::relation_search::IdealLatticeLllData& data,
        const silex::Ideal& ideal) noexcept {
    const slong degree = sflint::fmpz_mat_nrows(data.basis);
    if (degree <= 0 || sflint::fmpz_mat_ncols(data.basis) != degree ||
        sflint::fmpz_mat_nrows(data.scaled_gram) != degree ||
        sflint::fmpz_mat_ncols(data.scaled_gram) != degree ||
        sflint::fmpz_sgn(sflint::FmpzConstRef(data.gram_denominator)) <= 0 ||
        !same_lattice_as_ideal(sflint::FmpzMatConstRef(data.basis), ideal)) {
        return false;
    }
    for (slong i = 0; i < degree; ++i) {
        if (sflint::fmpz_sgn(sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(data.scaled_gram), i, i)) <= 0) {
            return false;
        }
        for (slong j = 0; j < degree; ++j) {
            if (!sflint::fmpz_equal(
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(data.scaled_gram), i,
                                j),
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(data.scaled_gram), j,
                                i))) {
                return false;
            }
        }
    }
    return true;
}

// The basis of the one-LLL route computed directly from the cached rounded
// rows R: U H with U the LLL transform of H R (delta 0.99) for the ideal's
// HNF basis H.
bool one_lll_basis_from_cache(
        sflint::FmpzMat& out,
        const silex::Ideal& ideal,
        const silex::detail::relation_search::IdealLatticeRoundedEmbeddingCache&
                cache) noexcept {
    const slong degree = sflint::fmpz_mat_nrows(cache.rounded_rows);
    if (degree <= 0) {
        return false;
    }
    sflint::FmpzMat hnf(degree, degree);
    if (!ideal.get_hnf(sflint::FmpzMatRef(hnf))) {
        return false;
    }
    sflint::FmpzMat embedded(degree, degree);
    sflint::FmpzMat transform(degree, degree);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(embedded),
                         sflint::FmpzMatConstRef(hnf),
                         sflint::FmpzMatConstRef(cache.rounded_rows));
    sflint::fmpz_mat_one(sflint::FmpzMatRef(transform));
    sflint::FmpzLll config(0.99, 0.51, Z_BASIS, APPROX);
    ::fmpz_lll(embedded.raw(), transform.raw(), config.raw());
    out = sflint::FmpzMat(degree, degree);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(out),
                         sflint::FmpzMatConstRef(transform),
                         sflint::FmpzMatConstRef(hnf));
    return true;
}

bool same_lll_data(
        const silex::detail::relation_search::IdealLatticeLllData& left,
        const silex::detail::relation_search::IdealLatticeLllData&
                right) noexcept {
    return sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left.basis),
                                  sflint::FmpzMatConstRef(right.basis)) &&
           sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left.scaled_gram),
                                  sflint::FmpzMatConstRef(right.scaled_gram)) &&
           sflint::fmpz_equal(sflint::FmpzConstRef(left.gram_denominator),
                              sflint::FmpzConstRef(right.gram_denominator));
}

int test_rounded_embedding_lll_contracts() {
    namespace rs = silex::detail::relation_search;
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const slong degree = order.degree();

    silex::Ideal prime_ideal(order);
    silex::Ideal large_ideal(order);
    assert(prime_ideal.is_defined() && large_ideal.is_defined());
    assert(degree_one_prime_ideal_above_five(prime_ideal, order));
    assert(build_large_ideal(large_ideal, prime_ideal, order, setup.field));

    rs::IdealLatticeRoundedEmbeddingCache cache;
    rs::IdealLatticeLllData prime_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            prime_data, prime_ideal, kReductionPrecision, cache));
    assert(rounded_embedding_data_has_contract(prime_data, prime_ideal));
    assert(cache.has_order && cache.has_rounded_rows);
    assert(sflint::fmpz_mat_nrows(cache.rounded_rows) == degree);
    assert(sflint::fmpz_mat_rank(sflint::FmpzMatConstRef(cache.rounded_rows)) ==
           degree);
    assert(cache.order_gram_precision >= kReductionPrecision);
    // The first basis passed the route check against the two-LLL basis.
    assert(cache.route_checked && !cache.use_two_lll_route);
    // The basis is the one-LLL route's.
    sflint::FmpzMat one_lll_basis(0, 0);
    assert(one_lll_basis_from_cache(one_lll_basis, prime_ideal, cache));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(one_lll_basis),
                                  sflint::FmpzMatConstRef(prime_data.basis)));

    // The cached per-order rows are reused unchanged for another ideal.
    sflint::FmpzMat cached_rows(degree, degree);
    sflint::fmpz_mat_set(sflint::FmpzMatRef(cached_rows),
                         sflint::FmpzMatConstRef(cache.rounded_rows));
    rs::IdealLatticeLllData large_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            large_data, large_ideal, kReductionPrecision, cache));
    assert(rounded_embedding_data_has_contract(large_data, large_ideal));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(cached_rows),
                                  sflint::FmpzMatConstRef(cache.rounded_rows)));

    // On x^3 + x + 200 the one-LLL basis of P^4, P above 5 of degree one,
    // differs from the two-LLL route's, so this shows the one-LLL result
    // was accepted.
    {
        FieldSetup skew = cubic_x3_plus_x_plus_200();
        silex::Ideal skew_prime(skew.maximal_order);
        silex::Ideal skew_square(skew.maximal_order);
        silex::Ideal skew_fourth(skew.maximal_order);
        assert(degree_one_prime_ideal_above_five(skew_prime,
                                                 skew.maximal_order));
        assert(skew_square.multiply(skew_prime, skew_prime));
        assert(skew_fourth.multiply(skew_square, skew_square));
        rs::IdealLatticeRoundedEmbeddingCache skew_cache;
        rs::IdealLatticeLllData skew_data;
        assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
                skew_data, skew_fourth, kReductionPrecision, skew_cache));
        assert(rounded_embedding_data_has_contract(skew_data, skew_fourth));
        assert(one_lll_basis_from_cache(one_lll_basis, skew_fourth,
                                        skew_cache));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(one_lll_basis),
                sflint::FmpzMatConstRef(skew_data.basis)));
    }

    // A cache filled for another order is rebuilt for the new order.
    FieldSetup quadratic = quadratic_x2_minus_6();
    silex::Ideal quadratic_prime(quadratic.maximal_order);
    assert(quadratic_prime.is_defined());
    assert(degree_one_prime_ideal_above_five(quadratic_prime,
                                             quadratic.maximal_order));
    rs::IdealLatticeLllData quadratic_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            quadratic_data, quadratic_prime, kReductionPrecision, cache));
    assert(rounded_embedding_data_has_contract(quadratic_data,
                                               quadratic_prime));
    assert(sflint::fmpz_mat_nrows(cache.rounded_rows) == 2);

    // Rounded rows that steer the LLL badly: with R = [[1, 0, 0], [M, 1, 0],
    // [0, 0, 1]] and M = 2^200, every basis that is short in the embedding
    // B R has a vector with first coordinate near -M times its second, far
    // above the T2 bounds.  The checks reject it and the result is the
    // two-LLL route's.
    rs::IdealLatticeRoundedEmbeddingCache steered_cache;
    rs::IdealLatticeLllData steered_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            steered_data, large_ideal, kReductionPrecision, steered_cache));
    sflint::fmpz_mat_one(sflint::FmpzMatRef(steered_cache.rounded_rows));
    sflint::fmpz_one(sflint::fmpz_mat_entry(steered_cache.rounded_rows, 1, 0));
    sflint::fmpz_mul_2exp(
            sflint::fmpz_mat_entry(steered_cache.rounded_rows, 1, 0),
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(steered_cache.rounded_rows), 1, 0),
            200);
    rs::IdealLatticeLllData fallback_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            steered_data, large_ideal, kReductionPrecision, steered_cache));
    assert(rs::build_ideal_lattice_lll_data(fallback_data, large_ideal,
                                            kReductionPrecision));
    assert(rounded_embedding_data_has_contract(steered_data, large_ideal));
    assert(same_lll_data(steered_data, fallback_data));

    // Deterministic failures leave the output untouched.
    prime_data.basis = sflint::FmpzMat(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(prime_data.basis, 0, 0), 11);
    assert(!rs::build_ideal_lattice_rounded_embedding_lll_data(
            prime_data, prime_ideal, 0, cache));
    silex::Ideal no_hnf_ideal(order);
    assert(no_hnf_ideal.is_defined() && !no_hnf_ideal.has_hnf());
    assert(!rs::build_ideal_lattice_rounded_embedding_lll_data(
            prime_data, no_hnf_ideal, kReductionPrecision, cache));
    assert(sflint::fmpz_mat_nrows(prime_data.basis) == 1);
    assert(sflint::fmpz_mat_ncols(prime_data.basis) == 2);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(prime_data.basis), 0, 0),
            11));
    return 0;
}

// On the real subfield of Q(zeta_136) the order rows rounded at the first
// full-rank scale steer the LLL to a basis that passes the T2 bounds but is
// far longer than the two-LLL basis: its T2 diagonal product exceeds the
// two-LLL one by far more than the factor 2^32 the route check allows.  The route check sends the order to the two-LLL
// route, for the checked ideal and every later one, and a cache for another
// order starts unchecked.
int test_rounded_embedding_route_check_for_badly_rounded_order() {
    namespace rs = silex::detail::relation_search;
    FieldSetup setup = real_cyclotomic_136();
    const silex::Order& order = setup.maximal_order;
    assert(order.degree() == 32);

    silex::Ideal first_prime(order);
    silex::Ideal second_prime(order);
    assert(first_prime.is_defined() && second_prime.is_defined());
    assert(degree_one_prime_ideal_above(first_prime, order, UWORD(137), 0));
    assert(degree_one_prime_ideal_above(second_prime, order, UWORD(137), 1));

    rs::IdealLatticeRoundedEmbeddingCache cache;
    rs::IdealLatticeLllData first_data;
    rs::IdealLatticeLllData first_expected;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            first_data, first_prime, kReductionPrecision, cache));
    assert(cache.has_rounded_rows && cache.route_checked &&
           cache.use_two_lll_route);
    assert(rs::build_ideal_lattice_lll_data(first_expected, first_prime,
                                            kReductionPrecision));
    assert(rounded_embedding_data_has_contract(first_data, first_prime));
    assert(same_lll_data(first_data, first_expected));
    sflint::FmpzMat one_lll_basis(0, 0);
    assert(one_lll_basis_from_cache(one_lll_basis, first_prime, cache));
    assert(!sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(one_lll_basis),
                                   sflint::FmpzMatConstRef(first_data.basis)));

    rs::IdealLatticeLllData second_data;
    rs::IdealLatticeLllData second_expected;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            second_data, second_prime, kReductionPrecision, cache));
    assert(rs::build_ideal_lattice_lll_data(second_expected, second_prime,
                                            kReductionPrecision));
    assert(rounded_embedding_data_has_contract(second_data, second_prime));
    assert(same_lll_data(second_data, second_expected));

    FieldSetup cubic = cubic_x3_minus_2();
    silex::Ideal cubic_prime(cubic.maximal_order);
    assert(degree_one_prime_ideal_above_five(cubic_prime,
                                             cubic.maximal_order));
    rs::IdealLatticeLllData cubic_data;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            cubic_data, cubic_prime, kReductionPrecision, cache));
    assert(cache.route_checked && !cache.use_two_lll_route);
    assert(one_lll_basis_from_cache(one_lll_basis, cubic_prime, cache));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(one_lll_basis),
                                  sflint::FmpzMatConstRef(cubic_data.basis)));
    return 0;
}

// Counts the route checks of the rounded-embedding reduction and the ideals
// reduced by the two-LLL route of a switched order.
struct RouteCheckCounter {
    long checks = 0;
    long two_lll_order = 0;
};

void route_check_callback(void* user,
                          silex::DiagnosticsModule,
                          silex::ProfileEvent event,
                          const char*,
                          const char* label) noexcept {
    if (event != silex::ProfileEvent::begin_scope || label == nullptr) {
        return;
    }
    RouteCheckCounter& counter = *static_cast<RouteCheckCounter*>(user);
    if (std::strcmp(label,
                    "class_group.ideal_lattice_reduction.rounded_embedding."
                    "route_check") == 0) {
        ++counter.checks;
    } else if (std::strcmp(label,
                           "class_group.ideal_lattice_reduction."
                           "rounded_embedding.two_lll_order") == 0) {
        ++counter.two_lll_order;
    }
}

// Runs the grh class/unit transaction on the maximal order of `setup` with
// the options of the class/unit benchmark tool, and copies the route-check
// state its class-group context kept.
void run_grh_route_check(bool& route_checked,
                         bool& use_two_lll_route,
                         RouteCheckCounter& counter,
                         const FieldSetup& setup) noexcept {
    const silex::Order& order = setup.maximal_order;
    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), order));
    silex::ClassGroupComputeOptions options;
    const slong n = order.degree();
    options.max_candidates = std::max<slong>(5000, n * n + 64);
    options.max_relations = 1000;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::class_group),
            route_check_callback, &counter);
    options.diagnostics = &diagnostics;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            class_group, order, sflint::FmpzConstRef(factor_base_bound),
            options, 128));
    assert(class_group.certification_status() ==
           silex::CertificationMode::grh);
    const silex::detail::relation_search::IdealLatticeRoundedEmbeddingCache*
            cache = silex::detail::ClassGroupRelationSearchAccess::
                    lll_rounded_embedding_cache(class_group);
    assert(cache != nullptr && cache->has_order &&
           cache->order.has_same_data(order));
    route_checked = cache->route_checked;
    use_two_lll_route = cache->use_two_lll_route;
}

// The route check runs once per order in a class-group context: the result
// for the first accepted ideal of the order is kept in the context and used
// by every later relation-route slice.  On the real subfield of Q(zeta_31)
// (benchmark row cyclotomic_real_31_deg15) the first ideal, the unit ideal
// of the first slice, passes the check (one-LLL T2 diagonal product about
// 2^13.3 times the two-LLL one, under 2^15), so no ideal of the computation
// takes the two-LLL route.  When each slice ran its own check, the first
// pivot ideals of four of the seven extension slices failed it (2^15.2 to
// 2^20.9), so those slices switched route and the relation set changed.  On the octic field of row
// multiquadratic_p5_m7_m11 the unit ideal fails the check (about 2^23.4
// against 2^8), and the order keeps the two-LLL route.
int test_rounded_embedding_route_check_once_per_context() {
    constexpr std::array<slong, 16> real_31 = {
            -1, -8, 28, 84, -126, -252, 210, 330,
            -165, -220, 66, 78, -13, -14, 1, 1};
    bool route_checked = false;
    bool use_two_lll_route = true;
    RouteCheckCounter counter;
    run_grh_route_check(route_checked, use_two_lll_route, counter,
                        monic_field(real_31.data(),
                                    static_cast<slong>(real_31.size())));
    assert(route_checked && !use_two_lll_route);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.checks == 1);
    assert(counter.two_lll_order == 0);
#endif

    constexpr std::array<slong, 9> p5_m7_m11 = {
            144, -300, 259, -196, 162, -48, 24, -2, 1};
    route_checked = false;
    use_two_lll_route = false;
    counter = RouteCheckCounter{};
    run_grh_route_check(route_checked, use_two_lll_route, counter,
                        monic_field(p5_m7_m11.data(),
                                    static_cast<slong>(p5_m7_m11.size())));
    assert(route_checked && use_two_lll_route);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.checks == 1);
    assert(counter.two_lll_order > 0);
#endif
    return 0;
}

// The route check compares P1 / 2^(n d1) with 2^n P2 / 2^(n d2), evaluated as
// P1 2^(n d2) <= P2 2^(n (d1 + 1)).  With n = 2, d1 = 3 and d2 = 1 the
// threshold is P1 <= 64 P2.  A dropped "+ 1" would put it at 16 P2, and
// swapped d1 and d2 at P2 / 4.
int test_rounded_embedding_route_check_boundary() {
    namespace rs = silex::detail::relation_search;
    const auto diagonal = [](ulong a, ulong b) {
        sflint::FmpzMat gram(2, 2);
        sflint::fmpz_mat_zero(sflint::FmpzMatRef(gram));
        sflint::fmpz_set_ui(sflint::fmpz_mat_entry(gram, 0, 0), a);
        sflint::fmpz_set_ui(sflint::fmpz_mat_entry(gram, 1, 1), b);
        return gram;
    };
    rs::IdealLatticeLllData two_lll;
    two_lll.scaled_gram = diagonal(1, 1);
    sflint::fmpz_set_ui(sflint::FmpzRef(two_lll.gram_denominator), UWORD(2));
    const auto matches = [&](ulong a, ulong b) {
        sflint::FmpzMat one = diagonal(a, b);
        return rs::ideal_lattice_rounded_embedding_basis_matches_two_lll(
                sflint::FmpzMatConstRef(one), 3, two_lll);
    };
    assert(matches(8, 8));    // ratio exactly 2^n
    assert(matches(64, 1));
    assert(!matches(65, 1));  // just above 2^n
    assert(!matches(9, 8));

    // Malformed input: a denominator that is not a power of two.
    sflint::fmpz_set_ui(sflint::FmpzRef(two_lll.gram_denominator), UWORD(3));
    assert(!matches(1, 1));
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

// Checks the scaled Gram matrix of a totally real power-basis lattice with
// basis coordinates `basis` against the exact trace form: the Minkowski Gram
// matrix of the power basis is Tr(theta^(i+j)), so `scaled_gram` must agree
// with 2^denominator_bits * basis * Tr * basis^T up to rounding, which is
// far below 2^-20 of the diagonal.  Power sums come from Newton's identities
// for the monic polynomial with the given ascending coefficients.
template <std::size_t N>
bool scaled_gram_matches_trace_form(
        const std::array<slong, N>& coefficients,
        const sflint::FmpzMat& basis,
        const sflint::FmpzMat& scaled_gram,
        slong denominator_bits) noexcept {
    const slong degree = static_cast<slong>(N) - 1;
    if (degree <= 0 || coefficients[N - 1] != 1 ||
        sflint::fmpz_mat_nrows(basis) != degree ||
        sflint::fmpz_mat_ncols(basis) != degree ||
        sflint::fmpz_mat_nrows(scaled_gram) != degree ||
        sflint::fmpz_mat_ncols(scaled_gram) != degree) {
        return false;
    }

    sflint::FmpzMat power_sums(1, 2 * degree - 1);
    sflint::Fmpz coefficient;
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(power_sums, 0, 0), degree);
    for (slong k = 1; k < 2 * degree - 1; ++k) {
        sflint::FmpzRef sum = sflint::fmpz_mat_entry(power_sums, 0, k);
        sflint::fmpz_zero(sum);
        if (k <= degree) {
            sflint::fmpz_set_si(sflint::FmpzRef(coefficient),
                                k * coefficients[degree - k]);
            sflint::fmpz_sub(sum, sflint::FmpzConstRef(sum.raw()),
                             sflint::FmpzConstRef(coefficient));
        }
        for (slong i = 1; i <= std::min(k - 1, degree); ++i) {
            sflint::fmpz_set_si(sflint::FmpzRef(coefficient),
                                coefficients[degree - i]);
            sflint::fmpz_submul(
                    sum, sflint::FmpzConstRef(coefficient),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(power_sums), 0, k - i));
        }
    }

    sflint::FmpzMat trace_form(degree, degree);
    for (slong i = 0; i < degree; ++i) {
        for (slong j = 0; j < degree; ++j) {
            sflint::fmpz_set(sflint::fmpz_mat_entry(trace_form, i, j),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(power_sums), 0,
                                     i + j));
        }
    }

    sflint::FmpzMat basis_transpose(degree, degree);
    sflint::FmpzMat product(degree, degree);
    sflint::FmpzMat expected(degree, degree);
    sflint::fmpz_mat_transpose(sflint::FmpzMatRef(basis_transpose),
                               sflint::FmpzMatConstRef(basis));
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(trace_form),
                         sflint::FmpzMatConstRef(basis_transpose));
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(expected),
                         sflint::FmpzMatConstRef(basis),
                         sflint::FmpzMatConstRef(product));

    sflint::Fmpz difference;
    sflint::Fmpz bound;
    for (slong i = 0; i < degree; ++i) {
        for (slong j = 0; j < degree; ++j) {
            // difference = scaled_gram - 2^bits * expected, times 2^20.
            sflint::fmpz_mul_2exp(
                    sflint::FmpzRef(difference),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(expected), i, j),
                    static_cast<ulong>(denominator_bits));
            sflint::fmpz_sub(sflint::FmpzRef(difference),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(scaled_gram), i,
                                     j),
                             sflint::FmpzConstRef(difference));
            sflint::fmpz_mul_2exp(sflint::FmpzRef(difference),
                                  sflint::FmpzConstRef(difference), 20);
            // bound = 2^bits * (expected_ii + expected_jj).
            sflint::fmpz_add(sflint::FmpzRef(bound),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(expected), i, i),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(expected), j,
                                     j));
            sflint::fmpz_mul_2exp(sflint::FmpzRef(bound),
                                  sflint::FmpzConstRef(bound),
                                  static_cast<ulong>(denominator_bits));
            if (sflint::fmpz_sgn(sflint::FmpzConstRef(bound)) <= 0 ||
                sflint::fmpz_cmpabs(sflint::FmpzConstRef(difference),
                                    sflint::FmpzConstRef(bound)) > 0) {
                return false;
            }
        }
    }
    return true;
}

// Real subfield of Q(zeta_37), degree 18, in its power basis.  The
// embeddings of theta^17 are near 2^17, so a fixed 64-bit relative
// embedding guard misses the absolute rounding radius 0.1 * 2^-precision at
// every precision; precision doubling never closes that gap.  The LLL data
// must build at the relation route's first precision, as the absolute
// tolerance of the reference minkowski_gram_mat_scaled allows.
int test_large_power_basis_embeddings_round_at_first_precision() {
    constexpr std::array<slong, 19> coefficients = {
            -1,  9,     45,   -120, -330, 462, 924, -792, -1287, 715,
            1001, -364, -455, 105,  120,  -16, -17, 1,    1};
    constexpr slong kFirstPrecision = 200;

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    const silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order order = silex::test::equation_order(field);
    const slong degree = order.degree();
    assert(degree == 18);

    silex::Ideal unit_ideal(order);
    assert(unit_ideal.is_defined());
    assert(unit_ideal.one());

    sflint::Fmpz expected_denominator;
    sflint::fmpz_one(sflint::FmpzRef(expected_denominator));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(expected_denominator),
                          sflint::FmpzConstRef(expected_denominator),
                          static_cast<ulong>(kFirstPrecision / 2));

    silex::detail::relation_search::IdealLatticeLllData data;
    assert(silex::detail::relation_search::build_ideal_lattice_lll_data(
            data, unit_ideal, kFirstPrecision));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(data.gram_denominator),
                              sflint::FmpzConstRef(expected_denominator)));
    // The unit ideal's HNF is the identity, so data.basis is its own
    // transform from that HNF.
    sflint::FmpzMat hnf(degree, degree);
    sflint::FmpzMat identity(degree, degree);
    assert(unit_ideal.get_hnf(sflint::FmpzMatRef(hnf)));
    sflint::fmpz_mat_one(sflint::FmpzMatRef(identity));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(hnf),
                                  sflint::FmpzMatConstRef(identity)));
    assert(reduction_matrices_have_expected_properties(
            unit_ideal, data.basis, data.basis, data.scaled_gram));
    assert(scaled_gram_matches_trace_form(coefficients, data.basis,
                                          data.scaled_gram,
                                          kFirstPrecision / 2));

    sflint::FmpzMat weights(1, degree);
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(weights));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(weights, 0, 0), 16);
    sflint::FmpzMat basis(degree, degree);
    sflint::FmpzMat transform(degree, degree);
    sflint::FmpzMat scaled_gram(degree, degree);
    sflint::Fmpz gram_denominator;
    assert(silex::detail::relation_search::
                   build_weighted_ideal_lattice_reduction(
                           basis, transform, scaled_gram, gram_denominator,
                           unit_ideal, sflint::FmpzMatConstRef(weights),
                           kFirstPrecision, nullptr, nullptr));
    assert(reduction_matrices_have_expected_properties(
            unit_ideal, basis, transform, scaled_gram));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(gram_denominator),
                              sflint::FmpzConstRef(expected_denominator)));
    return 0;
}

// Checks the T2 enumeration data of an ideal of a maximal order: the basis
// spans the ideal lattice (same HNF), the product of the Gram-Schmidt pivots
// on the quadratic-form diagonal is det(Gram) = |disc| N(I)^2, and the basis
// meets the LLL bounds the reduction accepts:
//   |b_1|^2 <= 2^floor((n + 1) / 2) (d^(1/n) + 1)  and
//   prod_i |b_i|^2 <= 2^(n (n - 1) / 2) d.
bool t2_enumeration_data_is_reduced_ideal_basis(
        const silex::Ideal& ideal,
        const silex::Order& order,
        const silex::detail::relation_search::FiniteIdealT2EnumerationData&
                data) noexcept {
    const slong degree = order.degree();
    if (degree <= 0 || sflint::fmpz_mat_nrows(data.basis) != degree ||
        sflint::fmpz_mat_ncols(data.basis) != degree ||
        data.quadratic_form_data.size() !=
                static_cast<std::size_t>(degree * degree) ||
        !std::isfinite(data.initial_bound_value) ||
        data.initial_bound_value <= 0.0) {
        return false;
    }

    sflint::FmpzMat ideal_hnf(degree, degree);
    sflint::FmpzMat basis_hnf(degree, degree);
    if (!ideal.get_hnf(sflint::FmpzMatRef(ideal_hnf))) {
        return false;
    }
    ::fmpz_mat_hnf(basis_hnf.raw(), data.basis.raw());
    ::fmpz_mat_hnf(ideal_hnf.raw(), ideal_hnf.raw());
    if (!sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(basis_hnf),
                                sflint::FmpzMatConstRef(ideal_hnf))) {
        return false;
    }

    sflint::Fmpz discriminant;
    sflint::Fmpz norm;
    sflint::Fmpz determinant;
    if (!order.discriminant(sflint::FmpzRef(discriminant)) ||
        !ideal.norm(sflint::FmpzRef(norm))) {
        return false;
    }
    sflint::fmpz_abs(sflint::FmpzRef(determinant),
                     sflint::FmpzConstRef(discriminant));
    sflint::fmpz_mul(sflint::FmpzRef(determinant),
                     sflint::FmpzConstRef(determinant),
                     sflint::FmpzConstRef(norm));
    sflint::fmpz_mul(sflint::FmpzRef(determinant),
                     sflint::FmpzConstRef(determinant),
                     sflint::FmpzConstRef(norm));
    const double log2_determinant =
            std::log2(sflint::fmpz_get_d(sflint::FmpzConstRef(determinant)));
    if (!std::isfinite(log2_determinant)) {
        return false;
    }

    // Row i of the quadratic form holds the Gram-Schmidt pivot B_i on the
    // diagonal and mu_{j,i} at column j > i, so
    // |b_j|^2 = B_j + sum_{i<j} mu_{j,i}^2 B_i.
    const auto entry = [&](slong row, slong column) {
        return data.quadratic_form_data[static_cast<std::size_t>(
                row * degree + column)];
    };
    double log2_pivot_product = 0.0;
    double log2_norm_product = 0.0;
    for (slong j = 0; j < degree; ++j) {
        const double pivot = entry(j, j);
        if (!std::isfinite(pivot) || pivot <= 0.0) {
            return false;
        }
        log2_pivot_product += std::log2(pivot);
        double squared_norm = pivot;
        for (slong i = 0; i < j; ++i) {
            squared_norm += entry(i, j) * entry(i, j) * entry(i, i);
        }
        log2_norm_product += std::log2(squared_norm);
    }
    constexpr double kLog2Tolerance = 1.0e-6;
    if (std::fabs(log2_pivot_product - log2_determinant) > kLog2Tolerance) {
        return false;
    }

    const double n = static_cast<double>(degree);
    const double log2_first_bound =
            static_cast<double>((degree + 1) / 2) +
            std::log2(std::exp2(log2_determinant / n) + 1.0);
    const double log2_product_bound =
            (n * (n - 1.0)) / 2.0 + log2_determinant;
    return std::log2(entry(0, 0)) <= log2_first_bound + kLog2Tolerance &&
           log2_norm_product <= log2_product_bound + kLog2Tolerance;
}

// A degree-10 field given by a polynomial with coefficients near 2^18, so
// the maximal-order basis
// embeddings reach about 2^34 while the ideal lattices have short vectors
// near 2^7.  Rounding those embeddings at the first full-rank scale (2^4)
// leaves an LLL transform whose rounding error dominates the reduced
// vectors: for a prime above 89 the double Gram-Schmidt of the resulting
// Gram matrix has a nonpositive last pivot at every working precision.  The
// T2 enumeration data must build for every prime above 89.
int test_skewed_basis_t2_enumeration_data_builds() {
    constexpr std::array<slong, 11> coefficients = {
            238831, -228566, -128567, 146666, -15113, -13352,
            2746,   334,     -97,     -2,     1};

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    const silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));
    assert(order.is_maximal());
    const slong degree = order.degree();
    assert(degree == 10);

    sflint::Fmpz p;
    sflint::fmpz_set_ui(sflint::FmpzRef(p), UWORD(89));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == degree);
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        silex::Ideal prime_ideal(order);
        assert(prime_ideal.is_defined());
        assert(prime->get_ideal(prime_ideal));

        silex::detail::relation_search::FiniteIdealT2EnumerationData data;
        assert(silex::detail::relation_search::
                       build_finite_ideal_t2_enumeration_data_with_retry(
                               data, prime_ideal, nullptr, nullptr));
        assert(t2_enumeration_data_is_reduced_ideal_basis(prime_ideal, order,
                                                          data));
    }
    return 0;
}

// The totally real octic x^8 - 20x^6 + 95x^4 - 100x^2 + 25 (LMFDB
// 8.8.82944000000.1).  For prime ideals above 5, 11 and 13 the basis reduced
// from the first full-rank rounding has quadratic-form data but misses the
// LLL bounds; the T2 enumeration data must then come from a basis that meets
// them.
int test_octic_t2_enumeration_data_meets_lll_bounds() {
    constexpr std::array<slong, 9> coefficients = {25, 0, -100, 0, 95,
                                                   0,  -20, 0,  1};

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    const silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));
    assert(order.is_maximal());
    assert(order.degree() == 8);

    constexpr std::array<ulong, 6> rational_primes = {2, 3, 5, 7, 11, 13};
    for (const ulong rational_prime : rational_primes) {
        sflint::Fmpz p;
        sflint::fmpz_set_ui(sflint::FmpzRef(p), rational_prime);
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order,
                                      sflint::FmpzConstRef(p)));
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            silex::Ideal prime_ideal(order);
            assert(prime->get_ideal(prime_ideal));

            silex::detail::relation_search::FiniteIdealT2EnumerationData data;
            assert(silex::detail::relation_search::
                           build_finite_ideal_t2_enumeration_data_with_retry(
                                   data, prime_ideal, nullptr, nullptr));
            assert(t2_enumeration_data_is_reduced_ideal_basis(prime_ideal,
                                                              order, data));
        }
    }
    return 0;
}

}  // namespace

// The integer LLL of the rounded-embedding route certifies FLINT's first
// phase with fmpz_mat_is_reduced instead of letting fmpz_lll run its own
// check; the result must be fmpz_lll's on every route.  On the octic field
// of row multiquadratic_m3_m5_m7 the one-LLL basis of a degree-one prime
// above 43 is certified on the first phase and then fails the T2 diagonal
// checks, so the ideal is reduced by the two-LLL route, as when the whole
// fmpz_lll ran.
int test_rounded_embedding_certified_lll_matches_fmpz_lll() {
    namespace rs = silex::detail::relation_search;
    constexpr std::array<slong, 9> m3_m5_m7 = {1201, 598, -524, -146, 69,
                                              2,    4,   -4,   1};
    FieldSetup setup = monic_field(m3_m5_m7.data(),
                                   static_cast<slong>(m3_m5_m7.size()));
    const silex::Order& order = setup.maximal_order;
    const slong degree = order.degree();
    assert(degree == 8);
    silex::Ideal prime(order);
    assert(degree_one_prime_ideal_above(prime, order, UWORD(43), 0));

    rs::IdealLatticeRoundedEmbeddingCache cache;
    rs::IdealLatticeLllData data;
    rs::IdealLatticeLllData expected;
    assert(rs::build_ideal_lattice_rounded_embedding_lll_data(
            data, prime, kReductionPrecision, cache));
    // The order data was formed and the one-LLL basis was rejected: an
    // accepted one would have run the route check.
    assert(cache.has_rounded_rows && !cache.route_checked &&
           !cache.use_two_lll_route);
    assert(rs::build_ideal_lattice_lll_data(expected, prime,
                                            kReductionPrecision));
    assert(rounded_embedding_data_has_contract(data, prime));
    assert(same_lll_data(data, expected));

    // On the rejected ideal's embedded basis H R, every route returns
    // fmpz_lll's matrix and transform.
    sflint::FmpzMat hnf(degree, degree);
    assert(prime.get_hnf(sflint::FmpzMatRef(hnf)));
    sflint::FmpzMat input(degree, degree);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(input),
                         sflint::FmpzMatConstRef(hnf),
                         sflint::FmpzMatConstRef(cache.rounded_rows));
    const auto check_routes = [](const sflint::FmpzMat& matrix,
                                 rs::IdealLatticeOneLllRoute route) {
        const slong rows = sflint::fmpz_mat_nrows(matrix);
        const slong cols = sflint::fmpz_mat_ncols(matrix);
        sflint::FmpzMat reference(rows, cols);
        sflint::FmpzMat reference_transform(rows, rows);
        sflint::fmpz_mat_set(sflint::FmpzMatRef(reference),
                             sflint::FmpzMatConstRef(matrix));
        sflint::fmpz_mat_one(sflint::FmpzMatRef(reference_transform));
        sflint::FmpzLll config(0.99, 0.51, Z_BASIS, APPROX);
        ::fmpz_lll(reference.raw(), reference_transform.raw(), config.raw());
        for (const auto failure :
             {rs::IdealLatticeOneLllTestFailure::none,
              rs::IdealLatticeOneLllTestFailure::reducer,
              rs::IdealLatticeOneLllTestFailure::certificate}) {
            sflint::FmpzMat reduced(rows, cols);
            sflint::FmpzMat transform(0, 0);
            sflint::fmpz_mat_set(sflint::FmpzMatRef(reduced),
                                 sflint::FmpzMatConstRef(matrix));
            const rs::IdealLatticeOneLllRoute expected_route =
                    route != rs::IdealLatticeOneLllRoute::certified
                            ? route
                    : failure == rs::IdealLatticeOneLllTestFailure::none
                            ? rs::IdealLatticeOneLllRoute::certified
                            : rs::IdealLatticeOneLllRoute::fallback;
            if (rs::ideal_lattice_rounded_embedding_integer_lll(
                        reduced, transform, failure) != expected_route ||
                !sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(reduced),
                                        sflint::FmpzMatConstRef(reference)) ||
                !sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(transform),
                        sflint::FmpzMatConstRef(reference_transform))) {
                return false;
            }
        }
        return true;
    };
    assert(check_routes(input, rs::IdealLatticeOneLllRoute::certified));

    // The one-LLL basis the route rejected is the fmpz_lll one, and it is
    // not the two-LLL basis the route returned.
    sflint::FmpzMat one_lll_basis(0, 0);
    assert(one_lll_basis_from_cache(one_lll_basis, prime, cache));
    assert(!sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(one_lll_basis),
                                   sflint::FmpzMatConstRef(data.basis)));

    // Entries above FLINT's 250-bit ULLL threshold, and a single row, go
    // to fmpz_lll directly.
    sflint::FmpzMat wide(degree, degree);
    ::fmpz_mat_scalar_mul_2exp(wide.raw(), input.raw(), 260);
    assert(check_routes(wide, rs::IdealLatticeOneLllRoute::ineligible));
    sflint::FmpzMat single(1, 3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(single, 0, 0), 6);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(single, 0, 2), -4);
    assert(check_routes(single, rs::IdealLatticeOneLllRoute::ineligible));
    return 0;
}

int main() {
    assert(test_ideal_lattice_result_contracts() == 0);
    assert(test_rounded_embedding_lll_contracts() == 0);
    assert(test_rounded_embedding_route_check_for_badly_rounded_order() == 0);
    assert(test_rounded_embedding_route_check_boundary() == 0);
    assert(test_rounded_embedding_route_check_once_per_context() == 0);
    assert(test_rounded_embedding_certified_lll_matches_fmpz_lll() == 0);
    assert(test_reduction_identities_and_mutation() == 0);
    assert(test_short_elements_and_reduction_matrices() == 0);
    assert(test_deterministic_failures_preserve_outputs() == 0);
    assert(test_large_power_basis_embeddings_round_at_first_precision() ==
           0);
    assert(test_skewed_basis_t2_enumeration_data_builds() == 0);
    assert(test_octic_t2_enumeration_data_meets_lll_bounds() == 0);
    return 0;
}
