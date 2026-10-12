#include <silex/diagnostics.hpp>
#include <silex/factored_element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include "class_group/class_group_internal.hpp"
#include "factored_element/compact_reconstruction_internal.hpp"
#include "order_unit/compact_reconstruction_bound_internal.hpp"
#include "test_support.hpp"

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>

#include <cassert>
#include <cstddef>
#include <cstring>
#include <span>
#include <vector>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField index_three_quartic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -5);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool basis_begins_with_one(const silex::Order& order) noexcept {
    sflint::FmpqMat basis(order.degree(), order.degree());
    if (!order.get_basis(sflint::FmpqMatRef(basis))) {
        return false;
    }
    for (slong column = 0; column < order.degree(); ++column) {
        const slong expected = column == 0 ? 1 : 0;
        if (!sflint::fmpq_equal_si(
                    sflint::fmpq_mat_entry(sflint::FmpqMatConstRef(basis), 0,
                                           column),
                    expected)) {
            return false;
        }
    }
    return true;
}

void exact_outputs_and_bound(
        std::vector<silex::Element>& outputs, sflint::Fmpz& bound,
        const silex::Order& order,
        std::span<const silex::FactoredElement> inputs) noexcept {
    assert(order.parent() != nullptr);
    sflint::fmpz_zero(sflint::FmpzRef(bound));
    outputs.clear();
    outputs.reserve(inputs.size());

    sflint::FmpqMat rational_coordinates(1, order.degree());
    sflint::FmpzMat integral_coordinates(1, order.degree());
    sflint::Fmpz absolute_value;
    for (const silex::FactoredElement& input : inputs) {
        outputs.emplace_back(*order.parent());
        assert(input.evaluate(outputs.back()));
        assert(order.coordinates(sflint::FmpqMatRef(rational_coordinates),
                                 outputs.back()));
        assert(fmpq_mat_get_fmpz_mat(integral_coordinates.raw(),
                                     rational_coordinates.raw()) != 0);
        for (slong column = 0; column < order.degree(); ++column) {
            fmpz_abs(absolute_value.raw(),
                     fmpz_mat_entry(integral_coordinates.raw(), 0, column));
            if (fmpz_cmp(absolute_value.raw(), bound.raw()) > 0) {
                fmpz_set(bound.raw(), absolute_value.raw());
            }
        }
    }
}

void assert_equal_outputs(const std::vector<silex::Element>& actual,
                          const std::vector<silex::Element>& expected) {
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) {
        assert(actual[i].equal(expected[i]));
    }
}

void assert_exact_crt_crossing(
        const silex::detail::BoundedCompactReconstructionReport& report,
        sflint::FmpzConstRef bound) noexcept {
    sflint::Fmpz twice_bound;
    sflint::Fmpz previous_modulus;
    fmpz_mul_2exp(twice_bound.raw(), bound.raw(), 1);
    assert(fmpz_cmp(report.centered_crt_modulus.raw(), twice_bound.raw()) > 0);
    assert(report.last_prime != 0);
    fmpz_divexact_ui(previous_modulus.raw(), report.centered_crt_modulus.raw(),
                     report.last_prime);
    assert(fmpz_cmp(previous_modulus.raw(), twice_bound.raw()) <= 0);
}

silex::FactoredElementSpan factored_span(
        const std::vector<silex::FactoredElement>& values) noexcept {
    return {values.data(), values.size()};
}

int test_shared_exponent_matrix_and_centered_crt() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::Element theta(field);
    silex::Element one(field);
    silex::Element numerator(field);
    silex::Element unit(field);
    silex::Element minus_one(field);
    assert(theta.gen());
    assert(one.one());
    assert(numerator.add(theta, one));
    assert(unit.scalar_div_si(numerator, 2));
    assert(minus_one.set_si(-1));
    assert(maximal.contains(unit));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    inputs.emplace_back(field);
    inputs.emplace_back(field);
    assert(inputs[0].push(unit, 257));
    assert(inputs[1].push(unit, -233));
    assert(inputs[1].push(minus_one, 1));
    assert(inputs[2].push(unit, 89));
    assert(inputs[2].push(minus_one, -7));

    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, maximal, inputs);
    assert(fmpz_bits(bound.raw()) > FLINT_BITS);

    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, maximal, inputs, sflint::FmpzConstRef(bound)));
    assert(report.status ==
           silex::detail::BoundedCompactReconstructionStatus::success);
    assert(report.factor_rows == 2);
    assert(report.primes_used >= 2);
    assert(report.prime_trials == report.primes_used);
    assert(report.denominator_prime_rejections == 0);
    assert(report.noninvertible_prime_rejections == 0);
    assert_exact_crt_crossing(report, sflint::FmpzConstRef(bound));
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_denominator_prime_exclusion() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::Element half(field);
    silex::Element two(field);
    assert(half.set_si_over_si(1, 2));
    assert(two.set_si(2));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(half, 1));
    assert(inputs[0].push(two, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, order, inputs);
    assert(sflint::fmpz_equal_si(bound, 1));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 1;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, order, inputs, sflint::FmpzConstRef(bound),
            options));
    assert(report.prime_trials == 2);
    assert(report.primes_used == 1);
    assert(report.denominator_prime_rejections == 1);
    assert(report.noninvertible_prime_rejections == 0);
    assert(report.last_prime == 3);
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_field_polynomial_denominator_prime_exclusion() {
    sflint::Fmpq half;
    sflint::FmpqPoly polynomial;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, half);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order order = silex::test::equation_order(field);
    silex::Element one(field);
    assert(one.one());

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(one, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, order, inputs);
    assert(sflint::fmpz_equal_si(bound, 1));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 1;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, order, inputs, sflint::FmpzConstRef(bound),
            options));
    assert(report.prime_trials == 2);
    assert(report.primes_used == 1);
    assert(report.field_polynomial_prime_rejections == 1);
    assert(report.order_basis_prime_rejections == 0);
    assert(report.denominator_prime_rejections == 0);
    assert(report.noninvertible_prime_rejections == 0);
    assert(report.last_prime == 3);
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_noninvertible_prime_exclusion() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::Element three(field);
    silex::Element six(field);
    assert(three.set_si(3));
    assert(six.set_si(6));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(three, -1));
    assert(inputs[0].push(six, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, order, inputs);
    assert(sflint::fmpz_equal_si(bound, 2));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 2;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, order, inputs, sflint::FmpzConstRef(bound),
            options));
    assert(report.prime_trials == 2);
    assert(report.primes_used == 1);
    assert(report.denominator_prime_rejections == 0);
    assert(report.noninvertible_prime_rejections == 1);
    assert(report.last_prime == 5);
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_order_basis_prime_exclusion() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::Element theta(field);
    silex::Element one(field);
    silex::Element numerator(field);
    silex::Element unit(field);
    assert(theta.gen());
    assert(one.one());
    assert(numerator.add(theta, one));
    assert(unit.scalar_div_si(numerator, 2));
    assert(maximal.contains(unit));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(unit, 1));
    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, maximal, inputs);
    assert(sflint::fmpz_equal_si(bound, 1));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 1;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, maximal, inputs, sflint::FmpzConstRef(bound),
            options));
    assert(report.prime_trials == 2);
    assert(report.primes_used == 1);
    assert(report.field_polynomial_prime_rejections == 0);
    assert(report.order_basis_prime_rejections == 1);
    assert(report.denominator_prime_rejections == 0);
    assert(report.noninvertible_prime_rejections == 0);
    assert(report.last_prime == 3);
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_bound_violation_fails_without_publication() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::Element two(field);
    assert(two.set_si(2));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(two, 1));

    sflint::Fmpz false_bound;
    sflint::fmpz_one(sflint::FmpzRef(false_bound));
    std::vector<silex::Element> output;
    output.emplace_back(field);
    assert(output[0].set_si(99));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 1;
    silex::detail::BoundedCompactReconstructionReport report;
    assert(!silex::detail::bounded_compact_reconstruct(
            report, output, order, inputs, sflint::FmpzConstRef(false_bound),
            options));
    assert(report.status == silex::detail::BoundedCompactReconstructionStatus::
                                    certified_bound_violated);
    assert(output.size() == 1);
    assert(output[0].equal_si(99));
    return 0;
}

int test_single_prime_centered_negative_lift() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::Element minus_two(field);
    assert(minus_two.set_si(-2));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(minus_two, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz bound;
    exact_outputs_and_bound(expected, bound, order, inputs);
    assert(sflint::fmpz_equal_si(bound, 2));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 3;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            report, actual, order, inputs, sflint::FmpzConstRef(bound),
            options));
    assert(report.primes_used == 1);
    assert(report.last_prime == 5);
    assert_exact_crt_crossing(report, sflint::FmpzConstRef(bound));
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_empty_batch() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    sflint::Fmpz zero_bound;
    sflint::fmpz_zero(sflint::FmpzRef(zero_bound));

    silex::detail::BoundedCompactReconstructionOptions options;
    options.prime_search_start = 1;
    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> output;
    std::span<const silex::FactoredElement> empty;
    assert(silex::detail::bounded_compact_reconstruct(
            report, output, order, empty, sflint::FmpzConstRef(zero_bound),
            options));
    assert(output.empty());
    assert(report.primes_used == 0);
    return 0;
}

int test_coordinate_bound_and_reconstruction() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(basis_begins_with_one(maximal));

    silex::Element theta(field);
    silex::Element one(field);
    silex::Element numerator(field);
    silex::Element unit(field);
    assert(theta.gen());
    assert(one.one());
    assert(numerator.add(theta, one));
    assert(unit.scalar_div_si(numerator, 2));
    assert(maximal.contains(unit));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(unit, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz exact_maximum;
    exact_outputs_and_bound(expected, exact_maximum, maximal, inputs);

    silex::EmbeddingContext embeddings(field);
    sflint::Fmpz source_bound;
    silex::detail::CompactCoordinateBoundReport bound_report;
    assert(silex::detail::compact_unit_coordinate_bound(
            bound_report, sflint::FmpzRef(source_bound), maximal,
            factored_span(inputs), embeddings));
    assert(bound_report.status ==
           silex::detail::CompactCoordinateBoundStatus::success);
    assert(bound_report.attempts == 1);
    assert(bound_report.precision == 128);
    assert(sflint::fmpz_equal(
            sflint::FmpzConstRef(source_bound),
            sflint::FmpzConstRef(bound_report.coordinate_bound)));
    assert(sflint::fmpz_cmp(sflint::FmpzConstRef(source_bound),
                            sflint::FmpzConstRef(exact_maximum)) >= 0);
    assert(sflint::arb_is_finite(bound_report.log2_minkowski_bound));
    assert(sflint::arb_is_finite(bound_report.log2_unit_bound));
    assert(sflint::arb_is_finite(bound_report.log2_coordinate_bound));

    silex::detail::BoundedCompactReconstructionReport reconstruction_report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            reconstruction_report, actual, maximal, inputs,
            sflint::FmpzConstRef(source_bound)));
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_coordinate_bound_rank_zero() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    assert(order.maximality_known());
    assert(order.is_maximal());

    silex::EmbeddingContext embeddings(field);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 99);
    silex::detail::CompactCoordinateBoundReport report;
    assert(silex::detail::compact_unit_coordinate_bound(
            report, sflint::FmpzRef(bound), order, {}, embeddings));
    assert(report.status ==
           silex::detail::CompactCoordinateBoundStatus::success);
    assert(report.bit_bound == 0);
    assert(sflint::fmpz_equal_si(bound, 1));
    assert(sflint::arb_is_zero(report.log2_minkowski_bound));
    assert(sflint::arb_is_zero(report.log2_unit_bound));
    assert(sflint::arb_is_zero(report.log2_coordinate_bound));
    return 0;
}

int test_coordinate_bound_rejects_nonmaximal_order() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    assert(equation.maximality_known());
    assert(!equation.is_maximal());

    silex::Element one(field);
    assert(one.one());
    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(one, 1));

    silex::EmbeddingContext embeddings(field);
    sflint::Fmpz output;
    sflint::fmpz_set_si(sflint::FmpzRef(output), 99);
    silex::detail::CompactCoordinateBoundReport report;
    assert(!silex::detail::compact_unit_coordinate_bound(
            report, sflint::FmpzRef(output), equation, factored_span(inputs),
            embeddings));
    assert(report.status ==
           silex::detail::CompactCoordinateBoundStatus::invalid_input);
    assert(sflint::fmpz_equal_si(output, 99));
    return 0;
}

int test_coordinate_bound_supports_reordered_basis() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::FmpqMat basis(2, 2);
    assert(maximal.get_basis(sflint::FmpqMatRef(basis)));
    for (slong column = 0; column < 2; ++column) {
        fmpq_swap(fmpq_mat_entry(basis.raw(), 0, column),
                  fmpq_mat_entry(basis.raw(), 1, column));
    }
    silex::Order reordered_basis =
            silex::Order::from_basis(field, sflint::FmpqMatConstRef(basis));
    assert(reordered_basis.is_defined());
    silex::Order reordered =
            silex::test::verified_maximal_order(reordered_basis);
    assert(!basis_begins_with_one(reordered));

    silex::Element theta(field);
    silex::Element one(field);
    silex::Element numerator(field);
    silex::Element unit(field);
    assert(theta.gen());
    assert(one.one());
    assert(numerator.add(theta, one));
    assert(unit.scalar_div_si(numerator, 2));
    assert(reordered.contains(unit));

    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    assert(inputs[0].push(unit, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz exact_maximum;
    exact_outputs_and_bound(expected, exact_maximum, reordered, inputs);

    silex::EmbeddingContext embeddings(field);
    sflint::Fmpz output;
    silex::detail::CompactCoordinateBoundReport report;
    assert(silex::detail::compact_unit_coordinate_bound(
            report, sflint::FmpzRef(output), reordered, factored_span(inputs),
            embeddings));
    assert(report.status ==
           silex::detail::CompactCoordinateBoundStatus::success);
    assert(sflint::fmpz_cmp(sflint::FmpzConstRef(output),
                            sflint::FmpzConstRef(exact_maximum)) >= 0);

    silex::detail::BoundedCompactReconstructionReport reconstruction_report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            reconstruction_report, actual, reordered, inputs,
            sflint::FmpzConstRef(output)));
    assert_equal_outputs(actual, expected);
    return 0;
}

int test_coordinate_bound_supports_index_three_maximal_basis() {
    silex::NumberField field = index_three_quartic_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Fmpz equation_index;
    assert(silex::order_index(sflint::FmpzRef(equation_index), equation,
                              maximal));
    assert(sflint::fmpz_equal_si(equation_index, 3));
    assert(!basis_begins_with_one(maximal));

    silex::Element one(field);
    assert(one.one());
    std::vector<silex::FactoredElement> inputs;
    inputs.emplace_back(field);
    inputs.emplace_back(field);
    assert(inputs[0].push(one, 1));
    assert(inputs[1].push(one, 1));

    std::vector<silex::Element> expected;
    sflint::Fmpz exact_maximum;
    exact_outputs_and_bound(expected, exact_maximum, maximal, inputs);

    silex::EmbeddingContext embeddings(field);
    sflint::Fmpz output;
    silex::detail::CompactCoordinateBoundReport report;
    assert(silex::detail::compact_unit_coordinate_bound(
            report, sflint::FmpzRef(output), maximal, factored_span(inputs),
            embeddings));
    assert(report.status ==
           silex::detail::CompactCoordinateBoundStatus::success);
    assert(sflint::fmpz_cmp(sflint::FmpzConstRef(output),
                            sflint::FmpzConstRef(exact_maximum)) >= 0);

    silex::detail::BoundedCompactReconstructionReport reconstruction_report;
    std::vector<silex::Element> actual;
    assert(silex::detail::bounded_compact_reconstruct(
            reconstruction_report, actual, maximal, inputs,
            sflint::FmpzConstRef(output)));
    assert_equal_outputs(actual, expected);
    return 0;
}

silex::NumberField monic_field(std::span<const slong> coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// Keys of the compact field-modulus cache: an entry is returned only for its
// own field data, prime and threshold, and entries stay sorted when a miss
// is stored before existing entries.
int test_compact_field_modulus_cache_keys() {
    const slong k_coefficients[] = {-1, -1, 0, 0, 1};
    const slong l_coefficients[] = {-2, 0, 0, 0, 1};
    const silex::NumberField k = monic_field(k_coefficients);
    const silex::NumberField k_copy = k;
    const silex::NumberField l = monic_field(l_coefficients);
    assert(k.is_defined() && l.is_defined() && !k.has_same_data(l));

    sflint::Fmpz p0;
    sflint::Fmpz p1;
    sflint::Fmpz p2;
    fmpz_one(p0.raw());
    fmpz_mul_2exp(p0.raw(), p0.raw(), 60);
    fmpz_nextprime(p0.raw(), p0.raw(), 1);
    fmpz_nextprime(p1.raw(), p0.raw(), 1);
    fmpz_nextprime(p2.raw(), p1.raw(), 1);

    sflint::FmpzModCtx ctx(p1.raw());
    sflint::FmpzModPoly modulus(ctx);
    assert(modulus.is_initialized());
    fmpz_mod_poly_set_coeff_si(modulus.raw(), 4, 1, ctx.raw());
    fmpz_mod_poly_set_coeff_si(modulus.raw(), 1, -1, ctx.raw());
    fmpz_mod_poly_set_coeff_si(modulus.raw(), 0, -1, ctx.raw());

    silex::detail::CompactFieldModulusCache cache;
    std::size_t position = 99;
    assert(silex::detail::compact_field_modulus_cache_find(
                   position, cache, k, sflint::FmpzConstRef(p1), 3) ==
           nullptr);
    assert(position == 0 && cache.field.has_same_data(k));
    silex::detail::compact_field_modulus_cache_store(
            cache, position, sflint::FmpzConstRef(p1), 3, true, modulus,
            ctx);

    // Hit for the same key, through a copy sharing the field data.
    const auto* hit = silex::detail::compact_field_modulus_cache_find(
            position, cache, k_copy, sflint::FmpzConstRef(p1), 3);
    assert(hit != nullptr && hit->usable && position == 0);
    assert(fmpz_equal(hit->prime.raw(), p1.raw()));
    sflint::FmpzPoly expected;
    fmpz_mod_poly_get_fmpz_poly(expected.raw(), modulus.raw(), ctx.raw());
    assert(fmpz_poly_equal(hit->modulus.raw(), expected.raw()));

    // Misses after that hit: another prime, another threshold.
    assert(silex::detail::compact_field_modulus_cache_find(
                   position, cache, k, sflint::FmpzConstRef(p2), 3) ==
           nullptr);
    assert(position == 1);
    silex::detail::compact_field_modulus_cache_store(
            cache, position, sflint::FmpzConstRef(p2), 3, false, modulus,
            ctx);
    assert(silex::detail::compact_field_modulus_cache_find(
                   position, cache, k, sflint::FmpzConstRef(p1), 4) ==
           nullptr);
    assert(position == 2);

    // A miss before the first entry is stored in order.
    assert(silex::detail::compact_field_modulus_cache_find(
                   position, cache, k, sflint::FmpzConstRef(p0), 3) ==
           nullptr);
    assert(position == 0);
    silex::detail::compact_field_modulus_cache_store(
            cache, position, sflint::FmpzConstRef(p0), 3, false, modulus,
            ctx);
    assert(cache.entries.size() == 3);
    const sflint::Fmpz* order[] = {&p0, &p1, &p2};
    for (std::size_t i = 0; i < cache.entries.size(); ++i) {
        assert(fmpz_equal(cache.entries[i].prime.raw(), order[i]->raw()));
    }
    hit = silex::detail::compact_field_modulus_cache_find(
            position, cache, k, sflint::FmpzConstRef(p2), 3);
    assert(hit != nullptr && !hit->usable && position == 2);
    hit = silex::detail::compact_field_modulus_cache_find(
            position, cache, k, sflint::FmpzConstRef(p1), 3);
    assert(hit != nullptr && hit->usable && position == 1);

    // Another field never sees these entries: the cache is reset.
    assert(silex::detail::compact_field_modulus_cache_find(
                   position, cache, l, sflint::FmpzConstRef(p1), 3) ==
           nullptr);
    assert(position == 0 && cache.entries.empty());
    assert(cache.field.has_same_data(l));
    return 0;
}

bool compact_power_root_squares_to(const silex::FactoredElement& root,
                                   const silex::FactoredElement& input) {
    const silex::NumberField* field = input.parent();
    silex::Element root_value(*field);
    silex::Element square(*field);
    silex::Element value(*field);
    return root.evaluate(root_value) &&
           square.multiply(root_value, root_value) &&
           input.evaluate(value) && square.equal(value);
}

silex::Element element_from_si(const silex::NumberField& field,
                               std::span<const slong> coefficients,
                               slong denominator) {
    sflint::FmpqPoly polynomial;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    ::fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), denominator);
    silex::Element element(field);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    return element;
}

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
void count_debug_failure(void* user,
                         silex::DiagnosticsModule,
                         silex::DebugLevel,
                         const char*,
                         const char*,
                         int,
                         const char*,
                         const char*) noexcept {
    ++*static_cast<slong*>(user);
}
#endif

// detail::compact_power_input_norm, the N(input) that
// compact_power_root_from_presentation passes to Element::is_power, equals
// the norm of the expanded input remainder_value * scale^n: bases of either
// norm sign, a non-integral base, negative exponents, scales other than one
// and the empty product.  Then unit compact power roots in the cubic field,
// one of norm -1, are found exactly with the debug cross-check of that norm
// against the expanded input's norm (diagnostics builds) silent.
int test_compact_power_input_norm() {
    const slong coefficients[] = {-1, -1, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const slong a_coefficients[] = {2, -1};
    const slong b_coefficients[] = {1, 1};
    const slong c_coefficients[] = {5, 0, -1};
    const slong d_coefficients[] = {-2, 1};
    const silex::Element a = element_from_si(field, a_coefficients, 1);
    const silex::Element b = element_from_si(field, b_coefficients, 3);
    const silex::Element c = element_from_si(field, c_coefficients, 1);
    const silex::Element d = element_from_si(field, d_coefficients, 1);

    silex::FactoredElement product(field);
    assert(product.push(a, 2) && product.push(b, -1) && product.push(c, 4) &&
           product.push(d, 1));
    silex::FactoredElement empty(field);
    assert(empty.one());
    for (const silex::FactoredElement* remainder : {&product, &empty}) {
        silex::Element value(field);
        assert(remainder->evaluate(value));
        for (slong n : {3, 5}) {
            for (slong scale_si : {1, 2, -3}) {
                sflint::Fmpz scale;
                sflint::Fmpz scale_power;
                fmpz_set_si(scale.raw(), scale_si);
                fmpz_pow_ui(scale_power.raw(), scale.raw(),
                            static_cast<ulong>(n));
                silex::Element scale_element(field);
                silex::Element input(field);
                sflint::Fmpq expected;
                sflint::Fmpq norm;
                assert(scale_element.set_fmpz(
                               sflint::FmpzConstRef(scale_power)) &&
                       input.multiply(value, scale_element) &&
                       input.norm(sflint::FmpqRef(expected)));
                assert(silex::detail::compact_power_input_norm(
                        norm, *remainder, sflint::FmpzConstRef(scale), n,
                        field.degree()));
                assert(fmpq_equal(norm.raw(), expected.raw()) != 0);
            }
        }
    }

    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::Element theta(field);
    silex::Element minus_theta(field);
    silex::Element theta_squared(field);
    assert(theta.gen() && minus_theta.negate(theta) &&
           theta_squared.multiply(theta, theta));
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    slong failures = 0;
#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    silex::diagnostics_set_debug_checks(
            diagnostics, silex::DebugLevel::expensive,
            silex::diagnostics_all_modules, count_debug_failure, &failures);
#endif
    constexpr slong precision = 128;
    // theta * (theta^2)^2 = theta^5 and (-theta) * (theta^2)^2 = (-theta)^5,
    // of norms 1 and -1.
    for (const silex::Element* unit : {&theta, &minus_theta}) {
        silex::FactoredElement element(field);
        assert(element.push(*unit, 1) && element.push(theta_squared, 2));
        bool is_power = false;
        silex::FactoredElement root(field);
        assert(silex::detail::compact_unit_presentation_power_root(
                is_power, root, order, element, 5, precision, precision,
                &diagnostics, nullptr));
        assert(is_power);
        silex::Element root_value(field);
        silex::Element power(field);
        silex::Element value(field);
        sflint::Fmpz five;
        fmpz_set_si(five.raw(), 5);
        assert(root.evaluate(root_value) &&
               power.pow_fmpz(root_value, sflint::FmpzConstRef(five)) &&
               element.evaluate(value) && power.equal(value));
    }
    assert(failures == 0);
    return 0;
}

// One field-modulus cache shared by several compact power-root calls: the
// first call fills it, a repeated call adds no entry, and every call gives
// the result of an uncached call.
int test_compact_field_modulus_cache_shared_across_calls() {
    // x^4 - x - 1 has degree above the threshold 3, so the auxiliary prime
    // walk rejects primes where it has an irreducible quartic factor.
    const slong coefficients[] = {-1, -1, 0, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));

    silex::Element x(field);
    silex::Element x_cubed(field);
    silex::Element x_plus_one(field);
    silex::Element one(field);
    sflint::Fmpz three;
    fmpz_set_si(three.raw(), 3);
    assert(x.gen() && one.one() && x_plus_one.add(x, one) &&
           x_cubed.pow_fmpz(x, sflint::FmpzConstRef(three)));

    // x * x^3 * (x + 1)^2 = x^12 and x * x^3 * (x + 1)^4 = x^20 are squares
    // of units with a nonempty formal leftover.
    silex::FactoredElement square(field);
    silex::FactoredElement other_square(field);
    assert(square.push(x, 1) && square.push(x_cubed, 1) &&
           square.push(x_plus_one, 2));
    assert(other_square.push(x, 1) && other_square.push(x_cubed, 1) &&
           other_square.push(x_plus_one, 4));

    constexpr slong precision = 128;
    silex::detail::CompactFieldModulusCache cache;
    bool is_power = false;
    silex::FactoredElement root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, square, 2, precision, precision, nullptr,
            &cache));
    assert(is_power && compact_power_root_squares_to(root, square));
    const std::size_t filled = cache.entries.size();
    assert(filled > 0 && cache.field.has_same_data(field));
    bool saw_rejected = false;
    for (const auto& entry : cache.entries) {
        saw_rejected = saw_rejected || !entry.usable;
    }
    assert(saw_rejected);

    silex::FactoredElement second_root(field);
    is_power = false;
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, second_root, order, square, 2, precision, precision,
            nullptr, &cache));
    assert(is_power && compact_power_root_squares_to(second_root, square));
    assert(cache.entries.size() == filled);

    silex::FactoredElement uncached_root(field);
    is_power = false;
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, uncached_root, order, square, 2, precision, precision,
            nullptr, nullptr));
    assert(is_power);
    silex::Element cached_value(field);
    silex::Element uncached_value(field);
    assert(second_root.evaluate(cached_value) &&
           uncached_root.evaluate(uncached_value) &&
           cached_value.equal(uncached_value));

    silex::FactoredElement other_root(field);
    is_power = false;
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, other_root, order, other_square, 2, precision,
            precision, nullptr, &cache));
    assert(is_power && compact_power_root_squares_to(other_root, other_square));
    assert(cache.entries.size() >= filled);
    return 0;
}

// A unit whose formal leftover is not a unit: in Q(sqrt 2) with
// eps = 1 + sqrt 2, gamma = eps^K sqrt 2 and delta = eps^e / 2 are not units,
// but gamma^2 delta = eps^(2K + e) is.  The formal split leaves the small
// non-unit delta, so the residual conjugate-log bound is below the square
// root of the original one.
struct ResidualNotUnitCounter {
    slong events = 0;
};

void residual_not_unit_callback(void* user,
                                silex::DiagnosticsModule,
                                silex::ProfileEvent event,
                                const char*,
                                const char* label) noexcept {
    if (event == silex::ProfileEvent::event && label != nullptr &&
        std::strcmp(label, "unit_group.compact_unit_residual_not_unit") == 0) {
        ++static_cast<ResidualNotUnitCounter*>(user)->events;
    }
}

int test_compact_unit_power_root_nonunit_small_leftover() {
    const slong coefficients[] = {-2, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::Element sqrt2(field);
    silex::Element eps(field);
    assert(sqrt2.gen() && eps.add_si(sqrt2, 1));
    sflint::Fmpz k;
    fmpz_set_si(k.raw(), 20);
    silex::Element eps_k(field);
    silex::Element gamma(field);
    silex::Element half_eps(field);
    silex::Element half_eps_sq(field);
    silex::Element eps_sq(field);
    assert(eps_k.pow_fmpz(eps, sflint::FmpzConstRef(k)) &&
           gamma.multiply(eps_k, sqrt2) &&
           half_eps.scalar_div_si(eps, 2) &&
           eps_sq.multiply(eps, eps) &&
           half_eps_sq.scalar_div_si(eps_sq, 2));

    constexpr slong precision = 128;
    // eps^41 is not a square (eps is fundamental of norm -1).
    silex::FactoredElement odd(field);
    assert(odd.push(gamma, 2) && odd.push(half_eps, 1));
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    ResidualNotUnitCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            residual_not_unit_callback, &counter);
    bool is_power = true;
    silex::FactoredElement root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, odd, 2, precision, precision,
            &diagnostics, nullptr));
    assert(!is_power);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The residual branch was considered and rejected for the non-unit
    // leftover.
    assert(counter.events == 1);
#endif

    // eps^42 = (eps^21)^2.
    silex::FactoredElement even(field);
    assert(even.push(gamma, 2) && even.push(half_eps_sq, 1));
    is_power = false;
    silex::FactoredElement even_root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, even_root, order, even, 2, precision, precision));
    assert(is_power && compact_power_root_squares_to(even_root, even));
    return 0;
}

// The residual-leftover unit gate rejects a leftover whose factored norm is
// not +-1 without expanding it, and expands and tests every other leftover.
// In Q(sqrt 2) with eps = 1 + sqrt 2 and gamma = eps^20 sqrt 2:
// gamma^2 (eps / 2) = eps^41 leaves eps / 2 (norm -1/4), eps^41 leaves eps
// (norm -1) and eps^41 * eps^3 = eps^44 leaves eps * eps^3 (norm 1), and the
// residual bound selects the residual branch for all three.
struct ResidualLeftoverGateCounter {
    slong norm_rejected = 0;
    slong unit_tests = 0;
    slong not_unit = 0;
    slong residual_branch = 0;
    slong original_branch = 0;
};

void residual_leftover_gate_callback(void* user,
                                     silex::DiagnosticsModule,
                                     silex::ProfileEvent event,
                                     const char*,
                                     const char* label) noexcept {
    if (label == nullptr) {
        return;
    }
    auto* counter = static_cast<ResidualLeftoverGateCounter*>(user);
    if (event == silex::ProfileEvent::begin_scope &&
        std::strcmp(label, "unit_group.compact_unit_residual_unit_test") ==
                0) {
        ++counter->unit_tests;
    }
    if (event != silex::ProfileEvent::event) {
        return;
    }
    if (std::strcmp(label,
                    "unit_group.compact_unit_residual_norm_not_unit") == 0) {
        ++counter->norm_rejected;
    } else if (std::strcmp(label,
                           "unit_group.compact_unit_residual_not_unit") ==
               0) {
        ++counter->not_unit;
    } else if (std::strcmp(label,
                           "unit_group.compact_unit_residual_branch") == 0) {
        ++counter->residual_branch;
    } else if (std::strcmp(label,
                           "unit_group.compact_unit_original_branch") == 0) {
        ++counter->original_branch;
    }
}

int test_compact_unit_residual_leftover_norm_prefilter() {
    const slong coefficients[] = {-2, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::Element sqrt2(field);
    silex::Element eps(field);
    assert(sqrt2.gen() && eps.add_si(sqrt2, 1));
    sflint::Fmpz k;
    fmpz_set_si(k.raw(), 20);
    silex::Element eps_k(field);
    silex::Element gamma(field);
    silex::Element half_eps(field);
    silex::Element eps_sq(field);
    silex::Element eps_cube(field);
    assert(eps_k.pow_fmpz(eps, sflint::FmpzConstRef(k)) &&
           gamma.multiply(eps_k, sqrt2) &&
           half_eps.scalar_div_si(eps, 2) &&
           eps_sq.multiply(eps, eps) && eps_cube.multiply(eps_sq, eps));

    constexpr slong precision = 128;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    ResidualLeftoverGateCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            residual_leftover_gate_callback, &counter);

    // Non-unit leftover eps / 2: rejected by its norm, never expanded.
    silex::FactoredElement nonunit(field);
    assert(nonunit.push(gamma, 2) && nonunit.push(half_eps, 1));
    bool is_power = true;
    silex::FactoredElement root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, nonunit, 2, precision, precision,
            &diagnostics, nullptr));
    assert(!is_power);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.norm_rejected == 1 && counter.unit_tests == 0 &&
           counter.not_unit == 1 && counter.residual_branch == 0);
#endif

    // Unit leftover eps: passes the norm filter, is expanded and tested, and
    // takes the residual branch; eps^41 is not a square.
    counter = ResidualLeftoverGateCounter{};
    silex::FactoredElement odd(field);
    assert(odd.push(eps, 41));
    is_power = true;
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, odd, 2, precision, precision,
            &diagnostics, nullptr));
    assert(!is_power);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.norm_rejected == 0 && counter.unit_tests == 1 &&
           counter.not_unit == 0 && counter.residual_branch == 1);
#endif

    // Unit leftover eps * eps^3: the same route, and eps^44 = (eps^22)^2.
    counter = ResidualLeftoverGateCounter{};
    silex::FactoredElement even(field);
    assert(even.push(eps, 41) && even.push(eps_cube, 1));
    is_power = false;
    silex::FactoredElement even_root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, even_root, order, even, 2, precision, precision,
            &diagnostics, nullptr));
    assert(is_power && compact_power_root_squares_to(even_root, even));
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.norm_rejected == 0 && counter.unit_tests == 1 &&
           counter.not_unit == 0 && counter.residual_branch == 1);
#endif
    return 0;
}

// A leftover of norm +-1 that is not a unit of a non-maximal order.  In
// Z[2 sqrt 2] (the equation order of x^2 - 8, theta = 2 sqrt 2), eps =
// 1 + sqrt 2 has norm -1 but lies outside the order, while eps^2 =
// 3 + theta is a unit of it.  With n = 3, eps^28 * (eps^2)^k splits as
// (eps^9)^3 * eps * (eps^2)^k.  For k = 1 the leftover is eps^3 =
// 7 + 5 sqrt 2, of norm -1 and outside the order, and eps^30 = (eps^10)^3;
// for k = 0 the leftover is eps and eps^28 is not a cube.  The residual
// bound selects the residual branch, the norm prefilter passes the
// leftover, and the exact test must reject it, so the input takes the
// original branch.  A "norm +-1 implies unit" shortcut would send the
// leftover to the unit-only reduction, whose invariant needs a unit of the
// order.
int test_compact_unit_residual_leftover_norm_one_nonmaximal() {
    const slong coefficients[] = {-8, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const silex::Order order = silex::test::equation_order(field);
    assert(!order.is_maximal());
    silex::Element theta(field);
    silex::Element eps(field);
    silex::Element eps_sq(field);
    assert(theta.gen() && eps.scalar_div_si(theta, 2) &&
           eps.add_si(eps, 1) && eps_sq.multiply(eps, eps));
    constexpr slong precision = 128;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    ResidualLeftoverGateCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            residual_leftover_gate_callback, &counter);

    silex::FactoredElement cube(field);
    assert(cube.push(eps, 28) && cube.push(eps_sq, 1));
    bool is_power = false;
    silex::FactoredElement root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, cube, 3, precision, precision,
            &diagnostics, nullptr));
    assert(is_power);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The norm prefilter passes the leftover, the exact test rejects it,
    // and the input takes the original branch.
    assert(counter.norm_rejected == 0 && counter.unit_tests == 1 &&
           counter.not_unit == 1 && counter.residual_branch == 0 &&
           counter.original_branch == 1);
#endif
    silex::Element root_value(field);
    silex::Element root_cube(field);
    silex::Element value(field);
    assert(root.evaluate(root_value) &&
           root_cube.multiply(root_value, root_value) &&
           root_cube.multiply(root_cube, root_value) && cube.evaluate(value) &&
           root_cube.equal(value));

    counter = ResidualLeftoverGateCounter{};
    silex::FactoredElement non_cube(field);
    assert(non_cube.push(eps, 28));
    is_power = true;
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, non_cube, 3, precision, precision,
            &diagnostics, nullptr));
    assert(!is_power);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.norm_rejected == 0 && counter.unit_tests == 1 &&
           counter.not_unit == 1 && counter.residual_branch == 0 &&
           counter.original_branch == 1);
#endif
    return 0;
}

// Log-cache events seen inside the residual log bound and inside the compact
// presentation of a unit power test.  The presentation embeds the bases of
// its input first and only later the multipliers it builds, so its leading
// hits (those before its first store) are base logs served by the cache.
struct SharedLogCacheCounters {
    slong residual_bound_depth = 0;
    slong presentation_depth = 0;
    slong residual_bound_hits = 0;
    slong residual_bound_stores = 0;
    slong presentation_leading_hits = 0;
    slong presentation_stores = 0;
};

void shared_log_cache_callback(void* user,
                               silex::DiagnosticsModule,
                               silex::ProfileEvent event,
                               const char*,
                               const char* label) noexcept {
    SharedLogCacheCounters* counters =
            static_cast<SharedLogCacheCounters*>(user);
    if (label == nullptr) {
        return;
    }
    const bool residual_bound =
            std::strcmp(label, "unit_group.compact_unit_residual_log_bound") ==
            0;
    const bool presentation =
            std::strcmp(label,
                        "unit_group.compact_unit_presentation_reduction") == 0;
    if (event == silex::ProfileEvent::begin_scope ||
        event == silex::ProfileEvent::end_scope) {
        const slong step = event == silex::ProfileEvent::begin_scope ? 1 : -1;
        if (residual_bound) {
            counters->residual_bound_depth += step;
        } else if (presentation) {
            counters->presentation_depth += step;
        }
        return;
    }
    const bool hit =
            std::strcmp(label, "element.logarithmic_embedding.cache_hit") == 0;
    const bool store =
            std::strcmp(label, "element.logarithmic_embedding.cache_store") ==
            0;
    if (counters->residual_bound_depth > 0) {
        counters->residual_bound_hits += hit ? 1 : 0;
        counters->residual_bound_stores += store ? 1 : 0;
    }
    if (counters->presentation_depth > 0) {
        if (hit && counters->presentation_stores == 0) {
            ++counters->presentation_leading_hits;
        }
        counters->presentation_stores += store ? 1 : 0;
    }
}

// A unit power test whose bases have huge coefficients: one conjugate of
// each base is tiny, so its logarithm cancels and the factored-log precision
// ladder climbs far above the requested 64 bits.  The two log bounds and the
// compact presentation share one embedding context, so the later embeddings
// take each base's log from the context's cache, stored at the top rung of
// the first ladder.  The root must still square to the input.
int test_compact_unit_power_root_reuses_high_precision_base_logs() {
    const slong coefficients[] = {-2, 0, 1};
    const silex::NumberField field = monic_field(coefficients);
    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));

    silex::Element x(field);
    silex::Element one(field);
    silex::Element unit(field);
    assert(x.gen() && one.one() && unit.add(x, one));
    const auto unit_power = [&](slong exponent) {
        silex::Element out(field);
        sflint::Fmpz e;
        fmpz_set_si(e.raw(), exponent);
        assert(out.pow_fmpz(unit, sflint::FmpzConstRef(e)));
        return out;
    };
    const silex::Element large = unit_power(401);
    const silex::Element small = unit_power(3);
    const silex::Element structural = unit_power(150);

    // (1 + x)^401 * (1 + x)^3 * ((1 + x)^150)^2 = ((1 + x)^352)^2, with the
    // formal leftover (1 + x)^401 * (1 + x)^3.
    silex::FactoredElement square(field);
    assert(square.push(large, 1) && square.push(small, 1) &&
           square.push(structural, 2));

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    SharedLogCacheCounters counters;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::element) |
                    silex::diagnostics_module_bit(
                            silex::DiagnosticsModule::unit_group),
            shared_log_cache_callback, &counters);

    constexpr slong precision = 128;
    bool is_power = false;
    silex::FactoredElement root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, root, order, square, 2, precision, precision,
            &diagnostics, nullptr));
    assert(is_power && compact_power_root_squares_to(root, square));
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The residual bound takes both leftover base logs from the context the
    // original bound filled and stores none.  The compact presentation shares
    // that context, so its first embeddings, those of the bases of its input
    // (at least the two leftover bases), are cache hits; with a context of
    // its own the first embedding would miss and store.
    assert(counters.residual_bound_depth == 0 &&
           counters.presentation_depth == 0);
    assert(counters.residual_bound_hits >= 2 &&
           counters.residual_bound_stores == 0);
    assert(counters.presentation_leading_hits >= 2);
#else
    assert(counters.residual_bound_hits == 0 &&
           counters.presentation_leading_hits == 0);
#endif

    // A unit that is not a square: (1 + x)^401 * (1 + x)^4.
    silex::FactoredElement non_square(field);
    assert(non_square.push(large, 1) && non_square.push(unit_power(4), 1));
    is_power = true;
    silex::FactoredElement no_root(field);
    assert(silex::detail::compact_unit_presentation_power_root(
            is_power, no_root, order, non_square, 2, precision, precision,
            nullptr, nullptr));
    assert(!is_power);
    return 0;
}

// Reference evaluation of the modular power products: one independent binary
// powering per (output, base) and a running product, the evaluation that the
// shared-table interleaved evaluation replaces.
void reference_power_products(
        sflint::NmodMat& out,
        const std::vector<silex::detail::ModularPowerProductFactor>& factors,
        const sflint::NmodPoly& modulus) {
    const ulong prime = modulus.raw()->mod.n;
    for (slong column = 0; column < nmod_mat_nrows(out.raw()); ++column) {
        sflint::NmodPoly value(prime);
        sflint::NmodPoly power(prime);
        sflint::NmodPoly product(prime);
        sflint::Fmpz magnitude;
        nmod_poly_one(value.raw());
        for (const auto& factor : factors) {
            const sflint::Fmpz& exponent =
                    (*factor.exponents)[static_cast<std::size_t>(column)];
            if (fmpz_is_zero(exponent.raw()) != 0) {
                continue;
            }
            const sflint::NmodPoly& base = fmpz_sgn(exponent.raw()) < 0
                                                   ? *factor.inverse
                                                   : *factor.base;
            fmpz_abs(magnitude.raw(), exponent.raw());
            if (fmpz_abs_fits_ui(magnitude.raw()) != 0) {
                nmod_poly_powmod_ui_binexp(power.raw(), base.raw(),
                                           fmpz_get_ui(magnitude.raw()),
                                           modulus.raw());
            } else {
                nmod_poly_powmod_fmpz_binexp(power.raw(), base.raw(),
                                             magnitude.raw(), modulus.raw());
            }
            nmod_poly_mulmod(product.raw(), value.raw(), power.raw(),
                             modulus.raw());
            value.swap(product);
        }
        for (slong i = 0; i < nmod_mat_ncols(out.raw()); ++i) {
            nmod_mat_entry(out.raw(), column, i) =
                    nmod_poly_get_coeff_ui(value.raw(), i);
        }
    }
}

struct PowerProductCase {
    std::vector<sflint::NmodPoly> bases;
    std::vector<sflint::NmodPoly> inverses;
    std::vector<std::vector<sflint::Fmpz>> exponents;
    std::vector<bool> has_inverse;
};

// Bases are random reduced residues that are invertible modulo the modulus
// when allow_negative is set; exponent bit sizes are drawn up to max_bits
// (zero included), with random signs when allow_negative is set.
void random_power_product_case(PowerProductCase& out, flint_rand_t state,
                               const sflint::NmodPoly& modulus,
                               std::size_t base_count,
                               std::size_t output_count,
                               flint_bitcnt_t max_bits, bool allow_negative) {
    const ulong prime = modulus.raw()->mod.n;
    const slong length = nmod_poly_length(modulus.raw());
    for (std::size_t k = 0; k < base_count; ++k) {
        out.bases.emplace_back(prime);
        out.inverses.emplace_back(prime);
        while (true) {
            nmod_poly_randtest(out.bases.back().raw(), state, length - 1);
            if (!allow_negative ||
                nmod_poly_invmod(out.inverses.back().raw(),
                                 out.bases.back().raw(),
                                 modulus.raw()) != 0) {
                break;
            }
        }
        out.exponents.emplace_back(output_count);
        bool negative = false;
        for (sflint::Fmpz& exponent : out.exponents.back()) {
            const flint_bitcnt_t bits = n_randint(state, max_bits + 1);
            if (bits != 0) {
                fmpz_randbits(exponent.raw(), state, bits);
            }
            if (!allow_negative) {
                fmpz_abs(exponent.raw(), exponent.raw());
            }
            negative = negative || fmpz_sgn(exponent.raw()) < 0;
        }
        out.has_inverse.push_back(negative);
    }
}

void assert_power_products_match_reference(const PowerProductCase& input,
                                           const sflint::NmodPoly& modulus,
                                           std::size_t output_count) {
    std::vector<silex::detail::ModularPowerProductFactor> factors;
    for (std::size_t k = 0; k < input.bases.size(); ++k) {
        factors.push_back({&input.bases[k],
                           input.has_inverse[k] ? &input.inverses[k] : nullptr,
                           &input.exponents[k]});
    }
    const ulong prime = modulus.raw()->mod.n;
    const slong degree = nmod_poly_degree(modulus.raw());
    sflint::NmodMat expected(static_cast<slong>(output_count), degree, prime);
    sflint::NmodMat actual(static_cast<slong>(output_count), degree, prime);
    reference_power_products(expected, factors, modulus);
    silex::detail::modular_power_products(actual, factors, modulus);
    assert(nmod_mat_equal(expected.raw(), actual.raw()) != 0);
}

void random_monic_modulus(sflint::NmodPoly& modulus, flint_rand_t state,
                          slong degree) {
    nmod_poly_randtest_monic(modulus.raw(), state, degree + 1);
}

int test_modular_power_products_match_independent_powering() {
    flint_rand_t state;
    flint_rand_init(state);
    const ulong primes[] = {UWORD(5), UWORD(101),
                            n_nextprime(UWORD(1) << 40, 1),
                            n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1)};
    struct Shape {
        slong degree;
        std::size_t bases;
        std::size_t outputs;
        flint_bitcnt_t max_bits;
    };
    // A single base, a single output, word-size and wide exponents, many
    // bases (the wide unit rebuild shape), and degree one.
    const Shape shapes[] = {{1, 1, 1, 3},   {1, 4, 3, 70},  {2, 1, 5, 1},
                            {3, 1, 1, 130}, {5, 7, 6, 64},  {8, 3, 9, 200},
                            {32, 59, 4, 73}, {6, 40, 12, 9}, {4, 12, 2, 66}};
    for (ulong prime : primes) {
        for (const Shape& shape : shapes) {
            for (int sign_mode = 0; sign_mode < 2; ++sign_mode) {
                sflint::NmodPoly modulus(prime);
                random_monic_modulus(modulus, state, shape.degree);
                PowerProductCase input;
                random_power_product_case(input, state, modulus, shape.bases,
                                          shape.outputs, shape.max_bits,
                                          sign_mode == 1);
                assert_power_products_match_reference(input, modulus,
                                                       shape.outputs);
            }
        }
    }

    // Edge cases on one modulus: an output whose exponents are all zero, a
    // zero base, the base one, exponents +-1, and magnitudes at the word
    // boundary (2^64 - 1 and 2^64) and of a single high bit.
    const ulong prime = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
    sflint::NmodPoly modulus(prime);
    random_monic_modulus(modulus, state, 6);
    PowerProductCase input;
    random_power_product_case(input, state, modulus, 4, 6, 8, true);
    for (auto& row : input.exponents) {
        fmpz_zero(row[0].raw());
    }
    fmpz_set_si(input.exponents[0][1].raw(), 1);
    fmpz_set_si(input.exponents[1][1].raw(), -1);
    fmpz_set_ui(input.exponents[0][2].raw(), UWORD_MAX);
    fmpz_set_ui(input.exponents[1][2].raw(), UWORD_MAX);
    fmpz_neg(input.exponents[1][2].raw(), input.exponents[1][2].raw());
    fmpz_one(input.exponents[2][3].raw());
    fmpz_mul_2exp(input.exponents[2][3].raw(), input.exponents[2][3].raw(),
                  FLINT_BITS);
    fmpz_neg(input.exponents[3][3].raw(), input.exponents[2][3].raw());
    fmpz_set_si(input.exponents[3][4].raw(), -64);
    input.has_inverse.assign(input.bases.size(), false);
    for (std::size_t k = 0; k < input.bases.size(); ++k) {
        for (const sflint::Fmpz& exponent : input.exponents[k]) {
            if (fmpz_sgn(exponent.raw()) < 0) {
                input.has_inverse[k] = true;
            }
        }
    }
    assert_power_products_match_reference(input, modulus, 6);

    PowerProductCase special;
    for (int k = 0; k < 2; ++k) {
        special.bases.emplace_back(prime);
        special.inverses.emplace_back(prime);
        special.exponents.emplace_back(3);
        special.has_inverse.push_back(false);
    }
    nmod_poly_zero(special.bases[0].raw());
    nmod_poly_one(special.bases[1].raw());
    fmpz_set_ui(special.exponents[0][0].raw(), 7);
    fmpz_set_ui(special.exponents[1][0].raw(), 5);
    fmpz_set_ui(special.exponents[1][1].raw(), UWORD_MAX);
    fmpz_set_ui(special.exponents[0][2].raw(), 1);
    assert_power_products_match_reference(special, modulus, 3);

    // No factors at all: every output is one.
    PowerProductCase empty;
    assert_power_products_match_reference(empty, modulus, 2);

    flint_rand_clear(state);
    return 0;
}

}  // namespace

int main() {
    assert(test_modular_power_products_match_independent_powering() == 0);
    assert(test_shared_exponent_matrix_and_centered_crt() == 0);
    assert(test_denominator_prime_exclusion() == 0);
    assert(test_field_polynomial_denominator_prime_exclusion() == 0);
    assert(test_noninvertible_prime_exclusion() == 0);
    assert(test_order_basis_prime_exclusion() == 0);
    assert(test_bound_violation_fails_without_publication() == 0);
    assert(test_single_prime_centered_negative_lift() == 0);
    assert(test_empty_batch() == 0);
    assert(test_coordinate_bound_and_reconstruction() == 0);
    assert(test_coordinate_bound_rank_zero() == 0);
    assert(test_coordinate_bound_rejects_nonmaximal_order() == 0);
    assert(test_coordinate_bound_supports_reordered_basis() == 0);
    assert(test_coordinate_bound_supports_index_three_maximal_basis() == 0);
    assert(test_compact_field_modulus_cache_keys() == 0);
    assert(test_compact_field_modulus_cache_shared_across_calls() == 0);
    assert(test_compact_unit_power_root_reuses_high_precision_base_logs() == 0);
    assert(test_compact_unit_power_root_nonunit_small_leftover() == 0);
    assert(test_compact_unit_residual_leftover_norm_prefilter() == 0);
    assert(test_compact_unit_residual_leftover_norm_one_nonmaximal() == 0);
    assert(test_compact_power_input_norm() == 0);
    return 0;
}
