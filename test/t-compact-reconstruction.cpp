#include <silex/factored_element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include "factored_element/compact_reconstruction_internal.hpp"
#include "order_unit/compact_reconstruction_bound_internal.hpp"
#include "test_support.hpp"

#include <cassert>
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
    silex::Order reordered =
            silex::Order::from_basis(field, sflint::FmpqMatConstRef(basis));
    assert(reordered.is_defined());
    reordered.set_maximality(true);
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

}  // namespace

int main() {
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
    return 0;
}
