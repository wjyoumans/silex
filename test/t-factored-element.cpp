#include <silex/factored_element.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>
#include <limits>
#include <utility>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

void set_rational(silex::Element& element,
                  slong numerator,
                  ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator, denominator);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool element_is_rational(const silex::Element& element,
                         slong numerator,
                         ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq expected;
    sflint::Fmpq actual;
    sflint::fmpq_set_si(expected, numerator, denominator);
    assert(element.get_fmpq_poly(sflint::FmpqPolyRef(polynomial)));
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual), polynomial, 0);
    return sflint::fmpq_poly_degree(polynomial) <= 0 &&
           sflint::fmpq_equal(actual, expected);
}

void set_one_plus_theta(silex::Element& element) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool has_radius_at_most_2exp(const arb_struct* value,
                             slong exponent) noexcept {
    return ::arb_is_finite(value) != 0 &&
           ::mag_cmp_2exp_si(arb_radref(value), exponent) <= 0;
}

}  // namespace

int main() {
    silex::NumberField field = degree_one_field();
    silex::NumberField other_field = degree_one_field();

    silex::Element two(field);
    silex::Element three(field);
    silex::Element four(field);
    silex::Element eight(field);
    silex::Element zero(field);
    silex::Element other_parent(other_field);
    assert(two.set_si(2));
    assert(three.set_si(3));
    assert(four.set_si(4));
    assert(eight.set_si(8));
    assert(zero.zero());
    assert(other_parent.set_si(2));

    silex::FactoredElement product(field);
    assert(product.is_defined());
    assert(product.parent() != nullptr &&
           product.parent()->has_same_data(field));
    assert(product.length() == 0);

    slong exponent = 3;
    assert(product.push(two, exponent));
    exponent = -1;
    assert(product.push(three, exponent));
    assert(product.length() == 2);

    silex::Element evaluated(field);
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement recursive_eval(field);
    assert(recursive_eval.push(two, 13));
    assert(recursive_eval.push(three, -11));
    assert(recursive_eval.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8192, 177147));

    silex::FactoredElement single(field);
    exponent = 2;
    assert(single.push(two, exponent));

    silex::FactoredElement combined(field);
    assert(combined.multiply(product, single));
    assert(combined.length() == 3);
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 3));

    assert(combined.divide(combined, single));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    assert(combined.invert(product));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 8));

    assert(combined.multiply(combined, product));
    assert(combined.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    exponent = -2;
    assert(combined.pow_si(product, exponent));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 9, 64));

    exponent = 0;
    assert(combined.pow_si(product, exponent));
    assert(combined.length() == 0);
    assert(combined.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    silex::FactoredElement copy(field);
    assert(copy.set(product));
    assert(copy.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement moved(std::move(copy));
    assert(moved.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    assert(!product.push(zero, exponent));
    assert(!product.push(other_parent, exponent));

    silex::FactoredElement other_product(other_field);
    assert(other_product.set_element(other_parent));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.multiply(product, other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.divide(product, other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.invert(other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.pow_si(other_product, exponent));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement root_input(field);
    assert(root_input.push(two, 6));
    assert(root_input.push(three, -4));
    assert(root_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    silex::FactoredElement root(field);
    assert(root.root_si(root_input, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 9));
    assert(combined.pow_si(root, 2));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    assert(root.root_si(root_input, 1));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    assert(root.root_si(root, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 9));

    silex::FactoredElement empty_root(field);
    assert(empty_root.root_si(empty_root, 7));
    assert(empty_root.length() == 0);
    assert(empty_root.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    assert(root.set_element(three));
    assert(!root.root_si(root_input, 5));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(root_input, 0));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(root_input, -2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(other_product, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));

    bool is_square = false;
    silex::FactoredElement square_input(field);
    silex::FactoredElement square_root(field);
    silex::FactoredElement square_power(field);
    assert(square_input.push(two, 4));
    assert(square_input.push(three, -2));
    assert(square_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));

    assert(square_input.is_square(is_square, square_input));
    assert(is_square);
    assert(square_power.pow_si(square_input, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));

    assert(square_input.one());
    assert(square_input.push(four, 1));
    assert(square_input.push(three, 2));
    assert(square_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 36, 1));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 36, 1));

    assert(square_input.one());
    assert(square_input.push(four, -1));
    assert(square_input.push(three, 2));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 2));
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 9, 4));

    assert(square_input.one());
    assert(square_input.push(two, 1));
    assert(square_root.set_element(four));
    is_square = true;
    assert(square_input.is_square(is_square, square_root));
    assert(!is_square);
    assert(square_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(square_input.one());
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_root.length() == 0);
    assert(!square_input.is_square(is_square, other_product));

    bool is_power = false;
    silex::FactoredElement power_input(field);
    silex::FactoredElement power_root(field);
    silex::FactoredElement power_check(field);
    assert(power_input.push(two, 6));
    assert(power_input.push(three, -3));
    assert(power_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));
    assert(power_input.is_power_si(is_power, power_root, 3,
                                   silex::FactoredRootStrategy::reduced));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));

    assert(power_input.is_power_si(is_power, power_input, 3));
    assert(is_power);
    assert(power_check.pow_si(power_input, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));

    assert(power_input.one());
    assert(power_input.push(four, 1));
    assert(power_input.is_power_si(is_power, power_root, 2));
    assert(is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 2, 1));

    assert(power_input.one());
    assert(power_input.push(two, 3));
    assert(power_input.push(three, 6));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));

    assert(power_input.one());
    assert(power_input.push(eight, -1));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 1, 8));

    assert(power_input.one());
    assert(power_input.push(two, 1));
    assert(power_root.set_element(four));
    is_power = true;
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(!is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(power_input.one());
    assert(power_input.is_power_si(is_power, power_input, 1));
    assert(is_power);
    assert(power_input.length() == 0);

    assert(power_root.set_element(four));
    is_power = true;
    assert(!power_input.is_power_si(is_power, power_root, 0));
    assert(is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    assert(!power_input.is_power_si(is_power, other_product, 3));

    silex::CompactElement compact(field);
    assert(compact.is_defined());
    assert(compact.parent() != nullptr &&
           compact.parent()->has_same_data(field));
    assert(compact.one(3));
    assert(compact.base() == 3);
    assert(compact.length() == 0);
    assert(compact.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    silex::FactoredElement compact_input(field);
    silex::FactoredElement compact_root(field);
    silex::FactoredElement compact_power(field);
    assert(compact_input.push(two, 3));
    assert(compact_input.push(three, 6));
    assert(compact_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(compact_root.length() == 2);
    slong compact_exp = 0;
    assert(compact_root.exponent(compact_exp, 0));
    assert(compact_exp == 1);
    assert(compact_root.exponent(compact_exp, 1));
    assert(compact_exp == 2);
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact.length() == 2);
    assert(compact.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));
    assert(compact.root_base(compact_root));
    assert(compact_power.pow_si(compact_root, 3));
    assert(compact_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));

    assert(compact_input.one());
    assert(compact_input.push(two, 4));
    assert(compact_input.push(four, 1));
    assert(!compact_root.root_si(compact_input, 3));
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    assert(compact_power.pow_si(compact_root, 3));
    assert(compact_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 1));

    assert(compact_input.one());
    assert(compact_input.push(two, 5));
    assert(compact_input.push(three, -4));
    assert(compact_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 81));
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 81));

    assert(compact_input.one());
    assert(compact_input.push(two, 1));
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact_root.set_element(four));
    assert(!compact.root_base(compact_root));
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    is_power = true;
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(!is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(compact_input.one());
    assert(!compact.set_factored_element(compact_input, 1));
    assert(!compact_input.is_power_si(is_power, compact_root, -3));
    assert(!compact.root_base(other_product));
    is_power = true;
    assert(!compact_input.is_power_si(
            is_power, compact_root, 3,
            static_cast<silex::FactoredRootStrategy>(99)));
    assert(is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    silex::FactoredElement overflow_input(field);
    assert(overflow_input.push(two, std::numeric_limits<slong>::min()));
    assert(combined.set_element(three));
    assert(!combined.invert(overflow_input));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!combined.pow_si(overflow_input, -1));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!combined.divide(product, overflow_input));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));

    assert(product.one());
    assert(product.length() == 0);
    assert(product.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    set_rational(two, 5, 2);
    assert(product.set_element(two));
    assert(product.length() == 1);
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 5, 2));
    assert(!product.set_element(zero));

    silex::EmbeddingContext rational_embeddings(field);
    silex::Element log_two(field);
    silex::Element log_four(field);
    silex::FactoredElement guarded_logs(field);
    sflint::ArbVec guarded_log_row(1);
    const slong large_even_exponent =
            (std::numeric_limits<slong>::max() / 2048) * 2;
    assert(log_two.set_si(2));
    assert(log_four.set_si(4));
    assert(guarded_logs.push(log_two, large_even_exponent));
    assert(guarded_logs.push(log_four, -large_even_exponent / 2));
    assert(guarded_logs.logarithmic_embedding(
            sflint::ArbVecRef(guarded_log_row), rational_embeddings,
            silex::LogEmbeddingMode::product, 32));
    assert(sflint::arb_contains_si(guarded_log_row.data() + 0, 0));
    assert(has_radius_at_most_2exp(guarded_log_row.data() + 0, -32));

    silex::NumberField quadratic = quadratic_field(2);
    silex::NumberField other_quadratic = quadratic_field(3);
    silex::EmbeddingContext embeddings(quadratic);
    silex::EmbeddingContext other_embeddings(other_quadratic);
    silex::Element unit_like(quadratic);
    silex::Element expanded(quadratic);
    sflint::ArbVec elt_logs(2);
    set_one_plus_theta(unit_like);

    silex::FactoredElement quadratic_square(quadratic);
    silex::FactoredElement quadratic_root(quadratic);
    silex::FactoredElement quadratic_power(quadratic);
    silex::Element quadratic_evaluated(quadratic);
    assert(expanded.multiply(unit_like, unit_like));
    assert(quadratic_square.push(expanded, 1));
    assert(quadratic_square.is_square(is_square, quadratic_root));
    assert(is_square);
    assert(quadratic_power.pow_si(quadratic_root, 2));
    assert(quadratic_power.evaluate(quadratic_evaluated));
    assert(quadratic_evaluated.equal(expanded));

    silex::NumberField cubic = cubic_field();
    silex::Element cubic_eight(cubic);
    silex::Element cubic_evaluated(cubic);
    silex::FactoredElement cubic_power(cubic);
    silex::FactoredElement cubic_root(cubic);
    silex::FactoredElement cubic_check(cubic);
    silex::FactoredElement cubic_compact_root(cubic);
    silex::Element cubic_theta(cubic);
    silex::Element cubic_shift(cubic);
    silex::Element cubic_cube(cubic);
    assert(cubic_eight.set_si(8));
    assert(cubic_power.push(cubic_eight, 1));
    assert(cubic_power.is_power_si(is_power, cubic_root, 3));
    assert(is_power);
    assert(cubic_root.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal_si(2));
    assert(cubic_check.pow_si(cubic_root, 3));
    assert(cubic_check.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal(cubic_eight));
    assert(cubic_theta.gen());
    assert(cubic_shift.add_si(cubic_theta, 1));
    sflint::Fmpz cubic_exponent;
    sflint::fmpz_set_ui(cubic_exponent, 3);
    assert(cubic_cube.pow_fmpz(cubic_shift,
                               sflint::FmpzConstRef(cubic_exponent)));
    assert(cubic_power.one());
    assert(cubic_power.push(cubic_cube, 1));
    is_power = false;
    assert(cubic_power.is_power_si(
            is_power, cubic_compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(cubic_check.pow_si(cubic_compact_root, 3));
    assert(cubic_check.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal(cubic_cube));

    silex::FactoredElement compact_log_source(quadratic);
    silex::CompactElement compact_logs(quadratic);
    silex::Element compact_log_value(quadratic);
    sflint::ArbVec compact_log_row(2);
    assert(compact_log_source.push(unit_like, 5));
    assert(compact_logs.set_factored_element(compact_log_source, 3));
    assert(compact_logs.evaluate(compact_log_value));
    assert(compact_logs.logarithmic_embedding(
            sflint::ArbVecRef(compact_log_row), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(silex::logarithmic_embedding(
            sflint::ArbVecRef(elt_logs), embeddings, compact_log_value,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_overlaps(
            sflint::ArbConstRef(compact_log_row.data() + 0),
            sflint::ArbConstRef(elt_logs.data() + 0)));
    assert(sflint::arb_overlaps(
            sflint::ArbConstRef(compact_log_row.data() + 1),
            sflint::ArbConstRef(elt_logs.data() + 1)));

    silex::FactoredElement empty_logs(quadratic);
    sflint::ArbVec fac_logs(2);
    assert(empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 0));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 0));

    exponent = 2;
    assert(empty_logs.push(unit_like, exponent));
    assert(empty_logs.evaluate(expanded));
    assert(empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(silex::logarithmic_embedding(
            sflint::ArbVecRef(elt_logs), embeddings, expanded,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_overlaps(sflint::ArbConstRef(fac_logs.data() + 0),
                                sflint::ArbConstRef(elt_logs.data() + 0)));
    assert(sflint::arb_overlaps(sflint::ArbConstRef(fac_logs.data() + 1),
                                sflint::ArbConstRef(elt_logs.data() + 1)));

    sflint::arb_set_si(fac_logs.data() + 0, 123);
    sflint::arb_set_si(fac_logs.data() + 1, 456);
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), other_embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            static_cast<silex::LogEmbeddingMode>(99), 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 0));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));

    return 0;
}
