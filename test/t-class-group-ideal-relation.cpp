#include <silex/class_group.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/prime_ideal.hpp>

#include "class_group/class_group_internal.hpp"
#include "class_group/lll_relation_search_internal.hpp"
#include "class_group/relation_search_internal.hpp"
#include "ideal/two_generator_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstring>
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

bool set_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

silex::PrimeIdeal prime_above(const silex::Order& order, slong p) noexcept {
    sflint::Fmpz prime;
    set_si(prime, p);
    silex::PrimeIdealList factors;
    assert(silex::decompose_prime(
            factors, order, sflint::FmpzConstRef(prime)));
    assert(factors.size() >= 1 && factors.at(0) != nullptr);
    silex::PrimeIdeal result(order);
    assert(result.set(*factors.at(0)));
    return result;
}

struct Fixture {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context{order};

    Fixture() noexcept {
        sflint::Fmpz bound;
        set_si(bound, 2);
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
        assert(context.factor_base() != nullptr);
    }
};

void assert_context_unchanged(const silex::ClassGroupContext& context,
                              slong relation_count,
                              bool had_presentation,
                              silex::CertificationMode certification) noexcept {
    assert(context.relation_count() == relation_count);
    assert(context.has_presentation() == had_presentation);
    assert(context.certification_status() == certification);
}

int test_direct_factor_base_witness() {
    Fixture fixture;
    const silex::FactorBase* base = fixture.context.factor_base();
    assert(base->length() == 1);

    silex::PrimeIdeal prime(fixture.order);
    silex::Ideal ideal(fixture.order);
    assert(base->prime(prime, 0));
    assert(prime.get_ideal(ideal));

    silex::FactoredElement multiplier(fixture.field);
    sflint::FmpzMat row(1, base->length());
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(row), 0, 0),
                        91);
    silex::detail::ClassGroupIdealRelationWitnessResult result;
    const slong relation_count = fixture.context.relation_count();
    const bool had_presentation = fixture.context.has_presentation();
    const auto certification = fixture.context.certification_status();
    assert(silex::detail::class_group_ideal_relation_witness(
            result, multiplier, sflint::FmpzMatRef(row), fixture.context,
            ideal));
    assert(result.status ==
           silex::detail::ClassGroupIdealRelationWitnessStatus::success);
    assert(result.stage ==
           silex::detail::ClassGroupIdealRelationWitnessStage::
                   direct_factorization);
    assert(!result.used_reduction);
    assert(silex::detail::verify_class_group_ideal_relation_witness(
            fixture.context, ideal, multiplier,
            sflint::FmpzMatConstRef(row)));
    assert_context_unchanged(fixture.context, relation_count,
                             had_presentation, certification);

    sflint::fmpz_add_ui(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(row), 0, 0),
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(row), 0, 0),
                        1);
    assert(!silex::detail::verify_class_group_ideal_relation_witness(
            fixture.context, ideal, multiplier,
            sflint::FmpzMatConstRef(row)));
    return 0;
}

int test_reduced_selected_prime_outside_factor_base() {
    silex::NumberField field = silex::test::quadratic_field(-23);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    set_si(bound, 2);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr && base->length() > 0);
    silex::PrimeIdeal selected = prime_above(order, 3);
    assert(!base->contains(selected));
    silex::Ideal ideal(order);
    assert(selected.get_ideal(ideal));

    silex::FactoredElement multiplier(field);
    sflint::FmpzMat row(1, base->length());
    silex::detail::ClassGroupIdealRelationWitnessResult result;
    const slong relation_count = context.relation_count();
    const auto certification = context.certification_status();
    assert(silex::detail::class_group_ideal_relation_witness(
            result, multiplier, sflint::FmpzMatRef(row), context, ideal));
    assert(result.stage ==
           silex::detail::ClassGroupIdealRelationWitnessStage::
                   reduced_factorization);
    assert(result.used_reduction);
    assert(silex::detail::verify_class_group_ideal_relation_witness(
            context, ideal, multiplier,
            sflint::FmpzMatConstRef(row)));
    assert_context_unchanged(context, relation_count, false, certification);
    return 0;
}

int test_continuation_witness() {
    silex::NumberField field = silex::test::quadratic_field(-23);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    set_si(bound, 2);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr && base->length() > 0);
    silex::OrderElement generator(order);
    silex::Ideal ideal(order);
    assert(generator.set_si(3));
    assert(ideal.set_principal(generator));

    silex::Element one(field);
    assert(one.set_si(1));
    silex::FactoredElement multiplier(field);
    sflint::FmpzMat row(1, base->length());
    silex::detail::ClassGroupIdealRelationWitnessOptions options;
    options.max_candidates = 16;
    options.max_candidates_per_ideal = 16;
    options.max_random_products = 0;
    silex::detail::ClassGroupIdealRelationWitnessResult result;
    assert(silex::detail::class_group_ideal_relation_continuation_witness(
            result, multiplier, sflint::FmpzMatRef(row), context,
            ideal, ideal, one, options));
    assert(result.stage ==
           silex::detail::ClassGroupIdealRelationWitnessStage::
                   continuation_enumeration);
    assert(result.candidates_tried > 0);
    assert(!result.used_random_product);
    assert(silex::detail::verify_class_group_ideal_relation_witness(
            context, ideal, multiplier,
            sflint::FmpzMatConstRef(row)));
    return 0;
}

int test_exhaustion_is_named_and_transactional() {
    Fixture fixture;
    const silex::FactorBase* base = fixture.context.factor_base();
    silex::PrimeIdeal selected = prime_above(fixture.order, 3);
    silex::Ideal ideal(fixture.order);
    assert(selected.get_ideal(ideal));

    silex::Element one(fixture.field);
    assert(one.set_si(1));
    silex::Element sentinel(fixture.field);
    assert(sentinel.set_si(7));
    silex::FactoredElement multiplier(fixture.field);
    assert(multiplier.set_element(sentinel));
    sflint::FmpzMat row(1, base->length());
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(row), 0, 0),
                        73);

    silex::detail::ClassGroupIdealRelationWitnessOptions options;
    options.max_candidates = 0;
    options.max_random_products = 0;
    silex::detail::ClassGroupIdealRelationWitnessResult result;
    assert(!silex::detail::class_group_ideal_relation_continuation_witness(
            result, multiplier, sflint::FmpzMatRef(row), fixture.context,
            ideal, ideal, one, options));
    assert(result.status ==
           silex::detail::ClassGroupIdealRelationWitnessStatus::exhausted);
    assert(result.stage ==
           silex::detail::ClassGroupIdealRelationWitnessStage::
                   search_exhaustion);
    assert(std::strcmp(
                   silex::detail::class_group_ideal_relation_witness_status_name(
                           result.status),
                   "exhausted") == 0);
    assert(std::strcmp(
                   silex::detail::class_group_ideal_relation_witness_stage_name(
                           result.stage),
                   "search_exhaustion") == 0);
    assert(sflint::fmpz_equal_si(sflint::fmpz_mat_entry(
                                         sflint::FmpzMatConstRef(row), 0, 0),
                                 73));
    silex::Element expanded(fixture.field);
    assert(multiplier.evaluate(expanded));
    assert(expanded.equal(sentinel));
    return 0;
}

int test_random_product_row_is_exact() {
    silex::NumberField field = silex::test::quadratic_field(-23);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    set_si(bound, 5);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr && base->length() > 0);

    silex::FactoredElement one(field);
    assert(one.one());
    for (slong sequence_index : {WORD(0), WORD(3)}) {
        silex::Ideal random_product(order);
        sflint::FmpzMat row(1, base->length());
        assert(silex::detail::relation_search::
                       build_random_factor_base_product(
                               random_product, sflint::FmpzMatRef(row),
                               *base, sequence_index,
                               UWORD(0x6a09e667f3bcc909), nullptr));
        assert(silex::detail::verify_class_group_ideal_relation_witness(
                context, random_product, one,
                sflint::FmpzMatConstRef(row)));
    }
    return 0;
}

bool ideal_hnf_equal(const silex::Ideal& left,
                     const silex::Ideal& right) noexcept {
    const slong n = left.degree();
    sflint::FmpzMat left_hnf(n, n);
    sflint::FmpzMat right_hnf(n, n);
    return left.get_hnf(sflint::FmpzMatRef(left_hnf)) &&
           right.get_hnf(sflint::FmpzMatRef(right_hnf)) &&
           sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left_hnf),
                                  sflint::FmpzMatConstRef(right_hnf));
}

// The random-product search multiplies by a factor-base prime through
// multiply_by_factor_base_prime, which finds the prime's two-element form
// once and reuses it.  Its product must be the HNF Ideal::multiply gives,
// whichever factor Ideal::multiply would search.  Dedekind's cubic
// x^3 - x^2 - 2x - 8 has 2 as a common index divisor: 2 splits into three
// primes of norm 2, none with a Kummer generator, so their forms come from
// the search alone.  Inputs include O, every prime, products whose minimum
// (2) is below the prime's, and a square of a prime above 2.
int test_factor_base_prime_form_product_matches_ideal_multiply() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -8);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation_order));

    std::vector<silex::PrimeIdeal> primes;
    std::vector<silex::Ideal> prime_ideals;
    for (slong p : {WORD(2), WORD(3), WORD(5), WORD(7), WORD(11),
                    WORD(13)}) {
        sflint::Fmpz rational_prime;
        set_si(rational_prime, p);
        silex::PrimeIdealList factors;
        assert(silex::decompose_prime(
                factors, order, sflint::FmpzConstRef(rational_prime)));
        for (slong i = 0; i < factors.size(); ++i) {
            assert(factors.at(i) != nullptr);
            primes.emplace_back(order);
            assert(primes.back().set(*factors.at(i)));
            prime_ideals.emplace_back(order);
            assert(primes.back().get_ideal(prime_ideals.back()));
        }
    }
    // The three primes above 2 come first.
    assert(prime_ideals.size() > 3);
    sflint::FmpzMat coordinates(1, order.degree());
    for (std::size_t i = 0; i < 3; ++i) {
        sflint::Fmpz p;
        assert(primes.at(i).rational_prime(sflint::FmpzRef(p)));
        assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(p), 2));
        assert(!primes.at(i).kummer_generator_coordinates(
                sflint::FmpzMatRef(coordinates)));
    }

    std::vector<silex::Ideal> inputs;
    inputs.emplace_back(order);
    assert(inputs.back().one());
    for (const silex::Ideal& prime_ideal : prime_ideals) {
        inputs.emplace_back(order);
        assert(inputs.back().set(prime_ideal));
    }
    inputs.emplace_back(order);
    assert(inputs.back().multiply(prime_ideals.at(0), prime_ideals.at(1)));
    inputs.emplace_back(order);
    assert(inputs.back().multiply(prime_ideals.at(0), prime_ideals.at(0)));
    silex::Ideal partial(order);
    assert(partial.multiply(prime_ideals.at(2), prime_ideals.at(3)));
    inputs.emplace_back(order);
    assert(inputs.back().multiply(partial, prime_ideals.back()));

    std::vector<silex::detail::RandomProductPrimeForm> forms;
    slong with_form = 0;
    for (std::size_t i = 0; i < prime_ideals.size(); ++i) {
        silex::detail::IdealTwoGeneratorForm form;
        if (silex::detail::ideal_two_generator_form(form,
                                                    prime_ideals.at(i))) {
            ++with_form;
        }
        // Twice: the first pass builds forms[i], the second reuses it.
        for (int pass = 0; pass < 2; ++pass) {
            for (const silex::Ideal& input : inputs) {
                silex::Ideal expected(order);
                silex::Ideal actual(order);
                assert(expected.multiply(input, prime_ideals.at(i)));
                assert(silex::detail::relation_search::
                               multiply_by_factor_base_prime(
                                       actual, input, prime_ideals.at(i),
                                       static_cast<slong>(i), forms));
                assert(actual.equal(expected));
                assert(ideal_hnf_equal(actual, expected));
            }
        }
        assert(forms.size() == i + 1);
    }
    // Every prime above 2 is nonscalar and has a form.
    for (std::size_t i = 0; i < 3; ++i) {
        assert(forms.at(i).has_form);
    }
    assert(with_form >= 3);

    // Index 0 now holding another prime rebuilds the entry.
    for (const silex::Ideal& input : inputs) {
        silex::Ideal expected(order);
        silex::Ideal actual(order);
        assert(expected.multiply(input, prime_ideals.at(1)));
        assert(silex::detail::relation_search::multiply_by_factor_base_prime(
                actual, input, prime_ideals.at(1), 0, forms));
        assert(ideal_hnf_equal(actual, expected));
    }
    sflint::FmpzMat key(order.degree(), order.degree());
    assert(prime_ideals.at(1).get_hnf(sflint::FmpzMatRef(key)));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(key),
                                  sflint::FmpzMatConstRef(forms.at(0).key)));

    // A form is used only for the ideal it was built for.
    assert(forms.at(1).has_form);
    sflint::FmpzMat product(order.degree(), order.degree());
    assert(!silex::detail::try_ideal_product_hnf_with_form(
            sflint::FmpzMatRef(product), inputs.at(1), prime_ideals.at(2),
            forms.at(1).form));
    return 0;
}

}  // namespace

int main() {
    assert(test_direct_factor_base_witness() == 0);
    assert(test_reduced_selected_prime_outside_factor_base() == 0);
    assert(test_continuation_witness() == 0);
    assert(test_exhaustion_is_named_and_transactional() == 0);
    assert(test_random_product_row_is_exact() == 0);
    assert(test_factor_base_prime_form_product_matches_ideal_multiply() == 0);
    return 0;
}
