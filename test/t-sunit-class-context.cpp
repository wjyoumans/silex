#include <silex/class_group.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>

#include "sunit/sunit_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct ProvenQuadraticFixture {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext class_group;
};

ProvenQuadraticFixture proven_quadratic(slong radicand) noexcept {
    ProvenQuadraticFixture out;
    out.field = silex::test::quadratic_field(radicand);
    silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    assert(out.order.maximal_order(equation));
    assert(out.order.is_maximal());

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 4096;
    options.max_relations = 256;
    options.requested_certification = silex::CertificationMode::proven;
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), out.order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            out.class_group, out.order, sflint::FmpzConstRef(bound), options,
            192));
    assert(out.class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(out.class_group.has_presentation());
    return out;
}

std::vector<silex::PrimeIdeal> first_prime_above(
        const silex::Order& order,
        slong rational_prime) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList decomposition;
    assert(silex::decompose_prime(
            decomposition, order, sflint::FmpzConstRef(p)));
    assert(decomposition.size() > 0 && decomposition.at(0) != nullptr);
    std::vector<silex::PrimeIdeal> out;
    out.emplace_back(order);
    assert(out.back().set(*decomposition.at(0)));
    return out;
}

bool selected_prime_product(
        silex::FractionalIdeal& out,
        const silex::detail::SUnitClassContext& context,
        sflint::FmpzMatConstRef row) noexcept {
    if (sflint::fmpz_mat_nrows(row) != 1 ||
        sflint::fmpz_mat_ncols(row) !=
                static_cast<slong>(context.selected_primes.size())) {
        return false;
    }
    silex::FractionalIdeal accumulator(context.order);
    silex::FractionalIdeal prime_ideal(context.order);
    silex::FractionalIdeal power(context.order);
    if (!accumulator.one()) {
        return false;
    }
    for (slong i = 0; i < static_cast<slong>(context.selected_primes.size());
         ++i) {
        sflint::FmpzConstRef exponent = sflint::fmpz_mat_entry(row, 0, i);
        if (sflint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!silex::detail::prime_to_fractional_ideal(
                    prime_ideal,
                    context.selected_primes[static_cast<std::size_t>(i)]) ||
            !power.pow_fmpz(prime_ideal, exponent) ||
            !accumulator.multiply(accumulator, power)) {
            return false;
        }
    }
    out.swap(accumulator);
    return true;
}

// Index k with selected prime j equal to factor-base prime Q_k, or -1.
// Test oracle: this deliberately duplicates
// silex::detail::selected_factor_base_indices rather than calling it, so the
// tests do not take their expected merge from the code under test; do not
// deduplicate.  Both match primes with the exact PrimeIdeal::equal, so this
// copy does not check the matching primitive itself; the valuations it is used
// with come independently from PrimeIdeal::valuation.
std::vector<slong> selected_factor_base_indices(
        const silex::detail::SUnitClassContext& context) noexcept {
    std::vector<slong> out(context.selected_primes.size(), -1);
    silex::PrimeIdeal prime(context.order);
    for (slong k = 0; k < context.factor_base.length(); ++k) {
        assert(context.factor_base.prime(prime, k));
        for (std::size_t j = 0; j < context.selected_primes.size(); ++j) {
            if (context.selected_primes[j].equal(prime)) {
                out[j] = k;
            }
        }
    }
    return out;
}

// Asserts that the factored valuations of `element` at every factor-base
// prime Q_k and every selected prime P_j equal the vectors derived in the
// ideal group: fb_part[k] + (s_part[j] when Q_k == P_j) at Q_k, and
// s_part[j] + (fb_part[k] when P_j == Q_k) at P_j.  Either row may be empty.
void assert_factored_valuations(
        const silex::detail::SUnitClassContext& context,
        const silex::FactoredElement& element,
        const std::vector<slong>& selected_fb_index,
        sflint::FmpzMatConstRef fb_part,
        sflint::FmpzMatConstRef s_part) noexcept {
    const slong fb_count = context.factor_base.length();
    const slong selected_count =
            static_cast<slong>(context.selected_primes.size());
    std::vector<sflint::Fmpz> expected_fb(static_cast<std::size_t>(fb_count));
    std::vector<sflint::Fmpz> expected_s(
            static_cast<std::size_t>(selected_count));
    for (slong k = 0; k < fb_count; ++k) {
        sflint::fmpz_set(sflint::FmpzRef(expected_fb[static_cast<std::size_t>(k)]),
                         sflint::fmpz_mat_entry(fb_part, 0, k));
    }
    for (slong j = 0; j < selected_count; ++j) {
        sflint::fmpz_set(sflint::FmpzRef(expected_s[static_cast<std::size_t>(j)]),
                         sflint::fmpz_mat_entry(s_part, 0, j));
    }
    for (slong j = 0; j < selected_count; ++j) {
        const slong k = selected_fb_index[static_cast<std::size_t>(j)];
        if (k < 0) {
            continue;
        }
        sflint::Fmpz merged;
        sflint::fmpz_add(sflint::FmpzRef(merged),
                         sflint::FmpzConstRef(
                                 expected_fb[static_cast<std::size_t>(k)]),
                         sflint::FmpzConstRef(
                                 expected_s[static_cast<std::size_t>(j)]));
        sflint::fmpz_set(sflint::FmpzRef(expected_fb[static_cast<std::size_t>(k)]),
                         sflint::FmpzConstRef(merged));
        sflint::fmpz_set(sflint::FmpzRef(expected_s[static_cast<std::size_t>(j)]),
                         sflint::FmpzConstRef(merged));
    }
    silex::PrimeIdeal prime(context.order);
    for (slong k = 0; k < fb_count; ++k) {
        slong valuation = 0;
        assert(context.factor_base.prime(prime, k));
        assert(prime.valuation(valuation, element));
        assert(sflint::fmpz_equal_si(
                expected_fb[static_cast<std::size_t>(k)], valuation));
    }
    for (slong j = 0; j < selected_count; ++j) {
        slong valuation = 0;
        assert(context.selected_primes[static_cast<std::size_t>(j)].valuation(
                valuation, element));
        assert(sflint::fmpz_equal_si(
                expected_s[static_cast<std::size_t>(j)], valuation));
    }
}

// Differential oracle for the class HNF witnesses: expansion (the former
// check) and factored valuations at every factor-base prime both match the
// witness's HNF row.
void assert_class_hnf_witnesses(
        const silex::detail::SUnitClassContext& context) noexcept {
    const slong fb_count = context.factor_base.length();
    sflint::FmpzMat row(1, fb_count);
    const std::vector<slong> no_selected(context.selected_primes.size(), -1);
    sflint::FmpzMat zero_s(1, static_cast<slong>(context.selected_primes.size()));
    for (std::size_t i = 0; i < context.class_hnf.witnesses.size(); ++i) {
        for (slong k = 0; k < fb_count; ++k) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, k),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(context.class_hnf.rows),
                            static_cast<slong>(i), k));
        }
        silex::FractionalIdeal expected(context.order);
        silex::FractionalIdeal principal(context.order);
        silex::Element value(*context.order.parent());
        assert(silex::detail::factor_base_row_ideal(
                expected, context.factor_base, sflint::FmpzMatConstRef(row)));
        assert(context.class_hnf.witnesses[i].evaluate(value));
        assert(principal.set_principal(value));
        assert(principal.equal(expected));
        assert_factored_valuations(context, context.class_hnf.witnesses[i],
                                   selected_factor_base_indices(context),
                                   sflint::FmpzMatConstRef(row),
                                   sflint::FmpzMatConstRef(zero_s));
    }
}

void assert_exact_context_identities(
        const silex::detail::SUnitClassContext& context) noexcept {
    const slong selected_count =
            static_cast<slong>(context.selected_primes.size());
    const slong generator_count = context.s_class_group.generator_count();
    const slong relation_count = context.s_class_group.relation_count();
    assert(context.defined);
    assert(context.source_class_certification ==
           silex::CertificationMode::proven);
    assert(context.s_class_proof_status == silex::ProofState::verified);
    assert(context.s_unit_mod_units_proof_status ==
           silex::ProofState::verified);
    assert(sflint::fmpz_mat_nrows(context.relation_kernel) == selected_count);
    assert(sflint::fmpz_mat_ncols(context.relation_kernel) == relation_count);
    assert(sflint::fmpz_mat_nrows(context.generator_coefficients) ==
           selected_count);
    assert(sflint::fmpz_mat_ncols(context.generator_coefficients) ==
           relation_count);
    assert(sflint::fmpz_mat_nrows(context.valuation_rows) == selected_count);
    assert(sflint::fmpz_mat_ncols(context.valuation_rows) == selected_count);
    assert(context.generators_mod_units.size() ==
           static_cast<std::size_t>(selected_count));

    sflint::FmpzMat zero(selected_count, generator_count);
    sflint::fmpz_mat_mul(
            sflint::FmpzMatRef(zero),
            sflint::FmpzMatConstRef(context.generator_coefficients),
            sflint::FmpzMatConstRef(context.augmented_relations));
    assert(::fmpz_mat_is_zero(zero.raw()) != 0);
    assert(selected_count == 0 ||
           sflint::fmpz_mat_rank(
                   sflint::FmpzMatConstRef(context.valuation_rows)) ==
                   selected_count);

    assert_class_hnf_witnesses(context);
    const std::vector<slong> selected_fb_index =
            selected_factor_base_indices(context);
    sflint::FmpzMat zero_fb(1, context.factor_base.length());
    sflint::FmpzMat row(1, selected_count);
    for (slong i = 0; i < selected_count; ++i) {
        for (slong j = 0; j < selected_count; ++j) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, j),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(context.valuation_rows),
                            i, j));
            slong valuation = 0;
            assert(context.selected_primes[static_cast<std::size_t>(j)].
                           valuation(
                                   valuation,
                                   context.generators_mod_units[
                                           static_cast<std::size_t>(i)]));
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(row), 0, j),
                    valuation));
        }
        silex::FractionalIdeal expected(context.order);
        silex::FractionalIdeal principal(context.order);
        silex::Element value(*context.order.parent());
        assert(selected_prime_product(expected, context,
                                      sflint::FmpzMatConstRef(row)));
        assert(context.generators_mod_units[static_cast<std::size_t>(i)].
                       evaluate(value));
        assert(principal.set_principal(value));
        assert(principal.equal(expected));
        assert_factored_valuations(
                context,
                context.generators_mod_units[static_cast<std::size_t>(i)],
                selected_fb_index, sflint::FmpzMatConstRef(zero_fb),
                sflint::FmpzMatConstRef(row));
    }

    const slong invariant_count = context.s_class_group.invariant_count();
    assert(context.s_class_invariant_ideals.size() ==
           static_cast<std::size_t>(invariant_count));
    assert(context.s_class_power_witnesses.size() ==
           static_cast<std::size_t>(invariant_count));
    assert(sflint::fmpz_mat_nrows(
                   context.s_class_power_selected_exponents) ==
           invariant_count);
    assert(sflint::fmpz_mat_ncols(
                   context.s_class_power_selected_exponents) ==
           selected_count);
    for (slong i = 0; i < invariant_count; ++i) {
        sflint::Fmpz invariant;
        assert(context.s_class_group.invariant(sflint::FmpzRef(invariant), i));
        silex::FractionalIdeal ideal_power(context.order);
        silex::FractionalIdeal selected_product(context.order);
        silex::FractionalIdeal expected(context.order);
        silex::FractionalIdeal principal(context.order);
        silex::Element witness_value(*context.order.parent());
        for (slong j = 0; j < selected_count; ++j) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, j),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(
                                    context.s_class_power_selected_exponents),
                            i, j));
        }
        assert(ideal_power.pow_fmpz(
                context.s_class_invariant_ideals[
                        static_cast<std::size_t>(i)],
                sflint::FmpzConstRef(invariant)));
        assert(selected_prime_product(selected_product, context,
                                      sflint::FmpzMatConstRef(row)));
        assert(expected.multiply(ideal_power, selected_product));
        assert(context.s_class_power_witnesses[static_cast<std::size_t>(i)].
                       evaluate(witness_value));
        assert(principal.set_principal(witness_value));
        assert(principal.equal(expected));

        // Derived valuation vector of the witness: d_i times the invariant
        // generator row at the factor-base primes, and the selected
        // exponents at S, merged where S meets the factor base.
        sflint::FmpzMat generator_row(1, context.factor_base.length());
        sflint::FmpzMat invariant_rows(invariant_count,
                                       context.factor_base.length());
        assert(context.s_class_group.invariant_generator_matrix(
                sflint::FmpzMatRef(invariant_rows)));
        for (slong k = 0; k < context.factor_base.length(); ++k) {
            sflint::fmpz_mul(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(generator_row),
                                           0, k),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(invariant_rows), i, k),
                    sflint::FmpzConstRef(invariant));
        }
        assert_factored_valuations(
                context,
                context.s_class_power_witnesses[static_cast<std::size_t>(i)],
                selected_fb_index, sflint::FmpzMatConstRef(generator_row),
                sflint::FmpzMatConstRef(row));
    }
}

slong selected_primes_in_factor_base(
        const silex::detail::SUnitClassContext& context) noexcept {
    slong count = 0;
    for (slong k : selected_factor_base_indices(context)) {
        count += k >= 0 ? 1 : 0;
    }
    return count;
}

void assert_group_order(const silex::FiniteAbelianGroup& group,
                        slong expected) noexcept {
    sflint::Fmpz order;
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, expected));
}

int test_empty_s_reproduces_class_group() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    const slong relations_before = fixture.class_group.relation_count();
    std::vector<silex::PrimeIdeal> selected;
    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert(result.success);
    assert_group_order(context.s_class_group, 2);
    assert(context.s_class_group.invariant_count() ==
           fixture.class_group.invariant_count());
    for (slong i = 0; i < fixture.class_group.invariant_count(); ++i) {
        sflint::Fmpz expected;
        sflint::Fmpz actual;
        assert(fixture.class_group.invariant(sflint::FmpzRef(expected), i));
        assert(context.s_class_group.invariant(sflint::FmpzRef(actual), i));
        assert(sflint::fmpz_equal(sflint::FmpzConstRef(actual),
                                  sflint::FmpzConstRef(expected)));
    }
    assert(context.selected_primes.empty());
    assert_exact_context_identities(context);
    assert(fixture.class_group.relation_count() == relations_before);
    assert(fixture.class_group.certification_status() ==
           silex::CertificationMode::proven);

    silex::ClassGroupContext unknown;
    silex::ClassGroupCandidateOptions unknown_options;
    unknown_options.max_candidates = 256;
    unknown_options.max_relations = 32;
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), fixture.order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    assert(unknown.compute_candidate(fixture.order,
                                     sflint::FmpzConstRef(bound),
                                     unknown_options));
    assert(unknown.has_presentation());
    assert(unknown.certification_status() ==
           silex::CertificationMode::unknown);
    assert(!silex::detail::build_sunit_class_context(
            result, context, unknown, selected));
    assert(result.stage ==
           silex::detail::SUnitClassBuildStage::input_validation);
    assert(context.defined);
    assert_group_order(context.s_class_group, 2);
    return 0;
}

int test_trivial_class_selected_prime() {
    ProvenQuadraticFixture fixture = proven_quadratic(2);
    assert(fixture.class_group.invariant_count() == 0);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert_group_order(context.s_class_group, 1);
    assert(sflint::fmpz_mat_nrows(context.valuation_rows) == 1);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(context.valuation_rows), 0, 0),
            1));
    assert_exact_context_identities(context);
    return 0;
}

int test_nonprincipal_selected_prime_kills_class_group() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::Ideal selected_ideal(fixture.order);
    silex::FractionalIdeal selected_fractional(fixture.order);
    assert(selected[0].get_ideal(selected_ideal));
    assert(selected_fractional.set_integral(selected_ideal));
    sflint::FmpzMat class_coordinates(
            1, fixture.class_group.invariant_count());
    assert(fixture.class_group.ideal_class_coordinates(
            sflint::FmpzMatRef(class_coordinates), selected_fractional));
    assert(!sflint::fmpz_is_zero(sflint::fmpz_mat_entry(
            sflint::FmpzMatConstRef(class_coordinates), 0, 0)));

    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert_group_order(context.s_class_group, 1);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(context.valuation_rows), 0, 0),
            2));
    assert_exact_context_identities(context);
    assert(selected_primes_in_factor_base(context) == 1);

    std::vector<silex::PrimeIdeal> duplicate;
    duplicate.emplace_back(fixture.order);
    duplicate.emplace_back(fixture.order);
    assert(duplicate[0].set(selected[0]));
    assert(duplicate[1].set(selected[0]));
    assert(!silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, duplicate));
    assert(result.stage ==
           silex::detail::SUnitClassBuildStage::input_validation);
    assert(context.defined);
    assert_group_order(context.s_class_group, 1);
    return 0;
}

int test_nontrivial_s_class_group() {
    ProvenQuadraticFixture fixture = proven_quadratic(-14);
    sflint::Fmpz ordinary_order;
    assert(fixture.class_group.order(sflint::FmpzRef(ordinary_order)));
    assert(sflint::fmpz_equal_si(ordinary_order, 4));
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert_group_order(context.s_class_group, 2);
    assert(context.s_class_group.invariant_count() == 1);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(context.valuation_rows), 0, 0),
            2));
    assert_exact_context_identities(context);
    return 0;
}

// S meeting the factor base in one prime and avoiding it in another: the
// derived valuation vectors must merge P_j with its equal factor-base prime.
int test_selected_primes_meeting_and_avoiding_factor_base() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    std::vector<silex::PrimeIdeal> large =
            first_prime_above(fixture.order, 1009);
    selected.emplace_back(fixture.order);
    assert(selected.back().set(large[0]));
    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert(selected_primes_in_factor_base(context) == 1);
    assert_group_order(context.s_class_group, 1);
    assert_exact_context_identities(context);
    return 0;
}

// Class number three, empty S and S meeting the factor base.
int test_class_number_three() {
    ProvenQuadraticFixture fixture = proven_quadratic(-23);
    std::vector<silex::PrimeIdeal> selected;
    silex::detail::SUnitClassContext context;
    silex::detail::SUnitClassBuildResult result;
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert_group_order(context.s_class_group, 3);
    assert_exact_context_identities(context);

    selected = first_prime_above(fixture.order, 2);
    assert(silex::detail::build_sunit_class_context(
            result, context, fixture.class_group, selected));
    assert_group_order(context.s_class_group, 1);
    assert_exact_context_identities(context);
    return 0;
}

// Composition and the factored valuation check fail closed on slong
// overflow instead of wrapping: FactoredElement::pow_si rejects an exponent
// product that does not fit, and PrimeIdeal::valuation rejects a factored
// valuation that does not fit, so the build stage that calls them fails.
int test_composition_overflow_fails_closed() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    std::vector<silex::PrimeIdeal> above_two =
            first_prime_above(fixture.order, 2);
    silex::Element two(fixture.field);
    assert(two.set_si(2));

    silex::FactoredElement base(fixture.field);
    silex::FactoredElement accumulator(fixture.field);
    assert(base.push(two, WORD_MAX / 2));
    assert(accumulator.one());
    sflint::Fmpz exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), 3);
    assert(!silex::detail::multiply_factored_element_power_fmpz(
            accumulator, base, sflint::FmpzConstRef(exponent)));
    assert(accumulator.length() == 0);

    silex::FactoredElement huge(fixture.field);
    assert(huge.push(two, WORD_MAX));
    slong valuation = 0;
    assert(!above_two[0].valuation(valuation, huge));
    return 0;
}

// When a selected prime P equals a factor-base prime Q_k, the derived
// valuation of an S-class witness at P adds d * generator_row_k to the
// selected exponent instead of treating FB and S as independent coordinates.
int test_invariant_witness_valuation_merges_factor_base_prime() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    silex::detail::SUnitClassContext context;
    context.order = fixture.order;
    assert(context.factor_base.set(*fixture.class_group.factor_base()));
    context.selected_primes = first_prime_above(fixture.order, 2);
    std::vector<slong> index;
    assert(silex::detail::selected_factor_base_indices(index, context));
    assert(index.size() == 1 && index[0] >= 0);

    // (2) = P^2 for the ramified prime P above 2.
    silex::Element two(fixture.field);
    silex::FactoredElement witness(fixture.field);
    assert(two.set_si(2));
    assert(witness.set_element(two));
    sflint::FmpzMat generator_row(1, context.factor_base.length());
    sflint::FmpzMat selected(1, 1);
    sflint::Fmpz invariant;
    sflint::fmpz_set_si(sflint::FmpzRef(invariant), 2);

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(generator_row), 0, index[0]),
                        1);
    assert(silex::detail::verify_s_class_invariant_witness(
            context, index, sflint::FmpzMatConstRef(generator_row),
            sflint::FmpzConstRef(invariant), witness,
            sflint::FmpzMatConstRef(selected)));

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(selected), 0, 0),
                        2);
    assert(!silex::detail::verify_s_class_invariant_witness(
            context, index, sflint::FmpzMatConstRef(generator_row),
            sflint::FmpzConstRef(invariant), witness,
            sflint::FmpzMatConstRef(selected)));

    sflint::fmpz_zero(sflint::fmpz_mat_entry(
            sflint::FmpzMatRef(generator_row), 0, index[0]));
    assert(silex::detail::verify_s_class_invariant_witness(
            context, index, sflint::FmpzMatConstRef(generator_row),
            sflint::FmpzConstRef(invariant), witness,
            sflint::FmpzMatConstRef(selected)));
    return 0;
}
}  // namespace

int main() {
    assert(test_empty_s_reproduces_class_group() == 0);
    assert(test_trivial_class_selected_prime() == 0);
    assert(test_nonprincipal_selected_prime_kills_class_group() == 0);
    assert(test_nontrivial_s_class_group() == 0);
    assert(test_selected_primes_meeting_and_avoiding_factor_base() == 0);
    assert(test_class_number_three() == 0);
    assert(test_composition_overflow_fails_closed() == 0);
    assert(test_invariant_witness_valuation_merges_factor_base_prime() == 0);
    return 0;
}
