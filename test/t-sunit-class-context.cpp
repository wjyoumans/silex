#include <silex/class_group.hpp>
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
    }
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

}  // namespace

int main() {
    assert(test_empty_s_reproduces_class_group() == 0);
    assert(test_trivial_class_selected_prime() == 0);
    assert(test_nonprincipal_selected_prime_kills_class_group() == 0);
    assert(test_nontrivial_s_class_group() == 0);
    return 0;
}
