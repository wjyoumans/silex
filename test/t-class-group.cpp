#include <silex/class_group.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>
#include <silex/residue_field.hpp>
#include <silex/zeta.hpp>

#include "class_group/class_group_certification_internal.hpp"
#include "class_group/class_group_internal.hpp"
#include "class_group/relation_factor_base_plan_internal.hpp"
#include "class_group/relation_saturation_internal.hpp"
#include "class_group/relation_search_internal.hpp"
#include "factor_base/factor_base_internal.hpp"
#include "ideal/ideal_access_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"
#include "zeta/zeta_internal.hpp"

#include <flint/nmod_mat.h>

#include <cassert>
#include <cstring>
#include <initializer_list>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_nothrow_default_constructible_v<
              silex::ClassGroupContext>);
static_assert(std::is_nothrow_move_constructible_v<
              silex::ClassGroupContext>);
static_assert(std::is_nothrow_move_assignable_v<silex::ClassGroupContext>);
static_assert(std::is_nothrow_destructible_v<silex::ClassGroupContext>);

template <typename Options>
concept CandidateOptionsAccepted = requires(
        silex::ClassGroupContext& context,
        const silex::Order& order,
        silex::flint::FmpzConstRef bound,
        const Options& options) {
    context.compute_candidate(order, bound, options);
};

template <typename Options>
concept HasRelationSearchPolicyKnob = requires(Options& options) {
    options.coordinate_search_radius;
} || requires(Options& options) {
    options.ideal_search_radius;
} || requires(Options& options) {
    options.target_relation_kernel_units;
} || requires(Options& options) {
    options.post_finite_refinement_phase_budget;
};

template <typename Options>
concept HasRequestedCertification = requires(Options& options) {
    options.requested_certification;
};

static_assert(CandidateOptionsAccepted<silex::ClassGroupCandidateOptions>);
static_assert(!CandidateOptionsAccepted<silex::ClassGroupComputeOptions>);
static_assert(!HasRelationSearchPolicyKnob<
              silex::ClassGroupCandidateOptions>);
static_assert(!HasRelationSearchPolicyKnob<silex::ClassGroupComputeOptions>);
static_assert(!HasRequestedCertification<silex::ClassGroupCandidateOptions>);
static_assert(HasRequestedCertification<silex::ClassGroupComputeOptions>);

namespace {
namespace sflint = silex::flint;
using CertificationAccess = silex::detail::ClassGroupCertificationAccess;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField shifted_quadratic_field(slong c1, slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField cubic_field(slong c1, slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quintic_field(slong c3,
                                 slong c2,
                                 slong c1,
                                 slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, c2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

bool class_group_bound_at_least_two(sflint::Fmpz& out,
                                    const silex::Order& order) noexcept {
    if (!silex::factor_base_class_group_bound(sflint::FmpzRef(out), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(out), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(out), 2);
    }
    return true;
}

silex::ClassGroupContext degree_one_class_group_context() noexcept {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    return context;
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool context_relation_rows_match_exact_generators(
        const silex::ClassGroupContext& context) noexcept {
    const silex::FactorBase* base = context.factor_base();
    const silex::Order* order = context.parent();
    if (base == nullptr || order == nullptr || order->parent() == nullptr ||
        context.generator_count() != base->length()) {
        return false;
    }

    sflint::FmpzMat rows(context.relation_count(), base->length());
    if (!context.relations(sflint::FmpzMatRef(rows))) {
        return false;
    }

    silex::Element generator(*order->parent());
    silex::Relation exact(*base);
    sflint::FmpzMat exact_row(1, base->length());
    for (slong i = 0; i < context.relation_count(); ++i) {
        if (!context.relation_generator(generator, i) ||
            !exact.set_generator(generator) ||
            !exact.exponents(sflint::FmpzMatRef(exact_row))) {
            return false;
        }
        for (slong j = 0; j < base->length(); ++j) {
            if (!sflint::fmpz_equal(
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(rows), i, j),
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(exact_row), 0, j))) {
                return false;
            }
        }
    }
    return true;
}

void assert_unknown_class_group_certification(
        const silex::ClassGroupContext& context) noexcept {
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(context.analytic_class_regulator_status() ==
           silex::ProofState::not_checked);
    assert(context.unit_proof_status() == silex::ProofState::not_checked);
    assert(context.regulator_proof_status() ==
           silex::ProofState::not_checked);
}

void assert_proven_class_group_certification(
        const silex::ClassGroupContext& context) noexcept {
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.unit_proof_status() == silex::ProofState::verified);
    assert(context.regulator_proof_status() ==
           silex::ProofState::verified);
}

int test_certification_metadata_defaults_and_invalidation() {
    silex::ClassGroupContext empty;
    assert_unknown_class_group_certification(empty);
    assert(empty.relation_saturation_record_count() == 0);

    sflint::Fmpz ell;
    silex::ProofState status = silex::ProofState::verified;
    assert(set_fmpz_si(ell, 77));
    assert(!empty.relation_saturation_record(
            sflint::FmpzRef(ell), status, 0));
    assert(sflint::fmpz_equal_si(ell, 77));
    assert(status == silex::ProofState::verified);

    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    assert_unknown_class_group_certification(context);
    assert(context.relation_saturation_record_count() == 0);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert_unknown_class_group_certification(context);
    assert(context.relation_saturation_record_count() == 0);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(!context.has_presentation());
    assert_unknown_class_group_certification(context);

    assert(context.publish_presentation());
    assert(context.has_presentation());
    assert_unknown_class_group_certification(context);
    assert(context.relation_saturation_record_count() == 0);
    assert(set_fmpz_si(ell, 77));
    status = silex::ProofState::verified;
    assert(!context.relation_saturation_record(
            sflint::FmpzRef(ell), status, 0));
    assert(sflint::fmpz_equal_si(ell, 77));
    assert(status == silex::ProofState::verified);

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(!context.has_presentation());
    assert_unknown_class_group_certification(context);
    assert(context.relation_saturation_record_count() == 0);

    return 0;
}

int test_factor_base_generation_metadata() {
    silex::NumberField field = quadratic_field(-5);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 128;
    options.max_relations = 16;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);

    sflint::Fmpz value;
    assert(context.factor_base_build_bound(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));
    assert(context.factor_base_generation_bound(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));
    assert(!context.factor_base_generation_checked_bound(
            sflint::FmpzRef(value)));
    assert(context.factor_base_generation_record_count() == 0);

    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_bound(
            sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));
    assert(context.factor_base_generation_record_count() == 1);

    silex::ProofState status = silex::ProofState::not_checked;
    assert(context.factor_base_generation_record(
            sflint::FmpzRef(value), status, 0));
    assert(sflint::fmpz_equal_si(value, 2));
    assert(status == silex::ProofState::verified);
    auto generation_record = context.factor_base_generation_record(0);
    assert(generation_record.has_value());
    assert(sflint::fmpz_equal_si(generation_record->p, 2));
    assert(generation_record->status == silex::ProofState::verified);
    assert(!context.factor_base_generation_record(
            sflint::FmpzRef(value), status, 1));
    assert(!context.factor_base_generation_record(1).has_value());

    assert(set_fmpz_si(bound, 3));
    assert(!context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::unavailable);
    assert(context.factor_base_generation_checked_bound(
            sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 3));
    assert(context.factor_base_generation_record_count() == 0);

    return 0;
}

int test_factor_base_honesty_receipt_requires_canonical_bound() {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;

    silex::NumberField field = quadratic_field(-5);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz build_bound;
    assert(set_fmpz_si(build_bound, 2));
    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 128;
    options.max_relations = 16;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(
            maximal_order, sflint::FmpzConstRef(build_bound), options));
    assert(context.has_presentation());
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);

    sflint::Fmpz coverage_bound;
    sflint::Fmpz checked_bound;
    assert(set_fmpz_si(coverage_bound, 1));
    assert(!CertificationAccess::record_factor_base_honesty_proof(
            context, sflint::FmpzConstRef(coverage_bound)));
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);
    assert(!context.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));

    assert(set_fmpz_si(coverage_bound, 3));
    assert(CertificationAccess::record_factor_base_honesty_proof(
            context, sflint::FmpzConstRef(coverage_bound)));
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));
    assert(sflint::fmpz_equal_si(checked_bound, 2));
    return 0;
}

int test_append_and_publish() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.has_factor_base());
    assert(context.generator_count() == 1);
    assert(!context.has_presentation());

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);

    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) == 0);

    silex::ClassGroupRelationSource source =
            silex::ClassGroupRelationSource::Unknown;
    assert(context.relation_source(source, 0));
    assert(source == silex::ClassGroupRelationSource::Supplied);
    assert(!context.relation_source(source, 1));
    assert(source == silex::ClassGroupRelationSource::Supplied);

    sflint::FmpzMat rows(1, 1);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, 2));
    auto owned_rows = context.relations();
    assert(owned_rows.has_value());
    assert(sflint::fmpz_mat_equal(*owned_rows, rows));
    assert(!context.has_presentation());

    assert(context.publish_presentation());
    assert(context.has_presentation());
    assert(context.invariant_count() == 1);

    sflint::Fmpz value;
    assert(context.invariant(sflint::FmpzRef(value), 0));
    assert(sflint::fmpz_equal_si(value, 2));
    auto owned_invariant = context.invariant(0);
    assert(owned_invariant.has_value());
    assert(sflint::fmpz_equal_si(*owned_invariant, 2));
    assert(context.order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));
    auto owned_order = context.order();
    assert(owned_order.has_value());
    assert(sflint::fmpz_equal_si(*owned_order, 2));
    assert(context.relation_kernel_unit_count() == 0);

    silex::FiniteAbelianGroup presentation;
    assert(context.presentation(presentation));
    auto owned_presentation = context.presentation();
    assert(owned_presentation.has_value());
    assert(presentation.order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));
    assert(owned_presentation->order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));

    silex::Element generator(field);
    assert(context.relation_generator(generator, 0));
    assert(generator.equal_si(4));

    silex::PrimeIdeal prime(order);
    assert(context.factor_base_prime(prime, 0));
    assert(prime.rational_prime(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 2));

    sflint::FmpzMat invariant_rows(1, 1);
    assert(context.invariant_generator_matrix(
            sflint::FmpzMatRef(invariant_rows)));
    auto owned_invariant_rows = context.invariant_generator_matrix();
    assert(owned_invariant_rows.has_value());
    assert(sflint::fmpz_mat_equal(*owned_invariant_rows, invariant_rows));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(invariant_rows, 0, 0), 1) ||
           sflint::fmpz_equal_si(
                   sflint::fmpz_mat_entry(invariant_rows, 0, 0), -1));

    silex::FractionalIdeal invariant_generator(order);
    sflint::FmpzMat factor_row(1, 1);
    sflint::FmpzMat coords(1, 1);
    assert(context.invariant_generator(invariant_generator, 0));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(factor_row),
                                         invariant_generator, *base));
    assert(sflint::fmpz_equal(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(factor_row), 0, 0),
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(invariant_rows), 0, 0)));
    assert(context.ideal_class_coordinates(sflint::FmpzMatRef(coords),
                                           invariant_generator));
    assert(mat_entry_is_si(coords, 0, 0, 1));

    silex::FractionalIdeal principal_relation(order);
    assert(principal_relation.set_principal(alpha));
    assert(context.ideal_class_coordinates(sflint::FmpzMatRef(coords),
                                           principal_relation));
    assert(mat_entry_is_si(coords, 0, 0, 0));

    sflint::FmpzMat bad_coords(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(bad_coords, 0, 0), 44);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(bad_coords, 0, 1), 45);
    assert(!context.ideal_class_coordinates(sflint::FmpzMatRef(bad_coords),
                                            invariant_generator));
    assert(mat_entry_is_si(bad_coords, 0, 0, 44));
    assert(mat_entry_is_si(bad_coords, 0, 1, 45));

    assert(alpha.set_si(3));
    silex::FractionalIdeal nonsmooth(order);
    assert(nonsmooth.set_principal(alpha));
    assert(!context.ideal_class_coordinates(sflint::FmpzMatRef(coords),
                                            nonsmooth));
    assert(mat_entry_is_si(coords, 0, 0, 0));

    silex::FractionalIdeal unset;
    assert(!context.invariant_generator(unset, 0));
    assert(!context.invariant_generator(invariant_generator, -1));
    assert(!context.invariant_generator(invariant_generator, 1));

    silex::FactoredElement power_witness(field);
    silex::Element witness_value(field);
    silex::FractionalIdeal witness_ideal(order);
    silex::FractionalIdeal expected_power(order);
    assert(context.invariant_generator_power_witness(power_witness, 0));
    assert(power_witness.evaluate(witness_value));
    assert(witness_ideal.set_principal(witness_value));
    assert(expected_power.pow_fmpz(invariant_generator,
                                  sflint::FmpzConstRef(value)));
    assert(witness_ideal.equal(expected_power));
    assert(!context.invariant_generator_power_witness(power_witness, -1));
    assert(!context.invariant_generator_power_witness(power_witness, 1));
    assert(!context.relation_kernel_unit(power_witness, 0));

    return 0;
}

int test_context_keeps_parent_order_alive() {
    silex::ClassGroupContext context = degree_one_class_group_context();
    assert(context.is_defined());
    assert(context.has_factor_base());
    assert(context.parent() != nullptr);
    assert(context.parent()->parent() != nullptr);
    assert(context.generator_count() == 1);

    silex::PrimeIdeal prime(*context.parent());
    assert(context.factor_base_prime(prime, 0));
    assert(silex::same_order_parent(prime.parent(), context.parent()));

    sflint::Fmpz rational_prime;
    assert(prime.rational_prime(sflint::FmpzRef(rational_prime)));
    assert(sflint::fmpz_equal_si(rational_prime, 2));

    return 0;
}

int test_append_rejects_wrong_factor_base() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    silex::FactorBase copied_base(order);
    assert(copied_base.build(sflint::FmpzConstRef(bound)));
    silex::Relation copied_relation(copied_base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(copied_relation.set_generator(alpha));
    assert(context.append_relation(copied_relation));
    assert(context.relation_count() == 1);

    silex::FactorBase other_base(order);
    assert(set_fmpz_si(bound, 3));
    assert(other_base.build(sflint::FmpzConstRef(bound)));
    silex::Relation wrong_relation(other_base);
    assert(wrong_relation.set_generator(alpha));
    assert(!context.append_relation(wrong_relation));
    assert(context.relation_count() == 1);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    assert(context.has_presentation());

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(!context.has_presentation());
    assert(context.relation_count() == 2);

    return 0;
}

// Under the keep_nonduplicate policy (the LLL relation route's) an exact
// repeat of a stored row is skipped and counted, while a row that differs
// only in sign is not a duplicate and is kept.
int test_keep_nonduplicate_skips_exact_repeat_only() {
    using Access = silex::detail::ClassGroupRelationSearchAccess;
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.generator_count() == 2);
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    bool skipped = true;

    // (1, 0) from -2: kept.
    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(Access::append_keep_nonduplicate_relation(skipped, context,
                                                     relation));
    assert(!skipped);
    assert(context.relation_count() == 1);
    assert(context.skipped_dependent_relation_count() == 0);

    // The same relation again: skipped once.
    assert(Access::append_keep_nonduplicate_relation(skipped, context,
                                                     relation));
    assert(skipped);
    assert(context.relation_count() == 1);
    assert(context.skipped_dependent_relation_count() == 1);

    // (-1, 0) from -1/2 differs only in sign: kept, though dependent.
    assert(alpha.set_si_over_si(-1, 2));
    assert(relation.set_generator(alpha));
    sflint::FmpzMat row(1, 2);
    assert(relation.exponents(sflint::FmpzMatRef(row)));
    assert(mat_entry_is_si(row, 0, 0, -1));
    assert(mat_entry_is_si(row, 0, 1, 0));
    skipped = true;
    assert(Access::append_keep_nonduplicate_relation(skipped, context,
                                                     relation));
    assert(!skipped);
    assert(context.relation_count() == 2);
    assert(context.skipped_dependent_relation_count() == 1);

    // Repeating the sign-flipped row is skipped in turn.
    assert(Access::append_keep_nonduplicate_relation(skipped, context,
                                                     relation));
    assert(skipped);
    assert(context.relation_count() == 2);
    assert(context.skipped_dependent_relation_count() == 2);
    return 0;
}

int test_set_relation_matrix_preserves_on_failure() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::FactorBase base2(order);
    silex::FactorBase base3(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base2.build(sflint::FmpzConstRef(bound)));
    assert(set_fmpz_si(bound, 3));
    assert(base3.build(sflint::FmpzConstRef(bound)));

    silex::Element alpha(field);
    silex::Relation relation2(base2);
    silex::RelationMatrix matrix2(base2);
    assert(alpha.set_si(4));
    assert(relation2.set_generator(alpha));
    assert(matrix2.append(relation2));

    silex::ClassGroupContext context(order);
    sflint::Fmpz order_out;
    assert(context.set_relation_matrix(matrix2));
    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() == 1);
    assert(context.factor_base_generation_status() ==
           silex::ProofState::not_checked);
    assert(!context.factor_base_build_bound(sflint::FmpzRef(order_out)));
    assert(!context.factor_base_generation_bound(sflint::FmpzRef(order_out)));
    assert(!context.factor_base_generation_checked_bound(
            sflint::FmpzRef(order_out)));
    assert(context.factor_base_generation_record_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 1);

    silex::ClassGroupRelationSource source =
            silex::ClassGroupRelationSource::Unknown;
    assert(context.relation_source(source, 0));
    assert(source == silex::ClassGroupRelationSource::Supplied);

    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    silex::Relation relation3(base3);
    silex::RelationMatrix rank_deficient(base3);
    assert(alpha.set_si(2));
    assert(relation3.set_generator(alpha));
    assert(rank_deficient.append(relation3));
    assert(!context.set_relation_matrix(rank_deficient));

    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() == 1);
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    sflint::FmpzMat rows(1, 1);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, 2));

    silex::NumberField other_field = degree_one_field();
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);
    silex::FactorBase other_base(other_order);
    assert(set_fmpz_si(bound, 2));
    assert(other_base.build(sflint::FmpzConstRef(bound)));
    silex::Relation other_relation(other_base);
    silex::RelationMatrix other_matrix(other_base);
    silex::Element other_alpha(other_field);
    assert(other_alpha.set_si(4));
    assert(other_relation.set_generator(other_alpha));
    assert(other_matrix.append(other_relation));
    assert(!context.set_relation_matrix(other_matrix));
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    return 0;
}

int test_relation_kernel_unit_witness() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    assert(context.invariant_count() == 0);
    assert(context.relation_count() == 2);
    assert(context.generator_count() == 1);
    assert(context.relation_kernel_unit_count() == 1);

    silex::FactoredElement unit_witness(field);
    silex::Element unit_value(field);
    silex::FractionalIdeal unit_ideal(order);
    sflint::FmpzMat row(1, 1);
    assert(context.relation_kernel_unit(unit_witness, 0));
    assert(unit_witness.evaluate(unit_value));
    assert(unit_ideal.set_principal(unit_value));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                         unit_ideal, *base));
    assert(mat_entry_is_si(row, 0, 0, 0));
    assert(!context.relation_kernel_unit(unit_witness, -1));
    assert(!context.relation_kernel_unit(unit_witness, 1));
    assert(!context.invariant_generator_power_witness(unit_witness, 0));

    return 0;
}

// relation_kernel_unit reads a reduced kernel cached with the presentation;
// appending a relation and publishing again must rebuild it, so every
// index of the new kernel yields a unit witness.
int test_relation_kernel_unit_follows_republished_presentation() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    // Exponent rows 8, 4, 2, 1: each strictly enlarges the relation
    // lattice, so the native policy keeps every one.
    const slong generators[] = {256, 16, 4, 2};
    silex::FactoredElement unit_witness(field);
    silex::Element unit_value(field);
    silex::FractionalIdeal unit_ideal(order);
    sflint::FmpzMat row(1, 1);
    for (slong count = 1; count <= 4; ++count) {
        assert(alpha.set_si(generators[count - 1]));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        assert(context.relation_count() == count);
        assert(context.relation_kernel_unit_count() == count - 1);
        for (slong i = 0; i < count - 1; ++i) {
            assert(context.relation_kernel_unit(unit_witness, i));
            assert(unit_witness.evaluate(unit_value));
            assert(unit_ideal.set_principal(unit_value));
            assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(row),
                                                 unit_ideal, *base));
            assert(mat_entry_is_si(row, 0, 0, 0));
        }
        assert(!context.relation_kernel_unit(unit_witness, count - 1));
    }
    return 0;
}

// An exponent outside slong is never expanded for a base that is not a
// root of unity: x^e would have height |e| h(x) > 2^63 h(x).  The power
// fails closed and leaves the accumulator unchanged.  A root of unity is
// raised to e modulo the order of the roots of unity instead.
int test_compact_power_fmpz_exponent_guard() {
    silex::NumberField field = silex::test::quadratic_field(-1);
    sflint::Fmpz huge;
    sflint::fmpz_set_ui(sflint::FmpzRef(huge), 1);
    sflint::fmpz_mul_2exp(sflint::FmpzRef(huge), sflint::FmpzConstRef(huge),
                          64);
    sflint::fmpz_add_ui(sflint::FmpzRef(huge), sflint::FmpzConstRef(huge), 3);

    silex::Element value(field);
    silex::Element expected(field);
    silex::FactoredElement base(field);
    silex::FactoredElement accumulator(field);
    assert(value.set_si(2));
    assert(base.set_element(value));
    assert(accumulator.one());
    assert(!silex::detail::multiply_factored_element_power_fmpz(
            accumulator, base, sflint::FmpzConstRef(huge)));
    assert(accumulator.length() == 0);

    // i^(2^64 + 3) = i^3 = -i and i^-(2^64 + 3) = i.
    assert(value.gen());
    assert(base.set_element(value));
    assert(silex::detail::multiply_factored_element_power_fmpz(
            accumulator, base, sflint::FmpzConstRef(huge)));
    assert(accumulator.evaluate(value));
    assert(expected.gen());
    assert(expected.negate(expected));
    assert(value.equal(expected));

    sflint::fmpz_neg(sflint::FmpzRef(huge), sflint::FmpzConstRef(huge));
    assert(accumulator.one());
    assert(silex::detail::multiply_factored_element_power_fmpz(
            accumulator, base, sflint::FmpzConstRef(huge)));
    assert(accumulator.evaluate(value));
    assert(expected.gen());
    assert(value.equal(expected));
    return 0;
}

int test_goal_publication_waits_for_kernel_target() {
    using RelationSearchAccess =
            silex::detail::ClassGroupRelationSearchAccess;
    assert(RelationSearchAccess::relation_kernel_row_count_from_dimensions(
                   0, 1, 0) == 1);
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    silex::detail::ClassUnitTransactionReport audit;
    audit.reset();
    silex::detail::ClassUnitTransactionContext run_context{
            audit,
            nullptr,
            {},
            {},
            false};
    silex::detail::ClassUnitTransactionAccess::set_run_context(
            context, &run_context);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_rank() == context.generator_count());
    assert(context.relation_kernel_unit_count() == 0);
    assert(!silex::detail::relation_search::publish_and_check_compute_goal(
            context, 1, silex::CertificationMode::unknown));
    assert(!context.has_presentation());

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(!context.has_presentation());
    assert(silex::detail::relation_search::publish_and_check_compute_goal(
            context, 1, silex::CertificationMode::unknown));
    assert(context.has_presentation());
    assert(context.relation_kernel_unit_count() == 1);
    assert(context.invariant_count() == 0);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context, nullptr);
    return 0;
}

int test_relation_sources_reset_with_factor_base() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(
            relation, silex::ClassGroupRelationSource::ClassGenerator));
    assert(context.relation_count() == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::ClassGenerator) == 1);

    assert(!context.append_relation(
            relation, static_cast<silex::ClassGroupRelationSource>(99)));
    assert(context.relation_count() == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::ClassGenerator) == 1);

    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.relation_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::ClassGenerator) == 0);

    silex::ClassGroupRelationSource source =
            silex::ClassGroupRelationSource::Unknown;
    assert(!context.relation_source(source, 0));
    assert(source == silex::ClassGroupRelationSource::Unknown);

    return 0;
}

int test_append_tracks_rank_and_skips_nonrefining_dependent_rows() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.generator_count() == 2);
    assert(context.relation_rank() == 0);
    assert(context.skipped_dependent_relation_count() == 0);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);

    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 1);
    assert(context.relation_rank() == 1);
    assert(context.skipped_dependent_relation_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 1);
    assert(!context.has_presentation());

    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 1);
    assert(context.relation_rank() == 1);
    assert(context.skipped_dependent_relation_count() == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 1);

    assert(alpha.set_si(-3));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 2);
    assert(context.relation_rank() == 2);
    assert(context.skipped_dependent_relation_count() == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 2);

    assert(context.publish_presentation());
    assert(context.has_presentation());
    assert(context.invariant_count() == 0);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));

    sflint::FmpzMat rows(2, 2);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, 1));
    assert(mat_entry_is_si(rows, 0, 1, 0));
    assert(mat_entry_is_si(rows, 1, 0, 0));
    assert(mat_entry_is_si(rows, 1, 1, 1));

    silex::Element generator(field);
    assert(context.relation_generator(generator, 0));
    assert(generator.equal_si(-2));
    assert(context.relation_generator(generator, 1));
    assert(generator.equal_si(-3));

    return 0;
}

int test_append_keeps_index_refining_dependent_rows() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);

    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 1);
    assert(context.relation_rank() == 1);
    assert(context.publish_presentation());

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.relation_count() == 2);
    assert(context.relation_rank() == 1);
    assert(context.skipped_dependent_relation_count() == 0);
    assert(!context.has_presentation());
    assert(context.publish_presentation());
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));
    assert(context.relation_kernel_unit_count() == 1);

    return 0;
}

int test_factor_base_prime_search_masks() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.generator_count() == 2);

    bool principal = true;
    bool covered = true;
    assert(!context.factor_base_prime_is_principal(principal, -1));
    assert(principal);
    assert(!context.factor_base_prime_is_hnf_covered(covered, 2));
    assert(covered);

    assert(context.factor_base_prime_is_principal(principal, 0));
    assert(!principal);
    assert(context.factor_base_prime_is_principal(principal, 1));
    assert(!principal);
    assert(context.factor_base_prime_is_hnf_covered(covered, 0));
    assert(!covered);
    assert(context.factor_base_prime_is_hnf_covered(covered, 1));
    assert(!covered);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);

    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.factor_base_prime_is_principal(principal, 0));
    assert(principal);
    assert(context.factor_base_prime_is_principal(principal, 1));
    assert(!principal);
    assert(context.factor_base_prime_is_hnf_covered(covered, 0));
    assert(covered);
    assert(context.factor_base_prime_is_hnf_covered(covered, 1));
    assert(!covered);

    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.skipped_dependent_relation_count() == 1);
    assert(context.factor_base_prime_is_principal(principal, 0));
    assert(principal);
    assert(context.factor_base_prime_is_hnf_covered(covered, 0));
    assert(covered);

    assert(alpha.set_si(-3));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.factor_base_prime_is_principal(principal, 0));
    assert(principal);
    assert(context.factor_base_prime_is_principal(principal, 1));
    assert(principal);
    assert(context.factor_base_prime_is_hnf_covered(covered, 0));
    assert(covered);
    assert(context.factor_base_prime_is_hnf_covered(covered, 1));
    assert(covered);

    return 0;
}

int test_compute_candidate_coordinate_search() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context;
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.compute_candidate(order, sflint::FmpzConstRef(bound)));
    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() == 1);
    assert(context.invariant_count() == 0);
    assert(context.relation_kernel_unit_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) == 1);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Supplied) == 0);

    silex::ClassGroupRelationSource source =
            silex::ClassGroupRelationSource::Unknown;
    assert(context.relation_source(source, 0));
    assert(source == silex::ClassGroupRelationSource::Search);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));

    sflint::FmpzMat rows(1, 1);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, -1) ||
           mat_entry_is_si(rows, 0, 0, 1));

    silex::Element generator(field);
    assert(context.relation_generator(generator, 0));
    assert(generator.equal_si(-2) || generator.equal_si(2));

    return 0;
}

int test_compute_candidate_skips_principal_prime_lattice_search() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 6;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                     options));
    assert(context.has_presentation());
    assert(context.generator_count() == 2);
    assert(context.relation_count() == 2);
    assert(context.invariant_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) == 2);

    sflint::FmpzMat rows(2, 2);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, -1) ||
           mat_entry_is_si(rows, 0, 0, 1));
    assert(mat_entry_is_si(rows, 0, 1, 0));
    assert(mat_entry_is_si(rows, 1, 0, 0));
    assert(mat_entry_is_si(rows, 1, 1, -1) ||
           mat_entry_is_si(rows, 1, 1, 1));

    return 0;
}

int test_compute_candidate_ideal_lattice_search() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                     options));
    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() == 1);
    assert(context.invariant_count() == 0);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) == 1);

    sflint::FmpzMat rows(1, 1);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, -1) ||
           mat_entry_is_si(rows, 0, 0, 1));

    silex::Element generator(field);
    assert(context.relation_generator(generator, 0));
    assert(generator.equal_si(-2) || generator.equal_si(2));

    silex::ClassGroupCandidateOptions limited = options;
    limited.max_candidates = 0;
    silex::ClassGroupContext limited_context;
    assert(!limited_context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                              limited));
    assert(!limited_context.has_presentation());

    return 0;
}

int test_compute_candidate_small_quadratic_policy() {
    silex::NumberField field = quadratic_field(-5);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 128;
    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() >= 2);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) ==
           context.relation_count());

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    sflint::FmpzMat rows(context.relation_count(), 1);
    assert(context.relations(sflint::FmpzMatRef(rows)));
    bool found_refining_relation = false;
    for (slong row = 0; row < context.relation_count(); ++row) {
        found_refining_relation = found_refining_relation ||
                mat_entry_is_si(rows, row, 0, 2) ||
                mat_entry_is_si(rows, row, 0, -2);
    }
    assert(found_refining_relation);

    return 0;
}

int test_compute_candidate_proven_quadratic() {
    silex::NumberField imaginary_field = quadratic_field(-5);
    silex::Order imaginary_equation;
    silex::Order imaginary_maximal(imaginary_field);
    imaginary_equation = silex::test::equation_order(imaginary_field);
    assert(imaginary_maximal.maximal_order(imaginary_equation));
    assert(imaginary_maximal.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), imaginary_maximal));
    assert(sflint::fmpz_equal_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 128;
    options.max_relations = 16;

    silex::ClassGroupContext existing;
    assert(existing.compute_candidate(imaginary_maximal,
                                      sflint::FmpzConstRef(bound), options));
    assert_unknown_class_group_certification(existing);
    assert(!existing.try_certify_quadratic(
            silex::CertificationMode::unknown));
    sflint::Fmpz order_out;
    assert(existing.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));
    assert(existing.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(existing.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(existing.try_certify_quadratic(
            silex::CertificationMode::proven));
    assert_proven_class_group_certification(existing);

    assert(existing.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));
    assert(existing.factor_base_generation_record_count() > 0);

    silex::ClassGroupContext proven_imaginary;
    silex::ClassGroupCandidateOptions proven_options;
    proven_options.max_candidates = 128;
    proven_options.max_relations = 16;
    assert(class_group_bound_at_least_two(bound, imaginary_maximal));
    assert(sflint::fmpz_equal_si(bound, 2));
    assert(proven_options.max_candidates == 128);
    assert(proven_options.max_relations == 16);
    assert(proven_imaginary.compute_candidate(
            imaginary_maximal, sflint::FmpzConstRef(bound), proven_options));
    assert_unknown_class_group_certification(proven_imaginary);
    assert(proven_imaginary.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(proven_imaginary.try_certify_quadratic(
            silex::CertificationMode::proven));
    assert_proven_class_group_certification(proven_imaginary);
    // The exact route (h_cand = h(D) = 2) stores a verified saturation
    // record at ell = 2.
    assert(proven_imaginary.relation_saturation_record_count() == 1);
    auto imaginary_record = proven_imaginary.relation_saturation_record(0);
    assert(imaginary_record.has_value());
    assert(sflint::fmpz_equal_si(imaginary_record->ell, 2));
    assert(imaginary_record->status == silex::ProofState::verified);
    assert(proven_imaginary.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));
    assert(proven_imaginary.relation_count() >= 2);

    bool has_refining_relation = false;
    sflint::FmpzMat rows(proven_imaginary.relation_count(),
                         proven_imaginary.generator_count());
    assert(proven_imaginary.relations(sflint::FmpzMatRef(rows)));
    for (slong i = 0; i < proven_imaginary.relation_count(); ++i) {
        if (mat_entry_is_si(rows, i, 0, 2) ||
            mat_entry_is_si(rows, i, 0, -2)) {
            has_refining_relation = true;
        }
    }
    assert(has_refining_relation);

    silex::NumberField real_field = quadratic_field(2);
    silex::Order real_order;
    real_order = silex::test::equation_order(real_field);
    assert(real_order.is_maximal());
    proven_options = silex::ClassGroupCandidateOptions{};
    proven_options.max_candidates = 256;
    proven_options.max_relations = 32;
    assert(class_group_bound_at_least_two(bound, real_order));
    assert(sflint::fmpz_equal_si(bound, 2));
    assert(proven_options.max_candidates == 256);
    assert(proven_options.max_relations == 32);

    silex::ClassGroupContext proven_real;
    assert(proven_real.compute_candidate(real_order,
                                         sflint::FmpzConstRef(bound),
                                         proven_options));
    assert_unknown_class_group_certification(proven_real);
    assert(proven_real.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(proven_real.try_certify_quadratic(
            silex::CertificationMode::proven));
    // h_cand = 1 with verified generation proves h = 1 without any units,
    // so the class group is proven but no unit or regulator proof is
    // claimed.
    assert(proven_real.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_real.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(proven_real.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(proven_real.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_real.unit_proof_status() == silex::ProofState::not_checked);
    assert(proven_real.regulator_proof_status() ==
           silex::ProofState::not_checked);
    assert(proven_real.relation_saturation_record_count() == 0);
    assert(proven_real.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_is_one(order_out));

    return 0;
}

// The silex::detail ideal shortcuts take caller-supplied generators, but
// they store only the HNF of the ideal those generators determine (see the
// comment above detail::multiply_integral_ideal_by_two_generator in
// src/ideal/ideal.cpp).  A "false" generator, here a pair a caller might take
// for a prime above 2, yields a true ideal that is equal to the public-API
// ideal with the same generators, so a proven class group's discrete
// logarithm sees the same ideal on both routes and the shortcuts leave the
// group's labels unchanged.
int test_detail_ideal_shortcuts_reach_labelled_dlog_as_true_ideals() {
    silex::NumberField field = quadratic_field(-5);
    silex::Order equation;
    silex::Order order(field);
    equation = silex::test::equation_order(field);
    assert(order.maximal_order(equation));

    sflint::Fmpz bound;
    assert(class_group_bound_at_least_two(bound, order));
    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 128;
    options.max_relations = 16;
    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                     options));
    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(context.try_certify_quadratic(silex::CertificationMode::proven));
    assert_proven_class_group_certification(context);
    sflint::Fmpz class_order;
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));

    silex::Element theta(field);
    silex::OrderElement theta_order(order);
    silex::OrderElement one(order);
    silex::OrderElement theta_plus_one(order);
    silex::OrderElement two_element(order);
    assert(theta.gen());
    assert(theta_order.set_element(theta));
    assert(one.set_si(1));
    assert(two_element.set_si(2));
    assert(theta_plus_one.add(theta_order, one));
    sflint::Fmpz two;
    sflint::fmpz_set_ui(two, 2);

    // Public route for a O + beta O: the sum of the two principal ideals.
    auto public_two_generator = [&](silex::Ideal& out,
                                    const silex::OrderElement& beta) {
        silex::Ideal scalar_ideal(order);
        silex::Ideal element_ideal(order);
        return scalar_ideal.set_principal(two_element) &&
               element_ideal.set_principal(beta) &&
               out.add(scalar_ideal, element_ideal);
    };
    auto class_coordinate = [&](const silex::Ideal& ideal, slong expected) {
        silex::FractionalIdeal fractional(order);
        silex::FractionalIdeal shortcut_fractional(order);
        sflint::FmpzMat coords(1, 1);
        sflint::FmpzMat shortcut_coords(1, 1);
        return fractional.set_integral(ideal) &&
               silex::detail::set_integral_ideal_known_hnf(
                       shortcut_fractional, ideal) &&
               shortcut_fractional.equal(fractional) &&
               context.ideal_class_coordinates(sflint::FmpzMatRef(coords),
                                               fractional) &&
               context.ideal_class_coordinates(
                       sflint::FmpzMatRef(shortcut_coords),
                       shortcut_fractional) &&
               mat_entry_is_si(coords, 0, 0, expected) &&
               mat_entry_is_si(shortcut_coords, 0, 0, expected);
    };

    // (2, 1 + theta) is the nonprincipal prime above 2.
    silex::Ideal prime_two(order);
    silex::Ideal public_prime_two(order);
    assert(silex::detail::set_known_two_generator_ideal(
            prime_two, sflint::FmpzConstRef(two), theta_plus_one));
    assert(public_two_generator(public_prime_two, theta_plus_one));
    assert(prime_two.equal(public_prime_two));
    assert(prime_two.known_principal_generator() == nullptr);
    assert(class_coordinate(prime_two, 1));

    // (2, theta) contains 2 and theta^2 = -5, so it is the unit ideal, not a
    // prime above 2; the shortcut stores O and the class is trivial.
    silex::Ideal false_prime(order);
    silex::Ideal public_false_prime(order);
    assert(silex::detail::set_known_two_generator_ideal(
            false_prime, sflint::FmpzConstRef(two), theta_order));
    assert(public_two_generator(public_false_prime, theta_order));
    assert(false_prime.equal(public_false_prime));
    assert(false_prime.is_one());
    assert(class_coordinate(false_prime, 0));

    // P2 (2, 1 + theta) = P2^2 = (2), and P2 (2, theta) = P2.
    silex::Ideal square(order);
    silex::Ideal public_square(order);
    assert(silex::detail::multiply_integral_ideal_by_two_generator(
            square, prime_two, sflint::FmpzConstRef(two), theta_plus_one));
    assert(public_square.multiply(public_prime_two, public_prime_two));
    assert(square.equal(public_square));
    assert(class_coordinate(square, 0));
    silex::Ideal times_unit_ideal(order);
    assert(silex::detail::multiply_integral_ideal_by_two_generator(
            times_unit_ideal, prime_two, sflint::FmpzConstRef(two),
            theta_order));
    assert(times_unit_ideal.equal(public_prime_two));
    assert(class_coordinate(times_unit_ideal, 1));

    // 2 P2 is in the class of P2; (1/2) P2 is not integral, so the element
    // shortcut refuses it and leaves its output unchanged.
    silex::Element two_field(field);
    silex::Element half(field);
    assert(two_field.set_si(2));
    assert(half.set_si_over_si(1, 2));
    silex::Ideal doubled(order);
    silex::Ideal public_doubled(order);
    silex::Ideal principal_two(order);
    assert(silex::detail::multiply_integral_ideal_by_element(
            doubled, prime_two, two_field));
    assert(principal_two.set_principal(two_element));
    assert(public_doubled.multiply(principal_two, public_prime_two));
    assert(doubled.equal(public_doubled));
    assert(class_coordinate(doubled, 1));
    assert(!silex::detail::multiply_integral_ideal_by_element(
            doubled, prime_two, half));
    assert(doubled.equal(public_doubled));

    assert_proven_class_group_certification(context);
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    return 0;
}

// A successful real-quadratic try_certify_quadratic, and
// try_certify_trivial_quotient, prove h = 1, so no prime needs a saturation
// record; both clear stale ell-records and the proof records behind them.
int test_certify_real_quadratic_clears_saturation_records() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(class_group_bound_at_least_two(bound, order));
    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz two;
    assert(set_fmpz_si(two, 2));
    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        assert(context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                         options));
        assert(context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(bound)));
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 80));

        // Saturation at every prime up to 2 leaves an ell = 2 record backed
        // by an ell = 2 proof record.
        sflint::Fmpz index_bound;
        sflint::Fmpz aux_bound;
        assert(set_fmpz_si(index_bound, 2));
        assert(set_fmpz_si(aux_bound, 31));
        assert(CertificationAccess::
                       try_prove_relation_saturation_index_bound_with_units(
                               context, units,
                               sflint::FmpzConstRef(index_bound),
                               sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 1);
        silex::ProofState status = silex::ProofState::not_checked;
        slong rank = -1;
        slong target_rank = -1;
        slong local_primes = -1;
        assert(CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(two), status, rank,
                target_rank, local_primes));
    };
    auto assert_cleared = [&](const silex::ClassGroupContext& context)
            noexcept {
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(context.relation_saturation_record_count() == 0);
        silex::ProofState status = silex::ProofState::not_checked;
        slong rank = -1;
        slong target_rank = -1;
        slong local_primes = -1;
        assert(!CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(two), status, rank,
                target_rank, local_primes));
    };

    silex::ClassGroupContext context;
    prepare_context(context);
    assert(context.try_certify_quadratic(silex::CertificationMode::proven));
    assert_cleared(context);

    silex::ClassGroupContext trivial_context;
    prepare_context(trivial_context);
    assert(trivial_context.try_certify_trivial_quotient(
            silex::CertificationMode::proven));
    assert_cleared(trivial_context);
    return 0;
}

int test_default_maximal_quadratic_proven_policy() {
    auto check_imaginary = [](slong radicand,
                              slong expected_order,
                              slong expected_invariant) noexcept {
        silex::NumberField field = quadratic_field(radicand);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(maximal.is_maximal());

        sflint::Fmpz bound;
        assert(class_group_bound_at_least_two(bound, maximal));
        silex::ClassGroupCandidateOptions options;
        silex::ClassGroupContext context;
        assert(context.compute_candidate(
                maximal, sflint::FmpzConstRef(bound), options));
        assert_unknown_class_group_certification(context);
        assert(context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(bound)));
        assert(context.try_certify_quadratic(
                silex::CertificationMode::proven));
        assert_proven_class_group_certification(context);

        sflint::Fmpz value;
        assert(context.order(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, expected_order));
        assert(context.invariant_count() == 1);
        assert(context.invariant(sflint::FmpzRef(value), 0));
        assert(sflint::fmpz_equal_si(value, expected_invariant));
    };

    check_imaginary(-14, 4, 4);
    check_imaginary(-47, 5, 5);

    return 0;
}

int test_nonmaximal_quadratic_candidate_rejection_preserves_output() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(!equation.is_maximal());
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    sflint::Fmpz bound;
    assert(class_group_bound_at_least_two(bound, maximal));
    silex::ClassGroupCandidateOptions options;

    silex::ClassGroupContext fresh;
    assert(!fresh.compute_candidate(
            equation, sflint::FmpzConstRef(bound), options));
    assert(!fresh.has_factor_base());
    assert(!fresh.has_presentation());

    silex::ClassGroupContext preserved;
    assert(preserved.compute_candidate(
            maximal, sflint::FmpzConstRef(bound), options));
    assert(preserved.check_factor_base_generation_bound(
            sflint::FmpzConstRef(bound)));
    assert(preserved.try_certify_quadratic(
            silex::CertificationMode::proven));
    assert_proven_class_group_certification(preserved);
    const slong relation_count_before = preserved.relation_count();
    const auto relations_before = preserved.relations();
    const auto generators_before = preserved.invariant_generator_matrix();
    const auto order_before = preserved.order();
    assert(relations_before.has_value());
    assert(generators_before.has_value());
    assert(order_before.has_value());

    assert(!preserved.compute_candidate(
            equation, sflint::FmpzConstRef(bound), options));
    assert(silex::same_order_parent(preserved.parent(), &maximal));
    assert(preserved.relation_count() == relation_count_before);
    assert_proven_class_group_certification(preserved);

    const auto relations_after = preserved.relations();
    const auto generators_after = preserved.invariant_generator_matrix();
    const auto order_after = preserved.order();
    assert(relations_after.has_value());
    assert(generators_after.has_value());
    assert(order_after.has_value());
    assert(sflint::fmpz_mat_equal(*relations_before, *relations_after));
    assert(sflint::fmpz_mat_equal(*generators_before, *generators_after));
    assert(sflint::fmpz_equal(*order_before, *order_after));

    return 0;
}

int test_compute_candidate_random_product_relation_cubic() {
    silex::NumberField field = cubic_field(-2, -5);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), maximal_order));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.factor_base() != nullptr);
    assert(context.factor_base()->length() >= 3);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::RandomProduct) > 0);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    return 0;
}

int test_one_large_prime_match_cubic() {
    silex::NumberField field = cubic_field(-6, -10);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));

    silex::ClassGroupContext context(maximal_order);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    sflint::FmpzMat coordinates(1, maximal_order.degree());
    silex::OrderElement order_element(maximal_order);
    silex::Element element(field);
    auto scan_until_large_prime_match = [&]() noexcept {
        const slong relation_matches_before =
                context.relation_source_count(
                        silex::ClassGroupRelationSource::LargePrimeMatch);
        const slong skipped_before =
                context.skipped_dependent_relation_count();
        for (slong a = -10; a <= 10; ++a) {
            for (slong b = -10; b <= 10; ++b) {
                for (slong c = -10; c <= 10; ++c) {
                    bool partial_throttle_exit = false;
                    sflint::fmpz_mat_zero(sflint::FmpzMatRef(coordinates));
                    sflint::fmpz_set_si(
                            sflint::fmpz_mat_entry(coordinates, 0, 0), a);
                    sflint::fmpz_set_si(
                            sflint::fmpz_mat_entry(coordinates, 0, 1), b);
                    sflint::fmpz_set_si(
                            sflint::fmpz_mat_entry(coordinates, 0, 2), c);
                    assert(order_element.set_coordinates(
                            sflint::FmpzMatConstRef(coordinates)));
                    assert(order_element.get_element(element));
                    assert(context.try_append_generator_relation(
                            partial_throttle_exit, element,
                            silex::ClassGroupRelationSource::Search));
                    if (partial_throttle_exit ||
                        context.relation_source_count(
                                silex::ClassGroupRelationSource::
                                        LargePrimeMatch) >
                                relation_matches_before ||
                        context.skipped_dependent_relation_count() >
                                skipped_before) {
                        return partial_throttle_exit;
                    }
                }
            }
        }
        return false;
    };

    assert(!scan_until_large_prime_match());
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::LargePrimeMatch) > 0);
    assert(context_relation_rows_match_exact_generators(context));
    assert(silex::detail::ClassGroupRelationSearchAccess::
                           direct_partial_residue_block_cache_size(context) >
            0);

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    const slong threshold =
            base->length() > 4 ? 2 * base->length() : 8;
    bool throttled = false;
    for (slong i = 0; i < threshold; ++i) {
        throttled = scan_until_large_prime_match();
        if (i + 1 < threshold) {
            assert(!throttled);
        }
    }
    assert(throttled);
    assert(context.skipped_dependent_relation_count() >= threshold);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(silex::detail::ClassGroupRelationSearchAccess::
                           direct_partial_residue_block_cache_size(context) ==
            0);
    return 0;
}

// The exact-valuation branch of the direct large-prime partial-relation
// row reuses the p-adic valuation of the element norm.  Its rows must match
// the per-prime element valuations of the public PrimeIdeal::valuation
// (which recomputes the norm), with the principal-ideal valuation as its
// fallback, whether or not the norm and coordinates are supplied, and
// whether the residue test runs word-size (nmod) or multi-limb (fmpz_mod).
int test_direct_partial_row_exact_valuation_cubic() {
    silex::NumberField field = cubic_field(-6, -10);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 30));
    silex::ClassGroupContext context(maximal_order);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);

    sflint::FmpzMat coordinates(1, maximal_order.degree());
    silex::OrderElement order_element(maximal_order);
    silex::Element element(field);
    sflint::Fmpq norm;
    sflint::FmpzMat row(0, 0);
    sflint::FmpzMat known_row(0, 0);
    sflint::FmpzMat multi_limb_row(0, 0);
    sflint::Fmpz large_prime;
    sflint::Fmpz known_large_prime;
    sflint::Fmpz multi_limb_large_prime;
    sflint::Fmpz leftover;
    sflint::Fmpz block_prime;
    sflint::Fmpz large_prime_bound;
    sflint::fmpz_mul(sflint::FmpzRef(large_prime_bound),
                     sflint::FmpzConstRef(bound), sflint::FmpzConstRef(bound));
    slong found = 0;
    slong exact_entries = 0;
    for (slong a = -8; a <= 8; ++a) {
        for (slong b = -8; b <= 8; ++b) {
            for (slong c = -8; c <= 8; ++c) {
                if (a == 0 && b == 0 && c == 0) {
                    continue;
                }
                sflint::fmpz_mat_zero(sflint::FmpzMatRef(coordinates));
                sflint::fmpz_set_si(
                        sflint::fmpz_mat_entry(coordinates, 0, 0), a);
                sflint::fmpz_set_si(
                        sflint::fmpz_mat_entry(coordinates, 0, 1), b);
                sflint::fmpz_set_si(
                        sflint::fmpz_mat_entry(coordinates, 0, 2), c);
                assert(order_element.set_coordinates(
                        sflint::FmpzMatConstRef(coordinates)));
                assert(order_element.get_element(element));
                const bool have_row = silex::detail::
                        direct_partial_relation_row_for_testing(
                                row, large_prime, element, *base,
                                sflint::FmpzConstRef(bound));
                assert(element.norm(sflint::FmpqRef(norm)));
                const bool have_known_row = silex::detail::
                        direct_partial_relation_row_for_testing(
                                known_row, known_large_prime, element,
                                *base, sflint::FmpzConstRef(bound),
                                norm.raw(), coordinates.raw());
                assert(have_row == have_known_row);
                const bool have_multi_limb_row = silex::detail::
                        direct_partial_relation_row_for_testing(
                                multi_limb_row, multi_limb_large_prime,
                                element, *base, sflint::FmpzConstRef(bound),
                                nullptr, nullptr, false);
                assert(have_row == have_multi_limb_row);

                // The row exists exactly when the norm is a factor-base
                // smooth part times one prime at most bound^2 (the
                // single large prime of reference Rel_add.jl).
                sflint::fmpz_abs(sflint::FmpzRef(leftover),
                                 sflint::fmpq_num_ref(
                                         sflint::FmpqConstRef(norm)));
                for (slong block = 0;
                     block < base->rational_prime_block_count(); ++block) {
                    slong block_length = 0;
                    assert(base->rational_prime_block_data(
                            sflint::FmpzRef(block_prime), block_length,
                            block));
                    while (!sflint::fmpz_is_zero(
                                   sflint::FmpzConstRef(leftover)) &&
                           sflint::fmpz_divisible(
                                   sflint::FmpzConstRef(leftover),
                                   sflint::FmpzConstRef(block_prime))) {
                        sflint::fmpz_divexact(
                                sflint::FmpzRef(leftover),
                                sflint::FmpzConstRef(leftover),
                                sflint::FmpzConstRef(block_prime));
                    }
                }
                const bool expect_row =
                        sflint::fmpz_cmp(sflint::FmpzConstRef(leftover),
                                         sflint::FmpzConstRef(
                                                 large_prime_bound)) <= 0 &&
                        sflint::fmpz_is_prime(
                                sflint::FmpzConstRef(leftover));
                assert(have_row == expect_row);
                if (!have_row) {
                    continue;
                }
                ++found;
                assert(sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(row),
                        sflint::FmpzMatConstRef(known_row)));
                assert(sflint::fmpz_equal(
                        sflint::FmpzConstRef(large_prime),
                        sflint::FmpzConstRef(known_large_prime)));
                assert(sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(row),
                        sflint::FmpzMatConstRef(multi_limb_row)));
                assert(sflint::fmpz_equal(
                        sflint::FmpzConstRef(large_prime),
                        sflint::FmpzConstRef(multi_limb_large_prime)));
                assert(sflint::fmpz_equal(
                        sflint::FmpzConstRef(large_prime),
                        sflint::FmpzConstRef(leftover)));
                assert(sflint::fmpz_mat_ncols(row) == base->length());

                silex::Ideal principal(maximal_order);
                bool have_principal = false;
                for (slong i = 0; i < base->length(); ++i) {
                    const silex::PrimeIdeal* prime = base->prime_at(i);
                    assert(prime != nullptr);
                    slong valuation = -1;
                    if (!prime->valuation(valuation, order_element)) {
                        if (!have_principal) {
                            assert(principal.set_principal(order_element));
                            have_principal = true;
                        }
                        assert(prime->valuation(valuation, principal));
                    }
                    assert(sflint::fmpz_equal_si(
                            sflint::fmpz_mat_entry(row, 0, i), valuation));
                    if (valuation >= 2) {
                        // The residue-degree shortcut writes only 0 or 1.
                        ++exact_entries;
                    }
                }
            }
        }
    }
    assert(found > 0);
    assert(exact_entries > 0);
    return 0;
}

int test_compute_candidate_post_finite_refinement_cubic() {
    silex::NumberField field = cubic_field(-4, -7);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), maximal_order));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext baseline;
    assert(baseline.compute_candidate(maximal_order,
                                      sflint::FmpzConstRef(bound), options));
    assert(baseline.has_presentation());

    sflint::Fmpz baseline_order;
    assert(baseline.order(sflint::FmpzRef(baseline_order)));
    assert(sflint::fmpz_equal_si(baseline_order, 1));

    silex::ClassGroupContext refined;
    assert(refined.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(refined.has_presentation());
    assert(refined.relation_count() >= baseline.relation_count());
    assert(refined.relation_source_count(
                   silex::ClassGroupRelationSource::RandomProduct) >=
           baseline.relation_source_count(
                   silex::ClassGroupRelationSource::RandomProduct));

    sflint::Fmpz order_out;
    assert(refined.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_sgn(sflint::FmpzConstRef(order_out)) > 0);

    return 0;
}

int test_try_certify_with_units_zeta_quadratic() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation;
    silex::Order maximal(field);
    equation = silex::test::equation_order(field);
    assert(maximal.maximal_order(equation));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 512;
    options.max_relations = 64;
    sflint::Fmpz bound;
    assert(class_group_bound_at_least_two(bound, maximal));

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal, sflint::FmpzConstRef(bound),
                                     options));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.analytic_class_regulator_status() ==
           silex::ProofState::not_checked);
    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 5));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::proven, 192));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(context.analytic_class_regulator_status() ==
           silex::ProofState::verified);
    // Index one against the unconditional L(1, chi) hR proves saturation
    // at every prime, with no ell-local records kept.
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 0);
    assert(context.unit_proof_status() == silex::ProofState::verified);
    assert(context.regulator_proof_status() == silex::ProofState::verified);
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 5));

    return 0;
}

int test_compute_candidate_relation_kernel_target_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.relation_kernel_unit_count() >= 6);
    assert(context.relation_count() >=
           context.generator_count() + 6);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));

    silex::FactoredElement unit(*maximal_order.parent());
    assert(context.relation_kernel_unit(unit, 0));

    return 0;
}

int test_compute_candidate_relation_kernel_target_shifted_quadratic() {
    silex::NumberField field = shifted_quadratic_field(-3, 1);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.relation_kernel_unit_count() >= 6);
    assert(context.relation_count() >=
           context.generator_count() + 6);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) ==
           context.relation_count());

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));

    return 0;
}

int test_compute_candidate_relation_kernel_target_quintic() {
    silex::NumberField field = quintic_field(0, 0, -1, -1);
    silex::Order equation_order;
    silex::Order maximal_order(field);
    equation_order = silex::test::equation_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), maximal_order));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(maximal_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.has_presentation());
    assert(context.relation_kernel_unit_count() >= 2);
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));

    return 0;
}

int test_relation_saturation_proof_with_units_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    assert(context.has_presentation());
    assert(context.invariant_count() == 0);

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 0);

    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(aux_bound, 31));
    assert(!context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_record_count() == 0);

    sflint::Fmpz required_bound;
    assert(context.factor_base_generation_bound(sflint::FmpzRef(required_bound)));
    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(context.unit_proof_status() == silex::ProofState::verified);
    assert(context.regulator_proof_status() == silex::ProofState::verified);

    // A failed ell-local proof leaves the certification metadata unchanged.
    assert(set_fmpz_si(aux_bound, 2));
    assert(!context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_record_count() == 0);
    assert(context.unit_proof_status() == silex::ProofState::verified);
    assert(context.regulator_proof_status() == silex::ProofState::verified);

    assert(set_fmpz_si(aux_bound, 31));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);

    sflint::Fmpz record_ell;
    silex::ProofState status = silex::ProofState::not_checked;
    assert(context.relation_saturation_record(
            sflint::FmpzRef(record_ell), status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 2));
    assert(status == silex::ProofState::verified);
    auto saturation_record = context.relation_saturation_record(0);
    assert(saturation_record.has_value());
    assert(sflint::fmpz_equal_si(saturation_record->ell, 2));
    assert(saturation_record->status == silex::ProofState::verified);
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);

    assert(set_fmpz_si(ell, 4));
    assert(!context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);

    return 0;
}

int test_relation_saturation_bounded_append_with_units_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    assert(context.has_presentation());
    assert(context.relation_count() == 1);

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 3));
    bool changed = false;
    bool saturated = false;
    assert(context.saturate_relations_bounded_with_units(
            changed, saturated, units, sflint::FmpzConstRef(aux_bound), 1,
            1));
    assert(changed);
    assert(saturated);
    assert(context.has_presentation());
    assert(context.relation_count() == 2);
    assert(context.relation_source_count(
                   silex::ClassGroupRelationSource::Saturation) == 1);

    silex::ClassGroupRelationSource source =
            silex::ClassGroupRelationSource::Unknown;
    assert(context.relation_source(source, 1));
    assert(source == silex::ClassGroupRelationSource::Saturation);
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 1));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);

    sflint::Fmpz record_ell;
    silex::ProofState status = silex::ProofState::not_checked;
    assert(context.relation_saturation_record(
            sflint::FmpzRef(record_ell), status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 2));
    assert(status == silex::ProofState::verified);
    auto saturation_record = context.relation_saturation_record(0);
    assert(saturation_record.has_value());
    assert(sflint::fmpz_equal_si(saturation_record->ell, 2));
    assert(saturation_record->status == silex::ProofState::verified);
    assert(context.unit_proof_status() == silex::ProofState::verified);
    assert(context.regulator_proof_status() == silex::ProofState::verified);

    return 0;
}

// Stacks residue characters until the reference stopping rule holds.
void grow_candidate_filter_to_stop(
        silex::detail::RelationSaturationCandidateFilter& filter) {
    while (filter.available && !filter.stopped) {
        const slong before = filter.stacked_prime_count;
        assert(silex::detail::relation_saturation_candidate_filter_grow(filter));
        // One degree-one prime per call: the cost gate stacks no more.
        assert(filter.stacked_prime_count == before + 1);
    }
}

slong candidate_filter_kernel_dimension(
        const silex::detail::RelationSaturationCandidateFilter& filter) {
    return filter.width - static_cast<slong>(filter.rows.size());
}

// The stacked residue-character filter of the class-relation saturation
// with units.  Over Q with the single relation -4 and the torsion generator
// -1 at ell = 2, only (-4) * (-1) = 4 is a square, so the common kernel of the
// quadratic characters is spanned by (1, 1); a row outside it is ruled out.
// The stack stops once the dimension has stayed the same for more than
// stable * width = 1.5 * 2 primes; a failed power test doubles stable and
// lets it grow to more than 3 * 2.
int test_relation_saturation_candidate_filter_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.free_rank() == 0);

    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    silex::detail::RelationSaturationCandidateFilter filter;
    assert(silex::detail::relation_saturation_candidate_filter_start(
            filter, context, units, true, sflint::FmpzConstRef(ell)));
    assert(filter.available);
    assert(filter.width == 2);
    assert(filter.stable == 1.5);
    assert(filter.rows.empty() && !filter.stopped);

    // Rows (3, 1) = (1, 1) mod 2 and (1, 0), (0, 1) as exponent vectors.
    sflint::FmpzMat rows(3, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 1), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 1, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 2, 1), 1);
    // An empty stack rules nothing out.
    for (slong i = 0; i < 3; ++i) {
        bool ruled_out = true;
        assert(silex::detail::relation_saturation_candidate_filter_rules_out(
                ruled_out, filter, sflint::FmpzMatConstRef(rows), i));
        assert(!ruled_out);
    }

    grow_candidate_filter_to_stop(filter);
    assert(filter.stopped);
    assert(candidate_filter_kernel_dimension(filter) == 1);
    assert(!silex::detail::relation_saturation_candidate_filter_kernel_empty(
            filter));
    assert(filter.stable_count == 4);
    const bool expected_ruled_out[] = {false, true, true};
    for (slong i = 0; i < 3; ++i) {
        bool ruled_out = !expected_ruled_out[i];
        assert(silex::detail::relation_saturation_candidate_filter_rules_out(
                ruled_out, filter, sflint::FmpzMatConstRef(rows), i));
        assert(ruled_out == expected_ruled_out[i]);
    }

    silex::detail::relation_saturation_candidate_filter_failed_power_test(
            filter);
    assert(filter.stable == 3.0);
    assert(!filter.stopped);
    grow_candidate_filter_to_stop(filter);
    assert(filter.stable_count == 7);
    assert(candidate_filter_kernel_dimension(filter) == 1);
    return 0;
}

// The count gate of the filter: over Q with the relation -4 and the torsion
// generator -1 at ell = 2 there are two inputs, so the first two candidates
// are power tested without starting the filter.  The third starts it and
// stacks one prime; screening on until the stack stops rules out (1, 0),
// since -4 is not a square.
int test_relation_saturation_candidate_filter_count_gate() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    silex::OrderUnitGroup units;
    assert(units.compute(order));

    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    sflint::FmpzMat rows(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
    silex::detail::RelationSaturationCandidateFilter filter;
    bool ruled_out = true;
    for (slong tests = 1; tests <= 2; ++tests) {
        silex::detail::relation_saturation_candidate_filter_screen(
                ruled_out, filter, context, units, true,
                sflint::FmpzConstRef(ell), sflint::FmpzMatConstRef(rows), 0);
        assert(!ruled_out);
        assert(!filter.started);
        assert(filter.stacked_prime_count == 0);
        assert(filter.power_tests == tests);
    }
    silex::detail::relation_saturation_candidate_filter_screen(
            ruled_out, filter, context, units, true,
            sflint::FmpzConstRef(ell), sflint::FmpzMatConstRef(rows), 0);
    assert(filter.started && filter.available);
    assert(filter.width == 2);
    assert(filter.stacked_prime_count == 1);
    while (!ruled_out && !filter.stopped) {
        silex::detail::relation_saturation_candidate_filter_screen(
                ruled_out, filter, context, units, true,
                sflint::FmpzConstRef(ell), sflint::FmpzMatConstRef(rows), 0);
    }
    assert(ruled_out);
    return 0;
}

// The round memo of the saturation with units: a kernel row run without an
// appended relation is settled for the rest of the round, matched by its
// exact entries in any later per-prime kernel, and the filter's lazy start
// keeps it.  A row with a negative entry is never settled.
int test_relation_saturation_round_settled_rows() {
    silex::detail::RelationSaturationCandidateFilter filter;
    sflint::FmpzMat first(2, 3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(first, 0, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(first, 0, 2), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(first, 1, 1), 1);
    assert(!silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(first), 0));
    silex::detail::relation_saturation_round_settle_row(
            filter, sflint::FmpzMatConstRef(first), 0);
    silex::detail::relation_saturation_round_settle_row(
            filter, sflint::FmpzMatConstRef(first), 0);
    assert(filter.settled_rows.size() == 1);
    assert(silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(first), 0));
    assert(!silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(first), 1));

    // The same exponent vector as another prime's kernel row.
    sflint::FmpzMat second(2, 3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(second, 0, 1), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(second, 1, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(second, 1, 2), 1);
    assert(!silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(second), 0));
    assert(silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(second), 1));
    // Only the exact entries match: (3, 0, 1) is not (1, 0, 1).
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(second, 1, 0), 3);
    assert(!silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(second), 1));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(second, 1, 0), -1);
    silex::detail::relation_saturation_round_settle_row(
            filter, sflint::FmpzMatConstRef(second), 1);
    assert(filter.settled_rows.size() == 1);
    assert(!silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(second), 1));

    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    silex::OrderUnitGroup units;
    assert(units.compute(order));
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    assert(silex::detail::relation_saturation_candidate_filter_start(
            filter, context, units, true, sflint::FmpzConstRef(ell)));
    assert(silex::detail::relation_saturation_round_row_settled(
            filter, sflint::FmpzMatConstRef(first), 0));
    return 0;
}

// Q(sqrt(-23)) (x^2 = x - 6) has h = 3.  The factor base up to 3 with the
// relations generated by 4, 3, x and 1 + x gives h_cand = 6, so the exact
// index h_cand / h is 2.
void prepare_minus_23_index_two_candidate(silex::ClassGroupContext& context,
                                          const silex::Order& order) noexcept {
    const silex::NumberField* field = order.parent();
    assert(field != nullptr);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 3));
    context = silex::ClassGroupContext(order);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    // Generators a + b x with x^2 = x - 6: 4, 3, x, 1 + x.
    const slong generators[][2] = {{4, 0}, {3, 0}, {0, 1}, {1, 1}};
    for (const auto& generator : generators) {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, generator[0]);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, generator[1]);
        silex::Element alpha(*field);
        assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        silex::Relation relation(*base);
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
    }
    assert(context.publish_presentation());
    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 6));
}

#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
struct SettledRowCounter {
    long candidates = 0;
    long settled = 0;
};

void settled_row_callback(void* user,
                          silex::DiagnosticsModule,
                          silex::ProfileEvent event,
                          const char*,
                          const char* label) noexcept {
    if (event != silex::ProfileEvent::event || label == nullptr) {
        return;
    }
    auto* counter = static_cast<SettledRowCounter*>(user);
    if (std::strcmp(label,
                    "class_group.saturation_relation_candidate_with_units") ==
        0) {
        ++counter->candidates;
    } else if (std::strcmp(label,
                           "class_group.saturation_row_settled_with_units") ==
               0) {
        ++counter->settled;
    }
}

// The scan of the saturation with units skips a kernel row the round has
// already run.  In Z[sqrt 2] with the single relation -2 = -(sqrt 2)^2 the
// scan at ell = 2 meets the non-squares of the stacked kernel (the unit
// 1 + sqrt 2, -1, -2 alone) at several auxiliary primes before the square
// -2 * (-1) appends, so the skip event fires.
int test_relation_saturation_scan_skips_settled_rows() {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::OrderUnitGroup units;
    assert(units.compute(order));
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    SettledRowCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::class_group),
            settled_row_callback, &counter);

    silex::ClassGroupContext context(order);
    sflint::Fmpz factor_base_bound;
    assert(set_fmpz_si(factor_base_bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(factor_base_bound)));
    silex::Relation relation(*context.factor_base());
    silex::Element alpha(field);
    assert(alpha.set_si(-2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    context.set_diagnostics(&diagnostics);

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 9));
    bool changed = false;
    bool saturated = false;
    assert(context.saturate_relations_bounded_with_units(
            changed, saturated, units, sflint::FmpzConstRef(aux_bound), 1,
            1));
    assert(changed);
    assert(counter.candidates > 0);
    assert(counter.settled > 0);
    return 0;
}
#endif

// On the Q(sqrt(-23)) candidate above (h_cand = 6, h = 3) the relations are
// already 3-saturated, so at ell = 3 (no torsion column: w = 2) no nonzero
// exponent vector mod 3 makes the relation generators a cube and the stacked
// kernel is empty.  At ell = 2 the generator 4 is a square, so (1, 0, 0, 0, 0)
// lies in the kernel and 3 (not a square in K) does not.
int test_relation_saturation_candidate_filter_minus_23() {
    silex::NumberField field = shifted_quadratic_field(-1, 6);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.free_rank() == 0);
    silex::ClassGroupContext context;
    prepare_minus_23_index_two_candidate(context, order);

    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 3));
    silex::detail::RelationSaturationCandidateFilter filter;
    assert(silex::detail::relation_saturation_candidate_filter_start(
            filter, context, units, false, sflint::FmpzConstRef(ell)));
    assert(filter.available);
    assert(filter.width == 4);
    grow_candidate_filter_to_stop(filter);
    assert(silex::detail::relation_saturation_candidate_filter_kernel_empty(
            filter));
    // An empty kernel ends the stack; a failed test does not reopen it.
    silex::detail::relation_saturation_candidate_filter_failed_power_test(
            filter);
    assert(filter.stopped && filter.stable == 1.5);

    assert(set_fmpz_si(ell, 2));
    silex::detail::RelationSaturationCandidateFilter square_filter;
    assert(silex::detail::relation_saturation_candidate_filter_start(
            square_filter, context, units, true, sflint::FmpzConstRef(ell)));
    assert(square_filter.available);
    assert(square_filter.width == 5);
    grow_candidate_filter_to_stop(square_filter);
    const slong dimension = candidate_filter_kernel_dimension(square_filter);
    assert(dimension >= 1 && dimension < 5);
    sflint::FmpzMat rows(2, 5);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 1, 1), 1);
    bool ruled_out = true;
    assert(silex::detail::relation_saturation_candidate_filter_rules_out(
            ruled_out, square_filter, sflint::FmpzMatConstRef(rows), 0));
    assert(!ruled_out);
    assert(silex::detail::relation_saturation_candidate_filter_rules_out(
            ruled_out, square_filter, sflint::FmpzMatConstRef(rows), 1));
    assert(ruled_out);
    return 0;
}

// The filter is not used when the defining polynomial is not monic integral:
// the auxiliary primes are then not known to avoid the index.  Q given by
// 2x + 1 has a maximal order and a unit group, but the filter reports itself
// unavailable.
int test_relation_saturation_candidate_filter_needs_nice_polynomial() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    sflint::FmpqMat basis(1, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 0, 0), 1, 1);
    silex::Order base_order = silex::Order::from_basis(
            field, sflint::FmpqMatConstRef(basis));
    assert(base_order.is_defined());
    silex::Order order(field);
    assert(order.maximal_order(base_order));
    assert(order.is_maximal());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    silex::ClassGroupContext context(order);

    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    silex::detail::RelationSaturationCandidateFilter filter;
    filter.available = true;
    assert(silex::detail::relation_saturation_candidate_filter_start(
            filter, context, units, true, sflint::FmpzConstRef(ell)));
    assert(filter.started);
    assert(!filter.available);
    assert(silex::detail::relation_saturation_candidate_filter_grow(filter));
    assert(filter.stacked_prime_count == 0);
    return 0;
}

// The word-size degree-one quotient-log row equals the exact fmpz row, and
// fails in the same cases.  In Q(theta), theta^2 = theta - 6, the prime
// q = 13 = 1 mod 3 has the degree-one root 9; the inputs have negative
// exponents and rational coefficients.
int test_degree_one_quotient_log_row_word_size_matches_general() {
    silex::NumberField field = shifted_quadratic_field(-1, 6);
    auto element = [&field](slong c0, slong c1, slong den) {
        sflint::FmpqPoly poly;
        sflint::fmpq_poly_zero(poly);
        sflint::fmpq_poly_set_coeff_si(poly, 0, c0);
        sflint::fmpq_poly_set_coeff_si(poly, 1, c1);
        fmpq_poly_scalar_div_si(poly.raw(), poly.raw(), den);
        silex::Element out(field);
        assert(out.set_fmpq_poly(sflint::FmpqPolyConstRef(poly)));
        return out;
    };
    auto factored = [&field](std::vector<silex::FactoredElement>& out,
                             const silex::Element& a, slong ea,
                             const silex::Element& b, slong eb) {
        silex::FactoredElement value(field);
        assert(value.one());
        assert(value.push(a, ea));
        assert(value.push(b, eb));
        out.push_back(std::move(value));
    };

    std::vector<silex::FactoredElement> inputs;
    factored(inputs, element(1, 1, 2), -3, element(3, 0, 1), 1);
    factored(inputs, element(0, 1, 1), 2, element(-5, 1, 1), -1);
    factored(inputs, element(7, 0, 5), 1, element(2, 3, 7), -2);

    sflint::Fmpz q;
    sflint::Fmpz root;
    assert(set_fmpz_si(q, 13));
    assert(set_fmpz_si(root, 9));
    sflint::FmpzMat general(0, 0);
    sflint::FmpzMat word(0, 0);
    assert(silex::detail::degree_one_quotient_log_row_general(
            general, inputs, sflint::FmpzConstRef(q),
            sflint::FmpzConstRef(root), 3));
    assert(silex::detail::degree_one_quotient_log_row_word_size(
            word, inputs, sflint::FmpzConstRef(q),
            sflint::FmpzConstRef(root), 3));
    assert(sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(word)) == 3);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(general),
                                  sflint::FmpzMatConstRef(word)));
    bool nonzero = false;
    for (slong j = 0; j < 3; ++j) {
        nonzero = nonzero ||
                  !sflint::fmpz_is_zero(sflint::fmpz_mat_entry(
                          sflint::FmpzMatConstRef(word), 0, j));
    }
    assert(nonzero);

    // q divides a coefficient denominator, or a factor vanishes at the root:
    // both paths fail.
    std::vector<silex::FactoredElement> bad_denominator;
    factored(bad_denominator, element(0, 1, 13), -1, element(3, 0, 1), 1);
    std::vector<silex::FactoredElement> zero_value;
    factored(zero_value, element(-9, 1, 1), -2, element(3, 0, 1), 1);
    for (const auto* bad : {&bad_denominator, &zero_value}) {
        assert(!silex::detail::degree_one_quotient_log_row_general(
                general, *bad, sflint::FmpzConstRef(q),
                sflint::FmpzConstRef(root), 3));
        assert(!silex::detail::degree_one_quotient_log_row_word_size(
                word, *bad, sflint::FmpzConstRef(q),
                sflint::FmpzConstRef(root), 3));
    }
    return 0;
}

// The saturation candidate search skips an auxiliary prime q that divides
// lc(f) * disc(f) when the defining polynomial f is not monic integral.  In
// K = Q(theta), 3 theta^2 = 7 q0^2 with q0 = next_prime(2^60), the roots of f
// modulo q0 do not give degree-one primes: f = 3 x^2 mod q0, and the root 0
// would read 21 = (3 theta / q0)^2 as the non-square 21 mod q0 (q0 is inert in
// Q(sqrt(21))).  The search must not use q0, so the true square 21 stays in
// the candidate space.
int test_saturation_candidates_skip_bad_prime_of_nonmonic_polynomial() {
    sflint::Fmpz q0;
    fmpz_one(q0.raw());
    fmpz_mul_2exp(q0.raw(), q0.raw(), 60);
    sflint::fmpz_nextprime(sflint::FmpzRef(q0), sflint::FmpzConstRef(q0));
    sflint::Fmpz small;
    assert(set_fmpz_si(small, 21));
    assert(fmpz_jacobi(small.raw(), q0.raw()) == -1);
    assert(set_fmpz_si(small, -1));
    assert(fmpz_jacobi(small.raw(), q0.raw()) == 1);

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::Fmpz c0;
    fmpz_mul(c0.raw(), q0.raw(), q0.raw());
    fmpz_mul_si(c0.raw(), c0.raw(), -7);
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, c0.raw());
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 3);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    // Z[s] with s = 3 theta / q0, s^2 = 21, in theta coordinates.
    sflint::FmpqMat basis(2, 2);
    fmpq_mat_zero(basis.raw());
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 0, 0), 1, 1);
    fmpz_set_ui(fmpq_mat_entry_num(basis.raw(), 1, 1), 3);
    fmpz_set(fmpq_mat_entry_den(basis.raw(), 1, 1), q0.raw());
    silex::Order input = silex::Order::from_basis(
            field, sflint::FmpqMatConstRef(basis));
    assert(input.is_defined());
    silex::Order order(field);
    assert(order.maximal_order(input));
    assert(order.is_maximal());

    silex::Element value(field);
    assert(value.set_si(21));
    silex::FactoredElement relation(field);
    assert(relation.one());
    assert(relation.push(value, 1));
    std::vector<silex::FactoredElement> relations;
    relations.push_back(std::move(relation));

    // Inputs (21, -1): the field torsion generator is appended for ell = 2.
    sflint::FmpzMat candidates(0, 0);
    slong auxiliary_prime_count = 0;
    assert(silex::detail::saturation_candidates_from_relations_from_prime(
            candidates, auxiliary_prime_count, relations, order, 2, 100.0,
            sflint::FmpzConstRef(q0), 1));
    assert(auxiliary_prime_count == 1);
    const slong rows = sflint::fmpz_mat_nrows(
            sflint::FmpzMatConstRef(candidates));
    const slong cols = sflint::fmpz_mat_ncols(
            sflint::FmpzMatConstRef(candidates));
    assert(rows == 2);

    // (1, 0) lies in the column span modulo 2.
    nmod_mat_t span;
    nmod_mat_t augmented;
    nmod_mat_init(span, rows, cols, 2);
    nmod_mat_init(augmented, rows, cols + 1, 2);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            const ulong entry = fmpz_fdiv_ui(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(candidates), i, j)
                            .raw(),
                    2);
            nmod_mat_set_entry(span, i, j, entry);
            nmod_mat_set_entry(augmented, i, j, entry);
        }
    }
    nmod_mat_set_entry(augmented, 0, cols, 1);
    nmod_mat_set_entry(augmented, 1, cols, 0);
    const slong span_rank = nmod_mat_rank(span);
    const slong augmented_rank = nmod_mat_rank(augmented);
    nmod_mat_clear(augmented);
    nmod_mat_clear(span);
    assert(span_rank == augmented_rank);
    return 0;
}

// A defining polynomial with non-trivial content, f = q0 (x^2 - 21), has the
// integer numerator f itself.  q0 divides its content, so f is zero modulo q0
// and the root search would abort on it; the search must reject q0 and go on
// to the next prime q1.  21 is a nonresidue modulo q1 and a residue modulo the
// prime q2 after it, so q1 has no degree-one prime and q2 has two: the search
// from q0 must match the search from q1 and not the one from q2.
int test_saturation_candidates_skip_content_prime_of_nonmonic_polynomial() {
    sflint::Fmpz q0;
    fmpz_one(q0.raw());
    fmpz_mul_2exp(q0.raw(), q0.raw(), 60);
    sflint::fmpz_nextprime(sflint::FmpzRef(q0), sflint::FmpzConstRef(q0));

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::Fmpz c0;
    fmpz_mul_si(c0.raw(), q0.raw(), -21);
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, c0.raw());
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 2, q0.raw());
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    // Maximal order Z[(1 + theta) / 2], theta^2 = 21.
    sflint::FmpqMat basis(2, 2);
    fmpq_mat_zero(basis.raw());
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 1, 0), 1, 2);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 1, 1), 1, 2);
    silex::Order input = silex::Order::from_basis(
            field, sflint::FmpqMatConstRef(basis));
    assert(input.is_defined());
    silex::Order order(field);
    assert(order.maximal_order(input));
    assert(order.is_maximal());

    silex::Element value(field);
    assert(value.set_si(21));
    silex::FactoredElement relation(field);
    assert(relation.one());
    assert(relation.push(value, 1));
    std::vector<silex::FactoredElement> relations;
    relations.push_back(std::move(relation));

    sflint::Fmpz q1;
    sflint::Fmpz q2;
    sflint::Fmpz twenty_one;
    sflint::fmpz_nextprime(sflint::FmpzRef(q1), sflint::FmpzConstRef(q0));
    sflint::fmpz_nextprime(sflint::FmpzRef(q2), sflint::FmpzConstRef(q1));
    fmpz_set_si(twenty_one.raw(), 21);
    assert(fmpz_jacobi(twenty_one.raw(), q1.raw()) == -1);
    assert(fmpz_jacobi(twenty_one.raw(), q2.raw()) == 1);

    sflint::FmpzMat candidates(0, 0);
    slong auxiliary_prime_count = 0;
    silex::detail::SaturationCandidateSearchResult search;
    assert(silex::detail::saturation_candidates_from_relations_from_prime(
            candidates, auxiliary_prime_count, relations, order, 2, 100.0,
            sflint::FmpzConstRef(q0), 1, &search));
    assert(auxiliary_prime_count == 1);
    assert(search.auxiliary_prime_count == 1);
    assert(search.degree_one_prime_count == 0);

    sflint::FmpzMat q1_candidates(0, 0);
    slong q1_count = 0;
    silex::detail::SaturationCandidateSearchResult q1_search;
    assert(silex::detail::saturation_candidates_from_relations_from_prime(
            q1_candidates, q1_count, relations, order, 2, 100.0,
            sflint::FmpzConstRef(q1), 1, &q1_search));
    assert(q1_count == 1 && q1_search.degree_one_prime_count == 0);
    assert(sflint::fmpz_mat_equal(candidates, q1_candidates));

    sflint::FmpzMat q2_candidates(0, 0);
    slong q2_count = 0;
    silex::detail::SaturationCandidateSearchResult q2_search;
    assert(silex::detail::saturation_candidates_from_relations_from_prime(
            q2_candidates, q2_count, relations, order, 2, 100.0,
            sflint::FmpzConstRef(q2), 1, &q2_search));
    assert(q2_count == 1 && q2_search.degree_one_prime_count == 2);
    return 0;
}

// Index-bounded saturation (the internal class/unit validation pass) clears
// the primes of a bound on h_cand / h.  Clearing them proves saturation only
// when the bound holds unconditionally; with a GRH-conditional bound the
// status stays below `verified` unless verified ell-records cover every
// p | h_cand.  Q(sqrt(-23)) has h = 3; the relations below give h_cand = 6,
// so the bound 2 = h_cand / h is cleared and h_cand becomes 3, which no
// ell-local proof covers.
int test_relation_saturation_index_bounded_status_needs_unconditional_bound() {
    silex::NumberField field = shifted_quadratic_field(-1, 6);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        prepare_minus_23_index_two_candidate(context, order);
    };

    sflint::Fmpz index_bound;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(index_bound, 2));
    assert(set_fmpz_si(aux_bound, 101));
    for (const bool unconditional : {false, true}) {
        silex::ClassGroupContext context;
        prepare_context(context);
        bool changed = false;
        bool saturated = false;
        assert(CertificationAccess::
                       saturate_relations_for_index_bound_with_units(
                               changed, saturated, context, units,
                               sflint::FmpzConstRef(index_bound),
                               unconditional,
                               sflint::FmpzConstRef(aux_bound), 4, 8));
        // The bounded search cleared ell = 2 from h_cand.
        assert(changed);
        assert(saturated);
        sflint::Fmpz order_out;
        assert(context.order(sflint::FmpzRef(order_out)));
        assert(sflint::fmpz_equal_si(order_out, 3));
        assert(context.relation_saturation_record_count() == 1);
        assert(context.relation_saturation_status() ==
               (unconditional ? silex::ProofState::verified
                              : silex::ProofState::unavailable));
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
    }
    return 0;
}

// The validation pass treats the exact imaginary-quadratic index as an
// unconditional bound exactly when factor-base generation is verified up to
// the required bound, whatever the request.  On the h_cand = 6, h = 3
// candidate above, clearing ell = 2 gives `verified` when generation is
// checked to the required bound, under a grh request as under a proven one,
// and `unavailable` when it is not checked, under either request.
int test_validation_exact_index_saturation_status_follows_generation() {
    silex::NumberField field = shifted_quadratic_field(-1, 6);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz index_bound;
    assert(set_fmpz_si(index_bound, 2));
    for (const silex::CertificationMode requested :
         {silex::CertificationMode::grh, silex::CertificationMode::proven}) {
        for (const bool check_generation : {false, true}) {
            silex::ClassGroupContext context;
            prepare_minus_23_index_two_candidate(context, order);
            if (check_generation) {
                sflint::Fmpz required_bound;
                assert(context.factor_base_generation_bound(
                        sflint::FmpzRef(required_bound)));
                assert(context.check_factor_base_generation_bound(
                        sflint::FmpzConstRef(required_bound)));
            }
            assert(context.factor_base_generation_checked_status() ==
                   (check_generation ? silex::ProofState::verified
                                     : silex::ProofState::not_checked));

            silex::ClassGroupComputeOptions options;
            options.requested_certification = requested;
            options.max_relations = 64;
            bool changed = false;
            assert(silex::detail::saturate_candidate_class_relations_with_units(
                    changed, context, units, options,
                    sflint::FmpzConstRef(index_bound), true));
            assert(changed);
            sflint::Fmpz order_out;
            assert(context.order(sflint::FmpzRef(order_out)));
            assert(sflint::fmpz_equal_si(order_out, 3));
            assert(context.relation_saturation_status() ==
                   (check_generation ? silex::ProofState::verified
                                     : silex::ProofState::unavailable));
            assert(context.certification_status() ==
                   silex::CertificationMode::unknown);
        }
    }
    return 0;
}

int test_relation_saturation_index_bound_with_units_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.is_defined());
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element alpha(field);
        assert(alpha.set_si(-2));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        assert(context.has_presentation());
        assert(context.invariant_count() == 0);

        sflint::Fmpz required_bound;
        assert(context.factor_base_generation_bound(
                sflint::FmpzRef(required_bound)));
        assert(context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(required_bound)));
    };

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    sflint::Fmpz index_bound;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 31));

    silex::ClassGroupContext invalid_context;
    prepare_context(invalid_context);
    assert(invalid_context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(set_fmpz_si(index_bound, 1));
    assert(!CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           invalid_context, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(invalid_context.relation_saturation_status() ==
           silex::ProofState::not_checked);

    // Over Q (h_cand = 1, w = 2, unit rank 0) the primes 3 and 5 below the
    // bound 6 have target rank 0: they divide neither h_cand nor w, so the
    // index argument proves them with nothing to check.  Each records a
    // verified ell-local proof with rank 0 and target 0, and ell = 2 is
    // proved through the torsion.
    silex::ClassGroupContext rank_zero_context;
    prepare_context(rank_zero_context);
    assert(rank_zero_context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(set_fmpz_si(index_bound, 6));
    assert(CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           rank_zero_context, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(rank_zero_context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(rank_zero_context.relation_saturation_record_count() == 3);
    assert(rank_zero_context.certification_status() ==
           silex::CertificationMode::proven);
    for (const slong small_ell : {3, 5}) {
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, small_ell));
        silex::ProofState proof_status = silex::ProofState::not_checked;
        slong rank = -1;
        slong target_rank = -1;
        slong local_primes = -1;
        assert(CertificationAccess::relation_saturation_proof_record(
                rank_zero_context, sflint::FmpzConstRef(ell), proof_status,
                rank, target_rank, local_primes));
        assert(proof_status == silex::ProofState::verified);
        assert(rank == 0);
        assert(target_rank == 0);
        assert(local_primes == 0);
    }

    silex::ClassGroupContext proven_context;
    prepare_context(proven_context);
    assert(proven_context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(set_fmpz_si(index_bound, 2));
    assert(CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           proven_context, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(proven_context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_context.relation_saturation_record_count() == 1);
    assert(proven_context.certification_status() ==
           silex::CertificationMode::proven);

    sflint::Fmpz record_ell;
    silex::ProofState status = silex::ProofState::not_checked;
    assert(proven_context.relation_saturation_record(
            sflint::FmpzRef(record_ell), status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 2));
    assert(status == silex::ProofState::verified);
    assert(proven_context.try_certify_with_units(
            units, silex::CertificationMode::proven, 80));

    return 0;
}

// Q(sqrt(-23)) has class number 3.  With the factor base of the primes above
// 2 and the relations 2 and ((3 + sqrt(-23))/2)^k, the candidate order is 3k.
// For k = 3 the candidate 9 is unsaturated at 3, so no sequence of saturation
// proofs that omits ell = 3 may certify it.
void prepare_minus_23_candidate(silex::ClassGroupContext& context,
                                const silex::Order& maximal,
                                slong beta_power,
                                slong expected_order,
                                bool check_generation = true) noexcept {
    const silex::NumberField* field = maximal.parent();
    assert(field != nullptr);
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                maximal));
    context = silex::ClassGroupContext(maximal);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Element two(*field);
    silex::Element half(*field);
    silex::Element beta(*field);
    silex::Element beta_power_element(*field);
    assert(two.set_si(2));
    assert(half.set_si_over_si(1, 2));
    assert(beta.gen());
    assert(beta.add_si(beta, 3));
    assert(beta.multiply(beta, half));
    sflint::Fmpz exponent;
    assert(set_fmpz_si(exponent, beta_power));
    assert(beta_power_element.pow_fmpz(beta, sflint::FmpzConstRef(exponent)));

    silex::Relation two_relation(*base);
    assert(two_relation.set_generator(two));
    assert(context.append_relation(two_relation));
    silex::Relation beta_relation(*base);
    assert(beta_relation.set_generator(beta_power_element));
    assert(context.append_relation(beta_relation));
    assert(context.publish_presentation());

    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, expected_order));
    if (!check_generation) {
        return;
    }

    sflint::Fmpz required_bound;
    assert(context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));
}

// A failed exact imaginary-quadratic certification leaves the
// factor-base generation check exactly as it found it.  With h_cand = 9 and
// the exact class number 3 the route runs the generation check and then
// rejects the candidate, so the check it ran must be rolled back.
int test_exact_order_certification_failure_restores_generation_check() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 3, 9, false);
    assert(context.factor_base_generation_status() ==
           silex::ProofState::verified);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);
    assert(context.factor_base_generation_record_count() == 0);

    sflint::Fmpz discriminant;
    sflint::Fmpz exact_order;
    sflint::Fmpz checked_bound;
    assert(set_fmpz_si(discriminant, -23));
    assert(set_fmpz_si(exact_order, 3));
    assert(!CertificationAccess::
                   try_certify_imaginary_quadratic_from_exact_order(
                           context, silex::CertificationMode::proven,
                           sflint::FmpzConstRef(discriminant),
                           sflint::FmpzConstRef(exact_order)));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.factor_base_generation_checked_status() ==
           silex::ProofState::not_checked);
    assert(!context.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));
    assert(sflint::fmpz_equal_si(checked_bound, 0));
    assert(context.factor_base_generation_record_count() == 0);
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_record_count() == 0);

    return 0;
}

int test_saturation_promotion_requires_every_class_order_prime() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    sflint::Fmpz index_bound;
    assert(set_fmpz_si(aux_bound, 200));

    // A failed proven request leaves every proof status untouched.
    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 3, 9);
    assert(!context.try_certify_with_units(
            units, silex::CertificationMode::proven, 128));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.unit_proof_status() == silex::ProofState::not_checked);
    assert(context.regulator_proof_status() ==
           silex::ProofState::not_checked);
    assert(context.analytic_class_regulator_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);

    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));
    assert(context.unit_proof_status() == silex::ProofState::verified);

    // The candidate is not saturated at 3.
    assert(set_fmpz_si(ell, 3));
    assert(!context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(context.relation_saturation_record_count() == 0);

    // Saturation at 2 is provable but does not cover 3 | h_cand = 9.
    assert(set_fmpz_si(ell, 2));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_record_count() == 1);
    assert(context.relation_saturation_status() !=
           silex::ProofState::verified);
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);

    assert(!context.try_certify_with_units(
            units, silex::CertificationMode::proven, 128));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.relation_saturation_record_count() == 1);
    assert(context.unit_proof_status() == silex::ProofState::verified);

    // A supplied index bound below 3 does not drop the prime 3 | h_cand.
    silex::ClassGroupContext bounded;
    prepare_minus_23_candidate(bounded, maximal, 3, 9);
    assert(bounded.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));
    assert(set_fmpz_si(index_bound, 2));
    assert(!CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           bounded, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(bounded.certification_status() ==
           silex::CertificationMode::unknown);
    assert(bounded.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(bounded.relation_saturation_record_count() == 0);
    assert(!bounded.try_certify_with_units(
            units, silex::CertificationMode::proven, 128));
    assert(bounded.certification_status() ==
           silex::CertificationMode::unknown);

    // The true class group h = 3 is certified by saturation at 2 and 3.
    silex::ClassGroupContext exact;
    prepare_minus_23_candidate(exact, maximal, 1, 3);
    assert(exact.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));
    assert(CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           exact, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(exact.certification_status() == silex::CertificationMode::proven);
    assert(exact.relation_saturation_status() == silex::ProofState::verified);
    assert(exact.relation_saturation_record_count() == 2);
    assert(exact.analytic_class_regulator_status() ==
           silex::ProofState::not_checked);

    return 0;
}

// T-010 follow-up: saturation completes ell-by-ell as each required prime is
// proven individually, not only through the one-shot index-bound helper
// exercised above.  Covers h_cand = 3 (ell = 2 then 3) and the k = 5 case
// h_cand = 15 = 3 * 5 (ell = 3 then 5).  Each low-level
// try_prove_relation_saturation_with_units call never promotes the class
// group itself (class_group.cpp:9937-9958); only a subsequent certification
// call, here the saturation-only CertificationAccess helper, does.
int test_relation_saturation_completes_prime_by_prime() {
    // h_cand = 3: Q(sqrt(-23)) has class number 3.
    {
        silex::NumberField field = quadratic_field(-23);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));

        silex::OrderUnitGroup units;
        assert(units.compute(maximal));
        assert(units.certification_status() ==
               silex::CertificationMode::proven);

        silex::ClassGroupContext context;
        prepare_minus_23_candidate(context, maximal, 1, 3);
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 128));
        assert(context.unit_proof_status() == silex::ProofState::verified);

        sflint::Fmpz aux_bound;
        assert(set_fmpz_si(aux_bound, 200));
        sflint::Fmpz ell;

        // ell = 2 does not divide h_cand = 3, but is itself provable; one
        // record is not enough to complete the proof.
        assert(set_fmpz_si(ell, 2));
        assert(context.try_prove_relation_saturation_with_units(
                units, sflint::FmpzConstRef(ell),
                sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 1);
        assert(context.relation_saturation_status() !=
               silex::ProofState::verified);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);

        // ell = 3 completes coverage of every prime dividing h_cand.
        assert(set_fmpz_si(ell, 3));
        assert(context.try_prove_relation_saturation_with_units(
                units, sflint::FmpzConstRef(ell),
                sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 2);
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        // Completing the per-ell records does not itself promote the class
        // group; only a certification call does.
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
        assert(CertificationAccess::try_prove_class_order_saturation_with_units(
                context, units, sflint::FmpzConstRef(aux_bound)));
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
    }

    // The k = 5 case, h_cand = 15 = 3 * 5: Q(sqrt(-239)) has class number 15
    // (cyclic).
    {
        silex::NumberField field = quadratic_field(-239);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));

        silex::ClassGroupCandidateOptions options;
        options.max_candidates = 20000;
        options.max_relations = 2000;
        sflint::Fmpz factor_base_bound;
        assert(class_group_bound_at_least_two(factor_base_bound, maximal));

        silex::ClassGroupContext context;
        assert(context.compute_candidate(
                maximal, sflint::FmpzConstRef(factor_base_bound), options));
        sflint::Fmpz order_out;
        assert(context.order(sflint::FmpzRef(order_out)));
        assert(sflint::fmpz_equal_si(order_out, 15));

        sflint::Fmpz required_bound;
        assert(context.factor_base_generation_bound(
                sflint::FmpzRef(required_bound)));
        assert(context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(required_bound)));

        silex::OrderUnitGroup units;
        assert(units.compute(maximal));
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 128));
        assert(context.unit_proof_status() == silex::ProofState::verified);

        sflint::Fmpz aux_bound;
        assert(set_fmpz_si(aux_bound, 200));
        sflint::Fmpz ell;

        // ell = 3 alone does not cover h_cand = 15 = 3 * 5.
        assert(set_fmpz_si(ell, 3));
        assert(context.try_prove_relation_saturation_with_units(
                units, sflint::FmpzConstRef(ell),
                sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 1);
        assert(context.relation_saturation_status() !=
               silex::ProofState::verified);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);

        // ell = 5 completes coverage.
        assert(set_fmpz_si(ell, 5));
        assert(context.try_prove_relation_saturation_with_units(
                units, sflint::FmpzConstRef(ell),
                sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 2);
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
        assert(CertificationAccess::try_prove_class_order_saturation_with_units(
                context, units, sflint::FmpzConstRef(aux_bound)));
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
    }

    return 0;
}

// The two saturation-proof entry points differ, and these differences are
// pinned here.  On Q(sqrt(-23)) with h_cand = 3:
// - The class-order path runs the factor-base generation check itself and
//   sets the unit and regulator proof states.  The index-bound path requires
//   both to be verified already and leaves them unchanged when they are not.
// - The class-order path proves only the primes dividing h_cand.  The
//   index-bound path also proves every prime up to its bound.
// - The class-order path keeps an ell that is already verified.  The
//   index-bound path proves it again by the ell-local proof.
// - A failed class-order attempt rolls back the states it set.
int test_saturation_proof_paths_differ() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz aux_bound;
    sflint::Fmpz index_bound;
    sflint::Fmpz ell;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(index_bound, 2));
    silex::ProofState proof_status = silex::ProofState::not_checked;
    slong rank = -1;
    slong target_rank = -1;
    slong local_primes = -1;

    // Preconditions and the states each path sets.
    {
        silex::ClassGroupContext context;
        prepare_minus_23_candidate(context, maximal, 1, 3, false);
        assert(context.unit_proof_status() == silex::ProofState::not_checked);
        assert(context.factor_base_generation_checked_status() ==
               silex::ProofState::not_checked);

        assert(!CertificationAccess::
                       try_prove_relation_saturation_index_bound_with_units(
                               context, units,
                               sflint::FmpzConstRef(index_bound),
                               sflint::FmpzConstRef(aux_bound)));
        assert(context.unit_proof_status() == silex::ProofState::not_checked);
        assert(context.regulator_proof_status() ==
               silex::ProofState::not_checked);
        assert(context.factor_base_generation_checked_status() ==
               silex::ProofState::not_checked);
        assert(context.relation_saturation_status() ==
               silex::ProofState::not_checked);
        assert(context.relation_saturation_record_count() == 0);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);

        assert(CertificationAccess::try_prove_class_order_saturation_with_units(
                context, units, sflint::FmpzConstRef(aux_bound)));
        assert(context.unit_proof_status() == silex::ProofState::verified);
        assert(context.regulator_proof_status() ==
               silex::ProofState::verified);
        assert(context.factor_base_generation_checked_status() ==
               silex::ProofState::verified);
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        // Only ell = 3 | h_cand; no prime below an index bound.
        assert(context.relation_saturation_record_count() == 1);
        assert(set_fmpz_si(ell, 2));
        assert(!CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(ell), proof_status, rank,
                target_rank, local_primes));
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
        assert(context.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
    }

    // The index-bound path, given the same prerequisites, also proves the
    // primes up to its bound.
    {
        silex::ClassGroupContext context;
        prepare_minus_23_candidate(context, maximal, 1, 3);
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 128));
        assert(context.unit_proof_status() == silex::ProofState::verified);
        assert(CertificationAccess::
                       try_prove_relation_saturation_index_bound_with_units(
                               context, units,
                               sflint::FmpzConstRef(index_bound),
                               sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(context.relation_saturation_record_count() == 2);
        assert(set_fmpz_si(ell, 2));
        assert(CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(ell), proof_status, rank,
                target_rank, local_primes));
        assert(proof_status == silex::ProofState::verified);
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
    }

    // A failed class-order attempt (h_cand = 9 is unsaturated at 3) rolls
    // back the generation check, the unit and regulator states, and the
    // records it wrote.
    {
        silex::ClassGroupContext context;
        prepare_minus_23_candidate(context, maximal, 3, 9, false);
        assert(!CertificationAccess::try_prove_class_order_saturation_with_units(
                context, units, sflint::FmpzConstRef(aux_bound)));
        assert(context.unit_proof_status() == silex::ProofState::not_checked);
        assert(context.regulator_proof_status() ==
               silex::ProofState::not_checked);
        assert(context.factor_base_generation_checked_status() ==
               silex::ProofState::not_checked);
        assert(context.relation_saturation_status() ==
               silex::ProofState::not_checked);
        assert(context.relation_saturation_record_count() == 0);
        assert(set_fmpz_si(ell, 3));
        assert(!CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(ell), proof_status, rank,
                target_rank, local_primes));
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
    }

    // An ell already verified by the exact class number (an exact-order
    // proof record with local_primes = 0) is kept by the class-order path
    // and proved again by the ell-local proof on the index-bound path.
    sflint::Fmpz discriminant;
    sflint::Fmpz exact_order;
    assert(set_fmpz_si(discriminant, -23));
    assert(set_fmpz_si(exact_order, 3));
    assert(set_fmpz_si(ell, 3));
    for (const bool index_bound_path : {false, true}) {
        silex::ClassGroupContext context;
        prepare_minus_23_candidate(context, maximal, 1, 3);
        assert(CertificationAccess::
                       try_certify_imaginary_quadratic_from_exact_order(
                               context, silex::CertificationMode::proven,
                               sflint::FmpzConstRef(discriminant),
                               sflint::FmpzConstRef(exact_order)));
        assert(CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(ell), proof_status, rank,
                target_rank, local_primes));
        assert(proof_status == silex::ProofState::verified);
        assert(local_primes == 0);

        if (index_bound_path) {
            assert(CertificationAccess::
                           try_prove_relation_saturation_index_bound_with_units(
                                   context, units,
                                   sflint::FmpzConstRef(index_bound),
                                   sflint::FmpzConstRef(aux_bound)));
        } else {
            assert(CertificationAccess::
                           try_prove_class_order_saturation_with_units(
                                   context, units,
                                   sflint::FmpzConstRef(aux_bound)));
        }
        assert(CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(ell), proof_status, rank,
                target_rank, local_primes));
        assert(proof_status == silex::ProofState::verified);
        assert(index_bound_path ? local_primes > 0 : local_primes == 0);
        assert(context.relation_saturation_record_count() ==
               (index_bound_path ? 2 : 1));
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
    }

    return 0;
}

// try_certify_class_unit_with_zeta keeps complete ell-local saturation
// records.  On Q(sqrt(-23)), h_cand = 3, a saturation proof at ell = 3
// covers every prime dividing h_cand; the later unconditional L(1, chi)
// index-one check then certifies `proven` and leaves those verified records
// and the `verified` status in place.
int test_class_unit_zeta_keeps_complete_saturation_records() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 1, 3);
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));

    sflint::Fmpz aux_bound;
    sflint::Fmpz ell;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(ell, 3));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);

    assert(context.try_certify_class_unit_with_zeta(units, 128));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.analytic_class_regulator_certification() ==
           silex::CertificationMode::proven);
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);
    sflint::Fmpz record_ell;
    silex::ProofState record_status = silex::ProofState::not_checked;
    assert(context.relation_saturation_record(sflint::FmpzRef(record_ell),
                                              record_status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 3));
    assert(record_status == silex::ProofState::verified);
    return 0;
}

// The unconditional analytic index-one route of a `proven`
// try_certify_with_units keeps complete ell-local saturation records.  On
// Q(sqrt(-23)), h_cand = 3, generation is checked below the required bound,
// which is enough for a saturation proof at ell = 3 but not for the early
// promotion, so try_certify_with_units falls through to the L(1, chi)
// index-one check.  That check raises the generation bound, certifies
// `proven`, and leaves the verified ell = 3 record in place.
int test_try_certify_with_units_analytic_keeps_complete_saturation_records() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 1, 3, false);
    sflint::Fmpz required_bound;
    assert(context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(sflint::fmpz_cmp_ui(sflint::FmpzConstRef(required_bound), 1) > 0);
    sflint::Fmpz low_bound;
    sflint::fmpz_sub_ui(sflint::FmpzRef(low_bound),
                        sflint::FmpzConstRef(required_bound), 1);
    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(low_bound)));
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));

    sflint::Fmpz aux_bound;
    sflint::Fmpz ell;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(ell, 3));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);

    assert(context.try_certify_with_units(
            units, silex::CertificationMode::proven, 128));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.analytic_class_regulator_certification() ==
           silex::CertificationMode::proven);
    sflint::Fmpz checked_bound;
    assert(context.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));
    assert(sflint::fmpz_cmp(sflint::FmpzConstRef(checked_bound),
                            sflint::FmpzConstRef(required_bound)) >= 0);
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);
    sflint::Fmpz record_ell;
    silex::ProofState record_status = silex::ProofState::not_checked;
    assert(context.relation_saturation_record(sflint::FmpzRef(record_ell),
                                              record_status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 3));
    assert(record_status == silex::ProofState::verified);
    return 0;
}

int test_relation_saturation_index_bound_checks_nondivisor_primes() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 512;
    options.max_relations = 64;
    sflint::Fmpz factor_base_bound;
    assert(class_group_bound_at_least_two(factor_base_bound, maximal));

    silex::ClassGroupContext context;
    assert(context.compute_candidate(
            maximal, sflint::FmpzConstRef(factor_base_bound), options));
    sflint::Fmpz order_out;
    assert(context.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 5));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::proven, 192));

    sflint::Fmpz ell_five;
    sflint::Fmpz index_bound;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(ell_five, 5));
    assert(set_fmpz_si(index_bound, 5));
    assert(set_fmpz_si(aux_bound, 101));
    // The only prime factor of B=5 is individually provable.
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell_five),
            sflint::FmpzConstRef(aux_bound)));
    // A bound I <= 5 also requires the primes 2 and 3 that do not divide
    // h_cand = 5.  ell = 2 divides w = 2 and is proved through the torsion;
    // ell = 3 divides neither h_cand nor w and the unit rank is 0, so its
    // target rank is 0 and it records a verified proof with rank 0 and
    // target 0.
    assert(CertificationAccess::
                   try_prove_relation_saturation_index_bound_with_units(
                           context, units, sflint::FmpzConstRef(index_bound),
                           sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 3);
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    sflint::Fmpz ell_three;
    assert(set_fmpz_si(ell_three, 3));
    silex::ProofState proof_status = silex::ProofState::not_checked;
    slong rank = -1;
    slong target_rank = -1;
    slong local_primes = -1;
    assert(CertificationAccess::relation_saturation_proof_record(
            context, sflint::FmpzConstRef(ell_three), proof_status, rank,
            target_rank, local_primes));
    assert(proof_status == silex::ProofState::verified);
    assert(rank == 0);
    assert(target_rank == 0);
    assert(local_primes == 0);

    return 0;
}

int test_analytic_class_unit_proof_requires_factor_base_generation() {
    silex::NumberField field = quadratic_field(210);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 10000;
    options.max_relations = 1000;

    sflint::Fmpz factor_base_bound;
    assert(set_fmpz_si(factor_base_bound, 10));
    silex::ClassGroupContext context;
    assert(context.compute_candidate(
            maximal, sflint::FmpzConstRef(factor_base_bound), options));
    assert(context.has_presentation());
    assert(context.factor_base_generation_status() ==
           silex::ProofState::unavailable);

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() ==
           silex::CertificationMode::proven);

    // Give the combined proof path a synthetic candidate product whose ratio
    // is exactly one.  It still may not publish the tentative class group
    // until the factor base is proved to generate the full class group.
    sflint::Arb analytic_hR;
    assert(units.class_regulator_product(
            sflint::ArbRef(analytic_hR), context, 192));
    // Pass the flag the route itself reports for this maximal quadratic
    // order, not a literal `true` (T-012): the quadratic L(1, chi) route can
    // fall back to Belabas-Friedman, so whether a held value is unconditional
    // depends on the route that produced it
    // (the `unconditional` output of
    // zeta_class_regulator_product_with_diagnostics, declared in
    // src/zeta/zeta_internal.hpp).
    sflint::Arb reported_hR;
    bool hr_unconditional = false;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(reported_hR), maximal, 192, nullptr, nullptr,
            nullptr, &hr_unconditional));
    assert(!CertificationAccess::try_certify_class_unit_with_units(
            context, units, sflint::ArbConstRef(analytic_hR), 192,
            hr_unconditional));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);
    assert(context.analytic_class_regulator_status() ==
           silex::ProofState::not_checked);
    assert(context.unit_proof_status() == silex::ProofState::not_checked);
    assert(context.regulator_proof_status() ==
           silex::ProofState::not_checked);

    return 0;
}

int test_relation_saturation_analytic_index_bound_with_units_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.is_defined());
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element alpha(field);
        assert(alpha.set_si(-2));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        assert(context.has_presentation());
        assert(context.invariant_count() == 0);
    };

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    // The flag the degree-one route itself reports (T-012: use the reported
    // route, not a literal `true`); degree one has no zeta computation and
    // is always unconditional (src/zeta/zeta.cpp:2401-2402).
    sflint::Arb reported_hR;
    bool hr_unconditional = false;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(reported_hR), order, 128, nullptr, nullptr,
            nullptr, &hr_unconditional));
    assert(hr_unconditional);

    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb analytic_hR;
    sflint::arb_set_fmpq(analytic_hR, half, 128);

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 31));

    silex::ClassGroupContext invalid_analytic_context;
    prepare_context(invalid_analytic_context);
    sflint::Arb zero;
    sflint::arb_zero(zero);
    assert(!CertificationAccess::try_analytic_index_bound_with_units(
            invalid_analytic_context, units, sflint::ArbConstRef(zero),
            sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
    assert(invalid_analytic_context.unit_proof_status() ==
           silex::ProofState::not_checked);
    assert(invalid_analytic_context.regulator_proof_status() ==
           silex::ProofState::not_checked);
    assert(invalid_analytic_context.relation_saturation_status() ==
           silex::ProofState::not_checked);

    silex::ClassGroupContext unavailable_context;
    prepare_context(unavailable_context);
    sflint::Fmpz required_bound;
    assert(unavailable_context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(unavailable_context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));
    assert(set_fmpz_si(aux_bound, 2));
    assert(!CertificationAccess::try_analytic_index_bound_with_units(
            unavailable_context, units, sflint::ArbConstRef(analytic_hR),
            sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
    assert(unavailable_context.unit_proof_status() ==
           silex::ProofState::not_checked);
    assert(unavailable_context.regulator_proof_status() ==
           silex::ProofState::not_checked);
    assert(unavailable_context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(unavailable_context.relation_saturation_record_count() == 0);
    assert(unavailable_context.certification_status() ==
           silex::CertificationMode::unknown);

    silex::ClassGroupContext index_one_context;
    prepare_context(index_one_context);
    assert(index_one_context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(index_one_context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));
    sflint::Arb exact_hR;
    sflint::arb_one(exact_hR);
    assert(CertificationAccess::try_analytic_index_bound_with_units(
            index_one_context, units, sflint::ArbConstRef(exact_hR),
            sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
    assert(index_one_context.unit_proof_status() ==
           silex::ProofState::verified);
    assert(index_one_context.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(index_one_context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(index_one_context.relation_saturation_record_count() == 0);
    assert(index_one_context.certification_status() ==
           silex::CertificationMode::proven);

    silex::ClassGroupContext nontrivial_bound_context;
    prepare_context(nontrivial_bound_context);
    assert(nontrivial_bound_context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(nontrivial_bound_context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));
    assert(set_fmpz_si(aux_bound, 31));
    assert(!CertificationAccess::try_analytic_index_bound_with_units(
            nontrivial_bound_context, units, sflint::ArbConstRef(analytic_hR),
            sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
    assert(nontrivial_bound_context.unit_proof_status() ==
           silex::ProofState::not_checked);
    assert(nontrivial_bound_context.regulator_proof_status() ==
           silex::ProofState::not_checked);
    assert(nontrivial_bound_context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(nontrivial_bound_context.relation_saturation_record_count() == 0);
    assert(nontrivial_bound_context.certification_status() ==
           silex::CertificationMode::unknown);

    return 0;
}

// In degree one hR = 1 exactly, so both Belabas-Friedman entry points prove
// the class/unit pair by the unconditional analytic index-one argument and
// report relation saturation as `verified` with no per-prime records, as
// try_certify_class_unit_with_zeta does on the same order.
int test_zeta_bf_degree_one_verifies_relation_saturation() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.degree() == 1);

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element alpha(field);
        assert(alpha.set_si(-2));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        sflint::Fmpz required_bound;
        assert(context.factor_base_generation_bound(
                sflint::FmpzRef(required_bound)));
        assert(context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(required_bound)));
        assert(context.relation_saturation_status() ==
               silex::ProofState::not_checked);
    };
    auto assert_proven_and_saturated =
            [](const silex::ClassGroupContext& context,
               const silex::OrderUnitGroup& units) noexcept {
                assert(context.certification_status() ==
                       silex::CertificationMode::proven);
                assert(context.analytic_class_regulator_certification() ==
                       silex::CertificationMode::proven);
                assert(context.relation_saturation_status() ==
                       silex::ProofState::verified);
                assert(context.relation_saturation_record_count() == 0);
                assert(units.certification_status() ==
                       silex::CertificationMode::proven);
            };

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    silex::ClassGroupContext zeta_context;
    prepare_context(zeta_context);
    assert(zeta_context.try_certify_class_unit_with_zeta(units, 64));
    assert_proven_and_saturated(zeta_context, units);

    silex::ClassGroupContext bf_context;
    prepare_context(bf_context);
    assert(bf_context.try_certify_class_unit_with_zeta_bf(units, 20000, 64));
    assert(bf_context.zeta_bf_proof_status() == silex::ProofState::verified);
    assert_proven_and_saturated(bf_context, units);

    sflint::Arb zeta_hR;
    sflint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;
    assert(silex::zeta_class_regulator_product_bf_audit(
            sflint::ArbRef(zeta_hR), sflint::ArbRef(error_bound), cutoff,
            work_precision, order, 20000, 64));
    silex::ClassGroupContext audit_context;
    prepare_context(audit_context);
    assert(CertificationAccess::try_certify_class_unit_with_bf_audit(
            audit_context, units, sflint::ArbConstRef(zeta_hR),
            sflint::ArbConstRef(error_bound), cutoff, 20000, 64,
            work_precision));
    assert(audit_context.zeta_bf_proof_status() ==
           silex::ProofState::verified);
    assert_proven_and_saturated(audit_context, units);
    return 0;
}

// The internal fallback try_analytic_index_bound_with_units_ keeps complete
// ell-local saturation records, like the other unconditional analytic
// index-one routes.  On Q(sqrt(-23)), h_cand = 3, a saturation proof at
// ell = 3 covers every prime dividing h_cand; the L(1, chi) index-one check
// then certifies `proven` and leaves the verified ell = 3 record in place.
int test_analytic_index_bound_keeps_complete_saturation_records() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 1, 3);
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));

    sflint::Fmpz aux_bound;
    sflint::Fmpz ell;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(ell, 3));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);

    sflint::Arb analytic_hR;
    bool hr_unconditional = false;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(analytic_hR), maximal, 128, nullptr, nullptr,
            nullptr, &hr_unconditional));
    assert(hr_unconditional);
    assert(CertificationAccess::try_analytic_index_bound_with_units(
            context, units, sflint::ArbConstRef(analytic_hR),
            sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);
    sflint::Fmpz record_ell;
    silex::ProofState record_status = silex::ProofState::not_checked;
    assert(context.relation_saturation_record(sflint::FmpzRef(record_ell),
                                              record_status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 3));
    assert(record_status == silex::ProofState::verified);
    return 0;
}

// A public generation check at a bound below the required generation bound
// does not establish that the factor base generates the class group, so it
// cannot stand in for the full check when a stored unconditional analytic
// record would promote.  check_factor_base_generation_bound(1) records a
// verified check at bound 1 once an earlier check has failed; the `proven`
// label reached afterwards must rest on a check that covers the required
// bound.  Q(sqrt(10)): h = 2.  A `grh` transaction records a
// Belabas-Friedman check, so the unconditional record (the quadratic
// L(1, chi) route) is stored here through the internal access.
int test_low_bound_generation_check_does_not_promote() {
    silex::NumberField field = quadratic_field(10);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), maximal_order));

    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            class_group, maximal_order,
            sflint::FmpzConstRef(factor_base_bound), options, 128));
    assert(class_group.certification_status() ==
           silex::CertificationMode::grh);
    assert(class_group.analytic_class_regulator_certification() ==
           silex::CertificationMode::grh);
    CertificationAccess::record_grh_acceptance_analytic_check(class_group,
                                                              true);
    assert(class_group.analytic_class_regulator_certification() ==
           silex::CertificationMode::proven);
    assert(class_group.certification_status() ==
           silex::CertificationMode::grh);

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    silex::ClassGroupComputeOptions proven_options = options;
    proven_options.requested_certification =
            silex::CertificationMode::proven;
    assert(proven_units.compute_with_class_group(
            proven_class_group, maximal_order,
            sflint::FmpzConstRef(factor_base_bound), proven_options, 128));
    assert(proven_units.certification_status() ==
           silex::CertificationMode::proven);

    sflint::Fmpz generation_bound;
    assert(class_group.factor_base_generation_bound(
            sflint::FmpzRef(generation_bound)));
    assert(sflint::fmpz_cmp_ui(sflint::FmpzConstRef(generation_bound), 1) >
           0);
    sflint::Fmpz build_bound;
    assert(class_group.factor_base_build_bound(sflint::FmpzRef(build_bound)));

    // A check beyond the build bound fails and clears the verified check;
    // a check at bound 1 then verifies trivially.
    sflint::Fmpz beyond_build;
    sflint::fmpz_add_ui(sflint::FmpzRef(beyond_build),
                        sflint::FmpzConstRef(build_bound), 1);
    assert(!class_group.check_factor_base_generation_bound(
            sflint::FmpzConstRef(beyond_build)));
    sflint::Fmpz one;
    sflint::fmpz_one(sflint::FmpzRef(one));
    assert(class_group.check_factor_base_generation_bound(
            sflint::FmpzConstRef(one)));
    sflint::Fmpz checked_bound;
    assert(class_group.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));
    assert(sflint::fmpz_is_one(checked_bound));

    // Precision 1 disables the L(1, chi) recompute (the zeta product is
    // unavailable there), so the stored record is the only analytic route,
    // and with a check at bound 1 it must not promote.  The class group
    // still becomes `proven` through the saturation route, which re-runs
    // the generation check at the required bound and leaves its ell = 2
    // record (h = 2); the stored record route leaves none.
    // try_certify_with_units offers no way to switch the saturation route
    // off, so the stored record route is seen only through these side
    // effects.
    {
        sflint::Arb precision_one_hR;
        assert(!silex::zeta_class_regulator_product(
                sflint::ArbRef(precision_one_hR), maximal_order, 1));
    }
    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    assert(class_group.relation_saturation_record_count() == 0);
    assert(class_group.try_certify_with_units(
            proven_units, silex::CertificationMode::proven, 1));
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.factor_base_generation_checked_bound(
            sflint::FmpzRef(checked_bound)));
    assert(sflint::fmpz_cmp(sflint::FmpzConstRef(checked_bound),
                            sflint::FmpzConstRef(generation_bound)) >= 0);
    assert(class_group.relation_saturation_record_count() == 1);
    sflint::Fmpz record_ell;
    silex::ProofState record_status = silex::ProofState::not_checked;
    assert(class_group.relation_saturation_record(
            sflint::FmpzRef(record_ell), record_status, 0));
    assert(sflint::fmpz_equal_si(record_ell, 2));
    assert(record_status == silex::ProofState::verified);
    return 0;
}

// A `grh` request takes the Belabas-Friedman residue in degree two as in
// every other degree; a `proven` request keeps the unconditional quadratic
// L(1, chi) route first.  On x^2 - 10007 (D = 40028, maximal order) the
// default route gives the unconditional value, and accept_grh_residue gives
// a Belabas-Friedman value with its audit data that overlaps it.  The cache
// recomputes a value taken that way for a request without the flag, and an
// unconditional value serves both.
int test_grh_request_takes_belabas_friedman_residue() {
    silex::NumberField field = quadratic_field(10007);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));

    sflint::Arb unconditional_hR;
    bool unconditional = false;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(unconditional_hR), maximal_order, 128, nullptr,
            nullptr, nullptr, &unconditional));
    assert(unconditional);

    sflint::Arb grh_hR;
    silex::detail::ZetaBfRouteAudit audit;
    unconditional = true;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(grh_hR), maximal_order, 128, nullptr, nullptr,
            nullptr, &unconditional, &audit, true));
    assert(!unconditional);
    assert(audit.cutoff > 0 && audit.requested_precision == 128);
    assert(::arb_overlaps(grh_hR.raw(), unconditional_hR.raw()) != 0);

    silex::detail::AnalyticClassRegulatorCache cache;
    assert(cache.ensure(maximal_order, 128, nullptr, nullptr, false, true));
    assert(!cache.value_unconditional());
    assert(cache.value_bf_audit() != nullptr);
    // A request without the flag does not reuse the GRH value.
    assert(cache.ensure(maximal_order, 128));
    assert(cache.value_unconditional());
    assert(cache.value_bf_audit() == nullptr);
    // The unconditional value now serves a grh request as well.
    assert(cache.ensure(maximal_order, 128, nullptr, nullptr, false, true));
    assert(cache.value_unconditional());
    return 0;
}

// The unconditional analytic index-one routes that do not keep a complete
// saturation proof clear the ell-local proof records together with the
// per-prime records, so no proof record is left beside the `verified`
// status.  On Q(sqrt(-23)), h_cand = 3, a saturation proof at ell = 2
// (ell | w = 2) leaves a verified ell = 2 record and proof record, which do
// not cover ell = 3 | h_cand; each route then certifies `proven` and clears
// both.
int test_analytic_index_one_routes_clear_saturation_proof_records() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::Fmpz aux_bound;
    sflint::Fmpz two;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(two, 2));
    auto has_proof_record = [&](const silex::ClassGroupContext& context)
            noexcept {
        silex::ProofState status = silex::ProofState::not_checked;
        slong rank = -1;
        slong target_rank = -1;
        slong local_primes = -1;
        return CertificationAccess::relation_saturation_proof_record(
                context, sflint::FmpzConstRef(two), status, rank,
                target_rank, local_primes);
    };
    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        prepare_minus_23_candidate(context, maximal, 1, 3);
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 128));
        assert(context.try_prove_relation_saturation_with_units(
                units, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(aux_bound)));
        assert(context.relation_saturation_record_count() == 1);
        assert(context.relation_saturation_status() !=
               silex::ProofState::verified);
        assert(has_proof_record(context));
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
    };
    auto assert_cleared = [&](const silex::ClassGroupContext& context)
            noexcept {
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
        assert(context.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(context.relation_saturation_record_count() == 0);
        assert(!has_proof_record(context));
    };

    // try_certify_with_units(proven): the L(1, chi) index-one check of
    // try_certify_analytic_class_regulator_.
    {
        silex::ClassGroupContext context;
        prepare_context(context);
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::proven, 128));
        assert_cleared(context);
    }

    // try_certify_with_units(proven): promotion from an unconditional
    // analytic check recorded earlier.
    {
        silex::ClassGroupContext context;
        prepare_context(context);
        CertificationAccess::record_grh_acceptance_analytic_check(context,
                                                                  true);
        assert(context.analytic_class_regulator_certification() ==
               silex::CertificationMode::proven);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
        assert(context.try_certify_with_units(
                units, silex::CertificationMode::proven, 128));
        assert_cleared(context);
    }

    // try_certify_class_unit_with_zeta.
    {
        silex::ClassGroupContext context;
        prepare_context(context);
        assert(context.try_certify_class_unit_with_zeta(units, 128));
        assert_cleared(context);
    }

    // The internal fallback try_analytic_index_bound_with_units_.
    {
        silex::ClassGroupContext context;
        prepare_context(context);
        sflint::Arb analytic_hR;
        bool hr_unconditional = false;
        assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
                sflint::ArbRef(analytic_hR), maximal, 128, nullptr, nullptr,
                nullptr, &hr_unconditional));
        assert(hr_unconditional);
        assert(CertificationAccess::try_analytic_index_bound_with_units(
                context, units, sflint::ArbConstRef(analytic_hR),
                sflint::FmpzConstRef(aux_bound), 128, hr_unconditional));
        assert_cleared(context);
    }
    return 0;
}

// The exact imaginary-quadratic route stores a proof at each ell | h only
// where none is verified yet: an ell-local proof already verified at
// ell | h is kept with its own data, and no second record is added.  On
// Q(sqrt(-23)), h = 3, the ell = 3 proof runs the ell-local test (target
// rank one), and try_certify_quadratic then certifies `proven` through the
// exact class number.
int test_exact_route_keeps_verified_ell_local_proof() {
    silex::NumberField field = quadratic_field(-23);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal));

    silex::ClassGroupContext context;
    prepare_minus_23_candidate(context, maximal, 1, 3);
    assert(context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 128));

    sflint::Fmpz aux_bound;
    sflint::Fmpz ell;
    assert(set_fmpz_si(aux_bound, 200));
    assert(set_fmpz_si(ell, 3));
    assert(context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.certification_status() ==
           silex::CertificationMode::unknown);

    silex::ProofState status = silex::ProofState::not_checked;
    slong rank = -1;
    slong target_rank = -1;
    slong local_primes = -1;
    assert(CertificationAccess::relation_saturation_proof_record(
            context, sflint::FmpzConstRef(ell), status, rank, target_rank,
            local_primes));
    assert(status == silex::ProofState::verified);
    assert(target_rank > 0);
    assert(local_primes > 0);

    assert(context.try_certify_quadratic(silex::CertificationMode::proven));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(context.relation_saturation_record_count() == 1);

    silex::ProofState status_after = silex::ProofState::not_checked;
    slong rank_after = -1;
    slong target_rank_after = -1;
    slong local_primes_after = -1;
    assert(CertificationAccess::relation_saturation_proof_record(
            context, sflint::FmpzConstRef(ell), status_after, rank_after,
            target_rank_after, local_primes_after));
    assert(status_after == silex::ProofState::verified);
    assert(rank_after == rank);
    assert(target_rank_after == target_rank);
    assert(local_primes_after == local_primes);
    return 0;
}

int test_compute_candidate_preserves_on_failure() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    sflint::Fmpz order_before;
    assert(context.order(sflint::FmpzRef(order_before)));
    assert(sflint::fmpz_equal_si(order_before, 2));

    silex::ClassGroupCandidateOptions bad_options;
    bad_options.max_candidates = 0;
    assert(!context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                      bad_options));
    assert(context.has_presentation());
    assert(context.generator_count() == 1);
    assert(context.relation_count() == 1);

    sflint::Fmpz order_after;
    assert(context.order(sflint::FmpzRef(order_after)));
    assert(sflint::fmpz_equal_si(order_after, 2));

    bad_options.max_candidates = -1;
    assert(!context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                      bad_options));
    assert(context.order(sflint::FmpzRef(order_after)));
    assert(sflint::fmpz_equal_si(order_after, 2));

    bad_options.max_candidates = WORD_MAX;
    bad_options.max_relations = 0;
    assert(!context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                      bad_options));
    assert(context.order(sflint::FmpzRef(order_after)));
    assert(sflint::fmpz_equal_si(order_after, 2));

    bad_options.max_relations = -1;
    assert(!context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                      bad_options));
    assert(context.order(sflint::FmpzRef(order_after)));
    assert(sflint::fmpz_equal_si(order_after, 2));

    return 0;
}

int test_move_and_swap() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext left(order);
    silex::ClassGroupContext right(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(left.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = left.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(left.append_relation(relation));
    assert(left.publish_presentation());
    assert(!right.has_factor_base());

    swap(left, right);
    assert(!left.has_factor_base());
    assert(right.has_presentation());

    silex::ClassGroupContext moved(std::move(right));
    assert(moved.has_presentation());
    assert(!right.has_presentation());

    sflint::Fmpz order_out;
    assert(moved.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    const silex::FactorBase* moved_base = moved.factor_base();
    assert(moved_base != nullptr);

    sflint::FmpzMat generator_rows(1, 1);
    assert(moved.invariant_generator_matrix(
            sflint::FmpzMatRef(generator_rows)));

    silex::FractionalIdeal class_generator(order);
    sflint::FmpzMat factor_row(1, 1);
    sflint::FmpzMat coords(1, 1);
    assert(moved.invariant_generator(class_generator, 0));
    assert(silex::ideal_factor_over_base(sflint::FmpzMatRef(factor_row),
                                         class_generator, *moved_base));
    assert(sflint::fmpz_equal(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(factor_row), 0, 0),
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(generator_rows), 0, 0)));
    assert(moved.ideal_class_coordinates(sflint::FmpzMatRef(coords),
                                         class_generator));
    assert(mat_entry_is_si(coords, 0, 0, 1));

    silex::FactoredElement power_witness(field);
    silex::Element witness_value(field);
    silex::FractionalIdeal witness_ideal(order);
    silex::FractionalIdeal expected_power(order);
    assert(moved.invariant_generator_power_witness(power_witness, 0));
    assert(power_witness.evaluate(witness_value));
    assert(witness_ideal.set_principal(witness_value));
    assert(expected_power.pow_fmpz(class_generator,
                                  sflint::FmpzConstRef(order_out)));
    assert(witness_ideal.equal(expected_power));

    assert(right.define(order));
    assert(right.is_defined());
    assert(!right.has_factor_base());

    silex::ClassGroupContext assigned;
    assigned = std::move(moved);
    assert(assigned.has_presentation());
    assert(!moved.has_presentation());
    assert(moved.define(order));
    assert(moved.is_defined());
    assert(!moved.has_factor_base());

    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.define(order));
    assert(assigned.is_defined());
    assert(!assigned.has_factor_base());

    return 0;
}

}  // namespace

int test_grh_provenance_records() {
    // GRH provenance of a grh class/unit result, and the plumbing of the
    // GRH generation record.  x^3 + x + 200 (|D| = 1080004): the grh
    // policy's GRH generation bound is the BDF bound 19 (the smallest
    // integer cutoff where the criterion holds), its base is built to 76
    // (19 doubled until the base has more than 20 prime ideals), and the
    // Minkowski-type bound is 295, so generation rests on GRH; completeness
    // rests on a Belabas-Friedman hR.
    silex::NumberField field = cubic_field(1, 200);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz minkowski;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(minkowski),
                                               maximal_order));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext context;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            context, maximal_order, sflint::FmpzConstRef(minkowski), options,
            128));
    assert(context.certification_status() == silex::CertificationMode::grh);
    assert(context.factor_base_generation_status() !=
           silex::ProofState::verified);

    sflint::Fmpz value;
    assert(context.factor_base_generation_certification() ==
           silex::CertificationMode::grh);
    assert(context.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::bdf);
    assert(context.factor_base_generation_certification_bound(
            sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 19));
    assert(context.class_unit_completeness_certification() ==
           silex::CertificationMode::grh);
    assert(context.class_unit_completeness_basis() ==
           silex::ClassUnitCompletenessBasis::belabas_friedman);

    // Base completeness, checked against prime decomposition directly:
    // every prime ideal of norm at most the build bound 76 is in the base,
    // including both primes of a split pair, inert primes and ramified
    // primes.  Under GRH these primes generate the class group because 76
    // is at least the GRH generation bound 19.
    const slong grh_build_bound = 76;
    assert(context.factor_base_build_bound(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, grh_build_bound));
    const silex::FactorBase* grh_base = context.factor_base();
    assert(grh_base != nullptr);
    slong checked_primes = 0;
    slong checked_higher_degree = 0;
    slong checked_ramified = 0;
    for (ulong p = 2; p <= static_cast<ulong>(grh_build_bound);
         p = n_nextprime(p, 1)) {
        sflint::Fmpz rational_prime;
        sflint::fmpz_set_ui(sflint::FmpzRef(rational_prime), p);
        silex::PrimeIdealList above;
        assert(silex::decompose_prime(above, maximal_order,
                                      sflint::FmpzConstRef(rational_prime)));
        assert(above.size() > 0);
        for (slong i = 0; i < above.size(); ++i) {
            const silex::PrimeIdeal* prime = above.at(i);
            assert(prime != nullptr);
            sflint::Fmpz norm;
            assert(prime->norm(sflint::FmpzRef(norm)));
            if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(norm),
                                    static_cast<ulong>(grh_build_bound)) >
                0) {
                continue;
            }
            assert(grh_base->contains(*prime));
            ++checked_primes;
            checked_higher_degree += prime->residue_degree() > 1 ? 1 : 0;
            checked_ramified += prime->ramification_index() > 1 ? 1 : 0;
        }
    }
    assert(checked_primes == grh_base->length());
    assert(checked_higher_degree > 0);
    assert(checked_ramified > 0);

    // Only the two GRH theorems are recordable, and the GRH bound must not
    // exceed the build bound (76 here).
    sflint::Fmpz grh_bound;
    assert(set_fmpz_si(grh_bound, 77));
    assert(!CertificationAccess::record_grh_factor_base_generation(
            context, sflint::FmpzConstRef(grh_bound),
            silex::ClassGroupGenerationBasis::bach));
    assert(set_fmpz_si(grh_bound, 15));
    assert(!CertificationAccess::record_grh_factor_base_generation(
            context, sflint::FmpzConstRef(grh_bound),
            silex::ClassGroupGenerationBasis::minkowski_type));
    assert(!CertificationAccess::record_grh_factor_base_generation(
            context, sflint::FmpzConstRef(grh_bound),
            silex::ClassGroupGenerationBasis::none));
    assert(context.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::bdf);
    // A Bach record (the theorem chosen when Bach's bound is the smaller)
    // is reported as such, with its bound.
    assert(CertificationAccess::record_grh_factor_base_generation(
            context, sflint::FmpzConstRef(grh_bound),
            silex::ClassGroupGenerationBasis::bach));
    assert(context.factor_base_generation_certification() ==
           silex::CertificationMode::grh);
    assert(context.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::bach);
    assert(context.factor_base_generation_certification_bound(
            sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 15));
    assert(context.factor_base_generation_status() !=
           silex::ProofState::verified);

    // The proven partner: unconditional generation and completeness.
    // Generation rests on Zimmert's twin-class bound
    // ceil(sqrt(1080004) / 3.9519...) = 263 (Korollar 1, certified enclosure
    // of F), below the Minkowski-type bound 295, with no hypothesis on the class of
    // the different (GB-kor1-generation).
    options.requested_certification = silex::CertificationMode::proven;
    silex::ClassGroupContext proven;
    silex::OrderUnitGroup proven_units;
    assert(proven_units.compute_with_class_group(
            proven, maximal_order, sflint::FmpzConstRef(minkowski), options,
            128));
    assert(proven.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven.factor_base_generation_certification() ==
           silex::CertificationMode::proven);
    assert(proven.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::zimmert_twin_class);
    assert(CertificationAccess::factor_base_generation_different_route(
                   proven) ==
           silex::detail::GenerationDifferentRoute::kor1_generation);
    assert(sflint::fmpz_equal_si(minkowski, 295));
    assert(proven.factor_base_generation_certification_bound(
            sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 263));
    assert(proven.class_unit_completeness_certification() ==
           silex::CertificationMode::proven);
    assert(proven.class_unit_completeness_basis() ==
           silex::ClassUnitCompletenessBasis::unconditional_certification);

    // A context without a presentation reports nothing.
    silex::ClassGroupContext empty;
    assert(empty.factor_base_generation_certification() ==
           silex::CertificationMode::unknown);
    assert(empty.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::none);
    assert(!empty.factor_base_generation_certification_bound(
            sflint::FmpzRef(value)));
    assert(empty.class_unit_completeness_certification() ==
           silex::CertificationMode::unknown);
    assert(empty.class_unit_completeness_basis() ==
           silex::ClassUnitCompletenessBasis::none);
    return 0;
}

struct GenerationRequirementLogCounts {
    // Factor-base builds that recorded the twin-class requirement.
    slong recorded = 0;
    // Builds that took it from the class/unit transaction instead of
    // evaluating it; recorded - kept counts the evaluations.
    slong kept = 0;
    // be_honest checks that used the recorded requirement.
    slong reused = 0;
};

void count_generation_requirement_logs(void* user,
                                       silex::DiagnosticsModule,
                                       silex::LogLevel,
                                       const char*,
                                       const char* message,
                                       const char*) noexcept {
    auto* counts = static_cast<GenerationRequirementLogCounts*>(user);
    if (message == nullptr) {
        return;
    }
    if (std::strcmp(message,
                    "generation requirement zimmert_twin_class, "
                    "different kor1_generation") == 0) {
        ++counts->recorded;
    } else if (std::strcmp(message,
                           "reusing the transaction generation "
                           "requirement") == 0) {
        ++counts->kept;
    } else if (std::strcmp(message,
                           "be_honest using the recorded generation "
                           "requirement") == 0) {
        ++counts->reused;
    }
}

void attach_generation_requirement_counter(
        silex::DiagnosticsContext& diagnostics,
        GenerationRequirementLogCounts& counts) {
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_logging(
            diagnostics, silex::LogLevel::detail,
            silex::diagnostics_module_bit(
                    silex::DiagnosticsModule::class_group),
            count_generation_requirement_logs, &counts);
}

// A proven class/unit run evaluates the twin-class generation requirement
// (Zimmert's G) once; the be_honest check reuses the requirement recorded
// with the factor base instead of evaluating G again.
int test_honesty_reuses_recorded_generation_requirement() {
    silex::NumberField field = cubic_field(1, 200);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz minkowski;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(minkowski),
                                               maximal_order));

    GenerationRequirementLogCounts counts;
    silex::DiagnosticsContext diagnostics;
    attach_generation_requirement_counter(diagnostics, counts);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    options.diagnostics = &diagnostics;
    silex::ClassGroupContext context;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            context, maximal_order, sflint::FmpzConstRef(minkowski), options,
            128));
    assert(context.certification_status() ==
           silex::CertificationMode::proven);
    assert(context.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::zimmert_twin_class);
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    assert(counts.recorded >= 1);
    assert(counts.recorded - counts.kept == 1);
    assert(counts.reused >= 1);
#else
    assert(counts.recorded == 0);
    assert(counts.kept == 0);
    assert(counts.reused == 0);
#endif

    using Access = silex::detail::ClassGroupRelationSearchAccess;
    // The recorded requirement is the one proven_generation_requirement
    // gives for the context's order.
    sflint::Fmpz expected;
    silex::detail::GenerationDifferentRoute route =
            silex::detail::GenerationDifferentRoute::none;
    assert(silex::detail::proven_generation_requirement(
            sflint::FmpzRef(expected), route, maximal_order));
    assert(route == silex::detail::GenerationDifferentRoute::kor1_generation);
    sflint::Fmpz recorded;
    assert(Access::recorded_generation_requirement(
            sflint::FmpzRef(recorded), context, maximal_order));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(recorded),
                              sflint::FmpzConstRef(expected)));

    // Another order, even the maximal order of the same field, is not the
    // context's order: nothing is returned and `out` is left unchanged.
    silex::Order other_order(field);
    assert(other_order.maximal_order(equation_order));
    sflint::Fmpz sentinel;
    sflint::fmpz_set_si(sflint::FmpzRef(sentinel), -7);
    sflint::Fmpz untouched;
    sflint::fmpz_set_si(sflint::FmpzRef(untouched), -7);
    assert(!Access::recorded_generation_requirement(
            sflint::FmpzRef(untouched), context, other_order));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(untouched),
                              sflint::FmpzConstRef(sentinel)));

    // Without a factor base nothing is recorded.
    silex::ClassGroupContext empty(maximal_order);
    assert(!Access::recorded_generation_requirement(
            sflint::FmpzRef(untouched), empty, maximal_order));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(untouched),
                              sflint::FmpzConstRef(sentinel)));

    // set_factor_base installs a base without recording a requirement.
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::ClassGroupContext installed(maximal_order);
    assert(installed.set_factor_base(*base));
    assert(installed.has_factor_base());
    assert(!Access::recorded_generation_requirement(
            sflint::FmpzRef(untouched), installed, maximal_order));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(untouched),
                              sflint::FmpzConstRef(sentinel)));
    return 0;
}

// Within one class/unit transaction a factor-base rebuild of the same order
// reuses the requirement of the first build instead of evaluating G again,
// while a build of another order (here another maximal order of the same
// field, with different data) evaluates it and replaces the kept value;
// without a transaction each build evaluates it.
int test_transaction_keeps_generation_requirement() {
    silex::NumberField field = cubic_field(1, 200);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz expected;
    silex::detail::GenerationDifferentRoute route =
            silex::detail::GenerationDifferentRoute::none;
    assert(silex::detail::proven_generation_requirement(
            sflint::FmpzRef(expected), route, maximal_order));
    assert(route == silex::detail::GenerationDifferentRoute::kor1_generation);
    sflint::Fmpz small_bound;
    sflint::Fmpz large_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(small_bound), 7);
    sflint::fmpz_set_si(sflint::FmpzRef(large_bound), 11);

    using Access = silex::detail::ClassGroupRelationSearchAccess;
    for (const bool attach_transaction : {false, true}) {
        GenerationRequirementLogCounts counts;
        silex::DiagnosticsContext diagnostics;
        attach_generation_requirement_counter(diagnostics, counts);
        silex::detail::ClassUnitTransactionReport audit;
        audit.reset();
        silex::detail::ClassUnitTransactionContext run_context{
                audit, nullptr, {}, {}, false};
        silex::ClassGroupContext context(maximal_order);
        context.set_diagnostics(&diagnostics);
        if (attach_transaction) {
            silex::detail::ClassUnitTransactionAccess::set_run_context(
                    context, &run_context);
        }
        assert(context.build_factor_base(sflint::FmpzConstRef(small_bound)));
        assert(context.build_factor_base(sflint::FmpzConstRef(large_bound)));
        sflint::Fmpz recorded;
        assert(Access::recorded_generation_requirement(
                sflint::FmpzRef(recorded), context, maximal_order));
        assert(sflint::fmpz_equal(sflint::FmpzConstRef(recorded),
                                  sflint::FmpzConstRef(expected)));
        if (attach_transaction) {
            assert(run_context.generation_requirement.order.has_same_data(
                    maximal_order));
        }

        // A build of another order does not reuse the kept requirement: it
        // evaluates G, gets the same value, and becomes the kept order.
        silex::Order other_order(field);
        assert(other_order.maximal_order(equation_order));
        assert(!other_order.has_same_data(maximal_order));
        silex::ClassGroupContext other_context(other_order);
        other_context.set_diagnostics(&diagnostics);
        if (attach_transaction) {
            silex::detail::ClassUnitTransactionAccess::set_run_context(
                    other_context, &run_context);
        }
        assert(other_context.build_factor_base(
                sflint::FmpzConstRef(small_bound)));
        sflint::Fmpz other_recorded;
        assert(Access::recorded_generation_requirement(
                sflint::FmpzRef(other_recorded), other_context, other_order));
        assert(sflint::fmpz_equal(sflint::FmpzConstRef(other_recorded),
                                  sflint::FmpzConstRef(expected)));
        if (attach_transaction) {
            assert(run_context.generation_requirement.order.has_same_data(
                    other_order));
            assert(!run_context.generation_requirement.order.has_same_data(
                    maximal_order));
            silex::detail::ClassUnitTransactionAccess::set_run_context(
                    other_context, nullptr);
            silex::detail::ClassUnitTransactionAccess::set_run_context(
                    context, nullptr);
        }
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
        assert(counts.recorded == 3);
        assert(counts.kept == (attach_transaction ? 1 : 0));
#else
        assert(counts.recorded == 0);
        assert(counts.kept == 0);
#endif
    }
    return 0;
}

// The GRH generation bound falls back to Bach's bound alone when the
// Belabas-Diaz y Diaz-Friedman criterion cannot be computed or decided, and
// the theorem that gave the bound reaches the grh record through the policy.
int test_grh_generation_bound_bach_fallback() {
    // 4x^3 + x + 100 defines the field of x^3 + x + 200 (its root is half a
    // root of x^3 + x + 200), so |D| = 1080004 for both presentations.  The
    // residue-degree route of the BDF criterion needs a monic integral
    // defining polynomial, so for the non-monic presentation BDF cannot be
    // computed and the bound is Bach's alone: ceil(12 log^2 1080004) = 2317
    // (Bach 1990, Theorem 4).  The monic presentation keeps the BDF bound 19.
    silex::NumberField monic_field = cubic_field(1, 200);
    silex::Order monic_equation = silex::test::equation_order(monic_field);
    silex::Order monic_maximal(monic_field);
    assert(monic_maximal.maximal_order(monic_equation));
    sflint::Fmpz bound;
    bool bach_selected = true;
    assert(silex::detail::grh_factor_base_bound_with_diagnostics(
            sflint::FmpzRef(bound), monic_maximal, nullptr, &bach_selected));
    assert(sflint::fmpz_equal_si(bound, 19));
    assert(!bach_selected);

    sflint::FmpqPoly nonmonic_polynomial;
    sflint::fmpq_poly_zero(nonmonic_polynomial);
    sflint::fmpq_poly_set_coeff_si(nonmonic_polynomial, 3, 4);
    sflint::fmpq_poly_set_coeff_si(nonmonic_polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(nonmonic_polynomial, 0, 100);
    silex::NumberField nonmonic_field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(nonmonic_polynomial));
    // Z[2 beta] = Z[alpha], with basis 1, 2 beta, 4 beta^2.
    sflint::FmpqMat alpha_basis(3, 3);
    fmpq_set_si(fmpq_mat_entry(alpha_basis.raw(), 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(alpha_basis.raw(), 1, 1), 2, 1);
    fmpq_set_si(fmpq_mat_entry(alpha_basis.raw(), 2, 2), 4, 1);
    silex::Order alpha_order = silex::Order::from_basis(
            nonmonic_field, sflint::FmpqMatConstRef(alpha_basis));
    assert(alpha_order.is_defined());
    silex::Order nonmonic_maximal(nonmonic_field);
    assert(nonmonic_maximal.maximal_order(alpha_order));
    sflint::Fmpz discriminant;
    assert(nonmonic_maximal.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, -1080004));
    bach_selected = false;
    assert(silex::detail::grh_factor_base_bound_with_diagnostics(
            sflint::FmpzRef(bound), nonmonic_maximal, nullptr,
            &bach_selected));
    assert(sflint::fmpz_equal_si(bound, 2317));
    assert(bach_selected);

    // A grh request on the non-monic presentation still fails closed, since
    // prime decomposition (and so the norm-bounded factor base) also needs a
    // monic integral defining polynomial.  The bound above now succeeds, so
    // the request still fails at policy selection, in the factor-base probe.
    // The failure publishes nothing.
    silex::FactorBase failed_probe_base;
    assert(!silex::detail::grh_factor_base_bound(bound, nonmonic_maximal,
                                                 nullptr, nullptr, nullptr,
                                                 &failed_probe_base));
    assert(!failed_probe_base.is_defined());
    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;
    {
        sflint::Fmpz requested_bound;
        sflint::fmpz_set_ui(sflint::FmpzRef(requested_bound), 500);
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        silex::detail::ClassUnitTransactionReport audit;
        assert(!silex::detail::compute_class_unit_transaction(
                units, class_group, nonmonic_maximal,
                sflint::FmpzConstRef(requested_bound), options, 128, audit));
        assert(audit.failure_stage ==
               silex::detail::ClassUnitStage::factor_base_bound);
        assert(!audit.policy.selected);
        assert(!audit.final_result_published);
        assert(audit.class_group_certification ==
               silex::CertificationMode::unknown);
        assert(audit.unit_group_certification ==
               silex::CertificationMode::unknown);
        assert(!class_group.has_factor_base());
        assert(!class_group.has_presentation());
        assert(!units.is_set());
    }

    // End to end through the policy: x^2 + x + 4 (D = -15, class group
    // Z/2).  The BDF bound is 5, the smallest integer cutoff where the
    // criterion holds, below Bach's bound ceil(6 log^2 15) = 45 (Bach 1990,
    // Theorem 4 and the remark following it), so the grh policy records
    // basis bdf with bound 5, and the class group's GRH generation record
    // carries it.  (Before the BDF search reached integer resolution it never
    // returned less than 50 and this field took Bach's bound; with the exact
    // cutoff no field small enough to test has Bach's bound below the BDF
    // bound, so the bach record is checked directly in
    // test_grh_provenance_records.)  The base is built to 80, the bound 5
    // doubled until the base has more than 20 prime ideals.  The base also
    // reaches the Minkowski-type bound of so small a field,
    // so the public generation basis is the unconditional one.
    silex::NumberField quadratic = shifted_quadratic_field(1, 4);
    silex::Order quadratic_equation = silex::test::equation_order(quadratic);
    silex::Order quadratic_maximal(quadratic);
    assert(quadratic_maximal.maximal_order(quadratic_equation));
    sflint::Fmpz minkowski;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(minkowski),
                                               quadratic_maximal));
    silex::ClassGroupGenerationBasis basis =
            silex::ClassGroupGenerationBasis::none;
    sflint::Fmpz generation_bound;
    assert(silex::detail::grh_factor_base_bound(
            bound, quadratic_maximal, nullptr, &generation_bound, &basis));
    assert(basis == silex::ClassGroupGenerationBasis::bdf);
    assert(sflint::fmpz_equal_si(generation_bound, 5));
    assert(sflint::fmpz_equal_si(bound, 80));
    // The probe base handed out with the bound is the norm-bounded base at
    // that bound, the one the transaction's strict search build takes over.
    silex::FactorBase probe_base;
    sflint::Fmpz probe_bound;
    assert(silex::detail::grh_factor_base_bound(probe_bound,
                                                quadratic_maximal, nullptr,
                                                nullptr, nullptr,
                                                &probe_base));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(probe_bound),
                              sflint::FmpzConstRef(bound)));
    silex::FactorBase rebuilt_base(quadratic_maximal);
    assert(rebuilt_base.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(probe_base.length() > 20);
    assert(probe_base.equal(rebuilt_base));

    // The transaction takes the probe base over only for a matching order
    // and bound; every other case yields nullptr (a fresh build).
    {
        silex::detail::ClassUnitTransactionReport takeover_audit;
        silex::detail::ClassUnitTransactionContext takeover{takeover_audit, nullptr, {}, {}, false};
        sflint::Fmpz other_bound;
        sflint::fmpz_set_ui(sflint::FmpzRef(other_bound), 81);
        silex::Order other_order(quadratic);
        assert(other_order.maximal_order(quadratic_equation));
        assert(silex::detail::grh_norm_bounded_base_for(
                       nullptr, quadratic_maximal,
                       sflint::FmpzConstRef(bound)) == nullptr);
        // Nothing stored yet.
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, quadratic_maximal,
                       sflint::FmpzConstRef(bound)) == nullptr);
        takeover.grh_norm_bounded_base = &probe_base;
        takeover.grh_norm_bounded_base_bound = &probe_bound;
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, quadratic_maximal,
                       sflint::FmpzConstRef(bound)) == &probe_base);
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, quadratic_maximal,
                       sflint::FmpzConstRef(other_bound)) == nullptr);
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, other_order,
                       sflint::FmpzConstRef(bound)) == nullptr);
        // A taken (cleared) base is not handed out again.
        silex::FactorBase taken;
        takeover.grh_norm_bounded_base = &taken;
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, quadratic_maximal,
                       sflint::FmpzConstRef(bound)) == nullptr);
        takeover.grh_norm_bounded_base = &probe_base;
        takeover.grh_norm_bounded_base_bound = nullptr;
        assert(silex::detail::grh_norm_bounded_base_for(
                       &takeover, quadratic_maximal,
                       sflint::FmpzConstRef(bound)) == nullptr);
    }

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport audit;
    assert(silex::detail::compute_class_unit_transaction(
            units, class_group, quadratic_maximal,
            sflint::FmpzConstRef(minkowski), options, 128, audit));
    assert(audit.final_result_published);
    assert(audit.policy.factor_base_bound ==
           silex::detail::NativeFactorBaseBoundStrategy::grh);
    assert(audit.policy.grh_generation_basis ==
           silex::ClassGroupGenerationBasis::bdf);
    assert(sflint::fmpz_equal_si(audit.policy.grh_generation_bound, 5));
    assert(class_group.certification_status() ==
           silex::CertificationMode::grh);
    assert(units.certification_status() == silex::CertificationMode::grh);
    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    assert(CertificationAccess::grh_factor_base_generation_covered(
            class_group));
    assert(CertificationAccess::grh_factor_base_generation_basis(
                   class_group) == silex::ClassGroupGenerationBasis::bdf);
    assert(class_group.factor_base_generation_basis() ==
           silex::ClassGroupGenerationBasis::minkowski_type);
    return 0;
}

// Proven generation on Zimmert's twin-class bound needs no hypothesis on
// the class of the different (GB-kor1-generation).  Two more presentations
// of the field of x^3 + x + 200 (class group Z/2, Minkowski-type bound 295,
// twin-class bound 263): x^3 + 25 x + 25000 (theta = 5 alpha, index 125,
// primes above 5 of norm 5) and x^3 + 361 x + 1371800 (theta = 19 alpha, a
// prime of norm 361 above 19, which formerly kept the Minkowski-type bound)
// both take requirement 263.  The class group is the same in every
// presentation.
int test_twin_class_generation_routes() {
    struct Row {
        slong c1;
        slong c0;
        slong requirement;
        silex::ClassGroupGenerationBasis basis;
        silex::detail::GenerationDifferentRoute route;
    };
    const Row rows[] = {
            {25, 25000, 263,
             silex::ClassGroupGenerationBasis::zimmert_twin_class,
             silex::detail::GenerationDifferentRoute::kor1_generation},
            {361, 1371800, 263,
             silex::ClassGroupGenerationBasis::zimmert_twin_class,
             silex::detail::GenerationDifferentRoute::kor1_generation},
    };
    for (const Row& row : rows) {
        silex::NumberField field = cubic_field(row.c1, row.c0);
        silex::Order equation_order = silex::test::equation_order(field);
        silex::Order maximal_order(field);
        assert(maximal_order.maximal_order(equation_order));
        sflint::Fmpz minkowski;
        assert(silex::factor_base_class_group_bound(
                sflint::FmpzRef(minkowski), maximal_order));
        assert(sflint::fmpz_equal_si(minkowski, 295));

        silex::ClassGroupComputeOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        options.zeta_bf_max_cutoff = 20000;
        options.requested_certification = silex::CertificationMode::proven;
        silex::ClassGroupContext proven;
        silex::OrderUnitGroup units;
        assert(units.compute_with_class_group(
                proven, maximal_order, sflint::FmpzConstRef(minkowski),
                options, 128));
        assert(proven.certification_status() ==
               silex::CertificationMode::proven);
        sflint::Fmpz value;
        assert(proven.order(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, 2));
        assert(proven.factor_base_generation_basis() == row.basis);
        assert(CertificationAccess::factor_base_generation_different_route(
                       proven) == row.route);
        assert(proven.factor_base_generation_bound(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, row.requirement));
        assert(proven.factor_base_generation_certification_bound(
                sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, row.requirement));
    }
    return 0;
}

// A `proven` imaginary quadratic request whose base up to the generation
// bound would exceed the relation module limit keeps its compact base and
// is certified by the structure certificate
// (ClassGroupGenerationBasis::structure_certificate).  The limit is lowered
// by the test hook so two small fields with non-cyclic class groups take
// that route: D = -4027 (Z/3 x Z/3, one ell = 3 combination step over two
// forms) and D = -20155 (Z/4 x Z/4).  The class group equals the one
// certified by the generation bound, and no generation bound is published.
int test_imaginary_quadratic_structure_certificate_route() {
    struct Row {
        slong radicand;
        slong discriminant;
        std::vector<slong> invariants;
    };
    const Row rows[] = {
            {-4027, -4027, {3, 3}},
            {-20155, -20155, {4, 4}},
    };
    for (const Row& row : rows) {
        silex::NumberField field = quadratic_field(row.radicand);
        silex::Order equation_order = silex::test::equation_order(field);
        silex::Order maximal_order(field);
        assert(maximal_order.maximal_order(equation_order));
        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(
                sflint::FmpzRef(bound), maximal_order));
        sflint::Fmpz discriminant;
        assert(set_fmpz_si(discriminant, row.discriminant));
        sflint::Fmpz class_number;
        fmpz_one(class_number.raw());
        for (slong invariant : row.invariants) {
            fmpz_mul_si(class_number.raw(), class_number.raw(), invariant);
        }

        silex::ClassGroupComputeOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        options.requested_certification = silex::CertificationMode::proven;

        silex::ClassGroupContext reference;
        silex::OrderUnitGroup reference_units;
        assert(reference_units.compute_with_class_group(
                reference, maximal_order, sflint::FmpzConstRef(bound),
                options, 128));
        assert(reference.certification_status() ==
               silex::CertificationMode::proven);
        assert(reference.factor_base_generation_basis() !=
               silex::ClassGroupGenerationBasis::structure_certificate);

        silex::detail::relation_search::
                set_quadratic_structure_certificate_column_limit_for_testing(
                        1);
        silex::ClassGroupContext forced;
        silex::OrderUnitGroup forced_units;
        const bool computed = forced_units.compute_with_class_group(
                forced, maximal_order, sflint::FmpzConstRef(bound), options,
                128);
        silex::detail::relation_search::
                set_quadratic_structure_certificate_column_limit_for_testing(
                        0);
        assert(computed);
        assert(forced.certification_status() ==
               silex::CertificationMode::proven);
        assert(forced_units.certification_status() ==
               silex::CertificationMode::proven);
        assert(forced.factor_base_generation_certification() ==
               silex::CertificationMode::proven);
        assert(forced.factor_base_generation_basis() ==
               silex::ClassGroupGenerationBasis::structure_certificate);
        assert(forced.factor_base_generation_status() !=
               silex::ProofState::verified);
        sflint::Fmpz value;
        assert(!forced.factor_base_generation_certification_bound(
                sflint::FmpzRef(value)));
        assert(forced.class_unit_completeness_certification() ==
               silex::CertificationMode::proven);
        assert(forced.class_unit_completeness_basis() ==
               silex::ClassUnitCompletenessBasis::unconditional_certification);
        assert(forced.relation_saturation_status() ==
               silex::ProofState::verified);

        assert(forced.invariant_count() ==
               static_cast<slong>(row.invariants.size()));
        assert(reference.invariant_count() == forced.invariant_count());
        for (slong j = 0; j < forced.invariant_count(); ++j) {
            assert(forced.invariant(sflint::FmpzRef(value), j));
            assert(sflint::fmpz_equal_si(
                    value, row.invariants[static_cast<std::size_t>(j)]));
            assert(reference.invariant(sflint::FmpzRef(value), j));
            assert(sflint::fmpz_equal_si(
                    value, row.invariants[static_cast<std::size_t>(j)]));
        }

        // The certificate itself, without recording: it passes against
        // the exact class number and fails against any other order.
        assert(CertificationAccess::imaginary_quadratic_structure_certificate(
                forced, sflint::FmpzConstRef(discriminant),
                sflint::FmpzConstRef(class_number), false));
        assert(CertificationAccess::imaginary_quadratic_structure_certificate(
                reference, sflint::FmpzConstRef(discriminant),
                sflint::FmpzConstRef(class_number), false));
        sflint::Fmpz doubled;
        fmpz_mul_ui(doubled.raw(), class_number.raw(), 2);
        assert(!CertificationAccess::imaginary_quadratic_structure_certificate(
                forced, sflint::FmpzConstRef(discriminant),
                sflint::FmpzConstRef(doubled), false));
    }
    return 0;
}

// The proven witness step for x^2 + 10000019 reduces binary forms above the
// fast reduction switch (more than 9000 bits).  Every published relation
// generator must factor to exactly its relation row.
int test_imaginary_quadratic_proven_large_witness_forms() {
    silex::NumberField field = shifted_quadratic_field(0, 10000019);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                maximal_order));

    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::proven;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, maximal_order,
                                          sflint::FmpzConstRef(bound),
                                          options, 128));
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_count() > 0);
    assert(context_relation_rows_match_exact_generators(class_group));
    return 0;
}

// A relation lattice of index h does not prove that the factor base
// generates the class group (the Cl(-199) example of the structure
// certificate's design).  D = -4363 has class group Z/9; 2, 3, 5 and 7 are
// inert and 11 = P Pbar splits with P of order 3: P^3 = (alpha),
// alpha = (31 + sqrt(-4363)) / 2 of norm 11^3.  The base up to 11 is
// {(2), (3), (5), (7), P, Pbar} (or its orientation), and the true
// relations (2), (3), (5), (7), (11) and alpha^3 (P^9) have index 9 = h, so
// the presentation is Z/9 generated by P.  Its order is h and g^9 is
// principal, but g^3 = P^3 is principal too, so the ell = 3 step rejects
// it.
int test_imaginary_quadratic_structure_certificate_rejects_subgroup() {
    silex::NumberField field = quadratic_field(-4363);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));

    silex::ClassGroupContext context(maximal_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 11));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    silex::Element theta(field);
    silex::Element shifted(field);
    silex::Element alpha(field);
    silex::Element alpha_cubed(field);
    sflint::Fmpz three;
    assert(set_fmpz_si(three, 3));
    assert(theta.gen() && shifted.add_si(theta, 31) &&
           alpha.scalar_div_si(shifted, 2) &&
           alpha_cubed.pow_fmpz(alpha, sflint::FmpzConstRef(three)));
    for (slong rational : {2, 3, 5, 7, 11}) {
        silex::Element generator(field);
        assert(generator.set_si(rational));
        assert(context.try_append_generator_relation(generator));
    }
    assert(context.try_append_generator_relation(alpha_cubed));
    assert(context.publish_presentation());

    sflint::Fmpz value;
    assert(context.order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 9));
    assert(context.invariant_count() == 1);
    sflint::Fmpz discriminant;
    assert(set_fmpz_si(discriminant, -4363));
    assert(!CertificationAccess::imaginary_quadratic_structure_certificate(
            context, sflint::FmpzConstRef(discriminant),
            sflint::FmpzConstRef(value), false));
    return 0;
}

// A factor base with more generators than the relation modules allow
// (2^14) is refused before their dense mod-p matrices are allocated, the
// same way on every machine: the context keeps no base, and an attached
// class/unit transaction names the reason factor_base_size_limit_exceeded
// until the next attempt clears it or a later base within the limit is
// built.
int test_factor_base_size_limit_fails_closed() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    // pi(200000) = 17984 primes of Q, one prime ideal each.
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 200000));
    silex::FactorBase base(order);
    assert(base.is_defined() &&
           base.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(bound)));
    assert(base.length() > (WORD(1) << 14));

    silex::ClassGroupContext detached(order);
    assert(!detached.set_factor_base(base));
    assert(!detached.has_factor_base());

    silex::detail::ClassUnitTransactionReport audit;
    silex::detail::ClassUnitTransactionContext transaction{
            audit, nullptr, {}, {}, false};
    silex::ClassGroupContext context(order);
    silex::detail::ClassUnitTransactionAccess::set_run_context(
            context, &transaction);
    assert(!context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(!context.has_factor_base());
    assert(audit.factor_base_size_limit_exceeded);
    assert(std::strcmp(
                   silex::detail::class_unit_computation_failure_reason(
                           audit),
                   "factor_base_size_limit_exceeded") == 0);
    silex::detail::ClassGroupRelationSearchAccess::
            clear_factor_base_honesty_failure(context);
    assert(!audit.factor_base_size_limit_exceeded);
    assert(std::strcmp(
                   silex::detail::class_unit_computation_failure_reason(
                           audit),
                   "class_unit_computation_failed") == 0);
    assert(!context.set_factor_base(base));
    assert(audit.factor_base_size_limit_exceeded);
    // A later base within the limit clears the flag, so a failure after it
    // is not reported as the earlier size-limit refusal.
    sflint::Fmpz small_bound;
    assert(set_fmpz_si(small_bound, 100));
    assert(context.build_factor_base(sflint::FmpzConstRef(small_bound)));
    assert(context.has_factor_base());
    assert(!audit.factor_base_size_limit_exceeded);
    assert(std::strcmp(
                   silex::detail::class_unit_computation_failure_reason(
                           audit),
                   "class_unit_computation_failed") == 0);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return 0;
}

// A grh class/unit transaction whose relation-search factor base exceeds
// the relation module size limit publishes the documented reason
// factor_base_size_limit_exceeded.  The relation search of a field of
// degree at least 2 builds its base on a fresh candidate context with no
// transaction attached, so the run's transaction is passed to the size
// check.  Q(sqrt(p)), p = 10^130 + 1771 (a prime, 3 mod 4): the GRH base has
// more than 2^14 generators.
int test_grh_relation_search_size_limit_reason() {
    sflint::Fmpz p;
    sflint::fmpz_set_ui(sflint::FmpzRef(p), 10);
    sflint::fmpz_pow_ui(sflint::FmpzRef(p), sflint::FmpzConstRef(p), 130);
    sflint::fmpz_add_ui(sflint::FmpzRef(p), sflint::FmpzConstRef(p), 1771);
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpz_neg(sflint::FmpzRef(p), sflint::FmpzConstRef(p));
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, p.raw());
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport audit;
    assert(!silex::detail::compute_class_unit_transaction(
            units, class_group, order, sflint::FmpzConstRef(bound), options,
            128, audit));
    assert(audit.policy.selected);
    assert(audit.factor_base_size_limit_exceeded);
    assert(audit.failure_reason != nullptr &&
           std::strcmp(audit.failure_reason,
                       "factor_base_size_limit_exceeded") == 0);
    assert(!audit.final_result_published);
    assert(!class_group.has_presentation());
    assert(!units.is_set());
    return 0;
}

// The two index routes of a rank-zero or quadratic class/unit pair compare
// the candidate with the same normalization as the value they test against.
// The analytic route divides the candidate product h_cand R_cand
// (OrderUnitGroup::class_regulator_product, without w) by the zeta value
// res w sqrt|D| / (2^r1 (2 pi)^r2) = h R (class_regulator_product_from_
// residue, src/zeta/zeta.cpp), so w enters only the analytic side, where
// the class number formula puts it.  The exact imaginary-quadratic route
// divides h_cand by h(D), the number of reduced primitive forms of the
// fundamental discriminant D (FLINT qfb_reduced_forms counts primitive
// forms only), which is the class number of the maximal order.  On proven
// pairs with w = 4, w = 6, h = 3, a real quadratic and a mixed-signature
// cubic field, the candidate product lies in the zeta value, the exact route
// gives index one, the analytic accessor gives its conservative bound 2, and
// the zeta value matches an independent h R.  Q(sqrt(-11)) = Q[x]/(x^2 + 11)
// has polynomial discriminant -44 (h(-44) = 3) but field discriminant -11
// (h = 1), so reading D from the equation order changes the exact index to 3
// and the zeta value by a factor of 2.
int test_class_unit_index_routes_share_normalization() {
    struct Case {
        silex::NumberField field;
        slong class_number;
        slong torsion;
        bool exact_quadratic;
    };
    std::vector<Case> cases;
    cases.push_back({quadratic_field(-1), 1, 4, true});
    cases.push_back({shifted_quadratic_field(-1, 1), 1, 6, true});
    cases.push_back({shifted_quadratic_field(-1, 6), 3, 2, true});
    cases.push_back({quadratic_field(10), 2, 2, false});
    cases.push_back({cubic_field(0, -2), 1, 2, false});
    cases.push_back({quadratic_field(-11), 1, 2, true});

    const slong precision = 256;
    for (Case& item : cases) {
        silex::Order equation = silex::test::equation_order(item.field);
        silex::Order order(item.field);
        assert(order.maximal_order(equation));
        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                   order));
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
        }
        silex::ClassGroupComputeOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        options.requested_certification = silex::CertificationMode::proven;
        silex::ClassGroupContext context;
        silex::OrderUnitGroup units;
        assert(units.compute_with_class_group(
                context, order, sflint::FmpzConstRef(bound), options,
                128));
        assert(context.certification_status() ==
               silex::CertificationMode::proven);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        sflint::Fmpz value;
        assert(context.order(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, item.class_number));
        assert(units.torsion_order(sflint::FmpzRef(value)));
        assert(sflint::fmpz_equal_si(value, item.torsion));

        sflint::Arb candidate;
        sflint::Arb analytic;
        assert(units.class_regulator_product(sflint::ArbRef(candidate),
                                             context, precision));
        assert(silex::zeta_class_regulator_product(
                sflint::ArbRef(analytic), order, precision));
        assert(::arb_overlaps(candidate.raw(), analytic.raw()) != 0);

        // Independent h R: h times log of a fundamental unit, 3 + sqrt(10)
        // for Q(sqrt(10)) and 1 + 2^(1/3) + 2^(2/3) for Q(2^(1/3)); R = 1
        // for the imaginary quadratic fields.
        sflint::Arb reference;
        if (item.exact_quadratic) {
            ::arb_set_si(reference.raw(), item.class_number);
        } else if (order.degree() == 2) {
            ::arb_sqrt_ui(reference.raw(), 10, precision);
            ::arb_add_ui(reference.raw(), reference.raw(), 3, precision);
            ::arb_log(reference.raw(), reference.raw(), precision);
            ::arb_mul_si(reference.raw(), reference.raw(),
                         item.class_number, precision);
        } else {
            sflint::Arb cube_root;
            ::arb_set_ui(cube_root.raw(), 2);
            ::arb_root_ui(cube_root.raw(), cube_root.raw(), 3, precision);
            ::arb_mul(reference.raw(), cube_root.raw(), cube_root.raw(),
                      precision);
            ::arb_add(reference.raw(), reference.raw(), cube_root.raw(),
                      precision);
            ::arb_add_ui(reference.raw(), reference.raw(), 1, precision);
            ::arb_log(reference.raw(), reference.raw(), precision);
            ::arb_mul_si(reference.raw(), reference.raw(),
                         item.class_number, precision);
        }
        assert(::arb_overlaps(analytic.raw(), reference.raw()) != 0);

        sflint::Fmpz index;
        assert(set_fmpz_si(index, 42));
        assert(units.class_regulator_index_bound(
                sflint::FmpzRef(index), context,
                sflint::ArbConstRef(analytic), precision));
        // The public accessor takes the ceiling of the upper endpoint of the
        // quotient enclosure and is conservative on purpose (order_unit.cpp),
        // so a true index one is reported as 2.  This pins that ceiling: a
        // quotient too large by a factor of 2 would give 3, and one too small
        // would fail closed.  The radius-independent check of the
        // normalization is the pair of overlap assertions above.
        assert(sflint::fmpz_equal_si(index, 2));

        if (item.exact_quadratic) {
            assert(set_fmpz_si(index, 42));
            assert(CertificationAccess::rank_zero_quadratic_class_index_bound(
                    sflint::FmpzRef(index), context, units,
                    silex::CertificationMode::proven));
            assert(sflint::fmpz_is_one(sflint::FmpzConstRef(index)));
            assert(context.certification_status() ==
                   silex::CertificationMode::proven);
        }
    }
    return 0;
}

// The exact imaginary-quadratic route has no clamp of its own: an h_cand
// below h is not divisible by h, so the route returns false and publishes
// nothing, and the analytic route fails closed on the same candidate
// (quotient h_cand / h = 1/2 < 1).  Q(sqrt(-14)) has h = 4.  The base up to
// 2 holds only the ramified prime above 2, of order 2 in the class group,
// and the relation (2) = p^2 gives h_cand = 2.  The GRH generation record
// for that base is written directly and is false: it stands for a wrong
// generation assumption, the case in which h_cand < h can occur.  A grh
// validation run then ends without accepting the pair.
int test_exact_quadratic_index_below_one_fails_closed() {
    silex::NumberField field = quadratic_field(-14);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    assert(units.certification_status() == silex::CertificationMode::proven);

    auto prepare = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element two(field);
        assert(two.set_si(2));
        assert(relation.set_generator(two));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        sflint::Fmpz class_order;
        assert(context.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 2));
        assert(CertificationAccess::record_grh_factor_base_generation(
                context, sflint::FmpzConstRef(bound),
                silex::ClassGroupGenerationBasis::bach));
        assert(CertificationAccess::grh_factor_base_generation_covered(
                context));
    };

    // Exact route: no index, nothing published.
    {
        silex::ClassGroupContext context;
        prepare(context);
        sflint::Fmpz index;
        assert(set_fmpz_si(index, 42));
        assert(!CertificationAccess::rank_zero_quadratic_class_index_bound(
                sflint::FmpzRef(index), context, units,
                silex::CertificationMode::grh));
        assert(sflint::fmpz_equal_si(index, 42));
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
        assert(context.class_unit_completeness_certification() ==
               silex::CertificationMode::unknown);

        // Analytic route on the same candidate: the zeta value is h = 4.
        sflint::Arb analytic;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(analytic),
                                                   order, 128));
        assert(::arb_contains_si(analytic.raw(), 4) != 0);
        assert(!units.class_regulator_index_bound(
                sflint::FmpzRef(index), context,
                sflint::ArbConstRef(analytic), 128));
        assert(sflint::fmpz_equal_si(index, 42));
    }

    // grh validation loop: neither route accepts the pair.
    {
        silex::ClassGroupContext context;
        prepare(context);
        silex::OrderUnitGroup loop_units;
        assert(loop_units.compute(order));
        silex::OrderUnitGroup scratch(order);
        silex::EmbeddingContext embeddings(field);
        silex::ClassGroupComputeOptions options;
        options.requested_certification = silex::CertificationMode::grh;
        options.max_relations = 64;
        silex::detail::AnalyticClassRegulatorCache analytic_cache;
        silex::detail::ValidateRefineSummary summary;
        assert(!silex::detail::try_validate_refine_loop(
                context, loop_units, scratch, order, options, embeddings,
                analytic_cache, summary, 128, false));
        assert(summary.outcome ==
               silex::detail::ValidateRefineOutcome::analytic_unavailable);
        assert(context.certification_status() ==
               silex::CertificationMode::unknown);
        assert(context.class_unit_completeness_certification() ==
               silex::CertificationMode::unknown);
        sflint::Fmpz class_order;
        assert(context.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 2));
    }
    return 0;
}

// End to end through the validation loop (try_validate_refine_loop): the
// exact imaginary-quadratic index is an unconditional saturation bound only
// when factor-base generation is checked to the required bound.  On the
// Q(sqrt(-23)) candidate with h_cand = 6, h = 3, a grh run assumes
// generation under GRH, takes the exact index 2, clears ell = 2 and then
// accepts index one; relation saturation ends `unavailable`, below
// `verified`.  The same candidate under a proven request checks generation
// first and ends `verified` and `proven`.  The grh generation record that
// the grh base build writes is written directly here for the hand-built
// base; primes of norm at most 3 do generate this class group (Minkowski
// bound (2 / pi) sqrt(23) < 3.1).
int test_validation_loop_exact_index_status_follows_request() {
    silex::NumberField field = shifted_quadratic_field(-1, 6);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    for (const silex::CertificationMode requested :
         {silex::CertificationMode::grh, silex::CertificationMode::proven}) {
        const bool grh = requested == silex::CertificationMode::grh;
        silex::ClassGroupContext context;
        prepare_minus_23_index_two_candidate(context, order);
        if (grh) {
            sflint::Fmpz grh_bound;
            assert(set_fmpz_si(grh_bound, 3));
            assert(CertificationAccess::record_grh_factor_base_generation(
                    context, sflint::FmpzConstRef(grh_bound),
                    silex::ClassGroupGenerationBasis::bach));
        }
        assert(context.factor_base_generation_checked_status() ==
               silex::ProofState::not_checked);

        silex::OrderUnitGroup units;
        assert(units.compute(order));
        silex::OrderUnitGroup scratch(order);
        silex::EmbeddingContext embeddings(field);
        silex::ClassGroupComputeOptions options;
        options.requested_certification = requested;
        options.max_relations = 64;
        silex::detail::AnalyticClassRegulatorCache analytic_cache;
        silex::detail::ValidateRefineSummary summary;
        assert(silex::detail::try_validate_refine_loop(
                context, units, scratch, order, options, embeddings,
                analytic_cache, summary, 128, false));
        assert(summary.outcome == silex::detail::ValidateRefineOutcome::proven);
        assert(summary.class_progress);
        sflint::Fmpz class_order;
        assert(context.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 3));

        if (grh) {
            // The caller publishes the grh labels; the saturation status
            // stays below `verified`.
            assert(context.certification_status() ==
                   silex::CertificationMode::unknown);
            assert(CertificationAccess::publish_grh_labels(context, units,
                                                           nullptr));
            assert(context.certification_status() ==
                   silex::CertificationMode::grh);
            assert(context.relation_saturation_status() ==
                   silex::ProofState::unavailable);
            assert(context.class_unit_completeness_certification() ==
                   silex::CertificationMode::proven);
            assert(context.class_unit_completeness_basis() ==
                   silex::ClassUnitCompletenessBasis::exact_class_number);
        } else {
            assert(context.certification_status() ==
                   silex::CertificationMode::proven);
            assert(context.relation_saturation_status() ==
                   silex::ProofState::verified);
            assert(units.certification_status() ==
                   silex::CertificationMode::proven);
        }
    }
    return 0;
}

// Asserts that two class/unit results report the same certification and
// provenance.
void assert_same_class_unit_provenance(const silex::ClassGroupContext& left,
                                       const silex::OrderUnitGroup& left_units,
                                       const silex::ClassGroupContext& right,
                                       const silex::OrderUnitGroup& right_units)
        noexcept {
    assert(left.certification_status() == right.certification_status());
    assert(left_units.certification_status() ==
           right_units.certification_status());
    assert(left.factor_base_generation_certification() ==
           right.factor_base_generation_certification());
    assert(left.factor_base_generation_basis() ==
           right.factor_base_generation_basis());
    sflint::Fmpz left_bound;
    sflint::Fmpz right_bound;
    const bool left_has_bound =
            left.factor_base_generation_certification_bound(
                    sflint::FmpzRef(left_bound));
    const bool right_has_bound =
            right.factor_base_generation_certification_bound(
                    sflint::FmpzRef(right_bound));
    assert(left_has_bound == right_has_bound);
    assert(!left_has_bound ||
           sflint::fmpz_equal(sflint::FmpzConstRef(left_bound),
                              sflint::FmpzConstRef(right_bound)));
    assert(left.factor_base_generation_checked_status() ==
           right.factor_base_generation_checked_status());
    assert(left.class_unit_completeness_certification() ==
           right.class_unit_completeness_certification());
    assert(left.class_unit_completeness_basis() ==
           right.class_unit_completeness_basis());
    assert(left.analytic_class_regulator_status() ==
           right.analytic_class_regulator_status());
    assert(left.analytic_class_regulator_certification() ==
           right.analytic_class_regulator_certification());
    assert(left.zeta_bf_proof_status() == right.zeta_bf_proof_status());
    assert(left.relation_saturation_status() ==
           right.relation_saturation_status());
    assert(left.unit_proof_status() == right.unit_proof_status());
    assert(left.regulator_proof_status() == right.regulator_proof_status());
}

// Reusing a context does not carry provenance from one result to the next.
// Q(sqrt(10)), h = 2.  The same context and unit group run grh, then
// proven, then grh again, and each result reports what a fresh run with that
// request reports.  A dependent relation appended after a grh acceptance
// leaves the record as it was (current behaviour).  A failed certification attempt
// inside a transaction (an index test against a wrong hR, which runs the
// generation check and sets the unit proof states before it fails) leaves
// every field as it was.
int test_class_unit_provenance_after_context_reuse() {
    silex::NumberField field = quadratic_field(10);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                               order));
    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;

    auto run = [&](silex::ClassGroupContext& context,
                   silex::OrderUnitGroup& units,
                   silex::CertificationMode requested) noexcept {
        options.requested_certification = requested;
        assert(units.compute_with_class_group(
                context, order, sflint::FmpzConstRef(bound), options, 128));
        assert(context.certification_status() == requested);
    };

    silex::ClassGroupContext fresh_grh;
    silex::OrderUnitGroup fresh_grh_units;
    run(fresh_grh, fresh_grh_units, silex::CertificationMode::grh);
    assert(fresh_grh.class_unit_completeness_certification() ==
           silex::CertificationMode::grh);
    assert(fresh_grh.class_unit_completeness_basis() ==
           silex::ClassUnitCompletenessBasis::belabas_friedman);
    silex::ClassGroupContext fresh_proven;
    silex::OrderUnitGroup fresh_proven_units;
    run(fresh_proven, fresh_proven_units, silex::CertificationMode::proven);
    assert(fresh_proven.class_unit_completeness_certification() ==
           silex::CertificationMode::proven);

    silex::ClassGroupContext reused;
    silex::OrderUnitGroup reused_units;
    run(reused, reused_units, silex::CertificationMode::grh);
    assert_same_class_unit_provenance(reused, reused_units, fresh_grh,
                                      fresh_grh_units);
    run(reused, reused_units, silex::CertificationMode::proven);
    assert_same_class_unit_provenance(reused, reused_units, fresh_proven,
                                      fresh_proven_units);
    run(reused, reused_units, silex::CertificationMode::grh);
    assert_same_class_unit_provenance(reused, reused_units, fresh_grh,
                                      fresh_grh_units);

    // A failed transaction restores every field, including the completeness
    // record.  The index test against twice the true hR fails (quotient
    // 1/2) after it has checked generation and set the unit proof states.
    assert(reused.factor_base_generation_status() ==
           silex::ProofState::verified);
    {
        sflint::Arb wrong_hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(wrong_hR),
                                                   order, 128));
        ::arb_mul_2exp_si(wrong_hR.raw(), wrong_hR.raw(), 1);
        sflint::Fmpz aux_bound;
        assert(set_fmpz_si(aux_bound, 31));
        assert(!CertificationAccess::try_analytic_index_bound_with_units(
                reused, fresh_proven_units, sflint::ArbConstRef(wrong_hR),
                sflint::FmpzConstRef(aux_bound), 128, true));
    }
    assert_same_class_unit_provenance(reused, reused_units, fresh_grh,
                                      fresh_grh_units);

    // A relation appended after the grh acceptance.  (2) is the square of
    // the ramified prime above 2, and the accepted module already holds
    // it: the append is skipped as dependent and returns before the
    // metadata reset, so the accepted record stands.  (A reset happens only
    // for a retained row, and an accepted module has full rank and no
    // index to refine.)  This pins current behaviour; it is not a
    // revert-fail test.
    {
        const silex::FactorBase* base = reused.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element two(field);
        assert(two.set_si(2));
        assert(relation.set_generator(two));
        assert(reused.append_relation(relation));
    }
    assert_same_class_unit_provenance(reused, reused_units, fresh_grh,
                                      fresh_grh_units);
    return 0;
}

// A real-quadratic pair accepted against a GRH-conditional hR records
// completeness `grh` with basis `belabas_friedman` although the order has
// the unconditional L(1, chi) route: the record follows the route that
// produced the value.  The value is seeded into the validation cache, as
// when the L(1, chi) evaluation fails and the zeta route falls back to
// Belabas-Friedman (a seeded value is never unconditional).  Q(sqrt(10)),
// h = 2.
int test_real_quadratic_conditional_hr_records_belabas_friedman() {
    silex::NumberField field = quadratic_field(10);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(silex::detail::zeta_unconditional_route_available(order));
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                               order));
    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;

    // An unlabelled candidate with generation checked to the required
    // bound, and proven units.
    silex::ClassGroupContext context;
    silex::ClassGroupCandidateOptions candidate_options;
    candidate_options.max_candidates = 5000;
    candidate_options.max_relations = 500;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(bound),
                                     candidate_options));
    assert(context.has_presentation() || context.publish_presentation());
    sflint::Fmpz class_order;
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    sflint::Fmpz required_bound;
    assert(context.factor_base_generation_bound(
            sflint::FmpzRef(required_bound)));
    assert(context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(required_bound)));

    silex::OrderUnitGroup units;
    assert(units.compute(order));
    silex::OrderUnitGroup scratch(order);
    silex::EmbeddingContext embeddings(field);
    options.requested_certification = silex::CertificationMode::grh;
    options.max_relations = 500;
    silex::detail::AnalyticClassRegulatorCache analytic_cache;
    sflint::Arb seeded;
    assert(silex::zeta_class_regulator_product(sflint::ArbRef(seeded), order,
                                               128));
    analytic_cache.seed(order, sflint::ArbConstRef(seeded), 128);
    assert(!analytic_cache.value_unconditional());
    silex::detail::ValidateRefineSummary summary;
    assert(silex::detail::try_validate_refine_loop(
            context, units, scratch, order, options, embeddings,
            analytic_cache, summary, 128, false));
    assert(summary.outcome == silex::detail::ValidateRefineOutcome::proven);
    // The completeness record is reported once the grh labels are
    // published.
    assert(CertificationAccess::publish_grh_labels(context, units, nullptr));
    assert(context.certification_status() == silex::CertificationMode::grh);
    assert(context.class_unit_completeness_certification() ==
           silex::CertificationMode::grh);
    assert(context.class_unit_completeness_basis() ==
           silex::ClassUnitCompletenessBasis::belabas_friedman);
    assert(context.analytic_class_regulator_certification() ==
           silex::CertificationMode::grh);
    return 0;
}

ulong relation_saturation_test_rank_mod(
        const std::vector<std::vector<ulong>>& columns,
        slong rows,
        ulong ell) {
    if (columns.empty()) {
        return 0;
    }
    nmod_mat_t matrix;
    nmod_mat_init(matrix, rows, static_cast<slong>(columns.size()), ell);
    for (std::size_t c = 0; c < columns.size(); ++c) {
        for (slong r = 0; r < rows; ++r) {
            nmod_mat_set_entry(matrix, r, static_cast<slong>(c),
                               columns[c][static_cast<std::size_t>(r)]);
        }
    }
    const slong rank = nmod_mat_rank(matrix);
    nmod_mat_clear(matrix);
    return static_cast<ulong>(rank);
}

// The ell-th-power character of a nonzero residue image; fails on zero.
bool relation_saturation_test_character(
        ulong& out,
        const silex::ResidueFieldElement& image,
        sflint::FmpzConstRef ell) {
    sflint::Fmpz log;
    if (image.is_zero() ||
        !image.quotient_log_mod_prime(sflint::FmpzRef(log), ell)) {
        return false;
    }
    out = fmpz_get_ui(log.raw());
    return true;
}

// Whether the last row k of the stack whose integral Gram-Schmidt data are
// `d` and `lambda` is size-reduced against the rows before it:
// |mu(k, j)| <= 1/2, that is 2 |lambda(k, j)| <= d_{j+1}, for all j < k.
bool relation_saturation_test_size_reduced(const sflint::FmpzMat& d,
                                           const sflint::FmpzMat& lambda,
                                           slong k) {
    for (slong j = 0; j < k; ++j) {
        sflint::Fmpz twice;
        fmpz_mul_2exp(twice.raw(),
                      sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(lambda),
                                             k, j)
                              .raw(),
                      1);
        fmpz_abs(twice.raw(), twice.raw());
        if (fmpz_cmp(twice.raw(),
                     sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(d), 0,
                                            j + 1)
                             .raw()) > 0) {
            return false;
        }
    }
    return true;
}

// Integral Gram-Schmidt data and size reduction (Cohen, GTM 138,
// Algorithm 2.6.7) on the basis b_0 = (3, 0), b_1 = (1, 3): d = (1, 9, 81),
// lambda(1, 0) = 3, and (10, 7) reduces to (10, 7) - 2 b_1 - 3 b_0 =
// (-1, 1), with mu = -1/3 and 1/3.
int test_relation_saturation_size_reduce_row() {
    sflint::FmpzMat basis(2, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(basis), 0, 0),
                        3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(basis), 1, 0),
                        1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(basis), 1, 1),
                        3);
    sflint::FmpzMat d(0, 0);
    sflint::FmpzMat lambda(0, 0);
    assert(silex::detail::integral_gram_schmidt(
            d, lambda, sflint::FmpzMatConstRef(basis)));
    assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(d)) == 1 &&
           sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(d)) == 3);
    const slong expected_d[3] = {1, 9, 81};
    for (slong j = 0; j < 3; ++j) {
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(d), 0, j),
                expected_d[j]));
    }
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(lambda), 1, 0), 3));

    sflint::FmpzMat row(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, 0),
                        10);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, 1),
                        7);
    assert(silex::detail::size_reduce_row_against_basis(
            sflint::FmpzMatRef(row), sflint::FmpzMatConstRef(basis),
            sflint::FmpzMatConstRef(d), sflint::FmpzMatConstRef(lambda)));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(row), 0, 0), -1));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(row), 0, 1), 1));

    // Dependent rows have no Gram-Schmidt data.
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(basis), 1, 0),
                        6);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(basis), 1, 1),
                        0);
    assert(!silex::detail::integral_gram_schmidt(
            d, lambda, sflint::FmpzMatConstRef(basis)));
    return 0;
}

// Replaces `kernel`, the presentation's exact relation-kernel basis, by the
// basis relation_saturation_beta_coefficients size-reduces against, after
// checking it: it spans the same lattice, is LLL-reduced with FLINT's
// default parameters (fmpz_mat_is_reduced), and is cached.  For more than
// one row and entries of at most 250 bits it comes from FLINT's
// double-precision first phase when that phase is certified, so the
// LLL-reduced basis read by relation_kernel_unit is left unbuilt;
// `expect_first_phase` asserts that route.  Otherwise it is a copy of that
// LLL-reduced basis, which is then built.
void relation_saturation_test_kernel_basis(
        sflint::FmpzMat& kernel,
        const silex::ClassGroupContext& context,
        bool expect_first_phase) {
    const slong rows = sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(kernel));
    const slong cols = sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(kernel));
    const slong raw_bits = ::fmpz_mat_max_bits(kernel.raw());
    const bool eligible = rows > 1 && raw_bits >= -250 && raw_bits <= 250;
    sflint::FmpzMat basis(0, 0);
    bool lll_basis_built = false;
    assert(silex::detail::relation_saturation_kernel_basis(
            basis, lll_basis_built, context));
    assert(eligible || lll_basis_built);
    assert(!expect_first_phase || (eligible && !lll_basis_built));
    assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(basis)) == rows &&
           sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(basis)) == cols);
    sflint::FmpzLll lll;
    assert(::fmpz_mat_is_reduced(basis.raw(), lll.raw()->delta,
                                 lll.raw()->eta));
    if (expect_first_phase) {
        // With FLINT's default routing a certified first phase is exactly
        // what fmpz_lll returns, so the basis equals it entry for entry.
        sflint::FmpzMat reference(rows, cols);
        ::fmpz_mat_set(reference.raw(), kernel.raw());
        ::fmpz_lll(reference.raw(), nullptr, lll.raw());
        assert(::fmpz_mat_equal(reference.raw(), basis.raw()));
    }
    sflint::FmpzMat basis_hnf(rows, cols);
    sflint::FmpzMat kernel_hnf(rows, cols);
    ::fmpz_mat_hnf(basis_hnf.raw(), basis.raw());
    ::fmpz_mat_hnf(kernel_hnf.raw(), kernel.raw());
    assert(::fmpz_mat_equal(basis_hnf.raw(), kernel_hnf.raw()));
    sflint::FmpzMat again(0, 0);
    bool built_again = false;
    assert(silex::detail::relation_saturation_kernel_basis(
            again, built_again, context));
    assert(::fmpz_mat_equal(again.raw(), basis.raw()));
    assert(built_again == lll_basis_built);
    kernel = std::move(basis);
}

// The saturation kernel basis from FLINT's certified first LLL phase.
// Q(sqrt(100000001)), x^2 - x - 25000000, in GRH mode: the relation kernel
// has more than one row and entries of at most 250 bits, so the basis the
// beta rows are size-reduced against comes from
// fmpz_lll_d_with_removal_knapsack certified by the exact certificate
// (lat::detail::lll_certificate_holds), without building the LLL-reduced
// basis of relation_kernel_unit.  With
// proven units, each beta coefficient row differs from the presentation's
// row by a kernel vector and is size-reduced against that basis.
int test_relation_saturation_kernel_basis_first_phase() {
    silex::NumberField field = shifted_quadratic_field(-1, -25000000);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), maximal_order));

    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext context;
    silex::OrderUnitGroup grh_units;
    assert(grh_units.compute_with_class_group(
            context, maximal_order, sflint::FmpzConstRef(factor_base_bound),
            options, 128));

    silex::FiniteAbelianGroup presentation;
    assert(context.presentation(presentation));
    const slong relation_count = context.relation_count();
    sflint::FmpzMat kernel(0, 0);
    {
        std::optional<sflint::FmpzMat> exact_kernel =
                presentation.relation_kernel_matrix();
        assert(exact_kernel.has_value());
        kernel = std::move(*exact_kernel);
    }
    const slong kernel_rows =
            sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(kernel));
    relation_saturation_test_kernel_basis(kernel, context, true);

    silex::OrderUnitGroup units(maximal_order);
    silex::EmbeddingContext embeddings(field);
    silex::ProofState unit_status = silex::ProofState::not_checked;
    bool units_changed = false;
    sflint::Fmpz unit_aux_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(unit_aux_bound), 1000);
    assert(units.prove_index_bound(unit_status, units_changed, grh_units, 1,
                                   sflint::FmpzConstRef(unit_aux_bound), 2,
                                   embeddings, 128));
    assert(unit_status == silex::ProofState::verified);
    assert(units.certification_status() == silex::CertificationMode::proven);

    sflint::FmpzMat raw(context.invariant_count(), relation_count);
    assert(presentation.invariant_generator_relation_matrix(
            sflint::FmpzMatRef(raw)));
    sflint::Fmpz invariant;
    slong checked_rows = 0;
    for (slong i = 0; i < context.invariant_count(); ++i) {
        assert(context.invariant(sflint::FmpzRef(invariant), i));
        sflint::Fmpz ell;
        sflint::fmpz_set_si(sflint::FmpzRef(ell), 2);
        while (!sflint::fmpz_divisible(sflint::FmpzConstRef(invariant),
                                       sflint::FmpzConstRef(ell))) {
            fmpz_nextprime(ell.raw(), ell.raw(), 1);
        }
        sflint::FmpzMat coefficients(0, 0);
        assert(silex::detail::relation_saturation_beta_coefficients(
                coefficients, context, units, sflint::FmpzConstRef(ell)));
        // Row b of `coefficients` is invariant index beta_indices[b].
        std::vector<slong> beta_indices;
        for (slong k = 0; k < context.invariant_count(); ++k) {
            sflint::Fmpz value;
            assert(context.invariant(sflint::FmpzRef(value), k));
            if (sflint::fmpz_divisible(sflint::FmpzConstRef(value),
                                       sflint::FmpzConstRef(ell))) {
                beta_indices.push_back(k);
            }
        }
        const slong beta_count = static_cast<slong>(beta_indices.size());
        assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(coefficients)) ==
               beta_count);
        for (slong b = 0; b < beta_count; ++b) {
            sflint::FmpzMat difference(kernel_rows + 1, relation_count);
            sflint::FmpzMat reduced_stack(kernel_rows + 1, relation_count);
            for (slong r = 0; r < kernel_rows; ++r) {
                for (slong j = 0; j < relation_count; ++j) {
                    fmpz_set(sflint::fmpz_mat_entry(
                                     sflint::FmpzMatRef(difference), r, j)
                                     .raw(),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(kernel), r, j)
                                     .raw());
                    fmpz_set(sflint::fmpz_mat_entry(
                                     sflint::FmpzMatRef(reduced_stack), r, j)
                                     .raw(),
                             sflint::fmpz_mat_entry(
                                     sflint::FmpzMatConstRef(kernel), r, j)
                                     .raw());
                }
            }
            for (slong j = 0; j < relation_count; ++j) {
                const sflint::FmpzConstRef reduced_entry =
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(coefficients), b, j);
                fmpz_set(sflint::fmpz_mat_entry(
                                 sflint::FmpzMatRef(reduced_stack),
                                 kernel_rows, j)
                                 .raw(),
                         reduced_entry.raw());
                fmpz_sub(sflint::fmpz_mat_entry(
                                 sflint::FmpzMatRef(difference),
                                 kernel_rows, j)
                                 .raw(),
                         sflint::fmpz_mat_entry(
                                 sflint::FmpzMatConstRef(raw),
                                 beta_indices[static_cast<std::size_t>(b)], j)
                                 .raw(),
                         reduced_entry.raw());
            }
            assert(fmpz_mat_rank(difference.raw()) == kernel_rows);
            sflint::FmpzMat d(0, 0);
            sflint::FmpzMat lambda(0, 0);
            assert(silex::detail::integral_gram_schmidt(
                    d, lambda, sflint::FmpzMatConstRef(reduced_stack)));
            assert(relation_saturation_test_size_reduced(d, lambda,
                                                         kernel_rows));
            ++checked_rows;
        }
    }
    assert(checked_rows > 0);
    // The size reduction left the LLL-reduced basis unbuilt.
    sflint::FmpzMat basis(0, 0);
    bool lll_basis_built = true;
    assert(silex::detail::relation_saturation_kernel_basis(
            basis, lll_basis_built, context));
    assert(!lll_basis_built);
    return 0;
}

// Checks `basis`, a relation-saturation kernel basis read back from
// `context`, against the presentation: it equals fmpz_lll of the exact
// relation-kernel basis entry for entry (so it spans the kernel), and the
// cached integral Gram-Schmidt data equal integral_gram_schmidt(basis).
void relation_saturation_test_kernel_basis_is_lll(
        const sflint::FmpzMat& basis,
        const silex::ClassGroupContext& context) {
    silex::FiniteAbelianGroup presentation;
    assert(context.presentation(presentation));
    std::optional<sflint::FmpzMat> kernel =
            presentation.relation_kernel_matrix();
    assert(kernel.has_value());
    const slong rows = sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(*kernel));
    assert(rows == context.relation_kernel_unit_count());
    sflint::FmpzLll lll;
    ::fmpz_lll(kernel->raw(), nullptr, lll.raw());
    assert(::fmpz_mat_equal(kernel->raw(), basis.raw()));

    sflint::FmpzMat d(0, 0);
    sflint::FmpzMat lambda(0, 0);
    assert(silex::detail::relation_saturation_kernel_gram(d, lambda, context));
    sflint::FmpzMat expected_d(0, 0);
    sflint::FmpzMat expected_lambda(0, 0);
    assert(silex::detail::integral_gram_schmidt(
            expected_d, expected_lambda, sflint::FmpzMatConstRef(basis)));
    assert(sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(d)) == rows + 1);
    assert(::fmpz_mat_equal(d.raw(), expected_d.raw()));
    assert(::fmpz_mat_equal(lambda.raw(), expected_lambda.raw()));
}

// Builds the relation-saturation kernel basis of `context` with `failure`
// injected into its first phase and checks the route: the phase ran and its
// result was used (none), or the injected failure replaced it once and the
// basis came from the fmpz_lll fallback, which builds the LLL-reduced basis
// of relation_kernel_unit.
sflint::FmpzMat relation_saturation_test_kernel_basis_with_failure(
        const silex::ClassGroupContext& context,
        silex::detail::CertifiedLllTestFailure failure) {
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
    const slong injected = TestHookAccess::certified_lll_failures_injected();
    TestHookAccess::set_certified_lll_failure(failure);
    sflint::FmpzMat basis(0, 0);
    bool lll_basis_built = false;
    const bool built = silex::detail::relation_saturation_kernel_basis(
            basis, lll_basis_built, context);
    TestHookAccess::set_certified_lll_failure(
            CertifiedLllTestFailure::none);
    const bool injecting = failure != CertifiedLllTestFailure::none;
    assert(built);
    assert(lll_basis_built == injecting);
    assert(TestHookAccess::certified_lll_failures_injected() ==
           injected + (injecting ? 1 : 0));
    relation_saturation_test_kernel_basis_is_lll(basis, context);
    return basis;
}

// Over Q with the factor base {2}, the relations 256, 16, 4, 2 have a
// relation kernel of rank 3 with small entries, so the first phase builds
// the saturation kernel basis.  When the phase reports failure or its
// result is not certified, the fallback builds the basis from fmpz_lll,
// equal to the basis the certified phase gives.  Publishing the
// presentation again resets the cached bases, so each case builds afresh.
int test_relation_saturation_kernel_basis_fallbacks() {
    using silex::detail::CertifiedLllTestFailure;
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    for (const slong generator : {256, 16, 4, 2}) {
        assert(alpha.set_si(generator));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
    }
    assert(context.relation_count() == 4);

    sflint::FmpzMat expected(0, 0);
    for (const CertifiedLllTestFailure failure :
         {CertifiedLllTestFailure::none,
          CertifiedLllTestFailure::first_phase,
          CertifiedLllTestFailure::certificate}) {
        assert(context.publish_presentation());
        assert(context.relation_kernel_unit_count() == 3);
        sflint::FmpzMat basis =
                relation_saturation_test_kernel_basis_with_failure(context,
                                                                   failure);
        if (failure == CertifiedLllTestFailure::none) {
            expected = std::move(basis);
        } else {
            assert(::fmpz_mat_equal(basis.raw(), expected.raw()));
        }
    }
    return 0;
}

// The saturation kernel basis is a copy of the LLL-reduced basis of
// relation_kernel_unit when that is already built, and both bases and the
// Gram data are rebuilt for the new kernel after the presentation changes.
// Over Q with the factor base {2}: the relations 256, 16, 4 have a kernel
// of rank 2, and adding 2 gives one of rank 3.
int test_relation_saturation_kernel_basis_reuse_and_invalidation() {
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    for (const slong generator : {256, 16, 4}) {
        assert(alpha.set_si(generator));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
    }
    assert(context.publish_presentation());
    assert(context.relation_kernel_unit_count() == 2);

    // relation_kernel_unit builds the LLL-reduced basis first.  With every
    // first-phase run failing, the saturation basis is still built and no
    // failure is injected: the phase did not run, and the basis is the
    // copy of that LLL-reduced basis.
    silex::FactoredElement unit_witness(field);
    assert(context.relation_kernel_unit(unit_witness, 0));
    const slong injected = TestHookAccess::certified_lll_failures_injected();
    TestHookAccess::set_certified_lll_failure(
            CertifiedLllTestFailure::first_phase);
    sflint::FmpzMat basis(0, 0);
    bool lll_basis_built = false;
    const bool built = silex::detail::relation_saturation_kernel_basis(
            basis, lll_basis_built, context);
    TestHookAccess::set_certified_lll_failure(
            CertifiedLllTestFailure::none);
    assert(built && lll_basis_built);
    assert(TestHookAccess::certified_lll_failures_injected() == injected);
    assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(basis)) == 2);
    relation_saturation_test_kernel_basis_is_lll(basis, context);

    // A new presentation resets both bases and the Gram data: the first
    // phase builds the rank-3 basis, and relation_kernel_unit's basis is
    // left unbuilt.
    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    assert(context.relation_kernel_unit_count() == 3);
    basis = relation_saturation_test_kernel_basis_with_failure(
            context, CertifiedLllTestFailure::none);
    assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(basis)) == 3);
    return 0;
}

// Forcing the fallbacks of the certified first LLL phase, including the
// saturation kernel basis, leaves a proven computation unchanged.
// Q(sqrt(100000001)), x^2 - x - 25000000: its relation kernel takes the
// certified first phase (see
// test_relation_saturation_kernel_basis_first_phase).  With the first phase
// failing, or its certificate failing, on every run at every call site
// (CertifiedLllTestSeams), the proven class group, its
// saturation records and the units are the same.
int test_relation_saturation_kernel_fallbacks_keep_proven_result() {
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
    silex::NumberField field = shifted_quadratic_field(-1, -25000000);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), maximal_order));

    struct ProvenResult {
        std::vector<sflint::Fmpz> invariants;
        std::vector<sflint::Fmpz> saturation_ells;
        std::vector<silex::ProofState> saturation_states;
        sflint::Arb regulator;
    };
    auto compute = [&](CertifiedLllTestFailure failure) {
        silex::ClassGroupComputeOptions options;
        options.zeta_bf_max_cutoff = 20000;
        options.requested_certification = silex::CertificationMode::proven;
        silex::ClassGroupContext context;
        silex::OrderUnitGroup units;
        const slong injected =
                TestHookAccess::certified_lll_failures_injected();
        TestHookAccess::set_certified_lll_failure(failure);
        const bool computed = units.compute_with_class_group(
                context, maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128);
        TestHookAccess::set_certified_lll_failure(
                CertifiedLllTestFailure::none);
        assert(computed);
        const slong injections =
                TestHookAccess::certified_lll_failures_injected() -
                injected;
        assert(failure == CertifiedLllTestFailure::none
                       ? injections == 0
                       : injections > 0);
        assert_proven_class_group_certification(context);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);

        ProvenResult result;
        for (slong i = 0; i < context.invariant_count(); ++i) {
            result.invariants.emplace_back();
            assert(context.invariant(
                    sflint::FmpzRef(result.invariants.back()), i));
        }
        assert(context.relation_saturation_record_count() > 0);
        for (slong i = 0; i < context.relation_saturation_record_count();
             ++i) {
            result.saturation_ells.emplace_back();
            silex::ProofState state = silex::ProofState::not_checked;
            assert(context.relation_saturation_record(
                    sflint::FmpzRef(result.saturation_ells.back()), state,
                    i));
            result.saturation_states.push_back(state);
        }
        assert(units.regulator(sflint::ArbRef(result.regulator)));
        return result;
    };

    const ProvenResult expected = compute(CertifiedLllTestFailure::none);
    assert(!expected.invariants.empty());
    for (const CertifiedLllTestFailure failure :
         {CertifiedLllTestFailure::first_phase,
          CertifiedLllTestFailure::certificate}) {
        const ProvenResult result = compute(failure);
        assert(result.invariants.size() == expected.invariants.size());
        for (std::size_t i = 0; i < expected.invariants.size(); ++i) {
            assert(sflint::fmpz_equal(
                    sflint::FmpzConstRef(result.invariants[i]),
                    sflint::FmpzConstRef(expected.invariants[i])));
        }
        assert(result.saturation_ells.size() ==
               expected.saturation_ells.size());
        for (std::size_t i = 0; i < expected.saturation_ells.size(); ++i) {
            assert(sflint::fmpz_equal(
                    sflint::FmpzConstRef(result.saturation_ells[i]),
                    sflint::FmpzConstRef(expected.saturation_ells[i])));
            assert(result.saturation_states[i] ==
                   expected.saturation_states[i]);
        }
        assert(::arb_equal(result.regulator.raw(),
                           expected.regulator.raw()));
    }
    return 0;
}

// The relation-saturation beta rows are size-reduced against an
// LLL-reduced relation-kernel basis.  Q(sqrt(100000005)), x^2 - x - 25000001:
// h = 212 with invariants [106, 2], w = 2, unit rank 1, regulator
// 18.19753721763815524762... (certified independently).  Each
// reduced coefficient row differs from the presentation's row for
// g_i^{d_i} by a relation-kernel vector and is size-reduced against the
// reduced kernel; the raw rows of this presentation are not.  A reduced
// beta_i therefore differs from the raw one by a unit, so at the same
// auxiliary primes the raw and reduced ell-th-power characters give the
// same rank modulo ell: for ell = 2 with the torsion row and for ell = 53
// without it (-1 is a 53rd power).  The raw characters are the linear
// combinations of the relation-generator characters.
int test_relation_saturation_beta_rows_size_reduced() {
    silex::NumberField field = shifted_quadratic_field(-1, -25000001);
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), maximal_order));

    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext context;
    silex::OrderUnitGroup grh_units;
    assert(grh_units.compute_with_class_group(
            context, maximal_order, sflint::FmpzConstRef(factor_base_bound),
            options, 128));
    sflint::Fmpz value;
    assert(context.order(sflint::FmpzRef(value)));
    assert(sflint::fmpz_equal_si(value, 212));

    silex::OrderUnitGroup units(maximal_order);
    silex::EmbeddingContext embeddings(field);
    silex::ProofState unit_status = silex::ProofState::not_checked;
    bool units_changed = false;
    sflint::Fmpz unit_aux_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(unit_aux_bound), 1000);
    assert(units.prove_index_bound(unit_status, units_changed, grh_units, 1,
                                   sflint::FmpzConstRef(unit_aux_bound), 2,
                                   embeddings, 128));
    assert(unit_status == silex::ProofState::verified);
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 1);
    silex::FactoredElement fundamental_unit(field);
    assert(units.free_generator(fundamental_unit, 0));
    silex::OrderElement torsion(maximal_order);
    assert(units.torsion_generator(torsion));

    silex::FiniteAbelianGroup presentation;
    assert(context.presentation(presentation));
    const slong relation_count = context.relation_count();
    sflint::FmpzMat raw(context.invariant_count(), relation_count);
    assert(presentation.invariant_generator_relation_matrix(
            sflint::FmpzMatRef(raw)));
    sflint::FmpzMat kernel(0, 0);
    {
        std::optional<sflint::FmpzMat> exact_kernel =
                presentation.relation_kernel_matrix();
        assert(exact_kernel.has_value());
        kernel = std::move(*exact_kernel);
    }
    const slong kernel_rows =
            sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(kernel));
    assert(kernel_rows > 1);
    // Size reduction below is checked against the basis the beta rows use.
    relation_saturation_test_kernel_basis(kernel, context, false);
    bool raw_size_reduced = true;

    std::vector<silex::Element> generators;
    generators.reserve(static_cast<std::size_t>(relation_count));
    for (slong j = 0; j < relation_count; ++j) {
        generators.emplace_back(field);
        assert(context.relation_generator(generators.back(), j));
    }

    for (const slong ell_value : {2, 53}) {
        sflint::Fmpz ell;
        sflint::fmpz_set_si(sflint::FmpzRef(ell), ell_value);
        const ulong ell_ui = static_cast<ulong>(ell_value);

        std::vector<silex::FactoredElement> reduced;
        assert(silex::detail::relation_saturation_beta_rows(
                reduced, context, units, sflint::FmpzConstRef(ell)));
        {
            // Units that are not proven are refused.
            std::vector<silex::FactoredElement> refused;
            assert(!silex::detail::relation_saturation_beta_rows(
                    refused, context, grh_units, sflint::FmpzConstRef(ell)));
        }
        std::vector<slong> beta_indices;
        for (slong i = 0; i < context.invariant_count(); ++i) {
            assert(context.invariant(sflint::FmpzRef(value), i));
            if (sflint::fmpz_divisible(sflint::FmpzConstRef(value),
                                       sflint::FmpzConstRef(ell))) {
                beta_indices.push_back(i);
            }
        }
        const slong beta_count = static_cast<slong>(beta_indices.size());
        assert(beta_count > 0 && beta_count == static_cast<slong>(reduced.size()));

        // Each reduced coefficient row differs from the raw row of
        // g_i^{d_i} by a relation-kernel vector (the integer kernel is
        // saturated, so membership in its rational span suffices) and is
        // size-reduced against that kernel basis: |mu(row, l)| <= 1/2.
        sflint::FmpzMat coefficients(0, 0);
        assert(silex::detail::relation_saturation_beta_coefficients(
                coefficients, context, units, sflint::FmpzConstRef(ell)));
        {
            sflint::FmpzMat refused(0, 0);
            assert(!silex::detail::relation_saturation_beta_coefficients(
                    refused, context, grh_units, sflint::FmpzConstRef(ell)));
        }
        assert(sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(coefficients)) ==
               beta_count);
        for (slong b = 0; b < beta_count; ++b) {
            sflint::FmpzMat difference(kernel_rows + 1, relation_count);
            sflint::FmpzMat raw_stack(kernel_rows + 1, relation_count);
            sflint::FmpzMat reduced_stack(kernel_rows + 1, relation_count);
            for (slong r = 0; r < kernel_rows; ++r) {
                for (slong j = 0; j < relation_count; ++j) {
                    const sflint::FmpzConstRef entry = sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(kernel), r, j);
                    for (sflint::FmpzMat* stack :
                         {&difference, &raw_stack, &reduced_stack}) {
                        sflint::fmpz_set(sflint::fmpz_mat_entry(
                                                 sflint::FmpzMatRef(*stack),
                                                 r, j),
                                         entry);
                    }
                }
            }
            for (slong j = 0; j < relation_count; ++j) {
                const sflint::FmpzConstRef raw_entry = sflint::fmpz_mat_entry(
                        sflint::FmpzMatConstRef(raw),
                        beta_indices[static_cast<std::size_t>(b)], j);
                const sflint::FmpzConstRef reduced_entry =
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(coefficients), b, j);
                sflint::fmpz_set(sflint::fmpz_mat_entry(
                                         sflint::FmpzMatRef(raw_stack),
                                         kernel_rows, j),
                                 raw_entry);
                sflint::fmpz_set(sflint::fmpz_mat_entry(
                                         sflint::FmpzMatRef(reduced_stack),
                                         kernel_rows, j),
                                 reduced_entry);
                fmpz_sub(sflint::fmpz_mat_entry(
                                 sflint::FmpzMatRef(difference),
                                 kernel_rows, j)
                                 .raw(),
                         raw_entry.raw(), reduced_entry.raw());
            }
            assert(fmpz_mat_rank(difference.raw()) == kernel_rows);
            sflint::FmpzMat d(0, 0);
            sflint::FmpzMat lambda(0, 0);
            assert(silex::detail::integral_gram_schmidt(
                    d, lambda, sflint::FmpzMatConstRef(reduced_stack)));
            assert(relation_saturation_test_size_reduced(d, lambda,
                                                         kernel_rows));
            assert(silex::detail::integral_gram_schmidt(
                    d, lambda, sflint::FmpzMatConstRef(raw_stack)));
            raw_size_reduced =
                    raw_size_reduced &&
                    relation_saturation_test_size_reduced(d, lambda,
                                                          kernel_rows);
        }
        const bool include_torsion = ell_value == 2;
        const slong rows = beta_count + 1 + (include_torsion ? 1 : 0);

        std::vector<std::vector<ulong>> raw_columns;
        std::vector<std::vector<ulong>> reduced_columns;
        ulong rank = 0;
        for (ulong q = 2; q <= 10000 && rank < static_cast<ulong>(rows);
             q = n_nextprime(q, 1)) {
            if ((q - 1) % ell_ui != 0) {
                continue;
            }
            sflint::Fmpz rational_prime;
            sflint::fmpz_set_ui(sflint::FmpzRef(rational_prime), q);
            silex::PrimeIdealList above;
            assert(silex::decompose_prime(
                    above, maximal_order,
                    sflint::FmpzConstRef(rational_prime)));
            for (slong k = 0; k < above.size() &&
                              rank < static_cast<ulong>(rows);
                 ++k) {
                const silex::PrimeIdeal* prime = above.at(k);
                assert(prime != nullptr);
                if (prime->residue_degree() != 1) {
                    continue;
                }
                silex::ResidueField residue_field(*prime);
                assert(residue_field.is_defined());
                silex::ResidueFieldElement image(residue_field);

                // Primes dividing a relation generator, a beta row or a unit
                // give no character; skip them.
                std::vector<ulong> relation_characters;
                bool coprime = true;
                for (const silex::Element& generator : generators) {
                    ulong character = 0;
                    coprime = image.set_element(generator) &&
                              relation_saturation_test_character(
                                      character, image,
                                      sflint::FmpzConstRef(ell));
                    if (!coprime) {
                        break;
                    }
                    relation_characters.push_back(character);
                }
                std::vector<ulong> raw_column;
                std::vector<ulong> reduced_column;
                for (slong b = 0; coprime && b < beta_count; ++b) {
                    sflint::Fmpz sum;
                    for (slong j = 0; j < relation_count; ++j) {
                        sflint::fmpz_addmul_ui(
                                sflint::FmpzRef(sum),
                                sflint::fmpz_mat_entry(
                                        sflint::FmpzMatConstRef(raw),
                                        beta_indices[static_cast<std::size_t>(
                                                b)],
                                        j),
                                relation_characters[static_cast<std::size_t>(
                                        j)]);
                    }
                    raw_column.push_back(fmpz_fdiv_ui(sum.raw(), ell_ui));
                    ulong character = 0;
                    coprime = image.set_factored_element(
                                      reduced[static_cast<std::size_t>(b)]) &&
                              relation_saturation_test_character(
                                      character, image,
                                      sflint::FmpzConstRef(ell));
                    reduced_column.push_back(character);
                }
                ulong unit_character = 0;
                coprime = coprime &&
                          image.set_factored_element(fundamental_unit) &&
                          relation_saturation_test_character(
                                  unit_character, image,
                                  sflint::FmpzConstRef(ell));
                ulong torsion_character = 0;
                coprime = coprime &&
                          (!include_torsion ||
                           (image.set_order_element(torsion) &&
                            relation_saturation_test_character(
                                    torsion_character, image,
                                    sflint::FmpzConstRef(ell))));
                if (!coprime) {
                    continue;
                }
                raw_column.push_back(unit_character);
                reduced_column.push_back(unit_character);
                if (include_torsion) {
                    raw_column.push_back(torsion_character);
                    reduced_column.push_back(torsion_character);
                }
                raw_columns.push_back(std::move(raw_column));
                reduced_columns.push_back(std::move(reduced_column));

                rank = relation_saturation_test_rank_mod(reduced_columns,
                                                         rows, ell_ui);
                assert(rank == relation_saturation_test_rank_mod(
                                       raw_columns, rows, ell_ui));
            }
        }
        assert(rank == static_cast<ulong>(rows));
        assert(raw_columns.size() > 1);
    }
    // The fixture exercises the reduction: some raw row is not size-reduced.
    assert(!raw_size_reduced);

    return 0;
}

int main() {
    test_certification_metadata_defaults_and_invalidation();
    test_factor_base_generation_metadata();
    test_factor_base_honesty_receipt_requires_canonical_bound();
    test_append_and_publish();
    test_context_keeps_parent_order_alive();
    test_append_rejects_wrong_factor_base();
    test_set_relation_matrix_preserves_on_failure();
    test_relation_kernel_unit_witness();
    test_relation_kernel_unit_follows_republished_presentation();
    test_compact_power_fmpz_exponent_guard();
    test_goal_publication_waits_for_kernel_target();
    test_relation_sources_reset_with_factor_base();
    test_append_tracks_rank_and_skips_nonrefining_dependent_rows();
    test_keep_nonduplicate_skips_exact_repeat_only();
    test_append_keeps_index_refining_dependent_rows();
    test_factor_base_prime_search_masks();
    test_compute_candidate_coordinate_search();
    test_compute_candidate_skips_principal_prime_lattice_search();
    test_compute_candidate_ideal_lattice_search();
    test_compute_candidate_small_quadratic_policy();
    test_compute_candidate_proven_quadratic();
    test_detail_ideal_shortcuts_reach_labelled_dlog_as_true_ideals();
    test_certify_real_quadratic_clears_saturation_records();
    test_default_maximal_quadratic_proven_policy();
    test_nonmaximal_quadratic_candidate_rejection_preserves_output();
    test_try_certify_with_units_zeta_quadratic();
    test_compute_candidate_random_product_relation_cubic();
    test_one_large_prime_match_cubic();
    test_direct_partial_row_exact_valuation_cubic();
    test_compute_candidate_post_finite_refinement_cubic();
    test_compute_candidate_relation_kernel_target_quadratic();
    test_compute_candidate_relation_kernel_target_shifted_quadratic();
    test_compute_candidate_relation_kernel_target_quintic();
    test_relation_saturation_proof_with_units_degree_one();
    test_relation_saturation_bounded_append_with_units_degree_one();
    test_relation_saturation_candidate_filter_degree_one();
    test_relation_saturation_candidate_filter_count_gate();
    test_relation_saturation_round_settled_rows();
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    test_relation_saturation_scan_skips_settled_rows();
#endif
    test_relation_saturation_index_bound_with_units_degree_one();
    test_relation_saturation_index_bound_checks_nondivisor_primes();
    test_relation_saturation_candidate_filter_minus_23();
    test_relation_saturation_candidate_filter_needs_nice_polynomial();
    test_degree_one_quotient_log_row_word_size_matches_general();
    test_saturation_candidates_skip_bad_prime_of_nonmonic_polynomial();
    test_saturation_candidates_skip_content_prime_of_nonmonic_polynomial();
    test_relation_saturation_index_bounded_status_needs_unconditional_bound();
    test_validation_exact_index_saturation_status_follows_generation();
    test_exact_order_certification_failure_restores_generation_check();
    test_saturation_promotion_requires_every_class_order_prime();
    test_relation_saturation_completes_prime_by_prime();
    test_saturation_proof_paths_differ();
    test_class_unit_zeta_keeps_complete_saturation_records();
    test_try_certify_with_units_analytic_keeps_complete_saturation_records();
    test_analytic_class_unit_proof_requires_factor_base_generation();
    test_relation_saturation_analytic_index_bound_with_units_degree_one();
    test_zeta_bf_degree_one_verifies_relation_saturation();
    test_analytic_index_bound_keeps_complete_saturation_records();
    test_analytic_index_one_routes_clear_saturation_proof_records();
    test_low_bound_generation_check_does_not_promote();
    test_relation_saturation_size_reduce_row();
    test_relation_saturation_beta_rows_size_reduced();
    test_relation_saturation_kernel_basis_first_phase();
    test_relation_saturation_kernel_basis_fallbacks();
    test_relation_saturation_kernel_basis_reuse_and_invalidation();
    test_relation_saturation_kernel_fallbacks_keep_proven_result();
    test_grh_request_takes_belabas_friedman_residue();
    test_exact_route_keeps_verified_ell_local_proof();
    test_compute_candidate_preserves_on_failure();
    test_move_and_swap();
    test_grh_provenance_records();
    test_honesty_reuses_recorded_generation_requirement();
    test_transaction_keeps_generation_requirement();
    test_grh_generation_bound_bach_fallback();
    test_twin_class_generation_routes();
    test_imaginary_quadratic_structure_certificate_route();
    test_imaginary_quadratic_proven_large_witness_forms();
    test_imaginary_quadratic_structure_certificate_rejects_subgroup();
    test_factor_base_size_limit_fails_closed();
    test_grh_relation_search_size_limit_reason();
    test_class_unit_index_routes_share_normalization();
    test_exact_quadratic_index_below_one_fails_closed();
    test_validation_loop_exact_index_status_follows_request();
    test_class_unit_provenance_after_context_reuse();
    test_real_quadratic_conditional_hr_records_belabas_friedman();
    return 0;
}
