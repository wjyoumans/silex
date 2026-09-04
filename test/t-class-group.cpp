#include <silex/class_group.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/relation.hpp>

#include "class_group/class_group_certification_internal.hpp"
#include "class_group/relation_search_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <type_traits>
#include <utility>

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
    assert(proven_imaginary.relation_saturation_record_count() == 0);
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
    assert_proven_class_group_certification(proven_real);
    assert(proven_real.relation_saturation_record_count() == 0);
    assert(proven_real.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_is_one(order_out));

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
    assert(context.relation_saturation_status() ==
           silex::ProofState::not_checked);
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

    assert(set_fmpz_si(aux_bound, 2));
    assert(!context.try_prove_relation_saturation_with_units(
            units, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::unavailable);
    assert(context.relation_saturation_record_count() == 0);

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
    assert(!invalid_context
                    .try_prove_relation_saturation_index_bound_with_units(
                            units, sflint::FmpzConstRef(index_bound),
                            sflint::FmpzConstRef(aux_bound)));
    assert(invalid_context.relation_saturation_status() ==
           silex::ProofState::not_checked);

    silex::ClassGroupContext unavailable_context;
    prepare_context(unavailable_context);
    assert(unavailable_context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(set_fmpz_si(index_bound, 6));
    assert(!unavailable_context
                    .try_prove_relation_saturation_index_bound_with_units(
                            units, sflint::FmpzConstRef(index_bound),
                            sflint::FmpzConstRef(aux_bound)));
    assert(unavailable_context.relation_saturation_status() ==
           silex::ProofState::unavailable);
    assert(unavailable_context.relation_saturation_record_count() == 0);
    assert(unavailable_context.certification_status() ==
           silex::CertificationMode::unknown);

    silex::ClassGroupContext proven_context;
    prepare_context(proven_context);
    assert(proven_context.try_certify_with_units(
            units, silex::CertificationMode::unknown, 80));
    assert(set_fmpz_si(index_bound, 2));
    assert(proven_context
                   .try_prove_relation_saturation_index_bound_with_units(
                           units, sflint::FmpzConstRef(index_bound),
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
    // A bound I <= 5 also requires excluding 2 and 3; ell=3 is unavailable.
    assert(!context.try_prove_relation_saturation_index_bound_with_units(
            units, sflint::FmpzConstRef(index_bound),
            sflint::FmpzConstRef(aux_bound)));
    assert(context.relation_saturation_status() ==
           silex::ProofState::unavailable);
    assert(context.relation_saturation_record_count() == 0);

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
    assert(!context.try_certify_class_unit_with_units(
            units, sflint::ArbConstRef(analytic_hR), 192));
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
    assert(!invalid_analytic_context.try_analytic_index_bound_with_units(
            units, sflint::ArbConstRef(zero),
            sflint::FmpzConstRef(aux_bound), 128));
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
    assert(!unavailable_context.try_analytic_index_bound_with_units(
            units, sflint::ArbConstRef(analytic_hR),
            sflint::FmpzConstRef(aux_bound), 128));
    assert(unavailable_context.unit_proof_status() ==
           silex::ProofState::verified);
    assert(unavailable_context.regulator_proof_status() ==
           silex::ProofState::verified);
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
    assert(index_one_context.try_analytic_index_bound_with_units(
            units, sflint::ArbConstRef(exact_hR),
            sflint::FmpzConstRef(aux_bound), 128));
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
    assert(!nontrivial_bound_context.try_analytic_index_bound_with_units(
            units, sflint::ArbConstRef(analytic_hR),
            sflint::FmpzConstRef(aux_bound), 128));
    assert(nontrivial_bound_context.unit_proof_status() ==
           silex::ProofState::verified);
    assert(nontrivial_bound_context.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(nontrivial_bound_context.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(nontrivial_bound_context.relation_saturation_record_count() == 0);
    assert(nontrivial_bound_context.certification_status() ==
           silex::CertificationMode::unknown);

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

int main() {
    test_certification_metadata_defaults_and_invalidation();
    test_factor_base_generation_metadata();
    test_factor_base_honesty_receipt_requires_canonical_bound();
    test_append_and_publish();
    test_context_keeps_parent_order_alive();
    test_append_rejects_wrong_factor_base();
    test_set_relation_matrix_preserves_on_failure();
    test_relation_kernel_unit_witness();
    test_goal_publication_waits_for_kernel_target();
    test_relation_sources_reset_with_factor_base();
    test_append_tracks_rank_and_skips_nonrefining_dependent_rows();
    test_append_keeps_index_refining_dependent_rows();
    test_factor_base_prime_search_masks();
    test_compute_candidate_coordinate_search();
    test_compute_candidate_skips_principal_prime_lattice_search();
    test_compute_candidate_ideal_lattice_search();
    test_compute_candidate_small_quadratic_policy();
    test_compute_candidate_proven_quadratic();
    test_default_maximal_quadratic_proven_policy();
    test_nonmaximal_quadratic_candidate_rejection_preserves_output();
    test_try_certify_with_units_zeta_quadratic();
    test_compute_candidate_random_product_relation_cubic();
    test_one_large_prime_match_cubic();
    test_compute_candidate_post_finite_refinement_cubic();
    test_compute_candidate_relation_kernel_target_quadratic();
    test_compute_candidate_relation_kernel_target_shifted_quadratic();
    test_compute_candidate_relation_kernel_target_quintic();
    test_relation_saturation_proof_with_units_degree_one();
    test_relation_saturation_bounded_append_with_units_degree_one();
    test_relation_saturation_index_bound_with_units_degree_one();
    test_relation_saturation_index_bound_checks_nondivisor_primes();
    test_analytic_class_unit_proof_requires_factor_base_generation();
    test_relation_saturation_analytic_index_bound_with_units_degree_one();
    test_compute_candidate_preserves_on_failure();
    test_move_and_swap();
    return 0;
}
