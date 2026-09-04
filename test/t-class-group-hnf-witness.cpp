#include <silex/class_group.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>
#include <silex/relation.hpp>

#include "class_group/class_group_internal.hpp"
#include "test_support.hpp"

#include <cassert>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField cubic_discriminant_81_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::ClassGroupContext relation_context(
        silex::NumberField& field,
        silex::Order& order) noexcept {
    field = degree_one_field();
    order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr && base->length() == 2);

    silex::Relation relation(*base);
    silex::Element generator(field);
    for (slong value : {WORD(4), WORD(9), WORD(6)}) {
        assert(generator.set_si(value));
        assert(relation.set_generator(generator));
        assert(context.append_relation(relation));
    }
    assert(context.relation_count() == 3);
    assert(context.relation_rank() == 2);
    return context;
}

void append_integer_relation(silex::ClassGroupContext& context,
                             silex::NumberField& field,
                             slong value) noexcept {
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element generator(field);
    assert(generator.set_si(value));
    assert(relation.set_generator(generator));
    assert(context.append_relation(relation));
}

void assert_finish_workspace_exact(
        const silex::ClassGroupContext& context,
        silex::detail::HnfFinishWorkspace& workspace) noexcept {
    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    assert(workspace.exact_valid);
    assert(workspace.processed_relation_count == relation_count);
    assert(workspace.generator_count == generator_count);
    assert(workspace.relation_rank == context.relation_rank());

    sflint::FmpzMat stored(relation_count, generator_count);
    sflint::FmpzMat transformed(relation_count, generator_count);
    sflint::FmpzMat expected_hnf(relation_count, generator_count);
    assert(context.relations(sflint::FmpzMatRef(stored)));
    sflint::fmpz_mat_mul(
            sflint::FmpzMatRef(transformed),
            sflint::FmpzMatConstRef(workspace.relation_transform),
            sflint::FmpzMatConstRef(stored));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(transformed),
            sflint::FmpzMatConstRef(workspace.hnf)));

    ::fmpz_mat_hnf(expected_hnf.raw(), stored.raw());
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(expected_hnf),
            sflint::FmpzMatConstRef(workspace.hnf)));

    sflint::Fmpz determinant;
    sflint::fmpz_mat_det(
            sflint::FmpzRef(determinant),
            sflint::FmpzMatConstRef(workspace.relation_transform));
    assert(sflint::fmpz_is_pm1(determinant));

    sflint::FmpzMat dependency_coefficients(0, 0);
    assert(silex::detail::hnf_finish_workspace_witness_coefficients(
            dependency_coefficients, workspace, context));
    sflint::FmpzMat zero_rows(
            relation_count - context.relation_rank(), generator_count);
    sflint::fmpz_mat_mul(
            sflint::FmpzMatRef(zero_rows),
            sflint::FmpzMatConstRef(dependency_coefficients),
            sflint::FmpzMatConstRef(stored));
    assert(::fmpz_mat_is_zero(zero_rows.raw()) != 0);
}

int test_witnessed_hnf_basis() {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext context = relation_context(field, order);
    const slong relation_count = context.relation_count();
    const bool had_presentation = context.has_presentation();
    const auto certification = context.certification_status();

    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(silex::detail::class_relation_witnessed_hnf_basis(
            basis, context));
    assert(sflint::fmpz_mat_nrows(basis.rows) == 2);
    assert(sflint::fmpz_mat_ncols(basis.rows) == 2);
    assert(sflint::fmpz_mat_nrows(basis.relation_coefficients) == 2);
    assert(sflint::fmpz_mat_ncols(basis.relation_coefficients) == 3);
    assert(basis.witnesses.size() == 2);

    sflint::FmpzMat stored(3, 2);
    sflint::FmpzMat product(2, 2);
    sflint::FmpzMat expected_hnf(3, 2);
    assert(context.relations(sflint::FmpzMatRef(stored)));
    ::fmpz_mat_hnf(expected_hnf.raw(), stored.raw());
    assert(::fmpz_mat_is_zero_row(expected_hnf.raw(), 2) != 0);
    for (slong i = 0; i < 2; ++i) {
        for (slong j = 0; j < 2; ++j) {
            assert(sflint::fmpz_equal(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(basis.rows), i, j),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(expected_hnf), i, j)));
        }
    }
    sflint::fmpz_mat_mul(
            sflint::FmpzMatRef(product),
            sflint::FmpzMatConstRef(basis.relation_coefficients),
            sflint::FmpzMatConstRef(stored));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(product),
                                  sflint::FmpzMatConstRef(basis.rows)));

    silex::Ideal one(order);
    sflint::FmpzMat row(1, 2);
    assert(one.one());
    for (slong i = 0; i < 2; ++i) {
        for (slong j = 0; j < 2; ++j) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, j),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(basis.rows), i, j));
        }
        assert(silex::detail::verify_class_group_ideal_relation_witness(
                context, one, basis.witnesses[static_cast<std::size_t>(i)],
                sflint::FmpzMatConstRef(row)));
    }

    assert(context.relation_count() == relation_count);
    assert(context.has_presentation() == had_presentation);
    assert(context.certification_status() == certification);
    return 0;
}

int test_zero_row_unit_consumer_shares_transform() {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext context = relation_context(field, order);
    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(silex::detail::class_relation_witnessed_hnf_basis(
            basis, context));

    sflint::FmpzMat zero_coefficients(0, 0);
    assert(silex::detail::hnf_unit_witness_coefficients(
            zero_coefficients, context));
    assert(sflint::fmpz_mat_nrows(zero_coefficients) == 1);
    assert(sflint::fmpz_mat_ncols(zero_coefficients) == 3);

    sflint::FmpzMat full_transform(3, 3);
    for (slong row = 0; row < 2; ++row) {
        for (slong col = 0; col < 3; ++col) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatRef(full_transform), row, col),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(
                                    basis.relation_coefficients),
                            row, col));
        }
    }
    for (slong col = 0; col < 3; ++col) {
        sflint::fmpz_set(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(full_transform), 2, col),
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatConstRef(zero_coefficients), 0, col));
    }

    sflint::FmpzMat stored(3, 2);
    sflint::FmpzMat transformed(3, 2);
    assert(context.relations(sflint::FmpzMatRef(stored)));
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(transformed),
                         sflint::FmpzMatConstRef(full_transform),
                         sflint::FmpzMatConstRef(stored));
    for (slong row = 0; row < 2; ++row) {
        for (slong col = 0; col < 2; ++col) {
            assert(sflint::fmpz_equal(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(transformed), row, col),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(basis.rows), row, col)));
        }
    }
    assert(::fmpz_mat_is_zero_row(transformed.raw(), 2) != 0);

    std::vector<silex::FactoredElement> units;
    assert(silex::detail::hnf_unit_witnesses(units, context));
    assert(units.size() == 1);
    silex::Ideal one(order);
    sflint::FmpzMat zero_row(1, 2);
    assert(one.one());
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(zero_row));
    assert(silex::detail::verify_class_group_ideal_relation_witness(
            context, one, units[0], sflint::FmpzMatConstRef(zero_row)));
    return 0;
}

int test_rank_deficient_failure_is_transactional() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    silex::detail::WitnessedClassRelationHnfBasis output;
    output.rows = sflint::FmpzMat(1, 1);
    output.relation_coefficients = sflint::FmpzMat(1, 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                                sflint::FmpzMatRef(output.rows), 0, 0),
                        71);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatRef(output.relation_coefficients), 0, 0),
            73);
    silex::FactoredElement sentinel(field);
    silex::Element five(field);
    assert(five.set_si(5));
    assert(sentinel.set_element(five));
    output.witnesses.push_back(std::move(sentinel));

    assert(!silex::detail::class_relation_witnessed_hnf_basis(
            output, context));
    assert(sflint::fmpz_equal_si(sflint::fmpz_mat_entry(
                                         sflint::FmpzMatConstRef(output.rows),
                                         0, 0),
                                 71));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(output.relation_coefficients),
                    0, 0),
            73));
    assert(output.witnesses.size() == 1);
    silex::Element expanded(field);
    assert(output.witnesses[0].evaluate(expanded));
    assert(expanded.equal(five));
    return 0;
}

int test_incremental_finish_workspace_lifecycle() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    auto* workspace =
            silex::detail::class_group_finish_workspace(context);
    assert(workspace != nullptr);
    append_integer_relation(context, field, 4);
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert_finish_workspace_exact(context, *workspace);
    assert(workspace->processed_relation_count == 1);
    assert(workspace->relation_rank == 1);

    append_integer_relation(context, field, 9);
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert_finish_workspace_exact(context, *workspace);
    assert(workspace->processed_relation_count == 2);
    assert(workspace->relation_rank == 2);

    workspace->relation_logs = sflint::ArbMat(2, 1);
    sflint::arb_set_si(
            sflint::arb_mat_entry_ref(workspace->relation_logs, 0, 0), 7);
    workspace->relation_log_precision = 160;
    workspace->logged_relation_count = 2;
    workspace->relation_logs_valid = true;
    workspace->unit_logs = sflint::ArbMat(0, 1);
    workspace->unit_log_precision = 160;
    workspace->unit_logs_valid = true;

    append_integer_relation(context, field, 6);
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert_finish_workspace_exact(context, *workspace);
    assert(workspace->relation_logs_valid);
    assert(workspace->relation_log_precision == 160);
    assert(workspace->logged_relation_count == 2);
    assert(sflint::arb_mat_nrows_value(workspace->relation_logs) == 2);
    assert(sflint::arb_contains_si(
            sflint::arb_mat_entry_ref(
                    sflint::ArbMatConstRef(workspace->relation_logs), 0, 0),
            7));
    assert(!workspace->unit_logs_valid);

    append_integer_relation(context, field, 12);
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert_finish_workspace_exact(context, *workspace);

    sflint::FmpzMat transform_before(
            context.relation_count(), context.relation_count());
    sflint::fmpz_mat_set(
            sflint::FmpzMatRef(transform_before),
            sflint::FmpzMatConstRef(workspace->relation_transform));
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(transform_before),
            sflint::FmpzMatConstRef(workspace->relation_transform)));

    const slong processed = workspace->processed_relation_count;
    silex::ClassGroupContext moved(std::move(context));
    assert(silex::detail::class_group_finish_workspace(context) ==
           nullptr);
    workspace = silex::detail::class_group_finish_workspace(moved);
    assert(workspace != nullptr && workspace->exact_valid);
    assert(workspace->processed_relation_count == processed);
    assert_finish_workspace_exact(moved, *workspace);

    silex::ClassGroupContext replacement(order);
    assert(replacement.build_factor_base(sflint::FmpzConstRef(bound)));
    auto* replacement_workspace =
            silex::detail::class_group_finish_workspace(replacement);
    assert(replacement_workspace != nullptr);
    assert(!replacement_workspace->exact_valid);
    auto* populated_workspace = workspace;
    moved.swap(replacement);
    assert(silex::detail::class_group_finish_workspace(moved) ==
           replacement_workspace);
    assert(silex::detail::class_group_finish_workspace(replacement) ==
           populated_workspace);
    assert_finish_workspace_exact(replacement, *populated_workspace);
    moved.swap(replacement);
    workspace = silex::detail::class_group_finish_workspace(moved);
    assert(workspace == populated_workspace);

    moved.clear();
    workspace = silex::detail::class_group_finish_workspace(moved);
    assert(workspace != nullptr);
    assert(!workspace->exact_valid);
    assert(workspace->processed_relation_count == 0);
    assert(sflint::fmpz_mat_nrows(workspace->hnf) == 0);
    assert(sflint::fmpz_mat_nrows(workspace->relation_transform) == 0);
    return 0;
}

int test_finish_workspace_requires_factor_base() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    auto* workspace =
            silex::detail::class_group_finish_workspace(context);
    assert(workspace != nullptr);
    assert(!silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert(!workspace->exact_valid);
    return 0;
}

int test_default_completion_populates_finish_workspace() {
    silex::NumberField field = cubic_discriminant_81_field();
    const silex::Order equation_order =
            silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz factor_base_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_base_bound), maximal_order));
    if (sflint::fmpz_cmp_ui(
                sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            class_group, maximal_order,
            sflint::FmpzConstRef(factor_base_bound), options, 128));

    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(sflint::FmpzConstRef(class_order)));
    assert(units.free_rank() == 2);
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(units.certification_status() ==
           silex::CertificationMode::proven);

    const auto* workspace =
            silex::detail::class_group_finish_workspace(class_group);
    assert(workspace != nullptr);
    assert(workspace->exact_valid);
    assert(workspace->relation_logs_valid);
    assert(workspace->unit_logs_valid);
    assert(workspace->processed_relation_count ==
           class_group.relation_count());
    assert(workspace->logged_relation_count ==
           class_group.relation_count());
    return 0;
}

}  // namespace

int main() {
    assert(test_witnessed_hnf_basis() == 0);
    assert(test_zero_row_unit_consumer_shares_transform() == 0);
    assert(test_rank_deficient_failure_is_transactional() == 0);
    assert(test_incremental_finish_workspace_lifecycle() == 0);
    assert(test_finish_workspace_requires_factor_base() == 0);
    assert(test_default_completion_populates_finish_workspace() == 0);
    return 0;
}
