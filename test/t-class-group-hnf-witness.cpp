#include <silex/archimedean.hpp>
#include <silex/class_group.hpp>
#include <silex/embedding.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order_unit.hpp>
#include <silex/relation.hpp>
#include <silex/zeta.hpp>

#include "class_group/class_group_internal.hpp"
#include "class_group/relation_search_driver_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/relation_unit_internal.hpp"
#include "relation/relation_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"

#include <cstdlib>
#include <cassert>
#include <cstdint>
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

// Q(sqrt 2) with factor base {P_2, 3} (bound 3).  The relation rows of the
// tests below are installed through RelationAccess::set_relation_from_known_row,
// which trusts its row, so only the per-generator ideal proof in
// class_relation_witnessed_hnf_basis can reject a wrong one.  These tests call
// that function directly: build_sunit_class_context needs a proven class-group
// presentation, and appending a relation to a proven context either drops the
// presentation or skips the row, so its class_hnf_basis stage cannot be
// reached with a corrupted row without a new test hook.
struct Sqrt2Fixture {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext context;
    silex::Element root;  // sqrt 2
    silex::Element unit;  // 1 + sqrt 2
};

void sqrt2_fixture(Sqrt2Fixture& out) noexcept {
    out.field = silex::test::quadratic_field(2);
    const silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    assert(out.order.maximal_order(equation));
    out.context = silex::ClassGroupContext(out.order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(out.context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(out.context.factor_base()->length() == 2);
    out.root = silex::Element(out.field);
    out.unit = silex::Element(out.field);
    assert(out.root.gen());
    silex::Element square(out.field);
    assert(square.multiply(out.root, out.root));
    assert(square.equal_si(2));
    assert(out.unit.add_si(out.root, 1));
}

sflint::FmpzMat honest_row(const silex::FactorBase& base,
                           const silex::Element& generator) noexcept {
    silex::Relation relation(base);
    assert(relation.set_generator(generator));
    sflint::FmpzMat row(1, base.length());
    assert(relation.exponents(sflint::FmpzMatRef(row)));
    return row;
}

void append_known_row(silex::ClassGroupContext& context,
                      const silex::Element& generator,
                      sflint::FmpzMatConstRef row) noexcept {
    const silex::FactorBase* base = context.factor_base();
    silex::Relation relation(*base);
    assert(silex::detail::RelationAccess::set_relation_from_known_row(
            relation, *base, generator, row));
    assert(context.append_relation(relation));
}

bool generator_ideal_matches_row(const silex::ClassGroupContext& context,
                                 const silex::Element& generator,
                                 sflint::FmpzMatConstRef row) noexcept {
    silex::FractionalIdeal principal(*context.parent());
    silex::FractionalIdeal expected(*context.parent());
    assert(principal.set_principal(generator));
    assert(silex::detail::factor_base_row_ideal(
            expected, *context.factor_base(), row));
    return principal.equal(expected);
}

// Differential oracle: every published witness generates its HNF row's
// ideal when expanded (the former expansion check), and its factored
// valuations at every factor-base prime equal that row.
void assert_witnesses_match_rows(
        const silex::ClassGroupContext& context,
        const silex::detail::WitnessedClassRelationHnfBasis& basis) noexcept {
    const silex::FactorBase* base = context.factor_base();
    const slong generator_count = base->length();
    silex::Ideal one(*context.parent());
    assert(one.one());
    sflint::FmpzMat row(1, generator_count);
    silex::PrimeIdeal prime(*context.parent());
    for (std::size_t i = 0; i < basis.witnesses.size(); ++i) {
        for (slong k = 0; k < generator_count; ++k) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(row), 0, k),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(basis.rows),
                            static_cast<slong>(i), k));
            slong valuation = 0;
            assert(base->prime(prime, k));
            assert(prime.valuation(valuation, basis.witnesses[i]));
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(row), 0, k),
                    valuation));
        }
        assert(silex::detail::verify_class_group_ideal_relation_witness(
                context, one, basis.witnesses[i],
                sflint::FmpzMatConstRef(row)));
    }
}

// A relation whose stored row is wrong at one factor-base prime must stop
// the witnessed basis, even though C * stored == rows still holds exactly.
int test_corrupted_known_row_rejected() {
    Sqrt2Fixture fixture;
    sqrt2_fixture(fixture);
    const silex::FactorBase& base = *fixture.context.factor_base();
    silex::Element three(fixture.field);
    assert(three.set_si(3));
    sflint::FmpzMat root_row = honest_row(base, fixture.root);
    sflint::FmpzMat bad_row = honest_row(base, three);
    ::fmpz_mat_add(bad_row.raw(), bad_row.raw(), root_row.raw());
    append_known_row(fixture.context, fixture.root,
                     sflint::FmpzMatConstRef(root_row));
    append_known_row(fixture.context, three,
                     sflint::FmpzMatConstRef(bad_row));
    assert(fixture.context.relation_rank() == 2);
    assert(!generator_ideal_matches_row(fixture.context, three,
                                        sflint::FmpzMatConstRef(bad_row)));

    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(!silex::detail::class_relation_witnessed_hnf_basis(
            basis, fixture.context));
    assert(basis.witnesses.empty());
    return 0;
}

// 15 = 3 * 5 has the factor-base valuations of 3, but (5) is a prime outside
// the factor base, so a check of factor-base valuations alone would accept
// the row of 3; the exact ideal proof must not.
int test_extra_prime_outside_factor_base_rejected() {
    Sqrt2Fixture fixture;
    sqrt2_fixture(fixture);
    const silex::FactorBase& base = *fixture.context.factor_base();
    silex::Element three(fixture.field);
    silex::Element fifteen(fixture.field);
    assert(three.set_si(3));
    assert(fifteen.set_si(15));
    sflint::FmpzMat root_row = honest_row(base, fixture.root);
    sflint::FmpzMat three_row = honest_row(base, three);
    silex::PrimeIdeal prime(fixture.order);
    for (slong k = 0; k < base.length(); ++k) {
        slong valuation = 0;
        assert(base.prime(prime, k));
        assert(prime.valuation(valuation, fifteen));
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(three_row), 0,
                                       k),
                valuation));
    }
    append_known_row(fixture.context, fixture.root,
                     sflint::FmpzMatConstRef(root_row));
    append_known_row(fixture.context, fifteen,
                     sflint::FmpzMatConstRef(three_row));
    assert(!generator_ideal_matches_row(fixture.context, fifteen,
                                        sflint::FmpzMatConstRef(three_row)));

    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(!silex::detail::class_relation_witnessed_hnf_basis(
            basis, fixture.context));
    return 0;
}

// A generator multiplied by a unit keeps its ideal, so its true row passes.
int test_unit_multiplied_generator_accepted() {
    Sqrt2Fixture fixture;
    sqrt2_fixture(fixture);
    const silex::FactorBase& base = *fixture.context.factor_base();
    silex::Element three(fixture.field);
    silex::Element cube(fixture.field);
    silex::Element root_unit(fixture.field);
    silex::Element three_unit(fixture.field);
    sflint::Fmpz exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), -3);
    assert(three.set_si(3));
    assert(cube.pow_fmpz(fixture.unit, sflint::FmpzConstRef(exponent)));
    assert(root_unit.multiply(fixture.root, fixture.unit));
    assert(three_unit.multiply(three, cube));
    sflint::FmpzMat root_row = honest_row(base, fixture.root);
    sflint::FmpzMat three_row = honest_row(base, three);
    append_known_row(fixture.context, root_unit,
                     sflint::FmpzMatConstRef(root_row));
    append_known_row(fixture.context, three_unit,
                     sflint::FmpzMatConstRef(three_row));

    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(silex::detail::class_relation_witnessed_hnf_basis(
            basis, fixture.context));
    assert(basis.witnesses.size() == 2);
    assert_witnesses_match_rows(fixture.context, basis);
    return 0;
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

    assert_witnesses_match_rows(context, basis);

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

    // The cached kernel bases are refreshed by a sync that adds relations.
    for (const bool finish : {true, false}) {
        const sflint::FmpzMat* const kernel =
                silex::detail::hnf_finish_workspace_kernel(
                        *workspace, context, finish);
        assert(kernel != nullptr);
        assert(sflint::fmpz_mat_nrows(*kernel) ==
               context.relation_count() - context.relation_rank());
        assert(sflint::fmpz_mat_ncols(*kernel) == context.relation_count());
    }
    assert(workspace->finish_kernel_valid && workspace->witness_kernel_valid);

    append_integer_relation(context, field, 12);
    assert(silex::detail::sync_hnf_finish_workspace(
            *workspace, context));
    assert(!workspace->finish_kernel_valid &&
           !workspace->witness_kernel_valid);
    for (const bool finish : {true, false}) {
        const sflint::FmpzMat* const kernel =
                silex::detail::hnf_finish_workspace_kernel(
                        *workspace, context, finish);
        assert(kernel != nullptr);
        assert(sflint::fmpz_mat_nrows(*kernel) ==
               context.relation_count() - context.relation_rank());
        assert(sflint::fmpz_mat_ncols(*kernel) == context.relation_count());
    }
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

// Two kernel bases span the same lattice (compared by row HNF).
void assert_same_row_lattice(const sflint::FmpzMat& a,
                             const sflint::FmpzMat& b) noexcept {
    assert(sflint::fmpz_mat_nrows(a) == sflint::fmpz_mat_nrows(b));
    assert(sflint::fmpz_mat_ncols(a) == sflint::fmpz_mat_ncols(b));
    sflint::FmpzMat a_hnf(sflint::fmpz_mat_nrows(a), sflint::fmpz_mat_ncols(a));
    sflint::FmpzMat b_hnf(sflint::fmpz_mat_nrows(b), sflint::fmpz_mat_ncols(b));
    ::fmpz_mat_hnf(a_hnf.raw(), a.raw());
    ::fmpz_mat_hnf(b_hnf.raw(), b.raw());
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(a_hnf),
                                  sflint::FmpzMatConstRef(b_hnf)));
}

// Appends the relation of the integer `value` under the search's keep
// policy, which stores a relation even when it leaves the relation lattice
// unchanged, so every call adds a row.
void append_kept_integer_relation(silex::ClassGroupContext& context,
                                  silex::NumberField& field,
                                  slong value) noexcept {
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element generator(field);
    assert(generator.set_si(value));
    assert(relation.set_generator(generator));
    const slong before = context.relation_count();
    bool skipped = true;
    assert(silex::detail::ClassGroupRelationSearchAccess::
                   append_keep_nonduplicate_relation(skipped, context,
                                                     relation));
    assert(!skipped);
    assert(context.relation_count() == before + 1);
}

sflint::FmpzMat copy_fmpz_mat(const sflint::FmpzMat& in) noexcept {
    sflint::FmpzMat out(sflint::fmpz_mat_nrows(in),
                        sflint::fmpz_mat_ncols(in));
    sflint::fmpz_mat_set(sflint::FmpzMatRef(out), sflint::FmpzMatConstRef(in));
    return out;
}

bool same_fmpz_mat(const sflint::FmpzMat& a, const sflint::FmpzMat& b) noexcept {
    return sflint::fmpz_mat_nrows(a) == sflint::fmpz_mat_nrows(b) &&
           sflint::fmpz_mat_ncols(a) == sflint::fmpz_mat_ncols(b) &&
           sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(a),
                                  sflint::FmpzMatConstRef(b));
}

// Independent model of the workspace's finish reduction.  `cached` is the
// reduced extension of the first m relations (0 x 0 for none).  The input
// has the cached rows in columns 0..g+m-1 of rows 0..m-1 and
// [2^c * relation_i | e_i] in row i >= m; the reduction is the first phase
// of fmpz_lll, fmpz_lll_d_with_removal_knapsack, followed by the exact
// kernel checks.  False when the phase or the checks fail.
bool model_finish_reduction(sflint::FmpzMat& rows_out,
                            sflint::FmpzMat& reduced_out,
                            const sflint::FmpzMat& cached,
                            const sflint::FmpzMat& relations,
                            slong rank) noexcept {
    const slong n = sflint::fmpz_mat_nrows(relations);
    const slong g = sflint::fmpz_mat_ncols(relations);
    const slong m = sflint::fmpz_mat_nrows(cached);
    assert(m < n);
    assert(m == 0 || sflint::fmpz_mat_ncols(cached) == g + m);
    sflint::FmpzMat input(n, g + n);
    for (slong i = 0; i < m; ++i) {
        for (slong j = 0; j < g + m; ++j) {
            ::fmpz_set(::fmpz_mat_entry(input.raw(), i, j),
                       ::fmpz_mat_entry(cached.raw(), i, j));
        }
    }
    for (slong i = m; i < n; ++i) {
        for (slong j = 0; j < g; ++j) {
            ::fmpz_mul_2exp(::fmpz_mat_entry(input.raw(), i, j),
                            ::fmpz_mat_entry(relations.raw(), i, j),
                            silex::detail::relation_kernel_lll_scale_bits);
        }
        ::fmpz_one(::fmpz_mat_entry(input.raw(), i, g + i));
    }
    sflint::FmpzLll lll;
    if (::fmpz_lll_d_with_removal_knapsack(input.raw(), nullptr, nullptr,
                                           lll.raw()) == -1 ||
        !silex::detail::checked_relation_kernel_rows(
                rows_out, sflint::FmpzMatConstRef(relations),
                sflint::FmpzMatConstRef(input), rank)) {
        return false;
    }
    reduced_out = std::move(input);
    return true;
}

sflint::FmpzMat context_relations(
        const silex::ClassGroupContext& context) noexcept {
    sflint::FmpzMat relations(context.relation_count(),
                              context.generator_count());
    assert(context.relations(sflint::FmpzMatRef(relations)));
    return relations;
}

// The workspace's finish kernel and kept extension are the model's
// reduction from scratch, which gives hnf_unit_finish_coefficients's rows.
void assert_finish_kernel_from_scratch(
        silex::ClassGroupContext& context,
        silex::detail::HnfFinishWorkspace& workspace) noexcept {
    const sflint::FmpzMat* kernel =
            silex::detail::hnf_finish_workspace_kernel(workspace, context,
                                                       true);
    assert(kernel != nullptr);
    const sflint::FmpzMat relations = context_relations(context);
    sflint::FmpzMat rows(0, 0);
    sflint::FmpzMat reduced(0, 0);
    assert(model_finish_reduction(rows, reduced, sflint::FmpzMat(0, 0),
                                  relations, context.relation_rank()));
    sflint::FmpzMat expected(0, 0);
    assert(silex::detail::hnf_unit_finish_coefficients(expected, context));
    assert(same_fmpz_mat(rows, expected));
    assert(same_fmpz_mat(*kernel, expected));
    assert(workspace.finish_extension_valid);
    assert(workspace.finish_extension_relation_count ==
           context.relation_count());
    assert(same_fmpz_mat(workspace.finish_extension, reduced));
}

// The workspace's finish kernel after a sync that appends relations starts
// its reduction from the reduced extension kept from the previous finish
// kernel; the kernel and the new kept extension are the model's reduction
// of the padded cached rows and the new relations, which here differ from
// the rows of the reduction from scratch but span the same kernel lattice.
// A kept extension that fails the exact checks falls back to the
// reduction from scratch.  A non-incremental sync and a workspace clear()
// drop the extension, and so does the final fallback.
int test_finish_workspace_kernel_extends_cached_reduction() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    auto* workspace = silex::detail::class_group_finish_workspace(context);
    assert(workspace != nullptr);

    for (const slong value : {24, 256, 32, 12}) {
        append_kept_integer_relation(context, field, value);
    }
    assert(context.relation_count() == 4);
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));
    assert(!workspace->finish_extension_valid);
    assert_finish_kernel_from_scratch(context, *workspace);

    // The sync drops the kernel but keeps the extension; the next kernel
    // extends it.
    const sflint::FmpzMat cached = copy_fmpz_mat(workspace->finish_extension);
    for (const slong value : {9, 96, 8}) {
        append_kept_integer_relation(context, field, value);
    }
    assert(context.relation_count() == 7);
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));
    assert(!workspace->finish_kernel_valid);
    assert(workspace->finish_extension_valid);
    assert(workspace->finish_extension_relation_count == 4);
    assert(same_fmpz_mat(workspace->finish_extension, cached));
    const sflint::FmpzMat* kernel =
            silex::detail::hnf_finish_workspace_kernel(*workspace, context,
                                                       true);
    assert(kernel != nullptr);
    {
        const sflint::FmpzMat relations = context_relations(context);
        const slong n = context.relation_count();
        const slong rank = context.relation_rank();
        sflint::FmpzMat rows(0, 0);
        sflint::FmpzMat reduced(0, 0);
        assert(model_finish_reduction(rows, reduced, cached, relations,
                                      rank));
        assert(sflint::fmpz_mat_nrows(rows) == n - rank);
        assert(same_fmpz_mat(*kernel, rows));
        assert(workspace->finish_extension_valid);
        assert(workspace->finish_extension_relation_count == n);
        assert(same_fmpz_mat(workspace->finish_extension, reduced));

        sflint::FmpzMat scratch(0, 0);
        assert(silex::detail::hnf_unit_finish_coefficients(scratch,
                                                           context));
        assert(!same_fmpz_mat(*kernel, scratch));
        assert_same_row_lattice(*kernel, scratch);
    }

    // A kept extension whose relation part no longer matches the
    // relations: the model shows that the extended input fails the exact
    // checks, so the kernel and the new extension are the reduction from
    // scratch.
    append_kept_integer_relation(context, field, 27);
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));
    const slong kept = workspace->finish_extension_relation_count;
    assert(workspace->finish_extension_valid &&
           kept == context.relation_count() - 1);
    for (slong i = 0; i < kept; ++i) {
        fmpz* entry = sflint::fmpz_mat_entry(
                              sflint::FmpzMatRef(workspace->finish_extension),
                              i, 0)
                              .raw();
        ::fmpz_add_ui(entry, entry, 1);
    }
    {
        sflint::FmpzMat rows(0, 0);
        sflint::FmpzMat reduced(0, 0);
        assert(!model_finish_reduction(rows, reduced,
                                       workspace->finish_extension,
                                       context_relations(context),
                                       context.relation_rank()));
    }
    assert_finish_kernel_from_scratch(context, *workspace);

    // A sync that rebuilds the exact state from the relations drops the
    // extension.  Only clear() drops the exact state through the API, and
    // it drops the extension too, so the test clears exact_valid alone to
    // reach the rebuild with an extension still kept.
    assert(workspace->finish_extension_valid);
    workspace->exact_valid = false;
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));
    assert(workspace->exact_valid);
    assert(!workspace->finish_extension_valid);
    assert(workspace->finish_extension_relation_count == 0);
    assert(sflint::fmpz_mat_nrows(workspace->finish_extension) == 0);
    assert_finish_kernel_from_scratch(context, *workspace);

    // clear(), as at the start of relation production, drops the
    // extension; the next sync and finish kernel start from scratch.
    workspace->clear();
    assert(!workspace->finish_extension_valid);
    assert(workspace->finish_extension_relation_count == 0);
    assert(sflint::fmpz_mat_nrows(workspace->finish_extension) == 0);
    append_kept_integer_relation(context, field, 16);
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));
    assert(!workspace->finish_extension_valid);
    assert_finish_kernel_from_scratch(context, *workspace);

    context.clear();
    workspace = silex::detail::class_group_finish_workspace(context);
    assert(workspace != nullptr);
    assert(!workspace->finish_extension_valid);
    assert(workspace->finish_extension_relation_count == 0);
    assert(sflint::fmpz_mat_nrows(workspace->finish_extension) == 0);
    return 0;
}

// The final fallback (here: an empty finish kernel) keeps no extension.
// n - rank never decreases as relations are appended, so an empty kernel
// never follows a kept extension through the API; the test plants one.
int test_finish_workspace_empty_kernel_drops_extension() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 3);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    auto* workspace = silex::detail::class_group_finish_workspace(context);
    assert(workspace != nullptr);
    for (const slong value : {4, 9}) {
        append_kept_integer_relation(context, field, value);
    }
    assert(context.relation_count() == 2 && context.relation_rank() == 2);
    assert(silex::detail::sync_hnf_finish_workspace(*workspace, context));

    const slong g = context.generator_count();
    workspace->finish_extension = sflint::FmpzMat(1, g + 1);
    ::fmpz_one(::fmpz_mat_entry(workspace->finish_extension.raw(), 0, g));
    workspace->finish_extension_relation_count = 1;
    workspace->finish_extension_valid = true;
    const sflint::FmpzMat* kernel =
            silex::detail::hnf_finish_workspace_kernel(*workspace, context,
                                                       true);
    assert(kernel != nullptr);
    assert(sflint::fmpz_mat_nrows(*kernel) == 0);
    assert(sflint::fmpz_mat_ncols(*kernel) == 2);
    assert(!workspace->finish_extension_valid);
    assert(workspace->finish_extension_relation_count == 0);
    assert(sflint::fmpz_mat_nrows(workspace->finish_extension) == 0);
    return 0;
}

// Local saturation compares the kernel relation's factor base with the
// context base once per residue-dlog kernel; a mismatch, an undefined relation,
// or a context without a factor base must still be rejected, and filling the
// relation from a candidate on the context base must keep it matched.
int test_saturation_kernel_relation_base_check() {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext context = relation_context(field, order);
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);

    silex::Relation relation(*base);
    assert(silex::detail::saturation_kernel_relation_base_matches(
            relation, context));
    silex::Relation candidate(*base);
    silex::Element generator(field);
    assert(generator.set_si(6));
    assert(candidate.set_generator(generator));
    assert(relation.set(candidate));
    assert(silex::detail::saturation_kernel_relation_base_matches(
            relation, context));

    silex::ClassGroupContext larger(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 5);
    assert(larger.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(larger.factor_base()->length() == 3);
    silex::Relation mismatched(*larger.factor_base());
    assert(mismatched.is_defined());
    assert(!silex::detail::saturation_kernel_relation_base_matches(
            mismatched, context));
    assert(!silex::detail::saturation_kernel_relation_base_matches(
            relation, larger));

    silex::Relation undefined;
    assert(!silex::detail::saturation_kernel_relation_base_matches(
            undefined, context));

    silex::ClassGroupContext no_base(order);
    assert(!silex::detail::saturation_kernel_relation_base_matches(
            relation, no_base));
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

    silex::detail::WitnessedClassRelationHnfBasis basis;
    assert(silex::detail::class_relation_witnessed_hnf_basis(
            basis, class_group));
    assert_witnesses_match_rows(class_group, basis);
    return 0;
}

// The analytic finish check hands the independent-unit pre-probe's
// unit-log matrix and pivot selection to the regulator product.  The
// product must return exactly what it returns when it rebuilds both, and
// must ignore a probe taken at another precision.
int test_preprobe_handoff_matches_rebuilt_product() {
    silex::NumberField field = cubic_discriminant_81_field();
    const silex::Order equation_order =
            silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));

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
    assert(units.free_rank() == 2);

    using silex::detail::RegulatorPivotOutcome;
    auto* workspace =
            silex::detail::class_group_finish_workspace(class_group);
    assert(workspace != nullptr);
    const slong precision = 160;
    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(precision));
    sflint::Arb analytic_hR;
    assert(silex::zeta_class_regulator_product(
            sflint::ArbRef(analytic_hR), maximal_order, precision));

    slong independent = -1;
    silex::detail::HnfUnitPivotProbe probe;
    assert(silex::detail::hnf_independent_unit_count(
                   independent, maximal_order, class_group, embeddings,
                   workspace, precision, &probe) ==
           RegulatorPivotOutcome::success);
    assert(independent == 2);
    assert(probe.valid && probe.precision == precision);
    assert(static_cast<slong>(probe.selected.size()) == independent);
    assert(probe.relation_count == class_group.relation_count());

    for (const slong product_precision : {precision, 2 * precision}) {
        sflint::Arb rebuilt;
        sflint::Arb handed;
        slong rebuilt_count = -1;
        slong handed_count = -1;
        const RegulatorPivotOutcome rebuilt_outcome =
                silex::detail::hnf_class_regulator_product(
                        sflint::ArbRef(rebuilt), rebuilt_count,
                        maximal_order, class_group, embeddings,
                        sflint::ArbConstRef(analytic_hR), workspace,
                        product_precision);
        const RegulatorPivotOutcome handed_outcome =
                silex::detail::hnf_class_regulator_product(
                        sflint::ArbRef(handed), handed_count,
                        maximal_order, class_group, embeddings,
                        sflint::ArbConstRef(analytic_hR), workspace,
                        product_precision, &probe);
        assert(rebuilt_outcome == RegulatorPivotOutcome::success);
        assert(handed_outcome == rebuilt_outcome);
        assert(rebuilt_count == 2 && handed_count == rebuilt_count);
        assert(::arb_equal(rebuilt.raw(), handed.raw()) != 0);
    }
    return 0;
}

// The analytic finish check keeps its embedding context between checks at
// one precision, rebuilds it when the precision changes, and releases it
// when the check returns finished.
int test_analytic_finish_embedding_reuse() {
    silex::NumberField field = cubic_discriminant_81_field();
    const silex::Order equation_order =
            silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));

    namespace rs = silex::detail::relation_search;
    rs::RelationSearchControlState route_state;
    const slong precision = route_state.analytic_finish_precision;
    silex::Element theta(field);
    assert(theta.gen());
    sflint::ArbVec logs(3);

    silex::EmbeddingContext* first =
            silex::detail::TestHookAccess::analytic_finish_embeddings(route_state, field,
                                                       precision);
    assert(first != nullptr && first->precision() == precision);
    assert(silex::logarithmic_embedding(
            sflint::ArbVecRef(logs), *first, theta,
            silex::LogEmbeddingMode::product, precision));
    std::uint64_t version = 0;
    assert(silex::detail::embedding_log_cache_version(version, *first));

    // Same precision: the context and its log cache are kept.
    silex::EmbeddingContext* again =
            silex::detail::TestHookAccess::analytic_finish_embeddings(route_state, field,
                                                       precision);
    assert(again == first && again->precision() == precision);
    std::uint64_t kept_version = 0;
    assert(silex::detail::embedding_log_cache_version(kept_version, *again));
    assert(kept_version == version);

    // Doubled precision: a fresh context, with no log cache yet.
    silex::EmbeddingContext* doubled =
            silex::detail::TestHookAccess::analytic_finish_embeddings(route_state, field,
                                                       2 * precision);
    assert(doubled != nullptr && doubled->precision() == 2 * precision);
    std::uint64_t rebuilt_version = 0;
    assert(!silex::detail::embedding_log_cache_version(rebuilt_version,
                                                       *doubled));

    // A finished check releases the kept context.
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

    rs::RelationSearchControlState check_state;
    assert(silex::detail::TestHookAccess::analytic_finish_embeddings(
                   check_state, field, check_state.analytic_finish_precision)
           != nullptr);
    assert(check_state.analytic_finish_embeddings.is_defined());
    bool goal_reached = true;
    const silex::detail::ClassGroupRelationOptions relation_options;
    assert(rs::apply_analytic_finish_check(
                   class_group, check_state, goal_reached, maximal_order,
                   relation_options, false, false,
                   silex::detail::class_group_finish_workspace(class_group),
                   nullptr) == rs::AnalyticFinishDecision::finished);
    assert(!check_state.analytic_finish_embeddings.is_defined());
    assert(check_state.analytic_finish_embedding_precision == 0);
    return 0;
}

int test_hnf_row_pivot_column_flags() {
    std::vector<char> covered;
    silex::fmpz_smat::HnfContext module;
    assert(module.reset(4, 7));
    assert(silex::detail::hnf_row_pivot_column_flags(covered, module, 4));
    assert(covered == std::vector<char>({0, 0, 0, 0}));

    // Rows (0, 2, 1, 0) and (0, 4, 0, 3) have HNF leading columns 1 and 2.
    sflint::FmpzMat rows(2, 4);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(rows), 0, 1),
                        2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(rows), 0, 2),
                        1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(rows), 1, 1),
                        4);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(rows), 1, 3),
                        3);
    bool independent = false;
    assert(module.add_fmpz_mat_row(&independent,
                                   sflint::FmpzMatConstRef(rows), 0));
    assert(independent);
    assert(module.add_fmpz_mat_row(&independent,
                                   sflint::FmpzMatConstRef(rows), 1));
    assert(independent);
    assert(silex::detail::hnf_row_pivot_column_flags(covered, module, 4));
    assert(covered == std::vector<char>({0, 1, 1, 0}));

    // A mismatched ambient dimension cannot extract the HNF rows.
    assert(!silex::detail::hnf_row_pivot_column_flags(covered, module, 3));
    return 0;
}

// Integer relation matrix with rank 2 and a 2-dimensional left kernel.
sflint::FmpzMat kernel_test_relations() noexcept {
    sflint::FmpzMat relations(4, 2);
    const slong entries[4][2] = {{2, 0}, {0, 3}, {4, 6}, {6, 9}};
    for (slong i = 0; i < 4; ++i) {
        for (slong j = 0; j < 2; ++j) {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(relations), i,
                                           j),
                    entries[i][j]);
        }
    }
    return relations;
}

sflint::FmpzMat transform_kernel_rows(sflint::FmpzMatConstRef transform,
                                      slong rank) noexcept {
    const slong n = sflint::fmpz_mat_nrows(transform);
    sflint::FmpzMat rows(n - rank, n);
    for (slong i = rank; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            sflint::fmpz_set(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(rows), i - rank,
                                           j),
                    sflint::fmpz_mat_entry(transform, i, j));
        }
    }
    return rows;
}

bool same_row_lattice(sflint::FmpzMatConstRef left,
                      sflint::FmpzMatConstRef right) noexcept {
    sflint::FmpzMat left_hnf(sflint::fmpz_mat_nrows(left),
                             sflint::fmpz_mat_ncols(left));
    sflint::FmpzMat right_hnf(sflint::fmpz_mat_nrows(right),
                              sflint::fmpz_mat_ncols(right));
    ::fmpz_mat_hnf(left_hnf.raw(), left.raw());
    ::fmpz_mat_hnf(right_hnf.raw(), right.raw());
    return sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left_hnf),
                                  sflint::FmpzMatConstRef(right_hnf));
}

bool zero_relation_part(sflint::FmpzMatConstRef extension, slong row,
                        slong generators) noexcept {
    for (slong j = 0; j < generators; ++j) {
        if (!sflint::fmpz_is_zero(sflint::fmpz_mat_entry(extension, row, j))) {
            return false;
        }
    }
    return true;
}

// The LLL kernel basis is accepted only after its two exact checks; a
// candidate failing either check, or a missing or misshapen candidate, gives
// the HNF transform rows.
int test_relation_kernel_basis_checks_and_fallback() {
    using silex::detail::RelationKernelBasisRoute;
    const sflint::FmpzMat relations = kernel_test_relations();
    const slong n = 4;
    const slong g = 2;
    const slong rank = 2;
    sflint::FmpzMat hnf(n, g);
    sflint::FmpzMat transform(n, n);
    ::fmpz_mat_hnf_transform(hnf.raw(), transform.raw(), relations.raw());
    const sflint::FmpzMat fallback =
            transform_kernel_rows(sflint::FmpzMatConstRef(transform), rank);

    // Accepted LLL route: n - rank rows in the left kernel, spanning the
    // same lattice as the HNF transform rows.
    sflint::FmpzMat basis(0, 0);
    RelationKernelBasisRoute route = RelationKernelBasisRoute::hnf_transform;
    assert(silex::detail::relation_kernel_basis(
            basis, &route, sflint::FmpzMatConstRef(relations), rank));
    assert(route == RelationKernelBasisRoute::lll_extended);
    assert(sflint::fmpz_mat_nrows(basis) == n - rank);
    assert(sflint::fmpz_mat_ncols(basis) == n);
    sflint::FmpzMat product(n - rank, g);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(basis),
                         sflint::FmpzMatConstRef(relations));
    assert(::fmpz_mat_is_zero(product.raw()) != 0);
    assert(same_row_lattice(sflint::FmpzMatConstRef(basis),
                            sflint::FmpzMatConstRef(fallback)));

    sflint::FmpzMat extension(0, 0);
    assert(silex::detail::lll_relation_kernel_extension(
            extension, sflint::FmpzMatConstRef(relations)));
    sflint::FmpzMat checked(0, 0);
    assert(silex::detail::checked_relation_kernel_rows(
            checked, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(extension), rank));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(checked),
                                  sflint::FmpzMatConstRef(basis)));

    const auto assert_fallback = [&](const sflint::FmpzMat* candidate) {
        sflint::FmpzMat out(0, 0);
        RelationKernelBasisRoute taken =
                RelationKernelBasisRoute::lll_extended;
        assert(silex::detail::relation_kernel_basis_from_extension(
                out, &taken, sflint::FmpzMatConstRef(relations), rank,
                candidate));
        assert(taken == RelationKernelBasisRoute::hnf_transform);
        assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(out),
                                      sflint::FmpzMatConstRef(fallback)));
    };

    // Count check: an extra row with a zero relation part gives n - rank + 1
    // candidate rows.
    sflint::FmpzMat extra_row(n, g + n);
    sflint::fmpz_mat_set(sflint::FmpzMatRef(extra_row),
                         sflint::FmpzMatConstRef(extension));
    slong nonkernel = -1;
    for (slong i = 0; i < n && nonkernel < 0; ++i) {
        if (!zero_relation_part(sflint::FmpzMatConstRef(extension), i, g)) {
            nonkernel = i;
        }
    }
    assert(nonkernel >= 0);
    for (slong j = 0; j < g; ++j) {
        sflint::fmpz_zero(sflint::fmpz_mat_entry(
                sflint::FmpzMatRef(extra_row), nonkernel, j));
    }
    assert(!silex::detail::checked_relation_kernel_rows(
            checked, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(extra_row), rank));
    assert_fallback(&extra_row);

    // Product check: a row whose relation part is zero but whose
    // coefficients are not in the kernel.
    sflint::FmpzMat outside_kernel(n, g + n);
    sflint::fmpz_mat_set(sflint::FmpzMatRef(outside_kernel),
                         sflint::FmpzMatConstRef(extension));
    slong kernel_row = -1;
    for (slong i = 0; i < n && kernel_row < 0; ++i) {
        if (zero_relation_part(sflint::FmpzMatConstRef(extension), i, g)) {
            kernel_row = i;
        }
    }
    assert(kernel_row >= 0);
    ::fmpz_add_ui(fmpz_mat_entry(outside_kernel.raw(), kernel_row, g),
                  fmpz_mat_entry(outside_kernel.raw(), kernel_row, g), 1);
    assert(!silex::detail::checked_relation_kernel_rows(
            checked, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(outside_kernel), rank));
    assert_fallback(&outside_kernel);

    // Missing count: the rank of a different matrix.
    assert(!silex::detail::checked_relation_kernel_rows(
            checked, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(extension), rank - 1));

    // Too large a rank: the zero-block has n - rank rows too few.
    assert(!silex::detail::checked_relation_kernel_rows(
            checked, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(extension), rank + 1));

    // Empty kernel (rank == n): no LLL, a 0-row basis by the transform route.
    {
        sflint::FmpzMat full(2, 2);
        ::fmpz_mat_one(full.raw());
        sflint::FmpzMat empty(5, 5);
        RelationKernelBasisRoute empty_route =
                RelationKernelBasisRoute::lll_extended;
        assert(silex::detail::relation_kernel_basis(
                empty, &empty_route, sflint::FmpzMatConstRef(full), 2));
        assert(empty_route == RelationKernelBasisRoute::hnf_transform);
        assert(sflint::fmpz_mat_nrows(empty) == 0);
        assert(sflint::fmpz_mat_ncols(empty) == 2);
    }

    // Misshapen or absent candidates.
    const sflint::FmpzMat misshapen(n, g + n - 1);
    assert_fallback(&misshapen);
    assert_fallback(nullptr);

    // The fallback computes the transform itself: its rows are those of
    // the HNF transform above.
    sflint::FmpzMat transform_rows(0, 0);
    assert(silex::detail::hnf_transform_kernel_rows(
            transform_rows, sflint::FmpzMatConstRef(relations), rank));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(transform_rows),
                                  sflint::FmpzMatConstRef(fallback)));

    // An inconsistent rank is a failure, not a fallback: out of range, or
    // (on the fallback) not the number of nonzero HNF rows.
    sflint::FmpzMat out(0, 0);
    assert(!silex::detail::relation_kernel_basis(
            out, nullptr, sflint::FmpzMatConstRef(relations), n + 1));
    assert(!silex::detail::relation_kernel_basis_from_extension(
            out, nullptr, sflint::FmpzMatConstRef(relations), rank - 1,
            nullptr));
    assert(!silex::detail::relation_kernel_basis_from_extension(
            out, nullptr, sflint::FmpzMatConstRef(relations), rank + 1,
            nullptr));
    return 0;
}

// The finish's kernel basis keeps FLINT's double-precision phase only
// after the two exact checks, and otherwise returns relation_kernel_basis's
// rows: when the phase status is -1, and when the checks reject the phase
// output.
int test_finish_relation_kernel_basis_fast_path_and_fallback() {
    using silex::detail::FinishKernelTestFailure;
    using silex::detail::RelationKernelBasisRoute;
    using silex::detail::TestHookAccess;
    const sflint::FmpzMat relations = kernel_test_relations();
    const slong n = 4;
    const slong g = 2;
    const slong rank = 2;

    sflint::FmpzMat shared(0, 0);
    RelationKernelBasisRoute shared_route =
            RelationKernelBasisRoute::hnf_transform;
    assert(silex::detail::relation_kernel_basis(
            shared, &shared_route, sflint::FmpzMatConstRef(relations), rank));
    assert(shared_route == RelationKernelBasisRoute::lll_extended);

    // The phase is kept: on this input FLINT's check accepts it, so the
    // rows are relation_kernel_basis's.
    sflint::FmpzMat basis(0, 0);
    RelationKernelBasisRoute route = RelationKernelBasisRoute::hnf_transform;
    assert(silex::detail::finish_relation_kernel_basis(
            basis, &route, sflint::FmpzMatConstRef(relations), rank));
    assert(route == RelationKernelBasisRoute::lll_double_extended);
    assert(sflint::fmpz_mat_nrows(basis) == n - rank);
    sflint::FmpzMat product(n - rank, g);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(basis),
                         sflint::FmpzMatConstRef(relations));
    assert(::fmpz_mat_is_zero(product.raw()) != 0);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(basis),
                                  sflint::FmpzMatConstRef(shared)));

    // The finish rows span the same lattice as an independent exact
    // oracle, the HNF transform rows (compared by row HNF).
    {
        sflint::FmpzMat oracle(0, 0);
        assert(silex::detail::hnf_transform_kernel_rows(
                oracle, sflint::FmpzMatConstRef(relations), rank));
        sflint::FmpzMat basis_hnf(n - rank, n);
        sflint::FmpzMat oracle_hnf(n - rank, n);
        ::fmpz_mat_hnf(basis_hnf.raw(), basis.raw());
        ::fmpz_mat_hnf(oracle_hnf.raw(), oracle.raw());
        assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(basis_hnf),
                                      sflint::FmpzMatConstRef(oracle_hnf)));
    }

    // Status -1 alone (the phase output left valid and passing both exact
    // checks) must still send the finish to the fmpz_lll route.
    {
        const slong injected_before =
                TestHookAccess::finish_kernel_failures_injected();
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::first_phase_status);
        sflint::FmpzMat out(0, 0);
        RelationKernelBasisRoute taken =
                RelationKernelBasisRoute::lll_double_extended;
        const bool ok = silex::detail::finish_relation_kernel_basis(
                out, &taken, sflint::FmpzMatConstRef(relations), rank);
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::none);
        assert(ok);
        assert(TestHookAccess::finish_kernel_failures_injected() ==
               injected_before + 1);
        assert(taken == RelationKernelBasisRoute::lll_extended);
        assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(out),
                                      sflint::FmpzMatConstRef(shared)));
    }

    // Entries over 250 bits skip the phase: the fmpz_lll route, nothing
    // injected.
    {
        sflint::FmpzMat big(0, 0);
        {
            sflint::FmpzMat scaled = kernel_test_relations();
            for (slong i = 0; i < n; ++i) {
                fmpz* e = sflint::fmpz_mat_entry(sflint::FmpzMatRef(scaled),
                                                 i, 1).raw();
                ::fmpz_mul_2exp(e, e, 300);
            }
            big = std::move(scaled);
        }
        const slong injected_before =
                TestHookAccess::finish_kernel_failures_injected();
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::kernel_rows);
        sflint::FmpzMat out(0, 0);
        RelationKernelBasisRoute taken =
                RelationKernelBasisRoute::lll_double_extended;
        const bool ok = silex::detail::finish_relation_kernel_basis(
                out, &taken, sflint::FmpzMatConstRef(big), rank);
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::none);
        assert(ok);
        assert(TestHookAccess::finish_kernel_failures_injected() ==
               injected_before);
        assert(taken == RelationKernelBasisRoute::lll_extended);
        assert(sflint::fmpz_mat_nrows(out) == n - rank);
    }

    // The kernel_rows injection zeroes the phase output: the fallback must
    // run fmpz_lll on the unchanged input, not on the phase output.
    for (const FinishKernelTestFailure failure :
         {FinishKernelTestFailure::first_phase,
          FinishKernelTestFailure::kernel_rows}) {
        const slong injected_before =
                TestHookAccess::finish_kernel_failures_injected();
        TestHookAccess::set_finish_kernel_failure(failure);
        sflint::FmpzMat out(0, 0);
        RelationKernelBasisRoute taken =
                RelationKernelBasisRoute::lll_double_extended;
        const bool ok = silex::detail::finish_relation_kernel_basis(
                out, &taken, sflint::FmpzMatConstRef(relations), rank);
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::none);
        assert(ok);
        assert(TestHookAccess::finish_kernel_failures_injected() ==
               injected_before + 1);
        assert(taken == RelationKernelBasisRoute::lll_extended);
        assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(out),
                                      sflint::FmpzMatConstRef(shared)));
    }

    // Empty kernel (rank == n): no phase, the transform route.
    {
        sflint::FmpzMat full(2, 2);
        ::fmpz_mat_one(full.raw());
        sflint::FmpzMat empty(5, 5);
        RelationKernelBasisRoute empty_route =
                RelationKernelBasisRoute::lll_extended;
        const slong injected_before =
                TestHookAccess::finish_kernel_failures_injected();
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::first_phase);
        const bool ok = silex::detail::finish_relation_kernel_basis(
                empty, &empty_route, sflint::FmpzMatConstRef(full), 2);
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::none);
        assert(ok);
        assert(TestHookAccess::finish_kernel_failures_injected() ==
               injected_before);
        assert(empty_route == RelationKernelBasisRoute::hnf_transform);
        assert(sflint::fmpz_mat_nrows(empty) == 0);
        assert(sflint::fmpz_mat_ncols(empty) == 2);
    }

    // An inconsistent rank is a failure, as for relation_kernel_basis.
    sflint::FmpzMat out(0, 0);
    assert(!silex::detail::finish_relation_kernel_basis(
            out, nullptr, sflint::FmpzMatConstRef(relations), n + 1));
    return 0;
}

// The HNF finish unit logs are the logs of the units composed from the
// witness coefficient rows, row by row.
int test_unit_logs_share_witness_rows() {
    silex::NumberField field = cubic_discriminant_81_field();
    const silex::Order equation_order =
            silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));

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
    assert(units.free_rank() == 2);

    const slong relation_count = class_group.relation_count();
    const slong rank = class_group.relation_rank();
    sflint::FmpzMat coefficients(0, 0);
    assert(silex::detail::hnf_unit_witness_coefficients(coefficients,
                                                        class_group));
    assert(sflint::fmpz_mat_nrows(coefficients) == relation_count - rank);
    assert(relation_count - rank > 0);

    // The witness rows are the shared kernel basis of the stored relations.
    sflint::FmpzMat stored(relation_count, class_group.generator_count());
    assert(class_group.relations(sflint::FmpzMatRef(stored)));
    assert(::fmpz_mat_rank(stored.raw()) == rank);
    sflint::FmpzMat shared(0, 0);
    silex::detail::RelationKernelBasisRoute route =
            silex::detail::RelationKernelBasisRoute::hnf_transform;
    assert(silex::detail::relation_kernel_basis(
            shared, &route, sflint::FmpzMatConstRef(stored), rank));
    assert(route == silex::detail::RelationKernelBasisRoute::lll_extended);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(shared),
                                  sflint::FmpzMatConstRef(coefficients)));

    const slong precision = 160;
    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(precision));
    const slong places = 3;
    sflint::ArbMat logs(relation_count - rank, places);
    assert(silex::detail::hnf_unit_log_matrix(
            logs, class_group, embeddings, nullptr, precision,
            silex::detail::HnfUnitLogKernel::witness));

    // The finish rows come from FLINT's double-precision phase without its
    // reducedness check; where that check accepts, as here, they are the
    // witness rows, and so are the finish logs.
    sflint::FmpzMat finish_coefficients(0, 0);
    assert(silex::detail::hnf_unit_finish_coefficients(finish_coefficients,
                                                       class_group));
    sflint::FmpzMat finish_rows(0, 0);
    silex::detail::RelationKernelBasisRoute finish_route =
            silex::detail::RelationKernelBasisRoute::hnf_transform;
    assert(silex::detail::finish_relation_kernel_basis(
            finish_rows, &finish_route, sflint::FmpzMatConstRef(stored),
            rank));
    assert(finish_route ==
           silex::detail::RelationKernelBasisRoute::lll_double_extended);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(finish_rows),
                                  sflint::FmpzMatConstRef(coefficients)));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(finish_coefficients),
            sflint::FmpzMatConstRef(coefficients)));
    sflint::ArbMat finish_logs(relation_count - rank, places);
    assert(silex::detail::hnf_unit_log_matrix(
            finish_logs, class_group, embeddings, nullptr, precision,
            silex::detail::HnfUnitLogKernel::finish));
    for (slong unit = 0; unit < relation_count - rank; ++unit) {
        for (slong place = 0; place < places; ++place) {
            assert(::arb_equal(arb_mat_entry(finish_logs.raw(), unit, place),
                               arb_mat_entry(logs.raw(), unit, place)) != 0);
        }
    }

    // The finish workspace takes its finish and witness rows, and so its
    // unit logs, from these reduced kernel bases, not from the unreduced
    // zero-row part of its HNF transform.
    {
        auto* workspace =
                silex::detail::class_group_finish_workspace(class_group);
        assert(workspace != nullptr);
        assert(silex::detail::sync_hnf_finish_workspace(*workspace,
                                                        class_group));
        const sflint::FmpzMat* const finish_kernel =
                silex::detail::hnf_finish_workspace_kernel(
                        *workspace, class_group, true);
        assert(finish_kernel != nullptr);
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(*finish_kernel),
                sflint::FmpzMatConstRef(finish_coefficients)));
        sflint::FmpzMat workspace_witness(0, 0);
        assert(silex::detail::hnf_finish_workspace_witness_coefficients(
                workspace_witness, *workspace, class_group));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(workspace_witness),
                sflint::FmpzMatConstRef(coefficients)));
        const silex::detail::HnfUnitLogKernel kinds[] = {
                silex::detail::HnfUnitLogKernel::finish,
                silex::detail::HnfUnitLogKernel::witness};
        for (const auto kind : kinds) {
            sflint::ArbMat with_workspace(relation_count - rank, places);
            assert(silex::detail::hnf_unit_log_matrix(
                    with_workspace, class_group, embeddings, workspace,
                    precision, kind));
            for (slong unit = 0; unit < relation_count - rank; ++unit) {
                for (slong place = 0; place < places; ++place) {
                    assert(::arb_equal(
                                   arb_mat_entry(with_workspace.raw(), unit,
                                                 place),
                                   arb_mat_entry(logs.raw(), unit, place)) !=
                           0);
                }
            }
        }
    }

    // Only the finish logs run the phase: with the checks forced to reject
    // it, the finish logs come from the fallback (the same rows here) and
    // the witness logs never reach the phase.
    {
        using silex::detail::FinishKernelTestFailure;
        using silex::detail::TestHookAccess;
        const slong injected_before =
                TestHookAccess::finish_kernel_failures_injected();
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::kernel_rows);
        sflint::ArbMat fallback_logs(relation_count - rank, places);
        const bool finish_ok = silex::detail::hnf_unit_log_matrix(
                fallback_logs, class_group, embeddings, nullptr, precision,
                silex::detail::HnfUnitLogKernel::finish);
        const slong injected_finish =
                TestHookAccess::finish_kernel_failures_injected();
        sflint::ArbMat witness_logs(relation_count - rank, places);
        const bool witness_ok = silex::detail::hnf_unit_log_matrix(
                witness_logs, class_group, embeddings, nullptr, precision,
                silex::detail::HnfUnitLogKernel::witness);
        TestHookAccess::set_finish_kernel_failure(
                FinishKernelTestFailure::none);
        assert(finish_ok && witness_ok);
        assert(injected_finish == injected_before + 1);
        assert(TestHookAccess::finish_kernel_failures_injected() ==
               injected_finish);
        for (slong unit = 0; unit < relation_count - rank; ++unit) {
            for (slong place = 0; place < places; ++place) {
                assert(::arb_equal(
                               arb_mat_entry(fallback_logs.raw(), unit,
                                             place),
                               arb_mat_entry(logs.raw(), unit, place)) != 0);
            }
        }
    }

    sflint::ArbVec row(places);
    silex::Element generator(field);
    for (slong unit = 0; unit < relation_count - rank; ++unit) {
        sflint::ArbVec expected(places);
        for (slong relation = 0; relation < relation_count; ++relation) {
            const auto exponent = sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(coefficients), unit, relation);
            if (sflint::fmpz_is_zero(exponent)) {
                continue;
            }
            assert(class_group.relation_generator(generator, relation));
            assert(silex::logarithmic_embedding(
                    sflint::ArbVecRef(row), embeddings, generator,
                    silex::LogEmbeddingMode::product, precision));
            for (slong place = 0; place < places; ++place) {
                ::arb_addmul_fmpz(expected.data() + place, row.data() + place,
                                  exponent.raw(), precision);
            }
        }
        for (slong place = 0; place < places; ++place) {
            assert(::arb_overlaps(expected.data() + place,
                                  arb_mat_entry(logs.raw(), unit, place)) !=
                   0);
        }
    }
    return 0;
}

sflint::FmpzMat integer_rows(slong rows, slong cols,
                             const std::vector<slong>& entries) {
    sflint::FmpzMat out(rows, cols);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            ::fmpz_set_si(::fmpz_mat_entry(out.raw(), i, j),
                          entries[static_cast<std::size_t>(i * cols + j)]);
        }
    }
    return out;
}

int test_unit_candidate_kernel_lll_matches_fmpz_lll() {
    sflint::FmpzLll lll;
    const double delta = lll.raw()->delta;
    const double eta = lll.raw()->eta;

    // Boundary cases of the exact test against FLINT's fmpz_mat_is_reduced:
    // mu = 1/2, mu = 51/100 (below the double 0.51), mu = 52/100, a
    // Lovasz ratio of exactly 99/100, a failing Lovasz pair, and a zero
    // Gram-Schmidt vector.
    struct Case {
        slong rows;
        slong cols;
        std::vector<slong> entries;
        bool reduced;
    };
    const Case cases[] = {
            {2, 2, {2, 0, 1, 2}, true},
            {2, 2, {100, 0, 51, 86}, true},
            {2, 2, {100, 0, 52, 86}, false},
            {2, 4, {10, 0, 0, 0, 0, 3, 3, 9}, true},
            {2, 2, {10, 0, 0, 9}, false},
            {2, 2, {1, 1, 0, 0}, false},
            {3, 3, {0, 0, 0, 1, 0, 0, 0, 1, 0}, false},
    };
    for (const Case& item : cases) {
        const sflint::FmpzMat basis =
                integer_rows(item.rows, item.cols, item.entries);
        assert(silex::lat::detail::integer_basis_is_lll_reduced(
                       sflint::FmpzMatConstRef(basis), delta, eta) ==
               item.reduced);
        assert((::fmpz_mat_is_reduced(basis.raw(), delta, eta) != 0) ==
               item.reduced);
    }

    // Equality cases of `<=` with dyadic parameters (exact in a double):
    // mu = 1/2 = eta, and a Lovasz ratio with equality, both against
    // fmpz_mat_is_reduced.  A `<` for `<=` slip fails the first and third.
    {
        struct DyadicCase {
            slong rows;
            slong cols;
            std::vector<slong> entries;
            double delta;
            double eta;
            bool reduced;
        };
        const DyadicCase dyadic[] = {
                // mu_{1,0} = 1/2 = eta; Lovasz holds strictly.
                {2, 2, {2, 0, 1, 2}, 0.75, 0.5, true},
                // mu = 1/2 > eta = 1/4.
                {2, 2, {2, 0, 1, 2}, 0.75, 0.25, false},
                // Lovasz equality: mu = 1/2, |b*_1|^2 = 1 and
                // (1/2 - 1/4) * 4 = 1.
                {2, 2, {2, 0, 1, 1}, 0.5, 0.5, true},
                // The same basis with delta just above 1/2 fails Lovasz.
                {2, 2, {2, 0, 1, 1}, 0.5078125, 0.5, false},
        };
        for (const DyadicCase& item : dyadic) {
            const sflint::FmpzMat basis =
                    integer_rows(item.rows, item.cols, item.entries);
            const bool flint_says =
                    ::fmpz_mat_is_reduced(basis.raw(), item.delta,
                                          item.eta) != 0;
            assert(flint_says == item.reduced);
            assert(silex::lat::detail::integer_basis_is_lll_reduced(
                           sflint::FmpzMatConstRef(basis), item.delta,
                           item.eta) == item.reduced);
        }
    }

    // The 250-bit guard in reduce_unit_candidate_kernel.  This basis has 254
    // bits; the bare first pass is accepted by the exact test but differs
    // from the basis fmpz_lll returns, so skipping the guard (taking the
    // fast path above 250 bits) returns a different basis.
    {
        sflint::FmpzMat basis(2, 2);
        ::fmpz_set_si(::fmpz_mat_entry(basis.raw(), 0, 0), 1);
        assert(::fmpz_set_str(
                       ::fmpz_mat_entry(basis.raw(), 0, 1),
                       "-374144419153988888124815421131069254946284848021503",
                       10) == 0);
        ::fmpz_set_si(::fmpz_mat_entry(basis.raw(), 1, 0), -1073741823);
        assert(::fmpz_set_str(
                       ::fmpz_mat_entry(basis.raw(), 1, 1),
                       "-28948022201489262187290378635446524851005985908420609"
                       "365455395862898437783552",
                       10) == 0);
        assert(std::labs(static_cast<long>(
                       ::fmpz_mat_max_bits(basis.raw()))) > 250);
        sflint::FmpzMat pass_basis(2, 2);
        ::fmpz_mat_set(pass_basis.raw(), basis.raw());
        sflint::FmpzMat pass_transform(2, 2);
        assert(::fmpz_lll_d_with_removal_knapsack(
                       pass_basis.raw(), pass_transform.raw(), nullptr,
                       lll.raw()) != -1);
        assert(silex::lat::detail::integer_basis_is_lll_reduced(
                sflint::FmpzMatConstRef(pass_basis), delta, eta));
        sflint::FmpzMat expected(2, 2);
        ::fmpz_mat_set(expected.raw(), basis.raw());
        sflint::FmpzMat expected_transform(2, 2);
        fmpz_lll(expected.raw(), expected_transform.raw(), lll.raw());
        assert(::fmpz_mat_equal(pass_basis.raw(), expected.raw()) == 0);
        silex::detail::reduce_unit_candidate_kernel(basis, nullptr);
        assert(::fmpz_mat_equal(basis.raw(), expected.raw()) != 0);
    }

    flint_rand_t state;
    flint_rand_init(state);
    slong reduced_seen = 0;
    slong fallback_seen = 0;
    slong fast_path_seen = 0;
    slong above_cutoff_seen = 0;
    for (slong iter = 0; iter < 3000; ++iter) {
        const slong rows = 2 + static_cast<slong>(n_randint(state, 9));
        const slong cols = rows + static_cast<slong>(n_randint(state, 8));
        flint_bitcnt_t bits = 1 + n_randint(state, 8);
        if (n_randint(state, 20) == 0) {
            bits = 240 + n_randint(state, 40);
        }
        sflint::FmpzMat input(rows, cols);
        ::fmpz_mat_randtest(input.raw(), state, bits);
        if (n_randint(state, 2) == 0) {
            sflint::FmpzMat transform(rows, rows);
            fmpz_lll(input.raw(), transform.raw(), lll.raw());
        }
        if (n_randint(state, 6) == 0) {
            const slong source = static_cast<slong>(n_randint(state, rows));
            const slong target = static_cast<slong>(n_randint(state, rows));
            for (slong j = 0; j < cols; ++j) {
                ::fmpz_add(::fmpz_mat_entry(input.raw(), target, j),
                           ::fmpz_mat_entry(input.raw(), target, j),
                           ::fmpz_mat_entry(input.raw(), source, j));
            }
        }

        const bool exact = silex::lat::detail::integer_basis_is_lll_reduced(
                sflint::FmpzMatConstRef(input), delta, eta);
        assert(exact == (::fmpz_mat_is_reduced(input.raw(), delta, eta) != 0));
        reduced_seen += exact ? 1 : 0;

        sflint::FmpzMat expected(rows, cols);
        ::fmpz_mat_set(expected.raw(), input.raw());
        sflint::FmpzMat expected_transform(rows, rows);
        fmpz_lll(expected.raw(), expected_transform.raw(), lll.raw());

        sflint::FmpzMat actual(rows, cols);
        ::fmpz_mat_set(actual.raw(), input.raw());
        silex::detail::reduce_unit_candidate_kernel(actual, nullptr);
        assert(::fmpz_mat_equal(actual.raw(), expected.raw()) != 0);

        // The same with each failure injected into the certified first
        // phase (CertifiedLllTestSeams): the fallback still returns
        // fmpz_lll's basis.
        for (const auto failure :
             {silex::detail::CertifiedLllTestFailure::first_phase,
              silex::detail::CertifiedLllTestFailure::certificate}) {
            sflint::FmpzMat injected(rows, cols);
            ::fmpz_mat_set(injected.raw(), input.raw());
            silex::detail::TestHookAccess::set_certified_lll_failure(failure);
            silex::detail::reduce_unit_candidate_kernel(injected, nullptr);
            silex::detail::TestHookAccess::set_certified_lll_failure(
                    silex::detail::CertifiedLllTestFailure::none);
            assert(::fmpz_mat_equal(injected.raw(), expected.raw()) != 0);
        }

        // The bare first pass, as reduce_unit_candidate_kernel runs it.  For
        // at most 250 bits an accepted pass is the fast path.  Above 250
        // bits fmpz_lll takes another route (see the fixed case below).
        {
            sflint::FmpzMat pass_basis(rows, cols);
            ::fmpz_mat_set(pass_basis.raw(), input.raw());
            sflint::FmpzMat pass_transform(rows, rows);
            const int pass = ::fmpz_lll_d_with_removal_knapsack(
                    pass_basis.raw(), pass_transform.raw(), nullptr,
                    lll.raw());
            if (pass != -1 &&
                silex::lat::detail::integer_basis_is_lll_reduced(
                        sflint::FmpzMatConstRef(pass_basis), delta, eta)) {
                const slong input_bits = ::fmpz_mat_max_bits(input.raw());
                if ((input_bits < 0 ? -input_bits : input_bits) <= 250) {
                    ++fast_path_seen;
                    // The shared helper certifies this phase output with
                    // fmpz_mat_is_reduced above the 16-bit cutoff.
                    const slong pass_bits =
                            ::fmpz_mat_max_bits(pass_basis.raw());
                    if ((pass_bits < 0 ? -pass_bits : pass_bits) >
                        silex::lat::detail::
                                certified_lll_integer_test_max_bits) {
                        ++above_cutoff_seen;
                    }
                }
            }
        }
        if (!silex::lat::detail::integer_basis_is_lll_reduced(
                    sflint::FmpzMatConstRef(expected), delta, eta)) {
            ++fallback_seen;
        }
    }
    flint_rand_clear(state);
    assert(reduced_seen > 0);
    assert(fallback_seen > 0);
    assert(fast_path_seen > 0);
    assert(above_cutoff_seen > 0);
    return 0;
}

// solve_unit_candidate_system returns what fmpz_mat_solve_dixon_den
// returns, fills the Dixon cache on the first basis, keeps it while the
// basis is unchanged, replaces it for a new basis, reaches the final
// reconstruction past the modulus bound, and reports a singular
// basis without caching it.
bool unit_candidate_system_matches_flint(
        silex::detail::RelationUnitCandidateBatchState& state,
        const sflint::FmpzMat& basis, const sflint::FmpzMat& rhs) {
    const slong n = ::fmpz_mat_nrows(basis.raw());
    const slong cols = ::fmpz_mat_ncols(rhs.raw());
    sflint::FmpzMat expected(n, cols);
    sflint::Fmpz expected_den;
    if (::fmpz_mat_solve_dixon_den(expected.raw(), expected_den.raw(),
                                   basis.raw(), rhs.raw()) == 0) {
        return false;
    }
    sflint::FmpzMat actual(n, cols);
    sflint::Fmpz actual_den;
    if (!silex::detail::solve_unit_candidate_system(
                actual, actual_den, state, sflint::FmpzMatConstRef(basis),
                sflint::FmpzMatConstRef(rhs))) {
        return false;
    }
    return ::fmpz_mat_equal(actual.raw(), expected.raw()) != 0 &&
           ::fmpz_equal(actual_den.raw(), expected_den.raw()) != 0 &&
           state.dixon_prime != 0 &&
           ::fmpz_mat_equal(state.dixon_basis.raw(), basis.raw()) != 0;
}

void fill_random(sflint::FmpzMat& matrix, ::flint_rand_t state,
                 flint_bitcnt_t bits) {
    for (slong i = 0; i < ::fmpz_mat_nrows(matrix.raw()); ++i) {
        for (slong j = 0; j < ::fmpz_mat_ncols(matrix.raw()); ++j) {
            ::fmpz_randtest(::fmpz_mat_entry(matrix.raw(), i, j), state,
                            bits);
        }
    }
}

int test_unit_candidate_system_dixon_cache() {
    ::flint_rand_t random;
    ::flint_rand_init(random);

    // Random nonsingular systems of several sizes; entry sizes up to 80 bits
    // make the lifting run past the first stabilization checks.
    const slong sizes[] = {1, 2, 6, 20, 40};
    const flint_bitcnt_t bits[] = {4, 20, 80};
    for (const slong n : sizes) {
        for (const flint_bitcnt_t b : bits) {
            sflint::FmpzMat basis(n, n);
            do {
                fill_random(basis, random, b);
            } while (::fmpz_mat_rank(basis.raw()) != n);
            silex::detail::RelationUnitCandidateBatchState state;
            assert(state.dixon_prime == 0);
            sflint::FmpzMat rhs(n, 3);
            fill_random(rhs, random, b);
            assert(unit_candidate_system_matches_flint(state, basis, rhs));
            const ulong prime = state.dixon_prime;
            const ulong* const inverse_entries =
                    state.dixon_inverse.raw()->entries;

            // Unchanged basis, new right-hand side: the cache is reused.  A
            // recomputed inverse is a new allocation swapped in while the
            // old one is still live, so an unchanged entries pointer shows
            // the cached inverse was kept.
            sflint::FmpzMat other_rhs(n, 5);
            fill_random(other_rhs, random, b + 7);
            assert(unit_candidate_system_matches_flint(state, basis,
                                                       other_rhs));
            assert(state.dixon_prime == prime);
            assert(state.dixon_inverse.raw()->entries == inverse_entries);

            // A new basis replaces the cached entry.
            sflint::FmpzMat next(n, n);
            do {
                fill_random(next, random, b);
            } while (::fmpz_mat_rank(next.raw()) != n ||
                     ::fmpz_mat_equal(next.raw(), basis.raw()) != 0);
            assert(unit_candidate_system_matches_flint(state, next, rhs));

            // Returning to the first basis still gives FLINT's answer.
            assert(unit_candidate_system_matches_flint(state, basis,
                                                       other_rhs));
        }
    }

    // The final reconstruction past the modulus bound.  For the 1 x 1
    // system a x = 1 with a = 2^300 + 1, fmpz_mat_solve_bound gives
    // D >= |a| and N = D, so the bound is 2 D^2 >= 2 a^2.  Rational
    // reconstruction of 1/a modulo m needs a <= floor(sqrt(m / 2)), that is
    // m >= 2 a^2, and the stabilization checks only run at odd moduli
    // p^i <= 2 D^2, which lie below 2 a^2 unless they fall in the window
    // [2 a^2, 2 D^2] of relative width about 2^-29 left by the det-bound
    // rounding.  So the checks at i = 1, 2, 3, 5, 8 fail and the solution
    // comes from the reconstruction after the loop.
    {
        sflint::FmpzMat basis(1, 1);
        ::fmpz_one(::fmpz_mat_entry(basis.raw(), 0, 0));
        ::fmpz_mul_2exp(::fmpz_mat_entry(basis.raw(), 0, 0),
                        ::fmpz_mat_entry(basis.raw(), 0, 0), 300);
        ::fmpz_add_ui(::fmpz_mat_entry(basis.raw(), 0, 0),
                      ::fmpz_mat_entry(basis.raw(), 0, 0), 1);
        sflint::FmpzMat rhs(1, 2);
        ::fmpz_one(::fmpz_mat_entry(rhs.raw(), 0, 0));
        ::fmpz_set_si(::fmpz_mat_entry(rhs.raw(), 0, 1), -1);
        silex::detail::RelationUnitCandidateBatchState state;
        assert(unit_candidate_system_matches_flint(state, basis, rhs));
    }

    // A singular basis fails and leaves no cache entry behind.
    {
        sflint::FmpzMat basis(3, 3);
        sflint::FmpzMat rhs(3, 1);
        for (slong j = 0; j < 3; ++j) {
            ::fmpz_set_si(::fmpz_mat_entry(basis.raw(), 0, j), j + 1);
            ::fmpz_set_si(::fmpz_mat_entry(basis.raw(), 1, j), 2 * (j + 1));
            ::fmpz_set_si(::fmpz_mat_entry(basis.raw(), 2, j), j * j - 4);
            ::fmpz_set_si(::fmpz_mat_entry(rhs.raw(), j, 0), j);
        }
        silex::detail::RelationUnitCandidateBatchState state;
        sflint::FmpzMat identity(3, 3);
        ::fmpz_mat_one(identity.raw());
        assert(unit_candidate_system_matches_flint(state, identity, rhs));
        sflint::FmpzMat solution(3, 1);
        sflint::Fmpz den;
        assert(!silex::detail::solve_unit_candidate_system(
                solution, den, state, sflint::FmpzMatConstRef(basis),
                sflint::FmpzMatConstRef(rhs)));
        assert(state.dixon_prime == 0);
    }

    ::flint_rand_clear(random);
    return 0;
}

}  // namespace

int main() {
    assert(test_witnessed_hnf_basis() == 0);
    assert(test_zero_row_unit_consumer_shares_transform() == 0);
    assert(test_rank_deficient_failure_is_transactional() == 0);
    assert(test_incremental_finish_workspace_lifecycle() == 0);
    assert(test_finish_workspace_requires_factor_base() == 0);
    assert(test_finish_workspace_kernel_extends_cached_reduction() == 0);
    assert(test_finish_workspace_empty_kernel_drops_extension() == 0);
    assert(test_saturation_kernel_relation_base_check() == 0);
    assert(test_default_completion_populates_finish_workspace() == 0);
    assert(test_corrupted_known_row_rejected() == 0);
    assert(test_extra_prime_outside_factor_base_rejected() == 0);
    assert(test_unit_multiplied_generator_accepted() == 0);
    assert(test_hnf_row_pivot_column_flags() == 0);
    assert(test_preprobe_handoff_matches_rebuilt_product() == 0);
    assert(test_analytic_finish_embedding_reuse() == 0);
    assert(test_relation_kernel_basis_checks_and_fallback() == 0);
    assert(test_finish_relation_kernel_basis_fast_path_and_fallback() == 0);
    assert(test_unit_logs_share_witness_rows() == 0);
    assert(test_unit_candidate_kernel_lll_matches_fmpz_lll() == 0);
    assert(test_unit_candidate_system_dixon_cache() == 0);
    return 0;
}
