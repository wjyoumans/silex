#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/relation.hpp>

#include "class_group/class_group_storage_internal.hpp"
#include "class_group/relation_admission_cache_internal.hpp"
#include "class_group/relation_factor_base_rebuild_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;
using RelationSearchAccess =
        silex::detail::ClassGroupRelationSearchAccess;

struct QuadraticFixture {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::DiagnosticsContext diagnostics;
    silex::ClassGroupContext context{order};

    explicit QuadraticFixture(ulong bound) noexcept {
        silex::diagnostics_context_init(diagnostics);
        context.set_diagnostics(&diagnostics);
        sflint::Fmpz bound_value;
        sflint::fmpz_set_ui(sflint::FmpzRef(bound_value), bound);
        assert(order.is_maximal());
        assert(RelationSearchAccess::build_relation_factor_base(
                context, sflint::FmpzConstRef(bound_value)));
    }
};

bool admit_generator(silex::ClassGroupContext& context,
                     silex::detail::RelationAdmissionCache& cache,
                     const silex::NumberField& field,
                     slong value) noexcept {
    const silex::FactorBase* base = context.factor_base();
    if (base == nullptr) {
        return false;
    }
    silex::Element generator(field);
    silex::Relation relation(*base);
    bool retained = false;
    return generator.set_si(value) && relation.set_generator(generator) &&
           silex::detail::try_admit_relation(
                   context, retained, cache, relation, false) &&
           retained;
}

bool append_generator(silex::ClassGroupContext& context,
                      const silex::NumberField& field,
                      slong value) noexcept {
    const silex::FactorBase* base = context.factor_base();
    if (base == nullptr) {
        return false;
    }
    silex::Element generator(field);
    silex::Relation relation(*base);
    return generator.set_si(value) && relation.set_generator(generator) &&
           context.append_relation(relation);
}

bool cache_equal(const silex::detail::RelationAdmissionCache& left,
                 const silex::detail::RelationAdmissionCache& right) noexcept {
    return left.target_relation_count == right.target_relation_count &&
           left.relation_surplus == right.relation_surplus &&
           left.relation_count_before_init ==
                   right.relation_count_before_init &&
           left.retained_relation_count == right.retained_relation_count &&
           left.trivial_relation_count == right.trivial_relation_count &&
           left.duplicate_relation_count == right.duplicate_relation_count &&
           left.skipped_dependent_relation_count ==
                   right.skipped_dependent_relation_count &&
           left.missing == right.missing &&
           left.modular_basis == right.modular_basis &&
           left.relation_row_hashes == right.relation_row_hashes;
}

bool context_generators_equal(
        const silex::ClassGroupContext& context,
        const std::vector<silex::Element>& expected) noexcept {
    if (context.relation_count() != static_cast<slong>(expected.size())) {
        return false;
    }
    const silex::Order* order = context.parent();
    if (order == nullptr || order->parent() == nullptr) {
        return false;
    }
    silex::Element actual(*order->parent());
    for (slong i = 0; i < context.relation_count(); ++i) {
        if (!context.relation_generator(actual, i) ||
            !actual.equal(expected[static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    return true;
}

std::vector<silex::Element> copy_context_generators(
        const silex::ClassGroupContext& context) {
    std::vector<silex::Element> generators;
    const silex::Order* order = context.parent();
    assert(order != nullptr);
    assert(order->parent() != nullptr);
    generators.reserve(static_cast<std::size_t>(context.relation_count()));
    for (slong i = 0; i < context.relation_count(); ++i) {
        silex::Element generator(*order->parent());
        assert(context.relation_generator(generator, i));
        generators.emplace_back(std::move(generator));
    }
    return generators;
}

bool matrix_entry_is(const sflint::FmpzMat& matrix,
                     slong row,
                     slong column,
                     slong expected) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(matrix), row, column),
            expected);
}

int test_success_replays_generators_in_stored_order() {
    QuadraticFixture fixture(3);
    assert(fixture.context.generator_count() == 1);

    silex::detail::RelationAdmissionCache cache;
    assert(silex::detail::initialize_relation_admission_cache(
            fixture.context, cache, fixture.order, 2));
    assert(admit_generator(fixture.context, cache, fixture.field, 8));
    assert(admit_generator(fixture.context, cache, fixture.field, 4));
    RelationSearchAccess::set_relation_kernel_units_target(
            fixture.context, 6);

    sflint::Fmpz next_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(next_bound), 7);
    assert(silex::detail::rebuild_relation_factor_base_and_replay(
            fixture.context, cache, fixture.order,
            sflint::FmpzConstRef(next_bound), 2));

    assert(fixture.context.diagnostics() == &fixture.diagnostics);
    assert(RelationSearchAccess::relation_kernel_units_target(
                   fixture.context) == 6);
    assert(fixture.context.generator_count() == 3);
    assert(fixture.context.relation_count() == 4);
    assert(fixture.context.relation_rank() == 2);
    assert(!fixture.context.has_presentation());
    assert(fixture.context.relation_source_count(
                   silex::ClassGroupRelationSource::Search) == 4);

    const slong expected_generators[] = {2, 7, 8, 4};
    silex::Element generator(fixture.field);
    for (slong i = 0; i < 4; ++i) {
        assert(fixture.context.relation_generator(generator, i));
        assert(generator.equal_si(expected_generators[i]));
    }

    sflint::FmpzMat rows(4, 3);
    assert(fixture.context.relations(sflint::FmpzMatRef(rows)));
    const slong expected_rows[4][3] = {
            {2, 0, 0},
            {0, 1, 1},
            {6, 0, 0},
            {4, 0, 0},
    };
    for (slong i = 0; i < 4; ++i) {
        for (slong j = 0; j < 3; ++j) {
            assert(matrix_entry_is(rows, i, j, expected_rows[i][j]));
        }
    }

    assert(cache.target_relation_count == 4);
    assert(cache.relation_surplus == 0);
    assert(cache.relation_count_before_init == 0);
    assert(cache.retained_relation_count == 4);
    assert(cache.trivial_relation_count == 2);
    assert(cache.duplicate_relation_count == 1);
    assert(cache.skipped_dependent_relation_count == 0);
    assert(cache.missing == 1);
    assert(cache.modular_basis.size() == 9);
    assert(cache.relation_row_hashes.size() == 4);
    return 0;
}

int test_failed_replay_preserves_context_and_cache() {
    QuadraticFixture fixture(7);
    assert(append_generator(fixture.context, fixture.field, 4));
    assert(append_generator(fixture.context, fixture.field, 7));
    RelationSearchAccess::set_relation_kernel_units_target(
            fixture.context, 6);

    silex::detail::RelationAdmissionCache cache;
    cache.target_relation_count = 91;
    cache.relation_surplus = 8;
    cache.relation_count_before_init = 7;
    cache.retained_relation_count = 6;
    cache.trivial_relation_count = 5;
    cache.duplicate_relation_count = 4;
    cache.skipped_dependent_relation_count = 3;
    cache.missing = 2;
    cache.modular_basis = {1, 2, 3, 4};
    cache.relation_row_hashes = {5, 6, 7};
    const silex::detail::RelationAdmissionCache cache_before = cache;

    const silex::FactorBase* live_base = fixture.context.factor_base();
    assert(live_base != nullptr);
    silex::FactorBase base_before(fixture.order);
    assert(base_before.set(*live_base));
    sflint::FmpzMat rows_before(
            fixture.context.relation_count(),
            fixture.context.generator_count());
    assert(fixture.context.relations(sflint::FmpzMatRef(rows_before)));
    const std::vector<silex::Element> generators_before =
            copy_context_generators(fixture.context);
    std::vector<silex::ClassGroupRelationSource> sources_before;
    for (slong i = 0; i < fixture.context.relation_count(); ++i) {
        silex::ClassGroupRelationSource source =
                silex::ClassGroupRelationSource::Search;
        assert(fixture.context.relation_source(source, i));
        sources_before.push_back(source);
    }
    const slong relation_count_before = fixture.context.relation_count();
    const slong relation_rank_before = fixture.context.relation_rank();
    const slong skipped_before =
            fixture.context.skipped_dependent_relation_count();
    const bool presentation_before = fixture.context.has_presentation();
    const silex::CertificationMode certification_before =
            fixture.context.certification_status();

    // The smaller candidate initializes successfully and accepts 4, but 7 is
    // not smooth over its relation base. This fails inside ordered replay.
    sflint::Fmpz next_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(next_bound), 3);
    assert(!silex::detail::rebuild_relation_factor_base_and_replay(
            fixture.context, cache, fixture.order,
            sflint::FmpzConstRef(next_bound), 0));

    assert(cache_equal(cache, cache_before));
    assert(fixture.context.diagnostics() == &fixture.diagnostics);
    assert(fixture.context.factor_base() != nullptr);
    assert(fixture.context.factor_base()->equal(base_before));
    assert(fixture.context.relation_count() == relation_count_before);
    assert(fixture.context.relation_rank() == relation_rank_before);
    assert(fixture.context.skipped_dependent_relation_count() ==
           skipped_before);
    assert(fixture.context.has_presentation() == presentation_before);
    assert(fixture.context.certification_status() == certification_before);
    assert(RelationSearchAccess::relation_kernel_units_target(
                   fixture.context) == 6);

    sflint::FmpzMat rows_after(
            fixture.context.relation_count(),
            fixture.context.generator_count());
    assert(fixture.context.relations(sflint::FmpzMatRef(rows_after)));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(rows_after),
            sflint::FmpzMatConstRef(rows_before)));
    assert(context_generators_equal(fixture.context, generators_before));
    for (slong i = 0; i < fixture.context.relation_count(); ++i) {
        silex::ClassGroupRelationSource source =
                silex::ClassGroupRelationSource::Search;
        assert(fixture.context.relation_source(source, i));
        assert(source == sources_before[static_cast<std::size_t>(i)]);
    }
    return 0;
}

}  // namespace

int main() {
    assert(test_success_replays_generators_in_stored_order() == 0);
    assert(test_failed_replay_preserves_context_and_cache() == 0);
    return 0;
}
