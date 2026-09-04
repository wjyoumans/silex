#include "relation_factor_base_rebuild_internal.hpp"

#include "relation_admission_cache_internal.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace silex::detail {

static_assert(std::is_nothrow_swappable_v<RelationAdmissionCache>);
static_assert(noexcept(std::declval<ClassGroupContext&>().swap(
        std::declval<ClassGroupContext&>())));

bool rebuild_relation_factor_base_and_replay(
        ClassGroupContext& context,
        RelationAdmissionCache& cache,
        const Order& order,
        flint::FmpzConstRef next_bound,
        slong add_need) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr || add_need < 0) {
        return false;
    }

    std::vector<Element> computed_generators;
    computed_generators.reserve(
            static_cast<std::size_t>(context.relation_count()));
    for (slong i = 0; i < context.relation_count(); ++i) {
        Element generator(*field);
        if (!generator.is_defined() ||
            !context.relation_generator(generator, i)) {
            return false;
        }
        computed_generators.emplace_back(std::move(generator));
    }

    ClassGroupContext candidate(order);
    candidate.set_diagnostics(context.diagnostics());
    const slong relation_kernel_units_target =
            ClassGroupRelationSearchAccess::relation_kernel_units_target(
                    context);
    ClassGroupRelationSearchAccess::set_relation_kernel_units_target(
            candidate, relation_kernel_units_target);
    if (!candidate.is_defined() ||
        !ClassGroupRelationSearchAccess::build_relation_factor_base(
                candidate, next_bound)) {
        return false;
    }
    ClassGroupRelationSearchAccess::set_relation_kernel_units_target(
            candidate, relation_kernel_units_target);

    RelationAdmissionCache next_cache;
    if (!initialize_relation_admission_cache(
                candidate, next_cache, order, add_need)) {
        return false;
    }

    const FactorBase* next_base = candidate.factor_base();
    if (next_base == nullptr) {
        return false;
    }
    for (const Element& generator : computed_generators) {
        Relation relation(*next_base);
        bool retained = false;
        if (!relation.is_defined() ||
            !relation.set_generator(generator, context.diagnostics()) ||
            !try_admit_relation(
                    candidate, retained, next_cache, relation, false)) {
            return false;
        }
    }
    next_cache.target_relation_count = candidate.relation_count();

    using std::swap;
    swap(cache, next_cache);
    context.swap(candidate);
    return true;
}

}  // namespace silex::detail
