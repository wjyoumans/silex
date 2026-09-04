#include "relation_completion_scheduler_internal.hpp"

#include "class_group_internal.hpp"

#include "relation_search_internal.hpp"
#include "../factor_base/factor_base_internal.hpp"

#include <silex/unit.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace silex::detail::relation_search {
namespace {

constexpr slong kMinimumSubfactorBaseSize = 3;
constexpr slong kSubfactorBaseMaxGrowthTrials = 3;
constexpr slong kSubfactorBaseDependencyMultiplier = 16;
constexpr slong kSubfactorBaseRotationDivisor = 10;

bool build_subfactor_base_exclusions(
        std::vector<char>& excluded_indices,
        const FactorBase& base) noexcept {
    excluded_indices.assign(static_cast<std::size_t>(base.length()), 0);
    for (slong block = 0; block < base.rational_prime_block_count(); ++block) {
        flint::Fmpz rational_prime;
        slong length = 0;
        bool complete = false;
        if (!base.rational_prime_block_data(
                    flint::FmpzRef(rational_prime), length, block) ||
            !detail::FactorBaseBlockAccess::rational_prime_block_is_complete(
                    complete, base, block)) {
            return false;
        }
        if (complete && length > 0) {
            slong index = -1;
            if (!base.rational_prime_block_index(index, block, length - 1)) {
                return false;
            }
            excluded_indices[static_cast<std::size_t>(index)] = 1;
        }
    }
    return true;
}

bool build_factor_base_norm_permutation(
        std::vector<slong>& permutation,
        const FactorBase& base) noexcept {
    struct NormIndex {
        slong index = 0;
        flint::Fmpz norm;
    };

    std::vector<NormIndex> entries;
    entries.reserve(static_cast<std::size_t>(base.length()));
    for (slong i = 0; i < base.length(); ++i) {
        const PrimeIdeal* prime = base.prime_at(i);
        NormIndex entry;
        entry.index = i;
        if (prime == nullptr || !prime->norm(flint::FmpzRef(entry.norm))) {
            return false;
        }
        entries.emplace_back(std::move(entry));
    }

    std::sort(entries.begin(), entries.end(),
              [](const NormIndex& left, const NormIndex& right) noexcept {
                  const int cmp = flint::fmpz_cmp(
                          flint::FmpzConstRef(left.norm),
                          flint::FmpzConstRef(right.norm));
                  return cmp < 0 || (cmp == 0 && left.index < right.index);
              });

    permutation.clear();
    permutation.reserve(entries.size());
    for (const NormIndex& entry : entries) {
        permutation.push_back(entry.index);
    }
    return true;
}

void set_subfactor_base_dependency_thresholds(
        SubfactorBaseSchedule& schedule) noexcept {
    const slong length =
            static_cast<slong>(schedule.subfactor_base.size());
    schedule.dependency_growth_threshold =
            length * kSubfactorBaseDependencyMultiplier;
    schedule.dependency_rotation_threshold =
            schedule.dependency_growth_threshold /
            kSubfactorBaseRotationDivisor;
}

}  // namespace

bool build_initial_subfactor_base(
        SubfactorBaseSchedule& schedule,
        const FactorBase& base,
        flint::FmpzConstRef product_bound) noexcept {
    schedule = SubfactorBaseSchedule{};
    std::vector<slong> norm_permutation;
    if (!build_factor_base_norm_permutation(norm_permutation, base) ||
        !build_subfactor_base_exclusions(schedule.excluded_indices, base)) {
        return false;
    }

    const double product_limit = flint::fmpz_get_d(product_bound);
    double product = 1.0;
    std::vector<slong> yes;
    std::vector<slong> no;
    std::size_t pos = 0;
    for (; pos < norm_permutation.size(); ++pos) {
        const slong index = norm_permutation[pos];
        if (index < 0 || index >= base.length()) {
            return false;
        }
        if (schedule.excluded_indices[static_cast<std::size_t>(index)] != 0) {
            no.push_back(index);
            continue;
        }

        const PrimeIdeal* prime = base.prime_at(index);
        flint::Fmpz norm;
        if (prime == nullptr || !prime->norm(flint::FmpzRef(norm))) {
            return false;
        }
        yes.push_back(index);
        product *= flint::fmpz_get_d(flint::FmpzConstRef(norm));
        if (static_cast<slong>(yes.size()) >=
                    kMinimumSubfactorBaseSize &&
            product > product_limit) {
            ++pos;
            break;
        }
    }

    schedule.permutation = yes;
    schedule.permutation.insert(schedule.permutation.end(), no.begin(),
                                no.end());
    for (; pos < norm_permutation.size(); ++pos) {
        schedule.permutation.push_back(norm_permutation[pos]);
    }
    schedule.subfactor_base = std::move(yes);
    set_subfactor_base_dependency_thresholds(schedule);
    return true;
}

bool subfactor_base_index_is_excluded(
        const SubfactorBaseSchedule& schedule,
        slong index) noexcept {
    return index < 0 ||
           index >= static_cast<slong>(schedule.excluded_indices.size()) ||
           schedule.excluded_indices[static_cast<std::size_t>(index)] != 0;
}

bool change_subfactor_base(
        SubfactorBaseSchedule& schedule,
        const FactorBase& base,
        const std::vector<slong>& useful,
        SubfactorBaseChange change) noexcept {
    const slong desired = change == SubfactorBaseChange::increase
            ? static_cast<slong>(schedule.subfactor_base.size()) + 1
            : static_cast<slong>(schedule.subfactor_base.size());
    if (desired <= 0 || base.length() <= 0 || schedule.permutation.empty() ||
        schedule.excluded_indices.size() !=
                static_cast<std::size_t>(base.length())) {
        return false;
    }

    std::vector<slong> next;
    std::vector<char> present(static_cast<std::size_t>(base.length()), 0);
    next.reserve(static_cast<std::size_t>(desired));
    std::size_t pos = 0;
    for (; pos < useful.size(); ++pos) {
        const slong index = useful[pos];
        if (index < 0 || index >= base.length() ||
            subfactor_base_index_is_excluded(schedule, index)) {
            continue;
        }
        next.push_back(index);
        present[static_cast<std::size_t>(index)] = 1;
        if (static_cast<slong>(next.size()) >= desired) {
            ++pos;
            break;
        }
    }

    for (; static_cast<slong>(next.size()) < desired &&
           pos < schedule.permutation.size(); ++pos) {
        const slong index = schedule.permutation[pos];
        if (index < 0 || index >= base.length() ||
            present[static_cast<std::size_t>(index)] != 0 ||
            subfactor_base_index_is_excluded(schedule, index)) {
            continue;
        }
        next.push_back(index);
        present[static_cast<std::size_t>(index)] = 1;
    }
    if (static_cast<slong>(next.size()) < desired) {
        return false;
    }

    if (next != schedule.subfactor_base) {
        schedule.subfactor_base = std::move(next);
        set_subfactor_base_dependency_thresholds(schedule);
    }
    return true;
}

SubfactorBaseAdvanceDecision advance_subfactor_base_schedule(
        SubfactorBaseSchedule& schedule,
        RelationCompletionState& completion,
        const FactorBase& base,
        const std::vector<slong>& useful,
        bool guarded_rebuild_available,
        bool rebuild_to_limit_available) noexcept {
    if (++completion.dependent_trials >
        schedule.dependency_growth_threshold) {
        if (++completion.subfactor_base_trials >
                    kSubfactorBaseMaxGrowthTrials &&
            guarded_rebuild_available) {
            return SubfactorBaseAdvanceDecision::
                    request_guarded_factor_base_rebuild;
        }
        if (!change_subfactor_base(schedule, base, useful,
                                   SubfactorBaseChange::increase)) {
            return rebuild_to_limit_available
                    ? SubfactorBaseAdvanceDecision::
                              request_factor_base_rebuild_to_limit
                    : SubfactorBaseAdvanceDecision::increase_failed;
        }
        completion.dependent_trials = 0;
    } else if (schedule.dependency_rotation_threshold > 0 &&
               completion.dependent_trials %
                               schedule.dependency_rotation_threshold ==
                       0) {
        if (!change_subfactor_base(schedule, base, useful,
                                   SubfactorBaseChange::rotate)) {
            return rebuild_to_limit_available
                    ? SubfactorBaseAdvanceDecision::
                              request_factor_base_rebuild_to_limit
                    : SubfactorBaseAdvanceDecision::rotation_failed;
        }
    }
    return SubfactorBaseAdvanceDecision::ready;
}

bool compute_relation_need(
        RelationNeedDecision& decision,
        ClassGroupContext& context,
        const FactorBase& base,
        slong target_relation_kernel_units) noexcept {
    decision = RelationNeedDecision{};

    decision.class_relation_defect =
            max_slong_value(0, base.length() - context.relation_rank());
    decision.need = decision.class_relation_defect;
    if (decision.class_relation_defect > 0) {
        decision.reason = RelationNeedReason::class_relation_defect;
    }

    if (decision.class_relation_defect == 0) {
        if (!context.publish_presentation()) {
            return false;
        }

        decision.relation_kernel_unit_defect = max_slong_value(
                0,
                target_relation_kernel_units -
                        context.relation_kernel_unit_count());
        if (decision.relation_kernel_unit_defect > 0) {
            decision.need = min_slong_value(
                    base.length(),
                    decision.need + decision.relation_kernel_unit_defect);
            decision.reason =
                    RelationNeedReason::relation_kernel_unit_defect;
        }
    }

    decision.goal_reached = decision.need == 0 &&
            context.has_presentation() &&
            context.relation_kernel_unit_count() >=
                    target_relation_kernel_units;
    return true;
}

bool relation_kernel_unit_target(
        slong& target,
        const Order& order,
        const ClassGroupRelationOptions& options) noexcept {
    if (order.parent() == nullptr || options.target_relation_kernel_units < 0) {
        return false;
    }

    slong rank = 0;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }
    target = max_slong_value(options.target_relation_kernel_units, rank);
    return true;
}

void update_relation_need(RelationCompletionState& completion,
                          slong next_need) noexcept {
    completion.need = next_need;
    if (completion.need != completion.old_need) {
        completion.dependent_trials = 0;
        completion.old_need = completion.need;
    }
}

slong relation_collection_target_count(
        const ClassGroupContext& context,
        const ClassGroupRelationOptions& options,
        const RelationCompletionState& completion) noexcept {
    if (completion.need <= 0) {
        return context.relation_count();
    }
    const slong remaining = options.max_relations - context.relation_count();
    if (remaining <= 0) {
        return context.relation_count();
    }
    const slong need = completion.need < remaining
            ? completion.need
            : remaining;
    return context.relation_count() + need;
}

}  // namespace silex::detail::relation_search
