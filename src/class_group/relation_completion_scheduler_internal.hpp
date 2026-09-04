#pragma once

#include <silex/class_group.hpp>

#include <vector>

namespace silex::detail {

struct SubfactorBaseSchedule {
    std::vector<slong> permutation;
    std::vector<slong> subfactor_base;
    std::vector<char> excluded_indices;
    slong dependency_growth_threshold = 0;
    slong dependency_rotation_threshold = 0;
};

struct RelationCompletionState {
    slong need = 0;
    slong old_need = 0;
    slong dependent_trials = 0;
    slong subfactor_base_trials = 0;
};

namespace relation_search {

inline constexpr slong kRelationSurplus = 5;

enum class RelationNeedReason {
    none,
    class_relation_defect,
    relation_kernel_unit_defect
};

struct RelationNeedDecision {
    slong need = 0;
    slong class_relation_defect = 0;
    slong relation_kernel_unit_defect = 0;
    bool goal_reached = false;
    RelationNeedReason reason = RelationNeedReason::none;
};

enum class SubfactorBaseChange {
    rotate,
    increase
};

enum class SubfactorBaseAdvanceDecision {
    ready,
    request_guarded_factor_base_rebuild,
    request_factor_base_rebuild_to_limit,
    increase_failed,
    rotation_failed
};

bool build_initial_subfactor_base(
        SubfactorBaseSchedule& schedule,
        const FactorBase& base,
        flint::FmpzConstRef product_bound) noexcept;

bool subfactor_base_index_is_excluded(
        const SubfactorBaseSchedule& schedule,
        slong index) noexcept;

bool change_subfactor_base(
        SubfactorBaseSchedule& schedule,
        const FactorBase& base,
        const std::vector<slong>& useful,
        SubfactorBaseChange change) noexcept;

SubfactorBaseAdvanceDecision advance_subfactor_base_schedule(
        SubfactorBaseSchedule& schedule,
        RelationCompletionState& completion,
        const FactorBase& base,
        const std::vector<slong>& useful,
        bool guarded_rebuild_available,
        bool rebuild_to_limit_available) noexcept;

bool compute_relation_need(
        RelationNeedDecision& decision,
        ClassGroupContext& context,
        const FactorBase& base,
        slong target_relation_kernel_units) noexcept;

bool relation_kernel_unit_target(
        slong& target,
        const Order& order,
        const ClassGroupRelationOptions& options) noexcept;

void update_relation_need(RelationCompletionState& completion,
                          slong next_need) noexcept;

slong relation_collection_target_count(
        const ClassGroupContext& context,
        const ClassGroupRelationOptions& options,
        const RelationCompletionState& completion) noexcept;

}  // namespace relation_search
}  // namespace silex::detail
