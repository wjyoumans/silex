#pragma once

#include "class_group_storage_internal.hpp"

#include <vector>

namespace silex::detail::relation_search {

struct NormPrefilter;

inline constexpr slong kAnalyticFinishInitialPrecision = 160;

using RelationSearchControlState = detail::RelationSearchControlState;

enum class AnalyticFinishDecision {
    incomplete,
    finished,
    old_cache_needs_more_relations,
    analytic_needs_more_relations,
    restart_factor_base,
    failed
};

void activate_finish_unit_log_rotation(
        RelationSearchControlState& route_state) noexcept;

bool build_useful_pivot_indices(
        std::vector<slong>& useful,
        ClassGroupContext& context,
        const FactorBase& base,
        const std::vector<slong>& permutation,
        slong need,
        bool finish_unit_log_rotation_active,
        slong finish_unit_log_rotation) noexcept;

AnalyticFinishDecision apply_analytic_finish_check(
        ClassGroupContext& context,
        RelationSearchControlState& route_state,
        bool& goal_reached,
        const Order& order,
        const ClassGroupRelationOptions& options,
        bool allow_rank_zero_regulator_reconstruction,
        bool factor_base_restart_available,
        HnfFinishWorkspace* finish_workspace,
        const DiagnosticsContext* diagnostics) noexcept;

AnalyticFinishDecision apply_honesty_check(
        ClassGroupContext& context,
        RelationSearchControlState& route_state,
        bool goal_reached,
        const Order& order,
        const SubfactorBaseSchedule* subfactor_base_schedule,
        flint::FmpzConstRef active_factor_base_bound,
        bool factor_base_restart_available,
        bool use_required_prime_factor_over_base,
        const DiagnosticsContext* diagnostics) noexcept;

bool collect_small_norm_relations(
        ClassGroupContext& context,
        const Order& order,
        NormPrefilter* norm_prefilter,
        slong target_relation_count,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        detail::RelationAdmissionCache& cache,
        const SubfactorBaseSchedule& subfactor_base_schedule,
        RelationSearchControlState& route_state,
        slong j0,
        bool apply_multiplier_filter,
        std::vector<slong>& small_multiplier,
        bool allow_silex_ideal_lattice_input) noexcept;

bool collect_random_relations(
        ClassGroupContext& context,
        const Order& order,
        NormPrefilter* norm_prefilter,
        slong target_relation_count,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        detail::RelationAdmissionCache& cache,
        SubfactorBaseSchedule& subfactor_base_schedule,
        RelationSearchControlState& route_state,
        bool factor_base_restart_available,
        bool factor_base_restart_to_max_available,
        const DiagnosticsContext* diagnostics) noexcept;

}  // namespace silex::detail::relation_search
