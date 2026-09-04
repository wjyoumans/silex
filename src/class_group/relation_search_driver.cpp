#include <silex/class_group.hpp>

#include <silex/embedding.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/lat.hpp>
#include <silex/order_element.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include "class_group_internal.hpp"
#include "class_group_certification_internal.hpp"
#include "factor_base_honesty_internal.hpp"
#include "ideal_t2_enumeration_internal.hpp"
#include "relation_admission_cache_internal.hpp"
#include "relation_search_driver_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"
#include "relation_factor_base_rebuild_internal.hpp"
#include "relation_factor_base_plan_internal.hpp"
#include "relation_search_internal.hpp"
#include "../ideal_factorization/ideal_factorization_internal.hpp"
#include "../order/order_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"
#include "../order_unit/order_unit_internal.hpp"
#include "../zeta/zeta_internal.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

namespace silex {

namespace detail {

class ClassGroupFinishAccess {
public:
    static bool relation_rows(flint::FmpzMatRef out,
                              const ClassGroupContext& context,
                              slong first_relation) noexcept {
        if (first_relation < 0 ||
            first_relation > context.relations_.length() ||
            flint::fmpz_mat_nrows(out) !=
                    context.relations_.length() - first_relation ||
            flint::fmpz_mat_ncols(out) != context.relations_.ncols()) {
            return false;
        }
        for (slong row = 0; row < flint::fmpz_mat_nrows(out); ++row) {
            flint::FmpzMatWindow target(
                    out, row, 0, row + 1,
                    flint::fmpz_mat_ncols(out));
            if (!context.relations_.row(target.ref(),
                                        first_relation + row)) {
                return false;
            }
        }
        return true;
    }

    static HnfFinishWorkspace* finish_workspace(
            ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_->hnf_finish_workspace;
    }

    static const HnfFinishWorkspace* finish_workspace(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_->hnf_finish_workspace;
    }

    static bool finish_torsion(
            const flint::Fmpz*& torsion_order,
            const Element*& torsion_generator,
            const ClassGroupContext& context) noexcept {
        torsion_order = nullptr;
        torsion_generator = nullptr;
        const Order* const order = context.parent();
        if (context.private_storage_ == nullptr || order == nullptr ||
            order->parent() == nullptr || !order->is_maximal() ||
            !context.private_storage_->finish_state.has_value()) {
            return false;
        }
        const RelationSearchControlState& route_state =
                context.private_storage_->finish_state
                        ->route_state;
        if (!route_state.analytic_finish_torsion_valid ||
            flint::fmpz_sgn(flint::FmpzConstRef(
                    route_state.analytic_finish_torsion_order)) <= 0 ||
            !route_state.analytic_finish_torsion_generator.has_parent(
                    *order->parent())) {
            return false;
        }
        torsion_order = &route_state.analytic_finish_torsion_order;
        torsion_generator =
                &route_state.analytic_finish_torsion_generator;
        return true;
    }

};

HnfFinishWorkspace* class_group_finish_workspace(
        ClassGroupContext& context) noexcept {
    return ClassGroupFinishAccess::finish_workspace(context);
}

const HnfFinishWorkspace* class_group_finish_workspace(
        const ClassGroupContext& context) noexcept {
    return ClassGroupFinishAccess::finish_workspace(context);
}

bool class_group_finish_torsion(
        const flint::Fmpz*& torsion_order,
        const Element*& torsion_generator,
        const ClassGroupContext& context) noexcept {
    return ClassGroupFinishAccess::finish_torsion(
            torsion_order, torsion_generator, context);
}

bool class_group_relation_rows(
        flint::FmpzMatRef out,
        const ClassGroupContext& context,
        slong first_relation) noexcept {
    return ClassGroupFinishAccess::relation_rows(out, context, first_relation);
}

}  // namespace detail

namespace detail::relation_search {

constexpr ulong kHonestyPhase = UWORD(7);
constexpr slong kRelationRandomExponentBits = 4;
constexpr slong kFailDivisor = 32;
constexpr slong kMinFail = 10;
struct RelationSearchResetState {
    slong candidates_tried = 0;
    slong analytic_finish_precision = kAnalyticFinishInitialPrecision;
    slong analytic_precision_doublings = 0;
    slong analytic_precision_inconclusive_count = 0;
    slong bounded_full_factor_base_restart_requests = 0;
    slong honesty_restart_requests = 0;
};

ulong next_relation_random_exponent(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) +
            UWORD(1442695040888963407);
    constexpr ulong mask =
            (UWORD(1) << kRelationRandomExponentBits) - 1;
    const ulong bits = static_cast<ulong>(8 * sizeof(ulong));
    const ulong shift = bits >
                    static_cast<ulong>(kRelationRandomExponentBits)
            ? bits - static_cast<ulong>(kRelationRandomExponentBits)
            : 0;
    return (state >> shift) & mask;
}

bool build_useful_pivot_indices(std::vector<slong>& useful,
                                     ClassGroupContext& context,
                                     const FactorBase& base,
                                     const std::vector<slong>& permutation,
                                     slong need,
                                     bool finish_unit_log_rotation_active,
                                     slong finish_unit_log_rotation) noexcept {
    if (!context.has_presentation()) {
        useful = permutation;
        return static_cast<slong>(useful.size()) == base.length();
    }

    std::vector<slong> nonprincipal;
    std::vector<char> hnf_covered;
    if (!build_nonprincipal_indices(nonprincipal, context, base) ||
        !build_hnf_covered_flags(hnf_covered, context, base) ||
        !build_uncovered_indices_from_flags(useful, hnf_covered,
                                            nonprincipal)) {
        return false;
    }
    if (!useful.empty()) {
        return true;
    }

    if (static_cast<slong>(permutation.size()) != base.length()) {
        return false;
    }

    if (finish_unit_log_rotation_active) {
        const slong length = static_cast<slong>(permutation.size());
        if (length <= 0) {
            return true;
        }
        const slong shift = length <= 1
                ? WORD(0)
                : finish_unit_log_rotation % length;
        useful.reserve(static_cast<std::size_t>(length));
        for (slong i = 0; i < length; ++i) {
            const slong pos = (i + shift) % length;
            const slong index = permutation[static_cast<std::size_t>(pos)];
            if (index < 0 || index >= base.length()) {
                return false;
            }
            useful.push_back(index);
        }
        return true;
    }

    const slong target = need <= 0
            ? WORD(1)
            : (need < base.length() ? need : base.length());
    useful.reserve(static_cast<std::size_t>(target));
    for (slong index : permutation) {
        if (index < 0 || index >= base.length()) {
            return false;
        }
        useful.push_back(index);
        if (static_cast<slong>(useful.size()) >= target) {
            break;
        }
    }
    return true;
}

void activate_finish_unit_log_rotation(
        RelationSearchControlState& route_state) noexcept {
    route_state.finish_unit_log_rotation_active = true;
    route_state.finish_unit_log_rotation = route_state.squash_index;
    ++route_state.squash_index;
}

void reset_route_after_factor_base_restart(
        RelationSearchControlState& route_state,
        const ClassGroupContext& context,
        const FactorBase& base) noexcept {
    RelationSearchResetState preserved;
    preserved.candidates_tried = route_state.candidates_tried;
    preserved.analytic_finish_precision =
            route_state.analytic_finish_precision;
    preserved.analytic_precision_doublings =
            route_state.analytic_precision_doublings;
    preserved.analytic_precision_inconclusive_count =
            route_state.analytic_precision_inconclusive_count;
    preserved.bounded_full_factor_base_restart_requests =
            route_state.bounded_full_factor_base_restart_requests;
    preserved.honesty_restart_requests = route_state.honesty_restart_requests;

    route_state = RelationSearchControlState{};
    route_state.candidates_tried = preserved.candidates_tried;
    route_state.accepted_relations = context.relation_count();
    route_state.fail_limit = base.length() + 1;
    route_state.analytic_finish_precision =
            preserved.analytic_finish_precision;
    route_state.analytic_precision_doublings =
            preserved.analytic_precision_doublings;
    route_state.analytic_precision_inconclusive_count =
            preserved.analytic_precision_inconclusive_count;
    route_state.bounded_full_factor_base_restart_requests =
            preserved.bounded_full_factor_base_restart_requests;
    route_state.honesty_restart_requests =
            preserved.honesty_restart_requests;
}

AnalyticFinishDecision request_honesty_factor_base_restart(
        RelationSearchControlState& route_state,
        const DiagnosticsContext* diagnostics) noexcept {
    route_state.bounded_full_factor_base_restart_pending = true;
    route_state.bounded_full_factor_base_restart_allow_past_half = true;
    ++route_state.bounded_full_factor_base_restart_requests;
    ++route_state.honesty_restart_requests;
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              "be_honest requested factor-base restart");
    return AnalyticFinishDecision::restart_factor_base;
}

AnalyticFinishDecision apply_honesty_check(
        ClassGroupContext& context,
        RelationSearchControlState& route_state,
        bool goal_reached,
        const Order& order,
        const SubfactorBaseSchedule* subfb_state,
        flint::FmpzConstRef active_factor_base_bound,
        bool factor_base_restart_available,
        bool use_required_prime_factor_over_base,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!goal_reached) {
        return AnalyticFinishDecision::incomplete;
    }
    flint::Fmpz installed_factor_base_bound;
    if (!context.factor_base_build_bound(
                flint::FmpzRef(installed_factor_base_bound)) ||
        !flint::fmpz_equal(
                flint::FmpzConstRef(installed_factor_base_bound),
                active_factor_base_bound)) {
        return AnalyticFinishDecision::failed;
    }
    if (route_state.honesty_checked) {
        return context.factor_base_generation_status() ==
                               ProofState::verified &&
                       context.factor_base_generation_checked_status() ==
                               ProofState::verified
                ? AnalyticFinishDecision::finished
                : AnalyticFinishDecision::failed;
    }

    flint::Fmpz required_bound;
    if (!factor_base_class_group_bound(flint::FmpzRef(required_bound),
                                       order)) {
        return AnalyticFinishDecision::failed;
    }
    if (flint::fmpz_cmp(flint::FmpzConstRef(required_bound),
                        active_factor_base_bound) <= 0) {
        if (!context.check_factor_base_generation_bound(
                    flint::FmpzConstRef(required_bound))) {
            return AnalyticFinishDecision::failed;
        }
        route_state.honesty_checked = true;
        return AnalyticFinishDecision::finished;
    }

    const FactorBase* base = context.factor_base();
    const ulong random_seed = base == nullptr
            ? UWORD(0)
            : relation_search_phase_seed(context, base->length(),
                                         kHonestyPhase,
                                         route_state.honesty_restart_requests);
    bool honest = false;
    if (base == nullptr ||
        !factor_base_honesty_check(
                honest, *base, active_factor_base_bound,
                flint::FmpzConstRef(required_bound), subfb_state,
                random_seed, use_required_prime_factor_over_base,
                route_state.analytic_finish_precision,
                diagnostics)) {
        return AnalyticFinishDecision::failed;
    }
    if (honest) {
        if (!detail::ClassGroupCertificationAccess::
                    record_factor_base_honesty_proof(
                            context,
                            flint::FmpzConstRef(required_bound))) {
            return AnalyticFinishDecision::failed;
        }
        route_state.honesty_checked = true;
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail, "be_honest check passed");
        return AnalyticFinishDecision::finished;
    }

    if (factor_base_restart_available) {
        return request_honesty_factor_base_restart(route_state,
                                                        diagnostics);
    }
    return AnalyticFinishDecision::failed;
}

template <class RunAutoSaturation>
AnalyticFinishDecision apply_old_cache_checkpoint(
        ClassGroupContext& context,
        RelationSearchControlState& route_state,
        bool& goal_reached,
        const ClassGroupRelationOptions& options,
        const DiagnosticsContext* diagnostics,
        RunAutoSaturation run_auto_saturation) noexcept {
    if (!goal_reached) {
        return AnalyticFinishDecision::incomplete;
    }

    const slong before = context.relation_count();
    if (!run_auto_saturation()) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "old_cache checkpoint saturation hook failed");
        return AnalyticFinishDecision::failed;
    }
    route_state.accepted_relations = context.relation_count();
    if (context.relation_count() != before) {
        route_state.checkpoint_relation_count = context.relation_count();
        return AnalyticFinishDecision::finished;
    }

    if (route_state.checkpoint_relation_count == context.relation_count() &&
        route_state.candidates_tried < options.max_candidates &&
        route_state.accepted_relations < options.max_relations) {
        goal_reached = false;
        ++route_state.checkpoint_extra_relation_requests;
        activate_finish_unit_log_rotation(route_state);
        route_state.finish_full_rank_relation_active = true;
        update_relation_need(route_state.completion, 1);
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "old_cache checkpoint requested one more relation");
        return AnalyticFinishDecision::old_cache_needs_more_relations;
    }

    route_state.checkpoint_relation_count = context.relation_count();
    return AnalyticFinishDecision::finished;
}

}  // namespace detail::relation_search

using namespace detail::relation_search;

bool ClassGroupContext::run_relation_search_route_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool emit_norm_prefilter_profile_event,
        bool clamp_relation_kernel_units_to_rank) noexcept {
    (void)emit_norm_prefilter_profile_event;
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    const bool native_quadratic_rank_zero_completion =
            class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_->audit.policy.selected &&
            class_unit_transaction_context_->audit.policy.degree == 2 &&
            class_unit_transaction_context_->audit.policy.unit_rank == 0 &&
            class_unit_transaction_context_->audit.policy.factor_base ==
                    detail::NativeFactorBaseStrategy::
                            relation_completion &&
            class_unit_transaction_context_->audit.policy.relations ==
                    detail::NativeRelationStrategy::relation_completion_table;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_route");
    if (private_storage_ == nullptr) {
        return false;
    }
    auto finish_workspace = [&]() noexcept
            -> detail::HnfFinishWorkspace* {
        return false && private_storage_ != nullptr
                ? &private_storage_->hnf_finish_workspace
                : nullptr;
    };
    if (detail::HnfFinishWorkspace* workspace =
                finish_workspace();
        workspace != nullptr) {
        workspace->clear();
    }

    slong relation_kernel_units_target = 0;
    if (clamp_relation_kernel_units_to_rank &&
        !relation_kernel_unit_target(relation_kernel_units_target,
                                     order, options)) {
        return false;
    }
    if (!clamp_relation_kernel_units_to_rank &&
        options.target_relation_kernel_units < 0) {
        return false;
    }
    if (!clamp_relation_kernel_units_to_rank) {
        relation_kernel_units_target =
                options.target_relation_kernel_units;
    }

    relation_kernel_units_target_ = relation_kernel_units_target;
    configure_partial_relations_(options);

    flint::Fmpz active_factor_base_bound;
    flint::fmpz_set(flint::FmpzRef(active_factor_base_bound),
                    factor_base_bound);
    flint::Fmpz restart_max_bound;
    const bool have_restart_max_bound =
            relation_factor_base_restart_limit(restart_max_bound, order);

    detail::SubfactorBaseSchedule subfb_state;
    {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.relation_completion."
                            "initial_subfactor_base");
        if (!build_initial_subfactor_base(
                    subfb_state, base_,
                    flint::FmpzConstRef(active_factor_base_bound))) {
            return false;
        }
    }
    if (subfb_state.subfactor_base.empty()) {
        return false;
    }

    detail::RelationAdmissionCache cache;
    if (relation_kernel_units_target >
        WORD_MAX - kRelationSurplus) {
        return false;
    }
    const slong add_need =
            kRelationSurplus + relation_kernel_units_target;
    if (!detail::initialize_relation_admission_cache(
                *this, cache, order, add_need)) {
        return false;
    }
    SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
              LogLevel::detail, "relation admission cache initialized");

    NormPrefilter norm_prefilter;
    {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.norm_prefilter");
        if (!build_norm_prefilter(
                    norm_prefilter, base_,
                    flint::FmpzConstRef(active_factor_base_bound), false)) {
            return false;
        }
    }

    RelationSearchControlState route_state;
    if (private_storage_->analytic_finish_precision > 0) {
        route_state.analytic_finish_precision = max_slong_value(
                route_state.analytic_finish_precision,
                private_storage_->analytic_finish_precision);
    }
    route_state.accepted_relations = relation_count();
    route_state.fail_limit = base_.length() + 1;
    std::vector<slong> small_multiplier(
            static_cast<std::size_t>(base_.length()), 0);

    bool goal_reached = false;
    auto need_with_cache_outstanding = [&](slong native_need) noexcept {
        const slong outstanding = max_slong_value(
                0, cache.target_relation_count - relation_count());
        return native_need < outstanding ? outstanding : native_need;
    };
    auto factor_base_restart_available = [&]() noexcept {
        if (!have_restart_max_bound) {
            return false;
        }
        flint::Fmpz next_bound;
        return next_relation_factor_base_bound(
                next_bound, flint::FmpzConstRef(active_factor_base_bound),
                flint::FmpzConstRef(restart_max_bound));
    };

    auto factor_base_restart_strict_half_available = [&]() noexcept {
        if (!have_restart_max_bound ||
            !flint::fmpz_fits_si(
                    flint::FmpzConstRef(active_factor_base_bound)) ||
            !flint::fmpz_fits_si(
                    flint::FmpzConstRef(restart_max_bound))) {
            return false;
        }
        const slong current_si = flint::fmpz_get_si(
                flint::FmpzConstRef(active_factor_base_bound));
        const slong max_si =
                flint::fmpz_get_si(flint::FmpzConstRef(restart_max_bound));
        if (current_si >= max_si / 2) {
            return false;
        }
        flint::Fmpz next_bound;
        return next_relation_factor_base_bound_to_limit(
                next_bound, flint::FmpzConstRef(active_factor_base_bound),
                flint::FmpzConstRef(restart_max_bound));
    };

    auto factor_base_restart_to_max_available = [&]() noexcept {
        if (!have_restart_max_bound) {
            return false;
        }
        flint::Fmpz next_bound;
        return next_relation_factor_base_bound_to_limit(
                next_bound, flint::FmpzConstRef(active_factor_base_bound),
                flint::FmpzConstRef(restart_max_bound));
    };

    auto perform_factor_base_restart = [&]() noexcept {
        if (!route_state.bounded_full_factor_base_restart_pending) {
            return true;
        }
        const bool allow_past_half =
                route_state.bounded_full_factor_base_restart_allow_past_half;
        route_state.bounded_full_factor_base_restart_pending = false;
        route_state.bounded_full_factor_base_restart_allow_past_half = false;

        flint::Fmpz next_bound;
        const bool have_next = have_restart_max_bound &&
                (allow_past_half
                         ? next_relation_factor_base_bound_to_limit(
                                   next_bound,
                                   flint::FmpzConstRef(
                                           active_factor_base_bound),
                                   flint::FmpzConstRef(
                                           restart_max_bound))
                         : next_relation_factor_base_bound(
                                   next_bound,
                                   flint::FmpzConstRef(
                                           active_factor_base_bound),
                                   flint::FmpzConstRef(
                                           restart_max_bound)));
        if (!have_next) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "factor-base restart skipped at bound cap");
            return true;
        }

        if (!detail::rebuild_relation_factor_base_and_replay(
                    *this, cache, order, flint::FmpzConstRef(next_bound),
                    add_need)) {
            return false;
        }
        if (detail::HnfFinishWorkspace* workspace =
                    finish_workspace();
            workspace != nullptr) {
            workspace->clear();
        }
        flint::fmpz_set(flint::FmpzRef(active_factor_base_bound),
                        flint::FmpzConstRef(next_bound));
        configure_partial_relations_(options);
        if (!build_initial_subfactor_base(
                    subfb_state, base_,
                    flint::FmpzConstRef(active_factor_base_bound)) ||
            subfb_state.subfactor_base.empty() ||
            !build_norm_prefilter(
                    norm_prefilter, base_,
                    flint::FmpzConstRef(active_factor_base_bound), false)) {
            return false;
        }
        small_multiplier.assign(static_cast<std::size_t>(base_.length()), 0);
        reset_route_after_factor_base_restart(route_state, *this, base_);
        goal_reached = false;
        SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "factor-base restart rebuilt and replayed relations");
        return true;
    };

    auto apply_finish_checkpoint = [&]() noexcept {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.finish_checkpoint");
        AnalyticFinishDecision decision = apply_old_cache_checkpoint(
                *this, route_state, goal_reached, options, active_diagnostics,
                [&]() noexcept {
                    return try_auto_relation_saturation_(
                            order,
                            flint::FmpzConstRef(active_factor_base_bound),
                            options);
                });
        if (decision != AnalyticFinishDecision::finished) {
            return decision;
        }

        decision = apply_analytic_finish_check(
                *this, route_state, goal_reached, order, options,
                false ||
                        native_quadratic_rank_zero_completion,
                factor_base_restart_available(),
                finish_workspace(),
                active_diagnostics);
        if (decision != AnalyticFinishDecision::finished) {
            return decision;
        }
        private_storage_->analytic_finish_precision =
                route_state.analytic_finish_precision;

        if (class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_
                    ->defer_relation_saturation_until_units &&
            options.target_relation_kernel_units > 0) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "be_honest deferred until class/unit validation");
            return AnalyticFinishDecision::finished;
        }

        decision = apply_honesty_check(
                *this, route_state, goal_reached, order, &subfb_state,
                flint::FmpzConstRef(active_factor_base_bound),
                factor_base_restart_to_max_available(),
                false,
                active_diagnostics);
        return decision;
    };

    auto handle_finish_checkpoint = [&](bool& restarted) noexcept {
        restarted = false;
        const AnalyticFinishDecision decision = apply_finish_checkpoint();
        if (decision == AnalyticFinishDecision::failed) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "finish checkpoint failed");
            return false;
        }
        if (decision == AnalyticFinishDecision::restart_factor_base) {
            if (!perform_factor_base_restart()) {
                return false;
            }
            restarted = true;
        }
        return true;
    };

    auto compute_search_relation_need = [&](RelationNeedDecision& decision) noexcept {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.relation_completion.compute_need");
        return compute_relation_need(
                decision, *this, base_, relation_kernel_units_target);
    };

    for (;;) {
        RelationNeedDecision need_decision;
        if (!compute_search_relation_need(need_decision)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "native need computation failed before small_norm");
            return false;
        }
        slong next_need =
                need_with_cache_outstanding(need_decision.need);
        goal_reached = need_decision.goal_reached && next_need <= 0;
        update_relation_need(route_state.completion, next_need);
        bool factor_base_restarted = false;
        if (!handle_finish_checkpoint(factor_base_restarted)) {
            return false;
        }
        if (factor_base_restarted) {
            continue;
        }
        if (goal_reached) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before small_norm: goal reached");
            break;
        }
        if (route_state.completion.need <= 0) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before small_norm: need exhausted");
            break;
        }
        if (route_state.candidates_tried >= options.max_candidates) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before small_norm: candidate cap");
            break;
        }
        if (route_state.accepted_relations >= options.max_relations) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before small_norm: relation cap");
            break;
        }

        cache.target_relation_count =
                relation_collection_target_count(
                        *this, options, route_state.completion);
        const bool have_published_presentation = has_presentation();
        if (have_published_presentation) {
            route_state.fail_limit = max_slong_value(
                    base_.length() / kFailDivisor, kMinFail);
        }
        const slong small_norm_relation_cap =
                cache.relation_count_before_init +
                2 * base_.length() +
                2 * (relation_kernel_units_target + 1) +
                kRelationSurplus;
        bool small_norm_attempted = false;
        // Strict reference treats Silex Ideal HNF as the route-neutral lattice
        // representation consumed by reference's finite_quadratic_form_ideal_search control flow.
        // Keep the public legacy selector's existing degree guard unchanged.
        const bool runtime_small_norm_available = false;
        const bool runtime_direct_t2_available =
                false || order.degree() <= 4;
        if (!runtime_small_norm_available &&
            !runtime_direct_t2_available) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before unsourced direct T2: "
                      "nf.zk column-HNF input unavailable");
            return false;
        }
        if (runtime_small_norm_available &&
            options.ideal_search_radius > 0 &&
            order.degree() > 2 &&
            (route_state.done_small <= base_.length() + 1 ||
             have_published_presentation) &&
            route_state.small_fail <= route_state.fail_limit &&
            relation_count() < small_norm_relation_cap) {
            small_norm_attempted = true;
            const slong before = relation_count();
            const slong small_norm_multiplier_index =
                    route_state.done_small % (base_.length() + 1);
            if (!collect_small_norm_relations(
                        *this, order, &norm_prefilter,
                        cache.target_relation_count,
                        relation_kernel_units_target, options, cache,
                        subfb_state, route_state,
                        small_norm_multiplier_index,
                        !have_published_presentation,
                        small_multiplier, false)) {
                SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "small_norm relation collection failed");
                return false;
            }
            if (!have_published_presentation && relation_count() != before) {
                route_state.small_fail = 0;
            } else {
                ++route_state.small_fail;
            }
            ++route_state.done_small;
        }

        if (!compute_search_relation_need(need_decision)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "native need computation failed before random_relation_search");
            return false;
        }
        next_need = need_with_cache_outstanding(need_decision.need);
        goal_reached = need_decision.goal_reached && next_need <= 0;
        update_relation_need(route_state.completion, next_need);
        factor_base_restarted = false;
        if (!handle_finish_checkpoint(factor_base_restarted)) {
            return false;
        }
        if (factor_base_restarted) {
            continue;
        }
        if (goal_reached) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before random_relation_search: goal reached");
            break;
        }
        if (small_norm_attempted) {
            continue;
        }
        if (route_state.completion.need <= 0) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before random_relation_search: need exhausted");
            break;
        }
        if (route_state.candidates_tried >= options.max_candidates) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before random_relation_search: candidate cap");
            break;
        }
        if (route_state.accepted_relations >= options.max_relations) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped before random_relation_search: relation cap");
            break;
        }

        cache.target_relation_count =
                relation_collection_target_count(
                        *this, options, route_state.completion);
        if (!collect_random_relations(
                    *this, order, &norm_prefilter,
                    cache.target_relation_count,
                    relation_kernel_units_target, options, cache,
                    subfb_state, route_state,
                    factor_base_restart_strict_half_available(),
                    factor_base_restart_to_max_available(),
                    active_diagnostics)) {
            return false;
        }
        if (route_state.bounded_full_factor_base_restart_pending) {
            if (!perform_factor_base_restart()) {
                return false;
            }
            continue;
        }

        if (!compute_search_relation_need(need_decision)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "native need computation failed after random_relation_search");
            return false;
        }
        next_need = need_with_cache_outstanding(need_decision.need);
        goal_reached = need_decision.goal_reached && next_need <= 0;
        update_relation_need(route_state.completion, next_need);
        factor_base_restarted = false;
        if (!handle_finish_checkpoint(factor_base_restarted)) {
            return false;
        }
        if (factor_base_restarted) {
            continue;
        }
        if (goal_reached) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "route stopped after random_relation_search: goal reached");
            break;
        }
    }

    detail::ClassGroupContextStorage& final_route_storage =
            *private_storage_;
    flint::Arb& cached_analytic_finish_product =
            final_route_storage.analytic_finish_product;
    if (goal_reached && route_state.analytic_finish_product_valid &&
        route_state.analytic_finish_product_precision > 0) {
        flint::arb_set(
                flint::ArbRef(cached_analytic_finish_product),
                flint::ArbConstRef(
                        route_state.analytic_finish_product));
        final_route_storage.analytic_finish_product_precision =
                route_state.analytic_finish_product_precision;
        final_route_storage.analytic_finish_product_valid = true;
    } else {
        flint::arb_zero(
                flint::ArbRef(cached_analytic_finish_product));
        final_route_storage.analytic_finish_product_precision = 0;
        final_route_storage.analytic_finish_product_valid = false;
    }

    return goal_reached;
}

}  // namespace silex
