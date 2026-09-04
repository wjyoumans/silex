#include <silex/class_group.hpp>

#include <silex/unit.hpp>

#include "class_group_internal.hpp"
#include "relation_admission_cache_internal.hpp"
#include "relation_search_driver_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_factor_base_plan_internal.hpp"
#include "relation_search_internal.hpp"
#include "../order/order_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"

#include <utility>
#include <vector>

namespace silex {
namespace {

// A standalone candidate is deliberately an unproven, subgroup-ready
// relation result.  Resource ceilings remain caller controls; relation-search
// scheduling is private and deterministic.  The fixed rank-zero/quadratic and
// rank-one surpluses retain the source-traced relation/unit search envelopes.
bool standalone_candidate_policy(
        detail::ClassGroupRelationOptions& out,
        const Order& order,
        const ClassGroupCandidateOptions& requested) noexcept {
    if (requested.max_candidates < 0 || requested.max_relations < 0 ||
        !detail::order_has_parented_basis(order)) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent()) || rank < 0) {
        return false;
    }

    out = detail::ClassGroupRelationOptions{};
    if (order.degree() == 2 && order.is_maximal()) {
        out.coordinate_search_radius = 4;
        out.ideal_search_radius = 2;
    } else if (order.degree() >= 3) {
        out.ideal_search_radius = 12;
    }
    constexpr slong kRankOneRelationSurplus = 5;
    out.target_relation_kernel_units = rank;
    if (rank == 0 && order.degree() == 2) {
        out.target_relation_kernel_units = 5;
    }
    if (rank == 1 && order.degree() != 3) {
        out.target_relation_kernel_units += kRankOneRelationSurplus;
    }
    out.max_candidates = requested.max_candidates;
    out.max_relations = requested.max_relations;
    out.requested_certification = CertificationMode::unknown;
    out.diagnostics = requested.diagnostics;
    return true;
}

}  // namespace

namespace detail::relation_search {
}  // namespace detail::relation_search

using namespace detail::relation_search;

bool ClassGroupContext::run_relation_production_route_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool emit_norm_prefilter_profile_event) noexcept {
    (void) emit_norm_prefilter_profile_event;
    const SpecializedRelationBackendStatus quadratic_status =
            run_maximal_imaginary_quadratic_relation_backend_(
                    order, factor_base_bound, options);
    if (quadratic_status ==
        SpecializedRelationBackendStatus::succeeded) {
        return true;
    }
    if (quadratic_status == SpecializedRelationBackendStatus::failed) {
        return false;
    }
    return run_relation_search_route_(order, factor_base_bound, options,
                                    emit_norm_prefilter_profile_event, false);
}

bool ClassGroupContext::run_relation_production_prepass_route_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool emit_norm_prefilter_profile_event) noexcept {
    (void) emit_norm_prefilter_profile_event;
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_completion_prepass_route");
    SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
              LogLevel::detail,
              "relation_completion relation prepass using native FLINT-backed "
              "relation primitives");
    if (private_storage_ == nullptr) {
        return false;
    }
    detail::ClassGroupContextStorage& execution_storage = *private_storage_;
    const bool native_hnf_finish =
            class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_->audit.policy.selected &&
            class_unit_transaction_context_->audit.policy.units ==
                    detail::NativeUnitStrategy::hnf;
    detail::HnfFinishWorkspace* const finish_workspace =
            native_hnf_finish
                    ? &execution_storage.hnf_finish_workspace
                    : nullptr;
    if (finish_workspace != nullptr) {
        finish_workspace->clear();
    }
    flint::Arb& cached_analytic_finish_product =
            execution_storage.analytic_finish_product;
    execution_storage.finish_state.reset();
    {
        SILEX_PROFILE_SCOPE(active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group.contract_scheduler_prepass");
        slong relation_kernel_units_target = 0;
        const FactorBase* base = factor_base();
        if (base != nullptr && base->is_defined() &&
            relation_kernel_unit_target(
                    relation_kernel_units_target, order, options) &&
            relation_kernel_units_target <=
                    WORD_MAX - kRelationSurplus) {
            relation_kernel_units_target_ =
                    relation_kernel_units_target;
            configure_partial_relations_(options);

            detail::SubfactorBaseSchedule subfactor_base_schedule;
            detail::RelationAdmissionCache cache;
            NormPrefilter norm_prefilter;
            const slong add_need =
                    kRelationSurplus +
                    relation_kernel_units_target;
            const bool prepass_ready =
                    build_initial_subfactor_base(
                            subfactor_base_schedule, *base,
                            factor_base_bound) &&
                    !subfactor_base_schedule.subfactor_base.empty() &&
                    detail::initialize_relation_admission_cache(
                            *this, cache, order, add_need) &&
                    build_norm_prefilter(norm_prefilter, *base,
                                         factor_base_bound, false);
            if (prepass_ready) {
                RelationSearchControlState route_state;
                route_state.accepted_relations = relation_count();
                route_state.fail_limit = base->length() + 1;
                std::vector<slong> small_multiplier(
                        static_cast<std::size_t>(base->length()), 0);
                auto store_live_state = [&]() noexcept {
                    detail::RelationFinishState next_state;
                    next_state.cache = cache;
                    next_state.subfactor_base_schedule =
                            subfactor_base_schedule;
                    next_state.route_state = std::move(route_state);
                    next_state.small_multiplier = small_multiplier;
                    flint::fmpz_set(
                            flint::FmpzRef(next_state.factor_base_bound),
                            factor_base_bound);
                    next_state.target_relation_kernel_units =
                            relation_kernel_units_target;
                    execution_storage.finish_state.emplace(
                            std::move(next_state));
                };

                bool goal_reached = false;
                auto compute_need =
                        [&](RelationNeedDecision& decision) noexcept {
                    return compute_relation_need(
                            decision, *this, *base,
                            relation_kernel_units_target);
                };
                auto update_need =
                        [&](RelationNeedDecision& decision) noexcept {
                    const slong outstanding = max_slong_value(
                            0, cache.target_relation_count -
                                       relation_count());
                    const slong next_need =
                            decision.need < outstanding
                                    ? outstanding
                                    : decision.need;
                    goal_reached = decision.goal_reached && next_need <= 0;
                    update_relation_need(route_state.completion, next_need);
                };

                RelationNeedDecision need_decision;
                bool small_norm_attempted = false;
                bool prepass_ok = compute_need(need_decision);
                if (prepass_ok) {
                    update_need(need_decision);
                }
                if (prepass_ok && !goal_reached &&
                    route_state.completion.need > 0 &&
                    route_state.candidates_tried < options.max_candidates &&
                    route_state.accepted_relations < options.max_relations) {
                    cache.target_relation_count =
                            relation_collection_target_count(
                                    *this, options, route_state.completion);
                    const slong small_norm_relation_cap =
                            cache.relation_count_before_init +
                            2 * base->length() +
                            2 * (relation_kernel_units_target + 1) +
                            kRelationSurplus;
                    if (options.ideal_search_radius > 0 &&
                        order.degree() > 2 &&
                        relation_count() < small_norm_relation_cap) {
                        small_norm_attempted = true;
                        const slong before = relation_count();
                        SILEX_PROFILE_SCOPE(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.contract_small_norm_prepass");
                        prepass_ok = collect_small_norm_relations(
                                *this, order, &norm_prefilter,
                                cache.target_relation_count,
                                relation_kernel_units_target, options,
                                cache, subfactor_base_schedule, route_state, 0,
                                true,
                                small_multiplier, true);
                        if (relation_count() != before) {
                            route_state.small_fail = 0;
                        } else {
                            ++route_state.small_fail;
                        }
                        ++route_state.done_small;
                    }
                }
                if (prepass_ok && compute_need(need_decision)) {
                    update_need(need_decision);
                }
                if (prepass_ok && !goal_reached && !small_norm_attempted &&
                    route_state.completion.need > 0 &&
                    route_state.candidates_tried < options.max_candidates &&
                    route_state.accepted_relations < options.max_relations) {
                    cache.target_relation_count =
                            relation_collection_target_count(
                                    *this, options, route_state.completion);
                    SILEX_PROFILE_SCOPE(
                            active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group.contract_random_relations_prepass");
                    prepass_ok = collect_random_relations(
                            *this, order, &norm_prefilter,
                            cache.target_relation_count,
                            relation_kernel_units_target, options, cache,
                            subfactor_base_schedule, route_state, false, false,
                            active_diagnostics);
                    if (prepass_ok && compute_need(need_decision)) {
                        update_need(need_decision);
                    }
                }
                slong unit_target_rank = -1;
                // Source relation-completion setup computes `invhr` before the
                // relation loop and then lets `reconstruct_regulator` validate the live
                // class/unit finish.  Run that existing finish check here so
                // an accepted `compute_inverse_residue` product is part of the reference
                // finish state before HNF
                // unit publication.  Class-group-only callers keep the old
                // fast prepass exit unless they explicitly request relation
                // kernel unit witnesses.
                const bool strict_finish_guard =
                        class_unit_transaction_context_ != nullptr &&
                        class_unit_transaction_context_
                                ->defer_relation_saturation_until_units;
                const bool opportunistic_finish_guard =
                        !strict_finish_guard &&
                        options.target_relation_kernel_units > 0;
                const bool run_finish_guard =
                        (strict_finish_guard ||
                         opportunistic_finish_guard) &&
                        order.parent() != nullptr &&
                        unit_rank(unit_target_rank, *order.parent()) &&
                        unit_target_rank > 0;
                while (prepass_ok && goal_reached &&
                       run_finish_guard) {
                    const bool goal_reached_before_finish = goal_reached;
                    AnalyticFinishDecision finish_decision =
                            apply_analytic_finish_check(
                                    *this, route_state, goal_reached, order,
                                    options, false, false, finish_workspace,
                                    active_diagnostics);
                    bool honesty_attempted = false;
                    if (finish_decision == AnalyticFinishDecision::finished) {
                        // The reference finish orders the compact-factor-base
                        // honesty proof after regulator reconstruction and
                        // before fundamental-unit extraction.  Share the
                        // ordinary relation route's proof checkpoint here so
                        // the prepass cannot publish an analytically finished
                        // class/unit state with unverified factor-base
                        // generation.
                        honesty_attempted = true;
                        finish_decision = apply_honesty_check(
                                *this, route_state, goal_reached, order,
                                &subfactor_base_schedule, factor_base_bound,
                                false, false, active_diagnostics);
                    }
                    if (finish_decision == AnalyticFinishDecision::finished) {
                        execution_storage.analytic_finish_precision =
                                route_state.analytic_finish_precision;
                        if (route_state.analytic_finish_product_valid &&
                            route_state.analytic_finish_product_precision >
                                    0) {
                            flint::arb_set(
                                    flint::ArbRef(
                                            cached_analytic_finish_product),
                                    flint::ArbConstRef(
                                            route_state
                                                    .analytic_finish_product));
                            execution_storage
                                    .analytic_finish_product_precision =
                                    route_state
                                            .analytic_finish_product_precision;
                            execution_storage
                                    .analytic_finish_product_valid = true;
                        } else {
                            flint::arb_zero(
                                    flint::ArbRef(
                                            cached_analytic_finish_product));
                            execution_storage
                                    .analytic_finish_product_precision = 0;
                            execution_storage
                                    .analytic_finish_product_valid = false;
                        }
                        store_live_state();
                        return true;
                    }
                    if (honesty_attempted) {
                        // Once the analytic finish reaches the honesty
                        // checkpoint, failure to establish the proof is a
                        // terminal prepass outcome.  In particular, do not
                        // reinterpret it as the opportunistic analytic
                        // fallback below and publish the unchecked state.
                        prepass_ok = false;
                        break;
                    }
                    if (opportunistic_finish_guard) {
                        // The non-deferred prepass is allowed to seed HNF
                        // unit publication with a reference compute_inverse_residue product
                        // only after reconstruct_regulator accepts it.  If reconstruct_regulator asks
                        // for more work, keep the previous prepass result and
                        // leave the product unpublished; the strict deferred
                        // path below still follows the source retry loop.
                        goal_reached = goal_reached_before_finish;
                        flint::arb_zero(
                                flint::ArbRef(
                                        cached_analytic_finish_product));
                        execution_storage
                                .analytic_finish_product_precision = 0;
                        execution_storage
                                .analytic_finish_product_valid = false;
                        SILEX_PROFILE_EVENT(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.contract_finish_deferred");
                        break;
                    }
                    if (finish_decision == AnalyticFinishDecision::failed ||
                        finish_decision ==
                                AnalyticFinishDecision::restart_factor_base) {
                        prepass_ok = false;
                        break;
                    }
                    if (route_state.completion.need <= 0 ||
                        route_state.candidates_tried >=
                                options.max_candidates ||
                        route_state.accepted_relations >=
                                options.max_relations) {
                        break;
                    }

                    cache.target_relation_count =
                            relation_collection_target_count(
                                    *this, options, route_state.completion);
                    bool small_norm_refinement_attempted = false;
                    const slong small_norm_relation_cap =
                            cache.relation_count_before_init +
                            2 * base->length() +
                            2 * (relation_kernel_units_target + 1) +
                            kRelationSurplus;
                    // Source relation-completion setup handles regulator reconstruction
                    // relation-defect requests by re-entering the
                    // small_norm branch first, with the full-rank "lie" and
                    // squash-index rotation active, before falling through to
                    // random_relation_search.  Native FLINT-backed small_norm contract input
                    // is available here, so mirror that high-level scheduling
                    // rather than adding random rows directly.
                    if (options.ideal_search_radius > 0 &&
                        order.degree() > 2 &&
                        (route_state.done_small <= base->length() + 1 ||
                         has_presentation()) &&
                        route_state.small_fail <= route_state.fail_limit &&
                        relation_count() < small_norm_relation_cap) {
                        small_norm_refinement_attempted = true;
                        const slong before = relation_count();
                        const slong small_norm_multiplier_index =
                                route_state.done_small %
                                (base->length() + 1);
                        SILEX_PROFILE_SCOPE(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.contract_analytic_small_norm");
                        prepass_ok = collect_small_norm_relations(
                                *this, order, &norm_prefilter,
                                cache.target_relation_count,
                                relation_kernel_units_target, options,
                                cache, subfactor_base_schedule, route_state,
                                small_norm_multiplier_index, false,
                                small_multiplier, true);
                        if (relation_count() != before) {
                            route_state.small_fail = 0;
                        } else {
                            ++route_state.small_fail;
                        }
                        ++route_state.done_small;
                    }
                    if (!small_norm_refinement_attempted) {
                        SILEX_PROFILE_SCOPE(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.contract_analytic_refinement");
                        prepass_ok = collect_random_relations(
                                *this, order, &norm_prefilter,
                                cache.target_relation_count,
                                relation_kernel_units_target, options,
                                cache, subfactor_base_schedule, route_state,
                                false, false,
                                active_diagnostics);
                    }
                    if (prepass_ok && compute_need(need_decision)) {
                        update_need(need_decision);
                    }
                }
                if (prepass_ok && goal_reached &&
                    (!run_finish_guard ||
                     opportunistic_finish_guard)) {
                    store_live_state();
                    return true;
                }
            }
        }
    }
    return false;
}

bool ClassGroupContext::extend_relation_kernel_units_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.extend_relation_kernel_units");
    if (!detail::order_has_parented_basis(order) || !order.is_maximal() ||
        !same_order_parent(parent(), &order) || !has_factor_base() ||
        options.coordinate_search_radius < 1 ||
        options.ideal_search_radius < 0 ||
        options.target_relation_kernel_units < 0 ||
        options.post_finite_refinement_phase_budget < 0 ||
        options.max_candidates < 0 || options.max_relations < 0 ||
        (options.relation_saturation_aux_prime_bound != 0 &&
         options.relation_saturation_aux_prime_bound < 2) ||
        options.relation_saturation_max_appends_per_ell < 0 ||
        !detail::valid_certification_request(options.requested_certification)) {
        return false;
    }

    return run_native_experimental_relation_route_(
            order, factor_base_bound, options,
            false);
}

bool ClassGroupContext::compute_relation_candidate_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.compute_candidate");
    SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::class_group,
                  VerboseLevel::progress,
                  "class-group candidate computation started");
    if (!detail::order_has_parented_basis(order) || !order.is_maximal() ||
        options.coordinate_search_radius < 1 ||
        options.ideal_search_radius < 0 ||
        options.target_relation_kernel_units < 0 ||
        options.post_finite_refinement_phase_budget < 0 ||
        options.max_candidates < 0 || options.max_relations < 0 ||
        (options.relation_saturation_aux_prime_bound != 0 &&
         options.relation_saturation_aux_prime_bound < 2) ||
        options.relation_saturation_max_appends_per_ell < 0 ||
        !detail::valid_certification_request(options.requested_certification)) {
        return false;
    }

    SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.build_factor_base");
    const bool paired_transaction =
            class_unit_transaction_context_ != nullptr;
    const detail::ClassUnitExecutionPolicy* const policy =
            paired_transaction
            ? &class_unit_transaction_context_->audit.policy
            : nullptr;
    const bool native_complete_factor_base =
            policy != nullptr &&
            policy->factor_base ==
                    detail::NativeFactorBaseStrategy::complete;
    const bool native_complete_factor_base_strategy_valid =
            native_complete_factor_base &&
            policy->relations ==
                    detail::NativeRelationStrategy::native_experimental;
    const bool native_relation_completion_strategy_valid =
            policy != nullptr &&
            policy->factor_base ==
                    detail::NativeFactorBaseStrategy::
                            relation_completion &&
            ((order.degree() == 2 &&
              policy->relations ==
                      detail::NativeRelationStrategy::
                              relation_completion_table) ||
             (order.degree() >= 3 && policy->relations ==
                      detail::NativeRelationStrategy::
                              relation_completion_prepass));
    if (paired_transaction &&
        (!class_unit_transaction_context_->audit.policy.selected ||
         (class_unit_transaction_context_->audit.policy.factor_base !=
                  detail::NativeFactorBaseStrategy::complete &&
          class_unit_transaction_context_->audit.policy.factor_base !=
                  detail::NativeFactorBaseStrategy::relation_completion) ||
         ((native_complete_factor_base &&
           !native_complete_factor_base_strategy_valid) ||
          (class_unit_transaction_context_->audit.policy.factor_base ==
                   detail::NativeFactorBaseStrategy::relation_completion &&
           !native_relation_completion_strategy_valid)))) {
        return false;
    }
    flint::Fmpz route_factor_base_bound;
    flint::fmpz_set(flint::FmpzRef(route_factor_base_bound),
                    factor_base_bound);
    bool retry_native_requested_bound = false;
    {
        SILEX_PROFILE_SCOPE(active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group.route_bound_select");
        constexpr ulong kNativeRelationCompletionStartMinPublicBound = 512;
        if (!paired_transaction &&
            order.degree() >= 5 &&
            flint::fmpz_cmp_ui(factor_base_bound,
                               kNativeRelationCompletionStartMinPublicBound) > 0) {
            detail::RelationFactorBasePlan plan;
            if (build_relation_factor_base_plan(
                        plan, order, active_diagnostics) &&
                flint::fmpz_cmp(
                        flint::FmpzConstRef(plan.working_bound),
                        factor_base_bound) < 0) {
                flint::fmpz_set(flint::FmpzRef(route_factor_base_bound),
                                flint::FmpzConstRef(plan.working_bound));
                retry_native_requested_bound = true;
            }
        }
        if (!flint::fmpz_equal(flint::FmpzConstRef(route_factor_base_bound),
                               factor_base_bound)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "relation completion initial relation factor-base bound "
                      "selected");
        }
    }

    const auto try_route_bound = [&](flint::FmpzConstRef bound) noexcept {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.route_bound_attempt");
        ClassGroupContext candidate(order);
        candidate.class_unit_transaction_context_ = class_unit_transaction_context_;
        candidate.set_diagnostics(active_diagnostics);
        bool factor_base_ok = false;
        {
            SILEX_PROFILE_SCOPE(active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route_bound_factor_base");
            if (candidate.is_defined()) {
                if (paired_transaction &&
                        class_unit_transaction_context_->audit.policy
                                        .factor_base ==
                                detail::NativeFactorBaseStrategy::
                                        relation_completion) {
                    factor_base_ok =
                            candidate.build_relation_factor_base_(bound);
                } else {
                    factor_base_ok = candidate.build_factor_base(bound);
                }
            }
        }
        if (!factor_base_ok) {
            return false;
        }
        bool search_ok = false;
        {
            SILEX_PROFILE_SCOPE(active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route_bound_search");
            search_ok = candidate.run_native_experimental_relation_route_(
                    order, bound, options, true);
        }
        if (!search_ok ||
            (options.requested_certification == CertificationMode::unknown &&
             candidate.certification_status() != CertificationMode::unknown)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "class-group candidate computation did not reach goal");
            return false;
        }

        {
            SILEX_PROFILE_SCOPE(active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route_bound_commit");
            swap(candidate);
            set_diagnostics(active_diagnostics);
        }
        return true;
    };

    bool route_ok =
            try_route_bound(flint::FmpzConstRef(route_factor_base_bound));
    if (!route_ok && retry_native_requested_bound) {
        SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "native relation route retrying requested factor-base bound");
        route_ok = try_route_bound(factor_base_bound);
    }
    if (!route_ok) {
        return false;
    }
    SILEX_DEBUG_CHECK(
            active_diagnostics, DiagnosticsModule::class_group,
            DebugLevel::cheap, "class-group relation source count",
            relations_.length() == static_cast<slong>(relation_sources_.size()));
    SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::class_group,
                  VerboseLevel::progress,
                  "class-group candidate computation finished");
    return true;
}

bool ClassGroupContext::compute_candidate(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupCandidateOptions& options) noexcept {
    detail::ClassGroupRelationOptions policy;
    return class_unit_transaction_context_ == nullptr &&
           standalone_candidate_policy(policy, order, options) &&
           compute_relation_candidate_(order, factor_base_bound, policy);
}

}  // namespace silex
