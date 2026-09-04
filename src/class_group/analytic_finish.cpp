#include <silex/class_group.hpp>

#include <silex/embedding.hpp>
#include <silex/order_element.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include "class_group_internal.hpp"
#include "ideal_t2_enumeration_internal.hpp"
#include "relation_search_driver_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"
#include "../order_unit/order_unit_internal.hpp"
#include "../zeta/zeta_internal.hpp"

#include <cstdio>
#include <utility>

namespace silex::detail::relation_search {

constexpr slong kAnalyticFinishMaxPrecision = 40960;
constexpr slong kPrecisionRestartPeriod = 8;

enum class AnalyticQuotientDecision {
    finished,
    relation_defect,
    precision_inconclusive
};

static AnalyticFinishDecision request_analytic_extra_relations(
        RelationSearchControlState& route_state,
        bool& goal_reached,
        slong need,
        bool full_rank_finish_request,
        const DiagnosticsContext* diagnostics,
        const char* message) noexcept {
#if !defined(SILEX_ENABLE_LOGGING) || !SILEX_ENABLE_LOGGING
    (void)message;
#endif
    if (need <= 0) {
        need = 1;
    }
    goal_reached = false;
    ++route_state.analytic_extra_relation_requests;
    activate_finish_unit_log_rotation(route_state);
    route_state.finish_full_rank_relation_active = full_rank_finish_request;
    update_relation_need(route_state.completion, need);
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              message);
    return AnalyticFinishDecision::analytic_needs_more_relations;
}

static slong next_finish_precision(slong precision) noexcept {
    if (precision <= 0) {
        return kAnalyticFinishInitialPrecision;
    }
    if (precision < 1280) {
        return precision <= kAnalyticFinishMaxPrecision / 2
                ? 2 * precision
                : kAnalyticFinishMaxPrecision;
    }
    if (precision <= kAnalyticFinishMaxPrecision / 3) {
        return 3 * precision / 2;
    }
    return kAnalyticFinishMaxPrecision;
}

static bool increase_analytic_finish_precision(
        RelationSearchControlState& route_state,
        const DiagnosticsContext* diagnostics) noexcept {
    if (route_state.analytic_finish_precision >=
        kAnalyticFinishMaxPrecision) {
        return false;
    }

    route_state.analytic_finish_precision = next_finish_precision(
            route_state.analytic_finish_precision);
    ++route_state.analytic_precision_doublings;
    ++route_state.analytic_precision_inconclusive_count;
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              "analytic finish check increased precision");
    return true;
}

static bool precision_restart_due(
        const RelationSearchControlState& route_state) noexcept {
    return route_state.analytic_precision_doublings > 0 &&
           route_state.analytic_precision_doublings %
                           kPrecisionRestartPeriod ==
                   kPrecisionRestartPeriod - 1;
}

static AnalyticFinishDecision request_factor_base_restart(
        RelationSearchControlState& route_state,
        const DiagnosticsContext* diagnostics) noexcept {
    route_state.bounded_full_factor_base_restart_pending = true;
    route_state.bounded_full_factor_base_restart_allow_past_half = false;
    ++route_state.bounded_full_factor_base_restart_requests;
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              "precision boundary requested factor-base restart");
    return AnalyticFinishDecision::restart_factor_base;
}

static RegulatorPivotFinishAction route_regulator_pivot_outcome(
        RelationSearchControlState& route_state,
        bool& goal_reached,
        RegulatorPivotFinishPath path,
        RegulatorPivotOutcome outcome,
        slong independent_count,
        slong target_rank,
        bool factor_base_restart_available,
        const DiagnosticsContext* diagnostics) noexcept {
    if (outcome == RegulatorPivotOutcome::invalid ||
        independent_count < 0 || target_rank < 0 ||
        independent_count > target_rank) {
        return RegulatorPivotFinishAction::failed;
    }
    if (outcome == RegulatorPivotOutcome::precision_inconclusive) {
        if (factor_base_restart_available &&
            precision_restart_due(route_state)) {
            (void)request_factor_base_restart(route_state, diagnostics);
            return RegulatorPivotFinishAction::restart_factor_base;
        }
        if (!increase_analytic_finish_precision(route_state,
                                                     diagnostics)) {
            return RegulatorPivotFinishAction::failed;
        }
        return RegulatorPivotFinishAction::retry_precision;
    }
    if (independent_count < target_rank) {
        // A certified pivot defect is certify_regulator_multiple's unit-column
        // defect on either live path.  The later bad_check quotient defect is
        // the distinct branch that activates full-rank relation scheduling.
        const bool full_rank_finish_request = false;
        (void)request_analytic_extra_relations(
                route_state, goal_reached, target_rank - independent_count,
                full_rank_finish_request, diagnostics,
                path == RegulatorPivotFinishPath::regulator_product
                        ? "reconstruct_regulator pivot rank defect requested more "
                          "relations"
                        : "certify_regulator_multiple HNF unit-column defect "
                          "requested more relations");
        return RegulatorPivotFinishAction::request_relations;
    }
    return RegulatorPivotFinishAction::proceed;
}

static bool route_control_equal(
        const RelationSearchControlState& left,
        const RelationSearchControlState& right) noexcept {
    return left.completion.need == right.completion.need &&
           left.completion.old_need == right.completion.old_need &&
           left.completion.dependent_trials ==
                   right.completion.dependent_trials &&
           left.completion.subfactor_base_trials ==
                   right.completion.subfactor_base_trials &&
           left.checkpoint_relation_count == right.checkpoint_relation_count &&
           left.done_small == right.done_small &&
           left.small_fail == right.small_fail &&
           left.fail_limit == right.fail_limit &&
           left.candidates_tried == right.candidates_tried &&
           left.accepted_relations == right.accepted_relations &&
           left.checkpoint_extra_relation_requests ==
                   right.checkpoint_extra_relation_requests &&
           left.analytic_extra_relation_requests ==
                   right.analytic_extra_relation_requests &&
           left.analytic_finish_precision == right.analytic_finish_precision &&
           left.analytic_finish_product_precision ==
                   right.analytic_finish_product_precision &&
           left.analytic_precision_doublings ==
                   right.analytic_precision_doublings &&
           left.analytic_precision_inconclusive_count ==
                   right.analytic_precision_inconclusive_count &&
           left.bounded_full_factor_base_restart_requests ==
                   right.bounded_full_factor_base_restart_requests &&
           left.honesty_restart_requests == right.honesty_restart_requests &&
           left.squash_index == right.squash_index &&
           left.finish_unit_log_rotation == right.finish_unit_log_rotation &&
           left.bounded_full_factor_base_restart_pending ==
                   right.bounded_full_factor_base_restart_pending &&
           left.bounded_full_factor_base_restart_allow_past_half ==
                   right.bounded_full_factor_base_restart_allow_past_half &&
           left.finish_unit_log_rotation_active ==
                   right.finish_unit_log_rotation_active &&
           left.finish_full_rank_relation_active ==
                   right.finish_full_rank_relation_active &&
           left.honesty_checked == right.honesty_checked &&
           left.analytic_finish_product_valid ==
                   right.analytic_finish_product_valid &&
           left.analytic_finish_torsion_valid ==
                   right.analytic_finish_torsion_valid;
}

RegulatorPivotFinishTestResult
regulator_pivot_finish_for_testing(
        RegulatorPivotFinishPath path,
        RegulatorPivotOutcome outcome,
        slong independent_count,
        slong target_rank,
        bool factor_base_restart_available,
        slong analytic_finish_precision,
        slong analytic_precision_doublings) noexcept {
    RelationSearchControlState route_state;
    route_state.completion.need = 17;
    route_state.completion.old_need = 19;
    route_state.completion.dependent_trials = 23;
    route_state.completion.subfactor_base_trials = 29;
    route_state.checkpoint_relation_count = 31;
    route_state.done_small = 37;
    route_state.small_fail = 41;
    route_state.fail_limit = 43;
    route_state.candidates_tried = 47;
    route_state.accepted_relations = 53;
    route_state.checkpoint_extra_relation_requests = 59;
    route_state.analytic_extra_relation_requests = 61;
    route_state.analytic_finish_precision = analytic_finish_precision;
    route_state.analytic_finish_product_precision = 67;
    route_state.analytic_precision_doublings = analytic_precision_doublings;
    route_state.analytic_precision_inconclusive_count = 71;
    route_state.bounded_full_factor_base_restart_requests = 73;
    route_state.honesty_restart_requests = 79;
    route_state.squash_index = 83;
    route_state.finish_unit_log_rotation = 89;
    route_state.bounded_full_factor_base_restart_pending = false;
    route_state.bounded_full_factor_base_restart_allow_past_half = true;
    route_state.finish_unit_log_rotation_active = false;
    route_state.finish_full_rank_relation_active = false;
    route_state.honesty_checked = true;
    route_state.analytic_finish_product_valid = true;
    route_state.analytic_finish_torsion_valid = true;
    const RelationSearchControlState before = std::move(route_state);
    // Reconstruct the move-only state's scalar portion without relying on a
    // copy operation; the testing seam intentionally observes route mutation.
    route_state.completion = before.completion;
    route_state.checkpoint_relation_count = before.checkpoint_relation_count;
    route_state.done_small = before.done_small;
    route_state.small_fail = before.small_fail;
    route_state.fail_limit = before.fail_limit;
    route_state.candidates_tried = before.candidates_tried;
    route_state.accepted_relations = before.accepted_relations;
    route_state.checkpoint_extra_relation_requests =
            before.checkpoint_extra_relation_requests;
    route_state.analytic_extra_relation_requests =
            before.analytic_extra_relation_requests;
    route_state.analytic_finish_precision = before.analytic_finish_precision;
    route_state.analytic_finish_product_precision =
            before.analytic_finish_product_precision;
    route_state.analytic_precision_doublings =
            before.analytic_precision_doublings;
    route_state.analytic_precision_inconclusive_count =
            before.analytic_precision_inconclusive_count;
    route_state.bounded_full_factor_base_restart_requests =
            before.bounded_full_factor_base_restart_requests;
    route_state.honesty_restart_requests = before.honesty_restart_requests;
    route_state.squash_index = before.squash_index;
    route_state.finish_unit_log_rotation = before.finish_unit_log_rotation;
    route_state.bounded_full_factor_base_restart_pending =
            before.bounded_full_factor_base_restart_pending;
    route_state.bounded_full_factor_base_restart_allow_past_half =
            before.bounded_full_factor_base_restart_allow_past_half;
    route_state.finish_unit_log_rotation_active =
            before.finish_unit_log_rotation_active;
    route_state.finish_full_rank_relation_active =
            before.finish_full_rank_relation_active;
    route_state.honesty_checked = before.honesty_checked;
    route_state.analytic_finish_product_valid =
            before.analytic_finish_product_valid;
    route_state.analytic_finish_torsion_valid =
            before.analytic_finish_torsion_valid;
    bool goal_reached = true;

    RegulatorPivotFinishTestResult result;
    result.goal_reached_before = goal_reached;
    result.relation_need_before = before.completion.need;
    result.completion_old_need_before = before.completion.old_need;
    result.completion_dependent_trials_before =
            before.completion.dependent_trials;
    result.completion_subfactor_base_trials_before =
            before.completion.subfactor_base_trials;
    result.analytic_extra_relation_requests_before =
            before.analytic_extra_relation_requests;
    result.analytic_finish_precision_before =
            before.analytic_finish_precision;
    result.analytic_precision_doublings_before =
            before.analytic_precision_doublings;
    result.analytic_precision_inconclusive_before =
            before.analytic_precision_inconclusive_count;
    result.factor_base_restart_requests_before =
            before.bounded_full_factor_base_restart_requests;
    result.factor_base_restart_pending_before =
            before.bounded_full_factor_base_restart_pending;
    result.factor_base_restart_allow_past_half_before =
            before.bounded_full_factor_base_restart_allow_past_half;
    result.finish_unit_log_rotation_before = before.finish_unit_log_rotation;
    result.squash_index_before = before.squash_index;
    result.candidates_tried_before = before.candidates_tried;
    result.accepted_relations_before = before.accepted_relations;
    result.finish_unit_log_rotation_active_before =
            before.finish_unit_log_rotation_active;
    result.finish_full_rank_relation_active_before =
            before.finish_full_rank_relation_active;

    result.action = route_regulator_pivot_outcome(
            route_state, goal_reached, path, outcome, independent_count,
            target_rank, factor_base_restart_available, nullptr);
    result.goal_reached_after = goal_reached;
    result.relation_need_after = route_state.completion.need;
    result.completion_old_need_after = route_state.completion.old_need;
    result.completion_dependent_trials_after =
            route_state.completion.dependent_trials;
    result.completion_subfactor_base_trials_after =
            route_state.completion.subfactor_base_trials;
    result.analytic_extra_relation_requests_after =
            route_state.analytic_extra_relation_requests;
    result.analytic_finish_precision_after =
            route_state.analytic_finish_precision;
    result.analytic_precision_doublings_after =
            route_state.analytic_precision_doublings;
    result.analytic_precision_inconclusive_after =
            route_state.analytic_precision_inconclusive_count;
    result.factor_base_restart_requests_after =
            route_state.bounded_full_factor_base_restart_requests;
    result.factor_base_restart_pending_after =
            route_state.bounded_full_factor_base_restart_pending;
    result.factor_base_restart_allow_past_half_after =
            route_state.bounded_full_factor_base_restart_allow_past_half;
    result.finish_unit_log_rotation_after = route_state.finish_unit_log_rotation;
    result.squash_index_after = route_state.squash_index;
    result.candidates_tried_after = route_state.candidates_tried;
    result.accepted_relations_after = route_state.accepted_relations;
    result.finish_unit_log_rotation_active_after =
            route_state.finish_unit_log_rotation_active;
    result.finish_full_rank_relation_active_after =
            route_state.finish_full_rank_relation_active;
    result.relation_control_state_unchanged =
            before.completion.need == route_state.completion.need &&
            before.completion.old_need == route_state.completion.old_need &&
            before.completion.dependent_trials ==
                    route_state.completion.dependent_trials &&
            before.completion.subfactor_base_trials ==
                    route_state.completion.subfactor_base_trials &&
            before.checkpoint_relation_count ==
                    route_state.checkpoint_relation_count &&
            before.done_small == route_state.done_small &&
            before.small_fail == route_state.small_fail &&
            before.fail_limit == route_state.fail_limit &&
            before.candidates_tried == route_state.candidates_tried &&
            before.accepted_relations == route_state.accepted_relations &&
            before.checkpoint_extra_relation_requests ==
                    route_state.checkpoint_extra_relation_requests &&
            before.analytic_extra_relation_requests ==
                    route_state.analytic_extra_relation_requests &&
            before.honesty_restart_requests ==
                    route_state.honesty_restart_requests &&
            before.squash_index == route_state.squash_index &&
            before.finish_unit_log_rotation == route_state.finish_unit_log_rotation &&
            before.finish_unit_log_rotation_active ==
                    route_state.finish_unit_log_rotation_active &&
            before.finish_full_rank_relation_active ==
                    route_state.finish_full_rank_relation_active &&
            before.honesty_checked == route_state.honesty_checked &&
            before.analytic_finish_product_valid ==
                    route_state.analytic_finish_product_valid &&
            before.analytic_finish_torsion_valid ==
                    route_state.analytic_finish_torsion_valid;
    result.all_control_state_unchanged =
            goal_reached && route_control_equal(before, route_state);
    return result;
}

static AnalyticQuotientDecision bad_check_quotient(flint::ArbConstRef quotient,
                                             slong precision) noexcept {
    if (precision <= 0 || !flint::arb_is_finite(quotient) ||
        !flint::arb_is_positive(quotient)) {
        return AnalyticQuotientDecision::precision_inconclusive;
    }

    flint::Fmpq low_q;
    flint::Fmpq high_q;
    flint::fmpq_set_si(low_q, 3, 4);
    flint::fmpq_set_si(high_q, 13, 10);
    flint::Arb low;
    flint::Arb high;
    flint::arb_set_fmpq(low, low_q, precision);
    flint::arb_set_fmpq(high, high_q, precision);
    flint::Arb value;
    flint::arb_set(flint::ArbRef(value), quotient);
    if (flint::arb_gt(value, high)) {
        return AnalyticQuotientDecision::relation_defect;
    }
    if (flint::arb_lt(value, low) ||
        flint::arb_contains_fmpq(flint::ArbConstRef(value), low_q) ||
        flint::arb_contains_fmpq(flint::ArbConstRef(value), high_q)) {
        return AnalyticQuotientDecision::precision_inconclusive;
    }
    return AnalyticQuotientDecision::finished;
}

static const char* bad_check_decision_name(AnalyticQuotientDecision decision)
        noexcept {
    switch (decision) {
        case AnalyticQuotientDecision::finished:
            return "finished";
        case AnalyticQuotientDecision::relation_defect:
            return "relation_defect";
        case AnalyticQuotientDecision::precision_inconclusive:
            return "precision_inconclusive";
    }
    return "unknown";
}

static void log_analytic_finish_check(
        const DiagnosticsContext* diagnostics,
        const ClassGroupContext& context,
        flint::ArbConstRef quotient,
        slong precision,
        AnalyticQuotientDecision decision) noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (!log_enabled(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::detail)) {
        return;
    }
    char detail[256];
    std::snprintf(detail, sizeof(detail),
                  "relations=%ld rank=%ld quotient_mid=%.17g "
                  "precision=%ld decision=%s",
                  static_cast<long>(context.relation_count()),
                  static_cast<long>(context.relation_rank()),
                  arb_midpoint_double(quotient),
                  static_cast<long>(precision),
                  bad_check_decision_name(decision));
    log_emit(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
             __func__, "analytic finish check", detail);
#else
    (void)diagnostics;
    (void)context;
    (void)quotient;
    (void)precision;
    (void)decision;
#endif
}

static bool analytic_finish_product(flint::ArbRef out,
                                  RelationSearchControlState& route_state,
                                  const Order& order,
                                  slong precision,
                                  const DiagnosticsContext* diagnostics,
                                  detail::ZetaBfResidueDegreeCache*
                                          residue_degree_cache)
        noexcept {
    if (route_state.analytic_finish_product_valid &&
        route_state.analytic_finish_product_precision == precision) {
        flint::arb_set(out, flint::ArbConstRef(
                                    route_state.analytic_finish_product));
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.analytic_finish_inverse_residue_product");
        flint::Fmpz torsion_order;
        Element torsion_generator;
        if (detail::class_regulator_product_estimate_with_diagnostics(
                    flint::ArbRef(route_state.analytic_finish_product), order,
                    precision, diagnostics, residue_degree_cache,
                    &torsion_order, &torsion_generator)) {
            route_state.analytic_finish_product_valid = true;
            route_state.analytic_finish_product_precision = precision;
            flint::fmpz_set(
                    flint::FmpzRef(
                            route_state.analytic_finish_torsion_order),
                    flint::FmpzConstRef(torsion_order));
            route_state.analytic_finish_torsion_generator =
                    std::move(torsion_generator);
            route_state.analytic_finish_torsion_valid =
                    flint::fmpz_sgn(flint::FmpzConstRef(
                            route_state.analytic_finish_torsion_order)) > 0 &&
                    route_state.analytic_finish_torsion_generator
                            .is_defined();
            flint::arb_set(
                    out,
                    flint::ArbConstRef(route_state.analytic_finish_product));
            return true;
        }
    }

    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.analytic_finish_zeta_product");
    if (!zeta_class_regulator_product(
                flint::ArbRef(route_state.analytic_finish_product), order,
                precision)) {
        route_state.analytic_finish_product_valid = false;
        route_state.analytic_finish_product_precision = 0;
        route_state.analytic_finish_torsion_generator.clear();
        route_state.analytic_finish_torsion_valid = false;
        return false;
    }
    route_state.analytic_finish_product_valid = true;
    route_state.analytic_finish_product_precision = precision;
    route_state.analytic_finish_torsion_generator.clear();
    route_state.analytic_finish_torsion_valid = false;
    flint::arb_set(out,
                   flint::ArbConstRef(route_state.analytic_finish_product));
    return true;
}

AnalyticFinishDecision apply_analytic_finish_check(
        ClassGroupContext& context,
        RelationSearchControlState& route_state,
        bool& goal_reached,
        const Order& order,
        const ClassGroupRelationOptions& options,
        bool allow_rank_zero_regulator_reconstruction,
        bool factor_base_restart_available,
        HnfFinishWorkspace* finish_workspace,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!goal_reached) {
        return AnalyticFinishDecision::incomplete;
    }
    if (route_state.candidates_tried >= options.max_candidates ||
        route_state.accepted_relations >= options.max_relations ||
        !context.has_presentation() || order.parent() == nullptr) {
        return AnalyticFinishDecision::finished;
    }

    for (;;) {
        const slong precision = route_state.analytic_finish_precision;
        EmbeddingContext embeddings(*order.parent());
        flint::Arb analytic_hR;
        flint::Arb candidate_hR;
        flint::Arb quotient;

        if (!embeddings.is_defined() || !embeddings.refine(precision)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "analytic finish embedding setup failed");
            return AnalyticFinishDecision::failed;
        }

        slong unit_target_rank = -1;
        slong independent_units = -1;
        if (!unit_rank(unit_target_rank, *order.parent())) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "analytic finish unit-rank setup failed");
            return AnalyticFinishDecision::failed;
        }
        const RegulatorPivotOutcome preprobe_outcome =
                detail::hnf_independent_unit_count(
                        independent_units, order, context, embeddings,
                        finish_workspace, precision);
        const RegulatorPivotFinishAction preprobe_action =
                route_regulator_pivot_outcome(
                        route_state, goal_reached,
                        RegulatorPivotFinishPath::independent_unit_preprobe,
                        preprobe_outcome, independent_units, unit_target_rank,
                        factor_base_restart_available, diagnostics);
        if (preprobe_action ==
            RegulatorPivotFinishAction::retry_precision) {
            continue;
        }
        if (preprobe_action ==
            RegulatorPivotFinishAction::restart_factor_base) {
            return AnalyticFinishDecision::restart_factor_base;
        }
        if (preprobe_action ==
            RegulatorPivotFinishAction::request_relations) {
            return AnalyticFinishDecision::analytic_needs_more_relations;
        }
        if (preprobe_action == RegulatorPivotFinishAction::failed) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "analytic finish HNF unit-column probe failed");
            return AnalyticFinishDecision::failed;
        }

        if (!analytic_finish_product(flint::ArbRef(analytic_hR),
                                          route_state, order, precision,
                                          diagnostics,
                                          detail::ClassUnitTransactionAccess::
                                                  relation_factor_base_plan_residue_degrees(
                                                          context,
                                                          order))) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "analytic finish zeta product failed");
            return AnalyticFinishDecision::failed;
        }

        slong product_independent_units = 0;
        const RegulatorPivotOutcome product_outcome =
                (unit_target_rank > 0 || allow_rank_zero_regulator_reconstruction)
                        ? detail::hnf_class_regulator_product(
                                  flint::ArbRef(candidate_hR),
                                  product_independent_units, order, context,
                                  embeddings, flint::ArbConstRef(analytic_hR),
                                  finish_workspace, precision)
                        : RegulatorPivotOutcome::precision_inconclusive;
        const RegulatorPivotFinishAction product_action =
                route_regulator_pivot_outcome(
                        route_state, goal_reached,
                        RegulatorPivotFinishPath::regulator_product,
                        product_outcome, product_independent_units,
                        unit_target_rank, factor_base_restart_available,
                        diagnostics);
        if (product_action ==
            RegulatorPivotFinishAction::retry_precision) {
            continue;
        }
        if (product_action ==
            RegulatorPivotFinishAction::restart_factor_base) {
            return AnalyticFinishDecision::restart_factor_base;
        }
        if (product_action ==
            RegulatorPivotFinishAction::request_relations) {
            return AnalyticFinishDecision::analytic_needs_more_relations;
        }
        if (product_action == RegulatorPivotFinishAction::failed) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "reconstruct_regulator pivot/product failed");
            return AnalyticFinishDecision::failed;
        }

        flint::arb_div(quotient, candidate_hR, analytic_hR, precision);
        const AnalyticQuotientDecision check = bad_check_quotient(
                flint::ArbConstRef(quotient), precision);
        log_analytic_finish_check(diagnostics, context,
                                       flint::ArbConstRef(quotient),
                                       precision, check);
        if (check == AnalyticQuotientDecision::finished) {
            return AnalyticFinishDecision::finished;
        }
        if (check == AnalyticQuotientDecision::relation_defect) {
            break;
        }
        if (factor_base_restart_available &&
            precision_restart_due(route_state)) {
            return request_factor_base_restart(route_state,
                                                    diagnostics);
        }
        if (!increase_analytic_finish_precision(route_state,
                                                     diagnostics)) {
            return AnalyticFinishDecision::failed;
        }
    }

    return request_analytic_extra_relations(
            route_state, goal_reached, 1, true, diagnostics,
            "reconstruct_regulator-style check requested one more relation");
}

}  // namespace silex::detail::relation_search
