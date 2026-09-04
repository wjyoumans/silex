#include <silex/order_unit.hpp>

#include "compute_internal.hpp"
#include "order_unit_internal.hpp"
#include "../class_group/class_group_certification_internal.hpp"
#include "../zeta/zeta_internal.hpp"

#include <silex/class_group.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include <string>

namespace silex::detail {

const char* validate_refine_outcome_name(ValidateRefineOutcome outcome)
        noexcept {
    switch (outcome) {
    case ValidateRefineOutcome::not_run: return "not_run";
    case ValidateRefineOutcome::proven: return "proven";
    case ValidateRefineOutcome::not_proven_request:
        return "not_proven_request";
    case ValidateRefineOutcome::analytic_unavailable:
        return "analytic_unavailable";
    case ValidateRefineOutcome::analytic_inconclusive:
        return "analytic_inconclusive";
    case ValidateRefineOutcome::index_one_publication_failed:
        return "index_one_publication_failed";
    case ValidateRefineOutcome::class_saturation_failed:
        return "class_saturation_failed";
    case ValidateRefineOutcome::no_progress: return "no_progress";
    case ValidateRefineOutcome::unit_refinement_failed:
        return "unit_refinement_failed";
    case ValidateRefineOutcome::local_rank_unavailable:
        return "local_rank_unavailable";
    case ValidateRefineOutcome::pass_cap_reached:
        return "pass_cap_reached";
    }
    return "unknown";
}

const char* validate_refine_outcome_profile_label(
        ValidateRefineOutcome outcome) noexcept {
    switch (outcome) {
    case ValidateRefineOutcome::not_run:
        return "unit_group.validation_outcome.not_run";
    case ValidateRefineOutcome::proven:
        return "unit_group.validation_outcome.proven";
    case ValidateRefineOutcome::not_proven_request:
        return "unit_group.validation_outcome.not_proven_request";
    case ValidateRefineOutcome::analytic_unavailable:
        return "unit_group.validation_outcome.analytic_unavailable";
    case ValidateRefineOutcome::analytic_inconclusive:
        return "unit_group.validation_outcome.analytic_inconclusive";
    case ValidateRefineOutcome::index_one_publication_failed:
        return "unit_group.validation_outcome.index_one_publication_failed";
    case ValidateRefineOutcome::class_saturation_failed:
        return "unit_group.validation_outcome.class_saturation_failed";
    case ValidateRefineOutcome::no_progress:
        return "unit_group.validation_outcome.no_progress";
    case ValidateRefineOutcome::unit_refinement_failed:
        return "unit_group.validation_outcome.unit_refinement_failed";
    case ValidateRefineOutcome::local_rank_unavailable:
        return "unit_group.validation_outcome.local_rank_unavailable";
    case ValidateRefineOutcome::pass_cap_reached:
        return "unit_group.validation_outcome.pass_cap_reached";
    }
    return "unit_group.validation_outcome.unknown";
}

bool analytic_class_regulator_product_for_validation(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return detail::zeta_class_regulator_product_with_diagnostics(
            out, order, precision, diagnostics, residue_degree_base,
            residue_degree_cache);
}

bool bf_class_regulator_product_for_validation(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return detail::zeta_class_regulator_product_bf_audit_with_diagnostics(
            out, error_bound, cutoff, work_precision, order, max_cutoff,
            precision, diagnostics, residue_degree_base,
            residue_degree_cache);
}

bool class_regulator_product_for_validation(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return detail::
            zeta_class_regulator_product_validation_with_diagnostics(
                    out, error_bound, cutoff, work_precision, order,
                    max_cutoff, precision, diagnostics, residue_degree_base,
                    residue_degree_cache);
}

slong relation_saturation_max_appends_per_ell(
        const ClassGroupComputeOptions& options) noexcept {
    return options.relation_saturation_max_appends_per_ell > 0
            ? options.relation_saturation_max_appends_per_ell
            : kComputeSatAuxTarget;
}

void relation_saturation_aux_bound(
        flint::Fmpz& out,
        const ClassGroupComputeOptions& options) noexcept {
    const slong bound = options.relation_saturation_aux_prime_bound >= 2
            ? options.relation_saturation_aux_prime_bound
            : kComputeSatAuxMax;
    flint::fmpz_set_si(flint::FmpzRef(out), bound);
}

bool relation_saturation_retry_aux_bound(
        flint::Fmpz& out,
        const ClassGroupComputeOptions& options) noexcept {
    relation_saturation_aux_bound(out, options);
    if (options.relation_saturation_aux_prime_bound >= 2 ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(out),
                           kComputeProofRetryAuxMax) >= 0) {
        return false;
    }

    flint::fmpz_set_si(flint::FmpzRef(out), kComputeProofRetryAuxMax);
    return true;
}

namespace {

bool validation_index_bound_from_product(
        flint::Fmpz& index_bound,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision) noexcept {
    const Order* const order = class_group.parent();
    slong expected_rank = -1;
    if (order == nullptr || order->parent() == nullptr ||
        !class_group.has_presentation() || !units.is_set() ||
        !same_order_parent(units.parent(), order) ||
        !unit_rank(expected_rank, *order->parent()) ||
        units.free_rank() != expected_rank) {
        return false;
    }

    flint::Arb candidate_class_regulator_product;
    if (!units.class_regulator_product(
                flint::ArbRef(candidate_class_regulator_product),
                class_group, precision)) {
        return false;
    }

    const bool integer_index_prerequisites =
            class_group.factor_base_generation_status() ==
                    ProofState::verified &&
            class_group.factor_base_generation_checked_status() ==
                    ProofState::verified;
    // In this validation context the quotient is an integer class/unit
    // index only after factor-base generation has been proved.  Keep the
    // public ratio accessor conservative, but recognize the source-complete
    // index-one case before taking an outward-rounded ceiling.
    if (integer_index_prerequisites &&
        class_regulator_index_is_one_from_candidate_product(
                flint::ArbConstRef(candidate_class_regulator_product),
                analytic_class_regulator_product, precision,
                units.diagnostics())) {
        flint::fmpz_one(flint::FmpzRef(index_bound));
        return true;
    }

    if (!class_regulator_index_bound_from_candidate_product(
                flint::FmpzRef(index_bound),
                flint::ArbConstRef(candidate_class_regulator_product),
                analytic_class_regulator_product, precision,
                units.diagnostics())) {
        return false;
    }
    return integer_index_prerequisites ||
           flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) > 0;
}

}  // namespace

bool analytic_index_bound_for_validation(
        flint::Fmpz& index_bound,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_analytic_index_bound");
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound.ensure_product");
        if (!analytic_cache.ensure(order, precision, units.diagnostics(),
                                   class_group.factor_base())) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound.index_bound");
        return validation_index_bound_from_product(
                index_bound, units, class_group, analytic_cache.value(),
                precision);
    }
}

bool class_unit_validation_estimate(
        flint::Fmpz& index_bound,
        flint::Arb& expected_regulator,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_estimate");
    if (!analytic_index_bound_for_validation(
                index_bound, units, class_group, analytic_cache, order,
                precision) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) < 0) {
        return false;
    }

    flint::Arb regulator;
    if (!units.regulator(flint::ArbRef(regulator)) ||
        !flint::arb_is_finite(regulator)) {
        return false;
    }

    flint::arb_abs(expected_regulator, regulator);
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) > 0) {
        flint::arb_div_fmpz(expected_regulator, expected_regulator,
                            flint::FmpzConstRef(index_bound), precision);
    }
    return flint::arb_is_finite(expected_regulator) &&
           flint::arb_is_positive(expected_regulator);
}

bool class_unit_bf_validation_estimate(
        flint::Fmpz& index_bound,
        flint::Arb& expected_regulator,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.bf_validation_estimate");
    if (max_cutoff == 0 || order.degree() <= 2) {
        return class_unit_validation_estimate(
                index_bound, expected_regulator, units, class_group,
                analytic_cache, order, precision);
    }

    const bool use_validation =
            analytic_cache.validation_enabled();
    if (use_validation) {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.bf_validation_estimate."
                "ensure_validation");
        if (!analytic_cache.ensure_validation(
                    order, max_cutoff, precision, diagnostics,
                    class_group.factor_base())) {
            return class_unit_validation_estimate(
                    index_bound, expected_regulator, units, class_group,
                    analytic_cache, order, precision);
        }
    } else {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.bf_validation_estimate.ensure_bf_audit");
        if (!analytic_cache.ensure_bf_audit(
                    order, max_cutoff, precision, diagnostics,
                    class_group.factor_base())) {
            return class_unit_validation_estimate(
                    index_bound, expected_regulator, units, class_group,
                    analytic_cache, order, precision);
        }
        analytic_cache.seed(order, analytic_cache.bf_value(),
                            analytic_cache.bf_work_precision());
    }

    const flint::ArbConstRef validation_value = use_validation
            ? analytic_cache.zeta_validation_value()
            : analytic_cache.bf_value();
    if (!validation_index_bound_from_product(
                index_bound, units, class_group, validation_value,
                precision) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) < 0) {
        return false;
    }

    flint::Arb regulator;
    if (!units.regulator(flint::ArbRef(regulator)) ||
        !flint::arb_is_finite(regulator)) {
        return false;
    }

    flint::arb_abs(expected_regulator, regulator);
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) > 0) {
        flint::arb_div_fmpz(expected_regulator, expected_regulator,
                            flint::FmpzConstRef(index_bound), precision);
    }
    return flint::arb_is_finite(expected_regulator) &&
           flint::arb_is_positive(expected_regulator);
}

bool expected_regulator_stop(
        const OrderUnitGroup& units,
        const flint::Arb& expected_regulator,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.expected_regulator_stop");
    if (precision <= 0 || !flint::arb_is_finite(expected_regulator) ||
        !flint::arb_is_positive(expected_regulator)) {
        return false;
    }

    flint::Arb regulator;
    flint::Arb half_regulator;
    if (!units.regulator(flint::ArbRef(regulator)) ||
        !flint::arb_is_finite(regulator)) {
        return false;
    }

    flint::arb_abs(regulator, regulator);
    flint::arb_div_ui(half_regulator, regulator, 2, precision);
    return flint::arb_gt(expected_regulator, half_regulator);
}

bool validate_refine_continue_with_bound(
        bool have_previous,
        bool progress_since_previous,
        flint::FmpzConstRef previous,
        flint::FmpzConstRef current) noexcept {
    return !have_previous ||
           flint::fmpz_cmp(current, previous) < 0 ||
           progress_since_previous;
}

bool validate_refine_bound_decreased(flint::FmpzConstRef current,
                                     flint::FmpzConstRef next) noexcept {
    return flint::fmpz_cmp_ui(next, 1) > 0 &&
           flint::fmpz_cmp(next, current) < 0;
}

bool validate_refine_has_progress(bool class_changed,
                                  bool next_check_available,
                                  flint::FmpzConstRef current,
                                  flint::FmpzConstRef next) noexcept {
    return class_changed ||
           (next_check_available &&
            (flint::fmpz_is_one(next) ||
             validate_refine_bound_decreased(current, next)));
}

bool validate_refine_extra_pass_allowed(
        const ValidateRefineSummary& summary,
        const ClassGroupComputeOptions& options,
        slong extra_passes_done,
        slong extra_pass_cap) noexcept {
    return options.requested_certification == CertificationMode::proven &&
           extra_passes_done >= 0 && extra_pass_cap > 0 &&
           extra_passes_done < extra_pass_cap &&
           (summary.class_progress || summary.unit_progress) &&
           flint::fmpz_cmp_ui(
                   flint::FmpzConstRef(summary.last_index_bound), 1) > 0 &&
           summary.outcome != ValidateRefineOutcome::proven;
}

bool validation_target_bump(slong& bump,
                            const ValidateRefineSummary& summary,
                            slong rank,
                            slong target,
                            slong max_relations,
                            bool requested_proven,
                            bool rank_complete) noexcept {
    bump = 0;
    if (!requested_proven || !rank_complete || rank <= 0 || target < 0 ||
        max_relations < target ||
        flint::fmpz_cmp_ui(
                flint::FmpzConstRef(summary.last_index_bound), 1) <= 0 ||
        !flint::fmpz_fits_si(
                flint::FmpzConstRef(summary.last_index_bound))) {
        return false;
    }

    if (summary.outcome != ValidateRefineOutcome::local_rank_unavailable &&
        summary.outcome != ValidateRefineOutcome::no_progress) {
        return false;
    }

    slong extra =
            min_slong(flint::fmpz_get_si(
                              flint::FmpzConstRef(summary.last_index_bound)) -
                              1,
                      max_relations - target);
    extra = min_slong(extra, kComputeAnalyticBumpStepMax);
    if (extra <= 0 || target + extra > max_relations) {
        return false;
    }

    bump = extra;
    return true;
}

bool missing_unit_target_bump(slong& bump,
                              slong rank,
                              slong unit_count,
                              slong target,
                              slong max_relations) noexcept {
    bump = 0;
    if (rank <= 1 || unit_count < 0 || unit_count >= rank || target < 0 ||
        max_relations < target) {
        return false;
    }

    // reference `relation_completion_parameters` adds the unit-rank defect `RU - 1 - zc` to the
    // current relation need.  Preserve the cheap low-target attempt first, but
    // after it fails follow reference `FindUnits.jl:find_candidates`: test
    // candidates from `max(10, r1 + r2 - 1)` relation rows instead of walking
    // one small unit-rank defect at a time.
    const slong batch =
            max_slong(kUnitCandidateBatchTarget, rank);
    slong extra = rank - unit_count;
    if (target < batch) {
        extra = max_slong(extra, batch - target);
    } else {
        extra = max_slong(extra, batch);
    }
    extra = min_slong(extra, max_relations - target);
    if (extra <= 0 || target + extra > max_relations) {
        return false;
    }

    bump = extra;
    return true;
}

bool class_group_order_is_one(const ClassGroupContext& class_group) noexcept {
    flint::Fmpz order;
    return class_group.order(flint::FmpzRef(order)) &&
           flint::fmpz_is_one(flint::FmpzConstRef(order));
}

slong initial_relation_kernel_target(slong rank,
                                     slong degree,
                                     slong lower_target) noexcept {
    if (rank <= 0) {
        return lower_target;
    }
    if (rank == 1 && degree == 3) {
        return lower_target;
    }
    if (rank <= 1) {
        return max_slong(lower_target, rank + kRankOneAdaptiveSurplus);
    }
    return lower_target;
}

bool saturate_candidate_class_relations_with_units(
        bool& changed,
        ClassGroupContext& class_group,
        const OrderUnitGroup& units,
        const ClassGroupComputeOptions& options,
        flint::FmpzConstRef index_bound,
        bool index_bound_is_exact) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.validation_class_relation_saturation");
    changed = false;
    if ((options.requested_certification != CertificationMode::grh &&
         options.requested_certification != CertificationMode::proven) ||
        flint::fmpz_cmp_ui(index_bound, 1) <= 0 || !units.is_set() ||
        !same_order_parent(units.parent(), class_group.parent()) ||
        !class_group.has_presentation() ||
        class_group.relation_count() >= options.max_relations) {
        return true;
    }

    flint::Fmpz aux_bound;
    relation_saturation_aux_bound(aux_bound, options);
    bool saturated = false;
    const slong remaining = options.max_relations - class_group.relation_count();
    const bool saturation_ok = index_bound_is_exact
            ? ClassGroupCertificationAccess::
                      saturate_relations_for_index_bound_with_units(
                              changed, saturated, class_group, units,
                              index_bound, flint::FmpzConstRef(aux_bound),
                              relation_saturation_max_appends_per_ell(options),
                              remaining)
            : class_group.saturate_relations_bounded_with_units(
                      changed, saturated, units,
                      flint::FmpzConstRef(aux_bound),
                      relation_saturation_max_appends_per_ell(options),
                      remaining);
    if (!saturation_ok) {
        return false;
    }

    if (changed) {
        flint::Fmpz required_bound;
        if (class_group.factor_base_generation_bound(
                    flint::FmpzRef(required_bound))) {
            (void) class_group.check_factor_base_generation_bound(
                    flint::FmpzConstRef(required_bound));
        }
    }

    return true;
}

bool try_certify_candidate_with_zeta(ClassGroupContext& class_group,
                                     OrderUnitGroup& units,
                                     const Order& order,
                                     const ClassGroupComputeOptions& options,
                                     AnalyticClassRegulatorCache& analytic_cache,
                                     slong precision) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_zeta_certify");
    if (options.requested_certification != CertificationMode::proven) {
        return false;
    }
    if (options.zeta_bf_max_cutoff != 0 && order.degree() > 2) {
        if (analytic_cache.validation_enabled()) {
            bool have_validation = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_zeta_validation");
                have_validation =
                        analytic_cache.ensure_validation(
                                order, options.zeta_bf_max_cutoff, precision,
                                diagnostics, class_group.factor_base());
            }
            if (have_validation) {
                flint::Fmpz required_bound;
                if (class_group.factor_base_generation_bound(
                            flint::FmpzRef(required_bound)) &&
                    class_group.check_factor_base_generation_bound(
                            flint::FmpzConstRef(required_bound)) &&
                    class_group.try_certify_class_unit_with_units(
                            units, analytic_cache.zeta_validation_value(),
                            precision)) {
                    return true;
                }
            }
            // An unavailable enclosure or an inconclusive index-one test must
            // not weaken certification.  Retain the former strict BF audit as
            // the fail-closed producer below.
        }

        bool has_bf_audit = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_zeta_bf_audit");
            has_bf_audit = analytic_cache.ensure_bf_audit(
                    order, options.zeta_bf_max_cutoff, precision,
                    diagnostics, class_group.factor_base());
        }
        if (!has_bf_audit) {
            return false;
        }

        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.validation_zeta_bf_certify");
        return ClassGroupCertificationAccess::
                try_certify_class_unit_with_bf_audit(
                        class_group, units, analytic_cache.bf_value(),
                        analytic_cache.bf_error_bound(),
                        analytic_cache.bf_cutoff(),
                        options.zeta_bf_max_cutoff, precision,
                        analytic_cache.bf_work_precision());
    }
    return class_group.try_certify_class_unit_with_zeta(units, precision);
}

bool try_validate_candidate_pair(ClassGroupContext& class_group,
                                 OrderUnitGroup& units,
                                 const Order& order,
                                 const ClassGroupComputeOptions& options,
                                 AnalyticClassRegulatorCache& analytic_cache,
                                 slong precision) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_candidate_pair");
    // Class-group-only imaginary-quadratic computation already certifies
    // through FLINT reduced forms.  Combined computation defers class proof
    // until its exact rank-zero unit context exists, so use the exact index
    // instead of a generic zeta bound; a nontrivial index still drives the
    // existing exact relation-saturation loop.
    if (order.degree() == 2) {
        flint::Fmpz exact_index_bound;
        if (ClassGroupCertificationAccess::
                    rank_zero_quadratic_class_index_bound(
                            flint::FmpzRef(exact_index_bound), class_group,
                            units, options.requested_certification)) {
            return flint::fmpz_is_one(
                    flint::FmpzConstRef(exact_index_bound));
        }
    }
    const bool prefer_zeta_bf =
            options.requested_certification == CertificationMode::proven &&
            options.zeta_bf_max_cutoff != 0 && order.degree() > 2 &&
            units.is_set() && units.free_rank() > 1;
    const bool prefer_validation =
            options.requested_certification == CertificationMode::proven &&
            options.zeta_bf_max_cutoff != 0 && order.degree() > 2 &&
            analytic_cache.validation_enabled();
    if ((prefer_zeta_bf || prefer_validation) &&
        try_certify_candidate_with_zeta(class_group, units, order, options,
                                        analytic_cache, precision)) {
        return true;
    }

    if (class_group.try_certify_with_units(
                units, options.requested_certification, precision,
                options.zeta_bf_max_cutoff)) {
        return true;
    }

    return try_certify_candidate_with_zeta(class_group, units, order, options,
                                           analytic_cache, precision);
}

bool try_validate_candidate_pair_after_progress(
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        const Order& order,
        const ClassGroupComputeOptions& options,
        AnalyticClassRegulatorCache& analytic_cache,
        slong precision) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_after_progress");
    if (options.requested_certification != CertificationMode::proven) {
        return false;
    }

    if (order.degree() == 2) {
        flint::Fmpz exact_index_bound;
        if (ClassGroupCertificationAccess::
                    rank_zero_quadratic_class_index_bound(
                            flint::FmpzRef(exact_index_bound), class_group,
                            units, options.requested_certification)) {
            return flint::fmpz_is_one(
                    flint::FmpzConstRef(exact_index_bound));
        }
    }

    if (!analytic_cache.ensure(order, precision, diagnostics,
                               class_group.factor_base())) {
        return false;
    }

    flint::Fmpz aux_bound;
    relation_saturation_aux_bound(aux_bound, options);
    if (class_group.try_analytic_index_bound_with_units(
                units, analytic_cache.value(), flint::FmpzConstRef(aux_bound),
                precision)) {
        return true;
    }

    flint::Fmpz retry_aux_bound;
    if (!relation_saturation_retry_aux_bound(retry_aux_bound, options)) {
        return false;
    }

    return class_group.try_analytic_index_bound_with_units(
            units, analytic_cache.value(), flint::FmpzConstRef(retry_aux_bound),
            precision);
}

const char* validation_recompute_profile_label(
        ValidationRecomputeCause cause) noexcept {
    switch (cause) {
    case ValidationRecomputeCause::initial_requested_proven:
        return "unit_group.validation_recompute.initial_requested_proven";
    case ValidationRecomputeCause::post_compute_initial_candidate:
        return "unit_group.validation_recompute.post_compute_initial_candidate";
    case ValidationRecomputeCause::missing_candidate_units:
        return "unit_group.validation_recompute.missing_candidate_units";
    case ValidationRecomputeCause::validation_loop_unit_refinement:
        return "unit_group.validation_recompute.validation_loop_unit_refinement";
    case ValidationRecomputeCause::proven_refinement_retry:
        return "unit_group.validation_recompute.proven_refinement_retry";
    }
    return "unit_group.validation_recompute.unknown";
}

bool validation_uses_class_relation_units(
        const ValidationUnitRefreshOptions& refresh_options) noexcept {
    return refresh_options.use_class_relation_units &&
           refresh_options.relation_unit_state != nullptr &&
           refresh_options.unit_add != nullptr &&
           (refresh_options.rank > 1 ||
            (refresh_options.require_source_units &&
             refresh_options.rank > 0));
}

bool recompute_units_from_class_context(OrderUnitGroup& units,
                                        const Order& order,
                                        const ClassGroupContext& class_group,
                                        EmbeddingContext& embeddings,
                                        AnalyticClassRegulatorCache&
                                                analytic_cache,
                                        bool prove_index_bound,
                                        ValidationRecomputeCause cause,
                                        slong precision,
                                        const ValidationUnitRefreshOptions&
                                                refresh_options) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_recompute_units");
    const char* const cause_label = validation_recompute_profile_label(cause);
    (void) cause_label;
    SILEX_PROFILE_EVENT(units.diagnostics(), DiagnosticsModule::unit_group,
                        cause_label);
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        cause_label);
    SILEX_PROFILE_EVENT(
            units.diagnostics(), DiagnosticsModule::unit_group,
            prove_index_bound
                    ? "unit_group.validation_recompute.proof_requested"
                    : "unit_group.validation_recompute.proof_skipped");
    if (validation_uses_class_relation_units(refresh_options)) {
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.class_relation_units");
        bool ready = false;
        slong improved = 0;
        if (!set_class_relation_units(
                    ready, improved, units, order, class_group, embeddings,
                    analytic_cache, *refresh_options.relation_unit_state,
                    precision, refresh_options.rank,
                    *refresh_options.unit_add,
                    refresh_options.validation_bf_max_cutoff)) {
            SILEX_PROFILE_EVENT(
                    units.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.result.failed");
            return false;
        }
        if (*refresh_options.unit_add == WORD_MAX) {
            SILEX_PROFILE_EVENT(
                    units.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.result.add_exhausted");
            return false;
        }
        ++(*refresh_options.unit_add);
        if (refresh_options.unit_improved != nullptr) {
            *refresh_options.unit_improved = improved;
        }
        if (!ready) {
            SILEX_PROFILE_EVENT(
                    units.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.result.incomplete");
            return false;
        }
        if (prove_index_bound) {
            // reference `_class_unit_group` tries finite-index saturation after
            // `_validate_class_unit_group` returns an index bound > 1.  The
            // reference class-relation unit refresh above supplies the tentative
            // full-rank unit context; reuse the native proof/saturation route
            // to certify or carry any accepted roots forward.
            OrderUnitGroup proven(order);
            proven.set_diagnostics(units.diagnostics());
            ProofState status = ProofState::not_checked;
            bool proof_changed = false;
            flint::Fmpz proof_aux_bound;
            flint::fmpz_set_si(flint::FmpzRef(proof_aux_bound),
                               kComputeProofAuxMax);
            if (!proven.is_defined() ||
                !proven.prove_index_bound(
                        status, proof_changed, units, kComputeSatAuxTarget,
                        flint::FmpzConstRef(proof_aux_bound),
                        kComputeSatMaxPasses, embeddings, precision)) {
                SILEX_PROFILE_EVENT(
                        units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_recompute.result.proof_unavailable");
                return true;
            }
            units.swap(proven);
            SILEX_PROFILE_EVENT(
                    units.diagnostics(), DiagnosticsModule::unit_group,
                    status == ProofState::verified
                            ? "unit_group.validation_recompute.result.proven"
                            : "unit_group.validation_recompute.result.proof_unavailable");
            return true;
        }
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.ready");
        return true;
    }

    flint::Fmpz aux_start;
    flint::Fmpz aux_max;
    bool changed = false;
    bool stable = false;
    flint::fmpz_set_si(flint::FmpzRef(aux_start), kComputeSatAuxStart);
    flint::fmpz_set_si(flint::FmpzRef(aux_max), kComputeSatAuxMax);

    OrderUnitGroup refined(order);
    refined.set_diagnostics(units.diagnostics());
    if (!refined.is_defined()) {
        SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "relation-kernel unit refinement output unavailable");
        return false;
    }
    if (!refined.set_relation_kernel_units_index_bounded_saturated(
                changed, stable, order, class_group, embeddings, precision,
                precision, kComputeSatAuxTarget,
                flint::FmpzConstRef(aux_start),
                flint::FmpzConstRef(aux_max), kComputeSatMaxPasses)) {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
        if (log_enabled(units.diagnostics(), DiagnosticsModule::unit_group,
                        LogLevel::detail)) {
            std::string message = "relation-kernel unit refinement failed";
            message += " relations=";
            message += std::to_string(static_cast<long long>(
                    class_group.relation_count()));
            message += " kernel_units=";
            message += std::to_string(static_cast<long long>(
                    class_group.relation_kernel_unit_count()));
            SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail, message.c_str());
        }
#endif
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.refinement_failed");
        return false;
    }

    if (!prove_index_bound) {
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.refined_unproved");
        units.swap(refined);
        return true;
    }

    OrderUnitGroup proven(order);
    proven.set_diagnostics(units.diagnostics());
    ProofState status = ProofState::not_checked;
    bool proof_changed = false;
    flint::Fmpz proof_aux_bound;
    flint::fmpz_set_si(flint::FmpzRef(proof_aux_bound), kComputeProofAuxMax);
    if (!proven.is_defined()) {
        SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "relation-kernel unit proof output unavailable");
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.proof_output_failed");
        return false;
    }

    if (!proven.prove_index_bound(
                status, proof_changed, refined, kComputeSatAuxTarget,
                flint::FmpzConstRef(proof_aux_bound),
                kComputeSatMaxPasses, embeddings, precision)) {
        // Match C `_nf_ord_unit_group_prove_index_bound_compute`: once the
        // relation-kernel subgroup exists, an unsuccessful bounded proof
        // attempt is proof-unavailable, not a hard failure of the class/unit
        // validation loop.
        SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "relation-kernel unit proof unavailable");
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.proof_unavailable");
        units.swap(refined);
        return true;
    }

    units.swap(proven);
    if (status != ProofState::verified) {
        // Source trace: C `_nf_ord_unit_group_prove_index_bound_compute`
        // publishes the bounded proof attempt output whenever the call
        // succeeds, even if the proof status is unavailable.  This preserves
        // accepted saturation roots for the next validation pass.
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_recompute.result.proof_unavailable");
        return true;
    }

    SILEX_PROFILE_EVENT(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_recompute.result.proven");
    return true;
}

bool complete_requested_proven_relation_saturation(
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        const Order& order,
        const ClassGroupComputeOptions& options,
        AnalyticClassRegulatorCache& analytic_cache,
        slong precision) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.validation_relation_saturation_completion");
    if (options.requested_certification != CertificationMode::proven ||
        class_group.relation_saturation_status() == ProofState::verified) {
        return true;
    }
    constexpr ulong kAnalyticProofSkipMinBfCutoff = 1024;
    if (class_group.certification_status() == CertificationMode::proven &&
        units.certification_status() == CertificationMode::proven &&
        class_group.analytic_class_regulator_status() ==
                ProofState::verified &&
        class_group.factor_base_generation_checked_status() !=
                ProofState::verified &&
        (analytic_cache.has_bf_audit(options.zeta_bf_max_cutoff, precision) ||
         options.zeta_bf_max_cutoff >= kAnalyticProofSkipMinBfCutoff)) {
        return true;
    }
    if (class_group.analytic_class_regulator_status() !=
                ProofState::verified &&
        class_group.zeta_bf_proof_status() != ProofState::verified) {
        return true;
    }
    if (class_group.certification_status() != CertificationMode::proven ||
        units.certification_status() != CertificationMode::proven) {
        return false;
    }
    return try_validate_candidate_pair_after_progress(
            class_group, units, order, options, analytic_cache, precision);
}

bool try_validate_refine_loop(ClassGroupContext& class_group,
                              OrderUnitGroup& units,
                              OrderUnitGroup& scratch_units,
                              const Order& order,
                              const ClassGroupComputeOptions& options,
                              EmbeddingContext& embeddings,
                              AnalyticClassRegulatorCache& analytic_cache,
                              ValidateRefineSummary& summary,
                              slong precision,
                              bool units_from_current_class_context,
                              const ValidationUnitRefreshOptions&
                                      refresh_options) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_refine_loop");
    summary.reset();
    flint::Fmpz index_bound;
    flint::Fmpz next_index_bound;
    flint::Fmpz previous_index_bound;
    bool have_previous = false;
    bool progress_since_previous = false;
    bool class_progress_since_unit_refresh = false;

    for (slong pass = 0; pass < kComputeSatMaxPasses; ++pass) {
        const bool requested_grh =
                options.requested_certification == CertificationMode::grh;
        if (!requested_grh &&
            try_validate_candidate_pair(class_group, units, order, options,
                                        analytic_cache, precision)) {
            summary.outcome = ValidateRefineOutcome::proven;
            return true;
        }

        if (!requested_grh &&
            options.requested_certification != CertificationMode::proven) {
            summary.outcome = ValidateRefineOutcome::not_proven_request;
            return false;
        }

        bool exact_quadratic_index = false;
        if (order.degree() == 2) {
            exact_quadratic_index = ClassGroupCertificationAccess::
                    rank_zero_quadratic_class_index_bound(
                            flint::FmpzRef(index_bound), class_group, units,
                            options.requested_certification);
        }
        if (!exact_quadratic_index &&
            !analytic_index_bound_for_validation(
                    index_bound, units, class_group, analytic_cache, order,
                    precision)) {
            summary.outcome = ValidateRefineOutcome::analytic_unavailable;
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(summary.last_index_bound),
                        flint::FmpzConstRef(index_bound));

        if (flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) <= 0) {
            if (exact_quadratic_index || requested_grh) {
                // reference `_class_unit_group` accepts analytic index one for a
                // GRH request and skips the later unconditional class/unit
                // proof passes.  The caller publishes the conditional labels.
                summary.outcome = ValidateRefineOutcome::proven;
                return true;
            }
            // Source trace: reference `Clgp.jl:_class_unit_group` finishes the
            // tentative class/unit loop when `_validate_class_unit_group`
            // returns index one.  Native validation computed the same
            // analytic index-one result above; try to publish through the
            // class/unit certification helper before reporting publication
            // failure.
            if (class_group.try_certify_class_unit_with_units(
                        units, analytic_cache.value(), precision)) {
                summary.outcome = ValidateRefineOutcome::proven;
                return true;
            }
            summary.outcome =
                    ValidateRefineOutcome::index_one_publication_failed;
            return false;
        }

        if (!validate_refine_continue_with_bound(
                    have_previous, progress_since_previous,
                    flint::FmpzConstRef(previous_index_bound),
                    flint::FmpzConstRef(index_bound))) {
            summary.outcome = ValidateRefineOutcome::no_progress;
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(previous_index_bound),
                        flint::FmpzConstRef(index_bound));
        have_previous = true;
        progress_since_previous = false;

        bool class_changed = false;
        if (!saturate_candidate_class_relations_with_units(
                    class_changed, class_group, units, options,
                    flint::FmpzConstRef(index_bound),
                    exact_quadratic_index)) {
            summary.outcome =
                    ValidateRefineOutcome::class_saturation_failed;
            return false;
        }
        if (class_changed) {
            progress_since_previous = true;
            summary.class_progress = true;
            class_progress_since_unit_refresh = true;
            if (!requested_grh &&
                try_validate_candidate_pair_after_progress(
                        class_group, units, order, options, analytic_cache,
                        precision)) {
                summary.outcome = ValidateRefineOutcome::proven;
                return true;
            }
            summary.outcome = ValidateRefineOutcome::local_rank_unavailable;
            continue;
        }

        const bool independently_proven_quadratic_unit =
                order.degree() == 2 && units.free_rank() == 1 &&
                !units_from_current_class_context &&
                units.certification_status() == CertificationMode::proven;
        if (independently_proven_quadratic_unit) {
            // A proven full unit group is independent of the tentative class
            // relation lattice.  A remaining index therefore asks for more
            // class relations; rebuilding relation-kernel units cannot refine
            // this unit group.
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.proven_units_skip");
            summary.outcome = ValidateRefineOutcome::no_progress;
            return false;
        }

        const bool defer_validation_unit_proof =
                options.zeta_bf_max_cutoff != 0 && order.degree() > 2;
        if (defer_validation_unit_proof) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.validation_loop_proof_deferred");
        }

        if (units_from_current_class_context &&
            !class_progress_since_unit_refresh) {
            // reference keeps the current UnitGrpCtx and only changes it through
            // saturation/root acceptance.  If this native validation pass has
            // not changed the class context since the units were built from it,
            // rebuilding the same relation-kernel units is redundant.
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_recompute.validation_loop_fresh_context_skip");
            summary.outcome = ValidateRefineOutcome::no_progress;
            return false;
        }

        if (validation_uses_class_relation_units(refresh_options) &&
            units.is_set() && !scratch_units.set(units)) {
            summary.outcome = ValidateRefineOutcome::unit_refinement_failed;
            return false;
        }

        if (!recompute_units_from_class_context(
                    scratch_units, order, class_group, embeddings,
                    analytic_cache,
                    !requested_grh && !defer_validation_unit_proof,
                    ValidationRecomputeCause::validation_loop_unit_refinement,
                    precision, refresh_options)) {
            summary.outcome =
                    ValidateRefineOutcome::unit_refinement_failed;
            return false;
        }

        bool next_check_available = analytic_index_bound_for_validation(
                next_index_bound, scratch_units, class_group, analytic_cache,
                order, precision);
        if (next_check_available) {
            flint::fmpz_set(flint::FmpzRef(summary.last_index_bound),
                            flint::FmpzConstRef(next_index_bound));
        }

        units.swap(scratch_units);
        class_progress_since_unit_refresh = false;
        const bool candidate_validated = !requested_grh &&
                (defer_validation_unit_proof
                         ? try_validate_candidate_pair(
                                   class_group, units, order, options,
                                   analytic_cache, precision)
                         : try_validate_candidate_pair_after_progress(
                                   class_group, units, order, options,
                                   analytic_cache, precision));
        if (candidate_validated) {
            summary.outcome = ValidateRefineOutcome::proven;
            return true;
        }

        if (next_check_available &&
            validate_refine_has_progress(
                    false, true, flint::FmpzConstRef(index_bound),
                    flint::FmpzConstRef(next_index_bound))) {
            progress_since_previous = true;
            summary.unit_progress = true;
            summary.outcome = ValidateRefineOutcome::local_rank_unavailable;
            continue;
        }

        summary.outcome = next_check_available
                ? ValidateRefineOutcome::no_progress
                : ValidateRefineOutcome::analytic_unavailable;
        return false;
    }

    summary.outcome = ValidateRefineOutcome::pass_cap_reached;
    return false;
}

bool try_validate_refine_extra_pass(slong& extra_passes_done,
                                    ClassGroupContext& class_group,
                                    OrderUnitGroup& units,
                                    OrderUnitGroup& scratch_units,
                                    const Order& order,
                                    const ClassGroupComputeOptions& options,
                                    EmbeddingContext& embeddings,
                                    AnalyticClassRegulatorCache& analytic_cache,
                                    ValidateRefineSummary& summary,
                                    slong precision,
                                    bool units_from_current_class_context,
                                    const ValidationUnitRefreshOptions&
                                            refresh_options) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_extra_pass");
    if (!validate_refine_extra_pass_allowed(
                summary, options, extra_passes_done,
                kComputeRevalidateMaxPasses)) {
        return false;
    }

    ++extra_passes_done;
    return try_validate_refine_loop(class_group, units, scratch_units, order,
                                    options, embeddings, analytic_cache,
                                    summary, precision,
                                    units_from_current_class_context,
                                    refresh_options);
}

}  // namespace silex::detail
