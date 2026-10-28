#include <silex/order_unit.hpp>

#include "compute_internal.hpp"
#include "order_unit_internal.hpp"
#include "../class_group/class_group_certification_internal.hpp"
#include "../zeta/zeta_internal.hpp"

#include <silex/class_group.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include <cstdio>
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
        ZetaBfResidueDegreeCache* residue_degree_cache,
        bool* unconditional,
        ZetaBfRouteAudit* bf_audit,
        bool accept_grh_residue) noexcept {
    return detail::zeta_class_regulator_product_with_diagnostics(
            out, order, precision, diagnostics, residue_degree_base,
            residue_degree_cache, unconditional, bf_audit,
            accept_grh_residue);
}

bool index_one_class_regulator_product_for_validation(
        flint::ArbRef out,
        ZetaBfRouteAudit& audit,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache,
        ulong log_radius_denominator) noexcept {
    return detail::zeta_class_regulator_product_index_one_with_diagnostics(
            out, audit, order, precision, diagnostics, residue_degree_base,
            residue_degree_cache, log_radius_denominator);
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

// `accept_grh_generation` lets GRH-conditional factor-base generation
// (record_grh_factor_base_generation) stand in for the unconditional
// generation check.  Only a grh request may pass true for an acceptance
// decision; that record is written only on a grh request's base.
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
std::string log_fmpz_text(flint::FmpzConstRef value) noexcept {
    char* const text = ::fmpz_get_str(nullptr, 10, value.raw());
    std::string out = text != nullptr ? text : "?";
    if (text != nullptr) {
        ::flint_free(text);
    }
    return out;
}

std::string log_arb_mid_rad(flint::ArbConstRef value) noexcept {
    char buffer[96];
    std::snprintf(buffer, sizeof(buffer), "mid=%.17g rad=%.6e",
                  ::arf_get_d(arb_midref(value.raw()), ARF_RND_NEAR),
                  ::mag_get_d(arb_radref(value.raw())));
    return std::string(buffer);
}
#endif

bool validation_integer_index_prerequisites(
        const ClassGroupContext& class_group,
        bool accept_grh_generation) noexcept {
    return ClassGroupCertificationAccess::
                   factor_base_generation_checked_to_required_bound(
                           class_group) ||
           (accept_grh_generation &&
            ClassGroupCertificationAccess::
                    grh_factor_base_generation_covered(class_group));
}

bool validation_index_bound_from_product(
        flint::Fmpz& index_bound,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision,
        bool accept_grh_generation) noexcept {
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

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (log_enabled(units.diagnostics(), DiagnosticsModule::unit_group,
                    LogLevel::detail)) {
        std::string message = "index_bound_product candidate_hR ";
        message += log_arb_mid_rad(
                flint::ArbConstRef(candidate_class_regulator_product));
        message += " analytic_hR ";
        message += log_arb_mid_rad(analytic_class_regulator_product);
        SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail, message.c_str());
    }
#endif

    const bool integer_index_prerequisites =
            validation_integer_index_prerequisites(class_group,
                                                   accept_grh_generation);
    // In this validation context the quotient is an integer class/unit
    // index only once the factor base generates the class group: proved
    // to the proven generation requirement (the Minkowski-type bound, or
    // Zimmert's twin-class bound where it applies), or, for a grh request,
    // assumed under GRH from a base containing every prime ideal up to the
    // GRH bound (Bach 1990, Thm 4 and the remark following it, p. 376;
    // Belabas-Diaz y Diaz-Friedman 2008, Cor 2.2).  Keep the public ratio
    // accessor conservative, but recognize the source-complete index-one
    // case before taking an outward-rounded ceiling.
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

enum class CoarseIndexOneOutcome {
    no_value,
    index_one,
    index_above_one,
    no_index_bound,
};

// The grh coarse index-one test of analytic_index_bound_for_validation on
// its own: the Belabas-Friedman hR accurate enough to recognize index one
// (ensure_index_one), compared with the candidate h R.  It never computes
// the default hR.  `coarse_index_bound`, when given, receives the coarse
// value's index bound on index_one and index_above_one.
CoarseIndexOneOutcome coarse_index_one_outcome(
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision,
        flint::Fmpz* coarse_index_bound_out = nullptr) noexcept {
    bool have_index_one_value = false;
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound."
                "ensure_index_one_product");
        have_index_one_value = analytic_cache.ensure_index_one(
                order, precision, units.diagnostics(),
                class_group.factor_base());
    }
    if (!have_index_one_value) {
        return CoarseIndexOneOutcome::no_value;
    }
    flint::Fmpz coarse_index_bound;
    if (!validation_index_bound_from_product(
                coarse_index_bound, units, class_group,
                analytic_cache.value(), precision, true)) {
        return CoarseIndexOneOutcome::no_index_bound;
    }
    if (coarse_index_bound_out != nullptr) {
        flint::fmpz_set(flint::FmpzRef(*coarse_index_bound_out),
                        flint::FmpzConstRef(coarse_index_bound));
    }
    return flint::fmpz_is_one(flint::FmpzConstRef(coarse_index_bound))
            ? CoarseIndexOneOutcome::index_one
            : CoarseIndexOneOutcome::index_above_one;
}

}  // namespace

bool validation_index_enclosure_has_single_integer(
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision) noexcept {
    if (precision <= 0 || !units.is_set() ||
        !flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product)) {
        return false;
    }
    flint::Arb candidate_class_regulator_product;
    if (!units.class_regulator_product(
                flint::ArbRef(candidate_class_regulator_product),
                class_group, precision) ||
        !flint::arb_is_finite(candidate_class_regulator_product) ||
        !flint::arb_is_positive(candidate_class_regulator_product)) {
        return false;
    }
    // The quotient of class_regulator_index_bound_from_candidate_product,
    // whose upper endpoint's ceiling is the index bound.
    flint::Arb quotient;
    ::arb_div(quotient.raw(), candidate_class_regulator_product.raw(),
              analytic_class_regulator_product.raw(), precision);
    if (!flint::arb_is_finite(quotient) || !flint::arb_is_positive(quotient)) {
        return false;
    }
    flint::Arf lower;
    flint::Arf upper;
    flint::arb_get_lbound_arf(lower, quotient, precision);
    flint::arb_get_ubound_arf(upper, quotient, precision);
    if (!flint::arf_is_finite(lower) || !flint::arf_is_finite(upper)) {
        return false;
    }
    flint::Fmpz least;
    flint::Fmpz greatest;
    flint::arf_get_fmpz(least, lower, ARF_RND_CEIL);
    flint::arf_get_fmpz(greatest, upper, ARF_RND_FLOOR);
    return flint::fmpz_equal(flint::FmpzConstRef(least),
                             flint::FmpzConstRef(greatest));
}

namespace {

// The refined step of the grh index bound, between the coarse index-one
// value and the default hR: a Belabas-Friedman hR with log radius below
// 1/q, q = max(8, 2 B) for the coarse value's index bound B
// (ensure_index_one_refined), used only when the enclosure of
// h_cand R_cand / hR holds a single integer k.  With the generation
// prerequisites the quotient is the class/unit index, an integer, and
// under GRH both values enclose it (Belabas-Friedman 2015, Theorem 7 for
// this one, Theorem 1 for the default route), so it is k, and the index
// bound this value gives, the ceiling of the upper endpoint, is the one the
// default hR gives whenever that hR's upper endpoint also lies in
// (k, k + 1) (or both shortcut to index one).  q only schedules the work:
// with true index k <= B and relative width about 2/q the enclosure is at
// most about one unit wide.  Returns false, leaving the default hR to
// decide, when the prerequisites fail, B is too large for q, the value is
// unavailable (no cutoff up to the default maximum reaches 1/q), or the
// enclosure holds two or more integers or none.
bool refined_index_bound_for_validation(
        flint::Fmpz& index_bound,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision,
        const flint::Fmpz& coarse_index_bound) noexcept {
    constexpr ulong kMinDenominator = 2 * kZetaBfIndexOneLogRadiusDenominator;
    if (!validation_integer_index_prerequisites(class_group, true) ||
        flint::fmpz_sgn(flint::FmpzConstRef(coarse_index_bound)) <= 0 ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(coarse_index_bound),
                           UWORD_MAX / 2) > 0) {
        return false;
    }
    const ulong denominator = FLINT_MAX(
            kMinDenominator,
            2 * flint::fmpz_get_ui(flint::FmpzConstRef(coarse_index_bound)));
    bool have_refined_value = false;
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound."
                "ensure_refined_product");
        have_refined_value = analytic_cache.ensure_index_one_refined(
                order, precision, denominator, units.diagnostics(),
                class_group.factor_base());
    }
    if (!have_refined_value) {
        return false;
    }
    if (!validation_index_enclosure_has_single_integer(
                units, class_group, analytic_cache.value(), precision)) {
        SILEX_PROFILE_EVENT(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound."
                "refined_product_ambiguous");
        return false;
    }
    return validation_index_bound_from_product(
            index_bound, units, class_group, analytic_cache.value(),
            precision, true);
}

}  // namespace

GrhIndexOneTest grh_validation_index_one_test(
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.grh_validation_index_one_test");
    switch (coarse_index_one_outcome(units, class_group, analytic_cache,
                                     order, precision)) {
    case CoarseIndexOneOutcome::index_one:
        return GrhIndexOneTest::index_one;
    case CoarseIndexOneOutcome::index_above_one:
        // Without the generation prerequisites the quotient is not an
        // integer index, and the coarse value never shows index one
        // (validation_index_bound_from_product).
        return validation_integer_index_prerequisites(class_group, true)
                ? GrhIndexOneTest::not_index_one
                : GrhIndexOneTest::unavailable;
    case CoarseIndexOneOutcome::no_value:
    case CoarseIndexOneOutcome::no_index_bound:
        break;
    }
    return GrhIndexOneTest::unavailable;
}

bool analytic_index_bound_for_validation(
        flint::Fmpz& index_bound,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision,
        bool accept_grh_generation) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_analytic_index_bound");
    // A grh request first tries a Belabas-Friedman hR that is only accurate
    // enough to recognize index one (log radius below 1/4; Belabas-Friedman
    // 2015, p. 1, footnote 1, and Theorem 7 for the cutoff).  Index one
    // from it is accepted as is.  When it shows an index above one, a
    // narrower value of the same route gives the bound if its enclosure
    // holds a single integer (refined_index_bound_for_validation).  Any
    // other outcome is recomputed below from the default hR, so every
    // failure comes from the same value as before.
    if (accept_grh_generation) {
        flint::Fmpz coarse_index_bound;
        const CoarseIndexOneOutcome coarse = coarse_index_one_outcome(
                units, class_group, analytic_cache, order, precision,
                &coarse_index_bound);
        if (coarse == CoarseIndexOneOutcome::index_one) {
            flint::fmpz_one(flint::FmpzRef(index_bound));
            return true;
        }
        if (coarse != CoarseIndexOneOutcome::no_value) {
            SILEX_PROFILE_EVENT(
                    units.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.validation_analytic_index_bound."
                    "index_one_product_refined");
        }
        if (coarse == CoarseIndexOneOutcome::index_above_one &&
            refined_index_bound_for_validation(
                    index_bound, units, class_group, analytic_cache, order,
                    precision, coarse_index_bound)) {
            return true;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound.ensure_product");
        // `accept_grh_generation` is true exactly for a grh request, which
        // also accepts a GRH-conditional hR: take the Belabas-Friedman
        // residue there instead of the quadratic L(1, chi) sum.
        if (!analytic_cache.ensure(order, precision, units.diagnostics(),
                                   class_group.factor_base(), false,
                                   accept_grh_generation)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                units.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound.index_bound");
        return validation_index_bound_from_product(
                index_bound, units, class_group, analytic_cache.value(),
                precision, accept_grh_generation);
    }
}

bool class_unit_validation_estimate(
        flint::Fmpz& index_bound,
        flint::Arb& expected_regulator,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision,
        bool accept_grh_generation) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.validation_estimate");
    // `accept_grh_generation` is true exactly for a grh request, whose base
    // alone can carry the GRH generation record.  An estimated index of one
    // can become an acceptance (unknown_candidate_validated), so the record
    // counts as generation only inside a grh request.
    if (!analytic_index_bound_for_validation(
                index_bound, units, class_group, analytic_cache, order,
                precision, accept_grh_generation) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) < 0) {
        return false;
    }
    return validation_expected_regulator(expected_regulator, units,
                                         index_bound, precision);
}

bool validation_expected_regulator(flint::Arb& expected_regulator,
                                   const OrderUnitGroup& units,
                                   const flint::Fmpz& index_bound,
                                   slong precision) noexcept {
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
        const DiagnosticsContext* diagnostics,
        bool accept_grh_generation) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.bf_validation_estimate");
    if (max_cutoff == 0 || order.degree() <= 2) {
        return class_unit_validation_estimate(
                index_bound, expected_regulator, units, class_group,
                analytic_cache, order, precision,
                accept_grh_generation);
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
                    analytic_cache, order, precision,
                    accept_grh_generation);
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
                    analytic_cache, order, precision,
                    accept_grh_generation);
        }
        analytic_cache.seed(order, analytic_cache.bf_value(),
                            analytic_cache.bf_work_precision());
    }

    const flint::ArbConstRef validation_value = use_validation
            ? analytic_cache.zeta_validation_value()
            : analytic_cache.bf_value();
    if (!validation_index_bound_from_product(
                index_bound, units, class_group, validation_value,
                precision, accept_grh_generation) ||
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

namespace {

// True when each of the two positive regulators is certified to be below
// twice the other.  When one unit lattice contains the other, the quotient
// of their regulators is the integer index, so the lattices are then equal.
bool regulators_within_factor_two(const flint::Arb& left,
                                  const flint::Arb& right) noexcept {
    if (!flint::arb_is_finite(left) || !flint::arb_is_finite(right) ||
        !flint::arb_is_positive(left) || !flint::arb_is_positive(right)) {
        return false;
    }
    flint::Arb twice_left;
    flint::Arb twice_right;
    flint::arb_mul_2exp_si(twice_left, left, 1);
    flint::arb_mul_2exp_si(twice_right, right, 1);
    return flint::arb_lt(left, twice_right) &&
           flint::arb_lt(right, twice_left);
}

}  // namespace

bool saturate_candidate_class_relations_with_units(
        bool& changed,
        ClassGroupContext& class_group,
        const OrderUnitGroup& units,
        const ClassGroupComputeOptions& options,
        flint::FmpzConstRef index_bound,
        bool index_bound_is_exact,
        RelationSaturationMemo* memo) noexcept {
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

    // Extension before the first bounded pass (see
    // RelationSaturationMemo::defer_until_extension_stalls).  While the
    // deferral is pending every bounded pass is deferred; the pending pass
    // runs only through run_deferred_class_relation_saturation, when an
    // extension slice fails or appends nothing, or before acceptance on the
    // relations of the last deferral.
    //
    // A stall alone (an extension slice that left the index bound at or
    // above the deferred bound) does not run the pass: that slice appended
    // relations, and a pass on them is deferred like the first one.  A stall
    // test that also asked for no relation appended since the deferral could
    // never hold, so there is none.  In compute_with_relation_class_group_,
    // extension_rounds counts only slices that end with more relations than
    // they started with; a slice that does not (a factor-base rebuild may
    // drop replayed relations) runs the pending pass and ends the deferral;
    // and outside extension slices the loop only appends relations to its
    // candidate class group.  So a new round since the deferral implies a
    // relation count above the one recorded at the deferral.
    if (memo != nullptr && !index_bound_is_exact &&
        memo->defer_until_extension_stalls && !memo->deferral_done) {
        // The bound of the pending pass: the current one after an
        // extension slice, else the smallest of the passes deferred since
        // that slice.
        if (!memo->deferred ||
            memo->extension_rounds != memo->deferred_extension_round ||
            flint::fmpz_cmp(index_bound,
                            flint::FmpzConstRef(
                                    memo->deferred_index_bound)) < 0) {
            flint::fmpz_set(flint::FmpzRef(memo->deferred_index_bound),
                            index_bound);
        }
        memo->deferred = true;
        memo->deferred_extension_round = memo->extension_rounds;
        memo->deferred_relation_count = class_group.relation_count();
        ++memo->deferrals;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.validation_class_relation_saturation."
                "deferred_to_extension");
        return true;
    }

    flint::Fmpz aux_bound;
    relation_saturation_aux_bound(aux_bound, options);
    // The exact imaginary-quadratic index h_cand / h bounds the saturation
    // index only once the factor base generates the class group.  The bound
    // is therefore unconditional exactly when generation is verified up to
    // the required Minkowski-type bound, whatever the request: a proven
    // request has checked it (rank_zero_quadratic_class_index_bound), and
    // a grh request may have checked it on another route; otherwise a grh
    // request assumes generation under GRH and the bound is conditional.
    // The other branch saturates at the primes of h_cand itself, an
    // unconditional bound; the analytic index bound only decides whether to
    // run it.
    const bool exact_index_bound_unconditional =
            ClassGroupCertificationAccess::
                    factor_base_generation_checked_to_required_bound(
                            class_group);
    bool saturated = false;
    const slong remaining = options.max_relations - class_group.relation_count();
    const slong max_appends_per_ell =
            relation_saturation_max_appends_per_ell(options);
    // The per-ell append limit of a pass that appends nothing.
    const slong attempt_limit = min_slong(max_appends_per_ell, remaining);
    if (memo != nullptr && !index_bound_is_exact && memo->valid) {
        flint::Fmpz class_order;
        flint::Arb regulator;
        if (class_group.order(flint::FmpzRef(class_order)) &&
            units.regulator(flint::ArbRef(regulator)) &&
            flint::fmpz_equal(flint::FmpzConstRef(class_order),
                              flint::FmpzConstRef(memo->class_order)) &&
            flint::fmpz_equal(flint::FmpzConstRef(aux_bound),
                              flint::FmpzConstRef(memo->aux_bound)) &&
            class_group.generator_count() == memo->generator_count &&
            max_appends_per_ell == memo->max_appends_per_ell &&
            attempt_limit == memo->attempt_limit &&
            regulators_within_factor_two(regulator, memo->regulator) &&
            ClassGroupCertificationAccess::
                    replay_unchanged_bounded_relation_saturation(
                            class_group, units,
                            flint::FmpzConstRef(class_order))) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_class_relation_saturation."
                    "memo_skip");
            ++memo->skips;
            return true;
        }
    }
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    const bool log_saturation = log_enabled(
            diagnostics, DiagnosticsModule::unit_group, LogLevel::detail);
    flint::Fmpz h_cand_before;
    const bool have_h_before = log_saturation
            && class_group.order(flint::FmpzRef(h_cand_before));
#endif
    const bool saturation_ok = index_bound_is_exact
            ? ClassGroupCertificationAccess::
                      saturate_relations_for_index_bound_with_units(
                              changed, saturated, class_group, units,
                              index_bound, exact_index_bound_unconditional,
                              flint::FmpzConstRef(aux_bound),
                              relation_saturation_max_appends_per_ell(options),
                              remaining)
            : class_group.saturate_relations_bounded_with_units(
                      changed, saturated, units,
                      flint::FmpzConstRef(aux_bound),
                      max_appends_per_ell, remaining);
    if (!saturation_ok) {
        return false;
    }
    if (memo != nullptr && !index_bound_is_exact) {
        memo->valid = !changed && !saturated &&
                class_group.order(flint::FmpzRef(memo->class_order)) &&
                units.regulator(flint::ArbRef(memo->regulator));
        if (memo->valid) {
            flint::fmpz_set(flint::FmpzRef(memo->aux_bound),
                            flint::FmpzConstRef(aux_bound));
            memo->generator_count = class_group.generator_count();
            memo->max_appends_per_ell = max_appends_per_ell;
            memo->attempt_limit = attempt_limit;
        }
    }

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (log_saturation) {
        flint::Fmpz h_cand_after;
        const bool have_h_after = class_group.order(flint::FmpzRef(h_cand_after));
        std::string message = "saturation changed=";
        message += changed ? "1" : "0";
        message += " h_cand_before=";
        message += have_h_before ? log_fmpz_text(h_cand_before) : "?";
        message += " h_cand_after=";
        message += have_h_after ? log_fmpz_text(h_cand_after) : "?";
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail, message.c_str());
    }
#endif

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

bool run_deferred_class_relation_saturation(
        bool& changed,
        ClassGroupContext& class_group,
        const OrderUnitGroup& units,
        const ClassGroupComputeOptions& options,
        RelationSaturationMemo& memo,
        bool before_acceptance) noexcept {
    changed = false;
    if (!memo.defer_until_extension_stalls || memo.deferral_done ||
        !memo.deferred) {
        return true;
    }
    memo.deferral_done = true;
    memo.deferred = false;
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    if (before_acceptance) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.validation_class_relation_saturation."
                "deferred_run_before_acceptance");
    } else {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.validation_class_relation_saturation."
                "deferred_run_after_failed_extension");
    }
    return saturate_candidate_class_relations_with_units(
            changed, class_group, units, options,
            flint::FmpzConstRef(memo.deferred_index_bound), false, &memo);
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
                    ClassGroupCertificationAccess::
                            try_certify_class_unit_with_units(
                                    class_group, units,
                                    analytic_cache.zeta_validation_value(),
                                    precision, false)) {
                    // The validation enclosure is a Belabas-Friedman value
                    // (GRH-conditional): the gate publishes `proven` only
                    // when saturation has already been proven.
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
    if (ClassGroupCertificationAccess::try_analytic_index_bound_with_units(
                class_group, units, analytic_cache.value(),
                flint::FmpzConstRef(aux_bound), precision,
                analytic_cache.value_unconditional())) {
        return true;
    }

    flint::Fmpz retry_aux_bound;
    if (!relation_saturation_retry_aux_bound(retry_aux_bound, options)) {
        return false;
    }

    return ClassGroupCertificationAccess::try_analytic_index_bound_with_units(
            class_group, units, analytic_cache.value(),
            flint::FmpzConstRef(retry_aux_bound), precision,
            analytic_cache.value_unconditional());
}

namespace {

bool try_prove_candidate_pair_by_saturation(
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        const Order& order,
        const ClassGroupComputeOptions& options,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    const DiagnosticsContext* diagnostics =
            options.diagnostics != nullptr ? options.diagnostics
                                           : class_group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.validation_saturation_proof");
    if (options.requested_certification != CertificationMode::proven ||
        !units.is_set() || !same_order_parent(units.parent(), &order) ||
        !class_group.has_presentation()) {
        return false;
    }

    // A GRH-conditional (Belabas-Friedman) analytic index-one result only
    // selects this candidate.  The proof is unconditional: the unit group is
    // proven from the regulator lower bound and unit saturation, and the
    // class group from relation saturation at every p | h_cand.  The proven
    // unit group is built aside and replaces `units` only after the class
    // proof succeeds, so a failure leaves both objects unchanged.
    // The class saturation proof scans auxiliary primes q = 1 mod ell with
    // no fixed bound, as source buch3.c:primecertify does (see
    // kComputeProofClassAuxMax).  The unit proof below keeps its bounded
    // pre-scan and its own full scan.
    flint::Fmpz aux_bound;
    flint::fmpz_set_ui(flint::FmpzRef(aux_bound), kComputeProofClassAuxMax);
    if (units.certification_status() != CertificationMode::proven) {
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
                    kComputeSatMaxPasses, embeddings, precision) ||
            status != ProofState::verified ||
            proven.certification_status() != CertificationMode::proven) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.validation_saturation_proof.units_unavailable");
            return false;
        }
        // The published unit group is proven only if its torsion is the
        // torsion Silex computes for the order, never a cached or supplied
        // value.  prove_index_bound checks that on its input before any
        // work and fails closed otherwise (event
        // unit_group.prove_index_bound.torsion_not_computed), and a proven
        // result carries that torsion, so it is not recomputed here.  The
        // check happens before the class proof, so a rejection leaves both
        // objects unchanged.
        if (!ClassGroupCertificationAccess::
                    try_prove_class_order_saturation_with_units(
                            class_group, proven,
                            flint::FmpzConstRef(aux_bound))) {
            return false;
        }
        units.swap(proven);
        return true;
    }

    // Units already labelled proven were labelled by a route that computed
    // or checked their torsion (the torsion invariant in
    // order_unit_internal.hpp), so the label is trusted here.
    return ClassGroupCertificationAccess::
            try_prove_class_order_saturation_with_units(
                    class_group, units, flint::FmpzConstRef(aux_bound));
}

}  // namespace

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
                    refresh_options.validation_bf_max_cutoff,
                    refresh_options.accept_grh_generation)) {
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
            // prove_index_bound checks the torsion before it labels the
            // result proven, so the swapped-in group is proven only with
            // computed torsion.
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

    // prove_index_bound checks the torsion of `refined` before it labels
    // the result proven; by the torsion invariant (order_unit_internal.hpp)
    // the check passes for this internally built group.
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
    // Only an unconditional analytic proof may stand in for relation
    // saturation; a Belabas-Friedman check assumes GRH.
    if (class_group.certification_status() == CertificationMode::proven &&
        units.certification_status() == CertificationMode::proven &&
        class_group.analytic_class_regulator_certification() ==
                CertificationMode::proven &&
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

        // Under a grh request the exact imaginary-quadratic index is
        // available only for |D| < 2^34; above it h comes from the GRH
        // relation lattice through the analytic index test below, without
        // listing reduced forms (rank_zero_quadratic_class_index_bound).
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
                    precision, requested_grh)) {
            summary.outcome = ValidateRefineOutcome::analytic_unavailable;
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(summary.last_index_bound),
                        flint::FmpzConstRef(index_bound));

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
        if (log_enabled(units.diagnostics(), DiagnosticsModule::unit_group,
                        LogLevel::detail)) {
            std::string message = "refine_pass pass=";
            message += std::to_string(static_cast<long long>(pass));
            message += " index_bound=";
            message += log_fmpz_text(flint::FmpzConstRef(index_bound));
            SILEX_LOG(units.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail, message.c_str());
        }
#endif

        if (flint::fmpz_cmp_ui(flint::FmpzConstRef(index_bound), 1) <= 0) {
            // A bounded pass deferred behind relation extension
            // (RelationSaturationMemo::defer_until_extension_stalls) on the
            // relations validated here would have run before this index-one
            // result without the deferral, and its relation_saturation_
            // status() and ell-records stay with the accepted pair; run it
            // now.  After any relation append since the deferral the pass's
            // records would have been cleared, so nothing runs then.
            RelationSaturationMemo* const saturation_memo =
                    refresh_options.relation_saturation_memo;
            if (requested_grh && saturation_memo != nullptr &&
                saturation_memo->deferred &&
                !saturation_memo->deferral_done &&
                saturation_memo->deferred_extension_round ==
                        saturation_memo->extension_rounds &&
                saturation_memo->deferred_relation_count ==
                        class_group.relation_count()) {
                bool class_changed = false;
                if (!run_deferred_class_relation_saturation(
                            class_changed, class_group, units, options,
                            *saturation_memo, true)) {
                    summary.outcome =
                            ValidateRefineOutcome::class_saturation_failed;
                    return false;
                }
                if (class_changed) {
                    progress_since_previous = true;
                    summary.class_progress = true;
                    class_progress_since_unit_refresh = true;
                    summary.outcome =
                            ValidateRefineOutcome::local_rank_unavailable;
                    continue;
                }
            }
            if (exact_quadratic_index || requested_grh) {
                // reference `_class_unit_group` accepts analytic index one for a
                // GRH request and skips the later unconditional class/unit
                // proof passes.  The caller publishes the conditional labels.
                // As in the reference, index one here needs factor-base
                // generation, which a grh request may take under GRH from a
                // base containing every prime ideal up to min(BDF, Bach)
                // (validation_index_bound_from_product), or unconditionally
                // when the base reaches the proven generation requirement
                // (the Minkowski-type bound, or Zimmert's twin-class bound
                // where it applies).  GRH then
                // enters through that generation bound and a
                // Belabas-Friedman hR.
                if (!exact_quadratic_index) {
                    // Record the analytic check that accepted the pair, with
                    // its own conditionality.  Informational only; the
                    // labels stay `grh`.  The exact imaginary-quadratic
                    // index uses no analytic value and records nothing.
                    ClassGroupCertificationAccess::
                            record_grh_acceptance_analytic_check(
                                    class_group,
                                    analytic_cache.value_unconditional());
                    // When that hR was a Belabas-Friedman value, also record
                    // its audit data.  A grh request never enables the
                    // validation enclosure (validation_bf_max_cutoff_for_
                    // options), so the value is either the coarse Theorem 7
                    // index-one value that was accepted as index one, or the
                    // default zeta route's value; a seeded value has no
                    // audit data and records none.
                    // The record cannot fail here: the presentation exists,
                    // both precisions are positive, and the error bound is
                    // finite because the Belabas-Friedman value it was added
                    // to is finite (bf_log_residue_cutoff).
                    if (const ZetaBfRouteAudit* bf_audit =
                                analytic_cache.value_bf_audit();
                        bf_audit != nullptr) {
                        static_cast<void>(ClassGroupCertificationAccess::
                                record_grh_acceptance_bf_audit(
                                        class_group,
                                        flint::ArbConstRef(
                                                bf_audit->error_bound),
                                        bf_audit->cutoff,
                                        bf_audit->max_cutoff,
                                        bf_audit->requested_precision,
                                        bf_audit->work_precision));
                    }
                }
                if (exact_quadratic_index && requested_grh) {
                    // Index one against the exact class number: the
                    // relations and units are complete unconditionally,
                    // relative to generation.  The labels stay `grh`.
                    ClassGroupCertificationAccess::
                            record_grh_acceptance_exact_class_number(
                                    class_group);
                }
                summary.outcome = ValidateRefineOutcome::proven;
                return true;
            }
            // Source trace: reference `Clgp.jl:_class_unit_group` finishes the
            // tentative class/unit loop when `_validate_class_unit_group`
            // returns index one.  Native validation computed the same
            // analytic index-one result above; try to publish through the
            // class/unit certification helper before reporting publication
            // failure.
            if (ClassGroupCertificationAccess::
                        try_certify_class_unit_with_units(
                                class_group, units, analytic_cache.value(),
                                precision,
                                analytic_cache.value_unconditional())) {
                summary.outcome = ValidateRefineOutcome::proven;
                return true;
            }
            // A GRH-conditional index one cannot publish `proven`; prove
            // the pair unconditionally by unit proof and saturation.
            if (try_prove_candidate_pair_by_saturation(
                        class_group, units, order, options, embeddings,
                        precision)) {
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
                    exact_quadratic_index,
                    refresh_options.relation_saturation_memo)) {
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
            // the rebuild (a fresh relation subset at a larger unit_add) is
            // deferred to the next outer iteration, matching Clgp.jl
            // _class_unit_group.
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

        // The refined index only decides whether this pass made progress.
        // Under a grh request in degree two it comes from the exact class
        // number where the pass check above does (|D| < 2^34), so no
        // unconditional L(1, chi) hR is evaluated just to measure progress;
        // above the cutoff, or when that route is unavailable, it comes from
        // the analytic index bound as before.
        bool next_check_available = false;
        if (requested_grh && order.degree() == 2) {
            next_check_available = ClassGroupCertificationAccess::
                    rank_zero_quadratic_class_index_bound(
                            flint::FmpzRef(next_index_bound), class_group,
                            scratch_units, options.requested_certification);
        }
        if (!next_check_available) {
            next_check_available = analytic_index_bound_for_validation(
                    next_index_bound, scratch_units, class_group,
                    analytic_cache, order, precision, requested_grh);
        }
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
