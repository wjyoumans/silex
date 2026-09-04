#pragma once

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/order_unit.hpp>

#include "relation_unit_internal.hpp"
#include "../zeta/zeta_internal.hpp"

namespace silex {

namespace detail {

inline constexpr slong kDefaultComputePrecision = 160;
inline constexpr slong kComputeSatAuxTarget = 1;
inline constexpr slong kComputeSatAuxStart = 2;
inline constexpr slong kComputeSatAuxMax = 31;
inline constexpr slong kComputeProofAuxMax = 1000;
inline constexpr slong kComputeProofRetryAuxMax = 127;
inline constexpr slong kComputeSatMaxPasses = 2;
inline constexpr slong kComputeRevalidateMaxPasses = 1;
inline constexpr slong kRankOneAdaptiveSurplus = 5;
inline constexpr slong kUnitCandidateBatchTarget = 10;
inline constexpr slong kComputeAnalyticBumpMaxPasses = 4;
inline constexpr slong kComputeAnalyticBumpStepMax = 4;
inline constexpr slong kComputePostFinitePhaseBudget = 64;

inline slong max_slong(slong left, slong right) noexcept {
    return left > right ? left : right;
}

inline slong min_slong(slong left, slong right) noexcept {
    return left < right ? left : right;
}

bool rank_zero_torsion(flint::FmpzRef best_order,
                       OrderElement& best_generator,
                       const Order& order) noexcept;

enum class ValidateRefineOutcome {
    not_run,
    proven,
    not_proven_request,
    analytic_unavailable,
    analytic_inconclusive,
    index_one_publication_failed,
    class_saturation_failed,
    no_progress,
    unit_refinement_failed,
    local_rank_unavailable,
    pass_cap_reached
};

struct ValidateRefineSummary {
    ValidateRefineOutcome outcome = ValidateRefineOutcome::not_run;
    flint::Fmpz last_index_bound;
    bool class_progress = false;
    bool unit_progress = false;

    void reset() noexcept {
        outcome = ValidateRefineOutcome::not_run;
        flint::fmpz_zero(flint::FmpzRef(last_index_bound));
        class_progress = false;
        unit_progress = false;
    }
};

const char* validate_refine_outcome_name(ValidateRefineOutcome outcome)
        noexcept;
const char* validate_refine_outcome_profile_label(
        ValidateRefineOutcome outcome) noexcept;

bool analytic_class_regulator_product_for_validation(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

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
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

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
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

class AnalyticClassRegulatorCache {
public:
    void configure_validation(ulong max_cutoff) noexcept {
        if (validation_max_cutoff_ != max_cutoff) {
            validation_active_ = false;
        }
        validation_enabled_ = max_cutoff != 0;
        validation_max_cutoff_ = max_cutoff;
    }

    bool validation_enabled() const noexcept {
        return validation_enabled_;
    }

    bool validation_active() const noexcept {
        return validation_active_;
    }

    bool ensure(const Order& order,
                slong precision,
                const DiagnosticsContext* diagnostics = nullptr,
                const FactorBase* residue_degree_base = nullptr,
                bool cache_bf_residue_degrees = false) noexcept {
        reset_for_order_if_needed_(order);
        if (validation_enabled_ && order.degree() > 2 &&
            ensure_validation(
                    order, validation_max_cutoff_, precision,
                    diagnostics, residue_degree_base)) {
            return true;
        }
        validation_active_ = false;
        if (valid_ > 0 && precision_ >= precision) {
            return true;
        }
        if (valid_ < 0 && precision_ == precision) {
            return false;
        }

        valid_ = analytic_class_regulator_product_for_validation(
                         flint::ArbRef(value_), order, precision, diagnostics,
                         residue_degree_base,
                         cache_bf_residue_degrees ? &bf_residue_degree_cache_
                                                  : nullptr)
                ? 1
                : -1;
        precision_ = precision;
        return valid_ > 0;
    }

    void seed(const Order& order,
              flint::ArbConstRef value,
              slong precision) noexcept {
        seed_order_(order);
        flint::arb_set(flint::ArbRef(value_), value);
        precision_ = precision;
        valid_ = 1;
        validation_active_ = false;
    }

    flint::ArbConstRef value() const noexcept {
        return validation_active_
                ? flint::ArbConstRef(validation_value_)
                : flint::ArbConstRef(value_);
    }

    bool ensure_validation(
            const Order& order,
            ulong max_cutoff,
            slong precision,
            const DiagnosticsContext* diagnostics,
            const FactorBase* residue_degree_base = nullptr) noexcept {
        if (max_cutoff == 0 || precision <= 0 || order.degree() <= 2) {
            return false;
        }
        reset_for_order_if_needed_(order);
        if (validation_valid_ != 0 &&
            validation_computed_max_cutoff_ == max_cutoff &&
            validation_precision_ == precision) {
            validation_active_ = validation_valid_ > 0;
            return validation_active_;
        }

        validation_computed_max_cutoff_ = max_cutoff;
        validation_precision_ = precision;
        validation_valid_ =
                class_regulator_product_for_validation(
                        flint::ArbRef(validation_value_),
                        flint::ArbRef(validation_error_bound_),
                        validation_cutoff_,
                        validation_work_precision_, order, max_cutoff,
                        precision, diagnostics, residue_degree_base,
                        &bf_residue_degree_cache_)
                ? 1
                : -1;
        validation_active_ = validation_valid_ > 0;
        return validation_active_;
    }

    bool has_validation(ulong max_cutoff,
                              slong precision) const noexcept {
        return validation_valid_ > 0 &&
               validation_computed_max_cutoff_ == max_cutoff &&
               validation_precision_ == precision;
    }

    flint::ArbConstRef zeta_validation_value() const noexcept {
        return flint::ArbConstRef(validation_value_);
    }

    flint::ArbConstRef validation_error_bound() const noexcept {
        return flint::ArbConstRef(validation_error_bound_);
    }

    ulong validation_cutoff() const noexcept {
        return validation_cutoff_;
    }

    slong validation_work_precision() const noexcept {
        return validation_work_precision_;
    }

    bool ensure_bf_audit(const Order& order,
                         ulong max_cutoff,
                         slong precision,
                         const DiagnosticsContext* diagnostics,
                         const FactorBase* residue_degree_base = nullptr)
            noexcept {
        if (max_cutoff == 0 || precision <= 0) {
            return false;
        }
        reset_for_order_if_needed_(order);
        if (bf_valid_ != 0 && bf_max_cutoff_ == max_cutoff &&
            bf_precision_ == precision) {
            return bf_valid_ > 0;
        }

        bf_max_cutoff_ = max_cutoff;
        bf_precision_ = precision;
        bf_valid_ = bf_class_regulator_product_for_validation(
                            flint::ArbRef(bf_value_),
                            flint::ArbRef(bf_error_bound_), bf_cutoff_,
                            bf_work_precision_, order, max_cutoff, precision,
                            diagnostics, residue_degree_base,
                            &bf_residue_degree_cache_)
                ? 1
                : -1;
        return bf_valid_ > 0;
    }

    bool has_bf_audit(ulong max_cutoff, slong precision) const noexcept {
        return bf_valid_ > 0 && bf_max_cutoff_ == max_cutoff &&
               bf_precision_ == precision;
    }

    flint::ArbConstRef bf_value() const noexcept {
        return flint::ArbConstRef(bf_value_);
    }

    flint::ArbConstRef bf_error_bound() const noexcept {
        return flint::ArbConstRef(bf_error_bound_);
    }

    bool seed_bf_residue_degree_cache(
            const Order& order,
            const ZetaBfResidueDegreeCache& source,
            const DiagnosticsContext* diagnostics) noexcept {
        if (source.entries.empty() || source.residue_degrees.empty()) {
            return false;
        }
        ulong previous_prime = 0;
        for (const ZetaBfResidueDegreeCacheEntry& entry : source.entries) {
            if (entry.p <= previous_prime || entry.length == 0 ||
                entry.offset > source.residue_degrees.size() ||
                entry.length >
                        source.residue_degrees.size() - entry.offset) {
                return false;
            }
            previous_prime = entry.p;
            for (std::size_t i = 0; i < entry.length; ++i) {
                const slong degree = source.residue_degrees[entry.offset + i];
                if (degree <= 0 || degree > order.degree()) {
                    return false;
                }
            }
        }
        reset_for_order_if_needed_(order);
        bf_residue_degree_cache_ = source;
        bf_residue_degree_cache_.lookup_hint = 0;
        bf_valid_ = 0;
        validation_valid_ = 0;
        validation_active_ = false;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.validation_analytic_index_bound."
                "prime_cache_from_relation_factor_base_plan");
        return true;
    }

    ulong bf_cutoff() const noexcept {
        return bf_cutoff_;
    }

    slong bf_work_precision() const noexcept {
        return bf_work_precision_;
    }

private:
    flint::Arb value_;
    slong precision_ = 0;
    int valid_ = 0;
    flint::Arb bf_value_;
    flint::Arb bf_error_bound_;
    ulong bf_cutoff_ = 0;
    ulong bf_max_cutoff_ = 0;
    slong bf_precision_ = 0;
    slong bf_work_precision_ = 0;
    int bf_valid_ = 0;
    flint::Arb validation_value_;
    flint::Arb validation_error_bound_;
    ulong validation_cutoff_ = 0;
    ulong validation_max_cutoff_ = 0;
    ulong validation_computed_max_cutoff_ = 0;
    slong validation_precision_ = 0;
    slong validation_work_precision_ = 0;
    int validation_valid_ = 0;
    bool validation_enabled_ = false;
    bool validation_active_ = false;
    ZetaBfResidueDegreeCache bf_residue_degree_cache_;
    Order bf_residue_degree_cache_order_;
    bool bf_residue_degree_cache_has_order_ = false;

    void reset_for_order_if_needed_(const Order& order) noexcept {
        if (bf_residue_degree_cache_has_order_ &&
            bf_residue_degree_cache_order_.has_same_data(order)) {
            return;
        }
        seed_order_(order);
        valid_ = 0;
        bf_valid_ = 0;
        validation_valid_ = 0;
        validation_active_ = false;
    }

    void seed_order_(const Order& order) noexcept {
        if (bf_residue_degree_cache_has_order_ &&
            bf_residue_degree_cache_order_.has_same_data(order)) {
            return;
        }
        bf_residue_degree_cache_.entries.clear();
        bf_residue_degree_cache_.residue_degrees.clear();
        bf_residue_degree_cache_.lookup_hint = 0;
        bf_residue_degree_cache_order_ = order;
        bf_residue_degree_cache_has_order_ = true;
        bf_valid_ = 0;
        validation_valid_ = 0;
        validation_active_ = false;
    }
};

bool class_group_order_is_one(const ClassGroupContext& class_group) noexcept;
slong initial_relation_kernel_target(slong rank,
                                     slong degree,
                                     slong lower_target) noexcept;
bool missing_unit_target_bump(slong& bump,
                              slong rank,
                              slong unit_count,
                              slong target,
                              slong max_relations) noexcept;
bool class_unit_validation_estimate(
        flint::Fmpz& index_bound,
        flint::Arb& expected_regulator,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision) noexcept;
bool class_unit_bf_validation_estimate(
        flint::Fmpz& index_bound,
        flint::Arb& expected_regulator,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;
bool expected_regulator_stop(
        const OrderUnitGroup& units,
        const flint::Arb& expected_regulator,
        slong precision) noexcept;
bool validation_target_bump(slong& bump,
                            const ValidateRefineSummary& summary,
                            slong rank,
                            slong target,
                            slong max_relations,
                            bool requested_proven,
                            bool rank_complete) noexcept;

enum class ValidationRecomputeCause {
    initial_requested_proven,
    post_compute_initial_candidate,
    missing_candidate_units,
    validation_loop_unit_refinement,
    proven_refinement_retry,
};

struct ValidationUnitRefreshOptions {
    RelationUnitExtractionState* relation_unit_state = nullptr;
    slong* unit_add = nullptr;
    slong* unit_improved = nullptr;
    slong rank = 0;
    ulong validation_bf_max_cutoff = 0;
    bool use_class_relation_units = false;
    bool require_source_units = false;
};

bool set_class_relation_units(bool& ready,
                                    slong& improved,
                                    OrderUnitGroup& out,
                                    const Order& order,
                                    const ClassGroupContext& class_group,
                                    EmbeddingContext& embeddings,
                                    AnalyticClassRegulatorCache& analytic_cache,
                                    RelationUnitExtractionState& extraction_state,
                                    slong precision,
                                    slong rank,
                                    slong outer_add,
                                    ulong validation_bf_max_cutoff) noexcept;
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
                                                refresh_options = {}) noexcept;
bool complete_requested_proven_relation_saturation(
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        const Order& order,
        const ClassGroupComputeOptions& options,
        AnalyticClassRegulatorCache& analytic_cache,
        slong precision) noexcept;
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
                                      refresh_options = {}) noexcept;
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
                                            refresh_options = {}) noexcept;

}  // namespace detail

}  // namespace silex
