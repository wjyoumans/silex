#pragma once

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/order_unit.hpp>

#include <atomic>
#include <vector>

#include "relation_unit_internal.hpp"
#include "../zeta/zeta_internal.hpp"

namespace silex {

namespace detail {

inline constexpr slong kDefaultComputePrecision = 160;
inline constexpr slong kComputeSatAuxTarget = 1;
inline constexpr slong kComputeSatAuxStart = 2;
inline constexpr slong kComputeSatAuxMax = 31;
// Upper end of the bounded auxiliary-prime pre-scan in prove_index_bound.  It
// is not a proof limit: an ell < kComputeProofAuxMax whose pre-scan leaves a
// nonempty character kernel continues with the unbounded reference `saturate!`
// scan, and every ell >= kComputeProofAuxMax uses that scan directly.
inline constexpr slong kComputeProofAuxMax = 1000;
// Upper end of the auxiliary primes q of the ell-local class-relation
// saturation proof in the proven paired transaction
// (prove_relation_saturation_dlog_ell_).  It is neither a proof limit nor a
// fixed resource bound.  The source bnfcertify local check (2.17.3
// src/basemath/buch3.c:primecertify, lines 1086-1146) walks every prime
// q = 1 mod p from u_forprime_arith_init(&T, 1, ULONG_MAX, 1, p) until the
// character matrix reaches full rank; the proof does the same over the whole
// machine-word range of q, and checks every q it uses exactly as before.
// The pair it proves passed the Belabas-Friedman index-one test, so under
// GRH its relations are saturated and the scan ends at full rank, as
// primecertify does on a bnfinit result.  Silex has no in-library time
// limit: only the caller's time or memory limit on the request (for example
// a process timeout) stops a longer scan, and a stopped request publishes
// nothing, so no result is labelled proven.  A scan that reaches the end of
// the range records the ell-local proof as unavailable and the request
// fails closed.
inline constexpr ulong kComputeProofClassAuxMax = UWORD_MAX;
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

// Relations per refill of the grh class/unit loop
// (compute_with_relation_class_group_): the reference
// `class_group_new_relations_via_lll` slice returns after 2 once h > 0
// (reference source Clgp/Main_LLL.jl).
inline constexpr slong kRelationRefillInitial = 2;

// The refill target of the next grh class/unit loop slice.  `repeated` says
// the validation just run ended with the same index bound (> 1) as the one
// before the previous refill, and that refill left the class order
// unchanged: the previous slice gave no measurable progress.  The target
// then doubles (a Silex rule, not in the reference: k stalled rounds of 2
// relations become about log2 k larger refills, each followed by one
// validation).  Any other round restarts at the reference 2.  Saturates at WORD_MAX; the route's
// max_relations cap still bounds every slice.
inline slong next_relation_refill_target(slong current,
                                         bool repeated) noexcept {
    if (!repeated || current < kRelationRefillInitial) {
        return kRelationRefillInitial;
    }
    return current > WORD_MAX / 2 ? WORD_MAX : 2 * current;
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

class TestHookAccess;

// Switches that only detail::TestHookAccess (test/test_hook_access.hpp) can
// set; compute.cpp reads them on cold failure paths.  All are off by default
// and no production code sets them, so a normal run never sees them.  Relaxed
// atomics so a reader never races with a test setting them.
class ComputeTestSeams {
public:
    // compute_with_class_group: after a failed candidate attempt, forget the
    // validation summary, as if validation never ran.
    static bool forget_validation_summary() noexcept {
        return forget_validation_summary_.load(std::memory_order_relaxed);
    }
    // compute_with_relation_class_group_: true once per arming, when the
    // continuation should skip its extend_tentative_relations_ call and act
    // as if it returned true without appending a relation.  Reading consumes
    // the arming, so a later continuation runs normally.
    static bool take_continuation_appends_nothing() noexcept {
        return continuation_appends_nothing_.exchange(
                false, std::memory_order_relaxed);
    }
    // word_quotient_log_resolve (relation_kernel.cpp): counts the calls that
    // took each discrete-log route, so a test can observe the route choice
    // (shared walk below kWordQuotientLogBabyGiantBits bits of ell, baby-step
    // giant-step from there on).  One relaxed increment per call, which is
    // once per residue_quotient_log_column; it never changes a result.
    static void count_word_quotient_log_walk() noexcept {
        word_quotient_log_walk_calls_.fetch_add(
                1, std::memory_order_relaxed);
    }
    static void count_word_quotient_log_baby_giant() noexcept {
        word_quotient_log_baby_giant_calls_.fetch_add(
                1, std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<unsigned long> word_quotient_log_walk_calls_{0};
    inline static std::atomic<unsigned long>
            word_quotient_log_baby_giant_calls_{0};
    inline static std::atomic<bool> forget_validation_summary_{false};
    inline static std::atomic<bool> continuation_appends_nothing_{false};
};

const char* validate_refine_outcome_name(ValidateRefineOutcome outcome)
        noexcept;
const char* validate_refine_outcome_profile_label(
        ValidateRefineOutcome outcome) noexcept;

// `unconditional` reports whether the product came from an unconditional
// route, and `bf_audit` receives the Belabas-Friedman audit data when it did
// not; `accept_grh_residue` takes the Belabas-Friedman residue for a grh
// request (see zeta_class_regulator_product_with_diagnostics).
bool analytic_class_regulator_product_for_validation(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        bool* unconditional = nullptr,
        ZetaBfRouteAudit* bf_audit = nullptr,
        bool accept_grh_residue = false) noexcept;

// The grh index-one hR (zeta_class_regulator_product_index_one_with_
// diagnostics): a Belabas-Friedman value with log radius below 1/4, enough
// for the analytic index-one test only.
bool index_one_class_regulator_product_for_validation(
        flint::ArbRef out,
        ZetaBfRouteAudit& audit,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept;

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

    // `accept_grh_residue` is set for a grh request, which accepts a
    // GRH-conditional hR: the default route then takes the Belabas-Friedman
    // residue instead of the quadratic L(1, chi) sum
    // (zeta_class_regulator_product_with_diagnostics).  A value computed that
    // way is not reused for a call without it, which keeps the unconditional
    // route for proven requests; an unconditional value serves both.
    bool ensure(const Order& order,
                slong precision,
                const DiagnosticsContext* diagnostics = nullptr,
                const FactorBase* residue_degree_base = nullptr,
                bool cache_bf_residue_degrees = false,
                bool accept_grh_residue = false) noexcept {
        reset_for_order_if_needed_(order);
        index_one_active_ = false;
        if (validation_enabled_ && order.degree() > 2 &&
            ensure_validation(
                    order, validation_max_cutoff_, precision,
                    diagnostics, residue_degree_base)) {
            return true;
        }
        validation_active_ = false;
        // The flag changes the route only where the default route has the
        // quadratic L(1, chi) sum; elsewhere the value is the same either
        // way, so it is recorded without the flag and serves every call.
        const bool route_accepts_grh_residue =
                accept_grh_residue &&
                zeta_unconditional_route_available(order);
        const bool value_serves_request =
                accept_grh_residue || value_unconditional_ ||
                !value_accept_grh_residue_;
        if (valid_ > 0 && precision_ >= precision && value_serves_request) {
            return true;
        }
        if (valid_ < 0 && precision_ == precision &&
            value_accept_grh_residue_ == route_accepts_grh_residue) {
            return false;
        }

        value_unconditional_ = false;
        value_bf_audit_valid_ = false;
        value_accept_grh_residue_ = route_accepts_grh_residue;
        // A grh request shares this order's cache with ensure_index_one,
        // which fills it over the same primes, so the default-route
        // Belabas-Friedman sum reuses its residue degrees and linear-factor
        // counts instead of recomputing them; both are exact per prime, so
        // the value is the same.
        const bool use_residue_degree_cache =
                cache_bf_residue_degrees || accept_grh_residue;
        valid_ = analytic_class_regulator_product_for_validation(
                         flint::ArbRef(value_), order, precision, diagnostics,
                         residue_degree_base,
                         use_residue_degree_cache ? &bf_residue_degree_cache_
                                                  : nullptr,
                         &value_unconditional_, &value_bf_audit_,
                         accept_grh_residue)
                ? 1
                : -1;
        value_bf_audit_valid_ = valid_ > 0 && !value_unconditional_ &&
                                order.degree() > 1;
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
        // A seeded value (a Belabas-Friedman audit or a route estimate) is
        // never an unconditional hR, and carries no default-route audit.
        // It serves every later request, as it is not a default-route
        // value that a call without accept_grh_residue would recompute.
        value_unconditional_ = false;
        value_bf_audit_valid_ = false;
        value_accept_grh_residue_ = false;
        validation_active_ = false;
        index_one_active_ = false;
    }

    flint::ArbConstRef value() const noexcept {
        if (validation_active_) {
            return flint::ArbConstRef(validation_value_);
        }
        return index_one_active_ ? flint::ArbConstRef(index_one_value_)
                                 : flint::ArbConstRef(value_);
    }

    // True only when value() is an unconditional hR: degree one or the
    // quadratic L(1, chi) route.  Belabas-Friedman values, including the
    // validation enclosure, assume GRH and report false.
    bool value_unconditional() const noexcept {
        return !validation_active_ && !index_one_active_ && valid_ > 0 &&
               value_unconditional_;
    }

    // The Belabas-Friedman audit data of value() when it came from the
    // default zeta route's Belabas-Friedman fallback or from the grh
    // index-one value (not the validation enclosure, not a seeded value);
    // nullptr otherwise.
    const ZetaBfRouteAudit* value_bf_audit() const noexcept {
        if (validation_active_) {
            return nullptr;
        }
        if (index_one_active_) {
            return &index_one_audit_;
        }
        return valid_ > 0 && value_bf_audit_valid_ ? &value_bf_audit_
                                                   : nullptr;
    }

    // Makes value() the grh index-one hR, a Belabas-Friedman value whose
    // log radius is below 1/4 (index_one_class_regulator_product_for_
    // validation), and returns true; the next ensure() or seed() switches
    // back to the default value.  Returns false, leaving value() as it is,
    // when the validation enclosure is enabled, when the default value is
    // already cached at `precision` or more (it is at least as narrow), or
    // when the index-one value is unavailable (degree one, or no cutoff up
    // to the default maximum reaches 1/4).  The caller must
    // use the default value from ensure() for anything other than
    // recognizing index one.
    bool ensure_index_one(const Order& order,
                          slong precision,
                          const DiagnosticsContext* diagnostics = nullptr,
                          const FactorBase* residue_degree_base =
                                  nullptr) noexcept {
        reset_for_order_if_needed_(order);
        index_one_active_ = false;
        if (validation_enabled_ || precision <= 0 ||
            (valid_ > 0 && precision_ >= precision)) {
            return false;
        }
        if (index_one_valid_ == 0 || index_one_precision_ != precision) {
            index_one_precision_ = precision;
            index_one_valid_ =
                    index_one_class_regulator_product_for_validation(
                            flint::ArbRef(index_one_value_),
                            index_one_audit_, order, precision, diagnostics,
                            residue_degree_base, &bf_residue_degree_cache_)
                    ? 1
                    : -1;
        }
        if (index_one_valid_ <= 0) {
            return false;
        }
        validation_active_ = false;
        index_one_active_ = true;
        return true;
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
        index_one_active_ = false;
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
        previous_prime = 0;
        for (const ZetaBfLinearFactorCountEntry& entry :
             source.linear_factor_counts) {
            if (entry.p <= previous_prime || entry.count < 0 ||
                entry.count > order.degree()) {
                return false;
            }
            previous_prime = entry.p;
        }
        reset_for_order_if_needed_(order);
        bf_residue_degree_cache_ = source;
        bf_residue_degree_cache_.lookup_hint = 0;
        bf_residue_degree_cache_.linear_factor_count_hint = 0;
        bf_valid_ = 0;
        validation_valid_ = 0;
        validation_active_ = false;
        index_one_valid_ = 0;
        index_one_active_ = false;
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
    bool value_unconditional_ = false;
    // Whether value_ came from a call with accept_grh_residue (see ensure).
    bool value_accept_grh_residue_ = false;
    ZetaBfRouteAudit value_bf_audit_;
    bool value_bf_audit_valid_ = false;
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
    flint::Arb index_one_value_;
    ZetaBfRouteAudit index_one_audit_;
    slong index_one_precision_ = 0;
    int index_one_valid_ = 0;
    bool index_one_active_ = false;
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
        index_one_valid_ = 0;
        index_one_active_ = false;
    }

    void seed_order_(const Order& order) noexcept {
        if (bf_residue_degree_cache_has_order_ &&
            bf_residue_degree_cache_order_.has_same_data(order)) {
            return;
        }
        bf_residue_degree_cache_.entries.clear();
        bf_residue_degree_cache_.residue_degrees.clear();
        bf_residue_degree_cache_.lookup_hint = 0;
        bf_residue_degree_cache_.linear_factor_counts.clear();
        bf_residue_degree_cache_.linear_factor_count_hint = 0;
        bf_residue_degree_cache_order_ = order;
        bf_residue_degree_cache_has_order_ = true;
        bf_valid_ = 0;
        validation_valid_ = 0;
        validation_active_ = false;
        index_one_valid_ = 0;
        index_one_active_ = false;
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
        slong precision,
        bool accept_grh_generation) noexcept;
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
        bool accept_grh_generation) noexcept;
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

// The last bounded class-relation saturation pass of grh validation that
// appended no relation and cleared no prime of h_cand, kept across the
// iterations of the class/unit loop (compute_with_relation_class_group_).
// While h_cand is exactly the recorded value, the factor base has the
// recorded number of generators and the candidate regulator is within a
// certified factor two of the recorded one, the pass is not repeated:
// saturate_candidate_class_relations_with_units replays its state instead
// (ClassGroupCertificationAccess::
// replay_unchanged_bounded_relation_saturation).  The skip is a heuristic
// over the same lattices and the same bounds (ell set, auxiliary-prime
// bound, append limit): it assumes that the pass would again append
// nothing, and that is empirical, not proven.  The pass power-tests kernel
// rows over the actual relation rows and unit generators, so appended rows
// that leave the lattices unchanged can change what it tests; and the
// factor-two regulator test proves the unit lattices equal only when one
// contains the other, which the units rebuilt in each iteration do not
// guarantee.  Identical outputs were observed on every grh row of the
// class/unit benchmark corpus.  This is a cost memo only: grh acceptance
// never rests on this saturation, which only adds relations when it finds
// them.
struct RelationSaturationMemo {
    bool valid = false;
    flint::Fmpz class_order;
    flint::Arb regulator;
    flint::Fmpz aux_bound;
    slong generator_count = 0;
    slong max_appends_per_ell = 0;
    slong attempt_limit = 0;
    // Number of passes skipped through the memo, for tests.
    slong skips = 0;
};

struct ValidationUnitRefreshOptions {
    RelationUnitExtractionState* relation_unit_state = nullptr;
    slong* unit_add = nullptr;
    slong* unit_improved = nullptr;
    slong rank = 0;
    ulong validation_bf_max_cutoff = 0;
    bool use_class_relation_units = false;
    bool require_source_units = false;
    // True exactly for a grh request: lets the validation estimates count
    // the GRH factor-base generation record.
    bool accept_grh_generation = false;
    // Set only by the grh class/unit loop; null runs every saturation pass.
    RelationSaturationMemo* relation_saturation_memo = nullptr;
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
                                    ulong validation_bf_max_cutoff,
                                    bool accept_grh_generation) noexcept;
// Replacement for the a priori denominator bound of add_dependent_unit
// (OrderUnitGroup::regulator_index_bound) when that bound does not fit in
// slong, the width the rational-reconstruction search takes.  Only a grh
// request supplies one.  The value is the class/unit index bound
// hR_cand / hR from class_unit_index_bound_for_improvement, which bounds
// the unit index [U:G] under GRH once the factor base generates the class
// group (see "Class groups and order units" in
// docs/reference/algorithms_and_sources.rst).  With the context pointers
// set, add_dependent_unit computes it for the current group on the first
// overflow and again after each adjoined unit; without them it uses
// `bound` as given when `usable` is set.  Every relation it finds is still
// checked exactly, so a bound that is too small only misses relations.
struct DependentUnitIndexFallback {
    const ClassGroupContext* class_group = nullptr;
    AnalyticClassRegulatorCache* analytic_cache = nullptr;
    const Order* order = nullptr;
    slong precision = 0;
    flint::Fmpz bound;
    bool computed = false;
    bool usable = false;
};

// add_dependent_unit, with `fallback` (nullable) replacing an a priori
// denominator bound that does not fit in slong.
bool add_dependent_unit(bool& changed,
                        OrderUnitGroup& group,
                        const FactoredElement& candidate,
                        EmbeddingContext& embeddings,
                        RelationUnitExtractionState& extraction_state,
                        slong precision,
                        DependentUnitIndexFallback* fallback) noexcept;

// The partial unit basis of the relation-unit search before it reaches full
// rank: the reference `add_unit!` over `_isindependent` (see "Class groups
// and order units" in docs/reference/algorithms_and_sources.rst).  A
// candidate is kept when compact_independent's rule shows it independent of
// the units already kept.  add_independent_unit first tries that rule on an
// incremental enclosure of the Gram determinant, with compact_independent's
// precision doubling and cap, and hands the candidate to
// compact_independent, unchanged, whenever the incremental test does not
// decide; see the soundness note in compute.cpp.  Wherever
// compact_independent decides, the decision is the same.  A candidate it
// dropped at the cap, or after a failed logarithmic embedding, may now be
// decided at higher precision (both need about 2^20 bits).
struct IncrementalUnitContext {
    const NumberField* field = nullptr;
    slong target_rank = 0;
    std::vector<FactoredElement> units;
    bool full_rank = false;

    // Cost memo of the incremental test at one working precision, never a
    // decision by itself.  rows[i] encloses the product-mode log embedding
    // of units[i].  When gram_units equals the number of kept units,
    // gram_inverse and gram_determinant enclose the inverse and the
    // determinant of the Gram matrix of those rows.
    struct LogState {
        slong precision = 0;
        std::vector<flint::ArbVec> rows;
        flint::ArbMat gram_inverse{0, 0};
        flint::Arb gram_determinant;
        slong gram_units = -1;
    };
    // One state per working precision the doubling loop has reached.
    std::vector<LogState> log_states;

    // For tests: candidates the incremental test decided, and candidates
    // passed to compact_independent.
    slong incremental_decisions = 0;
    slong fallback_decisions = 0;

    bool define(const NumberField& parent, slong rank) noexcept;
    // Appends `candidate` to `units` when it is independent of them and
    // sets `changed`.  An undecided or failed independence test leaves
    // `units` unchanged and returns true, as before.  Returns false on an
    // undefined context, a full-rank context, a non-positive precision, a
    // candidate from another field, or when copying the candidate fails.
    bool add_independent_unit(bool& changed,
                              const FactoredElement& candidate,
                              EmbeddingContext& embeddings,
                              slong precision,
                              const DiagnosticsContext* diagnostics) noexcept;
};

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
// The relation-saturation pass of class/unit validation against
// `index_bound`, a bound on h_cand / h.  With `index_bound_is_exact` it is
// the exact imaginary-quadratic index, which counts as unconditional only
// when factor-base generation is verified up to the required bound.
// A non-null `memo` lets the bounded pass (inexact `index_bound`) be
// skipped and replayed as described at RelationSaturationMemo.
bool saturate_candidate_class_relations_with_units(
        bool& changed,
        ClassGroupContext& class_group,
        const OrderUnitGroup& units,
        const ClassGroupComputeOptions& options,
        flint::FmpzConstRef index_bound,
        bool index_bound_is_exact,
        RelationSaturationMemo* memo = nullptr) noexcept;
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

// Defined in compute.cpp; declared here for tests that pin their behaviour.
bool round_arb_approx_fmpz(flint::FmpzRef out, const arb_struct* value,
                           slong precision) noexcept;
bool apply_unit_reduction(std::vector<FactoredElement>& values,
                          const OrderUnitGroup& group,
                          const flint::FmpzMat& coefficients) noexcept;
bool unknown_saturation_at_two(
        bool& validated,
        bool& presentation_rebuilt,
        RelationUnitExtractionState& extraction_state,
        ValidateRefineSummary& summary,
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong expected_rank,
        slong precision,
        const DiagnosticsContext* diagnostics,
        bool accept_grh_generation) noexcept;
bool unknown_small_index_saturation(
        bool& validated,
        bool& attempted,
        RelationUnitExtractionState& extraction_state,
        ValidateRefineSummary& summary,
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong expected_rank,
        slong precision,
        const DiagnosticsContext* diagnostics,
        bool accept_grh_generation,
        bool entered_after_rebuilt_presentation) noexcept;

}  // namespace detail

}  // namespace silex
