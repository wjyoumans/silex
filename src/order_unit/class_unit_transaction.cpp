#include "class_unit_transaction_internal.hpp"

#include "../class_group/relation_search_driver_internal.hpp"
#include "../class_group/relation_factor_base_plan_internal.hpp"

#include <silex/signature.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>
#include <vector>

namespace silex::detail {
namespace {

std::size_t stage_index(ClassUnitStage stage) noexcept {
    return static_cast<std::size_t>(stage);
}

double elapsed_milliseconds(
        std::chrono::steady_clock::time_point start) noexcept {
    const auto elapsed = std::chrono::steady_clock::now() - start;
    return std::chrono::duration<double, std::milli>(elapsed).count();
}

bool select_policy(
        ClassUnitExecutionPolicy& policy,
        flint::Fmpz& selected_bound,
        RelationFactorBasePlan& relation_factor_base_plan,
        const Order& order,
        flint::FmpzConstRef requested_bound,
        const ClassGroupComputeOptions& options) noexcept {
    policy = ClassUnitExecutionPolicy{};
    relation_factor_base_plan = RelationFactorBasePlan{};
    if ((options.requested_certification != CertificationMode::grh &&
         options.requested_certification != CertificationMode::proven) ||
        order.parent() == nullptr || !order.is_maximal() ||
        order.degree() <= 0) {
        return false;
    }

    Signature signature_value;
    flint::Fmpz discriminant;
    if (!signature(signature_value, *order.parent()) ||
        signature_value.degree() != order.degree() ||
        !order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }

    policy.degree = order.degree();
    policy.unit_rank = signature_value.r1() + signature_value.r2() - 1;
    policy.discriminant_bits =
            static_cast<slong>(fmpz_bits(discriminant.raw()));
    // Degree one has an exact native class/unit route.  In particular, do not
    // send a GRH request through the positive-discriminant Bach-bound/LLL
    // policy: its discriminant-one bound is not a relation-search input.  The
    // exact result is relabelled conservatively at the transaction boundary.
    if (options.requested_certification == CertificationMode::grh &&
        policy.degree != 1) {
        if (!grh_factor_base_bound(
                    selected_bound, order, options.diagnostics,
                    &policy.grh_generation_bound,
                    &policy.grh_generation_basis)) {
            return false;
        }
        policy.factor_base_bound =
                NativeFactorBaseBoundStrategy::grh;
        policy.factor_base = NativeFactorBaseStrategy::norm_bounded;
        policy.relations = NativeRelationStrategy::lll;
        policy.units = policy.unit_rank == 0
                ? NativeUnitStrategy::rank_zero
                : NativeUnitStrategy::class_relation;
        policy.continuation = NativeContinuationStrategy::lll;
        policy.validation = NativeValidationStrategy::shared_exact;
        flint::fmpz_set(
                flint::FmpzRef(policy.selected_factor_base_bound),
                flint::FmpzConstRef(selected_bound));
        policy.selected = true;
        return true;
    }

    policy.factor_base_bound = NativeFactorBaseBoundStrategy::requested;
    policy.factor_base = NativeFactorBaseStrategy::complete;
    policy.relations = NativeRelationStrategy::native_experimental;
    const bool quadratic_relation_completion_table =
            order.degree() == 2;
    if (policy.unit_rank == 0) {
        policy.units = NativeUnitStrategy::rank_zero;
        policy.continuation = NativeContinuationStrategy::native_target_growth;
    } else if (order.degree() == 2 && policy.unit_rank == 1) {
        policy.relations = NativeRelationStrategy::relation_completion_table;
        policy.units = NativeUnitStrategy::quadratic_rank_one;
        policy.continuation = NativeContinuationStrategy::native_target_growth;
    } else {
        policy.units = NativeUnitStrategy::class_relation;
        policy.continuation = NativeContinuationStrategy::lll;
    }
    policy.validation = NativeValidationStrategy::shared_exact;

    flint::fmpz_set(flint::FmpzRef(selected_bound), requested_bound);
    if (quadratic_relation_completion_table || order.degree() >= 3) {
        const bool specialized_imaginary_quadratic =
                order.degree() == 2 && policy.unit_rank == 0 &&
                order.parent()->backend_kind() ==
                        NumberFieldBackendKind::quadratic;
        const bool plan_built = specialized_imaginary_quadratic
                ? relation_search::
                          build_maximal_imaginary_quadratic_factor_base_plan(
                                  relation_factor_base_plan, order,
                                  options.diagnostics)
                : relation_search::build_relation_factor_base_plan(
                          relation_factor_base_plan, order,
                          options.diagnostics);
        if (!plan_built) {
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(selected_bound),
                        flint::FmpzConstRef(
                                relation_factor_base_plan.working_bound));
        policy.factor_base_bound =
                NativeFactorBaseBoundStrategy::relation_completion;
        policy.factor_base =
                NativeFactorBaseStrategy::relation_completion;
        policy.relations = order.degree() == 2
                ? NativeRelationStrategy::relation_completion_table
                : NativeRelationStrategy::relation_completion_prepass;
        if (order.degree() >= 3 &&
            native_hnf_unit_strategy_for_signature(
                    policy.degree, policy.unit_rank)) {
            policy.units = NativeUnitStrategy::hnf;
            policy.continuation =
                    NativeContinuationStrategy::regulator_multiple_reconstruction;
        }
    }
    flint::fmpz_set(flint::FmpzRef(policy.selected_factor_base_bound),
                    flint::FmpzConstRef(selected_bound));
    policy.selected = true;
    return true;
}

// LIMCMAX = (long)(4 log^2 |D|) of reference buch2.c:Buchall_param (2.17.3,
// line 3772), in doubles as there: the cap of the reference factor-base
// restarts.  It is a sizing heuristic, not a generation bound.
bool grh_factor_base_restart_cap(ulong& out,
                                 const Order& order) noexcept {
    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(discriminant),
                    flint::FmpzConstRef(discriminant));
    const double log_discriminant = fmpz_dlog(discriminant.raw());
    const double cap = 4.0 * log_discriminant * log_discriminant;
    if (!std::isfinite(cap) || cap < 1.0 ||
        cap >= static_cast<double>(WORD_MAX)) {
        return false;
    }
    out = static_cast<ulong>(cap);
    return true;
}

// The next factor-base bound after a failed grh candidate computation, as
// the reference (2.17.3): a restart only while LIMC < LIMCMAX / 2
// (buch2.c:Buchall_param, line 3953), then buch1.c:bnf_increase_LIMC (lines
// 185-192): double while LIMC <= LIMCMAX / 13.333, otherwise add
// max(1, LIMCMAX / 20), and clamp to LIMCMAX.
bool next_grh_factor_base_bound(flint::Fmpz& bound, ulong cap) noexcept {
    if (flint::fmpz_sgn(flint::FmpzConstRef(bound)) <= 0 ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(bound), cap / 2) >= 0) {
        return false;
    }
    ulong limc = flint::fmpz_get_ui(flint::FmpzConstRef(bound));
    if (static_cast<double>(limc) <= static_cast<double>(cap) / 13.333) {
        limc *= 2;
    } else {
        limc += FLINT_MAX(UWORD(1), cap / 20);
    }
    if (limc > cap) {
        limc = cap;
    }
    flint::fmpz_set_ui(flint::FmpzRef(bound), limc);
    return true;
}

}  // namespace

bool grh_factor_base_bound(
        flint::Fmpz& out,
        const Order& order,
        const DiagnosticsContext* diagnostics,
        flint::Fmpz* generation_bound,
        ClassGroupGenerationBasis* generation_basis) noexcept {
    bool bach_selected = false;
    if (!grh_factor_base_bound_with_diagnostics(
                flint::FmpzRef(out), order, diagnostics, &bach_selected)) {
        return false;
    }
    if (generation_basis != nullptr) {
        *generation_basis = bach_selected ? ClassGroupGenerationBasis::bach
                                          : ClassGroupGenerationBasis::bdf;
    }
    if (generation_bound != nullptr) {
        flint::fmpz_set(flint::FmpzRef(*generation_bound),
                        flint::FmpzConstRef(out));
    }
    // Reference buch2.c:Buchall_param (2.17.3, lines 3826-3827)
    // takes LIMC = max(LIMC2, nthideal(N)): the base contains at least
    // N = [K : Q] prime ideals.  nthideal (lines 604-641) is the norm of the
    // N-th smallest prime ideal over the primes that are not inert (an inert
    // p, of residue degree N, is principal).  Probe bounds double until the
    // norm-bounded base holds N such ideals; all ideals of norm at most the
    // probe bound are then in it, so the N-th smallest norm is exact.
    const slong degree = order.degree();
    if (degree <= 0) {
        return false;
    }
    flint::Fmpz probe_bound;
    flint::fmpz_set(flint::FmpzRef(probe_bound), flint::FmpzConstRef(out));
    std::vector<flint::Fmpz> norms;
    for (;;) {
        FactorBase probe(order);
        if (!probe.is_defined() || !probe.build_prime_ideal_norm_bounded(
                                           flint::FmpzConstRef(probe_bound))) {
            return false;
        }
        norms.clear();
        for (slong i = 0; i < probe.length(); ++i) {
            const PrimeIdeal* prime = probe.prime_at(i);
            if (prime == nullptr || prime->residue_degree() <= 0) {
                return false;
            }
            if (prime->residue_degree() >= degree) {
                continue;
            }
            flint::Fmpz norm;
            if (!prime->norm(flint::FmpzRef(norm))) {
                return false;
            }
            norms.push_back(std::move(norm));
        }
        if (norms.size() >= static_cast<std::size_t>(degree)) {
            break;
        }
        flint::fmpz_mul_2exp(flint::FmpzRef(probe_bound),
                             flint::FmpzConstRef(probe_bound), UWORD(1));
    }
    const auto nth = norms.begin() + (degree - 1);
    std::nth_element(norms.begin(), nth, norms.end(),
                     [](const flint::Fmpz& left, const flint::Fmpz& right) {
                         return flint::fmpz_cmp(flint::FmpzConstRef(left),
                                                flint::FmpzConstRef(right)) <
                                0;
                     });
    if (flint::fmpz_cmp(flint::FmpzConstRef(*nth),
                        flint::FmpzConstRef(out)) > 0) {
        flint::fmpz_set(flint::FmpzRef(out), flint::FmpzConstRef(*nth));
    }
    return true;
}

void ClassUnitTransactionAccess::set_run_context(
        ClassGroupContext& context,
        ClassUnitTransactionContext* run_context) noexcept {
    context.class_unit_transaction_context_ = run_context;
}

ZetaBfResidueDegreeCache*
ClassUnitTransactionAccess::relation_factor_base_plan_residue_degrees(
        ClassGroupContext& context,
        const Order& order) noexcept {
    ClassUnitTransactionContext* const run_context =
            context.class_unit_transaction_context_;
    return run_context != nullptr &&
                   run_context->relation_factor_base_plan != nullptr &&
                   same_order_parent(context.parent(), &order)
            ? silex::detail::relation_factor_base_plan_residue_degrees(
                      *run_context->relation_factor_base_plan, order)
            : nullptr;
}

bool ClassGroupRelationSearchAccess::defer_native_goal_publication(
        const ClassGroupContext& context) noexcept {
    const ClassUnitTransactionContext* const run_context =
            context.class_unit_transaction_context_;
    return run_context != nullptr &&
           !uses_class_unit_kernel(run_context);
}

void ClassGroupRelationSearchAccess::record_factor_base_honesty_failure(
        ClassGroupContext& context,
        FactorBaseHonestyFailure&& failure) noexcept {
    ClassUnitTransactionContext* const run_context =
            context.class_unit_transaction_context_;
    if (run_context != nullptr) {
        run_context->audit.factor_base_honesty_failure = std::move(failure);
    }
}

void ClassGroupRelationSearchAccess::clear_factor_base_honesty_failure(
        ClassGroupContext& context) noexcept {
    ClassUnitTransactionContext* const run_context =
            context.class_unit_transaction_context_;
    if (run_context != nullptr) {
        run_context->audit.factor_base_honesty_failure =
                FactorBaseHonestyFailure{};
    }
}

const char* class_unit_computation_failure_reason(
        const ClassUnitTransactionReport& audit) noexcept {
    return audit.factor_base_honesty_failure.recorded
            ? "factor_base_honesty_unwitnessed"
            : "class_unit_computation_failed";
}

void ClassUnitTransactionReport::reset() noexcept {
    failure_stage = ClassUnitStage::none;
    failure_reason = nullptr;
    timing_ms.fill(0.0);
    timing_recorded.fill(false);
    policy = ClassUnitExecutionPolicy{};
    class_group_certification = CertificationMode::unknown;
    unit_group_certification = CertificationMode::unknown;
    unit_proof_status = ProofState::not_checked;
    regulator_proof_status = ProofState::not_checked;
    final_result_published = false;
    factor_base_honesty_failure = FactorBaseHonestyFailure{};
}

void ClassUnitTransactionReport::record_timing(ClassUnitStage stage,
                                      double elapsed_ms) noexcept {
    const std::size_t index = stage_index(stage);
    if (index >= timing_ms.size()) {
        return;
    }
    timing_ms[index] += elapsed_ms;
    timing_recorded[index] = true;
}

const char* class_unit_stage_name(ClassUnitStage stage) noexcept {
    switch (stage) {
    case ClassUnitStage::none:
        return "none";
    case ClassUnitStage::factor_base_bound:
        return "factor_base_bound";
    case ClassUnitStage::factor_base_construction:
        return "factor_base_construction";
    case ClassUnitStage::relation_scheduling:
        return "relation_scheduling";
    case ClassUnitStage::relation_collection:
        return "relation_collection";
    case ClassUnitStage::matrix_hnf_publication:
        return "matrix_hnf_publication";
    case ClassUnitStage::class_group_publication:
        return "class_group_publication";
    case ClassUnitStage::unit_extraction:
        return "unit_extraction";
    case ClassUnitStage::regulator_computation:
        return "regulator_computation";
    case ClassUnitStage::saturation:
        return "saturation";
    case ClassUnitStage::analytic_bf_zeta_validation:
        return "analytic_bf_zeta_validation";
    case ClassUnitStage::relation_extension:
        return "relation_extension";
    case ClassUnitStage::certification:
        return "certification";
    case ClassUnitStage::final_result_publication:
        return "final_result_publication";
    case ClassUnitStage::total:
        return "total";
    case ClassUnitStage::count:
        break;
    }
    return "none";
}

bool native_hnf_unit_strategy_for_signature(
        slong degree,
        slong unit_rank) noexcept {
    // The measured source-backed crossover favors reference HNF for cubics,
    // rank-one quartics, totally real fields, and every installed
    // degree-six-or-larger family.
    return unit_rank > 0 &&
            (degree == 3 || (degree == 4 && unit_rank == 1) || degree >= 6 ||
             unit_rank == degree - 1);
}

bool compute_class_unit_transaction(
        OrderUnitGroup& units,
        ClassGroupContext& class_group,
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupComputeOptions& options,
        slong precision,
        ClassUnitTransactionReport& audit) noexcept {
    audit.reset();
    const auto total_start = std::chrono::steady_clock::now();
    if (options.requested_certification != CertificationMode::grh &&
        options.requested_certification != CertificationMode::proven) {
        audit.failure_stage = ClassUnitStage::certification;
        audit.failure_reason = "unsupported_certification_request";
        audit.record_timing(ClassUnitStage::total,
                            elapsed_milliseconds(total_start));
        return false;
    }

    ClassGroupComputeOptions transaction_options = options;
    flint::Fmpz transaction_factor_base_bound;
    RelationFactorBasePlan relation_factor_base_plan;
    const auto factor_base_bound_start =
            std::chrono::steady_clock::now();
    if (!select_policy(
                audit.policy, transaction_factor_base_bound,
                relation_factor_base_plan, order, factor_base_bound,
                transaction_options)) {
        audit.failure_stage = ClassUnitStage::factor_base_bound;
        audit.failure_reason = "execution_policy_selection_failed";
        audit.record_timing(
                ClassUnitStage::factor_base_bound,
                elapsed_milliseconds(factor_base_bound_start));
        audit.record_timing(ClassUnitStage::total,
                            elapsed_milliseconds(total_start));
        return false;
    }
    audit.record_timing(ClassUnitStage::factor_base_bound,
                        elapsed_milliseconds(factor_base_bound_start));

    ClassGroupContext candidate_class_group;
    OrderUnitGroup candidate_units;
    ulong restart_cap = 0;
    const bool restartable = audit.policy.factor_base_bound ==
                    NativeFactorBaseBoundStrategy::grh &&
            grh_factor_base_restart_cap(restart_cap, order);
    bool success = false;
    for (;;) {
        ClassUnitTransactionContext transaction_context{
                audit,
                options.requested_certification == CertificationMode::proven
                        ? &relation_factor_base_plan
                        : nullptr,
                {},
                {},
                false};
        ClassUnitTransactionAccess::set_run_context(
                candidate_class_group, &transaction_context);
        success = candidate_units.compute_with_class_group(
                candidate_class_group, order,
                flint::FmpzConstRef(transaction_factor_base_bound),
                transaction_options, precision);
        ClassUnitTransactionAccess::set_run_context(candidate_class_group,
                                                    nullptr);
        // grh policy: when the candidate computation fails, enlarge the
        // factor base on the reference schedule and start again from
        // scratch (reference buch2.c:Buchall_param, 2.17.3, `goto START`
        // when the relations do not suffice, line 3953, which restarts only
        // while LIMC < LIMCMAX / 2; the new LIMC is buch1.c:
        // bnf_increase_LIMC, lines 185-192).  The GRH generation bound is
        // unchanged: the larger base still contains every prime ideal up
        // to it.
        if (success || !restartable ||
            !next_grh_factor_base_bound(
                    transaction_factor_base_bound, restart_cap)) {
            break;
        }
        flint::fmpz_set(
                flint::FmpzRef(audit.policy.selected_factor_base_bound),
                flint::FmpzConstRef(transaction_factor_base_bound));
        ++audit.policy.grh_factor_base_restarts;
        candidate_class_group = ClassGroupContext();
        candidate_units = OrderUnitGroup();
    }
    if (!success) {
        audit.failure_stage = ClassUnitStage::total;
        audit.failure_reason = class_unit_computation_failure_reason(audit);
        audit.record_timing(ClassUnitStage::total,
                            elapsed_milliseconds(total_start));
        return false;
    }

    if (candidate_class_group.certification_status() !=
                options.requested_certification ||
        candidate_units.certification_status() !=
                options.requested_certification) {
        audit.failure_stage = ClassUnitStage::certification;
        audit.failure_reason = "certification_incomplete";
        audit.record_timing(ClassUnitStage::total,
                            elapsed_milliseconds(total_start));
        return false;
    }

    class_group.swap(candidate_class_group);
    units.swap(candidate_units);
    audit.class_group_certification = class_group.certification_status();
    audit.unit_group_certification = units.certification_status();
    audit.unit_proof_status = class_group.unit_proof_status();
    audit.regulator_proof_status = class_group.regulator_proof_status();
    audit.final_result_published = true;
    audit.record_timing(ClassUnitStage::final_result_publication, 0.0);
    audit.record_timing(ClassUnitStage::total,
                        elapsed_milliseconds(total_start));
    return true;
}

}  // namespace silex::detail
