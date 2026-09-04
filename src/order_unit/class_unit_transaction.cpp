#include "class_unit_transaction_internal.hpp"

#include "../class_group/relation_search_driver_internal.hpp"
#include "../class_group/relation_factor_base_plan_internal.hpp"

#include <silex/signature.hpp>

#include <chrono>

namespace silex::detail {
namespace {

constexpr slong kMinimumFactorBaseSize = 20;

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
                    selected_bound, order, options.diagnostics)) {
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

}  // namespace

bool grh_factor_base_bound(
        flint::Fmpz& out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!grh_factor_base_bound_with_diagnostics(
                flint::FmpzRef(out), order, diagnostics)) {
        return false;
    }
    for (;;) {
        FactorBase probe(order);
        if (!probe.is_defined() || !probe.build_prime_ideal_norm_bounded(
                                           flint::FmpzConstRef(out))) {
            return false;
        }
        for (slong i = 0; i < probe.length(); ++i) {
            const PrimeIdeal* prime = probe.prime_at(i);
            if (prime == nullptr || prime->residue_degree() <= 0) {
                return false;
            }
        }
        if (probe.length() > kMinimumFactorBaseSize) {
            return true;
        }
        flint::fmpz_mul_2exp(flint::FmpzRef(out), flint::FmpzConstRef(out),
                             UWORD(1));
    }
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
    const bool success = candidate_units.compute_with_class_group(
            candidate_class_group, order,
            flint::FmpzConstRef(transaction_factor_base_bound),
            transaction_options, precision);
    ClassUnitTransactionAccess::set_run_context(candidate_class_group,
                                                nullptr);
    if (!success) {
        audit.failure_stage = ClassUnitStage::total;
        audit.failure_reason = "class_unit_computation_failed";
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
