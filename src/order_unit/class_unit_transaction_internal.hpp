#pragma once

#include <array>
#include <cstddef>

#include <silex/class_group.hpp>
#include <silex/order_unit.hpp>

namespace silex::detail {

struct RelationFactorBasePlan;
struct ZetaBfResidueDegreeCache;

enum class ClassUnitStage : unsigned char {
    none = 0,
    factor_base_bound,
    factor_base_construction,
    relation_scheduling,
    relation_collection,
    matrix_hnf_publication,
    class_group_publication,
    unit_extraction,
    regulator_computation,
    saturation,
    analytic_bf_zeta_validation,
    relation_extension,
    certification,
    final_result_publication,
    total,
    count,
};

enum class NativeFactorBaseBoundStrategy : unsigned char {
    none = 0,
    requested = 1,
    relation_completion = 2,
    grh = 3,
};

enum class NativeFactorBaseStrategy : unsigned char {
    none = 0,
    complete = 1,
    relation_completion = 2,
    norm_bounded = 3,
};

enum class NativeRelationStrategy : unsigned char {
    none = 0,
    native_experimental = 1,
    relation_completion_prepass = 3,
    relation_completion_table = 4,
    lll = 5,
};

enum class NativeUnitStrategy : unsigned char {
    none = 0,
    rank_zero = 1,
    quadratic_rank_one = 2,
    class_relation = 4,
    hnf = 5,
};

enum class NativeContinuationStrategy : unsigned char {
    none = 0,
    native_target_growth = 1,
    lll = 2,
    regulator_multiple_reconstruction = 3,
};

enum class NativeValidationStrategy : unsigned char {
    none = 0,
    shared_exact = 1,
};

struct ClassUnitExecutionPolicy {
    bool selected = false;
    slong degree = 0;
    slong unit_rank = -1;
    slong discriminant_bits = 0;
    flint::Fmpz selected_factor_base_bound;
    NativeFactorBaseBoundStrategy factor_base_bound =
            NativeFactorBaseBoundStrategy::none;
    NativeFactorBaseStrategy factor_base = NativeFactorBaseStrategy::none;
    NativeRelationStrategy relations = NativeRelationStrategy::none;
    NativeUnitStrategy units = NativeUnitStrategy::none;
    NativeContinuationStrategy continuation =
            NativeContinuationStrategy::none;
    NativeValidationStrategy validation = NativeValidationStrategy::none;
};

inline constexpr std::size_t kClassUnitStageCount =
        static_cast<std::size_t>(ClassUnitStage::count);

struct ClassUnitTransactionReport {
    ClassUnitStage failure_stage = ClassUnitStage::none;
    const char* failure_reason = nullptr;
    std::array<double, kClassUnitStageCount> timing_ms{};
    std::array<bool, kClassUnitStageCount> timing_recorded{};
    ClassUnitExecutionPolicy policy{};
    CertificationMode class_group_certification =
            CertificationMode::unknown;
    CertificationMode unit_group_certification =
            CertificationMode::unknown;
    ProofState unit_proof_status = ProofState::not_checked;
    ProofState regulator_proof_status = ProofState::not_checked;
    bool final_result_published = false;

    void reset() noexcept;
    void record_timing(ClassUnitStage stage, double elapsed_ms) noexcept;
};

struct ClassUnitTransactionContext {
    ClassUnitTransactionReport& audit;
    RelationFactorBasePlan* relation_factor_base_plan = nullptr;
    flint::Fmpz imaginary_quadratic_exact_order_discriminant;
    flint::Fmpz imaginary_quadratic_exact_order;
    bool imaginary_quadratic_exact_order_valid = false;
    bool defer_relation_saturation_until_units = false;
};

inline bool uses_class_unit_kernel(
        const ClassUnitTransactionContext* context) noexcept {
    return context != nullptr &&
           context->audit.policy.selected &&
           context->audit.policy.relations ==
                   NativeRelationStrategy::lll;
}

class ClassUnitTransactionAccess {
public:
    static void set_run_context(
            ClassGroupContext& context,
            ClassUnitTransactionContext* run_context) noexcept;
    static ZetaBfResidueDegreeCache*
    relation_factor_base_plan_residue_degrees(
            ClassGroupContext& context,
            const Order& order) noexcept;
};

const char* class_unit_stage_name(ClassUnitStage stage) noexcept;
bool native_hnf_unit_strategy_for_signature(
        slong degree,
        slong unit_rank) noexcept;
bool grh_factor_base_bound(
        flint::Fmpz& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool compute_class_unit_transaction(
        OrderUnitGroup& units,
        ClassGroupContext& class_group,
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupComputeOptions& options,
        slong precision,
        ClassUnitTransactionReport& audit) noexcept;

}  // namespace silex::detail
