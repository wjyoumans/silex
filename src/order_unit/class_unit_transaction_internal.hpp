#pragma once

#include <array>
#include <cstddef>

#include <silex/class_group.hpp>
#include <silex/order_unit.hpp>

#include "../class_group/factor_base_honesty_failure_internal.hpp"
#include "../factor_base/generation_different_route_internal.hpp"

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
    // grh policy only: the GRH generation bound min(BDF, Bach), or Bach's
    // bound alone when the BDF criterion cannot be computed or decided,
    // before the selected bound is doubled to reach the minimum factor-base
    // size.
    // Under GRH every prime ideal of norm at most this bound generates the
    // class group; zero on every other policy.
    flint::Fmpz grh_generation_bound;
    // grh policy only: the theorem that gave grh_generation_bound, bdf
    // (Belabas-Diaz y Diaz-Friedman 2008, Cor 2.2) or bach (Bach 1990,
    // Thm 4); none on every other policy.
    ClassGroupGenerationBasis grh_generation_basis =
            ClassGroupGenerationBasis::none;
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
    // The unwitnessed required prime that ended the most recent class-group
    // candidate computation or extension, if the relation-search honesty
    // checkpoint ended it with no factor-base restart left.  Cleared when
    // each attempt starts; read only to name a failed transaction's reason.
    FactorBaseHonestyFailure factor_base_honesty_failure{};
    // Name (a static string) of the last unit-group validation outcome of
    // the candidate computation, or nullptr when validation never ran.  Read
    // only to describe a class_unit_computation_failed transaction.
    const char* last_validation_outcome = nullptr;
    // Whether the most recent class-group candidate attempt ended because its
    // factor base had more generators than the relation modules allow
    // (kRelationModuleMaxColumns, relation_module_limits_internal.hpp).
    // Cleared with factor_base_honesty_failure; read only to name a failed
    // transaction's reason.
    bool factor_base_size_limit_exceeded = false;

    void reset() noexcept;
    void record_timing(ClassUnitStage stage, double elapsed_ms) noexcept;
};

// proven_generation_requirement(bound, route, order) of the most recently
// evaluated order of a class/unit transaction, kept so that later
// factor-base builds of the same order reuse it; a build of another order
// evaluates and replaces it (ClassGroupContext::record_factor_base_generation_).
// The requirement depends only on the order, which `order` holds; nothing is
// kept while `order` is undefined.
struct TransactionGenerationRequirement {
    Order order;
    flint::Fmpz bound;
    GenerationDifferentRoute route = GenerationDifferentRoute::none;
};

struct ClassUnitTransactionContext {
    ClassUnitTransactionReport& audit;
    RelationFactorBasePlan* relation_factor_base_plan = nullptr;
    flint::Fmpz imaginary_quadratic_exact_order_discriminant;
    flint::Fmpz imaginary_quadratic_exact_order;
    bool imaginary_quadratic_exact_order_valid = false;
    bool defer_relation_saturation_until_units = false;
    TransactionGenerationRequirement generation_requirement{};
    // grh policy only: the norm-bounded factor base grh_factor_base_bound
    // built at the selected bound while choosing it, and that bound.  The
    // strict search-base build at that bound takes the base over instead of
    // decomposing every rational prime again (see
    // grh_norm_bounded_base_for); it is left empty once taken.
    FactorBase* grh_norm_bounded_base = nullptr;
    const flint::Fmpz* grh_norm_bounded_base_bound = nullptr;
};

inline bool uses_class_unit_kernel(
        const ClassUnitTransactionContext* context) noexcept {
    return context != nullptr &&
           context->audit.policy.selected &&
           context->audit.policy.relations ==
                   NativeRelationStrategy::lll;
}

// The base kept in `context` by grh_factor_base_bound when it can stand in
// for FactorBase::build_prime_ideal_norm_bounded(bound) on `order`: it is
// defined, belongs to the same order data, and was built at exactly
// `bound`.  nullptr otherwise.
FactorBase* grh_norm_bounded_base_for(
        ClassUnitTransactionContext* context,
        const Order& order,
        flint::FmpzConstRef bound) noexcept;

class ClassUnitTransactionAccess {
public:
    static void set_run_context(
            ClassGroupContext& context,
            ClassUnitTransactionContext* run_context) noexcept;
    // Records the validation outcome name (static string) in the run
    // context's report, if any.
    static void record_validation_outcome(ClassGroupContext& context,
                                          const char* outcome_name) noexcept;
    static ZetaBfResidueDegreeCache*
    relation_factor_base_plan_residue_degrees(
            ClassGroupContext& context,
            const Order& order) noexcept;
};

const char* class_unit_stage_name(ClassUnitStage stage) noexcept;
// The failure_reason of a transaction whose class/unit computation failed:
// factor_base_size_limit_exceeded when the last candidate attempt ended
// because its factor base exceeded the relation module size limit,
// factor_base_honesty_unwitnessed when it ended at an unwitnessed required
// prime with no factor-base restart left, otherwise
// class_unit_computation_failed.
const char* class_unit_computation_failure_reason(
        const ClassUnitTransactionReport& audit) noexcept;
bool native_hnf_unit_strategy_for_signature(
        slong degree,
        slong unit_rank) noexcept;
// `out` is the GRH generation bound min(BDF, Bach) (Bach's alone when the
// BDF criterion cannot be computed or decided), doubled until the
// norm-bounded factor base has more than the minimum number of primes.
// `generation_bound`, when given, receives the undoubled GRH bound, and
// `generation_basis`, when given, the theorem that gave it (bdf or bach).
// `norm_bounded_base`, when given, receives on success the factor base built
// by build_prime_ideal_norm_bounded(out) for that size test (every prime
// ideal in it checked to have positive residue degree); it is left
// unchanged on failure.
bool grh_factor_base_bound(
        flint::Fmpz& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr,
        flint::Fmpz* generation_bound = nullptr,
        ClassGroupGenerationBasis* generation_basis = nullptr,
        FactorBase* norm_bounded_base = nullptr) noexcept;

bool compute_class_unit_transaction(
        OrderUnitGroup& units,
        ClassGroupContext& class_group,
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupComputeOptions& options,
        slong precision,
        ClassUnitTransactionReport& audit) noexcept;

}  // namespace silex::detail
