#pragma once

#include <memory>
#include <vector>

#include <silex/class_group.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order_unit.hpp>

namespace silex {

class Signature;

namespace detail {

struct CompactFieldModulusCache;
struct HnfFinishWorkspace;
struct CompactRelationUnitGroupStorage;

class CompactRelationUnitGroup {
public:
    CompactRelationUnitGroup() noexcept;
    ~CompactRelationUnitGroup() noexcept;

    CompactRelationUnitGroup(const CompactRelationUnitGroup&) = delete;
    CompactRelationUnitGroup& operator=(
            const CompactRelationUnitGroup&) = delete;

    CompactRelationUnitGroup(CompactRelationUnitGroup&& other) noexcept;
    CompactRelationUnitGroup& operator=(
            CompactRelationUnitGroup&& other) noexcept;

    void swap(CompactRelationUnitGroup& other) noexcept;
    void clear() noexcept;
    bool is_set() const noexcept;
    slong free_rank() const noexcept;
    slong relation_generator_count() const noexcept;
    slong exponents_exceeding_slong() const noexcept;
    slong maximum_exponent_bits() const noexcept;
    bool exact_kernel_verified() const noexcept;
    bool regulator(flint::ArbRef out) const noexcept;
    bool class_regulator_index_bound(
            flint::FmpzRef out,
            const ClassGroupContext& class_group,
            flint::ArbConstRef analytic_class_regulator_product,
            slong precision) const noexcept;

private:
    friend bool set_hnf_compact_relation_units(
            CompactRelationUnitGroup& out,
            const Order& order,
            const ClassGroupContext& class_group,
            EmbeddingContext& embeddings,
            flint::ArbConstRef analytic_hR,
            HnfFinishWorkspace* workspace,
            slong precision) noexcept;

    std::unique_ptr<CompactRelationUnitGroupStorage> storage_;
};

enum class OrdinaryUnitCoordinateOutcome {
    unknown,
    not_unit,
    verified,
};

enum class OrdinaryUnitCoordinateStage {
    none,
    input_validation,
    exact_unit_test,
    precision_exhausted,
    exact_verification,
};

struct OrdinaryUnitCoordinates {
    flint::Fmpz torsion_exponent;
    flint::FmpzMat free_exponents{1, 0};
    slong work_precision = 0;
    bool defined = false;
};

struct OrdinaryUnitCoordinateResult {
    bool success = false;
    OrdinaryUnitCoordinateOutcome outcome =
            OrdinaryUnitCoordinateOutcome::unknown;
    OrdinaryUnitCoordinateStage stage = OrdinaryUnitCoordinateStage::none;
    slong work_precision = 0;
};

const char* ordinary_unit_coordinate_stage_name(
        OrdinaryUnitCoordinateStage stage) noexcept;
bool ordinary_unit_coordinates(
        OrdinaryUnitCoordinateResult& result,
        OrdinaryUnitCoordinates& out,
        const OrderUnitGroup& group,
        const Element& value,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept;
bool ordinary_unit_coordinates(
        OrdinaryUnitCoordinateResult& result,
        OrdinaryUnitCoordinates& out,
        const OrderUnitGroup& group,
        const FactoredElement& value,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept;

bool compute_power(Element& out,
                   const Element& base,
                   slong exponent) noexcept;
bool compact_places(slong& places, EmbeddingContext& embeddings) noexcept;
bool compact_log_matrix(flint::ArbMat& out,
                        EmbeddingContext& embeddings,
                        FactoredElementSpan generators,
                        slong precision) noexcept;
bool compact_regulator_from_log_matrix(flint::ArbRef out,
                                       const flint::ArbMat& logs,
                                       slong len,
                                       slong places,
                                       slong precision) noexcept;
bool arb_radius_lt_2exp(const flint::Arb& value, slong exponent) noexcept;
bool class_regulator_index_bound_from_candidate_product(
        flint::FmpzRef out,
        flint::ArbConstRef candidate_class_regulator_product,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;
bool class_regulator_index_is_one_from_candidate_product(
        flint::ArbConstRef candidate_class_regulator_product,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;
bool compact_independent(bool& independent,
                         EmbeddingContext& embeddings,
                         FactoredElementSpan generators,
                         slong precision) noexcept;
bool compact_independent_from_log_matrix(bool& independent,
                                         const flint::ArbMat& logs,
                                         slong len,
                                         slong places,
                                         slong precision) noexcept;
bool relation_kernel_rescale_log_matrix_for_testing(
        flint::FmpzMatRef rows,
        flint::ArbMatConstRef matrix,
        slong precision,
        slong word_bits) noexcept;
slong relation_kernel_native_word_bits() noexcept;
bool transformed_hnf_unit_coefficients_from_regulator_matrix(
        flint::FmpzMat& out,
        const ClassGroupContext& class_group,
        const flint::FmpzMat& witness_coefficients,
        const flint::FmpzMat& integer_coordinates,
        const DiagnosticsContext* diagnostics) noexcept;
bool reduce_relation_coefficients_by_log_lll(
        flint::FmpzMat& coefficients,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong& precision,
        const DiagnosticsContext* diagnostics) noexcept;
enum class RegulatorPivotOutcome {
    success,
    precision_inconclusive,
    invalid,
};
enum class RegulatorPivotFinishPath {
    independent_unit_preprobe,
    regulator_product,
};
enum class RegulatorPivotFinishAction {
    proceed,
    retry_precision,
    restart_factor_base,
    request_relations,
    failed,
};
struct RegulatorPivotFinishTestResult {
    RegulatorPivotFinishAction action =
            RegulatorPivotFinishAction::failed;
    bool goal_reached_before = true;
    bool goal_reached_after = true;
    slong relation_need_before = 0;
    slong relation_need_after = 0;
    slong analytic_extra_relation_requests_before = 0;
    slong analytic_extra_relation_requests_after = 0;
    slong analytic_finish_precision_before = 0;
    slong analytic_finish_precision_after = 0;
    slong analytic_precision_doublings_before = 0;
    slong analytic_precision_doublings_after = 0;
    slong analytic_precision_inconclusive_before = 0;
    slong analytic_precision_inconclusive_after = 0;
    slong factor_base_restart_requests_before = 0;
    slong factor_base_restart_requests_after = 0;
    bool factor_base_restart_pending_before = false;
    bool factor_base_restart_pending_after = false;
    bool factor_base_restart_allow_past_half_before = false;
    bool factor_base_restart_allow_past_half_after = false;
    bool relation_control_state_unchanged = false;
    bool all_control_state_unchanged = false;
    slong completion_old_need_before = 0;
    slong completion_old_need_after = 0;
    slong completion_dependent_trials_before = 0;
    slong completion_dependent_trials_after = 0;
    slong completion_subfactor_base_trials_before = 0;
    slong completion_subfactor_base_trials_after = 0;
    slong finish_unit_log_rotation_before = 0;
    slong finish_unit_log_rotation_after = 0;
    slong squash_index_before = 0;
    slong squash_index_after = 0;
    slong candidates_tried_before = 0;
    slong candidates_tried_after = 0;
    slong accepted_relations_before = 0;
    slong accepted_relations_after = 0;
    bool finish_unit_log_rotation_active_before = false;
    bool finish_unit_log_rotation_active_after = false;
    bool finish_full_rank_relation_active_before = false;
    bool finish_full_rank_relation_active_after = false;
};
RegulatorPivotOutcome regulator_pivot_unit_indices(
        std::vector<slong>& out,
        const flint::ArbMat& unit_logs,
        const Signature& sig,
        slong rank,
        slong places,
        slong precision,
        std::vector<slong>* pivot_rows_out = nullptr) noexcept;
RegulatorPivotOutcome regulator_pivot_rows_for_testing(
        std::vector<slong>& out,
        flint::ArbMatConstRef matrix,
        slong precision) noexcept;
RegulatorPivotOutcome regulator_pivot_unit_indices_for_testing(
        std::vector<slong>& out,
        std::vector<slong>& pivot_rows,
        const flint::ArbMat& unit_logs,
        const Signature& sig,
        slong rank,
        slong places,
        slong precision) noexcept;
bool regulator_ubound_lt(const flint::Arb& left,
                         const flint::Arb& right,
                         slong precision) noexcept;
bool evaluated_is_order_unit(const Order& order,
                             const Element& value) noexcept;

bool validate_relation_kernel_inputs(const Order& order,
                                     const ClassGroupContext& class_group,
                                     EmbeddingContext& embeddings,
                                     slong precision) noexcept;
bool relation_kernel_generator(FactoredElement& out,
                               const Order& order,
                               const ClassGroupContext& class_group,
                               slong index) noexcept;
bool relation_basis(std::vector<FactoredElement>& out,
                    const OrderUnitGroup& group,
                    const FactoredElement& root,
                    flint::FmpzMatConstRef relation,
                    flint::FmpzConstRef exponent) noexcept;
bool kernel_row_divisible(flint::FmpzMatConstRef kernel_rows,
                          slong row,
                          slong len,
                          flint::FmpzConstRef ell) noexcept;
bool kernel_row_root(bool& is_power,
                     FactoredElement& root,
                     const OrderUnitGroup& group,
                     flint::FmpzMatConstRef kernel_rows,
                     slong row,
                     flint::FmpzConstRef ell,
                     CompactFieldModulusCache* field_modulus_cache =
                             nullptr) noexcept;
bool residue_dlog_kernel(flint::FmpzMat& out,
                         const OrderUnitGroup& group,
                         PrimeIdealSpan primes,
                         flint::FmpzConstRef ell) noexcept;
bool residue_dlog_proof_kernel(flint::FmpzMat& out,
                               const OrderUnitGroup& group,
                               PrimeIdealSpan primes,
                               flint::FmpzConstRef ell) noexcept;
bool dlog_proof_rank(slong& out,
                     const OrderUnitGroup& group,
                     flint::FmpzConstRef ell) noexcept;
bool residue_dlog_proof_matrix(flint::FmpzMat& out,
                               const OrderUnitGroup& group,
                               const PrimeIdeal& prime,
                               flint::FmpzConstRef ell) noexcept;
bool residue_dlog_matrix_direct_degree_one(flint::FmpzMat& out,
                                           const OrderUnitGroup& group,
                                           const PrimeIdeal& prime,
                                           flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_column_direct_degree_one(
        flint::FmpzMat& out,
        const OrderUnitGroup& group,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_column_at_degree_one_root(
        flint::FmpzMat& out,
        const OrderUnitGroup& group,
        flint::FmpzConstRef p,
        flint::FmpzConstRef root,
        flint::FmpzConstRef ell) noexcept;
bool dlog_kernel_from_matrix(flint::FmpzMat& out,
                             const flint::FmpzMat& matrix,
                             flint::FmpzConstRef ell) noexcept;
bool saturation_prime_usable(const OrderUnitGroup& group,
                             const PrimeIdeal& prime,
                             flint::FmpzConstRef ell) noexcept;
bool saturation_prime_column(flint::FmpzMat& out,
                             const OrderUnitGroup& group,
                             const PrimeIdeal& prime,
                             flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_usable(const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_column(flint::FmpzMat& out,
                                   const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_rank_one_nonzero(bool& nonzero,
                                             const OrderUnitGroup& group,
                                             const PrimeIdeal& prime,
                                             flint::FmpzConstRef ell) noexcept;
bool rank_one_factored_image_nonzero_at_degree_one_root(
        bool& nonzero,
        const FactoredElement& generator,
        flint::FmpzConstRef p,
        flint::FmpzConstRef root,
        flint::FmpzConstRef ell) noexcept;
bool rank_one_factored_image_nonzero_at_degree_one_root_nmod(
        bool& nonzero,
        const FactoredElement& generator,
        ulong p,
        ulong root,
        flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_known_rank_one_free_generator_nonzero(
        bool& nonzero,
        const OrderUnitGroup& group,
        const FactoredElement& generator,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept;
bool saturation_proof_prime_known_rank_one_torsion_nonzero(
        bool& nonzero,
        const OrderUnitGroup& group,
        const OrderElement& torsion,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept;
bool copy_selected_primes(PrimeIdealList& out,
                          const Order& order,
                          const std::vector<PrimeIdeal>& selected) noexcept;
bool adjoin_verified_dependent_relation(bool& changed,
                                        OrderUnitGroup& out,
                                        const OrderUnitGroup& group,
                                        const FactoredElement& root,
                                        flint::FmpzMatConstRef rel,
                                        slong row,
                                        EmbeddingContext& embeddings,
                                        slong precision) noexcept;
bool dependent_relation_bounded_min_denominator(
        bool& recovered,
        FactoredElement& root,
        flint::FmpzMat& rel,
        flint::Fmpz& torsion_exp,
        const OrderUnitGroup& group,
        const FactoredElement& y,
        EmbeddingContext& embeddings,
        flint::FmpzConstRef denominator_bound,
        slong min_denominator,
        slong start_precision,
        slong max_precision,
        bool require_y_root,
        bool require_torsion_exponent) noexcept;
bool dependent_relation_bounded_with_inverse(
        bool& recovered,
        FactoredElement& root,
        flint::FmpzMat& rel,
        flint::Fmpz& torsion_exp,
        const OrderUnitGroup& group,
        const FactoredElement& y,
        EmbeddingContext& embeddings,
        const flint::ArbMat& inverse_cutoff,
        flint::FmpzConstRef denominator_bound,
        slong min_denominator,
        slong precision,
        bool require_y_root,
        bool require_torsion_exponent) noexcept;
bool set_initial_relation_kernel_units(OrderUnitGroup& out,
                                       const Order& order,
                                       const ClassGroupContext& class_group,
                                       EmbeddingContext& embeddings,
                                       slong precision) noexcept;
bool relation_kernel_independent_unit_count(slong& out,
                                            const Order& order,
                                            const ClassGroupContext& class_group,
                                            EmbeddingContext& embeddings,
                                            slong precision) noexcept;
bool relation_coefficients_log_matrix(
        flint::ArbMat& out,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        const flint::FmpzMat& coefficients,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;
bool reduced_regulator_from_coordinates(
        flint::ArbRef out,
        const flint::ArbMat& coordinates,
        flint::ArbConstRef regulator_multiple,
        flint::ArbConstRef z,
        flint::FmpzMat* integer_coordinates_out,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;
bool hnf_regulator_multiple_from_basis(
        flint::ArbRef out,
        const flint::ArbMat& basis,
        slong degree,
        slong precision) noexcept;
RegulatorPivotOutcome
hnf_class_regulator_product_and_unit_matrix(
        flint::ArbRef out,
        slong* independent_count_out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        flint::FmpzMat* unit_matrix_out,
        flint::Arb* regulator_out,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept;
bool factored_units_from_relation_coefficients(
        std::vector<FactoredElement>& out,
        const ClassGroupContext& class_group,
        const flint::FmpzMat& coefficients,
        const DiagnosticsContext* diagnostics) noexcept;
bool unit_log_row_sums_are_small(const flint::ArbMat& logs,
                                         slong precision) noexcept;
bool unit_log_row_sums_are_small(
        const std::vector<FactoredElement>& units,
        EmbeddingContext& embeddings,
        slong precision) noexcept;
bool unit_regulator_matches_reconstruction(
        flint::ArbConstRef actual_regulator,
        flint::ArbConstRef expected_regulator,
        slong precision) noexcept;
bool unit_regulator_matches_reconstruction(
        const OrderUnitGroup& units,
        flint::ArbConstRef expected_regulator,
        slong precision) noexcept;
bool publish_validated_hnf_units(
        OrderUnitGroup& out,
        const Order& order,
        std::vector<FactoredElement>& units,
        EmbeddingContext& embeddings,
        flint::ArbConstRef expected_regulator,
        slong rank,
        slong precision,
        const flint::Fmpz* cached_torsion_order,
        const OrderElement* cached_torsion_generator) noexcept;
bool hnf_unit_log_matrix(flint::ArbMat& out,
                              const ClassGroupContext& class_group,
                              EmbeddingContext& embeddings,
                              HnfFinishWorkspace* workspace,
                              slong precision) noexcept;
RegulatorPivotOutcome hnf_independent_unit_count(
        slong& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept;
bool set_hnf_units(OrderUnitGroup& out,
                        const Order& order,
                        const ClassGroupContext& class_group,
                        EmbeddingContext& embeddings,
                        flint::ArbConstRef analytic_hR,
                        HnfFinishWorkspace* workspace,
                        slong precision) noexcept;
bool set_hnf_compact_relation_units(
        CompactRelationUnitGroup& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept;
RegulatorPivotOutcome hnf_class_regulator_product(
        flint::ArbRef out,
        slong& independent_count,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept;

namespace relation_search {

RegulatorPivotFinishTestResult
regulator_pivot_finish_for_testing(
        RegulatorPivotFinishPath path,
        RegulatorPivotOutcome outcome,
        slong independent_count,
        slong target_rank,
        bool factor_base_restart_available,
        slong analytic_finish_precision,
        slong analytic_precision_doublings) noexcept;

}  // namespace relation_search

}  // namespace detail

}  // namespace silex
