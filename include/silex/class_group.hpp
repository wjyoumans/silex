#pragma once

#include <memory>
#include <optional>
#include <vector>

#include <flint/flint.h>

#include <silex/abelian_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/detail/class_relation_module_context.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/fmpz_smat.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>
#include <silex/status.hpp>

namespace silex {

class OrderUnitGroup;
class ClassGroupContext;
class FactoredElement;

namespace detail {
struct ClassUnitTransactionContext;
struct ClassGroupRelationOptions;
struct ClassGroupContextStorage;
class ClassGroupCertificationAccess;
class ClassGroupRelationAccess;
class ClassGroupFinishAccess;
class ClassGroupRelationSearchAccess;
class ClassUnitTransactionAccess;
}  // namespace detail

struct ClassGroupFactorBaseGenerationRecord {
    flint::Fmpz p;
    ProofState status = ProofState::not_checked;
};

struct ClassGroupRelationSaturationRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
};

struct ClassGroupZetaBfProofRecord {
    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
    flint::Arb error_bound;
};

// Resource ceilings for the deterministic standalone relation candidate.
// Successful candidates always retain CertificationMode::unknown.
struct ClassGroupCandidateOptions {
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    const DiagnosticsContext* diagnostics = nullptr;
};

// Proof resources for the atomic paired class/unit transaction.
struct ClassGroupComputeOptions {
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    slong relation_saturation_aux_prime_bound = 0;
    slong relation_saturation_max_appends_per_ell = 0;
    ulong zeta_bf_max_cutoff = 0;
    CertificationMode requested_certification = CertificationMode::unknown;
    const DiagnosticsContext* diagnostics = nullptr;
};

enum class ClassGroupRelationSource {
    Unknown = 0,
    Search = 1,
    RandomProduct = 2,
    Supplied = 3,
    Saturation = 4,
    LargePrimeMatch = 5,
    ClassGenerator = 6
};

class ClassGroupContext {
public:
    ClassGroupContext() noexcept;
    explicit ClassGroupContext(const Order& order) noexcept;
    ~ClassGroupContext() noexcept;

    ClassGroupContext(const ClassGroupContext&) = delete;
    ClassGroupContext& operator=(const ClassGroupContext&) = delete;

    ClassGroupContext(ClassGroupContext&& other) noexcept;
    ClassGroupContext& operator=(ClassGroupContext&& other) noexcept;

    void swap(ClassGroupContext& other) noexcept;
    void clear() noexcept;
    bool define(const Order& order) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    void set_diagnostics(const DiagnosticsContext* diagnostics) noexcept;
    const DiagnosticsContext* diagnostics() const noexcept;
    slong analytic_finish_precision() const noexcept;
    bool analytic_finish_product(flint::ArbRef out,
                                      slong& precision) const noexcept;

    bool has_factor_base() const noexcept;
    const FactorBase* factor_base() const noexcept;
    slong generator_count() const noexcept;
    bool factor_base_prime(PrimeIdeal& out, slong index) const noexcept;
    bool factor_base_prime_is_principal(bool& out,
                                        slong index) const noexcept;
    bool factor_base_prime_is_hnf_covered(bool& out,
                                          slong index) noexcept;
    bool build_factor_base(flint::FmpzConstRef bound) noexcept;
    bool set_factor_base(const FactorBase& base) noexcept;

    bool append_relation(const Relation& relation) noexcept;
    bool append_relation(const Relation& relation,
                         ClassGroupRelationSource source) noexcept;
    bool try_append_generator_relation(
            const Element& generator,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept;
    bool try_append_generator_relation(
            bool& partial_throttle_exit,
            const Element& generator,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept;
    // `norm` must be the exact field norm of `generator`.
    bool try_append_generator_relation_with_norm(
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpqConstRef norm,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept;
    slong relation_count() const noexcept;
    slong relation_rank() const noexcept;
    slong skipped_dependent_relation_count() const noexcept;
    bool relation_source(ClassGroupRelationSource& out,
                         slong index) const noexcept;
    slong relation_source_count(ClassGroupRelationSource source) const noexcept;
    bool relations(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::FmpzMat> relations() const noexcept;
    bool relation_generator(Element& out, slong index) const noexcept;
    const Element* relation_generator_at(slong index) const noexcept;

    bool set_relation_matrix(const RelationMatrix& matrix) noexcept;
    bool publish_presentation() noexcept;
    bool has_presentation() const noexcept;
    CertificationMode certification_status() const noexcept;
    ProofState factor_base_generation_status() const noexcept;
    ProofState factor_base_generation_checked_status() const noexcept;
    bool factor_base_build_bound(flint::FmpzRef out) const noexcept;
    bool factor_base_generation_bound(flint::FmpzRef out) const noexcept;
    bool factor_base_generation_checked_bound(
            flint::FmpzRef out) const noexcept;
    slong factor_base_generation_record_count() const noexcept;
    bool factor_base_generation_record(flint::FmpzRef p,
                                       ProofState& status,
                                       slong index) const noexcept;
    std::optional<ClassGroupFactorBaseGenerationRecord>
    factor_base_generation_record(slong index) const noexcept;
    bool check_factor_base_generation_bound(
            flint::FmpzConstRef required_bound) noexcept;
    ProofState relation_saturation_status() const noexcept;
    ProofState analytic_class_regulator_status() const noexcept;
    ProofState zeta_bf_proof_status() const noexcept;
    bool zeta_bf_proof_record(ulong& cutoff,
                              ulong& max_cutoff,
                              slong& requested_precision,
                              slong& work_precision,
                              flint::ArbRef error_bound) const noexcept;
    std::optional<ClassGroupZetaBfProofRecord>
    zeta_bf_proof_record() const noexcept;
    slong relation_saturation_record_count() const noexcept;
    bool relation_saturation_record(flint::FmpzRef ell,
                                    ProofState& status,
                                    slong index) const noexcept;
    std::optional<ClassGroupRelationSaturationRecord>
    relation_saturation_record(slong index) const noexcept;
    ProofState unit_proof_status() const noexcept;
    ProofState regulator_proof_status() const noexcept;
    bool try_certify_quadratic(CertificationMode requested) noexcept;
    bool try_certify_trivial_quotient(CertificationMode requested) noexcept;
    bool try_certify_with_units(const OrderUnitGroup& units,
                                CertificationMode requested,
                                slong precision) noexcept;
    bool try_certify_with_units(const OrderUnitGroup& units,
                                CertificationMode requested,
                                slong precision,
                                ulong zeta_bf_max_cutoff) noexcept;
    bool try_certify_class_unit_with_units(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            slong precision) noexcept;
    bool try_certify_class_unit_with_zeta(OrderUnitGroup& units,
                                          slong precision) noexcept;
    bool try_certify_class_unit_with_zeta_bf(OrderUnitGroup& units,
                                             ulong max_cutoff,
                                             slong precision) noexcept;
    bool try_prove_relation_saturation_with_units(
            const OrderUnitGroup& units,
            flint::FmpzConstRef ell,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool try_prove_relation_saturation_index_bound_with_units(
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool try_analytic_index_bound_with_units(
            const OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            flint::FmpzConstRef aux_prime_bound,
            slong precision) noexcept;
    bool saturate_relations_bounded_with_units(
            bool& changed,
            bool& saturated,
            const OrderUnitGroup& units,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;
    bool presentation(FiniteAbelianGroup& out) const noexcept;
    std::optional<FiniteAbelianGroup> presentation() const noexcept;
    slong invariant_count() const noexcept;
    bool invariant(flint::FmpzRef out, slong index) const noexcept;
    bool invariants(flint::FmpzVecRef out) const noexcept;
    bool order(flint::FmpzRef out) const noexcept;
    bool invariant_generator_matrix(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::Fmpz> invariant(slong index) const noexcept;
    std::optional<flint::Fmpz> order() const noexcept;
    std::optional<flint::FmpzMat> invariant_generator_matrix() const noexcept;
    bool invariant_generator(FractionalIdeal& out, slong index) const noexcept;
    bool ideal_class_coordinates(flint::FmpzMatRef out,
                                 const FractionalIdeal& ideal) const noexcept;
    bool invariant_generator_power_witness(FactoredElement& out,
                                           slong index) const noexcept;
    slong relation_kernel_unit_count() const noexcept;
    bool relation_kernel_unit(FactoredElement& out, slong index) const noexcept;
    bool compute_candidate(const Order& order,
                           flint::FmpzConstRef factor_base_bound,
                           const ClassGroupCandidateOptions& options = {})
            noexcept;

private:
    friend class OrderUnitGroup;
    friend class detail::ClassGroupCertificationAccess;
    friend class detail::ClassGroupRelationAccess;
    friend class detail::ClassGroupFinishAccess;
    friend class detail::ClassGroupRelationSearchAccess;
    friend class detail::ClassUnitTransactionAccess;

    enum class RelationAppendOutcome {
        none,
        rank,
        index,
        kernel,
        skipped_dependent
    };

    enum class DependentRelationPolicy {
        native,
        skip,
        keep,
        keep_nonduplicate
    };

    enum class SpecializedRelationBackendStatus {
        unavailable,
        succeeded,
        failed
    };

    bool prove_relation_saturation_dlog_ell_(
            const OrderUnitGroup& units,
            flint::FmpzConstRef ell,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool try_promote_proven_certification_() noexcept;
    bool try_certify_analytic_class_regulator_(
            const OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision) noexcept;
    bool record_analytic_class_unit_regulator_(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision) noexcept;
    bool try_certify_analytic_class_unit_regulator_(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision) noexcept;
    bool record_zeta_bf_audit_(flint::ArbConstRef error_bound,
                               ulong cutoff,
                               ulong max_cutoff,
                               slong requested_precision,
                               slong work_precision) noexcept;
    bool relation_row_refines_(const Relation& relation) noexcept;
    bool sync_row_module_checkpoint_() noexcept;
    bool append_relation_with_outcome_(RelationAppendOutcome& outcome,
                                       const Relation& relation,
                                       ClassGroupRelationSource source)
            noexcept;
    bool append_relation_with_outcome_(
            RelationAppendOutcome& outcome,
            const Relation& relation,
            ClassGroupRelationSource source,
            DependentRelationPolicy dependent_policy) noexcept;
    SpecializedRelationBackendStatus
    run_maximal_imaginary_quadratic_relation_backend_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool build_relation_factor_base_(flint::FmpzConstRef bound) noexcept;
    bool build_search_factor_base_(
            flint::FmpzConstRef bound,
            bool incomplete = false) noexcept;
    void reset_partial_relations_() noexcept;
    void configure_partial_relations_(
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool try_append_integral_generator_relation_(
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            ClassGroupRelationSource source) noexcept;
    bool try_partial_relation_(bool& partial_throttle_exit,
                               const Element& generator,
                               const fmpq* known_norm = nullptr,
                               const fmpz_mat_struct*
                                       known_integral_coordinates = nullptr,
                               const fmpz_poly_struct*
                                       known_integral_polynomial = nullptr)
            noexcept;
    bool saturate_local_once_(bool& changed,
                              const PrimeIdeal& prime,
                              flint::FmpzConstRef ell) noexcept;
    bool saturate_local_(bool& changed,
                         bool& index_cleared,
                         flint::FmpzConstRef ell,
                         flint::FmpzConstRef aux_prime_bound,
                         slong max_appends) noexcept;
    bool saturate_relations_bounded_(bool& changed,
                                     bool& saturated,
                                     flint::FmpzConstRef aux_prime_bound,
                                     slong max_appends_per_ell,
                                     slong max_appends_total) noexcept;
    bool try_auto_relation_saturation_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool* changed_out = nullptr) noexcept;
    bool run_native_experimental_relation_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_relation_search_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event,
            bool clamp_relation_kernel_units_to_rank) noexcept;
    bool run_relation_production_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_relation_production_prepass_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_lll_relation_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event,
            bool accept_tentative_presentation = false,
            bool extension_only = false) noexcept;
    bool pivot_info_(flint::FmpzRef h,
                           std::vector<slong>& pivots) noexcept;
    bool compute_tentative_candidate_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool extend_tentative_relations_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool extend_lll_relation_slice_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool saturate_local_once_with_units_(bool& changed,
                                         const OrderUnitGroup& units,
                                         bool include_torsion,
                                         const PrimeIdeal& prime,
                                         flint::FmpzConstRef ell) noexcept;
    bool saturate_local_with_units_(bool& changed,
                                    bool& index_cleared,
                                    const OrderUnitGroup& units,
                                    flint::FmpzConstRef ell,
                                    flint::FmpzConstRef aux_prime_bound,
                                    slong max_appends) noexcept;
    bool saturate_relations_bounded_for_index_with_units_(
            bool& changed,
            bool& saturated,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;
    bool extend_relation_kernel_units_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool compute_relation_candidate_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;

    void reset_certification_metadata_() noexcept;
    void reset_factor_base_generation_check_() noexcept;
    bool mark_factor_base_generation_check_(
            flint::FmpzConstRef checked_bound,
            ProofState status) noexcept;
    bool record_factor_base_generation_(
            flint::FmpzConstRef build_bound) noexcept;
    bool record_factor_base_honesty_proof_(
            flint::FmpzConstRef required_bound) noexcept;
    bool relation_saturation_proof_prereqs_verified_() const noexcept;
    bool record_relation_saturation_proof_(flint::FmpzConstRef ell,
                                           ProofState status,
                                           slong rank,
                                           slong target_rank,
                                           slong local_primes) noexcept;
    bool relation_saturation_proof_verified_(
            flint::FmpzConstRef ell) const noexcept;
    bool mark_relation_saturation_verified_(
            flint::FmpzConstRef ell) noexcept;
    bool mark_relation_saturation_(flint::FmpzConstRef ell,
                                   ProofState status) noexcept;
    bool complete_relation_saturation_proof_(
            flint::FmpzConstRef ell) noexcept;
    bool complete_relation_saturation_proof_(
            const std::vector<flint::Fmpz>& required_ells) noexcept;
    bool relation_saturation_proof_complete_() const noexcept;
    bool quadratic_completeness_verified_() const noexcept;
    bool ensure_private_storage_() noexcept;
    Order parent_;
    FactorBase base_;
    RelationMatrix relations_;
    std::vector<ClassGroupRelationSource> relation_sources_;
    std::vector<char> relation_basis_flags_;
    Relation generator_relation_scratch_;
    // Allocated transactionally by define(); null only while undefined.
    std::unique_ptr<detail::ClassGroupContextStorage> private_storage_;
    FiniteAbelianGroup quotient_;
    fmpz_smat::HnfContext row_module_;
    detail::ClassRelationModuleContext relation_module_;
    slong relation_rank_ = 0;
    bool row_module_synced_ = true;
    slong skipped_dependent_relations_ = 0;
    slong relation_kernel_units_target_ = 0;
    DependentRelationPolicy generator_relation_policy_ =
            DependentRelationPolicy::native;
    CertificationMode certification_ = CertificationMode::unknown;
    flint::Fmpz factor_base_build_bound_;
    flint::Fmpz factor_base_generation_bound_;
    flint::Fmpz factor_base_generation_checked_bound_;
    ProofState factor_base_generation_status_ = ProofState::not_checked;
    ProofState factor_base_generation_checked_status_ =
            ProofState::not_checked;
    ProofState relation_saturation_status_ = ProofState::not_checked;
    ProofState analytic_class_regulator_status_ = ProofState::not_checked;
    ProofState zeta_bf_status_ = ProofState::not_checked;
    ulong zeta_bf_cutoff_ = 0;
    ulong zeta_bf_max_cutoff_ = 0;
    slong zeta_bf_requested_precision_ = 0;
    slong zeta_bf_work_precision_ = 0;
    flint::Arb zeta_bf_error_bound_;
    const DiagnosticsContext* diagnostics_ = nullptr;
    // Driver-scoped and non-owning; intentionally not swapped with results.
    detail::ClassUnitTransactionContext* class_unit_transaction_context_ =
            nullptr;
    ProofState unit_proof_status_ = ProofState::not_checked;
    ProofState regulator_proof_status_ = ProofState::not_checked;
    bool has_base_ = false;
};

inline void swap(ClassGroupContext& left,
                 ClassGroupContext& right) noexcept {
    left.swap(right);
}

}  // namespace silex
