#pragma once

#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>

#include <vector>

#include "class_group_storage_internal.hpp"
#include "ideal_minkowski_embedding_internal.hpp"

namespace silex::detail {

struct RelationUnitCandidateBatchState;

struct ClassGroupRelationOptions {
    slong coordinate_search_radius = 2;
    slong ideal_search_radius = 1;
    slong target_relation_kernel_units = 0;
    slong post_finite_refinement_phase_budget = 0;
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    slong relation_saturation_aux_prime_bound = 0;
    slong relation_saturation_max_appends_per_ell = 0;
    ulong zeta_bf_max_cutoff = 0;
    CertificationMode requested_certification = CertificationMode::unknown;
    const DiagnosticsContext* diagnostics = nullptr;
};

inline ClassGroupRelationOptions class_group_relation_options(
        const Order& order,
        const ClassGroupComputeOptions& requested) noexcept {
    ClassGroupRelationOptions options;
    if (order.degree() >= 3) {
        options.ideal_search_radius = 12;
    }
    options.max_candidates = requested.max_candidates;
    options.max_relations = requested.max_relations;
    options.relation_saturation_aux_prime_bound =
            requested.relation_saturation_aux_prime_bound;
    options.relation_saturation_max_appends_per_ell =
            requested.relation_saturation_max_appends_per_ell;
    options.zeta_bf_max_cutoff = requested.zeta_bf_max_cutoff;
    options.requested_certification = requested.requested_certification;
    options.diagnostics = requested.diagnostics;
    return options;
}

bool imaginary_quadratic_class_number(
        flint::FmpzRef out,
        flint::FmpzConstRef discriminant) noexcept;

struct WitnessedClassRelationHnfBasis {
    flint::FmpzMat rows{0, 0};
    flint::FmpzMat relation_coefficients{0, 0};
    std::vector<FactoredElement> witnesses;
};

enum class ClassGroupIdealRelationWitnessStatus {
    not_started,
    success,
    invalid_input,
    arithmetic_failure,
    factorization_failure,
    verification_failure,
    exhausted,
};

enum class ClassGroupIdealRelationWitnessStage {
    none,
    input_validation,
    direct_factorization,
    ideal_reduction,
    reduced_factorization,
    continuation_setup,
    continuation_enumeration,
    continuation_factorization,
    random_product_construction,
    random_product_enumeration,
    exact_verification,
    search_exhaustion,
};

struct ClassGroupIdealRelationWitnessResult {
    ClassGroupIdealRelationWitnessStatus status =
            ClassGroupIdealRelationWitnessStatus::not_started;
    ClassGroupIdealRelationWitnessStage stage =
            ClassGroupIdealRelationWitnessStage::none;
    slong candidates_tried = 0;
    slong random_products_tried = 0;
    bool used_reduction = false;
    bool used_random_product = false;
};

struct ClassGroupIdealRelationWitnessOptions {
    slong reduction_precision = 200;
    slong max_candidates = 500;
    slong max_candidates_per_ideal = 64;
    slong max_random_products = 8;
    ulong random_seed = UWORD(0x6a09e667f3bcc909);
};

const char* class_group_ideal_relation_witness_status_name(
        ClassGroupIdealRelationWitnessStatus status) noexcept;
const char* class_group_ideal_relation_witness_stage_name(
        ClassGroupIdealRelationWitnessStage stage) noexcept;

bool verify_class_group_ideal_relation_witness(
        const ClassGroupContext& context,
        const Ideal& ideal,
        const FactoredElement& multiplier,
        flint::FmpzMatConstRef factor_base_row) noexcept;
bool factor_base_row_ideal(
        FractionalIdeal& out,
        const FactorBase& base,
        flint::FmpzMatConstRef factor_base_row) noexcept;
bool multiply_factored_element_power_fmpz(
        FactoredElement& accumulator,
        const FactoredElement& base,
        flint::FmpzConstRef exponent) noexcept;

bool class_group_ideal_relation_continuation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& original_ideal,
        const Ideal& reduced_ideal,
        const Element& reduction_multiplier,
        const ClassGroupIdealRelationWitnessOptions& options = {}) noexcept;

bool class_group_ideal_relation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& ideal,
        const ClassGroupIdealRelationWitnessOptions& options = {}) noexcept;

struct CompactFieldModulusCacheEntry {
    flint::Fmpz prime;
    flint::FmpzPoly modulus;
    slong threshold = 0;
    bool usable = false;
};

struct CompactFieldModulusCache {
    const NumberField* field = nullptr;
    std::vector<CompactFieldModulusCacheEntry> entries;
};

bool compact_first_reduction(
        Ideal& reduced,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool compact_presentation_reduction(
        Ideal& residual_ideal,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool compact_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool compact_unit_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const Order& order,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics = nullptr,
        CompactFieldModulusCache* field_modulus_cache = nullptr) noexcept;

inline bool valid_certification_request(CertificationMode mode) noexcept {
    return mode == CertificationMode::unknown ||
           mode == CertificationMode::proven;
}

inline bool valid_paired_certification_request(
        CertificationMode mode) noexcept {
    return mode == CertificationMode::grh ||
           mode == CertificationMode::proven;
}

inline bool order_has_quadratic_source_data(const Order& order) noexcept {
    flint::Fmpz radicand;
    return order.parent() != nullptr &&
           order.parent()->backend_kind() == NumberFieldBackendKind::quadratic &&
           order.parent()->quadratic_radicand(flint::FmpzRef(radicand));
}

inline bool order_supports_exact_quadratic_class_certificate(
        const Order& order) noexcept {
    if (order_has_quadratic_source_data(order)) {
        return true;
    }

    flint::Fmpz discriminant;
    return order.parent() != nullptr && order.degree() == 2 &&
           order.discriminant(flint::FmpzRef(discriminant)) &&
           flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) < 0;
}

inline bool prime_to_fractional_ideal(FractionalIdeal& out,
                                      const PrimeIdeal& prime) noexcept {
    if (!same_order_parent(out.parent(), prime.parent())) {
        return false;
    }

    Ideal integral(*out.parent());
    return integral.is_defined() && prime.get_ideal(integral) &&
           out.set_integral(integral);
}

bool unit_candidate_witnesses(std::vector<FactoredElement>& out,
                                    RelationUnitCandidateBatchState& state,
                                    const ClassGroupContext& context,
                                    slong add)
        noexcept;
HnfFinishWorkspace* class_group_finish_workspace(
        ClassGroupContext& context) noexcept;
const HnfFinishWorkspace* class_group_finish_workspace(
        const ClassGroupContext& context) noexcept;
bool class_group_finish_torsion(
        const flint::Fmpz*& torsion_order,
        const Element*& torsion_generator,
        const ClassGroupContext& context) noexcept;
bool class_group_relation_rows(
        flint::FmpzMatRef out,
        const ClassGroupContext& context,
        slong first_relation) noexcept;
bool sync_hnf_finish_workspace(
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept;
bool hnf_finish_workspace_witness_coefficients(
        flint::FmpzMat& out,
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept;
bool hnf_unit_witness_coefficients(flint::FmpzMat& out,
                                        const ClassGroupContext& context)
        noexcept;
bool hnf_unit_witnesses(std::vector<FactoredElement>& out,
                             const ClassGroupContext& context) noexcept;
bool class_relation_witnessed_hnf_basis(
        WitnessedClassRelationHnfBasis& out,
        const ClassGroupContext& context) noexcept;

}  // namespace silex::detail
