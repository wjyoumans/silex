#pragma once

#include "class_group_storage_internal.hpp"
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>

namespace silex::detail::relation_search {

struct OrderCoordinateElementConversion {
    const Order* order = nullptr;
    flint::FmpqMat basis{0, 0};
    flint::FmpqMat power_row{0, 0};
    flint::FmpqPoly polynomial;

    bool reset(const Order& next_order) noexcept;
    bool set(Element& out,
             flint::FmpzMatConstRef coordinates,
             flint::FmpzConstRef den) noexcept;
};

struct NormPrefilter {
    flint::Fmpz rational_prime_product;
    flint::Fmpz large_prime_norm_bound;
    flint::Fmpq norm_scratch;
    flint::Fmpz norm_resultant_scratch;
    flint::Fmpz norm_remainder;
    flint::Fmpz norm_gcd;
    flint::FmpzPoly field_polynomial;
    flint::FmpzPoly coordinate_polynomial;
    slong rejected_count = 0;
    bool valid = false;
    bool allow_large_prime = false;
    bool has_equation_order_norm_polynomial = false;
};

void coordinates_from_lattice_combination(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef coefficients,
        flint::FmpzMatConstRef basis) noexcept;

bool build_norm_prefilter(NormPrefilter& out,
                          const FactorBase& base,
                          flint::FmpzConstRef factor_base_bound,
                          bool allow_large_prime) noexcept;

bool reduced_basis_first_row_is_scalar_rational(
        bool& scalar,
        const Order& order,
        flint::FmpzMatConstRef basis) noexcept;

bool element_is_scalar_rational(bool& scalar,
                                const Element& element) noexcept;

bool try_coordinate_candidate_den(
        ClassGroupContext& context,
        const Order& order,
        flint::FmpzMat& coordinates,
        flint::FmpzConstRef den,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        bool& partial_throttle_exit,
        detail::RelationAdmissionCache* admission_cache = nullptr,
        bool random_relation = false,
        slong* factor_attempts = nullptr,
        slong max_factor_attempts = 0,
        bool* factor_attempt_limit_reached = nullptr,
        Relation* admission_scratch = nullptr,
        OrderCoordinateElementConversion* coordinate_conversion = nullptr,
        bool admission_scratch_base_verified = false) noexcept;

bool try_coordinate_candidate(
        ClassGroupContext& context,
        const Order& order,
        flint::FmpzMat& coordinates,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        bool& partial_throttle_exit,
        OrderCoordinateElementConversion* coordinate_conversion = nullptr)
        noexcept;

}  // namespace silex::detail::relation_search
