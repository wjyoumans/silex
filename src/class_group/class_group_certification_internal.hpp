#pragma once

#include <silex/class_group.hpp>

namespace silex::detail {

class ClassGroupCertificationAccess {
public:
    static bool exact_imaginary_quadratic_class_order_for_run(
            flint::FmpzRef out,
            ClassGroupContext& context,
            flint::FmpzConstRef discriminant) noexcept;

    static bool try_certify_imaginary_quadratic_from_exact_order(
            ClassGroupContext& context,
            CertificationMode requested,
            flint::FmpzConstRef discriminant,
            flint::FmpzConstRef exact_order) noexcept;

    static bool record_factor_base_honesty_proof(
            ClassGroupContext& context,
            flint::FmpzConstRef required_bound) noexcept;

    static bool rank_zero_quadratic_class_index_bound(
            flint::FmpzRef out,
            ClassGroupContext& context,
            OrderUnitGroup& units,
            CertificationMode requested) noexcept;

    static bool saturate_relations_for_index_bound_with_units(
            bool& changed,
            bool& saturated,
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;

    static bool try_certify_class_unit_with_bf_audit(
            ClassGroupContext& context,
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            flint::ArbConstRef error_bound,
            ulong cutoff,
            ulong max_cutoff,
            slong requested_precision,
            slong work_precision) noexcept;
};

}  // namespace silex::detail
