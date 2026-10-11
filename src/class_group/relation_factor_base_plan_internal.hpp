#pragma once

#include <silex/diagnostics.hpp>
#include <silex/order.hpp>

#include "../zeta/zeta_internal.hpp"

namespace silex::detail {

struct RelationFactorBasePlan {
    flint::Fmpz working_bound;
    ZetaBfResidueDegreeCache residue_degrees;
    Order order;
    bool valid = false;
    // Set by specialize_maximal_imaginary_quadratic_factor_base_plan when
    // the base up to the generation bound would exceed the relation module
    // size limit: the working bound then stays compact, and the `proven`
    // class group is certified by its structure certificate instead of the
    // generation bound (ClassGroupGenerationBasis::structure_certificate).
    bool imaginary_quadratic_structure_certificate = false;
};

ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        RelationFactorBasePlan& plan,
        const Order& order) noexcept;

const ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        const RelationFactorBasePlan& plan,
        const Order& order) noexcept;

namespace relation_search {

bool build_relation_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool build_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool specialize_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& plan,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool relation_factor_base_restart_limit(
        flint::Fmpz& out,
        const Order& order) noexcept;

bool next_relation_factor_base_bound_to_limit(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept;

bool next_relation_factor_base_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept;

bool quadratic_prime_is_bad(flint::FmpzConstRef discriminant,
                            flint::FmpzConstRef rational_prime) noexcept;

bool quadratic_nth_suitable_ideal_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef discriminant,
        slong count) noexcept;

bool quadratic_working_bound_floor(slong& out,
                                   flint::FmpzConstRef discriminant) noexcept;

// Test hook for the structure-certificate switch of
// specialize_maximal_imaginary_quadratic_factor_base_plan: a positive
// `columns` replaces kRelationModuleMaxColumns as the number of prime ideals
// up to the generation bound past which a `proven` imaginary quadratic plan
// keeps its compact base, so a test can take that route on a small field.
// Zero restores the default.
void set_quadratic_structure_certificate_column_limit_for_testing(
        slong columns) noexcept;

}  // namespace relation_search
}  // namespace silex::detail
