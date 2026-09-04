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

}  // namespace relation_search
}  // namespace silex::detail
