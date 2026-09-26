#pragma once

#include <silex/relation.hpp>

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>

namespace silex::detail {

// Noninstalled relation hooks for class-group relation search and tests.
//
// Every entry point here trusts caller-supplied data: an exact field norm,
// integral coordinates, or a complete exponent row.  None of them re-derives
// that data, so a wrong input installs a wrong relation.  They are therefore
// not part of the installed API; public callers use Relation::set_generator,
// which computes the norm and factorization itself.
class RelationAccess {
public:
    // Installs `generator` with the supplied exponent row unchecked.
    static bool set_relation_from_known_row(
            Relation& out,
            const FactorBase& base,
            const Element& generator,
            flint::FmpzMatConstRef row) noexcept;

    // `norm` must be the exact field norm of `alpha`; hot relation-search
    // loops use this to avoid recomputing a norm they already screened.
    static bool set_generator_with_norm(
            Relation& out,
            const Element& alpha,
            flint::FmpqConstRef norm,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static bool set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static bool set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static bool set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial = nullptr,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static bool factor_relation_row_from_integral_coordinates_and_norm(
            Relation& out,
            bool& handled,
            bool& smooth,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial = nullptr,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static bool commit_relation_generator_from_integral_coordinates(
            Relation& out,
            flint::FmpzMatConstRef integral_coordinates,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    static flint::FmpzMatConstRef pending_relation_exponents_ref(
            const Relation& relation) noexcept;
};

}  // namespace silex::detail
