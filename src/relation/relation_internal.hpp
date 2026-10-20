#pragma once

#include <silex/relation.hpp>

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include <vector>

namespace silex::detail {

class TestHookAccess;

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

private:
    friend class TestHookAccess;

    // Test hook: the residue-screen factorization of integral coordinates
    // with exact field norm `norm` over `base`.  `exponents` must be
    // 1 x base.length().  `handled` is false when the screen leaves the
    // element to the exact fallback; otherwise `smooth` tells whether the
    // row in `exponents` is a factorization.  With
    // `deduce_last_hit_valuation` false every residue-hit valuation is
    // computed, none is deduced from the norm; when
    // `deduced_hit_valuations` is given it receives the number of
    // valuations deduced.  For tests only.
    static bool residue_screen_row(
            flint::FmpzMatRef exponents,
            bool& handled,
            bool& smooth,
            const FactorBase& base,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            bool deduce_last_hit_valuation,
            slong* deduced_hit_valuations) noexcept;

    // Test hook: the linear-root residue test of the screen.  The 1 x n
    // `coordinates` are read as the polynomial sum c_i x^i and tested at
    // each entry of the 1 x k `roots` modulo p, both as coordinates and as
    // a polynomial, and both by the fmpz evaluation and by the production
    // dispatch, which takes the word-size path when 2 <= p < 2^64.  One
    // block scratch serves every root, as in the screen.  `results`
    // receives, per root, {coordinate fmpz, coordinate dispatch, polynomial
    // fmpz, polynomial dispatch} as 0/1; `word_path` tells whether p took
    // the word-size path.  For tests only.
    static bool linear_root_residue_paths(
            std::vector<int>& results,
            bool& word_path,
            flint::FmpzMatConstRef coordinates,
            flint::FmpzMatConstRef roots,
            flint::FmpzConstRef p) noexcept;
};

}  // namespace silex::detail
