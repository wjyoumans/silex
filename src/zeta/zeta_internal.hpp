#pragma once

#include <silex/diagnostics.hpp>
#include <silex/zeta.hpp>

#include <cstddef>
#include <vector>

namespace silex {
class Element;
class FactorBase;
}

namespace silex::detail {

struct ZetaBfResidueDegreeCacheEntry {
    ulong p = 0;
    std::size_t offset = 0;
    std::size_t length = 0;
};

// Caches, per prime `p`, the residue degrees of the primes of an order above
// `p` (used by the Belabas-Friedman per-prime term in `bf_term`).  Entries
// are keyed only by `p`: there is no order (or computation) identity stored
// alongside them.
//
// Decision (2026-09-26): the cache intentionally carries no order key.  Each
// order, and each computation that needs one, must own and pass its own
// `ZetaBfResidueDegreeCache` instance; callers must never reuse one cache
// across two different orders.  The two owning callers already follow this
// by pairing the cache with an explicit order and checking it before reuse:
//
// - `RelationFactorBasePlan`
//   (class_group/relation_factor_base_plan_internal.hpp) stores its `Order`
//   next to the cache, and `relation_factor_base_plan_residue_degrees`
//   returns the cache only when that order matches.
// - `AnalyticClassRegulatorCache` (order_unit/compute_internal.hpp) tracks
//   the order it last cached for and resets the cache whenever a different
//   order is seen (`reset_for_order_if_needed_`).
//
// The analytic-finish and validation paths (class_group/analytic_finish.cpp,
// order_unit/validation.cpp) only forward a cache pointer already owned this
// way; they do not own one themselves.  A caller that needs to cache residue
// degrees for more than one order must construct one
// `ZetaBfResidueDegreeCache` per order.
struct ZetaBfResidueDegreeCache {
    std::vector<ZetaBfResidueDegreeCacheEntry> entries;
    std::vector<slong> residue_degrees;
    std::size_t lookup_hint = 0;
};

// Audit data of the Belabas-Friedman evaluation behind a default-route
// zeta product: the truncation error bound added to the log residue (valid
// under GRH, Belabas-Friedman 2015, Theorem 1), the cutoff used, the
// maximum cutoff the default route allows, and the requested and working
// precisions.  The default route does not require the error bound to meet
// the precision target, so `cutoff` may equal `max_cutoff`.
struct ZetaBfRouteAudit {
    flint::Arb error_bound;
    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
};

// `out` is the GRH generation bound min(BDF, Bach), or Bach's bound alone
// when the BDF criterion cannot be computed or decided (for example for a
// non-monic defining polynomial).  `bach_selected`, when given, is set on
// success to whether Bach's bound gave `out`: strictly below the BDF bound,
// or the BDF bound unavailable; on a tie the BDF bound is reported.
bool grh_factor_base_bound_with_diagnostics(
        flint::FmpzRef out,
        const Order& order,
        const DiagnosticsContext* diagnostics,
        bool* bach_selected = nullptr) noexcept;

// `unconditional`, when given, is set to true exactly when the product came
// from an unconditional route: degree one, or the quadratic L(1, chi) route
// (explicit quadratic backend, maximal order, |D| fitting a ulong).  It is
// false when the product came from the Belabas-Friedman fallback, whose error
// bound assumes GRH for zeta_K and zeta_Q (Belabas-Friedman 2015, Theorem 1).
// `bf_audit`, when given, receives the audit data of that fallback; it is
// written only on success on the fallback route and left unchanged
// otherwise.
bool zeta_class_regulator_product_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        bool* unconditional = nullptr,
        ZetaBfRouteAudit* bf_audit = nullptr) noexcept;

// True when `zeta_class_regulator_product` has an unconditional route for
// `order`: degree one, or the quadratic L(1, chi) route preconditions.  The
// quadratic route can still fall back to Belabas-Friedman if the L-value
// evaluation fails, so a caller that holds the value should use the route
// reported by `zeta_class_regulator_product_with_diagnostics` instead.
bool zeta_unconditional_route_available(const Order& order) noexcept;

bool zeta_class_regulator_product_bf_audit_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_class_regulator_product_validation_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_bf_audit_cutoff_available(const Order& order,
                                    ulong max_cutoff,
                                    slong precision) noexcept;

bool class_regulator_product_estimate_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        flint::Fmpz* torsion_order = nullptr,
        Element* torsion_generator = nullptr) noexcept;

// The character of the quadratic L(1, chi) route: the Kronecker symbol
// (D / n), which for a fundamental discriminant D is the real primitive
// character modulo |D|.
int quadratic_character(flint::FmpzConstRef discriminant, ulong n) noexcept;

// L(1, chi_D) for a fundamental discriminant D with |D| >= 3 fitting a
// ulong, by the approximate functional equation with chi_D evaluated as a
// Kronecker symbol.  Returns false when D is outside that range or the
// enclosure is not finite; D is not checked to be fundamental.
bool quadratic_dirichlet_l1(flint::ArbRef out,
                            flint::FmpzConstRef discriminant,
                            slong precision) noexcept;

}  // namespace silex::detail
