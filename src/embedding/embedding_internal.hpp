#pragma once

#include <flint/acb.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>

namespace silex {

struct DiagnosticsContext;
class EmbeddingContext;

namespace detail {

enum class RootRefinementStatus {
    success,
    invalid_input,
    not_isolated,
    low_accuracy,
    real_root_validation,
    root_isolation,
    root_pairing,
    precision_cap
};

// Refinement of `previous` with FLINT acb_poly_find_roots, using the previous
// balls as initial approximations, into `roots` in the same place order.
// `roots` must start as a copy of `previous`.  The first working precision is
// `precision` plus a few guard bits; it doubles while an attempt fails, up to
// a fixed multiple of the first one, and then returns
// RootRefinementStatus::precision_cap so the caller can fall back to full
// isolation.  Success means every refined ball is isolated, has relative
// accuracy at least `precision`, passes acb_poly_validate_real_roots, and
// overlaps exactly one previous ball.
RootRefinementStatus refine_roots_from_initial(
        acb_ptr roots,
        acb_srcptr previous,
        const fmpq_poly_t polynomial,
        slong degree,
        slong num_real,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;

// Full root isolation with FLINT arb_fmpz_poly_complex_roots, returned in the
// place order of `previous` (real places first, then each complex place as
// its positive-imaginary root followed by the conjugate).
//
// `previous` must hold pairwise disjoint balls, each containing exactly one
// root of `numerator`, laid out as above with `num_real` real places.  The
// isolation starts at `precision` and doubles until every new ball overlaps
// exactly one previous ball.  That terminates: each root has positive
// distance from the closed previous balls that do not contain it, and the new
// radii shrink below that distance.  Returns false only for invalid input or
// when the working precision would overflow.  `attempts`, when not null,
// receives the number of isolations run.
bool isolate_roots_in_previous_order(acb_ptr roots,
                                     acb_srcptr previous,
                                     const fmpz_poly_t numerator,
                                     slong degree,
                                     slong num_real,
                                     slong precision,
                                     const DiagnosticsContext* diagnostics,
                                     slong* attempts = nullptr) noexcept;

// Noninstalled access to EmbeddingContext internals.
class EmbeddingContextAccess {
public:
    // Gives `out`, a defined context with no roots yet, a copy of the roots
    // stored in `source` for the same field (NumberField::has_same_data),
    // so its first refine() at or below source.precision() needs no root
    // isolation.  The copy keeps `source`'s place order and precision, and
    // `out` is left as `source` would be after the same refine() calls: a
    // request at or below that precision returns the copied balls, and a
    // higher one refines them (refine_roots_from_initial, then full
    // isolation in the same place order).  Nothing else is copied: `out`
    // keeps its own logarithmic-embedding cache, and its evaluation tables
    // are built from the copied roots on first use.
    //
    // Returns false and leaves `out` unchanged when `out` is undefined or
    // already has roots, `source` has none, the fields differ, or the copy
    // cannot be allocated.
    static bool seed_roots(EmbeddingContext& out,
                           const EmbeddingContext& source) noexcept;
};

}  // namespace detail
}  // namespace silex
