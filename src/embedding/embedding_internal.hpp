#pragma once

#include <flint/acb.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>

namespace silex {

struct DiagnosticsContext;

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

}  // namespace detail
}  // namespace silex
