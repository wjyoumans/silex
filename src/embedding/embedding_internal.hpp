#pragma once

#include <flint/acb.h>
#include <flint/fmpz_poly.h>

namespace silex {

struct DiagnosticsContext;

namespace detail {

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
