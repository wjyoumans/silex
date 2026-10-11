#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled access to Ideal internals.  set_hnf_direct is for tests only;
// set_known_ideal_hnf is for library code.
class IdealAccess {
public:
    // Stores the normalized HNF of the full-rank lattice spanned by the rows
    // of the n x n matrix `hnf`, as Ideal::set_hnf does, but without testing
    // that the lattice is an ideal of the parent order; builds with assertions
    // enabled (NDEBUG undefined) still assert that test.  Fails, leaving
    // `ideal` unchanged, where set_hnf fails for a reason other than that
    // test: an undefined `ideal`, a misshapen or singular matrix.  Only for a
    // lattice that is an ideal by construction, with the argument stated at
    // the call.
    static bool set_known_ideal_hnf(Ideal& ideal,
                                    flint::FmpzMatConstRef hnf) noexcept;


    // Stores `hnf` as the ideal's basis matrix without normalizing it or
    // checking that its row span is an ideal; only the n x n shape is
    // checked.  Every public construction route stores a full-rank upper
    // triangular HNF.  This hook exists so tests can build the other stored
    // shapes that reach the defensive lattice fallbacks of Ideal::contains.
    // Library code must not use it.
    static bool set_hnf_direct(Ideal& ideal,
                               flint::FmpzMatConstRef hnf) noexcept {
        return ideal.set_hnf_direct(hnf);
    }
};

}  // namespace silex::detail
