#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled access to Ideal internals, for tests only.
class IdealAccess {
public:
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
