#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled access to Ideal internals.  set_hnf_direct is for tests only;
// set_known_ideal_hnf is for library code.
class IdealAccess {
public:
    enum class HnfInput { unnormalized, normalized };
    enum class IdealTest { always, assert_only };

    // The one routine behind Ideal::set_hnf, set_known_ideal_hnf and the
    // private Ideal::set_known_ideal_rows.  `unnormalized` input is put in
    // normalized HNF (failing if singular or misshapen); `normalized` input
    // is only shape-checked and copied.  IdealTest::always rejects a lattice
    // that is not an ideal; assert_only only asserts it (NDEBUG removes it).
    // On failure `ideal` is unchanged.
    static bool store_hnf(Ideal& ideal, flint::FmpzMatConstRef hnf,
                          HnfInput input, IdealTest test) noexcept;

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


    // Implementations of detail::set_known_two_generator_ideal and
    // detail::set_integral_ideal_known_hnf, defined in ideal.cpp.
    static bool set_known_two_generator(
            Ideal& out, flint::FmpzConstRef scalar_generator,
            const OrderElement& element_generator) noexcept;
    static bool set_integral_known_hnf(FractionalIdeal& out,
                                       const Ideal& numerator) noexcept;

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

// Trust helpers, noninstalled.  set_known_two_generator_ideal stores the HNF
// of a O + beta O for a positive scalar a and an element beta of the same
// order, trusting that the lattice is an ideal.  set_integral_ideal_known_hnf
// sets `out` to the integral fractional ideal with the given numerator.
// Both return false, leaving `out` unchanged, on undefined or mismatched
// input.
inline bool set_known_two_generator_ideal(
        Ideal& out, flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator) noexcept {
    return IdealAccess::set_known_two_generator(out, scalar_generator,
                                                element_generator);
}

inline bool set_integral_ideal_known_hnf(FractionalIdeal& out,
                                         const Ideal& numerator) noexcept {
    return IdealAccess::set_integral_known_hnf(out, numerator);
}

}  // namespace silex::detail
