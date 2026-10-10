#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled inverse route for an integral ideal I of a maximal order O
// with p O contained in I and I != O, for a rational prime p (every prime
// ideal of O qualifies).  On success sets numerator to the integral ideal
// p I^{-1} and den to p; numerator / den is then exactly the normalized
// numerator and denominator that FractionalIdeal::invert produces for I.
//
// Returns false, leaving the outputs unchanged, when the order is not marked
// maximal, when the HNF of I does not have the form above, or when an
// internal consistency check fails; callers then use FractionalIdeal::invert.
// The result never depends on which of these inputs reach the fast route.
bool prime_ideal_inverse_num_den(Ideal& numerator,
                                 flint::FmpzRef den,
                                 const Ideal& ideal) noexcept;

}  // namespace silex::detail
