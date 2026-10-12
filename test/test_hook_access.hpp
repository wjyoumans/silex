#pragma once

// Test-only white-box hooks.
//
// detail::TestHookAccess is declared (incomplete) and befriended by the
// installed classes whose private state a test must force.  It is defined
// only here, in the test tree, with inline definitions, so libsilex contains
// none of these hooks and no production translation unit can call them:
// outside a test the class is incomplete.  tools/check-cpp-profile.py rejects
// any use of the class under src/ or include/ other than the forward
// declaration and the friend declaration.
//
// A hook here may break an invariant that every construction route keeps.
// Use it only to reach a state that tests cannot reach through the public or
// detail::*Access entry points, and say which invariant it skips.

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>

namespace silex::detail {

class TestHookAccess {
public:
    // Labels a set unit group proven without any check.  Tests use it to
    // build a proven-labelled group whose later rejection must come from the
    // class-group gate itself.  Does nothing for an unset group.
    static void mark_certification_proven(OrderUnitGroup& units) noexcept {
        units.mark_certification_proven_();
    }

    // Replaces the residue data of `prime` (which must already have prime
    // data) with `residue_polynomial` in the alpha convention, keeping its p,
    // e, f and ideal, without checking that p does not divide
    // [O : Z[alpha]].  No construction route records alpha residue data at
    // such a p; this lets tests reach the p | index refusal of
    // PrimeIdeal::reduce and residue_variable_numerator.  It does not check
    // that the degree of `residue_polynomial` matches the residue degree f of
    // the prime.  Returns false when `prime` has no prime data.
    static bool set_alpha_residue_polynomial(
            PrimeIdeal& prime,
            flint::FmpzPolyConstRef residue_polynomial) noexcept {
        if (!prime.has_prime_data()) {
            return false;
        }
        return prime.set_data(flint::FmpzConstRef(prime.p_), prime.e_,
                              prime.f_, prime.ideal_, residue_polynomial,
                              false);
    }
};

}  // namespace silex::detail
