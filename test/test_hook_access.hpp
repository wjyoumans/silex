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

#include "order_unit/compute_internal.hpp"

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

    // While set, compute_with_class_group clears its validation summary after
    // each failed candidate attempt, so a failure with compute_ok true looks
    // as if validation never ran (outcome not_run).  Skips the invariant that
    // the summary names the validation that ended the attempt; no real run
    // reaches not_run with compute_ok true.
    static void set_forget_validation_summary(bool on) noexcept {
        ComputeTestSeams::forget_validation_summary_.store(
                on, std::memory_order_relaxed);
    }

    // Arms (or disarms) one use: the next grh continuation in
    // compute_with_relation_class_group_ does not call
    // extend_tentative_relations_ and acts as if it returned true without
    // appending a relation.  Later continuations run normally.  Skips the invariant that a
    // successful extension grows relation_count().
    static void set_continuation_appends_nothing(bool on) noexcept {
        ComputeTestSeams::continuation_appends_nothing_.store(
                on, std::memory_order_relaxed);
    }

    // Calls of word_quotient_log_resolve that took the shared walk and the
    // baby-step giant-step route since the process started.  Read-only; a
    // test compares the counts before and after a call.
    static unsigned long word_quotient_log_walk_calls() noexcept {
        return ComputeTestSeams::word_quotient_log_walk_calls_.load(
                std::memory_order_relaxed);
    }
    static unsigned long word_quotient_log_baby_giant_calls() noexcept {
        return ComputeTestSeams::word_quotient_log_baby_giant_calls_.load(
                std::memory_order_relaxed);
    }
};

}  // namespace silex::detail
