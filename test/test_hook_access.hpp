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

#include "class_group/class_group_internal.hpp"
#include "class_group/relation_saturation_internal.hpp"
#include "element/element_internal.hpp"
#include "lat/certified_lll_internal.hpp"
#include "class_group/relation_factor_base_plan_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/precision_cap_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "sunit/sunit_internal.hpp"
#include "zeta/zeta_internal.hpp"
#include <silex/factor_base.hpp>

namespace silex::detail {

class TestHookAccess {
public:
    // Arms (n > 0) or disarms (0) the S-unit composition overflow: the n-th
    // composition call of the next builds (counting S-class invariant
    // witnesses first, then S-unit generators) is replaced by a power whose
    // exponent product overflows slong, so it fails.  Skips the invariant
    // that relation coefficients keep the exponent products in range.  The
    // countdown persists across builds until it fires; set 0 to disarm.
    static void set_sunit_composition_overflow_at(slong n) noexcept {
        SUnitTestSeams::composition_overflow_countdown_.store(
                n, std::memory_order_relaxed);
    }

    // Prime-ideal test entries, forwarded to private PrimeIdealAccess
    // members whose bodies use file-local helpers of prime_ideal.cpp.

    // Values `element` at `prime` by coordinate divisibility with the
    // anti-uniformizer from the HNF of P (the route for primes without a
    // Kummer quotient map), for every prime of a maximal order, so tests can
    // compare it with the power-containment fallback.
    static bool anti_uniformizer_valuation(slong& out,
                                           const PrimeIdeal& prime,
                                           const OrderElement& element,
                                           slong bound) noexcept {
        return PrimeIdealAccess::anti_uniformizer_valuation(out, prime,
                                                            element, bound);
    }

    // Splits O/(p, I_p) for a prime p dividing the index of the maximal
    // order, by bounded enumeration (`recursive_split` false) or the
    // recursive split (true); skips the separating-element search.
    static bool index_divisor_prime_ideals(std::vector<Ideal>& out,
                                           const Order& order,
                                           flint::FmpzConstRef p,
                                           bool recursive_split) noexcept {
        return PrimeIdealAccess::index_divisor_prime_ideals(
                out, order, p, recursive_split);
    }

    // The route decompose_prime(.., order, p, max_residue_degree) takes.
    static DecomposePrimeRoute decompose_prime_route(
            const Order& order,
            flint::FmpzConstRef p,
            slong max_residue_degree) noexcept {
        return PrimeIdealAccess::decompose_prime_route_for_max_degree(
                order, p, max_residue_degree);
    }

    // Runs the exact power-containment valuation fallback (the largest
    // i <= `bound` with `element` in P^i) directly, bypassing the norm,
    // Hensel and coordinate-divisibility routes.  Skips the routing that
    // production valuation does; the result is still the exact valuation.
    static bool valuation_by_power_containment(slong& out,
                                               const PrimeIdeal& prime,
                                               const OrderElement& element,
                                               slong bound) noexcept {
        return prime.valuation_by_power_containment(out, element, bound,
                                                    nullptr);
    }

    // A positive `cap` replaces the cap computed by precision_doubling_cap;
    // zero restores the default.
    static void set_precision_doubling_cap(slong cap) noexcept {
        PrecisionCapTestSeams::forced_cap_.store(
                cap > 0 ? cap : 0, std::memory_order_relaxed);
    }

    // A positive `max_aux_primes` replaces the per-scan guard of 2^16
    // examined primes q of the stable residue-character proof scan, and a
    // nonzero `max_q` replaces the machine-word upper end of the q range, so
    // a test can reach the guard or the end of the range quickly.  Zero
    // restores each default.
    static void set_stable_proof_scan_limits(slong max_aux_primes,
                                             ulong max_q) noexcept {
        StableProofTestSeams::max_aux_primes_.store(
                max_aux_primes > 0 ? max_aux_primes : 0,
                std::memory_order_relaxed);
        StableProofTestSeams::max_q_.store(max_q, std::memory_order_relaxed);
    }

    // A positive `columns` replaces kRelationModuleMaxColumns as the number
    // of prime ideals up to the generation bound past which a `proven`
    // imaginary quadratic plan keeps its compact base, so a test can take
    // the structure-certificate route on a small field.  Zero restores the
    // default.
    static void set_quadratic_structure_certificate_column_limit(
            slong columns) noexcept {
        QuadraticStructureCertificateTestSeams::column_limit_.store(
                columns > 0 ? columns : 0, std::memory_order_relaxed);
    }

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

    // Sets the failure injected into every later run of FLINT's certified
    // first LLL phase (CertifiedLllTestSeams, lat/certified_lll_internal.hpp),
    // at every call site, until it is set back to none.  The caller's
    // fallback then returns fmpz_lll's result, which every other route also
    // returns, so no invariant is skipped; it only forces the fallback where
    // the first phase would be certified.
    static void set_certified_lll_failure(
            CertifiedLllTestFailure failure) noexcept {
        CertifiedLllTestSeams::failure_.store(failure,
                                              std::memory_order_relaxed);
    }

    // The number of first-phase runs an injected failure has replaced since
    // the process started.
    static slong certified_lll_failures_injected() noexcept {
        return CertifiedLllTestSeams::injected_.load(
                std::memory_order_relaxed);
    }

    // Sets the failure injected into every later extension-route run of the
    // relation-saturation kernel basis
    // (SaturationKernelTestSeams) until it is set back to none.  The
    // fallback then builds the basis from fmpz_lll, which every other route
    // can also take, so no invariant is skipped; it only forces the
    // fallback where the route would be certified.
    static void set_saturation_kernel_failure(
            SaturationKernelTestFailure failure) noexcept {
        SaturationKernelTestSeams::failure_.store(failure,
                                                  std::memory_order_relaxed);
    }

    // Sends every later saturation kernel basis of more than one row
    // through the extension route (SaturationKernelTestSeams) until it is
    // set back to false.  That route is certified like the first phase and
    // falls back to fmpz_lll in the same way, so no invariant is skipped.
    static void set_saturation_kernel_force_extension(bool force) noexcept {
        SaturationKernelTestSeams::force_extension_route_.store(
                force, std::memory_order_relaxed);
    }

    // The number of extension-route runs an injected failure
    // has replaced since the process started.
    static slong saturation_kernel_failures_injected() noexcept {
        return SaturationKernelTestSeams::injected_.load(
                std::memory_order_relaxed);
    }

    // Arms (or disarms) one use: the next
    // compact_unit_presentation_power_root call returns
    // CompactPowerRootResult::error after its argument checks, as an
    // internal error would.  No valid call reaches that error: both
    // callers check the arguments first, and OrderUnitGroup::define
    // refuses an order without a basis, so the working ideal is defined.
    // So this skips no invariant a run keeps; it only lets a test reach
    // the callers' error handling.
    static void set_compact_unit_power_root_error(bool on) noexcept {
        CompactPowerRootTestSeams::unit_power_root_error_.store(
                on, std::memory_order_relaxed);
    }

    // Make the next call of the named certification route fail once, after
    // the route has written its certification state
    // (CertificationFailureTestSeams, class_group_internal.hpp).  Each
    // switch consumes itself.  A failure is a result every route already
    // returns; this only lets a test reach it with a valid fixture.
    static void arm_exact_order_publish_failure() noexcept {
        CertificationFailureTestSeams::exact_order_publish_failure_.store(
                true, std::memory_order_relaxed);
    }
    static void arm_exact_order_mismatch() noexcept {
        CertificationFailureTestSeams::exact_order_mismatch_.store(
                true, std::memory_order_relaxed);
    }
    static void arm_zeta_bf_record_failure() noexcept {
        CertificationFailureTestSeams::zeta_bf_record_failure_.store(
                true, std::memory_order_relaxed);
    }
    // True while any of the switches above is still armed, that is, the
    // armed route has not reached its seam.
    static bool certification_failure_pending() noexcept {
        return CertificationFailureTestSeams::exact_order_publish_failure_
                       .load(std::memory_order_relaxed) ||
               CertificationFailureTestSeams::exact_order_mismatch_.load(
                       std::memory_order_relaxed) ||
               CertificationFailureTestSeams::zeta_bf_record_failure_.load(
                       std::memory_order_relaxed);
    }
    static void disarm_certification_failures() noexcept {
        CertificationFailureTestSeams::exact_order_publish_failure_.store(
                false, std::memory_order_relaxed);
        CertificationFailureTestSeams::exact_order_mismatch_.store(
                false, std::memory_order_relaxed);
        CertificationFailureTestSeams::zeta_bf_record_failure_.store(
                false, std::memory_order_relaxed);
    }

    // The number of extension-route saturation kernel bases kept since the
    // process started.
    static slong saturation_kernel_extensions_kept() noexcept {
        return SaturationKernelTestSeams::extension_kept_.load(
                std::memory_order_relaxed);
    }

    // While set, pure_power_hensel_root returns "unsupported" in degree 10
    // and higher, as before the lift ran there, so Element::is_power leaves
    // a non-rational true power (exponent >= 2) undecided in such a field.
    // Skips no invariant: "unsupported" is a result every caller already
    // handles (a failed lift returns it too); it only lets a test reach the
    // paths after an undecided root with a small fixture.
    static void set_high_degree_power_lift_off(bool on) noexcept {
        ElementPowerTestSeams::high_degree_lift_off_.store(
                on, std::memory_order_relaxed);
    }

    // The number of power tests whose pure-power lift a SaturationLiftLimit
    // scope skipped since the process started.
    static slong saturation_lifts_skipped() noexcept {
        return SaturationLiftLimit::skipped_lifts_.load(
                std::memory_order_relaxed);
    }

    // FactorBase::build_lll_relation_base_impl(bound, true): the
    // incomplete-base LLL relation base.
    static bool build_lll_relation_base_incomplete(
            FactorBase& base, flint::FmpzConstRef bound) noexcept {
        return base.build_lll_relation_base_impl(bound, true);
    }

    // FactorBase::append_prime, which the builders call.  Fails, leaving
    // base unchanged, if base already holds `prime`.
    static bool append_prime(FactorBase& base,
                             const PrimeIdeal& prime) noexcept {
        return base.append_prime(prime);
    }

    // The maximal imaginary-quadratic builder that ClassGroupContext reaches
    // in an oriented quadratic class-unit transaction.
    static bool build_maximal_imaginary_quadratic(
            FactorBase& base, flint::FmpzConstRef bound) noexcept {
        return base.build_maximal_imaginary_quadratic_relation_base(bound,
                                                                    nullptr);
    }

    // The batched and per-prime linear factor counts of the zeta
    // Belabas-Friedman prime loop (see ZetaBfTestEntry).
    static bool zeta_bf_linear_factor_counts(
            std::vector<ulong>& primes,
            std::vector<slong>& batch_counts,
            std::vector<slong>& per_prime_counts,
            const Order& order,
            ulong first,
            ulong end) noexcept {
        return ZetaBfTestEntry::linear_factor_counts(
                primes, batch_counts, per_prime_counts, order, first, end);
    }
};

}  // namespace silex::detail
