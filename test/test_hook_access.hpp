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
#include "class_group/class_group_storage_internal.hpp"
#include "class_group/maximal_imaginary_quadratic_row_hnf_internal.hpp"
#include "class_group/factor_base_honesty_internal.hpp"
#include "class_group/ideal_t2_enumeration_internal.hpp"
#include "class_group/relation_saturation_internal.hpp"
#include "class_group/relation_search_driver_internal.hpp"
#include "element/element_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "lat/certified_lll_internal.hpp"
#include "class_group/relation_factor_base_plan_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/precision_cap_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "relation/relation_internal.hpp"
#include "sunit/sunit_internal.hpp"
#include "zeta/zeta_internal.hpp"
#include <silex/factor_base.hpp>

#include <vector>

namespace silex::detail {

class TestHookAccess {
public:
    // Arms or disarms the ideal-build failure of the unramified
    // index-divisor split (ideal_from_mod_rows, or pO for an inert p): while
    // armed every such build fails as an allocation failure would.  Skips
    // the invariant that the build only fails when an allocation fails.
    static void set_prime_ideal_build_failure(bool armed) noexcept {
        PrimeIdealTestSeams::ideal_build_failure_.store(
                armed, std::memory_order_relaxed);
    }

    // Arms or disarms the kernel_row_root fallback-error seam: while armed
    // the compact result counts as undecided and the factored fallback as an
    // argument error.  Skips the invariant that the fallback accepts the
    // arguments the compact routine accepted.
    static void set_kernel_row_root_fallback_error(bool armed) noexcept {
        KernelRowRootTestSeams::fallback_error_.store(
                armed, std::memory_order_relaxed);
    }

    // How many ideal builds the armed seam has failed since the process
    // started.
    static slong prime_ideal_build_failures_taken() noexcept {
        return PrimeIdealTestSeams::ideal_build_failures_taken_.load(
                std::memory_order_relaxed);
    }

    // Arms (delta != 0) or disarms (0) a one-shot perturbation: the next
    // IdealFactorization::factor adds delta to the first valuation v_P it
    // computes, so the norm accounting sum f_P v_P differs from v_p(N(I)).
    // Skips the invariant that PrimeIdeal::valuation is exact.  The
    // perturbation persists until a factor() call consumes it.
    static void set_ideal_factorization_valuation_perturbation(
            slong delta) noexcept {
        IdealFactorizationTestSeams::valuation_perturbation_.store(
                delta, std::memory_order_relaxed);
    }

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

    // The residue-screen factorization of integral coordinates with exact
    // field norm `norm` over `base`, forwarded to the private
    // RelationAccess member whose body uses file-local helpers of
    // relation.cpp.  `exponents` must be 1 x base.length().  `handled` is
    // false when the screen leaves the element to the exact fallback;
    // otherwise `smooth` tells whether the row is a factorization.  With
    // `deduce_last_hit_valuation` false every residue-hit valuation is
    // computed, none is deduced from the norm; `deduced_hit_valuations`, if
    // given, receives the number deduced.  Trusts the supplied norm.
    static bool residue_screen_row(flint::FmpzMatRef exponents,
                                   bool& handled,
                                   bool& smooth,
                                   const FactorBase& base,
                                   flint::FmpzMatConstRef integral_coordinates,
                                   flint::FmpqConstRef norm,
                                   bool deduce_last_hit_valuation,
                                   slong* deduced_hit_valuations) noexcept {
        return RelationAccess::residue_screen_row(
                exponents, handled, smooth, base, integral_coordinates, norm,
                deduce_last_hit_valuation, deduced_hit_valuations);
    }

    // The linear-root residue test by the fmpz evaluation and by the
    // production (word-size when p fits a limb) dispatch; see
    // RelationAccess::linear_root_residue_paths.
    static bool linear_root_residue_paths(std::vector<int>& results,
                                          bool& word_path,
                                          flint::FmpzMatConstRef coordinates,
                                          flint::FmpzMatConstRef roots,
                                          flint::FmpzConstRef p) noexcept {
        return RelationAccess::linear_root_residue_paths(
                results, word_path, coordinates, roots, p);
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

    // Replaces the slack the incremental independence test adds to its
    // target precision after an undecided ball (32 bits); `bits` may be
    // negative, so the target step stays undecided.  Skips no invariant:
    // the slack only chooses a precision to try, and the decisions are
    // still the ball tests'.  IncrementalTargetTestSeams::kOff restores the
    // default.
    static void set_incremental_target_slack(slong bits) noexcept {
        IncrementalTargetTestSeams::forced_slack_.store(
                bits, std::memory_order_relaxed);
    }

    // Forces the relation step of the incremental independence test to
    // find no candidate or to check a perturbed one, so it decides nothing
    // and the test goes on with the lower-bound rule alone.  Skips no
    // invariant: the step only adds dependent decisions, each proved by
    // the torsion check.  Mode::off restores the default.
    static void set_incremental_relation_mode(
            IncrementalRelationTestSeams::Mode mode) noexcept {
        IncrementalRelationTestSeams::forced_mode_.store(
                mode, std::memory_order_relaxed);
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
                              prime.f_, prime.integral_ideal(),
                              residue_polynomial, false);
    }

    // Whether `prime` still holds a deferred two-element generator, so its
    // HNF has not been built yet (PrimeIdeal::integral_ideal builds it).
    static bool prime_ideal_hnf_deferred(const PrimeIdeal& prime) noexcept {
        return prime.deferred_ideal_ != nullptr;
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

    // Sets the failure injected into every later double-precision phase of
    // the analytic finish's relation kernel basis (FinishKernelTestSeams,
    // class_group/class_group_internal.hpp) until it is set back to none.
    // The fallback then returns relation_kernel_basis's rows, which the
    // finish also takes whenever FLINT's reducedness check would accept the
    // phase, so no invariant is skipped; it only forces the fallback.
    static void set_finish_kernel_failure(
            FinishKernelTestFailure failure) noexcept {
        FinishKernelTestSeams::failure_.store(failure,
                                              std::memory_order_relaxed);
    }

    // The number of phase runs an injected failure has replaced since the
    // process started.
    static slong finish_kernel_failures_injected() noexcept {
        return FinishKernelTestSeams::injected_.load(
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
    // through the presentation's left transform first, whatever the
    // presentation's size (SaturationKernelTestSeams), until it is set back
    // to false.  It is the route those kernels take at or below the size
    // crossover and whenever the extension route fails, so no invariant is
    // skipped.  Overrides set_saturation_kernel_force_extension.
    static void set_saturation_kernel_force_transform(bool force) noexcept {
        SaturationKernelTestSeams::force_transform_route_.store(
                force, std::memory_order_relaxed);
    }

    // Sends every later saturation kernel basis of more than one row
    // through the extension route first, whatever the presentation's size
    // (SaturationKernelTestSeams), until it is set back to false.  It is
    // the route those kernels take above the size crossover, and every
    // failure of it falls back to the left transform, so no invariant is
    // skipped.
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

    // The number of beta coefficient sets solved from the extension
    // route's image rows since the process started.
    static slong saturation_beta_rows_solved() noexcept {
        return SaturationKernelTestSeams::beta_solved_.load(
                std::memory_order_relaxed);
    }

    // The number of extension-route saturation kernel bases kept from the
    // row-subset attempt since the process started.
    static slong saturation_kernel_subsets_kept() noexcept {
        return SaturationKernelTestSeams::subset_kept_.load(
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

    // The batched, per-prime and nmod linear factor counts of the zeta
    // Belabas-Friedman prime loop (see ZetaBfTestEntry).
    static bool zeta_bf_linear_factor_counts(
            std::vector<ulong>& primes,
            std::vector<slong>& batch_counts,
            std::vector<slong>& per_prime_counts,
            std::vector<slong>& nmod_counts,
            const Order& order,
            ulong first,
            ulong end) noexcept {
        return ZetaBfTestEntry::linear_factor_counts(
                primes, batch_counts, per_prime_counts, nmod_counts, order,
                first, end);
    }

    // The cyclotomic closed forms of those counts against the word kernel
    // and the nmod path (see ZetaBfTestEntry).
    static bool zeta_bf_cyclotomic_linear_factor_counts(
            ulong& cyclotomic_index,
            ulong& real_cyclotomic_index,
            std::size_t& batch_size,
            std::vector<ulong>& primes,
            std::vector<slong>& counts,
            std::vector<slong>& word_counts,
            std::vector<slong>& nmod_counts,
            const Order& order,
            ulong first,
            ulong end) noexcept {
        return ZetaBfTestEntry::cyclotomic_linear_factor_counts(
                cyclotomic_index, real_cyclotomic_index, batch_size, primes,
                counts, word_counts, nmod_counts, order, first, end);
    }

    // The cyclotomic closed forms of the residue degrees against the nmod
    // path (see ZetaBfTestEntry).
    static bool zeta_bf_cyclotomic_residue_degrees(
            ulong& cyclotomic_index,
            ulong& real_cyclotomic_index,
            std::vector<ulong>& primes,
            std::vector<std::vector<slong>>& degrees,
            std::vector<std::vector<slong>>& nmod_degrees,
            const Order& order,
            ulong first,
            ulong end) noexcept {
        return ZetaBfTestEntry::cyclotomic_residue_degrees(
                cyclotomic_index, real_cyclotomic_index, primes, degrees,
                nmod_degrees, order, first, end);
    }

    // The word kernel's two squarings and Frobenius powers (see
    // ZetaBfTestEntry).
    static bool zeta_bf_word_square(std::vector<ulong>& symmetric,
                                    std::vector<ulong>& full,
                                    const std::vector<ulong>& power) noexcept {
        return ZetaBfTestEntry::word_square(symmetric, full, power);
    }

    static bool zeta_bf_word_frobenius(std::vector<ulong>& symmetric,
                                       std::vector<ulong>& full,
                                       slong& symmetric_count,
                                       slong& full_count,
                                       const std::vector<ulong>& f,
                                       ulong p) noexcept {
        return ZetaBfTestEntry::word_frobenius(symmetric, full,
                                               symmetric_count, full_count, f,
                                               p);
    }

    // Zeta test entries (see ZetaTestEntry): the Belabas-Friedman Theorem 7
    // bound at a given cutoff, the GRH generation bound with a given stop
    // bound in place of Bach's, and the stepped incomplete gamma values and
    // blocks of the real-quadratic L(1, chi) evaluation.
    static bool zeta_bf_theorem7_error_bound(flint::ArbRef out,
                                             const Order& order,
                                             ulong cutoff,
                                             slong precision) noexcept {
        return ZetaTestEntry::bf_theorem7_error_bound(out, order, cutoff,
                                                      precision);
    }

    static bool zeta_grh_factor_base_bound_with_stop(
            flint::FmpzRef out,
            const Order& order,
            const flint::Fmpz& stop_bound,
            const DiagnosticsContext* diagnostics,
            bool* bach_selected) noexcept {
        return ZetaTestEntry::grh_factor_base_bound_with_stop(
                out, order, stop_bound, diagnostics, bach_selected);
    }

    static bool quadratic_afe_gamma_stepped(
            flint::ArbVec& out,
            QuadraticAfeGammaStepCount& count,
            ulong shift,
            ulong q,
            ulong first,
            ulong max_block,
            double target,
            slong prec,
            slong guard,
            QuadraticAfeStepFault fault = QuadraticAfeStepFault::none)
            noexcept {
        return ZetaTestEntry::quadratic_afe_gamma_stepped(
                out, count, shift, q, first, max_block, target, prec, guard,
                fault);
    }

    static bool quadratic_afe_gamma_step_block(
            QuadraticAfeGammaStepBlock& block,
            mag_struct* truncation,
            ulong shift,
            ulong q,
            ulong n,
            double target,
            slong prec) noexcept {
        return ZetaTestEntry::quadratic_afe_gamma_step_block(
                block, truncation, shift, q, n, target, prec);
    }

    // Class-group test entries.

    // The proven witness generator of the row ideal of a maximal imaginary
    // quadratic order (see MaximalImaginaryQuadraticRowTestEntry).
    static bool maximal_imaginary_quadratic_row_generator(
            Element& out,
            const Order& order,
            const FactorBase& base,
            flint::FmpzMatConstRef row,
            bool use_row_hnf_builder) noexcept {
        return MaximalImaginaryQuadraticRowTestEntry::row_generator(
                out, order, base, row, use_row_hnf_builder);
    }

    // The analytic finish embedding context kept in route_state, refined to
    // precision (see AnalyticFinishTestEntry).
    static EmbeddingContext* analytic_finish_embeddings(
            relation_search::RelationSearchControlState& route_state,
            const NumberField& field,
            slong precision) noexcept {
        return AnalyticFinishTestEntry::embeddings(
                route_state, field, precision);
    }

    // The direct large-prime partial relation factorization (see
    // DirectPartialRelationTestEntry).
    static bool direct_partial_relation_row(
            flint::FmpzMat& out_row,
            flint::Fmpz& out_rational_large_prime,
            const Element& generator,
            const FactorBase& base,
            flint::FmpzConstRef factor_base_bound,
            const fmpq* known_norm = nullptr,
            const fmpz_mat_struct* known_integral_coordinates = nullptr,
            bool word_size_residues = true,
            bool deduce_last_hit_valuation = true,
            slong* deduced_hit_valuations = nullptr) noexcept {
        return DirectPartialRelationTestEntry::row(
                out_row, out_rational_large_prime, generator, base,
                factor_base_bound, known_norm, known_integral_coordinates,
                word_size_residues, deduce_last_hit_valuation,
                deduced_hit_valuations);
    }

    // Factor-base honesty scan call counters (the test-only counters in
    // relation_search::FactorBaseHonestyCallCounts): the scans started per
    // route and escalation mode since the last reset.
    static void factor_base_honesty_call_counts_reset() noexcept {
        relation_search::FactorBaseHonestyCallCounts::reset();
    }
    static slong factor_base_honesty_call_count(
            bool direct_required_prime_witness,
            relation_search::FactorBaseWitnessEscalation escalation) noexcept {
        return relation_search::FactorBaseHonestyCallCounts::count(
                direct_required_prime_witness, escalation);
    }

    // Test-only T2 element-step cap of the factor-base honesty search (see
    // relation_search::FactorBaseHonestyStepCapSeams): a positive `steps`
    // replaces the stage-0 cap, zero restores the default; the hit counter
    // counts enumerations that stopped at the cap.
    static void set_factor_base_honesty_element_step_cap(
            slong steps) noexcept {
        relation_search::FactorBaseHonestyStepCapSeams::limit_.store(
                steps > 0 ? steps : 0, std::memory_order_relaxed);
    }
    static void factor_base_honesty_element_step_cap_hits_reset() noexcept {
        relation_search::FactorBaseHonestyStepCapSeams::hits_.store(
                0, std::memory_order_relaxed);
    }
    static slong factor_base_honesty_element_step_cap_hits() noexcept {
        return relation_search::FactorBaseHonestyStepCapSeams::hits_.load(
                std::memory_order_relaxed);
    }

    // Test-only forced misses of the factor-base honesty witness search (see
    // relation_search::FactorBaseHonestyForcedMissSeams): a positive
    // `stages` discards the witness found at each search stage below it (1
    // forces a stage-0 miss), zero restores the default; the miss counter
    // counts the discarded witnesses.
    static void set_factor_base_honesty_forced_miss_stages(
            slong stages) noexcept {
        relation_search::FactorBaseHonestyForcedMissSeams::stages_.store(
                stages > 0 ? stages : 0, std::memory_order_relaxed);
    }
    static void factor_base_honesty_forced_misses_reset() noexcept {
        relation_search::FactorBaseHonestyForcedMissSeams::misses_.store(
                0, std::memory_order_relaxed);
    }
    static slong factor_base_honesty_forced_misses() noexcept {
        return relation_search::FactorBaseHonestyForcedMissSeams::misses_.load(
                std::memory_order_relaxed);
    }

    // Test-only mode of the T2 double-row setup (see
    // relation_search::T2DoubleRowsTestSeams); `off` restores production.
    static void set_t2_double_rows_mode(
            relation_search::T2DoubleRowsTestSeams::Mode mode) noexcept {
        relation_search::T2DoubleRowsTestSeams::mode_.store(
                static_cast<int>(mode), std::memory_order_relaxed);
    }

    // Counts of analytic-bound searches that stopped at the cap (see
    // AnalyticBoundCapSeams).
    static void analytic_bound_cap_counts_reset() noexcept {
        for (auto& c : AnalyticBoundCapSeams::counts_) {
            c.store(0, std::memory_order_relaxed);
        }
    }
    static slong analytic_bound_cap_count(bool quadratic) noexcept {
        return AnalyticBoundCapSeams::counts_[quadratic ? 1 : 0].load(
                std::memory_order_relaxed);
    }

    // Relation-kernel log-matrix word rescaling (see
    // RelationKernelRescaleTestEntry).
    static bool relation_kernel_rescale_log_matrix(
            flint::FmpzMatRef rows,
            flint::ArbMatConstRef matrix,
            slong precision,
            slong word_bits) noexcept {
        return RelationKernelRescaleTestEntry::rescale_log_matrix(
                rows, matrix, precision, word_bits);
    }

    // Regulator pivot selection (see RegulatorPivotTestEntry).
    static RegulatorPivotOutcome regulator_pivot_rows(
            std::vector<slong>& out,
            flint::ArbMatConstRef matrix,
            slong precision) noexcept {
        return RegulatorPivotTestEntry::rows(out, matrix, precision);
    }
    static RegulatorPivotOutcome regulator_pivot_unit_indices(
            std::vector<slong>& out,
            std::vector<slong>& pivot_rows,
            const flint::ArbMat& unit_logs,
            const Signature& sig,
            slong rank,
            slong places,
            slong precision) noexcept {
        return RegulatorPivotTestEntry::unit_indices(
                out, pivot_rows, unit_logs, sig, rank, places, precision);
    }

    // Regulator pivot finish action (see RegulatorPivotFinishTestEntry).
    static RegulatorPivotFinishTestResult
    regulator_pivot_finish(
            RegulatorPivotFinishPath path,
            RegulatorPivotOutcome outcome,
            slong independent_count,
            slong target_rank,
            bool factor_base_restart_available,
            slong analytic_finish_precision,
            slong analytic_precision_doublings) noexcept {
        return relation_search::RegulatorPivotFinishTestEntry::finish(
                path, outcome, independent_count, target_rank,
                factor_base_restart_available, analytic_finish_precision,
                analytic_precision_doublings);
    }

    // Bounded best-denominator rational reconstruction (see
    // FmpqBestDenominatorTestEntry).
    static bool fmpq_best_denominator_bounded(
            flint::Fmpq& out,
            flint::FmpqConstRef value,
            flint::FmpzConstRef max_denominator) noexcept {
        return FmpqBestDenominatorTestEntry::bounded(out, value,
                                                     max_denominator);
    }

    // Word-size quotient logs of the degree-one residue column (see
    // WordQuotientLogsTestEntry).
    static bool word_quotient_logs(std::vector<ulong>& logs,
                                   std::vector<bool>& resolved,
                                   const std::vector<ulong>& images,
                                   ulong zeta,
                                   ulong p,
                                   ulong ell) noexcept {
        return WordQuotientLogsTestEntry::logs(logs, resolved, images, zeta,
                                               p, ell);
    }

    // Marks the finish workspace's exact state invalid while leaving its
    // kept finish extension in place, so the next sync rebuilds the exact
    // state with an extension still kept.  Skips the invariant that only
    // clear() drops exact_valid, and it drops the extension too.
    static void invalidate_finish_exact_state(
            HnfFinishWorkspace& workspace) noexcept {
        workspace.exact_valid = false;
    }

    // Adds one to the first entry of each of the first `rows` rows of the
    // workspace's kept finish extension, so its relation part no longer
    // matches the stored relations.  Skips the invariant that only the
    // finish phase writes the extension.
    static void perturb_finish_extension(HnfFinishWorkspace& workspace,
                                         slong rows) noexcept {
        for (slong i = 0; i < rows; ++i) {
            fmpz* entry = flint::fmpz_mat_entry(
                                  flint::FmpzMatRef(workspace.finish_extension),
                                  i, 0)
                                  .raw();
            ::fmpz_add_ui(entry, entry, 1);
        }
    }

    // Sets the context's relation rank to `rank` without syncing its row
    // module, as a stale or modular rank would be.  Skips the invariant
    // that the rank is the row module's rank.
    static void set_context_relation_rank(ClassGroupContext& context,
                                          slong rank) noexcept {
        context.relation_rank_ = rank;
    }
};

}  // namespace silex::detail
