#pragma once

#include <atomic>

#include "factor_base_honesty_failure_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"

namespace silex::detail::relation_search {

// Noninstalled predicate selection, independent of witness search strategy.
enum class FactorBaseWitnessPredicate {
    selected,
    full_factorization,
    order_element_direct,
};

struct FactorBaseWitnessAudit {
    bool used_reference = false;
    bool used_scalar_direct = false;
};

bool factor_base_scalar_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        flint::FmpzConstRef scalar,
        FactorBaseWitnessPredicate predicate = FactorBaseWitnessPredicate::selected,
        const DiagnosticsContext* diagnostics = nullptr,
        FactorBaseWitnessAudit* audit = nullptr) noexcept;

bool factor_base_principal_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const OrderElement& generator,
        FactorBaseWitnessPredicate predicate = FactorBaseWitnessPredicate::selected,
        const DiagnosticsContext* diagnostics = nullptr,
        FactorBaseWitnessAudit* audit = nullptr) noexcept;

// Whether the honesty witness search for a required prime may escalate past
// the stage-0 caps.  The reference buch2.c:be_honest gives up after
// 1 + maxtry_HONEST attempts and bnfinit enlarges the factor base, so a
// caller with a factor-base restart or another route to fall back on keeps
// that flow (none).  A caller for which an unwitnessed prime is a final
// failure asks for bounded escalation instead; callers select the mode with
// factor_base_witness_escalation().
enum class FactorBaseWitnessEscalation {
    none,
    bounded,
};

constexpr FactorBaseWitnessEscalation factor_base_witness_escalation(
        bool recovery_available) noexcept {
    return recovery_available ? FactorBaseWitnessEscalation::none
                              : FactorBaseWitnessEscalation::bounded;
}

// Counts the honesty scans started, by search route and escalation mode, so a
// route-level test can check which mode each production caller passed.  The
// routes differ only in the arguments they pass and a field's honest scans
// pass at stage 0, so the mode is not otherwise observable.  One relaxed
// increment per scan; it feeds no decision.
class FactorBaseHonestyCallCounts {
public:
    static void record(bool direct_required_prime_witness,
                       FactorBaseWitnessEscalation escalation) noexcept {
        counts_[direct_required_prime_witness ? 1 : 0]
               [escalation == FactorBaseWitnessEscalation::bounded ? 1 : 0]
                       .fetch_add(1, std::memory_order_relaxed);
    }
    static slong count(bool direct_required_prime_witness,
                       FactorBaseWitnessEscalation escalation) noexcept {
        return counts_[direct_required_prime_witness ? 1 : 0]
                      [escalation == FactorBaseWitnessEscalation::bounded ? 1
                                                                          : 0]
                              .load(std::memory_order_relaxed);
    }
    static void reset() noexcept {
        for (auto& row : counts_) {
            for (auto& entry : row) {
                entry.store(0, std::memory_order_relaxed);
            }
        }
    }

private:
    inline static std::atomic<slong> counts_[2][2]{};
};

// Under bounded escalation, when a required prime has no witness at the
// stage-0 caps, the witness search effort for that prime (lattice radius,
// twists, random tries, and the T2 factor-attempt and element-step caps)
// doubles at each of at most this many further stages before the prime is
// reported unwitnessed.  The doubling, the choice of caps that double, and
// the number of stages are Silex's own bounded choices with no upstream
// source.  The nearest reference schedule, buch2.c:SPLIT (used by
// bnftestprimes and isprincipal once the class group is known), keeps a fixed
// set of twisting directions and instead widens the number of Vbase primes in
// each random product, doubling its try limit as that number grows and then
// removing the limit; the Silex random tries keep the sub-factor base fixed.
// Only search effort changes: every witness found at any stage passes the
// same exact check.
inline constexpr slong kFactorBaseHonestyEscalationStages = 3;

struct FactorBaseHonestyScanAudit {
    slong rational_prime_checks = 0;
    slong checks_at_or_below_active_bound = 0;
    slong witness_targets = 0;
    slong witnessed_targets = 0;
    // Targets witnessed only after the stage-0 caps were raised, and the
    // highest escalation stage any target's search reached, whether or not
    // that search found a witness.
    slong escalated_witnessed_targets = 0;
    slong max_search_stage = 0;
    ulong final_random_state = 0;
    // Set when the scan stops at a required prime with no witness.
    FactorBaseHonestyFailure unwitnessed;
};

bool factor_base_honesty_primitive_part(
        Ideal& ideal,
        Element& back_multiplier) noexcept;

bool factor_base_honesty_reduce_large_ideal(
        Ideal& ideal,
        Element& back_multiplier,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;

// Sets `out` to the discriminant D of `order` and returns true when `order`
// is the maximal order of a quadratic field whose defining polynomial is
// monic integral (the p-independent decompose_prime preconditions), on either
// number-field backend.  The honesty scan then skips
// every rational prime p with Kronecker symbol (D/p) = -1 without building
// its decomposition.  For a quadratic field of discriminant D, p is inert
// exactly when the Kronecker symbol (D/p) is -1, p = 2 included (Cohen,
// GTM 138, Proposition 5.1.4(2); Corollary 4.8.12).  Then
// pO_K is the only prime above p (e = 1, f = 2), and
// select_factor_base_proof_targets returns no target for it at any bound,
// so the per-prime check is honest without a witness search; the reference
// buch2.c:be_honest and bnftestprimes likewise leave nothing to test for an
// inert p.  The reference quadratic class group code reads the splitting of
// p from kronecker(D, p) in the same way (buch1.c:nthidealquad).
bool factor_base_honesty_quadratic_discriminant(
        flint::FmpzRef out,
        const Order& order) noexcept;

// use_direct_required_prime_witness selects the capped T2 search instead of
// the scalar/lattice search.  Both classify each candidate with
// factor_base_principal_witness under the given predicate, and both treat a
// null subfactor_base_schedule, or a schedule entry outside the base, as no
// further twists rather than as a failed scan.
bool factor_base_honesty_check(
        bool& honest,
        const FactorBase& base,
        flint::FmpzConstRef active_bound,
        flint::FmpzConstRef required_bound,
        const SubfactorBaseSchedule* subfactor_base_schedule,
        ulong random_seed,
        bool use_direct_required_prime_witness,
        slong ideal_reduction_precision,
        FactorBaseWitnessEscalation escalation,
        const DiagnosticsContext* diagnostics,
        FactorBaseHonestyScanAudit* audit = nullptr,
        FactorBaseWitnessPredicate predicate =
                FactorBaseWitnessPredicate::selected) noexcept;

}  // namespace silex::detail::relation_search
