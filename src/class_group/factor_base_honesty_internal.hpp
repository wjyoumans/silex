#pragma once

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

// The predicate selector affects only the legacy scalar/lattice search;
// use_direct_required_prime_witness independently selects the T2 search.
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
