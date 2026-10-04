#pragma once

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

// When a required prime has no witness at the stage-0 caps, the witness
// search effort for that prime (lattice radius, twists, random tries, and the
// T2 factor-attempt and element-step caps) doubles at each of at most this
// many further stages before the prime is reported unwitnessed.  The
// reference buch2.c:be_honest stops after stage 0 (maxtry_HONEST) and
// enlarges the factor base; bnftestprimes calls SPLIT, which doubles its
// random-try limit and widens its twisting set in stages without a final
// cap.  Silex keeps the stages bounded.  Only search effort changes: every
// witness found at any stage passes the same exact check.
inline constexpr slong kFactorBaseHonestyEscalationStages = 3;

struct FactorBaseHonestyScanAudit {
    slong rational_prime_checks = 0;
    slong checks_at_or_below_active_bound = 0;
    slong witness_targets = 0;
    slong witnessed_targets = 0;
    // Targets witnessed only after the stage-0 caps were raised, and the
    // highest escalation stage any target reached.
    slong escalated_witnessed_targets = 0;
    slong max_witness_stage = 0;
    ulong final_random_state = 0;
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
        const DiagnosticsContext* diagnostics,
        FactorBaseHonestyScanAudit* audit = nullptr,
        FactorBaseWitnessPredicate predicate =
                FactorBaseWitnessPredicate::selected) noexcept;

}  // namespace silex::detail::relation_search
