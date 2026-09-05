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
};

bool factor_base_principal_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const OrderElement& generator,
        FactorBaseWitnessPredicate predicate = FactorBaseWitnessPredicate::selected,
        const DiagnosticsContext* diagnostics = nullptr,
        FactorBaseWitnessAudit* audit = nullptr) noexcept;

struct FactorBaseHonestyScanAudit {
    slong rational_prime_checks = 0;
    slong checks_at_or_below_active_bound = 0;
};

bool factor_base_honesty_primitive_part(
        Ideal& ideal,
        Element& back_multiplier) noexcept;

bool factor_base_honesty_reduce_large_ideal(
        Ideal& ideal,
        Element& back_multiplier,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept;

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
        FactorBaseHonestyScanAudit* audit = nullptr) noexcept;

}  // namespace silex::detail::relation_search
