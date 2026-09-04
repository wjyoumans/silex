#include <silex/class_group.hpp>

#include "class_group_internal.hpp"

#include "ideal_minkowski_embedding_internal.hpp"
#include "relation_search_driver_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"
#include "relation_search_internal.hpp"

#include <cstddef>
#include <vector>

namespace silex::detail::relation_search {

constexpr ulong kRandomRelationPhase = UWORD(6);
constexpr slong kRandomRelationsPerIdeal = 1;

static bool integral_ideal_is_scalar(const Ideal& ideal, bool& scalar) noexcept {
    scalar = false;
    if (!ideal.has_hnf()) {
        return false;
    }

    const slong n = ideal.degree();
    flint::FmpzMat hnf(n, n);
    if (!ideal.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            const fmpz* entry = fmpz_mat_entry(hnf.raw(), i, j);
            if (i == j) {
                if (fmpz_equal(entry, fmpz_mat_entry(hnf.raw(), 0, 0)) == 0) {
                    return true;
                }
            } else if (fmpz_is_zero(entry) == 0) {
                return true;
            }
        }
    }

    scalar = true;
    return true;
}

static bool build_random_product_ideal(
        Ideal& out,
        std::vector<slong>& exponents,
        ulong& random_state,
        PrimeIdeal& prime,
        Ideal& factor,
        Ideal& product,
        const FactorBase& base,
        const SubfactorBaseSchedule& state) noexcept {
    if (state.subfactor_base.empty()) {
        return true;
    }

    exponents.assign(state.subfactor_base.size(), 0);
    for (;;) {
        if (!out.one()) {
            return false;
        }

        bool nonzero = false;
        for (std::size_t i = 0; i < state.subfactor_base.size(); ++i) {
            const slong exponent =
                    static_cast<slong>(
                            next_relation_random_exponent(random_state));
            exponents[i] = exponent;
            if (exponent == 0) {
                continue;
            }

            nonzero = true;
            const slong index = state.subfactor_base[i];
            if (!base.prime(prime, index) || !prime.get_ideal(factor)) {
                return false;
            }
            for (slong k = 0; k < exponent; ++k) {
                if (!product.multiply(out, factor)) {
                    return false;
                }
                out.swap(product);
            }
        }

        bool scalar = false;
        if (!integral_ideal_is_scalar(out, scalar)) {
            return false;
        }
        if (nonzero && !scalar) {
            return true;
        }
    }
}

static void request_dependency_factor_base_restart(
        RelationSearchControlState& route_state,
        bool allow_past_half,
        const DiagnosticsContext* diagnostics) noexcept {
    route_state.bounded_full_factor_base_restart_pending = true;
    route_state.bounded_full_factor_base_restart_allow_past_half =
            allow_past_half;
    ++route_state.bounded_full_factor_base_restart_requests;
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              "dependency counters requested factor-base restart");
}

bool collect_random_relations(
        ClassGroupContext& context,
        const Order& order,
        NormPrefilter* norm_prefilter,
        slong target_relation_count,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        detail::RelationAdmissionCache& cache,
        SubfactorBaseSchedule& subfb_state,
        RelationSearchControlState& route_state,
        bool factor_base_restart_available,
        bool factor_base_restart_to_max_available,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.random_relations");
    const FactorBase* base = context.factor_base();
    if (base == nullptr || subfb_state.subfactor_base.empty() ||
        target_relation_count <= context.relation_count() ||
        options.ideal_search_radius <= 0 ||
        options.max_candidates <= 0 ||
        options.max_relations <= context.relation_count()) {
        return true;
    }

    std::vector<slong> useful;
    if (!build_useful_pivot_indices(useful, context, *base,
                                         subfb_state.permutation,
                                         route_state.completion.need,
                                         route_state
                                                 .finish_unit_log_rotation_active,
                                         route_state.finish_unit_log_rotation)) {
        return false;
    }
    if (useful.empty()) {
        return true;
    }

    const SubfactorBaseAdvanceDecision subfactor_base_decision =
            advance_subfactor_base_schedule(
                    subfb_state, route_state.completion, *base, useful,
                    factor_base_restart_available,
                    factor_base_restart_to_max_available);
    switch (subfactor_base_decision) {
    case SubfactorBaseAdvanceDecision::ready:
        break;
    case SubfactorBaseAdvanceDecision::
            request_guarded_factor_base_rebuild:
        request_dependency_factor_base_restart(route_state, false,
                                                    diagnostics);
        return true;
    case SubfactorBaseAdvanceDecision::
            request_factor_base_rebuild_to_limit:
        request_dependency_factor_base_restart(route_state, true,
                                                    diagnostics);
        return true;
    case SubfactorBaseAdvanceDecision::increase_failed:
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "random_relation_search subFB increase failed at bound cap");
        return false;
    case SubfactorBaseAdvanceDecision::rotation_failed:
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "random_relation_search subFB rotate failed at bound cap");
        return false;
    }

    const slong max_relations =
            target_relation_count < options.max_relations
                    ? target_relation_count
                    : options.max_relations;
    ulong random_state = relation_search_phase_seed(
            context, base->length(), kRandomRelationPhase,
            route_state.completion.dependent_trials);
    std::vector<slong> exponents;

    PrimeIdeal prime(order);
    PrimeIdeal pivot_prime(order);
    Ideal random_product(order);
    Ideal factor(order);
    Ideal random_product_step(order);
    Ideal pivot_ideal(order);
    Ideal search_ideal(order);
    flint::Fmpz denominator;
    fmpz_one(denominator.raw());
    if (!prime.is_defined() || !pivot_prime.is_defined() ||
        !random_product.is_defined() || !factor.is_defined() ||
        !random_product_step.is_defined() || !pivot_ideal.is_defined() ||
        !search_ideal.is_defined()) {
        return false;
    }
    detail::OrderMinkowskiEmbeddingCache t2_embedding_cache;

    if (subfb_state.subfactor_base.empty()) {
        return true;
    }
    if (!build_random_product_ideal(
                random_product, exponents, random_state, prime, factor,
                random_product_step, *base, subfb_state)) {
        return false;
    }

    auto stop_random_pass = [&]() noexcept {
        return context.relation_count() >= target_relation_count ||
               route_state.candidates_tried >= options.max_candidates ||
               route_state.accepted_relations >= max_relations;
    };
    auto build_search_ideal = [&](slong index) noexcept {
        return base->prime(pivot_prime, index) &&
               pivot_prime.get_ideal(pivot_ideal) &&
               search_ideal.multiply(random_product, pivot_ideal);
    };

    for (slong index : useful) {
        if (stop_random_pass()) {
            return true;
        }

        bool ideal_goal_reached = false;
        if (!build_search_ideal(index) ||
            !collect_private_direct_t2_relations_den(
                    context, search_ideal,
                    flint::FmpzConstRef(denominator),
                    norm_prefilter, ClassGroupRelationSource::RandomProduct,
                    target_relation_kernel_units,
                    options.max_candidates, max_relations,
                    kRandomRelationsPerIdeal,
                    route_state.candidates_tried,
                    route_state.accepted_relations,
                    options.requested_certification, ideal_goal_reached,
                    cache, true, &t2_embedding_cache)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "random_relation_search direct T2 relation search failed");
            return false;
        }
    }

    return true;
}

}  // namespace silex::detail::relation_search
