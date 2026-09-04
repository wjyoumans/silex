#include <silex/class_group.hpp>

#include "class_group_internal.hpp"

#include "ideal_minkowski_embedding_internal.hpp"
#include "relation_search_driver_internal.hpp"
#include "relation_admission_cache_internal.hpp"
#include "relation_search_internal.hpp"

#include <cstddef>
#include <cstdio>
#include <vector>

namespace silex::detail::relation_search {

// Source-derived caps limit retained relations per searched ideal in both the
// small-norm and random-product searches.
constexpr slong kSmallNormRelationsPerIdeal = 4;

static bool build_full_rank_log_indices(std::vector<slong>& indices,
                                      const ClassGroupContext& context,
                                      const FactorBase& base,
                                      const std::vector<slong>& permutation)
        noexcept {
    const slong target = context.relation_rank() < base.length()
            ? context.relation_rank()
            : base.length();
    if (target < 0 ||
        static_cast<slong>(permutation.size()) != base.length() ||
        static_cast<slong>(permutation.size()) < target) {
        return false;
    }

    indices.clear();
    indices.reserve(static_cast<std::size_t>(target));
    for (slong i = 0; i < target; ++i) {
        const slong index = permutation[static_cast<std::size_t>(i)];
        if (index < 0 || index >= base.length()) {
            return false;
        }
        indices.push_back(index);
    }
    return true;
}

static bool floor_log_fmpz(slong& out,
                    flint::FmpzConstRef bound,
                    flint::FmpzConstRef base) noexcept {
    out = 0;
    if (flint::fmpz_cmp_ui(bound, 1) < 0 ||
        flint::fmpz_cmp_ui(base, 2) < 0) {
        return false;
    }

    flint::Fmpz power;
    flint::Fmpz next;
    flint::fmpz_one(flint::FmpzRef(power));
    for (;;) {
        flint::fmpz_mul(flint::FmpzRef(next),
                        flint::FmpzConstRef(power), base);
        if (flint::fmpz_cmp(flint::FmpzConstRef(next), bound) > 0) {
            return true;
        }
        if (out == WORD_MAX) {
            return false;
        }
        power.swap(next);
        ++out;
    }
}

static bool ideal_pow_ui(Ideal& out,
                  const Ideal& input,
                  slong exponent,
                  Ideal& product) noexcept {
    if (exponent < 0 || !out.one()) {
        return false;
    }

    for (slong i = 0; i < exponent; ++i) {
        if (!product.multiply(out, input)) {
            return false;
        }
        out.swap(product);
    }
    return true;
}

static bool build_small_norm_multiplier(Ideal& multiplier,
                                      slong& exponent,
                                      PrimeIdeal& multiplier_prime,
                                      Ideal& multiplier_base,
                                      Ideal& product,
                                      const FactorBase& base,
                                      slong j0) noexcept {
    exponent = 0;
    if (j0 <= 0) {
        return multiplier.one();
    }

    const slong multiplier_index = j0 - 1;
    if (multiplier_index < 0 || multiplier_index >= base.length() ||
        base.length() <= 0 ||
        !base.prime(multiplier_prime, multiplier_index) ||
        !multiplier_prime.get_ideal(multiplier_base)) {
        return false;
    }

    const PrimeIdeal* last_prime = base.prime_at(base.length() - 1);
    flint::Fmpz last_norm;
    flint::Fmpz multiplier_norm;
    flint::Fmpz bound;
    if (last_prime == nullptr ||
        !last_prime->norm(flint::FmpzRef(last_norm)) ||
        !multiplier_prime.norm(flint::FmpzRef(multiplier_norm))) {
        return false;
    }
    flint::fmpz_mul(flint::FmpzRef(bound),
                    flint::FmpzConstRef(last_norm),
                    flint::FmpzConstRef(last_norm));
    if (!floor_log_fmpz(exponent, flint::FmpzConstRef(bound),
                        flint::FmpzConstRef(multiplier_norm))) {
        return false;
    }

    return ideal_pow_ui(multiplier, multiplier_base, exponent, product);
}

static void filter_small_norm_multiplier_useful(
        std::vector<slong>& useful,
        std::vector<slong>& small_multiplier,
        slong j0) noexcept {
    if (j0 <= 0 ||
        small_multiplier.size() < useful.size()) {
        return;
    }

    const slong j0_index = j0 - 1;
    if (j0_index < 0 ||
        j0_index >= static_cast<slong>(small_multiplier.size())) {
        useful.clear();
        return;
    }

    const slong previous_multiplier = small_multiplier[
            static_cast<std::size_t>(j0_index)];
    std::vector<slong> filtered;
    filtered.reserve(useful.size());
    for (slong index : useful) {
        if (index < 0 ||
            index >= static_cast<slong>(small_multiplier.size())) {
            continue;
        }
        if (index + 1 > previous_multiplier) {
            small_multiplier[static_cast<std::size_t>(index)] = j0;
            filtered.push_back(index);
        }
    }
    useful.swap(filtered);
}

static bool small_norm_skip_trivial_multiplier_relation(
        const PrimeIdeal& prime,
        slong multiplier_exponent,
        slong degree) noexcept {
    const slong ramification = prime.ramification_index();
    const slong residue_degree = prime.residue_degree();
    return ramification > 0 &&
           ((multiplier_exponent + 1) % ramification == 0) &&
           ramification * residue_degree == degree;
}

class FullRankLogScope {
public:
    FullRankLogScope() noexcept = default;
    ~FullRankLogScope() noexcept {
        if (active_) {
            end_selected_pivot_recollection(*cache_, indices_, n_);
        }
    }

    FullRankLogScope(const FullRankLogScope&) = delete;
    FullRankLogScope& operator=(const FullRankLogScope&) = delete;

    bool begin(detail::RelationAdmissionCache& cache,
               const std::vector<slong>& indices,
               slong n) noexcept {
        if (!begin_selected_pivot_recollection(cache, indices, n)) {
            return false;
        }
        cache_ = &cache;
        indices_ = indices;
        n_ = n;
        active_ = true;
        return true;
    }

private:
    detail::RelationAdmissionCache* cache_ = nullptr;
    std::vector<slong> indices_;
    slong n_ = 0;
    bool active_ = false;
};

static void log_small_norm_trace(const DiagnosticsContext* diagnostics,
                               const char* phase,
                               slong j0,
                               bool full_rank_log_selection,
                               bool apply_multiplier_filter,
                               slong useful_count,
                               slong target_relation_count,
                               const detail::RelationAdmissionCache& cache,
                               const ClassGroupContext& context,
                               const RelationSearchControlState& route_state)
        noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (!log_enabled(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::trace)) {
        return;
    }
    char detail[256];
    std::snprintf(detail, sizeof(detail),
                  "phase=%s j0=%ld lie=%d filter=%d useful=%ld "
                  "target=%ld relations=%ld rank=%ld missing=%ld relsup=%ld "
                  "done_small=%ld need=%ld",
                  phase, static_cast<long>(j0), full_rank_log_selection ? 1 : 0,
                  apply_multiplier_filter ? 1 : 0,
                  static_cast<long>(useful_count),
                  static_cast<long>(target_relation_count),
                  static_cast<long>(context.relation_count()),
                  static_cast<long>(context.relation_rank()),
                  static_cast<long>(cache.missing),
                  static_cast<long>(cache.relation_surplus),
                  static_cast<long>(route_state.done_small),
                  static_cast<long>(route_state.completion.need));
    log_emit(diagnostics, DiagnosticsModule::class_group, LogLevel::trace,
             __func__, "small_norm pass", detail);
#else
    (void)diagnostics;
    (void)phase;
    (void)j0;
    (void)full_rank_log_selection;
    (void)apply_multiplier_filter;
    (void)useful_count;
    (void)target_relation_count;
    (void)cache;
    (void)context;
    (void)route_state;
#endif
}

bool collect_small_norm_relations(
        ClassGroupContext& context,
        const Order& order,
        NormPrefilter* norm_prefilter,
        slong target_relation_count,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        detail::RelationAdmissionCache& cache,
        const SubfactorBaseSchedule& subfb_state,
        RelationSearchControlState& route_state,
        slong j0,
        bool apply_multiplier_filter,
        std::vector<slong>& small_multiplier,
        bool allow_silex_ideal_lattice_input
) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    const FactorBase* base = context.factor_base();
    if (base == nullptr || target_relation_count <= context.relation_count() ||
        options.ideal_search_radius <= 0 ||
        options.max_candidates <= 0 ||
        options.max_relations <= context.relation_count()) {
        return true;
    }

    std::vector<slong> useful;
    FullRankLogScope log_scope;
    bool full_rank_log_selection = route_state.finish_full_rank_relation_active &&
            context.has_presentation() && route_state.done_small % 2 != 0;
    if (full_rank_log_selection) {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.small_norm.useful_indices");
            if (!build_full_rank_log_indices(
                        useful, context, *base, subfb_state.permutation)) {
                return false;
            }
        }
        if (!useful.empty()) {
            if (!log_scope.begin(cache, useful, base->length())) {
                return false;
            }
            const slong lie_target =
                    context.relation_count() +
                    static_cast<slong>(useful.size());
            target_relation_count = lie_target < options.max_relations
                    ? lie_target
                    : options.max_relations;
        }
    } else {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.small_norm.useful_indices");
            if (!build_useful_pivot_indices(
                        useful, context, *base, subfb_state.permutation,
                        route_state.completion.need,
                        route_state.finish_unit_log_rotation_active,
                        route_state.finish_unit_log_rotation)) {
                return false;
            }
        }
        if (apply_multiplier_filter) {
            filter_small_norm_multiplier_useful(useful,
                                                     small_multiplier, j0);
        }
    }
    log_small_norm_trace(context.diagnostics(), "begin", j0,
                              full_rank_log_selection, apply_multiplier_filter,
                              static_cast<slong>(useful.size()),
                              target_relation_count, cache, context,
                              route_state);
    if (useful.empty()) {
        return true;
    }

    if (!allow_silex_ideal_lattice_input) {
        SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "small_norm stopped: Silex ideal-lattice input "
                  "unavailable for this route");
        return false;
    }

    PrimeIdeal prime(order);
    PrimeIdeal multiplier_prime(order);
    Ideal ideal(order);
    Ideal multiplier(order);
    Ideal multiplier_base(order);
    Ideal multiplier_product(order);
    Ideal search_ideal(order);
    if (!prime.is_defined() || !multiplier_prime.is_defined() ||
        !ideal.is_defined() || !multiplier.is_defined() ||
        !multiplier_base.is_defined() || !multiplier_product.is_defined() ||
        !search_ideal.is_defined()) {
        return false;
    }
    slong multiplier_exponent = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.small_norm.multiplier");
        if (!build_small_norm_multiplier(
                    multiplier, multiplier_exponent, multiplier_prime,
                    multiplier_base, multiplier_product, *base, j0)) {
            return false;
        }
    }

    const slong max_relations =
            target_relation_count < options.max_relations
                    ? target_relation_count
                    : options.max_relations;
    flint::Fmpz denominator;
    fmpz_one(denominator.raw());
    detail::OrderMinkowskiEmbeddingCache t2_embedding_cache;

    auto stop_small_norm_pass = [&]() noexcept {
        return context.relation_count() >= target_relation_count ||
               route_state.candidates_tried >= options.max_candidates ||
               route_state.accepted_relations >= max_relations;
    };
    auto build_search_ideal = [&](slong index, bool& skip) noexcept {
        skip = false;
        if (!base->prime(prime, index)) {
            return false;
        }
        if (j0 > 0 && index == j0 - 1 &&
            small_norm_skip_trivial_multiplier_relation(
                    prime, multiplier_exponent, order.degree())) {
            skip = true;
            return true;
        }
        // reference small_norm enumerates p0 * P_j (or P_j) as an ideal lattice.
        // Silex Ideal multiplication and HNF provide that route-neutral
        // lattice in the active maximal-order basis.
        return prime.get_ideal(ideal) &&
               (j0 > 0 ? search_ideal.multiply(multiplier, ideal)
                       : search_ideal.set(ideal));
    };

    for (std::size_t pos = useful.size(); pos > 0; --pos) {
        if (stop_small_norm_pass()) {
            return true;
        }

        const slong useful_index = useful[pos - 1];

        bool skip = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.small_norm.build_search_ideal");
            if (!build_search_ideal(useful_index, skip)) {
                return false;
            }
        }
        if (skip) {
            continue;
        }

        const slong relations_before_ideal = context.relation_count();
        const slong accepted_before_ideal = route_state.accepted_relations;
        const slong candidates_before_ideal = route_state.candidates_tried;
        if (log_enabled(diagnostics, DiagnosticsModule::class_group,
                        LogLevel::trace)) {
            char detail[256];
            std::snprintf(
                    detail, sizeof(detail),
                    "j0=%ld useful_index=%ld index=%ld useful_pos=%ld "
                    "useful_count=%ld target=%ld relations=%ld accepted=%ld "
                    "candidates=%ld",
                    static_cast<long>(j0),
                    static_cast<long>(useful_index),
                    static_cast<long>(useful_index + 1),
                    static_cast<long>(pos),
                    static_cast<long>(useful.size()),
                    static_cast<long>(target_relation_count),
                    static_cast<long>(relations_before_ideal),
                    static_cast<long>(accepted_before_ideal),
                    static_cast<long>(candidates_before_ideal));
            log_emit(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::trace, __func__,
                     "small_norm search ideal begin", detail);
        }
        bool ideal_goal_reached = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.small_norm.collect_ideal");
            if (!collect_private_direct_t2_relations_den(
                        context, search_ideal,
                        flint::FmpzConstRef(denominator), norm_prefilter,
                        ClassGroupRelationSource::Search,
                        target_relation_kernel_units,
                        options.max_candidates, max_relations,
                        kSmallNormRelationsPerIdeal,
                        route_state.candidates_tried,
                        route_state.accepted_relations,
                        options.requested_certification, ideal_goal_reached,
                        cache, false, &t2_embedding_cache)) {
                return false;
            }
        }
        if (log_enabled(diagnostics, DiagnosticsModule::class_group,
                        LogLevel::trace)) {
            char detail[256];
            std::snprintf(
                    detail, sizeof(detail),
                    "j0=%ld useful_index=%ld index=%ld relations_delta=%ld "
                    "accepted_delta=%ld candidates_delta=%ld relations=%ld "
                    "accepted=%ld candidates=%ld goal=%d",
                    static_cast<long>(j0),
                    static_cast<long>(useful_index),
                    static_cast<long>(useful_index + 1),
                    static_cast<long>(context.relation_count() -
                                      relations_before_ideal),
                    static_cast<long>(route_state.accepted_relations -
                                      accepted_before_ideal),
                    static_cast<long>(route_state.candidates_tried -
                                      candidates_before_ideal),
                    static_cast<long>(context.relation_count()),
                    static_cast<long>(route_state.accepted_relations),
                    static_cast<long>(route_state.candidates_tried),
                    ideal_goal_reached ? 1 : 0);
            log_emit(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::trace, __func__,
                     "small_norm search ideal end", detail);
        }
    }

    log_small_norm_trace(context.diagnostics(), "end", j0,
                              full_rank_log_selection, apply_multiplier_filter,
                              static_cast<slong>(useful.size()),
                              target_relation_count, cache, context,
                              route_state);
    return true;
}

}  // namespace silex::detail::relation_search
