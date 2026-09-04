#include <silex/class_group.hpp>

#include <silex/lat.hpp>

#include "ideal_t2_enumeration_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {

bool publish_and_check_compute_goal(ClassGroupContext& context,
                                    slong target_relation_kernel_units,
                                    CertificationMode requested_certification)
        noexcept {
    if (context.relation_rank() < context.generator_count()) {
        SILEX_PROFILE_EVENT(
                context.diagnostics(), DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_publish.skip_rank_deficient");
        return false;
    }
    const bool defer_presentation =
            detail::ClassGroupRelationSearchAccess::
                    defer_native_goal_publication(context) &&
            requested_certification != CertificationMode::proven;
    if (!defer_presentation) {
        if (!context.publish_presentation()) {
            return false;
        }

        if (requested_certification == CertificationMode::proven) {
            return context.try_certify_quadratic(requested_certification);
        }

        return context.relation_kernel_unit_count() >=
                target_relation_kernel_units;
    }
    // Both source baselines keep the full-rank relation module current while
    // deferring group/SNF publication.  For a full-rank m-by-n relation
    // matrix, m - n is the exact relation-kernel dimension.
    if (!detail::ClassGroupRelationSearchAccess::
                sync_row_module_checkpoint(context)) {
        return false;
    }
    if (detail::ClassGroupRelationSearchAccess::relation_kernel_row_count(
                context) <
                target_relation_kernel_units) {
        return false;
    }
    if (!context.publish_presentation()) {
        return false;
    }
    return true;
}

struct IdealLatticeVisitContext {
    ClassGroupContext* context = nullptr;
    const Order* order = nullptr;
    const flint::FmpzMat* basis = nullptr;
    const fmpz* denominator = nullptr;
    NormPrefilter* norm_prefilter = nullptr;
    slong max_candidates = 0;
    slong max_relations = 0;
    slong target_relation_kernel_units = 0;
    slong* candidates_tried = nullptr;
    slong* accepted_relations = nullptr;
    bool* goal_reached = nullptr;
    slong skipped_before = 0;
    CertificationMode requested_certification = CertificationMode::unknown;
    ClassGroupRelationSource source = ClassGroupRelationSource::Search;
    detail::RelationAdmissionCache* admission_cache = nullptr;
    bool random_relation = false;
    slong factor_attempts = 0;
    bool factor_attempt_limit_reached = false;
    bool failed = false;
};

int ideal_lattice_visit(const fmpz_mat_t coefficients, void* user) {
    IdealLatticeVisitContext* visit =
            static_cast<IdealLatticeVisitContext*>(user);
    if (visit == nullptr || visit->context == nullptr ||
        visit->order == nullptr || visit->basis == nullptr ||
        visit->candidates_tried == nullptr ||
        visit->accepted_relations == nullptr ||
        visit->goal_reached == nullptr) {
        return 0;
    }

    if (*visit->goal_reached ||
        *visit->candidates_tried >= visit->max_candidates ||
        *visit->accepted_relations >= visit->max_relations) {
        return 0;
    }

    if (visit->admission_cache != nullptr &&
        !fmpz_mat_single_row_is_primitive(
                flint::FmpzMatConstRef(coefficients))) {
        return 1;
    }

    flint::FmpzMat coordinates(1, visit->order->degree());
    coordinates_from_lattice_combination(
            coordinates, flint::FmpzMatConstRef(coefficients),
            flint::FmpzMatConstRef(*visit->basis));
    bool partial_throttle_exit = false;
    if (visit->denominator == nullptr ||
        !try_coordinate_candidate_den(
                *visit->context, *visit->order, coordinates,
                flint::FmpzConstRef(visit->denominator),
                visit->norm_prefilter, visit->source,
                visit->target_relation_kernel_units,
                visit->max_candidates, visit->max_relations,
                *visit->candidates_tried, *visit->accepted_relations,
                visit->requested_certification, *visit->goal_reached,
                partial_throttle_exit, visit->admission_cache,
                visit->random_relation, &visit->factor_attempts,
                kMaxFactorAttempts,
                &visit->factor_attempt_limit_reached)) {
        visit->failed = true;
        return 0;
    }

    return visit->factor_attempt_limit_reached ||
                   *visit->goal_reached ||
                   *visit->candidates_tried >= visit->max_candidates ||
                   *visit->accepted_relations >= visit->max_relations ||
                   visit->context->skipped_dependent_relation_count() >
                           visit->skipped_before
            ? 0
            : 1;
}

bool ideal_lattice_shell_recurse(ClassGroupContext& context,
                                 const Order& order,
                                 flint::FmpzMat& basis,
                                 flint::FmpzMat& coefficients,
                                 slong shell,
                                 slong column,
                                 bool used_shell,
                                 flint::FmpzConstRef den,
                                 NormPrefilter* norm_prefilter,
                                 ulong restart_seed,
                                 ClassGroupRelationSource source,
                                 slong target_relation_kernel_units,
                                 slong max_candidates,
                                 slong max_relations,
                                 slong& candidates_tried,
                                 slong& accepted_relations,
                                 CertificationMode requested_certification,
                                 bool& goal_reached,
                                 detail::RelationAdmissionCache* admission_cache =
                                         nullptr,
                                 bool random_relation = false,
                                 slong* factor_attempts = nullptr,
                                 bool* factor_attempt_limit_reached =
                                         nullptr) noexcept {
    if ((factor_attempt_limit_reached != nullptr &&
         *factor_attempt_limit_reached) ||
        goal_reached ||
        candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    if (column == flint::fmpz_mat_ncols(coefficients)) {
        if (!used_shell) {
            return true;
        }

        if (admission_cache != nullptr &&
            !fmpz_mat_single_row_is_primitive(
                    flint::FmpzMatConstRef(coefficients))) {
            return true;
        }

        flint::FmpzMat coordinates(1, order.degree());
        coordinates_from_lattice_combination(
                coordinates, flint::FmpzMatConstRef(coefficients),
                flint::FmpzMatConstRef(basis));
        bool partial_throttle_exit = false;
        return try_coordinate_candidate_den(
                context, order, coordinates, den, norm_prefilter, source,
                target_relation_kernel_units, max_candidates, max_relations,
                candidates_tried, accepted_relations,
                requested_certification, goal_reached,
                partial_throttle_exit, admission_cache,
                random_relation, factor_attempts,
                kMaxFactorAttempts,
                factor_attempt_limit_reached);
    }

    const slong count = 2 * shell + 1;
    const slong start =
            static_cast<slong>(restart_seed % static_cast<ulong>(count));
    for (slong k = 0; k < count; ++k) {
        const slong c = -shell + ((start + k) % count);
        flint::fmpz_set_si(flint::fmpz_mat_entry(coefficients, 0, column), c);
        if (!ideal_lattice_shell_recurse(context, order, basis, coefficients,
                                         shell, column + 1,
                                         used_shell || c == shell ||
                                                 c == -shell,
                                         den, norm_prefilter, restart_seed,
                                         source,
                                         target_relation_kernel_units,
                                         max_candidates, max_relations,
                                         candidates_tried, accepted_relations,
                                         requested_certification,
                                         goal_reached, admission_cache,
                                         random_relation,
                                         factor_attempts,
                                         factor_attempt_limit_reached)) {
            return false;
        }
        if ((factor_attempt_limit_reached != nullptr &&
             *factor_attempt_limit_reached) ||
            goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }
    }

    return true;
}

bool collect_integral_ideal_lattice_relations_den(
        ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzConstRef den,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        ulong restart_seed,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::RelationAdmissionCache* admission_cache,
        bool random_relation,
        slong* factor_attempts,
        bool* factor_attempt_limit_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || ideal_search_radius <= 0 ||
        fmpz_sgn(den.raw()) <= 0) {
        return false;
    }

    flint::FmpzMat hnf(order->degree(), order->degree());
    lat::Lat lattice(order->degree());
    lat::Lat reduced(order->degree());
    if (!ideal.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }

    bool reduced_lattice_cached = false;
    if (reduced_lattice_cache != nullptr) {
        // The full HNF equality check makes hits exact; a collision-prone
        // fingerprint is not used to decide reuse.
        for (const detail::ReducedIdealLatticeCacheEntry& entry :
             reduced_lattice_cache->entries) {
            if (fmpz_mat_equal(entry.ideal_hnf.raw(), hnf.raw()) == 0) {
                continue;
            }
            SILEX_PROFILE_EVENT(
                    context.diagnostics(), DiagnosticsModule::class_group,
                    "class_group.native_post_finite_reduced_lattice_cache_hit");
            if (!reduced.set_basis(
                        flint::FmpzMatConstRef(entry.reduced_basis))) {
                return false;
            }
            reduced_lattice_cached = true;
            break;
        }
        if (!reduced_lattice_cached) {
            SILEX_PROFILE_EVENT(
                    context.diagnostics(), DiagnosticsModule::class_group,
                    "class_group.native_post_finite_reduced_lattice_cache_miss");
        }
    }

    if ((!reduced_lattice_cached &&
         (!lattice.set_basis(flint::FmpzMatConstRef(hnf)) ||
          !lattice.lll_reduce(reduced))) ||
        reduced.nrows() <= 0) {
        return false;
    }

    if (!reduced_lattice_cached && reduced_lattice_cache != nullptr) {
        detail::ReducedIdealLatticeCacheEntry entry(
                flint::FmpzMatConstRef(hnf), reduced.basis_ref());
        if (reduced_lattice_cache->entries.size() <
            detail::ReducedIdealLatticeCache::capacity) {
            reduced_lattice_cache->entries.push_back(std::move(entry));
        } else {
            SILEX_PROFILE_EVENT(
                    context.diagnostics(), DiagnosticsModule::class_group,
                    "class_group.native_post_finite_reduced_lattice_cache_evict");
            reduced_lattice_cache
                    ->entries[reduced_lattice_cache->next_eviction] =
                    std::move(entry);
            reduced_lattice_cache->next_eviction =
                    (reduced_lattice_cache->next_eviction + 1) %
                    detail::ReducedIdealLatticeCache::capacity;
        }
    }

    flint::FmpzMat basis(reduced.nrows(), order->degree());
    flint::FmpzMat coefficients(1, reduced.nrows());
    flint::Fmpz max_diag;
    flint::Fmpz gram_ii;
    flint::Arb bound_sq;
    if (!reduced.get_basis(flint::FmpzMatRef(basis))) {
        return false;
    }

    flint::fmpz_zero(flint::FmpzRef(max_diag));
    for (slong i = 0; i < flint::fmpz_mat_nrows(basis); ++i) {
        flint::fmpz_zero(flint::FmpzRef(gram_ii));
        for (slong k = 0; k < flint::fmpz_mat_ncols(basis); ++k) {
            flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(flint::FmpzMatConstRef(basis), i, k);
            flint::fmpz_addmul(flint::FmpzRef(gram_ii),
                               entry, entry);
        }
        if (flint::fmpz_cmp(flint::FmpzConstRef(gram_ii),
                            flint::FmpzConstRef(max_diag)) > 0) {
            flint::fmpz_set(flint::FmpzRef(max_diag),
                            flint::FmpzConstRef(gram_ii));
        }
    }

    flint::fmpz_mul_ui(flint::FmpzRef(max_diag),
                       flint::FmpzConstRef(max_diag),
                       static_cast<ulong>(reduced.nrows()));

    IdealLatticeVisitContext visit;
    visit.context = &context;
    visit.order = order;
    visit.basis = &basis;
    visit.denominator = den.raw();
    visit.norm_prefilter = norm_prefilter;
    visit.max_candidates = max_candidates;
    visit.max_relations = max_relations;
    visit.target_relation_kernel_units = target_relation_kernel_units;
    visit.candidates_tried = &candidates_tried;
    visit.accepted_relations = &accepted_relations;
    visit.goal_reached = &goal_reached;
    visit.requested_certification = requested_certification;
    visit.source = source;
    visit.admission_cache = admission_cache;
    visit.random_relation = random_relation;
    if (factor_attempts != nullptr) {
        visit.factor_attempts = *factor_attempts;
    }
    if (factor_attempt_limit_reached != nullptr) {
        visit.factor_attempt_limit_reached =
                *factor_attempt_limit_reached;
    }
    auto store_factor_attempt_state = [&]() noexcept {
        if (factor_attempts != nullptr) {
            *factor_attempts = visit.factor_attempts;
        }
        if (factor_attempt_limit_reached != nullptr) {
            *factor_attempt_limit_reached =
                    visit.factor_attempt_limit_reached;
        }
    };
    if (visit.factor_attempt_limit_reached) {
        store_factor_attempt_state();
        return true;
    }

    bool enum_ok = false;
    const slong accepted_before = accepted_relations;
    const slong skipped_before = context.skipped_dependent_relation_count();
    visit.skipped_before = skipped_before;
    for (slong tries = 0; tries < 4; ++tries) {
        flint::arb_set_fmpz(bound_sq, max_diag);
        enum_ok = reduced.enum_short_vectors_arb(
                flint::ArbConstRef(bound_sq), ideal_search_radius, 64,
                ideal_lattice_visit, &visit);
        if (!enum_ok && !visit.failed && !goal_reached &&
            candidates_tried < max_candidates &&
            accepted_relations < max_relations) {
            enum_ok = reduced.enum_short_vectors_arb(
                    flint::ArbConstRef(bound_sq), ideal_search_radius, 256,
                    ideal_lattice_visit, &visit);
        }

        if (visit.failed || goal_reached ||
            visit.factor_attempt_limit_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations ||
            context.skipped_dependent_relation_count() > skipped_before ||
            accepted_relations > accepted_before) {
            break;
        }

        flint::fmpz_mul_ui(flint::FmpzRef(max_diag),
                           flint::FmpzConstRef(max_diag), UWORD(2));
    }

    if (visit.failed) {
        store_factor_attempt_state();
        return false;
    }

    if (!enum_ok && !goal_reached &&
        !visit.factor_attempt_limit_reached &&
        candidates_tried < max_candidates &&
        accepted_relations < max_relations) {
        for (slong shell = 1; shell <= ideal_search_radius; ++shell) {
            if (!ideal_lattice_shell_recurse(
                        context, *order, basis, coefficients, shell, 0, false,
                        den, norm_prefilter, restart_seed, source,
                        target_relation_kernel_units,
                        max_candidates, max_relations, candidates_tried,
                        accepted_relations, requested_certification,
                        goal_reached, admission_cache,
                        random_relation, &visit.factor_attempts,
                        &visit.factor_attempt_limit_reached)) {
                store_factor_attempt_state();
                return false;
            }
            if (visit.factor_attempt_limit_reached ||
                goal_reached ||
                candidates_tried >= max_candidates ||
                accepted_relations >= max_relations) {
                break;
            }
        }
    }

    store_factor_attempt_state();
    return true;
}

bool collect_finite_quadratic_form_relations_den(
        ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzConstRef den,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong max_relations_per_ideal,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::RelationAdmissionCache& admission_cache,
        bool random_relation,
        slong* factor_attempts_inout = nullptr,
        bool* factor_attempt_limit_reached_out = nullptr,
        detail::OrderMinkowskiEmbeddingCache* t2_embedding_cache = nullptr)
        noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.finite_ideal_search");
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() ||
        flint::fmpz_sgn(den) <= 0) {
        return false;
    }
    if (goal_reached) {
        return true;
    }


    FiniteIdealT2EnumerationData t2_context;
    OrderCoordinateElementConversion coordinate_conversion;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_ideal_search.setup");
        if (!build_finite_ideal_t2_enumeration_data_with_retry(
                    t2_context, ideal, diagnostics, t2_embedding_cache) ||
            !coordinate_conversion.reset(*order)) {
            SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "direct T2 enumeration context build failed");
            return false;
        }
    }

    FiniteQuadraticFormEnumerationContext enumeration;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_ideal_search.enumeration_init");
        if (!enumeration.reset(t2_context.quadratic_form_data, order->degree())) {
            SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "direct T2 enumeration reset failed");
            return false;
        }
    }

    bool skip_first_scalar = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_ideal_search.scalar_first_row");
        if (!reduced_basis_first_row_is_scalar_rational(
                    skip_first_scalar, *order,
                    flint::FmpzMatConstRef(t2_context.basis))) {
            SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "direct T2 scalar-first-row check failed");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_ideal_search.enumeration_start");
        if (!enumeration.start(t2_context.initial_bound_value,
                               kMaxElementSteps,
                               skip_first_scalar)) {
            SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "direct T2 enumeration start failed");
            return false;
        }
    }
    flint::FmpzMat coefficients(1, order->degree());
    flint::FmpzMat coordinates(1, order->degree());
    const FactorBase* base = context.factor_base();
    Relation admission_scratch;
    if (base == nullptr || !admission_scratch.define(*base)) {
        SILEX_LOG(context.diagnostics(), DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "direct T2 scratch relation allocation failed");
        return false;
    }
    slong factor_attempts =
            factor_attempts_inout == nullptr
                    ? 0
                    : *factor_attempts_inout;
    bool factor_attempt_limit_reached =
            factor_attempt_limit_reached_out != nullptr &&
            *factor_attempt_limit_reached_out;
    auto store_factor_attempt_state = [&]() noexcept {
        if (factor_attempts_inout != nullptr) {
            *factor_attempts_inout = factor_attempts;
        }
        if (factor_attempt_limit_reached_out != nullptr) {
            *factor_attempt_limit_reached_out =
                    factor_attempt_limit_reached;
        }
    };
    const slong accepted_before = accepted_relations;
    if (factor_attempt_limit_reached) {
        store_factor_attempt_state();
        return true;
    }

    for (;;) {
        bool has_next = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_ideal_search.enumeration_next");
            has_next = enumeration.next();
        }
        if (!has_next) {
            break;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_ideal_search.current_row");
            if (!enumeration.current_row(flint::FmpzMatRef(coefficients))) {
                store_factor_attempt_state();
                SILEX_LOG(context.diagnostics(),
                          DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "direct T2 current row read failed");
                return false;
            }
        }
        bool primitive = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_ideal_search.primitive_check");
            primitive = fmpz_mat_single_row_is_primitive(
                    flint::FmpzMatConstRef(coefficients));
        }
        if (!primitive) {
            continue;
        }

        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_ideal_search.coordinates");
            coordinates_from_lattice_combination(
                    coordinates, flint::FmpzMatConstRef(coefficients),
                    flint::FmpzMatConstRef(t2_context.basis));
        }
        bool partial_throttle_exit = false;
        bool candidate_goal_reached = false;
        if (!try_coordinate_candidate_den(
                    context, *order, coordinates, den,
                    norm_prefilter, source,
                    target_relation_kernel_units,
                    max_candidates, max_relations,
                    candidates_tried, accepted_relations,
                    requested_certification, candidate_goal_reached,
                    partial_throttle_exit, &admission_cache,
                    random_relation, &factor_attempts,
                    kMaxFactorAttempts,
                    &factor_attempt_limit_reached,
                    &admission_scratch, &coordinate_conversion,
                    true)) {
            store_factor_attempt_state();
            return false;
        }

        if (factor_attempt_limit_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations ||
            (max_relations_per_ideal > 0 &&
             accepted_relations - accepted_before >=
                     max_relations_per_ideal)) {
            store_factor_attempt_state();
            return true;
        }
    }

    store_factor_attempt_state();
    return true;
}

bool collect_private_direct_t2_relations_den(
        ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzConstRef den,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong max_relations_per_ideal,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::RelationAdmissionCache& admission_cache,
        bool random_relation,
        detail::OrderMinkowskiEmbeddingCache* t2_embedding_cache = nullptr)
        noexcept {
    slong factor_attempts = 0;
    bool factor_attempt_limit_reached = false;
    return collect_finite_quadratic_form_relations_den(
            context, ideal, den, norm_prefilter, source,
            target_relation_kernel_units, max_candidates, max_relations,
            max_relations_per_ideal, candidates_tried, accepted_relations,
            requested_certification, goal_reached, admission_cache,
            random_relation,
            &factor_attempts, &factor_attempt_limit_reached,
            t2_embedding_cache);
}

bool collect_integral_ideal_lattice_relations(
        ClassGroupContext& context,
        const Ideal& ideal,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        ulong restart_seed,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache) noexcept {
    flint::Fmpz den;
    fmpz_one(den.raw());
    return collect_integral_ideal_lattice_relations_den(
            context, ideal, flint::FmpzConstRef(den), ideal_search_radius,
            norm_prefilter, restart_seed, source, target_relation_kernel_units,
            max_candidates, max_relations, candidates_tried,
            accepted_relations, requested_certification, goal_reached,
            nullptr, false, nullptr, nullptr, reduced_lattice_cache);
}


ulong relation_search_phase_seed(const ClassGroupContext& context,
                                 slong factor_base_length,
                                 ulong phase,
                                 slong phase_restarts) noexcept;



ulong relation_search_phase_seed(const ClassGroupContext& context,
                                 slong factor_base_length,
                                 ulong phase,
                                 slong phase_restarts) noexcept {
    ulong seed = static_cast<ulong>(factor_base_length);
    seed ^= static_cast<ulong>(context.relation_count()) +
            UWORD(0x9e3779b97f4a7c15);
    seed ^= static_cast<ulong>(context.relation_rank()) << 1;
    seed ^= phase << 7;
    seed ^= static_cast<ulong>(phase_restarts) << 13;
    return seed;
}

bool build_nonprincipal_indices(std::vector<slong>& nonprincipal,
                                ClassGroupContext& context,
                                const FactorBase& base) noexcept {
    nonprincipal.clear();
    nonprincipal.reserve(static_cast<std::size_t>(base.length()));
    for (slong i = 0; i < base.length(); ++i) {
        bool is_principal = false;
        if (!context.factor_base_prime_is_principal(is_principal, i)) {
            return false;
        }
        if (!is_principal) {
            nonprincipal.push_back(i);
        }
    }
    return true;
}

bool build_hnf_covered_flags(std::vector<char>& covered,
                             ClassGroupContext& context,
                             const FactorBase& base) noexcept {
    covered.assign(static_cast<std::size_t>(base.length()), 0);
    for (slong i = 0; i < base.length(); ++i) {
        bool is_covered = false;
        if (!context.factor_base_prime_is_hnf_covered(is_covered, i)) {
            return false;
        }
        covered[static_cast<std::size_t>(i)] = is_covered ? 1 : 0;
    }
    return true;
}

bool build_uncovered_indices_from_flags(
        std::vector<slong>& uncovered,
        const std::vector<char>& hnf_covered,
        const std::vector<slong>& nonprincipal) noexcept {
    uncovered.clear();
    uncovered.reserve(hnf_covered.size());
    for (std::size_t i = 0; i < hnf_covered.size(); ++i) {
        if (hnf_covered[i] == 0) {
            uncovered.push_back(static_cast<slong>(i));
        }
    }
    if (uncovered.empty()) {
        uncovered = nonprincipal;
    }
    return true;
}


slong max_slong_value(slong left, slong right) noexcept {
    return left < right ? right : left;
}

slong min_slong_value(slong left, slong right) noexcept {
    return left < right ? left : right;
}


}  // namespace detail::relation_search

}  // namespace silex
