#include <silex/class_group.hpp>

#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>

#include "class_group_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {

constexpr ulong kCoordsPhase = UWORD(1);
constexpr ulong kRandomProductPhase = UWORD(2);
constexpr ulong kIdealLatticePhase = UWORD(3);
constexpr ulong kProductIdealPhase = UWORD(4);
constexpr ulong kPivotRandomPhase = UWORD(5);
constexpr slong kSearchSubfbDefault = 3;
constexpr slong kSearchSubfbMax = 5;
constexpr slong kPostFiniteMinPhaseBudget = 4;
constexpr slong kPostFiniteMinTotalBudget = 16;
constexpr slong kPostFiniteMaxTotalBudget = 64;

bool multiply_ideal_by_factor_base_index(
        Ideal& out,
        PrimeIdeal& prime,
        Ideal& generic_factor,
        const Ideal& input,
        const FactorBase& base,
        slong index,
        const DiagnosticsContext* diagnostics) noexcept;

bool collect_coordinate_relations(ClassGroupContext& context,
                                  const Order& order,
                                  slong search_radius,
                                  NormPrefilter* norm_prefilter,
                                  OrderCoordinateElementConversion*
                                          coordinate_conversion,
                                  slong target_relation_kernel_units,
                                  slong max_candidates,
                                  slong max_relations,
                                  flint::FmpzMat& coordinates,
                                  ulong restart_seed,
                                  slong column,
                                  slong& candidates_tried,
                                  slong& accepted_relations,
                                  CertificationMode requested_certification,
                                  bool& goal_reached,
                                  bool& partial_throttle_exit) noexcept {
    if (goal_reached || partial_throttle_exit ||
        candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    if (column == order.degree()) {
        return try_coordinate_candidate(context, order, coordinates,
                                        norm_prefilter,
                                        ClassGroupRelationSource::Search,
                                        target_relation_kernel_units,
                                        max_candidates, max_relations,
                                        candidates_tried, accepted_relations,
                                        requested_certification,
                                        goal_reached, partial_throttle_exit,
                                        coordinate_conversion);
    }

    const slong count = 2 * search_radius + 1;
    const slong start =
            static_cast<slong>(restart_seed % static_cast<ulong>(count));
    for (slong k = 0; k < count; ++k) {
        const slong c = -search_radius + ((start + k) % count);
        flint::fmpz_set_si(flint::fmpz_mat_entry(coordinates, 0, column), c);
        if (!collect_coordinate_relations(context, order, search_radius,
                                          norm_prefilter,
                                          coordinate_conversion,
                                          target_relation_kernel_units,
                                          max_candidates, max_relations,
                                          coordinates, restart_seed,
                                          column + 1,
                                          candidates_tried, accepted_relations,
                                          requested_certification,
                                          goal_reached,
                                          partial_throttle_exit)) {
            return false;
        }
        if (goal_reached || partial_throttle_exit ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }
    }

    return true;
}

bool collect_prime_ideal_lattice_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong phase_restart,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || ideal_search_radius <= 0) {
        return true;
    }

    PrimeIdeal prime(order);
    Ideal ideal(order);
    if (!prime.is_defined() || !ideal.is_defined()) {
        return false;
    }

    const ulong seed = relation_search_phase_seed(
            context, base->length(), kIdealLatticePhase, phase_restart);
    slong start = phase_restart % base->length();
    if (start < 0) {
        start += base->length();
    }
    for (slong i = 0; i < base->length(); ++i) {
        const slong index = (start + i) % base->length();
        bool is_principal = false;
        if (goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }
        if (!context.factor_base_prime_is_principal(is_principal, index)) {
            return false;
        }
        if (is_principal) {
            continue;
        }
        if (!base->prime(prime, index) || !prime.get_ideal(ideal) ||
            !collect_integral_ideal_lattice_relations(
                    context, ideal, ideal_search_radius, norm_prefilter, seed,
                    ClassGroupRelationSource::Search,
                    target_relation_kernel_units, max_candidates,
                    max_relations, candidates_tried, accepted_relations,
                    requested_certification,
                    goal_reached, reduced_lattice_cache)) {
            return false;
        }
    }

    return true;
}

bool collect_product_ideal_lattice_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong phase_restart,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || ideal_search_radius <= 0) {
        return true;
    }

    PrimeIdeal left_prime(order);
    PrimeIdeal right_prime(order);
    Ideal left_ideal(order);
    Ideal right_ideal(order);
    Ideal product(order);
    if (!left_prime.is_defined() || !right_prime.is_defined() ||
        !left_ideal.is_defined() || !right_ideal.is_defined() ||
        !product.is_defined()) {
        return false;
    }

    const ulong seed = relation_search_phase_seed(
            context, base->length(), kProductIdealPhase, phase_restart);
    slong outer_start = phase_restart % base->length();
    if (outer_start < 0) {
        outer_start += base->length();
    }
    for (slong i = 0; i < base->length(); ++i) {
        const slong left_index = (outer_start + i) % base->length();
        bool left_principal = false;
        if (goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }
        if (!context.factor_base_prime_is_principal(left_principal,
                                                    left_index)) {
            return false;
        }
        if (left_principal) {
            continue;
        }
        if (!base->prime(left_prime, left_index) ||
            !left_prime.get_ideal(left_ideal)) {
            return false;
        }

        slong inner_start = (phase_restart + left_index) % base->length();
        if (inner_start < 0) {
            inner_start += base->length();
        }
        for (slong j = 0; j < base->length(); ++j) {
            const slong right_index = (inner_start + j) % base->length();
            if (right_index < left_index) {
                continue;
            }
            bool right_principal = false;
            if (goal_reached ||
                candidates_tried >= max_candidates ||
                accepted_relations >= max_relations) {
                return true;
            }
            if (!context.factor_base_prime_is_principal(right_principal,
                                                        right_index)) {
                return false;
            }
            if (right_principal) {
                continue;
            }
            if (!multiply_ideal_by_factor_base_index(
                        product, right_prime, right_ideal, left_ideal, *base,
                        right_index, context.diagnostics()) ||
                !collect_integral_ideal_lattice_relations(
                        context, product, ideal_search_radius, norm_prefilter,
                        seed,
                        ClassGroupRelationSource::Search,
                        target_relation_kernel_units, max_candidates,
                        max_relations, candidates_tried, accepted_relations,
                        requested_certification,
                        goal_reached, reduced_lattice_cache)) {
                return false;
            }
        }
    }

    return true;
}

ulong random_product_signature(const std::vector<slong>& indices,
                               slong len) noexcept {
    ulong h = UWORD(1469598103934665603);
    for (slong i = 0; i < len; ++i) {
        h ^= static_cast<ulong>(indices[static_cast<std::size_t>(i)] + 1);
        h *= UWORD(1099511628211);
        h ^= static_cast<ulong>(len);
        h *= UWORD(1099511628211);
    }
    return h;
}

void random_product_select(std::vector<slong>& indices,
                           slong& len,
                           const std::vector<slong>& nonprincipal,
                           ulong attempt,
                           ulong seed) noexcept {
    const slong nonprincipal_len =
            static_cast<slong>(nonprincipal.size());
    if (nonprincipal_len <= 0) {
        len = 0;
        return;
    }

    if (attempt == 0) {
        indices[0] = nonprincipal[0];
        indices[1] = nonprincipal[0];
        len = 2;
        return;
    }

    ulong state = seed ^ (UWORD(0x9e3779b97f4a7c15) + attempt);
    len = (state & 1) != 0 && nonprincipal_len >= 3 ? 3 : 2;
    if (len > nonprincipal_len + 1) {
        len = nonprincipal_len >= 3 ? 3 : 2;
    }

    for (slong i = 0; i < len; ++i) {
        state = state * UWORD(6364136223846793005) +
                UWORD(1442695040888963407);
        indices[static_cast<std::size_t>(i)] =
                nonprincipal[static_cast<std::size_t>(
                        state % static_cast<ulong>(nonprincipal_len))];
    }

    for (slong i = 0; i + 1 < len; ++i) {
        for (slong j = i + 1; j < len; ++j) {
            if (indices[static_cast<std::size_t>(j)] <
                indices[static_cast<std::size_t>(i)]) {
                const slong t = indices[static_cast<std::size_t>(i)];
                indices[static_cast<std::size_t>(i)] =
                        indices[static_cast<std::size_t>(j)];
                indices[static_cast<std::size_t>(j)] = t;
            }
        }
    }
}

slong pivot_subfb_rotate_threshold(slong subfb_len) noexcept {
    const slong width = subfb_len < 2 ? 2 : subfb_len;
    const slong threshold = ((width - 1) * 16) / 10;
    return threshold < 1 ? 1 : threshold;
}

slong pivot_subfb_widen_threshold(slong subfb_len) noexcept {
    const slong width = subfb_len < 2 ? 2 : subfb_len;
    return (width - 1) * 16;
}

slong build_subfb(std::vector<slong>& subfb,
                  const std::vector<slong>& nonprincipal,
                  slong rotation,
                  slong cap) noexcept {
    const slong np_len = static_cast<slong>(nonprincipal.size());
    cap = cap < kSearchSubfbMax ? cap : kSearchSubfbMax;
    const slong subfb_len = cap < np_len ? cap : np_len;
    subfb.clear();
    if (subfb_len <= 0 || np_len <= 0) {
        return 0;
    }

    slong rot = rotation % np_len;
    if (rot < 0) {
        rot += np_len;
    }

    for (slong i = 0; i < subfb_len; ++i) {
        subfb.push_back(
                nonprincipal[static_cast<std::size_t>((rot + i) % np_len)]);
    }
    return subfb_len;
}

slong pivot_start_from_window_seed(ulong seed,
                                   slong nonprincipal_len,
                                   slong uncovered_len) noexcept {
    if (nonprincipal_len <= 0 || uncovered_len <= 0) {
        return 0;
    }

    const slong window_start = static_cast<slong>(
            seed % static_cast<ulong>(nonprincipal_len));
    return window_start % uncovered_len;
}

slong build_scored_subfb(std::vector<slong>& subfb,
                         const std::vector<slong>& nonprincipal,
                         const std::vector<char>& hnf_covered,
                         const std::vector<slong>& accepted_score,
                         const std::vector<slong>& skipped_score,
                         const std::vector<slong>& rejected_score,
                         slong rotation,
                         slong cap,
                         bool& used_scoring) noexcept {
    used_scoring = false;
    const slong np_len = static_cast<slong>(nonprincipal.size());
    const slong score_len = static_cast<slong>(accepted_score.size());
    if (np_len <= 0 || score_len != static_cast<slong>(skipped_score.size()) ||
        score_len != static_cast<slong>(rejected_score.size()) ||
        score_len != static_cast<slong>(hnf_covered.size())) {
        return build_subfb(subfb, nonprincipal, rotation, cap);
    }

    cap = cap < kSearchSubfbMax ? cap : kSearchSubfbMax;
    const slong subfb_len = cap < np_len ? cap : np_len;
    if (subfb_len <= 0) {
        subfb.clear();
        return 0;
    }

    bool have_signal = false;
    for (slong index : nonprincipal) {
        if (index < 0 || index >= score_len) {
            return build_subfb(subfb, nonprincipal, rotation, cap);
        }
        if (accepted_score[static_cast<std::size_t>(index)] != 0 ||
            skipped_score[static_cast<std::size_t>(index)] != 0 ||
            rejected_score[static_cast<std::size_t>(index)] != 0) {
            have_signal = true;
            break;
        }
    }
    if (!have_signal) {
        return build_subfb(subfb, nonprincipal, rotation, cap);
    }

    slong rot = rotation % np_len;
    if (rot < 0) {
        rot += np_len;
    }

    std::vector<slong> selected;
    std::vector<slong> selected_scores;
    selected.reserve(static_cast<std::size_t>(subfb_len));
    selected_scores.reserve(static_cast<std::size_t>(subfb_len));

    for (slong i = 0; i < np_len; ++i) {
        const slong index =
                nonprincipal[static_cast<std::size_t>((rot + i) % np_len)];
        slong score =
                64 * accepted_score[static_cast<std::size_t>(index)] -
                8 * skipped_score[static_cast<std::size_t>(index)] -
                2 * rejected_score[static_cast<std::size_t>(index)];
        if (hnf_covered[static_cast<std::size_t>(index)] == 0) {
            ++score;
        }

        if (static_cast<slong>(selected.size()) < subfb_len) {
            selected.push_back(index);
            selected_scores.push_back(score);
        } else if (score <= selected_scores.back()) {
            continue;
        } else {
            selected.back() = index;
            selected_scores.back() = score;
        }

        std::size_t pos = selected.size() - 1;
        while (pos > 0 && selected_scores[pos] > selected_scores[pos - 1]) {
            const slong score_tmp = selected_scores[pos - 1];
            selected_scores[pos - 1] = selected_scores[pos];
            selected_scores[pos] = score_tmp;

            const slong index_tmp = selected[pos - 1];
            selected[pos - 1] = selected[pos];
            selected[pos] = index_tmp;
            --pos;
        }
    }

    subfb = std::move(selected);
    used_scoring = true;
    return static_cast<slong>(subfb.size());
}


void pivot_score_update(std::vector<slong>& accepted_score,
                        std::vector<slong>& skipped_score,
                        std::vector<slong>& rejected_score,
                        slong pivot_idx,
                        const std::vector<slong>& product_indices,
                        slong product_len,
                        slong accepted_delta,
                        slong skipped_delta,
                        slong rejected_delta) noexcept {
    const slong score_len = static_cast<slong>(accepted_score.size());
    if (pivot_idx < 0 || pivot_idx >= score_len ||
        score_len != static_cast<slong>(skipped_score.size()) ||
        score_len != static_cast<slong>(rejected_score.size())) {
        return;
    }

    accepted_score[static_cast<std::size_t>(pivot_idx)] += accepted_delta;
    skipped_score[static_cast<std::size_t>(pivot_idx)] += skipped_delta;
    rejected_score[static_cast<std::size_t>(pivot_idx)] += rejected_delta;

    for (slong i = 0; i < product_len; ++i) {
        const slong index = product_indices[static_cast<std::size_t>(i)];
        if (index < 0 || index >= score_len) {
            continue;
        }
        accepted_score[static_cast<std::size_t>(index)] += accepted_delta;
        skipped_score[static_cast<std::size_t>(index)] += skipped_delta;
        rejected_score[static_cast<std::size_t>(index)] += rejected_delta;
    }
}

bool pivot_subfb_adapt(std::vector<slong>& subfb,
                       slong& subfb_len,
                       slong& subfb_cap,
                       const std::vector<slong>& nonprincipal,
                       const std::vector<char>& hnf_covered,
                       const std::vector<slong>& accepted_score,
                       const std::vector<slong>& skipped_score,
                       const std::vector<slong>& rejected_score,
                       slong& subfb_rotation,
                       slong& dep_no_progress,
                       slong& rotate_threshold,
                       slong& widen_threshold,
                       ulong& seed,
                       ulong base_seed) noexcept {
    if (dep_no_progress < rotate_threshold) {
        return false;
    }

    if (dep_no_progress >= widen_threshold && subfb_cap < kSearchSubfbMax) {
        subfb_cap = kSearchSubfbMax;
    }

    ++subfb_rotation;
    bool used_scoring = false;
    subfb_len = build_scored_subfb(subfb, nonprincipal, hnf_covered,
                                   accepted_score, skipped_score,
                                   rejected_score,
                                   subfb_rotation, subfb_cap, used_scoring);
    (void)used_scoring;
    rotate_threshold = pivot_subfb_rotate_threshold(subfb_len);
    widen_threshold = pivot_subfb_widen_threshold(subfb_len);
    seed = base_seed ^ (static_cast<ulong>(subfb_rotation) << 17);
    dep_no_progress = 0;
    return true;
}

bool set_ideal_from_factor_base_index(Ideal& out,
                                      PrimeIdeal& prime,
                                      const FactorBase& base,
                                      slong index) noexcept {
    return base.prime(prime, index) && prime.get_ideal(out);
}

bool multiply_ideal_by_factor_base_index(
        Ideal& out,
        PrimeIdeal& prime,
        Ideal& generic_factor,
        const Ideal& input,
        const FactorBase& base,
        slong index,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = input.parent();
    if (order == nullptr || !base.prime(prime, index)) {
        return false;
    }

    flint::FmpzMat generator_coordinates(1, order->degree());
    if (!prime.kummer_generator_coordinates(
                flint::FmpzMatRef(generator_coordinates))) {
        return prime.get_ideal(generic_factor) &&
               out.multiply(input, generic_factor);
    }

    flint::Fmpz p;
    OrderElement generator(*order);
    if (!prime.rational_prime(flint::FmpzRef(p)) ||
        !generator.is_defined() ||
        !generator.set_coordinates(
                flint::FmpzMatConstRef(generator_coordinates))) {
        return false;
    }

    return detail::multiply_integral_ideal_by_two_generator(
            out, input, flint::FmpzConstRef(p), generator, diagnostics);
}

bool build_product_ideal_from_indices(Ideal& out,
                                      PrimeIdeal& prime,
                                      Ideal& factor,
                                      Ideal& next_product,
                                      const FactorBase& base,
                                      const std::vector<slong>& indices,
                                      slong len,
                                      const DiagnosticsContext* diagnostics)
        noexcept {
    if (len <= 0) {
        return false;
    }
    if (!set_ideal_from_factor_base_index(out, prime, base, indices[0])) {
        return false;
    }
    for (slong i = 1; i < len; ++i) {
        if (!multiply_ideal_by_factor_base_index(
                    next_product, prime, factor, out, base,
                    indices[static_cast<std::size_t>(i)], diagnostics)) {
            return false;
        }
        out.swap(next_product);
    }
    return true;
}

bool collect_pivot_random_ideal_lattice_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong phase_restart,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || order.degree() <= 2 || base->length() < 2 ||
        ideal_search_radius <= 0 || max_candidates <= 0) {
        return true;
    }

    std::vector<slong> nonprincipal;
    std::vector<slong> uncovered;
    std::vector<char> hnf_covered;
    std::vector<slong> subfb;
    std::vector<slong> accepted_score(
            static_cast<std::size_t>(base->length()), 0);
    std::vector<slong> skipped_score(
            static_cast<std::size_t>(base->length()), 0);
    std::vector<slong> rejected_score(
            static_cast<std::size_t>(base->length()), 0);
    if (!build_nonprincipal_indices(nonprincipal, context, *base) ||
        !build_hnf_covered_flags(hnf_covered, context, *base) ||
        !build_uncovered_indices_from_flags(uncovered, hnf_covered,
                                            nonprincipal)) {
        return false;
    }
    if (nonprincipal.empty() || uncovered.empty()) {
        return true;
    }

    slong subfb_cap = kSearchSubfbDefault;
    slong subfb_rotation = phase_restart;
    slong subfb_len =
            build_subfb(subfb, nonprincipal, subfb_rotation, subfb_cap);
    if (subfb_len < 1) {
        return true;
    }

    slong uncov_len = static_cast<slong>(uncovered.size());
    const slong phase_limit = uncov_len < 4 ? 4 : uncov_len;
    const slong attempt_limit = max_candidates < phase_limit
            ? max_candidates
            : phase_limit;
    ulong base_seed = relation_search_phase_seed(
            context, base->length(), kPivotRandomPhase, phase_restart);
    ulong seed = base_seed;
    slong dep_no_progress = 0;
    slong rotate_threshold = pivot_subfb_rotate_threshold(subfb_len);
    slong widen_threshold = pivot_subfb_widen_threshold(subfb_len);
    slong pivot_start = pivot_start_from_window_seed(
            seed, static_cast<slong>(nonprincipal.size()), uncov_len);
    slong mask_relation_count = context.relation_count();
    std::vector<slong> indices(kSearchSubfbMax, 0);

    PrimeIdeal prime(order);
    Ideal base_product(order);
    Ideal factor(order);
    Ideal next_product(order);
    Ideal product(order);
    if (!prime.is_defined() || !base_product.is_defined() ||
        !factor.is_defined() ||
        !next_product.is_defined() || !product.is_defined()) {
        return false;
    }

    for (slong attempt = 0; attempt < attempt_limit; ++attempt) {
        if (goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }

        if (context.relation_count() != mask_relation_count) {
            if (!build_nonprincipal_indices(nonprincipal, context, *base) ||
                !build_hnf_covered_flags(hnf_covered, context, *base) ||
                !build_uncovered_indices_from_flags(uncovered, hnf_covered,
                                                    nonprincipal)) {
                return false;
            }
            if (nonprincipal.empty() || uncovered.empty()) {
                return true;
            }
            uncov_len = static_cast<slong>(uncovered.size());
            subfb_len = build_subfb(subfb, nonprincipal, subfb_rotation,
                                    subfb_cap);
            if (subfb_len < 1) {
                return true;
            }
            rotate_threshold = pivot_subfb_rotate_threshold(subfb_len);
            widen_threshold = pivot_subfb_widen_threshold(subfb_len);
            base_seed = relation_search_phase_seed(
                    context, base->length(), kPivotRandomPhase,
                    phase_restart);
            seed = base_seed;
            pivot_start = pivot_start_from_window_seed(
                    seed, static_cast<slong>(nonprincipal.size()), uncov_len);
            mask_relation_count = context.relation_count();
        }

        slong product_len = 0;
        random_product_select(indices, product_len, subfb,
                              static_cast<ulong>(attempt), seed);
        if (product_len < 2) {
            continue;
        }
        if (!build_product_ideal_from_indices(base_product, prime, factor,
                                              next_product, *base, indices,
                                              product_len,
                                              context.diagnostics())) {
            return false;
        }

        const slong accepted_product_before = context.relation_count();
        for (slong j = 0; j < uncov_len; ++j) {
            if (goal_reached ||
                candidates_tried >= max_candidates ||
                accepted_relations >= max_relations) {
                return true;
            }

            const slong pivot_idx =
                    uncovered[static_cast<std::size_t>(
                            (pivot_start + j) % uncov_len)];
            const slong accepted_before = context.relation_count();
            const slong skipped_before =
                    context.skipped_dependent_relation_count();
            const slong rejected_before = norm_prefilter == nullptr
                    ? 0
                    : norm_prefilter->rejected_count;
            if (!multiply_ideal_by_factor_base_index(
                        product, prime, factor, base_product, *base,
                        pivot_idx, context.diagnostics()) ||
                !collect_integral_ideal_lattice_relations(
                        context, product, ideal_search_radius, norm_prefilter,
                        seed,
                        ClassGroupRelationSource::RandomProduct,
                        target_relation_kernel_units, max_candidates,
                        max_relations, candidates_tried, accepted_relations,
                        requested_certification,
                        goal_reached, reduced_lattice_cache)) {
                return false;
            }
            const slong accepted_after = context.relation_count();
            const slong skipped_after =
                    context.skipped_dependent_relation_count();
            const slong rejected_after = norm_prefilter == nullptr
                    ? 0
                    : norm_prefilter->rejected_count;
            pivot_score_update(accepted_score, skipped_score, rejected_score,
                               pivot_idx,
                               indices, product_len,
                               accepted_after - accepted_before,
                               skipped_after - skipped_before,
                               rejected_after - rejected_before);
            if (accepted_after > accepted_before) {
                dep_no_progress = 0;
            } else if (skipped_after > skipped_before) {
                dep_no_progress += skipped_after - skipped_before;
            }
            pivot_subfb_adapt(subfb, subfb_len, subfb_cap, nonprincipal,
                              hnf_covered, accepted_score, skipped_score,
                              rejected_score,
                              subfb_rotation, dep_no_progress,
                              rotate_threshold, widen_threshold, seed,
                              base_seed);
        }

        if (context.relation_count() == accepted_product_before) {
            ++dep_no_progress;
            pivot_subfb_adapt(subfb, subfb_len, subfb_cap, nonprincipal,
                              hnf_covered, accepted_score, skipped_score,
                              rejected_score,
                              subfb_rotation, dep_no_progress,
                              rotate_threshold, widen_threshold, seed,
                              base_seed);
        }

        pivot_start = (pivot_start + 1) % uncov_len;
    }

    return true;
}

bool collect_random_product_ideal_lattice_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong phase_restart,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || order.degree() <= 2 || base->length() < 2 ||
        ideal_search_radius <= 0 || max_candidates <= 0) {
        return true;
    }

    std::vector<slong> nonprincipal;
    if (!build_nonprincipal_indices(nonprincipal, context, *base)) {
        return false;
    }
    if (nonprincipal.empty()) {
        return true;
    }

    ulong seed = relation_search_phase_seed(
            context, base->length(), kRandomProductPhase, phase_restart);
    slong pivot_start =
            static_cast<slong>(seed % static_cast<ulong>(nonprincipal.size()));
    std::vector<slong> pivot_nonprincipal(nonprincipal.size(), 0);
    for (std::size_t i = 0; i < nonprincipal.size(); ++i) {
        pivot_nonprincipal[i] =
                nonprincipal[(static_cast<std::size_t>(pivot_start) + i) %
                             nonprincipal.size()];
    }
    slong mask_relation_count = context.relation_count();

    const slong phase_limit = base->length() > WORD_MAX / 4
            ? WORD_MAX
            : 4 * base->length();
    const slong bounded_phase_limit = phase_limit < 8 ? 8 : phase_limit;
    const slong attempt_limit = max_candidates < bounded_phase_limit
            ? max_candidates
            : bounded_phase_limit;
    std::vector<ulong> seen;
    seen.reserve(static_cast<std::size_t>(attempt_limit));
    std::vector<slong> indices(3, 0);

    PrimeIdeal prime(order);
    Ideal product(order);
    Ideal factor(order);
    Ideal next_product(order);
    if (!prime.is_defined() || !product.is_defined() ||
        !factor.is_defined() || !next_product.is_defined()) {
        return false;
    }

    for (slong attempt = 0; attempt < attempt_limit; ++attempt) {
        if (goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }

        if (context.relation_count() != mask_relation_count) {
            if (!build_nonprincipal_indices(nonprincipal, context, *base)) {
                return false;
            }
            if (nonprincipal.empty()) {
                return true;
            }
            pivot_nonprincipal = nonprincipal;
            seen.clear();
            seed = relation_search_phase_seed(
                    context, base->length(), kRandomProductPhase,
                    phase_restart);
            mask_relation_count = context.relation_count();
        }

        slong product_len = 0;
        random_product_select(indices, product_len, pivot_nonprincipal,
                              static_cast<ulong>(attempt), seed);
        if (product_len < 2) {
            continue;
        }

        const ulong signature = random_product_signature(indices, product_len);
        bool duplicate = false;
        for (ulong previous : seen) {
            if (previous == signature) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            continue;
        }
        seen.push_back(signature);

        if (!build_product_ideal_from_indices(product, prime, factor,
                                              next_product, *base, indices,
                                              product_len,
                                              context.diagnostics())) {
            return false;
        }

        if (!collect_integral_ideal_lattice_relations(
                    context, product, ideal_search_radius, norm_prefilter,
                    seed,
                    ClassGroupRelationSource::RandomProduct,
                    target_relation_kernel_units, max_candidates,
                    max_relations, candidates_tried, accepted_relations,
                    requested_certification,
                    goal_reached, reduced_lattice_cache)) {
            return false;
        }
    }

    return true;
}

bool target_indices_contain(const std::vector<slong>& indices,
                            slong value) noexcept {
    for (slong index : indices) {
        if (index == value) {
            return true;
        }
    }
    return false;
}

bool build_class_generator_target_indices(
        std::vector<slong>& targets,
        ClassGroupContext& context,
        const FactorBase& base) noexcept {
    targets.clear();
    targets.reserve(kSearchSubfbMax);

    for (slong i = 0; i < base.length() &&
                      static_cast<slong>(targets.size()) < kSearchSubfbMax;
         ++i) {
        bool is_principal = false;
        bool is_covered = false;
        if (!context.factor_base_prime_is_principal(is_principal, i) ||
            !context.factor_base_prime_is_hnf_covered(is_covered, i)) {
            return false;
        }
        if (!is_principal && !is_covered) {
            targets.push_back(i);
        }
    }

    for (slong i = 0; i < base.length() &&
                      static_cast<slong>(targets.size()) < kSearchSubfbMax;
         ++i) {
        bool is_principal = false;
        if (!context.factor_base_prime_is_principal(is_principal, i)) {
            return false;
        }
        if (!is_principal && !target_indices_contain(targets, i)) {
            targets.push_back(i);
        }
    }

    return true;
}

bool multiply_by_factor_base_prime(FractionalIdeal& out,
                                   const FractionalIdeal& input,
                                   PrimeIdeal& prime,
                                   FractionalIdeal& prime_ideal,
                                   FractionalIdeal& product,
                                   const FactorBase& base,
                                   slong index) noexcept {
    return base.prime(prime, index) &&
           detail::prime_to_fractional_ideal(prime_ideal, prime) &&
           product.multiply(input, prime_ideal) &&
           out.set(product);
}

bool collect_class_generator_target_relations(
        ClassGroupContext& context,
        const FractionalIdeal& target,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        bool& order_changed,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    if (candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    const Order* order = target.parent();
    if (order == nullptr || !context.has_presentation()) {
        return false;
    }

    flint::Fmpz before;
    flint::Fmpz after;
    FractionalIdeal inverse(*order);
    Ideal numerator(*order);
    flint::Fmpz den;
    bool refinement_goal_reached = false;
    if (!context.order(flint::FmpzRef(before)) ||
        !inverse.is_defined() || !numerator.is_defined() ||
        !inverse.invert(target) ||
        !inverse.get_integral_den(numerator, flint::FmpzRef(den)) ||
        !collect_integral_ideal_lattice_relations_den(
                context, numerator, flint::FmpzConstRef(den),
                ideal_search_radius, norm_prefilter, 0,
                ClassGroupRelationSource::ClassGenerator, WORD_MAX,
                max_candidates, max_relations, candidates_tried,
                accepted_relations, CertificationMode::unknown,
                refinement_goal_reached, nullptr, false, nullptr, nullptr,
                reduced_lattice_cache) ||
        !context.order(flint::FmpzRef(after))) {
        return false;
    }

    order_changed =
            fmpz_cmp(before.raw(), after.raw()) != 0;
    return true;
}

bool collect_class_generator_refinement_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || context.invariant_count() <= 0 ||
        base->length() <= 0 || candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    std::vector<slong> targets;
    if (!build_class_generator_target_indices(targets, context, *base)) {
        return false;
    }

    FractionalIdeal generator(order);
    FractionalIdeal target(order);
    FractionalIdeal pair_target(order);
    FractionalIdeal prime_ideal(order);
    FractionalIdeal product(order);
    PrimeIdeal prime(order);
    if (!generator.is_defined() || !target.is_defined() ||
        !pair_target.is_defined() || !prime_ideal.is_defined() ||
        !product.is_defined() || !prime.is_defined()) {
        return false;
    }

    for (slong i = 0; i < context.invariant_count(); ++i) {
        if (candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }

        bool order_changed = false;
        if (!context.invariant_generator(generator, i) ||
            !collect_class_generator_target_relations(
                    context, generator, ideal_search_radius, norm_prefilter,
                    max_candidates, max_relations, candidates_tried,
                    accepted_relations, order_changed,
                    reduced_lattice_cache)) {
            return false;
        }
        if (order_changed) {
            return true;
        }

        for (std::size_t left_pos = 0; left_pos < targets.size();
             ++left_pos) {
            if (candidates_tried >= max_candidates ||
                accepted_relations >= max_relations) {
                return true;
            }
            const slong left_index = targets[left_pos];
            if (!multiply_by_factor_base_prime(
                        target, generator, prime, prime_ideal, product, *base,
                        left_index) ||
                !collect_class_generator_target_relations(
                        context, target, ideal_search_radius, norm_prefilter,
                        max_candidates, max_relations, candidates_tried,
                        accepted_relations, order_changed,
                        reduced_lattice_cache)) {
                return false;
            }
            if (order_changed) {
                return true;
            }

            for (std::size_t right_pos = left_pos; right_pos < targets.size();
                 ++right_pos) {
                if (candidates_tried >= max_candidates ||
                    accepted_relations >= max_relations) {
                    return true;
                }
                const slong right_index = targets[right_pos];
                if (!multiply_by_factor_base_prime(
                            pair_target, target, prime, prime_ideal, product,
                            *base, right_index) ||
                    !collect_class_generator_target_relations(
                            context, pair_target, ideal_search_radius,
                            norm_prefilter, max_candidates, max_relations,
                            candidates_tried, accepted_relations,
                            order_changed, reduced_lattice_cache)) {
                    return false;
                }
                if (order_changed) {
                    return true;
                }
            }
        }
    }

    return true;
}


template <class RunPhase>
bool run_restartable_relation_phase(const Order& order,
                                    slong restart_limit,
                                    slong max_candidates,
                                    slong max_relations,
                                    slong& candidates_tried,
                                    slong& accepted_relations,
                                    bool& goal_reached,
                                    RunPhase run_phase,
                                    const bool* phase_exit = nullptr)
        noexcept {
    if (order.degree() <= 2) {
        restart_limit = 0;
    }

    slong accepted_before = accepted_relations;
    for (slong restart = 0; restart <= restart_limit; ++restart) {
        if (goal_reached ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }

        if (!run_phase(restart)) {
            return false;
        }

        if (phase_exit != nullptr && *phase_exit) {
            return true;
        }

        if (goal_reached ||
            accepted_relations > accepted_before ||
            candidates_tried >= max_candidates ||
            accepted_relations >= max_relations) {
            return true;
        }

        if (restart == restart_limit) {
            return true;
        }

        accepted_before = accepted_relations;
    }

    return true;
}

bool collect_post_finite_refinement_relations(
        ClassGroupContext& context,
        const Order& order,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        slong phase_budget_option,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations) noexcept {
    const FactorBase* base = context.factor_base();
    if (base == nullptr || order.degree() <= 2 || ideal_search_radius <= 0 ||
        !context.has_presentation() || candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    slong phase_budget = max_slong_value(kPostFiniteMinPhaseBudget,
                                         base->length());
    phase_budget = max_slong_value(phase_budget, phase_budget_option);
    slong total_budget = base->length() > WORD_MAX / 4
            ? kPostFiniteMaxTotalBudget
            : max_slong_value(kPostFiniteMinTotalBudget, 4 * base->length());
    total_budget = min_slong_value(kPostFiniteMaxTotalBudget, total_budget);

    const slong total_candidate_limit = max_candidates - candidates_tried <
                    total_budget
            ? max_candidates
            : candidates_tried + total_budget;
    const slong refinement_target = WORD_MAX;
    bool refinement_goal_reached = false;
    bool pivot_random_ran = false;
    detail::ReducedIdealLatticeCache* reduced_lattice_cache =
            phase_budget_option <= 0
            ? nullptr
            : detail::ClassGroupRelationSearchAccess::
                      native_post_finite_reduced_lattices(context);

    auto phase_limit = [&]() noexcept {
        if (total_candidate_limit - candidates_tried <= phase_budget) {
            return total_candidate_limit;
        }
        return candidates_tried + phase_budget;
    };

    if (!collect_class_generator_refinement_relations(
                context, order, ideal_search_radius, norm_prefilter,
                phase_limit(),
                max_relations, candidates_tried, accepted_relations,
                reduced_lattice_cache)) {
        return false;
    }

    if (candidates_tried < total_candidate_limit && base->length() >= 3) {
        pivot_random_ran = true;
        const slong cap = phase_limit();
        if (!run_restartable_relation_phase(
                    order, 1, cap, max_relations,
                    candidates_tried, accepted_relations,
                    refinement_goal_reached,
                    [&](slong restart) noexcept {
                        return collect_pivot_random_ideal_lattice_relations(
                                context, order, ideal_search_radius,
                                norm_prefilter, restart, refinement_target,
                                cap, max_relations, candidates_tried,
                                accepted_relations,
                                CertificationMode::unknown,
                                refinement_goal_reached,
                                reduced_lattice_cache);
                    })) {
            return false;
        }
    }

    if (candidates_tried < total_candidate_limit &&
        accepted_relations < max_relations) {
        const slong cap = phase_limit();
        if (!run_restartable_relation_phase(
                    order, 1, cap, max_relations,
                    candidates_tried, accepted_relations,
                    refinement_goal_reached,
                    [&](slong restart) noexcept {
                        return collect_random_product_ideal_lattice_relations(
                                context, order, ideal_search_radius,
                                norm_prefilter, restart, refinement_target,
                                cap, max_relations, candidates_tried,
                                accepted_relations,
                                CertificationMode::unknown,
                                refinement_goal_reached,
                                reduced_lattice_cache);
                    })) {
            return false;
        }
    }

    if (candidates_tried < total_candidate_limit &&
        accepted_relations < max_relations && !pivot_random_ran) {
        const slong cap = phase_limit();
        if (!run_restartable_relation_phase(
                    order, 1, cap, max_relations,
                    candidates_tried, accepted_relations,
                    refinement_goal_reached,
                    [&](slong restart) noexcept {
                        return collect_pivot_random_ideal_lattice_relations(
                                context, order, ideal_search_radius,
                                norm_prefilter, restart, refinement_target,
                                cap, max_relations, candidates_tried,
                                accepted_relations,
                                CertificationMode::unknown,
                                refinement_goal_reached,
                                reduced_lattice_cache);
                    })) {
            return false;
        }
    }

    if (candidates_tried < total_candidate_limit &&
        accepted_relations < max_relations) {
        const slong cap = phase_limit();
        if (!run_restartable_relation_phase(
                    order, 1, cap, max_relations,
                    candidates_tried, accepted_relations,
                    refinement_goal_reached,
                    [&](slong restart) noexcept {
                        return collect_prime_ideal_lattice_relations(
                                context, order, ideal_search_radius,
                                norm_prefilter, restart, refinement_target,
                                cap, max_relations, candidates_tried,
                                accepted_relations,
                                CertificationMode::unknown,
                                refinement_goal_reached,
                                reduced_lattice_cache);
                    })) {
            return false;
        }
    }

    if (candidates_tried < total_candidate_limit &&
        accepted_relations < max_relations) {
        const slong cap = phase_limit();
        if (!run_restartable_relation_phase(
                    order, 1, cap, max_relations,
                    candidates_tried, accepted_relations,
                    refinement_goal_reached,
                    [&](slong restart) noexcept {
                        return collect_product_ideal_lattice_relations(
                                context, order, ideal_search_radius,
                                norm_prefilter, restart, refinement_target,
                                cap, max_relations, candidates_tried,
                                accepted_relations,
                                CertificationMode::unknown,
                                refinement_goal_reached,
                                reduced_lattice_cache);
                    })) {
            return false;
        }
    }

    return true;
}

}  // namespace detail::relation_search

using namespace detail::relation_search;

bool ClassGroupContext::run_native_experimental_relation_route_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool emit_norm_prefilter_profile_event) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;

    const bool run_relation_completion_table =
            class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_->audit.policy.relations ==
                    detail::NativeRelationStrategy::relation_completion_table;
    if (run_relation_completion_table) {
        return run_relation_production_route_(
                order, factor_base_bound, options,
                emit_norm_prefilter_profile_event);
    }

    const bool run_relation_completion_prepass =
            class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_->audit.policy.relations ==
                    detail::NativeRelationStrategy::relation_completion_prepass;
    if (run_relation_completion_prepass) {
        return run_relation_production_prepass_route_(
                order, factor_base_bound, options,
                emit_norm_prefilter_profile_event);
    }

    relation_kernel_units_target_ = options.target_relation_kernel_units;
    configure_partial_relations_(options);
    if (private_storage_ == nullptr) {
        return false;
    }

    NormPrefilter norm_prefilter;
    if (emit_norm_prefilter_profile_event) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.build_norm_prefilter");
    }
    if (!build_norm_prefilter(norm_prefilter, base_, factor_base_bound,
                              private_storage_->use_partial_relations)) {
        return false;
    }

    slong candidates_tried = 0;
    slong accepted_relations = relation_count();
    bool goal_reached = has_presentation() &&
            relation_kernel_unit_count() >= options.target_relation_kernel_units;
    bool partial_throttle_exit = false;
    bool pivot_random_ran = false;
    const bool run_early_pivot_random =
            order.degree() > 2 && base_.length() >= 3;

    flint::FmpzMat coordinates(1, order.degree());
    // The order basis is fixed for the entire coordinate box.  Retain the
    // exact coordinate-to-element map just as the reference T2 path does.
    OrderCoordinateElementConversion coordinate_conversion;
    bool search_ok = true;
    if (search_ok && !goal_reached && !partial_throttle_exit) {
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_phase");
        search_ok = coordinate_conversion.reset(order) &&
                run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    const ulong seed = relation_search_phase_seed(
                            *this, base_.length(), kCoordsPhase, restart);
                    return collect_coordinate_relations(
                            *this, order, options.coordinate_search_radius,
                            &norm_prefilter,
                            &coordinate_conversion,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            coordinates, seed, 0, candidates_tried,
                            accepted_relations,
                            options.requested_certification, goal_reached,
                            partial_throttle_exit);
                },
                &partial_throttle_exit);
    }
    if (search_ok && !goal_reached && run_early_pivot_random) {
        pivot_random_ran = true;
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_random_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_random_phase");
        search_ok = run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    return collect_pivot_random_ideal_lattice_relations(
                            *this, order, options.ideal_search_radius,
                            &norm_prefilter, restart,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            candidates_tried, accepted_relations,
                            options.requested_certification, goal_reached);
                });
    }
    if (search_ok && !goal_reached) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.random_product_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.random_product_phase");
        search_ok = run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    return collect_random_product_ideal_lattice_relations(
                            *this, order, options.ideal_search_radius,
                            &norm_prefilter, restart,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            candidates_tried, accepted_relations,
                            options.requested_certification, goal_reached);
                });
    }
    if (search_ok && !goal_reached && !pivot_random_ran) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_random_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_random_phase");
        search_ok = run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    return collect_pivot_random_ideal_lattice_relations(
                            *this, order, options.ideal_search_radius,
                            &norm_prefilter, restart,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            candidates_tried, accepted_relations,
                            options.requested_certification, goal_reached);
                });
    }
    if (search_ok && !goal_reached) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.prime_ideal_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.prime_ideal_phase");
        search_ok = run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    return collect_prime_ideal_lattice_relations(
                            *this, order, options.ideal_search_radius,
                            &norm_prefilter, restart,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            candidates_tried, accepted_relations,
                            options.requested_certification, goal_reached);
                });
    }
    if (search_ok && !goal_reached) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.product_ideal_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.product_ideal_phase");
        search_ok = run_restartable_relation_phase(
                order, 1, options.max_candidates,
                options.max_relations, candidates_tried, accepted_relations,
                goal_reached,
                [&](slong restart) noexcept {
                    return collect_product_ideal_lattice_relations(
                            *this, order, options.ideal_search_radius,
                            &norm_prefilter, restart,
                            options.target_relation_kernel_units,
                            options.max_candidates, options.max_relations,
                            candidates_tried, accepted_relations,
                            options.requested_certification, goal_reached);
                });
    }
    bool auto_saturation_changed = false;
    if (search_ok && !try_auto_relation_saturation_(
                             order, factor_base_bound, options,
                             &auto_saturation_changed)) {
        search_ok = false;
    }
    if (search_ok) {
        accepted_relations = relation_count();
        if (!goal_reached && has_presentation() &&
            relation_kernel_unit_count() >=
                    options.target_relation_kernel_units) {
            goal_reached =
                    options.requested_certification !=
                            CertificationMode::proven ||
                    try_certify_quadratic(
                            options.requested_certification);
        }
    }
    const bool target_blocks_default_refinement =
            options.post_finite_refinement_phase_budget <= 0 &&
            options.target_relation_kernel_units > 0 &&
            relation_kernel_unit_count() >=
                    options.target_relation_kernel_units;
    bool post_finite_refinement_ran = false;
    if (search_ok && goal_reached &&
        options.requested_certification == CertificationMode::unknown &&
        !target_blocks_default_refinement) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.post_finite_refinement_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.post_finite_refinement_phase");
        post_finite_refinement_ran = true;
        search_ok = collect_post_finite_refinement_relations(
                *this, order, options.ideal_search_radius,
                &norm_prefilter, options.post_finite_refinement_phase_budget,
                options.max_candidates, options.max_relations,
                candidates_tried, accepted_relations);
    }
    if (search_ok &&
        (auto_saturation_changed || post_finite_refinement_ran) &&
        !try_auto_relation_saturation_(
                order, factor_base_bound, options)) {
        search_ok = false;
    }
    if (search_ok) {
        accepted_relations = relation_count();
        if (!goal_reached && has_presentation() &&
            relation_kernel_unit_count() >=
                    options.target_relation_kernel_units) {
            goal_reached =
                    options.requested_certification !=
                            CertificationMode::proven ||
                    try_certify_quadratic(
                            options.requested_certification);
        }
    }

    return search_ok && goal_reached;
}

}  // namespace silex
