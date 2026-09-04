#include <silex/class_group.hpp>
#include <silex/ideal.hpp>
#include <silex/ideal_factorization.hpp>

#include "class_group_internal.hpp"
#include "lll_relation_search_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"
#include "../order/order_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"

#include <flint/fmpz.h>

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {

constexpr slong kSmallLllInitialPrecision = 200;
constexpr ulong kLllRoutePhase = UWORD(8);

bool small_lll_selection_bound(flint::Fmpz& out,
                                     const Ideal& ideal,
                                     const Order& order,
                                     const IdealLatticeLllData& context)
        noexcept {
    const slong degree = order.degree();
    if (degree <= 0 ||
        flint::fmpz_sgn(flint::FmpzConstRef(context.gram_denominator)) <= 0) {
        return false;
    }

    flint::Fmpz discriminant;
    flint::Fmpz ideal_norm;
    flint::Fmpz product;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        !ideal.norm(flint::FmpzRef(ideal_norm)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(ideal_norm)) <= 0) {
        return false;
    }

    flint::fmpz_abs(flint::FmpzRef(discriminant),
                    flint::FmpzConstRef(discriminant));
    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(ideal_norm),
                    flint::FmpzConstRef(ideal_norm));
    flint::fmpz_mul(flint::FmpzRef(product), flint::FmpzConstRef(product),
                    flint::FmpzConstRef(discriminant));
    fmpz_root(out.raw(), product.raw(), static_cast<ulong>(degree));
    flint::fmpz_mul(flint::FmpzRef(out), flint::FmpzConstRef(out),
                    flint::FmpzConstRef(context.gram_denominator));
    return flint::fmpz_sgn(flint::FmpzConstRef(out)) > 0;
}

bool select_small_lll_rows(
        std::vector<slong>& selected,
        const IdealLatticeLllData& context,
        const Ideal& ideal,
        const Order& order) noexcept {
    selected.clear();
    const slong degree = order.degree();
    if (degree <= 0 ||
        flint::fmpz_mat_nrows(context.scaled_gram) != degree ||
        flint::fmpz_mat_ncols(context.scaled_gram) != degree) {
        return false;
    }

    slong target = degree / 4;
    if (target < 2) {
        target = degree;
    }
    if (target <= 0 || target > degree) {
        return false;
    }

    flint::Fmpz bound;
    if (!small_lll_selection_bound(bound, ideal, order, context)) {
        return false;
    }

    auto select_with_bound = [&]() noexcept {
        selected.clear();
        for (slong row = 0; row < degree; ++row) {
            if (flint::fmpz_cmp(
                        flint::fmpz_mat_entry(
                                flint::FmpzMatConstRef(context.scaled_gram),
                                row, row),
                        flint::FmpzConstRef(bound)) < 0) {
                selected.push_back(row);
            }
        }
    };

    select_with_bound();
    while (static_cast<slong>(selected.size()) < target) {
        flint::fmpz_mul_2exp(flint::FmpzRef(bound),
                             flint::FmpzConstRef(bound), UWORD(1));
        select_with_bound();
    }
    return true;
}

bool set_small_lll_combination(flint::FmpzMat& out,
                                     flint::FmpzMatConstRef basis,
                                     slong left_row,
                                     slong right_row,
                                     slong sign) noexcept {
    const slong columns = flint::fmpz_mat_ncols(basis);
    if (flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != columns ||
        left_row < 0 || left_row >= flint::fmpz_mat_nrows(basis) ||
        right_row < 0 || right_row >= flint::fmpz_mat_nrows(basis) ||
        (sign != 1 && sign != -1)) {
        return false;
    }

    for (slong column = 0; column < columns; ++column) {
        flint::FmpzRef entry = flint::fmpz_mat_entry(out, 0, column);
        flint::fmpz_set(entry, flint::fmpz_mat_entry(basis, left_row, column));
        if (sign > 0) {
            flint::fmpz_add(entry, flint::FmpzConstRef(entry.raw()),
                            flint::fmpz_mat_entry(basis, right_row, column));
        } else {
            flint::fmpz_sub(entry, flint::FmpzConstRef(entry.raw()),
                            flint::fmpz_mat_entry(basis, right_row, column));
        }
    }
    return true;
}

bool set_small_lll_next_candidate(
        flint::FmpzMat& out,
        slong& count_after_next,
        const IdealLatticeLllData& context,
        const std::vector<slong>& selected) noexcept {
    const slong selected_count = static_cast<slong>(selected.size());
    if (selected_count <= 0) {
        return false;
    }

    if (count_after_next < selected_count) {
        if (!copy_ideal_basis_row_coordinates(
                    out, flint::FmpzMatConstRef(context.basis),
                    selected[static_cast<std::size_t>(count_after_next)])) {
            return false;
        }
        ++count_after_next;
        return true;
    }

    if (count_after_next > selected_count * selected_count) {
        return false;
    }

    slong c = count_after_next - selected_count;
    slong sign = 1;
    if (c > selected_count * (selected_count - 1) / 2) {
        c -= selected_count * (selected_count - 1) / 2;
        sign = -1;
    }

    slong left = 1;
    while (c + left - 1 >= selected_count) {
        c -= selected_count - left;
        ++left;
    }
    const slong right = c + left - 1;
    if (left < 1 || left > selected_count ||
        right < 0 || right >= selected_count) {
        return false;
    }

    if (!set_small_lll_combination(
                out, flint::FmpzMatConstRef(context.basis),
                selected[static_cast<std::size_t>(left - 1)],
                selected[static_cast<std::size_t>(right)], sign)) {
        return false;
    }
    ++count_after_next;
    return true;
}

bool collect_order_small_lll_relations(
        ClassGroupContext& context,
        const Order& order,
        const Ideal& order_ideal,
        NormPrefilter* norm_prefilter,
        slong target_relation_kernel_units,
        const detail::ClassGroupRelationOptions& options,
        slong& candidates_tried,
        slong& accepted_relations,
        bool& goal_reached,
        bool& partial_throttle_exit,
        ClassGroupRelationSource source,
        slong max_good_relations) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.small_lll.collect_order");
    if (goal_reached || partial_throttle_exit ||
        candidates_tried >= options.max_candidates ||
        accepted_relations >= options.max_relations) {
        return true;
    }

    IdealLatticeLllData t2_context;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.small_lll.build_t2_context");
        if (!build_ideal_lattice_lll_data(
                    t2_context, order_ideal, kSmallLllInitialPrecision)) {
            return false;
        }
    }

    std::vector<slong> selected;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.small_lll.select_rows");
        if (!select_small_lll_rows(selected, t2_context, order_ideal,
                                         order)) {
            return false;
        }
    }


    flint::FmpzMat coordinates(1, order.degree());
    slong count = 0;
    slong good_relations = 0;
    const slong processed_limit =
            static_cast<slong>(selected.size() * selected.size());
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.small_lll.candidate_loop");
    while (set_small_lll_next_candidate(coordinates, count,
                                              t2_context, selected)) {
        // reference single_env checks I.cnt > length(I.b)^2 after fetching the
        // next candidate; mirror that boundary before exact relation testing.
        if (count > processed_limit) {
            break;
        }
        const slong relation_count_before = context.relation_count();
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.small_lll.candidate_try");
            if (!try_coordinate_candidate(
                        context, order, coordinates, norm_prefilter,
                        source,
                        target_relation_kernel_units,
                        options.max_candidates, options.max_relations,
                        candidates_tried, accepted_relations,
                        options.requested_certification, goal_reached,
                        partial_throttle_exit)) {
                return false;
            }
        }
        if (context.relation_count() > relation_count_before) {
            good_relations += context.relation_count() - relation_count_before;
            if (max_good_relations > -1 &&
                max_good_relations < good_relations) {
                break;
            }
        }
        if (goal_reached || partial_throttle_exit ||
            candidates_tried >= options.max_candidates ||
            accepted_relations >= options.max_relations) {
            break;
        }
    }

    return true;
}

ulong random_next(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) +
            UWORD(1442695040888963407);
    return state;
}

slong random_index(ulong& state, slong length) noexcept {
    if (length <= 0) {
        return -1;
    }
    return static_cast<slong>(
            random_next(state) % static_cast<ulong>(length));
}

bool indices_contain(const std::vector<slong>& indices,
                     slong value) noexcept {
    for (slong index : indices) {
        if (index == value) {
            return true;
        }
    }
    return false;
}

void append_unique_index(std::vector<slong>& indices, slong value) {
    if (!indices_contain(indices, value)) {
        indices.push_back(value);
    }
}

void normalize_index_set(std::vector<slong>& indices) {
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()),
                  indices.end());
}

bool same_index_set(std::vector<slong> left,
                    std::vector<slong> right) {
    normalize_index_set(left);
    normalize_index_set(right);
    return left == right;
}

bool supplement_pivots_if_empty(std::vector<slong>& pivots,
                                      const FactorBase& base,
                                      ulong& random_state) noexcept {
    if (!pivots.empty()) {
        return true;
    }
    if (base.length() <= 0) {
        return false;
    }
    for (slong i = 0; i < 5; ++i) {
        const slong index = random_index(random_state, base.length());
        if (index < 0) {
            return false;
        }
        append_unique_index(pivots, index);
    }
    normalize_index_set(pivots);
    return true;
}

bool build_random_base_indices(std::vector<slong>& out,
                                     const FactorBase& base,
                                     slong rand_exp,
                                     ulong& random_state) noexcept {
    out.clear();
    const slong length = base.length();
    if (base.parent() == nullptr || length <= 0 || rand_exp < 1) {
        return false;
    }

    const slong half_start = length / 2 + 1;
    const slong tail_start = length - 10 * (1 + rand_exp / 3);
    slong start = half_start > tail_start ? half_start : tail_start;
    if (start < 1) {
        start = 1;
    }
    slong target = length - start + 1;
    if (target < 1) {
        target = 1;
    }
    out.reserve(static_cast<std::size_t>(target));
    out.push_back(length - 1);

    PrimeIdeal prime(*base.parent());
    if (!prime.is_defined()) {
        return false;
    }
    std::vector<flint::Fmpz> minimums;
    flint::Fmpz rational_prime;
    if (!base.prime(prime, length - 1) ||
        !prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }
    minimums.push_back(std::move(rational_prime));

    for (slong i = length - 1;
         i > 0 && static_cast<slong>(out.size()) < target;) {
        --i;
        flint::Fmpz p;
        if (!base.prime(prime, i) ||
            !prime.rational_prime(flint::FmpzRef(p))) {
            return false;
        }
        bool repeated_minimum = false;
        for (const flint::Fmpz& seen : minimums) {
            if (flint::fmpz_equal(flint::FmpzConstRef(seen),
                                  flint::FmpzConstRef(p))) {
                repeated_minimum = true;
                break;
            }
        }
        if (prime.residue_degree() > 1 || repeated_minimum) {
            continue;
        }
        out.push_back(i);
        minimums.push_back(std::move(p));
    }

    while (static_cast<slong>(out.size()) < target) {
        const slong index = random_index(random_state, length);
        if (index < 0) {
            return false;
        }
        if (!indices_contain(out, index)) {
            out.push_back(index);
        }
    }
    return true;
}

class RandomIdealSearchEnvironment {
public:
    bool reset(const Order& order,
               const FactorBase& factor_base,
               const std::vector<slong>& indices,
               flint::FmpzConstRef lower_bound,
               flint::FmpzConstRef upper_bound,
               ulong seed) noexcept {
        clear();
        if (indices.empty() || factor_base.parent() == nullptr ||
            !factor_base.parent()->has_same_data(order) ||
            flint::fmpz_sgn(lower_bound) < 0 ||
            flint::fmpz_sgn(upper_bound) < 0) {
            return false;
        }

        order_ = &order;
        random_state_ = seed;
        flint::fmpz_set(flint::FmpzRef(lower_bound_), lower_bound);
        flint::fmpz_set(flint::FmpzRef(upper_bound_), upper_bound);
        if (!random_ideal_.define(order) || !random_ideal_.one()) {
            return false;
        }

        PrimeIdeal prime(order);
        if (!prime.is_defined()) {
            return false;
        }
        base_.reserve(indices.size());
        base_fractional_.reserve(indices.size());
        inverse_base_.reserve(indices.size());
        base_indices_.reserve(indices.size());
        exponents_.assign(indices.size(), 0);
        for (slong index : indices) {
            if (index < 0 || index >= factor_base.length()) {
                return false;
            }
            base_.emplace_back(order);
            base_fractional_.emplace_back(order);
            inverse_base_.emplace_back(order);
            if (!base_.back().is_defined() ||
                !base_fractional_.back().is_defined() ||
                !inverse_base_.back().is_defined() ||
                !factor_base.prime(prime, index) ||
                !prime.get_ideal(base_.back()) ||
                !base_fractional_.back().set_integral(base_.back()) ||
                !inverse_base_.back().invert(base_fractional_.back())) {
                return false;
            }
            base_indices_.push_back(index);
        }

        flint::Fmpz norm;
        if (!random_ideal_.norm(flint::FmpzRef(norm))) {
            return false;
        }
        while (flint::fmpz_cmp(flint::FmpzConstRef(norm),
                               lower_bound) <= 0) {
            const slong index = random_index(
                    random_state_, static_cast<slong>(base_.size()));
            if (index < 0 || !multiply_random_by_base(index)) {
                return false;
            }
            ++exponents_[static_cast<std::size_t>(index)];
            if (!random_ideal_.norm(flint::FmpzRef(norm))) {
                return false;
            }
        }
        return true;
    }

    bool reset_from_state(const Order& order,
                          const FactorBase& factor_base,
                          detail::RandomIdealSearchState& state)
            noexcept {
        clear();
        if (!state.defined || state.base_indices.empty() ||
            state.base_indices.size() != state.exponents.size() ||
            factor_base.parent() == nullptr ||
            !factor_base.parent()->has_same_data(order) ||
            flint::fmpz_sgn(flint::FmpzConstRef(state.lower_bound)) < 0 ||
            flint::fmpz_sgn(flint::FmpzConstRef(state.upper_bound)) < 0) {
            return false;
        }

        const bool has_cached_environment =
                state.random_ideal.is_defined() &&
                state.base_ideals.size() == state.base_indices.size() &&
                state.base_fractional_ideals.size() ==
                        state.base_indices.size() &&
                state.inverse_base_ideals.size() ==
                        state.base_indices.size();

        order_ = &order;
        random_state_ = state.random_state;
        flint::fmpz_set(flint::FmpzRef(lower_bound_),
                        flint::FmpzConstRef(state.lower_bound));
        flint::fmpz_set(flint::FmpzRef(upper_bound_),
                        flint::FmpzConstRef(state.upper_bound));
        if (has_cached_environment) {
            base_indices_.swap(state.base_indices);
            exponents_.swap(state.exponents);
            base_.swap(state.base_ideals);
            base_fractional_.swap(state.base_fractional_ideals);
            inverse_base_.swap(state.inverse_base_ideals);
            random_ideal_.swap(state.random_ideal);
            state.defined = false;
            if (base_.empty() || base_indices_.empty() ||
                base_.size() != base_indices_.size() ||
                base_fractional_.size() != base_indices_.size() ||
                inverse_base_.size() != base_indices_.size() ||
                exponents_.size() != base_indices_.size() ||
                random_ideal_.parent() == nullptr ||
                !random_ideal_.parent()->has_same_data(order)) {
                clear();
                return false;
            }
            for (slong i = 0; i < static_cast<slong>(base_indices_.size());
                 ++i) {
                if (base_indices_[static_cast<std::size_t>(i)] < 0 ||
                    base_indices_[static_cast<std::size_t>(i)] >=
                            factor_base.length() ||
                    exponents_[static_cast<std::size_t>(i)] < 0) {
                    clear();
                    return false;
                }
            }
            return true;
        }

        if (!random_ideal_.define(order) || !random_ideal_.one()) {
            return false;
        }

        PrimeIdeal prime(order);
        if (!prime.is_defined()) {
            return false;
        }
        base_.reserve(state.base_indices.size());
        base_fractional_.reserve(state.base_indices.size());
        inverse_base_.reserve(state.base_indices.size());
        base_indices_.reserve(state.base_indices.size());
        exponents_ = state.exponents;
        for (slong index : state.base_indices) {
            if (index < 0 || index >= factor_base.length()) {
                return false;
            }
            base_.emplace_back(order);
            base_fractional_.emplace_back(order);
            inverse_base_.emplace_back(order);
            if (!base_.back().is_defined() ||
                !base_fractional_.back().is_defined() ||
                !inverse_base_.back().is_defined() ||
                !factor_base.prime(prime, index) ||
                !prime.get_ideal(base_.back()) ||
                !base_fractional_.back().set_integral(base_.back()) ||
                !inverse_base_.back().invert(base_fractional_.back())) {
                return false;
            }
            base_indices_.push_back(index);
        }

        if (state.random_ideal.is_defined()) {
            if (!random_ideal_.set(state.random_ideal)) {
                return false;
            }
        } else {
            for (slong i = 0; i < static_cast<slong>(exponents_.size());
                 ++i) {
                const slong exponent =
                        exponents_[static_cast<std::size_t>(i)];
                if (exponent < 0) {
                    return false;
                }
                for (slong k = 0; k < exponent; ++k) {
                    if (!multiply_random_by_base(i)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool save_to_state(detail::RandomIdealSearchState& state)
            noexcept {
        if (order_ == nullptr || base_indices_.empty() ||
            base_indices_.size() != exponents_.size()) {
            return false;
        }
        state.defined = false;
        state.base_indices.clear();
        state.exponents.clear();
        state.base_ideals.clear();
        state.base_fractional_ideals.clear();
        state.inverse_base_ideals.clear();
        state.random_ideal.clear();
        state.base_indices.swap(base_indices_);
        state.exponents.swap(exponents_);
        state.base_ideals.swap(base_);
        state.base_fractional_ideals.swap(base_fractional_);
        state.inverse_base_ideals.swap(inverse_base_);
        state.random_ideal.swap(random_ideal_);
        flint::fmpz_set(flint::FmpzRef(state.lower_bound),
                        flint::FmpzConstRef(lower_bound_));
        flint::fmpz_set(flint::FmpzRef(state.upper_bound),
                        flint::FmpzConstRef(upper_bound_));
        state.random_state = random_state_;
        if (!state.random_ideal.is_defined() ||
            state.random_ideal.parent() == nullptr ||
            !state.random_ideal.parent()->has_same_data(*order_)) {
            return false;
        }
        state.defined = true;
        return true;
    }

    bool next(Ideal& out) noexcept {
        if (order_ == nullptr || !out.is_defined() ||
            !out.parent()->has_same_data(*order_)) {
            return false;
        }

        flint::Fmpz norm;
        if (!random_ideal_.norm(flint::FmpzRef(norm))) {
            return false;
        }

        slong delta = 0;
        if (flint::fmpz_cmp(flint::FmpzConstRef(norm),
                            flint::FmpzConstRef(upper_bound_)) >= 0) {
            delta = -1;
        } else if (flint::fmpz_cmp(flint::FmpzConstRef(norm),
                                   flint::FmpzConstRef(lower_bound_)) <= 0) {
            delta = 1;
        } else {
            delta = (random_next(random_state_) & UWORD(1)) == 0
                    ? WORD(-1)
                    : WORD(1);
        }

        slong index = -1;
        if (delta > 0) {
            index = random_index(random_state_,
                                       static_cast<slong>(base_.size()));
        } else {
            slong nonzero = 0;
            for (slong exponent : exponents_) {
                if (exponent != 0) {
                    ++nonzero;
                }
            }
            if (nonzero == 0) {
                return out.set(random_ideal_);
            }
            slong chosen = random_index(random_state_, nonzero);
            if (chosen < 0) {
                return false;
            }
            for (slong i = 0; i < static_cast<slong>(exponents_.size());
                 ++i) {
                if (exponents_[static_cast<std::size_t>(i)] == 0) {
                    continue;
                }
                if (chosen == 0) {
                    index = i;
                    break;
                }
                --chosen;
            }
        }
        if (index < 0) {
            return false;
        }

        exponents_[static_cast<std::size_t>(index)] += delta;
        if (delta > 0) {
            if (!multiply_random_by_base(index)) {
                return false;
            }
        } else if (!multiply_random_by_inverse_base(index)) {
            return false;
        }

        return out.set(random_ideal_);
    }

    bool extend_by_factor(flint::FmpzConstRef factor) noexcept {
        if (order_ == nullptr || flint::fmpz_sgn(factor) <= 0) {
            return false;
        }

        flint::fmpz_mul(flint::FmpzRef(lower_bound_),
                        flint::FmpzConstRef(lower_bound_), factor);
        flint::fmpz_mul(flint::FmpzRef(upper_bound_),
                        flint::FmpzConstRef(lower_bound_), factor);

        flint::Fmpz norm;
        if (!random_ideal_.norm(flint::FmpzRef(norm))) {
            return false;
        }
        while (flint::fmpz_cmp(flint::FmpzConstRef(norm),
                               flint::FmpzConstRef(lower_bound_)) < 0) {
            const slong index = random_index(
                    random_state_, static_cast<slong>(base_.size()));
            if (index < 0 || !multiply_random_by_base(index)) {
                return false;
            }
            ++exponents_[static_cast<std::size_t>(index)];
            if (!random_ideal_.norm(flint::FmpzRef(norm))) {
                return false;
            }
        }
        return true;
    }

    slong base_length() const noexcept {
        return static_cast<slong>(base_.size());
    }

private:
    void clear() noexcept {
        order_ = nullptr;
        base_.clear();
        base_fractional_.clear();
        inverse_base_.clear();
        base_indices_.clear();
        exponents_.clear();
        random_ideal_.clear();
        flint::fmpz_zero(flint::FmpzRef(lower_bound_));
        flint::fmpz_zero(flint::FmpzRef(upper_bound_));
        random_state_ = 0;
    }

    bool multiply_random_by_base(slong index) noexcept {
        if (order_ == nullptr || index < 0 ||
            index >= static_cast<slong>(base_.size())) {
            return false;
        }
        Ideal product(*order_);
        if (!product.is_defined() ||
            !product.multiply(random_ideal_,
                              base_[static_cast<std::size_t>(index)])) {
            return false;
        }
        random_ideal_.swap(product);
        return true;
    }

    bool multiply_random_by_inverse_base(slong index) noexcept {
        if (order_ == nullptr || index < 0 ||
            index >= static_cast<slong>(inverse_base_.size())) {
            return false;
        }
        FractionalIdeal random_fractional(*order_);
        FractionalIdeal product(*order_);
        Ideal numerator(*order_);
        flint::Fmpz den;
        if (!random_fractional.is_defined() || !product.is_defined() ||
            !numerator.is_defined() ||
            !random_fractional.set_integral(random_ideal_) ||
            !product.multiply(
                    random_fractional,
                    inverse_base_[static_cast<std::size_t>(index)]) ||
            !product.get_integral_den(numerator, flint::FmpzRef(den))) {
            return false;
        }
        random_ideal_.swap(numerator);
        return true;
    }

    const Order* order_ = nullptr;
    std::vector<Ideal> base_;
    std::vector<FractionalIdeal> base_fractional_;
    std::vector<FractionalIdeal> inverse_base_;
    std::vector<slong> base_indices_;
    std::vector<slong> exponents_;
    Ideal random_ideal_;
    flint::Fmpz lower_bound_;
    flint::Fmpz upper_bound_;
    ulong random_state_ = 0;
};

bool build_random_factor_base_product(
        Ideal& out,
        flint::FmpzMatRef factor_base_row,
        const FactorBase& base,
        slong sequence_index,
        ulong seed,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*order) || base.length() <= 0 ||
        sequence_index < 0 || flint::fmpz_mat_nrows(factor_base_row) != 1 ||
        flint::fmpz_mat_ncols(factor_base_row) != base.length()) {
        return false;
    }

    ulong random_state = seed;
    std::vector<slong> random_base_indices;
    flint::Fmpz discriminant;
    flint::Fmpz lower_bound;
    if (!build_random_base_indices(
                random_base_indices, base, WORD(1), random_state) ||
        !order->discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(discriminant),
                    flint::FmpzConstRef(discriminant));
    fmpz_sqrt(lower_bound.raw(), discriminant.raw());

    RandomIdealSearchEnvironment environment;
    Ideal candidate(*order);
    if (!candidate.is_defined() ||
        !environment.reset(*order, base, random_base_indices,
                           flint::FmpzConstRef(lower_bound),
                           flint::FmpzConstRef(discriminant), random_state)) {
        return false;
    }
    for (slong i = 0; i <= sequence_index; ++i) {
        if (!environment.next(candidate)) {
            return false;
        }
    }

    flint::FmpzMat candidate_row(1, base.length());
    if (!ideal_factor_over_base(flint::FmpzMatRef(candidate_row), candidate,
                                base, diagnostics)) {
        return false;
    }
    out.swap(candidate);
    flint::fmpz_mat_set(factor_base_row,
                        flint::FmpzMatConstRef(candidate_row));
    return true;
}

bool factor_base_prime_power_ideal(Ideal& out,
                                   const FactorBase& base,
                                   PrimeIdeal& prime,
                                   Ideal& factor,
                                   Ideal& product,
                                   slong index,
                                   slong exponent) noexcept {
    if (index < 0 || index >= base.length() || exponent < 0) {
        return false;
    }
    if (exponent == 0) {
        return out.one();
    }
    if (!base.prime(prime, index) || !prime.get_ideal(factor) ||
        !out.set(factor)) {
        return false;
    }
    for (slong i = 1; i < exponent; ++i) {
        if (!product.multiply(out, factor)) {
            return false;
        }
        out.swap(product);
    }
    return true;
}

bool index_is_coprime_to_h(slong exponent,
                                 flint::FmpzConstRef h) noexcept {
    if (exponent <= 0 || flint::fmpz_sgn(h) <= 0) {
        return true;
    }
    const ulong reduced = flint::fmpz_fdiv_ui(h, static_cast<ulong>(exponent));
    return std::gcd(static_cast<ulong>(exponent), reduced) == UWORD(1);
}

bool collect_pivot_lll_ideal(
        ClassGroupContext& context,
        const Order& order,
        const FactorBase& base,
        slong pivot_index,
        slong rand_exp,
        RandomIdealSearchEnvironment& random_env,
        NormPrefilter* norm_prefilter,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        slong& candidates_tried,
        slong& accepted_relations,
        bool& goal_reached,
        bool& partial_throttle_exit) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.pivot_lll_ideal");
    if (goal_reached || partial_throttle_exit ||
        candidates_tried >= options.max_candidates ||
        accepted_relations >= options.max_relations ||
        base.length() <= 0) {
        return true;
    }
    if (pivot_index < 0 || pivot_index >= base.length()) {
        return false;
    }

    PrimeIdeal prime(order);
    Ideal random_ideal(order);
    Ideal pivot_power(order);
    Ideal pivot_factor(order);
    Ideal pivot_product(order);
    Ideal search_ideal(order);
    if (!prime.is_defined() || !random_ideal.is_defined() ||
        !pivot_power.is_defined() || !pivot_factor.is_defined() ||
        !pivot_product.is_defined() || !search_ideal.is_defined()) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_lll_ideal.random_next");
        if (!random_env.next(random_ideal)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.pivot_lll_ideal.pivot_power");
        if (!factor_base_prime_power_ideal(
                    pivot_power, base, prime, pivot_factor, pivot_product,
                    pivot_index, rand_exp)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.pivot_lll_ideal.multiply");
        if (!search_ideal.multiply(random_ideal, pivot_power)) {
            return false;
        }
    }

    const slong max_good = context.relation_rank() > 0 &&
                    base.length() < 2 * context.relation_rank()
            ? WORD(1)
            : WORD(-1);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.pivot_lll_ideal.collect_small_lll");
        if (!collect_order_small_lll_relations(
                    context, order, search_ideal, norm_prefilter,
                    target_relation_kernel_units, options,
                    candidates_tried, accepted_relations, goal_reached,
                    partial_throttle_exit,
                    ClassGroupRelationSource::RandomProduct, max_good)) {
            return false;
        }
    }

    return true;
}


}  // namespace detail::relation_search

using namespace detail::relation_search;

bool ClassGroupContext::pivot_info_(flint::FmpzRef h,
                                          std::vector<slong>& pivots)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.pivot_info");
    pivots.clear();
    if (!has_factor_base() || h.raw() == nullptr) {
        return false;
    }

    const slong n = generator_count();
    flint::fmpz_zero(h);
    if (relation_rank() < n) {
        for (slong index = 0; index < n; ++index) {
            bool covered = false;
            if (!factor_base_prime_is_hnf_covered(covered, index)) {
                return false;
            }
            if (!covered) {
                pivots.push_back(index);
            }
        }
        return true;
    }

    if (!sync_row_module_checkpoint_() || !row_module_.full_rank_index(h)) {
        return false;
    }

    flint::FmpzMat hnf(n, n);
    if (!row_module_.get_hnf_rows(flint::FmpzMatRef(hnf))) {
        return false;
    }
    for (slong row = 0; row < n; ++row) {
        slong pivot_col = -1;
        for (slong col = 0; col < n; ++col) {
            if (!flint::fmpz_is_zero(
                        flint::fmpz_mat_entry(
                                flint::FmpzMatConstRef(hnf), row, col))) {
                pivot_col = col;
                break;
            }
        }
        if (pivot_col < 0) {
            return false;
        }
        if (!flint::fmpz_is_one(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(hnf), row, pivot_col))) {
            pivots.push_back(pivot_col);
        }
    }
    return true;
}


bool ClassGroupContext::run_lll_relation_route_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool emit_norm_prefilter_profile_event,
        bool accept_tentative_presentation,
        bool extension_only) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_route");
    if (private_storage_ == nullptr) {
        return false;
    }
    detail::ClassGroupContextStorage& route_storage = *private_storage_;
    struct GeneratorRelationPolicyGuard {
        ClassGroupContext& context;
        ClassGroupContext::DependentRelationPolicy previous;

        GeneratorRelationPolicyGuard(
                ClassGroupContext& context_value,
                ClassGroupContext::DependentRelationPolicy policy) noexcept
            : context(context_value),
              previous(context_value.generator_relation_policy_) {
            context.generator_relation_policy_ = policy;
        }

        ~GeneratorRelationPolicyGuard() noexcept {
            context.generator_relation_policy_ = previous;
        }
    } generator_relation_policy_guard(
            *this, DependentRelationPolicy::keep_nonduplicate);

    relation_kernel_units_target_ = options.target_relation_kernel_units;
    configure_partial_relations_(options);
    // reference's LLL route stops at `single_env` /
    // `class_group_new_relations_via_lll` source boundaries, not when the
    // native presentation target first becomes available.
    constexpr slong source_loop_target_units = WORD_MAX;

    NormPrefilter norm_prefilter;
    if (emit_norm_prefilter_profile_event) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.build_norm_prefilter");
    }
    {
        SILEX_PROFILE_SCOPE(
                active_diagnostics, DiagnosticsModule::class_group,
                "class_group.route.build_norm_prefilter");
        if (!build_norm_prefilter(norm_prefilter, base_, factor_base_bound,
                                  route_storage.use_partial_relations)) {
            return false;
        }
    }

    Ideal order_ideal(order);
    {
        SILEX_PROFILE_SCOPE(active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group.route.order_ideal_one");
        if (!order_ideal.is_defined() || !order_ideal.one()) {
            return false;
        }
    }

    slong candidates_tried = 0;
    slong accepted_relations = relation_count();
    bool goal_reached = !extension_only && has_presentation() &&
            relation_kernel_unit_count() >= options.target_relation_kernel_units;
    bool partial_throttle_exit = false;
    bool search_ok = true;
    if (!extension_only && !goal_reached) {
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.order_small_lll_phase");
        SILEX_PROFILE_SCOPE(active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group.route.order_small_lll");
        search_ok = collect_order_small_lll_relations(
                *this, order, order_ideal, &norm_prefilter,
                source_loop_target_units, options,
                candidates_tried, accepted_relations, goal_reached,
                partial_throttle_exit);
    }
    if (search_ok && (extension_only || !goal_reached) &&
        !partial_throttle_exit) {
        flint::Fmpz h;
        std::vector<slong> pivots;
        {
            SILEX_PROFILE_SCOPE(active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route.initial_pivot_info");
            if (!pivot_info_(flint::FmpzRef(h), pivots)) {
                return false;
            }
        }
        if (!extension_only &&
            flint::fmpz_is_one(flint::FmpzConstRef(h))) {
            SILEX_PROFILE_SCOPE(active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route.initial_publish");
            goal_reached = publish_and_check_compute_goal(
                    *this, options.target_relation_kernel_units,
                    options.requested_certification);
        }
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.new_relations_via_lll");
        if (extension_only || !goal_reached) {
            ulong random_state = relation_search_phase_seed(
                    *this, base_.length(), kLllRoutePhase, WORD(0));
            {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.supplement_pivots");
                if (!supplement_pivots_if_empty(
                            pivots, base_, random_state)) {
                    return false;
                }
            }

            slong rand_exp = 1;
            {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.select_rand_exp");
                if (flint::fmpz_sgn(flint::FmpzConstRef(h)) > 0) {
                    ++rand_exp;
                    while (!index_is_coprime_to_h(
                            rand_exp, flint::FmpzConstRef(h))) {
                        ++rand_exp;
                    }
                }
            }

            std::vector<slong> random_base_indices;
            {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.random_base_indices");
                if (!build_random_base_indices(
                            random_base_indices, base_, rand_exp,
                            random_state)) {
                    return false;
                }
            }

            flint::Fmpz discriminant;
            flint::Fmpz lower_bound;
            {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.discriminant_bounds");
                if (!order.discriminant(flint::FmpzRef(discriminant))) {
                    return false;
                }
                flint::fmpz_abs(flint::FmpzRef(discriminant),
                                flint::FmpzConstRef(discriminant));
                fmpz_sqrt(lower_bound.raw(), discriminant.raw());
            }

            RandomIdealSearchEnvironment random_env;
            {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.random_env_setup");
                if (route_storage.random_ideal_search_state.defined) {
                    flint::Fmpz extension_factor;
                    flint::fmpz_set_ui(flint::FmpzRef(extension_factor), 2);
                    {
                        SILEX_PROFILE_SCOPE(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.route.random_env_reset_from_state");
                        if (!random_env.reset_from_state(
                                    order, base_,
                                    route_storage.random_ideal_search_state)) {
                            return false;
                        }
                    }
                    {
                      SILEX_PROFILE_SCOPE(
                          active_diagnostics, DiagnosticsModule::class_group,
                          "class_group.route.random_env_extend_factor");
                      if (!random_env.extend_by_factor(
                              flint::FmpzConstRef(extension_factor))) {
                        return false;
                      }
                    }
                } else {
                  SILEX_PROFILE_SCOPE(
                      active_diagnostics, DiagnosticsModule::class_group,
                      "class_group.route.random_env_initial_reset");
                  if (!random_env.reset(order, base_, random_base_indices,
                                        flint::FmpzConstRef(lower_bound),
                                        flint::FmpzConstRef(discriminant),
                                        random_state)) {
                    return false;
                  }
                }
            }


            auto route_budget_exhausted = [&]() noexcept {
              return (!extension_only && goal_reached) ||
                     partial_throttle_exit ||
                     candidates_tried >= options.max_candidates ||
                     accepted_relations >= options.max_relations;
            };

            bool source_loop_returned = false;
            {
              SILEX_PROFILE_SCOPE(active_diagnostics,
                                  DiagnosticsModule::class_group,
                                  "class_group.route.outer_loop");
              while (!route_budget_exhausted() && !source_loop_returned) {
                const slong relation_start = relation_count();
                while (!route_budget_exhausted() &&
                       relation_count() - relation_start < 2 &&
                       !source_loop_returned) {
                  std::vector<slong> pivot_order;
                  {
                    SILEX_PROFILE_SCOPE(active_diagnostics,
                                        DiagnosticsModule::class_group,
                                        "class_group.route.pivot_order");
                    pivot_order = pivots;
                    const bool sort_reverse =
                        10 * relation_rank() < 9 * base_.length();
                    std::sort(pivot_order.begin(), pivot_order.end());
                    if (sort_reverse) {
                      std::reverse(pivot_order.begin(), pivot_order.end());
                    }
                  }

                  for (slong pivot_index : pivot_order) {
                    if (route_budget_exhausted()) {
                      break;
                    }
                    {
                      SILEX_PROFILE_SCOPE(
                          active_diagnostics, DiagnosticsModule::class_group,
                          "class_group.route.collect_pivot");
                      if (!collect_pivot_lll_ideal(
                              *this, order, base_, pivot_index, rand_exp,
                              random_env, &norm_prefilter,
                              source_loop_target_units, options,
                              candidates_tried, accepted_relations,
                              goal_reached, partial_throttle_exit)) {
                        search_ok = false;
                        break;
                      }
                    }
                    if (route_budget_exhausted()) {
                      break;
                    }

                    if (flint::fmpz_is_zero(flint::FmpzConstRef(h)) &&
                        relation_rank() == base_.length()) {
                      std::vector<slong> refreshed_pivots;
                      {
                        SILEX_PROFILE_SCOPE(
                            active_diagnostics, DiagnosticsModule::class_group,
                            "class_group.route.rank_full_pivot_refresh");
                        if (!pivot_info_(flint::FmpzRef(h),
                                               refreshed_pivots)) {
                          return false;
                        }
                      }
                      if (flint::fmpz_sgn(flint::FmpzConstRef(h)) > 0) {
                        break;
                      }
                    }

                    const slong relation_delta =
                        relation_count() - relation_start;
                    if (flint::fmpz_sgn(flint::FmpzConstRef(h)) > 0 &&
                        relation_delta > 1) {
                      source_loop_returned = true;
                      break;
                    }
                    if (flint::fmpz_is_zero(flint::FmpzConstRef(h)) &&
                        relation_delta >
                            static_cast<slong>(pivots.size()) + 2) {
                      break;
                    }
                  }
                  if (!search_ok || pivot_order.empty()) {
                    break;
                  }
                }
                if (!search_ok || route_budget_exhausted() ||
                    source_loop_returned) {
                  break;
                }

                std::vector<slong> pivots_new;
                {
                  SILEX_PROFILE_SCOPE(
                      active_diagnostics, DiagnosticsModule::class_group,
                      "class_group.route.eval_pivot_info");
                  if (!pivot_info_(flint::FmpzRef(h), pivots_new)) {
                    return false;
                  }
                  normalize_index_set(pivots_new);
                }

                if (same_index_set(pivots_new, pivots)) {
                  {
                    SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.same_pivot_update");
                    if (flint::fmpz_sgn(flint::FmpzConstRef(h)) > 0) {
                      if (!random_env.extend_by_factor(
                              flint::FmpzConstRef(lower_bound))) {
                        return false;
                      }
                    }
                    ++rand_exp;
                    if (rand_exp > 13) {
                      rand_exp = 13;
                    }
                    if (flint::fmpz_sgn(flint::FmpzConstRef(h)) > 0) {
                      while (!index_is_coprime_to_h(
                          rand_exp, flint::FmpzConstRef(h))) {
                        ++rand_exp;
                      }
                    }
                  }
                }
                pivots.swap(pivots_new);
                if (flint::fmpz_is_one(flint::FmpzConstRef(h))) {
                  if (extension_only) {
                    source_loop_returned = true;
                    break;
                  }
                  {
                    SILEX_PROFILE_SCOPE(active_diagnostics,
                                        DiagnosticsModule::class_group,
                                        "class_group.route.loop_publish");
                    goal_reached = publish_and_check_compute_goal(
                        *this, options.target_relation_kernel_units,
                        options.requested_certification);
                  }
                  if (goal_reached) {
                    break;
                  }
                  {
                    SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.route.loop_supplement_pivots");
                    if (!supplement_pivots_if_empty(pivots, base_,
                                                          random_state)) {
                      return false;
                    }
                  }
                }
              }
            }
            {
                SILEX_PROFILE_SCOPE(active_diagnostics,
                                    DiagnosticsModule::class_group,
                                    "class_group.route.save_random_env");
                if (!random_env.save_to_state(
                            route_storage.random_ideal_search_state)) {
                    return false;
                }
            }
        }
    }

    if (search_ok) {
        accepted_relations = relation_count();
        if (!extension_only && !goal_reached && has_presentation() &&
            relation_kernel_unit_count() >=
                    options.target_relation_kernel_units) {
            goal_reached = true;
        }
    }

    SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
              LogLevel::detail,
              goal_reached
                      ? "route deterministic LLL slice reached goal"
                      : (accept_tentative_presentation && has_presentation())
                                ? "route deterministic LLL slice kept "
                                  "tentative presentation"
                                : "route deterministic LLL slice "
                                  "incomplete");
    return search_ok &&
           (goal_reached ||
            (accept_tentative_presentation && has_presentation()));
}


bool ClassGroupContext::compute_tentative_candidate_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.compute_tentative_candidate");
    if (!detail::order_has_parented_basis(order) || !order.is_maximal() ||
        options.coordinate_search_radius < 1 ||
        options.ideal_search_radius < 0 ||
        options.target_relation_kernel_units < 0 ||
        options.post_finite_refinement_phase_budget < 0 ||
        options.max_candidates < 0 || options.max_relations < 0 ||
        (options.relation_saturation_aux_prime_bound != 0 &&
         options.relation_saturation_aux_prime_bound < 2) ||
        options.relation_saturation_max_appends_per_ell < 0 ||
        !detail::valid_certification_request(options.requested_certification) ||
        (options.requested_certification == CertificationMode::proven &&
         !detail::order_has_quadratic_source_data(order))) {
        return false;
    }

    ClassGroupContext candidate(order);
    candidate.set_diagnostics(active_diagnostics);
    const bool strict_transaction =
            detail::uses_class_unit_kernel(
                    class_unit_transaction_context_);
    SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.build_factor_base");
    if (!candidate.is_defined() ||
        !candidate.build_search_factor_base_(factor_base_bound,
                                            strict_transaction)) {
        return false;
    }

    detail::ClassGroupRelationOptions local_options = options;
    local_options.diagnostics = active_diagnostics;
    // reference class_group_ctx builds a tentative class-group presentation first;
    // unit-kernel continuation belongs to the _class_unit_group loop.
    local_options.target_relation_kernel_units = 0;

    if (!candidate.run_lll_relation_route_(
                order, factor_base_bound, local_options, true, true) ||
        !candidate.has_presentation()) {
        return false;
    }

    swap(candidate);
    set_diagnostics(active_diagnostics);
    return true;
}

bool ClassGroupContext::extend_tentative_relations_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.extend_tentative_relations");
    if (!detail::order_has_parented_basis(order) || !order.is_maximal() ||
        !same_order_parent(parent(), &order) || !has_factor_base() ||
        !has_presentation() || options.coordinate_search_radius < 1 ||
        options.ideal_search_radius < 0 ||
        options.target_relation_kernel_units < 0 ||
        options.post_finite_refinement_phase_budget < 0 ||
        options.max_candidates < 0 || options.max_relations < 0 ||
        relation_count() >= options.max_relations ||
        (options.relation_saturation_aux_prime_bound != 0 &&
         options.relation_saturation_aux_prime_bound < 2) ||
        options.relation_saturation_max_appends_per_ell < 0 ||
        !detail::valid_certification_request(options.requested_certification)) {
        return false;
    }

    detail::ClassGroupRelationOptions local_options = options;
    local_options.diagnostics = active_diagnostics;
    return run_lll_relation_route_(order, factor_base_bound, local_options,
                                     false, true, true) &&
           has_presentation();
}

bool ClassGroupContext::extend_lll_relation_slice_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.extend_lll_relation_slice");
    if (!detail::order_has_parented_basis(order) || !order.is_maximal() ||
        !same_order_parent(parent(), &order) || !has_factor_base() ||
        !has_presentation() || options.coordinate_search_radius < 1 ||
        options.ideal_search_radius < 0 ||
        options.target_relation_kernel_units < 0 ||
        options.post_finite_refinement_phase_budget < 0 ||
        options.max_candidates < 0 || options.max_relations < 0 ||
        relation_count() >= options.max_relations ||
        (options.relation_saturation_aux_prime_bound != 0 &&
         options.relation_saturation_aux_prime_bound < 2) ||
        options.relation_saturation_max_appends_per_ell < 0 ||
        !detail::valid_certification_request(options.requested_certification)) {
        return false;
    }

    detail::ClassGroupRelationOptions local_options = options;
    local_options.diagnostics = active_diagnostics;
    // Source trace: reference `Clgp/Main_LLL.jl:
    // class_group_new_relations_via_lll` is invoked from
    // `_class_unit_group` independently of the relation-kernel target.  The
    // ported route uses `extension_only=true`, which keeps the same source
    // loop boundary and passes `WORD_MAX` to relation insertion so dependent
    // relation rows can become additional unit-kernel witnesses.
    return run_lll_relation_route_(order, factor_base_bound, local_options,
                                     false, true, true) &&
           has_presentation();
}


}  // namespace silex
