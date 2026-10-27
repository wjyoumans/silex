#include <silex/class_group.hpp>

#include <silex/flint/fmpz_lll.hpp>

#include "../lat/certified_lll_internal.hpp"
#include "class_group_internal.hpp"
#include "ideal_minkowski_embedding_internal.hpp"
#include "ideal_t2_enumeration_internal.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {

constexpr slong kIdealLatticeRoundingGuardBits = 64;
// First scale (in bits) at which the embedding rows are rounded for LLL; the
// scale doubles until the rounded rows have full rank, as the source rounds
// its embedding matrix (see "Class groups and order units" in
// docs/reference/algorithms_and_sources.rst).
constexpr slong kT2EnumerationInitialRoundingScale = 4;

bool fmpz_mat_single_row_is_primitive(flint::FmpzMatConstRef row) noexcept {
    if (flint::fmpz_mat_nrows(row) != 1) {
        return false;
    }

    flint::Fmpz content;
    flint::fmpz_zero(flint::FmpzRef(content));
    for (slong column = 0; column < flint::fmpz_mat_ncols(row); ++column) {
        flint::fmpz_gcd(flint::FmpzRef(content),
                        flint::FmpzConstRef(content),
                        flint::fmpz_mat_entry(row, 0, column));
        if (flint::fmpz_is_one(content)) {
            return true;
        }
    }

    return false;
}

static bool slong_from_double(slong& out, double value) noexcept {
    if (!std::isfinite(value) ||
        value < static_cast<double>(std::numeric_limits<slong>::min()) ||
        value > static_cast<double>(std::numeric_limits<slong>::max())) {
        return false;
    }
    out = static_cast<slong>(value);
    return true;
}

static bool quadratic_form_data_from_double_gram(std::vector<double>& out,
                                   std::vector<double>& gram,
                                   slong dimension) noexcept {
    if (dimension <= 0 ||
        static_cast<slong>(gram.size()) != dimension * dimension) {
        return false;
    }

    auto entry = [&](slong row, slong column) -> double& {
        return gram[static_cast<std::size_t>(row * dimension + column)];
    };

    for (slong i = 0; i < dimension - 1; ++i) {
        if (!std::isfinite(entry(i, i)) || entry(i, i) <= 0.0) {
            return false;
        }
        for (slong j = i + 1; j < dimension; ++j) {
            entry(j, i) = entry(i, j);
            entry(i, j) = entry(i, j) / entry(i, i);
            if (!std::isfinite(entry(i, j))) {
                return false;
            }
        }
        for (slong k = i + 1; k < dimension; ++k) {
            for (slong l = i + 1; l < dimension; ++l) {
                entry(k, l) -= entry(k, i) * entry(i, l);
                if (!std::isfinite(entry(k, l))) {
                    return false;
                }
            }
        }
    }
    for (slong j = 0; j < dimension - 1; ++j) {
        if (!std::isfinite(entry(j, j)) || entry(j, j) <= 0.0) {
            return false;
        }
        for (slong i = j + 1; i < dimension; ++i) {
            entry(i, j) = 0.0;
        }
    }
    if (!std::isfinite(entry(dimension - 1, dimension - 1)) ||
        entry(dimension - 1, dimension - 1) <= 0.0) {
        return false;
    }

    out = gram;
    return true;
}

static bool quadratic_form_data_from_arb_gram(std::vector<double>& out,
                                const flint::ArbMat& gram) noexcept;
static bool finite_quadratic_form_volume_bound(
        double& out,
        const std::vector<double>& gram_schmidt_diagonal,
        slong target_factor_attempts) noexcept;

static double quadratic_form_entry(const std::vector<double>& quadratic_form_data,
                     slong dimension,
                     slong row,
                     slong column) noexcept {
    return quadratic_form_data[static_cast<std::size_t>(row * dimension + column)];
}

static bool initial_bound_from_quadratic_form(double& out,
                                     const std::vector<double>& quadratic_form_data,
                                     slong dimension,
                                     slong target_factor_attempts) noexcept {
    if (dimension <= 0 ||
        static_cast<slong>(quadratic_form_data.size()) != dimension * dimension) {
        return false;
    }

    const double first_diagonal = quadratic_form_entry(quadratic_form_data, dimension, 0, 0);
    if (!std::isfinite(first_diagonal) || first_diagonal <= 0.0) {
        return false;
    }

    double second_vector_norm = first_diagonal;
    if (dimension > 1) {
        const double second_diagonal =
                quadratic_form_entry(quadratic_form_data, dimension, 1, 1);
        const double first_second_mu =
                quadratic_form_entry(quadratic_form_data, dimension, 0, 1);
        second_vector_norm =
                second_diagonal +
                first_diagonal * first_second_mu * first_second_mu;
    }
    if (!std::isfinite(second_vector_norm) || second_vector_norm <= 0.0) {
        return false;
    }

    std::vector<double> gram_schmidt_diagonal;
    gram_schmidt_diagonal.reserve(static_cast<std::size_t>(dimension));
    for (slong i = 0; i < dimension; ++i) {
        const double diagonal = quadratic_form_entry(quadratic_form_data, dimension, i, i);
        if (!std::isfinite(diagonal) || diagonal <= 0.0) {
            return false;
        }
        gram_schmidt_diagonal.push_back(diagonal);
    }

    double volume_bound = 0.0;
    if (!finite_quadratic_form_volume_bound(volume_bound,
                                        gram_schmidt_diagonal,
                                        target_factor_attempts)) {
        return false;
    }
    out = std::max(2.0 * second_vector_norm, volume_bound);
    return std::isfinite(out) && out > 0.0;
}

bool FiniteQuadraticFormEnumerationContext::reset(
        const std::vector<double>& quadratic_form_data, slong dimension) noexcept {
    if (dimension <= 0 ||
        static_cast<slong>(quadratic_form_data.size()) != dimension * dimension) {
        clear();
        return false;
    }
    for (slong i = 0; i < dimension; ++i) {
        const double diagonal = quadratic_form_entry(quadratic_form_data, dimension, i, i);
        if (!std::isfinite(diagonal) || diagonal <= 0.0) {
            clear();
            return false;
        }
        for (slong j = i + 1; j < dimension; ++j) {
            if (!std::isfinite(quadratic_form_entry(quadratic_form_data, dimension, i, j))) {
                clear();
                return false;
            }
        }
    }
    dimension_ = dimension;
    quadratic_form_data_ = quadratic_form_data;
    coefficients_.assign(static_cast<std::size_t>(dimension_), 0);
    increments_.assign(static_cast<std::size_t>(dimension_), 1);
    lengths_.assign(static_cast<std::size_t>(dimension_), 0.0);
    tails_.assign(static_cast<std::size_t>(dimension_), 0.0);
    started_ = false;
    exhausted_ = false;
    need_outer_step_ = false;
    element_steps_ = 0;
    element_step_limit_reached_ = false;
    return true;
}

bool FiniteQuadraticFormEnumerationContext::start(double bound,
                                              slong max_element_steps,
                                              bool skip_first_scalar) noexcept {
    if (dimension_ <= 0 || !std::isfinite(bound) || bound <= 0.0 ||
        max_element_steps <= 0) {
        return false;
    }
    bound_ = bound;
    max_element_steps_ = max_element_steps;
    skip_first_ = skip_first_scalar ? 1 : 0;
    std::fill(coefficients_.begin(), coefficients_.end(), 0);
    std::fill(increments_.begin(), increments_.end(), 1);
    std::fill(lengths_.begin(), lengths_.end(), 0.0);
    std::fill(tails_.begin(), tails_.end(), 0.0);
    level_ = dimension_ - 1;
    started_ = true;
    exhausted_ = false;
    need_outer_step_ = false;
    element_steps_ = 0;
    element_step_limit_reached_ = false;
    return true;
}

bool FiniteQuadraticFormEnumerationContext::next() noexcept {
    if (!started_ || exhausted_) {
        return false;
    }
    if (need_outer_step_) {
        if (!step(level_)) {
            exhausted_ = true;
            return false;
        }
        need_outer_step_ = false;
    }

    while (true) {
        do {
            bool skip_coordinate_test = false;
            if (level_ > 0) {
                const slong next_level = level_ - 1;
                tails_[static_cast<std::size_t>(next_level)] = 0.0;
                for (slong j = level_; j < dimension_; ++j) {
                    tails_[static_cast<std::size_t>(next_level)] +=
                            q(next_level, j) *
                            static_cast<double>(
                                    coefficients_[static_cast<std::size_t>(j)]);
                }
                const double p = static_cast<double>(
                                         coefficients_[static_cast<std::size_t>(
                                                 level_)]) +
                                 tails_[static_cast<std::size_t>(level_)];
                lengths_[static_cast<std::size_t>(next_level)] =
                        lengths_[static_cast<std::size_t>(level_)] +
                        p * p * v(level_);
                if (next_level < skip_first_ && lengths_[0] == 0.0) {
                    skip_coordinate_test = true;
                }
                const double centered = std::floor(
                        -tails_[static_cast<std::size_t>(next_level)] + 0.5);
                if (!slong_from_double(
                            coefficients_[static_cast<std::size_t>(next_level)],
                            centered)) {
                    exhausted_ = true;
                    return false;
                }
                level_ = next_level;
            }

            while (true) {
                if (!skip_coordinate_test) {
                    if (++element_steps_ > max_element_steps_) {
                        element_step_limit_reached_ = true;
                        exhausted_ = true;
                        return false;
                    }
                    double p = static_cast<double>(
                                       coefficients_[static_cast<std::size_t>(
                                               level_)]) +
                               tails_[static_cast<std::size_t>(level_)];
                    if (lengths_[static_cast<std::size_t>(level_)] +
                                p * p * v(level_) <=
                        bound_) {
                        break;
                    }

                    if (!step(level_)) {
                        exhausted_ = true;
                        return false;
                    }
                    p = static_cast<double>(
                                coefficients_[static_cast<std::size_t>(
                                        level_)]) +
                        tails_[static_cast<std::size_t>(level_)];
                    if (lengths_[static_cast<std::size_t>(level_)] +
                                p * p * v(level_) <=
                        bound_) {
                        break;
                    }
                }

                skip_coordinate_test = false;
                increments_[static_cast<std::size_t>(level_)] = 1;
                ++level_;
                if (level_ >= dimension_) {
                    exhausted_ = true;
                    return false;
                }
                if (!step(level_)) {
                    exhausted_ = true;
                    return false;
                }
            }
        } while (level_ > 0);

        need_outer_step_ = true;
        return true;
    }
}

bool FiniteQuadraticFormEnumerationContext::current_row(
        flint::FmpzMatRef out) const noexcept {
    if (dimension_ <= 0 || flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != dimension_) {
        return false;
    }
    for (slong i = 0; i < dimension_; ++i) {
        flint::fmpz_set_si(flint::fmpz_mat_entry(out, 0, i),
                           coefficients_[static_cast<std::size_t>(i)]);
    }
    return true;
}

slong FiniteQuadraticFormEnumerationContext::element_steps() const noexcept {
    return element_steps_;
}

bool FiniteQuadraticFormEnumerationContext::element_step_limit_reached()
        const noexcept {
    return element_step_limit_reached_;
}

void FiniteQuadraticFormEnumerationContext::clear() noexcept {
    dimension_ = 0;
    quadratic_form_data_.clear();
    coefficients_.clear();
    increments_.clear();
    lengths_.clear();
    tails_.clear();
    started_ = false;
    exhausted_ = true;
    need_outer_step_ = false;
    element_steps_ = 0;
    element_step_limit_reached_ = false;
}

double FiniteQuadraticFormEnumerationContext::q(slong row,
                                            slong column) const noexcept {
    return quadratic_form_entry(quadratic_form_data_, dimension_, row, column);
}

double FiniteQuadraticFormEnumerationContext::v(slong index) const noexcept {
    return q(index, index);
}

bool FiniteQuadraticFormEnumerationContext::step(slong index) noexcept {
    slong& coefficient = coefficients_[static_cast<std::size_t>(index)];
    if (lengths_[static_cast<std::size_t>(index)] == 0.0) {
        if (coefficient == std::numeric_limits<slong>::max()) {
            return false;
        }
        ++coefficient;
        return true;
    }

    slong& increment = increments_[static_cast<std::size_t>(index)];
    if ((increment > 0 &&
         coefficient > std::numeric_limits<slong>::max() - increment) ||
        (increment < 0 &&
         coefficient < std::numeric_limits<slong>::min() - increment)) {
        return false;
    }
    coefficient += increment;
    increment = increment > 0 ? -1 - increment : 1 - increment;
    return true;
}

double arb_midpoint_double(flint::ArbConstRef value) noexcept {
    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value.raw()));
    if (!flint::arf_is_finite(midpoint)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double result = flint::arf_get_d(midpoint, ARF_RND_NEAR);
    return std::isfinite(result)
            ? result
            : std::numeric_limits<double>::quiet_NaN();
}

static bool quadratic_form_data_from_arb_gram(std::vector<double>& out,
                                const flint::ArbMat& gram) noexcept {
    const slong dimension = flint::arb_mat_nrows_value(gram);
    if (dimension <= 0 ||
        dimension != flint::arb_mat_ncols_value(gram)) {
        return false;
    }

    std::vector<double> double_gram(
            static_cast<std::size_t>(dimension * dimension), 0.0);
    for (slong i = 0; i < dimension; ++i) {
        for (slong j = 0; j < dimension; ++j) {
            const double value = arb_midpoint_double(
                    flint::arb_mat_entry_ref(gram, i, j));
            if (!std::isfinite(value)) {
                return false;
            }
            double_gram[static_cast<std::size_t>(i * dimension + j)] = value;
        }
    }
    return quadratic_form_data_from_double_gram(out, double_gram, dimension);
}

// The reference ideal LLL accepts a result only when it meets the LLL
// guarantees of the ideal lattice, with d = det(Gram) = |disc| N(I)^2:
//   |b_1|^2 <= 2^floor((n + 1) / 2) (d^(1/n) + 1)  and
//   prod_i |b_i|^2 <= 2^(n (n - 1) / 2) d,
// and otherwise redoes the reduction at a higher precision (see "Class groups
// and order units" in docs/reference/algorithms_and_sources.rst).
// Here |b_i|^2 is the Gram diagonal and log2_determinant is log2 d, taken
// from the exact |disc| N(I)^2 as the source does
// (ideal_t2_lattice_log2_determinant).  When an input cannot be read as a
// finite (positive) double the test is not decided and the basis is accepted
// as before.
static bool finite_enumeration_diagonal_within_lll_bounds(
        const std::vector<double>& diagonal,
        double log2_determinant) noexcept {
    const slong dimension = static_cast<slong>(diagonal.size());
    if (dimension <= 0 || !std::isfinite(log2_determinant)) {
        return true;
    }

    double log2_diagonal_product = 0.0;
    double log2_first_diagonal = 0.0;
    for (slong i = 0; i < dimension; ++i) {
        const double entry = diagonal[static_cast<std::size_t>(i)];
        if (!std::isfinite(entry) || entry <= 0.0) {
            return true;
        }
        log2_diagonal_product += std::log2(entry);
        if (i == 0) {
            log2_first_diagonal = std::log2(entry);
        }
    }

    const double n = static_cast<double>(dimension);
    const double first_bound =
            static_cast<double>((dimension + 1) / 2) +
            std::log2(std::exp2(log2_determinant / n) + 1.0);
    const double product_bound =
            (n * (n - 1.0)) / 2.0 + log2_determinant;
    if (!std::isfinite(first_bound) || !std::isfinite(product_bound)) {
        return true;
    }
    return log2_first_diagonal <= first_bound &&
           log2_diagonal_product <= product_bound;
}

static bool finite_enumeration_basis_within_lll_bounds(
        const flint::ArbMat& gram,
        double log2_determinant) noexcept {
    const slong dimension = flint::arb_mat_nrows_value(gram);
    if (dimension <= 0 || dimension != flint::arb_mat_ncols_value(gram)) {
        return true;
    }
    std::vector<double> diagonal(static_cast<std::size_t>(dimension));
    for (slong i = 0; i < dimension; ++i) {
        diagonal[static_cast<std::size_t>(i)] = arb_midpoint_double(
                flint::arb_mat_entry_ref(gram, i, i));
    }
    return finite_enumeration_diagonal_within_lll_bounds(diagonal,
                                                         log2_determinant);
}

// log2 of d = det(Gram) for the T2 Gram matrix of any basis of the ideal
// lattice.  The weighted Minkowski rows are the T2 form, so for the HNF basis
// det(Gram) = |disc(O)| det(HNF)^2 = |disc(O)| N(I)^2, and a unimodular LLL
// transform keeps it.  The source takes this exact value in its ideal LLL
// acceptance test (see "Class groups and order units" in
// docs/reference/algorithms_and_sources.rst).  NaN when the discriminant or
// the norm cannot be read or is zero, so the bound test is not decided.
static double ideal_t2_lattice_log2_determinant(const Ideal& ideal,
                                                const Order& order) noexcept {
    flint::Fmpz discriminant;
    flint::Fmpz norm;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        !ideal.norm(flint::FmpzRef(norm)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(norm))) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    // |x| = |m| 2^e with |m| in [1/2, 1), so log2 |x| = log2 |m| + e without
    // overflowing a double.
    const auto log2_abs = [](flint::FmpzConstRef value) noexcept {
        slong exponent = 0;
        const double mantissa = flint::fmpz_get_d_2exp(exponent, value);
        return std::log2(std::fabs(mantissa)) + static_cast<double>(exponent);
    };
    return log2_abs(flint::FmpzConstRef(discriminant)) +
           2.0 * log2_abs(flint::FmpzConstRef(norm));
}

static bool round_scaled_midpoint_to_integer(flint::FmpzRef out,
                                           flint::ArbConstRef value,
                                           slong scale) noexcept {
    if (scale < 0) {
        return false;
    }

    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value.raw()));
    if (!flint::arf_is_finite(midpoint)) {
        return false;
    }
    ::arf_mul_2exp_si(midpoint.raw(), midpoint.raw(), scale);
    flint::arf_get_fmpz(out, midpoint, ARF_RND_NEAR);
    return true;
}

static bool round_embedding_rows_full_rank(flint::FmpzMat& out,
                                       slong& used_scale,
                                       const flint::ArbMat& embedding_rows,
                                       slong min_scale,
                                       slong max_scale) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (rows <= 0 || rows != cols ||
        min_scale < kT2EnumerationInitialRoundingScale ||
        max_scale < min_scale ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != cols) {
        return false;
    }

    for (slong scale = min_scale; scale <= max_scale;) {
        for (slong i = 0; i < rows; ++i) {
            for (slong j = 0; j < cols; ++j) {
                if (!round_scaled_midpoint_to_integer(
                            flint::fmpz_mat_entry(out, i, j),
                            flint::arb_mat_entry_ref(embedding_rows, i, j),
                            scale)) {
                    return false;
                }
            }
        }
        if (flint::fmpz_mat_rank(flint::FmpzMatConstRef(out)) == rows) {
            used_scale = scale;
            return true;
        }
        if (scale > max_scale / 2) {
            break;
        }
        scale *= 2;
    }
    return false;
}

static bool finite_enumeration_basis_lll_reduce(flint::FmpzMat& reduced_basis,
                              flint::FmpzMat& transform,
                              slong& used_rounding_scale,
                              const flint::FmpzMat& basis,
                              const flint::ArbMat& embedding_rows,
                              slong precision,
                              slong min_rounding_scale,
                              const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = flint::fmpz_mat_nrows(basis);
    if (precision <= 0 || degree != flint::fmpz_mat_ncols(basis) ||
        flint::fmpz_mat_nrows(reduced_basis) != degree ||
        flint::fmpz_mat_ncols(reduced_basis) != degree ||
        flint::fmpz_mat_nrows(transform) != degree ||
        flint::fmpz_mat_ncols(transform) != degree) {
        return false;
    }

    flint::FmpzMat scaled_rows(degree, degree);
    flint::FmpzMat reduced_scaled_rows(degree, degree);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.finite_enumeration_basis_lll.round_embedding_rows");
        if (!round_embedding_rows_full_rank(scaled_rows, used_rounding_scale,
                                            embedding_rows,
                                            min_rounding_scale, precision)) {
            return false;
        }
    }

    flint::fmpz_mat_set(flint::FmpzMatRef(reduced_scaled_rows),
                        flint::FmpzMatConstRef(scaled_rows));
    flint::fmpz_mat_one(flint::FmpzMatRef(transform));
    flint::FmpzLll config(0.99, 0.51, Z_BASIS, APPROX);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_basis_lll.flint_lll");
        // FLINT's first LLL phase with an exact reducedness certificate
        // (lat::detail::certified_fmpz_lll, the Lat::lll_reduce route,
        // docs/reference/algorithms_and_sources.rst), which skips FLINT's
        // multiprecision check of the phase output.  The rows and transform
        // are what fmpz_lll returns on every route: a failed certificate
        // restores both and runs fmpz_lll, and rows with entries above
        // FLINT's 250-bit ULLL threshold go to fmpz_lll directly.
        const lat::detail::CertifiedLllRoute route =
                lat::detail::certified_fmpz_lll(reduced_scaled_rows.raw(),
                                                transform.raw(),
                                                flint::FmpzLllConstRef(config));
        if (route == lat::detail::CertifiedLllRoute::fallback) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_basis_lll.certified_lll_fallback");
        } else if (route == lat::detail::CertifiedLllRoute::ineligible) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_basis_lll.certified_lll_ineligible");
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.finite_enumeration_basis_lll.apply_transform");
        flint::fmpz_mat_mul(flint::FmpzMatRef(reduced_basis),
                            flint::FmpzMatConstRef(transform),
                            flint::FmpzMatConstRef(basis));
    }
    return true;
}

// The LLL reductions of one T2 setup at one precision, by the minimum
// rounding scale they started from.  finite_enumeration_basis_lll_reduce is
// a deterministic function of that scale once the ideal basis, its embedding
// rows and the precision are fixed, so when the double_rows route falls back
// to the Arb route (T2EnumerationGram) the Arb scale loop reuses the
// reductions the double_rows loop already made instead of repeating them;
// its result is the one it would compute.
struct T2LllReduction {
    slong min_rounding_scale = 0;
    bool reduced = false;
    slong rounding_scale = 0;
    flint::FmpzMat basis{0, 0};
    flint::FmpzMat transform{0, 0};
};
using T2LllReductions = std::vector<T2LllReduction>;

static bool finite_enumeration_basis_lll_reduce_reusing(
        flint::FmpzMat& reduced_basis,
        flint::FmpzMat& transform,
        slong& used_rounding_scale,
        const flint::FmpzMat& basis,
        const flint::ArbMat& embedding_rows,
        slong precision,
        slong min_rounding_scale,
        T2LllReductions* reductions,
        const DiagnosticsContext* diagnostics) noexcept {
    if (reductions != nullptr) {
        for (const T2LllReduction& reduction : *reductions) {
            if (reduction.min_rounding_scale != min_rounding_scale) {
                continue;
            }
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.lll_reduction_reused");
            if (!reduction.reduced) {
                return false;
            }
            flint::fmpz_mat_set(flint::FmpzMatRef(reduced_basis),
                                flint::FmpzMatConstRef(reduction.basis));
            flint::fmpz_mat_set(flint::FmpzMatRef(transform),
                                flint::FmpzMatConstRef(reduction.transform));
            used_rounding_scale = reduction.rounding_scale;
            return true;
        }
    }
    const bool reduced = finite_enumeration_basis_lll_reduce(
            reduced_basis, transform, used_rounding_scale, basis,
            embedding_rows, precision, min_rounding_scale, diagnostics);
    if (reductions != nullptr) {
        T2LllReduction reduction;
        reduction.min_rounding_scale = min_rounding_scale;
        reduction.reduced = reduced;
        if (reduced) {
            reduction.rounding_scale = used_rounding_scale;
            reduction.basis = flint::FmpzMat(flint::fmpz_mat_nrows(reduced_basis),
                                             flint::fmpz_mat_ncols(reduced_basis));
            reduction.transform = flint::FmpzMat(flint::fmpz_mat_nrows(transform),
                                                 flint::fmpz_mat_ncols(transform));
            flint::fmpz_mat_set(flint::FmpzMatRef(reduction.basis),
                                flint::FmpzMatConstRef(reduced_basis));
            flint::fmpz_mat_set(flint::FmpzMatRef(reduction.transform),
                                flint::FmpzMatConstRef(transform));
        }
        reductions->push_back(std::move(reduction));
    }
    return reduced;
}

static double euclidean_ball_volume(slong dimension) noexcept {
    if (dimension <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double pi = std::acos(-1.0);
    double volume = (dimension % 2 == 0) ? 1.0 : 2.0;
    for (slong n = dimension; n >= 2; n -= 2) {
        volume *= (2.0 * pi) / static_cast<double>(n);
    }
    return volume;
}

static bool finite_quadratic_form_volume_bound(
        double& out,
        const std::vector<double>& gram_schmidt_diagonal,
        slong target_factor_attempts) noexcept {
    const slong dimension =
            static_cast<slong>(gram_schmidt_diagonal.size());
    if (dimension <= 0) {
        return false;
    }

    const double ball_volume = euclidean_ball_volume(dimension);
    if (!std::isfinite(ball_volume) || ball_volume <= 0.0) {
        return false;
    }

    if (target_factor_attempts <= 0) {
        return false;
    }
    const double target =
            (4.0 * static_cast<double>(target_factor_attempts)) /
                          ball_volume;
    if (!std::isfinite(target) || target <= 0.0) {
        return false;
    }

    double product = 1.0;
    for (slong i = 0; i < dimension; ++i) {
        const double diagonal =
                gram_schmidt_diagonal[static_cast<std::size_t>(i)];
        if (!std::isfinite(diagonal) || diagonal <= 0.0) {
            return false;
        }
        product *= diagonal;
        const double candidate =
                std::pow(target * target * product,
                         1.0 / static_cast<double>(i + 1));
        if (!std::isfinite(candidate) || candidate <= 0.0) {
            return false;
        }
        if (i + 1 == dimension ||
            candidate < gram_schmidt_diagonal[
                                static_cast<std::size_t>(i + 1)]) {
            out = candidate;
            return true;
        }
    }

    return false;
}

// The Gram matrix of the embedding rows, G[i][j] = <row_i, row_j>.  G is
// symmetric, so only the lower triangle j <= i is computed, and each entry is
// copied to G[j][i].  Each lower entry is one arb_dot of two contiguous rows.
// Like every Arb operation, arb_dot returns a ball that contains the exact
// dot product of the input balls, so G stays a rigorous enclosure.  Its
// midpoint comes from FLINT's dot-product accumulation rather than a sum of
// rounded products, so an entry can differ from the arb_mul/arb_add sum by
// rounding noise.  The readers of G use the diagonal and the entry midpoints.
static bool build_embedding_gram(flint::ArbMat& out,
                          const flint::ArbMat& embedding_rows,
                          slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (precision <= 0 || flint::arb_mat_nrows_value(out) != rows ||
        flint::arb_mat_ncols_value(out) != rows) {
        return false;
    }

    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j <= i; ++j) {
            flint::ArbRef entry = flint::arb_mat_entry_ref(out, i, j);
            ::arb_dot(entry.raw(), nullptr, 0,
                      flint::arb_mat_entry_ref(embedding_rows, i, 0).raw(), 1,
                      flint::arb_mat_entry_ref(embedding_rows, j, 0).raw(), 1,
                      cols, precision);
            if (j != i) {
                flint::arb_set(flint::arb_mat_entry_ref(out, j, i),
                               flint::ArbConstRef(entry.raw()));
            }
        }
    }
    return true;
}

bool arb_mat_midpoints_double(std::vector<double>& out,
                              const flint::ArbMat& rows) noexcept {
    const slong row_count = flint::arb_mat_nrows_value(rows);
    const slong col_count = flint::arb_mat_ncols_value(rows);
    out.assign(static_cast<std::size_t>(row_count * col_count), 0.0);
    for (slong i = 0; i < row_count; ++i) {
        for (slong j = 0; j < col_count; ++j) {
            const double value =
                    arb_midpoint_double(flint::arb_mat_entry_ref(rows, i, j));
            if (!std::isfinite(value)) {
                return false;
            }
            out[static_cast<std::size_t>(i * col_count + j)] = value;
        }
    }
    return true;
}

// The double rows are the exact reduced basis (order coordinates) times the
// order's Minkowski rows, so the cancellation in a row is bounded by the
// conditioning of the order basis rather than by the size of the LLL
// transform.  Each basis entry is read exactly (|b| <= 2^53).  G is
// symmetric: the lower triangle is computed and copied, as in
// build_embedding_gram.
bool reduced_basis_t2_gram_double(std::vector<double>& gram,
                                  flint::FmpzMatConstRef reduced_basis,
                                  const std::vector<double>& order_rows) noexcept {
    const slong degree = flint::fmpz_mat_nrows(reduced_basis);
    if (degree <= 0 || flint::fmpz_mat_ncols(reduced_basis) != degree ||
        order_rows.size() != static_cast<std::size_t>(degree * degree)) {
        return false;
    }
    const auto index = [degree](slong i, slong j) noexcept {
        return static_cast<std::size_t>(i * degree + j);
    };

    std::vector<double> basis(static_cast<std::size_t>(degree * degree));
    for (slong i = 0; i < degree; ++i) {
        for (slong j = 0; j < degree; ++j) {
            const flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(reduced_basis, i, j);
            if (::fmpz_bits(entry.raw()) > 53) {
                return false;
            }
            basis[index(i, j)] = flint::fmpz_get_d(entry);
        }
    }

    std::vector<double> rows(static_cast<std::size_t>(degree * degree), 0.0);
    for (slong i = 0; i < degree; ++i) {
        for (slong k = 0; k < degree; ++k) {
            const double coefficient = basis[index(i, k)];
            if (coefficient == 0.0) {
                continue;
            }
            for (slong j = 0; j < degree; ++j) {
                rows[index(i, j)] += coefficient * order_rows[index(k, j)];
            }
        }
    }
    if (T2DoubleRowsTestSeams::mode() ==
        T2DoubleRowsTestSeams::Mode::nonfinite_rows) {
        rows[0] = std::numeric_limits<double>::quiet_NaN();
    }
    for (const double value : rows) {
        if (!std::isfinite(value)) {
            return false;
        }
    }

    gram.assign(static_cast<std::size_t>(degree * degree), 0.0);
    for (slong i = 0; i < degree; ++i) {
        for (slong j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (slong k = 0; k < degree; ++k) {
                sum += rows[index(i, k)] * rows[index(j, k)];
            }
            if (!std::isfinite(sum)) {
                return false;
            }
            gram[index(i, j)] = sum;
            gram[index(j, i)] = sum;
        }
    }
    return true;
}

// Reduces the ideal basis from the first full-rank rounding at or above
// min_rounding_scale and builds the Gram matrix of the reduced rows.  The
// LLL bound test takes d from the exact |disc| N(I)^2 and positive
// definiteness is read from the double Gram-Schmidt of the Gram midpoints
// (quadratic_form_data_from_arb_gram), so no Arb Cholesky factor is formed.
// The reduced basis is the LLL transform times the ideal
// HNF, so it spans the same lattice at every scale.
static bool reduce_finite_enumeration_basis_from_scale(
        flint::FmpzMat& basis,
        flint::FmpzMat& transform,
        slong& rounding_scale,
        flint::ArbMat& embedding_rows,
        flint::ArbMat& gram,
        const flint::FmpzMat& original_basis,
        const flint::ArbMat& original_embedding_rows,
        slong precision,
        slong embedding_precision,
        slong min_rounding_scale,
        T2LllReductions* reductions,
        const DiagnosticsContext* diagnostics) noexcept {
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.lll_reduce");
        if (!finite_enumeration_basis_lll_reduce_reusing(
                    basis, transform, rounding_scale, original_basis,
                    original_embedding_rows, precision, min_rounding_scale,
                    reductions, diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail, "T2 context LLL reduction failed");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.finite_enumeration_setup.reduced_embedding_rows");
        if (!detail::multiply_integer_arb_matrices(
                    embedding_rows, flint::FmpzMatConstRef(transform),
                    original_embedding_rows, embedding_precision)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context reduced embedding rows failed");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.gram");
        if (!build_embedding_gram(gram, embedding_rows, precision)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail, "T2 context Gram build failed");
            return false;
        }
    }
    return true;
}

enum class T2BasisSelection { failed, fallback, accepted };

// The rounding-scale loop of one T2 enumeration setup, started at
// start_rounding_scale (see build_finite_ideal_t2_enumeration_data_at_precision
// for the rule).  On `accepted`, basis and quadratic_form_data hold the
// reduced basis accepted at accepted_rounding_scale; on `fallback`, they hold
// the first reduced basis that had quadratic-form data.
static T2BasisSelection select_finite_enumeration_basis(
        flint::FmpzMat& basis,
        std::vector<double>& quadratic_form_data,
        slong& accepted_rounding_scale,
        flint::FmpzMat& transform,
        flint::ArbMat& embedding_rows,
        flint::ArbMat& gram,
        const flint::FmpzMat& original_basis,
        const flint::ArbMat& original_embedding_rows,
        double log2_determinant,
        slong precision,
        slong embedding_precision,
        slong start_rounding_scale,
        T2LllReductions* reductions,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = flint::fmpz_mat_nrows(original_basis);
    flint::FmpzMat fallback_basis(degree, degree);
    std::vector<double> fallback_quadratic_form_data;
    bool have_fallback = false;
    slong min_rounding_scale = start_rounding_scale;
    for (;;) {
        slong rounding_scale = 0;
        if (!reduce_finite_enumeration_basis_from_scale(
                    basis, transform, rounding_scale, embedding_rows, gram,
                    original_basis, original_embedding_rows,
                    precision, embedding_precision, min_rounding_scale,
                    reductions, diagnostics)) {
            if (!have_fallback) {
                return T2BasisSelection::failed;
            }
            basis.swap(fallback_basis);
            quadratic_form_data.swap(fallback_quadratic_form_data);
            return T2BasisSelection::fallback;
        }
        bool have_quadratic_form = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.quadratic_form_data");
            const bool within_lll_bounds =
                    finite_enumeration_basis_within_lll_bounds(gram,
                                                               log2_determinant);
            have_quadratic_form =
                    quadratic_form_data_from_arb_gram(quadratic_form_data,
                                                      gram);
            if (have_quadratic_form && within_lll_bounds) {
                accepted_rounding_scale = rounding_scale;
                return T2BasisSelection::accepted;
            }
            if (!have_quadratic_form) {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "T2 context quadratic_form_data data failed");
            } else {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "T2 context reduced basis outside LLL bounds");
            }
        }
        if (have_quadratic_form && !have_fallback) {
            fallback_basis.swap(basis);
            fallback_quadratic_form_data.swap(quadratic_form_data);
            have_fallback = true;
        }
        if (rounding_scale > precision / 2) {
            if (!have_fallback) {
                return T2BasisSelection::failed;
            }
            basis.swap(fallback_basis);
            quadratic_form_data.swap(fallback_quadratic_form_data);
            return T2BasisSelection::fallback;
        }
        min_rounding_scale = 2 * rounding_scale;
    }
}

// Whether the Gram-Schmidt pivots of double quadratic-form data reproduce
// d = det(Gram) = |disc(O)| N(I)^2, the exact value the LLL-bound test takes
// (ideal_t2_lattice_log2_determinant): the pivots are the squared
// Gram-Schmidt norms, whose product is det(Gram) for any basis of the
// lattice.  The double rows of reduced_basis_t2_gram_double are sums of
// products |b_ik| |o_kj|, and when those terms are many orders larger than
// the reduced row they cancel to, rounding leaves errors comparable to the
// row itself: the Gram is then still positive definite and can pass the
// pivot and LLL-bound tests, but it is not the T2 form of the lattice (on
// multiquadratic_m23_m31_m47 the pivot product was 2^17 times d and varied
// with the rounding scale).  The product of the pivots is a consequence of
// every Gram entry, so this compares the computed data with an invariant
// known exactly and falls back to the Arb setup on a mismatch.
//
// Tolerance: 2^-20 in log2 d (a relative determinant error of about 7e-7).
// With rounding errors only, the discrepancy is a small multiple of
// n 2^-53 times the conditioning of the reduced basis; on the proven corpus
// rows measured it stayed below 3e-8 (cyclotomic_real_41_deg20, degree 20)
// and below 1e-12 on most rows, so data with no cancellation does not fall
// back.  On multiquadratic_m23_m31_m47 no setup was below 1e-6 and most were
// above one bit, so they all take the Arb route.  An accepted discrepancy
// changes the volume-based initial enumeration bound
// (finite_quadratic_form_volume_bound) by a factor of at most 2^(2^-20 / n)
// when that bound uses all n pivots.
// When d is not available (NaN) the data cannot be checked and is not used.
constexpr double kT2DoubleRowsLog2DeterminantTolerance = 0x1p-20;

static bool double_gram_schmidt_reproduces_determinant(
        const std::vector<double>& quadratic_form_data,
        slong degree,
        double log2_determinant) noexcept {
    if (degree <= 0 || !std::isfinite(log2_determinant) ||
        quadratic_form_data.size() !=
                static_cast<std::size_t>(degree * degree)) {
        return false;
    }
    double log2_pivot_product = 0.0;
    for (slong i = 0; i < degree; ++i) {
        const double pivot =
                quadratic_form_data[static_cast<std::size_t>(i * degree + i)];
        if (!std::isfinite(pivot) || pivot <= 0.0) {
            return false;
        }
        log2_pivot_product += std::log2(pivot);
    }
    return std::fabs(log2_pivot_product - log2_determinant) <=
           kT2DoubleRowsLog2DeterminantTolerance;
}

// The rounding-scale loop of select_finite_enumeration_basis with the
// reduced rows and Gram formed in double (T2EnumerationGram::double_rows).
// The rounding and the LLL are those of the Arb loop, so a scale gives the
// same transform and reduced basis.  Returns true with the basis accepted at
// accepted_rounding_scale.  Returns false, for the caller to run the Arb
// loop, when the LLL step fails, the double rows or Gram are not available
// (reduced_basis_t2_gram_double), the double Gram-Schmidt pivot test fails,
// the pivots do not reproduce d (double_gram_schmidt_reproduces_determinant;
// determinant_mismatch is then set and the caller goes straight to the Arb
// route), or no scale from start_rounding_scale on passes the LLL-bound test.
static bool select_finite_enumeration_basis_double_rows(
        flint::FmpzMat& basis,
        std::vector<double>& quadratic_form_data,
        slong& accepted_rounding_scale,
        flint::FmpzMat& transform,
        const flint::FmpzMat& original_basis,
        const flint::ArbMat& original_embedding_rows,
        const std::vector<double>& order_rows,
        double log2_determinant,
        slong precision,
        slong start_rounding_scale,
        T2LllReductions& reductions,
        bool& determinant_mismatch,
        const DiagnosticsContext* diagnostics) noexcept {
    determinant_mismatch = false;
    const slong degree = flint::fmpz_mat_nrows(original_basis);
    std::vector<double> gram;
    std::vector<double> diagonal(static_cast<std::size_t>(degree));
    slong min_rounding_scale = start_rounding_scale;
    for (;;) {
        slong rounding_scale = 0;
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.finite_enumeration_setup.lll_reduce");
            if (!finite_enumeration_basis_lll_reduce_reusing(
                        basis, transform, rounding_scale, original_basis,
                        original_embedding_rows, precision, min_rounding_scale,
                        &reductions, diagnostics)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.double_rows_gram");
            if (!reduced_basis_t2_gram_double(gram,
                                              flint::FmpzMatConstRef(basis),
                                              order_rows)) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.finite_enumeration_setup.double_rows_unavailable");
                return false;
            }
            if (!quadratic_form_data_from_double_gram(quadratic_form_data,
                                                      gram, degree)) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.finite_enumeration_setup.double_rows_pivot_failed");
                return false;
            }
        }
        if (!double_gram_schmidt_reproduces_determinant(
                    quadratic_form_data, degree, log2_determinant)) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.double_rows_determinant_mismatch");
            determinant_mismatch = true;
            return false;
        }
        for (slong i = 0; i < degree; ++i) {
            diagonal[static_cast<std::size_t>(i)] =
                    gram[static_cast<std::size_t>(i * degree + i)];
        }
        if (finite_enumeration_diagonal_within_lll_bounds(diagonal,
                                                          log2_determinant)) {
            accepted_rounding_scale = rounding_scale;
            return true;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.finite_enumeration_setup.double_rows_outside_lll_bounds");
        if (rounding_scale > precision / 2) {
            return false;
        }
        min_rounding_scale = 2 * rounding_scale;
    }
}

static bool build_finite_ideal_t2_enumeration_data_at_precision(
        FiniteIdealT2EnumerationData& out,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache,
        slong target_factor_attempts,
        bool double_rows,
        bool& used_double_rows)
        noexcept {
    used_double_rows = false;
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || precision <= 0) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail, "T2 context initial guard failed");
        return false;
    }

    // The source finite quadratic-form search enumerates the ideal lattice after
    // pr_hnf/ideal multiplication.  Silex stores the same lattice in the
    // active maximal-order basis; exact nf.zk coordinates are a private
    // representation, not a route contract.
    const slong degree = order->degree();
    const slong embedding_precision =
            precision + kIdealLatticeRoundingGuardBits;
    flint::FmpzMat original_basis(degree, degree);
    flint::FmpzMat basis(degree, degree);
    flint::FmpzMat transform(degree, degree);
    flint::ArbMat original_embedding_rows(degree, degree);
    flint::ArbMat embedding_rows(degree, degree);
    flint::ArbMat gram(degree, degree);
    std::vector<double> quadratic_form_data;

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.ideal_hnf");
        if (!ideal.get_hnf(flint::FmpzMatRef(original_basis))) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context ideal HNF read failed");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.finite_enumeration_setup.original_embedding_rows");
        if (!detail::build_ideal_minkowski_embedding_rows(
                    original_embedding_rows,
                    flint::FmpzMatConstRef(original_basis), *order,
                    embedding_precision, embedding_cache, diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context original embedding rows failed");
            return false;
        }
    }
    // Rounding embedding rows at the first full-rank scale follows the source,
    // which rounds the embeddings of a T2-reduced integral basis.  Silex
    // rounds the ideal's rows in the active order basis, which need not be
    // reduced.  When that basis is far from reduced, the rounding error
    // multiplied by the LLL transform can dominate the reduced vectors: the
    // result misses the LLL bounds, or the double Gram-Schmidt below fails.
    // Both read only the transform and the Gram midpoints, so a higher working
    // precision with the same transform gives the same result.  In either case
    // the same lattice is reduced again from the next finer full-rank
    // rounding, as the reference ideal LLL redoes the reduction at a higher
    // precision.  When a finer attempt fails, or no finer scale is left, the
    // first basis that had quadratic-form data is kept.  So a basis that a
    // single reduction from the first scale accepted before is still the
    // result when no finer basis is accepted.  The ideal lattice and the
    // bound rule applied to the result do not change.
    //
    // The scale loop starts at the smallest scale accepted so far for an ideal
    // of this order (the embedding cache's hint) rather than at the first
    // scale, skipping the coarser attempts that the ideals before it needed
    // to reject.  From the hint on, the attempts are those of the full loop.
    // An ideal whose full loop would accept a coarser scale instead gets the
    // basis reduced from the hint's rounding: the same lattice, accepted by
    // the same LLL bound test.  When no scale at or above the hint is
    // accepted, the loop runs again from the first scale, so the fallback and
    // failure results are those of the full loop.
    slong start_rounding_scale = kT2EnumerationInitialRoundingScale;
    if (embedding_cache != nullptr &&
        embedding_cache->t2_rounding_scale_hint > start_rounding_scale) {
        start_rounding_scale = embedding_cache->t2_rounding_scale_hint;
    }
    const double log2_determinant =
            ideal_t2_lattice_log2_determinant(ideal, *order);
    slong accepted_rounding_scale = 0;
    double initial_enumeration_bound = 0.0;
    T2LllReductions reductions;
    if (double_rows) {
        // The double_rows route (T2EnumerationGram) runs the same two scale
        // loops as the Arb route below, from the hint and then from the first
        // scale.  When neither accepts a basis it falls back to the whole Arb
        // route, so it never fails where the Arb route would not.
        std::vector<double> order_rows;
        bool have_order_rows = false;
        if (embedding_cache != nullptr &&
            embedding_cache->has_order_embedding_rows &&
            embedding_cache->has_order &&
            embedding_cache->order.has_same_data(*order) &&
            embedding_cache->embedding_precision == embedding_precision) {
            have_order_rows = arb_mat_midpoints_double(
                    order_rows, embedding_cache->order_embedding_rows);
        } else {
            flint::FmpzMat identity(degree, degree);
            flint::fmpz_mat_one(flint::FmpzMatRef(identity));
            have_order_rows =
                    detail::build_ideal_minkowski_embedding_rows(
                            embedding_rows, flint::FmpzMatConstRef(identity),
                            *order, embedding_precision, nullptr,
                            diagnostics) &&
                    arb_mat_midpoints_double(order_rows, embedding_rows);
        }
        if (have_order_rows) {
            bool determinant_mismatch = false;
            used_double_rows = select_finite_enumeration_basis_double_rows(
                    basis, quadratic_form_data, accepted_rounding_scale,
                    transform, original_basis, original_embedding_rows,
                    order_rows, log2_determinant, precision,
                    start_rounding_scale, reductions, determinant_mismatch,
                    diagnostics);
            if (!used_double_rows && !determinant_mismatch &&
                start_rounding_scale != kT2EnumerationInitialRoundingScale) {
                used_double_rows = select_finite_enumeration_basis_double_rows(
                        basis, quadratic_form_data, accepted_rounding_scale,
                        transform, original_basis, original_embedding_rows,
                        order_rows, log2_determinant, precision,
                        kT2EnumerationInitialRoundingScale, reductions,
                        determinant_mismatch, diagnostics);
            }
        }
        if (used_double_rows) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.initial_bound");
            used_double_rows = initial_bound_from_quadratic_form(
                    initial_enumeration_bound, quadratic_form_data, degree,
                    target_factor_attempts);
        }
        if (!used_double_rows) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.finite_enumeration_setup.double_rows_arb_fallback");
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context double rows not accepted; Arb route");
        }
    }
    T2BasisSelection selection = T2BasisSelection::accepted;
    if (!used_double_rows) {
        selection = select_finite_enumeration_basis(
                basis, quadratic_form_data, accepted_rounding_scale,
                transform, embedding_rows, gram, original_basis,
                original_embedding_rows, log2_determinant, precision,
                embedding_precision, start_rounding_scale,
                double_rows ? &reductions : nullptr, diagnostics);
    }
    if (!used_double_rows && selection != T2BasisSelection::accepted &&
        start_rounding_scale != kT2EnumerationInitialRoundingScale) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "T2 context rounding scale hint not accepted; full loop");
        selection = select_finite_enumeration_basis(
                basis, quadratic_form_data, accepted_rounding_scale,
                transform, embedding_rows, gram, original_basis,
                original_embedding_rows, log2_determinant, precision,
                embedding_precision, kT2EnumerationInitialRoundingScale,
                double_rows ? &reductions : nullptr, diagnostics);
    }
    if (selection == T2BasisSelection::failed) {
        return false;
    }
    if (selection == T2BasisSelection::fallback) {
        accepted_rounding_scale = 0;
    } else if (embedding_cache != nullptr &&
               (embedding_cache->t2_rounding_scale_hint == 0 ||
                accepted_rounding_scale <
                        embedding_cache->t2_rounding_scale_hint)) {
        // A running minimum: a later ideal accepted only at a finer scale
        // does not raise the hint.
        embedding_cache->t2_rounding_scale_hint = accepted_rounding_scale;
    }

    if (!used_double_rows) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.initial_bound");
        if (!initial_bound_from_quadratic_form(initial_enumeration_bound,
                                             quadratic_form_data, degree,
                                             target_factor_attempts)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context initial bound computation failed");
            return false;
        }
    }

    out.basis = std::move(basis);
    out.quadratic_form_data = std::move(quadratic_form_data);
    out.initial_bound_value = initial_enumeration_bound;
    out.rounding_scale = accepted_rounding_scale;
    return true;
}

bool build_finite_ideal_t2_enumeration_data_with_retry(
        FiniteIdealT2EnumerationData& out,
        const Ideal& ideal,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache,
        slong target_factor_attempts,
        slong initial_precision,
        T2EnumerationGram gram)
        noexcept {
    if (initial_precision <= 0 ||
        initial_precision > kT2EnumerationMaxPrecision) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail, "T2 context initial precision invalid");
        return false;
    }
    for (slong precision = initial_precision;
         precision <= kT2EnumerationMaxPrecision;) {
        // The double_rows route runs at the start precision only; the
        // precision retries take the Arb route.
        bool used_double_rows = false;
        if (build_finite_ideal_t2_enumeration_data_at_precision(
                    out, ideal, precision, diagnostics, embedding_cache,
                    target_factor_attempts,
                    gram == T2EnumerationGram::double_rows &&
                            precision == initial_precision,
                    used_double_rows)) {
            out.double_rows = used_double_rows;
            return true;
        }
        if (precision == kT2EnumerationMaxPrecision) {
            break;
        }
        slong next_precision = (6 * precision + 4) / 5;
        if (next_precision <= precision) {
            next_precision = precision + 1;
        }
        if (next_precision > kT2EnumerationMaxPrecision) {
            next_precision = kT2EnumerationMaxPrecision;
        }
        precision = next_precision;
    }
    SILEX_LOG(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
              "T2 context precision retries exhausted");
    return false;
}

}  // namespace detail::relation_search
}  // namespace silex
