#include <silex/class_group.hpp>

#include <silex/flint/fmpz_lll.hpp>

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

constexpr slong kT2EnumerationInitialPrecision = 200;
constexpr slong kT2EnumerationMaxPrecision = 1000;
constexpr slong kIdealLatticeRoundingGuardBits = 64;

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
        const std::vector<double>& gram_schmidt_diagonal) noexcept;

static double quadratic_form_entry(const std::vector<double>& quadratic_form_data,
                     slong dimension,
                     slong row,
                     slong column) noexcept {
    return quadratic_form_data[static_cast<std::size_t>(row * dimension + column)];
}

static bool initial_bound_from_quadratic_form(double& out,
                                     const std::vector<double>& quadratic_form_data,
                                     slong dimension) noexcept {
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
                                        gram_schmidt_diagonal)) {
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
                                       const flint::ArbMat& embedding_rows,
                                       slong max_scale) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (rows <= 0 || rows != cols || max_scale < 4 ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != cols) {
        return false;
    }

    for (slong scale = 4; scale <= max_scale;) {
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
                              const flint::FmpzMat& basis,
                              const flint::ArbMat& embedding_rows,
                              slong precision,
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
        if (!round_embedding_rows_full_rank(scaled_rows, embedding_rows,
                                               precision)) {
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
        fmpz_lll(reduced_scaled_rows.raw(), transform.raw(), config.raw());
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
        const std::vector<double>& gram_schmidt_diagonal) noexcept {
    const slong dimension =
            static_cast<slong>(gram_schmidt_diagonal.size());
    if (dimension <= 0) {
        return false;
    }

    const double ball_volume = euclidean_ball_volume(dimension);
    if (!std::isfinite(ball_volume) || ball_volume <= 0.0) {
        return false;
    }

    const double target = (4.0 * static_cast<double>(kMaxFactorAttempts)) /
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

static bool build_embedding_gram(flint::ArbMat& out,
                          const flint::ArbMat& embedding_rows,
                          slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (precision <= 0 || flint::arb_mat_nrows_value(out) != rows ||
        flint::arb_mat_ncols_value(out) != rows) {
        return false;
    }

    flint::Arb sum;
    flint::Arb term;
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < rows; ++j) {
            flint::arb_zero(sum);
            for (slong k = 0; k < cols; ++k) {
                flint::arb_mul(term,
                               flint::arb_mat_entry_ref(embedding_rows, i, k),
                               flint::arb_mat_entry_ref(embedding_rows, j, k),
                               precision);
                flint::arb_add(sum, sum, term, precision);
            }
            flint::arb_set(flint::arb_mat_entry_ref(out, i, j),
                           flint::ArbConstRef(sum));
        }
    }
    return true;
}

static bool build_finite_ideal_t2_enumeration_data_at_precision(
        FiniteIdealT2EnumerationData& out,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache)
        noexcept {
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
    flint::ArbMat cholesky(degree, degree);
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
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.lll_reduce");
        if (!finite_enumeration_basis_lll_reduce(basis, transform, original_basis,
                                      original_embedding_rows, precision,
                                      diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context LLL reduction failed");
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
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.cholesky");
        if (flint::arb_mat_cho(cholesky, gram, precision) == 0) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail, "T2 context Cholesky failed");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.quadratic_form_data");
        if (!quadratic_form_data_from_arb_gram(quadratic_form_data, gram)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context quadratic_form_data data failed");
            return false;
        }
    }

    double initial_enumeration_bound = 0.0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.finite_enumeration_setup.initial_bound");
        if (!initial_bound_from_quadratic_form(initial_enumeration_bound,
                                             quadratic_form_data, degree)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "T2 context initial bound computation failed");
            return false;
        }
    }

    out.basis = std::move(basis);
    out.quadratic_form_data = std::move(quadratic_form_data);
    out.initial_bound_value = initial_enumeration_bound;
    return true;
}

bool build_finite_ideal_t2_enumeration_data_with_retry(
        FiniteIdealT2EnumerationData& out,
        const Ideal& ideal,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache)
        noexcept {
    for (slong precision = kT2EnumerationInitialPrecision;
         precision <= kT2EnumerationMaxPrecision;) {
        if (build_finite_ideal_t2_enumeration_data_at_precision(
                    out, ideal, precision, diagnostics, embedding_cache)) {
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
