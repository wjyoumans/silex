#include <silex/class_group.hpp>

#include <silex/flint/fmpz_lll.hpp>

#include "class_group_internal.hpp"
#include "ideal_lattice_lll_internal.hpp"
#include "ideal_minkowski_embedding_internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace silex {
namespace detail::relation_search {

constexpr slong kIdealLatticeLllMaxPrecision = WORD(1) << 18;
constexpr slong kIdealLatticeRoundingGuardBits = 64;

// reference implementation NfOrd/NfOrd.jl:minkowski_gram_mat_scaled rounds
// minkowski_matrix(O, prec + 2), whose entries meet the *absolute* tolerance
// 2^-(prec + 2) (minkowski_map raises its working precision until they do).
// Silex embeds at the relative precision `precision + 64`, which suffices
// unless the embedded entries are large (high-degree power bases, weights).
// Intentional deviation in the precision schedule only: when the rounding
// check rejects an entry, the next attempt keeps `precision` and raises the
// embedding precision by the measured shortfall.  It raises again at the
// same precision only while the measured shortfall at least halves, so a
// first shortfall s allows at most about log2(s) + 1 raises; otherwise
// precision doubles as before, keeping the raise.  Calls in which
// every attempt passes the rounding check run exactly as before.
enum class IdealLatticeReductionStatus {
    success,
    retry_precision,
    retry_embedding_precision,
    failed
};

class IdealLatticePrecisionSchedule {
  public:
    explicit IdealLatticePrecisionSchedule(slong precision) noexcept
        : active_precision_(precision) {}

    slong active_precision() const noexcept { return active_precision_; }

    slong embedding_precision() const noexcept {
        return active_precision_ + kIdealLatticeRoundingGuardBits +
               extra_embedding_bits_;
    }

    // Moves to the next attempt after a retry status.  Returns false when
    // the precision cap is reached.
    bool advance(IdealLatticeReductionStatus status,
                 slong embedding_shortfall_bits) noexcept {
        if (status == IdealLatticeReductionStatus::retry_embedding_precision &&
            embedding_shortfall_bits > 0 &&
            (last_shortfall_bits_ == 0 ||
             2 * embedding_shortfall_bits <= last_shortfall_bits_) &&
            embedding_shortfall_bits <=
                    kIdealLatticeLllMaxPrecision - extra_embedding_bits_) {
            extra_embedding_bits_ += embedding_shortfall_bits;
            last_shortfall_bits_ = embedding_shortfall_bits;
            return true;
        }
        if (active_precision_ > kIdealLatticeLllMaxPrecision / 2) {
            return false;
        }
        active_precision_ *= 2;
        last_shortfall_bits_ = 0;
        return true;
    }

  private:
    slong active_precision_;
    slong extra_embedding_bits_ = 0;
    // Shortfall that triggered the last raise at this precision; 0 if none.
    slong last_shortfall_bits_ = 0;
};

static bool ideal_lattice_round_scaled_arb_to_fmpz(
        flint::FmpzRef out,
        flint::ArbConstRef value,
        slong scale,
        slong precision) noexcept {
    flint::Arb scaled;
    flint::Arb radius;
    flint::Arf radius_upper;
    flint::Arf midpoint;
    flint::arb_mul_2exp_si(scaled, value.raw(), scale);
    flint::arb_get_rad_arb(radius, scaled);
    flint::arb_get_ubound_arf(radius_upper, radius, precision);
    if (!flint::arf_is_finite(radius_upper)) {
        return false;
    }

    const double radius_value =
            flint::arf_get_d(radius_upper, ARF_RND_UP);
    if (!std::isfinite(radius_value) || radius_value > 0.1) {
        return false;
    }

    flint::arf_set(midpoint, arb_midref(scaled.raw()));
    if (!flint::arf_is_finite(midpoint)) {
        return false;
    }
    flint::arf_get_fmpz(out, midpoint, ARF_RND_NEAR);
    return true;
}

// After ideal_lattice_scaled_embedding_rows rejects a radius, measures how
// many more embedding bits bring every scaled radius above 0.1 below
// 2^-4 < 0.1, assuming the absolute radius shrinks by one bit per bit of
// embedding precision.  Returns false, leaving the plain precision
// doubling in charge, when a ball is not finite or no shortfall is found.
static bool ideal_lattice_rounding_shortfall_bits(
        slong& out,
        const flint::ArbMat& embedding_rows,
        slong scale,
        slong precision) noexcept {
    out = 0;
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    flint::Arb scaled;
    flint::Arb radius;
    flint::Arf radius_upper;
    flint::Arf midpoint;
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            flint::arb_mul_2exp_si(
                    scaled,
                    flint::arb_mat_entry_ref(embedding_rows, i, j).raw(),
                    scale);
            flint::arb_get_rad_arb(radius, scaled);
            flint::arb_get_ubound_arf(radius_upper, radius, precision);
            flint::arf_set(midpoint, arb_midref(scaled.raw()));
            if (!flint::arf_is_finite(radius_upper) ||
                !flint::arf_is_finite(midpoint)) {
                out = 0;
                return false;
            }
            const double radius_value =
                    flint::arf_get_d(radius_upper, ARF_RND_UP);
            if (std::isfinite(radius_value) && radius_value <= 0.1) {
                continue;
            }
            // |radius| < 2^exponent, with exponent >= -3 here.
            const slong exponent =
                    ::arf_abs_bound_lt_2exp_si(radius_upper.raw());
            if (exponent >= kIdealLatticeLllMaxPrecision) {
                out = 0;
                return false;
            }
            out = std::max(out, exponent + 4);
        }
    }
    return out > 0;
}

static IdealLatticeReductionStatus ideal_lattice_rounding_retry_status(
        slong& embedding_shortfall_bits,
        const flint::ArbMat& embedding_rows,
        slong scale,
        slong precision) noexcept {
    return ideal_lattice_rounding_shortfall_bits(
                   embedding_shortfall_bits, embedding_rows, scale,
                   precision)
                   ? IdealLatticeReductionStatus::retry_embedding_precision
                   : IdealLatticeReductionStatus::retry_precision;
}

static bool ideal_lattice_scaled_embedding_rows(
        flint::FmpzMat& out,
        const flint::ArbMat& embedding_rows,
        slong scale,
        slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (precision <= 0 || flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != cols) {
        return false;
    }

    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            if (!ideal_lattice_round_scaled_arb_to_fmpz(
                        flint::fmpz_mat_entry(out, i, j),
                        flint::arb_mat_entry_ref(embedding_rows, i, j),
                        scale, precision)) {
                return false;
            }
        }
    }
    return true;
}
static bool ideal_lattice_scaled_minkowski_gram_numerator_without_diagonal(
        flint::FmpzMat& out,
        const flint::ArbMat& embedding_rows,
        slong scale,
        slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (scale < 0 || precision <= 0 ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != rows) {
        return false;
    }

    flint::FmpzMat scaled_rows(rows, cols);
    if (!ideal_lattice_scaled_embedding_rows(scaled_rows, embedding_rows, scale,
                                     precision)) {
        return false;
    }

    flint::fmpz_mat_gram(flint::FmpzMatRef(out),
                         flint::FmpzMatConstRef(scaled_rows));
    flint::fmpz_mat_scalar_tdiv_q_2exp(flint::FmpzMatRef(out),
                                       flint::FmpzMatConstRef(out),
                                       static_cast<ulong>(scale));
    return true;
}

static bool ideal_lattice_scaled_minkowski_gram_numerator(
        flint::FmpzMat& out,
        const flint::ArbMat& embedding_rows,
        slong scale,
        slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (!ideal_lattice_scaled_minkowski_gram_numerator_without_diagonal(
                out, embedding_rows, scale, precision)) {
        return false;
    }
    if (rows <= 0 || cols <= 0 ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != rows ||
        flint::fmpz_mat_nrows(out) != flint::fmpz_mat_ncols(out)) {
        return false;
    }
    for (slong i = 0; i < rows; ++i) {
        flint::FmpzRef diagonal = flint::fmpz_mat_entry(out, i, i);
        flint::fmpz_add_ui(diagonal,
                           flint::FmpzConstRef(diagonal.raw()),
                           static_cast<ulong>(rows));
    }
    return true;
}

static bool gram_transform_matches(flint::FmpzMatConstRef reduced_gram,
                            flint::FmpzMatConstRef original_gram,
                            flint::FmpzMatConstRef transform) noexcept {
    const slong n = flint::fmpz_mat_nrows(original_gram);
    if (n != flint::fmpz_mat_ncols(original_gram) ||
        flint::fmpz_mat_nrows(reduced_gram) != n ||
        flint::fmpz_mat_ncols(reduced_gram) != n ||
        flint::fmpz_mat_nrows(transform) != n ||
        flint::fmpz_mat_ncols(transform) != n) {
        return false;
    }

    flint::FmpzMat transform_transpose(n, n);
    flint::FmpzMat tmp(n, n);
    flint::FmpzMat product(n, n);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transform_transpose),
                              transform);
    flint::fmpz_mat_mul(flint::FmpzMatRef(tmp), transform, original_gram);
    flint::fmpz_mat_mul(flint::FmpzMatRef(product),
                        flint::FmpzMatConstRef(tmp),
                        flint::FmpzMatConstRef(transform_transpose));
    return flint::fmpz_mat_equal(flint::FmpzMatConstRef(product),
                                 reduced_gram);
}

bool ideal_lattice_weight_vector_shape_is_valid(flint::FmpzMatConstRef weights,
                                                slong degree) noexcept {
    return degree > 0 && flint::fmpz_mat_nrows(weights) == 1 &&
           flint::fmpz_mat_ncols(weights) == degree;
}

bool ideal_lattice_weight_vector_is_zero(
        flint::FmpzMatConstRef weights) noexcept {
    if (flint::fmpz_mat_nrows(weights) != 1) {
        return false;
    }
    for (slong column = 0; column < flint::fmpz_mat_ncols(weights);
         ++column) {
        if (!flint::fmpz_is_zero(flint::fmpz_mat_entry(weights, 0, column))) {
            return false;
        }
    }
    return true;
}

static bool ideal_lattice_nonnegative_weight_sum_slong(
        slong& out,
        flint::FmpzMatConstRef weights) noexcept {
    if (flint::fmpz_mat_nrows(weights) != 1) {
        return false;
    }

    flint::Fmpz sum;
    for (slong column = 0; column < flint::fmpz_mat_ncols(weights);
         ++column) {
        flint::fmpz_add(flint::FmpzRef(sum), flint::FmpzConstRef(sum),
                        flint::fmpz_mat_entry(weights, 0, column));
    }
    if (flint::fmpz_sgn(flint::FmpzConstRef(sum)) < 0) {
        out = 0;
        return true;
    }
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(sum))) {
        return false;
    }
    out = flint::fmpz_get_si(flint::FmpzConstRef(sum));
    return true;
}

static bool apply_power_two_column_weights(
        flint::ArbMat& embedding_rows,
        flint::FmpzMatConstRef weights) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (!ideal_lattice_weight_vector_shape_is_valid(weights, cols)) {
        return false;
    }

    for (slong column = 0; column < cols; ++column) {
        flint::FmpzConstRef weight =
                flint::fmpz_mat_entry(weights, 0, column);
        if (!flint::fmpz_fits_si(weight)) {
            return false;
        }
        const slong shift = flint::fmpz_get_si(weight);
        if (shift == 0) {
            continue;
        }
        for (slong row = 0; row < rows; ++row) {
            arb_struct* entry =
                    flint::arb_mat_entry_ref(embedding_rows, row, column)
                            .raw();
            ::arb_mul_2exp_si(entry, entry, shift);
        }
    }
    return true;
}

static flint_bitcnt_t fmpz_mat_max_abs_bits(flint::FmpzMatConstRef matrix) noexcept {
    flint_bitcnt_t max_bits = 0;
    flint::Fmpz absolute_entry;
    for (slong i = 0; i < flint::fmpz_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < flint::fmpz_mat_ncols(matrix); ++j) {
            flint::fmpz_abs(flint::FmpzRef(absolute_entry),
                            flint::fmpz_mat_entry(matrix, i, j));
            max_bits = std::max(
                    max_bits,
                    ::fmpz_bits(flint::FmpzConstRef(absolute_entry).raw()));
        }
    }
    return max_bits;
}

static bool fmpz_mat_diagonal_product(flint::Fmpz& out,
                               flint::FmpzMatConstRef matrix) noexcept {
    const slong degree = flint::fmpz_mat_nrows(matrix);
    if (degree <= 0 || degree != flint::fmpz_mat_ncols(matrix)) {
        return false;
    }

    flint::fmpz_one(flint::FmpzRef(out));
    for (slong i = 0; i < degree; ++i) {
        flint::fmpz_mul(flint::FmpzRef(out), flint::FmpzConstRef(out),
                        flint::fmpz_mat_entry(matrix, i, i));
    }
    return true;
}

static bool ideal_lattice_t2_diagonal_checks_pass(
        flint::FmpzMatConstRef reduced_gram,
        const Ideal& ideal,
        const Order& order,
        slong denominator_bits,
        slong positive_weight_sum) noexcept;

static bool ideal_lattice_low_precision_checks_pass(
        flint::FmpzMatConstRef reduced_gram,
        flint::FmpzMatConstRef gram_transform,
        const Ideal& ideal,
        const Order& order,
        slong denominator_bits,
        slong positive_weight_sum = 0) noexcept {
    const slong degree = order.degree();
    if (degree <= 0 || denominator_bits < 0 || positive_weight_sum < 0 ||
        flint::fmpz_mat_nrows(reduced_gram) != degree ||
        flint::fmpz_mat_ncols(reduced_gram) != degree ||
        flint::fmpz_mat_nrows(gram_transform) != degree ||
        flint::fmpz_mat_ncols(gram_transform) != degree) {
        return false;
    }

    // reference implementation NfOrd/LLL.jl:_lll throws LowPrecisionLLL when the Gram
    // transform grows beyond nbits(maximum(abs, t)) > div(prec, 2), where
    // `prec` is already the denominator exponent after the half shift.
    const flint_bitcnt_t transform_bits =
            fmpz_mat_max_abs_bits(gram_transform);
    if (transform_bits >
        static_cast<flint_bitcnt_t>(denominator_bits / 2)) {
        return false;
    }

    return ideal_lattice_t2_diagonal_checks_pass(
            reduced_gram, ideal, order, denominator_bits,
            positive_weight_sum);
}

// The T2 diagonal part of ideal_lattice_low_precision_checks_pass, shared
// with the rounded-embedding route, which has no Gram transform.
static bool ideal_lattice_t2_diagonal_checks_pass(
        flint::FmpzMatConstRef reduced_gram,
        const Ideal& ideal,
        const Order& order,
        slong denominator_bits,
        slong positive_weight_sum) noexcept {
    const slong degree = order.degree();
    if (degree <= 0 || denominator_bits < 0 || positive_weight_sum < 0 ||
        flint::fmpz_mat_nrows(reduced_gram) != degree ||
        flint::fmpz_mat_ncols(reduced_gram) != degree) {
        return false;
    }

    flint::Fmpz discriminant;
    flint::Fmpz ideal_norm;
    flint::Fmpz source_disc;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        !ideal.norm(flint::FmpzRef(ideal_norm)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(ideal_norm)) <= 0) {
        return false;
    }

    // reference implementation NfOrd/LLL.jl:_lll uses
    //   disc = abs(discriminant(order(A))) * norm(A)^2 *
    //          den^(2*n) * 2^(2*sv)
    // The native weighted path forms the embedded rows directly in Arb, so
    // the rational order-basis denominator is already represented before
    // rounding.  The source `v = 0` path has sv = 0.
    flint::fmpz_abs(flint::FmpzRef(discriminant),
                    flint::FmpzConstRef(discriminant));
    flint::fmpz_mul(flint::FmpzRef(source_disc),
                    flint::FmpzConstRef(ideal_norm),
                    flint::FmpzConstRef(ideal_norm));
    flint::fmpz_mul(flint::FmpzRef(source_disc),
                    flint::FmpzConstRef(source_disc),
                    flint::FmpzConstRef(discriminant));
    if (positive_weight_sum >
        std::numeric_limits<slong>::max() / 2) {
        return false;
    }
    flint::fmpz_mul_2exp(
            flint::FmpzRef(source_disc),
            flint::FmpzConstRef(source_disc),
            static_cast<ulong>(2 * positive_weight_sum));
    if (flint::fmpz_sgn(flint::FmpzConstRef(source_disc)) <= 0) {
        return false;
    }

    flint::Fmpz first_diagonal_bound;
    ::fmpz_root(first_diagonal_bound.raw(), source_disc.raw(),
                static_cast<ulong>(degree));
    flint::fmpz_add_ui(flint::FmpzRef(first_diagonal_bound),
                       flint::FmpzConstRef(first_diagonal_bound),
                       UWORD(1));
    flint::fmpz_mul_2exp(flint::FmpzRef(first_diagonal_bound),
                         flint::FmpzConstRef(first_diagonal_bound),
                         static_cast<ulong>((degree + 1) / 2));
    flint::fmpz_mul_2exp(flint::FmpzRef(first_diagonal_bound),
                         flint::FmpzConstRef(first_diagonal_bound),
                         static_cast<ulong>(denominator_bits));
    if (flint::fmpz_cmp(
                flint::fmpz_mat_entry(reduced_gram, 0, 0),
                flint::FmpzConstRef(first_diagonal_bound)) > 0) {
        return false;
    }

    flint::Fmpz diagonal_product;
    flint::Fmpz product_bound;
    if (!fmpz_mat_diagonal_product(diagonal_product, reduced_gram)) {
        return false;
    }
    flint::fmpz_set(flint::FmpzRef(product_bound),
                    flint::FmpzConstRef(source_disc));
    flint::fmpz_mul_2exp(
            flint::FmpzRef(product_bound),
            flint::FmpzConstRef(product_bound),
            static_cast<ulong>((degree * (degree - 1)) / 2));
    flint::fmpz_mul_2exp(
            flint::FmpzRef(product_bound),
            flint::FmpzConstRef(product_bound),
            static_cast<ulong>(degree * denominator_bits));
    return flint::fmpz_cmp(flint::FmpzConstRef(diagonal_product),
                           flint::FmpzConstRef(product_bound)) <= 0;
}

// Runs certified_fmpz_lll on `matrix` with `transform` set to the
// identity of the row count.
static lat::detail::CertifiedLllRoute ideal_lattice_certified_lll(
        flint::FmpzMat& matrix,
        flint::FmpzMat& transform,
        const flint::FmpzLll& config) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    if (flint::fmpz_mat_nrows(transform) != rows ||
        flint::fmpz_mat_ncols(transform) != rows) {
        transform = flint::FmpzMat(rows, rows);
    }
    flint::fmpz_mat_one(flint::FmpzMatRef(transform));
    return lat::detail::certified_fmpz_lll(matrix.raw(), transform.raw(),
                                           flint::FmpzLllConstRef(config));
}

// The two-LLL route's fmpz_lll calls, through FLINT's certified first
// phase (lat::detail::certified_fmpz_lll, the Lat::lll_reduce route,
// docs/reference/algorithms_and_sources.rst), which skips FLINT's
// multiprecision pre-check of the phase output.  The returned matrices are
// fmpz_lll's on every route.
lat::detail::CertifiedLllRoute ideal_lattice_pre_reduction_lll(
        flint::FmpzMat& basis,
        flint::FmpzMat& transform) noexcept {
    // reference implementation NfOrd/LLL.jl:_lll first applies
    // lll_with_transform(basis_matrix(A), LLLContext(0.5, 0.51)).
    const flint::FmpzLll config(0.5, 0.51, Z_BASIS, APPROX);
    return ideal_lattice_certified_lll(basis, transform, config);
}

lat::detail::CertifiedLllRoute ideal_lattice_gram_lll(
        flint::FmpzMat& gram,
        flint::FmpzMat& transform) noexcept {
    const flint::FmpzLll config(0.99, 0.51, GRAM, APPROX);
    return ideal_lattice_certified_lll(gram, transform, config);
}

static IdealLatticeReductionStatus
ideal_lattice_gram_lll_reduce(flint::FmpzMat& reduced_basis,
                              flint::FmpzMat& transform,
                              flint::FmpzMat& reduced_gram,
                              const flint::FmpzMat& basis,
                              const Ideal& ideal,
                              const Order& order,
                              slong embedding_precision,
                              slong precision,
                              slong& embedding_shortfall_bits) noexcept;

static IdealLatticeReductionStatus
weighted_ideal_lattice_gram_lll_reduce(flint::FmpzMat& reduced_basis,
                                       flint::FmpzMat& transform,
                                       flint::FmpzMat& reduced_gram,
                                       const flint::FmpzMat& basis,
                                       const Ideal& ideal,
                                       const Order& order,
                                       flint::FmpzMatConstRef weights,
                                       slong embedding_precision,
                                       slong precision,
                                       slong& embedding_shortfall_bits,
                                       const DiagnosticsContext* diagnostics,
                                       detail::OrderMinkowskiEmbeddingCache* cache)
        noexcept {
    const slong degree = flint::fmpz_mat_nrows(basis);
    if (precision <= 0 || degree != flint::fmpz_mat_ncols(basis) ||
        embedding_precision <= 0 || order.degree() != degree ||
        !ideal_lattice_weight_vector_shape_is_valid(weights, degree) ||
        flint::fmpz_mat_nrows(reduced_basis) != degree ||
        flint::fmpz_mat_ncols(reduced_basis) != degree ||
        flint::fmpz_mat_nrows(transform) != degree ||
        flint::fmpz_mat_ncols(transform) != degree ||
        flint::fmpz_mat_nrows(reduced_gram) != degree ||
        flint::fmpz_mat_ncols(reduced_gram) != degree) {
        return IdealLatticeReductionStatus::failed;
    }
    if (ideal_lattice_weight_vector_is_zero(weights)) {
        return ideal_lattice_gram_lll_reduce(
                reduced_basis, transform, reduced_gram, basis, ideal, order,
                embedding_precision, precision, embedding_shortfall_bits);
    }

    flint::FmpzMat pre_reduced_basis(degree, degree);
    flint::FmpzMat pre_transform(degree, degree);
    flint::FmpzMat pre_transform_check(degree, degree);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.pre_lll");
        flint::fmpz_mat_set(flint::FmpzMatRef(pre_reduced_basis),
                            flint::FmpzMatConstRef(basis));
        // reference implementation NfOrd/LLL.jl:_lll first applies
        // lll_with_transform(basis_matrix(A), LLLContext(0.5, 0.51)).
        (void)ideal_lattice_pre_reduction_lll(pre_reduced_basis,
                                              pre_transform);
        flint::fmpz_mat_mul(flint::FmpzMatRef(pre_transform_check),
                            flint::FmpzMatConstRef(pre_transform),
                            flint::FmpzMatConstRef(basis));
        if (!flint::fmpz_mat_equal(
                    flint::FmpzMatConstRef(pre_transform_check),
                    flint::FmpzMatConstRef(pre_reduced_basis))) {
            return IdealLatticeReductionStatus::failed;
        }
    }

    slong positive_weight_sum = 0;
    if (!ideal_lattice_nonnegative_weight_sum_slong(positive_weight_sum,
                                                    weights)) {
        return IdealLatticeReductionStatus::failed;
    }

    const slong lll_denominator_bits = precision / 2;
    flint::ArbMat weighted_embedding_rows(degree, degree);
    flint::FmpzMat original_gram(degree, degree);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.embedding_rows");
        if (!detail::build_ideal_minkowski_embedding_rows(
                    weighted_embedding_rows,
                    flint::FmpzMatConstRef(pre_reduced_basis), order,
                    embedding_precision, cache, diagnostics)) {
            return IdealLatticeReductionStatus::retry_precision;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.apply_weights");
        if (!apply_power_two_column_weights(weighted_embedding_rows,
                                            weights)) {
            return IdealLatticeReductionStatus::retry_precision;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.gram_round");
        if (!ideal_lattice_scaled_minkowski_gram_numerator_without_diagonal(
                    original_gram, weighted_embedding_rows, precision,
                    precision)) {
            return ideal_lattice_rounding_retry_status(
                    embedding_shortfall_bits, weighted_embedding_rows,
                    precision, precision);
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.scale_diagonal");
        flint::fmpz_mat_scalar_tdiv_q_2exp(
                flint::FmpzMatRef(original_gram),
                flint::FmpzMatConstRef(original_gram),
                static_cast<ulong>(lll_denominator_bits));
        for (slong i = 0; i < degree; ++i) {
            flint::FmpzRef diagonal =
                    flint::fmpz_mat_entry(original_gram, i, i);
            flint::fmpz_add_ui(diagonal,
                               flint::FmpzConstRef(diagonal.raw()),
                               static_cast<ulong>(degree));
        }
    }

    flint::FmpzMat gram_transform(degree, degree);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.gram_lll");
        flint::fmpz_mat_set(flint::FmpzMatRef(reduced_gram),
                            flint::FmpzMatConstRef(original_gram));
        (void)ideal_lattice_gram_lll(reduced_gram, gram_transform);

        if (!gram_transform_matches(flint::FmpzMatConstRef(reduced_gram),
                                    flint::FmpzMatConstRef(original_gram),
                                    flint::FmpzMatConstRef(gram_transform))) {
            return IdealLatticeReductionStatus::failed;
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.precision_check");
        if (!ideal_lattice_low_precision_checks_pass(
                    flint::FmpzMatConstRef(reduced_gram),
                    flint::FmpzMatConstRef(gram_transform), ideal, order,
                    lll_denominator_bits, positive_weight_sum)) {
            return IdealLatticeReductionStatus::retry_precision;
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.weighted.combine");
        flint::fmpz_mat_mul(flint::FmpzMatRef(transform),
                            flint::FmpzMatConstRef(gram_transform),
                            flint::FmpzMatConstRef(pre_transform));
        flint::fmpz_mat_mul(flint::FmpzMatRef(reduced_basis),
                            flint::FmpzMatConstRef(gram_transform),
                            flint::FmpzMatConstRef(pre_reduced_basis));
        flint::fmpz_mat_mul(flint::FmpzMatRef(pre_transform_check),
                            flint::FmpzMatConstRef(transform),
                            flint::FmpzMatConstRef(basis));
        if (!flint::fmpz_mat_equal(
                    flint::FmpzMatConstRef(pre_transform_check),
                    flint::FmpzMatConstRef(reduced_basis))) {
            return IdealLatticeReductionStatus::failed;
        }
    }

    return IdealLatticeReductionStatus::success;
}

static IdealLatticeReductionStatus
ideal_lattice_gram_lll_reduce(flint::FmpzMat& reduced_basis,
                              flint::FmpzMat& transform,
                              flint::FmpzMat& reduced_gram,
                              const flint::FmpzMat& basis,
                              const Ideal& ideal,
                              const Order& order,
                              slong embedding_precision,
                              slong precision,
                              slong& embedding_shortfall_bits) noexcept {
    const slong degree = flint::fmpz_mat_nrows(basis);
    if (precision <= 0 || degree != flint::fmpz_mat_ncols(basis) ||
        embedding_precision <= 0 || order.degree() != degree ||
        flint::fmpz_mat_nrows(reduced_basis) != degree ||
        flint::fmpz_mat_ncols(reduced_basis) != degree ||
        flint::fmpz_mat_nrows(transform) != degree ||
        flint::fmpz_mat_ncols(transform) != degree ||
        flint::fmpz_mat_nrows(reduced_gram) != degree ||
        flint::fmpz_mat_ncols(reduced_gram) != degree) {
        return IdealLatticeReductionStatus::failed;
    }

    flint::FmpzMat pre_reduced_basis(degree, degree);
    flint::FmpzMat pre_transform(degree, degree);
    flint::FmpzMat pre_transform_check(degree, degree);
    flint::fmpz_mat_set(flint::FmpzMatRef(pre_reduced_basis),
                        flint::FmpzMatConstRef(basis));
    // reference implementation NfOrd/LLL.jl:_lll first applies
    // lll_with_transform(basis_matrix(A), LLLContext(0.5, 0.51)).
    (void)ideal_lattice_pre_reduction_lll(pre_reduced_basis, pre_transform);
    flint::fmpz_mat_mul(flint::FmpzMatRef(pre_transform_check),
                        flint::FmpzMatConstRef(pre_transform),
                        flint::FmpzMatConstRef(basis));
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(pre_transform_check),
                               flint::FmpzMatConstRef(pre_reduced_basis))) {
        return IdealLatticeReductionStatus::failed;
    }

    const slong lll_denominator_bits = precision / 2;
    flint::FmpzMat order_basis(degree, degree);
    flint::ArbMat order_embedding_rows(degree, degree);
    flint::FmpzMat order_gram(degree, degree);
    flint::FmpzMat original_gram(degree, degree);
    flint::FmpzMat pre_reduced_basis_transpose(degree, degree);
    flint::FmpzMat tmp_gram(degree, degree);
    flint::fmpz_mat_one(flint::FmpzMatRef(order_basis));
    if (!detail::build_ideal_minkowski_embedding_rows(
                order_embedding_rows,
                flint::FmpzMatConstRef(order_basis), order,
                embedding_precision)) {
        return IdealLatticeReductionStatus::retry_precision;
    }
    if (!ideal_lattice_scaled_minkowski_gram_numerator(
                order_gram, order_embedding_rows, precision, precision)) {
        return ideal_lattice_rounding_retry_status(
                embedding_shortfall_bits, order_embedding_rows, precision,
                precision);
    }
    flint::fmpz_mat_transpose(
            flint::FmpzMatRef(pre_reduced_basis_transpose),
            flint::FmpzMatConstRef(pre_reduced_basis));
    flint::fmpz_mat_mul(flint::FmpzMatRef(tmp_gram),
                        flint::FmpzMatConstRef(order_gram),
                        flint::FmpzMatConstRef(
                                pre_reduced_basis_transpose));
    flint::fmpz_mat_mul(flint::FmpzMatRef(original_gram),
                        flint::FmpzMatConstRef(pre_reduced_basis),
                        flint::FmpzMatConstRef(tmp_gram));
    flint::fmpz_mat_scalar_tdiv_q_2exp(
            flint::FmpzMatRef(original_gram),
            flint::FmpzMatConstRef(original_gram),
            static_cast<ulong>(lll_denominator_bits));
    for (slong i = 0; i < degree; ++i) {
        flint::FmpzRef diagonal =
                flint::fmpz_mat_entry(original_gram, i, i);
        flint::fmpz_add_ui(diagonal,
                           flint::FmpzConstRef(diagonal.raw()),
                           static_cast<ulong>(degree));
    }

    flint::fmpz_mat_set(flint::FmpzMatRef(reduced_gram),
                        flint::FmpzMatConstRef(original_gram));
    flint::FmpzMat gram_transform(degree, degree);
    (void)ideal_lattice_gram_lll(reduced_gram, gram_transform);

    if (!gram_transform_matches(flint::FmpzMatConstRef(reduced_gram),
                                flint::FmpzMatConstRef(original_gram),
                                flint::FmpzMatConstRef(gram_transform))) {
        return IdealLatticeReductionStatus::failed;
    }

    if (!ideal_lattice_low_precision_checks_pass(
                flint::FmpzMatConstRef(reduced_gram),
                flint::FmpzMatConstRef(gram_transform), ideal, order,
                lll_denominator_bits)) {
        return IdealLatticeReductionStatus::retry_precision;
    }

    flint::fmpz_mat_mul(flint::FmpzMatRef(transform),
                        flint::FmpzMatConstRef(gram_transform),
                        flint::FmpzMatConstRef(pre_transform));
    flint::fmpz_mat_mul(flint::FmpzMatRef(reduced_basis),
                        flint::FmpzMatConstRef(gram_transform),
                        flint::FmpzMatConstRef(pre_reduced_basis));
    flint::fmpz_mat_mul(flint::FmpzMatRef(pre_transform_check),
                        flint::FmpzMatConstRef(transform),
                        flint::FmpzMatConstRef(basis));
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(pre_transform_check),
                               flint::FmpzMatConstRef(reduced_basis))) {
        return IdealLatticeReductionStatus::failed;
    }

    return IdealLatticeReductionStatus::success;
}
bool copy_ideal_basis_row_coordinates(flint::FmpzMat& out,
                                      flint::FmpzMatConstRef basis,
                                      slong row) noexcept {
    if (flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(basis) ||
        row < 0 || row >= flint::fmpz_mat_nrows(basis)) {
        return false;
    }

    for (slong column = 0; column < flint::fmpz_mat_ncols(basis); ++column) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, 0, column),
                        flint::fmpz_mat_entry(basis, row, column));
    }
    return true;
}
bool build_ideal_lattice_lll_data(
        IdealLatticeLllData& out,
        const Ideal& ideal,
        slong precision) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || precision <= 0) {
        return false;
    }

    const slong degree = order->degree();
    flint::FmpzMat original_basis(degree, degree);
    if (!ideal.get_hnf(flint::FmpzMatRef(original_basis))) {
        return false;
    }

    if (degree > kIdealLatticeLllMaxPrecision / 4) {
        return false;
    }
    slong active_precision = std::max(precision, 4 * degree);
    if (active_precision > kIdealLatticeLllMaxPrecision) {
        return false;
    }

    IdealLatticePrecisionSchedule schedule(active_precision);
    while (schedule.active_precision() <= kIdealLatticeLllMaxPrecision) {
        active_precision = schedule.active_precision();
        flint::FmpzMat basis(degree, degree);
        flint::FmpzMat transform(degree, degree);
        flint::FmpzMat scaled_gram(degree, degree);
        flint::Fmpz gram_denominator;

        slong embedding_shortfall_bits = 0;
        const IdealLatticeReductionStatus reduction_status =
                ideal_lattice_gram_lll_reduce(
                        basis, transform, scaled_gram, original_basis, ideal,
                        *order, schedule.embedding_precision(),
                        active_precision, embedding_shortfall_bits);
        if (reduction_status == IdealLatticeReductionStatus::retry_precision ||
            reduction_status ==
                    IdealLatticeReductionStatus::retry_embedding_precision) {
            if (!schedule.advance(reduction_status,
                                  embedding_shortfall_bits)) {
                return false;
            }
            continue;
        }
        if (reduction_status != IdealLatticeReductionStatus::success) {
            return false;
        }

        flint::fmpz_one(flint::FmpzRef(gram_denominator));
        flint::fmpz_mul_2exp(flint::FmpzRef(gram_denominator),
                             flint::FmpzConstRef(gram_denominator),
                             static_cast<ulong>(active_precision / 2));

        out.basis = std::move(basis);
        out.scaled_gram = std::move(scaled_gram);
        out.gram_denominator = std::move(gram_denominator);
        return true;
    }

    return false;
}

// The reference nf_get_roundG (2.17.3) is the matrix stored by
// base1.c:1851, RM_round_maxrank(G) (base4.c:2631-2640): round 2^e G for
// e = 4, 8, 16, ... until the rounded matrix has full rank.  Silex rounds
// the midpoints of the order's Minkowski rows (rows indexed by the order
// basis, the transpose of the reference G), so the rounded rows R satisfy B R ~ 2^e B E for
// an ideal basis B in order coordinates.  Any integer matrix is a valid LLL
// input here; R only steers the reduction, and the caller checks the
// result exactly.
static bool ideal_lattice_round_embedding_rows_full_rank(
        flint::FmpzMat& out,
        const flint::ArbMat& embedding_rows,
        slong max_scale) noexcept {
    const slong rows = flint::arb_mat_nrows_value(embedding_rows);
    const slong cols = flint::arb_mat_ncols_value(embedding_rows);
    if (rows <= 0 || rows != cols || max_scale < 4) {
        return false;
    }

    flint::FmpzMat rounded(rows, cols);
    flint::Arf midpoint;
    for (slong scale = 4; scale <= max_scale; scale *= 2) {
        for (slong i = 0; i < rows; ++i) {
            for (slong j = 0; j < cols; ++j) {
                flint::arf_set(
                        midpoint,
                        arb_midref(flint::arb_mat_entry_ref(embedding_rows, i,
                                                            j)
                                           .raw()));
                if (!flint::arf_is_finite(midpoint)) {
                    return false;
                }
                ::arf_mul_2exp_si(midpoint.raw(), midpoint.raw(), scale);
                flint::arf_get_fmpz(flint::fmpz_mat_entry(rounded, i, j),
                                    midpoint, ARF_RND_NEAR);
            }
        }
        if (flint::fmpz_mat_rank(flint::FmpzMatConstRef(rounded)) == rows) {
            out = std::move(rounded);
            return true;
        }
    }
    return false;
}

// Fills the per-order cache: the order T2 Gram numerator follows the same
// precision schedule as the order Gram of ideal_lattice_gram_lll_reduce, so
// at a given starting precision both routes use the same order Gram.  The
// rounded rows come from the first embedding whose Gram rounds.
static bool ideal_lattice_fill_rounded_embedding_cache(
        IdealLatticeRoundedEmbeddingCache& cache,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = order.degree();
    if (degree <= 0 || precision <= 0) {
        return false;
    }
    if (!cache.has_order || !cache.order.has_same_data(order)) {
        cache.order = order;
        cache.has_order = true;
        cache.has_rounded_rows = false;
        cache.rounded_rows = flint::FmpzMat(0, 0);
        cache.order_gram_requested_precision = 0;
        cache.order_gram_precision = 0;
        cache.order_gram = flint::FmpzMat(0, 0);
        cache.route_checked = false;
        cache.use_two_lll_route = false;
    }
    if (cache.has_rounded_rows &&
        cache.order_gram_requested_precision == precision) {
        return true;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.ideal_lattice_reduction.rounded_embedding.fill_cache");
    flint::FmpzMat order_basis(degree, degree);
    flint::fmpz_mat_one(flint::FmpzMatRef(order_basis));
    IdealLatticePrecisionSchedule schedule(precision);
    while (schedule.active_precision() <= kIdealLatticeLllMaxPrecision) {
        const slong active_precision = schedule.active_precision();
        flint::ArbMat order_embedding_rows(degree, degree);
        flint::FmpzMat order_gram(degree, degree);
        IdealLatticeReductionStatus status =
                IdealLatticeReductionStatus::retry_precision;
        slong embedding_shortfall_bits = 0;
        if (detail::build_ideal_minkowski_embedding_rows(
                    order_embedding_rows,
                    flint::FmpzMatConstRef(order_basis), order,
                    schedule.embedding_precision())) {
            if (ideal_lattice_scaled_minkowski_gram_numerator(
                        order_gram, order_embedding_rows, active_precision,
                        active_precision)) {
                status = IdealLatticeReductionStatus::success;
            } else {
                status = ideal_lattice_rounding_retry_status(
                        embedding_shortfall_bits, order_embedding_rows,
                        active_precision, active_precision);
            }
        }
        if (status != IdealLatticeReductionStatus::success) {
            if (!schedule.advance(status, embedding_shortfall_bits)) {
                return false;
            }
            continue;
        }

        if (!cache.has_rounded_rows) {
            flint::FmpzMat rounded_rows(0, 0);
            if (!ideal_lattice_round_embedding_rows_full_rank(
                        rounded_rows, order_embedding_rows,
                        schedule.embedding_precision())) {
                return false;
            }
            cache.rounded_rows = std::move(rounded_rows);
            cache.has_rounded_rows = true;
        }
        cache.order_gram = std::move(order_gram);
        cache.order_gram_precision = active_precision;
        cache.order_gram_requested_precision = precision;
        return true;
    }
    return false;
}

// The integer LLL of the rounded-embedding route, with fmpz_lll's result,
// through FLINT's certified first phase
// (lat::detail::certified_fmpz_lll, the Lat::lll_reduce route,
// docs/reference/algorithms_and_sources.rst), which skips FLINT's
// inconclusive multiprecision pre-check that dominated the cost of this
// LLL.
lat::detail::CertifiedLllRoute ideal_lattice_rounded_embedding_integer_lll(
        flint::FmpzMat& embedded,
        flint::FmpzMat& transform) noexcept {
    const flint::FmpzLll config(0.99, 0.51, Z_BASIS, APPROX);
    return ideal_lattice_certified_lll(embedded, transform, config);
}

bool ideal_lattice_symmetric_congruence_product(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef basis,
        flint::FmpzMatConstRef gram) noexcept {
    const slong rows = flint::fmpz_mat_nrows(basis);
    const slong cols = flint::fmpz_mat_ncols(basis);
    if (flint::fmpz_mat_nrows(gram) != cols ||
        flint::fmpz_mat_ncols(gram) != cols ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != rows) {
        return false;
    }
    if (rows == 0) {
        return true;
    }
    if (cols == 0) {
        flint::fmpz_mat_zero(flint::FmpzMatRef(out));
        return true;
    }
    flint::FmpzMat basis_gram(rows, cols);
    flint::fmpz_mat_mul(flint::FmpzMatRef(basis_gram), basis, gram);
    const fmpz_mat_struct* b = basis.raw();
    const fmpz_mat_struct* c = basis_gram.raw();
    fmpz_mat_struct* o = out.raw();
    for (slong i = 0; i < rows; ++i) {
        for (slong j = i; j < rows; ++j) {
            ::_fmpz_vec_dot(::fmpz_mat_entry(o, i, j),
                            ::fmpz_mat_entry(b, i, 0),
                            ::fmpz_mat_entry(c, j, 0), cols);
        }
        for (slong j = 0; j < i; ++j) {
            ::fmpz_set(::fmpz_mat_entry(o, i, j), ::fmpz_mat_entry(o, j, i));
        }
    }
    return true;
}

// The reference 2.17.3 buch2.c:2461-2462 (Fincke_Pohst_ideal):
//   u = ZM_lll(ZM_mul(G0, I), 0.99, LLL_IM); ideal = ZM_mul(I, u);
// one integer LLL of the ideal basis in the rounded embedding.  The
// reduced basis is checked exactly against the LLL output, and its scaled
// T2 Gram is formed exactly from the cached order Gram as in
// ideal_lattice_gram_lll_reduce.  Returns retry_precision when the T2
// diagonal checks reject the result.
static IdealLatticeReductionStatus ideal_lattice_rounded_embedding_lll_reduce(
        flint::FmpzMat& reduced_basis,
        flint::FmpzMat& reduced_gram,
        const flint::FmpzMat& basis,
        const Ideal& ideal,
        const Order& order,
        const IdealLatticeRoundedEmbeddingCache& cache,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = flint::fmpz_mat_nrows(basis);
    if (degree <= 0 || degree != flint::fmpz_mat_ncols(basis) ||
        order.degree() != degree || !cache.has_rounded_rows ||
        cache.order_gram_precision <= 0 ||
        flint::fmpz_mat_nrows(cache.rounded_rows) != degree ||
        flint::fmpz_mat_ncols(cache.rounded_rows) != degree ||
        flint::fmpz_mat_nrows(cache.order_gram) != degree ||
        flint::fmpz_mat_ncols(cache.order_gram) != degree ||
        flint::fmpz_mat_nrows(reduced_basis) != degree ||
        flint::fmpz_mat_ncols(reduced_basis) != degree ||
        flint::fmpz_mat_nrows(reduced_gram) != degree ||
        flint::fmpz_mat_ncols(reduced_gram) != degree) {
        return IdealLatticeReductionStatus::failed;
    }

    flint::FmpzMat embedded(degree, degree);
    flint::FmpzMat transform(degree, degree);
    flint::FmpzMat check(degree, degree);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.rounded_embedding.lll");
        flint::fmpz_mat_mul(flint::FmpzMatRef(embedded),
                            flint::FmpzMatConstRef(basis),
                            flint::FmpzMatConstRef(cache.rounded_rows));
        ideal_lattice_rounded_embedding_integer_lll(embedded, transform);
        flint::fmpz_mat_mul(flint::FmpzMatRef(reduced_basis),
                            flint::FmpzMatConstRef(transform),
                            flint::FmpzMatConstRef(basis));
        flint::fmpz_mat_mul(flint::FmpzMatRef(check),
                            flint::FmpzMatConstRef(reduced_basis),
                            flint::FmpzMatConstRef(cache.rounded_rows));
        if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(check),
                                   flint::FmpzMatConstRef(embedded))) {
            return IdealLatticeReductionStatus::failed;
        }
    }

    const slong lll_denominator_bits = cache.order_gram_precision / 2;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.rounded_embedding.gram");
        // The cached order Gram is symmetric (fmpz_mat_gram, an entrywise
        // truncation and a diagonal shift), so the exact product
        // B * (G * B^T) is symmetric and its upper triangle determines it.
        if (!ideal_lattice_symmetric_congruence_product(
                    reduced_gram, flint::FmpzMatConstRef(reduced_basis),
                    flint::FmpzMatConstRef(cache.order_gram))) {
            return IdealLatticeReductionStatus::failed;
        }
        flint::fmpz_mat_scalar_tdiv_q_2exp(
                flint::FmpzMatRef(reduced_gram),
                flint::FmpzMatConstRef(reduced_gram),
                static_cast<ulong>(lll_denominator_bits));
        for (slong i = 0; i < degree; ++i) {
            flint::FmpzRef diagonal =
                    flint::fmpz_mat_entry(reduced_gram, i, i);
            flint::fmpz_add_ui(diagonal,
                               flint::FmpzConstRef(diagonal.raw()),
                               static_cast<ulong>(degree));
        }
    }

    if (!ideal_lattice_t2_diagonal_checks_pass(
                flint::FmpzMatConstRef(reduced_gram), ideal, order,
                lll_denominator_bits, 0)) {
        return IdealLatticeReductionStatus::retry_precision;
    }
    return IdealLatticeReductionStatus::success;
}

// Silex quality gate, not reference behavior.  The reference rounds the
// embeddings of a T2-reduced integral basis; Silex rounds those of its
// order basis, which may be badly conditioned (for example a power basis
// whose rows span many binary orders of magnitude).  Rounding at the first
// full-rank scale then loses the short vectors, and the one-LLL basis,
// though it passes the T2 diagonal bounds, is much longer than the two-LLL
// basis, so the relation search finds few relations per ideal.  Both are
// 0.99-LLL bases of the same lattice, but LLL does not bound the ratio of
// their T2 diagonal products by 2^n: the factor 2^degree is an empirical
// Silex tolerance with no theoretical backing.  Returns true when the
// one-LLL product (denominator 2^one_denominator_bits) is at most
// 2^degree times the two-LLL product, i.e. at most a factor 2 per T2
// diagonal entry on average.  Returns false on malformed inputs (a
// non-power-of-two denominator, an overflow guard, or a failed diagonal
// product); the caller then uses the two-LLL route, which is always valid.
bool ideal_lattice_rounded_embedding_basis_matches_two_lll(
        flint::FmpzMatConstRef one_lll_gram,
        slong one_denominator_bits,
        const IdealLatticeLllData& two_lll) noexcept {
    const slong degree = flint::fmpz_mat_nrows(one_lll_gram);
    if (degree <= 0 || one_denominator_bits < 0 ||
        flint::fmpz_mat_nrows(two_lll.scaled_gram) != degree ||
        flint::fmpz_sgn(flint::FmpzConstRef(two_lll.gram_denominator)) <= 0) {
        return false;
    }
    const flint::FmpzConstRef two_denominator(two_lll.gram_denominator);
    const flint_bitcnt_t two_denominator_bits =
            ::fmpz_val2(two_denominator.raw());
    if (::fmpz_bits(two_denominator.raw()) != two_denominator_bits + 1 ||
        one_denominator_bits >
                std::numeric_limits<slong>::max() / (degree + 1) ||
        static_cast<ulong>(two_denominator_bits) >
                static_cast<ulong>(std::numeric_limits<slong>::max() /
                                   (degree + 1))) {
        return false;
    }

    flint::Fmpz one_product;
    flint::Fmpz two_product;
    if (!fmpz_mat_diagonal_product(one_product, one_lll_gram) ||
        !fmpz_mat_diagonal_product(
                two_product,
                flint::FmpzMatConstRef(two_lll.scaled_gram))) {
        return false;
    }
    // one_product / 2^(degree d1) <= 2^degree two_product / 2^(degree d2).
    flint::fmpz_mul_2exp(
            flint::FmpzRef(one_product), flint::FmpzConstRef(one_product),
            static_cast<ulong>(degree) * two_denominator_bits);
    flint::fmpz_mul_2exp(
            flint::FmpzRef(two_product), flint::FmpzConstRef(two_product),
            static_cast<ulong>(degree) *
                    static_cast<ulong>(one_denominator_bits + 1));
    return flint::fmpz_cmp(flint::FmpzConstRef(one_product),
                           flint::FmpzConstRef(two_product)) <= 0;
}

bool build_ideal_lattice_rounded_embedding_lll_data(
        IdealLatticeLllData& out,
        const Ideal& ideal,
        slong precision,
        IdealLatticeRoundedEmbeddingCache& cache,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || precision <= 0) {
        return false;
    }

    const slong degree = order->degree();
    if (degree <= 0 || degree > kIdealLatticeLllMaxPrecision / 4) {
        return false;
    }
    const slong active_precision = std::max(precision, 4 * degree);
    if (active_precision > kIdealLatticeLllMaxPrecision) {
        return false;
    }

    flint::FmpzMat original_basis(degree, degree);
    if (!ideal.get_hnf(flint::FmpzMatRef(original_basis))) {
        return false;
    }

    if (cache.use_two_lll_route && cache.has_order &&
        cache.order.has_same_data(*order)) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.rounded_embedding."
                "two_lll_order");
        return build_ideal_lattice_lll_data(out, ideal, precision);
    }

    if (ideal_lattice_fill_rounded_embedding_cache(cache, *order,
                                                   active_precision,
                                                   diagnostics)) {
        flint::FmpzMat basis(degree, degree);
        flint::FmpzMat scaled_gram(degree, degree);
        if (ideal_lattice_rounded_embedding_lll_reduce(
                    basis, scaled_gram, original_basis, ideal, *order, cache,
                    diagnostics) == IdealLatticeReductionStatus::success) {
            if (!cache.route_checked) {
                // The first accepted one-LLL basis of the order is compared
                // with the two-LLL basis of the same ideal; when that fails,
                // this and every later ideal of the order take the two-LLL
                // route.  When the two-LLL route fails, the check waits for
                // the next ideal.
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction."
                        "rounded_embedding.route_check");
                IdealLatticeLllData two_lll_data;
                if (build_ideal_lattice_lll_data(two_lll_data, ideal,
                                                 precision)) {
                    cache.route_checked = true;
                    if (!ideal_lattice_rounded_embedding_basis_matches_two_lll(
                                flint::FmpzMatConstRef(scaled_gram),
                                cache.order_gram_precision / 2,
                                two_lll_data)) {
                        cache.use_two_lll_route = true;
                        out = std::move(two_lll_data);
                        return true;
                    }
                }
            }
            flint::Fmpz gram_denominator;
            flint::fmpz_one(flint::FmpzRef(gram_denominator));
            flint::fmpz_mul_2exp(
                    flint::FmpzRef(gram_denominator),
                    flint::FmpzConstRef(gram_denominator),
                    static_cast<ulong>(cache.order_gram_precision / 2));
            out.basis = std::move(basis);
            out.scaled_gram = std::move(scaled_gram);
            out.gram_denominator = std::move(gram_denominator);
            return true;
        }
    }

    // The T2 checks rejected the one-LLL basis, or the order data could not
    // be formed: use the two-LLL route with its own precision schedule.
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.ideal_lattice_reduction.rounded_embedding.fallback");
    return build_ideal_lattice_lll_data(out, ideal, precision);
}

bool build_weighted_ideal_lattice_reduction(
        flint::FmpzMat& basis,
        flint::FmpzMat& transform,
        flint::FmpzMat& scaled_gram,
        flint::Fmpz& gram_denominator,
        const Ideal& ideal,
        flint::FmpzMatConstRef weights,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr,
        detail::OrderMinkowskiEmbeddingCache* cache = nullptr) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || precision <= 0) {
        return false;
    }

    const slong degree = order->degree();
    if (!ideal_lattice_weight_vector_shape_is_valid(weights, degree) ||
        flint::fmpz_mat_nrows(basis) != degree ||
        flint::fmpz_mat_ncols(basis) != degree ||
        flint::fmpz_mat_nrows(transform) != degree ||
        flint::fmpz_mat_ncols(transform) != degree ||
        flint::fmpz_mat_nrows(scaled_gram) != degree ||
        flint::fmpz_mat_ncols(scaled_gram) != degree) {
        return false;
    }

    flint::FmpzMat original_basis(degree, degree);
    if (!ideal.get_hnf(flint::FmpzMatRef(original_basis))) {
        return false;
    }

    if (degree > kIdealLatticeLllMaxPrecision / 4) {
        return false;
    }
    slong active_precision = std::max(precision, 4 * degree);
    if (active_precision > kIdealLatticeLllMaxPrecision) {
        return false;
    }

    IdealLatticePrecisionSchedule schedule(active_precision);
    while (schedule.active_precision() <= kIdealLatticeLllMaxPrecision) {
        active_precision = schedule.active_precision();
        slong embedding_shortfall_bits = 0;
        const IdealLatticeReductionStatus reduction_status =
                weighted_ideal_lattice_gram_lll_reduce(
                        basis, transform, scaled_gram, original_basis, ideal,
                        *order, weights, schedule.embedding_precision(),
                        active_precision, embedding_shortfall_bits,
                        diagnostics, cache);
        if (reduction_status == IdealLatticeReductionStatus::retry_precision ||
            reduction_status ==
                    IdealLatticeReductionStatus::retry_embedding_precision) {
            if (!schedule.advance(reduction_status,
                                  embedding_shortfall_bits)) {
                return false;
            }
            continue;
        }
        if (reduction_status != IdealLatticeReductionStatus::success) {
            return false;
        }

        flint::fmpz_one(flint::FmpzRef(gram_denominator));
        flint::fmpz_mul_2exp(flint::FmpzRef(gram_denominator),
                             flint::FmpzConstRef(gram_denominator),
                             static_cast<ulong>(active_precision / 2));
        return true;
    }

    return false;
}

}  // namespace detail::relation_search
}  // namespace silex
