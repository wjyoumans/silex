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

enum class IdealLatticeReductionStatus {
    success,
    retry_precision,
    failed
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

static IdealLatticeReductionStatus
ideal_lattice_gram_lll_reduce(flint::FmpzMat& reduced_basis,
                              flint::FmpzMat& transform,
                              flint::FmpzMat& reduced_gram,
                              const flint::FmpzMat& basis,
                              const Ideal& ideal,
                              const Order& order,
                              slong embedding_precision,
                              slong precision) noexcept;

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
                embedding_precision, precision);
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
        flint::fmpz_mat_one(flint::FmpzMatRef(pre_transform));

        // reference implementation NfOrd/LLL.jl:_lll first applies
        // lll_with_transform(basis_matrix(A), LLLContext(0.5, 0.51)).
        flint::FmpzLll pre_config(0.5, 0.51, Z_BASIS, APPROX);
        fmpz_lll(pre_reduced_basis.raw(), pre_transform.raw(),
                 pre_config.raw());
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
            return IdealLatticeReductionStatus::retry_precision;
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
        flint::fmpz_mat_one(flint::FmpzMatRef(gram_transform));
        flint::FmpzLll config(0.99, 0.51, GRAM, APPROX);
        fmpz_lll(reduced_gram.raw(), gram_transform.raw(), config.raw());

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
                              slong precision) noexcept {
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
    flint::fmpz_mat_one(flint::FmpzMatRef(pre_transform));

    // reference implementation NfOrd/LLL.jl:_lll first applies
    // lll_with_transform(basis_matrix(A), LLLContext(0.5, 0.51)).
    flint::FmpzLll pre_config(0.5, 0.51, Z_BASIS, APPROX);
    fmpz_lll(pre_reduced_basis.raw(), pre_transform.raw(),
             pre_config.raw());
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
                embedding_precision) ||
        !ideal_lattice_scaled_minkowski_gram_numerator(
                order_gram, order_embedding_rows, precision, precision)) {
        return IdealLatticeReductionStatus::retry_precision;
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
    flint::fmpz_mat_one(flint::FmpzMatRef(gram_transform));
    flint::FmpzLll config(0.99, 0.51, GRAM, APPROX);
    fmpz_lll(reduced_gram.raw(), gram_transform.raw(), config.raw());

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

    while (active_precision <= kIdealLatticeLllMaxPrecision) {
        const slong embedding_precision =
                active_precision + kIdealLatticeRoundingGuardBits;
        flint::FmpzMat basis(degree, degree);
        flint::FmpzMat transform(degree, degree);
        flint::FmpzMat scaled_gram(degree, degree);
        flint::Fmpz gram_denominator;

        const IdealLatticeReductionStatus reduction_status =
                ideal_lattice_gram_lll_reduce(
                        basis, transform, scaled_gram, original_basis, ideal,
                        *order, embedding_precision, active_precision);
        if (reduction_status == IdealLatticeReductionStatus::retry_precision) {
            if (active_precision > kIdealLatticeLllMaxPrecision / 2) {
                return false;
            }
            active_precision *= 2;
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

    while (active_precision <= kIdealLatticeLllMaxPrecision) {
        const slong embedding_precision =
                active_precision + kIdealLatticeRoundingGuardBits;
        const IdealLatticeReductionStatus reduction_status =
                weighted_ideal_lattice_gram_lll_reduce(
                        basis, transform, scaled_gram, original_basis, ideal,
                        *order, weights, embedding_precision,
                        active_precision, diagnostics, cache);
        if (reduction_status == IdealLatticeReductionStatus::retry_precision) {
            if (active_precision > kIdealLatticeLllMaxPrecision / 2) {
                return false;
            }
            active_precision *= 2;
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
