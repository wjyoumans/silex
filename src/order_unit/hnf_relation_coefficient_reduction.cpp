#include "order_unit_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>

#include <silex/flint/arf.hpp>
#include <silex/flint/fmpz_lll.hpp>

namespace silex::detail {

namespace {

slong trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 && ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), i) != 0) {
            return -1;
        }
    }
    return rows;
}

}  // namespace

bool arf_2exp_valuation_slong(slong& out, const arf_struct* value) noexcept {
    if (arf_is_zero(value) != 0) {
        return false;
    }

    flint::Fmpz mantissa;
    flint::Fmpz exponent;
    arf_get_fmpz_2exp(mantissa.raw(), exponent.raw(), value);
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(exponent)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(mantissa))) {
        return false;
    }

    const slong exponent_si =
            flint::fmpz_get_si(flint::FmpzConstRef(exponent));
    const ulong mantissa_v2 = fmpz_val2(mantissa.raw());
    if (mantissa_v2 >
        static_cast<ulong>(std::numeric_limits<slong>::max())) {
        return false;
    }
    const slong mantissa_v2_si = static_cast<slong>(mantissa_v2);
    if (exponent_si > std::numeric_limits<slong>::max() - mantissa_v2_si) {
        return false;
    }

    out = exponent_si + mantissa_v2_si;
    return true;
}

slong relation_kernel_native_word_bits() noexcept {
#ifdef _WIN64
    using SelectedSignedWord = long long;
    using SelectedUnsignedWord = unsigned long long;
#else
    using SelectedSignedWord = long;
    using SelectedUnsignedWord = unsigned long;
#endif
    static_assert(std::is_same_v<slong, SelectedSignedWord>);
    static_assert(std::is_same_v<ulong, SelectedUnsignedWord>);
    constexpr auto word_bits =
            std::numeric_limits<SelectedUnsignedWord>::digits;
    static_assert(word_bits == 32 || word_bits == 64);
    return static_cast<slong>(word_bits);
}

slong round_positive_word_bits(slong bits, slong word_bits) noexcept {
    if ((word_bits != 32 && word_bits != 64) || bits <= 0 ||
        bits > std::numeric_limits<slong>::max() - (word_bits - 1)) {
        return -1;
    }
    return ((bits + (word_bits - 1)) / word_bits) * word_bits;
}

bool rescale_arb_matrix_to_int(flint::FmpzMat& out,
                                    const flint::ArbMat& matrix,
                                    slong precision,
                                    slong word_bits) noexcept {
    // reference lll(real matrix) calls polarit2.c:RgM_rescale_to_int before
    // integer LLL.  For t_REAL entries, reference uses represented real precision:
    // e = expo(c) + 1 - bit_prec(c), plus trailing zero bits in the
    // represented limbs, then grndtoi(x * 2^(-emin)).  Arb is still not reference
    // t_REAL, but rounding midpoints to the local word-rounded represented
    // precision keeps this relation-kernel boundary source-shaped.
    const slong rows = flint::arb_mat_nrows_value(matrix);
    const slong cols = flint::arb_mat_ncols_value(matrix);
    if ((word_bits != 32 && word_bits != 64) || precision <= 0 ||
        flint::fmpz_mat_nrows(out) != rows ||
        flint::fmpz_mat_ncols(out) != cols) {
        return false;
    }

    const slong represented_precision =
            round_positive_word_bits(precision, word_bits);
    if (represented_precision <= 0) {
        return false;
    }

    std::vector<flint::Arf> represented(
            static_cast<std::size_t>(rows * cols));
    slong emin = std::numeric_limits<slong>::max();
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            const arb_struct* entry = arb_mat_entry(matrix.raw(), i, j);
            if (arb_is_finite(entry) == 0) {
                return false;
            }
            flint::Arf& rounded =
                    represented[static_cast<std::size_t>(i * cols + j)];
            ::arf_set_round(rounded.raw(), arb_midref(entry),
                            represented_precision, ARF_RND_NEAR);
            if (!flint::arf_is_finite(rounded)) {
                return false;
            }
            if (arf_is_zero(rounded.raw()) != 0) {
                continue;
            }
            slong entry_emin = 0;
            if (!arf_2exp_valuation_slong(entry_emin, rounded.raw())) {
                return false;
            }
            emin = std::min(emin, entry_emin);
        }
    }

    if (emin == std::numeric_limits<slong>::max()) {
        flint::fmpz_mat_zero(flint::FmpzMatRef(out));
        return true;
    }

    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            arf_get_fmpz_fixed_si(
                    flint::fmpz_mat_entry(out, i, j).raw(),
                    represented[static_cast<std::size_t>(i * cols + j)].raw(),
                    emin);
        }
    }

    return true;
}

bool relation_kernel_rescale_log_matrix_for_testing(
        flint::FmpzMatRef rows,
        flint::ArbMatConstRef matrix,
        slong precision,
        slong word_bits) noexcept {
    const slong row_count = flint::arb_mat_nrows_value(matrix);
    const slong column_count = flint::arb_mat_ncols_value(matrix);
    if (round_positive_word_bits(precision, word_bits) <= 0 ||
        flint::fmpz_mat_nrows(rows) != row_count ||
        flint::fmpz_mat_ncols(rows) != column_count) {
        return false;
    }

    flint::ArbMat owned_matrix(row_count, column_count);
    flint::FmpzMat owned_rows(row_count, column_count);
    flint::arb_mat_set(flint::ArbMatRef(owned_matrix), matrix);
    if (!rescale_arb_matrix_to_int(owned_rows, owned_matrix, precision,
                                        word_bits)) {
        return false;
    }
    flint::fmpz_mat_set(rows, flint::FmpzMatConstRef(owned_rows));
    return true;
}

bool transform_relation_coefficients_by_rows(
        flint::FmpzMat& coefficients,
        flint::FmpzMatConstRef transform) noexcept;

bool reduce_relation_coefficients_by_log_lll(
        flint::FmpzMat& coefficients,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong& precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong rank = flint::fmpz_mat_nrows(coefficients);
    if (rank <= 0 ||
        flint::fmpz_mat_ncols(coefficients) != class_group.relation_count()) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "exact coefficient log LLL failed: invalid shape");
        return false;
    }

    slong places = 0;
    if (!compact_places(places, embeddings)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "exact coefficient log LLL failed: compact places");
        return false;
    }

    slong work_precision = precision;
    for (slong attempt = 0; attempt < 4; ++attempt) {
        if (!embeddings.refine(work_precision, diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "exact coefficient log LLL failed: embedding "
                      "refine");
            return false;
        }

        flint::ArbMat logs(rank, places);
        flint::FmpzMat scaled(rank, places);
        flint::FmpzMat reduced(rank, places);
        flint::FmpzMat transform(rank, rank);
        const bool have_logs =
                relation_coefficients_log_matrix(
                        logs, class_group, embeddings, coefficients,
                        work_precision, diagnostics);
        const bool rounded =
                have_logs &&
                rescale_arb_matrix_to_int(
                        scaled, logs, work_precision,
                        relation_kernel_native_word_bits());
        if (!have_logs) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "exact coefficient log LLL failed: log matrix");
        } else if (!rounded) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "exact coefficient log LLL failed: integer "
                      "rescale");
        }

        if (rounded) {
            fmpz_mat_set(reduced.raw(), scaled.raw());
            fmpz_mat_one(transform.raw());
            flint::FmpzLll lll;
            fmpz_lll(reduced.raw(), transform.raw(), lll.raw());
            const slong reduced_rank =
                    trimmed_nonzero_rows(flint::FmpzMatConstRef(reduced));
            if (reduced_rank == rank) {
                if (fmpz_mat_is_one(transform.raw()) == 0 &&
                    !transform_relation_coefficients_by_rows(
                            coefficients,
                            flint::FmpzMatConstRef(transform))) {
                    SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "exact coefficient log LLL failed: "
                              "transform coefficients");
                    return false;
                }
                precision = work_precision;
                return true;
            }
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "exact coefficient log LLL failed: reduced rank "
                      "defect");
        }

        if (work_precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        work_precision *= 2;
    }
    return false;
}

bool transform_relation_coefficients_by_rows(
        flint::FmpzMat& coefficients,
        flint::FmpzMatConstRef transform) noexcept {
    const slong rank = flint::fmpz_mat_nrows(coefficients);
    const slong relation_count = flint::fmpz_mat_ncols(coefficients);
    if (rank <= 0 || fmpz_mat_nrows(transform.raw()) != rank ||
        fmpz_mat_ncols(transform.raw()) != rank) {
        return false;
    }

    flint::FmpzMat transformed(rank, relation_count);
    fmpz_mat_mul(transformed.raw(), transform.raw(), coefficients.raw());
    coefficients = std::move(transformed);
    return true;
}

}  // namespace silex::detail
