#include "order_unit_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>

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

bool column_lattice_hnf(flint::FmpzMat& out,
                        const flint::FmpzMat& matrix) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    if (rows <= 0 || cols < 0) {
        return false;
    }
    if (cols == 0) {
        flint::FmpzMat empty(0, rows);
        out = std::move(empty);
        return true;
    }

    flint::FmpzMat transpose(cols, rows);
    flint::FmpzMat hnf(cols, rows);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transpose),
                              flint::FmpzMatConstRef(matrix));
    ::fmpz_mat_hnf(hnf.raw(), transpose.raw());

    const slong rank = trimmed_nonzero_rows(flint::FmpzMatConstRef(hnf));
    if (rank < 0) {
        return false;
    }

    flint::FmpzMat trimmed(rank, rows);
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < rows; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(trimmed, i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(hnf), i, j));
        }
    }
    out = std::move(trimmed);
    return true;
}

bool column_lattice_hnf_equal(const flint::FmpzMat& left,
                              const flint::FmpzMat& right) noexcept {
    return flint::fmpz_mat_nrows(left) == flint::fmpz_mat_nrows(right) &&
           flint::fmpz_mat_ncols(left) == flint::fmpz_mat_ncols(right) &&
           flint::fmpz_mat_equal(flint::FmpzMatConstRef(left),
                                 flint::FmpzMatConstRef(right));
}

bool select_columns(flint::FmpzMat& out,
                    const flint::FmpzMat& matrix,
                    const std::vector<slong>& columns) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    if (rows < 0 || cols < 0 ||
        columns.size() >
                static_cast<std::size_t>(std::numeric_limits<slong>::max())) {
        return false;
    }

    flint::FmpzMat selected(rows, static_cast<slong>(columns.size()));
    for (slong j = 0; j < static_cast<slong>(columns.size()); ++j) {
        const slong source_col = columns[static_cast<std::size_t>(j)];
        if (source_col < 0 || source_col >= cols) {
            return false;
        }
        for (slong i = 0; i < rows; ++i) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(selected, i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(matrix), i, source_col));
        }
    }
    out = std::move(selected);
    return true;
}

bool extract_full_lattice_column_indices(
        std::vector<slong>& out,
        bool& extracted,
        const flint::FmpzMat& matrix) noexcept {
    out.clear();
    extracted = false;
    const slong columns = flint::fmpz_mat_ncols(matrix);
    if (columns < 0) {
        return false;
    }
    if (columns + 1 < 200) {
        return true;
    }

    flint::FmpzMat full_hnf(0, 0);
    if (!column_lattice_hnf(full_hnf, matrix)) {
        return false;
    }
    flint::FmpzMat current_hnf(0, flint::fmpz_mat_nrows(matrix));

    // The source full-lattice extraction uses 1-based columns and lg(x) as
    // column_count + 1.  Keep the same variables here to mirror its block
    // splitting and elimination decisions.
    const slong l = columns + 1;
    slong dj = 1;
    for (slong j = 1; j < l;) {
        const slong previous_size = static_cast<slong>(out.size());
        const slong block = std::min(dj, l - j);
        for (slong k = 0; k < block; ++k) {
            out.push_back(j + k - 1);
        }

        flint::FmpzMat selected(0, 0);
        flint::FmpzMat hnf(0, 0);
        if (!select_columns(selected, matrix, out) ||
            !column_lattice_hnf(hnf, selected)) {
            return false;
        }

        if (column_lattice_hnf_equal(hnf, current_hnf)) {
            out.resize(static_cast<std::size_t>(previous_size));
            j += block;
            if (j >= l) {
                break;
            }
            dj = block;
            if (dj > std::numeric_limits<slong>::max() / 2) {
                return false;
            }
            dj <<= 1;
            if (j + dj >= l) {
                dj = (l - j) >> 1;
                if (dj == 0) {
                    dj = 1;
                }
            }
        } else if (block > 1) {
            out.resize(static_cast<std::size_t>(previous_size));
            dj = block >> 1;
            if (dj == 0) {
                dj = 1;
            }
        } else {
            if (column_lattice_hnf_equal(hnf, full_hnf)) {
                break;
            }
            current_hnf = std::move(hnf);
            ++j;
        }
    }

    extracted = true;
    return true;
}

}  // namespace

bool transformed_hnf_unit_coefficients_from_regulator_matrix(
        flint::FmpzMat& out,
        const ClassGroupContext& class_group,
        const flint::FmpzMat& witness_coefficients,
        const flint::FmpzMat& integer_coordinates,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong rank = flint::fmpz_mat_nrows(integer_coordinates);
    const slong columns = flint::fmpz_mat_ncols(integer_coordinates);
    const slong witness_count = flint::fmpz_mat_nrows(witness_coefficients);
    const Order* order = class_group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (rank <= 0 || columns <= 0 ||
        columns != witness_count ||
        flint::fmpz_mat_ncols(witness_coefficients) !=
                class_group.relation_count() ||
        field == nullptr) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "transform witnesses failed: invalid dimensions");
        return false;
    }
    bool extracted = false;
    std::vector<slong> selected_indices;
    if (!extract_full_lattice_column_indices(
                selected_indices, extracted, integer_coordinates)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "transform witnesses failed: extract_full_lattice");
        return false;
    }

    flint::FmpzMat selected_matrix(0, 0);
    const flint::FmpzMat* unit_matrix = &integer_coordinates;
    const std::vector<slong>* witness_indices = nullptr;
    if (extracted) {
        if (selected_indices.empty()) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "transform witnesses failed: empty extracted lattice");
            return false;
        }
        if (!select_columns(selected_matrix, integer_coordinates,
                            selected_indices)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "transform witnesses failed: select columns");
            return false;
        }
        unit_matrix = &selected_matrix;
        witness_indices = &selected_indices;
    }

    const slong active_columns = flint::fmpz_mat_ncols(*unit_matrix);
    if (active_columns <= 0) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "transform witnesses failed: no active columns");
        return false;
    }

    flint::FmpzMat transposed(active_columns, rank);
    flint::FmpzMat reduced(active_columns, rank);
    flint::FmpzMat transform(active_columns, active_columns);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transposed),
                              flint::FmpzMatConstRef(*unit_matrix));
    fmpz_mat_set(reduced.raw(), transposed.raw());
    fmpz_mat_one(transform.raw());
    flint::FmpzLll lll;
    fmpz_lll(reduced.raw(), transform.raw(), lll.raw());

    flint::FmpzMat transformed(rank, class_group.relation_count());
    slong output_row = 0;
    for (slong row = 0; row < active_columns &&
                        output_row < rank;
         ++row) {
        if (::fmpz_mat_is_zero_row(reduced.raw(), row) != 0) {
            continue;
        }

        for (slong col = 0; col < active_columns; ++col) {
            const slong witness_index =
                    witness_indices == nullptr
                            ? col
                            : (*witness_indices)[static_cast<std::size_t>(
                                      col)];
            const flint::FmpzConstRef exponent =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(transform), row, col);
            if (flint::fmpz_is_zero(exponent)) {
                continue;
            }
            for (slong rel = 0; rel < class_group.relation_count(); ++rel) {
                ::fmpz_addmul(
                        ::fmpz_mat_entry(transformed.raw(), output_row, rel),
                        exponent.raw(),
                        ::fmpz_mat_entry(witness_coefficients.raw(),
                                         witness_index, rel));
            }
        }
        ++output_row;
    }

    if (output_row != rank) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "transform witnesses failed: transformed rank defect");
        return false;
    }
    out = std::move(transformed);
    return true;
}

}  // namespace silex::detail
