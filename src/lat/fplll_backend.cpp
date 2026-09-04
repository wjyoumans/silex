#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
#include <gmp.h>
#endif

#include "fplll_backend_internal.hpp"

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/flint.h>

#include <limits>

#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
#include <fplll/bkz.h>
#include <fplll/wrapper.h>
#endif

namespace silex::lat::detail {

namespace {

bool valid_output_dimensions(flint::FmpzMatRef reduced,
                             flint::FmpzMatRef transform,
                             flint::FmpzMatConstRef input) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    return rows >= 0 && cols >= 0 &&
           fmpz_mat_nrows(reduced.raw()) == rows &&
           fmpz_mat_ncols(reduced.raw()) == cols &&
           fmpz_mat_nrows(transform.raw()) == rows &&
           fmpz_mat_ncols(transform.raw()) == rows;
}

bool valid_column_image_output_dimensions(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        slong rank) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    return rows >= 0 && cols >= 0 && rank >= 0 &&
           fmpz_mat_nrows(reduced.raw()) == rows &&
           fmpz_mat_ncols(reduced.raw()) == rank &&
           fmpz_mat_nrows(transform.raw()) == cols &&
           fmpz_mat_ncols(transform.raw()) == rank;
}

bool dimensions_fit_int(flint::FmpzMatConstRef input) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    return rows <= std::numeric_limits<int>::max() &&
           cols <= std::numeric_limits<int>::max();
}

#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
void set_fplll_entry_from_fmpz(fplll::Z_NR<mpz_t>& out,
                               const fmpz* value) noexcept {
    fmpz_get_mpz(out.get_data(), value);
}

void set_fmpz_from_fplll_entry(fmpz* out,
                               const fplll::Z_NR<mpz_t>& value) noexcept {
    fmpz_set_mpz(out, value.get_data());
}

void copy_to_fplll(fplll::ZZ_mat<mpz_t>& out,
                   flint::FmpzMatConstRef input) noexcept {
    const int rows = static_cast<int>(fmpz_mat_nrows(input.raw()));
    const int cols = static_cast<int>(fmpz_mat_ncols(input.raw()));
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            set_fplll_entry_from_fmpz(
                    out(row, col),
                    fmpz_mat_entry(input.raw(), row, col));
        }
    }
}

void copy_from_fplll(flint::FmpzMatRef out,
                     const fplll::ZZ_mat<mpz_t>& input) noexcept {
    const int rows = static_cast<int>(fmpz_mat_nrows(out.raw()));
    const int cols = static_cast<int>(fmpz_mat_ncols(out.raw()));
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            set_fmpz_from_fplll_entry(
                    fmpz_mat_entry(out.raw(), row, col), input(row, col));
        }
    }
}
#endif

}  // namespace

FplllLllResult fplll_row_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double delta) noexcept {
    if (!valid_output_dimensions(reduced, transform, input) ||
        !(delta > 0.25 && delta < 1.0)) {
        return {FplllBackendStatus::invalid_dimensions, 0};
    }
    if (!dimensions_fit_int(input)) {
        return {FplllBackendStatus::dimension_overflow, 0};
    }

#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
    const int rows = static_cast<int>(fmpz_mat_nrows(input.raw()));
    const int cols = static_cast<int>(fmpz_mat_ncols(input.raw()));
    fplll::ZZ_mat<mpz_t> basis(rows, cols);
    fplll::ZZ_mat<mpz_t> row_transform;
    row_transform.gen_identity(rows);
    copy_to_fplll(basis, input);

    const int backend_status =
            fplll::lll_reduction(basis, row_transform, delta, 0.51,
                                 fplll::LM_WRAPPER, fplll::FT_DEFAULT, 0,
                                 fplll::LLL_DEFAULT);
    if (backend_status != 0) {
        return {FplllBackendStatus::reduction_failed, backend_status};
    }

    copy_from_fplll(reduced, basis);
    copy_from_fplll(transform, row_transform);
    return {FplllBackendStatus::success, 0};
#else
    (void)reduced;
    (void)transform;
    (void)input;
    return {FplllBackendStatus::unavailable, 0};
#endif
}

FplllLllResult fplll_row_bkz_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        int block_size,
        int max_loops) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    if (!valid_output_dimensions(reduced, transform, input) ||
        block_size < 2 || block_size > rows || max_loops <= 0) {
        return {FplllBackendStatus::invalid_dimensions, 0};
    }
    if (!dimensions_fit_int(input)) {
        return {FplllBackendStatus::dimension_overflow, 0};
    }

#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
    const int row_count = static_cast<int>(rows);
    const int cols = static_cast<int>(fmpz_mat_ncols(input.raw()));
    fplll::ZZ_mat<mpz_t> basis(row_count, cols);
    fplll::ZZ_mat<mpz_t> row_transform;
    row_transform.gen_identity(row_count);
    copy_to_fplll(basis, input);

    fplll::vector<fplll::Strategy> strategies;
    fplll::BKZParam params(
            block_size, strategies, fplll::LLL_DEF_DELTA,
            fplll::BKZ_DEFAULT | fplll::BKZ_MAX_LOOPS, max_loops);
    const int backend_status = fplll::bkz_reduction(
            &basis, &row_transform, params, fplll::FT_DEFAULT, 0);
    // A loop limit is the expected successful stop for a bounded tour.
    if (backend_status != fplll::RED_SUCCESS &&
        backend_status != fplll::RED_BKZ_LOOPS_LIMIT) {
        return {FplllBackendStatus::reduction_failed, backend_status};
    }

    copy_from_fplll(reduced, basis);
    copy_from_fplll(transform, row_transform);
    return {FplllBackendStatus::success, backend_status};
#else
    (void)reduced;
    (void)transform;
    (void)input;
    (void)block_size;
    (void)max_loops;
    return {FplllBackendStatus::unavailable, 0};
#endif
}

FplllLllResult fplll_column_image_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double delta) noexcept {
    const slong rank = fmpz_mat_rank(input.raw());
    if (!valid_column_image_output_dimensions(
                reduced, transform, input, rank) ||
        !(delta > 0.25 && delta < 1.0)) {
        return {FplllBackendStatus::invalid_dimensions, 0};
    }

    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    if (rank == 0) {
        fmpz_mat_zero(reduced.raw());
        fmpz_mat_zero(transform.raw());
        return {FplllBackendStatus::success, 0};
    }

    flint::FmpzMat transpose(cols, rows);
    flint::FmpzMat reduced_rows(cols, rows);
    flint::FmpzMat row_transform(cols, cols);
    fmpz_mat_transpose(transpose.raw(), input.raw());

    const FplllLllResult row_result = fplll_row_lll_transform(
            flint::FmpzMatRef(reduced_rows),
            flint::FmpzMatRef(row_transform),
            flint::FmpzMatConstRef(transpose), delta);
    if (row_result.status != FplllBackendStatus::success) {
        return row_result;
    }

    fmpz_mat_zero(reduced.raw());
    fmpz_mat_zero(transform.raw());
    slong selected = 0;
    for (slong source_row = 0; source_row < cols && selected < rank;
         ++source_row) {
        if (fmpz_mat_is_zero_row(reduced_rows.raw(), source_row) != 0) {
            continue;
        }
        for (slong row = 0; row < rows; ++row) {
            fmpz_set(fmpz_mat_entry(reduced.raw(), row, selected),
                     fmpz_mat_entry(reduced_rows.raw(), source_row, row));
        }
        for (slong row = 0; row < cols; ++row) {
            fmpz_set(fmpz_mat_entry(transform.raw(), row, selected),
                     fmpz_mat_entry(row_transform.raw(), source_row, row));
        }
        ++selected;
    }
    if (selected != rank) {
        return {FplllBackendStatus::reduction_failed, 0};
    }

    flint::FmpzMat expected(rows, rank);
    fmpz_mat_mul(expected.raw(), input.raw(), transform.raw());
    if (fmpz_mat_equal(expected.raw(), reduced.raw()) == 0 ||
        fmpz_mat_rank(reduced.raw()) != rank) {
        return {FplllBackendStatus::reduction_failed, 0};
    }
    return {FplllBackendStatus::success, 0};
}

}  // namespace silex::lat::detail
