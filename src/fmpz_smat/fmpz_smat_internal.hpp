#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex::detail {

enum class SquareAbsDetRoute { diagonal_product, determinant };

// Sets out to |det(matrix)| for a square matrix. An upper triangular matrix
// takes the product of its diagonal; any other square matrix falls back to
// fmpz_mat_det. A non-square matrix is also passed to fmpz_mat_det, which
// rejects it. The returned route is for tests.
SquareAbsDetRoute fmpz_mat_square_abs_det(flint::FmpzRef out,
        flint::FmpzMatConstRef matrix) noexcept;

// Folds rows [hnf_rows, nrows) of `matrix` into the canonical row HNF held in
// its first hnf_rows rows, in place.  The first hnf_rows rows must be a
// canonical row HNF with no zero rows (as fmpz_mat_hnf returns them).  On
// success the first `rank` rows are the canonical row HNF of all rows (the
// fmpz_mat_hnf normalization) and the remaining rows are zero.  Returns
// false, with `matrix` unchanged, when hnf_rows is out of range or the first
// hnf_rows rows do not have positive pivots in strictly increasing columns.
bool fmpz_mat_hnf_add_rows(slong& rank, flint::FmpzMatRef matrix,
                           slong hnf_rows) noexcept;

}  // namespace silex::detail
