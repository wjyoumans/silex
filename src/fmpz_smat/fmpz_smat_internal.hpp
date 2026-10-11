#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex::detail {

enum class SquareAbsDetRoute { diagonal_product, determinant };

// Sets out to |det(matrix)| for a square matrix. An upper triangular matrix
// takes the product of its diagonal; any other square matrix falls back to
// fmpz_mat_det. The returned route is for tests.
SquareAbsDetRoute fmpz_mat_square_abs_det(flint::FmpzRef out,
        flint::FmpzMatConstRef matrix) noexcept;

}  // namespace silex::detail
