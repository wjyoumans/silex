#pragma once

#include <flint/fmpz_mat.h>

namespace silex::detail {

// Returns true when `basis` is square, upper triangular and has a nonzero
// diagonal, the shape `Ideal` stores for a full-rank HNF.
bool is_full_rank_upper_triangular(const fmpz_mat_t basis) noexcept;

// Exact membership of `vector` (length n) in the row lattice of the n x n
// basis B, which must pass `is_full_rank_upper_triangular`.  `vector` is used
// as the residual and is overwritten.
bool triangular_basis_contains(const fmpz_mat_t basis, fmpz* vector) noexcept;

}  // namespace silex::detail
