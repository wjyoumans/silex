#pragma once

#include <flint/fmpz_mat.h>

namespace silex::detail {

// Returns true when `basis` is square, upper triangular and has a nonzero
// diagonal, the shape `Ideal` stores for a full-rank HNF.
bool is_full_rank_upper_triangular(const fmpz_mat_t basis) noexcept;

// Exact membership of `vector` (length n) in the row lattice of the n x n
// basis B, which must pass `is_full_rank_upper_triangular`.  `vector` may be
// used as the residual and overwritten.  When |det B| fits in a word the test
// runs on residues modulo |det B| and leaves `vector` unchanged.
bool triangular_basis_contains(const fmpz_mat_t basis, fmpz* vector) noexcept;

// The same test on a read-only `vector`, copied into a residual only when
// |det B| does not fit in a word.
bool triangular_basis_contains_row(const fmpz_mat_t basis,
                                   const fmpz* vector) noexcept;

}  // namespace silex::detail
