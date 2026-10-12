#pragma once

#include <silex/flint/fmpz_mat.hpp>

#include "certified_lll_internal.hpp"

namespace silex::lat::detail {

// Input must be normalized, full row rank, and distinct from reduced.
// Both matrices have the same dimensions.  Sets `reduced` to the basis
// fmpz_lll returns for `input` (default context, zero transform).  For at
// most 14 rows and 14 columns it takes certified_fmpz_lll; larger inputs
// call fmpz_lll directly and report ineligible.
CertifiedLllRoute reduce_normalized_basis(flint::FmpzMatRef reduced,
        flint::FmpzMatConstRef input) noexcept;

}  // namespace silex::lat::detail
