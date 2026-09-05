#pragma once

#include <silex/flint/fmpz_mat.hpp>

namespace silex::lat::detail {

enum class LllRoute { certified, fallback, ineligible };
enum class LllTestFailure { none, reducer, certification };

// Input must be normalized, full row rank, and distinct from reduced.
// Both matrices have the same dimensions. Failure injection is test-only.
LllRoute reduce_normalized_basis(flint::FmpzMatRef reduced,
        flint::FmpzMatConstRef input,
        LllTestFailure failure = LllTestFailure::none) noexcept;

}  // namespace silex::lat::detail
