#pragma once

#include <silex/flint/fmpz_mat.hpp>

namespace silex::lat::detail {

enum class FplllBackendStatus {
    success,
    unavailable,
    invalid_dimensions,
    dimension_overflow,
    conversion_failed,
    reduction_failed
};

struct FplllLllResult {
    FplllBackendStatus status = FplllBackendStatus::invalid_dimensions;
    int backend_status = 0;
};

FplllLllResult fplll_row_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double delta = 0.99) noexcept;

FplllLllResult fplll_row_bkz_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        int block_size,
        int max_loops = 1) noexcept;

FplllLllResult fplll_column_image_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double delta = 0.99) noexcept;

}  // namespace silex::lat::detail
