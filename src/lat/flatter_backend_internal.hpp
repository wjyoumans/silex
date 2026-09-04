#pragma once

#include <silex/flint/fmpz_mat.hpp>

namespace silex::lat::detail {

enum class FlatterBackendStatus {
    success,
    unavailable,
    invalid_dimensions,
    dimension_overflow,
    conversion_failed,
    transform_unavailable,
    reduction_failed
};

struct FlatterLllResult {
    FlatterBackendStatus status = FlatterBackendStatus::invalid_dimensions;
    slong rank = 0;
};

FlatterLllResult flatter_column_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double root_hermite_factor = 1.02,
        unsigned int max_threads = 1) noexcept;

}  // namespace silex::lat::detail
