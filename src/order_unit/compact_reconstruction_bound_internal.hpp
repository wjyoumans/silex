#pragma once

#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>
#include <silex/order_unit.hpp>

namespace silex::detail {

enum class CompactCoordinateBoundStatus {
    success,
    invalid_input,
    unsupported_order_basis,
    precision_exhausted,
    bound_overflow,
};

struct CompactCoordinateBoundOptions {
    slong start_precision = 128;
    slong max_precision = 4096;
};

struct CompactCoordinateBoundReport {
    CompactCoordinateBoundStatus status =
            CompactCoordinateBoundStatus::invalid_input;
    slong precision = 0;
    slong attempts = 0;
    ulong bit_bound = 0;
    flint::Fmpz coordinate_bound;
    flint::Arb log2_minkowski_bound;
    flint::Arb log2_unit_bound;
    flint::Arb log2_coordinate_bound;
};

// The caller certifies that generators are the complete fundamental units for
// order and evaluate to elements of order. This private producer only ports
// reference's coordinate-bound calculation; it does not verify those two facts by
// expanding the compact elements.
bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options = {}) noexcept;

}  // namespace silex::detail
