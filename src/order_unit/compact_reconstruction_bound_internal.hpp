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
    // The reference expression: log2Mbound (Hadamard estimate) plus
    // log2fubound, plus the basis normalization scale.
    flint::Arb log2_minkowski_bound;
    flint::Arb log2_unit_bound;
    flint::Arb log2_coordinate_bound;
    // The quantity log2Mbound bounds, log2 of the max-row-sum norm of the
    // inverse split real/imaginary embedding matrix, and the same total with
    // it in place of log2_minkowski_bound.  Indeterminate when the inverse
    // was not certified at `precision`.
    flint::Arb log2_inverse_norm_bound;
    flint::Arb log2_inverse_norm_coordinate_bound;
    // Ceiling of log2_coordinate_bound.  bit_bound is the ceiling of
    // log2_inverse_norm_coordinate_bound when the inverse was certified and
    // that ceiling is stable (inverse_norm_bound_used), and this ceiling
    // otherwise.
    ulong hadamard_bit_bound = 0;
    bool inverse_norm_bound_used = false;
};

// The caller certifies that generators are units of order and evaluate to
// elements of order; order must be a maximal order with a known basis
// (checked: otherwise the bound is refused as invalid_input).  Exactly
// r1 + r2 - 1 generators are required (a full-rank subgroup; fewer are
// refused).  They need not be the complete fundamental units: the bound uses
// only each generator's embeddings, and membership in order makes the bounded
// coordinates integral, which reconstruction needs.  So it also applies to a
// finite-index subgroup in the middle of extraction.
// This private producer only ports reference's coordinate-bound calculation;
// it does not verify those facts by expanding the compact elements.
bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options = {}) noexcept;

struct WideFactoredUnits;

// The same bound for a batch with arbitrary-size exponents.  The bound uses
// only the units' log embeddings (reference log2fubound) and the order basis
// (log2Mbound), so the representation of the units does not matter; the logs
// are formed as exponent matrix times base logs.
bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, const WideFactoredUnits& generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options = {}) noexcept;

}  // namespace silex::detail
