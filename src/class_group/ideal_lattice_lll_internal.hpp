#pragma once

#include <silex/class_group.hpp>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

}  // namespace silex::detail

namespace silex::detail::relation_search {

struct IdealLatticeLllData {
    flint::FmpzMat basis{0, 0};
    flint::FmpzMat scaled_gram{0, 0};
    flint::Fmpz gram_denominator;
};

// Per-order data for the reference-style per-ideal reduction of
// build_ideal_lattice_rounded_embedding_lll_data: the integer-rounded
// Minkowski embedding rows of the order basis (reference nf_get_roundG) and
// the scaled order T2 Gram numerator at `order_gram_precision`, the precision
// the build_ideal_lattice_lll_data schedule reaches from the starting
// precision `order_gram_requested_precision` (0 while no Gram is cached).
// The cache is keyed on the order's data and is filled on first use.
// `route_checked` records that the first one-LLL basis accepted for the
// order was compared with the two-LLL basis of the same ideal, and
// `use_two_lll_route` that the comparison failed, so every further ideal of
// the order is reduced by build_ideal_lattice_lll_data.
struct IdealLatticeRoundedEmbeddingCache {
    Order order;
    bool has_order = false;
    bool has_rounded_rows = false;
    flint::FmpzMat rounded_rows{0, 0};
    slong order_gram_requested_precision = 0;
    slong order_gram_precision = 0;
    flint::FmpzMat order_gram{0, 0};
    bool route_checked = false;
    bool use_two_lll_route = false;
};

bool copy_ideal_basis_row_coordinates(flint::FmpzMat& out,
                                      flint::FmpzMatConstRef basis,
                                      slong row) noexcept;

bool build_ideal_lattice_lll_data(IdealLatticeLllData& out,
                                  const Ideal& ideal,
                                  slong precision) noexcept;

enum class IdealLatticeOneLllRoute { certified, fallback, ineligible };
enum class IdealLatticeOneLllTestFailure { none, reducer, certificate };

// The integer LLL (delta 0.99, eta 0.51) of
// build_ideal_lattice_rounded_embedding_lll_data: replaces `embedded` by the
// matrix fmpz_lll returns for it and sets `transform` (resized to a square
// matrix of the row count if needed) to fmpz_lll's transform from the
// identity.  certified: FLINT's first phase, certified by
// fmpz_mat_is_reduced; fallback: the phase failed and fmpz_lll ran;
// ineligible: fewer than two rows or entries above 250 bits, fmpz_lll ran.
// Failure injection is test-only: it zeroes the phase output, and
// `reducer` also reads the phase status as -1.
IdealLatticeOneLllRoute ideal_lattice_rounded_embedding_integer_lll(
        flint::FmpzMat& embedded,
        flint::FmpzMat& transform,
        IdealLatticeOneLllTestFailure failure =
                IdealLatticeOneLllTestFailure::none) noexcept;

// One integer LLL per ideal on the ideal basis times the cached rounded
// order embedding (reference buch2.c:Fincke_Pohst_ideal), with the same output contract
// as build_ideal_lattice_lll_data.  Falls back to build_ideal_lattice_lll_data
// when the result fails that route's T2 reduction-quality checks, and for
// every ideal of an order whose first accepted one-LLL basis has a T2
// diagonal product above 2^degree times that of the two-LLL basis.
bool build_ideal_lattice_rounded_embedding_lll_data(
        IdealLatticeLllData& out,
        const Ideal& ideal,
        slong precision,
        IdealLatticeRoundedEmbeddingCache& cache,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

// The route check of build_ideal_lattice_rounded_embedding_lll_data, exposed
// for tests: true when the T2 diagonal product of `one_lll_gram` (denominator
// 2^one_denominator_bits) is at most 2^degree times that of `two_lll`.  False
// on malformed inputs, which sends the order to the two-LLL route.
bool ideal_lattice_rounded_embedding_basis_matches_two_lll(
        flint::FmpzMatConstRef one_lll_gram,
        slong one_denominator_bits,
        const IdealLatticeLllData& two_lll) noexcept;

bool ideal_lattice_weight_vector_shape_is_valid(flint::FmpzMatConstRef weights,
                                                slong degree) noexcept;

bool ideal_lattice_weight_vector_is_zero(
        flint::FmpzMatConstRef weights) noexcept;

bool build_weighted_ideal_lattice_reduction(
        flint::FmpzMat& basis, flint::FmpzMat& transform,
        flint::FmpzMat& scaled_gram, flint::Fmpz& gram_denominator,
        const Ideal& ideal, flint::FmpzMatConstRef weights, slong precision,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* cache) noexcept;

}  // namespace silex::detail::relation_search
