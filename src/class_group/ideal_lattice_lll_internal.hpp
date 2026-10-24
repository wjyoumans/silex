#pragma once

#include <silex/class_group.hpp>

#include "../lat/certified_lll_internal.hpp"

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

// Sets `out` to basis * gram * basis^T for a symmetric `gram`, computing
// basis * gram with one matrix product and only the upper triangle of the
// symmetric result by dot products, then mirroring it.  Exact, so equal
// entry for entry to basis * (gram * basis^T) when gram is symmetric.
// Returns false on a shape mismatch (out must be rows(basis) square, gram
// cols(basis) square).  `out` must not alias `basis` (it may alias `gram`).
bool ideal_lattice_symmetric_congruence_product(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef basis,
        flint::FmpzMatConstRef gram) noexcept;

bool copy_ideal_basis_row_coordinates(flint::FmpzMat& out,
                                      flint::FmpzMatConstRef basis,
                                      slong row) noexcept;

bool build_ideal_lattice_lll_data(IdealLatticeLllData& out,
                                  const Ideal& ideal,
                                  slong precision) noexcept;

// The two LLLs of the two-LLL route (ideal_lattice_gram_lll_reduce and its
// weighted form), each through lat::detail::certified_fmpz_lll, whose
// route they return.  Each replaces its first argument by the matrix
// fmpz_lll returns for it and sets `transform` (resized to a square matrix
// of the row count if needed) to fmpz_lll's transform from the identity.
// ideal_lattice_pre_reduction_lll: the Z_BASIS pre-reduction of the ideal
// basis, delta 0.5, eta 0.51.  ideal_lattice_gram_lll: the GRAM reduction
// of the scaled T2 Gram matrix, delta 0.99, eta 0.51.
lat::detail::CertifiedLllRoute ideal_lattice_pre_reduction_lll(
        flint::FmpzMat& basis,
        flint::FmpzMat& transform) noexcept;
lat::detail::CertifiedLllRoute ideal_lattice_gram_lll(
        flint::FmpzMat& gram,
        flint::FmpzMat& transform) noexcept;

// The integer LLL (delta 0.99, eta 0.51) of
// build_ideal_lattice_rounded_embedding_lll_data: replaces `embedded` by the
// matrix fmpz_lll returns for it and sets `transform` (resized to a square
// matrix of the row count if needed) to fmpz_lll's transform from the
// identity, through lat::detail::certified_fmpz_lll, whose route it returns.
lat::detail::CertifiedLllRoute ideal_lattice_rounded_embedding_integer_lll(
        flint::FmpzMat& embedded,
        flint::FmpzMat& transform) noexcept;

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
