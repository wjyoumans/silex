#pragma once

#include <silex/class_group.hpp>
#include <silex/flint/arb_mat.hpp>

#include "../lat/certified_lll_internal.hpp"

#include <vector>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

}  // namespace silex::detail

namespace silex::detail::relation_search {

struct IdealLatticeLllData {
    flint::FmpzMat basis{0, 0};
    flint::FmpzMat scaled_gram{0, 0};
    flint::Fmpz gram_denominator;
};

// One entry of the two-LLL order data (ideal_lattice_two_lll_order_gram):
// the order's Minkowski rows at `embedding_precision` and the scaled order
// T2 Gram numerator formed from them at `precision`.
struct IdealLatticeTwoLllOrderGram {
    slong embedding_precision = 0;
    slong precision = 0;
    flint::ArbMat source_rows{0, 0};
    flint::FmpzMat order_gram{0, 0};
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
// `two_lll_order_grams` holds the order data of that two-LLL route (see
// IdealLatticeTwoLllOrderGram), cleared with the rest when the order changes.
// `two_lll_warm_precision` is the precision at which the last two-LLL build
// of the order (build_ideal_lattice_lll_data) succeeded, 0 if none: later
// builds of the order start their schedule there instead of at the lower
// requested precision, skipping attempts that failed the rounding checks
// for an earlier ideal.  Every result is still accepted by the same checks
// at the precision it was built at; an ideal that would have passed at a
// lower precision gets the (equally valid) reduced basis of the warm one.
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
    slong two_lll_warm_precision = 0;
    std::vector<IdealLatticeTwoLllOrderGram> two_lll_order_grams;
};

// Keys `cache` on `order`: when the cache holds no order or another order's
// data, every cached field is cleared and the order stored.  A cache already
// keyed on `order` is left as is.
void bind_ideal_lattice_rounded_embedding_cache(
        IdealLatticeRoundedEmbeddingCache& cache,
        const Order& order) noexcept;

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

// The two-LLL route (reference NfOrd/LLL.jl:_lll).  With a non-null `cache`
// whose order is the ideal's, the order T2 Gram numerator of each attempt is
// taken from `cache->two_lll_order_grams` instead of being recomputed; see
// ideal_lattice_two_lll_order_gram.  Such a cache also warm-starts the
// precision schedule at `cache->two_lll_warm_precision` when that exceeds
// the starting precision, and records the precision of a successful build
// there.  The result equals that of a build without a cache started at the
// same precision.
bool build_ideal_lattice_lll_data(
        IdealLatticeLllData& out,
        const Ideal& ideal,
        slong precision,
        IdealLatticeRoundedEmbeddingCache* cache = nullptr) noexcept;

// The order data of one two-LLL attempt at (embedding_precision, precision):
// the order's Minkowski rows from detail::order_minkowski_embedding_rows at
// `embedding_precision`, multiplied by the identity basis at that precision,
// and the scaled order T2 Gram numerator of those rows at scale and
// precision `precision`.  That Gram is a function of the order rows and the
// two precisions, so an entry of `cache->two_lll_order_grams` whose
// precisions match and whose `source_rows` equal the order's current rows
// as balls (midpoint and radius) holds exactly the Gram a recomputation
// gives, and is returned; otherwise the Gram is computed and, on success,
// stored (at most kIdealLatticeTwoLllOrderGramEntries entries, the oldest
// replaced first).  A null `cache`, or one for another order, recomputes.
// Returns false when the rows cannot be formed or do not round at
// `precision`; exposed for tests.
inline constexpr slong kIdealLatticeTwoLllOrderGramEntries = 4;
bool ideal_lattice_two_lll_order_gram(
        flint::FmpzMat& out,
        const Order& order,
        slong embedding_precision,
        slong precision,
        IdealLatticeRoundedEmbeddingCache* cache) noexcept;

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
