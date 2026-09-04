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

bool copy_ideal_basis_row_coordinates(flint::FmpzMat& out,
                                      flint::FmpzMatConstRef basis,
                                      slong row) noexcept;

bool build_ideal_lattice_lll_data(IdealLatticeLllData& out,
                                  const Ideal& ideal,
                                  slong precision) noexcept;

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
