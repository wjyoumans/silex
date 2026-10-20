#pragma once

#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order.hpp>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache {
    Order order;
    slong embedding_precision = 0;
    flint::ArbMat order_embedding_rows{0, 0};
    bool has_order = false;
    bool has_order_embedding_rows = false;
    // Smallest rounding scale at which a T2 enumeration basis of an ideal of
    // `order` was accepted (ideal_t2_enumeration.cpp), or 0 when none was.
    // The next ideal's scale loop starts there.  It lowers only when the rerun
    // from scale 4 accepts a coarser scale.  Cleared when `order` changes.
    slong t2_rounding_scale_hint = 0;
};

bool multiply_integer_arb_matrices(flint::ArbMat& out,
                                   flint::FmpzMatConstRef left,
                                   const flint::ArbMat& right,
                                   slong precision) noexcept;

bool build_ideal_minkowski_embedding_rows(
        flint::ArbMat& out,
        flint::FmpzMatConstRef basis,
        const Order& order,
        slong precision,
        OrderMinkowskiEmbeddingCache* cache = nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace silex::detail
