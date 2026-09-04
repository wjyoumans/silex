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
