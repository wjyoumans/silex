#include "ideal_minkowski_embedding_internal.hpp"

#include <silex/diagnostics.hpp>

namespace silex::detail {

bool multiply_integer_arb_matrices(flint::ArbMat& out,
                                   flint::FmpzMatConstRef left,
                                   const flint::ArbMat& right,
                                   slong precision) noexcept {
    const slong rows = flint::fmpz_mat_nrows(left);
    const slong middle = flint::fmpz_mat_ncols(left);
    const slong cols = flint::arb_mat_ncols_value(right);
    if (precision <= 0 || flint::arb_mat_nrows_value(right) != middle ||
        flint::arb_mat_nrows_value(out) != rows ||
        flint::arb_mat_ncols_value(out) != cols) {
        return false;
    }

    for (slong row = 0; row < rows; ++row) {
        for (slong column = 0; column < cols; ++column) {
            flint::ArbRef out_entry =
                    flint::arb_mat_entry_ref(out, row, column);
            flint::arb_zero(out_entry);
            for (slong k = 0; k < middle; ++k) {
                flint::arb_addmul_fmpz(
                        out_entry,
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(right), k, column),
                        flint::fmpz_mat_entry(left, row, k),
                        precision);
            }
        }
    }
    return true;
}

bool build_ideal_minkowski_embedding_rows(
        flint::ArbMat& out,
        flint::FmpzMatConstRef basis,
        const Order& order,
        slong precision,
        OrderMinkowskiEmbeddingCache* cache,
        const DiagnosticsContext* diagnostics) noexcept {
    const NumberField* field = order.parent();
    const slong degree = order.degree();
    if (field == nullptr || precision <= 0 ||
        flint::fmpz_mat_nrows(basis) != degree ||
        flint::fmpz_mat_ncols(basis) != degree ||
        flint::arb_mat_nrows_value(out) != degree ||
        flint::arb_mat_ncols_value(out) != degree) {
        return false;
    }

    flint::ArbMat local_order_embedding_rows(0, 0);
    flint::ArbMat* order_embedding_rows = nullptr;
    const bool cache_hit =
            cache != nullptr && cache->has_order_embedding_rows &&
            cache->has_order && cache->order.has_same_data(order) &&
            cache->embedding_precision == precision &&
            flint::arb_mat_nrows_value(cache->order_embedding_rows) == degree &&
            flint::arb_mat_ncols_value(cache->order_embedding_rows) == degree;
    if (cache_hit) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_minkowski_embedding.cache_hit");
        order_embedding_rows = &cache->order_embedding_rows;
    } else {
        if (cache == nullptr) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.ideal_minkowski_embedding.cache_unavailable");
        } else if (!cache->has_order_embedding_rows) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.ideal_minkowski_embedding.cache_empty");
        } else if (!cache->has_order ||
                   !cache->order.has_same_data(order)) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.ideal_minkowski_embedding.cache_order_miss");
        } else if (cache->embedding_precision != precision) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.ideal_minkowski_embedding.cache_precision_miss");
        } else {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.ideal_minkowski_embedding.cache_shape_miss");
        }
        if (cache != nullptr) {
            cache->has_order_embedding_rows = false;
            cache->order_embedding_rows = flint::ArbMat(degree, degree);
            order_embedding_rows = &cache->order_embedding_rows;
        } else {
            local_order_embedding_rows = flint::ArbMat(degree, degree);
            order_embedding_rows = &local_order_embedding_rows;
        }

        if (!order_minkowski_embedding_rows(
                    flint::ArbMatRef(*order_embedding_rows), order, precision,
                    diagnostics)) {
            return false;
        }
        if (cache != nullptr) {
            cache->order = order;
            cache->embedding_precision = precision;
            cache->has_order = true;
            cache->has_order_embedding_rows = true;
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_minkowski_embedding.multiply");
        return multiply_integer_arb_matrices(
                out, basis, *order_embedding_rows, precision);
    }
}

}  // namespace silex::detail
