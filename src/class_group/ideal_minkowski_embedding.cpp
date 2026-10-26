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

    if (rows == 0 || cols == 0) {
        return true;
    }
    if (middle == 0) {
        for (slong row = 0; row < rows; ++row) {
            for (slong column = 0; column < cols; ++column) {
                flint::arb_zero(flint::arb_mat_entry_ref(out, row, column));
            }
        }
        return true;
    }

    // Each entry is one arb_dot_fmpz of a row of `left` (contiguous) with a
    // column of `right`, read with the row stride of `right` (an Arb matrix
    // stores its rows at a uniform stride).  Like every Arb operation,
    // arb_dot_fmpz returns a ball that contains the exact sum
    // sum_k left[row][k] * x_k for every choice of x_k in the balls
    // right[k][column], so the product stays a rigorous enclosure.  Its
    // midpoint comes from FLINT's dot-product accumulation with one final
    // rounding rather than from a sum of rounded arb_addmul_fmpz steps, so an
    // entry can differ from that sum by rounding noise.
    const arb_struct* right_first = arb_mat_entry(right.raw(), 0, 0);
    const slong right_step =
            middle > 1 ? static_cast<slong>(
                                 arb_mat_entry(right.raw(), 1, 0) - right_first)
                       : 1;
    for (slong row = 0; row < rows; ++row) {
        const fmpz* left_row = flint::fmpz_mat_entry(left, row, 0).raw();
        for (slong column = 0; column < cols; ++column) {
            ::arb_dot_fmpz(
                    flint::arb_mat_entry_ref(out, row, column).raw(), nullptr,
                    0, right_first + column, right_step, left_row, 1, middle,
                    precision);
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
            if (!cache->has_order || !cache->order.has_same_data(order)) {
                cache->t2_rounding_scale_hint = 0;
            }
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
