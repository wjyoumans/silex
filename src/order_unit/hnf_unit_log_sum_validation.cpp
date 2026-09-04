#include "order_unit_internal.hpp"

namespace silex::detail {

namespace {

bool arb_abs_upper_lt_2exp(const arb_struct* value,
                           slong exponent,
                           slong precision) noexcept {
    if (precision <= 0 || ::arb_is_finite(value) == 0) {
        return false;
    }

    flint::Arf upper;
    ::arb_get_abs_ubound_arf(upper.raw(), value, precision);
    return flint::arf_is_finite(upper) &&
           ::arf_cmpabs_2exp_si(upper.raw(), exponent) < 0;
}

}  // namespace

bool unit_log_row_sums_are_small(const flint::ArbMat& logs,
                                         slong precision) noexcept {
    const slong rank = flint::arb_mat_nrows_value(logs);
    const slong places = flint::arb_mat_ncols_value(logs);
    if (rank <= 0 || places <= 0 || precision <= 0) {
        return false;
    }

    for (slong row = 0; row < rank; ++row) {
        flint::Arb sum;
        flint::arb_zero(sum);
        for (slong place = 0; place < places; ++place) {
            flint::arb_add(sum, sum, arb_mat_entry(logs.raw(), row, place),
                           precision);
        }
        // reference validate_unit_log_sums rejects gexpo(sum(real_i(x))) > -10.
        // Since exponent -10 means |sum| < 2^-9, require that bound.
        if (!arb_abs_upper_lt_2exp(sum.raw(), -9, precision)) {
            return false;
        }
    }

    return true;
}

bool unit_log_row_sums_are_small(
        const std::vector<FactoredElement>& units,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    const slong rank = static_cast<slong>(units.size());
    slong places = 0;
    if (rank <= 0 || precision <= 0 ||
        !compact_places(places, embeddings)) {
        return false;
    }

    flint::ArbMat logs(rank, places);
    return compact_log_matrix(
                   logs, embeddings,
                   FactoredElementSpan(units.data(), units.size()),
                   precision) &&
           unit_log_row_sums_are_small(logs, precision);
}

}  // namespace silex::detail
