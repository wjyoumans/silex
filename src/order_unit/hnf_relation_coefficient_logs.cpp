#include "order_unit_internal.hpp"

#include <silex/archimedean.hpp>

namespace silex::detail {

bool relation_coefficients_log_matrix(
        flint::ArbMat& out,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        const flint::FmpzMat& coefficients,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = class_group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = flint::fmpz_mat_nrows(coefficients);
    const slong relation_count = class_group.relation_count();
    slong places = 0;
    if (field == nullptr || precision <= 0 ||
        flint::fmpz_mat_ncols(coefficients) != relation_count ||
        !compact_places(places, embeddings) ||
        flint::arb_mat_nrows_value(out) != rank ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    for (slong row = 0; row < rank; ++row) {
        for (slong place = 0; place < places; ++place) {
            flint::arb_zero(flint::arb_mat_entry_ref(out, row, place));
        }
    }

    Element generator(*field);
    flint::ArbVec relation_logs(places);
    if (!generator.is_defined()) {
        return false;
    }
    for (slong rel = 0; rel < relation_count; ++rel) {
        bool used = false;
        for (slong row = 0; row < rank; ++row) {
            if (!flint::fmpz_is_zero(
                        flint::fmpz_mat_entry(
                                flint::FmpzMatConstRef(coefficients),
                                row, rel))) {
                used = true;
                break;
            }
        }
        if (!used) {
            continue;
        }

        if (!class_group.relation_generator(generator, rel) ||
            !silex::logarithmic_embedding(
                    flint::ArbVecRef(relation_logs), embeddings, generator,
                    LogEmbeddingMode::product, precision, diagnostics)) {
            return false;
        }

        for (slong row = 0; row < rank; ++row) {
            const flint::FmpzConstRef exponent =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(coefficients), row, rel);
            if (flint::fmpz_is_zero(exponent)) {
                continue;
            }
            for (slong place = 0; place < places; ++place) {
                flint::arb_addmul_fmpz(
                        flint::arb_mat_entry_ref(out, row, place),
                        flint::ArbConstRef(relation_logs.data() + place),
                        exponent, precision);
            }
        }
    }
    return true;
}

}  // namespace silex::detail
