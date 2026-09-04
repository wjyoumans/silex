#include "order_unit_internal.hpp"

namespace silex::detail {

namespace {

bool relation_product_from_coefficients(
        FactoredElement& out,
        const ClassGroupContext& class_group,
        flint::FmpzMatConstRef coefficients,
        slong row,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = class_group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field) || row < 0 ||
        row >= flint::fmpz_mat_nrows(coefficients) ||
        flint::fmpz_mat_ncols(coefficients) !=
                class_group.relation_count()) {
        return false;
    }

    FactoredElement factored(*field);
    Element generator(*field);
    bool all_exponents_fit = true;
    if (!factored.is_defined() || !generator.is_defined() ||
        !factored.one()) {
        return false;
    }
    for (slong i = 0; i < class_group.relation_count(); ++i) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(coefficients, row, i);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!class_group.relation_generator(generator, i)) {
            return false;
        }
        if (!flint::fmpz_fits_si(exponent)) {
            all_exponents_fit = false;
            break;
        }
        if (!factored.push(generator, flint::fmpz_get_si(exponent))) {
            return false;
        }
    }
    if (all_exponents_fit) {
        factored.normalize();
        out.swap(factored);
        return true;
    }

    SILEX_LOG(diagnostics, DiagnosticsModule::unit_group, LogLevel::detail,
              "transform witnesses failed: combined exponent exceeds slong");
    return false;
}

}  // namespace

bool factored_units_from_relation_coefficients(
        std::vector<FactoredElement>& out,
        const ClassGroupContext& class_group,
        const flint::FmpzMat& coefficients,
        const DiagnosticsContext* diagnostics) noexcept {
    out.clear();
    const Order* order = class_group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = flint::fmpz_mat_nrows(coefficients);
    if (field == nullptr || rank < 0 ||
        flint::fmpz_mat_ncols(coefficients) != class_group.relation_count()) {
        return false;
    }

    out.reserve(static_cast<std::size_t>(rank));
    for (slong row = 0; row < rank; ++row) {
        out.emplace_back(*field);
        FactoredElement& unit = out.back();
        if (!unit.is_defined() ||
            !relation_product_from_coefficients(
                    unit, class_group, flint::FmpzMatConstRef(coefficients),
                    row, diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "transform witnesses failed: build combined witness");
            return false;
        }
        unit.normalize();
    }
    return true;
}

}  // namespace silex::detail
