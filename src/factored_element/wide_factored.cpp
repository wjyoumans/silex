#include "wide_factored_internal.hpp"

#include <silex/order_unit.hpp>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace silex::detail {
namespace {

bool copy_base(std::vector<Element>& out, const Element& input) noexcept {
    if (input.parent() == nullptr) {
        return false;
    }
    out.emplace_back(*input.parent());
    return out.back().is_defined() && out.back().set(input);
}

}  // namespace

void WideFactoredUnitsDelete::operator()(
        WideFactoredUnits* units) const noexcept {
    delete units;
}

bool wide_factored_units_copy(WideFactoredUnits& out,
                              const WideFactoredUnits& input) noexcept {
    WideFactoredUnits copy;
    copy.bases.reserve(input.bases.size());
    for (const Element& base : input.bases) {
        if (!copy_base(copy.bases, base)) {
            return false;
        }
    }
    flint::FmpzMat exponents(flint::fmpz_mat_nrows(input.exponents),
                             flint::fmpz_mat_ncols(input.exponents));
    fmpz_mat_set(exponents.raw(), input.exponents.raw());
    copy.exponents.swap(exponents);
    out.bases.swap(copy.bases);
    out.exponents.swap(copy.exponents);
    return true;
}

bool wide_factored_units_exceed_slong(const WideFactoredUnits& units) noexcept {
    const slong rows = flint::fmpz_mat_nrows(units.exponents);
    const slong columns = flint::fmpz_mat_ncols(units.exponents);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < columns; ++j) {
            if (fmpz_fits_si(fmpz_mat_entry(units.exponents.raw(), i, j)) ==
                0) {
                return true;
            }
        }
    }
    return false;
}

bool wide_factored_units_transform(WideFactoredUnits& out,
                                   std::span<const FactoredElement> units,
                                   flint::FmpzMatConstRef transform) noexcept {
    if (units.empty() ||
        units.size() >
                static_cast<std::size_t>(std::numeric_limits<slong>::max()) ||
        units.front().parent() == nullptr) {
        return false;
    }
    const slong rank = static_cast<slong>(units.size());
    if (flint::fmpz_mat_nrows(transform) != rank ||
        flint::fmpz_mat_ncols(transform) != rank) {
        return false;
    }
    const NumberField& field = *units.front().parent();

    // Merge equal bases as reference FacElem does: column k of `merged`
    // collects the exponents of distinct base k over every input unit.
    std::vector<const Element*> distinct;
    std::vector<std::vector<std::pair<slong, slong>>> entries(units.size());
    for (std::size_t row = 0; row < units.size(); ++row) {
        const FactoredElement& unit = units[row];
        if (!unit.is_defined() || unit.parent() == nullptr ||
            !unit.parent()->has_same_data(field)) {
            return false;
        }
        for (const auto& entry : unit.factors()) {
            if (entry.exponent == 0) {
                continue;
            }
            if (!entry.factor.is_defined() ||
                !entry.factor.has_parent(field)) {
                return false;
            }
            std::size_t column = 0;
            while (column < distinct.size() &&
                   !distinct[column]->equal(entry.factor)) {
                ++column;
            }
            if (column == distinct.size()) {
                distinct.push_back(&entry.factor);
            }
            entries[row].emplace_back(static_cast<slong>(column),
                                      entry.exponent);
        }
    }

    const slong columns = static_cast<slong>(distinct.size());
    flint::FmpzMat merged(rank, columns);
    for (slong row = 0; row < rank; ++row) {
        for (const auto& [column, exponent] :
             entries[static_cast<std::size_t>(row)]) {
            fmpz* target = fmpz_mat_entry(merged.raw(), row, column);
            if (exponent >= 0) {
                fmpz_add_ui(target, target, static_cast<ulong>(exponent));
            } else {
                // -(exponent + 1) + 1 avoids negating WORD_MIN.
                fmpz_sub_ui(target, target,
                            static_cast<ulong>(-(exponent + 1)));
                fmpz_sub_ui(target, target, 1);
            }
        }
    }

    flint::FmpzMat product(rank, columns);
    fmpz_mat_mul(product.raw(), transform.raw(), merged.raw());

    std::vector<slong> used;
    for (slong column = 0; column < columns; ++column) {
        for (slong row = 0; row < rank; ++row) {
            if (!fmpz_is_zero(fmpz_mat_entry(product.raw(), row, column))) {
                used.push_back(column);
                break;
            }
        }
    }

    WideFactoredUnits result;
    result.bases.reserve(used.size());
    flint::FmpzMat exponents(rank, static_cast<slong>(used.size()));
    for (std::size_t k = 0; k < used.size(); ++k) {
        if (!copy_base(result.bases,
                       *distinct[static_cast<std::size_t>(used[k])])) {
            return false;
        }
        for (slong row = 0; row < rank; ++row) {
            fmpz_set(fmpz_mat_entry(exponents.raw(), row,
                                    static_cast<slong>(k)),
                     fmpz_mat_entry(product.raw(), row, used[k]));
        }
    }
    result.exponents.swap(exponents);
    out.bases.swap(result.bases);
    out.exponents.swap(result.exponents);
    return true;
}

}  // namespace silex::detail
