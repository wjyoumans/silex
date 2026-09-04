#pragma once

#include <silex/diagnostics.hpp>
#include <silex/zeta.hpp>

#include <cstddef>
#include <vector>

namespace silex {
class Element;
class FactorBase;
}

namespace silex::detail {

struct ZetaBfResidueDegreeCacheEntry {
    ulong p = 0;
    std::size_t offset = 0;
    std::size_t length = 0;
};

struct ZetaBfResidueDegreeCache {
    std::vector<ZetaBfResidueDegreeCacheEntry> entries;
    std::vector<slong> residue_degrees;
    std::size_t lookup_hint = 0;
};

bool grh_factor_base_bound_with_diagnostics(
        flint::FmpzRef out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept;

bool zeta_class_regulator_product_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_class_regulator_product_bf_audit_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_class_regulator_product_validation_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_bf_audit_cutoff_available(const Order& order,
                                    ulong max_cutoff,
                                    slong precision) noexcept;

bool class_regulator_product_estimate_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        flint::Fmpz* torsion_order = nullptr,
        Element* torsion_generator = nullptr) noexcept;

}  // namespace silex::detail
