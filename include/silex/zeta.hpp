#pragma once

#include <flint/flint.h>

#include <silex/flint/arb.hpp>
#include <silex/order.hpp>

#include <optional>

namespace silex {

struct ZetaBfAuditResult {
    flint::Arb value;
    flint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;
};

bool zeta_residue(flint::ArbRef out,
                  const Order& order,
                  slong precision) noexcept;

bool zeta_log_residue(flint::ArbRef out,
                      const Order& order,
                      slong precision) noexcept;

bool zeta_log_residue_bf(flint::ArbRef out,
                         const Order& order,
                         ulong max_cutoff,
                         slong precision) noexcept;

bool zeta_log_residue_bf_audit(flint::ArbRef out,
                               flint::ArbRef error_bound,
                               ulong& cutoff,
                               slong& work_precision,
                               const Order& order,
                               ulong max_cutoff,
                               slong precision) noexcept;
std::optional<ZetaBfAuditResult> zeta_log_residue_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept;

bool zeta_residue_bf(flint::ArbRef out,
                     const Order& order,
                     ulong max_cutoff,
                     slong precision) noexcept;

bool zeta_residue_bf_audit(flint::ArbRef out,
                           flint::ArbRef error_bound,
                           ulong& cutoff,
                           slong& work_precision,
                           const Order& order,
                           ulong max_cutoff,
                           slong precision) noexcept;
std::optional<ZetaBfAuditResult> zeta_residue_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept;

bool zeta_class_regulator_product(flint::ArbRef out,
                                  const Order& order,
                                  slong precision) noexcept;

bool zeta_class_regulator_product_bf(flint::ArbRef out,
                                     const Order& order,
                                     ulong max_cutoff,
                                     slong precision) noexcept;

bool zeta_class_regulator_product_bf_audit(flint::ArbRef out,
                                           flint::ArbRef error_bound,
                                           ulong& cutoff,
                                           slong& work_precision,
                                           const Order& order,
                                           ulong max_cutoff,
                                           slong precision) noexcept;
std::optional<ZetaBfAuditResult> zeta_class_regulator_product_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept;

}  // namespace silex
