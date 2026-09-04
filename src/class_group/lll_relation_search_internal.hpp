#pragma once

#include "ideal_lattice_lll_internal.hpp"

namespace silex::detail::relation_search {

bool select_small_lll_rows(
        std::vector<slong>& selected,
        const IdealLatticeLllData& context,
        const Ideal& ideal,
        const Order& order) noexcept;

bool set_small_lll_next_candidate(
        flint::FmpzMat& out,
        slong& count_after_next,
        const IdealLatticeLllData& context,
        const std::vector<slong>& selected) noexcept;

bool build_random_factor_base_product(
        Ideal& out,
        flint::FmpzMatRef factor_base_row,
        const FactorBase& base,
        slong sequence_index,
        ulong seed,
        const DiagnosticsContext* diagnostics) noexcept;

}  // namespace silex::detail::relation_search
