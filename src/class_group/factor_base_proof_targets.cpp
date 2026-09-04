#include "factor_base_proof_targets_internal.hpp"

#include <cstddef>
#include <vector>

namespace silex::detail {

bool select_factor_base_proof_targets(
        std::vector<slong>& out,
        const PrimeIdealList& decomposition,
        flint::FmpzConstRef rational_prime,
        flint::FmpzConstRef bound) noexcept {
    if (flint::fmpz_sgn(rational_prime) <= 0 ||
        flint::fmpz_sgn(bound) < 0) {
        return false;
    }

    std::vector<slong> selected;
    selected.reserve(static_cast<std::size_t>(decomposition.size()));
    flint::Fmpz norm;
    for (slong i = 0; i < decomposition.size(); ++i) {
        const PrimeIdeal* prime = decomposition.at(i);
        if (prime == nullptr || prime->residue_degree() <= 0) {
            return false;
        }

        flint::fmpz_pow_ui(
                flint::FmpzRef(norm), rational_prime,
                static_cast<ulong>(prime->residue_degree()));
        if (flint::fmpz_cmp(flint::FmpzConstRef(norm), bound) <= 0) {
            selected.push_back(i);
        }
    }

    // The principal (p) relation alone discharges a final e = 1 factor when
    // every other factor is also a proof target.  Without an independently
    // established subgroup witness, retain every selected factor in a
    // norm-truncated decomposition.
    if (selected.size() == static_cast<std::size_t>(decomposition.size()) &&
        selected.size() > 1) {
        const PrimeIdeal* last = decomposition.at(selected.back());
        if (last != nullptr && last->ramification_index() == 1) {
            selected.pop_back();
        }
    }

    out.swap(selected);
    return true;
}

}  // namespace silex::detail
