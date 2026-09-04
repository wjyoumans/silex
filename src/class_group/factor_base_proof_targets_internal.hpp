#pragma once

#include <silex/prime_ideal.hpp>

#include <vector>

namespace silex::detail {

// `decomposition` must be the complete, unfiltered decomposition above
// `rational_prime`.  The selector preserves its order and publishes `out`
// only after every retained entry has been validated.
bool select_factor_base_proof_targets(
        std::vector<slong>& out,
        const PrimeIdealList& decomposition,
        flint::FmpzConstRef rational_prime,
        flint::FmpzConstRef bound) noexcept;

}  // namespace silex::detail
