#pragma once

#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled experiment/test controls. The public operation uses these defaults.
struct TwoGeneratorSearchOptions {
    slong random_trials = 32;
    bool basis_candidates = true;
};

struct TwoGeneratorSearchReport {
    slong basis_trials = 0;
    slong random_trials = 0;
    bool verified = false;
};

// Least positive integer in the ideal whose basis is the rows of the n x n
// matrix h, given the 1 x n order coordinates of 1: the least m > 0 with
// m x integral, where x h = coordinates. ideal_minimum uses forward
// substitution when h is upper triangular with a nonzero diagonal (every
// stored ideal HNF) and the generic solve otherwise; ideal_minimum_generic is
// the fmpz_mat_solve route, kept as the reference. Both return false for a
// singular h; ideal_minimum also rejects mismatched shapes.
bool ideal_minimum(flint::Fmpz& minimum, flint::FmpzMatConstRef coordinates,
                   flint::FmpzMatConstRef h) noexcept;
bool ideal_minimum_generic(flint::Fmpz& minimum,
                           flint::FmpzMatConstRef coordinates,
                           flint::FmpzMatConstRef h) noexcept;

// On failure, scalar and element remain unchanged. Search exhaustion is not a
// claim that the ideal has no two-generator presentation.
bool ideal_two_generator(flint::FmpzRef scalar, OrderElement& element,
                         const Ideal& ideal,
                         const TwoGeneratorSearchOptions& options = {},
                         TwoGeneratorSearchReport* report = nullptr) noexcept;

// Produces a verified product HNF, or leaves out unchanged for generic fallback.
bool try_ideal_product_hnf(flint::FmpzMatRef out,
                           const Ideal& left, const Ideal& right,
                           const TwoGeneratorSearchOptions& options = {},
                           TwoGeneratorSearchReport* report = nullptr) noexcept;

}  // namespace silex::detail
