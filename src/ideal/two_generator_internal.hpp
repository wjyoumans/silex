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
