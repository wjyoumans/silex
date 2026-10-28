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

// out = left * right, the product fmpz_mat_mul returns, computed by
// accumulating right's rows over the nonzero entries of left only.  The ideal
// product forms H M with H an ideal HNF, which in large-degree fields is
// mostly zero, so this costs (nonzeros of H) * n multiply-adds instead of a
// dense n^3 product.  Shapes must agree; out must not alias left or right.
void zero_skipping_mul(flint::FmpzMatRef out, flint::FmpzMatConstRef left,
                       flint::FmpzMatConstRef right) noexcept;

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

// A verified two-element presentation of a nonscalar integral ideal I of a
// maximal order, the one try_ideal_product_hnf searches for when I is the
// factor it extracts from: `hnf` is the HNF of I, I = content K with K
// primitive, K = minimum O + beta O with `minimum` the least positive
// integer in K, and `multiplication` is the full (unreduced) multiplication
// matrix of beta.
struct IdealTwoGeneratorForm {
    flint::FmpzMat hnf{0, 0};
    flint::Fmpz content;
    flint::Fmpz minimum;
    flint::FmpzMat multiplication{0, 0};
};

// Runs the search try_ideal_product_hnf runs on a factor (extract_two_generator
// on the primitive part, with these options) once, so that later products
// with the ideal can reuse it.  Fails, leaving out unchanged, when the order
// is not maximal, the ideal is scalar (its products take the scalar branch,
// which needs no form) or the search is exhausted.
bool ideal_two_generator_form(IdealTwoGeneratorForm& out, const Ideal& ideal,
                              const TwoGeneratorSearchOptions& options = {},
                              TwoGeneratorSearchReport* report = nullptr) noexcept;

// The HNF try_ideal_product_hnf(out, left, right) produces, with the
// two-element form of right given instead of searched for.  right_form must
// come from ideal_two_generator_form(right): its stored HNF is compared with
// right's, and a mismatch fails.  The scalar and principal-generator branches
// of try_ideal_product_hnf, which search nothing, are taken through it.
// Leaves out unchanged on failure.
bool try_ideal_product_hnf_with_form(flint::FmpzMatRef out, const Ideal& left,
                                     const Ideal& right,
                                     const IdealTwoGeneratorForm& right_form) noexcept;

// The HNF of the product of two ideals of one order, from the n^2 basis
// products formed with the order's cached multiplication table and reduced by
// fmpz_mat_hnf_modular_eldiv modulo the product of the minima of the primitive
// parts (reference `Arithmetic.jl:mul_gen`).  Equals the exact HNF of those
// products; validity does not depend on maximality.  Leaves out unchanged on
// failure.
bool basis_product_hnf(flint::FmpzMatRef out, const Ideal& left,
                       const Ideal& right) noexcept;

}  // namespace silex::detail
