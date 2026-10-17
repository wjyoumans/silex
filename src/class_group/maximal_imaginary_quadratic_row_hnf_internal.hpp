#pragma once

#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order.hpp>

namespace silex::detail {

// Route that builds the row ideal of the maximal imaginary quadratic witness
// (src/class_group/maximal_imaginary_quadratic_relation_backend.cpp).
enum class QuadraticRowHnfRoute {
    // Hensel lifts of the prime roots combined by CRT into the primitive
    // part A Z + (omega - s) Z, times its content; no fallback.
    binary_form,
    // Ideal powers and products of the factor-base primes (the fallback).
    ideal_product,
};

// Test access to the witness row ideal: for a maximal imaginary quadratic
// order, a factor base the relation backend accepts, and a 1 x base.length()
// exponent row, sets hnf to the 2 x 2 HNF in the order basis (1, omega) of
// denominator * prod P_i^(row_i), norm to its norm, and denominator to the
// positive integer that the witness divides by.  Returns false when the
// factor base is not accepted or when the selected route fails.
bool maximal_imaginary_quadratic_row_hnf(flint::FmpzMat& hnf,
                                         flint::Fmpz& norm,
                                         flint::Fmpz& denominator,
                                         const Order& order,
                                         const FactorBase& base,
                                         flint::FmpzMatConstRef row,
                                         QuadraticRowHnfRoute route) noexcept;

class TestHookAccess;

// Test entry of the relation backend; only detail::TestHookAccess
// (test/test_hook_access.hpp) can call it.
class MaximalImaginaryQuadraticRowTestEntry {
    friend class TestHookAccess;

    // The proven witness generator of the row ideal, as the relation
    // backend computes it for a principal row (out generates
    // prod P_i^(row_i)).  With use_row_hnf_builder false the binary-form
    // builder is skipped, as when it declines, so the row HNF comes from the
    // Ideal route (the fallback).  Returns false when the factor base is not
    // accepted or the witness step fails.
    static bool row_generator(
            Element& out,
            const Order& order,
            const FactorBase& base,
            flint::FmpzMatConstRef row,
            bool use_row_hnf_builder) noexcept;
};

}  // namespace silex::detail
