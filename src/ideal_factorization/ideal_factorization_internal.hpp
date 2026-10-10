#pragma once

#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex {

class FactorBase;
class Ideal;
class OrderElement;
class PrimeIdeal;
struct DiagnosticsContext;

namespace detail {

enum class OneLargePrimeFactorStatus {
    no_candidate,
    found,
};

bool order_element_factor_over_base_with_one_large_prime(
        OneLargePrimeFactorStatus& status,
        flint::FmpzMatRef exponents,
        PrimeIdeal& large_prime,
        const OrderElement& element,
        flint::FmpqConstRef norm,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool ideal_factor_over_base_with_required_prime(
        bool& matches,
        const Ideal& ideal,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

// Records which path classified the element; for tests, independent of
// profiling builds.
struct RequiredPrimeDirectAudit {
    bool rejected_by_cofactor = false;
};

bool order_element_factor_over_base_with_required_prime(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics = nullptr,
        RequiredPrimeDirectAudit* audit = nullptr) noexcept;

// Per-search data for applying the cofactor rejection of
// order_element_factor_over_base_with_required_prime to many candidates of
// one order, factor base and required prime, given by their order-basis
// coordinates, without building an element. Ready only when every
// precondition of the direct path holds, every factor-base rational prime
// block is readable, and the degree is 2 to 4; otherwise callers keep the
// direct path. Scratch members make it single-threaded.
struct RequiredPrimeNormPrefilter {
    bool ready = false;
    slong degree = 0;
    // Product of the factor-base rational primes and the required prime's
    // rational prime.
    flint::Fmpz smooth_product;
    // Row i * degree + j holds the order coordinates of w_i * w_j.
    flint::FmpzMat multiplication_table{0, 0};
    flint::FmpzMat multiplication_matrix{0, 0};
    flint::Fmpz norm;
    flint::Fmpz cofactor;
    flint::Fmpz divisor;
};

bool required_prime_norm_prefilter_init(
        RequiredPrimeNormPrefilter& prefilter,
        const FactorBase& base,
        const PrimeIdeal& required_prime) noexcept;

// The exact norm of the order element with these coordinates (a 1 x degree
// row), from the multiplication table. Requires a ready prefilter.
bool required_prime_norm_prefilter_norm(
        flint::FmpzRef out,
        RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept;

// Sets rejected exactly when the direct path would reject the element by its
// norm cofactor: the norm is nonzero and not smooth over the primes of
// smooth_product. A false rejection leaves the decision to the direct path.
bool required_prime_norm_prefilter_rejects(
        bool& rejected,
        RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept;

// Evaluation success is separate from classification. Requires a nonzero
// integral scalar and a known maximal order; input objects are not modified.
bool scalar_factor_over_base_with_required_prime(
        bool& matches,
        flint::FmpzConstRef scalar,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace detail
}  // namespace silex
