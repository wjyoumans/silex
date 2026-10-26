#pragma once

#include <cstddef>
#include <vector>

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/prime_ideal.hpp>

namespace silex::detail {

// Returns the image of the NumberField defining generator, even when a direct
// maximal-quadratic prime stores its residue polynomial in the integral omega
// basis used by the order.
bool degree_one_prime_root_mod_p(flint::Fmpz& root,
                                 flint::Fmpz& p,
                                 const PrimeIdeal& prime) noexcept;
bool evaluate_fmpq_poly_at_degree_one_prime(
        flint::Fmpz& out,
        flint::FmpqPolyConstRef polynomial,
        flint::FmpzConstRef root,
        flint::FmpzConstRef p) noexcept;
// Smallest primitive root modulo the prime p, in fmpz arithmetic for any p.
// Returns false when p is not prime.
bool primitive_root_mod_prime(flint::Fmpz& out,
                              flint::FmpzConstRef p) noexcept;
// The same smallest primitive root for a prime p below the word size,
// factoring p - 1 once.  The caller guarantees that p is prime; returns 0
// only for p < 2.
ulong primitive_root_mod_prime_ui(ulong p) noexcept;
bool quotient_log_mod_prime_setup(flint::Fmpz& cofactor,
                                  flint::Fmpz& quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept;
bool quotient_log_mod_prime_apply(flint::FmpzRef out,
                                  flint::FmpzConstRef value,
                                  flint::FmpzConstRef cofactor,
                                  flint::FmpzConstRef quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept;
bool factored_value_at_degree_one_root(flint::Fmpz& out,
                                       const FactoredElement& element,
                                       flint::FmpzConstRef p,
                                       flint::FmpzConstRef root) noexcept;
// factored_value_at_degree_one_root in word arithmetic for a prime p below
// the word size: fails exactly where it does and otherwise gives the same
// residue.  Each factor is evaluated from its stored form num(x) / den in
// place (no polynomial copy), and the per-coefficient route of
// evaluate_fmpq_poly_at_degree_one_prime is taken only for a factor whose
// stored denominator is divisible by p.
bool factored_value_at_degree_one_root_nmod(ulong& out,
                                           const FactoredElement& element,
                                           ulong p,
                                           ulong root) noexcept;
// The reductions mod a word-size prime p that
// factored_value_at_degree_one_root_nmod makes for each factor of a
// FactoredElement (its stored denominator and numerator coefficients), made
// once so that the element can be evaluated at several roots modulo p, such
// as the degree-one primes over one rational prime: the reductions depend on
// p only.  After reduce(element, p), value_at_degree_one_root(out, root) is
// factored_value_at_degree_one_root_nmod(out, element, p, root): it fails
// exactly where that does and otherwise gives the same residue.  The
// reduction borrows `element`, which must outlive it unchanged.
class FactoredElementNmodReduction {
public:
    void reduce(const FactoredElement& element, ulong p) noexcept;
    bool value_at_degree_one_root(ulong& out, ulong root) const noexcept;

private:
    // A factor with nonzero exponent: its exponent, its index in
    // element.factors(), its stored denominator mod p and, when that is
    // nonzero, the `length` reduced numerator coefficients starting at
    // `offset` in coefficients_, lowest degree first.
    struct Factor {
        slong exponent = 0;
        std::size_t index = 0;
        std::size_t offset = 0;
        std::size_t length = 0;
        ulong denominator = 0;
    };

    const FactoredElement* element_ = nullptr;
    ulong p_ = 0;
    ulong p_inverse_ = 0;
    nmod_t mod_ = {};
    bool defined_ = false;
    std::vector<Factor> factors_;
    std::vector<ulong> coefficients_;
};
bool element_value_at_degree_one_root(flint::Fmpz& out,
                                      const Element& element,
                                      flint::FmpzConstRef p,
                                      flint::FmpzConstRef root) noexcept;

}  // namespace silex::detail
