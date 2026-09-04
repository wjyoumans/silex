#pragma once

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
bool factored_value_at_degree_one_root_nmod(ulong& out,
                                           const FactoredElement& element,
                                           ulong p,
                                           ulong root) noexcept;
bool element_value_at_degree_one_root(flint::Fmpz& out,
                                      const Element& element,
                                      flint::FmpzConstRef p,
                                      flint::FmpzConstRef root) noexcept;

}  // namespace silex::detail
