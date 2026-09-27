#pragma once

#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>

namespace silex {
class Element;
class NumberField;
}  // namespace silex

namespace silex::detail {

// Search step of roots_of_unity (see src/unit/roots_of_unity.cpp), exposed
// for tests.  `bound` must be a multiple of the number w of roots of unity in
// `field`; the function enforces this by failing unless `bound` is a multiple
// of the reduced good-prime bound, itself a proven multiple of w.  Succeeds
// only when it finds a root of unity of exact order `bound`, which then
// proves w = bound; writes order = bound and that generator.  Fails closed,
// leaving the outputs unchanged, when it cannot attain the bound, including
// when the bound is a strict multiple of w, and when the precondition fails.
bool roots_of_unity_for_bound(flint::FmpzRef order,
                              Element& generator,
                              const NumberField& field,
                              ulong bound) noexcept;

// Rigorous lower bounds for the regulator R of a number field with r1 real
// and r2 complex places and w roots of unity.  Each function writes an Arb
// enclosure of its bound, or returns false when the bound does not apply or
// its evaluation is not finite.  unit_lower_regulator_bound takes the maximum
// of the lower endpoints; see docs/reference/algorithms_and_sources.rst.

// Friedman 1989, Theorem B (p. 599): R >= 0.2052 for every number field.
bool regulator_lower_bound_friedman_minimum(flint::ArbRef out,
                                            slong precision) noexcept;

// Zimmert 1981, Korollar (i) (p. 375): 2R/w >= 0.04 exp(0.46 r1 + 0.1 r2),
// that is R >= 0.02 w exp(0.46 r1 + 0.1 r2).
bool regulator_lower_bound_zimmert_corollary(flint::ArbRef out,
                                             slong r1,
                                             slong r2,
                                             flint::FmpzConstRef w,
                                             slong precision) noexcept;

// Friedman 1989, Corollary on p. 620: R/w > 0.0031 exp(0.241 n + 0.497 r1)
// with n = r1 + 2 r2.
bool regulator_lower_bound_friedman_corollary(flint::ArbRef out,
                                              slong r1,
                                              slong r2,
                                              flint::FmpzConstRef w,
                                              slong precision) noexcept;

// Zimmert 1981, Satz 3 (p. 374), at the exact rational
// gamma = gamma_numerator / gamma_denominator > 0; the result bounds R, that
// is w times the Satz 3 bound for R/w.
bool regulator_lower_bound_zimmert_satz3(flint::ArbRef out,
                                         slong r1,
                                         slong r2,
                                         flint::FmpzConstRef w,
                                         slong gamma_numerator,
                                         slong gamma_denominator,
                                         slong precision) noexcept;

// Friedman 1989, Table 6 (p. 621), for the signatures of positive unit rank
// that the table lists; false for any other signature.
bool regulator_lower_bound_friedman_table6(flint::ArbRef out,
                                           slong r1,
                                           slong r2,
                                           slong precision) noexcept;

// The maximum of the lower endpoints of all the bounds above, with Satz 3
// evaluated at a fixed set of rational gamma values; an exact positive
// point.
bool regulator_lower_bound_from_signature(flint::ArbRef out,
                                          slong r1,
                                          slong r2,
                                          flint::FmpzConstRef w,
                                          slong precision) noexcept;

}  // namespace silex::detail
