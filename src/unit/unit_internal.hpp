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

// The index n when the defining polynomial of `field` is exactly the
// cyclotomic polynomial Phi_n (integral, with Phi_n's coefficients), and 0
// otherwise, including for translates Phi_n(+-x + c) with c != 0, scalar
// multiples and undefined fields.  Candidates are the n with phi(n) equal to
// the degree, compared with FLINT's fmpz_poly_cyclotomic.
ulong cyclotomic_defining_polynomial_index(const NumberField& field) noexcept;

// The proven multiple of w that roots_of_unity searches against, for the
// monic integral model of the defining polynomial of `field`, exposed for
// tests.  With use_cyclotomic_shortcut, a cyclotomic model polynomial Phi_n
// gives w = lcm(2, n) directly; without it the reduced good-prime bound is
// always computed.  False when the bound is not computed.
bool proven_root_of_unity_bound(ulong& bound,
                                const NumberField& field,
                                bool use_cyclotomic_shortcut) noexcept;

// Per-thread counts of the roots-of-unity memo in
// src/unit/roots_of_unity.cpp, which stores w and a generator per defining
// polynomial and backend.  Diagnostic only: the memo never changes a
// computed value.
struct RootsOfUnityCacheCounters {
    slong hits = 0;
    slong misses = 0;
};
RootsOfUnityCacheCounters roots_of_unity_cache_counters() noexcept;

// Empties the roots-of-unity memo on the calling thread, so a test can
// compare a cold computation with a memoized one.  Diagnostic only.
void clear_roots_of_unity_cache() noexcept;

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

// Pohst-Zassenhaus 1989, (6.22) (p. 366), which excerpts it from Pohst 1977
// and 1978 (Remak's method, Remak 1932): for a field K of degree n = r1 + 2 r2 with no proper subfield other than Q,
// unit rank r = r1 + r2 - 1 and |d_K| > c,
//   R >= ((3 log^2(|d_K| / c) / ((n-1) n (n+1) - 6 r2))^r
//         * 2^r2 / (n gamma_r^r))^(1/2),
// with c = 4^floor(n/2) when r2 = 0 and n <= 11 and c = n^n otherwise, and
// gamma_r^r Hermite's constant to the power r: the exact values of
// Pohst-Zassenhaus (3.35a) (p. 198) for r <= 8, else Blichfeldt's upper
// bound (2/pi)^r Gamma(2 + r/2)^2, (3.35b) (p. 199).  Silex certifies that
// K has no proper subfield other than Q only by n being prime, so the
// function returns false unless n is prime, r >= 1 and |d_K| > c.
// `field_discriminant` must be the discriminant of K (of its maximal order).
bool regulator_lower_bound_remak_pohst(flint::ArbRef out,
                                       slong r1,
                                       slong r2,
                                       flint::FmpzConstRef field_discriminant,
                                       slong precision) noexcept;

// unit_lower_regulator_bound, raised to the lower endpoint of
// regulator_lower_bound_remak_pohst when `field_discriminant` is non-null
// (it must then be the discriminant of `field`) and that bound applies.
// With a null discriminant this is exactly unit_lower_regulator_bound.
bool field_lower_regulator_bound(flint::ArbRef out,
                                 const NumberField& field,
                                 const flint::Fmpz* field_discriminant,
                                 slong precision) noexcept;

}  // namespace silex::detail
