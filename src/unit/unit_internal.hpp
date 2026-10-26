#pragma once

#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>

#include <flint/fmpz_poly.h>

#include <vector>

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
// otherwise, including for translates Phi_n(+-x + c) that are not themselves a
// cyclotomic polynomial, scalar multiples and undefined fields.  It reads the
// raw defining polynomial, not the monic integral model that the root bound
// uses, so it returns 0 for 3*Phi_7 where the root bound still gets 14.  Candidates are the n with phi(n) equal to
// the degree, compared with FLINT's fmpz_poly_cyclotomic.
ulong cyclotomic_defining_polynomial_index(const NumberField& field) noexcept;

// The index n >= 3 when `polynomial` is exactly psi_n, the minimal
// polynomial of 2 cos(2 pi / n) over Q (FLINT fmpz_poly_cos_minpoly), which
// defines the maximal real subfield Q(zeta_n + zeta_n^-1) of Q(zeta_n) and
// has degree phi(n) / 2; 0 otherwise, including for translates, scalar
// multiples and null or constant input.  The candidate n is read from
// x^d psi(x + 1/x), d = deg(psi), which is Phi_n exactly when psi = psi_n
// (its roots are then the zeta^(+-k), k coprime to n), and confirmed by an
// exact comparison with fmpz_poly_cos_minpoly(n).  Degree one is never
// matched (psi_3 = x + 1 and psi_6 = x - 1 are also Phi_2 and Phi_1).
ulong real_cyclotomic_polynomial_index(
        const fmpz_poly_struct* polynomial) noexcept;

// real_cyclotomic_polynomial_index of the raw defining polynomial of
// `field` when it is integral, and 0 otherwise, read like
// cyclotomic_defining_polynomial_index (no monic integral model is taken).
ulong real_cyclotomic_defining_polynomial_index(
        const NumberField& field) noexcept;

// The defining polynomial f of a field when it is recognized exactly as Phi_n
// (cyclotomic_index = n, from cyclotomic_defining_polynomial_index) or else
// as psi_n (real_cyclotomic_index = n, from
// real_cyclotomic_defining_polynomial_index); both are 0 otherwise.
struct CyclotomicDefiningPolynomial {
    ulong cyclotomic_index = 0;
    ulong real_cyclotomic_index = 0;
};

CyclotomicDefiningPolynomial cyclotomic_defining_polynomial(
        const NumberField& field) noexcept;

// Closed form of the test "f modulo the prime p is a product of distinct
// linear factors" for a recognized f = Phi_n or psi_n and p not dividing n:
// writes the answer to `splits` and returns true.  Returns false, leaving
// `splits` unchanged, when f is not recognized or p divides n; the caller
// then tests the reduction of f itself.  `prime` must be prime.  For p not
// dividing n, Phi_n splits so exactly when p = 1 mod n and psi_n exactly
// when p = +-1 mod n (see roots_of_unity.cpp for the sources).
bool cyclotomic_defining_polynomial_splits(
        bool& splits,
        const CyclotomicDefiningPolynomial& polynomial,
        ulong prime) noexcept;

// The cyclotomic units of source nfcyclotomicunits (buch2.c) when the
// defining polynomial of `field` is exactly Phi_m
// (cyclotomic_defining_polynomial_index), and none otherwise, including
// for fields that merely contain a root of unity.  With w = lcm(2, m) and
// z a primitive n-th root of unity, n = w or w / 2 so that n != 2 (mod 4)
// (z = x, or z = -x when m = 2 (mod 4)): nothing when w <= 6; otherwise
// (Z^a - 1) / (Z - 1) for each p^k || n, Z = z^(n / p^k),
// 1 < a <= (p^k - 1) / 2, p not dividing a; and, when n is not a prime
// power, z^a - 1 for 1 <= a < n / 2 with a = 1 or a not dividing n and
// n / gcd(a, n) not a prime power.  When n is a prime power the source
// states the list is independent, and it has phi(n) / 2 - 1 elements, the
// unit rank; otherwise it may be dependent or longer than the unit rank.  Every element is
// checked to have norm +-1 (each is integral by construction).  Returns
// false, with `out` empty, when an operation or a norm check fails.
bool cyclotomic_field_units(std::vector<Element>& out,
                            const NumberField& field) noexcept;

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
