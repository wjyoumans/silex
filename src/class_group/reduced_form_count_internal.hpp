#pragma once

#include <flint/flint.h>

namespace silex::detail {

// True when D = -abs_discriminant is a fundamental discriminant: either
// D = 1 mod 4 and |D| is squarefree, or D = 4m with m = 2 or 3 mod 4 and m
// squarefree (Cohen, GTM 138, Definition 5.1.2).
bool imaginary_quadratic_discriminant_is_fundamental(
        ulong abs_discriminant) noexcept;

// The number of reduced forms of a fundamental discriminant D < 0, counted
// without building them; this is h(D) (Cohen, GTM 138, Proposition 5.3.3).
// Every reduced form (a, b, c) has a <= sqrt(|D| / 3) (Lemma 5.3.4(1)), and
// for fundamental D every form is primitive (pp. 225 and 232).  For each
// such a the forms are the roots b of b^2 = D mod 4a, with b defined mod 2a
// (Cohen, p. 253) and c = (b^2 - D) / (4a):
//
// - when 4a^2 < |D|, every b with -a < b <= a gives a reduced form with
//   c > a (Lemma 5.3.4(2)), so a contributes R(4a) / 2, where R(m) is the
//   number of roots of x^2 = D mod m.  R(4a) is the product of R(p^e) over
//   the prime powers p^e exactly dividing 4a (the Chinese remainder theorem,
//   Cohen, Theorem 1.3.9 and Corollary 1.3.10), and each R(p^e) is the root
//   count of FLINT `n_sqrtmod_primepow`;
// - when 4a^2 >= |D|, the roots mod 4a are the Chinese-remainder
//   combinations of the roots mod the prime powers exactly dividing 4a, kept
//   from the `n_sqrtmod_primepow` calls above (4a is factored by a sieve over
//   the band), and are kept when (a, b, c) is reduced: |b| <= a <= c, and
//   b >= 0 when |b| = a or a = c (Cohen, Definition 5.3.2).
//
// Sets `out` and returns true when D = -abs_discriminant is fundamental and
// |D| < 2^62.  Returns false, leaving `out` unspecified, otherwise, or when
// an internal consistency check fails: the band factorisation does not
// rebuild a, the number of band roots disagrees with the product R(4a), or
// FLINT returns a root set that a fundamental D does not allow.  On a
// non-fundamental D the sum would also count imprimitive forms, so it is
// refused.  Memory is about 5 bytes per a <= sqrt(|D| / 3); the caller
// bounds |D|.
bool imaginary_quadratic_fundamental_reduced_form_count(
        ulong& out,
        ulong abs_discriminant) noexcept;

}  // namespace silex::detail
