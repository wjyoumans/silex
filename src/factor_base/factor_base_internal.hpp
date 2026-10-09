#pragma once

#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>

#include "generation_different_route_internal.hpp"

namespace silex::detail {

class FactorBaseBlockAccess {
public:
    static bool rational_prime_block_is_complete(
            bool& complete,
            const FactorBase& base,
            slong block_index) noexcept;
};

// Proven generation bound used by factor_base_class_group_bound for degree
// n = real_places + 2 * complex_pairs >= 3: every ideal class of a field
// with this absolute discriminant and signature contains an integral ideal
// of norm at most out.  It is the smallest of Zimmert's bound (n <= 20),
// Minkowski's bound with (4/pi)^r2, and the exact integer Minkowski form
// with 2^r2.  Fails, leaving out unchanged, if n < 3, the signature is
// negative, or abs_discriminant is not positive.
bool generation_bound(flint::FmpzRef out,
                      flint::FmpzConstRef abs_discriminant,
                      slong real_places,
                      slong complex_pairs) noexcept;

// Zimmert's bound alone (Zimmert 1981, Satz 2 with alpha from Bemerkung 1):
// the ceiling of a rigorous upper bound for sqrt|d| exp(-L(r1, r2, gamma))
// with the tabulated gamma.  Fails, leaving out unchanged, outside
// 3 <= n <= 20 or on the inputs generation_bound rejects.
bool zimmert_generation_bound(flint::FmpzRef out,
                              flint::FmpzConstRef abs_discriminant,
                              slong real_places,
                              slong complex_pairs) noexcept;

// Encloses (F_{a,b}(gamma) + n log pi) / 2 for gamma = gamma_hundredths/100,
// a = r1 + r2, b = r2 and n = r1 + 2 r2, where F_{a,b} is Zimmert's series
// (Zimmert 1981, p. 376) enclosed as by zimmert_f_enclosure.  Sets *terms
// and *order, when not null, to the L and K used.  Fails for gamma outside
// (0, 3) or r1 + r2 < 1.
bool zimmert_twin_class_log_lower_bound(flint::ArbRef out,
                                        slong real_places,
                                        slong complex_pairs,
                                        ulong gamma_hundredths,
                                        slong* terms = nullptr,
                                        slong* order = nullptr) noexcept;

// Certified enclosure of Zimmert's F_{a,b}(gamma), gamma =
// gamma_numerator/gamma_denominator > 0 (integers allowed), by part (iv) of
// GB-zimmert-F-tail-bound (silex-research,
// research/proofs/generation/GB-zimmert-F-tail-bound.md at commit 06191ff)
// and the evaluation procedure of that proof file, each tail bounded by
// 2^-128.  With terms = 0 on entry it chooses L and K as that proof file
// describes; with terms > 0 it uses L = terms and K = order, failing unless
// L > gamma/2, K >= 2 and (L, 0, K) is admissible.  On success sets terms
// and order to the L and K used.  Fails for a < 1, b < 0 or gamma <= 0.
bool zimmert_f_enclosure(flint::ArbRef out,
                         slong& terms,
                         slong& order,
                         slong a,
                         slong b,
                         slong gamma_numerator,
                         slong gamma_denominator) noexcept;

// Gamma (in hundredths) used for the twin-class bound of this signature,
// or 0 if the signature has none (n < 3, n > 20, or a negative signature).
ulong zimmert_twin_class_gamma_hundredths(slong real_places,
                                          slong complex_pairs) noexcept;

// The precomputed twin-class table entry of this signature: the exact
// decimal rational E = significand / 10^scale, stated to be at least
// exp(-c) for c the lower endpoint of zimmert_twin_class_log_lower_bound at
// gamma_hundredths.  The table is generated offline by
// tools/zimmert_twin_class_table.cpp and verified by
// test/t-generation-bound.cpp; nothing here re-checks it.  Fails, leaving
// the outputs unchanged, for a signature without an entry.
bool zimmert_twin_class_table_entry(flint::FmpzRef significand,
                                    slong& scale,
                                    ulong& gamma_hundredths,
                                    slong real_places,
                                    slong complex_pairs) noexcept;

// Twin-class generation bound for 3 <= n <= 20: the smaller of
// generation_bound and ceil(E sqrt|d|), E the precomputed table entry above,
// an exact upper bound for exp(-(F_{r1+r2,r2} + n log pi)/2) at the
// tabulated gamma (Zimmert 1981, Korollar 1).  For
// every ideal class K, K or D K^(-1) (D the class of the different) contains
// an integral ideal of norm at most out, and the prime ideals of norm at most
// out generate the class group with no hypothesis on D (GB-kor1-generation;
// see proven_generation_requirement).  It does not promise a small ideal in
// every class, so it is deliberately not part of generation_bound or
// factor_base_class_group_bound, whose contract is that every class
// contains a small ideal.  Fails, leaving out
// unchanged, outside 3 <= n <= 20 or on the inputs generation_bound
// rejects.
bool twin_class_generation_bound(flint::FmpzRef out,
                                 flint::FmpzConstRef abs_discriminant,
                                 slong real_places,
                                 slong complex_pairs) noexcept;

// The bound to which a proven request must verify factor-base generation
// for the maximal order `order`, and how it treats the different.  It is
// factor_base_class_group_bound unless Zimmert's twin-class bound B (Zimmert
// 1981, Korollar 1; 3 <= n <= 20) is strictly smaller, in which case it is B
// with route kor1_generation, for any defining polynomial: by
// GB-kor1-generation (silex-research, commit 8926120) B is a generation
// bound with no hypothesis on the class of the different.  Every prime
// ideal of norm at most `bound` then generates the class group.
// Fails, leaving `bound` and `route` unchanged, when
// factor_base_class_group_bound fails.  Never part of
// factor_base_class_group_bound.
bool proven_generation_requirement(flint::FmpzRef bound,
                                   GenerationDifferentRoute& route,
                                   const Order& order) noexcept;

}  // namespace silex::detail
