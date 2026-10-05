#pragma once

#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>

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

// Encloses (G_{m,j}(gamma) + n log pi) / 2 for gamma = gamma_hundredths/100,
// where G_{m,j} is the left side of Zimmert's Hilfssatz (Zimmert 1981,
// p. 378) for the function H of the proof of Satz 4 (p. 379) with a fixed
// number m of factors and j = ceil(1/2 + gamma).  Fails for an integer gamma,
// gamma outside (0, 3), or r1 + r2 < 1.
bool zimmert_twin_class_log_lower_bound(flint::ArbRef out,
                                        slong real_places,
                                        slong complex_pairs,
                                        ulong gamma_hundredths) noexcept;

// Twin-class generation bound for 3 <= n <= 20: the smaller of
// generation_bound and the ceiling of a rigorous upper bound for
// sqrt|d| exp(-(G_{m,j} + n log pi)/2) at the tabulated gamma (Zimmert 1981,
// Korollar 1, in the finite-m form of its proof).  For every ideal class K,
// K or D K^(-1) (D the class of the different) contains an integral ideal of
// norm at most out.  So the prime ideals of norm at most out generate the
// class group only if the class of the different lies in the subgroup they
// generate; callers must establish that separately.  It is deliberately not
// part of generation_bound or factor_base_class_group_bound, whose contract
// is that every class contains a small ideal.  Fails, leaving out
// unchanged, outside 3 <= n <= 20 or on the inputs generation_bound
// rejects.
bool twin_class_generation_bound(flint::FmpzRef out,
                                 flint::FmpzConstRef abs_discriminant,
                                 slong real_places,
                                 slong complex_pairs) noexcept;

// How a proven generation requirement handles the class of the different.
// - none: the requirement is factor_base_class_group_bound (every class
//   contains a small ideal), so the different plays no part.
// - index_one: the defining polynomial f is monic integral and
//   O_K = Z[theta], so the different is f'(theta) O_K, a principal ideal.
// - conductor_primes: f is monic integral and [O_K : Z[theta]] > 1; the
//   requirement also covers every prime ideal above each p dividing the
//   index, which includes every prime dividing the conductor F of Z[theta],
//   and f'(theta) O_K = F D makes [D] = [F]^(-1).
enum class GenerationDifferentRoute {
    none = 0,
    index_one = 1,
    conductor_primes = 2
};

// The bound to which a proven request must verify factor-base generation
// for the maximal order `order`, and how it treats the different.  It is
// factor_base_class_group_bound unless Zimmert's twin-class bound (Zimmert
// 1981, Korollar 1, finite-m form; 3 <= n <= 20) gives a strictly smaller
// requirement with the different settled exactly:
// - index_one: the requirement is the twin-class bound B;
// - conductor_primes: it is max(B, largest norm of a prime ideal above any
//   p | [O_K : Z[theta]]), so every prime ideal up to it generating the
//   class group also puts the conductor primes in the generated subgroup.
// Every prime ideal of norm at most `bound` then generates the class group.
// Fails, leaving `bound` and `route` unchanged, when
// factor_base_class_group_bound fails.  Never part of
// factor_base_class_group_bound.
bool proven_generation_requirement(flint::FmpzRef bound,
                                   GenerationDifferentRoute& route,
                                   const Order& order) noexcept;

}  // namespace silex::detail
