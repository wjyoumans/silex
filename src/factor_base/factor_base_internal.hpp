#pragma once

#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>
#include <silex/prime_ideal.hpp>

#include <vector>

#include "generation_different_route_internal.hpp"
#include "../zeta/zeta_internal.hpp"

namespace silex::detail {

// Decompositions of the rational primes 2, 3, 5, ..., in that order, kept
// across norm-bounded builds of one maximal order: entry k is
// decompose_prime of the k-th rational prime, or, when
// skipped_smallest_residue_degree[k] > 0, an empty placeholder for a prime
// that a build skipped without decomposing it.  A prime is skipped only
// when residue_degrees (the defining-polynomial route, under its own gate)
// gives its residue degrees and p^f exceeds the bound for the smallest of
// them, f, so that no prime ideal above p has norm at most the bound.  The
// placeholder keeps that f: a later build at a bound of at least p^f
// decomposes p then and fills the entry.  `order` is a copy of the parent
// order of the first build, for which residue_degrees is prepared; a later
// build over an order with other data fails.
struct PrimeDecompositionCache {
    std::vector<PrimeIdealList> decompositions;
    std::vector<slong> skipped_smallest_residue_degree;
    DefiningPolynomialResidueDegrees residue_degrees;
    Order order;
};

class FactorBaseBlockAccess {
public:
    // Builds exactly the base base.build_prime_ideal_norm_bounded(bound)
    // builds, taking the decomposition of each rational prime from `cache`
    // when the cache already holds it and appending the decompositions it
    // computes, or a placeholder for a prime it skips (see
    // PrimeDecompositionCache).  The cache must come from builds over base's parent order
    // (a cached decomposition of another order fails the build).  On
    // failure base is unchanged; the cache keeps every decomposition it
    // holds and may have gained some.  The order check runs before the
    // prime loop: a cache that already names another order fails even when
    // no cached decomposition would be read.  Once the base, bound and
    // parent-order checks pass, the first build over an empty cache stores
    // the parent order before any prime is decomposed, so a later failure
    // still leaves the cache bound to that order.  Use one cache per order
    // (the class-unit transaction does).  `diagnostics` only feeds profile
    // scopes; it does not change the build.
    static bool build_prime_ideal_norm_bounded_cached(
            FactorBase& base,
            flint::FmpzConstRef bound,
            PrimeDecompositionCache& cache,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;
    // The body of FactorBase::build_prime_ideal_norm_bounded_impl, with
    // the optional decomposition cache of the overload above.
    static bool build_prime_ideal_norm_bounded(
            FactorBase& base,
            flint::FmpzConstRef bound,
            bool include_inert_primes,
            slong degree_limit,
            bool retain_one_split_prime,
            PrimeDecompositionCache* cache,
            const DiagnosticsContext* diagnostics) noexcept;
    static bool rational_prime_block_is_complete(
            bool& complete,
            const FactorBase& base,
            slong block_index) noexcept;
    // The ordering step of the LLL relation base: sorts the primes by norm,
    // largest first, ties kept in their current order, and rebuilds the
    // rational-prime blocks.  Applied to the result of
    // build_prime_ideal_norm_bounded(bound) it gives exactly the base
    // build_lll_relation_base_impl(bound, true) builds.  Fails, possibly
    // leaving base unchanged, if a norm cannot be computed or a copy fails.
    static bool order_lll_relation_base_by_norm(FactorBase& base) noexcept;
};

// Proven generation bound used by factor_base_class_group_bound for degree
// n = real_places + 2 * complex_pairs >= 3: every ideal class of a field
// with this absolute discriminant and signature contains an integral ideal
// of norm at most out.  It is the smallest of Zimmert's bound (Satz 2),
// Minkowski's bound with (4/pi)^r2, and the exact integer Minkowski form
// with 2^r2.  Fails, leaving out unchanged, if n < 3, the signature is
// negative, or abs_discriminant is not positive.
bool generation_bound(flint::FmpzRef out,
                      flint::FmpzConstRef abs_discriminant,
                      slong real_places,
                      slong complex_pairs) noexcept;

// Zimmert's bound alone (Zimmert 1981, Satz 2 with alpha from Bemerkung 1):
// the ceiling of a rigorous upper bound for sqrt|d| exp(-L(r1, r2, gamma))
// with gamma = zimmert_satz2_gamma_hundredths(r1, r2) / 100.  Fails,
// leaving out unchanged, on the inputs generation_bound rejects.
bool zimmert_generation_bound(flint::FmpzRef out,
                              flint::FmpzConstRef abs_discriminant,
                              slong real_places,
                              slong complex_pairs) noexcept;

// Gamma (in hundredths) used for Zimmert's Satz 2 bound of this signature:
// a hand-checked table for n <= 20, a table generated by
// tools/zimmert_gamma_table.cpp for 21 <= n <= 64, and
// zimmert_satz2_gamma_search above.  0 for n < 3 or a negative signature.
ulong zimmert_satz2_gamma_hundredths(slong real_places,
                                     slong complex_pairs) noexcept;

// The deterministic gamma search behind the tables: the k in [1, 300]
// maximising the midpoint of zimmert_satz2_log_lower_bound(r1, r2, k), the
// smallest such k on ties.  0 for n < 3 or a negative signature.
ulong zimmert_satz2_gamma_search(slong real_places,
                                 slong complex_pairs) noexcept;

// Encloses the right side L of Zimmert 1981, Satz 2 (p. 372) at
// gamma = gamma_hundredths/100, alpha from Bemerkung 1 (p. 373), with Arb at
// 128 bits.  Fails, leaving out unchanged, for n < 3, a negative signature,
// gamma_hundredths = 0, or an enclosure that is not finite.
bool zimmert_satz2_log_lower_bound(flint::ArbRef out,
                                   slong real_places,
                                   slong complex_pairs,
                                   ulong gamma_hundredths) noexcept;

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
// and the evaluation procedure of that proof file.  With terms <= 0 on
// entry it chooses L and K as that proof file describes, each tail bounded
// by 2^-128, failing when the proof's search finds no admissible pair
// meeting that bound.  The search starts at the proof's lower estimate for
// L (failing at once if that estimate is not finite or exceeds 2^16,
// kZimmertMaxTailTerms), grows L by about a quarter each step up to 2^16,
// tries K up to 64 (kZimmertMaxTailOrder) and stops at the first K that is
// not admissible for a given L; it is not exhaustive, so an admissible
// pair meeting the bound may exist and not be found.  With terms > 0 it uses L = terms and K = order (the ball is still
// valid but may be wider), failing unless 2 <= K <= 64, L <= 2^16,
// L > gamma/2, and both (L, 0, K) and (L, 1/2, K) are admissible (the
// second follows from the first, since admissibility only needs
// N = 2(L + 1 + delta) large enough).  On success sets terms and order to
// the L and K used.  In either case fails for a < 1, b < 0, a numerator or
// denominator of gamma that is <= 0 or exceeds 2^16, or when an internal
// check fails: a ball that is not finite (or, for the M_Phi denominator,
// not positive), or the Phi_0 = 0, Phi_1 = -c^2/4 check of the proof's
// Statement.
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
