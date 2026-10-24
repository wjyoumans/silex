#include "class_group/relation_candidate_internal.hpp"

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/flint/nmod_poly_factor.hpp>

#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace {
namespace sflint = silex::flint;
namespace relation_search = silex::detail::relation_search;

ulong next_random(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) + UWORD(1442695040888963407);
    return state >> 11;
}

// Cyclotomic Phi_m for prime m: 1 + x + ... + x^(m - 1).
void set_prime_cyclotomic(sflint::FmpzPoly& out, slong m) noexcept {
    fmpz_poly_zero(out.raw());
    for (slong i = 0; i < m; ++i) {
        fmpz_poly_set_coeff_si(out.raw(), i, 1);
    }
}

// Degree-17 subfield of Q(zeta_103): totally real and abelian, so primes
// splitting completely have density 1/17 and norms take both signs.
void set_real_abelian_degree_17(sflint::FmpzPoly& out) noexcept {
    constexpr slong coefficients[] = {
            2837,   25763,  75722, 58507, -102976, -198615,
            -50288, 107759, 74838, -11031, -23311, -3653,
            2579,   763,    -105,  -48,   1,       1};
    fmpz_poly_zero(out.raw());
    slong i = 0;
    for (const slong coefficient : coefficients) {
        fmpz_poly_set_coeff_si(out.raw(), i, coefficient);
        ++i;
    }
}

void set_random_polynomial(sflint::FmpzPoly& out,
                           ulong& state,
                           slong length,
                           flint_bitcnt_t bits) noexcept {
    fmpz_poly_zero(out.raw());
    sflint::Fmpz coefficient;
    for (slong i = 0; i < length; ++i) {
        fmpz_zero(coefficient.raw());
        for (flint_bitcnt_t done = 0; done < bits; done += 32) {
            fmpz_mul_2exp(coefficient.raw(), coefficient.raw(), 32);
            fmpz_add_ui(coefficient.raw(), coefficient.raw(),
                        next_random(state) & UWORD(0xffffffff));
        }
        if (bits > 0) {
            fmpz_fdiv_r_2exp(coefficient.raw(), coefficient.raw(), bits);
        }
        if ((next_random(state) & 1) != 0) {
            fmpz_neg(coefficient.raw(), coefficient.raw());
        }
        fmpz_poly_set_coeff_fmpz(out.raw(), i, coefficient.raw());
    }
    // Keep the requested length so the bound sees the full degree.
    if (length > 0 && fmpz_poly_length(out.raw()) < length) {
        fmpz_poly_set_coeff_si(out.raw(), length - 1, 1);
    }
}

// Compares the split-prime value with fmpz_poly_resultant over random
// coordinate polynomials; returns the number of handled candidates or -1
// on a mismatch.  This catches sign and CRT errors and grossly short
// bounds, but random candidates lie many bits below the Hadamard bound, so
// a bound short by a few bits would still pass: the bound's correctness
// rests on the Hadamard argument in the source, not on these tests.
slong compare_with_resultant(relation_search::SplitPrimeResultant& evaluator,
                             const sflint::FmpzPoly& field,
                             ulong seed,
                             slong count,
                             flint_bitcnt_t max_bits) noexcept {
    const slong degree = fmpz_poly_degree(field.raw());
    ulong state = seed;
    sflint::FmpzPoly coordinate;
    sflint::Fmpz expected;
    sflint::Fmpz actual;
    slong handled_count = 0;
    for (slong k = 0; k < count; ++k) {
        const slong length =
                1 + static_cast<slong>(next_random(state) %
                                       static_cast<ulong>(degree));
        const flint_bitcnt_t bits = static_cast<flint_bitcnt_t>(
                next_random(state) % (max_bits + 1));
        set_random_polynomial(coordinate, state, length, bits);
        fmpz_poly_resultant(expected.raw(), field.raw(), coordinate.raw());
        bool handled = false;
        evaluator.resultant(handled, actual.raw(), coordinate.raw());
        if (handled) {
            ++handled_count;
            if (!fmpz_equal(expected.raw(), actual.raw())) {
                std::fprintf(stderr,
                             "split-prime resultant mismatch: degree %ld, "
                             "candidate %ld\n",
                             static_cast<long>(degree), static_cast<long>(k));
                return -1;
            }
        }
    }
    return handled_count;
}

int test_cyclotomic_matches_resultant() {
    // Phi_17 (degree 16, abelian) and Phi_41 (degree 40): every candidate
    // with small coefficients is covered and equals the exact resultant.
    for (const slong m : {17, 41}) {
        sflint::FmpzPoly field;
        set_prime_cyclotomic(field, m);
        relation_search::SplitPrimeResultant evaluator;
        evaluator.reset(field.raw());
        if (!evaluator.enabled) {
            return 1;
        }
        constexpr slong count = 200;
        if (compare_with_resultant(evaluator, field,
                                   static_cast<ulong>(m), count, 12) !=
            count) {
            return 1;
        }
        // Large coefficients: some bounds exceed the split primes, so the
        // search runs out and those candidates fall back; the ones that are
        // covered still match.
        constexpr slong large_count = 40;
        const slong large_handled = compare_with_resultant(
                evaluator, field, static_cast<ulong>(m) + 1, large_count,
                200);
        if (large_handled <= 0 || large_handled >= large_count ||
            !evaluator.search_exhausted) {
            return 1;
        }
    }
    return 0;
}

int test_real_abelian_matches_resultant() {
    sflint::FmpzPoly field;
    set_real_abelian_degree_17(field);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());
    if (!evaluator.enabled) {
        return 1;
    }
    constexpr slong count = 200;
    if (compare_with_resultant(evaluator, field, 103, count, 12) != count) {
        return 1;
    }

    // Res(f, -1) = (-1)^17 needs one prime and a negative representative.
    sflint::FmpzPoly coordinate;
    sflint::Fmpz actual;
    bool handled = false;
    fmpz_poly_set_si(coordinate.raw(), -1);
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (!handled || !fmpz_equal_si(actual.raw(), -1)) {
        return 1;
    }
    // Res(f, x) = -f(0) for odd degree.
    fmpz_poly_zero(coordinate.raw());
    fmpz_poly_set_coeff_si(coordinate.raw(), 1, 1);
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    return handled && fmpz_equal_si(actual.raw(), -2837) ? 0 : 1;
}

int test_constant_and_zero_polynomials() {
    sflint::FmpzPoly field;
    set_prime_cyclotomic(field, 17);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());

    sflint::FmpzPoly coordinate;
    sflint::Fmpz actual;
    bool handled = false;
    // Res(f, c) = c^16 for a constant c.
    fmpz_poly_set_si(coordinate.raw(), -3);
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    sflint::Fmpz expected;
    fmpz_set_ui(expected.raw(), 43046721);  // 3^16
    if (!handled || !fmpz_equal(actual.raw(), expected.raw())) {
        return 1;
    }
    // Res(f, x - 1) = (-1)^16 f(1) = 17.
    fmpz_poly_zero(coordinate.raw());
    fmpz_poly_set_coeff_si(coordinate.raw(), 0, -1);
    fmpz_poly_set_coeff_si(coordinate.raw(), 1, 1);
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (!handled || !fmpz_equal_si(actual.raw(), 17)) {
        return 1;
    }
    // The zero polynomial is left to fmpz_poly_resultant.
    fmpz_poly_zero(coordinate.raw());
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    return handled ? 1 : 0;
}

int test_disabled_and_exhausted_cases() {
    sflint::FmpzPoly coordinate;
    sflint::Fmpz actual;
    bool handled = true;
    fmpz_poly_set_coeff_si(coordinate.raw(), 0, 2);
    fmpz_poly_set_coeff_si(coordinate.raw(), 1, 1);

    // Below the degree cutoff the evaluator stays disabled.
    sflint::FmpzPoly small_field;
    set_prime_cyclotomic(small_field, 13);
    relation_search::SplitPrimeResultant small_evaluator;
    small_evaluator.reset(small_field.raw());
    small_evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (small_evaluator.enabled || handled) {
        return 1;
    }

    // x^16 - x - 1 has Galois group S_16, so no tested prime splits it
    // completely and the bounded search gives up.
    sflint::FmpzPoly field;
    fmpz_poly_set_coeff_si(field.raw(), 0, -1);
    fmpz_poly_set_coeff_si(field.raw(), 1, -1);
    fmpz_poly_set_coeff_si(field.raw(), 16, 1);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());
    handled = true;
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (!evaluator.enabled || handled || !evaluator.search_exhausted ||
        evaluator.tested_primes > 4 * 16) {
        return 1;
    }
    // Later candidates do not restart the search.
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    return handled || evaluator.tested_primes > 4 * 16 ? 1 : 0;
}

int test_reset_keeps_split_primes_for_same_field() {
    sflint::FmpzPoly field;
    set_prime_cyclotomic(field, 17);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());
    if (compare_with_resultant(evaluator, field, 5, 20, 40) != 20) {
        return 1;
    }
    const std::size_t found = evaluator.moduli.size();
    const slong tested = evaluator.tested_primes;
    if (found == 0) {
        return 1;
    }

    // Rebuilding the prefilter for the same f keeps the primes found and
    // does not search again for them.
    evaluator.reset(field.raw());
    if (!evaluator.enabled || evaluator.moduli.size() != found ||
        evaluator.tested_primes != tested) {
        return 1;
    }
    if (compare_with_resultant(evaluator, field, 5, 20, 40) != 20 ||
        evaluator.moduli.size() != found ||
        evaluator.tested_primes != tested) {
        return 1;
    }

    // A different f starts a fresh search.
    sflint::FmpzPoly other;
    set_prime_cyclotomic(other, 41);
    evaluator.reset(other.raw());
    if (!evaluator.enabled || !evaluator.moduli.empty() ||
        evaluator.tested_primes != 0 ||
        compare_with_resultant(evaluator, other, 7, 20, 12) != 20) {
        return 1;
    }
    // No field polynomial disables the evaluator.
    evaluator.reset(nullptr);
    return evaluator.enabled || !evaluator.moduli.empty() ? 1 : 0;
}


// Res(f, a) for constants and monomials c x^j has |Res| = |c|^n, where the
// AM-GM inequality |N(a)|^2 <= (T2(a) / n)^n is an equality, so these
// candidates sit within a fraction of a bit of the T2 bound.  Sweeping
// |c| across 2^k - 1, 2^k and 2^k + 1 places |Res| on both sides of every
// modulus size the bound can choose.  Counts the handled calls and, among
// them, the values c for which c x^(n-1) used no more split primes than
// the constant c: the two have the same T2 bound up to rounding, while
// the Hadamard bound of c x^(n-1) is larger by ||f||_2^(n-1).
struct TightCandidateCounts {
    int handled = 0;
    int pairs = 0;
    int pairs_without_extra_primes = 0;
};

int compare_tight_candidates(relation_search::SplitPrimeResultant& evaluator,
                             const sflint::FmpzPoly& field,
                             TightCandidateCounts& counts) noexcept {
    const slong degree = fmpz_poly_degree(field.raw());
    sflint::FmpzPoly coordinate;
    sflint::Fmpz c;
    sflint::Fmpz expected;
    sflint::Fmpz actual;
    counts = TightCandidateCounts{};
    for (slong k = 1; k <= 40; ++k) {
        for (const slong offset : {-1, 0, 1}) {
            fmpz_one(c.raw());
            fmpz_mul_2exp(c.raw(), c.raw(), static_cast<ulong>(k));
            fmpz_add_si(c.raw(), c.raw(), offset);
            slong constant_prime_count = 0;
            for (const slong power : {static_cast<slong>(0), degree / 2,
                                      degree - 1}) {
                for (const int sign : {1, -1}) {
                    fmpz_poly_zero(coordinate.raw());
                    fmpz_poly_set_coeff_fmpz(coordinate.raw(), power,
                                             c.raw());
                    if (sign < 0) {
                        fmpz_poly_neg(coordinate.raw(), coordinate.raw());
                    }
                    fmpz_poly_resultant(expected.raw(), field.raw(),
                                        coordinate.raw());
                    bool handled = false;
                    evaluator.resultant(handled, actual.raw(),
                                        coordinate.raw());
                    if (handled &&
                        !fmpz_equal(expected.raw(), actual.raw())) {
                        std::fprintf(stderr,
                                     "split-prime resultant mismatch on a "
                                     "tight candidate: degree %ld, k %ld\n",
                                     static_cast<long>(degree),
                                     static_cast<long>(k));
                        return 1;
                    }
                    if (!handled) {
                        continue;
                    }
                    ++counts.handled;
                    if (sign < 0) {
                        continue;
                    }
                    if (power == 0) {
                        constant_prime_count = evaluator.last_prime_count;
                    } else if (power == degree - 1 &&
                               constant_prime_count > 0) {
                        ++counts.pairs;
                        if (evaluator.last_prime_count <=
                            constant_prime_count) {
                            ++counts.pairs_without_extra_primes;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

// Each field makes 720 calls; at least half must be handled so the values
// are actually compared.  On Phi_17 and Phi_41 the T2 bound must be the
// one selected for most c x^(n-1).
int test_t2_bound_matches_resultant_on_tight_candidates() {
    constexpr int kMinimumHandled = 360;
    for (const slong m : {17, 41}) {
        sflint::FmpzPoly field;
        set_prime_cyclotomic(field, m);
        relation_search::SplitPrimeResultant evaluator;
        evaluator.reset(field.raw());
        TightCandidateCounts counts;
        if (compare_tight_candidates(evaluator, field, counts) != 0 ||
            counts.handled < kMinimumHandled || counts.pairs == 0 ||
            2 * counts.pairs_without_extra_primes < counts.pairs) {
            return 1;
        }
    }
    sflint::FmpzPoly field;
    set_real_abelian_degree_17(field);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());
    TightCandidateCounts counts;
    return compare_tight_candidates(evaluator, field, counts) != 0 ||
                           counts.handled < kMinimumHandled
                   ? 1
                   : 0;
}

int test_t2_bound_for_prime_cyclotomic() {
    // For Phi_41 in the power basis 1, ..., x^39,
    // T2(a) = 41 ||a||^2 - (a(1))^2, so the largest eigenvalue of the T2
    // Gram matrix is 41 and the certified ratio must lie in
    // [41/40, 41/40 (1 + 1/16)].
    sflint::FmpzPoly field;
    set_prime_cyclotomic(field, 41);
    relation_search::SplitPrimeResultant evaluator;
    evaluator.reset(field.raw());

    // Res(Phi_41, 1 + x) = Phi_41(-1) = 1.  The Hadamard bound
    // 41^(1/2) 2^20 needs two split primes; the T2 bound
    // (41/40 * 2)^20 < 2^21 needs one.
    sflint::FmpzPoly coordinate;
    sflint::Fmpz actual;
    fmpz_poly_set_coeff_si(coordinate.raw(), 0, 1);
    fmpz_poly_set_coeff_si(coordinate.raw(), 1, 1);
    bool handled = false;
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (!handled || !fmpz_is_one(actual.raw()) ||
        evaluator.last_prime_count != 1 ||
        evaluator.t2_bound_state !=
                relation_search::SplitPrimeResultant::T2BoundState::
                        available) {
        return 1;
    }

    arf_t lower;
    arf_t upper;
    arf_init(lower);
    arf_init(upper);
    arf_set_ui(lower, 41);
    arf_div_ui(lower, lower, 40, 64, ARF_RND_DOWN);
    arf_set_ui(upper, 41 * 17);
    arf_div_ui(upper, upper, 40 * 16, 64, ARF_RND_UP);
    const bool in_range =
            arf_cmp(evaluator.t2_ratio_bound.raw(), lower) >= 0 &&
            arf_cmp(evaluator.t2_ratio_bound.raw(), upper) <= 0;
    arf_clear(lower);
    arf_clear(upper);
    if (!in_range) {
        return 1;
    }

    // The minimal polynomial of 2 cos(2 pi / 41) has power sums
    // p_38 > 20 ||f||_2^2 (its roots lie near +-2), so the T2 bound cannot
    // beat the Hadamard bound and is not computed; values stay exact.
    // Resetting the Phi_41 evaluator to this field discards its bound.
    constexpr slong real_coefficients[] = {
            1,    -10,  -55, 165,  495,  -792, -1716,
            1716, 3003, -2002, -3003, 1365, 1820, -560,
            -680, 136,  153, -18,  -19,  1,    1};
    sflint::FmpzPoly real_field;
    slong i = 0;
    for (const slong coefficient : real_coefficients) {
        fmpz_poly_set_coeff_si(real_field.raw(), i, coefficient);
        ++i;
    }
    evaluator.reset(real_field.raw());
    if (evaluator.t2_bound_state !=
        relation_search::SplitPrimeResultant::T2BoundState::not_computed) {
        return 1;
    }
    return compare_with_resultant(evaluator, real_field, 11, 50, 12) == 50 &&
                           evaluator.t2_bound_state ==
                                   relation_search::SplitPrimeResultant::
                                           T2BoundState::unavailable
                   ? 0
                   : 1;
}

// The split-prime search without any screen: walks the primes from 2^27
// with the evaluator's stopping rule (at most 64 split primes; stop once
// 4 deg(f) (found + 1) primes have been tested) and keeps p when f mod p
// is a product of deg(f) distinct monic linear factors, read from the
// roots with multiplicities rather than from x^p mod f.
struct ReferenceSearch {
    std::vector<ulong> split_primes;
    // The roots modulo each split prime, sorted.
    std::vector<std::vector<ulong>> roots;
    slong tested_primes = 0;
    ulong last_tested_prime = 0;
};

void reference_split_prime_search(ReferenceSearch& out,
                                  const sflint::FmpzPoly& field) noexcept {
    const slong degree = fmpz_poly_degree(field.raw());
    out.split_primes.clear();
    out.roots.clear();
    out.tested_primes = 0;
    out.last_tested_prime = UWORD(1) << 27;
    for (;;) {
        const auto found = static_cast<slong>(out.split_primes.size());
        if (found >= 64 || out.tested_primes >= 4 * degree * (found + 1)) {
            return;
        }
        const ulong p = n_nextprime(out.last_tested_prime, 1);
        if (p >= (UWORD(1) << 28)) {
            return;
        }
        out.last_tested_prime = p;
        ++out.tested_primes;
        sflint::NmodPoly reduced(p);
        fmpz_poly_get_nmod_poly(reduced.raw(), field.raw());
        sflint::NmodPolyFactor roots;
        nmod_poly_roots(roots.raw(), reduced.raw(), 1);
        bool distinct = roots.raw()->num == degree;
        for (slong i = 0; distinct && i < roots.raw()->num; ++i) {
            distinct = roots.raw()->exp[i] == 1;
        }
        if (distinct) {
            out.split_primes.push_back(p);
            std::vector<ulong> found_roots;
            for (slong i = 0; i < roots.raw()->num; ++i) {
                const nmod_poly_struct* factor = roots.raw()->p + i;
                found_roots.push_back(
                        nmod_neg(nmod_poly_get_coeff_ui(factor, 0),
                                 factor->mod));
            }
            std::sort(found_roots.begin(), found_roots.end());
            out.roots.push_back(found_roots);
        }
    }
}

// Runs the evaluator's search to the end with one candidate whose bound
// needs more than the 64 split primes kept, then compares the primes
// found, their order, the number of primes tested and the set of roots
// modulo each prime (r_i = root_powers[i n + 1]) with the reference.
int compare_split_prime_search(
        relation_search::SplitPrimeResultant& evaluator,
        const sflint::FmpzPoly& field) noexcept {
    evaluator.reset(field.raw());
    if (!evaluator.enabled) {
        return 1;
    }
    sflint::FmpzPoly coordinate;
    ulong state = 709;
    set_random_polynomial(coordinate, state, fmpz_poly_degree(field.raw()),
                          4096);
    sflint::Fmpz actual;
    bool handled = true;
    evaluator.resultant(handled, actual.raw(), coordinate.raw());
    if (handled || !evaluator.search_exhausted) {
        return 1;
    }
    ReferenceSearch reference;
    reference_split_prime_search(reference, field);
    if (reference.split_primes.empty() ||
        evaluator.moduli.size() != reference.split_primes.size() ||
        evaluator.tested_primes != reference.tested_primes ||
        evaluator.last_tested_prime != reference.last_tested_prime) {
        std::fprintf(stderr,
                     "split-prime search differs: found %zu vs %zu, "
                     "tested %ld vs %ld\n",
                     evaluator.moduli.size(), reference.split_primes.size(),
                     static_cast<long>(evaluator.tested_primes),
                     static_cast<long>(reference.tested_primes));
        return 1;
    }
    const auto n = static_cast<std::size_t>(evaluator.degree);
    if (evaluator.root_powers.size() != reference.split_primes.size() * n * n) {
        return 1;
    }
    for (std::size_t i = 0; i < reference.split_primes.size(); ++i) {
        if (evaluator.moduli[i].n != reference.split_primes[i]) {
            return 1;
        }
        std::vector<ulong> roots(n);
        for (std::size_t j = 0; j < n; ++j) {
            roots[j] = evaluator.root_powers[(i * n + j) * n + 1];
        }
        std::sort(roots.begin(), roots.end());
        if (roots != reference.roots[i]) {
            std::fprintf(stderr, "split-prime roots differ at prime %lu\n",
                         static_cast<unsigned long>(reference.split_primes[i]));
            return 1;
        }
    }
    return 0;
}

int compare_split_prime_search(const sflint::FmpzPoly& field) noexcept {
    relation_search::SplitPrimeResultant evaluator;
    return compare_split_prime_search(evaluator, field);
}

// The discriminant screen in the split-prime search skips only primes
// that cannot split f into distinct linear factors, so the search finds
// the same primes in the same order after the same number of tested
// primes.  The fields cover a non-square discriminant (Phi_17 and the
// degree-17 real abelian field), a square one (Phi_32 = x^16 + 1, discriminant
// 2^64), and a first tested prime q = 2^27 + 29 dividing the discriminant:
// x (x - q) (x - 1) ... (x - 14) splits over Q (square discriminant) and
// x (x - q) Phi_17 does not (non-square), and modulo q both have the double
// root 0.
int test_screened_search_matches_unscreened_search() {
    const ulong q = n_nextprime(UWORD(1) << 27, 1);
    if (q != (UWORD(1) << 27) + 29) {
        return 1;
    }
    sflint::FmpzPoly field;
    set_prime_cyclotomic(field, 17);
    if (compare_split_prime_search(field) != 0) {
        return 1;
    }
    set_real_abelian_degree_17(field);
    if (compare_split_prime_search(field) != 0) {
        return 1;
    }
    fmpz_poly_zero(field.raw());
    fmpz_poly_set_coeff_si(field.raw(), 0, 1);
    fmpz_poly_set_coeff_si(field.raw(), 16, 1);
    if (compare_split_prime_search(field) != 0) {
        return 1;
    }

    // A negative discriminant: (x^2 + x + 2) (x - 1) ... (x - 14) has one
    // complex pair (r2 = 1), and disc(x^2 + x + 2) = -7 is a square modulo
    // primes of both residues modulo 4, so a sign mistake in the screen
    // would skip primes that split it.
    {
        sflint::FmpzPoly negative;
        fmpz_poly_set_coeff_si(negative.raw(), 0, 2);
        fmpz_poly_set_coeff_si(negative.raw(), 1, 1);
        fmpz_poly_set_coeff_si(negative.raw(), 2, 1);
        sflint::FmpzPoly factor;
        fmpz_poly_set_coeff_si(factor.raw(), 1, 1);
        for (slong a = 1; a <= 14; ++a) {
            fmpz_poly_set_coeff_si(factor.raw(), 0, -a);
            fmpz_poly_mul(negative.raw(), negative.raw(), factor.raw());
        }
        sflint::Fmpz disc;
        fmpz_poly_discriminant(disc.raw(), negative.raw());
        if (fmpz_sgn(disc.raw()) >= 0 ||
            compare_split_prime_search(negative) != 0) {
            return 1;
        }
    }

    // x (x - q), then x (x - q) (x - 1) ... (x - 14).
    fmpz_poly_zero(field.raw());
    fmpz_poly_set_coeff_si(field.raw(), 1, -static_cast<slong>(q));
    fmpz_poly_set_coeff_si(field.raw(), 2, 1);
    sflint::FmpzPoly linear;
    fmpz_poly_set_coeff_si(linear.raw(), 1, 1);
    sflint::FmpzPoly split_field;
    fmpz_poly_set(split_field.raw(), field.raw());
    for (slong a = 1; a <= 14; ++a) {
        fmpz_poly_set_coeff_si(linear.raw(), 0, -a);
        fmpz_poly_mul(split_field.raw(), split_field.raw(), linear.raw());
    }
    sflint::FmpzPoly cyclotomic;
    set_prime_cyclotomic(cyclotomic, 17);
    fmpz_poly_mul(field.raw(), field.raw(), cyclotomic.raw());
    sflint::Fmpz discriminant;
    for (const sflint::FmpzPoly* f : {&split_field, &field}) {
        fmpz_poly_discriminant(discriminant.raw(), f->raw());
        if (fmpz_fdiv_ui(discriminant.raw(), q) != 0 ||
            compare_split_prime_search(*f) != 0) {
            return 1;
        }
    }
    // The fully split field finds its 64 primes among the first 65 primes
    // tested: every prime except q splits it.
    ReferenceSearch reference;
    reference_split_prime_search(reference, split_field);
    return reference.split_primes.size() == 64 &&
                           reference.tested_primes == 65 &&
                           reference.split_primes.front() != q
                   ? 0
                   : 1;
}

// The closed-form split test for Phi_n (p = 1 mod n) and psi_n, the
// minimal polynomial of 2 cos(2 pi / n) (p = +-1 mod n), against the
// generic search: the same split primes in the same order, the same
// number of primes tested, the same root sets, and resultant values equal
// to fmpz_poly_resultant.  The fields cover prime, prime-power, odd
// composite and n = 2 mod 4 indices: Phi_17, Phi_19, Phi_32 = x^16 + 1,
// Phi_34, Phi_60, psi_37 (the cyclotomic_real_37_deg18 polynomial), psi_51,
// psi_64 and psi_120.  The degree-17 subfield of Q(zeta_103) and
// Phi_17(x + 1) take the generic test.
int test_closed_form_split_primes_match_reference() {
    using ClosedForm = relation_search::SplitPrimeResultant::SplitClosedForm;
    struct Case {
        bool real;
        ulong n;
    };
    const Case cases[] = {{false, 17}, {false, 19}, {false, 32},
                          {false, 34}, {false, 60}, {true, 37},
                          {true, 51},  {true, 64},  {true, 120}};
    for (const Case& c : cases) {
        sflint::FmpzPoly field;
        if (c.real) {
            fmpz_poly_cos_minpoly(field.raw(), c.n);
        } else {
            fmpz_poly_cyclotomic(field.raw(), c.n);
        }
        relation_search::SplitPrimeResultant evaluator;
        const int search_differs = compare_split_prime_search(evaluator, field);
        if (search_differs != 0 ||
            evaluator.split_closed_form !=
                    (c.real ? ClosedForm::real_cyclotomic
                            : ClosedForm::cyclotomic) ||
            evaluator.split_closed_form_index != c.n ||
            evaluator.closed_form_split_primes !=
                    static_cast<slong>(evaluator.moduli.size())) {
            std::fprintf(stderr, "closed form differs for %s_%lu\n",
                         c.real ? "psi" : "Phi",
                         static_cast<unsigned long>(c.n));
            return 1;
        }
        relation_search::SplitPrimeResultant values;
        values.reset(field.raw());
        if (compare_with_resultant(values, field, 211 + c.n, 40, 64) <= 0 ||
            values.closed_form_split_primes == 0) {
            return 1;
        }
    }
    // cyclotomic_real_37_deg18 is psi_37.
    {
        constexpr slong coefficients[] = {-1,  9,    45,  -120, -330,
                                          462, 924,  -792, -1287, 715,
                                          1001, -364, -455, 105,  120,
                                          -16, -17,  1,    1};
        sflint::FmpzPoly row;
        slong i = 0;
        for (const slong coefficient : coefficients) {
            fmpz_poly_set_coeff_si(row.raw(), i, coefficient);
            ++i;
        }
        sflint::FmpzPoly psi;
        fmpz_poly_cos_minpoly(psi.raw(), 37);
        if (fmpz_poly_equal(row.raw(), psi.raw()) == 0) {
            return 1;
        }
    }
    sflint::FmpzPoly field;
    set_real_abelian_degree_17(field);
    sflint::FmpzPoly translated;
    set_prime_cyclotomic(translated, 17);
    sflint::FmpzPoly shift;
    fmpz_poly_set_coeff_si(shift.raw(), 0, 1);
    fmpz_poly_set_coeff_si(shift.raw(), 1, 1);
    fmpz_poly_compose(translated.raw(), translated.raw(), shift.raw());
    for (const sflint::FmpzPoly* f : {&field, &translated}) {
        relation_search::SplitPrimeResultant evaluator;
        if (compare_split_prime_search(evaluator, *f) != 0 ||
            evaluator.split_closed_form != ClosedForm::none ||
            evaluator.closed_form_split_primes != 0) {
            return 1;
        }
    }
    return 0;
}

}  // namespace

// The closed-form degree-2 resultant equals fmpz_poly_resultant on random
// f of degree 2 (monic and non-monic, either sign of the leading
// coefficient) and g of length 0, 1 or 2, and, for g / d with an integer
// denominator d, Res(f, g) / d^2 equals fmpq_poly_resultant(f, g / d).
int test_degree_two_resultant_matches_flint() {
    ulong state = 2026;
    sflint::FmpzPoly field;
    sflint::FmpzPoly coordinate;
    sflint::Fmpz expected;
    sflint::Fmpz actual;
    sflint::Fmpz den;
    fmpq_poly_t rational_field;
    fmpq_poly_t rational_coordinate;
    fmpq_t rational_expected;
    fmpq_t rational_actual;
    fmpq_poly_init(rational_field);
    fmpq_poly_init(rational_coordinate);
    fmpq_init(rational_expected);
    fmpq_init(rational_actual);
    int failures = 0;
    for (slong k = 0; k < 4000 && failures == 0; ++k) {
        const flint_bitcnt_t field_bits =
                static_cast<flint_bitcnt_t>(next_random(state) % 161);
        set_random_polynomial(field, state, 3, field_bits);
        if (k % 4 == 0) {
            fmpz_poly_set_coeff_si(field.raw(), 2, 1);
        }
        const slong length = static_cast<slong>(next_random(state) % 3);
        const flint_bitcnt_t bits =
                static_cast<flint_bitcnt_t>(next_random(state) % 161);
        set_random_polynomial(coordinate, state, length, bits);
        fmpz_poly_resultant(expected.raw(), field.raw(), coordinate.raw());
        bool handled = false;
        relation_search::degree_two_resultant(handled, actual.raw(),
                                              field.raw(), coordinate.raw());
        if (!handled || !fmpz_equal(expected.raw(), actual.raw())) {
            std::fprintf(stderr,
                         "degree-2 resultant mismatch: candidate %ld\n",
                         static_cast<long>(k));
            ++failures;
            break;
        }

        fmpz_set_ui(den.raw(), 1 + (next_random(state) & UWORD(0xffffffffff)));
        fmpq_poly_set_fmpz_poly(rational_field, field.raw());
        fmpq_poly_set_fmpz_poly(rational_coordinate, coordinate.raw());
        fmpq_poly_scalar_div_fmpz(rational_coordinate, rational_coordinate,
                                  den.raw());
        fmpq_poly_resultant(rational_expected, rational_field,
                            rational_coordinate);
        fmpz_mul(den.raw(), den.raw(), den.raw());
        fmpq_set_fmpz_frac(rational_actual, actual.raw(), den.raw());
        if (!fmpq_equal(rational_expected, rational_actual)) {
            std::fprintf(stderr,
                         "degree-2 resultant denominator mismatch: "
                         "candidate %ld\n",
                         static_cast<long>(k));
            ++failures;
        }
    }
    fmpq_clear(rational_actual);
    fmpq_clear(rational_expected);
    fmpq_poly_clear(rational_coordinate);
    fmpq_poly_clear(rational_field);
    if (failures != 0) {
        return 1;
    }

    // Res(-3x^2 + 5x + 7, 2) = 4, not lc(f) * 4.
    fmpz_poly_zero(field.raw());
    fmpz_poly_set_coeff_si(field.raw(), 0, 7);
    fmpz_poly_set_coeff_si(field.raw(), 1, 5);
    fmpz_poly_set_coeff_si(field.raw(), 2, -3);
    fmpz_poly_set_si(coordinate.raw(), 2);
    bool handled = false;
    relation_search::degree_two_resultant(handled, actual.raw(), field.raw(),
                                          coordinate.raw());
    if (!handled || !fmpz_equal_si(actual.raw(), 4)) {
        return 1;
    }

    // Other degrees are left to the caller and leave out unchanged.
    fmpz_set_si(actual.raw(), 99);
    fmpz_poly_set_coeff_si(coordinate.raw(), 2, 1);
    relation_search::degree_two_resultant(handled, actual.raw(), field.raw(),
                                          coordinate.raw());
    if (handled || !fmpz_equal_si(actual.raw(), 99)) {
        return 1;
    }
    fmpz_poly_set_si(coordinate.raw(), 2);
    fmpz_poly_set_coeff_si(field.raw(), 3, 1);
    relation_search::degree_two_resultant(handled, actual.raw(), field.raw(),
                                          coordinate.raw());
    return handled || !fmpz_equal_si(actual.raw(), 99) ? 1 : 0;
}

// Sets out to row left of basis, or to row left + sign * row right.
void set_row_combination_polynomial(sflint::FmpzPoly& out,
                                    const sflint::FmpzMat& basis,
                                    slong left,
                                    slong right,
                                    slong sign) noexcept {
    fmpz_poly_zero(out.raw());
    sflint::Fmpz coefficient;
    for (slong j = 0; j < fmpz_mat_ncols(basis.raw()); ++j) {
        fmpz_set(coefficient.raw(), fmpz_mat_entry(basis.raw(), left, j));
        if (right >= 0 && sign > 0) {
            fmpz_add(coefficient.raw(), coefficient.raw(),
                     fmpz_mat_entry(basis.raw(), right, j));
        } else if (right >= 0) {
            fmpz_sub(coefficient.raw(), coefficient.raw(),
                     fmpz_mat_entry(basis.raw(), right, j));
        }
        fmpz_poly_set_coeff_fmpz(out.raw(), j, coefficient.raw());
    }
}

// Rows and sums and differences of two rows of a basis, in the small-LLL
// candidate order: the per-row values give the value of the direct
// evaluation and of fmpz_poly_resultant, with the same split primes.  The
// values are filled lazily, once per (row, prime): a row first reached by
// a candidate needing few primes is extended when a later candidate needs
// more, and is reused, not refilled, by every candidate covered by the
// primes already filled.  A description that does not match the
// polynomial, or one with no basis bound, is not used, and a description
// applies to one call only.
int test_row_combinations_match_direct() {
    for (const slong m : {17, 41}) {
        sflint::FmpzPoly field;
        set_prime_cyclotomic(field, m);
        const slong n = m - 1;
        relation_search::SplitPrimeResultant cached;
        relation_search::SplitPrimeResultant direct;
        cached.reset(field.raw());
        direct.reset(field.raw());
        if (!cached.enabled || !direct.enabled) {
            return 1;
        }

        // Small rows first, then rows whose bounds need more primes.
        constexpr slong rows = 6;
        constexpr flint_bitcnt_t row_bits[rows] = {2, 3, 3, 10, 16, 24};
        sflint::FmpzMat basis(rows, n);
        sflint::FmpzPoly row;
        ulong state = static_cast<ulong>(m) * 7;
        for (slong r = 0; r < rows; ++r) {
            set_random_polynomial(row, state, n, row_bits[r]);
            for (slong j = 0; j < n; ++j) {
                fmpz_poly_get_coeff_fmpz(fmpz_mat_entry(basis.raw(), r, j),
                                         row.raw(), j);
            }
        }

        cached.bind_row_basis(basis.raw());
        std::vector<slong> filled(static_cast<std::size_t>(rows), 0);
        slong expected_fills = 0;
        slong reused = 0;
        slong extended = 0;
        slong handled_count = 0;
        sflint::FmpzPoly candidate;
        sflint::Fmpz expected;
        sflint::Fmpz via_rows;
        sflint::Fmpz via_direct;
        const auto check = [&](slong left, slong right, slong sign) -> bool {
            set_row_combination_polynomial(candidate, basis, left, right,
                                           sign);
            fmpz_poly_resultant(expected.raw(), field.raw(),
                                candidate.raw());
            bool rows_handled = false;
            bool direct_handled = false;
            cached.set_row_combination(left, right, sign);
            cached.resultant(rows_handled, via_rows.raw(), candidate.raw());
            direct.resultant(direct_handled, via_direct.raw(),
                             candidate.raw());
            if (fmpz_poly_is_zero(candidate.raw()) != 0) {
                // b_l - b_l: neither path handles the zero polynomial.
                return !rows_handled && !direct_handled &&
                       cached.row_value_fills == expected_fills;
            }
            if (!rows_handled || !direct_handled ||
                !cached.last_used_row_values ||
                cached.last_prime_count != direct.last_prime_count ||
                !fmpz_equal(via_rows.raw(), expected.raw()) ||
                !fmpz_equal(via_direct.raw(), expected.raw())) {
                std::fprintf(stderr,
                             "row-combination resultant mismatch: degree "
                             "%ld, rows %ld %ld sign %ld\n",
                             static_cast<long>(n), static_cast<long>(left),
                             static_cast<long>(right),
                             static_cast<long>(sign));
                return false;
            }
            ++handled_count;
            const slong primes = cached.last_prime_count;
            // b_l + b_l reaches one row.
            for (const slong r : {left, right == left ? -1 : right}) {
                if (r < 0) {
                    continue;
                }
                slong& done = filled[static_cast<std::size_t>(r)];
                if (done >= primes) {
                    ++reused;
                } else {
                    if (done > 0) {
                        ++extended;
                    }
                    expected_fills += primes - done;
                    done = primes;
                }
            }
            return cached.row_value_fills == expected_fills;
        };

        for (slong left = 0; left < rows; ++left) {
            if (!check(left, -1, 0)) {
                return 1;
            }
        }
        for (const slong sign : {1, -1}) {
            for (slong left = 0; left < rows; ++left) {
                for (slong right = left; right < rows; ++right) {
                    if (!check(left, right, sign)) {
                        return 1;
                    }
                }
            }
        }
        if (handled_count != rows + rows * rows || reused <= 0 ||
            extended <= 0) {
            std::fprintf(stderr,
                         "row-combination coverage: degree %ld, handled "
                         "%ld, reused %ld, extended %ld\n",
                         static_cast<long>(n),
                         static_cast<long>(handled_count),
                         static_cast<long>(reused),
                         static_cast<long>(extended));
            return 1;
        }

        // A description that does not match the polynomial is not used.
        set_row_combination_polynomial(candidate, basis, 0, 1, -1);
        fmpz_poly_resultant(expected.raw(), field.raw(), candidate.raw());
        bool handled = false;
        cached.set_row_combination(0, 1, 1);
        cached.resultant(handled, via_rows.raw(), candidate.raw());
        if (!handled || cached.last_used_row_values ||
            !fmpz_equal(via_rows.raw(), expected.raw()) ||
            cached.row_value_fills != expected_fills) {
            return 1;
        }
        // A matching description applies to the next call only.
        cached.set_row_combination(0, 1, -1);
        cached.resultant(handled, via_rows.raw(), candidate.raw());
        if (!handled || !cached.last_used_row_values) {
            return 1;
        }
        cached.resultant(handled, via_rows.raw(), candidate.raw());
        if (!handled || cached.last_used_row_values ||
            !fmpz_equal(via_rows.raw(), expected.raw())) {
            return 1;
        }
        // With no basis bound the values are not used.
        cached.bind_row_basis(nullptr);
        cached.set_row_combination(0, 1, -1);
        cached.resultant(handled, via_rows.raw(), candidate.raw());
        if (!handled || cached.last_used_row_values ||
            !fmpz_equal(via_rows.raw(), expected.raw()) ||
            cached.row_value_fills != 0) {
            return 1;
        }
    }
    return 0;
}

// Row values are dropped when they stop describing the bound rows: a
// rebind to a different basis refills the rows from the new entries, and
// a reset to a different f (deg g = deg f, so the description still
// fits) while a basis is bound refills them for the new roots.  Stale
// values from either would give a wrong norm with a matching description.
// An unhandled call also clears last_used_row_values.
int test_row_values_invalidated() {
    constexpr slong m = 17;
    constexpr slong n = m - 1;
    constexpr slong rows = 4;
    sflint::FmpzPoly field;
    set_prime_cyclotomic(field, m);
    relation_search::SplitPrimeResultant cached;
    cached.reset(field.raw());
    if (!cached.enabled) {
        return 1;
    }

    // Two bases whose rows differ in every row.
    ulong state = 12345;
    sflint::FmpzMat first(rows, n);
    sflint::FmpzMat second(rows, n);
    sflint::FmpzPoly row;
    for (sflint::FmpzMat* basis : {&first, &second}) {
        for (slong r = 0; r < rows; ++r) {
            set_random_polynomial(row, state, n, 8);
            for (slong j = 0; j < n; ++j) {
                fmpz_poly_get_coeff_fmpz(
                        fmpz_mat_entry(basis->raw(), r, j), row.raw(), j);
            }
        }
    }
    for (slong r = 0; r < rows; ++r) {
        bool differs = false;
        for (slong j = 0; j < n && !differs; ++j) {
            differs = fmpz_equal(fmpz_mat_entry(first.raw(), r, j),
                                 fmpz_mat_entry(second.raw(), r, j)) == 0;
        }
        if (!differs) {
            return 1;
        }
    }

    sflint::FmpzPoly candidate;
    sflint::Fmpz expected;
    sflint::Fmpz actual;
    // Every row and row sum of basis, from row values, against
    // fmpz_poly_resultant(f, .); each row is filled once per prime.
    const auto check_basis = [&](const sflint::FmpzPoly& f,
                                 const sflint::FmpzMat& basis,
                                 slong fills_before) -> bool {
        std::vector<slong> filled(static_cast<std::size_t>(rows), 0);
        slong expected_fills = fills_before;
        for (slong left = 0; left < rows; ++left) {
            // Row left alone, then row left plus the row before it.
            for (const slong right : {slong{-1}, left - 1}) {
                set_row_combination_polynomial(candidate, basis, left,
                                               right, 1);
                fmpz_poly_resultant(expected.raw(), f.raw(),
                                    candidate.raw());
                bool handled = false;
                cached.set_row_combination(left, right, 1);
                cached.resultant(handled, actual.raw(), candidate.raw());
                if (!handled || !cached.last_used_row_values ||
                    !fmpz_equal(actual.raw(), expected.raw())) {
                    std::fprintf(stderr,
                                 "stale row values: rows %ld %ld\n",
                                 static_cast<long>(left),
                                 static_cast<long>(right));
                    return false;
                }
                const slong primes = cached.last_prime_count;
                for (const slong r : {left, right}) {
                    if (r < 0) {
                        continue;
                    }
                    slong& done = filled[static_cast<std::size_t>(r)];
                    if (done < primes) {
                        expected_fills += primes - done;
                        done = primes;
                    }
                }
                if (cached.row_value_fills != expected_fills) {
                    std::fprintf(stderr,
                                 "row value fills %ld, expected %ld\n",
                                 static_cast<long>(cached.row_value_fills),
                                 static_cast<long>(expected_fills));
                    return false;
                }
            }
        }
        return true;
    };

    cached.bind_row_basis(first.raw());
    if (!check_basis(field, first, 0)) {
        return 1;
    }
    // An unhandled call (the zero polynomial) does not report row values.
    fmpz_poly_zero(candidate.raw());
    bool handled = true;
    cached.resultant(handled, actual.raw(), candidate.raw());
    if (handled || cached.last_used_row_values) {
        return 1;
    }

    // Rebinding a different basis restarts the fills.
    cached.bind_row_basis(second.raw());
    if (cached.row_value_fills != 0 || !check_basis(field, second, 0)) {
        return 1;
    }

    // Reset to g(x) = f(x + 1) with the basis still bound: the same split
    // primes, in the same order, with shifted roots.
    sflint::FmpzPoly shifted;
    sflint::Fmpz one;
    fmpz_one(one.raw());
    fmpz_poly_taylor_shift(shifted.raw(), field.raw(), one.raw());
    cached.reset(shifted.raw());
    if (!cached.enabled ||
        !check_basis(shifted, second, cached.row_value_fills)) {
        return 1;
    }
    cached.bind_row_basis(nullptr);
    return 0;
}

int main() {
    return test_cyclotomic_matches_resultant() != 0 ||
                   test_screened_search_matches_unscreened_search() != 0 ||
                   test_closed_form_split_primes_match_reference() != 0 ||
                   test_real_abelian_matches_resultant() != 0 ||
                   test_constant_and_zero_polynomials() != 0 ||
                   test_disabled_and_exhausted_cases() != 0 ||
                   test_reset_keeps_split_primes_for_same_field() != 0 ||
                   test_t2_bound_matches_resultant_on_tight_candidates() !=
                           0 ||
                   test_t2_bound_for_prime_cyclotomic() != 0 ||
                   test_degree_two_resultant_matches_flint() != 0 ||
                   test_row_combinations_match_direct() != 0 ||
                   test_row_values_invalidated() != 0
            ? 1
            : 0;
}
