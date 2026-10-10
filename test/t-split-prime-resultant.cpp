#include "class_group/relation_candidate_internal.hpp"

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include <flint/fmpz_poly.h>

#include <cstddef>
#include <cstdio>

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

}  // namespace

int main() {
    return test_cyclotomic_matches_resultant() != 0 ||
                   test_real_abelian_matches_resultant() != 0 ||
                   test_constant_and_zero_polynomials() != 0 ||
                   test_disabled_and_exhausted_cases() != 0 ||
                   test_reset_keeps_split_primes_for_same_field() != 0
            ? 1
            : 0;
}
