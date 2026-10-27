#include "reduced_form_count_internal.hpp"

#include <flint/ulong_extras.h>

#include <cstdint>
#include <vector>

namespace silex::detail {
namespace {

// Owns a root array allocated by FLINT `n_sqrtmod_primepow` or `n_sqrtmodn`.
struct FlintRootArray {
    FlintRootArray() noexcept = default;
    FlintRootArray(const FlintRootArray&) = delete;
    FlintRootArray& operator=(const FlintRootArray&) = delete;
    ~FlintRootArray() noexcept {
        if (roots != nullptr) {
            ::flint_free(roots);
        }
    }

    ulong* roots = nullptr;
};

// R(p^e): the number of roots of x^2 = -abs_discriminant mod p^e, from FLINT
// `n_sqrtmod_primepow`.  p^e divides 4a <= 4 sqrt(|D| / 3), so it fits a
// word.
ulong prime_power_root_count(ulong abs_discriminant, ulong p,
                             slong exponent) noexcept {
    const ulong modulus = ::n_pow(p, static_cast<ulong>(exponent));
    FlintRootArray square_roots;
    const slong count = ::n_sqrtmod_primepow(
            &square_roots.roots,
            ::n_negmod(abs_discriminant % modulus, modulus), p, exponent);
    return static_cast<ulong>(count);
}

// The exponent of a prime power exactly dividing a word fits this bound.
constexpr int kMaxExponent = FLINT_BITS;

}  // namespace

bool imaginary_quadratic_discriminant_is_fundamental(
        ulong abs_discriminant) noexcept {
    // D = 1 mod 4, that is |D| = 3 mod 4: |D| squarefree.
    if ((abs_discriminant & 3) == 3) {
        return ::n_is_squarefree(abs_discriminant) != 0;
    }
    // D = 4m: m = -|D| / 4 = 2 or 3 mod 4, that is |D| / 4 = 2 or 1 mod 4,
    // and m squarefree.
    if ((abs_discriminant & 3) == 0) {
        const ulong quarter = abs_discriminant >> 2;
        const ulong residue = quarter & 3;
        return (residue == 1 || residue == 2) &&
               ::n_is_squarefree(quarter) != 0;
    }
    return false;
}

bool imaginary_quadratic_fundamental_reduced_form_count(
        ulong& out,
        ulong abs_discriminant) noexcept {
    const ulong n = abs_discriminant;
    if (n >= (UWORD(1) << 62) ||
        !imaginary_quadratic_discriminant_is_fundamental(n)) {
        return false;
    }
    // Every reduced form has a <= sqrt(|D| / 3) (Cohen, GTM 138,
    // Lemma 5.3.4(1)); floor(sqrt(floor(n / 3))) = floor(sqrt(n / 3)).
    const ulong a_limit = ::n_sqrt(n / 3);

    // root_count[a] = R(4a) = prod over p^e || 4a of R(p^e) (Chinese
    // remainder theorem, Cohen, GTM 138, Theorem 1.3.9 and Corollary 1.3.10),
    // built by a sieve over the primes p <= a_limit.  R(p^e) depends only on
    // (p, e), so it is computed once per pair.  Each R(p^e) is at most 4
    // (p = 2) or 2 (odd p), and a <= a_limit has fewer than 32 distinct
    // prime factors, so the product fits 32 bits.
    std::vector<std::uint32_t> root_count(a_limit + 1, 0);

    // p = 2: 4a is exactly divisible by 2^(v_2(a) + 2).
    {
        ulong count_at[kMaxExponent + 2] = {};
        bool known[kMaxExponent + 2] = {};
        for (ulong a = 1; a <= a_limit; ++a) {
            const int exponent = 2 + static_cast<int>(::flint_ctz(a));
            if (!known[exponent]) {
                count_at[exponent] =
                        prime_power_root_count(n, 2, exponent);
                known[exponent] = true;
            }
            root_count[a] = static_cast<std::uint32_t>(count_at[exponent]);
        }
    }

    // Odd primes p <= a_limit, found by a sieve of Eratosthenes over the odd
    // integers.
    std::vector<unsigned char> composite(a_limit + 1, 0);
    for (ulong p = 3; p <= a_limit; p += 2) {
        if (composite[p] != 0) {
            continue;
        }
        if (p <= a_limit / p) {
            for (ulong multiple = p * p; multiple <= a_limit;
                 multiple += 2 * p) {
                composite[multiple] = 1;
            }
        }
        ulong count_at[kMaxExponent + 1] = {};
        bool known[kMaxExponent + 1] = {};
        for (ulong a = p; a <= a_limit; a += p) {
            if (root_count[a] == 0) {
                continue;
            }
            ulong cofactor = a / p;
            int exponent = 1;
            while (cofactor % p == 0) {
                cofactor /= p;
                ++exponent;
            }
            if (!known[exponent]) {
                count_at[exponent] =
                        prime_power_root_count(n, p, exponent);
                known[exponent] = true;
            }
            root_count[a] *= static_cast<std::uint32_t>(count_at[exponent]);
        }
    }

    ulong count = 0;
    for (ulong a = 1; a <= a_limit; ++a) {
        const ulong roots = root_count[a];
        if (roots == 0) {
            continue;
        }
        const ulong four_a = 4 * a;
        if (a * four_a < n) {
            // 4a^2 < |D|: every b with -a < b <= a gives a reduced form with
            // c > a (Cohen, GTM 138, Lemma 5.3.4(2)).  b is defined mod 2a
            // (Cohen, p. 253), and b and b + 2a are both roots mod 4a, so
            // the window holds exactly half of the R(4a) roots.
            count += roots / 2;
            continue;
        }

        // 4a^2 >= |D|: list the roots b mod 4a and keep those that give a
        // reduced form (a, b, c), c = (b^2 + |D|) / (4a): |b| <= a <= c, and
        // b >= 0 when |b| = a or a = c (Cohen, GTM 138, Definition 5.3.2).
        n_factor_t factors;
        ::n_factor_init(&factors);
        ::n_factor(&factors, four_a, 1);
        FlintRootArray square_roots;
        const slong listed = ::n_sqrtmodn(
                &square_roots.roots, ::n_negmod(n % four_a, four_a),
                &factors);
        if (listed < 0 || static_cast<ulong>(listed) != roots) {
            return false;
        }
        for (slong i = 0; i < listed; ++i) {
            // The representative of b in (-a, a]; the other residues mod
            // 4a are b + 2a and lie outside it.
            const ulong root = square_roots.roots[i];
            ulong abs_b = 0;
            bool negative = false;
            if (root <= a) {
                abs_b = root;
            } else if (root > 3 * a) {
                abs_b = four_a - root;
                negative = true;
            } else {
                continue;
            }
            const ulong c = (abs_b * abs_b + n) / four_a;
            if (c >= a && (!negative || c != a)) {
                ++count;
            }
        }
    }

    out = count;
    return true;
}

}  // namespace silex::detail
