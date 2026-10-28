#include "reduced_form_count_internal.hpp"

#include <flint/ulong_extras.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace silex::detail {
namespace {

// Owns a root array allocated by FLINT `n_sqrtmod_primepow`.
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

// For a fundamental D, R(2^e) is at most 4 and R(p^e) at most 2 for odd p.
constexpr slong kMaxPrimePowerRoots = 4;

// The roots of x^2 = -abs_discriminant mod p^e, from FLINT
// `n_sqrtmod_primepow`; their number is R(p^e).
struct PrimePowerRoots {
    ulong count = 0;
    ulong roots[kMaxPrimePowerRoots] = {};
};

// Fills `out` with the roots of x^2 = -abs_discriminant mod p^e.  p^e
// divides 4a <= 4 sqrt(|D| / 3), so it fits a word.  Returns false when FLINT
// reports more roots than a fundamental D allows.
bool prime_power_roots(PrimePowerRoots& out, ulong abs_discriminant, ulong p,
                       slong exponent) noexcept {
    const ulong modulus = ::n_pow(p, static_cast<ulong>(exponent));
    FlintRootArray square_roots;
    const slong count = ::n_sqrtmod_primepow(
            &square_roots.roots,
            ::n_negmod(abs_discriminant % modulus, modulus), p, exponent);
    if (count < 0 || count > kMaxPrimePowerRoots) {
        return false;
    }
    out.count = static_cast<ulong>(count);
    for (slong i = 0; i < count; ++i) {
        out.roots[i] = square_roots.roots[i];
    }
    return true;
}

// The exponent of a prime power exactly dividing a word fits this bound.
constexpr int kMaxExponent = FLINT_BITS;

// The roots mod an odd prime power q = p^e < 2^31 are stored in one word:
// kNoRoot when there is none, 0 when 0 is the only root, and otherwise one
// root r != 0 of the pair {r, q - r} (the roots are closed under negation, and
// r = -r mod q only for r = 0).  kUnknown marks a pair not computed yet.
constexpr std::uint32_t kUnknown = 0xFFFFFFFFU;
constexpr std::uint32_t kNoRoot = 0xFFFFFFFEU;

std::uint32_t odd_root_count(std::uint32_t root) noexcept {
    if (root == kNoRoot) {
        return 0;
    }
    return root == 0 ? 1 : 2;
}

// Encodes the roots mod an odd p^e as above; returns false when they are not
// {}, {0} or {r, q - r} with r != 0.
bool encode_odd_roots(std::uint32_t& out, const PrimePowerRoots& roots,
                      ulong modulus) noexcept {
    if (roots.count == 0) {
        out = kNoRoot;
        return true;
    }
    const ulong root = roots.roots[0];
    if (roots.count == 1) {
        out = 0;
        return root == 0;
    }
    if (roots.count != 2 || root == 0 || root >= modulus ||
        roots.roots[1] != modulus - root) {
        return false;
    }
    out = static_cast<std::uint32_t>(root);
    return true;
}

// Exponents e >= 2 occur only for p <= sqrt(a_limit) < 2^16, so e < 32.
constexpr ulong kOddExponentSlots = 32;

// The roots mod p^e for the odd primes p <= a_limit, indexed by the position
// of p in `primes`; exponent 1 for every prime, exponents >= 2 for the primes
// with p^2 <= a_limit, which come first in `primes`.  Each pair (p, e) is
// computed once, by FLINT `n_sqrtmod_primepow`.
struct OddRootCache {
    ulong abs_discriminant = 0;
    std::vector<std::uint32_t> primes;
    std::vector<std::uint32_t> first_power;
    std::vector<std::uint32_t> higher_powers;

    bool root(std::uint32_t& out, std::size_t index, int exponent) noexcept {
        std::uint32_t* slot = nullptr;
        if (exponent == 1) {
            slot = &first_power[index];
        } else {
            if (exponent < 2 ||
                static_cast<ulong>(exponent) >= kOddExponentSlots) {
                return false;
            }
            const std::size_t position =
                    index * kOddExponentSlots +
                    static_cast<std::size_t>(exponent);
            if (position >= higher_powers.size()) {
                return false;
            }
            slot = &higher_powers[position];
        }
        if (*slot == kUnknown) {
            const ulong p = primes[index];
            PrimePowerRoots roots;
            if (!prime_power_roots(roots, abs_discriminant, p, exponent) ||
                !encode_odd_roots(*slot, roots,
                                  ::n_pow(p, static_cast<ulong>(exponent)))) {
                return false;
            }
        }
        out = *slot;
        return true;
    }
};

// A band a <= a_limit < 2^31 has at most 8 distinct odd prime factors
// (3 * 5 * ... * 29 > 2^31), so 4a has at most 4 * 2^8 square roots of D.
constexpr std::size_t kMaxOddFactors = 8;
constexpr std::size_t kMaxBandRoots = std::size_t{4} << kMaxOddFactors;
constexpr std::uint32_t kExponentBits = 5;
constexpr std::uint32_t kExponentMask = (1U << kExponentBits) - 1;

// Band chunks hold at least this many a, and the band takes at most
// kMaxBandChunks chunks, so the band factorisation sieve visits every odd
// prime at most kMaxBandChunks times.
constexpr ulong kMinBandChunk = UWORD(1) << 16;
constexpr ulong kMaxBandChunks = 16;

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
    // Lemma 5.3.4(1)); floor(sqrt(floor(n / 3))) = floor(sqrt(n / 3)), and
    // n < 2^62 gives a_limit < 2^31.
    const ulong a_limit = ::n_sqrt(n / 3);

    // root_count[a] = R(4a) = prod over p^e || 4a of R(p^e) (Chinese
    // remainder theorem, Cohen, GTM 138, Theorem 1.3.9 and Corollary 1.3.10),
    // built by a sieve over the primes p <= a_limit.  R(p^e) depends only on
    // (p, e), so its roots are computed once per pair and kept for the band
    // below.  Each R(p^e) is at most 4 (p = 2) or 2 (odd p), and a <= a_limit
    // has fewer than 32 distinct prime factors, so the product fits 32 bits.
    std::vector<std::uint32_t> root_count(a_limit + 1, 0);

    // p = 2: 4a is exactly divisible by 2^(v_2(a) + 2).
    PrimePowerRoots two_power_roots[kMaxExponent + 2];
    bool two_power_known[kMaxExponent + 2] = {};
    for (ulong a = 1; a <= a_limit; ++a) {
        const int exponent = 2 + static_cast<int>(::flint_ctz(a));
        if (!two_power_known[exponent]) {
            if (!prime_power_roots(two_power_roots[exponent], n, 2,
                                   exponent)) {
                return false;
            }
            two_power_known[exponent] = true;
        }
        root_count[a] =
                static_cast<std::uint32_t>(two_power_roots[exponent].count);
    }

    // Odd primes p <= a_limit, found by a sieve of Eratosthenes over the odd
    // integers (bit (x - 1) / 2 is set for an odd composite x).
    OddRootCache odd_roots;
    odd_roots.abs_discriminant = n;
    std::vector<std::uint64_t> odd_composite((a_limit >> 7) + 1, 0);
    for (ulong p = 3; p <= a_limit; p += 2) {
        const ulong bit = p >> 1;
        if (((odd_composite[bit >> 6] >> (bit & 63)) & 1) != 0) {
            continue;
        }
        if (p <= a_limit / p) {
            for (ulong multiple = p * p; multiple <= a_limit;
                 multiple += 2 * p) {
                const ulong multiple_bit = multiple >> 1;
                odd_composite[multiple_bit >> 6] |= UWORD(1)
                                                    << (multiple_bit & 63);
            }
        }
        const std::size_t index = odd_roots.primes.size();
        odd_roots.primes.push_back(static_cast<std::uint32_t>(p));
        odd_roots.first_power.push_back(kUnknown);
        if (p <= a_limit / p) {
            odd_roots.higher_powers.resize((index + 1) * kOddExponentSlots,
                                           kUnknown);
        }
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
            std::uint32_t root = 0;
            if (!odd_roots.root(root, index, exponent)) {
                return false;
            }
            root_count[a] *= odd_root_count(root);
        }
    }
    odd_composite = std::vector<std::uint64_t>();

    // The band is a >= band_start, where 4a^2 >= |D|; below it 4a^2 < |D|.
    ulong band_start = ::n_sqrt(n) / 2;
    while (4 * band_start * band_start < n) {
        ++band_start;
    }

    ulong count = 0;
    for (ulong a = 1; a <= a_limit && a < band_start; ++a) {
        // 4a^2 < |D|: every b with -a < b <= a gives a reduced form with
        // c > a (Cohen, GTM 138, Lemma 5.3.4(2)).  b is defined mod 2a
        // (Cohen, p. 253), and b and b + 2a are both roots mod 4a, so the
        // window holds exactly half of the R(4a) roots.
        count += root_count[a] / 2;
    }

    // 4a^2 >= |D|: list the roots b mod 4a and keep those that give a
    // reduced form (a, b, c), c = (b^2 + |D|) / (4a): |b| <= a <= c, and
    // b >= 0 when |b| = a or a = c (Cohen, GTM 138, Definition 5.3.2).  The
    // roots mod 4a are the CRT combinations of the cached roots mod the prime
    // powers exactly dividing 4a (Cohen, Theorem 1.3.9); the band is factored
    // in chunks by a sieve over the odd primes.
    if (band_start <= a_limit) {
        const ulong width = a_limit - band_start + 1;
        ulong chunk = (width + kMaxBandChunks - 1) / kMaxBandChunks;
        if (chunk < kMinBandChunk) {
            chunk = width < kMinBandChunk ? width : kMinBandChunk;
        }
        std::vector<std::uint32_t> odd_rest(chunk);
        std::vector<unsigned char> factor_count(chunk);
        std::vector<std::uint32_t> factors(chunk * kMaxOddFactors);
        std::vector<ulong> band_roots(kMaxBandRoots);
        const std::vector<std::uint32_t>& primes = odd_roots.primes;
        for (ulong low = band_start; low <= a_limit; low += chunk) {
            const ulong high =
                    a_limit - low < chunk ? a_limit : low + chunk - 1;
            for (ulong a = low; a <= high; ++a) {
                const ulong j = a - low;
                odd_rest[j] = root_count[a] == 0
                                      ? 0
                                      : static_cast<std::uint32_t>(
                                                a >> ::flint_ctz(a));
                factor_count[j] = 0;
            }
            for (std::size_t index = 0;
                 index < primes.size() && primes[index] <= high; ++index) {
                const ulong p = primes[index];
                for (ulong a = (low + p - 1) / p * p; a <= high; a += p) {
                    const ulong j = a - low;
                    ulong rest = odd_rest[j];
                    if (rest == 0) {
                        continue;
                    }
                    std::uint32_t exponent = 0;
                    do {
                        rest /= p;
                        ++exponent;
                    } while (rest % p == 0);
                    odd_rest[j] = static_cast<std::uint32_t>(rest);
                    if (factor_count[j] == kMaxOddFactors) {
                        return false;
                    }
                    factors[j * kMaxOddFactors + factor_count[j]] =
                            (static_cast<std::uint32_t>(index)
                             << kExponentBits) |
                            exponent;
                    ++factor_count[j];
                }
            }

            for (ulong a = low; a <= high; ++a) {
                const ulong j = a - low;
                const ulong roots = root_count[a];
                if (roots == 0) {
                    continue;
                }
                if (odd_rest[j] != 1) {
                    return false;
                }
                const ulong four_a = 4 * a;
                const int two_exponent = 2 + static_cast<int>(::flint_ctz(a));
                const PrimePowerRoots& two_roots =
                        two_power_roots[two_exponent];
                if (!two_power_known[two_exponent]) {
                    return false;
                }
                // The idempotent for the modulus q of a component is
                // (4a / q) * ((4a / q)^-1 mod q), so the component with root
                // r adds (4a / q) * (r * (4a / q)^-1 mod q) < 4a.
                const ulong two_mask = (UWORD(1) << two_exponent) - 1;
                const ulong odd_part = four_a >> two_exponent;
                const ulong two_inverse =
                        ::n_invmod(odd_part & two_mask, two_mask + 1);
                ulong listed = 0;
                for (ulong i = 0; i < two_roots.count; ++i) {
                    band_roots[listed++] =
                            odd_part *
                            ((two_roots.roots[i] * two_inverse) & two_mask);
                }
                for (std::size_t f = 0; f < factor_count[j]; ++f) {
                    const std::uint32_t entry = factors[j * kMaxOddFactors + f];
                    const std::size_t index = entry >> kExponentBits;
                    const int exponent =
                            static_cast<int>(entry & kExponentMask);
                    std::uint32_t root = 0;
                    if (!odd_roots.root(root, index, exponent) ||
                        root == kNoRoot) {
                        return false;
                    }
                    if (root == 0) {
                        continue;
                    }
                    const ulong modulus =
                            exponent == 1
                                    ? primes[index]
                                    : ::n_pow(primes[index],
                                              static_cast<ulong>(exponent));
                    const ulong cofactor = four_a / modulus;
                    const ulong term =
                            cofactor *
                            ((root * ::n_invmod(cofactor % modulus, modulus)) %
                             modulus);
                    if (2 * listed > kMaxBandRoots) {
                        return false;
                    }
                    for (ulong i = 0; i < listed; ++i) {
                        const ulong x = band_roots[i];
                        band_roots[i] =
                                x >= four_a - term ? x - (four_a - term)
                                                   : x + term;
                        band_roots[listed + i] =
                                x >= term ? x - term : x + (four_a - term);
                    }
                    listed *= 2;
                }
                if (listed != roots) {
                    return false;
                }
                for (ulong i = 0; i < listed; ++i) {
                    // The representative of b in (-a, a]; the other residues
                    // mod 4a are b + 2a and lie outside it.
                    const ulong root = band_roots[i];
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
        }
    }

    out = count;
    return true;
}

}  // namespace silex::detail
