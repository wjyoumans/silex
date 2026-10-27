#pragma once

#include <flint/flint.h>
#include <flint/longlong.h>

#include <array>
#include <cstddef>

namespace silex::detail {

// Word-size kernel for x^e mod f, f monic of degree d over F_p for a prime p
// below the word size, the quantity FLINT's nmod_poly_powmod_x_ui_preinv
// computes.  It serves the batched root counts of the zeta BF sums (e = p,
// zeta.cpp) and the Rabin root split of the proven stable scan
// (e = (p - 1) / 2, saturation.cpp), for d at most kWordFrobeniusMaxDegree
// and p below word_frobenius_prime_limit(d).  Coefficients are reduced
// lazily, once per squaring, by Barrett reduction (Brent and Zimmermann,
// Modern Computer Arithmetic, 2010, Section 2.4.1).
//
// Word-size overflow.  Every stored residue is at most M = p - 1.  The
// unreduced sums the kernel and its callers form are:
//   - square[s] = sum of power[i] power[j] over i + j = s, 0 <= i, j < d:
//     at most min(s + 1, 2d - 1 - s) <= d products, each at most M^2.
//     The kernel adds power[i]^2 and 2 power[i] power[j] for i < j
//     (word_square); these nonnegative terms sum to the same integer,
//     so the sum and this bound are unchanged, and every partial sum is
//     at most the sum;
//   - for j < d, square[j] + sum over d <= k <= 2d - 2 of
//     reduce(square[k]) reduction[k - d][j]: at most (j + 1) + (d - 1)
//     <= 2d - 1 products, each at most M^2 (j = d - 1 attains 2d - 1);
//   - in the sparse fold (word_x_power_mod_sparse), each coefficient of
//     the square plus its fold products: at most 2d - 1 products, each at
//     most M^2 (counted there);
//   - a + top reduction[0][j] in the reduction table and the
//     multiplication by x: at most M + M^2 <= 2 M^2 (M >= 1);
//   - single products in word_mulmod and the Euclidean steps of the
//     callers: M^2.
// So every intermediate value is at most (2d - 1) M^2 once d >= 2, and the
// kernel is exact when (2d - 1) M^2 <= 2^64 - 1, which word_reduce then
// reduces correctly (its quotient estimate holds for every 64-bit value).
// The exponent e only selects the squarings and multiplications by x, so
// the bound does not depend on it.  p < L gives M <= L - 2, and the limits
// below satisfy the condition with M = L - 1: L = 2^30 for d <= 8
// (2d - 1 <= 15 < 16 = 2^64 / 2^60), and L = 2^29 for 9 <= d <= 20
// (2d - 1 <= 39 < 64 = 2^64 / 2^58; at d = 20, 39 (2^29 - 1)^2 < 39 * 2^58
// < 2^63.3).  The bound is checked at compile time for every instantiated
// degree (word_frobenius_sums_fit).  At d = 9 the limit 2^30 would fail:
// 17 (2^30 - 2)^2 > 2^64.  Primes at or above the limit, and higher
// degrees, are left to the callers' FLINT paths.
constexpr slong kWordFrobeniusMaxDegree = WORD(20);

constexpr ulong word_frobenius_prime_limit(slong degree) noexcept {
    return degree <= WORD(8) ? UWORD(1) << 30 : UWORD(1) << 29;
}

// Whether (2 degree - 1) (L - 1)^2 <= 2^64 - 1 for the prime limit L of
// `degree`, the overflow condition above, with M bounded by L - 1.
constexpr bool word_frobenius_sums_fit(slong degree) noexcept {
    const ulong residue = word_frobenius_prime_limit(degree) - 1;
    return degree >= 2 && residue <= UWORD(0xffffffff) &&
           residue * residue <=
                   UWORD_MAX / static_cast<ulong>(2 * degree - 1);
}

static_assert(word_frobenius_sums_fit(kWordFrobeniusMaxDegree));

static_assert(FLINT_BITS == 64);

struct WordModulus {
    ulong n = 0;
    // floor((2^64 - 1) / n): the quotient estimate is at most one short,
    // since inverse >= (2^64 - n) / n gives value * inverse / 2^64 >
    // value / n - 1 for value < 2^64.
    ulong inverse = 0;
};

inline WordModulus word_modulus(ulong n) noexcept {
    return WordModulus{n, UWORD_MAX / n};
}

inline ulong word_reduce(ulong value, const WordModulus& modulus) noexcept {
    ulong quotient = 0;
    ulong low = 0;
    umul_ppmm(quotient, low, value, modulus.inverse);
    (void)low;
    const ulong remainder = value - quotient * modulus.n;
    return remainder >= modulus.n ? remainder - modulus.n : remainder;
}

inline ulong word_mulmod(ulong a, ulong b, const WordModulus& modulus)
        noexcept {
    return word_reduce(a * b, modulus);
}

// How the word kernel squares x^(2^k ...) mod f: `symmetric` forms each
// cross product power[i] power[j], i < j, once and doubles it,
// d (d + 1) / 2 products; `full` is the plain i, j double loop, d^2
// products, kept as the reference the tests compare against.
enum class WordSquare { symmetric, full };

// square[s] = sum of power[i] power[j] over i + j = s, 0 <= i, j < d, as an
// unreduced integer.  Both orders add the same nonnegative products, so the
// result is the same integer, and every partial sum is at most it.
template <slong Degree, WordSquare Square>
void word_square(
        std::array<ulong, 2 * static_cast<std::size_t>(Degree) - 1>& square,
        const std::array<ulong, static_cast<std::size_t>(Degree)>& power)
        noexcept {
    constexpr std::size_t kDegree = static_cast<std::size_t>(Degree);
    square.fill(0);
    if constexpr (Square == WordSquare::full) {
        for (std::size_t i = 0; i < kDegree; ++i) {
            for (std::size_t j = 0; j < kDegree; ++j) {
                square[i + j] += power[i] * power[j];
            }
        }
    } else {
        for (std::size_t i = 0; i < kDegree; ++i) {
            const ulong twice = power[i] + power[i];
            square[2 * i] += power[i] * power[i];
            for (std::size_t j = i + 1; j < kDegree; ++j) {
                square[i + j] += twice * power[j];
            }
        }
    }
}

// How the word kernel reduces each square modulo f.  `dense` folds every
// coefficient of degree d to 2d - 2 through a table of x^(d + k) mod f,
// d (d - 1) products per squaring after the d (d - 2) products that build
// the table; `sparse` divides by f from the top coefficient down, adding
// one product per nonzero low coefficient of f, (d - 1) t products for t
// nonzero coefficients among f_0, ..., f_(d - 1), with no table; and
// `automatic` takes `sparse` when t <= word_sparse_fold_limit(d) and
// `dense` otherwise.  Both compute the remainder of division by the monic
// f, which is unique, from the same exact integers modulo p, and reduce
// every coefficient to [0, p), so they return the same residues; the tests
// compare them.
enum class WordFold { automatic, dense, sparse };

// The largest number t of nonzero low coefficients of f for which
// `automatic` takes the sparse fold.  The sparse fold forms at most as many
// products as the dense one for every t <= d and builds no table, but its
// top-down order chains each reduction on the folds above it, which costs
// time where the dense fold's independent sums are cheap.  The limits are
// set below the crossover measured for every degree 2 to 20, with primes
// below 20000 and e = p as in the zeta BF sums: up to degree 7 only the
// leading-term-only f (t = 0) gains; from degree 8 to 18 the sparse fold
// is faster up to at least t = d / 2; at degrees 19 and 20 it is faster
// for every t.
constexpr std::size_t word_sparse_fold_limit(slong degree) noexcept {
    if (degree <= WORD(7)) {
        return 0;
    }
    if (degree <= WORD(18)) {
        return static_cast<std::size_t>(degree) / 2;
    }
    return static_cast<std::size_t>(degree);
}

// x^exponent mod f by left-to-right binary powering from x, with the
// dense fold: negated[j] = -f_j mod p.
template <slong Degree, WordSquare Square>
void word_x_power_mod_dense(
        std::array<ulong, static_cast<std::size_t>(Degree)>& power,
        const std::array<ulong, static_cast<std::size_t>(Degree)>& negated,
        ulong exponent,
        const WordModulus& modulus) noexcept {
    constexpr std::size_t kDegree = static_cast<std::size_t>(Degree);

    // reduction[k][j]: coefficient of x^j in x^(Degree + k) mod f.
    std::array<std::array<ulong, kDegree>, kDegree - 1> reduction{};
    reduction[0] = negated;
    for (std::size_t k = 1; k + 1 < kDegree; ++k) {
        const ulong top = reduction[k - 1][kDegree - 1];
        reduction[k][0] = word_mulmod(top, reduction[0][0], modulus);
        for (std::size_t j = 1; j < kDegree; ++j) {
            reduction[k][j] = word_reduce(
                    reduction[k - 1][j - 1] + top * reduction[0][j],
                    modulus);
        }
    }

    power.fill(0);
    power[1] = 1;
    for (slong bit = static_cast<slong>(FLINT_BIT_COUNT(exponent)) - 2;
         bit >= 0; --bit) {
        std::array<ulong, 2 * kDegree - 1> square;
        word_square<Degree, Square>(square, power);
        for (std::size_t k = kDegree; k < 2 * kDegree - 1; ++k) {
            square[k] = word_reduce(square[k], modulus);
        }
        for (std::size_t j = 0; j < kDegree; ++j) {
            ulong sum = square[j];
            for (std::size_t k = kDegree; k < 2 * kDegree - 1; ++k) {
                sum += square[k] * reduction[k - kDegree][j];
            }
            power[j] = word_reduce(sum, modulus);
        }
        if (((exponent >> bit) & UWORD(1)) != 0) {
            const ulong top = power[kDegree - 1];
            for (std::size_t j = kDegree - 1; j >= 1; --j) {
                power[j] = word_reduce(
                        power[j - 1] + top * reduction[0][j], modulus);
            }
            power[0] = word_mulmod(top, reduction[0][0], modulus);
        }
    }
}

// x^exponent mod f by left-to-right binary powering from x, with the
// sparse fold: negated[j] = -f_j mod p, nonzero exactly at the `count`
// increasing indices in `terms`.  The fold is division with remainder by
// the monic f (von zur Gathen and Gerhard, Modern Computer Algebra, 3rd
// ed., 2013, Algorithm 2.5), from coefficient 2d - 2 down to d: the
// coefficient c of x^k, reduced, is replaced by c x^(k - d) (x^d - f),
// that is c negated[i] added at k - d + i for each term i.  Every term
// lands below k, so each coefficient has received all its additions when
// it is reduced.  The unreduced sums stay within the dense bound of the
// overflow comment above: coefficient j < d holds its at most min(j + 1,
// d) squaring products plus one product for each term i with
// max(0, j - d + 2) <= i <= j, at most j + 1 for j <= d - 2 and d - 1 for
// j = d - 1, so at most 2d - 1 products of residues; coefficient
// j >= d holds at most (2d - 1 - j) + (2d - 2 - j) <= 2d - 3.  The
// multiplication by x adds top negated[i] to a reduced residue, at most
// M + M^2, as in the dense path.
template <slong Degree, WordSquare Square>
void word_x_power_mod_sparse(
        std::array<ulong, static_cast<std::size_t>(Degree)>& power,
        const std::array<ulong, static_cast<std::size_t>(Degree)>& negated,
        const std::array<std::size_t, static_cast<std::size_t>(Degree)>&
                terms,
        std::size_t count,
        ulong exponent,
        const WordModulus& modulus) noexcept {
    constexpr std::size_t kDegree = static_cast<std::size_t>(Degree);

    power.fill(0);
    power[1] = 1;
    for (slong bit = static_cast<slong>(FLINT_BIT_COUNT(exponent)) - 2;
         bit >= 0; --bit) {
        std::array<ulong, 2 * kDegree - 1> square;
        word_square<Degree, Square>(square, power);
        for (std::size_t k = 2 * kDegree - 2; k >= kDegree; --k) {
            const ulong top = word_reduce(square[k], modulus);
            const std::size_t shift = k - kDegree;
            for (std::size_t t = 0; t < count; ++t) {
                const std::size_t i = terms[t];
                square[shift + i] += top * negated[i];
            }
        }
        for (std::size_t j = 0; j < kDegree; ++j) {
            power[j] = word_reduce(square[j], modulus);
        }
        if (((exponent >> bit) & UWORD(1)) != 0) {
            const ulong top = power[kDegree - 1];
            for (std::size_t j = kDegree - 1; j >= 1; --j) {
                power[j] = power[j - 1];
            }
            power[0] = 0;
            for (std::size_t t = 0; t < count; ++t) {
                const std::size_t i = terms[t];
                power[i] = word_reduce(power[i] + top * negated[i], modulus);
            }
        }
    }
}

// `power` receives x^exponent mod f, exponent >= 1, for the monic `f` of
// degree Degree with coefficients reduced modulo the prime modulus.n.
template <slong Degree,
          WordSquare Square,
          WordFold Fold = WordFold::automatic>
void word_x_power_mod(std::array<ulong, static_cast<std::size_t>(Degree)>&
                              power,
                      const ulong* f,
                      ulong exponent,
                      const WordModulus& modulus) noexcept {
    static_assert(Degree >= 2 && Degree <= kWordFrobeniusMaxDegree &&
                  word_frobenius_sums_fit(Degree));
    constexpr std::size_t kDegree = static_cast<std::size_t>(Degree);
    const ulong n = modulus.n;

    // negated[j] = -f_j mod p, the coefficient of x^j in x^Degree mod f.
    std::array<ulong, kDegree> negated;
    std::array<std::size_t, kDegree> terms{};
    std::size_t count = 0;
    for (std::size_t j = 0; j < kDegree; ++j) {
        negated[j] = f[j] == 0 ? 0 : n - f[j];
        if (negated[j] != 0) {
            terms[count++] = j;
        }
    }
    const bool sparse =
            Fold == WordFold::sparse ||
            (Fold == WordFold::automatic &&
             count <= word_sparse_fold_limit(Degree));
    if (sparse) {
        word_x_power_mod_sparse<Degree, Square>(power, negated, terms, count,
                                                exponent, modulus);
    } else {
        word_x_power_mod_dense<Degree, Square>(power, negated, exponent,
                                               modulus);
    }
}

template <slong Degree>
bool word_x_power_mod_from(ulong* power,
                           const ulong* f,
                           slong degree,
                           ulong exponent,
                           const WordModulus& modulus) noexcept {
    if constexpr (Degree > kWordFrobeniusMaxDegree) {
        return false;
    } else {
        if (degree != Degree) {
            return word_x_power_mod_from<Degree + 1>(power, f, degree,
                                                     exponent, modulus);
        }
        std::array<ulong, static_cast<std::size_t>(Degree)> result;
        word_x_power_mod<Degree, WordSquare::symmetric>(result, f, exponent,
                                                        modulus);
        for (std::size_t j = 0; j < static_cast<std::size_t>(Degree); ++j) {
            power[j] = result[j];
        }
        return true;
    }
}

// `power` (`degree` entries) receives x^exponent mod f, for `f` with the
// `degree` + 1 coefficients reduced modulo the prime p, when the
// word kernel covers the case: 2 <= degree <= kWordFrobeniusMaxDegree,
// p < word_frobenius_prime_limit(degree), exponent >= 1 and f[degree] == 1.
// Returns false, leaving `power` unchanged, otherwise: the kernel reads f as
// monic and never looks at f[degree], so a non-monic f is declined.
inline bool word_x_power_mod_dispatch(ulong* power,
                                      const ulong* f,
                                      slong degree,
                                      ulong exponent,
                                      ulong p) noexcept {
    if (degree < 2 || degree > kWordFrobeniusMaxDegree || exponent == 0 ||
        p < 2 || p >= word_frobenius_prime_limit(degree) ||
        f[degree] != 1) {
        return false;
    }
    return word_x_power_mod_from<2>(power, f, degree, exponent,
                                    word_modulus(p));
}

}  // namespace silex::detail
