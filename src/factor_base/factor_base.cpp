#include <silex/factor_base.hpp>

#include "factor_base_internal.hpp"
#include "../prime_ideal/prime_ideal_internal.hpp"

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arf.hpp>
#include <silex/signature.hpp>

#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>

#include <algorithm>
#include <utility>

namespace silex {
namespace {

void floor_sqrt_div_ui(flint::Fmpz& out,
                       flint::FmpzConstRef value,
                       ulong denominator) noexcept {
    flint::Fmpz quotient;
    fmpz_fdiv_q_ui(quotient.raw(), value.raw(), denominator);
    fmpz_sqrt(out.raw(), quotient.raw());
}

void ceil_sqrt(flint::Fmpz& out, flint::FmpzConstRef value) noexcept {
    flint::Fmpz remainder;
    fmpz_sqrtrem(out.raw(), remainder.raw(), value.raw());
    if (fmpz_is_zero(remainder.raw()) == 0) {
        fmpz_add_ui(out.raw(), out.raw(), 1);
    }
}

// Exact integer form of Minkowski's bound with 2^r2 in place of (4/pi)^r2:
// ceil(n! * ceil(sqrt|d|) * 2^r2 / n^n).  Since 4/pi < 2 it is never below
// the bounds computed below, and it is the fallback if their Arb evaluation
// does not produce a finite enclosure.
void generic_minkowski_bound(flint::Fmpz& out,
                             flint::FmpzConstRef abs_discriminant,
                             slong degree,
                             slong complex_pairs) noexcept {
    flint::Fmpz numerator;
    flint::Fmpz denominator;
    flint::Fmpz sqrt_discriminant;

    ceil_sqrt(sqrt_discriminant, abs_discriminant);
    fmpz_fac_ui(numerator.raw(), static_cast<ulong>(degree));
    fmpz_mul(numerator.raw(), numerator.raw(), sqrt_discriminant.raw());
    fmpz_mul_2exp(numerator.raw(), numerator.raw(),
                  static_cast<ulong>(complex_pairs));
    fmpz_set_ui(denominator.raw(), static_cast<ulong>(degree));
    fmpz_pow_ui(denominator.raw(), denominator.raw(),
                static_cast<ulong>(degree));
    fmpz_cdiv_q(out.raw(), numerator.raw(), denominator.raw());
    if (fmpz_is_zero(out.raw()) != 0) {
        fmpz_one(out.raw());
    }
}

constexpr slong kGenerationBoundPrecision = 128;
constexpr slong kZimmertMaxDegree = 20;

// Parameter gamma of Zimmert's Satz 2, in hundredths, indexed by
// [n - 3][r2] for 3 <= n <= 20.  Source: H. Zimmert, "Ideale kleiner Norm in
// Idealklassen und eine Regulatorabschaetzung", Invent. Math. 62 (1981),
// 367-380.  Satz 2 (p. 372) holds for every gamma > alpha > 0, so gamma
// affects only the size of the bound, never its validity.  Each entry
// maximises the right side of Satz 2 over gamma in steps of 1/100, with alpha
// from Bemerkung 1 (p. 373), as Zimmert did for Tabelle 1 (p. 368, explained
// on pp. 373-374).  For every signature in Tabelle 1 the entry equals
// Zimmert's printed gamma.
constexpr unsigned char kZimmertGammaHundredths[18][11] = {
    {156, 184},                                           // n = 3
    {118, 134, 154},                                      // n = 4
    {96, 107, 120},                                       // n = 5
    {83, 90, 99, 109},                                    // n = 6
    {73, 79, 85, 92},                                     // n = 7
    {66, 70, 75, 81, 87},                                 // n = 8
    {60, 64, 68, 72, 77},                                 // n = 9
    {56, 59, 62, 65, 69, 74},                             // n = 10
    {52, 55, 57, 60, 63, 67},                             // n = 11
    {49, 51, 54, 56, 59, 62, 65},                         // n = 12
    {47, 49, 50, 52, 55, 57, 60},                         // n = 13
    {45, 46, 48, 50, 51, 54, 56, 59},                     // n = 14
    {43, 44, 45, 47, 49, 50, 53, 55},                     // n = 15
    {41, 42, 43, 45, 46, 48, 50, 52, 54},                 // n = 16
    {39, 40, 42, 43, 44, 46, 47, 49, 51},                 // n = 17
    {38, 39, 40, 41, 42, 43, 45, 46, 48, 50},             // n = 18
    {37, 38, 38, 39, 41, 42, 43, 44, 46, 47},             // n = 19
    {36, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46},         // n = 20
};

// Sets out to the ceiling of the upper endpoint of value.  Fails if the
// enclosure is not finite.
bool ceil_upper_endpoint(flint::Fmpz& out, const flint::Arb& value) noexcept {
    if (arb_is_finite(value.raw()) == 0) {
        return false;
    }
    flint::Arf upper;
    arb_get_ubound_arf(upper.raw(), value.raw(), kGenerationBoundPrecision);
    if (arf_is_finite(upper.raw()) == 0) {
        return false;
    }
    arf_get_fmpz(out.raw(), upper.raw(), ARF_RND_CEIL);
    return true;
}

// Encloses the right side of Zimmert 1981, Satz 2 (p. 372):
//
//   log(d^(1/2) / N(a)) >=
//       r1 (-psi((1+g)/2) - log Gamma(1/2+g) + log Gamma(1+g) + (1/2) log pi)
//     + r2 (-2 psi(1+g) + 2 log 2 + log(1/2+g) + log pi)
//     - 2/(g-a) - log[(1+1/a) (1+1/g)^(-2) (1+1/(2g-a))^(-1)],
//
// valid for every g > a > 0, where every ideal class contains an integral
// ideal a with this norm bound.  Here g = gamma_hundredths/100 and
// a = g - g(g+1)/sqrt(1+3g+3g^2), Zimmert's choice of alpha in Bemerkung 1
// (p. 373); 0 < a < g for every g > 0.  Ball arithmetic encloses the value at
// the exact g and a, so the lower endpoint is a valid lower bound.
bool zimmert_log_lower_bound(flint::Arb& out,
                             slong r1,
                             slong r2,
                             ulong gamma_hundredths) noexcept {
    const slong prec = kGenerationBoundPrecision;
    flint::Arb g;
    flint::Arb alpha;
    flint::Arb x;
    flint::Arb y;
    flint::Arb log_pi;
    flint::Arb real_term;
    flint::Arb complex_term;
    flint::Arb total;

    arb_set_ui(g.raw(), gamma_hundredths);
    arb_div_ui(g.raw(), g.raw(), 100, prec);

    // alpha = g - g(g+1)/sqrt(1+3g+3g^2)
    arb_mul_ui(x.raw(), g.raw(), 3, prec);
    arb_add_ui(y.raw(), g.raw(), 1, prec);
    arb_mul(x.raw(), x.raw(), y.raw(), prec);
    arb_add_ui(x.raw(), x.raw(), 1, prec);
    arb_sqrt(x.raw(), x.raw(), prec);
    arb_mul(y.raw(), y.raw(), g.raw(), prec);
    arb_div(y.raw(), y.raw(), x.raw(), prec);
    arb_sub(alpha.raw(), g.raw(), y.raw(), prec);
    arb_sub(x.raw(), g.raw(), alpha.raw(), prec);
    if (arb_is_positive(alpha.raw()) == 0 || arb_is_positive(x.raw()) == 0) {
        return false;
    }

    arb_const_pi(log_pi.raw(), prec);
    arb_log(log_pi.raw(), log_pi.raw(), prec);

    // Real places: -psi((1+g)/2) - log Gamma(1/2+g) + log Gamma(1+g)
    //              + (1/2) log pi
    arb_add_ui(x.raw(), g.raw(), 1, prec);
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_neg(real_term.raw(), y.raw());
    arb_one(x.raw());
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_add(x.raw(), x.raw(), g.raw(), prec);
    arb_lgamma(y.raw(), x.raw(), prec);
    arb_sub(real_term.raw(), real_term.raw(), y.raw(), prec);
    arb_add_ui(x.raw(), g.raw(), 1, prec);
    arb_lgamma(y.raw(), x.raw(), prec);
    arb_add(real_term.raw(), real_term.raw(), y.raw(), prec);
    arb_mul_2exp_si(y.raw(), log_pi.raw(), -1);
    arb_add(real_term.raw(), real_term.raw(), y.raw(), prec);

    // Complex places: -2 psi(1+g) + 2 log 2 + log(1/2+g) + log pi
    arb_add_ui(x.raw(), g.raw(), 1, prec);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_mul_2exp_si(complex_term.raw(), y.raw(), 1);
    arb_neg(complex_term.raw(), complex_term.raw());
    arb_const_log2(y.raw(), prec);
    arb_mul_2exp_si(y.raw(), y.raw(), 1);
    arb_add(complex_term.raw(), complex_term.raw(), y.raw(), prec);
    arb_one(x.raw());
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_add(x.raw(), x.raw(), g.raw(), prec);
    arb_log(y.raw(), x.raw(), prec);
    arb_add(complex_term.raw(), complex_term.raw(), y.raw(), prec);
    arb_add(complex_term.raw(), complex_term.raw(), log_pi.raw(), prec);

    arb_mul_si(total.raw(), real_term.raw(), r1, prec);
    arb_mul_si(y.raw(), complex_term.raw(), r2, prec);
    arb_add(total.raw(), total.raw(), y.raw(), prec);

    // - 2/(g-a)
    arb_sub(x.raw(), g.raw(), alpha.raw(), prec);
    arb_ui_div(y.raw(), 2, x.raw(), prec);
    arb_sub(total.raw(), total.raw(), y.raw(), prec);

    // - log(1+1/a) + 2 log(1+1/g) + log(1+1/(2g-a))
    arb_inv(x.raw(), alpha.raw(), prec);
    arb_log1p(y.raw(), x.raw(), prec);
    arb_sub(total.raw(), total.raw(), y.raw(), prec);
    arb_inv(x.raw(), g.raw(), prec);
    arb_log1p(y.raw(), x.raw(), prec);
    arb_mul_2exp_si(y.raw(), y.raw(), 1);
    arb_add(total.raw(), total.raw(), y.raw(), prec);
    arb_mul_2exp_si(x.raw(), g.raw(), 1);
    arb_sub(x.raw(), x.raw(), alpha.raw(), prec);
    arb_inv(x.raw(), x.raw(), prec);
    arb_log1p(y.raw(), x.raw(), prec);
    arb_add(total.raw(), total.raw(), y.raw(), prec);

    if (arb_is_finite(total.raw()) == 0) {
        return false;
    }
    arb_swap(out.raw(), total.raw());
    return true;
}

// Zimmert 1981, Satz 2: every ideal class contains an integral ideal of norm
// at most sqrt|d| * exp(-L), with L the right side above.  Sets out to the
// ceiling of a rigorous upper bound for sqrt|d| * exp(-L).
bool zimmert_bound(flint::Fmpz& out,
                   flint::FmpzConstRef abs_discriminant,
                   slong degree,
                   slong real_places,
                   slong complex_pairs) noexcept {
    if (degree < 3 || degree > kZimmertMaxDegree || real_places < 0 ||
        complex_pairs < 0 || real_places + 2 * complex_pairs != degree) {
        return false;
    }
    const ulong gamma =
        kZimmertGammaHundredths[degree - 3][complex_pairs];
    flint::Arb log_lower;
    if (gamma == 0 ||
        !zimmert_log_lower_bound(log_lower, real_places, complex_pairs,
                                 gamma)) {
        return false;
    }
    flint::Arb value;
    flint::Arb root;
    arb_neg(log_lower.raw(), log_lower.raw());
    arb_exp(value.raw(), log_lower.raw(), kGenerationBoundPrecision);
    arb_set_fmpz(root.raw(), abs_discriminant.raw());
    arb_sqrt(root.raw(), root.raw(), kGenerationBoundPrecision);
    arb_mul(value.raw(), value.raw(), root.raw(), kGenerationBoundPrecision);
    return ceil_upper_endpoint(out, value);
}

// Parameter gamma of Zimmert's Korollar 1 (twin classes), in hundredths,
// indexed by [n - 3][r2] for 3 <= n <= 20.  Korollar 1 (p. 376) holds for
// every gamma > 0, and the finite-m inequality evaluated below for every
// gamma > 0 that is not an integer, so gamma affects only the size of the
// bound, never its validity.  Each entry maximises
// zimmert_twin_class_log_lower_bound over 0 < gamma < 3 in steps of 1/100,
// integers excluded (test/zimmert_twin_class_oracle.py --gammas).
constexpr unsigned char kZimmertTwinClassGammaHundredths[18][11] = {
    {153, 182},                                           // n = 3
    {116, 133, 153},                                      // n = 4
    {95, 106, 119},                                       // n = 5
    {82, 89, 98, 109},                                    // n = 6
    {72, 78, 85, 92},                                     // n = 7
    {65, 70, 75, 81, 87},                                 // n = 8
    {60, 63, 67, 72, 77},                                 // n = 9
    {56, 58, 62, 65, 69, 74},                             // n = 10
    {52, 54, 57, 60, 63, 67},                             // n = 11
    {49, 51, 53, 56, 59, 62, 65},                         // n = 12
    {46, 48, 50, 52, 55, 57, 60},                         // n = 13
    {44, 46, 48, 49, 51, 54, 56, 59},                     // n = 14
    {42, 44, 45, 47, 49, 50, 53, 55},                     // n = 15
    {41, 42, 43, 45, 46, 48, 50, 52, 54},                 // n = 16
    {39, 40, 41, 43, 44, 45, 47, 49, 51},                 // n = 17
    {38, 39, 40, 41, 42, 43, 45, 46, 48, 50},             // n = 18
    {36, 37, 38, 39, 40, 42, 43, 44, 46, 47},             // n = 19
    {35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 47},         // n = 20
};

// Number m of factors in Zimmert's Q_a and Q_b (proof of Satz 4, p. 379).
// Every m >= 1 gives a valid inequality; G_{m,j} increases to F_{a,b} as m
// grows, and at m = 128 it is within 1e-4 of F for the tabulated gammas.
constexpr slong kZimmertTwinClassTerms = 128;

// Sets out to psi(x).  For x not positive, uses psi(x) = psi(x + k) -
// sum_{i=0}^{k-1} 1/(x + i) with k = floor(-x) + 1, so arb_digamma is only
// evaluated at a positive argument.  Fails if x + i contains zero for some i
// or the result is not finite.
bool digamma_by_shift(flint::Arb& out,
                      const flint::Arb& x,
                      slong prec) noexcept {
    flint::Arb shifted;
    flint::Arb term;
    flint::Arb correction;
    arb_set(shifted.raw(), x.raw());
    arb_zero(correction.raw());
    if (arb_is_positive(x.raw()) == 0) {
        flint::Fmpz k;
        arf_get_fmpz(k.raw(), arb_midref(x.raw()), ARF_RND_FLOOR);
        fmpz_neg(k.raw(), k.raw());
        fmpz_add_ui(k.raw(), k.raw(), 1);
        if (fmpz_sgn(k.raw()) <= 0 || fmpz_cmp_ui(k.raw(), 64) > 0) {
            return false;
        }
        const ulong steps = fmpz_get_ui(k.raw());
        for (ulong i = 0; i < steps; ++i) {
            if (arb_contains_zero(shifted.raw()) != 0) {
                return false;
            }
            arb_inv(term.raw(), shifted.raw(), prec);
            arb_add(correction.raw(), correction.raw(), term.raw(), prec);
            arb_add_ui(shifted.raw(), shifted.raw(), 1, prec);
        }
        if (arb_is_positive(shifted.raw()) == 0) {
            return false;
        }
    }
    arb_digamma(term.raw(), shifted.raw(), prec);
    arb_sub(out.raw(), term.raw(), correction.raw(), prec);
    return arb_is_finite(out.raw()) != 0;
}

// Sets out to (scale/D) (psi(y + 1/2) - psi(y)) - 1/u - 1/v with y =
// (numerator)/D, where numerator, u and v are given as balls; v may be
// omitted (null).  Every argument of psi here is positive.
bool zimmert_twin_class_term(flint::Arb& out,
                             const flint::Arb& numerator,
                             const flint::Arb& d,
                             ulong scale,
                             const flint::Arb& u,
                             const flint::Arb* v,
                             slong prec) noexcept {
    flint::Arb y;
    flint::Arb z;
    flint::Arb value;
    arb_div(y.raw(), numerator.raw(), d.raw(), prec);
    arb_one(z.raw());
    arb_mul_2exp_si(z.raw(), z.raw(), -1);
    arb_add(z.raw(), z.raw(), y.raw(), prec);
    if (arb_is_positive(y.raw()) == 0) {
        return false;
    }
    arb_digamma(value.raw(), z.raw(), prec);
    arb_digamma(z.raw(), y.raw(), prec);
    arb_sub(value.raw(), value.raw(), z.raw(), prec);
    arb_mul_ui(value.raw(), value.raw(), scale, prec);
    arb_div(value.raw(), value.raw(), d.raw(), prec);
    if (arb_contains_zero(u.raw()) != 0) {
        return false;
    }
    arb_inv(z.raw(), u.raw(), prec);
    arb_sub(value.raw(), value.raw(), z.raw(), prec);
    if (v != nullptr) {
        if (arb_contains_zero(v->raw()) != 0) {
            return false;
        }
        arb_inv(z.raw(), v->raw(), prec);
        arb_sub(value.raw(), value.raw(), z.raw(), prec);
    }
    arb_swap(out.raw(), value.raw());
    return true;
}

// Encloses (G_{m,j}(g) + n log pi) / 2, where G_{m,j} is the left side of
// Zimmert's Hilfssatz (p. 378) for the function H = -Q_a^a Q_b^b U V of the
// proof of Satz 4 (p. 379), with a = r1 + r2, b = r2, g = gamma/100,
// m = kZimmertTwinClassTerms and j = ceil(1/2 + g), the smallest natural
// number with j >= 1/2 + g.  For every such m and j the proof establishes
// G_{m,j}(g) <= 2 log A + f'/f(1+g) + g'/g(1+g), and Satz 4 is the limit
// m -> infinity; so, as in the derivation of Korollar 1 (p. 376), for every
// ideal class K there are integral ideals a in K and b in D K^(-1) (D the
// class of the different) with log(d / (N(a) N(b))) >= G_{m,j}(g) + n log pi.
// With D = 2(1 + 2g), y_l = (2l - 1 + g)/D, y'_l = (2l + g)/D and
// delta(y) = psi(y + 1/2) - psi(y), the log-derivatives of Zimmert's
// definitions give
//
//   G = a [sum_{l=1}^{m} ((4/D) delta(y_l) - 1/(2l-2-g) - 1/(2l-1+g))
//          + sum_{l=m+1}^{m+j} ((2/D) delta(y_l) - 1/(2l-2-g))]
//     + b [sum_{l=1}^{m} ((4/D) delta(y'_l) - 1/(2l-1-g) - 1/(2l+g))
//          + sum_{l=m+1}^{m+j} ((2/D) delta(y'_l) - 1/(2l-1-g))]
//     + (4/D) [psi(g/D) + psi((1+g)/D) - psi((1+3g)/D) - psi((2+3g)/D)]
//     - (a/2) [psi(1/2 + g/2) + psi(-g/2)] - (b/2) [psi(1 + g/2) + psi(1/2 - g/2)].
//
// This is the inequality the proof establishes for finite m, not the limit
// F_{a,b} of p. 376; no tail estimate is used.  Fails for an integer g
// (Zimmert's terms are singular there) or g outside (0, 3).
bool zimmert_twin_class_log_lower_bound(flint::Arb& out,
                                        slong r1,
                                        slong r2,
                                        ulong gamma_hundredths) noexcept {
    const slong prec = kGenerationBoundPrecision;
    if (r1 < 0 || r2 < 0 || r1 + r2 < 1 || gamma_hundredths == 0 ||
        gamma_hundredths >= 300 || gamma_hundredths % 100 == 0) {
        return false;
    }
    const slong a = r1 + r2;
    const slong b = r2;
    const slong m = kZimmertTwinClassTerms;
    const slong j = static_cast<slong>((gamma_hundredths + 149) / 100);

    flint::Arb g;
    flint::Arb d;
    flint::Arb numerator;
    flint::Arb u;
    flint::Arb v;
    flint::Arb term;
    flint::Arb real_sum;
    flint::Arb complex_sum;
    flint::Arb x;
    flint::Arb y;
    flint::Arb total;

    arb_set_ui(g.raw(), gamma_hundredths);
    arb_div_ui(g.raw(), g.raw(), 100, prec);
    arb_mul_2exp_si(d.raw(), g.raw(), 1);
    arb_add_ui(d.raw(), d.raw(), 1, prec);
    arb_mul_2exp_si(d.raw(), d.raw(), 1);

    arb_zero(real_sum.raw());
    arb_zero(complex_sum.raw());
    for (slong l = 1; l <= m + j; ++l) {
        const bool full = l <= m;
        const ulong scale = full ? 4 : 2;
        // Real places: y_l = (2l - 1 + g)/D, u = 2l - 2 - g, v = 2l - 1 + g.
        arb_add_si(numerator.raw(), g.raw(), 2 * l - 1, prec);
        arb_neg(u.raw(), g.raw());
        arb_add_si(u.raw(), u.raw(), 2 * l - 2, prec);
        if (!zimmert_twin_class_term(term, numerator, d, scale, u,
                                     full ? &numerator : nullptr, prec)) {
            return false;
        }
        arb_add(real_sum.raw(), real_sum.raw(), term.raw(), prec);
        if (b == 0) {
            continue;
        }
        // Complex places: y'_l = (2l + g)/D, u = 2l - 1 - g, v = 2l + g.
        arb_add_si(numerator.raw(), g.raw(), 2 * l, prec);
        arb_neg(u.raw(), g.raw());
        arb_add_si(u.raw(), u.raw(), 2 * l - 1, prec);
        if (!zimmert_twin_class_term(term, numerator, d, scale, u,
                                     full ? &numerator : nullptr, prec)) {
            return false;
        }
        arb_add(complex_sum.raw(), complex_sum.raw(), term.raw(), prec);
    }

    // Real-place constant: psi(1/2 + g/2) + psi(-g/2), weight -a/2.
    flint::Arb constant;
    arb_mul_2exp_si(x.raw(), g.raw(), -1);
    arb_neg(x.raw(), x.raw());
    if (!digamma_by_shift(constant, x, prec)) {
        return false;
    }
    arb_one(x.raw());
    arb_add(x.raw(), x.raw(), g.raw(), prec);
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_add(constant.raw(), constant.raw(), y.raw(), prec);
    arb_mul_2exp_si(constant.raw(), constant.raw(), -1);
    arb_sub(real_sum.raw(), real_sum.raw(), constant.raw(), prec);

    // Complex-place constant: psi(1 + g/2) + psi(1/2 - g/2), weight -b/2.
    if (b != 0) {
        arb_mul_2exp_si(x.raw(), g.raw(), -1);
        arb_add_ui(y.raw(), x.raw(), 1, prec);
        arb_digamma(constant.raw(), y.raw(), prec);
        arb_neg(x.raw(), x.raw());
        arb_one(y.raw());
        arb_mul_2exp_si(y.raw(), y.raw(), -1);
        arb_add(x.raw(), x.raw(), y.raw(), prec);
        if (!digamma_by_shift(y, x, prec)) {
            return false;
        }
        arb_add(constant.raw(), constant.raw(), y.raw(), prec);
        arb_mul_2exp_si(constant.raw(), constant.raw(), -1);
        arb_sub(complex_sum.raw(), complex_sum.raw(), constant.raw(), prec);
    }

    arb_mul_si(total.raw(), real_sum.raw(), a, prec);
    arb_mul_si(y.raw(), complex_sum.raw(), b, prec);
    arb_add(total.raw(), total.raw(), y.raw(), prec);

    // U V part: (4/D) [psi(g/D) + psi((1+g)/D) - psi((1+3g)/D)
    //                  - psi((2+3g)/D)]; every argument is positive.
    flint::Arb uv;
    arb_div(x.raw(), g.raw(), d.raw(), prec);
    arb_digamma(uv.raw(), x.raw(), prec);
    arb_add_ui(x.raw(), g.raw(), 1, prec);
    arb_div(x.raw(), x.raw(), d.raw(), prec);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_add(uv.raw(), uv.raw(), y.raw(), prec);
    arb_mul_ui(x.raw(), g.raw(), 3, prec);
    arb_add_ui(x.raw(), x.raw(), 1, prec);
    arb_div(x.raw(), x.raw(), d.raw(), prec);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_sub(uv.raw(), uv.raw(), y.raw(), prec);
    arb_mul_ui(x.raw(), g.raw(), 3, prec);
    arb_add_ui(x.raw(), x.raw(), 2, prec);
    arb_div(x.raw(), x.raw(), d.raw(), prec);
    arb_digamma(y.raw(), x.raw(), prec);
    arb_sub(uv.raw(), uv.raw(), y.raw(), prec);
    arb_mul_ui(uv.raw(), uv.raw(), 4, prec);
    arb_div(uv.raw(), uv.raw(), d.raw(), prec);
    arb_add(total.raw(), total.raw(), uv.raw(), prec);

    // (G + n log pi) / 2
    arb_const_pi(x.raw(), prec);
    arb_log(x.raw(), x.raw(), prec);
    arb_mul_si(x.raw(), x.raw(), r1 + 2 * r2, prec);
    arb_add(total.raw(), total.raw(), x.raw(), prec);
    arb_mul_2exp_si(total.raw(), total.raw(), -1);

    if (arb_is_finite(total.raw()) == 0) {
        return false;
    }
    arb_swap(out.raw(), total.raw());
    return true;
}

// Zimmert 1981, Korollar 1 (twin classes), in the finite-m form above: for
// every ideal class K, K or D K^(-1) contains an integral ideal of norm at
// most sqrt|d| exp(-c), with c = (G_{m,j}(g) + n log pi)/2 at the tabulated
// gamma.  Sets out to the ceiling of a rigorous upper bound for it.  This is
// not a bound on the smallest ideal of every class; it proves generation
// only together with the class of the different (see the callers).
bool zimmert_twin_class_bound(flint::Fmpz& out,
                              flint::FmpzConstRef abs_discriminant,
                              slong degree,
                              slong real_places,
                              slong complex_pairs) noexcept {
    if (degree < 3 || degree > kZimmertMaxDegree || real_places < 0 ||
        complex_pairs < 0 || real_places + 2 * complex_pairs != degree) {
        return false;
    }
    const ulong gamma =
        kZimmertTwinClassGammaHundredths[degree - 3][complex_pairs];
    flint::Arb log_lower;
    if (gamma == 0 ||
        !zimmert_twin_class_log_lower_bound(log_lower, real_places,
                                            complex_pairs, gamma)) {
        return false;
    }
    flint::Arb value;
    flint::Arb root;
    arb_neg(log_lower.raw(), log_lower.raw());
    arb_exp(value.raw(), log_lower.raw(), kGenerationBoundPrecision);
    arb_set_fmpz(root.raw(), abs_discriminant.raw());
    arb_sqrt(root.raw(), root.raw(), kGenerationBoundPrecision);
    arb_mul(value.raw(), value.raw(), root.raw(), kGenerationBoundPrecision);
    return ceil_upper_endpoint(out, value);
}


// Minkowski's bound n! n^(-n) (4/pi)^r2 sqrt|d| (Lang, Algebraic Number
// Theory, p. 119, Th. 4, as quoted by Zimmert 1981, p. 367).  Sets out to the
// ceiling of a rigorous upper bound for it.
bool minkowski_bound(flint::Fmpz& out,
                     flint::FmpzConstRef abs_discriminant,
                     slong degree,
                     slong complex_pairs) noexcept {
    const slong prec = kGenerationBoundPrecision;
    flint::Arb value;
    flint::Arb x;
    flint::Fmpz integer;

    arb_set_fmpz(value.raw(), abs_discriminant.raw());
    arb_sqrt(value.raw(), value.raw(), prec);
    fmpz_fac_ui(integer.raw(), static_cast<ulong>(degree));
    arb_mul_fmpz(value.raw(), value.raw(), integer.raw(), prec);
    fmpz_set_ui(integer.raw(), static_cast<ulong>(degree));
    fmpz_pow_ui(integer.raw(), integer.raw(), static_cast<ulong>(degree));
    arb_div_fmpz(value.raw(), value.raw(), integer.raw(), prec);
    arb_const_pi(x.raw(), prec);
    arb_ui_div(x.raw(), 4, x.raw(), prec);
    arb_pow_ui(x.raw(), x.raw(), static_cast<ulong>(complex_pairs), prec);
    arb_mul(value.raw(), value.raw(), x.raw(), prec);
    return ceil_upper_endpoint(out, value);
}

// Proven bound for degree >= 3: the smallest of Zimmert's bound (degree at
// most 20), Minkowski's bound with (4/pi)^r2, and the exact integer
// Minkowski form with 2^r2.  Each is a theorem, so the minimum is too.
void generation_bound(flint::Fmpz& out,
                      flint::FmpzConstRef abs_discriminant,
                      slong degree,
                      slong real_places,
                      slong complex_pairs) noexcept {
    generic_minkowski_bound(out, abs_discriminant, degree, complex_pairs);
    flint::Fmpz candidate;
    if (minkowski_bound(candidate, abs_discriminant, degree, complex_pairs) &&
        fmpz_cmp(candidate.raw(), out.raw()) < 0) {
        fmpz_swap(out.raw(), candidate.raw());
    }
    if (degree <= kZimmertMaxDegree &&
        zimmert_bound(candidate, abs_discriminant, degree, real_places,
                      complex_pairs) &&
        fmpz_cmp(candidate.raw(), out.raw()) < 0) {
        fmpz_swap(out.raw(), candidate.raw());
    }
    if (fmpz_cmp_ui(out.raw(), 1) < 0) {
        fmpz_one(out.raw());
    }
}

bool prime_ideal_norm_at_most_bound(const PrimeIdeal& prime,
                                    flint::FmpzConstRef rational_prime,
                                    flint::FmpzConstRef bound) noexcept {
    const slong residue_degree = prime.residue_degree();
    if (residue_degree <= 0) {
        return false;
    }

    flint::Fmpz norm;
    fmpz_pow_ui(norm.raw(), rational_prime.raw(),
                static_cast<ulong>(residue_degree));
    return fmpz_cmp(norm.raw(), bound.raw()) <= 0;
}

bool is_inert_prime_decomposition(const PrimeIdealList& decomposed,
                                  slong degree) noexcept {
    if (degree <= 0 || decomposed.size() != 1) {
        return false;
    }

    const PrimeIdeal* prime = decomposed.at(0);
    return prime != nullptr &&
           prime->ramification_index() == 1 &&
           prime->residue_degree() == degree;
}

bool generation_bound_inputs_valid(flint::FmpzConstRef abs_discriminant,
                                   slong real_places,
                                   slong complex_pairs) noexcept {
    return real_places >= 0 && complex_pairs >= 0 &&
           real_places + 2 * complex_pairs >= 3 &&
           fmpz_sgn(abs_discriminant.raw()) > 0;
}

}  // namespace

bool detail::generation_bound(flint::FmpzRef out,
                              flint::FmpzConstRef abs_discriminant,
                              slong real_places,
                              slong complex_pairs) noexcept {
    if (!generation_bound_inputs_valid(abs_discriminant, real_places,
                                       complex_pairs)) {
        return false;
    }
    flint::Fmpz bound;
    silex::generation_bound(bound, abs_discriminant,
                            real_places + 2 * complex_pairs, real_places,
                            complex_pairs);
    fmpz_set(out.raw(), bound.raw());
    return true;
}

bool detail::zimmert_generation_bound(flint::FmpzRef out,
                                      flint::FmpzConstRef abs_discriminant,
                                      slong real_places,
                                      slong complex_pairs) noexcept {
    if (!generation_bound_inputs_valid(abs_discriminant, real_places,
                                       complex_pairs)) {
        return false;
    }
    flint::Fmpz bound;
    if (!zimmert_bound(bound, abs_discriminant,
                       real_places + 2 * complex_pairs, real_places,
                       complex_pairs)) {
        return false;
    }
    fmpz_set(out.raw(), bound.raw());
    return true;
}

bool detail::zimmert_twin_class_log_lower_bound(
        flint::ArbRef out,
        slong real_places,
        slong complex_pairs,
        ulong gamma_hundredths) noexcept {
    flint::Arb value;
    if (!silex::zimmert_twin_class_log_lower_bound(
                value, real_places, complex_pairs, gamma_hundredths)) {
        return false;
    }
    arb_swap(out.raw(), value.raw());
    return true;
}

bool detail::twin_class_generation_bound(
        flint::FmpzRef out,
        flint::FmpzConstRef abs_discriminant,
        slong real_places,
        slong complex_pairs) noexcept {
    if (!generation_bound_inputs_valid(abs_discriminant, real_places,
                                       complex_pairs)) {
        return false;
    }
    const slong degree = real_places + 2 * complex_pairs;
    flint::Fmpz bound;
    if (!zimmert_twin_class_bound(bound, abs_discriminant, degree,
                                  real_places, complex_pairs)) {
        return false;
    }
    flint::Fmpz classical;
    silex::generation_bound(classical, abs_discriminant, degree, real_places,
                            complex_pairs);
    if (fmpz_cmp(classical.raw(), bound.raw()) < 0) {
        fmpz_swap(bound.raw(), classical.raw());
    }
    if (fmpz_cmp_ui(bound.raw(), 1) < 0) {
        fmpz_one(bound.raw());
    }
    fmpz_set(out.raw(), bound.raw());
    return true;
}

FactorBase::FactorBase(const Order& parent) noexcept {
    define(parent);
}

FactorBase::~FactorBase() noexcept = default;

FactorBase::FactorBase(FactorBase&& other) noexcept {
    swap(other);
}

FactorBase& FactorBase::operator=(FactorBase&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FactorBase::swap(FactorBase& other) noexcept {
    parent_.swap(other.parent_);
    primes_.swap(other.primes_);
    blocks_.swap(other.blocks_);
    std::swap(complete_rational_prime_blocks_,
              other.complete_rational_prime_blocks_);
}

void FactorBase::clear() noexcept {
    primes_.clear();
    blocks_.clear();
    parent_.clear();
    complete_rational_prime_blocks_ = false;
}

bool FactorBase::define(const Order& parent) noexcept {
    if (!parent.has_basis()) {
        return false;
    }

    clear();
    parent_ = parent;
    return true;
}

bool FactorBase::set(const FactorBase& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    FactorBase copy(other.parent_);
    if (!copy.is_defined()) {
        return false;
    }

    for (const PrimeIdeal& prime : other.primes_) {
        if (!copy.append_prime(prime)) {
            return false;
        }
    }
    if (copy.blocks_.size() != other.blocks_.size()) {
        return false;
    }
    for (std::size_t i = 0; i < copy.blocks_.size(); ++i) {
        copy.blocks_[i].complete = other.blocks_[i].complete;
    }
    copy.complete_rational_prime_blocks_ =
            other.complete_rational_prime_blocks_;

    swap(copy);
    return true;
}

bool FactorBase::is_defined() const noexcept {
    return parent_.has_basis();
}

const Order* FactorBase::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong FactorBase::length() const noexcept {
    return is_defined() ? static_cast<slong>(primes_.size()) : 0;
}

slong FactorBase::rational_prime_block_count() const noexcept {
    return is_defined() ? static_cast<slong>(blocks_.size()) : 0;
}

bool FactorBase::rational_prime_blocks_are_complete() const noexcept {
    return is_defined() && complete_rational_prime_blocks_;
}

bool detail::FactorBaseBlockAccess::rational_prime_block_is_complete(
        bool& complete,
        const FactorBase& base,
        slong block_index) noexcept {
    complete = false;
    if (!base.is_defined() || block_index < 0 ||
        block_index >= base.rational_prime_block_count()) {
        return false;
    }
    complete = base.complete_rational_prime_blocks_ ||
               base.blocks_[static_cast<std::size_t>(block_index)].complete;
    return true;
}

bool FactorBase::build(flint::FmpzConstRef bound) noexcept {
    return build(bound, nullptr);
}

bool FactorBase::build(flint::FmpzConstRef bound,
                       const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "factor_base.build_complete");
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0 ||
        !parent_.is_maximal()) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList decomposed;
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "factor_base.decompose_prime");
            if (!decompose_prime(decomposed, parent_,
                                 flint::FmpzConstRef(p), 0,
                                 diagnostics)) {
                return false;
            }
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "factor_base.append_prime_block");
            for (slong i = 0; i < decomposed.size(); ++i) {
                const PrimeIdeal* prime = decomposed.at(i);
                if (prime == nullptr || !candidate.append_prime(*prime)) {
                    return false;
                }
            }
        }

        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    candidate.complete_rational_prime_blocks_ = true;
    swap(candidate);
    return true;
}

bool FactorBase::build_prime_ideal_norm_bounded(
        flint::FmpzConstRef bound) noexcept {
    return build_prime_ideal_norm_bounded_impl(bound, true, 0);
}

bool FactorBase::build_relation_completion_base(
        flint::FmpzConstRef bound) noexcept {
    return build_prime_ideal_norm_bounded_impl(bound, false, 0);
}

bool FactorBase::build_lll_relation_base(
        flint::FmpzConstRef bound) noexcept {
    return build_lll_relation_base_impl(bound, false);
}

bool FactorBase::build_maximal_imaginary_quadratic_relation_base(
        flint::FmpzConstRef bound,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "factor_base.maximal_imaginary_quadratic.build");
    flint::Fmpz discriminant;
    if (parent_.degree() != 2 || !parent_.is_maximal() ||
        !parent_.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0) {
        return false;
    }
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdeal retained(parent_);
        detail::RetainedQuadraticPrimeKind kind =
                detail::RetainedQuadraticPrimeKind::inert;
        bool retained_ready = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "factor_base.maximal_imaginary_quadratic.classify_prime");
            retained_ready = retained.is_defined() &&
                    detail::MaximalQuadraticPrimeAccess::
                            set_first_degree_one_prime(
                                    retained, kind, parent_,
                                    flint::FmpzConstRef(p), diagnostics);
        }
        if (!retained_ready) {
            // The direct path has not published candidate. Preserve the
            // established generic decomposition route as the atomic fallback.
            return build_prime_ideal_norm_bounded_impl(bound, false, 0, true);
        }
        if (kind == detail::RetainedQuadraticPrimeKind::inert) {
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "factor_base.maximal_imaginary_quadratic.append_prime");
            if (!candidate.append_prime(std::move(retained))) {
                return false;
            }

            bool found = false;
            const slong block = candidate.block_position(
                    flint::FmpzConstRef(p), found);
            if (!found) {
                return false;
            }
            candidate.blocks_[static_cast<std::size_t>(block)].complete =
                    kind == detail::RetainedQuadraticPrimeKind::ramified;
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    swap(candidate);
    return true;
}

bool FactorBase::build_lll_relation_base_impl(
        flint::FmpzConstRef bound,
        bool incomplete) noexcept {
    if (incomplete
                ? !build_prime_ideal_norm_bounded_impl(
                          bound, true, 0)
                : !build(bound)) {
        return false;
    }

    struct NormIndex {
        slong index = 0;
        flint::Fmpz norm;
    };

    std::vector<NormIndex> order;
    order.reserve(primes_.size());
    for (slong i = 0; i < static_cast<slong>(primes_.size()); ++i) {
        NormIndex entry;
        entry.index = i;
        if (!primes_[static_cast<std::size_t>(i)].norm(
                    flint::FmpzRef(entry.norm))) {
            return false;
        }
        order.emplace_back(std::move(entry));
    }

    std::sort(order.begin(), order.end(),
              [](const NormIndex& left, const NormIndex& right) noexcept {
                  const int cmp = flint::fmpz_cmp(
                          flint::FmpzConstRef(left.norm),
                          flint::FmpzConstRef(right.norm));
                  if (cmp != 0) {
                      return cmp > 0;
                  }
                  return left.index < right.index;
              });

    std::vector<PrimeIdeal> sorted;
    sorted.reserve(primes_.size());
    for (const NormIndex& entry : order) {
        PrimeIdeal copy(parent_);
        if (!copy.is_defined() ||
            !copy.set(primes_[static_cast<std::size_t>(entry.index)])) {
            return false;
        }
        sorted.emplace_back(std::move(copy));
    }

    primes_.swap(sorted);
    return rebuild_blocks();
}

bool FactorBase::build_prime_ideal_norm_bounded_impl(
        flint::FmpzConstRef bound,
        bool include_inert_primes,
        slong degree_limit,
        bool retain_one_split_prime) noexcept {
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0 ||
        !parent_.is_maximal() || degree_limit < 0) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList decomposed;
        if (!decompose_prime(decomposed, parent_,
                             flint::FmpzConstRef(p))) {
            return false;
        }
        // The source factor-base generator skips inert rational primes for the relation
        // factor base, even when p^degree is still below the tested bound.
        if (!include_inert_primes &&
            is_inert_prime_decomposition(decomposed, parent_.degree())) {
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        slong retained = 0;
        for (slong i = 0; i < decomposed.size(); ++i) {
            const PrimeIdeal* prime = decomposed.at(i);
            if (prime == nullptr) {
                return false;
            }
            if (degree_limit > 0 &&
                prime->residue_degree() > degree_limit) {
                continue;
            }
            if (!prime_ideal_norm_at_most_bound(
                        *prime, flint::FmpzConstRef(p), bound)) {
                continue;
            }
            if (!candidate.append_prime(*prime)) {
                return false;
            }
            ++retained;
            // reference's quadratic FBquad stores one prime form for a split
            // rational prime; its inverse represents the conjugate class.
            if (retain_one_split_prime && decomposed.size() == 2) {
                break;
            }
        }
        if (retained > 0) {
            // The source factor-base generator marks LV[p] when every prime in the cached
            // decomposition survives the relation-base filters.  Retain that
            // fact here so init_rel and subFBgen need not decompose p again.
            bool found = false;
            const slong block = candidate.block_position(
                    flint::FmpzConstRef(p), found);
            if (!found) {
                return false;
            }
            candidate.blocks_[static_cast<std::size_t>(block)].complete =
                    retained == decomposed.size();
        }

        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    swap(candidate);
    return true;
}

bool FactorBase::prime(PrimeIdeal& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        !same_order_parent(out.parent(), &parent_)) {
        return false;
    }
    return out.set(primes_[static_cast<std::size_t>(index)]);
}

const PrimeIdeal* FactorBase::prime_at(slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length()) {
        return nullptr;
    }
    return &primes_[static_cast<std::size_t>(index)];
}

bool FactorBase::rational_prime_block_data(
        flint::FmpzRef rational_prime,
        slong& length,
        slong block_index) const noexcept {
    length = 0;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    fmpz_set(rational_prime.raw(), block.rational_prime.raw());
    length = block.length;
    return true;
}

bool FactorBase::rational_prime_block_index(slong& index,
                                            slong block_index,
                                            slong offset) const noexcept {
    index = -1;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    if (offset < 0 || offset >= static_cast<slong>(block.indices.size())) {
        return false;
    }

    index = block.indices[static_cast<std::size_t>(offset)];
    return true;
}

bool FactorBase::rational_prime_block(flint::FmpzRef rational_prime,
                                      slong& start,
                                      slong& length,
                                      slong block_index) const noexcept {
    start = 0;
    length = 0;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    if (!block.contiguous) {
        return false;
    }
    fmpz_set(rational_prime.raw(), block.rational_prime.raw());
    start = block.start;
    length = block.length;
    return true;
}

bool FactorBase::rational_prime_block_index_for_prime(
        slong& block_index,
        flint::FmpzConstRef rational_prime) const noexcept {
    block_index = -1;
    if (!is_defined()) {
        return false;
    }

    bool found = false;
    const slong position = block_position(rational_prime, found);
    if (!found) {
        return false;
    }
    block_index = position;
    return true;
}

bool FactorBase::equal(const FactorBase& other) const noexcept {
    if (!is_defined() || !other.is_defined() ||
        !same_order_parent(parent(), other.parent()) ||
        length() != other.length() ||
        rational_prime_block_count() != other.rational_prime_block_count()) {
        return false;
    }

    for (slong i = 0; i < rational_prime_block_count(); ++i) {
        const PrimeBlock& left = blocks_[static_cast<std::size_t>(i)];
        const PrimeBlock& right = other.blocks_[static_cast<std::size_t>(i)];
        if (!flint::fmpz_equal(left.rational_prime, right.rational_prime) ||
            left.start != right.start || left.length != right.length ||
            left.contiguous != right.contiguous ||
            left.indices != right.indices) {
            return false;
        }
    }

    for (slong i = 0; i < length(); ++i) {
        if (!primes_[static_cast<std::size_t>(i)].equal(
                    other.primes_[static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    return true;
}

slong FactorBase::index(const PrimeIdeal& prime) const noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return -1;
    }

    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return -1;
    }

    bool found = false;
    const slong block =
            block_position(flint::FmpzConstRef(rational_prime), found);
    if (!found) {
        return -1;
    }

    const PrimeBlock& prime_block = blocks_[static_cast<std::size_t>(block)];
    for (slong index : prime_block.indices) {
        if (primes_[static_cast<std::size_t>(index)].equal(prime)) {
            return index;
        }
    }
    return -1;
}

bool FactorBase::append_prime(const PrimeIdeal& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return false;
    }

    PrimeIdeal copy(parent_);
    flint::Fmpz rational_prime;
    if (!copy.is_defined() || !copy.set(prime) ||
        !copy.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    const slong index = static_cast<slong>(primes_.size());
    primes_.emplace_back(std::move(copy));
    return add_block_index(flint::FmpzConstRef(rational_prime), index);
}

bool FactorBase::append_prime(PrimeIdeal&& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return false;
    }

    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    const slong index = static_cast<slong>(primes_.size());
    primes_.emplace_back(std::move(prime));
    return add_block_index(flint::FmpzConstRef(rational_prime), index);
}

slong FactorBase::block_position(
        flint::FmpzConstRef rational_prime,
        bool& found) const noexcept {
    found = false;
    slong lo = 0;
    slong hi = static_cast<slong>(blocks_.size());
    while (lo < hi) {
        const slong mid = lo + (hi - lo) / 2;
        if (fmpz_cmp(blocks_[static_cast<std::size_t>(mid)]
                             .rational_prime.raw(),
                     rational_prime.raw()) < 0) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    if (lo < static_cast<slong>(blocks_.size()) &&
        fmpz_equal(blocks_[static_cast<std::size_t>(lo)]
                           .rational_prime.raw(),
                   rational_prime.raw()) != 0) {
        found = true;
    }
    return lo;
}

void FactorBase::refresh_block(PrimeBlock& block) noexcept {
    block.length = static_cast<slong>(block.indices.size());
    block.start = block.indices.empty() ? 0 : block.indices.front();
    block.contiguous = true;
    for (slong i = 0; i < static_cast<slong>(block.indices.size()); ++i) {
        if (block.indices[static_cast<std::size_t>(i)] != block.start + i) {
            block.contiguous = false;
            break;
        }
    }
}

bool FactorBase::add_block_index(flint::FmpzConstRef rational_prime,
                                 slong prime_index) noexcept {
    if (prime_index < 0 || prime_index >= static_cast<slong>(primes_.size())) {
        return false;
    }

    bool found = false;
    const slong pos = block_position(rational_prime, found);
    if (found) {
        PrimeBlock& block = blocks_[static_cast<std::size_t>(pos)];
        block.indices.push_back(prime_index);
        refresh_block(block);
        return true;
    }

    PrimeBlock block;
    fmpz_set(block.rational_prime.raw(), rational_prime.raw());
    block.indices.push_back(prime_index);
    refresh_block(block);
    blocks_.insert(blocks_.begin() + pos, std::move(block));
    return true;
}

bool FactorBase::rebuild_blocks() noexcept {
    std::vector<PrimeBlock> previous_blocks = std::move(blocks_);
    blocks_.clear();
    flint::Fmpz rational_prime;
    for (slong i = 0; i < static_cast<slong>(primes_.size()); ++i) {
        if (!primes_[static_cast<std::size_t>(i)].rational_prime(
                    flint::FmpzRef(rational_prime)) ||
            !add_block_index(flint::FmpzConstRef(rational_prime), i)) {
            return false;
        }
    }
    if (blocks_.size() != previous_blocks.size()) {
        return false;
    }
    for (std::size_t i = 0; i < blocks_.size(); ++i) {
        if (!flint::fmpz_equal(blocks_[i].rational_prime,
                               previous_blocks[i].rational_prime)) {
            return false;
        }
        blocks_[i].complete = previous_blocks[i].complete;
    }
    return true;
}

bool FactorBase::contains(const PrimeIdeal& prime) const noexcept {
    return index(prime) >= 0;
}

bool factor_base_class_group_bound(flint::FmpzRef out,
                                   const Order& order) noexcept {
    if (!order.is_maximal()) {
        return false;
    }

    const slong degree = order.degree();
    flint::Fmpz bound;
    if (degree == 1) {
        fmpz_one(bound.raw());
        fmpz_set(out.raw(), bound.raw());
        return true;
    }

    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        fmpz_is_zero(discriminant.raw()) != 0 ||
        order.parent() == nullptr) {
        return false;
    }

    Signature sig;
    if (!sig.compute(*order.parent())) {
        return false;
    }

    fmpz_abs(abs_discriminant.raw(), discriminant.raw());
    if (degree == 2 && fmpz_sgn(discriminant.raw()) < 0) {
        floor_sqrt_div_ui(bound, flint::FmpzConstRef(abs_discriminant), 3);
    } else if (degree == 2) {
        fmpz_sqrt(bound.raw(), abs_discriminant.raw());
        fmpz_fdiv_q_2exp(bound.raw(), bound.raw(), 1);
    } else {
        generation_bound(bound,
                         flint::FmpzConstRef(abs_discriminant),
                         degree,
                         sig.r1(),
                         sig.r2());
    }

    fmpz_set(out.raw(), bound.raw());
    return true;
}

namespace {

// True when the stored defining polynomial of `field` is monic with integer
// coefficients, so that its root theta is an algebraic integer and
// Z[theta] is the equation order.
bool field_polynomial_is_monic_integral(const NumberField& field) noexcept {
    const nf_struct* raw = field.raw_flint_field();
    if (raw == nullptr) {
        return false;
    }
    const slong length = fmpq_poly_length(raw->pol);
    return length >= 2 && fmpz_is_one(fmpq_poly_denref(raw->pol)) != 0 &&
           fmpz_is_one(fmpq_poly_numref(raw->pol) + length - 1) != 0;
}

// The largest norm of a prime ideal of the maximal order `order` above a
// rational prime dividing `index`.
bool largest_prime_norm_above(flint::Fmpz& out,
                              const Order& order,
                              flint::FmpzConstRef index) noexcept {
    flint::FmpzFactor factorization;
    flint::fmpz_factor(flint::FmpzFactorRef(factorization), index);
    const slong factor_count =
            flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization));
    flint::Fmpz largest;
    flint::Fmpz p;
    flint::Fmpz norm;
    for (slong i = 0; i < factor_count; ++i) {
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(p), flint::FmpzFactorConstRef(factorization),
                i);
        PrimeIdealList primes;
        if (!decompose_prime(primes, order, flint::FmpzConstRef(p)) ||
            primes.size() < 1) {
            return false;
        }
        for (slong k = 0; k < primes.size(); ++k) {
            const PrimeIdeal* prime = primes.at(k);
            if (prime == nullptr || !prime->norm(flint::FmpzRef(norm))) {
                return false;
            }
            if (fmpz_cmp(norm.raw(), largest.raw()) > 0) {
                fmpz_set(largest.raw(), norm.raw());
            }
        }
    }
    fmpz_swap(out.raw(), largest.raw());
    return true;
}

}  // namespace

// Zimmert 1981, Korollar 1 gives, for every class K, an integral ideal of
// norm at most B in K or in D K^(-1).  With every prime ideal of norm at
// most the requirement in the generated subgroup H, K lies in H as soon as
// [D] does.  For theta with monic minimal polynomial f in Z[x], Cohen,
// Advanced Topics in Computational Number Theory (GTM 193), Lemma 10.1.20
// and the proof of Corollary 10.1.21 (pp. 484-485) give
// f'(theta) O_K = F D for the conductor F = {x : x O_K in Z[theta]}, so
// [D] = [F]^(-1).  F = O_K when the index is 1 (Cohen, A Course in
// Computational Algebraic Number Theory, GTM 138, Exercise 4.30), and
// otherwise every prime dividing F lies above some p | [O_K : Z[theta]],
// because the index lies in F.  Anything else falls back to the bound of
// factor_base_class_group_bound (Satz 2 and Minkowski).
bool detail::proven_generation_requirement(
        flint::FmpzRef bound,
        GenerationDifferentRoute& route,
        const Order& order) noexcept {
    flint::Fmpz classical;
    if (!factor_base_class_group_bound(flint::FmpzRef(classical), order)) {
        return false;
    }
    flint::Fmpz chosen;
    fmpz_set(chosen.raw(), classical.raw());
    GenerationDifferentRoute chosen_route = GenerationDifferentRoute::none;

    const slong degree = order.degree();
    const NumberField* field = order.parent();
    flint::Fmpz discriminant;
    Signature sig;
    flint::Fmpz twin;
    if (degree >= 3 && degree <= kZimmertMaxDegree && field != nullptr &&
        field_polynomial_is_monic_integral(*field) &&
        order.discriminant(flint::FmpzRef(discriminant)) &&
        fmpz_is_zero(discriminant.raw()) == 0 && sig.compute(*field)) {
        fmpz_abs(discriminant.raw(), discriminant.raw());
        if (zimmert_twin_class_bound(twin,
                                     flint::FmpzConstRef(discriminant),
                                     degree, sig.r1(), sig.r2()) &&
            fmpz_cmp(twin.raw(), classical.raw()) < 0) {
            if (fmpz_cmp_ui(twin.raw(), 1) < 0) {
                fmpz_one(twin.raw());
            }
            const Order equation = Order::equation_order(*field);
            flint::Fmpz index;
            flint::Fmpz conductor_norm;
            if (equation.is_defined() &&
                order_index(flint::FmpzRef(index), equation, order)) {
                if (fmpz_is_one(index.raw()) != 0) {
                    fmpz_set(chosen.raw(), twin.raw());
                    chosen_route = GenerationDifferentRoute::index_one;
                } else if (largest_prime_norm_above(
                                   conductor_norm, order,
                                   flint::FmpzConstRef(index))) {
                    if (fmpz_cmp(conductor_norm.raw(), twin.raw()) < 0) {
                        fmpz_set(conductor_norm.raw(), twin.raw());
                    }
                    if (fmpz_cmp(conductor_norm.raw(), classical.raw()) < 0) {
                        fmpz_swap(chosen.raw(), conductor_norm.raw());
                        chosen_route =
                                GenerationDifferentRoute::conductor_primes;
                    }
                }
            }
        }
    }

    fmpz_set(bound.raw(), chosen.raw());
    route = chosen_route;
    return true;
}

}  // namespace silex
