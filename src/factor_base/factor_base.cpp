#include <silex/factor_base.hpp>

#include "factor_base_internal.hpp"
#include "../prime_ideal/prime_ideal_internal.hpp"

#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/signature.hpp>

#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <utility>

namespace silex {
namespace {

// Source of FactorBase content identities.  One counter for every base keeps
// an identity from being reused by a later, different base; zero is never
// issued.
std::atomic<std::uint64_t> next_factor_base_identity{1};

std::uint64_t fresh_factor_base_identity() noexcept {
    return next_factor_base_identity.fetch_add(1, std::memory_order_relaxed);
}

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
// every gamma > 0, and the certified enclosure of F evaluated below
// (GB-zimmert-F-tail-bound) holds for every gamma > 0, so gamma affects only
// the size of the bound, never its validity.  Each entry maximises F over
// 0 < gamma < 3 in steps of 1/100 (integers excluded in the oracle,
// test/zimmert_twin_class_oracle.py --gammas).
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

// Certified enclosure of Zimmert's F_{a,b}(gamma) (Zimmert 1981, p. 376)
// by the theorem GB-zimmert-F-tail-bound of the silex-research repository,
// research/proofs/generation/GB-zimmert-F-tail-bound.md at commit 06191ff,
// following the evaluation procedure stated there ("Evaluation with a
// certified enclosure").  Notation as in the proof: c = 1 + 2 gamma,
// h = 2c, rho = min(1/2, 3/(2c)), sigma = rho/2,
//
//   T(t) = (2/c)(psi((2t + 3 gamma)/h) - psi((2t - 1 + gamma)/h))
//          - 1/(2t - 2 - gamma) - 1/(2t - 1 + gamma),
//
// and for N = 2(L + 1 + delta), delta in {0, 1/2}, the main term
// mu(N, K) = sum_{k=1}^{K-1} Phi_k k!/N^(k+1) and the error bound
//
//   E(N, K) = 2 M_Phi K!/(rho^K N^(K+1))
//             + 3 e^(-(N-2-gamma) sigma)/((N - 2 - gamma)(1 - e^(-2 sigma)))
//             + M_Phi e^(-N sigma)/(N - (K - 1)/sigma),
//
// valid when (L, delta, K) is admissible: N > 2 + gamma and
// N > (K - 1)/sigma.  Part (i): |sum_{l>L} T(l + delta) - mu(N, K)| <=
// E(N, K).  Part (iv): for integers L >= 1 with L > gamma/2, K >= 2 and
// (L, 0, K) admissible, F_{a,b}(gamma) lies within a E(2L+2, K) +
// b E(2L+3, K) of S_L(a, b; gamma) + a mu(2L+2, K) + b mu(2L+3, K), where
//
//   S_L = a (H_a(L) - psi(1/2 + gamma/2)/2 - psi(L - gamma/2)/2)
//       + b (H_b(L) - psi(1 + gamma/2)/2 - psi(L + 1/2 - gamma/2)/2) + C,
//   H_a(L) = sum_{l=1}^{L} [(2/c)(psi((2l+3gamma)/h) - psi((2l-1+gamma)/h))
//                           - 1/(2l-1+gamma)],
//   H_b(L) = sum_{l=1}^{L} [(2/c)(psi((2l+1+3gamma)/h) - psi((2l+gamma)/h))
//                           - 1/(2l+gamma)],
//   C = -4/gamma + (2/c)(psi((1+gamma)/h) - psi((1+3gamma)/h)
//                        + psi((2+5gamma)/h) - psi((2+3gamma)/h)).
//
// Every psi argument and denominator is an exact positive rational, so
// arb_digamma is only evaluated at positive arguments; integer gamma needs
// no special case (part (iii)).

// Target 2^-kZimmertTailBits for each tail error E (proof, "Choosing L and
// K for a target error 2^-p per tail").
constexpr slong kZimmertTailBits = kGenerationBoundPrecision;
// Largest Taylor order K tried before L is increased.
constexpr slong kZimmertMaxTailOrder = 64;
// Largest L tried.
constexpr slong kZimmertMaxTailTerms = slong{1} << 16;

// Sets out to num/den (den > 0).
void set_rational(flint::Fmpq& out, slong num, slong den) noexcept {
    fmpq_set_si(out.raw(), num, static_cast<ulong>(den));
}

// Multiplies x by num/den (den > 0).
void scale_rational(flint::Fmpq& x, slong num, slong den) noexcept {
    flint::Fmpq ratio;
    set_rational(ratio, num, den);
    fmpq_mul(x.raw(), x.raw(), ratio.raw());
}

// The gamma-dependent constants of GB-zimmert-F-tail-bound: gamma, c, rho,
// sigma (exact rationals), M_Phi (a ball holding the upper endpoint of its
// enclosure, exactly) and the Taylor coefficients Phi_k of
// Phi(s) = phi(s)/(1 - e^(-2s)), computed exactly by the recurrences of the
// Statement for k < order.  The q_k, phi_k and the factorial series
// coefficients are kept to extend the recurrences.
struct ZimmertTailConstants {
    // gamma = numerator/denominator in lowest terms.
    slong numerator = 0;
    slong denominator = 1;
    flint::Fmpq gamma;
    flint::Fmpq c;
    flint::Fmpq rho;
    flint::Fmpq sigma;
    flint::Arb m_phi;
    // Largest k with q_k, phi_k and the series coefficients set (-1 for
    // none); Phi_0, ..., Phi_{computed-1} are set.
    slong computed = -1;
    flint::Fmpq q[kZimmertMaxTailOrder + 1];
    flint::Fmpq phi_small[kZimmertMaxTailOrder + 1];  // phi_k
    flint::Fmpq exp_c[kZimmertMaxTailOrder + 1];      // (-c)^i/i!
    flint::Fmpq exp_a[kZimmertMaxTailOrder + 1];      // (1-gamma)^i/i!
    flint::Fmpq exp_b[kZimmertMaxTailOrder + 1];      // (2+gamma)^i/i!
    flint::Fmpq d[kZimmertMaxTailOrder + 1];          // (-1)^k 2^(k+1)/(k+1)!
    flint::Fmpq phi_big[kZimmertMaxTailOrder];        // Phi_k
};

// Computes q_k, phi_k and the series coefficients for k <= index, and
// Phi_k for k < index, extending what is already set.
void zimmert_tail_extend(ZimmertTailConstants& t, slong index) noexcept {
    flint::Fmpq term;
    flint::Fmpq factor;
    for (slong k = t.computed + 1; k <= index; ++k) {
        if (k == 0) {
            fmpq_one(t.exp_c[0].raw());
            fmpq_one(t.exp_a[0].raw());
            fmpq_one(t.exp_b[0].raw());
            set_rational(t.q[0], 1, 2);
            fmpq_set_si(t.d[0].raw(), 2, 1);
        } else {
            // (-c)^k/k!, (1-gamma)^k/k!, (2+gamma)^k/k!.
            fmpq_neg(factor.raw(), t.c.raw());
            fmpq_mul(t.exp_c[k].raw(), t.exp_c[k - 1].raw(), factor.raw());
            scale_rational(t.exp_c[k], 1, k);
            fmpq_one(factor.raw());
            fmpq_sub(factor.raw(), factor.raw(), t.gamma.raw());
            fmpq_mul(t.exp_a[k].raw(), t.exp_a[k - 1].raw(), factor.raw());
            scale_rational(t.exp_a[k], 1, k);
            fmpq_add_si(factor.raw(), t.gamma.raw(), 2);
            fmpq_mul(t.exp_b[k].raw(), t.exp_b[k - 1].raw(), factor.raw());
            scale_rational(t.exp_b[k], 1, k);
            // d_k = d_{k-1} (-2)/(k+1).
            fmpq_set(t.d[k].raw(), t.d[k - 1].raw());
            scale_rational(t.d[k], -2, k + 1);
            // q_k = -(1/2) sum_{i=1}^{k} ((-c)^i/i!) q_{k-i}.
            fmpq_zero(t.q[k].raw());
            for (slong i = 1; i <= k; ++i) {
                fmpq_mul(term.raw(), t.exp_c[i].raw(), t.q[k - i].raw());
                fmpq_add(t.q[k].raw(), t.q[k].raw(), term.raw());
            }
            scale_rational(t.q[k], -1, 2);
        }
        // phi_k = sum_{i=0}^{k} ((1-gamma)^i/i!) g_{k-i} - (2+gamma)^k/k!,
        // with g_0 = 1 and g_j = 4 q_j for j >= 1.
        fmpq_neg(t.phi_small[k].raw(), t.exp_b[k].raw());
        for (slong i = 0; i <= k; ++i) {
            if (i == k) {
                fmpq_set(factor.raw(), t.exp_a[i].raw());  // g_0 = 1
            } else {
                fmpq_mul_si(factor.raw(), t.q[k - i].raw(), 4);
                fmpq_mul(factor.raw(), factor.raw(), t.exp_a[i].raw());
            }
            fmpq_add(t.phi_small[k].raw(), t.phi_small[k].raw(),
                     factor.raw());
        }
        // Phi_{k-1} = (1/2)(phi_k - sum_{i=1}^{k-1} d_i Phi_{k-1-i}).
        if (k >= 1) {
            const slong j = k - 1;
            fmpq_set(t.phi_big[j].raw(), t.phi_small[k].raw());
            for (slong i = 1; i <= j; ++i) {
                fmpq_mul(term.raw(), t.d[i].raw(), t.phi_big[j - i].raw());
                fmpq_sub(t.phi_big[j].raw(), t.phi_big[j].raw(), term.raw());
            }
            fmpq_div_2exp(t.phi_big[j].raw(), t.phi_big[j].raw(), 1);
        }
        t.computed = k;
    }
}

// Sets up the constants for gamma = numerator/denominator > 0 and checks
// Phi_0 = 0 and Phi_1 = -c^2/4 (proof, Statement).  Fails if M_Phi is not
// a finite ball with a positive denominator.
bool zimmert_tail_setup(ZimmertTailConstants& t,
                        slong numerator,
                        slong denominator,
                        slong prec) noexcept {
    if (numerator <= 0 || denominator <= 0) {
        return false;
    }
    set_rational(t.gamma, numerator, denominator);
    t.numerator = fmpz_get_si(fmpq_numref(t.gamma.raw()));
    t.denominator = fmpz_get_si(fmpq_denref(t.gamma.raw()));
    fmpq_mul_2exp(t.c.raw(), t.gamma.raw(), 1);
    fmpq_add_si(t.c.raw(), t.c.raw(), 1);
    // rho = min(1/2, 3/(2c)), sigma = rho/2.
    flint::Fmpq other;
    fmpq_inv(other.raw(), t.c.raw());
    fmpq_mul_si(other.raw(), other.raw(), 3);
    fmpq_div_2exp(other.raw(), other.raw(), 1);
    set_rational(t.rho, 1, 2);
    if (fmpq_cmp(other.raw(), t.rho.raw()) < 0) {
        fmpq_set(t.rho.raw(), other.raw());
    }
    fmpq_div_2exp(t.sigma.raw(), t.rho.raw(), 1);

    // M_Phi = [e^(|1-gamma| rho)(1 + 2 tan(c rho/2)) + e^((2+gamma) rho)]
    //         / (4 rho + 1 - e^(2 rho)).
    flint::Fmpq exact;
    flint::Arb x;
    flint::Arb y;
    flint::Arb value;
    fmpq_one(exact.raw());
    fmpq_sub(exact.raw(), exact.raw(), t.gamma.raw());
    fmpq_abs(exact.raw(), exact.raw());
    fmpq_mul(exact.raw(), exact.raw(), t.rho.raw());
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_exp(value.raw(), x.raw(), prec);
    fmpq_mul(exact.raw(), t.c.raw(), t.rho.raw());
    fmpq_div_2exp(exact.raw(), exact.raw(), 1);
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_tan(y.raw(), x.raw(), prec);
    arb_mul_2exp_si(y.raw(), y.raw(), 1);
    arb_add_ui(y.raw(), y.raw(), 1, prec);
    arb_mul(value.raw(), value.raw(), y.raw(), prec);
    fmpq_add_si(exact.raw(), t.gamma.raw(), 2);
    fmpq_mul(exact.raw(), exact.raw(), t.rho.raw());
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_exp(y.raw(), x.raw(), prec);
    arb_add(value.raw(), value.raw(), y.raw(), prec);
    fmpq_mul_2exp(exact.raw(), t.rho.raw(), 1);
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_exp(y.raw(), x.raw(), prec);
    arb_neg(y.raw(), y.raw());
    fmpq_mul_2exp(exact.raw(), t.rho.raw(), 2);
    fmpq_add_si(exact.raw(), exact.raw(), 1);
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_add(y.raw(), y.raw(), x.raw(), prec);
    if (arb_is_positive(y.raw()) == 0) {
        return false;
    }
    arb_div(value.raw(), value.raw(), y.raw(), prec);
    if (arb_is_finite(value.raw()) == 0) {
        return false;
    }
    flint::Arf upper;
    arb_get_ubound_arf(upper.raw(), value.raw(), prec);
    arb_set_arf(t.m_phi.raw(), upper.raw());

    t.computed = -1;
    zimmert_tail_extend(t, 2);
    fmpq_mul(exact.raw(), t.c.raw(), t.c.raw());
    scale_rational(exact, -1, 4);
    return fmpq_is_zero(t.phi_big[0].raw()) != 0 &&
           fmpq_equal(t.phi_big[1].raw(), exact.raw()) != 0;
}

// Checks that N > 2 + gamma and N > (K - 1)/sigma (admissibility, for the
// N = 2(L + 1 + delta) given).
bool zimmert_tail_admissible(const ZimmertTailConstants& t,
                             slong n,
                             slong order) noexcept {
    flint::Fmpq bound;
    fmpq_add_si(bound.raw(), t.gamma.raw(), 2);
    if (fmpq_cmp_si(bound.raw(), n) >= 0) {
        return false;
    }
    fmpq_set_si(bound.raw(), order - 1, 1);
    fmpq_div(bound.raw(), bound.raw(), t.sigma.raw());
    return fmpq_cmp_si(bound.raw(), n) < 0;
}

// Sets out to the exact ball mu(N, K) = sum_{k=1}^{K-1} Phi_k k!/N^(k+1).
void zimmert_tail_mean(flint::Arb& out,
                       ZimmertTailConstants& t,
                       slong n,
                       slong order,
                       slong prec) noexcept {
    zimmert_tail_extend(t, order);
    flint::Fmpq sum;
    flint::Fmpq term;
    flint::Fmpz weight;
    flint::Fmpz power;
    for (slong k = 1; k < order; ++k) {
        fmpz_fac_ui(weight.raw(), static_cast<ulong>(k));
        fmpz_set_si(power.raw(), n);
        fmpz_pow_ui(power.raw(), power.raw(), static_cast<ulong>(k + 1));
        fmpq_mul_fmpz(term.raw(), t.phi_big[k].raw(), weight.raw());
        fmpq_div_fmpz(term.raw(), term.raw(), power.raw());
        fmpq_add(sum.raw(), sum.raw(), term.raw());
    }
    arb_set_fmpq(out.raw(), sum.raw(), prec);
}

// Sets out to a ball for E(N, K); (N, K) must be admissible.  The caller
// uses its upper endpoint.
void zimmert_tail_error(flint::Arb& out,
                        const ZimmertTailConstants& t,
                        slong n,
                        slong order,
                        slong prec) noexcept {
    flint::Fmpq exact;
    flint::Fmpq other;
    flint::Fmpz integer;
    flint::Arb x;
    flint::Arb y;
    flint::Arb total;
    // 2 M_Phi K!/(rho^K N^(K+1)).
    fmpz_fac_ui(integer.raw(), static_cast<ulong>(order));
    fmpq_one(exact.raw());
    fmpq_mul_fmpz(exact.raw(), exact.raw(), integer.raw());
    fmpq_pow_si(other.raw(), t.rho.raw(), order);
    fmpq_div(exact.raw(), exact.raw(), other.raw());
    fmpz_set_si(integer.raw(), n);
    fmpz_pow_ui(integer.raw(), integer.raw(), static_cast<ulong>(order + 1));
    fmpq_div_fmpz(exact.raw(), exact.raw(), integer.raw());
    fmpq_mul_2exp(exact.raw(), exact.raw(), 1);
    arb_set_fmpq(total.raw(), exact.raw(), prec);
    arb_mul(total.raw(), total.raw(), t.m_phi.raw(), prec);
    // 3 e^(-(N-2-gamma) sigma)/((N - 2 - gamma)(1 - e^(-2 sigma))).
    fmpq_set_si(exact.raw(), n - 2, 1);
    fmpq_sub(exact.raw(), exact.raw(), t.gamma.raw());  // N - 2 - gamma
    fmpq_mul(other.raw(), exact.raw(), t.sigma.raw());
    fmpq_neg(other.raw(), other.raw());
    arb_set_fmpq(x.raw(), other.raw(), prec);
    arb_exp(x.raw(), x.raw(), prec);
    arb_mul_ui(x.raw(), x.raw(), 3, prec);
    arb_set_fmpq(y.raw(), exact.raw(), prec);
    arb_div(x.raw(), x.raw(), y.raw(), prec);
    fmpq_mul_2exp(other.raw(), t.sigma.raw(), 1);
    fmpq_neg(other.raw(), other.raw());
    arb_set_fmpq(y.raw(), other.raw(), prec);
    arb_expm1(y.raw(), y.raw(), prec);
    arb_neg(y.raw(), y.raw());  // 1 - e^(-2 sigma)
    arb_div(x.raw(), x.raw(), y.raw(), prec);
    arb_add(total.raw(), total.raw(), x.raw(), prec);
    // M_Phi e^(-N sigma)/(N - (K - 1)/sigma).
    fmpq_mul_si(other.raw(), t.sigma.raw(), -n);
    arb_set_fmpq(x.raw(), other.raw(), prec);
    arb_exp(x.raw(), x.raw(), prec);
    arb_mul(x.raw(), x.raw(), t.m_phi.raw(), prec);
    fmpq_set_si(exact.raw(), order - 1, 1);
    fmpq_div(exact.raw(), exact.raw(), t.sigma.raw());
    fmpq_neg(exact.raw(), exact.raw());
    fmpq_add_si(exact.raw(), exact.raw(), n);
    arb_set_fmpq(y.raw(), exact.raw(), prec);
    arb_div(x.raw(), x.raw(), y.raw(), prec);
    arb_add(out.raw(), total.raw(), x.raw(), prec);
}

// True if the upper endpoint of E(N, K) is at most 2^-bits.
bool zimmert_tail_error_at_most(const ZimmertTailConstants& t,
                                slong n,
                                slong order,
                                slong bits,
                                slong prec) noexcept {
    flint::Arb error;
    flint::Arf upper;
    zimmert_tail_error(error, t, n, order, prec);
    if (arb_is_finite(error.raw()) == 0) {
        return false;
    }
    arb_get_ubound_arf(upper.raw(), error.raw(), prec);
    return arf_cmp_2exp_si(upper.raw(), -bits) <= 0;
}

// Chooses (L, K) as in the proof's "Choosing L and K": L >= max(1,
// floor(gamma/2) + 1) with (2L + 2) sigma >= (p + 5) log 2 + log(3 + M_Phi),
// then the least K >= 2 with K - 1 < (2L + 2) sigma and E(2L+2, K) <= 2^-p,
// increasing L if no such K <= kZimmertMaxTailOrder exists.  The bound is
// rechecked for the result, and (L, 0, K) is admissible.
bool zimmert_tail_parameters(slong& terms,
                             slong& order,
                             ZimmertTailConstants& t,
                             slong bits,
                             slong prec) noexcept {
    flint::Arb x;
    flint::Arb y;
    arb_add_ui(x.raw(), t.m_phi.raw(), 3, prec);
    arb_log(x.raw(), x.raw(), prec);
    arb_const_log2(y.raw(), prec);
    arb_mul_si(y.raw(), y.raw(), bits + 5, prec);
    arb_add(x.raw(), x.raw(), y.raw(), prec);
    arb_set_fmpq(y.raw(), t.sigma.raw(), prec);
    arb_mul_2exp_si(y.raw(), y.raw(), 1);
    arb_div(x.raw(), x.raw(), y.raw(), prec);  // L + 1 >= this
    flint::Arf upper;
    flint::Fmpz integer;
    arb_get_ubound_arf(upper.raw(), x.raw(), prec);
    if (arf_is_finite(upper.raw()) == 0 ||
        arf_cmp_si(upper.raw(), kZimmertMaxTailTerms) > 0) {
        return false;
    }
    arf_get_fmpz(integer.raw(), upper.raw(), ARF_RND_CEIL);
    slong l = fmpz_get_si(integer.raw()) - 1;
    const slong least = t.numerator / (2 * t.denominator) + 1;
    l = std::max(l, std::max(least, slong{1}));
    while (l <= kZimmertMaxTailTerms) {
        const slong n = 2 * l + 2;
        for (slong k = 2; k <= kZimmertMaxTailOrder; ++k) {
            if (!zimmert_tail_admissible(t, n, k)) {
                break;
            }
            if (zimmert_tail_error_at_most(t, n, k, bits, prec)) {
                terms = l;
                order = k;
                return true;
            }
        }
        l += l / 4 + 1;
    }
    return false;
}

// Sets out to psi(x) for an exact rational x = num/den > 0.
bool digamma_at_rational(flint::Arb& out,
                         slong num,
                         slong den,
                         slong prec) noexcept {
    if (num <= 0 || den <= 0) {
        return false;
    }
    flint::Fmpq exact;
    flint::Arb x;
    set_rational(exact, num, den);
    arb_set_fmpq(x.raw(), exact.raw(), prec);
    arb_digamma(out.raw(), x.raw(), prec);
    return arb_is_finite(out.raw()) != 0;
}

// Adds weight * (psi(x1) - psi(x0)) - 1/v to sum, for x0 = n0/den,
// x1 = n1/den and v = vn/vd, all exact positive rationals.
bool zimmert_head_term(flint::Arb& sum,
                       const flint::Arb& weight,
                       slong n1,
                       slong n0,
                       slong den,
                       slong vn,
                       slong vd,
                       slong prec) noexcept {
    flint::Arb x;
    flint::Arb y;
    if (vn <= 0 || !digamma_at_rational(x, n1, den, prec) ||
        !digamma_at_rational(y, n0, den, prec)) {
        return false;
    }
    arb_sub(x.raw(), x.raw(), y.raw(), prec);
    arb_mul(x.raw(), x.raw(), weight.raw(), prec);
    arb_add(sum.raw(), sum.raw(), x.raw(), prec);
    flint::Fmpq exact;
    set_rational(exact, vd, vn);
    arb_set_fmpq(y.raw(), exact.raw(), prec);
    arb_sub(sum.raw(), sum.raw(), y.raw(), prec);
    return true;
}

// Encloses F_{a,b}(gamma) for gamma = numerator/denominator > 0 by part
// (iv) of GB-zimmert-F-tail-bound and the evaluation procedure of its proof
// file (steps 1-5), with the tails bounded by 2^-kZimmertTailBits each.
// Sets terms and order to the L and K used.  If terms > 0 on entry, uses
// the given L and K instead, failing unless 2 <= K <= kZimmertMaxTailOrder,
// L <= kZimmertMaxTailTerms, L > gamma/2, and (L, 0, K) and (L, 1/2, K) are
// admissible.
bool zimmert_f_enclosure(flint::Arb& out,
                         slong& terms,
                         slong& order,
                         slong a,
                         slong b,
                         slong numerator,
                         slong denominator) noexcept {
    if (a < 1 || b < 0 || numerator <= 0 || denominator <= 0 ||
        numerator > kZimmertMaxTailTerms ||
        denominator > kZimmertMaxTailTerms) {
        return false;
    }
    const slong base_prec = kGenerationBoundPrecision + 10;
    ZimmertTailConstants t;
    if (!zimmert_tail_setup(t, numerator, denominator, base_prec)) {
        return false;
    }
    // Step 1: L and K, with (L, 0, K) and (L, 1/2, K) admissible.
    slong l = terms;
    slong k = order;
    if (l > 0) {
        if (k < 2 || k > kZimmertMaxTailOrder || l > kZimmertMaxTailTerms ||
            l <= t.numerator / (2 * t.denominator)) {
            return false;
        }
    } else if (!zimmert_tail_parameters(l, k, t, kZimmertTailBits,
                                        base_prec)) {
        return false;
    }
    if (!zimmert_tail_admissible(t, 2 * l + 2, k) ||
        !zimmert_tail_admissible(t, 2 * l + 3, k)) {
        return false;
    }
    // Working precision: target + 10 + log2 L (proof, after step 5).
    const slong prec = base_prec + static_cast<slong>(FLINT_BIT_COUNT(
                                           static_cast<ulong>(l)));
    // gamma = p/q; every psi argument below is an integer over hq = 2(q + 2p)
    // (that is, over h q) or over 2q.
    const slong p = t.numerator;
    const slong q = t.denominator;
    const slong hq = 2 * (q + 2 * p);
    flint::Arb two_over_c;  // 2/c = 2q/(q + 2p)
    {
        flint::Fmpq exact;
        set_rational(exact, 2 * q, q + 2 * p);
        arb_set_fmpq(two_over_c.raw(), exact.raw(), prec);
    }

    // Step 4: S_L(a, b; gamma).
    flint::Arb head_a;
    flint::Arb head_b;
    flint::Arb x;
    for (slong j = 1; j <= l; ++j) {
        // (2/c)(psi((2j+3g)/h) - psi((2j-1+g)/h)) - 1/(2j-1+g).
        if (!zimmert_head_term(head_a, two_over_c, 2 * j * q + 3 * p,
                               (2 * j - 1) * q + p, hq, (2 * j - 1) * q + p,
                               q, prec)) {
            return false;
        }
        if (b == 0) {
            continue;
        }
        // (2/c)(psi((2j+1+3g)/h) - psi((2j+g)/h)) - 1/(2j+g).
        if (!zimmert_head_term(head_b, two_over_c, (2 * j + 1) * q + 3 * p,
                               2 * j * q + p, hq, 2 * j * q + p, q, prec)) {
            return false;
        }
    }
    // a-part: - psi(1/2 + g/2)/2 - psi(L - g/2)/2.
    if (!digamma_at_rational(x, q + p, 2 * q, prec)) {
        return false;
    }
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_sub(head_a.raw(), head_a.raw(), x.raw(), prec);
    if (!digamma_at_rational(x, 2 * l * q - p, 2 * q, prec)) {
        return false;
    }
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_sub(head_a.raw(), head_a.raw(), x.raw(), prec);
    if (b != 0) {
        // b-part: - psi(1 + g/2)/2 - psi(L + 1/2 - g/2)/2.
        if (!digamma_at_rational(x, 2 * q + p, 2 * q, prec)) {
            return false;
        }
        arb_mul_2exp_si(x.raw(), x.raw(), -1);
        arb_sub(head_b.raw(), head_b.raw(), x.raw(), prec);
        if (!digamma_at_rational(x, (2 * l + 1) * q - p, 2 * q, prec)) {
            return false;
        }
        arb_mul_2exp_si(x.raw(), x.raw(), -1);
        arb_sub(head_b.raw(), head_b.raw(), x.raw(), prec);
    }
    // C(gamma) = -4/g + (2/c)(psi((1+g)/h) - psi((1+3g)/h)
    //                         + psi((2+5g)/h) - psi((2+3g)/h)).
    flint::Arb constant;
    flint::Arb y;
    if (!digamma_at_rational(constant, q + p, hq, prec) ||
        !digamma_at_rational(x, q + 3 * p, hq, prec)) {
        return false;
    }
    arb_sub(constant.raw(), constant.raw(), x.raw(), prec);
    if (!digamma_at_rational(x, 2 * q + 5 * p, hq, prec) ||
        !digamma_at_rational(y, 2 * q + 3 * p, hq, prec)) {
        return false;
    }
    arb_add(constant.raw(), constant.raw(), x.raw(), prec);
    arb_sub(constant.raw(), constant.raw(), y.raw(), prec);
    arb_mul(constant.raw(), constant.raw(), two_over_c.raw(), prec);
    {
        flint::Fmpq exact;
        set_rational(exact, 4 * q, p);
        arb_set_fmpq(x.raw(), exact.raw(), prec);
    }
    arb_sub(constant.raw(), constant.raw(), x.raw(), prec);

    // Steps 2, 3 and 5: add a mu(2L+2, K) + b mu(2L+3, K) and widen by
    // a E(2L+2, K) + b E(2L+3, K), from the upper endpoints of E.
    flint::Arb error;
    flint::Arb bound;
    flint::Arf upper;
    zimmert_tail_mean(x, t, 2 * l + 2, k, prec);
    arb_add(head_a.raw(), head_a.raw(), x.raw(), prec);
    zimmert_tail_error(bound, t, 2 * l + 2, k, prec);
    if (arb_is_finite(bound.raw()) == 0) {
        return false;
    }
    arb_get_ubound_arf(upper.raw(), bound.raw(), prec);
    arb_set_arf(bound.raw(), upper.raw());
    arb_mul_si(error.raw(), bound.raw(), a, prec);
    arb_mul_si(head_a.raw(), head_a.raw(), a, prec);
    if (b != 0) {
        zimmert_tail_mean(x, t, 2 * l + 3, k, prec);
        arb_add(head_b.raw(), head_b.raw(), x.raw(), prec);
        zimmert_tail_error(bound, t, 2 * l + 3, k, prec);
        if (arb_is_finite(bound.raw()) == 0) {
            return false;
        }
        arb_get_ubound_arf(upper.raw(), bound.raw(), prec);
        arb_set_arf(bound.raw(), upper.raw());
        arb_mul_si(bound.raw(), bound.raw(), b, prec);
        arb_add(error.raw(), error.raw(), bound.raw(), prec);
        arb_mul_si(head_b.raw(), head_b.raw(), b, prec);
        arb_add(head_a.raw(), head_a.raw(), head_b.raw(), prec);
    }
    arb_add(head_a.raw(), head_a.raw(), constant.raw(), prec);
    arb_add_error(head_a.raw(), error.raw());
    if (arb_is_finite(head_a.raw()) == 0) {
        return false;
    }
    arb_swap(out.raw(), head_a.raw());
    terms = l;
    order = k;
    return true;
}

// Encloses (F_{a,b}(g) + n log pi) / 2 with a = r1 + r2, b = r2,
// g = gamma_hundredths/100 and n = r1 + 2 r2, F_{a,b} enclosed by
// zimmert_f_enclosure.  Zimmert 1981, Korollar 1 (p. 376): for every ideal
// class K there are integral ideals a in K and b in D K^(-1) (D the class
// of the different) with N(a) N(b) <= d exp(-(F_{a,b}(g) + n log pi)).
// Fails for g outside (0, 3) or r1 + r2 < 1.  Not evaluated at run time:
// it generates and verifies kZimmertTwinClassTable below.
bool zimmert_twin_class_log_lower_bound(flint::Arb& out,
                                        slong& terms,
                                        slong& order,
                                        slong r1,
                                        slong r2,
                                        ulong gamma_hundredths) noexcept {
    const slong prec = kGenerationBoundPrecision;
    if (r1 < 0 || r2 < 0 || r1 + r2 < 1 || gamma_hundredths == 0 ||
        gamma_hundredths >= 300) {
        return false;
    }
    flint::Arb total;
    flint::Arb x;
    terms = 0;
    order = 0;
    if (!zimmert_f_enclosure(total, terms, order, r1 + r2, r2,
                             static_cast<slong>(gamma_hundredths), 100)) {
        return false;
    }
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


// Upper bound significand / 10^scale (an exact decimal rational) for
// exp(-c), c = (F_{a,b}(g) + n log pi)/2 at the tabulated gamma.  A null
// significand marks a signature without a table entry.
struct ZimmertTwinClassEntry {
    ulong gamma_hundredths;
    const char* significand;
    slong scale;
};

// kZimmertTwinClassTable[n - 3][r2]: precomputed, not evaluated at run time.
// Generated offline by tools/zimmert_twin_class_table.cpp (which says how to
// regenerate it) from zimmert_twin_class_log_lower_bound, and verified entry
// by entry against it by test/t-generation-bound.cpp.
#include "zimmert_twin_class_table.inc"

// Zimmert 1981, Korollar 1 (twin classes): for every ideal class K, K or
// D K^(-1) contains an integral ideal of norm at most sqrt|d| exp(-c), with
// c = (F_{a,b}(g) + n log pi)/2 at the tabulated gamma.  Sets out to ceil(E sqrt|d|) for the table entry E >= exp(-c),
// computed exactly in integers: with E = s / 10^k, ceil(s sqrt|d|) is the
// integer square root of s^2 |d| rounded up, and
// ceil(E sqrt|d|) = ceil(ceil(s sqrt|d|) / 10^k).  This is not a bound on
// the smallest ideal of every class; it is a generation bound by
// GB-kor1-generation (see proven_generation_requirement).
bool zimmert_twin_class_bound(flint::Fmpz& out,
                              flint::FmpzConstRef abs_discriminant,
                              slong degree,
                              slong real_places,
                              slong complex_pairs) noexcept {
    if (degree < 3 || degree > kZimmertMaxDegree || real_places < 0 ||
        complex_pairs < 0 || real_places + 2 * complex_pairs != degree ||
        fmpz_sgn(abs_discriminant.raw()) < 0) {
        return false;
    }
    const ZimmertTwinClassEntry& entry =
        kZimmertTwinClassTable[degree - 3][complex_pairs];
    flint::Fmpz significand;
    if (entry.significand == nullptr || entry.scale < 0 ||
        fmpz_set_str(significand.raw(), entry.significand, 10) != 0 ||
        fmpz_sgn(significand.raw()) <= 0) {
        return false;
    }
    flint::Fmpz square;
    flint::Fmpz root;
    flint::Fmpz power;
    fmpz_mul(square.raw(), significand.raw(), significand.raw());
    fmpz_mul(square.raw(), square.raw(), abs_discriminant.raw());
    ceil_sqrt(root, flint::FmpzConstRef(square));
    fmpz_ui_pow_ui(power.raw(), 10, static_cast<ulong>(entry.scale));
    fmpz_cdiv_q(out.raw(), root.raw(), power.raw());
    return true;
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

// Two primes of one order with the same rational prime are the same prime
// ideal exactly when their integral ideals are equal; e and f are then equal
// too.  Ideal::equal compares HNFs, which are canonical, so this does not
// depend on the residue data a construction recorded (PrimeIdeal::equal also
// compares the residue polynomial, its variable convention and whether one
// was recorded at all).
bool same_prime_ideal(const PrimeIdeal& left,
                      const PrimeIdeal& right) noexcept {
    return detail::PrimeIdealAccess::ideal(left).equal(
            detail::PrimeIdealAccess::ideal(right));
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
        ulong gamma_hundredths,
        slong* terms,
        slong* order) noexcept {
    flint::Arb value;
    slong l = 0;
    slong k = 0;
    if (!silex::zimmert_twin_class_log_lower_bound(
                value, l, k, real_places, complex_pairs, gamma_hundredths)) {
        return false;
    }
    arb_swap(out.raw(), value.raw());
    if (terms != nullptr) {
        *terms = l;
    }
    if (order != nullptr) {
        *order = k;
    }
    return true;
}

bool detail::zimmert_f_enclosure(flint::ArbRef out,
                                 slong& terms,
                                 slong& order,
                                 slong a,
                                 slong b,
                                 slong gamma_numerator,
                                 slong gamma_denominator) noexcept {
    flint::Arb value;
    slong l = terms;
    slong k = order;
    if (!silex::zimmert_f_enclosure(value, l, k, a, b, gamma_numerator,
                                    gamma_denominator)) {
        return false;
    }
    arb_swap(out.raw(), value.raw());
    terms = l;
    order = k;
    return true;
}

ulong detail::zimmert_twin_class_gamma_hundredths(
        slong real_places,
        slong complex_pairs) noexcept {
    const slong degree = real_places + 2 * complex_pairs;
    if (real_places < 0 || complex_pairs < 0 || degree < 3 ||
        degree > kZimmertMaxDegree) {
        return 0;
    }
    return kZimmertTwinClassGammaHundredths[degree - 3][complex_pairs];
}

bool detail::zimmert_twin_class_table_entry(flint::FmpzRef significand,
                                            slong& scale,
                                            ulong& gamma_hundredths,
                                            slong real_places,
                                            slong complex_pairs) noexcept {
    const slong degree = real_places + 2 * complex_pairs;
    if (real_places < 0 || complex_pairs < 0 || degree < 3 ||
        degree > kZimmertMaxDegree) {
        return false;
    }
    const ZimmertTwinClassEntry& entry =
        kZimmertTwinClassTable[degree - 3][complex_pairs];
    flint::Fmpz value;
    if (entry.significand == nullptr ||
        fmpz_set_str(value.raw(), entry.significand, 10) != 0) {
        return false;
    }
    fmpz_swap(significand.raw(), value.raw());
    scale = entry.scale;
    gamma_hundredths = entry.gamma_hundredths;
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
    std::swap(identity_, other.identity_);
}

void FactorBase::clear() noexcept {
    primes_.clear();
    blocks_.clear();
    parent_.clear();
    complete_rational_prime_blocks_ = false;
    identity_ = 0;
}

bool FactorBase::define(const Order& parent) noexcept {
    if (!parent.has_basis()) {
        return false;
    }

    clear();
    parent_ = parent;
    identity_ = fresh_factor_base_identity();
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

    // The copy holds the same parent, primes and blocks as other, so it
    // carries other's identity until either side changes.
    swap(copy);
    identity_ = other.identity_;
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
    return detail::FactorBaseBlockAccess::order_lll_relation_base_by_norm(
            *this);
}

bool detail::FactorBaseBlockAccess::append_prime_for_testing(
        FactorBase& base,
        const PrimeIdeal& prime) noexcept {
    return base.append_prime(prime);
}

bool detail::FactorBaseBlockAccess::order_lll_relation_base_by_norm(
        FactorBase& base) noexcept {
    struct NormIndex {
        slong index = 0;
        flint::Fmpz norm;
    };

    std::vector<NormIndex> order;
    order.reserve(base.primes_.size());
    for (slong i = 0; i < static_cast<slong>(base.primes_.size()); ++i) {
        NormIndex entry;
        entry.index = i;
        if (!base.primes_[static_cast<std::size_t>(i)].norm(
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

    base.identity_ = fresh_factor_base_identity();
    std::vector<PrimeIdeal> sorted;
    sorted.reserve(base.primes_.size());
    for (const NormIndex& entry : order) {
        PrimeIdeal copy(base.parent_);
        if (!copy.is_defined() ||
            !copy.set(base.primes_[static_cast<std::size_t>(entry.index)])) {
            return false;
        }
        sorted.emplace_back(std::move(copy));
    }

    base.primes_.swap(sorted);
    return base.rebuild_blocks();
}

bool FactorBase::build_prime_ideal_norm_bounded_impl(
        flint::FmpzConstRef bound,
        bool include_inert_primes,
        slong degree_limit,
        bool retain_one_split_prime) noexcept {
    return detail::FactorBaseBlockAccess::build_prime_ideal_norm_bounded(
            *this, bound, include_inert_primes, degree_limit,
            retain_one_split_prime, nullptr, nullptr);
}

bool detail::FactorBaseBlockAccess::build_prime_ideal_norm_bounded_cached(
        FactorBase& base,
        flint::FmpzConstRef bound,
        PrimeDecompositionCache& cache,
        const DiagnosticsContext* diagnostics) noexcept {
    return build_prime_ideal_norm_bounded(base, bound, true, 0, false,
                                          &cache, diagnostics);
}

bool detail::FactorBaseBlockAccess::build_prime_ideal_norm_bounded(
        FactorBase& base,
        flint::FmpzConstRef bound,
        bool include_inert_primes,
        slong degree_limit,
        bool retain_one_split_prime,
        PrimeDecompositionCache* cache,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "factor_base.norm_bounded.build");
    const Order& parent = base.parent_;
    if (!base.is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0 ||
        !parent.is_maximal() || degree_limit < 0) {
        return false;
    }

    FactorBase candidate(parent);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    PrimeIdealList local;
    std::size_t prime_position = 0;
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        // Entry k of the cache is the decomposition of the k-th rational
        // prime, so a build at a larger bound reuses every decomposition an
        // earlier build of the same order computed and decomposes only the
        // primes beyond it.
        const PrimeIdealList* decomposition = nullptr;
        if (cache != nullptr &&
            prime_position < cache->decompositions.size()) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "factor_base.norm_bounded.decomposition_cache_hit");
            decomposition = &cache->decompositions[prime_position];
            const PrimeIdeal* first = decomposition->at(0);
            if (first == nullptr ||
                !same_order_parent(first->parent(), &parent)) {
                return false;
            }
        } else {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "factor_base.norm_bounded.decompose_prime");
            PrimeIdealList computed;
            if (!decompose_prime(computed, parent,
                                 flint::FmpzConstRef(p))) {
                return false;
            }
            if (cache != nullptr) {
                cache->decompositions.emplace_back(std::move(computed));
                decomposition = &cache->decompositions.back();
            } else {
                local.swap(computed);
                decomposition = &local;
            }
        }
        ++prime_position;
        const PrimeIdealList& decomposed = *decomposition;
        // The source factor-base generator skips inert rational primes for the relation
        // factor base, even when p^degree is still below the tested bound.
        if (!include_inert_primes &&
            is_inert_prime_decomposition(decomposed, parent.degree())) {
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

    base.swap(candidate);
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
    // A shared identity means one base is an unchanged copy of the other
    // (or the same object), whose contents compare equal below.
    if (identity_ != 0 && identity_ == other.identity_ && is_defined() &&
        other.is_defined()) {
        return true;
    }
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

    // Every prime in the block lies over rational_prime, so the ideal
    // comparison decides equality.
    const PrimeBlock& prime_block = blocks_[static_cast<std::size_t>(block)];
    for (slong index : prime_block.indices) {
        if (same_prime_ideal(primes_[static_cast<std::size_t>(index)],
                             prime)) {
            return index;
        }
    }
    return -1;
}

// The base holds each prime ideal at most once, and append_prime fails
// without changing the base when asked to repeat one.  The completeness
// checks rely on this: a block is complete when the sum of e f over its
// primes is the degree (ideal_factor_over_base), and the builders mark a
// block complete when it retains every prime of the decomposition, which is
// only sound for distinct primes.  The builders cannot repeat a prime: each
// rational prime is decomposed once, decompose_prime lists the distinct
// primes above it, and set() copies a base that already holds distinct
// primes.
bool FactorBase::append_prime(const PrimeIdeal& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data() || index(prime) >= 0) {
        return false;
    }

    PrimeIdeal copy(parent_);
    flint::Fmpz rational_prime;
    if (!copy.is_defined() || !copy.set(prime) ||
        !copy.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    identity_ = fresh_factor_base_identity();
    const slong index = static_cast<slong>(primes_.size());
    primes_.emplace_back(std::move(copy));
    return add_block_index(flint::FmpzConstRef(rational_prime), index);
}

bool FactorBase::append_prime(PrimeIdeal&& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data() || index(prime) >= 0) {
        return false;
    }

    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    identity_ = fresh_factor_base_identity();
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

    identity_ = fresh_factor_base_identity();
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
    identity_ = fresh_factor_base_identity();
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

// Zimmert 1981, Korollar 1 gives, for every class K, an integral ideal of
// norm at most B(gamma) in K or in D K^(-1), D the class of the different.
// GB-kor1-generation (silex-research,
// research/proofs/generation/GB-kor1-generation.md at commit 8926120) proves
// that the prime ideals of norm at most B(gamma) generate the class group
// with no hypothesis on [D], for every number field and every gamma > 0, and
// that any certified B' >= B(gamma) may replace B(gamma).  The twin-class
// bound ceil(E sqrt|d|) is such a B': every table entry E is at least the
// upper endpoint of Silex's own certified enclosure of
// exp(-(F + n log pi)/2) at the tabulated gamma
// (GB-zimmert-F-tail-bound, research/proofs/generation/ at commit 06191ff),
// and at least the independent c_up of silex-research's
// src/silex_research/zimmert.py (checked for all 117 entries by
// tools/check_zimmert_twin_class_table.py).  So the
// requirement is the twin-class bound whenever it is below
// factor_base_class_group_bound, whatever the defining polynomial.
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
            fmpz_swap(chosen.raw(), twin.raw());
            chosen_route = GenerationDifferentRoute::kor1_generation;
        }
    }

    fmpz_set(bound.raw(), chosen.raw());
    route = chosen_route;
    return true;
}

}  // namespace silex
