#include <silex/class_group.hpp>

#include <silex/order_element.hpp>
#include <silex/relation.hpp>

#include "relation_admission_cache_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"
#include "../relation/relation_internal.hpp"
#include "../unit/unit_internal.hpp"

#include <silex/flint/acb_vec.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/nmod_poly_factor.hpp>

#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/mag.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {

namespace {

// Measured speed threshold, not from a source, and revisable: below this
// degree the per-row cost of the exact prefilter resultants is a few
// milliseconds on the grh corpus, so the split-prime search does not pay
// for itself; those fields keep fmpz_poly_resultant only.  Either path
// gives the same exact resultant, so the value never depends on it.
constexpr slong kSplitPrimeResultantMinDegree = 16;

// Split primes are taken from [2^27, 2^28).  For p < 2^28,
// 255 * (p - 1)^2 < 2^64, so up to 255 products of reduced residues
// accumulate in one word before a reduction.
constexpr ulong kSplitPrimeSearchStart = UWORD(1) << 27;
constexpr ulong kSplitPrimeSearchEnd = UWORD(1) << 28;
constexpr slong kSplitPrimeDotChunk = 255;

// By Chebotarev, primes splitting completely in f have density
// 1/|Gal(f)| <= 1/deg(f).  The search stops once
// kSplitPrimeSearchFactor * deg(f) primes per split prime found (plus one)
// have been tested, which bounds the work spent on fields where such
// primes are rare.  Stopping only routes candidates to the exact
// fmpz_poly_resultant.
constexpr slong kSplitPrimeSearchFactor = 4;
// Measured speed threshold, not from a source, and revisable: at most this
// many split primes are kept; candidates whose bound needs more fall back
// to the exact fmpz_poly_resultant.
constexpr std::size_t kMaxSplitPrimes = 64;

// Precisions of the root enclosures and of the Gram matrix and Cholesky
// checks behind the T2 bound.  Numerical choices, not from a source: when
// they are too low for a field, the Cholesky checks fail and the bound
// stays at trace(G), which is still certified.
constexpr slong kT2BoundRootPrecision = 32;
constexpr slong kT2BoundPrecision = 64;
// Power-iteration steps for the estimate of lambda_max(G) and the relative
// margin 2^-7 added to it before the Cholesky check; the search stops once
// the certified upper end is within a factor 1 + 2^-6 of the lower end, or
// after kT2BoundBisectionSteps checks.  A factor 1 + 2^-6 in lambda costs at
// most n / 2 * log2(1 + 2^-6) < n / 87 bits of the resultant bound.
constexpr int kT2BoundPowerSteps = 32;
constexpr slong kT2BoundMarginExponent = -7;
constexpr slong kT2BoundStopExponent = -6;
constexpr int kT2BoundBisectionSteps = 12;

// Approximates lambda_max(G) 2^-scale_exponent by power iteration in double
// precision on the midpoints of G 2^-scale_exponent (Golub and Van Loan,
// Matrix Computations, Section 8.2.1: the Rayleigh quotient of the iterate
// converges to the largest eigenvalue of a symmetric positive semidefinite
// matrix).  Only an estimate: it picks the first lambda to certify.
double t2_gram_power_estimate(const arb_mat_t gram,
                              slong degree,
                              slong scale_exponent) noexcept {
    const auto n = static_cast<std::size_t>(degree);
    std::vector<double> matrix(n * n);
    flint::Arf scaled;
    for (slong j = 0; j < degree; ++j) {
        for (slong k = 0; k < degree; ++k) {
            arf_mul_2exp_si(scaled.raw(),
                            arb_midref(arb_mat_entry(gram, j, k)),
                            -scale_exponent);
            matrix[static_cast<std::size_t>(j) * n +
                   static_cast<std::size_t>(k)] =
                    arf_get_d(scaled.raw(), ARF_RND_NEAR);
        }
    }
    // A fixed pseudo-random start, so the estimate is deterministic and the
    // start is not orthogonal to the top eigenspace in structured cases
    // (the all-ones vector is an eigenvector of Phi_p's G for eigenvalue 1).
    std::vector<double> vector(n);
    std::vector<double> image(n);
    ulong state = 1;
    for (std::size_t j = 0; j < n; ++j) {
        state = state * UWORD(6364136223846793005) +
                UWORD(1442695040888963407);
        vector[j] = 0.5 + static_cast<double>(state >> 11) * 0x1p-53;
    }
    double estimate = 0.0;
    for (int step = 0; step < kT2BoundPowerSteps; ++step) {
        double vector_norm = 0.0;
        double rayleigh = 0.0;
        double image_norm = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                sum += matrix[j * n + k] * vector[k];
            }
            image[j] = sum;
            vector_norm += vector[j] * vector[j];
            rayleigh += vector[j] * sum;
            image_norm += sum * sum;
        }
        if (!(vector_norm > 0.0) || !(image_norm > 0.0)) {
            return 0.0;
        }
        estimate = rayleigh / vector_norm;
        const double inverse = 1.0 / std::sqrt(image_norm);
        for (std::size_t j = 0; j < n; ++j) {
            vector[j] = image[j] * inverse;
        }
    }
    return estimate;
}

// Returns whether lambda I - G is certainly positive definite, which gives
// lambda_max(G) < lambda for the symmetric matrix G inside the balls.
bool t2_gram_bound_certified(flint::ArbMat& shifted,
                             flint::ArbMat& cholesky,
                             const arb_mat_t gram,
                             const arf_t lambda,
                             slong degree) noexcept {
    arb_mat_neg(shifted.raw(), gram);
    for (slong j = 0; j < degree; ++j) {
        arb_add_arf(arb_mat_entry(shifted.raw(), j, j),
                    arb_mat_entry(shifted.raw(), j, j), lambda,
                    kT2BoundPrecision);
    }
    return arb_mat_cho(cholesky.raw(), shifted.raw(), kT2BoundPrecision) != 0;
}

// Sets out to a certified upper bound for lambda_max(G) / n, where
// G[j][k] = Re sum_i r_i^j conj(r_i)^k over the n complex roots r_i of the
// squarefree f, so that T2(a) = a^T G a for the coefficient vector of a.
// G = A^T A + B^T B with A = Re(r_i^j) and B = Im(r_i^j) is positive
// semidefinite, so lambda_max(G) <= trace(G); each smaller lambda for which
// FLINT arb_mat_cho certifies lambda I - G positive definite gives
// lambda_max(G) < lambda.  The first lambda tried is the power-iteration
// estimate plus a margin; a geometric bisection between the lower end (the
// largest diagonal entry or the estimate) and the certified upper end
// follows only while they are far apart.  Returns false when no finite
// bound is found.
bool certified_t2_ratio_bound(flint::Arf& out,
                              const fmpz_poly_struct* field_polynomial,
                              slong degree) noexcept {
    const slong prec = kT2BoundPrecision;
    if (degree <= 0 || fmpz_poly_is_squarefree(field_polynomial) == 0) {
        return false;
    }
    // lambda_max(G) >= G[j][j] = sum_i |r_i|^(2j) >= |p_(2j)|, the exact
    // power sums of the monic f.  When some |p_(2j)| >= n ||f||_2^2, then
    // lambda / n >= ||f||_2^2 and the real-valued T2 bound is never below
    // the real-valued Hadamard bound, so the roots are not computed.  (The
    // rounded bit counts can still differ by a few bits; that costs speed
    // only.)
    {
        flint::FmpzPoly power_sums;
        flint::Fmpz field_norm_squared;
        flint::Fmpz threshold;
        fmpz_poly_power_sums(power_sums.raw(), field_polynomial,
                             2 * degree - 1);
        _fmpz_vec_dot(field_norm_squared.raw(), field_polynomial->coeffs,
                      field_polynomial->coeffs, field_polynomial->length);
        fmpz_mul_si(threshold.raw(), field_norm_squared.raw(), degree);
        for (slong j = 0; 2 * j < power_sums.raw()->length; ++j) {
            if (fmpz_cmpabs(power_sums.raw()->coeffs + 2 * j,
                            threshold.raw()) >= 0) {
                return false;
            }
        }
    }
    flint::AcbVec roots(degree);
    arb_fmpz_poly_complex_roots(roots.data(), field_polynomial, 0,
                                kT2BoundRootPrecision);

    flint::ArbMat real_powers(degree, degree);
    flint::ArbMat imag_powers(degree, degree);
    {
        acb_t power;
        acb_init(power);
        for (slong i = 0; i < degree; ++i) {
            acb_one(power);
            for (slong j = 0; j < degree; ++j) {
                arb_set(arb_mat_entry(real_powers.raw(), i, j),
                        acb_realref(power));
                arb_set(arb_mat_entry(imag_powers.raw(), i, j),
                        acb_imagref(power));
                acb_mul(power, power, roots.data() + i, prec);
            }
        }
        acb_clear(power);
    }
    flint::ArbMat transpose(degree, degree);
    flint::ArbMat gram(degree, degree);
    flint::ArbMat term(degree, degree);
    arb_mat_transpose(transpose.raw(), real_powers.raw());
    arb_mat_mul(gram.raw(), transpose.raw(), real_powers.raw(), prec);
    arb_mat_transpose(transpose.raw(), imag_powers.raw());
    arb_mat_mul(term.raw(), transpose.raw(), imag_powers.raw(), prec);
    arb_mat_add(gram.raw(), gram.raw(), term.raw(), prec);

    flint::Arb trace;
    flint::Arf high;
    flint::Arf low;
    flint::Arf scratch;
    arb_mat_trace(trace.raw(), gram.raw(), prec);
    arb_get_ubound_arf(high.raw(), trace.raw(), prec);
    if (arf_is_finite(high.raw()) == 0 || arf_sgn(high.raw()) <= 0) {
        return false;
    }
    arf_zero(low.raw());
    for (slong j = 0; j < degree; ++j) {
        arb_get_lbound_arf(scratch.raw(), arb_mat_entry(gram.raw(), j, j),
                           prec);
        if (arf_cmp(scratch.raw(), low.raw()) > 0) {
            arf_set(low.raw(), scratch.raw());
        }
    }

    flint::ArbMat shifted(degree, degree);
    flint::ArbMat cholesky(degree, degree);
    if (fmpz_fits_si(ARF_EXPREF(high.raw())) != 0) {
        const slong scale_exponent = fmpz_get_si(ARF_EXPREF(high.raw()));
        const double estimate =
                t2_gram_power_estimate(gram.raw(), degree, scale_exponent);
        if (estimate > 0.0 && estimate <= 1.0) {
            flint::Arf candidate;
            arf_set_d(scratch.raw(), estimate);
            arf_mul_2exp_si(scratch.raw(), scratch.raw(), scale_exponent);
            arf_mul_2exp_si(candidate.raw(), scratch.raw(),
                            kT2BoundMarginExponent);
            arf_add(candidate.raw(), candidate.raw(), scratch.raw(), prec,
                    ARF_RND_UP);
            if (arf_cmp(candidate.raw(), high.raw()) < 0) {
                if (t2_gram_bound_certified(shifted, cholesky, gram.raw(),
                                            candidate.raw(), degree)) {
                    arf_set(high.raw(), candidate.raw());
                    // The estimate only steers the search below.
                    if (arf_cmp(scratch.raw(), low.raw()) > 0) {
                        arf_set(low.raw(), scratch.raw());
                    }
                } else if (arf_cmp(candidate.raw(), low.raw()) > 0) {
                    arf_set(low.raw(), candidate.raw());
                }
            }
        }
    }

    if (arf_sgn(low.raw()) > 0 && arf_is_finite(low.raw()) != 0) {
        flint::Arf middle;
        for (int step = 0; step < kT2BoundBisectionSteps; ++step) {
            arf_mul_2exp_si(scratch.raw(), low.raw(), kT2BoundStopExponent);
            arf_add(scratch.raw(), scratch.raw(), low.raw(), prec,
                    ARF_RND_UP);
            if (arf_cmp(high.raw(), scratch.raw()) <= 0) {
                break;
            }
            arf_mul(middle.raw(), low.raw(), high.raw(), prec, ARF_RND_NEAR);
            arf_sqrt(middle.raw(), middle.raw(), prec, ARF_RND_NEAR);
            if (t2_gram_bound_certified(shifted, cholesky, gram.raw(),
                                        middle.raw(), degree)) {
                arf_set(high.raw(), middle.raw());
            } else {
                arf_set(low.raw(), middle.raw());
            }
        }
    }

    arf_div_ui(out.raw(), high.raw(), static_cast<ulong>(degree), prec,
               ARF_RND_UP);
    return true;
}

// Sets out so that |Res(f, a)| < 2^out, from
// |Res(f, a)|^2 <= (ratio ||a||_2^2)^n with ratio >= lambda_max(G) / n.
// mag_t arithmetic rounds every step upward.  Returns false when the bound
// does not fit a word.
bool t2_resultant_half_bits(flint_bitcnt_t& out,
                            const fmpz_t norm_squared,
                            const arf_t ratio,
                            slong degree) noexcept {
    mag_t bound;
    mag_t factor;
    mag_init(bound);
    mag_init(factor);
    mag_set_fmpz(bound, norm_squared);
    arf_get_mag(factor, ratio);
    mag_mul(bound, bound, factor);
    mag_pow_ui(bound, bound, static_cast<ulong>(degree));
    bool ok = true;
    if (mag_is_zero(bound) != 0) {
        out = 0;
    } else if (mag_is_finite(bound) == 0 ||
               fmpz_fits_si(MAG_EXPREF(bound)) == 0) {
        ok = false;
    } else {
        // bound < 2^exponent: a mag_t is m 2^(e - MAG_BITS) with
        // m < 2^MAG_BITS; the loop only guards that representation.
        slong exponent = fmpz_get_si(MAG_EXPREF(bound));
        while (mag_cmp_2exp_si(bound, exponent) >= 0) {
            ++exponent;
        }
        // |Res|^2 < 2^exponent, so |Res| < 2^ceil(exponent / 2).
        out = exponent <= 0 ? 0
                            : static_cast<flint_bitcnt_t>((exponent + 1) / 2);
    }
    mag_clear(bound);
    mag_clear(factor);
    return ok;
}

// Not from a source, and revisable: the bases a tried for a primitive n-th
// root before a closed-form split prime falls back to the generic test.
// Each try succeeds with probability about phi(n) / (2n) or more, and a
// fallback finds the same split prime, so the cap affects speed only.
constexpr ulong kClosedFormRootAttempts = 64;

// f(x) mod p by Horner's rule on the reduced coefficients c_0, ..., c_n.
ulong evaluate_reduced(const std::vector<ulong>& coefficients,
                       ulong x,
                       nmod_t mod) noexcept {
    ulong value = 0;
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        value = nmod_add(nmod_mul(value, x, mod), coefficients[i - 1], mod);
    }
    return value;
}

// V_e(a) mod p for the Lucas sequence V_0 = 2, V_1 = a,
// V_(k+1) = a V_k - V_(k-1) (P = a, Q = 1), so V_k = alpha^k + alpha^-k
// for a root alpha of x^2 - a x + 1, by the binary ladder on
// (V_k, V_(k+1)) with V_2k = V_k^2 - 2 and V_(2k+1) = V_k V_(k+1) - a
// (Crandall and Pomerance, Prime Numbers: A Computational Perspective,
// 2nd ed., 2005, Section 3.6.1).  e >= 1 and p odd.
ulong lucas_v(ulong a, ulong e, nmod_t mod) noexcept {
    const ulong two = nmod_set_ui(2, mod);
    ulong low = two;
    ulong high = a;
    for (slong bit = static_cast<slong>(FLINT_BIT_COUNT(e)) - 1; bit >= 0;
         --bit) {
        const ulong cross = nmod_sub(nmod_mul(low, high, mod), a, mod);
        if (((e >> bit) & 1) != 0) {
            low = cross;
            high = nmod_sub(nmod_mul(high, high, mod), two, mod);
        } else {
            high = cross;
            low = nmod_sub(nmod_mul(low, low, mod), two, mod);
        }
    }
    return low;
}

}  // namespace

void SplitPrimeResultant::reset(
        const fmpz_poly_struct* monic_field_polynomial) noexcept {
    // Rebuilding the prefilter for the same f keeps the split primes
    // already found and the search state, so the search is not repeated.
    if (enabled && monic_field_polynomial != nullptr &&
        fmpz_poly_equal(field_polynomial.raw(), monic_field_polynomial) !=
                0) {
        return;
    }
    enabled = false;
    search_exhausted = false;
    degree = 0;
    field_norm_squared_bits = 0;
    last_tested_prime = kSplitPrimeSearchStart;
    tested_primes = 0;
    last_prime_count = 0;
    t2_bound_state = T2BoundState::not_computed;
    discriminant_computed = false;
    split_closed_form = SplitClosedForm::not_checked;
    split_closed_form_index = 0;
    closed_form_split_primes = 0;
    moduli.clear();
    root_powers.clear();
    residues.clear();
    // The row values belong to the split primes just dropped.
    for (std::vector<ulong>& values : row_values) {
        values.clear();
    }
    if (monic_field_polynomial == nullptr) {
        return;
    }

    const slong field_degree = fmpz_poly_degree(monic_field_polynomial);
    if (field_degree < kSplitPrimeResultantMinDegree ||
        fmpz_is_one(fmpz_poly_lead(monic_field_polynomial)) == 0) {
        return;
    }

    fmpz_poly_set(field_polynomial.raw(), monic_field_polynomial);
    _fmpz_vec_dot(coordinate_norm_squared.raw(),
                  monic_field_polynomial->coeffs,
                  monic_field_polynomial->coeffs,
                  monic_field_polynomial->length);
    field_norm_squared_bits = fmpz_bits(coordinate_norm_squared.raw());
    degree = field_degree;
    residues.assign(static_cast<std::size_t>(field_degree), 0);
    enabled = true;
}

bool SplitPrimeResultant::add_next_split_prime() noexcept {
    const auto n = static_cast<std::size_t>(degree);
    while (!search_exhausted) {
        const auto found = static_cast<slong>(moduli.size());
        if (moduli.size() >= kMaxSplitPrimes ||
            tested_primes >=
                    kSplitPrimeSearchFactor * degree * (found + 1)) {
            search_exhausted = true;
            break;
        }
        const ulong p = n_nextprime(last_tested_prime, 1);
        if (p >= kSplitPrimeSearchEnd) {
            search_exhausted = true;
            break;
        }
        last_tested_prime = p;
        ++tested_primes;

        // Closed form for f = Phi_n or f = psi_n, n < 2^27 <= p, so p does
        // not divide n.  Modulo such p, the roots of Phi_n are the elements
        // of exact order n (the definition of the cyclotomic polynomial
        // Q_n), which all lie in F_p exactly when n | p - 1, so f splits
        // into distinct linear factors exactly when p = 1 mod n (Lidl and Niederreiter, Finite Fields, 2nd ed.,
        // 1997, Theorem 2.47(ii); Washington, Introduction to Cyclotomic
        // Fields, 2nd ed., 1997, Theorem 2.13).  Since
        // x^d psi_n(x + 1/x) = Phi_n(x), d = phi(n) / 2, the roots of psi_n are zeta^k + zeta^-k for zeta
        // of exact order n and k coprime to n, 0 < k < n / 2, distinct for
        // n >= 3.  zeta + zeta^-1 lies in F_p exactly when zeta^p is zeta
        // or zeta^-1, that is when n | p - 1 or n | p + 1, so psi_n splits
        // exactly when p = +-1 mod n.  Other primes fail the generic test
        // below, so skipping them leaves the split primes found, their
        // order and tested_primes unchanged.  The roots are built from a
        // primitive root of unity and checked exactly; if that fails, the
        // generic test decides the prime.
        if (split_closed_form == SplitClosedForm::not_checked) {
            split_closed_form = SplitClosedForm::none;
            const ulong cyclotomic =
                    fmpz_poly_is_cyclotomic(field_polynomial.raw());
            const ulong real_cyclotomic =
                    cyclotomic != 0
                            ? 0
                            : detail::real_cyclotomic_polynomial_index(
                                      field_polynomial.raw());
            if (cyclotomic != 0 && cyclotomic < kSplitPrimeSearchStart) {
                split_closed_form = SplitClosedForm::cyclotomic;
                split_closed_form_index = cyclotomic;
            } else if (real_cyclotomic != 0 &&
                       real_cyclotomic < kSplitPrimeSearchStart) {
                split_closed_form = SplitClosedForm::real_cyclotomic;
                split_closed_form_index = real_cyclotomic;
            }
        }
        if (split_closed_form == SplitClosedForm::cyclotomic ||
            split_closed_form == SplitClosedForm::real_cyclotomic) {
            const ulong residue = p % split_closed_form_index;
            if (residue != 1 &&
                (split_closed_form != SplitClosedForm::real_cyclotomic ||
                 residue != split_closed_form_index - 1)) {
                continue;
            }
            nmod_t mod;
            nmod_init(&mod, p);
            if (add_closed_form_roots(mod)) {
                moduli.push_back(mod);
                ++closed_form_split_primes;
                return true;
            }
        }

        // f splits into distinct linear factors modulo p exactly when
        // f divides x^p - x = prod_{a in F_p} (x - a) in F_p[x] (Fermat).  The
        // roots a_1, ..., a_n of f mod p are then distinct elements of F_p and
        // disc(f) mod p = prod_{i<j} (a_i - a_j)^2 is a nonzero square, so
        // the Jacobi symbol (disc(f) / p) is +1.  Here disc(f) is the
        // discriminant of the monic integer polynomial f itself (FLINT's
        // fmpz_poly_discriminant: (-1)^(n(n-1)/2) Res(f, f') / lc(f) =
        // lc(f)^(2n-2) prod_{i<j} (alpha_i - alpha_j)^2 over the roots
        // alpha_i of f), not the field discriminant; it is a polynomial in
        // the coefficients, so for monic f its reduction mod p is the
        // discriminant of f mod p.  Every odd prime p with
        // (disc(f) / p) = -1, or with p | disc(f) (symbol 0: f mod p has a
        // repeated root, whatever the reason, ramification or an index
        // divisor), therefore fails the test below, and skipping it leaves
        // the split primes found and their order unchanged.  The prime
        // still counts in tested_primes, so the search stops where it did
        // without the screen.  The screen only saves work when disc(f) is
        // not a square: then about half of the tested primes have symbol
        // -1.  The tested primes are odd, as the Jacobi symbol requires.
        if (!discriminant_computed) {
            fmpz_poly_discriminant(discriminant.raw(),
                                   field_polynomial.raw());
            discriminant_computed = true;
        }
        if (n_jacobi_unsigned(fmpz_fdiv_ui(discriminant.raw(), p), p) != 1) {
            continue;
        }

        // deg(f) >= 2 here, so x mod f = x.
        flint::NmodPoly reduced(p);
        flint::NmodPoly reversed_inverse(p);
        flint::NmodPoly frobenius(p);
        fmpz_poly_get_nmod_poly(reduced.raw(), field_polynomial.raw());
        nmod_poly_reverse(reversed_inverse.raw(), reduced.raw(), degree + 1);
        nmod_poly_inv_series(reversed_inverse.raw(), reversed_inverse.raw(),
                             degree + 1);
        nmod_poly_powmod_x_ui_preinv(frobenius.raw(), p, reduced.raw(),
                                     reversed_inverse.raw());
        if (nmod_poly_length(frobenius.raw()) != 2 ||
            nmod_poly_get_coeff_ui(frobenius.raw(), 0) != 0 ||
            nmod_poly_get_coeff_ui(frobenius.raw(), 1) != 1) {
            continue;
        }

        flint::NmodPolyFactor roots;
        nmod_poly_roots(roots.raw(), reduced.raw(), 0);
        if (roots.raw()->num != degree) {
            continue;
        }

        nmod_t mod;
        nmod_init(&mod, p);
        const std::size_t base = root_powers.size();
        root_powers.resize(base + n * n);
        bool linear = true;
        for (std::size_t i = 0; i < n && linear; ++i) {
            const nmod_poly_struct* factor =
                    roots.raw()->p + static_cast<slong>(i);
            if (nmod_poly_length(factor) != 2 ||
                nmod_poly_get_coeff_ui(factor, 1) != 1) {
                linear = false;
                break;
            }
            const ulong root =
                    nmod_neg(nmod_poly_get_coeff_ui(factor, 0), mod);
            ulong power = 1;
            for (std::size_t j = 0; j < n; ++j) {
                root_powers[base + i * n + j] = power;
                power = nmod_mul(power, root, mod);
            }
        }
        if (!linear) {
            root_powers.resize(base);
            continue;
        }
        moduli.push_back(mod);
        return true;
    }
    return false;
}

bool SplitPrimeResultant::add_closed_form_roots(nmod_t mod) noexcept {
    const ulong p = mod.n;
    const ulong index = split_closed_form_index;
    const bool real = split_closed_form == SplitClosedForm::real_cyclotomic;
    const auto n = static_cast<std::size_t>(degree);
    reduced_coefficients.resize(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        reduced_coefficients[i] = fmpz_fdiv_ui(
                field_polynomial.raw()->coeffs + static_cast<slong>(i), p);
    }

    // A first root.  When n | p - 1, zeta = a^((p - 1) / n) has exact order
    // n exactly when Phi_n(zeta) = 0, or, for psi_n, when
    // psi_n(zeta + zeta^-1) = 0 (the roots zeta^k + zeta^-k determine
    // zeta^(+-k)).  When n | p + 1 (psi_n only), alpha = (a + sqrt(a^2 - 4)) / 2
    // lies in F_(p^2) with alpha^(p+1) = 1 when a^2 - 4 is not a square, and
    // y = V_e(a) = beta + beta^-1 with beta = alpha^e, e = (p + 1) / n, is a
    // root of psi_n when beta has exact order n.  Every candidate is
    // checked exactly against f, which also rejects the bases a where
    // these conditions fail.
    const bool minus_one = (p - 1) % index != 0;
    const ulong exponent = minus_one ? (p + 1) / index : (p - 1) / index;
    ulong first = 0;
    bool found = false;
    // a = 2 gives alpha = 1, so the n | p + 1 search starts at a = 3.
    const ulong first_base = minus_one ? 3 : 2;
    for (ulong a = first_base; a < first_base + kClosedFormRootAttempts &&
                               !found;
         ++a) {
        ulong candidate = 0;
        if (minus_one) {
            candidate = lucas_v(a, exponent, mod);
        } else {
            candidate = nmod_pow_ui(a, exponent, mod);
            if (real) {
                candidate = nmod_add(candidate, nmod_inv(candidate, mod),
                                     mod);
            }
        }
        if (evaluate_reduced(reduced_coefficients, candidate, mod) == 0) {
            first = candidate;
            found = true;
        }
    }
    if (!found) {
        return false;
    }

    // The roots in increasing k coprime to n: zeta^k, 0 < k < n, for
    // Phi_n, and V_k(y) = beta^k + beta^-k, 0 < k < n / 2, for psi_n
    // (V_0 = 2, V_1 = y, V_(k+1) = y V_k - V_(k-1)).  Each is checked to be
    // a root of f, and the n roots to be distinct, so f = prod (x - r_i)
    // modulo p.
    const std::size_t base = root_powers.size();
    root_powers.resize(base + n * n);
    std::size_t count = 0;
    bool ok = true;
    ulong previous = nmod_set_ui(2, mod);
    ulong current = real ? first : 1;
    const ulong limit = real ? (index + 1) / 2 : index;
    for (ulong k = 1; k < limit && ok; ++k) {
        if (real) {
            if (k > 1) {
                const ulong next = nmod_sub(nmod_mul(first, current, mod),
                                            previous, mod);
                previous = current;
                current = next;
            }
        } else {
            current = nmod_mul(current, first, mod);
        }
        if (n_gcd(k, index) != 1) {
            continue;
        }
        if (count >= n ||
            evaluate_reduced(reduced_coefficients, current, mod) != 0) {
            ok = false;
            break;
        }
        ulong power = 1;
        for (std::size_t j = 0; j < n; ++j) {
            root_powers[base + count * n + j] = power;
            power = nmod_mul(power, current, mod);
        }
        ++count;
    }
    if (ok && count == n) {
        // root_powers[base + i n + 1] is r_i (n >= 2).
        std::vector<ulong> roots(n);
        for (std::size_t i = 0; i < n; ++i) {
            roots[i] = root_powers[base + i * n + 1];
        }
        std::sort(roots.begin(), roots.end());
        ok = std::adjacent_find(roots.begin(), roots.end()) == roots.end();
    }
    if (!ok || count != n) {
        root_powers.resize(base);
        return false;
    }
    return true;
}

void degree_two_resultant(bool& handled,
                          fmpz_t out,
                          const fmpz_poly_struct* f,
                          const fmpz_poly_struct* g) noexcept {
    handled = false;
    if (f == nullptr || g == nullptr || f->length != 3 || g->length > 2) {
        return;
    }
    handled = true;
    if (g->length == 0) {
        fmpz_zero(out);
        return;
    }
    const fmpz* g0 = g->coeffs;
    if (g->length == 1) {
        fmpz_mul(out, g0, g0);
        return;
    }
    const fmpz* g1 = g->coeffs + 1;
    const fmpz* f0 = f->coeffs;
    const fmpz* f1 = f->coeffs + 1;
    const fmpz* f2 = f->coeffs + 2;
    // g0 (f2 g0 - f1 g1) + (f0 g1) g1.
    flint::Fmpz t;
    fmpz_mul(t.raw(), f2, g0);
    fmpz_submul(t.raw(), f1, g1);
    fmpz_mul(out, t.raw(), g0);
    fmpz_mul(t.raw(), f0, g1);
    fmpz_addmul(out, t.raw(), g1);
}

namespace {

// a(r) mod p = sum_j r^j a_j for the powers r^j of one root and the
// residues a_j mod p, j < length.  Each product is below 2^56 (p < 2^28),
// so a chunk of kSplitPrimeDotChunk of them fits in one word.
ulong split_prime_root_evaluation(const ulong* powers,
                                  const ulong* residues,
                                  std::size_t length,
                                  nmod_t mod) noexcept {
    ulong evaluation = 0;
    for (std::size_t start = 0; start < length;
         start += static_cast<std::size_t>(kSplitPrimeDotChunk)) {
        std::size_t end =
                start + static_cast<std::size_t>(kSplitPrimeDotChunk);
        if (end > length) {
            end = length;
        }
        ulong accumulator = 0;
        for (std::size_t j = start; j < end; ++j) {
            accumulator += powers[j] * residues[j];
        }
        evaluation = nmod_add(evaluation,
                              n_mod2_preinv(accumulator, mod.n, mod.ninv),
                              mod);
    }
    return evaluation;
}

}  // namespace

void SplitPrimeResultant::bind_row_basis(
        const fmpz_mat_struct* basis) noexcept {
    row_basis = basis;
    const std::size_t rows =
            basis == nullptr ? 0
                             : static_cast<std::size_t>(fmpz_mat_nrows(basis));
    if (row_values.size() < rows) {
        row_values.resize(rows);
    }
    for (std::vector<ulong>& values : row_values) {
        values.clear();
    }
    row_value_fills = 0;
    clear_row_combination();
}

void SplitPrimeResultant::set_row_combination(slong left,
                                              slong right,
                                              slong sign) noexcept {
    row_candidate_set = true;
    row_candidate_left = left;
    row_candidate_right = right;
    row_candidate_sign = sign;
}

void SplitPrimeResultant::clear_row_combination() noexcept {
    row_candidate_set = false;
    row_candidate_left = -1;
    row_candidate_right = -1;
    row_candidate_sign = 0;
}

bool SplitPrimeResultant::row_combination_matches(
        const fmpz_poly_struct* coordinate_polynomial) noexcept {
    if (!row_candidate_set || row_basis == nullptr ||
        fmpz_mat_ncols(row_basis) != degree) {
        return false;
    }
    const slong rows = fmpz_mat_nrows(row_basis);
    const slong left = row_candidate_left;
    const slong right = row_candidate_right;
    if (left < 0 || left >= rows ||
        static_cast<std::size_t>(rows) > row_values.size() ||
        (right >= 0 && (right >= rows || (row_candidate_sign != 1 &&
                                          row_candidate_sign != -1)))) {
        return false;
    }
    // The polynomial is the described combination, coefficient by
    // coefficient, so the row values evaluate exactly this polynomial.
    for (slong j = 0; j < degree; ++j) {
        const fmpz* expected = fmpz_mat_entry(row_basis, left, j);
        if (right >= 0) {
            if (row_candidate_sign > 0) {
                fmpz_add(row_scratch.raw(), expected,
                         fmpz_mat_entry(row_basis, right, j));
            } else {
                fmpz_sub(row_scratch.raw(), expected,
                         fmpz_mat_entry(row_basis, right, j));
            }
            expected = row_scratch.raw();
        }
        if (j < coordinate_polynomial->length) {
            if (fmpz_equal(expected, coordinate_polynomial->coeffs + j) ==
                0) {
                return false;
            }
        } else if (fmpz_is_zero(expected) == 0) {
            return false;
        }
    }
    return true;
}

const ulong* SplitPrimeResultant::row_values_at(
        slong row,
        std::size_t prime_index) noexcept {
    const auto n = static_cast<std::size_t>(degree);
    std::vector<ulong>& values = row_values[static_cast<std::size_t>(row)];
    // Primes are filled in order; a call reaches prime k only after it
    // has used primes 0, ..., k - 1.
    while (values.size() < (prime_index + 1) * n) {
        const std::size_t prime = values.size() / n;
        const nmod_t mod = moduli[prime];
        for (std::size_t j = 0; j < n; ++j) {
            residues[j] = fmpz_fdiv_ui(
                    fmpz_mat_entry(row_basis, row, static_cast<slong>(j)),
                    mod.n);
        }
        const ulong* powers = root_powers.data() + prime * n * n;
        const std::size_t base = values.size();
        values.resize(base + n);
        for (std::size_t i = 0; i < n; ++i) {
            values[base + i] = split_prime_root_evaluation(
                    powers + i * n, residues.data(), n, mod);
        }
        ++row_value_fills;
    }
    return values.data() + prime_index * n;
}

void SplitPrimeResultant::resultant(
        bool& handled,
        fmpz_t out,
        const fmpz_poly_struct* coordinate_polynomial) noexcept {
    handled = false;
    last_used_row_values = false;
    // The row description applies to this call only.
    const bool row_described = row_candidate_set;
    const slong row_left = row_candidate_left;
    const slong row_right = row_candidate_right;
    const slong row_sign = row_candidate_sign;
    if (!enabled || coordinate_polynomial == nullptr) {
        clear_row_combination();
        return;
    }
    const bool use_rows =
            row_described && row_combination_matches(coordinate_polynomial);
    clear_row_combination();
    const slong length = coordinate_polynomial->length;
    if (length <= 0 || length > degree) {
        return;
    }

    // Hadamard bound on the Sylvester matrix:
    // |Res(f, a)| <= ||f||_2^(len(a) - 1) ||a||_2^deg(f) < 2^half_bits.
    _fmpz_vec_dot(coordinate_norm_squared.raw(),
                  coordinate_polynomial->coeffs,
                  coordinate_polynomial->coeffs, length);
    const flint_bitcnt_t doubled_bits =
            static_cast<flint_bitcnt_t>(length - 1) *
                    field_norm_squared_bits +
            static_cast<flint_bitcnt_t>(degree) *
                    fmpz_bits(coordinate_norm_squared.raw());
    flint_bitcnt_t half_bits = (doubled_bits + 1) / 2;
    // The T2 bound, when smaller (see SplitPrimeResultant).  It is set up
    // only once a split prime exists, so fields without one never pay
    // for it; the loop below needs that first prime anyway.
    if (moduli.empty() && !add_next_split_prime()) {
        return;
    }
    if (t2_bound_state == T2BoundState::not_computed) {
        t2_bound_state = certified_t2_ratio_bound(t2_ratio_bound,
                                                  field_polynomial.raw(),
                                                  degree)
                                 ? T2BoundState::available
                                 : T2BoundState::unavailable;
    }
    flint_bitcnt_t t2_half_bits = 0;
    if (t2_bound_state == T2BoundState::available &&
        t2_resultant_half_bits(t2_half_bits, coordinate_norm_squared.raw(),
                               t2_ratio_bound.raw(), degree) &&
        t2_half_bits < half_bits) {
        half_bits = t2_half_bits;
    }
    // The symmetric residue is the integer once the modulus exceeds twice
    // the bound: fmpz_bits(modulus) >= half_bits + 2 gives
    // modulus >= 2^(half_bits + 1).
    const flint_bitcnt_t modulus_bits = half_bits + 2;

    const auto n = static_cast<std::size_t>(degree);
    const auto used_length = static_cast<std::size_t>(length);
    fmpz_one(modulus.raw());
    fmpz_zero(out);
    std::size_t used = 0;
    while (fmpz_bits(modulus.raw()) < modulus_bits) {
        if (used == moduli.size() && !add_next_split_prime()) {
            return;
        }
        const nmod_t mod = moduli[used];

        // Res(f, a) = prod_i a(r_i) mod p for monic f = prod_i (x - r_i).
        ulong value = 1;
        if (use_rows) {
            row_values_at(row_left, used);
            if (row_right >= 0) {
                row_values_at(row_right, used);
            }
            const ulong* left_values =
                    row_values[static_cast<std::size_t>(row_left)].data() +
                    used * n;
            const ulong* right_values =
                    row_right < 0
                            ? nullptr
                            : row_values[static_cast<std::size_t>(row_right)]
                                              .data() +
                                      used * n;
            for (std::size_t i = 0; i < n && value != 0; ++i) {
                ulong evaluation = left_values[i];
                if (right_values != nullptr) {
                    evaluation = row_sign > 0
                                         ? nmod_add(evaluation,
                                                    right_values[i], mod)
                                         : nmod_sub(evaluation,
                                                    right_values[i], mod);
                }
                value = nmod_mul(value, evaluation, mod);
            }
        } else {
            const ulong* powers = root_powers.data() + used * n * n;
            for (std::size_t j = 0; j < used_length; ++j) {
                residues[j] = fmpz_fdiv_ui(
                        coordinate_polynomial->coeffs + static_cast<slong>(j),
                        mod.n);
            }
            for (std::size_t i = 0; i < n && value != 0; ++i) {
                value = nmod_mul(value,
                                 split_prime_root_evaluation(
                                         powers + i * n, residues.data(),
                                         used_length, mod),
                                 mod);
            }
        }

        if (used == 0) {
            fmpz_set_ui(out, value);
            if (value > mod.n / 2) {
                fmpz_sub_ui(out, out, mod.n);
            }
        } else {
            fmpz_CRT_ui(crt_scratch.raw(), out, modulus.raw(), value, mod.n,
                        1);
            fmpz_swap(out, crt_scratch.raw());
        }
        fmpz_mul_ui(modulus.raw(), modulus.raw(), mod.n);
        ++used;
    }
    last_prime_count = static_cast<slong>(used);
    last_used_row_values = use_rows;
    handled = true;
}

bool OrderCoordinateElementConversion::reset(
        const Order& next_order) noexcept {
    if (next_order.parent() == nullptr || next_order.degree() <= 0) {
        return false;
    }

    const slong degree = next_order.degree();
    flint::FmpqMat next_basis(degree, degree);
    flint::FmpqMat next_power_row(1, degree);
    if (!next_order.get_basis(flint::FmpqMatRef(next_basis))) {
        return false;
    }

    order = &next_order;
    basis = std::move(next_basis);
    power_row = std::move(next_power_row);
    return true;
}

bool OrderCoordinateElementConversion::set(
        Element& out,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den) noexcept {
    if (order == nullptr || order->parent() == nullptr ||
        out.parent() == nullptr ||
        !out.parent()->has_same_data(*order->parent()) ||
        fmpz_sgn(den.raw()) <= 0 ||
        flint::fmpz_mat_nrows(coordinates) != 1 ||
        flint::fmpz_mat_ncols(coordinates) != order->degree()) {
        return false;
    }

    ::fmpq_mat_mul_r_fmpz_mat(power_row.raw(), coordinates.raw(),
                              basis.raw());
    flint::fmpq_poly_zero(flint::FmpqPolyRef(polynomial));
    for (slong col = 0; col < order->degree(); ++col) {
        flint::fmpq_poly_set_coeff_fmpq(
                flint::FmpqPolyRef(polynomial), col,
                flint::fmpq_mat_entry(
                        flint::FmpqMatConstRef(power_row), 0, col));
    }
    if (fmpz_is_one(den.raw()) == 0) {
        flint::fmpq_poly_scalar_div_fmpz(
                flint::FmpqPolyRef(polynomial),
                flint::FmpqPolyConstRef(polynomial), den);
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

void coordinates_from_lattice_combination(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef coefficients,
        flint::FmpzMatConstRef basis) noexcept {
    flint::fmpz_mat_zero(flint::FmpzMatRef(out));

    for (slong i = 0; i < flint::fmpz_mat_nrows(basis); ++i) {
        for (slong j = 0; j < flint::fmpz_mat_ncols(basis); ++j) {
            flint::fmpz_addmul(flint::fmpz_mat_entry(out, 0, j),
                               flint::fmpz_mat_entry(coefficients, 0, i),
                               flint::fmpz_mat_entry(basis, i, j));
        }
    }
}

bool fmpq_poly_to_monic_integral_fmpz_poly(
        flint::FmpzPoly& out,
        const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree <= 0) {
        return false;
    }

    flint::Fmpq coefficient;
    fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial, degree);
    if (fmpq_is_one(coefficient.raw()) == 0) {
        return false;
    }

    fmpz_poly_zero(out.raw());
    flint::Fmpz integer_coefficient;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coefficient.raw())) == 0) {
            fmpz_poly_zero(out.raw());
            return false;
        }
        fmpz_set(integer_coefficient.raw(), fmpq_numref(coefficient.raw()));
        fmpz_poly_set_coeff_fmpz(out.raw(), i, integer_coefficient.raw());
    }
    return true;
}

void integral_coordinates_to_polynomial(
        flint::FmpzPoly& out,
        flint::FmpzMatConstRef coordinates) noexcept {
    fmpz_poly_zero(out.raw());
    const slong cols = flint::fmpz_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        fmpz_poly_set_coeff_fmpz(
                out.raw(), j,
                flint::fmpz_mat_entry(coordinates, 0, j).raw());
    }
}

bool build_norm_prefilter(NormPrefilter& out,
                          const FactorBase& base,
                          flint::FmpzConstRef factor_base_bound,
                          bool allow_large_prime) noexcept {
    out.valid = false;
    out.allow_large_prime = allow_large_prime;
    out.rejected_count = 0;
    out.has_equation_order_norm_polynomial = false;
    if (!base.is_defined() || base.length() <= 0) {
        return true;
    }

    const Order* order = base.parent();
    if (order == nullptr) {
        return false;
    }

    flint::fmpz_one(flint::FmpzRef(out.rational_prime_product));
    flint::fmpz_mul(flint::FmpzRef(out.large_prime_norm_bound),
                    factor_base_bound, factor_base_bound);

    flint::Fmpz rational_prime;
    for (slong block = 0; block < base.rational_prime_block_count();
         ++block) {
        slong block_length = 0;
        if (!base.rational_prime_block_data(
                    flint::FmpzRef(rational_prime), block_length, block)) {
            return false;
        }
        if (block_length <= 0) {
            continue;
        }
        flint::fmpz_mul(flint::FmpzRef(out.rational_prime_product),
                        flint::FmpzConstRef(out.rational_prime_product),
                        flint::FmpzConstRef(rational_prime));
    }

    const NumberField* field = order->parent();
    const nf_struct* raw_field = field == nullptr
                                         ? nullptr
                                         : field->raw_flint_field();
    out.has_equation_order_norm_polynomial =
            order->is_equation_order() && raw_field != nullptr &&
            fmpq_poly_to_monic_integral_fmpz_poly(out.field_polynomial,
                                                  raw_field->pol);
    out.split_prime_resultant.reset(
            out.has_equation_order_norm_polynomial
                    ? out.field_polynomial.raw()
                    : nullptr);

    out.valid = true;
    return true;
}

bool norm_prefilter_accepts_current_norm(
        bool& accepts,
        NormPrefilter& prefilter,
        const DiagnosticsContext* diagnostics) noexcept {
    accepts = true;
    if (!flint::fmpz_is_one(flint::fmpq_den_ref(
                flint::FmpqConstRef(prefilter.norm_scratch)))) {
        return true;
    }

    flint::fmpz_abs(flint::FmpzRef(prefilter.norm_remainder),
                    flint::fmpq_num_ref(
                            flint::FmpqConstRef(prefilter.norm_scratch)));
    if (flint::fmpz_cmp_ui(
                flint::FmpzConstRef(prefilter.norm_remainder), 1) <= 0) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.factor_base_gcd");
        flint::fmpz_gcd(flint::FmpzRef(prefilter.norm_gcd),
                        flint::FmpzConstRef(prefilter.rational_prime_product),
                        flint::FmpzConstRef(prefilter.norm_remainder));
        while (!flint::fmpz_is_one(flint::FmpzConstRef(prefilter.norm_gcd))) {
            flint::fmpz_divexact(
                    flint::FmpzRef(prefilter.norm_remainder),
                    flint::FmpzConstRef(prefilter.norm_remainder),
                    flint::FmpzConstRef(prefilter.norm_gcd));
            flint::fmpz_gcd(flint::FmpzRef(prefilter.norm_gcd),
                            flint::FmpzConstRef(prefilter.norm_gcd),
                            flint::FmpzConstRef(prefilter.norm_remainder));
        }
    }

    if (flint::fmpz_is_one(flint::FmpzConstRef(prefilter.norm_remainder))) {
        return true;
    }
    if (prefilter.allow_large_prime) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.large_prime_screen");
        if (flint::fmpz_cmp(flint::FmpzConstRef(prefilter.norm_remainder),
                            flint::FmpzConstRef(
                                    prefilter.large_prime_norm_bound)) <= 0 &&
            flint::fmpz_is_prime(
                    flint::FmpzConstRef(prefilter.norm_remainder))) {
            return true;
        }
    }

    accepts = false;
    ++prefilter.rejected_count;
    return true;
}

bool norm_prefilter_accepts(bool& accepts,
                            NormPrefilter* prefilter,
                            const Element& alpha,
                            const DiagnosticsContext* diagnostics) noexcept {
    accepts = true;
    if (prefilter == nullptr || !prefilter->valid) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_prefilter.norm");
        if (!alpha.norm(flint::FmpqRef(prefilter->norm_scratch))) {
            return false;
        }
    }
    return norm_prefilter_accepts_current_norm(accepts, *prefilter,
                                               diagnostics);
}

bool norm_prefilter_accepts_integral_equation_coordinates(
        bool& handled,
        bool& accepts,
        NormPrefilter* prefilter,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den,
        const DiagnosticsContext* diagnostics) noexcept {
    handled = false;
    accepts = true;
    if (prefilter == nullptr || !prefilter->valid ||
        !prefilter->has_equation_order_norm_polynomial ||
        fmpz_is_one(den.raw()) == 0 ||
        flint::fmpz_mat_nrows(coordinates) != 1) {
        return true;
    }
    if (flint::fmpz_mat_ncols(coordinates) !=
        fmpz_poly_degree(prefilter->field_polynomial.raw())) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.coordinate_norm");
        integral_coordinates_to_polynomial(prefilter->coordinate_polynomial,
                                           coordinates);
        // Quadratic fields: the closed-form Sylvester determinant.
        bool resultant_handled = false;
        degree_two_resultant(resultant_handled,
                             prefilter->norm_resultant_scratch.raw(),
                             prefilter->field_polynomial.raw(),
                             prefilter->coordinate_polynomial.raw());
        if (!resultant_handled) {
            prefilter->split_prime_resultant.resultant(
                    resultant_handled,
                    prefilter->norm_resultant_scratch.raw(),
                    prefilter->coordinate_polynomial.raw());
        }
        if (!resultant_handled) {
            fmpz_poly_resultant(prefilter->norm_resultant_scratch.raw(),
                                prefilter->field_polynomial.raw(),
                                prefilter->coordinate_polynomial.raw());
        }
        flint::fmpq_set_fmpz(
                flint::FmpqRef(prefilter->norm_scratch),
                flint::FmpzConstRef(prefilter->norm_resultant_scratch));
    }
    handled = true;
    return norm_prefilter_accepts_current_norm(accepts, *prefilter,
                                               diagnostics);
}

bool try_generator_relation(ClassGroupContext& context,
                            const Element& alpha,
                            const fmpq* known_norm,
                            const fmpz_mat_struct* known_integral_coordinates,
                            const flint::FmpzPoly*
                                    known_integral_coordinate_polynomial,
                            ClassGroupRelationSource source,
                            slong& accepted_relations,
                            bool& partial_throttle_exit,
                            bool& relation_appended) noexcept {
    const slong relation_count_before = context.relation_count();
    bool candidate_partial_throttle = false;
    relation_appended = false;
    bool ok = false;
    if (known_norm != nullptr && known_integral_coordinates != nullptr) {
        // reference Rel_add.jl and reference factorgen keep the exact norm/factor data
        // produced during admission.  Keep this prefilter data on the same
        // private path instead of reconstructing it for relation materialization.
        ok = detail::ClassGroupRelationSearchAccess::
                try_append_integral_generator_relation(
                        context, candidate_partial_throttle, alpha,
                        flint::FmpzMatConstRef(known_integral_coordinates),
                        flint::FmpqConstRef(known_norm),
                        known_integral_coordinate_polynomial, source);
    } else if (known_norm != nullptr) {
        ok = detail::ClassGroupRelationSearchAccess::
                try_append_generator_relation_with_norm(
                        context, candidate_partial_throttle, alpha,
                        flint::FmpqConstRef(known_norm), source);
    } else {
        ok = context.try_append_generator_relation(
                candidate_partial_throttle, alpha, source);
    }
    if (!ok) {
        return false;
    }
    if (context.relation_count() > relation_count_before) {
        ++accepted_relations;
        relation_appended = true;
    }
    if (candidate_partial_throttle &&
        source == ClassGroupRelationSource::Search) {
        partial_throttle_exit = true;
    }
    return true;
}

bool element_from_order_coordinates_den(
        Element& out,
        const Order& order,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den) noexcept {
    if (out.parent() == nullptr || order.parent() == nullptr ||
        !out.parent()->has_same_data(*order.parent()) ||
        fmpz_sgn(den.raw()) <= 0) {
        return false;
    }

    OrderElement order_element(order);
    Element candidate(*order.parent());
    if (!order_element.is_defined() || !candidate.is_defined() ||
        !order_element.set_coordinates(coordinates) ||
        !order_element.get_element(candidate)) {
        return false;
    }

    if (fmpz_is_one(den.raw()) == 0) {
        flint::FmpqPoly polynomial;
        if (!candidate.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
            return false;
        }
        flint::fmpq_poly_scalar_div_fmpz(
                polynomial, polynomial, den);
        if (!candidate.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial))) {
            return false;
        }
    }

    out.swap(candidate);
    return true;
}

bool element_is_scalar_rational(bool& scalar, const Element& element) noexcept {
    scalar = false;
    flint::FmpqPoly polynomial;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }

    scalar = flint::fmpq_poly_degree(
                     flint::FmpqPolyConstRef(polynomial)) <= 0;
    return true;
}

bool reduced_basis_first_row_is_scalar_rational(
        bool& scalar,
        const Order& order,
        flint::FmpzMatConstRef basis) noexcept {
    scalar = false;
    const NumberField* field = order.parent();
    const slong degree = order.degree();
    if (field == nullptr || degree <= 0 ||
        flint::fmpz_mat_nrows(basis) != degree ||
        flint::fmpz_mat_ncols(basis) != degree) {
        return false;
    }

    flint::FmpzMat row(1, degree);
    for (slong j = 0; j < degree; ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(row, 0, j),
                        flint::fmpz_mat_entry(basis, 0, j));
    }

    flint::Fmpz denominator;
    flint::fmpz_one(flint::FmpzRef(denominator));
    Element element(*field);
    return element.is_defined() &&
           element_from_order_coordinates_den(
                   element, order, flint::FmpzMatConstRef(row),
                   flint::FmpzConstRef(denominator)) &&
           element_is_scalar_rational(scalar, element);
}

bool try_coordinate_candidate_den(ClassGroupContext& context,
                                  const Order& order,
                                  flint::FmpzMat& coordinates,
                                  flint::FmpzConstRef den,
                                  NormPrefilter* norm_prefilter,
                                  ClassGroupRelationSource source,
                                  slong target_relation_kernel_units,
                                  slong max_candidates,
                                  slong max_relations,
                                  slong& candidates_tried,
                                  slong& accepted_relations,
                                  CertificationMode requested_certification,
                                  bool& goal_reached,
                                  bool& partial_throttle_exit,
                                  detail::RelationAdmissionCache* admission_cache,
                                  bool random_relation,
                                  slong* factor_attempts,
                                  slong max_factor_attempts,
                                  bool* factor_attempt_limit_reached,
                                  Relation* admission_scratch,
                                  OrderCoordinateElementConversion*
                                          coordinate_conversion,
                                  bool admission_scratch_base_verified)
        noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate");
    if (candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    bool scalar_check_done = false;
    if (admission_cache != nullptr && order.is_equation_order() &&
        flint::fmpz_mat_nrows(flint::FmpzMatConstRef(coordinates)) == 1 &&
        flint::fmpz_mat_ncols(flint::FmpzMatConstRef(coordinates)) ==
                order.degree()) {
        bool scalar = true;
        for (slong col = 1; col < order.degree(); ++col) {
            if (!flint::fmpz_is_zero(
                        flint::fmpz_mat_entry(coordinates, 0, col))) {
                scalar = false;
                break;
            }
        }
        scalar_check_done = true;
        if (scalar) {
            return true;
        }
    }

    Element alpha(*order.parent());
    bool alpha_ready = false;
    auto ensure_alpha = [&]() noexcept -> bool {
        if (alpha_ready) {
            return true;
        }
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_element");
        const bool converted =
                alpha.is_defined() &&
                (coordinate_conversion == nullptr
                         ? element_from_order_coordinates_den(
                                   alpha, order,
                                   flint::FmpzMatConstRef(coordinates), den)
                         : coordinate_conversion->set(
                                   alpha,
                                   flint::FmpzMatConstRef(coordinates), den));
        if (!converted) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate conversion failed");
            return false;
        }
        alpha_ready = true;
        return true;
    };

    if (admission_cache != nullptr && !scalar_check_done) {
        bool scalar = false;
        if (!ensure_alpha() || !element_is_scalar_rational(scalar, alpha)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate scalar check failed");
            return false;
        }
        if (scalar) {
            return true;
        }
    }

    if (admission_cache != nullptr && factor_attempts != nullptr &&
        max_factor_attempts > 0) {
        // The source finite quadratic-form ideal search increments try_factor before
        // factorgen, so candidates rejected by the norm/factorability gate
        // still consume the per-ideal factor-attempt budget.
        if (*factor_attempts >= max_factor_attempts) {
            if (factor_attempt_limit_reached != nullptr) {
                *factor_attempt_limit_reached = true;
            }
            return true;
        }
        ++*factor_attempts;
    }

    bool prefilter_accepts = true;
    const fmpq* known_norm = nullptr;
    const fmpz_mat_struct* known_integral_coordinates = nullptr;
    const flint::FmpzPoly* known_integral_coordinate_polynomial = nullptr;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_prefilter");
        bool norm_prefilter_handled = false;
        if (!norm_prefilter_accepts_integral_equation_coordinates(
                    norm_prefilter_handled, prefilter_accepts,
                    norm_prefilter,
                    flint::FmpzMatConstRef(coordinates), den,
                    diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate coordinate norm prefilter failed");
            return false;
        }
        if (!norm_prefilter_handled &&
            (!ensure_alpha() ||
             !norm_prefilter_accepts(prefilter_accepts, norm_prefilter,
                                     alpha, diagnostics))) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate norm prefilter failed");
            return false;
        }
        if (norm_prefilter != nullptr && norm_prefilter->valid) {
            known_norm = norm_prefilter->norm_scratch.raw();
            if (flint::fmpz_is_one(den)) {
                known_integral_coordinates =
                        flint::FmpzMatConstRef(coordinates).raw();
            }
            if (norm_prefilter_handled) {
                known_integral_coordinate_polynomial =
                        &norm_prefilter->coordinate_polynomial;
            }
        }
    }
    if (!prefilter_accepts) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                            "class_group.candidate_prefilter_reject");
        return true;
    }
    const bool can_defer_integral_alpha =
            admission_cache != nullptr && known_norm != nullptr &&
            flint::fmpz_is_one(den) && order.is_equation_order() &&
            flint::fmpz_mat_nrows(flint::FmpzMatConstRef(coordinates)) == 1 &&
            flint::fmpz_mat_ncols(flint::FmpzMatConstRef(coordinates)) ==
                    order.degree();
    if (!can_defer_integral_alpha && !ensure_alpha()) {
        return false;
    }

    ++candidates_tried;
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.candidate_relation");
    bool candidate_partial_throttle = false;
    bool relation_appended = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_try_relation");
        if (admission_cache != nullptr) {
            const FactorBase* base = context.factor_base();
            if (base == nullptr) {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "relation admission missing factor base");
                return false;
            }
            Relation local_relation;
            Relation* relation = admission_scratch;
            if (relation == nullptr) {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.admission_scratch_define");
                if (!local_relation.define(*base)) {
                    SILEX_LOG(
                            diagnostics, DiagnosticsModule::class_group,
                            LogLevel::detail,
                            "relation admission scratch allocation failed");
                    return false;
                }
                relation = &local_relation;
            } else if (!relation->is_defined() ||
                       relation->factor_base() == nullptr ||
                       (admission_scratch_base_verified
                                ? relation->length() != base->length()
                                : !relation->factor_base()->equal(*base))) {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "relation admission scratch mismatch");
                return false;
            }
            bool relation_set = false;
            bool deferred_integral_relation_row = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate_set_relation");
                if (can_defer_integral_alpha) {
                    bool handled = false;
                    bool smooth = false;
                    if (!detail::RelationAccess::
                                factor_relation_row_from_integral_coordinates_and_norm(
                                        *relation, handled, smooth,
                                        flint::FmpzMatConstRef(coordinates),
                                        flint::FmpqConstRef(known_norm),
                                        known_integral_coordinate_polynomial,
                                        diagnostics)) {
                        return false;
                    }
                    if (handled) {
                        relation_set = smooth;
                        deferred_integral_relation_row = smooth;
                    } else {
                        if (!ensure_alpha()) {
                            return false;
                        }
                        relation_set =
                                detail::RelationAccess::
                                        set_relation_from_integral_coordinates_and_norm(
                                        *relation, alpha,
                                        flint::FmpzMatConstRef(coordinates),
                                        flint::FmpqConstRef(known_norm),
                                        diagnostics);
                    }
                } else if (known_norm != nullptr &&
                           flint::fmpz_is_one(den)) {
                    relation_set =
                            detail::RelationAccess::
                                    set_relation_from_integral_coordinates_and_norm(
                                    *relation, alpha,
                                    flint::FmpzMatConstRef(coordinates),
                                    flint::FmpqConstRef(known_norm),
                                    diagnostics);
                } else {
                    relation_set = known_norm == nullptr
                            ? relation->set_generator(alpha, diagnostics)
                            : detail::RelationAccess::set_generator_with_norm(
                                      *relation, alpha,
                                      flint::FmpqConstRef(known_norm),
                                      diagnostics);
                }
            }
            if (!relation_set) {
                return true;
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate_admit_relation");
                bool admitted = false;
                if (deferred_integral_relation_row) {
                    admitted = detail::try_admit_deferred_integral_relation(
                            context, relation_appended, *admission_cache,
                            *relation,
                            flint::FmpzMatConstRef(coordinates),
                            random_relation,
                            admission_scratch_base_verified);
                } else {
                    admitted = detail::try_admit_relation(
                            context, relation_appended, *admission_cache,
                            *relation, random_relation,
                            admission_scratch_base_verified);
                }
                if (!admitted) {
                    SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                              LogLevel::detail,
                              "relation admission failed");
                    return false;
                }
            }
            if (relation_appended) {
                ++accepted_relations;
            }
        } else {
            if (!try_generator_relation(
                        context, alpha, known_norm,
                        known_integral_coordinates,
                        known_integral_coordinate_polynomial, source,
                        accepted_relations, candidate_partial_throttle,
                        relation_appended)) {
                return false;
            }
        }
    }

    if (relation_appended) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_publish");
        if (admission_cache != nullptr) {
            goal_reached = context.relation_count() >= max_relations;
        } else {
            goal_reached = publish_and_check_compute_goal(
                    context, target_relation_kernel_units,
                    requested_certification);
        }
    }
    if (admission_cache == nullptr && candidate_partial_throttle &&
        !goal_reached) {
        partial_throttle_exit = true;
    }
    return true;
}

bool try_coordinate_candidate(ClassGroupContext& context,
                              const Order& order,
                              flint::FmpzMat& coordinates,
                              NormPrefilter* norm_prefilter,
                              ClassGroupRelationSource source,
                              slong target_relation_kernel_units,
                              slong max_candidates,
                              slong max_relations,
                              slong& candidates_tried,
                              slong& accepted_relations,
                              CertificationMode requested_certification,
                              bool& goal_reached,
                              bool& partial_throttle_exit,
                              OrderCoordinateElementConversion*
                                      coordinate_conversion) noexcept {
    flint::Fmpz den;
    fmpz_one(den.raw());
    return try_coordinate_candidate_den(
            context, order, coordinates, flint::FmpzConstRef(den),
            norm_prefilter, source, target_relation_kernel_units,
            max_candidates, max_relations, candidates_tried,
            accepted_relations, requested_certification, goal_reached,
            partial_throttle_exit, nullptr, false, nullptr, 0, nullptr,
            nullptr, coordinate_conversion);
}

}  // namespace detail::relation_search
}  // namespace silex
