#include <silex/class_group.hpp>

#include <silex/order_element.hpp>
#include <silex/relation.hpp>

#include "relation_admission_cache_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"
#include "../relation/relation_internal.hpp"

#include <silex/flint/acb_vec.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/nmod_poly_factor.hpp>

#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/mag.h>
#include <flint/ulong_extras.h>

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
    moduli.clear();
    root_powers.clear();
    residues.clear();
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

        // f splits into distinct linear factors modulo p exactly when
        // f divides x^p - x in F_p[x]; deg(f) >= 2 here, so x mod f = x.
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

void SplitPrimeResultant::resultant(
        bool& handled,
        fmpz_t out,
        const fmpz_poly_struct* coordinate_polynomial) noexcept {
    handled = false;
    if (!enabled || coordinate_polynomial == nullptr) {
        return;
    }
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
    // The T2 bound, when smaller (see SplitPrimeResultant).
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
        const ulong* powers = root_powers.data() + used * n * n;
        for (std::size_t j = 0; j < used_length; ++j) {
            residues[j] = fmpz_fdiv_ui(
                    coordinate_polynomial->coeffs + static_cast<slong>(j),
                    mod.n);
        }

        // Res(f, a) = prod_i a(r_i) mod p for monic f = prod_i (x - r_i).
        ulong value = 1;
        for (std::size_t i = 0; i < n && value != 0; ++i) {
            const ulong* row = powers + i * n;
            ulong evaluation = 0;
            for (std::size_t start = 0; start < used_length;
                 start += static_cast<std::size_t>(kSplitPrimeDotChunk)) {
                std::size_t end =
                        start + static_cast<std::size_t>(kSplitPrimeDotChunk);
                if (end > used_length) {
                    end = used_length;
                }
                ulong accumulator = 0;
                for (std::size_t j = start; j < end; ++j) {
                    accumulator += row[j] * residues[j];
                }
                evaluation = nmod_add(
                        evaluation,
                        n_mod2_preinv(accumulator, mod.n, mod.ninv), mod);
            }
            value = nmod_mul(value, evaluation, mod);
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
        bool split_prime_handled = false;
        prefilter->split_prime_resultant.resultant(
                split_prime_handled,
                prefilter->norm_resultant_scratch.raw(),
                prefilter->coordinate_polynomial.raw());
        if (!split_prime_handled) {
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
