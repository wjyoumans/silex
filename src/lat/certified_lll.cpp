#include "certified_lll_internal.hpp"

#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpz.hpp>

#include <flint/arb.h>
#include <flint/arb_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

#include <cmath>

namespace silex::lat::detail {

namespace {

// The precision fmpz_lll_wrapper_with_removal_knapsack passes to
// fmpz_lll_is_reduced_with_removal (wrapper_with_removal_knapsack.c).
constexpr flint_bitcnt_t kFlintWrapperCheckPrecision = 120;

// Writes the finite positive double `value` exactly as
// `mantissa * 2^-shift` with `shift >= 0`.
bool double_as_dyadic(flint::Fmpz& mantissa,
                      ulong& shift,
                      double value) noexcept {
    if (!std::isfinite(value) || !(value > 0.0)) {
        return false;
    }
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    // `fraction * 2^53` is an integer below 2^53, so the conversion is exact.
    ::fmpz_set_d(mantissa.raw(), std::ldexp(fraction, 53));
    const slong scaled_shift = 53 - static_cast<slong>(exponent);
    if (scaled_shift < 0) {
        ::fmpz_mul_2exp(mantissa.raw(), mantissa.raw(),
                        static_cast<ulong>(-scaled_shift));
        shift = 0;
    } else {
        shift = static_cast<ulong>(scaled_shift);
    }
    return true;
}

}  // namespace

bool gram_interval_screen_proves_lll_reduced(flint::FmpzMatConstRef gram_matrix,
                                             double delta,
                                             double eta) noexcept {
    // Why true implies integer_gram_is_lll_reduced(gram, delta, eta):
    //
    // For a Gram matrix G of linearly independent b_0, ..., b_{n-1}, the
    // Gram-Schmidt data follow from G alone (the Cholesky factorization of
    // the Gram matrix: Nguyen and Stehle, "An LLL algorithm with quadratic
    // complexity", SIAM J. Comput. 39 (2009), Section 3.2, Figure 4, the
    // Cholesky factorization algorithm, whose s_i^(i) is the sum below
    // with j = i): with r_ij = <b_i, b*_j>,
    //   r_ij = G_ij - sum_{k<j} mu_jk r_ik   (j <= i),
    //   mu_ij = r_ij / r_jj                  (j < i),
    //   |b*_i|^2 = r_ii.
    // In the notation of integer_gram_is_lll_reduced, r_ii = d_{i+1} / d_i
    // and mu_ij = lambda_ij / d_{j+1}; this holds row by row as long as the
    // earlier leading minors are positive, which the loop below has proved
    // by then.  Its three conditions are therefore, for every row i:
    //   (a) r_ii > 0                     (d_{i+1} > 0),
    //   (b) |mu_ij| <= eta for j < i      (|lambda_ij| <= eta d_{j+1}),
    //   (c) r_ii - (delta - mu_{i,i-1}^2) r_{i-1,i-1} >= 0 for i >= 1
    //       (the Lovasz test, divided by d_i d_{i-1} > 0).
    //
    // The loop evaluates the same recurrence in FLINT's ball arithmetic.
    // Every arb operation returns a ball containing the exact result for
    // every choice of exact values in its input balls (arb_set_round_fmpz,
    // arb_dot, arb_div, arb_sqr, arb_sub, arb_mul, arb_add; arb_set_d is
    // exact for a finite double, so delta and eta are the exact rationals
    // of the doubles, as in the exact test).  By induction over the
    // recurrence each ball contains its exact quantity, because every
    // divisor r_jj was proved positive (a) before any later row divides
    // by it.  arb_is_positive, arb_le and arb_is_nonnegative hold only
    // when the predicate holds for every point of the balls, so for the
    // exact quantities.  If a leading minor is not positive, the first such
    // row has exact r_ii <= 0 in its ball and (a) fails.  Hence a true
    // return proves (a), (b) and (c) for every row, which is exactly what
    // integer_gram_is_lll_reduced checks; any undecided ball returns false
    // and leaves the decision to the exact recurrence.  Like the exact
    // test, the loop reads only the lower triangle and the diagonal.
    const slong rows = flint::fmpz_mat_nrows(gram_matrix);
    if (rows < 2 || flint::fmpz_mat_ncols(gram_matrix) != rows ||
        !std::isfinite(delta) || !(delta > 0.0) || !std::isfinite(eta) ||
        !(eta > 0.0)) {
        return false;
    }
    constexpr slong prec = gram_interval_screen_precision;
    const fmpz_mat_struct* gram = gram_matrix.raw();
    // Strict lower triangle: mu_ij; diagonal: r_ii = |b*_i|^2.
    flint::ArbMat mu(rows, rows);
    // r_ij of the current row i.
    flint::ArbVec r(rows);
    flint::Arb entry;
    flint::Arb bound;
    flint::Arb delta_ball;
    flint::Arb eta_ball;
    ::arb_set_d(delta_ball.raw(), delta);
    ::arb_set_d(eta_ball.raw(), eta);
    arb_mat_struct* m = mu.raw();
    arb_ptr row_r = r.data();
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j <= i; ++j) {
            ::arb_set_round_fmpz(entry.raw(), ::fmpz_mat_entry(gram, i, j),
                                 prec);
            // r_ij = G_ij - sum_{k<j} mu_jk r_ik.
            ::arb_dot(row_r + j, entry.raw(), 1, arb_mat_entry(m, j, 0), 1,
                      row_r, 1, j, prec);
            if (j < i) {
                arb_ptr mu_ij = arb_mat_entry(m, i, j);
                ::arb_div(mu_ij, row_r + j, arb_mat_entry(m, j, j), prec);
                ::arb_abs(bound.raw(), mu_ij);
                if (!::arb_le(bound.raw(), eta_ball.raw())) {  // (b)
                    return false;
                }
            }
        }
        arb_ptr r_ii = arb_mat_entry(m, i, i);
        ::arb_set(r_ii, row_r + i);
        if (!::arb_is_positive(r_ii)) {  // (a)
            return false;
        }
        if (i >= 1) {
            // (mu_{i,i-1}^2 - delta) r_{i-1,i-1} + r_ii >= 0.
            ::arb_sqr(bound.raw(), arb_mat_entry(m, i, i - 1), prec);
            ::arb_sub(bound.raw(), bound.raw(), delta_ball.raw(), prec);
            ::arb_mul(bound.raw(), bound.raw(),
                      arb_mat_entry(m, i - 1, i - 1), prec);
            ::arb_add(bound.raw(), bound.raw(), r_ii, prec);
            if (!::arb_is_nonnegative(bound.raw())) {  // (c)
                return false;
            }
        }
    }
    return true;
}

bool integer_gram_is_lll_reduced(flint::FmpzMatConstRef gram_matrix,
                                 double delta,
                                 double eta) noexcept {
    // Exact form of FLINT's `fmpz_mat_is_reduced_gram(A, delta, eta)` on a
    // positive definite Gram matrix, and so of
    // `fmpz_mat_is_reduced(B, delta, eta)` (`gr_mat_is_row_lll_reduced_naive`)
    // on a basis B with Gram matrix A = B B^T: every Gram-Schmidt norm is
    // nonzero, |mu_ij| <= eta for j < i, and
    // (delta - mu_{i,i-1}^2) |b*_{i-1}|^2 <= |b*_i|^2, with `delta` and
    // `eta` read as the exact rationals of the doubles.  It uses the
    // integral Gram-Schmidt recurrence of Cohen, A Course in Computational
    // Algebraic Number Theory, Algorithm 2.6.7, step 2: d_0 = 1,
    // d_i = det of the leading i x i Gram minor, lambda_ij = d_j mu_ij, so
    // |b*_i|^2 = d_i / d_{i-1}.  Then |mu_ij| <= eta is
    // |lambda_ij| <= eta d_j, and the Lovasz test, multiplied by
    // d_{i-1} d_{i-2} > 0, is delta d_{i-1}^2 <= d_i d_{i-2} + lambda^2.
    // Only the lower triangle and the diagonal are read, as in
    // fmpz_mat_is_reduced_gram.  Each d_i is checked positive before it is
    // divided by, so no input divides by zero.  Matrices with fewer than
    // two rows, or not square, return false (callers fall back).
    const slong rows = flint::fmpz_mat_nrows(gram_matrix);
    if (rows < 2 || flint::fmpz_mat_ncols(gram_matrix) != rows) {
        return false;
    }
    flint::Fmpz eta_mantissa;
    flint::Fmpz delta_mantissa;
    ulong eta_shift = 0;
    ulong delta_shift = 0;
    if (!double_as_dyadic(eta_mantissa, eta_shift, eta) ||
        !double_as_dyadic(delta_mantissa, delta_shift, delta)) {
        return false;
    }
    // A sufficient test in ball arithmetic, much cheaper than the exact
    // recurrence below; it returns true only when it proves the conditions
    // this test checks, so the result is unchanged either way.
    if (gram_interval_screen_proves_lll_reduced(gram_matrix, delta, eta)) {
        return true;
    }

    const fmpz_mat_struct* gram = gram_matrix.raw();
    flint::FmpzMat lambda(rows, rows);
    // determinants(0, i) holds d_i; d_0 = 1.
    flint::FmpzMat determinants(1, rows + 1);
    ::fmpz_one(::fmpz_mat_entry(determinants.raw(), 0, 0));

    flint::Fmpz u;
    flint::Fmpz left;
    flint::Fmpz right;
    for (slong k = 0; k < rows; ++k) {
        for (slong j = 0; j <= k; ++j) {
            ::fmpz_set(u.raw(), ::fmpz_mat_entry(gram, k, j));
            for (slong i = 0; i < j; ++i) {
                ::fmpz_mul(u.raw(), u.raw(),
                           ::fmpz_mat_entry(determinants.raw(), 0, i + 1));
                ::fmpz_submul(u.raw(), ::fmpz_mat_entry(lambda.raw(), k, i),
                              ::fmpz_mat_entry(lambda.raw(), j, i));
                ::fmpz_divexact(u.raw(), u.raw(),
                                ::fmpz_mat_entry(determinants.raw(), 0, i));
            }
            ::fmpz_set(j < k ? ::fmpz_mat_entry(lambda.raw(), k, j)
                             : ::fmpz_mat_entry(determinants.raw(), 0, k + 1),
                       u.raw());
        }
        const fmpz* d_k = ::fmpz_mat_entry(determinants.raw(), 0, k + 1);
        if (::fmpz_sgn(d_k) <= 0) {
            return false;
        }
        for (slong j = 0; j < k; ++j) {
            ::fmpz_abs(left.raw(), ::fmpz_mat_entry(lambda.raw(), k, j));
            ::fmpz_mul_2exp(left.raw(), left.raw(), eta_shift);
            ::fmpz_mul(right.raw(), eta_mantissa.raw(),
                       ::fmpz_mat_entry(determinants.raw(), 0, j + 1));
            if (::fmpz_cmp(left.raw(), right.raw()) > 0) {
                return false;
            }
        }
        if (k >= 1) {
            const fmpz* d_previous =
                    ::fmpz_mat_entry(determinants.raw(), 0, k);
            const fmpz* lambda_previous =
                    ::fmpz_mat_entry(lambda.raw(), k, k - 1);
            ::fmpz_mul(left.raw(), d_previous, d_previous);
            ::fmpz_mul(left.raw(), left.raw(), delta_mantissa.raw());
            ::fmpz_mul(right.raw(), d_k,
                       ::fmpz_mat_entry(determinants.raw(), 0, k - 1));
            ::fmpz_addmul(right.raw(), lambda_previous, lambda_previous);
            ::fmpz_mul_2exp(right.raw(), right.raw(), delta_shift);
            if (::fmpz_cmp(left.raw(), right.raw()) > 0) {
                return false;
            }
        }
    }
    return true;
}

bool integer_basis_is_lll_reduced(flint::FmpzMatConstRef basis,
                                  double delta,
                                  double eta) noexcept {
    // The Gram matrix B B^T of a basis is positive semidefinite, and
    // positive definite exactly when integer_gram_is_lll_reduced can
    // return true, so this decides fmpz_mat_is_reduced(B, delta, eta).
    const slong rows = flint::fmpz_mat_nrows(basis);
    const slong cols = flint::fmpz_mat_ncols(basis);
    if (rows < 2) {
        return false;
    }
    const fmpz_mat_struct* b = basis.raw();
    flint::FmpzMat gram(rows, rows);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j <= i; ++j) {
            ::_fmpz_vec_dot(::fmpz_mat_entry(gram.raw(), i, j),
                            ::fmpz_mat_entry(b, i, 0),
                            ::fmpz_mat_entry(b, j, 0), cols);
        }
    }
    return integer_gram_is_lll_reduced(flint::FmpzMatConstRef(gram), delta,
                                       eta);
}

bool lll_certificate_holds(flint::FmpzMatConstRef basis,
                           flint::FmpzLllConstRef config) noexcept {
    const fmpz_lll_struct* fl = config.raw();
    const slong bits = ::fmpz_mat_max_bits(basis.raw());
    if (fl->rt == GRAM) {
        const slong rows = flint::fmpz_mat_nrows(basis);
        if (rows * (bits < 0 ? -bits : bits) <=
            certified_lll_gram_integer_test_max_size) {
            return integer_gram_is_lll_reduced(basis, fl->delta, fl->eta);
        }
        // FLINT's own acceptance test of the phase output
        // (wrapper_with_removal_knapsack.c, through
        // fmpz_lll_is_reduced_with_removal with gs_B NULL), so the helper
        // accepts exactly when fmpz_lll does.
        return rows >= 2 &&
               ::fmpz_lll_is_reduced(basis.raw(), fl,
                                     kFlintWrapperCheckPrecision) != 0;
    }
    if ((bits < 0 ? -bits : bits) <= certified_lll_integer_test_max_bits) {
        return integer_basis_is_lll_reduced(basis, fl->delta, fl->eta);
    }
    return flint::fmpz_mat_nrows(basis) >= 2 &&
           ::fmpz_mat_is_reduced(basis.raw(), fl->delta, fl->eta) != 0;
}

namespace {

// FLINT's ULLL truncation threshold (lll_with_removal_ulll.c): fmpz_lll
// reduces a Z_BASIS matrix whose entries have at most this many bits with
// one call of its default wrapper (wrapper_with_removal_knapsack.c).  A
// GRAM matrix is reduced with one wrapper call at every size.
constexpr slong kFlintUlllTruncationBits = 250;

bool first_phase_eligible(const fmpz_mat_struct* basis,
                          const fmpz_lll_struct* fl) noexcept {
    if (fl->gt != APPROX || ::fmpz_mat_nrows(basis) < 2) {
        return false;
    }
    if (fl->rt == GRAM) {
        return ::fmpz_mat_ncols(basis) == ::fmpz_mat_nrows(basis);
    }
    if (fl->rt != Z_BASIS) {
        return false;
    }
    const slong bits = ::fmpz_mat_max_bits(basis);
    return bits >= -kFlintUlllTruncationBits &&
           bits <= kFlintUlllTruncationBits;
}

}  // namespace

CertifiedLllRoute lll_first_phase_certified(fmpz_mat_struct* basis,
                                            fmpz_mat_struct* transform,
                                            flint::FmpzLllConstRef config)
        noexcept {
    if (!first_phase_eligible(basis, config.raw())) {
        return CertifiedLllRoute::ineligible;
    }
    // The first phase of fmpz_lll's default wrapper.  Its return value is a
    // status, not a certificate (nor a rank, with no removal bound).
    int status = ::fmpz_lll_d_with_removal_knapsack(basis, transform, nullptr,
                                                    config.raw());
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::CertifiedLllTestSeams;
    const CertifiedLllTestFailure failure = CertifiedLllTestSeams::failure();
    if (failure != CertifiedLllTestFailure::none) {
        // Test-only; see CertifiedLllTestSeams.
        ::fmpz_mat_zero(basis);
        if (failure == CertifiedLllTestFailure::first_phase) {
            status = -1;
        }
        CertifiedLllTestSeams::note_injected();
    }
    if (status != -1 &&
        lll_certificate_holds(flint::FmpzMatConstRef(basis), config)) {
        return CertifiedLllRoute::certified;
    }
    return CertifiedLllRoute::fallback;
}

CertifiedLllRoute certified_fmpz_lll(fmpz_mat_struct* basis,
                                     fmpz_mat_struct* transform,
                                     flint::FmpzLllConstRef config) noexcept {
    CertifiedLllRoute route = CertifiedLllRoute::ineligible;
    if (first_phase_eligible(basis, config.raw())) {
        flint::FmpzMat input(::fmpz_mat_nrows(basis), ::fmpz_mat_ncols(basis));
        ::fmpz_mat_set(input.raw(), basis);
        flint::FmpzMat input_transform(0, 0);
        if (transform != nullptr) {
            input_transform = flint::FmpzMat(::fmpz_mat_nrows(transform),
                                             ::fmpz_mat_ncols(transform));
            ::fmpz_mat_set(input_transform.raw(), transform);
        }
        route = lll_first_phase_certified(basis, transform, config);
        if (route == CertifiedLllRoute::certified) {
            return route;
        }
        ::fmpz_mat_set(basis, input.raw());
        if (transform != nullptr) {
            ::fmpz_mat_set(transform, input_transform.raw());
        }
    }
    ::fmpz_lll(basis, transform, config.raw());
    return route;
}

}  // namespace silex::lat::detail
