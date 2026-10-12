#pragma once

#include <silex/flint/fmpz.hpp>

namespace silex::detail {

// A 2 x 2 integer matrix [[m11, m12], [m21, m22]] acting on binary quadratic
// forms by substitution: (q o M)(x, y) = q(m11 x + m12 y, m21 x + m22 y).
struct QuadraticFormTransform {
    flint::Fmpz m11;
    flint::Fmpz m12;
    flint::Fmpz m21;
    flint::Fmpz m22;
};

// Reference implementation `Qfb.c:QFBRED_LIMIT`.  A positive definite form
// (a, b, c) of discriminant D is reduced by the fast route when
// 2 * qfb_maxexpi(a, b, c) - expi(D) exceeds this limit.
inline constexpr slong quadratic_form_fast_reduction_limit = 9000;

// True when (a, b, c) is positive definite (a > 0, b^2 - 4ac < 0) and lies
// above `quadratic_form_fast_reduction_limit`.
bool quadratic_form_uses_fast_reduction(flint::FmpzConstRef a,
                                        flint::FmpzConstRef b,
                                        flint::FmpzConstRef c) noexcept;

// Port of the reference implementation `Qfb.c:pqfbred_rec(Q, 0, &U)` for a
// form with nonnegative coefficients a, b, c and negative discriminant: sets
// (out_a, out_b, out_c) = (a, b, c) o transform with transform in SL(2, Z).
// The output is close to reduced, not necessarily reduced.  Returns false when
// the input has a negative coefficient or a nonnegative discriminant, or when
// an internal step that the reference reports as a bug is reached.
bool positive_quadratic_form_partial_reduce(
        flint::Fmpz& out_a,
        flint::Fmpz& out_b,
        flint::Fmpz& out_c,
        QuadraticFormTransform& transform,
        flint::FmpzConstRef a,
        flint::FmpzConstRef b,
        flint::FmpzConstRef c) noexcept;

// Port of the reference implementation `Qfb.c:qfbredsl2_imag_basecase`: sets
// the reduced form (reduced_a, reduced_b, reduced_c) of the positive definite
// form (a, b, c) and the first column (u1, v1) of the SL(2, Z) transform to
// it.  Returns false unless a > 0 and c > 0.
bool quadratic_form_reduce_basecase(flint::Fmpz& reduced_a,
                                    flint::Fmpz& reduced_b,
                                    flint::Fmpz& reduced_c,
                                    flint::Fmpz& u1,
                                    flint::Fmpz& v1,
                                    flint::FmpzConstRef a,
                                    flint::FmpzConstRef b,
                                    flint::FmpzConstRef c) noexcept;

// Same contract as `quadratic_form_reduce_basecase`.  Forms above
// `quadratic_form_fast_reduction_limit` are first brought close to reduced by
// `positive_quadratic_form_partial_reduce`, as the reference implementation
// `Qfb.c:qfbred_imag_av` and `qfbredsl2_imag` do; the reduced form is the
// same, and (u1, v1) can be the negative of the basecase vector.
bool quadratic_form_reduce(flint::Fmpz& reduced_a,
                           flint::Fmpz& reduced_b,
                           flint::Fmpz& reduced_c,
                           flint::Fmpz& u1,
                           flint::Fmpz& v1,
                           flint::FmpzConstRef a,
                           flint::FmpzConstRef b,
                           flint::FmpzConstRef c) noexcept;

// Sets (u1, v1) to a shortest vector of the positive definite form (a, b, c)
// and shortest_value to its value, through `quadratic_form_reduce`.
bool quadratic_form_shortest_vector(flint::Fmpz& u1,
                                    flint::Fmpz& v1,
                                    flint::Fmpz& shortest_value,
                                    flint::FmpzConstRef a,
                                    flint::FmpzConstRef b,
                                    flint::FmpzConstRef c) noexcept;

}  // namespace silex::detail
