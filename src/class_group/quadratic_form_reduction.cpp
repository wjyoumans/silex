#include "quadratic_form_reduction_internal.hpp"

#include <flint/fmpz.h>

#include <algorithm>

namespace silex::detail {
namespace {

void quadratic_form_reduce_basis_step(flint::Fmpz& b,
                                      flint::Fmpz& c,
                                      const flint::Fmpz& a,
                                      const flint::Fmpz& u1,
                                      flint::Fmpz& u2) noexcept {
    flint::Fmpz quotient;
    flint::Fmpz remainder;
    flint::Fmpz twice_a;
    flint::Fmpz half_sum;
    flint::Fmpz change;
    flint::fmpz_mul_2exp(flint::FmpzRef(twice_a),
                         flint::FmpzConstRef(a), 1);
    flint::fmpz_tdiv_qr(
            flint::FmpzRef(quotient), flint::FmpzRef(remainder),
            flint::FmpzConstRef(b), flint::FmpzConstRef(twice_a));
    const int remainder_comparison = flint::fmpz_cmpabs(
            flint::FmpzConstRef(remainder), flint::FmpzConstRef(a));
    if (flint::fmpz_sgn(flint::FmpzConstRef(b)) >= 0) {
        if (remainder_comparison > 0) {
            flint::fmpz_add_ui(flint::FmpzRef(quotient),
                               flint::FmpzConstRef(quotient), 1);
            flint::fmpz_sub(flint::FmpzRef(remainder),
                            flint::FmpzConstRef(remainder),
                            flint::FmpzConstRef(twice_a));
        }
    } else if (remainder_comparison >= 0) {
        flint::fmpz_sub_ui(flint::FmpzRef(quotient),
                           flint::FmpzConstRef(quotient), 1);
        flint::fmpz_add(flint::FmpzRef(remainder),
                        flint::FmpzConstRef(remainder),
                        flint::FmpzConstRef(twice_a));
    }
    flint::fmpz_add(flint::FmpzRef(half_sum), flint::FmpzConstRef(b),
                    flint::FmpzConstRef(remainder));
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(half_sum),
                            flint::FmpzConstRef(half_sum), 1);
    flint::fmpz_mul(flint::FmpzRef(change),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(half_sum));
    flint::fmpz_sub(flint::FmpzRef(c), flint::FmpzConstRef(c),
                    flint::FmpzConstRef(change));
    flint::fmpz_set(flint::FmpzRef(b), flint::FmpzConstRef(remainder));
    flint::fmpz_mul(flint::FmpzRef(change),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(u1));
    flint::fmpz_sub(flint::FmpzRef(u2), flint::FmpzConstRef(u2),
                    flint::FmpzConstRef(change));
}

// Fast reduction of positive forms: reference implementation
// `src/basemath/Qfb.c`, lines 655-843 (`pqfbred_1`, `is_minimal`,
// `pqfbred_iter_1`, `pqfbred_basecase`, `qfb_maxexpi`, `qfb_minexpi`,
// `qfb3_SL2_apply`, `pqfbred_rec`).  The source cites A. Schoenhage, "Fast reduction and composition of binary quadratic
// forms", ISSAC 1991, pp. 128-133.  A form carries its discriminant d, as a
// reference `t_QFB` does; truncated forms (`qfb3`) carry their own b^2 - 4ac.

struct PositiveForm {
    flint::Fmpz a;
    flint::Fmpz b;
    flint::Fmpz c;
    flint::Fmpz d;
};

// Reference `expi(0)` is -HIGHEXPOBIT; any value far below every bit count
// works.
constexpr slong expi_of_zero = -(WORD(1) << 60);

// Reference `expi`: the exponent of the leading bit of |x|.
slong leading_bit_exponent(const fmpz_t x) noexcept {
    if (fmpz_is_zero(x)) {
        return expi_of_zero;
    }
    return static_cast<slong>(fmpz_bits(x)) - 1;
}

// Reference `lti2n`: x < 2^m.
bool less_than_power_of_two(const fmpz_t x, slong m) noexcept {
    return fmpz_sgn(x) < 0 || leading_bit_exponent(x) < m;
}

slong form_max_expi(const fmpz_t a, const fmpz_t b, const fmpz_t c) noexcept {
    return 1 + std::max(leading_bit_exponent(a),
                        std::max(leading_bit_exponent(b),
                                 leading_bit_exponent(c)));
}

slong form_max_expi(const PositiveForm& q) noexcept {
    return form_max_expi(q.a.raw(), q.b.raw(), q.c.raw());
}

slong form_min_expi(const PositiveForm& q) noexcept {
    const slong m = std::min(leading_bit_exponent(q.a.raw()),
                             std::min(leading_bit_exponent(q.b.raw()),
                                      leading_bit_exponent(q.c.raw())));
    return m < 0 ? 0 : m;
}

// Reference `truedivii`: the Euclidean quotient (remainder in [0, |y|)).
void euclidean_quotient(fmpz_t q, const fmpz_t x, const fmpz_t y) noexcept {
    if (fmpz_sgn(y) > 0) {
        fmpz_fdiv_q(q, x, y);
        return;
    }
    flint::Fmpz negated;
    fmpz_neg(negated.raw(), y);
    fmpz_fdiv_q(q, x, negated.raw());
    fmpz_neg(q, q);
}

void transform_set_identity(QuadraticFormTransform& u) noexcept {
    fmpz_one(u.m11.raw());
    fmpz_zero(u.m12.raw());
    fmpz_zero(u.m21.raw());
    fmpz_one(u.m22.raw());
}

// Reference `ZM2_mul`: u = u * v.
void transform_mul_right(QuadraticFormTransform& u,
                         const QuadraticFormTransform& v) noexcept {
    flint::Fmpz n11;
    flint::Fmpz n12;
    flint::Fmpz n21;
    flint::Fmpz n22;
    fmpz_fmma(n11.raw(), u.m11.raw(), v.m11.raw(), u.m12.raw(), v.m21.raw());
    fmpz_fmma(n12.raw(), u.m11.raw(), v.m12.raw(), u.m12.raw(), v.m22.raw());
    fmpz_fmma(n21.raw(), u.m21.raw(), v.m11.raw(), u.m22.raw(), v.m21.raw());
    fmpz_fmma(n22.raw(), u.m21.raw(), v.m12.raw(), u.m22.raw(), v.m22.raw());
    u.m11.swap(n11);
    u.m12.swap(n12);
    u.m21.swap(n21);
    u.m22.swap(n22);
}

// Reference `pqfbred_1`: one normalization step against the smaller outer
// coefficient.  Returns false where the reference raises a bug error.
bool positive_form_step(PositiveForm& q,
                        slong m,
                        QuadraticFormTransform& u) noexcept {
    const bool use_a = fmpz_cmpabs(q.a.raw(), q.c.raw()) < 0;
    flint::Fmpz& lead = use_a ? q.a : q.c;
    flint::Fmpz& other = use_a ? q.c : q.a;
    flint::Fmpz r2;
    flint::Fmpz r;
    flint::Fmpz t;
    flint::Fmpz scratch;
    fmpz_mul_2exp(r2.raw(), lead.raw(), static_cast<ulong>(m + 2));
    fmpz_add(r2.raw(), r2.raw(), q.d.raw());
    const slong e2 = leading_bit_exponent(r2.raw());
    const slong r_exponent =
            fmpz_sgn(r2.raw()) < 0 || e2 < 2 * m + 2 ? m + 1 : e2 >> 1;
    fmpz_one(r.raw());
    fmpz_mul_2exp(r.raw(), r.raw(), static_cast<ulong>(r_exponent));
    fmpz_sub(t.raw(), q.b.raw(), r.raw());
    fmpz_mul_2exp(scratch.raw(), lead.raw(), 1);
    euclidean_quotient(t.raw(), t.raw(), scratch.raw());
    if (fmpz_is_zero(t.raw())) {
        return false;
    }
    // other = other - b t + (lead t) t;  b = b - 2 (lead t).
    flint::Fmpz lead_t;
    fmpz_mul(lead_t.raw(), lead.raw(), t.raw());
    fmpz_submul(other.raw(), q.b.raw(), t.raw());
    fmpz_addmul(other.raw(), lead_t.raw(), t.raw());
    fmpz_mul_2exp(lead_t.raw(), lead_t.raw(), 1);
    fmpz_sub(q.b.raw(), q.b.raw(), lead_t.raw());
    if (use_a) {
        fmpz_submul(u.m12.raw(), u.m11.raw(), t.raw());
        fmpz_submul(u.m22.raw(), u.m21.raw(), t.raw());
    } else {
        fmpz_submul(u.m11.raw(), u.m12.raw(), t.raw());
        fmpz_submul(u.m21.raw(), u.m22.raw(), t.raw());
    }
    return true;
}

// Reference `is_minimal`.
bool positive_form_is_minimal(const PositiveForm& q, slong m) noexcept {
    flint::Fmpz scratch;
    fmpz_sub(scratch.raw(), q.a.raw(), q.b.raw());
    fmpz_add(scratch.raw(), scratch.raw(), q.c.raw());
    if (less_than_power_of_two(scratch.raw(), m)) {
        return true;
    }
    fmpz_mul_2exp(scratch.raw(), q.a.raw(), 1);
    fmpz_sub(scratch.raw(), q.b.raw(), scratch.raw());
    if (!less_than_power_of_two(scratch.raw(), m + 1)) {
        return false;
    }
    fmpz_mul_2exp(scratch.raw(), q.c.raw(), 1);
    fmpz_sub(scratch.raw(), q.b.raw(), scratch.raw());
    return less_than_power_of_two(scratch.raw(), m + 1);
}

// Reference `pqfbred_iter_1`.
bool positive_form_iterate(PositiveForm& q,
                           slong m,
                           QuadraticFormTransform& u) noexcept {
    while (!positive_form_is_minimal(q, m)) {
        if (!positive_form_step(q, m, u)) {
            return false;
        }
    }
    return true;
}

// Reference `qfb3_SL2_apply`: (a, b, c) o M, written into the triple of `out`.
void positive_form_apply(PositiveForm& out,
                         const PositiveForm& q,
                         const QuadraticFormTransform& transform) noexcept {
    const fmpz* x = transform.m11.raw();
    const fmpz* y = transform.m21.raw();
    const fmpz* z = transform.m12.raw();
    const fmpz* t = transform.m22.raw();
    flint::Fmpz first;
    flint::Fmpz second;
    flint::Fmpz scratch;
    // A = x (a x + b y) + c y^2.
    fmpz_fmma(first.raw(), q.a.raw(), x, q.b.raw(), y);
    fmpz_mul(out.a.raw(), x, first.raw());
    fmpz_mul(scratch.raw(), y, y);
    fmpz_addmul(out.a.raw(), q.c.raw(), scratch.raw());
    // B = x (2 a z + b t) + y (2 c t + b z).
    fmpz_mul_2exp(scratch.raw(), q.a.raw(), 1);
    fmpz_fmma(first.raw(), scratch.raw(), z, q.b.raw(), t);
    fmpz_mul_2exp(scratch.raw(), q.c.raw(), 1);
    fmpz_fmma(second.raw(), scratch.raw(), t, q.b.raw(), z);
    fmpz_fmma(out.b.raw(), x, first.raw(), y, second.raw());
    // C = z (a z + b t) + c t^2.
    fmpz_fmma(first.raw(), q.a.raw(), z, q.b.raw(), t);
    fmpz_mul(out.c.raw(), z, first.raw());
    fmpz_mul(scratch.raw(), t, t);
    fmpz_addmul(out.c.raw(), q.c.raw(), scratch.raw());
}

// Reference `pqfbred_rec(Q, m, &U)`: reduces q in place and sets u.
bool positive_form_reduce_recursive(PositiveForm& q,
                                    slong m,
                                    QuadraticFormTransform& u) noexcept {
    const slong n = form_max_expi(q) - m;
    if (n < 170) {
        // pqfbred_basecase.
        transform_set_identity(u);
        return positive_form_iterate(q, m, u);
    }

    if (form_min_expi(q) <= m + 2) {
        transform_set_identity(u);
    } else {
        slong mm = m;
        slong p = 0;
        PositiveForm low;
        PositiveForm reduced;
        if (m > n) {
            mm = n;
            p = m + 1 - n;
            const ulong shift = static_cast<ulong>(p);
            // remi2n and shifti(., -p) truncate toward zero.
            fmpz_tdiv_r_2exp(low.a.raw(), q.a.raw(), shift);
            fmpz_tdiv_r_2exp(low.b.raw(), q.b.raw(), shift);
            fmpz_tdiv_r_2exp(low.c.raw(), q.c.raw(), shift);
            fmpz_tdiv_q_2exp(reduced.a.raw(), q.a.raw(), shift);
            fmpz_tdiv_q_2exp(reduced.b.raw(), q.b.raw(), shift);
            fmpz_tdiv_q_2exp(reduced.c.raw(), q.c.raw(), shift);
            // qfb3: the truncated form carries b^2 - 4ac.
            fmpz_mul(reduced.d.raw(), reduced.a.raw(), reduced.c.raw());
            fmpz_mul_2exp(reduced.d.raw(), reduced.d.raw(), 2);
            fmpz_neg(reduced.d.raw(), reduced.d.raw());
            fmpz_addmul(reduced.d.raw(), reduced.b.raw(), reduced.b.raw());
        } else {
            fmpz_set(reduced.a.raw(), q.a.raw());
            fmpz_set(reduced.b.raw(), q.b.raw());
            fmpz_set(reduced.c.raw(), q.c.raw());
            fmpz_set(reduced.d.raw(), q.d.raw());
        }
        const slong h = mm + (n >> 1);
        if (form_min_expi(reduced) <= h) {
            transform_set_identity(u);
        } else if (!positive_form_reduce_recursive(reduced, h, u)) {
            return false;
        }
        bool going_to_r8 = false;
        while (form_max_expi(reduced) > h) {
            if (positive_form_is_minimal(reduced, mm)) {
                going_to_r8 = true;
                break;
            }
            if (!positive_form_step(reduced, mm, u)) {
                return false;
            }
        }
        if (!going_to_r8) {
            QuadraticFormTransform v;
            if (!positive_form_reduce_recursive(reduced, mm, v)) {
                return false;
            }
            transform_mul_right(u, v);
        }
        if (p > 0) {
            PositiveForm low_image;
            positive_form_apply(low_image, low, u);
            const ulong shift = static_cast<ulong>(p);
            fmpz_mul_2exp(q.a.raw(), reduced.a.raw(), shift);
            fmpz_add(q.a.raw(), q.a.raw(), low_image.a.raw());
            fmpz_mul_2exp(q.b.raw(), reduced.b.raw(), shift);
            fmpz_add(q.b.raw(), q.b.raw(), low_image.b.raw());
            fmpz_mul_2exp(q.c.raw(), reduced.c.raw(), shift);
            fmpz_add(q.c.raw(), q.c.raw(), low_image.c.raw());
        } else {
            q.a.swap(reduced.a);
            q.b.swap(reduced.b);
            q.c.swap(reduced.c);
        }
    }
    return positive_form_iterate(q, m, u);
}

bool quadratic_form_discriminant_above_limit(flint::Fmpz& discriminant,
                                             flint::FmpzConstRef a,
                                             flint::FmpzConstRef b,
                                             flint::FmpzConstRef c) noexcept {
    fmpz_mul(discriminant.raw(), a.raw(), c.raw());
    fmpz_mul_2exp(discriminant.raw(), discriminant.raw(), 2);
    fmpz_neg(discriminant.raw(), discriminant.raw());
    fmpz_addmul(discriminant.raw(), b.raw(), b.raw());
    if (fmpz_sgn(a.raw()) <= 0 || fmpz_sgn(discriminant.raw()) >= 0) {
        return false;
    }
    // Reference `qfbred_imag_av`: 2*qfb_maxexpi(Q) - expi(D) > QFBRED_LIMIT.
    return 2 * form_max_expi(a.raw(), b.raw(), c.raw()) -
                   leading_bit_exponent(discriminant.raw()) >
           quadratic_form_fast_reduction_limit;
}

}  // namespace

bool quadratic_form_uses_fast_reduction(flint::FmpzConstRef a,
                                        flint::FmpzConstRef b,
                                        flint::FmpzConstRef c) noexcept {
    flint::Fmpz discriminant;
    return quadratic_form_discriminant_above_limit(discriminant, a, b, c);
}

bool positive_quadratic_form_partial_reduce(
        flint::Fmpz& out_a,
        flint::Fmpz& out_b,
        flint::Fmpz& out_c,
        QuadraticFormTransform& transform,
        flint::FmpzConstRef a,
        flint::FmpzConstRef b,
        flint::FmpzConstRef c) noexcept {
    if (fmpz_sgn(a.raw()) < 0 || fmpz_sgn(b.raw()) < 0 ||
        fmpz_sgn(c.raw()) < 0) {
        return false;
    }
    PositiveForm q;
    fmpz_set(q.a.raw(), a.raw());
    fmpz_set(q.b.raw(), b.raw());
    fmpz_set(q.c.raw(), c.raw());
    fmpz_mul(q.d.raw(), a.raw(), c.raw());
    fmpz_mul_2exp(q.d.raw(), q.d.raw(), 2);
    fmpz_neg(q.d.raw(), q.d.raw());
    fmpz_addmul(q.d.raw(), b.raw(), b.raw());
    if (fmpz_sgn(q.d.raw()) >= 0 ||
        !positive_form_reduce_recursive(q, 0, transform)) {
        return false;
    }
    out_a.swap(q.a);
    out_b.swap(q.b);
    out_c.swap(q.c);
    return true;
}

bool quadratic_form_reduce_basecase(flint::Fmpz& reduced_a,
                                    flint::Fmpz& reduced_b,
                                    flint::Fmpz& reduced_c,
                                    flint::Fmpz& u1,
                                    flint::Fmpz& v1,
                                    flint::FmpzConstRef input_a,
                                    flint::FmpzConstRef input_b,
                                    flint::FmpzConstRef input_c) noexcept {
    if (flint::fmpz_sgn(input_a) <= 0 || flint::fmpz_sgn(input_c) <= 0) {
        return false;
    }

    // Source: reference implementation `Qfb.c:qfbredsl2_imag_basecase`.  The first
    // column of its returned SL(2,Z) transform is a shortest vector for the
    // positive definite input form.
    flint::Fmpz a;
    flint::Fmpz b;
    flint::Fmpz c;
    flint::Fmpz u2;
    flint::fmpz_set(flint::FmpzRef(a), input_a);
    flint::fmpz_set(flint::FmpzRef(b), input_b);
    flint::fmpz_set(flint::FmpzRef(c), input_c);
    flint::fmpz_one(flint::FmpzRef(u1));
    flint::fmpz_zero(flint::FmpzRef(u2));

    int comparison = flint::fmpz_cmpabs(
            flint::FmpzConstRef(a), flint::FmpzConstRef(b));
    if (comparison < 0) {
        quadratic_form_reduce_basis_step(b, c, a, u1, u2);
    } else if (comparison == 0 && flint::fmpz_sgn(flint::FmpzConstRef(b)) < 0) {
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::fmpz_one(flint::FmpzRef(u2));
    }

    for (;;) {
        comparison = flint::fmpz_cmpabs(
                flint::FmpzConstRef(a), flint::FmpzConstRef(c));
        if (comparison <= 0) {
            break;
        }
        a.swap(c);
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::Fmpz old_u1;
        flint::fmpz_set(flint::FmpzRef(old_u1), flint::FmpzConstRef(u1));
        flint::fmpz_set(flint::FmpzRef(u1), flint::FmpzConstRef(u2));
        flint::fmpz_neg(flint::FmpzRef(u2), flint::FmpzConstRef(old_u1));
        quadratic_form_reduce_basis_step(b, c, a, u1, u2);
    }
    if (comparison == 0 && flint::fmpz_sgn(flint::FmpzConstRef(b)) < 0) {
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::Fmpz old_u1;
        flint::fmpz_set(flint::FmpzRef(old_u1), flint::FmpzConstRef(u1));
        flint::fmpz_set(flint::FmpzRef(u1), flint::FmpzConstRef(u2));
        flint::fmpz_neg(flint::FmpzRef(u2), flint::FmpzConstRef(old_u1));
    }

    flint::Fmpz z;
    flint::Fmpz numerator;
    flint::Fmpz product;
    flint::fmpz_sub(flint::FmpzRef(z), flint::FmpzConstRef(b), input_b);
    if (flint::fmpz_fdiv_ui(flint::FmpzConstRef(z), 2) != 0) {
        return false;
    }
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(z), flint::FmpzConstRef(z), 1);
    flint::fmpz_mul(flint::FmpzRef(numerator), flint::FmpzConstRef(z),
                    flint::FmpzConstRef(u1));
    flint::fmpz_mul(flint::FmpzRef(product), flint::FmpzConstRef(a),
                    flint::FmpzConstRef(u2));
    flint::fmpz_sub(flint::FmpzRef(numerator),
                    flint::FmpzConstRef(numerator),
                    flint::FmpzConstRef(product));
    if (!flint::fmpz_divisible(flint::FmpzConstRef(numerator), input_c)) {
        return false;
    }
    flint::fmpz_divexact(flint::FmpzRef(v1),
                         flint::FmpzConstRef(numerator), input_c);
    reduced_a.swap(a);
    reduced_b.swap(b);
    reduced_c.swap(c);
    return true;
}

bool quadratic_form_reduce(flint::Fmpz& reduced_a,
                           flint::Fmpz& reduced_b,
                           flint::Fmpz& reduced_c,
                           flint::Fmpz& u1,
                           flint::Fmpz& v1,
                           flint::FmpzConstRef a,
                           flint::FmpzConstRef b,
                           flint::FmpzConstRef c) noexcept {
    if (flint::fmpz_sgn(a) <= 0 || flint::fmpz_sgn(c) <= 0) {
        return false;
    }
    flint::Fmpz discriminant;
    if (!quadratic_form_discriminant_above_limit(discriminant, a, b, c)) {
        return quadratic_form_reduce_basecase(reduced_a, reduced_b, reduced_c,
                                              u1, v1, a, b, c);
    }

    // Reference implementation `Qfb.c:qfbred_imag_av` (lines 991-1002) and
    // `qfbredsl2_imag` (lines 886-907): run pqfbred_rec on (a, |b|, c), then
    // the basecase on its output.  For b < 0 the reference multiplies U by
    // S = diag(1, -1) on the left, which gives a transform of determinant -1
    // and so the reduced form of (a, -b, c).  Here the result of pqfbred_rec
    // for (a, -b, c) = (a, b, c) o S is conjugated instead: with
    // (a, -b, c) o U = (A, B, C), (a, b, c) o (S U S) = (A, -B, C), and S U S
    // is U with its off-diagonal entries negated.
    const bool negative_b = flint::fmpz_sgn(b) < 0;
    flint::Fmpz absolute_b;
    flint::fmpz_abs(flint::FmpzRef(absolute_b), b);
    flint::Fmpz partial_a;
    flint::Fmpz partial_b;
    flint::Fmpz partial_c;
    QuadraticFormTransform transform;
    flint::Fmpz w1;
    flint::Fmpz w2;
    if (!positive_quadratic_form_partial_reduce(
                partial_a, partial_b, partial_c, transform, a,
                flint::FmpzConstRef(absolute_b), c)) {
        return quadratic_form_reduce_basecase(reduced_a, reduced_b, reduced_c,
                                              u1, v1, a, b, c);
    }
    if (negative_b) {
        fmpz_neg(partial_b.raw(), partial_b.raw());
        fmpz_neg(transform.m12.raw(), transform.m12.raw());
        fmpz_neg(transform.m21.raw(), transform.m21.raw());
    }
    if (!quadratic_form_reduce_basecase(
                reduced_a, reduced_b, reduced_c, w1, w2,
                flint::FmpzConstRef(partial_a), flint::FmpzConstRef(partial_b),
                flint::FmpzConstRef(partial_c))) {
        return quadratic_form_reduce_basecase(reduced_a, reduced_b, reduced_c,
                                              u1, v1, a, b, c);
    }
    // First column of transform * W.
    fmpz_fmma(u1.raw(), transform.m11.raw(), w1.raw(), transform.m12.raw(),
              w2.raw());
    fmpz_fmma(v1.raw(), transform.m21.raw(), w1.raw(), transform.m22.raw(),
              w2.raw());
    // The caller publishes (u1, v1) as a generator of norm reduced_a, so the
    // fast path's transform is checked directly: unless the input form takes
    // the value reduced_a at (u1, v1), the basecase reduction is used.
    flint::Fmpz value;
    flint::Fmpz partial;
    fmpz_fmma(partial.raw(), a.raw(), u1.raw(), b.raw(), v1.raw());
    fmpz_mul(value.raw(), u1.raw(), partial.raw());
    fmpz_mul(partial.raw(), v1.raw(), v1.raw());
    fmpz_addmul(value.raw(), c.raw(), partial.raw());
    if (!fmpz_equal(value.raw(), reduced_a.raw())) {
        return quadratic_form_reduce_basecase(reduced_a, reduced_b, reduced_c,
                                              u1, v1, a, b, c);
    }
    return true;
}

bool quadratic_form_shortest_vector(flint::Fmpz& u1,
                                    flint::Fmpz& v1,
                                    flint::Fmpz& shortest_value,
                                    flint::FmpzConstRef a,
                                    flint::FmpzConstRef b,
                                    flint::FmpzConstRef c) noexcept {
    flint::Fmpz reduced_b;
    flint::Fmpz reduced_c;
    return quadratic_form_reduce(shortest_value, reduced_b, reduced_c, u1, v1,
                                 a, b, c);
}

}  // namespace silex::detail
