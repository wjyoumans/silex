#include "class_group/quadratic_form_reduction_internal.hpp"

#include <silex/flint/fmpz.hpp>

#include <flint/fmpz.h>

#include <cstdio>

namespace {
namespace sflint = silex::flint;
namespace detail = silex::detail;

int failures = 0;

void check(bool condition, const char* message, int index) noexcept {
    if (!condition) {
        std::fprintf(stderr, "FAIL [%d]: %s\n", index, message);
        ++failures;
    }
}

ulong next_random(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) + UWORD(1442695040888963407);
    return state >> 11;
}

struct Form {
    sflint::Fmpz a;
    sflint::Fmpz b;
    sflint::Fmpz c;
};

// out = q o M, with (q o M)(x, y) = q(m11 x + m12 y, m21 x + m22 y).
void apply_transform(Form& out,
                     const Form& q,
                     const fmpz_t m11,
                     const fmpz_t m12,
                     const fmpz_t m21,
                     const fmpz_t m22) noexcept {
    sflint::Fmpz first;
    sflint::Fmpz second;
    Form result;
    // A = q(m11, m21).
    fmpz_mul(first.raw(), q.a.raw(), m11);
    fmpz_addmul(first.raw(), q.b.raw(), m21);
    fmpz_mul(result.a.raw(), first.raw(), m11);
    fmpz_mul(second.raw(), m21, m21);
    fmpz_addmul(result.a.raw(), q.c.raw(), second.raw());
    // C = q(m12, m22).
    fmpz_mul(first.raw(), q.a.raw(), m12);
    fmpz_addmul(first.raw(), q.b.raw(), m22);
    fmpz_mul(result.c.raw(), first.raw(), m12);
    fmpz_mul(second.raw(), m22, m22);
    fmpz_addmul(result.c.raw(), q.c.raw(), second.raw());
    // B = 2 a m11 m12 + b (m11 m22 + m12 m21) + 2 c m21 m22.
    fmpz_mul(first.raw(), m11, m12);
    fmpz_mul(first.raw(), first.raw(), q.a.raw());
    fmpz_mul(second.raw(), m21, m22);
    fmpz_addmul(first.raw(), second.raw(), q.c.raw());
    fmpz_mul_2exp(result.b.raw(), first.raw(), 1);
    fmpz_fmma(first.raw(), m11, m22, m12, m21);
    fmpz_addmul(result.b.raw(), first.raw(), q.b.raw());
    out.a.swap(result.a);
    out.b.swap(result.b);
    out.c.swap(result.c);
}

void form_value(fmpz_t out, const Form& q, const fmpz_t x, const fmpz_t y) {
    sflint::Fmpz scratch;
    fmpz_mul(out, q.a.raw(), x);
    fmpz_addmul(out, q.b.raw(), y);
    fmpz_mul(out, out, x);
    fmpz_mul(scratch.raw(), y, y);
    fmpz_addmul(out, q.c.raw(), scratch.raw());
}

// A positive definite form scale * (3, 1, k) o M, where M is a product of
// `steps` random matrices [[t, -1], [1, 0]] with 51-bit t.  The scale mimics
// the content N(I) of an ideal norm form.  The primitive part has
// discriminant 1 - 12k < -4, so its only automorphisms are +-1, and its
// reduced form (3, 1, k) is not properly equivalent to (3, -1, k).
void random_form(Form& out,
                 ulong& state,
                 slong steps,
                 slong scale_bits) noexcept {
    sflint::Fmpz scale;
    fmpz_one(scale.raw());
    for (slong bits = 0; bits < scale_bits; bits += 62) {
        fmpz_mul_2exp(scale.raw(), scale.raw(), 62);
        fmpz_add_ui(scale.raw(), scale.raw(), next_random(state) >> 2);
    }
    Form q;
    fmpz_mul_ui(q.a.raw(), scale.raw(), 3);
    fmpz_set(q.b.raw(), scale.raw());
    fmpz_mul_ui(q.c.raw(), scale.raw(), 1000003 + (next_random(state) % 1000));
    sflint::Fmpz t;
    sflint::Fmpz minus_one;
    sflint::Fmpz one;
    sflint::Fmpz zero;
    fmpz_set_si(minus_one.raw(), -1);
    fmpz_one(one.raw());
    for (slong i = 0; i < steps; ++i) {
        fmpz_set_ui(t.raw(), (next_random(state) >> 2) + 1);
        apply_transform(q, q, t.raw(), minus_one.raw(), one.raw(), zero.raw());
    }
    out.a.swap(q.a);
    out.b.swap(q.b);
    out.c.swap(q.c);
}

bool vectors_equal_up_to_sign(const fmpz_t x0,
                              const fmpz_t y0,
                              const fmpz_t x1,
                              const fmpz_t y1,
                              bool& negated) noexcept {
    if (fmpz_equal(x0, x1) && fmpz_equal(y0, y1)) {
        negated = false;
        return true;
    }
    sflint::Fmpz nx;
    sflint::Fmpz ny;
    fmpz_neg(nx.raw(), x1);
    fmpz_neg(ny.raw(), y1);
    negated = true;
    return fmpz_equal(x0, nx.raw()) && fmpz_equal(y0, ny.raw());
}

bool fast(const Form& q) noexcept {
    return detail::quadratic_form_uses_fast_reduction(
            sflint::FmpzConstRef(q.a), sflint::FmpzConstRef(q.b),
            sflint::FmpzConstRef(q.c));
}

// The dispatching reduction against the basecase: same reduced form, the
// same first column up to sign, and a first column of value A.  Below the
// switch the first column is identical.
void check_against_basecase(const Form& q, int index, int& negated_count) {
    Form basecase;
    Form reduced;
    sflint::Fmpz x0;
    sflint::Fmpz y0;
    sflint::Fmpz x1;
    sflint::Fmpz y1;
    sflint::Fmpz shortest;
    sflint::Fmpz value;
    const bool basecase_ok = detail::quadratic_form_reduce_basecase(
            basecase.a, basecase.b, basecase.c, x0, y0,
            sflint::FmpzConstRef(q.a), sflint::FmpzConstRef(q.b),
            sflint::FmpzConstRef(q.c));
    const bool reduce_ok = detail::quadratic_form_reduce(
            reduced.a, reduced.b, reduced.c, x1, y1,
            sflint::FmpzConstRef(q.a), sflint::FmpzConstRef(q.b),
            sflint::FmpzConstRef(q.c));
    check(basecase_ok && reduce_ok, "reduction succeeds", index);
    if (!basecase_ok || !reduce_ok) {
        return;
    }
    check(fmpz_equal(basecase.a.raw(), reduced.a.raw()) &&
                  fmpz_equal(basecase.b.raw(), reduced.b.raw()) &&
                  fmpz_equal(basecase.c.raw(), reduced.c.raw()),
          "reduced form matches the basecase", index);
    bool negated = false;
    check(vectors_equal_up_to_sign(x0.raw(), y0.raw(), x1.raw(), y1.raw(),
                                   negated),
          "first column matches the basecase up to sign", index);
    if (negated) {
        ++negated_count;
        check(fast(q), "a sign change happens only above the switch", index);
    }
    form_value(value.raw(), q, x1.raw(), y1.raw());
    check(fmpz_equal(value.raw(), reduced.a.raw()),
          "the first column represents the reduced leading coefficient",
          index);
    check(detail::quadratic_form_shortest_vector(
                  x0, y0, shortest, sflint::FmpzConstRef(q.a),
                  sflint::FmpzConstRef(q.b), sflint::FmpzConstRef(q.c)) &&
                  fmpz_equal(x0.raw(), x1.raw()) &&
                  fmpz_equal(y0.raw(), y1.raw()) &&
                  fmpz_equal(shortest.raw(), reduced.a.raw()),
          "shortest vector agrees with the dispatching reduction", index);
}

// pqfbred_rec output: (a, b, c) o U with det U = 1 and positive outer
// coefficients, far below the input size.
void check_partial_reduction(const Form& q, int index) {
    Form input;
    fmpz_set(input.a.raw(), q.a.raw());
    fmpz_abs(input.b.raw(), q.b.raw());
    fmpz_set(input.c.raw(), q.c.raw());
    Form partial;
    Form image;
    detail::QuadraticFormTransform u;
    const bool ok = detail::positive_quadratic_form_partial_reduce(
            partial.a, partial.b, partial.c, u,
            sflint::FmpzConstRef(input.a), sflint::FmpzConstRef(input.b),
            sflint::FmpzConstRef(input.c));
    check(ok, "partial reduction succeeds", index);
    if (!ok) {
        return;
    }
    sflint::Fmpz det;
    fmpz_fmms(det.raw(), u.m11.raw(), u.m22.raw(), u.m12.raw(), u.m21.raw());
    check(fmpz_is_one(det.raw()), "transform has determinant one", index);
    apply_transform(image, input, u.m11.raw(), u.m12.raw(), u.m21.raw(),
                    u.m22.raw());
    check(fmpz_equal(image.a.raw(), partial.a.raw()) &&
                  fmpz_equal(image.b.raw(), partial.b.raw()) &&
                  fmpz_equal(image.c.raw(), partial.c.raw()),
          "output is the input transformed by U", index);
    check(fmpz_sgn(partial.a.raw()) > 0 && fmpz_sgn(partial.c.raw()) > 0,
          "output outer coefficients are positive", index);
    check(!fast(partial), "output is below the switch", index);
}

void test_threshold() {
    Form small;
    fmpz_set_ui(small.a.raw(), 5);
    fmpz_set_ui(small.b.raw(), 3);
    fmpz_set_ui(small.c.raw(), 7);
    check(!fast(small), "small form uses the basecase", 0);
    Form indefinite;
    fmpz_set_ui(indefinite.a.raw(), 1);
    fmpz_one(indefinite.b.raw());
    fmpz_mul_2exp(indefinite.b.raw(), indefinite.b.raw(), 20000);
    fmpz_set_ui(indefinite.c.raw(), 1);
    check(!fast(indefinite), "indefinite form uses the basecase", 0);
    // (1, 0, 2^k): 2 * (k + 1) - (k + 2) = k, so the switch is at k = 9000.
    Form edge;
    fmpz_one(edge.a.raw());
    fmpz_zero(edge.b.raw());
    fmpz_one(edge.c.raw());
    fmpz_mul_2exp(edge.c.raw(), edge.c.raw(), 9000);
    check(!fast(edge), "form at the switch uses the basecase", 0);
    fmpz_mul_2exp(edge.c.raw(), edge.c.raw(), 1);
    check(fast(edge), "form just above the switch uses the fast route", 0);
    int negated = 0;
    check_against_basecase(edge, 0, negated);
}

void test_below_switch_is_unchanged() {
    ulong state = 17;
    int negated = 0;
    for (int i = 0; i < 40; ++i) {
        Form q;
        random_form(q, state, 2 + i % 20, 64 * (i % 5));
        check(!fast(q), "small random form uses the basecase", 100 + i);
        check_against_basecase(q, 100 + i, negated);
        fmpz_neg(q.b.raw(), q.b.raw());
        check_against_basecase(q, 100 + i, negated);
    }
    check(negated == 0, "no sign change below the switch", 100);
}

void test_above_switch() {
    ulong state = 4242;
    int negated = 0;
    int compared = 0;
    int sign_sensitive = 0;
    const slong steps[] = {60, 100, 220, 400};
    const slong scale_bits[] = {0, 512, 3000};
    int index = 1000;
    for (slong step_count : steps) {
        for (slong bits : scale_bits) {
            for (int repeat = 0; repeat < 3; ++repeat, ++index) {
                Form q;
                random_form(q, state, step_count, bits);
                check(fast(q), "large random form uses the fast route",
                      index);
                check_partial_reduction(q, index);
                check_against_basecase(q, index, negated);
                Form flipped;
                fmpz_set(flipped.a.raw(), q.a.raw());
                fmpz_neg(flipped.b.raw(), q.b.raw());
                fmpz_set(flipped.c.raw(), q.c.raw());
                check_against_basecase(flipped, index, negated);
                compared += 2;
                // A reduced form (A, B, C) with 0 < B < A < C is not
                // properly equivalent to (A, -B, C): the b < 0 branch must
                // not return the class of (a, -b, c).
                Form reduced;
                sflint::Fmpz x;
                sflint::Fmpz y;
                if (detail::quadratic_form_reduce_basecase(
                            reduced.a, reduced.b, reduced.c, x, y,
                            sflint::FmpzConstRef(q.a),
                            sflint::FmpzConstRef(q.b),
                            sflint::FmpzConstRef(q.c)) &&
                    fmpz_sgn(reduced.b.raw()) != 0 &&
                    fmpz_cmp(reduced.b.raw(), reduced.a.raw()) != 0 &&
                    fmpz_cmp(reduced.a.raw(), reduced.c.raw()) != 0) {
                    ++sign_sensitive;
                }
            }
        }
    }
    check(compared == 72, "all large cases compared", 1000);
    check(sign_sensitive > 0, "some case separates (a, b, c) and (a, -b, c)",
          1000);
    std::printf("above switch: %d forms, %d first columns negated\n", compared,
                negated);
}

}  // namespace

int main() {
    test_threshold();
    test_below_switch_is_unchanged();
    test_above_switch();
    return failures == 0 ? 0 : 1;
}
