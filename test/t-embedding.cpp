#include <silex/diagnostics.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <flint/acb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include <cassert>
#include <cstring>
#include <utility>

namespace sflint = silex::flint;

namespace {

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// White-box helper: returns a defined field whose FLINT context has been
// replaced by one for the non-squarefree polynomial x^2.  NumberField
// construction rejects x^2, so this is the only way to keep the internal
// non-squarefree failure paths covered.  The field starts as the generic
// (non-quadratic-backend) field defined by x^2 - x - 1 so that no quadratic
// backend data disagrees with the installed polynomial.
silex::NumberField nonsquarefree_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    assert(field.backend_kind() == silex::NumberFieldBackendKind::generic);

    sflint::FmpqPoly square;
    sflint::fmpq_poly_set_coeff_si(square, 2, 1);
    nf_struct* raw = field.raw_flint_field();
    assert(raw != nullptr);
    nf_clear(raw);
    nf_init(raw, square.raw());
    assert(field.degree() == 2);
    return field;
}

void check_square_rejected() noexcept {
    sflint::FmpqPoly square;
    sflint::fmpq_poly_set_coeff_si(square, 2, 1);
    assert(!silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(square)).is_defined());
}

bool contains_si(const acb_t value, slong expected) noexcept {
    return sflint::arb_contains_si(sflint::acb_realref_ptr(value), expected) &&
           sflint::arb_contains_zero(sflint::acb_imagref_ptr(value));
}

bool contains_si(const sflint::Acb& value, slong expected) noexcept {
    return sflint::arb_contains_si(
                   sflint::acb_realref_ptr(value), expected) &&
           sflint::arb_contains_zero(sflint::acb_imagref_ptr(value));
}

bool satisfies_x2_minus(const acb_t root, slong a, slong precision) noexcept {
    sflint::Acb value;
    sflint::acb_mul(value, root, root, precision);
    sflint::acb_sub_si(value, value, a, precision);
    return sflint::acb_contains_zero(value);
}

bool satisfies_x2_minus(
        const sflint::Acb& root, slong a, slong precision) noexcept {
    sflint::Acb value;
    sflint::acb_mul(
            value, sflint::AcbConstRef(root), sflint::AcbConstRef(root), precision);
    sflint::acb_sub_si(value, value, a, precision);
    return sflint::acb_contains_zero(value);
}

bool satisfies_x3_minus(const acb_t root, slong a, slong precision) noexcept {
    sflint::Acb value;
    sflint::acb_pow_ui(value, root, 3, precision);
    sflint::acb_sub_si(value, value, a, precision);
    return sflint::acb_contains_zero(value);
}

bool satisfies_x3_minus(
        const sflint::Acb& root, slong a, slong precision) noexcept {
    sflint::Acb value;
    sflint::acb_pow_ui(value, sflint::AcbConstRef(root), 3, precision);
    sflint::acb_sub_si(value, value, a, precision);
    return sflint::acb_contains_zero(value);
}

// Degree-4 field with two conjugate pairs whose imaginary parts agree to
// about 300 bits:
//
//     f = (x^2 + 2x + 2)(x^2 - 8x + 17) - 2^-300 x
//       = x^4 - 6x^3 + 3x^2 + (18 - 2^-300)x + 34.
//
// The unperturbed roots are A = -1 +/- i and B = 4 +/- i.  To first order the
// perturbation moves a root r by 2^-300 r / f0'(r), which gives
// |Im A| - |Im B| ~ +4.9e-92 (about 2^-302), so the pair with the smaller
// real part has the larger modulus of its imaginary part.
//
// FLINT's _acb_vec_sort_pretty (acb/vec_sort_pretty.c, acb_cmp_pretty) orders
// by |Im| and falls back to the real part only when the |Im| difference ball
// contains zero.  With roots accurate to about 64 bits the difference is
// unresolved and A sorts first; with roots accurate to 512 bits it is resolved
// and B sorts first.  The pre-T-022 refine path re-sorted the refined roots
// with that comparator, so a 64 -> 512 refine swapped the two complex places.
// test_place_order_stable_across_refine checks the comparator flip directly,
// and this test failed on the pre-T-022 embedding code.
silex::NumberField close_imaginary_pairs_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -6);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 34);

    fmpq_t linear;
    fmpq_init(linear);
    fmpz_one(fmpq_denref(linear));
    fmpz_mul_2exp(fmpq_denref(linear), fmpq_denref(linear), 300);
    fmpz_mul_si(fmpq_numref(linear), fmpq_denref(linear), 18);
    fmpz_sub_ui(fmpq_numref(linear), fmpq_numref(linear), 1);
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, linear);
    fmpq_clear(linear);

    return field_by_polynomial(polynomial);
}

// Real quadratic field with two roots 1 +/- sqrt(3) 2^-500.  Rounding the
// constant coefficient 1 - 3 * 2^-1000 to fewer than about 1000 bits merges
// the roots, so refining from 64 to 128 bits cannot isolate them from the
// previous approximations until the working precision reaches 2048 bits
// (checked with acb_poly_find_roots on FLINT 3.6.0).  That is beyond the 8x
// refine cap, so the refine must fall back to full isolation.
silex::NumberField close_real_roots_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -2);

    fmpq_t constant;
    fmpq_init(constant);
    fmpz_one(fmpq_denref(constant));
    fmpz_mul_2exp(fmpq_denref(constant), fmpq_denref(constant), 1000);
    fmpz_sub_ui(fmpq_numref(constant), fmpq_denref(constant), 3);
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), 0, constant);
    fmpq_clear(constant);

    return field_by_polynomial(polynomial);
}

struct RefineProfileCounters {
    slong find_roots = 0;
    slong full_roots_fallback = 0;
    slong events = 0;
};

void refine_profile_callback(void* user,
                             silex::DiagnosticsModule,
                             silex::ProfileEvent event,
                             const char*,
                             const char* label) noexcept {
    RefineProfileCounters* counters =
            static_cast<RefineProfileCounters*>(user);
    ++counters->events;
    if (label == nullptr) {
        return;
    }
    if (event == silex::ProfileEvent::begin_scope &&
        std::strcmp(label,
                    "element.embedding.refine_from_initial.find_roots") == 0) {
        ++counters->find_roots;
    }
    if (event == silex::ProfileEvent::event &&
        std::strcmp(label, "element.embedding.refine.full_roots_fallback") ==
                0) {
        ++counters->full_roots_fallback;
    }
}

int test_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.is_defined());
    assert(embeddings.degree() == 1);
    assert(embeddings.refine(64));
    assert(embeddings.is_set());
    assert(embeddings.precision() == 64);
    assert(embeddings.refine(32));
    assert(embeddings.precision() == 64);

    silex::Element value(field);
    assert(value.set_si(7));

    sflint::Acb out;
    assert(embeddings.evaluate(sflint::AcbRef(out), value, 0, 64));
    assert(contains_si(out, 7));

    sflint::acb_set_si(out, -123);
    assert(!embeddings.evaluate(sflint::AcbRef(out), value, 1, 64));
    assert(sflint::acb_equal_si(out, -123));
    assert(!embeddings.evaluate(sflint::AcbRef(out), value, -1, 64));
    assert(sflint::acb_equal_si(out, -123));
    assert(!embeddings.evaluate(sflint::AcbRef(out), value, 0, 0));
    assert(sflint::acb_equal_si(out, -123));

    sflint::AcbVec all(1);
    assert(embeddings.evaluate_all(sflint::AcbVecRef(all), value, 64));
    assert(contains_si(all.data() + 0, 7));
    return 0;
}

int test_quadratic_real() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(128));
    assert(embeddings.signature().r1() == 2);
    assert(embeddings.signature().r2() == 0);

    silex::Element theta(field);
    assert(theta.gen());

    sflint::AcbVec values(2);
    assert(embeddings.evaluate_all(sflint::AcbVecRef(values), theta, 128));
    assert(satisfies_x2_minus(values.data() + 0, 2, 128));
    assert(satisfies_x2_minus(values.data() + 1, 2, 128));
    assert(sflint::arb_contains_zero(acb_imagref(values.data() + 0)));
    assert(sflint::arb_contains_zero(acb_imagref(values.data() + 1)));
    return 0;
}

int test_quadratic_complex() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -1);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(128));
    assert(embeddings.signature().r1() == 0);
    assert(embeddings.signature().r2() == 1);

    sflint::Acb root0;
    sflint::Acb root1;
    assert(embeddings.get_root(sflint::AcbRef(root0), 0));
    assert(embeddings.get_root(sflint::AcbRef(root1), 1));
    assert(satisfies_x2_minus(root0, -1, 128));
    assert(satisfies_x2_minus(root1, -1, 128));
    assert(sflint::arb_contains_zero(sflint::acb_realref_ptr(root0)));
    assert(sflint::arb_contains_zero(sflint::acb_realref_ptr(root1)));
    assert(sflint::arb_is_positive(sflint::acb_imagref_ptr(root0)));
    assert(sflint::arb_is_negative(sflint::acb_imagref_ptr(root1)));
    return 0;
}

int test_cubic_trace_norm() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Element theta(field);
    assert(theta.gen());

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(160));
    assert(embeddings.signature().r1() == 1);
    assert(embeddings.signature().r2() == 1);

    sflint::AcbVec values(3);
    assert(embeddings.evaluate_all(sflint::AcbVecRef(values), theta, 160));
    for (slong i = 0; i < 3; ++i) {
        assert(satisfies_x3_minus(values.data() + i, 2, 160));
    }
    assert(sflint::arb_contains_zero(
            sflint::acb_imagref_ptr(values.data() + 0)));
    assert(sflint::arb_is_positive(
            sflint::acb_imagref_ptr(values.data() + 1)));
    assert(sflint::arb_is_negative(
            sflint::acb_imagref_ptr(values.data() + 2)));

    sflint::Acb sum;
    sflint::acb_add(sum, values.data() + 0, values.data() + 1, 160);
    sflint::acb_add(sum, sum, values.data() + 2, 160);

    sflint::Fmpq trace;
    assert(theta.trace(sflint::FmpqRef(trace)));
    assert(sflint::acb_contains_fmpq(sum, trace));

    sflint::Acb product;
    sflint::acb_mul(product, values.data() + 0, values.data() + 1, 160);
    sflint::acb_mul(product, product, values.data() + 2, 160);

    sflint::Fmpq norm;
    assert(theta.norm(sflint::FmpqRef(norm)));
    assert(sflint::acb_contains_fmpq(product, norm));

    silex::Element theta_plus_one(field);
    silex::Element one(field);
    assert(one.set_si(1));
    assert(theta_plus_one.add(theta, one));
    assert(embeddings.evaluate_all(
            sflint::AcbVecRef(values), theta_plus_one, 160));
    for (slong i = 0; i < 3; ++i) {
        sflint::Acb shifted;
        ::acb_sub_si(shifted.raw(), values.data() + i, 1, 160);
        assert(satisfies_x3_minus(shifted, 2, 160));
    }

    assert(embeddings.evaluate_all(
            sflint::AcbVecRef(values), theta_plus_one, 256));
    for (slong i = 0; i < 3; ++i) {
        sflint::Acb shifted;
        ::acb_sub_si(shifted.raw(), values.data() + i, 1, 256);
        assert(satisfies_x3_minus(shifted, 2, 256));
    }
    return 0;
}

int test_undefined_field_failure() {
    check_square_rejected();

    silex::NumberField undefined;
    silex::EmbeddingContext embeddings(undefined);
    assert(!embeddings.is_defined());
    assert(!embeddings.refine(64));
    assert(!embeddings.is_set());

    silex::Element theta(undefined);
    assert(!theta.gen());

    sflint::Acb out;
    sflint::acb_set_si(out, -123);
    assert(!embeddings.evaluate(sflint::AcbRef(out), theta, 0, 64));
    assert(sflint::acb_equal_si(out, -123));

    sflint::AcbVec values(2);
    sflint::acb_set_si(values.data() + 0, -55);
    sflint::acb_set_si(values.data() + 1, -66);
    assert(!embeddings.evaluate_all(sflint::AcbVecRef(values), theta, 64));
    assert(sflint::acb_equal_si(values.data() + 0, -55));
    assert(sflint::acb_equal_si(values.data() + 1, -66));

    assert(!embeddings.get_root(sflint::AcbRef(out), 0));
    assert(sflint::acb_equal_si(out, -123));
    return 0;
}

int test_failure_preserves_output() {
    silex::NumberField field = nonsquarefree_field();

    silex::EmbeddingContext embeddings(field);
    assert(!embeddings.refine(64));
    assert(!embeddings.is_set());

    silex::Element theta(field);
    assert(theta.gen());

    sflint::Acb out;
    sflint::acb_set_si(out, -123);
    assert(!embeddings.evaluate(sflint::AcbRef(out), theta, 0, 64));
    assert(sflint::acb_equal_si(out, -123));

    sflint::AcbVec values(2);
    sflint::acb_set_si(values.data() + 0, -55);
    sflint::acb_set_si(values.data() + 1, -66);
    assert(!embeddings.evaluate_all(sflint::AcbVecRef(values), theta, 64));
    assert(sflint::acb_equal_si(values.data() + 0, -55));
    assert(sflint::acb_equal_si(values.data() + 1, -66));

    assert(!embeddings.get_root(sflint::AcbRef(out), 0));
    assert(sflint::acb_equal_si(out, -123));
    return 0;
}

int test_move_swap_and_clear() {
    sflint::FmpqPoly quadratic_polynomial;
    poly_x2_minus(quadratic_polynomial, 2);

    silex::NumberField quadratic = field_by_polynomial(quadratic_polynomial);

    silex::EmbeddingContext quadratic_embeddings(quadratic);
    assert(quadratic_embeddings.refine(128));

    silex::EmbeddingContext moved(std::move(quadratic_embeddings));
    assert(moved.is_defined());
    assert(moved.parent() != nullptr &&
           moved.parent()->has_same_data(quadratic));
    assert(moved.degree() == 2);
    assert(moved.is_set());
    assert(moved.precision() == 128);
    assert(!quadratic_embeddings.is_defined());

    sflint::Acb root;
    assert(moved.get_root(sflint::AcbRef(root), 0));
    assert(satisfies_x2_minus(root, 2, 128));

    sflint::FmpqPoly cubic_polynomial;
    poly_x3_minus(cubic_polynomial, 2);

    silex::NumberField cubic = field_by_polynomial(cubic_polynomial);

    silex::EmbeddingContext cubic_embeddings(cubic);
    assert(cubic_embeddings.refine(160));

    swap(moved, cubic_embeddings);
    assert(moved.parent() != nullptr &&
           moved.parent()->has_same_data(cubic));
    assert(moved.degree() == 3);
    assert(moved.precision() == 160);
    assert(cubic_embeddings.parent() != nullptr &&
           cubic_embeddings.parent()->has_same_data(quadratic));
    assert(cubic_embeddings.degree() == 2);
    assert(cubic_embeddings.precision() == 128);

    assert(moved.get_root(sflint::AcbRef(root), 0));
    assert(satisfies_x3_minus(root, 2, 160));
    assert(cubic_embeddings.get_root(sflint::AcbRef(root), 0));
    assert(satisfies_x2_minus(root, 2, 128));

    silex::EmbeddingContext assigned;
    assigned = std::move(moved);
    assert(assigned.parent() != nullptr &&
           assigned.parent()->has_same_data(cubic));
    assert(assigned.degree() == 3);
    assert(assigned.is_set());
    assert(!moved.is_defined());

    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(assigned.degree() == 0);
    assert(!assigned.is_set());
    assert(assigned.precision() == 0);
    assert(assigned.signature().degree() == 0);
    return 0;
}

int test_define_failure_preserves_context() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -1);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(96));

    silex::NumberField undefined;
    assert(!embeddings.define(undefined));
    assert(embeddings.parent() != nullptr &&
           embeddings.parent()->has_same_data(field));
    assert(embeddings.degree() == 2);
    assert(embeddings.is_set());
    assert(embeddings.precision() == 96);

    sflint::Acb root;
    assert(embeddings.get_root(sflint::AcbRef(root), 0));
    assert(satisfies_x2_minus(root, -1, 96));
    return 0;
}

int test_signature_before_refine() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);
    silex::NumberField cubic = field_by_polynomial(polynomial);

    silex::EmbeddingContext embeddings(cubic);
    assert(embeddings.is_defined());
    assert(!embeddings.is_set());
    assert(embeddings.signature().r1() == 1);
    assert(embeddings.signature().r2() == 1);
    assert(embeddings.refine(64));
    assert(embeddings.signature().r1() == 1);
    assert(embeddings.signature().r2() == 1);

    silex::EmbeddingContext quartic(close_imaginary_pairs_field());
    assert(!quartic.is_set());
    assert(quartic.signature().r1() == 0);
    assert(quartic.signature().r2() == 2);

    silex::EmbeddingContext deferred;
    assert(deferred.signature().degree() == 0);
    assert(deferred.define(cubic));
    assert(!deferred.is_set());
    assert(deferred.signature().r1() == 1);
    assert(deferred.signature().r2() == 1);

    // The signature is recomputed for the new parent on redefinition.
    poly_x2_minus(polynomial, -1);
    assert(deferred.define(field_by_polynomial(polynomial)));
    assert(deferred.signature().r1() == 0);
    assert(deferred.signature().r2() == 1);

    silex::EmbeddingContext undefined;
    assert(undefined.signature().degree() == 0);
    return 0;
}

int test_place_order_stable_across_refine() {
    silex::NumberField field = close_imaginary_pairs_field();
    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(64));

    sflint::AcbVec low(4);
    for (slong i = 0; i < 4; ++i) {
        assert(embeddings.get_root(sflint::AcbRef(low.data() + i), i));
    }

    assert(embeddings.refine(512));
    assert(embeddings.precision() == 512);
    assert(embeddings.signature().r1() == 0);
    assert(embeddings.signature().r2() == 2);

    sflint::AcbVec high(4);
    for (slong i = 0; i < 4; ++i) {
        assert(embeddings.get_root(sflint::AcbRef(high.data() + i), i));
        assert(::acb_rel_accuracy_bits(high.data() + i) >= 512);
        // Place i still refers to the same complex root.
        assert(::acb_overlaps(high.data() + i, low.data() + i) != 0);
        for (slong j = 0; j < 4; ++j) {
            if (j != i) {
                assert(::acb_overlaps(high.data() + i, low.data() + j) == 0);
            }
        }
    }
    // Complex places keep the positive-imaginary root first and its conjugate
    // second.
    for (slong pair = 0; pair < 2; ++pair) {
        acb_srcptr positive = high.data() + 2 * pair;
        acb_srcptr negative = high.data() + 2 * pair + 1;
        assert(::arb_is_positive(acb_imagref(positive)) != 0);
        sflint::Acb conjugate;
        ::acb_conj(conjugate.raw(), positive);
        assert(::acb_equal(conjugate.raw(), negative) != 0);
    }

    // The input really is adversarial for the naive comparator: sorting the
    // 64-bit and 512-bit root vectors with _acb_vec_sort_pretty puts
    // different conjugate pairs first.
    sflint::AcbVec low_sorted(4);
    sflint::AcbVec high_sorted(4);
    ::_acb_vec_set(low_sorted.data(), low.data(), 4);
    ::_acb_vec_set(high_sorted.data(), high.data(), 4);
    ::_acb_vec_sort_pretty(low_sorted.data(), 4);
    ::_acb_vec_sort_pretty(high_sorted.data(), 4);
    assert(::arb_is_negative(acb_realref(low_sorted.data() + 0)) != 0);
    assert(::arb_is_positive(acb_realref(high_sorted.data() + 0)) != 0);

    // A further refine keeps the order as well.
    assert(embeddings.refine(2048));
    for (slong i = 0; i < 4; ++i) {
        sflint::Acb root;
        assert(embeddings.get_root(sflint::AcbRef(root), i));
        assert(::acb_overlaps(root.raw(), low.data() + i) != 0);
    }
    return 0;
}

int test_refine_precision_cap() {
    silex::NumberField field = close_real_roots_field();
    silex::EmbeddingContext embeddings(field);
    assert(embeddings.signature().r1() == 2);
    assert(embeddings.refine(64));

    sflint::AcbVec low(2);
    for (slong i = 0; i < 2; ++i) {
        assert(embeddings.get_root(sflint::AcbRef(low.data() + i), i));
    }
    assert(::acb_overlaps(low.data() + 0, low.data() + 1) == 0);

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    RefineProfileCounters counters;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::element),
            refine_profile_callback, &counters);

    assert(embeddings.refine(128, &diagnostics));
    assert(embeddings.precision() == 128);
    for (slong i = 0; i < 2; ++i) {
        sflint::Acb root;
        assert(embeddings.get_root(sflint::AcbRef(root), i));
        assert(::acb_rel_accuracy_bits(root.raw()) >= 128);
        assert(::arb_is_zero(acb_imagref(root.raw())) != 0);
        assert(::acb_overlaps(root.raw(), low.data() + i) != 0);
        assert(::acb_overlaps(root.raw(), low.data() + (1 - i)) == 0);
    }

#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // Working precisions 128, 256, 512 and 1024 = 8 * 128 are tried from the
    // previous roots; the uncapped loop would go on to 2048.
    assert(counters.events > 0);
    assert(counters.find_roots == 4);
    assert(counters.full_roots_fallback == 1);
#else
    assert(counters.events == 0);
#endif
    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_quadratic_real() == 0);
    assert(test_quadratic_complex() == 0);
    assert(test_cubic_trace_norm() == 0);
    assert(test_undefined_field_failure() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_move_swap_and_clear() == 0);
    assert(test_define_failure_preserves_context() == 0);
    assert(test_signature_before_refine() == 0);
    assert(test_place_order_stable_across_refine() == 0);
    assert(test_refine_precision_cap() == 0);
    return 0;
}
