#include <silex/embedding.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>
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

int test_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);

    silex::NumberField field = field_by_polynomial(polynomial);

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

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_quadratic_real() == 0);
    assert(test_quadratic_complex() == 0);
    assert(test_cubic_trace_norm() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_move_swap_and_clear() == 0);
    assert(test_define_failure_preserves_context() == 0);
    return 0;
}
