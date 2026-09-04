#include <silex/archimedean.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>

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

void element_one_plus_theta(silex::Element& element) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool contains_si(const arb_t value, slong expected) noexcept {
    return sflint::arb_contains_si(value, expected);
}

bool contains_si(sflint::ArbRef value, slong expected) noexcept {
    return sflint::arb_contains_si(value, expected);
}

bool contains_si(const sflint::Arb& value, slong expected) noexcept {
    return sflint::arb_contains_si(value, expected);
}

bool contains_abs_norm(sflint::ArbConstRef value,
                       const silex::Element& element) noexcept {
    sflint::Fmpq norm;
    assert(element.norm(sflint::FmpqRef(norm)));
    sflint::fmpq_abs(norm, norm);
    return sflint::arb_contains_fmpq(value, norm);
}

bool log_sum_overlaps_abs_norm(sflint::ArbVec& logs,
                               const silex::Element& element,
                               slong precision) noexcept {
    sflint::Arb sum;
    for (slong i = 0; i < logs.length(); ++i) {
        sflint::arb_add(sum, sum, logs.data() + i, precision);
    }

    sflint::Fmpq norm;
    sflint::Arb expected;
    assert(element.norm(sflint::FmpqRef(norm)));
    sflint::fmpq_abs(norm, norm);
    sflint::arb_set_fmpq(expected, norm, precision);
    sflint::arb_log(expected, expected, precision);
    return sflint::arb_overlaps(sum, expected);
}

void arb_mat_row_square_sum(sflint::Arb& sum,
                            const sflint::ArbMat& matrix,
                            slong precision) noexcept {
    sflint::Arb term;
    sflint::arb_zero(sum);
    for (slong i = 0; i < sflint::arb_mat_ncols_value(matrix); ++i) {
        const auto entry = sflint::arb_mat_entry_ref(matrix, 0, i);
        sflint::arb_mul(term, entry, entry, precision);
        sflint::arb_add(sum, sum, term, precision);
    }
}

int test_absolute_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::EmbeddingContext embeddings(field);
    silex::Element value(field);
    assert(value.set_si(-7));

    sflint::Arb out;
    sflint::arb_set_si(out, 123);
    assert(archimedean_absolute(sflint::ArbRef(out), embeddings, value,
                                0, silex::ArchAbsMode::plain, 80));
    assert(contains_si(out, 7));

    assert(archimedean_absolute(sflint::ArbRef(out), embeddings, value,
                                0, silex::ArchAbsMode::product, 80));
    assert(contains_si(out, 7));

    sflint::arb_set_si(out, 123);
    assert(!archimedean_absolute(sflint::ArbRef(out), embeddings, value,
                                 1, silex::ArchAbsMode::plain, 80));
    assert(contains_si(out, 123));
    assert(!archimedean_absolute(sflint::ArbRef(out), embeddings, value,
                                 0, silex::ArchAbsMode::plain, 0));
    assert(contains_si(out, 123));
    return 0;
}

int test_absolute_quadratic_and_cubic() {
    sflint::FmpqPoly polynomial;
    sflint::Arb out;
    sflint::Arb product;

    poly_x2_minus(polynomial, -1);
    silex::NumberField complex_quadratic = field_by_polynomial(polynomial);
    silex::EmbeddingContext complex_embeddings(complex_quadratic);
    silex::Element one_plus_i(complex_quadratic);
    element_one_plus_theta(one_plus_i);

    assert(archimedean_absolute(sflint::ArbRef(out), complex_embeddings,
                                one_plus_i, 0, silex::ArchAbsMode::plain, 128));
    sflint::arb_mul(product, out, out, 128);
    assert(contains_si(product, 2));

    assert(archimedean_absolute(sflint::ArbRef(out), complex_embeddings,
                                one_plus_i, 0, silex::ArchAbsMode::product, 128));
    assert(contains_si(out, 2));
    assert(contains_abs_norm(out, one_plus_i));

    poly_x3_minus(polynomial, 2);
    silex::NumberField cubic = field_by_polynomial(polynomial);
    silex::EmbeddingContext cubic_embeddings(cubic);
    silex::Element theta(cubic);
    assert(theta.gen());

    sflint::arb_one(product);
    assert(archimedean_absolute(sflint::ArbRef(out), cubic_embeddings,
                                theta, 0, silex::ArchAbsMode::product, 160));
    sflint::arb_mul(product, product, out, 160);
    assert(archimedean_absolute(sflint::ArbRef(out), cubic_embeddings,
                                theta, 1, silex::ArchAbsMode::product, 160));
    sflint::arb_mul(product, product, out, 160);
    assert(contains_abs_norm(product, theta));
    assert(!archimedean_absolute(sflint::ArbRef(out), cubic_embeddings,
                                 theta, 2, silex::ArchAbsMode::product, 160));
    return 0;
}

int test_log_embedding() {
    sflint::FmpqPoly polynomial;

    poly_x(polynomial);
    silex::NumberField degree_one = field_by_polynomial(polynomial);
    silex::EmbeddingContext embeddings(degree_one);
    silex::Element value(degree_one);
    assert(value.set_si(-7));

    sflint::ArbVec logs1(1);
    assert(logarithmic_embedding(sflint::ArbVecRef(logs1), embeddings,
                                 value, silex::LogEmbeddingMode::plain, 80));
    sflint::Arb expected;
    sflint::arb_log_ui(expected, 7, 80);
    assert(sflint::arb_overlaps(logs1.data() + 0, expected));

    poly_x2_minus(polynomial, -1);
    silex::NumberField complex_quadratic = field_by_polynomial(polynomial);
    silex::EmbeddingContext complex_embeddings(complex_quadratic);
    silex::Element one_plus_i(complex_quadratic);
    element_one_plus_theta(one_plus_i);

    sflint::ArbVec logs2(1);
    assert(logarithmic_embedding(sflint::ArbVecRef(logs2),
                                 complex_embeddings, one_plus_i,
                                 silex::LogEmbeddingMode::plain, 128));
    sflint::arb_mul_2exp_si(expected, logs2.data() + 0, 1);
    assert(logarithmic_embedding(sflint::ArbVecRef(logs2),
                                 complex_embeddings, one_plus_i,
                                 silex::LogEmbeddingMode::product, 128));
    assert(sflint::arb_overlaps(logs2.data() + 0, expected));
    assert(log_sum_overlaps_abs_norm(logs2, one_plus_i, 128));

    poly_x3_minus(polynomial, 2);
    silex::NumberField cubic = field_by_polynomial(polynomial);
    silex::EmbeddingContext cubic_embeddings(cubic);
    silex::Element theta(cubic);
    assert(theta.gen());
    sflint::ArbVec logs_cubic(2);
    assert(logarithmic_embedding(sflint::ArbVecRef(logs_cubic),
                                 cubic_embeddings, theta,
                                 silex::LogEmbeddingMode::product, 160));
    assert(log_sum_overlaps_abs_norm(logs_cubic, theta, 160));

    assert(theta.zero());
    sflint::arb_set_si(logs_cubic.data() + 0, 123);
    sflint::arb_set_si(logs_cubic.data() + 1, 456);
    assert(!logarithmic_embedding(sflint::ArbVecRef(logs_cubic),
                                  cubic_embeddings, theta,
                                  silex::LogEmbeddingMode::product, 160));
    assert(contains_si(logs_cubic.data() + 0, 123));
    assert(contains_si(logs_cubic.data() + 1, 456));

    assert(value.set_si(7));
    sflint::arb_set_si(logs1.data() + 0, 123);
    assert(!logarithmic_embedding(sflint::ArbVecRef(logs1), embeddings,
                                  value, silex::LogEmbeddingMode::product, 0));
    assert(contains_si(logs1.data() + 0, 123));

    assert(!logarithmic_embedding(
            sflint::ArbVecRef(logs1), embeddings, value,
            static_cast<silex::LogEmbeddingMode>(99), 80));
    assert(contains_si(logs1.data() + 0, 123));

    poly_x2_minus(polynomial, 2);
    silex::NumberField real_quadratic = field_by_polynomial(polynomial);
    silex::EmbeddingContext real_embeddings(real_quadratic);
    silex::Element real_zero(real_quadratic);
    assert(real_zero.zero());
    sflint::ArbVec logs_real(2);
    sflint::arb_set_si(logs_real.data() + 0, 123);
    sflint::arb_set_si(logs_real.data() + 1, 456);
    assert(!logarithmic_embedding(sflint::ArbVecRef(logs_real),
                                  real_embeddings, real_zero,
                                  silex::LogEmbeddingMode::product, 128));
    assert(contains_si(logs_real.data() + 0, 123));
    assert(contains_si(logs_real.data() + 1, 456));
    return 0;
}

int test_minkowski_embedding() {
    sflint::FmpqPoly polynomial;
    sflint::Arb sum;

    poly_x(polynomial);
    silex::NumberField degree_one = field_by_polynomial(polynomial);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::Element value(degree_one);
    assert(value.set_si(-7));
    sflint::ArbMat row1(1, 1);

    assert(minkowski_embedding(sflint::ArbMatRef(row1),
                               degree_one_embeddings, value,
                               silex::MinkowskiEmbeddingMode::plain, 80));
    assert(contains_si(sflint::arb_mat_entry_ref(row1, 0, 0), -7));

    poly_x2_minus(polynomial, 2);
    silex::NumberField real_quadratic = field_by_polynomial(polynomial);
    silex::EmbeddingContext real_embeddings(real_quadratic);
    silex::Element theta(real_quadratic);
    assert(theta.gen());
    sflint::ArbMat row2(1, 2);
    assert(minkowski_embedding(sflint::ArbMatRef(row2), real_embeddings,
                               theta, silex::MinkowskiEmbeddingMode::plain, 128));
    arb_mat_row_square_sum(sum, row2, 128);
    assert(contains_si(sum, 4));

    poly_x2_minus(polynomial, -1);
    silex::NumberField complex_quadratic = field_by_polynomial(polynomial);
    silex::EmbeddingContext complex_embeddings(complex_quadratic);
    silex::Element one_plus_i(complex_quadratic);
    element_one_plus_theta(one_plus_i);
    assert(minkowski_embedding(sflint::ArbMatRef(row2), complex_embeddings,
                               one_plus_i, silex::MinkowskiEmbeddingMode::plain,
                               128));
    arb_mat_row_square_sum(sum, row2, 128);
    assert(contains_si(sum, 2));

    assert(minkowski_embedding(sflint::ArbMatRef(row2), complex_embeddings,
                               one_plus_i, silex::MinkowskiEmbeddingMode::weighted,
                               128));
    arb_mat_row_square_sum(sum, row2, 128);
    assert(contains_si(sum, 4));

    sflint::arb_set_si(sflint::arb_mat_entry_ref(row1, 0, 0), 123);
    assert(!minkowski_embedding(sflint::ArbMatRef(row1), complex_embeddings,
                                one_plus_i, silex::MinkowskiEmbeddingMode::plain,
                                128));
    assert(contains_si(sflint::arb_mat_entry_ref(row1, 0, 0), 123));

    assert(minkowski_embedding(sflint::ArbMatRef(row1),
                               degree_one_embeddings, value,
                               silex::MinkowskiEmbeddingMode::plain, 80));
    sflint::arb_set_si(sflint::arb_mat_entry_ref(row1, 0, 0), 123);
    assert(!minkowski_embedding(sflint::ArbMatRef(row1),
                                degree_one_embeddings, value,
                                silex::MinkowskiEmbeddingMode::plain, 0));
    assert(contains_si(sflint::arb_mat_entry_ref(row1, 0, 0), 123));

    assert(!minkowski_embedding(
            sflint::ArbMatRef(row1), degree_one_embeddings, value,
            static_cast<silex::MinkowskiEmbeddingMode>(99), 80));
    assert(contains_si(sflint::arb_mat_entry_ref(row1, 0, 0), 123));
    return 0;
}

int test_nonsquarefree_failure() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::EmbeddingContext embeddings(field);
    silex::Element theta(field);
    assert(theta.gen());

    sflint::Arb out;
    sflint::arb_set_si(out, 123);
    assert(!archimedean_absolute(sflint::ArbRef(out), embeddings, theta,
                                 0, silex::ArchAbsMode::plain, 80));
    assert(contains_si(out, 123));

    sflint::ArbVec logs(1);
    sflint::arb_set_si(logs.data() + 0, 123);
    assert(!logarithmic_embedding(sflint::ArbVecRef(logs), embeddings,
                                  theta, silex::LogEmbeddingMode::plain, 80));
    assert(contains_si(logs.data() + 0, 123));

    sflint::ArbMat row(1, 2);
    sflint::arb_set_si(sflint::arb_mat_entry_ref(row, 0, 0), 123);
    sflint::arb_set_si(sflint::arb_mat_entry_ref(row, 0, 1), 456);
    assert(!minkowski_embedding(sflint::ArbMatRef(row), embeddings, theta,
                                silex::MinkowskiEmbeddingMode::plain, 80));
    assert(contains_si(sflint::arb_mat_entry_ref(row, 0, 0), 123));
    assert(contains_si(sflint::arb_mat_entry_ref(row, 0, 1), 456));
    return 0;
}

}  // namespace

int main() {
    assert(test_absolute_degree_one() == 0);
    assert(test_absolute_quadratic_and_cubic() == 0);
    assert(test_log_embedding() == 0);
    assert(test_minkowski_embedding() == 0);
    assert(test_nonsquarefree_failure() == 0);
    return 0;
}
