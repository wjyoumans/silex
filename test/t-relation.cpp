#include <silex/abelian_group.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_element.hpp>
#include <silex/relation.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

bool set_rational(silex::Element& element,
                  slong numerator,
                  ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator, denominator);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    return element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool set_linear_rational(silex::Element& element,
                         slong constant_num,
                         ulong constant_den,
                         slong linear_num,
                         ulong linear_den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, constant_num, constant_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    sflint::fmpq_set_si(coefficient, linear_num, linear_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool relation_row_is_si(const silex::Relation& relation,
                        slong index,
                        slong value) noexcept {
    sflint::FmpzMat row(1, relation.length());
    return relation.exponents(sflint::FmpzMatRef(row)) &&
           mat_entry_is_si(row, 0, index, value);
}

bool relation_row_equals(const silex::Relation& relation,
                         const slong* expected,
                         slong length) noexcept {
    if (relation.length() != length) {
        return false;
    }
    sflint::FmpzMat row(1, length);
    if (!relation.exponents(sflint::FmpzMatRef(row))) {
        return false;
    }
    for (slong i = 0; i < length; ++i) {
        if (!mat_entry_is_si(row, 0, i, expected[i])) {
            return false;
        }
    }
    return true;
}

bool set_generator_from_order_coordinates(silex::Element& out,
                                          const silex::Order& order,
                                          slong c0,
                                          slong c1) noexcept {
    sflint::FmpzMat coordinates(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, 0), c0);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, 1), c1);

    silex::OrderElement order_element(order);
    return order_element.set_coordinates(sflint::FmpzMatConstRef(coordinates)) &&
           order_element.get_element(out);
}

bool set_element_from_order_coordinates(silex::Element& out,
                                        const silex::Order& order,
                                        const sflint::FmpzMat& coordinates)
        noexcept {
    silex::OrderElement order_element(order);
    return order_element.set_coordinates(
                   sflint::FmpzMatConstRef(coordinates)) &&
           order_element.get_element(out);
}

int test_degree_one_relation() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::Relation relation(base);
    silex::Relation copy(base);
    silex::Relation other(base);
    silex::Element alpha(field);
    silex::Element generator(field);
    sflint::FmpzMat row(1, 1);

    assert(relation.is_defined());
    assert(!relation.is_set());
    assert(!relation.exponents(sflint::FmpzMatRef(row)));

    assert(alpha.set_si(1));
    assert(relation.set_generator(alpha));
    assert(relation.exponents(sflint::FmpzMatRef(row)));
    assert(mat_entry_is_si(row, 0, 0, 0));
    assert(relation.generator(generator));
    assert(generator.equal_si(1));

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 1));
    {
        const sflint::FmpzMatConstRef row_view = relation.exponents_ref();
        assert(sflint::fmpz_mat_nrows(row_view) == 1);
        assert(sflint::fmpz_mat_ncols(row_view) == 1);
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(row_view, 0, 0), 1));
    }

    assert(copy.set(relation));
    assert(relation_row_is_si(copy, 0, 1));
    assert(copy.generator(generator));
    assert(generator.equal_si(2));

    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 2));

    assert(set_rational(alpha, 1, 2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -1));

    swap(relation, other);
    assert(!relation.is_set());
    assert(other.is_set());
    assert(relation_row_is_si(other, 0, -1));

    return 0;
}

int test_failure_preserves_relation() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Relation relation(base);
    silex::Element alpha(field);
    silex::Element generator(field);

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 1));

    assert(alpha.set_si(3));
    assert(!relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 1));
    assert(relation.generator(generator));
    assert(generator.equal_si(2));

    assert(set_rational(alpha, 3, 2));
    assert(!relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 1));

    assert(alpha.set_si(0));
    assert(!relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 1));

    assert(set_rational(alpha, 1, 2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -1));

    assert(set_rational(alpha, 3, 2));
    assert(!relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -1));

    return 0;
}

int test_known_row_reuses_defined_relation() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Relation relation(base);
    silex::Element alpha(field);
    silex::Element generator(field);
    sflint::FmpzMat row(1, 1);

    assert(alpha.set_si(2));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    assert(silex::detail::set_relation_from_known_row(
            relation, base, alpha, sflint::FmpzMatConstRef(row)));
    assert(relation.generator(generator));
    assert(generator.equal_si(2));
    assert(relation_row_is_si(relation, 0, 1));

    assert(alpha.set_si(4));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 2);
    assert(silex::detail::set_relation_from_known_row(
            relation, base, alpha, sflint::FmpzMatConstRef(row)));
    assert(relation.generator(generator));
    assert(generator.equal_si(4));
    assert(relation_row_is_si(relation, 0, 2));

    return 0;
}

int test_set_generator_with_norm() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::Relation relation(base);
    silex::Element alpha(field);
    sflint::Fmpq norm;

    assert(alpha.set_si(4));
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(relation.set_generator_with_norm(alpha,
                                            sflint::FmpqConstRef(norm)));
    assert(relation_row_is_si(relation, 0, 2));

    assert(alpha.set_si(3));
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(!relation.set_generator_with_norm(alpha,
                                             sflint::FmpqConstRef(norm)));
    assert(relation_row_is_si(relation, 0, 2));

    return 0;
}

int test_set_generator_with_integral_coordinates_and_norm() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());
    assert(order.is_equation_order());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::Relation normal(base);
    silex::Relation coordinate_relation(base);
    silex::Relation polynomial_relation(base);
    silex::Element alpha(field);
    sflint::Fmpq norm;
    sflint::FmpzPoly coordinate_polynomial;
    sflint::FmpzMat coordinates(1, 2);
    sflint::FmpzMat normal_row(1, base.length());
    sflint::FmpzMat coordinate_row(1, base.length());
    sflint::FmpzMat polynomial_row(1, base.length());

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, 0), 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, 1), 0);
    assert(set_element_from_order_coordinates(alpha, order, coordinates));
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(normal.set_generator_with_norm(alpha,
                                          sflint::FmpqConstRef(norm)));
    assert(silex::detail::set_relation_from_integral_coordinates_and_norm(
            coordinate_relation, alpha,
            sflint::FmpzMatConstRef(coordinates),
            sflint::FmpqConstRef(norm)));
    sflint::fmpz_poly_set_coeff_si(coordinate_polynomial, 0, 2);
    assert(silex::detail::set_relation_from_integral_coordinates_and_norm(
            polynomial_relation, alpha,
            sflint::FmpzMatConstRef(coordinates),
            sflint::FmpqConstRef(norm), &coordinate_polynomial));
    assert(normal.exponents(sflint::FmpzMatRef(normal_row)));
    assert(coordinate_relation.exponents(
            sflint::FmpzMatRef(coordinate_row)));
    assert(polynomial_relation.exponents(
            sflint::FmpzMatRef(polynomial_row)));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(normal_row),
            sflint::FmpzMatConstRef(coordinate_row)));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(normal_row),
            sflint::FmpzMatConstRef(polynomial_row)));
    assert(relation_row_is_si(coordinate_relation, 0, 2));

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, 0), 3);
    assert(set_element_from_order_coordinates(alpha, order, coordinates));
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(!silex::detail::set_relation_from_integral_coordinates_and_norm(
            coordinate_relation, alpha,
            sflint::FmpzMatConstRef(coordinates),
            sflint::FmpqConstRef(norm)));
    assert(relation_row_is_si(coordinate_relation, 0, 2));

    sflint::FmpzMat bad_coordinates(2, 2);
    assert(!silex::detail::set_relation_from_integral_coordinates_and_norm(
            coordinate_relation, alpha,
            sflint::FmpzMatConstRef(bad_coordinates),
            sflint::FmpqConstRef(norm)));
    assert(relation_row_is_si(coordinate_relation, 0, 2));

    return 0;
}

int test_quadratic_relation() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::Relation relation(base);
    silex::Element alpha(field);

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, 2));

    assert(set_rational(alpha, 1, 2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -2));

    return 0;
}

int test_repeated_residue_valuation_relation() {
    {
        silex::NumberField field = degree_one_field();
        silex::Order order = silex::test::equation_order(field);

        silex::FactorBase base(order);
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        assert(base.build(sflint::FmpzConstRef(bound)));
        assert(base.length() == 1);

        silex::Relation relation(base);
        silex::Element alpha(field);

        assert(alpha.set_si(8));
        assert(relation.set_generator(alpha));
        assert(relation_row_is_si(relation, 0, 3));
    }

    {
        silex::NumberField field = quadratic_field(2);
        silex::Order order = silex::test::equation_order(field);
        assert(order.is_maximal());

        silex::FactorBase base(order);
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        assert(base.build(sflint::FmpzConstRef(bound)));
        assert(base.length() == 1);

        silex::Relation relation(base);
        silex::Element alpha(field);

        assert(alpha.set_si(2));
        assert(relation.set_generator(alpha));
        assert(relation_row_is_si(relation, 0, 2));
    }

    return 0;
}

int test_nonintegral_equation_relation_denominator_split() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());
    assert(order.is_equation_order());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 1);

    silex::Relation relation(base);
    silex::Element alpha(field);

    assert(set_linear_rational(alpha, 1, 2, 1, 2));
    assert(relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -2));

    assert(set_rational(alpha, 1, 3));
    assert(!relation.set_generator(alpha));
    assert(relation_row_is_si(relation, 0, -2));

    return 0;
}

int test_half_integral_quadratic_relation_fallback() {
    {
        silex::NumberField field = quadratic_field(-47);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(maximal.is_maximal());
        assert(!maximal.is_equation_order());

        silex::FactorBase base(maximal);
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 3));
        assert(base.build(sflint::FmpzConstRef(bound)));
        assert(base.length() == 4);

        silex::Relation relation(base);
        silex::Element alpha(field);

        const slong row0[] = {0, 2, 0, 1};
        assert(set_generator_from_order_coordinates(alpha, maximal, 1, -1));
        assert(relation.set_generator(alpha));
        assert(relation_row_equals(relation, row0, 4));

        const slong row1[] = {1, 0, 0, 2};
        assert(set_generator_from_order_coordinates(alpha, maximal, 2, 1));
        assert(relation.set_generator(alpha));
        assert(relation_row_equals(relation, row1, 4));
    }

    {
        silex::NumberField field = quadratic_field(5);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(maximal.is_maximal());
        assert(!maximal.is_equation_order());

        silex::FactorBase base(maximal);
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        assert(base.build(sflint::FmpzConstRef(bound)));
        assert(base.length() == 1);

        silex::Relation relation(base);
        silex::Element alpha(field);

        const slong row0[] = {2};
        assert(set_generator_from_order_coordinates(alpha, maximal, 4, 4));
        assert(relation.set_generator(alpha));
        assert(relation_row_equals(relation, row0, 1));

        const slong row1[] = {1};
        assert(set_generator_from_order_coordinates(alpha, maximal, 4, -2));
        assert(relation.set_generator(alpha));
        assert(relation_row_equals(relation, row1, 1));
    }

    return 0;
}

int test_generic_cubic_relation() {
    silex::NumberField field = cubic_field();
    silex::Order equation = silex::test::equation_order(field);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::FactorBase base(maximal);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Relation relation(base);
    silex::Element alpha(field);
    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));

    sflint::FmpzMat row(1, base.length());
    assert(relation.exponents(sflint::FmpzMatRef(row)));

    silex::PrimeIdeal prime(maximal);
    sflint::Fmpz p;
    for (slong i = 0; i < base.length(); ++i) {
        assert(base.prime(prime, i));
        assert(prime.rational_prime(sflint::FmpzRef(p)));
        const slong expected = sflint::fmpz_equal_si(p, 2)
                                     ? prime.ramification_index()
                                     : 0;
        assert(mat_entry_is_si(row, 0, i, expected));
    }

    return 0;
}

int test_relation_matrix() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    silex::FactorBase copied_base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(copied_base.set(base));

    silex::Relation relation(base);
    silex::Relation copied_base_relation(copied_base);
    silex::RelationMatrix matrix(base);
    silex::RelationMatrix copy(base);
    silex::Element alpha(field);
    silex::Element generator(field);
    sflint::FmpzMat row(1, 1);
    sflint::FmpzMat rows(2, 1);
    slong first_nonzero = -1;

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(matrix.length() == 1);
    assert(matrix.row_first_nonzero(first_nonzero, 0));
    assert(first_nonzero == 0);
    assert(matrix.row(sflint::FmpzMatRef(row), 0));
    assert(mat_entry_is_si(row, 0, 0, 1));
    assert(matrix.row_equal(sflint::FmpzMatConstRef(row), 0));
    sflint::fmpz_zero(sflint::fmpz_mat_entry(row, 0, 0));
    assert(!matrix.row_equal(sflint::FmpzMatConstRef(row), 0));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    assert(matrix.generator(generator, 0));
    assert(generator.equal_si(2));
    const silex::Element* stored_generator = matrix.generator_at(0);
    assert(stored_generator != nullptr);
    assert(stored_generator->equal_si(2));
    assert(matrix.generator_at(-1) == nullptr);
    assert(matrix.generator_at(1) == nullptr);

    assert(copy.set(matrix));
    assert(copy.length() == 1);
    assert(copy.row(sflint::FmpzMatRef(row), 0));
    assert(mat_entry_is_si(row, 0, 0, 1));
    assert(copy.generator(generator, 0));
    assert(generator.equal_si(2));

    assert(set_rational(alpha, 1, 2));
    assert(relation.set_generator(alpha));
    assert(matrix.generator(generator, 0));
    assert(!generator.equal(alpha));
    assert(copy.generator(generator, 0));
    assert(!generator.equal(alpha));
    assert(matrix.append(relation));
    assert(matrix.length() == 2);
    assert(matrix.rows(sflint::FmpzMatRef(rows)));
    assert(mat_entry_is_si(rows, 0, 0, 1));
    assert(mat_entry_is_si(rows, 1, 0, -1));
    assert(matrix.row(sflint::FmpzMatRef(row), 1));
    assert(matrix.row_equal(sflint::FmpzMatConstRef(row), 1));
    assert(!matrix.row_equal(sflint::FmpzMatConstRef(rows), 0));

    assert(alpha.set_si(1));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(matrix.length() == 3);
    assert(matrix.row_first_nonzero(first_nonzero, 2));
    assert(first_nonzero == matrix.ncols());
    sflint::FmpzMat rows_with_zero(3, 1);
    assert(matrix.rows(sflint::FmpzMatRef(rows_with_zero)));
    assert(mat_entry_is_si(rows_with_zero, 2, 0, 0));
    sflint::FmpzMat zero_row(1, matrix.ncols());
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(zero_row));
    assert(matrix.row_equal(sflint::FmpzMatConstRef(zero_row), 2));

    assert(alpha.set_si(4));
    assert(copied_base_relation.set_generator(alpha));
    assert(matrix.append(copied_base_relation));
    assert(matrix.length() == 4);

    silex::FactorBase larger_base(order);
    assert(set_fmpz_si(bound, 3));
    assert(larger_base.build(sflint::FmpzConstRef(bound)));
    silex::Relation larger_base_relation(larger_base);
    assert(larger_base_relation.set_generator(alpha));
    assert(!matrix.append(larger_base_relation));
    assert(matrix.length() == 4);

    return 0;
}

int test_relation_matrix_to_abelian_group() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base2(order);
    silex::FactorBase base3(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base2.build(sflint::FmpzConstRef(bound)));
    assert(set_fmpz_si(bound, 3));
    assert(base3.build(sflint::FmpzConstRef(bound)));

    silex::Relation relation(base2);
    silex::RelationMatrix matrix(base2);
    silex::FiniteAbelianGroup group;
    silex::Element alpha(field);
    sflint::Fmpz order_out;
    sflint::FmpzMat rel5(1, 1);

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(matrix.to_abelian_group(group));
    assert(group.invariant_count() == 0);
    assert(group.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_is_one(order_out));

    matrix = silex::RelationMatrix(base2);
    assert(matrix.is_defined());
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(matrix.to_abelian_group(group));
    assert(group.invariant_count() == 1);
    assert(group.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 2));

    assert(alpha.set_si(2));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(alpha.set_si(8));
    assert(relation.set_generator(alpha));
    assert(matrix.append(relation));
    assert(matrix.length() == 3);
    assert(matrix.to_abelian_group(group));
    assert(group.invariant_count() == 0);
    assert(group.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_is_one(order_out));

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel5, 0, 0), 5);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(rel5)));
    silex::Relation relation3(base3);
    silex::RelationMatrix rank_deficient(base3);
    assert(alpha.set_si(2));
    assert(relation3.set_generator(alpha));
    assert(rank_deficient.append(relation3));
    assert(!rank_deficient.to_abelian_group(group));
    assert(group.order(sflint::FmpzRef(order_out)));
    assert(sflint::fmpz_equal_si(order_out, 5));

    return 0;
}

int test_cross_parent_copy_swap() {
    silex::NumberField left_field = degree_one_field();
    silex::NumberField right_field = degree_one_field();
    silex::Order left_order = silex::test::equation_order(left_field);
    silex::Order right_order = silex::test::equation_order(right_field);

    silex::FactorBase left_base(left_order);
    silex::FactorBase right_base(right_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(left_base.build(sflint::FmpzConstRef(bound)));
    assert(right_base.build(sflint::FmpzConstRef(bound)));

    silex::Relation left(left_base);
    silex::Relation right(right_base);
    silex::Element left_alpha(left_field);
    silex::Element right_alpha(right_field);
    silex::Element left_generator(left_field);
    silex::Element right_generator(right_field);
    assert(left_alpha.set_si(2));
    assert(right_alpha.set_si(4));
    assert(left.set_generator(left_alpha));
    assert(right.set_generator(right_alpha));

    assert(left.set(right));
    assert(left.factor_base() != nullptr && left.factor_base()->equal(right_base));
    assert(left.generator(right_generator));
    assert(right_generator.equal_si(4));
    assert(!left.generator(left_generator));

    swap(left, right);
    assert(right.factor_base() != nullptr && right.factor_base()->equal(right_base));
    assert(right.generator(right_generator));
    assert(right_generator.equal_si(4));

    return 0;
}

silex::Relation local_relation() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Element alpha(field);
    assert(alpha.set_si(2));
    silex::Relation relation(base);
    assert(relation.set_generator(alpha));
    return relation;
}

silex::RelationMatrix local_relation_matrix() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(base.build(sflint::FmpzConstRef(bound)));

    silex::Element alpha(field);
    assert(alpha.set_si(2));
    silex::Relation relation(base);
    assert(relation.set_generator(alpha));

    silex::RelationMatrix matrix(base);
    assert(matrix.append(relation));
    return matrix;
}

int test_keeps_factor_base_alive() {
    silex::Relation relation = local_relation();
    assert(relation.is_defined());
    assert(relation.is_set());
    assert(relation.factor_base() != nullptr);
    assert(relation.parent() != nullptr);
    assert(relation.parent()->parent() != nullptr);
    assert(relation.length() == 1);
    assert(relation_row_is_si(relation, 0, 1));

    silex::Element generator(*relation.parent()->parent());
    assert(relation.generator(generator));
    assert(generator.equal_si(2));

    silex::RelationMatrix matrix = local_relation_matrix();
    assert(matrix.is_defined());
    assert(matrix.factor_base() != nullptr);
    assert(matrix.parent() != nullptr);
    assert(matrix.parent()->parent() != nullptr);
    assert(matrix.length() == 1);
    assert(matrix.ncols() == 1);

    sflint::FmpzMat row(1, 1);
    assert(matrix.row(sflint::FmpzMatRef(row), 0));
    assert(mat_entry_is_si(row, 0, 0, 1));

    silex::Element matrix_generator(*matrix.parent()->parent());
    assert(matrix.generator(matrix_generator, 0));
    assert(matrix_generator.equal_si(2));
    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one_relation() == 0);
    assert(test_failure_preserves_relation() == 0);
    assert(test_known_row_reuses_defined_relation() == 0);
    assert(test_set_generator_with_norm() == 0);
    assert(test_set_generator_with_integral_coordinates_and_norm() == 0);
    assert(test_quadratic_relation() == 0);
    assert(test_repeated_residue_valuation_relation() == 0);
    assert(test_nonintegral_equation_relation_denominator_split() == 0);
    assert(test_half_integral_quadratic_relation_fallback() == 0);
    assert(test_generic_cubic_relation() == 0);
    assert(test_relation_matrix() == 0);
    assert(test_relation_matrix_to_abelian_group() == 0);
    assert(test_cross_parent_copy_swap() == 0);
    assert(test_keeps_factor_base_alive() == 0);
    return 0;
}
