#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include "order/order_internal.hpp"
#include "test_support.hpp"

#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_vec.h>
#include <flint/flint.h>

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

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

void mat_entry_si(sflint::FmpqMat& matrix,
                  slong row,
                  slong col,
                  slong value) noexcept {
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(matrix, row, col), value, 1);
}

void mat_entry_frac_si(sflint::FmpqMat& matrix,
                       slong row,
                       slong col,
                       slong num,
                       slong den) noexcept {
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(matrix, row, col), num, den);
}

bool fmpq_entry_is_si(const sflint::FmpqMat& matrix,
                      slong row,
                      slong col,
                      slong value) noexcept {
    return sflint::fmpq_equal_si(
            sflint::fmpq_mat_entry(matrix, row, col), value);
}

bool fmpq_entry_is_frac_si(const sflint::FmpqMat& matrix,
        slong row,
        slong col,
        slong numerator,
        ulong denominator) noexcept {
    sflint::Fmpq expected;
    sflint::fmpq_set_si(expected, numerator, denominator);
    return sflint::fmpq_equal(sflint::fmpq_mat_entry(matrix, row, col),
            sflint::FmpqConstRef(expected));
}

bool fmpz_mat_entry_is_si(const sflint::FmpzMat& matrix,
                          slong row,
                          slong col,
                          slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

void element_from_coeff(silex::Element& element,
                        slong index,
                        slong num,
                        slong den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_set_si(coeff, num, den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, index, coeff);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

int test_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order = silex::test::equation_order(field);
    assert(order.has_basis());
    assert(order.is_equation_order());
    assert(order.maximality_known());
    assert(order.is_maximal());

    sflint::FmpzMat trace(1, 1);
    assert(order.trace_matrix(sflint::FmpzMatRef(trace)));
    assert(fmpz_mat_entry_is_si(trace, 0, 0, 1));

    sflint::Fmpz disc;
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 1));

    silex::Element value(field);
    assert(value.set_si(7));
    assert(order.contains(value));

    sflint::FmpqMat coords(1, 1);
    assert(order.coordinates(sflint::FmpqMatRef(coords), value));
    assert(fmpq_entry_is_si(coords, 0, 0, 7));
    auto owned_basis = order.basis();
    assert(owned_basis.has_value());
    assert(fmpq_mat_nrows(owned_basis->raw()) == 1);
    assert(fmpq_mat_ncols(owned_basis->raw()) == 1);
    assert(fmpq_entry_is_si(*owned_basis, 0, 0, 1));
    auto owned_coords = order.coordinates(value);
    assert(owned_coords.has_value());
    assert(fmpq_mat_nrows(owned_coords->raw()) == 1);
    assert(fmpq_mat_ncols(owned_coords->raw()) == 1);
    assert(fmpq_entry_is_si(*owned_coords, 0, 0, 7));

    element_from_coeff(value, 0, 1, 2);
    assert(!order.contains(value));
    assert(order.coordinates(value).has_value());

    silex::Order maximal(field);
    sflint::Fmpz idx;
    assert(maximal.maximal_order(order));
    assert(order_index(sflint::FmpzRef(idx), order, maximal));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(order_index(sflint::FmpzRef(idx), maximal, order));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 1));
    assert(maximal.maximality_known());
    assert(maximal.is_maximal());

    sflint::Fmpz p;
    sflint::fmpz_set_ui(p, 7);
    silex::Order pmax(field);
    assert(pmax.pmaximal_overorder(order, sflint::FmpzConstRef(p)));
    assert(order_index(sflint::FmpzRef(idx), order, pmax));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(order_index(sflint::FmpzRef(idx), pmax, order));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(pmax.maximality_known());
    assert(pmax.is_maximal());
    return 0;
}

int test_quadratic_equation_order() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order = silex::test::equation_order(field);

    sflint::FmpzMat trace(2, 2);
    assert(order.trace_matrix(sflint::FmpzMatRef(trace)));
    assert(fmpz_mat_entry_is_si(trace, 0, 0, 2));
    assert(fmpz_mat_entry_is_si(trace, 0, 1, 0));
    assert(fmpz_mat_entry_is_si(trace, 1, 0, 0));
    assert(fmpz_mat_entry_is_si(trace, 1, 1, 4));

    sflint::Fmpz disc;
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));

    silex::Element theta(field);
    assert(theta.gen());
    sflint::FmpqMat coords(1, 2);
    assert(order.coordinates(sflint::FmpqMatRef(coords), theta));
    assert(fmpq_entry_is_si(coords, 0, 0, 0));
    assert(fmpq_entry_is_si(coords, 0, 1, 1));

    element_from_coeff(theta, 1, 1, 2);
    assert(!order.contains(theta));
    return 0;
}

int test_cubic_equation_order() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order = silex::test::equation_order(field);

    sflint::FmpzMat trace(3, 3);
    assert(order.trace_matrix(sflint::FmpzMatRef(trace)));
    assert(fmpz_mat_entry_is_si(trace, 0, 0, 3));
    assert(fmpz_mat_entry_is_si(trace, 0, 1, 0));
    assert(fmpz_mat_entry_is_si(trace, 0, 2, 0));
    assert(fmpz_mat_entry_is_si(trace, 1, 1, 0));
    assert(fmpz_mat_entry_is_si(trace, 1, 2, 6));
    assert(fmpz_mat_entry_is_si(trace, 2, 1, 6));
    assert(fmpz_mat_entry_is_si(trace, 2, 2, 0));

    sflint::Fmpz disc;
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -108));

    silex::Element value(field);
    element_from_coeff(value, 3, 1, 1);
    sflint::FmpqMat coords(1, 3);
    assert(order.coordinates(sflint::FmpqMatRef(coords), value));
    assert(fmpq_entry_is_si(coords, 0, 0, 2));
    assert(fmpq_entry_is_si(coords, 0, 1, 0));
    assert(fmpq_entry_is_si(coords, 0, 2, 0));
    return 0;
}

int test_explicit_basis_index_and_table() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order equation = silex::test::equation_order(field);

    silex::Order suborder(field);
    sflint::FmpqMat basis(2, 2);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_si(basis, 1, 1, 2);
    assert(suborder.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(suborder.maximality_known());
    sflint::Fmpz conductor;
    assert(suborder.quadratic_conductor(sflint::FmpzRef(conductor)));
    assert(sflint::fmpz_equal_si(conductor, 2));
    assert(!suborder.is_maximal());

    sflint::FmpqMat got(2, 2);
    assert(suborder.get_basis(sflint::FmpqMatRef(got)));
    assert(sflint::fmpq_mat_equal(got, basis));

    sflint::FmpzMat trace(2, 2);
    assert(suborder.trace_matrix(sflint::FmpzMatRef(trace)));
    assert(fmpz_mat_entry_is_si(trace, 0, 0, 2));
    assert(fmpz_mat_entry_is_si(trace, 1, 1, 16));

    sflint::Fmpz disc;
    assert(suborder.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    sflint::Fmpz idx;
    assert(order_index(sflint::FmpzRef(idx), suborder, equation));
    assert(sflint::fmpz_equal_si(idx, 2));

    sflint::fmpz_set_si(idx, 99);
    assert(!order_index(sflint::FmpzRef(idx), equation, suborder));
    assert(sflint::fmpz_equal_si(idx, 99));

    silex::Element two_theta(field);
    element_from_coeff(two_theta, 1, 2, 1);
    assert(suborder.contains(two_theta));
    sflint::FmpqMat coords(1, 2);
    assert(suborder.coordinates(sflint::FmpqMatRef(coords), two_theta));
    assert(fmpq_entry_is_si(coords, 0, 0, 0));
    assert(fmpq_entry_is_si(coords, 0, 1, 1));

    silex::Element theta(field);
    assert(theta.gen());
    assert(!suborder.contains(theta));

    sflint::FmpzMat table(4, 2);
    assert(equation.multiplication_table(sflint::FmpzMatRef(table)));
    assert(fmpz_mat_entry_is_si(table, 0, 0, 1));
    assert(fmpz_mat_entry_is_si(table, 1, 1, 1));
    assert(fmpz_mat_entry_is_si(table, 2, 1, 1));
    assert(fmpz_mat_entry_is_si(table, 3, 0, 2));

    assert(suborder.multiplication_table(sflint::FmpzMatRef(table)));
    assert(fmpz_mat_entry_is_si(table, 0, 0, 1));
    assert(fmpz_mat_entry_is_si(table, 1, 1, 1));
    assert(fmpz_mat_entry_is_si(table, 2, 1, 1));
    assert(fmpz_mat_entry_is_si(table, 3, 0, 8));
    return 0;
}

int test_undefined_order_quadratic_conductor() {
    silex::Order undefined;
    assert(!undefined.is_defined());
    sflint::Fmpz conductor;
    sflint::fmpz_set_si(sflint::FmpzRef(conductor), 17);
    assert(!undefined.quadratic_conductor(sflint::FmpzRef(conductor)));
    assert(sflint::fmpz_equal_si(conductor, 17));
    return 0;
}

int test_failure_preserves_order() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_si(basis, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    sflint::Fmpz disc;
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    sflint::FmpqMat wrong_shape(1, 2);
    assert(!order.set_basis(sflint::FmpqMatConstRef(wrong_shape)));
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    sflint::fmpq_mat_zero(basis);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_si(basis, 1, 0, 1);
    assert(!order.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    sflint::fmpq_mat_zero(basis);
    mat_entry_si(basis, 0, 0, 2);
    mat_entry_si(basis, 1, 1, 1);
    assert(!order.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    sflint::fmpq_mat_zero(basis);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(!order.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(order.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 32));

    return 0;
}

int test_copy_and_swap() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order source = silex::test::equation_order(field);

    silex::Order copy(field);
    assert(copy.set(source));
    sflint::Fmpz disc;
    assert(copy.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));

    silex::Order other(field);
    swap(copy, other);
    assert(other.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));
    return 0;
}

int test_move_clear_and_redefine() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order source = silex::test::equation_order(field);

    sflint::Fmpz disc;
    assert(source.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));

    silex::Order moved(std::move(source));
    assert(moved.is_defined());
    assert(moved.parent() != nullptr && moved.parent()->has_same_data(field));
    assert(moved.has_basis());
    assert(moved.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));
    assert(!source.is_defined());
    assert(source.parent() == nullptr);
    assert(source.degree() == 0);
    sflint::Fmpz moved_from_conductor;
    sflint::fmpz_set_si(sflint::FmpzRef(moved_from_conductor), 17);
    assert(!source.quadratic_conductor(
            sflint::FmpzRef(moved_from_conductor)));
    assert(sflint::fmpz_equal_si(moved_from_conductor, 17));

    silex::Order assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(assigned.parent() != nullptr &&
           assigned.parent()->has_same_data(field));
    assert(assigned.has_basis());
    assert(assigned.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));
    assert(!moved.is_defined());

    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(assigned.degree() == 0);
    assert(!assigned.has_basis());

    assert(assigned.define(field));
    assert(assigned.is_defined());
    assert(assigned.parent() != nullptr &&
           assigned.parent()->has_same_data(field));
    assert(!assigned.has_basis());
    assert(assigned.set(source));
    assert(!assigned.is_defined());

    assert(source.define_equation_order(field));
    assert(assigned.set(source));
    assert(assigned.has_basis());
    assert(assigned.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));
    return 0;
}

int test_quadratic_maximal_orders() {
    silex::NumberField field = quadratic_field(5);

    silex::Order equation = silex::test::equation_order(field);
    assert(equation.maximality_known());
    assert(!equation.is_maximal());

    silex::Order maximal(field);
    sflint::Fmpz idx;
    sflint::Fmpz disc;
    assert(maximal.maximal_order(equation));
    assert(order_index(sflint::FmpzRef(idx), equation, maximal));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 5));
    assert(maximal.maximality_known());
    assert(maximal.is_maximal());

    silex::Element alpha(field);
    sflint::FmpqPoly alpha_poly;
    sflint::fmpq_poly_set_coeff_si(alpha_poly, 0, 1);
    sflint::fmpq_poly_set_coeff_si(alpha_poly, 1, 1);
    sflint::fmpq_poly_scalar_div_ui(alpha_poly, alpha_poly, 2);
    assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(alpha_poly)));
    assert(maximal.contains(alpha));

    silex::Order again(field);
    assert(again.maximal_order(maximal));
    assert(order_index(sflint::FmpzRef(idx), maximal, again));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(order_index(sflint::FmpzRef(idx), again, maximal));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(again.is_maximal());

    sflint::Fmpz p;
    sflint::fmpz_set_ui(p, 2);
    silex::Order pmax(field);
    assert(pmax.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
    assert(order_index(sflint::FmpzRef(idx), equation, pmax));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(pmax.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 5));
    assert(pmax.is_maximal());

    sflint::fmpz_set_ui(p, 3);
    assert(pmax.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
    assert(order_index(sflint::FmpzRef(idx), equation, pmax));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(pmax.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 20));
    assert(pmax.maximality_known());
    assert(!pmax.is_maximal());

    silex::Order explicit_max(field);
    sflint::FmpqMat basis(2, 2);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_frac_si(basis, 1, 0, 1, 2);
    mat_entry_frac_si(basis, 1, 1, 1, 2);
    assert(explicit_max.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(explicit_max.maximality_known());
    assert(explicit_max.is_maximal());

    sflint::fmpz_set_ui(p, 2);
    assert(pmax.pmaximal_overorder(explicit_max, sflint::FmpzConstRef(p)));
    assert(order_index(sflint::FmpzRef(idx), explicit_max, pmax));
    assert(sflint::fmpz_equal_si(idx, 1));
    sflint::fmpz_set_ui(p, 5);
    assert(pmax.pmaximal_overorder(explicit_max, sflint::FmpzConstRef(p)));
    assert(order_index(sflint::FmpzRef(idx), explicit_max, pmax));
    assert(sflint::fmpz_equal_si(idx, 1));

    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_equation = silex::test::equation_order(sqrt2);
    assert(sqrt2_equation.is_maximal());
    silex::Order sqrt2_max(sqrt2);
    assert(sqrt2_max.maximal_order(sqrt2_equation));
    assert(order_index(sflint::FmpzRef(idx), sqrt2_equation, sqrt2_max));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(sqrt2_max.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 8));
    assert(sqrt2_max.is_maximal());

    return 0;
}

int test_polynomial_quadratic_and_generic_global_maximal_orders() {
    sflint::FmpqPoly quadratic_polynomial;
    poly_x2_minus(quadratic_polynomial, 5);
    silex::NumberField polynomial_quadratic =
            field_by_polynomial(quadratic_polynomial);
    assert(polynomial_quadratic.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    silex::Order quadratic_order =
            silex::test::equation_order(polynomial_quadratic);

    silex::Order maximal_quadratic(polynomial_quadratic);
    sflint::Fmpz idx;
    sflint::Fmpz disc;
    assert(maximal_quadratic.maximal_order(quadratic_order));
    assert(order_index(sflint::FmpzRef(idx),
                       quadratic_order, maximal_quadratic));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(maximal_quadratic.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 5));
    assert(maximal_quadratic.maximality_known());
    assert(maximal_quadratic.is_maximal());
    sflint::FmpqMat maximal_basis(2, 2);
    assert(maximal_quadratic.get_basis(sflint::FmpqMatRef(maximal_basis)));
    assert(fmpq_entry_is_si(maximal_basis, 0, 0, 1));
    assert(fmpq_entry_is_si(maximal_basis, 0, 1, 0));
    assert(fmpq_entry_is_frac_si(maximal_basis, 1, 0, 1, 2));
    assert(fmpq_entry_is_frac_si(maximal_basis, 1, 1, 1, 2));

    silex::Element alpha(polynomial_quadratic);
    sflint::FmpqPoly alpha_poly;
    sflint::fmpq_poly_set_coeff_si(alpha_poly, 0, 1);
    sflint::fmpq_poly_set_coeff_si(alpha_poly, 1, 1);
    sflint::fmpq_poly_scalar_div_ui(alpha_poly, alpha_poly, 2);
    assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(alpha_poly)));
    assert(maximal_quadratic.contains(alpha));

    sflint::FmpqPoly imaginary_polynomial;
    poly_x2_minus(imaginary_polynomial, -47);
    silex::NumberField imaginary_quadratic =
            field_by_polynomial(imaginary_polynomial);
    assert(imaginary_quadratic.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    silex::Order imaginary_equation =
            silex::test::equation_order(imaginary_quadratic);
    silex::Order imaginary_maximal(imaginary_quadratic);
    assert(imaginary_maximal.maximal_order(imaginary_equation));
    assert(order_index(
            sflint::FmpzRef(idx), imaginary_equation, imaginary_maximal));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(imaginary_maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -47));
    assert(imaginary_maximal.get_basis(sflint::FmpqMatRef(maximal_basis)));
    assert(fmpq_entry_is_si(maximal_basis, 0, 0, 1));
    assert(fmpq_entry_is_si(maximal_basis, 0, 1, 0));
    assert(fmpq_entry_is_frac_si(maximal_basis, 1, 0, 1, 2));
    assert(fmpq_entry_is_frac_si(maximal_basis, 1, 1, 1, 2));

    sflint::FmpqPoly nonsquarefree_polynomial;
    poly_x2_minus(nonsquarefree_polynomial, 12);
    silex::NumberField generic_quadratic =
            field_by_polynomial(nonsquarefree_polynomial);
    assert(generic_quadratic.backend_kind() ==
            silex::NumberFieldBackendKind::generic);
    silex::Order generic_quadratic_order =
            silex::test::equation_order(generic_quadratic);
    silex::Order generic_quadratic_maximal(generic_quadratic);
    assert(generic_quadratic_maximal.maximal_order(generic_quadratic_order));
    assert(order_index(sflint::FmpzRef(idx),
            generic_quadratic_order,
            generic_quadratic_maximal));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(generic_quadratic_maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 12));
    assert(generic_quadratic_maximal.maximality_known());
    assert(generic_quadratic_maximal.is_maximal());

    silex::Element half_theta(generic_quadratic);
    sflint::fmpq_poly_zero(alpha_poly);
    sflint::fmpq_poly_set_coeff_si(alpha_poly, 1, 1);
    sflint::fmpq_poly_scalar_div_ui(alpha_poly, alpha_poly, 2);
    assert(half_theta.set_fmpq_poly(sflint::FmpqPolyConstRef(alpha_poly)));
    assert(generic_quadratic_maximal.contains(half_theta));

    sflint::FmpqPoly cubic_polynomial;
    poly_x3_minus(cubic_polynomial, 2);
    silex::NumberField cubic = field_by_polynomial(cubic_polynomial);
    silex::Order cubic_order = silex::test::equation_order(cubic);

    silex::Order maximal_cubic(cubic);
    assert(maximal_cubic.maximal_order(cubic_order));
    assert(order_index(sflint::FmpzRef(idx), cubic_order, maximal_cubic));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(maximal_cubic.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -108));
    assert(maximal_cubic.maximality_known());
    assert(maximal_cubic.is_maximal());

    return 0;
}

int test_maximal_order_failure_preserves_output() {
    silex::NumberField quadratic = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(quadratic);

    sflint::FmpqPoly generic_polynomial;
    poly_x2_minus(generic_polynomial, 12);
    silex::NumberField generic_field = field_by_polynomial(generic_polynomial);
    silex::Order generic_order = silex::test::equation_order(generic_field);

    silex::Order out = silex::test::equation_order(quadratic);
    assert(out.maximality_known());
    assert(!out.is_maximal());

    sflint::Fmpz p;
    sflint::fmpz_set_ui(p, 1);
    assert(!out.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
    assert(out.maximality_known());
    assert(!out.is_maximal());

    sflint::fmpz_set_ui(p, 2);
    assert(!out.pmaximal_overorder(out, sflint::FmpzConstRef(p)));
    assert(out.maximality_known());
    assert(!out.is_maximal());
    assert(!out.maximal_order(out));
    assert(out.maximality_known());
    assert(!out.is_maximal());

    silex::Order wrong_parent(generic_field);
    assert(!wrong_parent.maximal_order(equation));
    assert(!wrong_parent.has_basis());
    assert(!wrong_parent.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
    assert(!wrong_parent.has_basis());

    silex::Order generic_out = silex::test::equation_order(generic_field);
    assert(!generic_out.maximality_known());
    assert(generic_out.maximal_order(generic_order));
    assert(generic_out.maximality_known());
    assert(generic_out.is_maximal());
    assert(generic_out.pmaximal_overorder(generic_order,
                                          sflint::FmpzConstRef(p)));
    sflint::Fmpz idx;
    sflint::Fmpz disc;
    assert(order_index(sflint::FmpzRef(idx), generic_order, generic_out));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(generic_out.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, 12));
    assert(!generic_out.maximality_known());
    assert(!generic_out.is_maximal());

    return 0;
}

int test_dedekind_criterion_equation_orders() {
    // x^3 - 10 has discriminant -2700 = -2^2 3^3 5^2.  It is Eisenstein at 2
    // and 5, so Z[theta] is 2- and 5-maximal.  At 3, x^3 - 10 = (x - 1)^3 and
    // (f - (x - 1)^3)/3 = x^2 - x - 3 = x(x - 1) mod 3, so Dedekind's gcd is
    // x - 1 and Z[theta] is not 3-maximal: since 10 = 1 mod 9,
    // (1 + theta + theta^2)/3 is integral, [O_K : Z[theta]] = 3 and
    // disc(O_K) = -300.
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 10);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order equation = silex::test::equation_order(field);
    assert(equation.is_equation_order());

    sflint::FmpqPoly third_poly;
    sflint::Fmpq third;
    sflint::fmpq_set_si(third, 1, 3);
    for (slong i = 0; i < 3; ++i) {
        sflint::fmpq_poly_set_coeff_fmpq(third_poly, i, third);
    }
    silex::Element third_element(field);
    assert(third_element.set_fmpq_poly(sflint::FmpqPolyConstRef(third_poly)));
    assert(!equation.contains(third_element));

    sflint::Fmpz p;
    sflint::Fmpz idx;
    sflint::Fmpz disc;
    for (slong prime = 2; prime <= 5; prime += 3) {
        sflint::fmpz_set_si(p, prime);
        silex::Order pmax(field);
        assert(pmax.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
        assert(pmax.is_equation_order());
        assert(!pmax.maximality_known());
        assert(order_index(sflint::FmpzRef(idx), equation, pmax));
        assert(sflint::fmpz_equal_si(idx, 1));
    }

    sflint::fmpz_set_si(p, 3);
    silex::Order pmax3(field);
    assert(pmax3.pmaximal_overorder(equation, sflint::FmpzConstRef(p)));
    assert(!pmax3.is_equation_order());
    assert(pmax3.contains(third_element));
    assert(order_index(sflint::FmpzRef(idx), equation, pmax3));
    assert(sflint::fmpz_equal_si(idx, 3));
    assert(pmax3.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -300));

    // The global route enlarges at 3 and keeps the enlarged order at 5.
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.contains(third_element));
    assert(order_index(sflint::FmpzRef(idx), equation, maximal));
    assert(sflint::fmpz_equal_si(idx, 3));
    assert(order_index(sflint::FmpzRef(idx), pmax3, maximal));
    assert(sflint::fmpz_equal_si(idx, 1));
    assert(maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -300));
    assert(maximal.is_maximal());

    // x^3 - 2 (discriminant -108 = -2^2 3^3) passes Dedekind at 2 and 3, so
    // the maximal order is the input equation order itself.
    poly_x3_minus(polynomial, 2);
    silex::NumberField pure_cubic = field_by_polynomial(polynomial);
    silex::Order pure_equation = silex::test::equation_order(pure_cubic);
    silex::Order pure_maximal(pure_cubic);
    assert(pure_maximal.maximal_order(pure_equation));
    assert(pure_maximal.is_equation_order());
    assert(pure_maximal.is_maximal());
    sflint::fmpz_set_si(p, 3);
    silex::Order pure_pmax(pure_cubic);
    assert(pure_pmax.pmaximal_overorder(pure_equation,
                                        sflint::FmpzConstRef(p)));
    assert(pure_pmax.is_equation_order());

    // Dedekind's x^3 - x^2 - 2x - 8 has discriminant -2012 = -2^2 * 503.
    // Modulo 2 it is x^2 (x + 1) (mixed multiplicities), the criterion fails
    // at 2, [O_K : Z[theta]] = 2 and disc(O_K) = -503.
    sflint::FmpqPoly dedekind_poly;
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 0, -8);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 1, -2);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 2, -1);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 3, 1);
    silex::NumberField dedekind_field = field_by_polynomial(dedekind_poly);
    silex::Order dedekind_equation = silex::test::equation_order(dedekind_field);
    silex::Order dedekind_maximal(dedekind_field);
    assert(dedekind_maximal.maximal_order(dedekind_equation));
    assert(!dedekind_maximal.is_equation_order());
    assert(order_index(sflint::FmpzRef(idx), dedekind_equation,
                       dedekind_maximal));
    assert(sflint::fmpz_equal_si(idx, 2));
    assert(dedekind_maximal.discriminant(sflint::FmpzRef(disc)));
    assert(sflint::fmpz_equal_si(disc, -503));

    return 0;
}

// Prime decomposition reads b_i * b_j at row i * n + j of
// Order::multiplication_table.  Check that table entry by entry against the
// products of the basis elements formed one at a time with OrderElement.
bool table_matches_basis_products(const silex::Order& input) noexcept {
    silex::Order order = input;
    const slong n = order.degree();
    sflint::FmpzMat table(n * n, n);
    if (!order.multiplication_table(sflint::FmpzMatRef(table))) {
        return false;
    }

    silex::OrderElement left(order);
    silex::OrderElement right(order);
    silex::OrderElement product(order);
    sflint::FmpzMat left_coords(1, n);
    sflint::FmpzMat right_coords(1, n);
    sflint::FmpzMat coords(1, n);
    for (slong i = 0; i < n; ++i) {
        sflint::fmpz_mat_zero(sflint::FmpzMatRef(left_coords));
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(left_coords, 0, i), 1);
        if (!left.set_coordinates(sflint::FmpzMatConstRef(left_coords))) {
            return false;
        }
        for (slong j = 0; j < n; ++j) {
            sflint::fmpz_mat_zero(sflint::FmpzMatRef(right_coords));
            sflint::fmpz_set_si(sflint::fmpz_mat_entry(right_coords, 0, j), 1);
            if (!right.set_coordinates(
                        sflint::FmpzMatConstRef(right_coords)) ||
                !product.multiply(left, right) ||
                !product.get_coordinates(sflint::FmpzMatRef(coords))) {
                return false;
            }
            for (slong k = 0; k < n; ++k) {
                if (::fmpz_equal(::fmpz_mat_entry(table.raw(), i * n + j, k),
                                 ::fmpz_mat_entry(coords.raw(), 0, k)) == 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

int test_maximal_order_table_matches_basis_products() {
    // Dedekind's x^3 - x^2 - 2x - 8: O_K is not monogenic and is not the
    // equation order.
    sflint::FmpqPoly dedekind_poly;
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 0, -8);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 1, -2);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 2, -1);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 3, 1);
    silex::NumberField dedekind_field = field_by_polynomial(dedekind_poly);
    silex::Order dedekind_maximal(dedekind_field);
    assert(dedekind_maximal.maximal_order(
            silex::test::equation_order(dedekind_field)));
    assert(!dedekind_maximal.is_equation_order());
    assert(table_matches_basis_products(dedekind_maximal));
    // A second read comes from the memoised table on the shared order data.
    assert(table_matches_basis_products(dedekind_maximal));

    // x^3 - 16 = (2 * 2^(1/3))^3 - 16: O_K = Z[2^(1/3)] is monogenic but is
    // not the equation order of the defining polynomial.
    sflint::FmpqPoly pure_poly;
    poly_x3_minus(pure_poly, 16);
    silex::NumberField pure_field = field_by_polynomial(pure_poly);
    silex::Order pure_maximal(pure_field);
    assert(pure_maximal.maximal_order(
            silex::test::equation_order(pure_field)));
    assert(!pure_maximal.is_equation_order());
    assert(table_matches_basis_products(pure_maximal));

    // Equation order of x^3 - 2 (monogenic, power basis).
    poly_x3_minus(pure_poly, 2);
    silex::NumberField power_field = field_by_polynomial(pure_poly);
    assert(table_matches_basis_products(
            silex::test::equation_order(power_field)));

    // Z + 2 O_K of x^3 - 2: neither maximal nor an equation order.
    silex::Order sub(power_field);
    sflint::FmpqMat sub_basis(3, 3);
    mat_entry_si(sub_basis, 0, 0, 1);
    mat_entry_si(sub_basis, 1, 1, 2);
    mat_entry_si(sub_basis, 2, 2, 2);
    assert(sub.set_basis(sflint::FmpqMatConstRef(sub_basis)));
    assert(!sub.is_equation_order());
    assert(table_matches_basis_products(sub));
    return 0;
}

silex::Order local_equation_order() {
    silex::NumberField field = quadratic_field(2);
    return silex::Order::equation_order(field);
}

int test_order_keeps_parent_alive() {
    silex::Order order = local_equation_order();
    assert(order.is_defined());
    assert(order.parent() != nullptr);
    assert(order.parent()->is_defined());
    assert(order.has_basis());

    silex::Element generator(*order.parent());
    assert(generator.gen());
    assert(order.contains(generator));
    return 0;
}

// Maximality is computed, never asserted.  Order has no accessible
// member that marks or withdraws maximality; each concept must be false.
template <typename T>
concept can_set_maximality = requires(T& order, bool value) {
    order.set_maximality(value);
};
template <typename T>
concept can_record_maximality = requires(T& order, bool value) {
    order.record_maximality(value);
};
template <typename T>
concept can_clear_maximality = requires(T& order) {
    order.clear_maximality();
};
static_assert(!can_set_maximality<silex::Order>);
static_assert(!can_record_maximality<silex::Order>);
static_assert(!can_clear_maximality<silex::Order>);

// Background: a wrong maximal label on Z[sqrt(5)] could let a unit-index
// proof accept 2 + sqrt(5) as fundamental.  This test checks only the
// maximality label: through the public API the equation order stays
// non-maximal in every handle and copy, and the only maximal label comes
// from computing the maximal order, which is strictly larger.
int test_maximality_is_only_computed() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    assert(equation.maximality_known());
    assert(!equation.is_maximal());

    silex::Order shared = equation;
    silex::Order copied;
    assert(copied.set(equation));
    assert(!shared.is_maximal());
    assert(!copied.is_maximal());

    // The same lattice through an explicit, non-canonical basis carries no
    // maximality label.
    sflint::FmpqMat basis(2, 2);
    mat_entry_si(basis, 0, 1, 1);
    mat_entry_si(basis, 1, 0, 1);
    silex::Order permuted =
            silex::Order::from_basis(field, sflint::FmpqMatConstRef(basis));
    assert(permuted.is_defined());
    assert(!permuted.is_maximal());

    silex::Order maximal(field);
    assert(maximal.maximal_order(permuted));
    assert(maximal.is_maximal());
    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));
    assert(!equation.is_maximal());
    assert(!permuted.is_maximal());

    // A quadratic basis in conductor form carries computed maximality.
    sflint::FmpqMat conductor_form(2, 2);
    mat_entry_si(conductor_form, 0, 0, 1);
    mat_entry_frac_si(conductor_form, 1, 0, 1, 2);
    mat_entry_frac_si(conductor_form, 1, 1, 1, 2);
    silex::Order recognised = silex::Order::from_basis(
            field, sflint::FmpqMatConstRef(conductor_form));
    assert(recognised.maximality_known());
    assert(recognised.is_maximal());

    // In degree one only the identity basis carries the flag.
    sflint::FmpqPoly linear;
    poly_x(linear);
    silex::NumberField rational = field_by_polynomial(linear);
    sflint::FmpqMat negated(1, 1);
    mat_entry_si(negated, 0, 0, -1);
    silex::Order negated_order = silex::Order::from_basis(
            rational, sflint::FmpqMatConstRef(negated));
    assert(negated_order.is_defined());
    assert(!negated_order.maximality_known());
    return 0;
}

// Ideal of `order` whose basis rows, in coordinates of `order`, are `rows`
// (row-major, degree x degree).
silex::Ideal ideal_from_rows(const silex::Order& order,
                             const slong* rows) noexcept {
    const slong n = order.degree();
    sflint::FmpzMat matrix(n, n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            ::fmpz_set_si(fmpz_mat_entry(matrix.raw(), i, j), rows[i * n + j]);
        }
    }
    silex::Ideal out(order);
    assert(out.set_hnf(sflint::FmpzMatConstRef(matrix)));
    return out;
}

// f * O as an ideal of O.
silex::Ideal scalar_ideal(const silex::Order& order, slong scalar) noexcept {
    const slong n = order.degree();
    sflint::FmpzMat matrix(n, n);
    for (slong i = 0; i < n; ++i) {
        ::fmpz_set_si(fmpz_mat_entry(matrix.raw(), i, i), scalar);
    }
    silex::Ideal out(order);
    assert(out.set_hnf(sflint::FmpzMatConstRef(matrix)));
    return out;
}

// For Z + f*O_K the conductor in O_K is f*O_K.  Check the general conductor
// against Order::quadratic_conductor and against the contraction of f*O_K.
void check_quadratic_conductor(slong radicand, slong conductor) noexcept {
    silex::NumberField field = quadratic_field(radicand);
    sflint::Fmpz f;
    sflint::fmpz_set_si(sflint::FmpzRef(f), conductor);
    silex::Order order =
            silex::Order::quadratic_order(field, sflint::FmpzConstRef(f));
    assert(order.is_defined());
    sflint::Fmpz recorded;
    assert(order.quadratic_conductor(sflint::FmpzRef(recorded)));
    assert(sflint::fmpz_equal_si(recorded, conductor));

    silex::Order maximal(field);
    assert(maximal.maximal_order(order));

    silex::Ideal general(order);
    assert(silex::detail::order_conductor(general, order, maximal));
    assert(general.parent()->has_same_data(order));

    silex::Ideal expected(order);
    assert(silex::detail::ideal_contract(
            expected, scalar_ideal(maximal, sflint::fmpz_get_si(recorded)),
            order));
    assert(general.equal(expected));

    // [R : f O_K] = [O_K : f O_K] / [O_K : R] = f^2 / f.
    sflint::Fmpz norm;
    assert(general.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, conductor));
}

int test_general_conductor_quadratic() {
    check_quadratic_conductor(5, 2);   // Z[sqrt(5)]: f = 2 O_K
    check_quadratic_conductor(2, 3);   // Z[3 sqrt(2)]: f = 3 O_K
    check_quadratic_conductor(5, 1);
    check_quadratic_conductor(5, 6);
    check_quadratic_conductor(2, 1);
    check_quadratic_conductor(-1, 4);
    check_quadratic_conductor(-3, 2);
    check_quadratic_conductor(-3, 5);
    check_quadratic_conductor(13, 3);
    check_quadratic_conductor(-7, 6);

    // Z[sqrt(5)] is the equation order; f = 2 O_K there as well.
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    silex::Ideal conductor(equation);
    assert(silex::detail::order_conductor(conductor, equation, maximal));
    silex::Ideal expected(equation);
    assert(silex::detail::ideal_contract(expected, scalar_ideal(maximal, 2),
                                         equation));
    assert(conductor.equal(expected));

    // The conductor of an order in itself is the unit ideal, and contraction
    // to the ideal's own order is a copy.
    silex::Ideal self(maximal);
    assert(silex::detail::order_conductor(self, maximal, maximal));
    assert(self.is_one());
    silex::Ideal copied(maximal);
    silex::Ideal two = scalar_ideal(maximal, 2);
    assert(silex::detail::ideal_contract(copied, two, maximal));
    assert(copied.equal(two));

    // Contraction of ideals of O_K that are not inside R = Z[sqrt(5)], so the
    // intersection A ∩ R is a proper sublattice of A.  The maximal order has
    // basis (1, (1+sqrt(5))/2); a+b*sqrt(5) lies in an ideal A of O_K iff
    // the division below stays integral.  By hand (basis (1, sqrt(5)) of R):
    //   O_K            -> R (unit ideal)
    //   sqrt(5) O_K    -> a+b sqrt5 in it iff 5 | a: rows (5,0),(0,1)
    //   2 O_K          -> iff a = b mod 2: rows (1,1),(0,2)
    //   3 O_K          -> iff 3 | a and 3 | b: 3R, rows (3,0),(0,3), norm 9
    sflint::FmpqMat basis(2, 2);
    assert(maximal.get_basis(sflint::FmpqMatRef(basis)));
    assert(fmpq_entry_is_si(basis, 0, 0, 1));
    assert(fmpq_entry_is_si(basis, 0, 1, 0));
    assert(fmpq_entry_is_frac_si(basis, 1, 0, 1, 2));
    assert(fmpq_entry_is_frac_si(basis, 1, 1, 1, 2));
    // sqrt(5) = -1 + 2 phi and sqrt(5) phi = 2 + phi, with phi = (1+sqrt5)/2;
    // the HNF of rows (-1,2),(2,1) is (1,3),(0,5).
    const slong sqrt5_rows[4] = {1, 3, 0, 5};
    silex::Ideal contracted(equation);
    assert(silex::detail::ideal_contract(
            contracted, ideal_from_rows(maximal, sqrt5_rows), equation));
    const slong sqrt5_expected[4] = {5, 0, 0, 1};
    assert(contracted.equal(ideal_from_rows(equation, sqrt5_expected)));
    assert(silex::detail::ideal_contract(
            contracted, scalar_ideal(maximal, 1), equation));
    assert(contracted.is_one());
    const slong two_expected[4] = {1, 1, 0, 2};
    assert(silex::detail::ideal_contract(
            contracted, scalar_ideal(maximal, 2), equation));
    assert(contracted.equal(ideal_from_rows(equation, two_expected)));
    assert(silex::detail::ideal_contract(
            contracted, scalar_ideal(maximal, 3), equation));
    assert(contracted.equal(scalar_ideal(equation, 3)));
    sflint::Fmpz contracted_norm;
    assert(contracted.norm(sflint::FmpzRef(contracted_norm)));
    assert(sflint::fmpz_equal_si(contracted_norm, 9));
    return 0;
}

// Cubic equation orders against the external reference GP 2.17.4.  For
// Z[alpha] the conductor is f'(alpha) times the inverse different, and the
// gp values below were computed with:
//   f(p)={my(nf=nfinit(p),F,B);
//     F=idealdiv(nf,deriv(nf.pol),nf.diff);
//     B=matrix(3,3,i,j,polcoeff(lift(nfbasistoalg(nf,F[,j])),i-1));
//     print(F," ",mathnf(B));}
//   f(x^3-x-216)  -> F = [2,1,0;0,1,0;0,0,2] on nf.zk,
//                    mathnf(B) = [2,1,1;0,1,0;0,0,1]
//   f(x^3-54)     -> F = 9 * identity on nf.zk (F = 9 O_K),
//                    mathnf(B) = [9,0,0;0,3,0;0,0,1]
// mathnf(B) holds the conductor's generators as columns in power-basis
// coordinates, which are the coordinates of the equation order; the rows
// below are its columns.
int test_general_conductor_cubic() {
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, -216);
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));

        silex::Ideal conductor(equation);
        assert(silex::detail::order_conductor(conductor, equation, maximal));
        const slong gp_rows[9] = {2, 0, 0, 1, 1, 0, 1, 0, 1};
        assert(conductor.equal(ideal_from_rows(equation, gp_rows)));
        sflint::Fmpz norm;
        assert(conductor.norm(sflint::FmpzRef(norm)));
        assert(sflint::fmpz_equal_si(norm, 2));
    }
    {
        sflint::FmpqPoly polynomial;
        poly_x3_minus(polynomial, 54);
        silex::NumberField field = field_by_polynomial(polynomial);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));

        silex::Ideal conductor(equation);
        assert(silex::detail::order_conductor(conductor, equation, maximal));
        const slong gp_rows[9] = {9, 0, 0, 0, 3, 0, 0, 0, 1};
        assert(conductor.equal(ideal_from_rows(equation, gp_rows)));

        silex::Ideal contracted(equation);
        assert(silex::detail::ideal_contract(
                contracted, scalar_ideal(maximal, 9), equation));
        assert(conductor.equal(contracted));
    }
    return 0;
}

int test_conductor_and_contraction_failures() {
    silex::NumberField field = quadratic_field(5);
    sflint::Fmpz two;
    sflint::fmpz_set_si(sflint::FmpzRef(two), 2);
    silex::Order order =
            silex::Order::quadratic_order(field, sflint::FmpzConstRef(two));
    silex::Order maximal(field);
    assert(maximal.maximal_order(order));

    // O_K is not contained in Z[sqrt(5)]: both routines fail and leave the
    // output unchanged.
    silex::Ideal out = scalar_ideal(maximal, 3);
    silex::Ideal before = scalar_ideal(maximal, 3);
    assert(!silex::detail::order_conductor(out, maximal, order));
    assert(out.equal(before));
    silex::Ideal two_order = scalar_ideal(order, 2);
    assert(!silex::detail::ideal_contract(out, two_order, maximal));
    assert(out.equal(before));

    // The output must already be defined over the order it is an ideal of.
    silex::Ideal wrong_parent = scalar_ideal(maximal, 3);
    assert(!silex::detail::order_conductor(wrong_parent, order, maximal));
    assert(wrong_parent.equal(before));
    silex::Ideal undefined;
    assert(!silex::detail::order_conductor(undefined, order, maximal));
    assert(!undefined.is_defined());

    // Orders of different fields.
    silex::NumberField other = quadratic_field(2);
    silex::Order other_order = silex::test::equation_order(other);
    silex::Ideal other_out(other_order);
    assert(!silex::detail::order_conductor(other_out, other_order, maximal));
    silex::Ideal other_before = scalar_ideal(other_order, 3);
    other_out = scalar_ideal(other_order, 3);
    assert(!silex::detail::order_conductor(other_out, other_order, maximal));
    assert(other_out.equal(other_before));
    assert(!silex::detail::ideal_contract(other_out, before, other_order));
    assert(other_out.equal(other_before));

    // ideal_contract: output over the wrong order, or not defined.
    silex::Ideal contract_wrong = scalar_ideal(maximal, 3);
    assert(!silex::detail::ideal_contract(contract_wrong, two_order, order));
    assert(contract_wrong.equal(before));
    silex::Ideal contract_undefined;
    assert(!silex::detail::ideal_contract(contract_undefined, two_order,
                                          order));
    assert(!contract_undefined.is_defined());
    return 0;
}

// The cached index [O : Z[theta]] equals a fresh order_index computation, on
// a cache hit, through handle copies and Order::set, and a new basis on the
// same handle (or on a copy) is not answered from the old order's cache.
bool cached_index_is(const silex::Order& order, slong expected) noexcept {
    sflint::Fmpz cached;
    sflint::Fmpz fresh;
    silex::Order equation = silex::Order::equation_order(*order.parent());
    return silex::detail::OrderAccess::equation_order_index(
                   sflint::FmpzRef(cached), order) &&
           order_index(sflint::FmpzRef(fresh), equation, order) &&
           fmpz_equal(cached.raw(), fresh.raw()) != 0 &&
           sflint::fmpz_equal_si(cached, expected);
}

int test_cached_equation_order_index() {
    // Q(cbrt 54) = Q(cbrt 2): disc(x^3 - 54) = -78732 = -108 * 27^2, so the
    // maximal order has index 27 over Z[cbrt 54].
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 54);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    assert(cached_index_is(equation, 1));
    assert(cached_index_is(maximal, 27));
    assert(cached_index_is(maximal, 27));  // cache hit

    silex::Order shared = maximal;  // same data, same cache
    assert(shared.has_same_data(maximal));
    assert(cached_index_is(shared, 27));

    silex::Order copied;
    assert(copied.set(maximal));  // deep copy carries the cached value
    assert(!copied.has_same_data(maximal));
    assert(cached_index_is(copied, 27));

    // A new basis on a handle whose order has a cached index.
    sflint::FmpqMat identity(3, 3);
    fmpq_mat_one(identity.raw());
    assert(copied.set_basis(sflint::FmpqMatConstRef(identity)));
    assert(cached_index_is(copied, 1));
    assert(cached_index_is(maximal, 27));

    sflint::FmpqMat maximal_basis(3, 3);
    assert(maximal.get_basis(sflint::FmpqMatRef(maximal_basis)));
    silex::Order rebuilt = silex::Order::from_basis(
            field, sflint::FmpqMatConstRef(maximal_basis));
    assert(rebuilt.is_defined());
    assert(cached_index_is(rebuilt, 27));
    assert(shared.set_basis(sflint::FmpqMatConstRef(identity)));
    assert(cached_index_is(shared, 1));
    assert(cached_index_is(maximal, 27));

    // An order that does not contain Z[theta] has no index; nothing is
    // cached, and a repeated call fails the same way.
    silex::NumberField quadratic = quadratic_field(5);
    sflint::Fmpz four;
    sflint::fmpz_set_si(sflint::FmpzRef(four), 4);
    silex::Order small =
            silex::Order::quadratic_order(quadratic, sflint::FmpzConstRef(four));
    assert(small.is_defined());
    sflint::Fmpz index;
    assert(!silex::detail::OrderAccess::equation_order_index(
            sflint::FmpzRef(index), small));
    assert(!silex::detail::OrderAccess::equation_order_index(
            sflint::FmpzRef(index), small));

    silex::Order undefined;
    assert(!silex::detail::OrderAccess::equation_order_index(
            sflint::FmpzRef(index), undefined));
    return 0;
}

// The noninstalled view of the cached multiplication table: the same entries
// the public copy returns, the same storage on a later call and through a
// handle sharing the order data, a separate table for a deep copy, and no
// table for an order without a basis.
int test_cached_multiplication_table_view() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order suborder(field);
    sflint::FmpqMat basis(2, 2);
    mat_entry_si(basis, 0, 0, 1);
    mat_entry_si(basis, 1, 1, 2);
    assert(suborder.set_basis(sflint::FmpqMatConstRef(basis)));

    // View first (fills the cache), then the public copy.
    const sflint::FmpzMat* view =
            silex::detail::OrderAccess::multiplication_table(suborder);
    assert(view != nullptr);
    assert(fmpz_mat_nrows(view->raw()) == 4 && fmpz_mat_ncols(view->raw()) == 2);
    assert(fmpz_mat_entry_is_si(*view, 3, 0, 8));
    sflint::FmpzMat table(4, 2);
    assert(suborder.multiplication_table(sflint::FmpzMatRef(table)));
    assert(fmpz_mat_equal(table.raw(), view->raw()) != 0);
    assert(silex::detail::OrderAccess::multiplication_table(suborder) == view);

    silex::Order shared = suborder;
    assert(silex::detail::OrderAccess::multiplication_table(shared) == view);

    silex::Order copied;
    assert(copied.set(suborder));
    const sflint::FmpzMat* copied_view =
            silex::detail::OrderAccess::multiplication_table(copied);
    assert(copied_view != nullptr && copied_view != view);
    assert(fmpz_mat_equal(copied_view->raw(), view->raw()) != 0);

    // Public copy first, then the view, on a fresh order.
    silex::Order equation = silex::test::equation_order(field);
    assert(equation.multiplication_table(sflint::FmpzMatRef(table)));
    const sflint::FmpzMat* equation_view =
            silex::detail::OrderAccess::multiplication_table(equation);
    assert(equation_view != nullptr);
    assert(fmpz_mat_equal(table.raw(), equation_view->raw()) != 0);
    assert(fmpz_mat_entry_is_si(*equation_view, 3, 0, 2));

    silex::Order undefined;
    assert(silex::detail::OrderAccess::multiplication_table(undefined) ==
           nullptr);
    silex::Order no_basis(field);
    assert(silex::detail::OrderAccess::multiplication_table(no_basis) ==
           nullptr);
    return 0;
}


// Reference power in O/pO: the product reduces every table entry and every
// partial sum modulo p, as Order::p_radical did before its products moved to
// a table reduced once.  Same right-to-left binary powering.
void reference_coords_pow_mod(sflint::FmpzVec& out,
                              const sflint::FmpzVec& base,
                              const fmpz_t exponent,
                              const sflint::FmpzVec& one,
                              const sflint::FmpzMat& table,
                              const sflint::FmpzModCtx& ctx) noexcept {
    const slong n = out.length();
    sflint::FmpzVec acc(n);
    sflint::FmpzVec b(n);
    sflint::FmpzVec product(n);
    sflint::Fmpz e;
    sflint::Fmpz ab;
    sflint::Fmpz coeff;
    sflint::Fmpz term;
    _fmpz_vec_set(acc.data(), one.data(), n);
    _fmpz_vec_set(b.data(), base.data(), n);
    fmpz_set(e.raw(), exponent);
    auto mul = [&](sflint::FmpzVec& c, const sflint::FmpzVec& x,
                   const sflint::FmpzVec& y) {
        _fmpz_vec_zero(product.data(), n);
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < n; ++j) {
                fmpz_mod_mul(ab.raw(), x.data() + i, y.data() + j, ctx.raw());
                for (slong k = 0; k < n; ++k) {
                    fmpz_mod_set_fmpz(coeff.raw(),
                                      fmpz_mat_entry(table.raw(), i * n + j, k),
                                      ctx.raw());
                    fmpz_mod_mul(term.raw(), ab.raw(), coeff.raw(), ctx.raw());
                    fmpz_mod_add(product.data() + k, product.data() + k,
                                 term.raw(), ctx.raw());
                }
            }
        }
        _fmpz_vec_set(c.data(), product.data(), n);
    };
    while (fmpz_is_zero(e.raw()) == 0) {
        if (fmpz_is_odd(e.raw()) != 0) {
            mul(acc, acc, b);
        }
        fmpz_fdiv_q_2exp(e.raw(), e.raw(), 1);
        if (fmpz_is_zero(e.raw()) == 0) {
            mul(b, b, b);
        }
    }
    _fmpz_vec_set(out.data(), acc.data(), n);
}

// Checks the word and fmpz product paths of the p_radical power against the
// reference for one modulus, on the given table and random bases with some
// zero coordinates.
void check_coords_pow_paths(const sflint::FmpzMat& table,
                            const sflint::FmpzVec& one_coords,
                            const fmpz_t p,
                            flint_rand_t state) noexcept {
    const slong n = fmpz_mat_ncols(table.raw());
    sflint::FmpzModCtx ctx(p);
    sflint::FmpzVec one(n);
    sflint::FmpzVec base(n);
    sflint::FmpzVec word_out(n);
    sflint::FmpzVec fmpz_out(n);
    sflint::FmpzVec expected(n);
    for (slong i = 0; i < n; ++i) {
        fmpz_mod(one.data() + i, one_coords.data() + i, p);
    }

    sflint::Fmpz exponents[5];
    fmpz_zero(exponents[0].raw());
    fmpz_one(exponents[1].raw());
    fmpz_set_ui(exponents[2].raw(), 2);
    // The p_radical exponent: the least power of p above n.
    fmpz_set(exponents[3].raw(), p);
    while (fmpz_cmp_si(exponents[3].raw(), n) <= 0) {
        fmpz_mul(exponents[3].raw(), exponents[3].raw(), p);
    }
    fmpz_one(exponents[4].raw());
    fmpz_mul_2exp(exponents[4].raw(), exponents[4].raw(), 70);
    fmpz_add_ui(exponents[4].raw(), exponents[4].raw(), 12345);

    for (int trial = 0; trial < 6; ++trial) {
        for (slong i = 0; i < n; ++i) {
            if (trial > 0 && n_randint(state, 4) == 0) {
                fmpz_zero(base.data() + i);
            } else if (trial == 0) {
                fmpz_sub_ui(base.data() + i, p, 1);
            } else {
                fmpz_randm(base.data() + i, state, p);
            }
        }
        for (const sflint::Fmpz& e : exponents) {
            reference_coords_pow_mod(expected, base, e.raw(), one, table, ctx);
            assert(silex::detail::order_coords_pow_mod_p(
                    word_out, base, sflint::FmpzConstRef(e), one, table,
                    sflint::FmpzConstRef(p), true));
            assert(silex::detail::order_coords_pow_mod_p(
                    fmpz_out, base, sflint::FmpzConstRef(e), one, table,
                    sflint::FmpzConstRef(p), false));
            assert(_fmpz_vec_equal(word_out.data(), expected.data(), n) != 0);
            assert(_fmpz_vec_equal(fmpz_out.data(), expected.data(), n) != 0);
        }
    }
}

int test_p_radical_power_word_and_fmpz_paths() {
    flint_rand_t state;
    flint_rand_init(state);

    // Dedekind's x^3 - x^2 - 2x - 8: a non-equation maximal order, whose
    // multiplication table p_radical uses.
    sflint::FmpqPoly dedekind_poly;
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 0, -8);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 1, -2);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 2, -1);
    sflint::fmpq_poly_set_coeff_si(dedekind_poly, 3, 1);
    silex::NumberField field = field_by_polynomial(dedekind_poly);
    silex::Order maximal(field);
    assert(maximal.maximal_order(silex::test::equation_order(field)));
    assert(!maximal.is_equation_order());
    sflint::FmpzMat order_table(9, 3);
    assert(maximal.multiplication_table(sflint::FmpzMatRef(order_table)));
    sflint::FmpzVec order_one(3);
    {
        sflint::FmpqMat basis(3, 3);
        sflint::FmpqMat inverse(3, 3);
        assert(maximal.get_basis(sflint::FmpqMatRef(basis)));
        assert(fmpq_mat_inv(inverse.raw(), basis.raw()) != 0);
        // 1 is the first row of the inverse basis (power-basis vector e_0).
        for (slong j = 0; j < 3; ++j) {
            assert(fmpz_is_one(fmpq_mat_entry_den(inverse.raw(), 0, j)) != 0);
            fmpz_set(order_one.data() + j,
                     fmpq_mat_entry_num(inverse.raw(), 0, j));
        }
    }

    // A synthetic 4 x 4 table with signed entries of up to 200 bits, so the
    // table reduction and the unreduced accumulation both see entries far
    // outside [0, p).
    const slong n = 4;
    sflint::FmpzMat wide_table(n * n, n);
    for (slong row = 0; row < n * n; ++row) {
        for (slong k = 0; k < n; ++k) {
            fmpz_randtest(fmpz_mat_entry(wide_table.raw(), row, k), state,
                          200);
        }
    }
    sflint::FmpzVec wide_one(n);
    for (slong i = 0; i < n; ++i) {
        fmpz_randtest(wide_one.data() + i, state, 80);
    }

    sflint::Fmpz p;
    const ulong small_primes[] = {2, 3, 5, 101, 503};
    for (ulong small : small_primes) {
        fmpz_set_ui(p.raw(), small);
        check_coords_pow_paths(order_table, order_one, p.raw(), state);
        check_coords_pow_paths(wide_table, wide_one, p.raw(), state);
    }
    // 2^61 - 1, below a full word; 2^64 - 59, the largest prime below 2^64
    // (nmod norm 0); and the least prime above 2^64, where the word request
    // falls back to fmpz residues.
    fmpz_one(p.raw());
    fmpz_mul_2exp(p.raw(), p.raw(), 61);
    fmpz_sub_ui(p.raw(), p.raw(), 1);
    check_coords_pow_paths(order_table, order_one, p.raw(), state);
    check_coords_pow_paths(wide_table, wide_one, p.raw(), state);
    fmpz_set_ui(p.raw(), UWORD_MAX);
    fmpz_sub_ui(p.raw(), p.raw(), 58);
    assert(fmpz_is_prime(p.raw()) == 1);
    check_coords_pow_paths(order_table, order_one, p.raw(), state);
    check_coords_pow_paths(wide_table, wide_one, p.raw(), state);
    fmpz_set_ui(p.raw(), UWORD_MAX);
    fmpz_nextprime(p.raw(), p.raw(), 1);
    assert(fmpz_abs_fits_ui(p.raw()) == 0);
    check_coords_pow_paths(order_table, order_one, p.raw(), state);
    check_coords_pow_paths(wide_table, wide_one, p.raw(), state);

    // Rejected inputs leave the output unchanged.
    sflint::FmpzVec out(3);
    sflint::FmpzVec base(3);
    fmpz_set_ui(out.data(), 7);
    fmpz_set_ui(p.raw(), 5);
    sflint::Fmpz e;
    fmpz_set_ui(e.raw(), 3);
    fmpz_set_ui(base.data(), 5);
    assert(!silex::detail::order_coords_pow_mod_p(
            out, base, sflint::FmpzConstRef(e), order_one, order_table,
            sflint::FmpzConstRef(p), true));
    fmpz_zero(base.data());
    fmpz_set_si(e.raw(), -1);
    assert(!silex::detail::order_coords_pow_mod_p(
            out, base, sflint::FmpzConstRef(e), order_one, order_table,
            sflint::FmpzConstRef(p), true));
    fmpz_set_ui(e.raw(), 3);
    assert(!silex::detail::order_coords_pow_mod_p(
            out, base, sflint::FmpzConstRef(e), order_one, wide_table,
            sflint::FmpzConstRef(p), true));
    // Matching vector lengths (3) but a table of the wrong row count.
    sflint::FmpzMat short_table(6, 3);
    assert(!silex::detail::order_coords_pow_mod_p(
            out, base, sflint::FmpzConstRef(e), order_one, short_table,
            sflint::FmpzConstRef(p), true));
    assert(fmpz_equal_ui(out.data(), 7) != 0);

    // p_radical of the non-equation order at word-size and multiword
    // primes: 503 ramifies (503 = P^2 Q, radical PQ of norm 503^2), and the
    // primes near 2^64 are unramified (radical pO, norm p^3).
    silex::Ideal radical(maximal);
    sflint::Fmpz norm;
    sflint::Fmpz expected_norm;
    fmpz_set_ui(p.raw(), 503);
    assert(maximal.p_radical(radical, sflint::FmpzConstRef(p)));
    assert(radical.norm(sflint::FmpzRef(norm)));
    assert(fmpz_equal_ui(norm.raw(), 503 * 503) != 0);
    fmpz_set_ui(p.raw(), UWORD_MAX);
    fmpz_sub_ui(p.raw(), p.raw(), 58);
    for (int round = 0; round < 2; ++round) {
        assert(maximal.p_radical(radical, sflint::FmpzConstRef(p)));
        assert(radical.norm(sflint::FmpzRef(norm)));
        fmpz_pow_ui(expected_norm.raw(), p.raw(), 3);
        assert(fmpz_equal(norm.raw(), expected_norm.raw()) != 0);
        fmpz_set_ui(p.raw(), UWORD_MAX);
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    flint_rand_clear(state);
    return 0;
}
}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_quadratic_equation_order() == 0);
    assert(test_cubic_equation_order() == 0);
    assert(test_explicit_basis_index_and_table() == 0);
    assert(test_undefined_order_quadratic_conductor() == 0);
    assert(test_failure_preserves_order() == 0);
    assert(test_copy_and_swap() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_quadratic_maximal_orders() == 0);
    assert(test_polynomial_quadratic_and_generic_global_maximal_orders() == 0);
    assert(test_maximal_order_failure_preserves_output() == 0);
    assert(test_dedekind_criterion_equation_orders() == 0);
    assert(test_order_keeps_parent_alive() == 0);
    assert(test_maximality_is_only_computed() == 0);
    assert(test_general_conductor_quadratic() == 0);
    assert(test_general_conductor_cubic() == 0);
    assert(test_conductor_and_contraction_failures() == 0);
    assert(test_cached_equation_order_index() == 0);
    assert(test_maximal_order_table_matches_basis_products() == 0);
    assert(test_cached_multiplication_table_view() == 0);
    assert(test_p_radical_power_word_and_fmpz_paths() == 0);
    return 0;
}
