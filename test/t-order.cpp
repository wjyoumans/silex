#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order.hpp>

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
    suborder.set_maximality(false);
    assert(suborder.maximality_known());
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
    out.set_maximality(false);

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
    generic_out.set_maximality(false);
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

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_quadratic_equation_order() == 0);
    assert(test_cubic_equation_order() == 0);
    assert(test_explicit_basis_index_and_table() == 0);
    assert(test_failure_preserves_order() == 0);
    assert(test_copy_and_swap() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_quadratic_maximal_orders() == 0);
    assert(test_polynomial_quadratic_and_generic_global_maximal_orders() == 0);
    assert(test_maximal_order_failure_preserves_output() == 0);
    assert(test_order_keeps_parent_alive() == 0);
    return 0;
}
