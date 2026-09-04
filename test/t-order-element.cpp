#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_element.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

bool same_parent(const silex::Order* parent,
                 const silex::Order& order) noexcept {
    return parent != nullptr && parent->has_same_data(order);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

void fmpq_mat_entry_si(sflint::FmpqMat& matrix,
                       slong row,
                       slong col,
                       slong value) noexcept {
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(matrix, row, col), value, 1);
}

void fmpz_mat_entry_si(sflint::FmpzMat& matrix,
                       slong row,
                       slong col,
                       slong value) noexcept {
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool fmpz_entry_is_si(const sflint::FmpzMat& matrix,
                      slong row,
                      slong col,
                      slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

void element_from_coeff(silex::Element& element,
                        slong index,
                        slong value) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, index, value);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

int test_basic_nonmaximal_coordinates() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -1);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    fmpq_mat_entry_si(basis, 0, 0, 1);
    fmpq_mat_entry_si(basis, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::OrderElement a(order);
    silex::OrderElement b(order);
    assert(a.is_defined());
    assert(same_parent(a.parent(), order));
    assert(a.degree() == 2);
    assert(a.equal_si(0));

    assert(a.one());
    assert(a.equal_si(1));

    silex::Element theta(field);
    assert(theta.gen());
    assert(!a.set_element(theta));
    assert(a.equal_si(1));

    silex::Element two_theta(field);
    element_from_coeff(two_theta, 1, 2);
    assert(a.set_element(two_theta));

    sflint::FmpzMat row(1, 2);
    assert(a.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 0));
    assert(fmpz_entry_is_si(row, 0, 1, 1));
    auto owned_coords = a.coordinates();
    assert(owned_coords.has_value());
    assert(fmpz_entry_is_si(*owned_coords, 0, 0, 0));
    assert(fmpz_entry_is_si(*owned_coords, 0, 1, 1));

    silex::Element out(field);
    assert(a.get_element(out));
    assert(out.equal(two_theta));

    fmpz_mat_entry_si(row, 0, 0, 1);
    fmpz_mat_entry_si(row, 0, 1, 1);
    assert(a.set_coordinates(sflint::FmpzMatConstRef(row)));
    assert(a.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 1));
    assert(fmpz_entry_is_si(row, 0, 1, 1));

    fmpz_mat_entry_si(row, 0, 0, 1);
    fmpz_mat_entry_si(row, 0, 1, -1);
    assert(b.set_coordinates(sflint::FmpzMatConstRef(row)));
    assert(!a.equal(b));

    assert(b.set(a));
    assert(a.equal(b));
    return 0;
}

int test_arithmetic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, -1);

    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order(field);
    sflint::FmpqMat basis(2, 2);
    fmpq_mat_entry_si(basis, 0, 0, 1);
    fmpq_mat_entry_si(basis, 1, 1, 2);
    assert(order.set_basis(sflint::FmpqMatConstRef(basis)));

    silex::OrderElement a(order);
    silex::OrderElement b(order);
    silex::OrderElement c(order);
    sflint::FmpzMat row(1, 2);

    fmpz_mat_entry_si(row, 0, 0, 1);
    fmpz_mat_entry_si(row, 0, 1, 1);
    assert(a.set_coordinates(sflint::FmpzMatConstRef(row)));
    fmpz_mat_entry_si(row, 0, 0, 1);
    fmpz_mat_entry_si(row, 0, 1, -1);
    assert(b.set_coordinates(sflint::FmpzMatConstRef(row)));

    assert(c.add(a, b));
    assert(c.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 2));
    assert(fmpz_entry_is_si(row, 0, 1, 0));

    assert(c.subtract(a, b));
    assert(c.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 0));
    assert(fmpz_entry_is_si(row, 0, 1, 2));

    assert(c.negate(a));
    assert(c.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, -1));
    assert(fmpz_entry_is_si(row, 0, 1, -1));

    assert(c.multiply(a, b));
    assert(c.equal_si(5));
    silex::Element out(field);
    assert(c.get_element(out));
    assert(out.equal_si(5));
    return 0;
}

int test_integer_and_failure_preservation() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::OrderElement value(order);
    assert(value.set_si(7));
    assert(value.equal_si(7));

    sflint::FmpzMat row(1, 2);
    assert(value.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 7));
    assert(fmpz_entry_is_si(row, 0, 1, 0));

    sflint::FmpzMat wrong_shape(2, 1);
    assert(!value.set_coordinates(sflint::FmpzMatConstRef(wrong_shape)));
    assert(value.get_coordinates(sflint::FmpzMatRef(row)));
    assert(fmpz_entry_is_si(row, 0, 0, 7));
    assert(fmpz_entry_is_si(row, 0, 1, 0));

    sflint::FmpqPoly other_polynomial;
    poly_x2_minus(other_polynomial, 3);
    silex::NumberField other_field = field_by_polynomial(other_polynomial);
    silex::Element other_element(other_field);
    assert(other_element.set_si(11));
    assert(!value.set_element(other_element));
    assert(value.equal_si(7));

    silex::Element wrong_output(other_field);
    assert(!value.get_element(wrong_output));
    assert(value.equal_si(7));
    return 0;
}

int test_arithmetic_parent_mismatch_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    sflint::FmpqPoly other_polynomial;
    poly_x2_minus(other_polynomial, 3);
    silex::NumberField other_field = field_by_polynomial(other_polynomial);
    silex::Order other_order = silex::test::equation_order(other_field);

    silex::OrderElement left(order);
    silex::OrderElement out(order);
    silex::OrderElement other(other_order);
    assert(left.set_si(2));
    assert(out.set_si(99));
    assert(other.set_si(5));

    assert(!out.add(left, other));
    assert(out.equal_si(99));
    assert(!out.subtract(left, other));
    assert(out.equal_si(99));
    assert(!out.multiply(left, other));
    assert(out.equal_si(99));
    assert(!out.negate(other));
    assert(out.equal_si(99));
    return 0;
}

int test_undefined_order_rejected() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order empty(field);
    assert(!empty.has_basis());

    silex::OrderElement value;
    assert(!value.define(empty));
    assert(!value.is_defined());
    return 0;
}

int test_swap() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order = silex::test::equation_order(field);

    silex::OrderElement left(order);
    silex::OrderElement right(order);
    assert(left.set_si(3));
    assert(right.set_si(9));

    swap(left, right);
    assert(left.equal_si(9));
    assert(right.equal_si(3));
    return 0;
}

int test_move_clear_and_redefine() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);

    silex::Order order = silex::test::equation_order(field);

    silex::OrderElement source(order);
    assert(source.set_si(7));

    silex::OrderElement moved(std::move(source));
    assert(moved.is_defined());
    assert(same_parent(moved.parent(), order));
    assert(moved.degree() == 2);
    assert(moved.equal_si(7));
    assert(!source.is_defined());
    assert(source.parent() == nullptr);
    assert(source.degree() == 0);

    silex::OrderElement assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(same_parent(assigned.parent(), order));
    assert(assigned.equal_si(7));
    assert(!moved.is_defined());
    assert(moved.parent() == nullptr);
    assert(moved.degree() == 0);

    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(assigned.degree() == 0);
    assert(!assigned.set_si(1));

    assert(assigned.define(order));
    assert(assigned.is_defined());
    assert(same_parent(assigned.parent(), order));
    assert(assigned.equal_si(0));
    assert(assigned.set_si(11));
    assert(assigned.equal_si(11));
    return 0;
}

silex::OrderElement local_order_element() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::OrderElement value(order);
    assert(value.set_si(13));
    return value;
}

int test_keeps_parent_order_alive() {
    silex::OrderElement value = local_order_element();
    assert(value.is_defined());
    assert(value.parent() != nullptr);
    assert(value.parent()->parent() != nullptr);
    assert(value.degree() == 2);
    assert(value.equal_si(13));

    silex::Element lifted(*value.parent()->parent());
    assert(value.get_element(lifted));
    assert(lifted.equal_si(13));
    return 0;
}

}  // namespace

int main() {
    assert(test_basic_nonmaximal_coordinates() == 0);
    assert(test_arithmetic() == 0);
    assert(test_integer_and_failure_preservation() == 0);
    assert(test_arithmetic_parent_mismatch_preserves_output() == 0);
    assert(test_undefined_order_rejected() == 0);
    assert(test_swap() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    return 0;
}
