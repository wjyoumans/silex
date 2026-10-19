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

// Checks set_coordinates against the element `coeffs` (numerators over
// `den`) in the power basis, and that get_coordinates returns the input row.
void assert_set_coordinates_gives(const silex::Order& order,
                                  const sflint::FmpzMat& row,
                                  const slong* numerators,
                                  ulong den) noexcept {
    const slong degree = order.degree();
    silex::OrderElement value(order);
    assert(value.set_coordinates(sflint::FmpzMatConstRef(row)));

    sflint::FmpqPoly expected_polynomial;
    for (slong j = 0; j < degree; ++j) {
        sflint::Fmpq coefficient;
        sflint::fmpq_set_si(coefficient, numerators[j], den);
        sflint::fmpq_poly_set_coeff_fmpq(expected_polynomial, j,
                                         sflint::FmpqConstRef(coefficient));
    }
    silex::Element expected(*order.parent());
    assert(expected.set_fmpq_poly(
            sflint::FmpqPolyConstRef(expected_polynomial)));
    silex::Element actual(*order.parent());
    assert(value.get_element(actual));
    assert(actual.equal(expected));

    sflint::FmpzMat back(1, degree);
    assert(value.get_coordinates(sflint::FmpzMatRef(back)));
    const sflint::FmpzMat& back_view = back;
    for (slong j = 0; j < degree; ++j) {
        assert(sflint::fmpz_equal(sflint::fmpz_mat_entry(back_view, 0, j),
                                  sflint::fmpz_mat_entry(row, 0, j)));
    }
}

// set_coordinates on an identity basis (equation order, built by the factory
// or through set_basis) reads the coordinates as power-basis coefficients;
// a basis with denominators, including one whose numerators alone form the
// identity, goes through the basis multiply.
int test_set_coordinates_identity_and_general_basis() {
    sflint::FmpqPoly cubic;
    sflint::fmpq_poly_set_coeff_si(cubic, 3, 1);
    sflint::fmpq_poly_set_coeff_si(cubic, 1, 1);
    sflint::fmpq_poly_set_coeff_si(cubic, 0, -3);
    silex::NumberField cubic_field = field_by_polynomial(cubic);

    sflint::FmpzMat cubic_row(1, 3);
    fmpz_mat_entry_si(cubic_row, 0, 0, 3);
    fmpz_mat_entry_si(cubic_row, 0, 1, -1);
    fmpz_mat_entry_si(cubic_row, 0, 2, 2);
    const slong cubic_expected[] = {3, -1, 2};

    silex::Order factory = silex::test::equation_order(cubic_field);
    assert(factory.is_equation_order());
    assert_set_coordinates_gives(factory, cubic_row, cubic_expected, 1);

    silex::Order explicit_identity(cubic_field);
    sflint::FmpqMat identity(3, 3);
    for (slong i = 0; i < 3; ++i) {
        fmpq_mat_entry_si(identity, i, i, 1);
    }
    assert(explicit_identity.set_basis(sflint::FmpqMatConstRef(identity)));
    assert(explicit_identity.is_equation_order());
    assert_set_coordinates_gives(explicit_identity, cubic_row, cubic_expected,
                                 1);

    // x^2 - 8: theta / 2 is integral, so diag(1, 1/2) is an order basis whose
    // numerators form the identity matrix.
    sflint::FmpqPoly quadratic;
    poly_x2_minus(quadratic, 8);
    silex::NumberField quadratic_field = field_by_polynomial(quadratic);
    silex::Order half(quadratic_field);
    sflint::FmpqMat half_basis(2, 2);
    fmpq_mat_entry_si(half_basis, 0, 0, 1);
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(half_basis, 1, 1), 1, 2);
    assert(half.set_basis(sflint::FmpqMatConstRef(half_basis)));
    assert(!half.is_equation_order());
    sflint::FmpzMat quadratic_row(1, 2);
    fmpz_mat_entry_si(quadratic_row, 0, 0, 3);
    fmpz_mat_entry_si(quadratic_row, 0, 1, 4);
    const slong half_expected[] = {3, 2};
    assert_set_coordinates_gives(half, quadratic_row, half_expected, 1);

    // Order::set copies orders whose identity test has already run.
    silex::Order factory_copy(cubic_field);
    assert(factory_copy.set(factory));
    assert_set_coordinates_gives(factory_copy, cubic_row, cubic_expected, 1);
    silex::Order half_copy(quadratic_field);
    assert(half_copy.set(half));
    assert_set_coordinates_gives(half_copy, quadratic_row, half_expected, 1);

    // x^2 - 5 with basis 1, (1 + theta) / 2: off-diagonal denominators.
    sflint::FmpqPoly golden;
    poly_x2_minus(golden, 5);
    silex::NumberField golden_field = field_by_polynomial(golden);
    silex::Order maximal(golden_field);
    sflint::FmpqMat golden_basis(2, 2);
    fmpq_mat_entry_si(golden_basis, 0, 0, 1);
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(golden_basis, 1, 0), 1, 2);
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(golden_basis, 1, 1), 1, 2);
    assert(maximal.set_basis(sflint::FmpqMatConstRef(golden_basis)));
    fmpz_mat_entry_si(quadratic_row, 0, 0, 1);
    fmpz_mat_entry_si(quadratic_row, 0, 1, 1);
    const slong golden_expected[] = {3, 1};
    assert_set_coordinates_gives(maximal, quadratic_row, golden_expected, 2);
    return 0;
}

// is_equation_order must agree with the uncached identity test on the
// order's own basis (false when no basis is set).
void assert_equation_order_matches_uncached(const silex::Order& order) {
    bool expected = false;
    if (order.has_basis()) {
        sflint::FmpqMat basis(order.degree(), order.degree());
        assert(order.get_basis(sflint::FmpqMatRef(basis)));
        expected = fmpq_mat_is_one(basis.raw()) != 0;
    }
    assert(order.is_equation_order() == expected);
}

int test_is_equation_order_matches_uncached() {
    sflint::FmpqPoly cubic;
    sflint::fmpq_poly_set_coeff_si(cubic, 3, 1);
    sflint::fmpq_poly_set_coeff_si(cubic, 1, 1);
    sflint::fmpq_poly_set_coeff_si(cubic, 0, -3);
    silex::NumberField cubic_field = field_by_polynomial(cubic);

    silex::Order bare(cubic_field);
    assert_equation_order_matches_uncached(bare);

    silex::Order factory = silex::test::equation_order(cubic_field);
    // A copy made before the flag is read on the source, then after.
    silex::Order early(cubic_field);
    assert(early.set(factory));
    assert_equation_order_matches_uncached(factory);
    silex::Order late(cubic_field);
    assert(late.set(factory));
    assert_equation_order_matches_uncached(early);
    assert_equation_order_matches_uncached(late);

    silex::Order identity_order(cubic_field);
    sflint::FmpqMat identity(3, 3);
    for (slong i = 0; i < 3; ++i) {
        fmpq_mat_entry_si(identity, i, i, 1);
    }
    assert(identity_order.set_basis(sflint::FmpqMatConstRef(identity)));
    assert_equation_order_matches_uncached(identity_order);

    sflint::FmpqPoly quadratic;
    poly_x2_minus(quadratic, 8);
    silex::NumberField quadratic_field = field_by_polynomial(quadratic);
    silex::Order half(quadratic_field);
    sflint::FmpqMat half_basis(2, 2);
    fmpq_mat_entry_si(half_basis, 0, 0, 1);
    sflint::fmpq_set_si(sflint::fmpq_mat_entry(half_basis, 1, 1), 1, 2);
    assert(half.set_basis(sflint::FmpqMatConstRef(half_basis)));
    silex::Order half_early(quadratic_field);
    assert(half_early.set(half));
    assert_equation_order_matches_uncached(half);
    assert(!half.is_equation_order());
    silex::Order half_late(quadratic_field);
    assert(half_late.set(half));
    assert_equation_order_matches_uncached(half_early);
    assert_equation_order_matches_uncached(half_late);
    assert(!half_late.is_equation_order());
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
    assert(test_is_equation_order_matches_uncached() == 0);
    assert(test_arithmetic() == 0);
    assert(test_integer_and_failure_preservation() == 0);
    assert(test_arithmetic_parent_mismatch_preserves_output() == 0);
    assert(test_undefined_order_rejected() == 0);
    assert(test_swap() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    assert(test_set_coordinates_identity_and_general_basis() == 0);
    return 0;
}
