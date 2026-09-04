#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/hom.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -radicand);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

void set_linear(silex::Element& element,
                slong constant,
                slong coefficient) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, coefficient);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool element_is_linear(const silex::Element& element,
                       slong constant,
                       slong coefficient) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq actual;
    sflint::Fmpq expected;
    assert(element.get_fmpq_poly(sflint::FmpqPolyRef(polynomial)));

    sflint::fmpq_set_si(expected, constant, 1);
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual),
                                     polynomial, 0);
    if (!sflint::fmpq_equal(actual, expected)) {
        return false;
    }

    sflint::fmpq_set_si(expected, coefficient, 1);
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual),
                                     polynomial, 1);
    return sflint::fmpq_equal(actual, expected) &&
           sflint::fmpq_poly_degree(polynomial) <= 1;
}

void set_matrix_entry_si(sflint::FmpqMat& matrix,
                         slong row,
                         slong col,
                         slong value) noexcept {
    fmpq_set_si(fmpq_mat_entry(matrix.raw(), row, col), value, 1);
}

bool fmpz_mat_entry_is_si(const sflint::FmpzMat& matrix,
                          slong row,
                          slong col,
                          slong value) noexcept {
    return fmpz_cmp_si(fmpz_mat_entry(matrix.raw(), row, col), value) == 0;
}

silex::OrderHom local_order_hom() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::Element generator(field);
    assert(generator.gen());
    silex::FieldHom identity(field, field);
    assert(identity.set_generator_image(generator));

    silex::OrderHom hom(order, order);
    assert(hom.set_field_homomorphism(identity));
    return hom;
}

}  // namespace

int main() {
    silex::NumberField rational = degree_one_field();
    silex::Element rational_gen(rational);
    silex::Element rational_x(rational);
    silex::Element rational_y(rational);
    assert(rational_gen.gen());

    silex::FieldHom empty;
    assert(!empty.is_defined());
    assert(empty.domain() == nullptr);
    assert(empty.codomain() == nullptr);

    silex::FieldHom unset(rational, rational);
    assert(unset.is_defined());
    assert(unset.domain() != nullptr &&
           unset.domain()->has_same_data(rational));
    assert(unset.codomain() != nullptr &&
           unset.codomain()->has_same_data(rational));
    assert(!unset.has_generator_image());
    assert(!unset.generator_image(rational_y));
    assert(!unset.apply(rational_y, rational_gen));
    assert(!unset.is_identity());
    assert(!unset.is_isomorphism());

    assert(unset.set_generator_image(rational_gen));
    assert(unset.has_generator_image());
    assert(unset.is_identity());
    assert(unset.is_isomorphism());
    assert(rational_x.set_si(5));
    assert(unset.apply(rational_y, rational_x));
    assert(rational_y.equal_si(5));

    silex::NumberField quadratic = quadratic_field(2);
    silex::Element theta(quadratic);
    silex::Element minus_theta(quadratic);
    silex::Element x(quadratic);
    silex::Element y(quadratic);
    silex::Element sentinel(quadratic);
    assert(theta.gen());
    assert(minus_theta.negate(theta));

    silex::FieldHom identity(quadratic, quadratic);
    assert(identity.set_generator_image(theta));
    assert(identity.is_identity());
    assert(identity.is_isomorphism());
    set_linear(x, 3, 4);
    assert(identity.apply(y, x));
    assert(element_is_linear(y, 3, 4));

    silex::FieldHom conjugation(quadratic, quadratic);
    assert(conjugation.set_generator_image(minus_theta));
    assert(!conjugation.is_identity());
    assert(conjugation.is_isomorphism());
    assert(conjugation.apply(y, x));
    assert(element_is_linear(y, 3, -4));
    assert(y.set(theta));
    assert(conjugation.apply(y, y));
    assert(y.equal(minus_theta));

    assert(rational_y.set_si(7));
    assert(!conjugation.generator_image(rational_y));
    assert(rational_y.equal_si(7));
    assert(conjugation.generator_image(y));
    assert(y.equal(minus_theta));

    silex::FieldHom copy(quadratic, quadratic);
    assert(copy.set(identity));
    assert(copy.is_identity());
    assert(copy.set(conjugation));
    assert(copy.generator_image(y));
    assert(y.equal(minus_theta));

    silex::FieldHom moved(std::move(copy));
    assert(moved.generator_image(y));
    assert(y.equal(minus_theta));

    assert(copy.set(identity));
    swap(copy, moved);
    assert(copy.generator_image(y));
    assert(y.equal(minus_theta));
    assert(moved.is_identity());

    assert(conjugation.generator_image(y));
    assert(y.equal(minus_theta));
    assert(!conjugation.set_generator_image(rational_gen));
    assert(conjugation.generator_image(y));
    assert(y.equal(minus_theta));

    silex::Element one(quadratic);
    assert(one.one());
    assert(!conjugation.set_generator_image(one));
    assert(conjugation.generator_image(y));
    assert(y.equal(minus_theta));

    silex::NumberField other_quadratic = quadratic_field(2);
    silex::Element theta_other(other_quadratic);
    silex::Element x_other(other_quadratic);
    silex::Element y_other(other_quadratic);
    assert(theta_other.gen());

    silex::FieldHom cross(quadratic, other_quadratic);
    assert(cross.set_generator_image(theta_other));
    assert(!cross.is_identity());
    assert(cross.is_isomorphism());
    set_linear(x, 1, 2);
    assert(cross.apply(y_other, x));
    assert(element_is_linear(y_other, 1, 2));
    assert(!cross.apply(y, x));
    assert(!cross.apply(y_other, x_other));

    silex::Order rational_order = silex::test::equation_order(rational);

    silex::OrderHom empty_order_hom;
    assert(!empty_order_hom.is_defined());
    assert(empty_order_hom.source_order() == nullptr);
    assert(empty_order_hom.target_order() == nullptr);

    silex::OrderHom unset_order_hom(rational_order, rational_order);
    sflint::FmpzMat rational_image_matrix(1, 1);
    assert(unset_order_hom.is_defined());
    assert(silex::same_order_parent(unset_order_hom.source_order(),
                                    &rational_order));
    assert(silex::same_order_parent(unset_order_hom.target_order(),
                                    &rational_order));
    assert(!unset_order_hom.has_field_homomorphism());
    assert(!unset_order_hom.apply(rational_y, rational_x));
    assert(!unset_order_hom.image_matrix(
            sflint::FmpzMatRef(rational_image_matrix)));

    assert(unset_order_hom.set_field_homomorphism(unset));
    assert(unset_order_hom.has_field_homomorphism());
    assert(unset_order_hom.image_matrix(
            sflint::FmpzMatRef(rational_image_matrix)));
    assert(fmpz_mat_entry_is_si(rational_image_matrix, 0, 0, 1));
    assert(rational_x.set_si(7));
    assert(unset_order_hom.apply(rational_y, rational_x));
    assert(rational_y.equal_si(7));

    silex::Order equation_order = silex::test::equation_order(quadratic);

    silex::OrderHom order_identity(equation_order, equation_order);
    assert(order_identity.set_field_homomorphism(identity));
    sflint::FmpzMat order_matrix(2, 2);
    assert(order_identity.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 0, 0, 1));
    assert(fmpz_mat_entry_is_si(order_matrix, 0, 1, 0));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 0, 0));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, 1));
    set_linear(x, 3, 4);
    assert(order_identity.apply(y, x));
    assert(y.equal(x));

    silex::OrderHom order_conjugation(equation_order, equation_order);
    assert(order_conjugation.set_field_homomorphism(conjugation));
    assert(order_conjugation.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 0, 0, 1));
    assert(fmpz_mat_entry_is_si(order_matrix, 0, 1, 0));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 0, 0));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, -1));
    assert(order_conjugation.apply(y, x));
    assert(element_is_linear(y, 3, -4));

    sflint::FmpzMat owned_order_matrix(0, 0);
    assert(order_conjugation.image_matrix(owned_order_matrix));
    assert(sflint::fmpz_mat_nrows(owned_order_matrix) == 2);
    assert(sflint::fmpz_mat_ncols(owned_order_matrix) == 2);
    assert(fmpz_mat_entry_is_si(owned_order_matrix, 1, 1, -1));

    silex::Order narrow_order(quadratic);
    sflint::FmpqMat narrow_basis(2, 2);
    set_matrix_entry_si(narrow_basis, 0, 0, 1);
    set_matrix_entry_si(narrow_basis, 1, 1, 2);
    assert(narrow_order.set_basis(sflint::FmpqMatConstRef(narrow_basis)));

    silex::OrderHom rejected(equation_order, narrow_order);
    assert(!rejected.set_field_homomorphism(identity));
    assert(!rejected.has_field_homomorphism());

    silex::OrderHom preserve(equation_order, equation_order);
    assert(preserve.set_field_homomorphism(identity));
    assert(!preserve.set_field_homomorphism(cross));
    assert(preserve.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 0, 0, 1));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, 1));

    silex::OrderHom copied_order_hom(equation_order, equation_order);
    assert(copied_order_hom.set(order_conjugation));
    assert(copied_order_hom.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, -1));

    silex::OrderHom moved_order_hom(std::move(copied_order_hom));
    assert(moved_order_hom.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, -1));

    assert(copied_order_hom.set(order_identity));
    swap(copied_order_hom, moved_order_hom);
    assert(copied_order_hom.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, -1));
    assert(moved_order_hom.image_matrix(sflint::FmpzMatRef(order_matrix)));
    assert(fmpz_mat_entry_is_si(order_matrix, 1, 1, 1));

    silex::OrderHom local_hom = local_order_hom();
    assert(local_hom.is_defined());
    assert(local_hom.has_field_homomorphism());
    assert(local_hom.source_order() != nullptr);
    assert(local_hom.target_order() != nullptr);
    assert(local_hom.source_order()->parent() != nullptr);
    assert(local_hom.target_order()->parent() != nullptr);
    sflint::FmpzMat local_matrix(1, 1);
    assert(local_hom.image_matrix(sflint::FmpzMatRef(local_matrix)));
    assert(fmpz_mat_entry_is_si(local_matrix, 0, 0, 1));

    silex::Element local_input(*local_hom.source_order()->parent());
    silex::Element local_output(*local_hom.target_order()->parent());
    assert(local_input.set_si(9));
    assert(local_hom.apply(local_output, local_input));
    assert(local_output.equal_si(9));

    return 0;
}
