#include <silex/aut.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/number_field.hpp>

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

silex::NumberField quadratic_field(slong constant,
                                   slong linear = 0,
                                   slong leading = 1) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, linear);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, leading);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

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
                                     polynomial,
                                     0);
    if (!sflint::fmpq_equal(actual, expected)) {
        return false;
    }

    sflint::fmpq_set_si(expected, coefficient, 1);
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual),
                                     polynomial,
                                     1);
    return sflint::fmpq_equal(actual, expected) &&
           sflint::fmpq_poly_degree(polynomial) <= 1;
}

}  // namespace

int main() {
    silex::FieldAutomorphism empty;
    assert(!empty.is_defined());
    assert(empty.parent() == nullptr);
    assert(!empty.has_homomorphism());

    silex::NumberField rational = degree_one_field();
    silex::Element rational_x(rational);
    silex::Element rational_y(rational);
    silex::FieldHom rational_hom;

    silex::FieldAutomorphism unset(rational);
    assert(unset.is_defined());
    assert(unset.parent() != nullptr &&
           unset.parent()->has_same_data(rational));
    assert(!unset.has_homomorphism());
    assert(!unset.apply(rational_y, rational_x));
    assert(!unset.homomorphism(rational_hom));
    assert(!unset.is_identity());

    assert(unset.set_identity());
    assert(unset.has_homomorphism());
    assert(unset.is_identity());
    assert(unset.homomorphism(rational_hom));
    assert(rational_hom.is_identity());
    assert(rational_x.set_si(9));
    assert(unset.apply(rational_y, rational_x));
    assert(rational_y.equal_si(9));
    assert(!unset.set_quadratic_conjugation());
    assert(unset.is_identity());

    silex::NumberField quadratic = quadratic_field(-2);
    silex::Element theta(quadratic);
    silex::Element x(quadratic);
    silex::Element y(quadratic);
    silex::Element expected(quadratic);
    assert(theta.gen());

    silex::FieldAutomorphism identity(quadratic);
    assert(identity.set_identity());
    assert(identity.is_identity());
    set_linear(x, 3, 4);
    assert(identity.apply(y, x));
    assert(y.equal(x));

    silex::FieldAutomorphism conjugation(quadratic);
    assert(conjugation.set_quadratic_conjugation());
    assert(!conjugation.is_identity());
    assert(conjugation.apply(y, theta));
    assert(expected.negate(theta));
    assert(y.equal(expected));
    assert(conjugation.apply(y, x));
    assert(element_is_linear(y, 3, -4));
    assert(y.set(theta));
    assert(conjugation.apply(y, y));
    assert(y.equal(expected));

    silex::FieldHom extracted;
    assert(conjugation.homomorphism(extracted));
    assert(extracted.apply(expected, x));
    assert(conjugation.apply(y, x));
    assert(expected.equal(y));

    silex::Element s(quadratic);
    silex::Element p(quadratic);
    silex::Element ax(quadratic);
    silex::Element ay(quadratic);
    silex::Element as(quadratic);
    silex::Element ap(quadratic);
    silex::Element check(quadratic);
    set_linear(x, 1, 2);
    set_linear(y, 3, -1);
    assert(s.add(x, y));
    assert(p.multiply(x, y));
    assert(conjugation.apply(ax, x));
    assert(conjugation.apply(ay, y));
    assert(conjugation.apply(as, s));
    assert(conjugation.apply(ap, p));
    assert(check.add(ax, ay));
    assert(as.equal(check));
    assert(check.multiply(ax, ay));
    assert(ap.equal(check));

    silex::NumberField shifted = quadratic_field(-1, 1);
    silex::Element shifted_theta(shifted);
    silex::Element shifted_x(shifted);
    silex::Element shifted_y(shifted);
    assert(shifted_theta.gen());

    silex::FieldAutomorphism shifted_conjugation(shifted);
    assert(shifted_conjugation.set_quadratic_conjugation());
    assert(shifted_conjugation.apply(shifted_y, shifted_theta));
    assert(element_is_linear(shifted_y, -1, -1));
    set_linear(shifted_x, 2, 3);
    assert(shifted_conjugation.apply(shifted_y, shifted_x));
    assert(element_is_linear(shifted_y, -1, -3));

    silex::NumberField cubic = cubic_field();
    silex::Element cubic_x(cubic);
    silex::Element cubic_y(cubic);
    silex::FieldAutomorphism cubic_identity(cubic);
    assert(cubic_identity.set_identity());
    assert(cubic_identity.is_identity());
    assert(!cubic_identity.set_quadratic_conjugation());
    assert(cubic_identity.is_identity());
    set_linear(cubic_x, 1, 2);
    assert(cubic_identity.apply(cubic_y, cubic_x));
    assert(cubic_y.equal(cubic_x));

    silex::FieldAutomorphism copy(quadratic);
    assert(copy.set(identity));
    assert(copy.is_identity());
    assert(copy.set(conjugation));
    assert(!copy.is_identity());
    assert(copy.apply(y, theta));
    assert(element_is_linear(y, 0, -1));
    assert(expected.negate(theta));
    assert(y.equal(expected));

    silex::FieldAutomorphism moved(std::move(copy));
    assert(moved.apply(y, theta));
    assert(y.equal(expected));

    assert(copy.set(identity));
    swap(copy, moved);
    assert(!copy.is_identity());
    assert(moved.is_identity());

    silex::Element rational_out(rational);
    assert(!conjugation.apply(rational_out, theta));
    assert(!conjugation.apply(y, rational_x));

    return 0;
}
