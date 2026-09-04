#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include <cassert>
#include <iostream>

namespace sflint = silex::flint;

namespace {

bool element_theta_over_two(silex::Element& out) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_set_si(sflint::FmpqRef(coeff), 1, 2);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coeff);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool fmpq_is_si(const sflint::Fmpq& value, slong num, slong den) noexcept {
    sflint::Fmpq expected;
    sflint::fmpq_set_si(sflint::FmpqRef(expected), num, den);
    return sflint::fmpq_equal(value, expected);
}

}  // namespace

int main() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 2);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());

    sflint::Fmpz discriminant;
    assert(order.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, 8));

    silex::Element theta_element(field);
    assert(theta_element.gen());

    silex::OrderElement theta(order);
    assert(theta.set_element(theta_element));

    silex::Ideal theta_ideal(order);
    sflint::Fmpz ideal_norm;
    assert(theta_ideal.set_principal(theta));
    assert(theta_ideal.contains(theta));
    assert(theta_ideal.norm(sflint::FmpzRef(ideal_norm)));
    assert(sflint::fmpz_equal_si(ideal_norm, 2));

    silex::Element theta_half(field);
    assert(element_theta_over_two(theta_half));

    silex::FractionalIdeal fractional(order);
    sflint::Fmpq fractional_norm;
    assert(fractional.set_principal(theta_half));
    assert(fractional.contains(theta_half));
    assert(fractional.norm(sflint::FmpqRef(fractional_norm)));
    assert(fmpq_is_si(fractional_norm, 1, 2));

    std::cout << "K = Q(theta), theta^2 - 2 = 0\n";
    std::cout << "degree = " << field.degree() << "\n";
    std::cout << "disc(Z[theta]) = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(discriminant))
              << "\n";
    std::cout << "Norm((theta)) = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(ideal_norm))
              << "\n";
    std::cout << "Norm((theta/2)) = "
              << sflint::fmpz_get_si(sflint::fmpq_num_ref(fractional_norm))
              << "/"
              << sflint::fmpz_get_si(sflint::fmpq_den_ref(fractional_norm))
              << "\n";
    return 0;
}
