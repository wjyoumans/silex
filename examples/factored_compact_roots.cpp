#include <silex/factored_element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/number_field.hpp>

#include <cassert>

namespace {
namespace sflint = silex::flint;

silex::NumberField rational_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    silex::NumberField field =
            silex::NumberField::by_polynomial(
                    sflint::FmpqPolyConstRef(polynomial));
    assert(field.is_defined());
    return field;
}

bool element_is_rational(const silex::Element& element,
                         slong numerator,
                         ulong denominator) noexcept {
    sflint::Fmpq expected;
    sflint::Fmpq actual;
    sflint::fmpq_set_si(expected, numerator, denominator);

    auto polynomial = element.to_fmpq_poly();
    assert(polynomial.has_value());
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual), *polynomial, 0);
    return sflint::fmpq_poly_degree(*polynomial) <= 0 &&
           sflint::fmpq_equal(actual, expected);
}

}  // namespace

int main() {
    silex::NumberField field = rational_field();
    silex::Element two(field);
    silex::Element three(field);
    silex::Element value(field);
    assert(two.set_si(2));
    assert(three.set_si(3));

    silex::FactoredElement product(field);
    assert(product.push(two, 3));
    assert(product.push(three, 6));
    assert(product.evaluate(value));
    assert(element_is_rational(value, 5832, 1));

    silex::CompactElement compact(field);
    assert(compact.set_factored_element(product, 3));
    assert(compact.evaluate(value));
    assert(element_is_rational(value, 5832, 1));

    silex::FactoredElement compact_root(field);
    silex::FactoredElement check(field);
    bool is_power = false;
    assert(product.is_power_si(is_power, compact_root, 3,
                               silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(check.pow_si(compact_root, 3));
    assert(check.evaluate(value));
    assert(element_is_rational(value, 5832, 1));

    silex::FactoredElement exact_root(field);
    assert(exact_root.root_si(product, 3));
    assert(exact_root.evaluate(value));
    assert(element_is_rational(value, 18, 1));

    return 0;
}
