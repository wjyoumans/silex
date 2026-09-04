#include <silex/order_unit.hpp>

#include <cassert>
#include <iostream>

namespace sflint = silex::flint;

int main() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 5);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order maximal_order(field);
    silex::Order equation_order = silex::Order::equation_order(field);
    assert(equation_order.is_defined());
    assert(maximal_order.maximal_order(equation_order));

    silex::OrderUnitGroup units;
    assert(units.compute(maximal_order));

    auto torsion_order = units.torsion_order();
    auto regulator = units.regulator();
    assert(torsion_order.has_value());
    assert(regulator.has_value());

    std::cout << "O = maximal order of Q(sqrt(5))\n";
    std::cout << "free rank = " << units.free_rank() << "\n";
    std::cout << "torsion order = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*torsion_order))
              << "\n";
    std::cout << "regulator certified = "
              << sflint::arb_is_positive(*regulator) << "\n";
    std::cout << "full group certified = "
              << (units.certification_status() ==
                  silex::CertificationMode::proven)
              << "\n";
    return 0;
}
