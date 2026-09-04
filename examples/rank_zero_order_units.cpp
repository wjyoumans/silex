#include <silex/order_unit.hpp>

#include <cassert>
#include <iostream>

namespace {
namespace sflint = silex::flint;

int run() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), -1);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    auto torsion_order = units.torsion_order();
    assert(torsion_order.has_value());

    std::cout << "O = Z[i]\n";
    std::cout << "free rank = " << units.free_rank() << "\n";
    std::cout << "torsion order = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*torsion_order))
              << "\n";
    std::cout << "certified = "
              << (units.certification_status() ==
                  silex::CertificationMode::proven)
              << "\n";
    return 0;
}

}  // namespace

int main() {
    return run();
}
