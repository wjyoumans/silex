#include <silex/order_unit.hpp>

#include <cassert>
#include <iostream>
#include <utility>

namespace sflint = silex::flint;

int main() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 2);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());

    silex::Element theta(field);
    silex::Element epsilon(field);
    assert(theta.gen());
    assert(epsilon.add_si(theta, 1));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup units;
    assert(units.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    auto torsion_order = units.torsion_order();
    auto regulator = units.regulator();
    assert(torsion_order.has_value());
    assert(regulator.has_value());

    std::cout << "O = Z[sqrt(2)]\n";
    std::cout << "free rank = " << units.free_rank() << "\n";
    std::cout << "torsion order = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*torsion_order))
              << "\n";
    std::cout << "subgroup regulator certified = "
              << sflint::arb_is_positive(*regulator) << "\n";
    std::cout << "full group certified = "
              << (units.certification_status() ==
                  silex::CertificationMode::proven)
              << "\n";
    return 0;
}
