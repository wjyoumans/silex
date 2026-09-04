#include <silex/class_group.hpp>
#include <silex/order_unit.hpp>

#include <cassert>
#include <iostream>

namespace sflint = silex::flint;

int main() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 2);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());
    assert(order.is_maximal());

    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 2);

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext class_group;
    assert(class_group.compute_candidate(order, sflint::FmpzConstRef(bound),
                                         options));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup units;
    assert(units.set_relation_kernel_units(order, class_group, embeddings,
                                           160));

    auto torsion_order = units.torsion_order();
    auto regulator = units.regulator();
    assert(torsion_order.has_value());
    assert(regulator.has_value());

    std::cout << "O = Z[sqrt(2)]\n";
    std::cout << "relation-kernel units = "
              << class_group.relation_kernel_unit_count() << "\n";
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
