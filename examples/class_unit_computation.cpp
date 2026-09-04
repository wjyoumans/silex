#include <silex/class_group.hpp>
#include <silex/order_unit.hpp>

#include <cassert>
#include <iostream>

namespace {
namespace sflint = silex::flint;

const char* certification_name(silex::CertificationMode mode) noexcept {
    switch (mode) {
    case silex::CertificationMode::unknown:
        return "unknown";
    case silex::CertificationMode::heuristic:
        return "heuristic";
    case silex::CertificationMode::grh:
        return "grh";
    case silex::CertificationMode::proven:
        return "proven";
    }
    return "unknown";
}

const char* proof_state_name(silex::ProofState state) noexcept {
    switch (state) {
    case silex::ProofState::not_checked:
        return "not_checked";
    case silex::ProofState::unavailable:
        return "unavailable";
    case silex::ProofState::verified:
        return "verified";
    }
    return "not_checked";
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
    assert(order.is_maximal());

    sflint::Fmpz factor_base_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(factor_base_bound), 2);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            class_group, order, sflint::FmpzConstRef(factor_base_bound),
            options, 160));

    auto class_order = class_group.order();
    auto torsion_order = units.torsion_order();
    auto regulator = units.regulator();
    assert(class_order.has_value());
    assert(torsion_order.has_value());
    assert(regulator.has_value());

    std::cout << "O = Z[sqrt(2)]\n";
    std::cout << "class group generators = "
              << class_group.generator_count() << "\n";
    std::cout << "class group order = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*class_order))
              << "\n";
    std::cout << "relation-kernel units = "
              << class_group.relation_kernel_unit_count() << "\n";
    std::cout << "unit free rank = " << units.free_rank() << "\n";
    std::cout << "torsion order = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*torsion_order))
              << "\n";
    std::cout << "regulator certified = "
              << sflint::arb_is_positive(*regulator) << "\n";
    std::cout << "class group certification = "
              << certification_name(class_group.certification_status())
              << "\n";
    std::cout << "relation saturation proof = "
              << proof_state_name(class_group.relation_saturation_status())
              << "\n";
    std::cout << "BF zeta proof = "
              << proof_state_name(class_group.zeta_bf_proof_status())
              << "\n";
    std::cout << "unit proof = "
              << proof_state_name(class_group.unit_proof_status())
              << "\n";
    std::cout << "regulator proof = "
              << proof_state_name(class_group.regulator_proof_status())
              << "\n";
    std::cout << "order-unit certification = "
              << certification_name(units.certification_status())
              << "\n";
    std::cout << "local unit proof records = "
              << units.unit_proof_record_count()
              << "\n";
    return 0;
}
