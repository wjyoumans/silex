#pragma once

#include <silex/sunit.hpp>

#include <vector>

namespace silex::detail {

struct SClassGroupStorage {
    Order order;
    std::vector<PrimeIdeal> selected_primes;
    FiniteAbelianGroup group;
    std::vector<FractionalIdeal> invariant_ideals;
    std::vector<FactoredElement> power_witnesses;
    flint::FmpzMat power_selected_exponents{0, 0};
    CertificationMode source_class_certification =
            CertificationMode::unknown;
    ProofState proof_status = ProofState::not_checked;
    bool defined = false;
};

struct SUnitGroupStorage {
    explicit SUnitGroupStorage(const Order& source_order) noexcept
        : order(source_order),
          ordinary_units(source_order) {
    }

    Order order;
    std::vector<PrimeIdeal> selected_primes;
    OrderUnitGroup ordinary_units;
    std::vector<FactoredElement> nonunit_generators;
    flint::FmpzMat nonunit_valuations{0, 0};
    flint::Arb regulator;
    slong regulator_precision = 0;
    CertificationMode source_class_certification =
            CertificationMode::unknown;
    CertificationMode source_unit_certification =
            CertificationMode::unknown;
    ProofState source_relation_saturation = ProofState::not_checked;
    ProofState source_unit_proof = ProofState::not_checked;
    ProofState source_regulator_proof = ProofState::not_checked;
    ProofState proof_status = ProofState::not_checked;
    ProofState regulator_proof_status = ProofState::not_checked;
    const DiagnosticsContext* diagnostics = nullptr;
    bool defined = false;
};

}  // namespace silex::detail
