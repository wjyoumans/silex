#pragma once

#include "../class_group/class_group_internal.hpp"

#include <silex/abelian_group.hpp>

#include <vector>

namespace silex::detail {

enum class SUnitClassBuildStage {
    none,
    input_validation,
    class_hnf_basis,
    selected_ideal_relation,
    s_class_presentation,
    s_class_invariant_witnesses,
    relation_kernel,
    valuation_hnf,
    generator_composition,
};

struct SUnitClassBuildResult {
    bool success = false;
    SUnitClassBuildStage stage = SUnitClassBuildStage::none;
    slong selected_index = -1;
    ClassGroupIdealRelationWitnessResult selected_relation;
};

struct SUnitClassContext {
    Order order;
    FactorBase factor_base;
    std::vector<PrimeIdeal> selected_primes;

    WitnessedClassRelationHnfBasis class_hnf;
    flint::FmpzMat selected_relation_rows{0, 0};
    std::vector<FactoredElement> selected_multipliers;
    std::vector<ClassGroupIdealRelationWitnessResult>
            selected_relation_results;
    flint::FmpzMat augmented_relations{0, 0};

    FiniteAbelianGroup s_class_group;
    std::vector<FractionalIdeal> s_class_invariant_ideals;
    std::vector<FactoredElement> s_class_power_witnesses;
    flint::FmpzMat s_class_power_selected_exponents{0, 0};

    flint::FmpzMat relation_kernel{0, 0};
    flint::FmpzMat generator_coefficients{0, 0};
    flint::FmpzMat valuation_rows{0, 0};
    std::vector<FactoredElement> generators_mod_units;

    CertificationMode source_class_certification =
            CertificationMode::unknown;
    ProofState s_class_proof_status = ProofState::not_checked;
    ProofState s_unit_mod_units_proof_status = ProofState::not_checked;
    bool defined = false;
};

const char* sunit_class_build_stage_name(SUnitClassBuildStage stage) noexcept;

bool build_sunit_class_context(
        SUnitClassBuildResult& result,
        SUnitClassContext& out,
        const ClassGroupContext& class_group,
        const std::vector<PrimeIdeal>& selected_primes,
        const ClassGroupIdealRelationWitnessOptions& witness_options = {})
        noexcept;

}  // namespace silex::detail
