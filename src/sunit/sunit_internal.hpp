#pragma once

#include "../class_group/class_group_internal.hpp"

#include <silex/abelian_group.hpp>

#include <atomic>
#include <vector>

namespace silex::detail {

class TestHookAccess;

// Switch that only detail::TestHookAccess (test/test_hook_access.hpp) can
// set.  Off (0) by default and never set by production code; one relaxed
// atomic load per composition call.
class SUnitTestSeams {
public:
    // Counts down one composition call (an S-class invariant witness or an
    // S-unit generator) per call while armed with n > 0; returns true on the
    // call that reaches zero, which should then fail as an slong exponent
    // overflow would.  Returns false when not armed.
    static bool take_composition_overflow() noexcept {
        slong n = composition_overflow_countdown_.load(
                std::memory_order_relaxed);
        while (n > 0) {
            if (composition_overflow_countdown_.compare_exchange_weak(
                        n, n - 1, std::memory_order_relaxed)) {
                return n == 1;
            }
        }
        return false;
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<slong> composition_overflow_countdown_{0};
};

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

// For each selected prime, the index of the equal factor-base prime, or -1.
bool selected_factor_base_indices(std::vector<slong>& out,
                                  const SUnitClassContext& context) noexcept;

// Factored valuation check of an S-class invariant witness at the selected
// primes: v_{P_j}(witness) == selected_exponents_j + invariant *
// generator_row_k when P_j equals factor-base prime Q_k, and
// selected_exponents_j otherwise.
bool verify_s_class_invariant_witness(
        const SUnitClassContext& context,
        const std::vector<slong>& selected_fb_index,
        flint::FmpzMatConstRef generator_row,
        flint::FmpzConstRef invariant,
        const FactoredElement& witness,
        flint::FmpzMatConstRef selected_exponents,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool build_sunit_class_context(
        SUnitClassBuildResult& result,
        SUnitClassContext& out,
        const ClassGroupContext& class_group,
        const std::vector<PrimeIdeal>& selected_primes,
        const ClassGroupIdealRelationWitnessOptions& witness_options = {})
        noexcept;

}  // namespace silex::detail
