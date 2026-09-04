#include <silex/class_group.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/unit.hpp>

#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

void diagnostic_log_callback(void*,
                             silex::DiagnosticsModule module,
                             silex::LogLevel level,
                             const char*,
                             const char* message,
                             const char* detail) noexcept {
    std::cerr << "["
              << silex::diagnostics_module_name(module) << ":"
              << silex::log_level_name(level) << "] " << message;
    if (detail != nullptr) {
        std::cerr << ": " << detail;
    }
    std::cerr << "\n";
}

void count_verbose_callback(void* user,
                            silex::DiagnosticsModule,
                            silex::VerboseLevel,
                            const char*,
                            const char*,
                            const char*) noexcept {
    if (user != nullptr) {
        ++*static_cast<std::size_t*>(user);
    }
}

struct OptionalFmpzSnapshot {
    bool available = false;
    sflint::Fmpz value;
};

struct OptionalArbSnapshot {
    bool available = false;
    sflint::Arb value;
};

struct FactoredElementSnapshot {
    bool defined = false;
    silex::NumberField parent;
    std::vector<silex::Element> factors;
    std::vector<slong> exponents;
};

struct PairPublicSnapshot {
    bool class_defined = false;
    bool class_parent_available = false;
    silex::Order class_parent;
    const silex::DiagnosticsContext* class_diagnostics = nullptr;
    slong analytic_finish_precision = 0;
    slong analytic_product_precision = 0;
    OptionalArbSnapshot analytic_product;

    bool factor_base_available = false;
    silex::FactorBase factor_base;
    bool factor_base_blocks_complete = false;
    std::vector<bool> factor_base_principal_available;
    std::vector<bool> factor_base_principal;
    slong generator_count = 0;

    slong relation_count = 0;
    slong relation_rank = 0;
    slong skipped_dependent_relation_count = 0;
    bool relations_available = false;
    sflint::FmpzMat relations{0, 0};
    std::vector<silex::ClassGroupRelationSource> relation_sources;
    std::vector<bool> relation_generator_available;
    std::vector<silex::Element> relation_generators;

    bool has_presentation = false;
    slong presentation_relation_count = 0;
    slong presentation_generator_count = 0;
    slong invariant_count = 0;
    slong relation_kernel_count = 0;
    sflint::FmpzMat presentation_relations{0, 0};
    std::vector<sflint::Fmpz> invariants;
    OptionalFmpzSnapshot class_order;
    sflint::FmpzMat invariant_generator_matrix{0, 0};
    std::vector<silex::FractionalIdeal> invariant_generators;

    silex::CertificationMode class_certification =
            silex::CertificationMode::unknown;
    silex::ProofState factor_base_generation_status =
            silex::ProofState::not_checked;
    silex::ProofState factor_base_generation_checked_status =
            silex::ProofState::not_checked;
    OptionalFmpzSnapshot factor_base_build_bound;
    OptionalFmpzSnapshot factor_base_generation_bound;
    OptionalFmpzSnapshot factor_base_generation_checked_bound;
    std::vector<silex::ClassGroupFactorBaseGenerationRecord>
            factor_base_generation_records;
    silex::ProofState relation_saturation_status =
            silex::ProofState::not_checked;
    std::vector<silex::ClassGroupRelationSaturationRecord>
            relation_saturation_records;
    silex::ProofState analytic_class_regulator_status =
            silex::ProofState::not_checked;
    silex::ProofState zeta_bf_status = silex::ProofState::not_checked;
    bool zeta_bf_record_available = false;
    ulong zeta_bf_cutoff = 0;
    ulong zeta_bf_max_cutoff = 0;
    slong zeta_bf_requested_precision = 0;
    slong zeta_bf_work_precision = 0;
    sflint::Arb zeta_bf_error_bound;
    silex::ProofState unit_proof_status = silex::ProofState::not_checked;
    silex::ProofState regulator_proof_status =
            silex::ProofState::not_checked;

    bool units_defined = false;
    bool units_parent_available = false;
    silex::Order units_parent;
    const silex::DiagnosticsContext* units_diagnostics = nullptr;
    bool units_set = false;
    slong free_rank = -1;
    silex::CertificationMode units_certification =
            silex::CertificationMode::unknown;
    OptionalFmpzSnapshot torsion_order;
    bool torsion_generator_available = false;
    silex::OrderElement torsion_generator;
    std::vector<FactoredElementSnapshot> free_generators;
    OptionalArbSnapshot regulator;
    std::vector<silex::OrderUnitProofRecord> unit_proof_records;
};

bool capture_factored_element(FactoredElementSnapshot& out,
                              const silex::FactoredElement& value) noexcept {
    out.defined = value.is_defined();
    if (!out.defined) {
        return true;
    }
    if (value.parent() == nullptr) {
        return false;
    }

    out.parent = *value.parent();
    out.factors.reserve(static_cast<std::size_t>(value.length()));
    out.exponents.reserve(static_cast<std::size_t>(value.length()));
    for (slong i = 0; i < value.length(); ++i) {
        const silex::Element* factor = value.factor(i);
        slong exponent = 0;
        if (factor == nullptr || factor->parent() == nullptr ||
            !value.exponent(exponent, i)) {
            return false;
        }
        silex::Element copy(*factor->parent());
        if (!copy.is_defined() || !copy.set(*factor)) {
            return false;
        }
        out.factors.emplace_back(std::move(copy));
        out.exponents.push_back(exponent);
    }
    return true;
}

bool capture_optional_bound(
        OptionalFmpzSnapshot& out,
        bool (silex::ClassGroupContext::*accessor)(sflint::FmpzRef)
                const noexcept,
        const silex::ClassGroupContext& class_group) noexcept {
    out.available = (class_group.*accessor)(sflint::FmpzRef(out.value));
    return true;
}

bool capture_pair_public_snapshot(
        PairPublicSnapshot& out,
        const silex::ClassGroupContext& class_group,
        const silex::OrderUnitGroup& units) noexcept {
    out.class_defined = class_group.is_defined();
    out.class_parent_available = class_group.parent() != nullptr;
    if (out.class_parent_available) {
        out.class_parent = *class_group.parent();
    }
    out.class_diagnostics = class_group.diagnostics();
    out.analytic_finish_precision =
            class_group.analytic_finish_precision();
    out.analytic_product.available =
            class_group.analytic_finish_product(
                    sflint::ArbRef(out.analytic_product.value),
                    out.analytic_product_precision);

    out.factor_base_available = class_group.has_factor_base();
    const silex::FactorBase* factor_base = class_group.factor_base();
    if (out.factor_base_available) {
        if (factor_base == nullptr || !out.factor_base.set(*factor_base)) {
            return false;
        }
        out.factor_base_blocks_complete =
                factor_base->rational_prime_blocks_are_complete();
    } else if (factor_base != nullptr) {
        return false;
    }

    out.generator_count = class_group.generator_count();
    out.factor_base_principal_available.reserve(
            static_cast<std::size_t>(out.generator_count));
    out.factor_base_principal.reserve(
            static_cast<std::size_t>(out.generator_count));
    for (slong i = 0; i < out.generator_count; ++i) {
        bool principal = false;
        const bool available =
                class_group.factor_base_prime_is_principal(principal, i);
        out.factor_base_principal_available.push_back(available);
        out.factor_base_principal.push_back(principal);
    }

    out.relation_count = class_group.relation_count();
    out.relation_rank = class_group.relation_rank();
    out.skipped_dependent_relation_count =
            class_group.skipped_dependent_relation_count();
    out.relations =
            sflint::FmpzMat(out.relation_count, out.generator_count);
    out.relations_available =
            class_group.relations(sflint::FmpzMatRef(out.relations));
    out.relation_sources.reserve(
            static_cast<std::size_t>(out.relation_count));
    out.relation_generator_available.reserve(
            static_cast<std::size_t>(out.relation_count));
    out.relation_generators.reserve(
            static_cast<std::size_t>(out.relation_count));
    for (slong i = 0; i < out.relation_count; ++i) {
        silex::ClassGroupRelationSource source =
                silex::ClassGroupRelationSource::Unknown;
        if (!class_group.relation_source(source, i)) {
            return false;
        }
        out.relation_sources.push_back(source);

        const silex::Element* generator =
                class_group.relation_generator_at(i);
        const bool available = generator != nullptr &&
                generator->parent() != nullptr;
        out.relation_generator_available.push_back(available);
        if (!available) {
            out.relation_generators.emplace_back();
            continue;
        }
        silex::Element copy(*generator->parent());
        if (!copy.is_defined() || !copy.set(*generator)) {
            return false;
        }
        out.relation_generators.emplace_back(std::move(copy));
    }

    out.has_presentation = class_group.has_presentation();
    silex::FiniteAbelianGroup presentation;
    const bool presentation_available =
            class_group.presentation(presentation);
    if (presentation_available != out.has_presentation) {
        return false;
    }
    if (presentation_available) {
        out.presentation_relation_count = presentation.relation_count();
        out.presentation_generator_count = presentation.generator_count();
        out.invariant_count = presentation.invariant_count();
        out.relation_kernel_count = presentation.relation_kernel_count();

        out.presentation_relations = sflint::FmpzMat(
                out.presentation_relation_count,
                out.presentation_generator_count);
        out.invariant_generator_matrix = sflint::FmpzMat(
                out.invariant_count, out.presentation_generator_count);
        if (!presentation.relations(
                    sflint::FmpzMatRef(out.presentation_relations)) ||
            !presentation.invariant_generator_matrix(
                    sflint::FmpzMatRef(out.invariant_generator_matrix))) {
            return false;
        }

        out.invariants.reserve(static_cast<std::size_t>(out.invariant_count));
        out.invariant_generators.reserve(
                static_cast<std::size_t>(out.invariant_count));
        for (slong i = 0; i < out.invariant_count; ++i) {
            sflint::Fmpz invariant;
            if (!presentation.invariant(sflint::FmpzRef(invariant), i) ||
                class_group.parent() == nullptr ||
                class_group.parent()->parent() == nullptr) {
                return false;
            }
            out.invariants.emplace_back(std::move(invariant));

            silex::FractionalIdeal generator(*class_group.parent());
            if (!generator.is_defined() ||
                !class_group.invariant_generator(generator, i)) {
                return false;
            }
            out.invariant_generators.emplace_back(std::move(generator));
        }
        out.class_order.available =
                presentation.order(sflint::FmpzRef(out.class_order.value));
        if (!out.class_order.available) {
            return false;
        }
    }

    out.class_certification = class_group.certification_status();
    out.factor_base_generation_status =
            class_group.factor_base_generation_status();
    out.factor_base_generation_checked_status =
            class_group.factor_base_generation_checked_status();
    capture_optional_bound(
            out.factor_base_build_bound,
            &silex::ClassGroupContext::factor_base_build_bound,
            class_group);
    capture_optional_bound(
            out.factor_base_generation_bound,
            &silex::ClassGroupContext::factor_base_generation_bound,
            class_group);
    capture_optional_bound(
            out.factor_base_generation_checked_bound,
            &silex::ClassGroupContext::factor_base_generation_checked_bound,
            class_group);

    const slong factor_base_record_count =
            class_group.factor_base_generation_record_count();
    if (factor_base_record_count < 0) {
        return false;
    }
    out.factor_base_generation_records.reserve(
            static_cast<std::size_t>(factor_base_record_count));
    for (slong i = 0; i < factor_base_record_count; ++i) {
        auto record = class_group.factor_base_generation_record(i);
        if (!record.has_value()) {
            return false;
        }
        out.factor_base_generation_records.emplace_back(
                std::move(*record));
    }

    out.relation_saturation_status =
            class_group.relation_saturation_status();
    const slong saturation_record_count =
            class_group.relation_saturation_record_count();
    if (saturation_record_count < 0) {
        return false;
    }
    out.relation_saturation_records.reserve(
            static_cast<std::size_t>(saturation_record_count));
    for (slong i = 0; i < saturation_record_count; ++i) {
        auto record = class_group.relation_saturation_record(i);
        if (!record.has_value()) {
            return false;
        }
        out.relation_saturation_records.emplace_back(std::move(*record));
    }

    out.analytic_class_regulator_status =
            class_group.analytic_class_regulator_status();
    out.zeta_bf_status = class_group.zeta_bf_proof_status();
    out.zeta_bf_record_available = class_group.zeta_bf_proof_record(
            out.zeta_bf_cutoff, out.zeta_bf_max_cutoff,
            out.zeta_bf_requested_precision, out.zeta_bf_work_precision,
            sflint::ArbRef(out.zeta_bf_error_bound));
    out.unit_proof_status = class_group.unit_proof_status();
    out.regulator_proof_status = class_group.regulator_proof_status();

    out.units_defined = units.is_defined();
    out.units_parent_available = units.parent() != nullptr;
    if (out.units_parent_available) {
        out.units_parent = *units.parent();
    }
    out.units_diagnostics = units.diagnostics();
    out.units_set = units.is_set();
    out.free_rank = units.free_rank();
    out.units_certification = units.certification_status();
    out.torsion_order.available =
            units.torsion_order(sflint::FmpzRef(out.torsion_order.value));
    out.regulator.available =
            units.regulator(sflint::ArbRef(out.regulator.value));

    if (out.units_parent_available) {
        silex::OrderElement torsion(*units.parent());
        out.torsion_generator_available = units.torsion_generator(torsion);
        if (out.torsion_generator_available) {
            out.torsion_generator = std::move(torsion);
        }
    }

    if (out.free_rank >= 0) {
        if (units.parent() == nullptr || units.parent()->parent() == nullptr) {
            return false;
        }
        out.free_generators.reserve(
                static_cast<std::size_t>(out.free_rank));
        for (slong i = 0; i < out.free_rank; ++i) {
            silex::FactoredElement generator(*units.parent()->parent());
            FactoredElementSnapshot snapshot;
            if (!generator.is_defined() ||
                !units.free_generator(generator, i) ||
                !capture_factored_element(snapshot, generator)) {
                return false;
            }
            out.free_generators.emplace_back(std::move(snapshot));
        }
    }

    const slong unit_record_count = units.unit_proof_record_count();
    if (unit_record_count < 0) {
        return false;
    }
    out.unit_proof_records.reserve(
            static_cast<std::size_t>(unit_record_count));
    for (slong i = 0; i < unit_record_count; ++i) {
        auto record = units.unit_proof_record(i);
        if (!record.has_value()) {
            return false;
        }
        out.unit_proof_records.emplace_back(std::move(*record));
    }
    return true;
}

bool same_fmpz(const sflint::Fmpz& left,
               const sflint::Fmpz& right) noexcept {
    return sflint::fmpz_equal(sflint::FmpzConstRef(left),
                              sflint::FmpzConstRef(right));
}

bool same_optional_fmpz(const OptionalFmpzSnapshot& left,
                        const OptionalFmpzSnapshot& right) noexcept {
    return left.available == right.available &&
           (!left.available || same_fmpz(left.value, right.value));
}

bool same_arb(const sflint::Arb& left,
              const sflint::Arb& right) noexcept {
    return ::arb_equal(left.raw(), right.raw()) != 0;
}

bool same_optional_arb(const OptionalArbSnapshot& left,
                       const OptionalArbSnapshot& right) noexcept {
    return left.available == right.available &&
           (!left.available || same_arb(left.value, right.value));
}

bool same_matrix(const sflint::FmpzMat& left,
                 const sflint::FmpzMat& right) noexcept {
    return sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left),
                                  sflint::FmpzMatConstRef(right));
}

bool same_elements(const std::vector<bool>& left_available,
                   const std::vector<silex::Element>& left,
                   const std::vector<bool>& right_available,
                   const std::vector<silex::Element>& right) noexcept {
    if (left_available != right_available || left.size() != right.size() ||
        left_available.size() != left.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left_available[i] && !left[i].equal(right[i])) {
            return false;
        }
    }
    return true;
}

bool same_factored_element(const FactoredElementSnapshot& left,
                           const FactoredElementSnapshot& right) noexcept {
    if (left.defined != right.defined ||
        left.factors.size() != right.factors.size() ||
        left.exponents != right.exponents) {
        return false;
    }
    if (!left.defined) {
        return true;
    }
    if (!left.parent.has_same_data(right.parent)) {
        return false;
    }
    for (std::size_t i = 0; i < left.factors.size(); ++i) {
        if (!left.factors[i].equal(right.factors[i])) {
            return false;
        }
    }
    return true;
}

bool same_fmpz_vector(const std::vector<sflint::Fmpz>& left,
                      const std::vector<sflint::Fmpz>& right) noexcept {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!same_fmpz(left[i], right[i])) {
            return false;
        }
    }
    return true;
}

bool same_invariant_generators(
        const std::vector<silex::FractionalIdeal>& left,
        const std::vector<silex::FractionalIdeal>& right) noexcept {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!left[i].equal(right[i])) {
            return false;
        }
    }
    return true;
}

bool same_factor_base_generation_records(
        const std::vector<silex::ClassGroupFactorBaseGenerationRecord>& left,
        const std::vector<silex::ClassGroupFactorBaseGenerationRecord>& right)
        noexcept {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!same_fmpz(left[i].p, right[i].p) ||
            left[i].status != right[i].status) {
            return false;
        }
    }
    return true;
}

bool same_saturation_records(
        const std::vector<silex::ClassGroupRelationSaturationRecord>& left,
        const std::vector<silex::ClassGroupRelationSaturationRecord>& right)
        noexcept {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!same_fmpz(left[i].ell, right[i].ell) ||
            left[i].status != right[i].status) {
            return false;
        }
    }
    return true;
}

bool same_unit_proof_records(
        const std::vector<silex::OrderUnitProofRecord>& left,
        const std::vector<silex::OrderUnitProofRecord>& right) noexcept {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!same_fmpz(left[i].ell, right[i].ell) ||
            left[i].status != right[i].status ||
            !same_fmpz(left[i].aux_prime_bound,
                       right[i].aux_prime_bound) ||
            left[i].local_primes != right[i].local_primes ||
            left[i].changed != right[i].changed) {
            return false;
        }
    }
    return true;
}

bool same_pair_public_snapshot(const PairPublicSnapshot& left,
                               const PairPublicSnapshot& right) noexcept {
    if (left.class_defined != right.class_defined ||
        left.class_parent_available != right.class_parent_available ||
        (left.class_parent_available &&
         !left.class_parent.has_same_data(right.class_parent)) ||
        left.class_diagnostics != right.class_diagnostics ||
        left.analytic_finish_precision != right.analytic_finish_precision ||
        left.analytic_product_precision != right.analytic_product_precision ||
        !same_optional_arb(left.analytic_product,
                           right.analytic_product) ||
        left.factor_base_available != right.factor_base_available ||
        (left.factor_base_available &&
         !left.factor_base.equal(right.factor_base)) ||
        left.factor_base_blocks_complete !=
                right.factor_base_blocks_complete ||
        left.factor_base_principal_available !=
                right.factor_base_principal_available ||
        left.factor_base_principal != right.factor_base_principal ||
        left.generator_count != right.generator_count ||
        left.relation_count != right.relation_count ||
        left.relation_rank != right.relation_rank ||
        left.skipped_dependent_relation_count !=
                right.skipped_dependent_relation_count ||
        left.relations_available != right.relations_available ||
        (left.relations_available &&
         !same_matrix(left.relations, right.relations)) ||
        left.relation_sources != right.relation_sources ||
        !same_elements(left.relation_generator_available,
                       left.relation_generators,
                       right.relation_generator_available,
                       right.relation_generators) ||
        left.has_presentation != right.has_presentation ||
        left.presentation_relation_count !=
                right.presentation_relation_count ||
        left.presentation_generator_count !=
                right.presentation_generator_count ||
        left.invariant_count != right.invariant_count ||
        left.relation_kernel_count != right.relation_kernel_count ||
        !same_matrix(left.presentation_relations,
                     right.presentation_relations) ||
        !same_fmpz_vector(left.invariants, right.invariants) ||
        !same_optional_fmpz(left.class_order, right.class_order) ||
        !same_matrix(left.invariant_generator_matrix,
                     right.invariant_generator_matrix) ||
        !same_invariant_generators(left.invariant_generators,
                                   right.invariant_generators) ||
        left.class_certification != right.class_certification ||
        left.factor_base_generation_status !=
                right.factor_base_generation_status ||
        left.factor_base_generation_checked_status !=
                right.factor_base_generation_checked_status ||
        !same_optional_fmpz(left.factor_base_build_bound,
                            right.factor_base_build_bound) ||
        !same_optional_fmpz(left.factor_base_generation_bound,
                            right.factor_base_generation_bound) ||
        !same_optional_fmpz(left.factor_base_generation_checked_bound,
                            right.factor_base_generation_checked_bound) ||
        !same_factor_base_generation_records(
                left.factor_base_generation_records,
                right.factor_base_generation_records) ||
        left.relation_saturation_status !=
                right.relation_saturation_status ||
        !same_saturation_records(left.relation_saturation_records,
                                 right.relation_saturation_records) ||
        left.analytic_class_regulator_status !=
                right.analytic_class_regulator_status ||
        left.zeta_bf_status != right.zeta_bf_status ||
        left.zeta_bf_record_available != right.zeta_bf_record_available ||
        (left.zeta_bf_record_available &&
         (left.zeta_bf_cutoff != right.zeta_bf_cutoff ||
          left.zeta_bf_max_cutoff != right.zeta_bf_max_cutoff ||
          left.zeta_bf_requested_precision !=
                  right.zeta_bf_requested_precision ||
          left.zeta_bf_work_precision != right.zeta_bf_work_precision ||
          !same_arb(left.zeta_bf_error_bound,
                    right.zeta_bf_error_bound))) ||
        left.unit_proof_status != right.unit_proof_status ||
        left.regulator_proof_status != right.regulator_proof_status ||
        left.units_defined != right.units_defined ||
        left.units_parent_available != right.units_parent_available ||
        (left.units_parent_available &&
         !left.units_parent.has_same_data(right.units_parent)) ||
        left.units_diagnostics != right.units_diagnostics ||
        left.units_set != right.units_set ||
        left.free_rank != right.free_rank ||
        left.units_certification != right.units_certification ||
        !same_optional_fmpz(left.torsion_order, right.torsion_order) ||
        left.torsion_generator_available !=
                right.torsion_generator_available ||
        (left.torsion_generator_available &&
         !left.torsion_generator.equal(right.torsion_generator)) ||
        left.free_generators.size() != right.free_generators.size() ||
        !same_optional_arb(left.regulator, right.regulator) ||
        !same_unit_proof_records(left.unit_proof_records,
                                 right.unit_proof_records)) {
        return false;
    }
    for (std::size_t i = 0; i < left.free_generators.size(); ++i) {
        if (!same_factored_element(left.free_generators[i],
                                   right.free_generators[i])) {
            return false;
        }
    }
    return true;
}

bool snapshot_is_wholly_unset(const PairPublicSnapshot& snapshot) noexcept {
    return !snapshot.class_defined && !snapshot.class_parent_available &&
           snapshot.class_diagnostics == nullptr &&
           snapshot.analytic_finish_precision == 0 &&
           !snapshot.analytic_product.available &&
           snapshot.analytic_product_precision == 0 &&
           !snapshot.factor_base_available &&
           !snapshot.factor_base_blocks_complete &&
           snapshot.factor_base_principal_available.empty() &&
           snapshot.factor_base_principal.empty() &&
           snapshot.generator_count == 0 && snapshot.relation_count == 0 &&
           snapshot.relation_rank == 0 &&
           snapshot.skipped_dependent_relation_count == 0 &&
           !snapshot.relations_available && snapshot.relation_sources.empty() &&
           snapshot.relation_generator_available.empty() &&
           snapshot.relation_generators.empty() &&
           !snapshot.has_presentation &&
           snapshot.presentation_relation_count == 0 &&
           snapshot.presentation_generator_count == 0 &&
           snapshot.invariant_count == 0 &&
           snapshot.relation_kernel_count == 0 &&
           snapshot.invariants.empty() && !snapshot.class_order.available &&
           snapshot.invariant_generators.empty() &&
           snapshot.class_certification ==
                   silex::CertificationMode::unknown &&
           snapshot.factor_base_generation_status ==
                   silex::ProofState::not_checked &&
           snapshot.factor_base_generation_checked_status ==
                   silex::ProofState::not_checked &&
           !snapshot.factor_base_build_bound.available &&
           !snapshot.factor_base_generation_bound.available &&
           !snapshot.factor_base_generation_checked_bound.available &&
           snapshot.factor_base_generation_records.empty() &&
           snapshot.relation_saturation_status ==
                   silex::ProofState::not_checked &&
           snapshot.relation_saturation_records.empty() &&
           snapshot.analytic_class_regulator_status ==
                   silex::ProofState::not_checked &&
           snapshot.zeta_bf_status == silex::ProofState::not_checked &&
           !snapshot.zeta_bf_record_available &&
           snapshot.unit_proof_status == silex::ProofState::not_checked &&
           snapshot.regulator_proof_status ==
                   silex::ProofState::not_checked &&
           !snapshot.units_defined && !snapshot.units_parent_available &&
           snapshot.units_diagnostics == nullptr && !snapshot.units_set &&
           snapshot.free_rank == -1 &&
           snapshot.units_certification ==
                   silex::CertificationMode::unknown &&
           !snapshot.torsion_order.available &&
           !snapshot.torsion_generator_available &&
           snapshot.free_generators.empty() && !snapshot.regulator.available &&
           snapshot.unit_proof_records.empty();
}

struct FieldSetup {
    silex::NumberField field;
    silex::Order equation_order;
    silex::Order maximal_order;
};

bool set_monic_polynomial(sflint::FmpqPoly& out,
                          const slong* coefficients,
                          slong degree) noexcept {
    if (coefficients == nullptr || degree < 1) {
        return false;
    }

    sflint::fmpq_poly_zero(out);
    sflint::fmpq_poly_set_coeff_si(out, degree, 1);
    for (slong i = 0; i < degree; ++i) {
        if (coefficients[i] != 0) {
            sflint::fmpq_poly_set_coeff_si(out, i, coefficients[i]);
        }
    }
    return true;
}

FieldSetup setup_from_coefficients(const slong* coefficients,
                                   slong degree) noexcept {
    sflint::FmpqPoly polynomial;
    assert(set_monic_polynomial(polynomial, coefficients, degree));

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    setup.equation_order = silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(setup.equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

bool configure_matrix_options(silex::ClassGroupComputeOptions& options,
                              sflint::Fmpz& factor_base_bound,
                              const silex::Order& order,
                              silex::CertificationMode requested)
        noexcept {
    options = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }

    if (order.degree() >= 3) {
        options.max_candidates = 5000;
        options.max_relations = order.degree() >= 5 ? 1000 : 500;
    }
    options.requested_certification = requested;
    options.zeta_bf_max_cutoff = 20000;
    return true;
}

bool configure_candidate_options(
        silex::ClassGroupCandidateOptions& options,
        sflint::Fmpz& factor_base_bound,
        const silex::Order& order) noexcept {
    options = silex::ClassGroupCandidateOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }
    if (order.degree() >= 3) {
        options.max_candidates = 5000;
        options.max_relations = order.degree() >= 5 ? 1000 : 500;
    }
    return true;
}

bool check_class_group_candidate(const char* name,
                                 const FieldSetup& setup,
                                 slong expected_candidate_order,
                                 slong min_relations) noexcept {
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupCandidateOptions options;
    if (!configure_candidate_options(
                options, factor_base_bound, setup.maximal_order)) {
        std::cerr << name << ": options unavailable\n";
        return false;
    }

    silex::ClassGroupContext context;
    if (!context.compute_candidate(setup.maximal_order,
                                   sflint::FmpzConstRef(factor_base_bound),
                                   options)) {
        std::cerr << name << ": class group computation failed\n";
        return false;
    }
    if (!context.has_presentation()) {
        std::cerr << name << ": missing presentation\n";
        return false;
    }
    if (context.relation_count() < min_relations) {
        std::cerr << name << ": too few relations\n";
        return false;
    }

    sflint::Fmpz class_order;
    if (!context.order(sflint::FmpzRef(class_order)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(class_order),
                               expected_candidate_order)) {
        std::cerr << name << ": unexpected class order "
                  << sflint::fmpz_get_si(
                             sflint::FmpzConstRef(class_order))
                  << " (expected " << expected_candidate_order << ")\n";
        return false;
    }
    return true;
}

bool order_element_has_exact_order(const silex::OrderElement& generator,
                                   slong expected_order) noexcept {
    const silex::Order* const order = generator.parent();
    if (order == nullptr || order->parent() == nullptr ||
        expected_order <= 0) {
        return false;
    }

    silex::Element value(*order->parent());
    silex::Element power(*order->parent());
    if (!value.is_defined() || !power.is_defined() ||
        !generator.get_element(value)) {
        return false;
    }
    for (slong exponent = 1; exponent <= expected_order; ++exponent) {
        sflint::Fmpz exponent_value;
        sflint::fmpz_set_si(sflint::FmpzRef(exponent_value), exponent);
        if (!power.pow_fmpz(
                    value, sflint::FmpzConstRef(exponent_value)) ||
            (exponent < expected_order && power.equal_si(1))) {
            return false;
        }
    }
    return power.equal_si(1);
}

bool check_class_unit_pair(const char* name,
                           const FieldSetup& setup,
                           silex::CertificationMode requested,
                           slong expected_class_order,
                           slong expected_unit_rank,
                           const slong* expected_invariants = nullptr,
                           slong expected_invariant_count = -1,
                           slong expected_torsion_order = 0,
                           bool expect_exact_proof_metadata = false,
                           bool expect_honesty_checkpoint = false)
        noexcept {
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_matrix_options(options, factor_base_bound,
                                  setup.maximal_order, requested)) {
        std::cerr << name << ": options unavailable\n";
        return false;
    }

    silex::DiagnosticsContext diagnostics;
    if (std::getenv("SILEX_TEST_CLASS_UNIT_LOG") != nullptr) {
        silex::diagnostics_context_init(diagnostics);
        const auto modules =
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::class_group) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group);
        silex::diagnostics_set_logging(
                diagnostics, silex::LogLevel::detail, modules,
                diagnostic_log_callback, nullptr);
        options.diagnostics = &diagnostics;
    }

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    if (!units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128)) {
        std::cerr << name << ": class/unit computation failed\n";
        std::cerr << name << ": partial class presentation="
                  << (class_group.has_presentation() ? 1 : 0)
                  << " relations=" << class_group.relation_count()
                  << " kernel_units="
                  << class_group.relation_kernel_unit_count()
                  << " unit_set=" << (units.is_set() ? 1 : 0)
                  << " unit_rank=" << units.free_rank() << "\n";
        return false;
    }
    if (!class_group.has_presentation()) {
        std::cerr << name << ": missing class presentation\n";
        return false;
    }
    if (!units.is_set()) {
        std::cerr << name << ": unit group unset\n";
        return false;
    }
    if (class_group.parent() == nullptr ||
        !class_group.parent()->has_same_data(setup.maximal_order) ||
        units.parent() == nullptr ||
        !units.parent()->has_same_data(setup.maximal_order)) {
        std::cerr << name << ": result parent mismatch\n";
        return false;
    }
    if (units.free_rank() != expected_unit_rank) {
        std::cerr << name << ": unexpected unit rank\n";
        return false;
    }

    sflint::Fmpz class_order;
    if (!class_group.order(sflint::FmpzRef(class_order)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(class_order),
                               expected_class_order)) {
        std::cerr << name << ": unexpected class order\n";
        return false;
    }
    if (expected_invariant_count >= 0) {
        if (class_group.invariant_count() != expected_invariant_count) {
            std::cerr << name << ": unexpected invariant count\n";
            return false;
        }
        for (slong i = 0; i < expected_invariant_count; ++i) {
            sflint::Fmpz invariant;
            if (expected_invariants == nullptr ||
                !class_group.invariant(sflint::FmpzRef(invariant), i) ||
                !sflint::fmpz_equal_si(
                        sflint::FmpzConstRef(invariant),
                        expected_invariants[i])) {
                std::cerr << name << ": unexpected invariant\n";
                return false;
            }
        }
    }

    if (class_group.certification_status() != requested ||
        units.certification_status() != requested) {
        std::cerr << name << ": wrong certification labels\n";
        return false;
    }
    if (expect_honesty_checkpoint) {
        sflint::Fmpz build_bound;
        sflint::Fmpz generation_bound;
        sflint::Fmpz checked_bound;
        if (class_group.factor_base_generation_status() !=
                    silex::ProofState::verified ||
            class_group.factor_base_generation_checked_status() !=
                    silex::ProofState::verified ||
            !class_group.factor_base_build_bound(
                    sflint::FmpzRef(build_bound)) ||
            !class_group.factor_base_generation_bound(
                    sflint::FmpzRef(generation_bound)) ||
            !class_group.factor_base_generation_checked_bound(
                    sflint::FmpzRef(checked_bound)) ||
            sflint::fmpz_cmp(sflint::FmpzConstRef(build_bound),
                             sflint::FmpzConstRef(generation_bound)) >= 0 ||
            sflint::fmpz_cmp(sflint::FmpzConstRef(checked_bound),
                             sflint::FmpzConstRef(generation_bound)) < 0) {
            std::cerr << name << ": missing compact factor-base honesty "
                      << "receipt\n";
            return false;
        }
    }
    if (expected_torsion_order > 0) {
        sflint::Fmpz torsion_order;
        silex::OrderElement torsion_generator(setup.maximal_order);
        sflint::Arb regulator;
        if (!units.torsion_order(sflint::FmpzRef(torsion_order)) ||
            !sflint::fmpz_equal_si(
                    sflint::FmpzConstRef(torsion_order),
                    expected_torsion_order) ||
            !units.torsion_generator(torsion_generator) ||
            !order_element_has_exact_order(
                    torsion_generator, expected_torsion_order) ||
            !units.regulator(sflint::ArbRef(regulator)) ||
            !sflint::arb_is_one(regulator)) {
            std::cerr << name << ": unexpected torsion subgroup\n";
            return false;
        }
    }
    if (expect_exact_proof_metadata) {
        if (class_group.factor_base_generation_status() !=
                    silex::ProofState::verified ||
            class_group.relation_saturation_status() !=
                    silex::ProofState::verified ||
            class_group.analytic_class_regulator_status() !=
                    silex::ProofState::not_checked ||
            class_group.zeta_bf_proof_status() !=
                    silex::ProofState::not_checked ||
            class_group.unit_proof_status() != silex::ProofState::verified ||
            class_group.regulator_proof_status() !=
                    silex::ProofState::verified ||
            class_group.relation_source_count(
                    silex::ClassGroupRelationSource::Saturation) != 0 ||
            class_group.relation_saturation_record_count() != 0) {
            std::cerr << name << ": unexpected exact proof metadata\n";
            return false;
        }
    }
    return true;
}

bool check_exact_edge_failure_atomicity(
        const char* name,
        const FieldSetup& setup,
        silex::CertificationMode requested) noexcept {
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions productive_options;
    if (!configure_matrix_options(
                productive_options, factor_base_bound, setup.maximal_order,
                requested)) {
        std::cerr << name << ": options unavailable\n";
        return false;
    }

    silex::ClassGroupContext populated_class_group;
    silex::OrderUnitGroup populated_units;
    if (!populated_units.compute_with_class_group(
                populated_class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), productive_options,
                128)) {
        std::cerr << name << ": productive pair unavailable\n";
        return false;
    }
    PairPublicSnapshot populated_snapshot;
    if (!capture_pair_public_snapshot(
                populated_snapshot, populated_class_group,
                populated_units)) {
        std::cerr << name << ": productive snapshot unavailable\n";
        return false;
    }

    auto check_failure = [&](const char* suffix,
                             silex::ClassGroupComputeOptions options,
                             slong precision) noexcept {
        silex::ClassGroupContext fresh_class_group;
        silex::OrderUnitGroup fresh_units;
        PairPublicSnapshot fresh_before;
        PairPublicSnapshot fresh_after;
        if (!capture_pair_public_snapshot(
                    fresh_before, fresh_class_group, fresh_units) ||
            fresh_units.compute_with_class_group(
                    fresh_class_group, setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound), options,
                    precision) ||
            !capture_pair_public_snapshot(
                    fresh_after, fresh_class_group, fresh_units) ||
            !snapshot_is_wholly_unset(fresh_after) ||
            !same_pair_public_snapshot(fresh_before, fresh_after)) {
            std::cerr << name << ": " << suffix
                      << " fresh atomicity mismatch\n";
            return false;
        }

        PairPublicSnapshot before;
        PairPublicSnapshot after;
        if (!capture_pair_public_snapshot(
                    before, populated_class_group, populated_units) ||
            !same_pair_public_snapshot(populated_snapshot, before) ||
            populated_units.compute_with_class_group(
                    populated_class_group, setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound), options,
                    precision) ||
            !capture_pair_public_snapshot(
                    after, populated_class_group, populated_units) ||
            !same_pair_public_snapshot(before, after)) {
            std::cerr << name << ": " << suffix
                      << " populated atomicity mismatch\n";
            return false;
        }
        return true;
    };

    silex::ClassGroupComputeOptions exhausted_options = productive_options;
    exhausted_options.max_candidates = 0;
    exhausted_options.max_relations = 0;
    return check_failure("zero-resource", exhausted_options, 128) &&
           check_failure("invalid-precision", productive_options, 0);
}

bool check_candidate_boundary_and_proven_pair(
        const char* name,
        const FieldSetup& setup,
        slong expected_candidate_order,
        slong expected_class_order,
        slong expected_unit_rank,
        const slong* expected_invariants = nullptr,
        slong expected_invariant_count = -1) noexcept {
    if (!check_class_group_candidate(
                name, setup, expected_candidate_order, 1)) {
        return false;
    }

    return check_class_unit_pair(
            name, setup, silex::CertificationMode::proven,
            expected_class_order, expected_unit_rank, expected_invariants,
            expected_invariant_count, 0, false, true);
}

bool check_zero_resource_failure_atomicity(
        const char* name,
        const FieldSetup& setup,
        silex::CertificationMode requested) noexcept {
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_matrix_options(
                options, factor_base_bound, setup.maximal_order, requested)) {
        std::cerr << name << ": options unavailable\n";
        return false;
    }
    options.max_candidates = 0;
    options.max_relations = 0;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    PairPublicSnapshot before;
    PairPublicSnapshot after;
    if (!capture_pair_public_snapshot(before, class_group, units) ||
        units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128) ||
        !capture_pair_public_snapshot(after, class_group, units) ||
        !snapshot_is_wholly_unset(after) ||
        !same_pair_public_snapshot(before, after)) {
        std::cerr << name << ": zero-resource atomicity mismatch\n";
        return false;
    }
    return true;
}

int test_public_paired_certification_request_validation() {
    const slong cubic_disc81[] = {1, -3, 0};
    FieldSetup setup = setup_from_coefficients(cubic_disc81, 3);

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions productive_options;
    if (!configure_matrix_options(
                productive_options, factor_base_bound, setup.maximal_order,
                silex::CertificationMode::proven)) {
        std::cerr << "paired certification validation: options unavailable\n";
        return 1;
    }

    std::size_t verbose_callback_count = 0;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_verbose(
            diagnostics, silex::VerboseLevel::trace,
            silex::diagnostics_all_modules, count_verbose_callback,
            &verbose_callback_count);
    productive_options.diagnostics = &diagnostics;

    silex::ClassGroupContext populated_class_group;
    silex::OrderUnitGroup populated_units;
    if (!populated_units.compute_with_class_group(
                populated_class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), productive_options,
                128) ||
        populated_class_group.certification_status() !=
                silex::CertificationMode::proven ||
        populated_units.certification_status() !=
                silex::CertificationMode::proven) {
        std::cerr << "paired certification validation: productive proven "
                     "pair unavailable\n";
        return 1;
    }

    PairPublicSnapshot productive_snapshot;
    if (!capture_pair_public_snapshot(
                productive_snapshot, populated_class_group,
                populated_units)) {
        std::cerr << "paired certification validation: productive snapshot "
                     "unavailable\n";
        return 1;
    }

    auto check_failed_mode = [&](const char* name,
                                 silex::CertificationMode requested,
                                 bool exhaust_resources,
                                 bool expect_verbose_activity) noexcept {
        silex::ClassGroupComputeOptions options = productive_options;
        options.requested_certification = requested;
        if (exhaust_resources) {
            options.max_candidates = 0;
            options.max_relations = 0;
        }

        silex::ClassGroupContext fresh_class_group;
        silex::OrderUnitGroup fresh_units;
        PairPublicSnapshot fresh_before;
        if (!capture_pair_public_snapshot(
                    fresh_before, fresh_class_group, fresh_units) ||
            !snapshot_is_wholly_unset(fresh_before)) {
            std::cerr << name << ": fresh pre-state is not wholly unset\n";
            return false;
        }

        verbose_callback_count = 0;
        const bool fresh_result = fresh_units.compute_with_class_group(
                fresh_class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128);
        const std::size_t fresh_callback_count = verbose_callback_count;
        PairPublicSnapshot fresh_after;
        if (fresh_result ||
            !capture_pair_public_snapshot(
                    fresh_after, fresh_class_group, fresh_units) ||
            !snapshot_is_wholly_unset(fresh_after) ||
            !same_pair_public_snapshot(fresh_before, fresh_after) ||
            (expect_verbose_activity
                     ? fresh_callback_count == 0
                     : fresh_callback_count != 0)) {
            std::cerr << name << ": fresh failure contract mismatch\n";
            return false;
        }

        PairPublicSnapshot populated_before;
        if (!capture_pair_public_snapshot(
                    populated_before, populated_class_group,
                    populated_units) ||
            !same_pair_public_snapshot(
                    productive_snapshot, populated_before)) {
            std::cerr << name << ": populated pre-state drifted\n";
            return false;
        }

        verbose_callback_count = 0;
        const bool populated_result =
                populated_units.compute_with_class_group(
                        populated_class_group, setup.maximal_order,
                        sflint::FmpzConstRef(factor_base_bound), options,
                        128);
        const std::size_t populated_callback_count = verbose_callback_count;
        PairPublicSnapshot populated_after;
        if (populated_result ||
            !capture_pair_public_snapshot(
                    populated_after, populated_class_group,
                    populated_units) ||
            !same_pair_public_snapshot(
                    populated_before, populated_after) ||
            (expect_verbose_activity
                     ? populated_callback_count == 0
                     : populated_callback_count != 0)) {
            std::cerr << name << ": populated failure contract mismatch\n";
            return false;
        }
        return true;
    };

    if (!check_failed_mode(
                "paired proven resource exhaustion",
                silex::CertificationMode::proven, true, true) ||
        !check_failed_mode(
                "paired heuristic request rejection",
                silex::CertificationMode::heuristic, false, false) ||
        !check_failed_mode(
                "paired GRH resource exhaustion",
                silex::CertificationMode::grh, true, true)) {
        return 1;
    }

    const slong imaginary_quadratic_47[] = {47, 0};
    FieldSetup grh_setup =
            setup_from_coefficients(imaginary_quadratic_47, 2);
    sflint::Fmpz grh_requested_bound;
    silex::ClassGroupComputeOptions grh_options;
    if (!configure_matrix_options(
                grh_options, grh_requested_bound, grh_setup.maximal_order,
                silex::CertificationMode::grh)) {
        std::cerr << "paired GRH computation: options unavailable\n";
        return 1;
    }
    silex::ClassGroupContext grh_class_group;
    silex::OrderUnitGroup grh_units;
    if (!grh_units.compute_with_class_group(
                grh_class_group, grh_setup.maximal_order,
                sflint::FmpzConstRef(grh_requested_bound), grh_options,
                128)) {
        std::cerr << "paired GRH computation: computation failed\n";
        return 1;
    }
    if (grh_class_group.certification_status() !=
                silex::CertificationMode::grh ||
        grh_units.certification_status() != silex::CertificationMode::grh) {
        std::cerr << "paired GRH computation: wrong certification labels\n";
        return 1;
    }

    sflint::Fmpz grh_class_order;
    sflint::Fmpz grh_factor_base_bound;
    const silex::FactorBase* const grh_factor_base =
            grh_class_group.factor_base();
    if (!grh_class_group.order(sflint::FmpzRef(grh_class_order)) ||
        sflint::fmpz_cmp_ui(
                sflint::FmpzConstRef(grh_class_order), 5) != 0 ||
        grh_units.free_rank() != 0 || grh_factor_base == nullptr ||
        !grh_class_group.factor_base_build_bound(
                sflint::FmpzRef(grh_factor_base_bound))) {
        std::cerr << "paired GRH computation: unexpected result data\n";
        return 1;
    }

    silex::FactorBase expected_grh_factor_base(grh_setup.maximal_order);
    if (!expected_grh_factor_base.build_prime_ideal_norm_bounded(
                sflint::FmpzConstRef(grh_factor_base_bound)) ||
        expected_grh_factor_base.length() != grh_factor_base->length()) {
        std::cerr << "paired GRH computation: factor base shape mismatch\n";
        return 1;
    }
    silex::PrimeIdeal prime(grh_setup.maximal_order);
    for (slong i = 0; i < grh_factor_base->length(); ++i) {
        if (!grh_factor_base->prime(prime, i) ||
            !expected_grh_factor_base.contains(prime)) {
            std::cerr << "paired GRH computation: factor base membership "
                         "mismatch\n";
            return 1;
        }
    }

    auto check_grh_pair = [&](const char* name,
                              const slong* coefficients,
                              slong degree,
                              ulong expected_class_order,
                              slong expected_unit_rank) noexcept {
        FieldSetup local_setup =
                setup_from_coefficients(coefficients, degree);
        sflint::Fmpz local_bound;
        silex::ClassGroupComputeOptions local_options;
        if (!configure_matrix_options(
                    local_options, local_bound, local_setup.maximal_order,
                    silex::CertificationMode::grh)) {
            std::cerr << name << ": options unavailable\n";
            return false;
        }

        silex::ClassGroupContext local_class_group;
        silex::OrderUnitGroup local_units;
        sflint::Fmpz local_class_order;
        if (!local_units.compute_with_class_group(
                    local_class_group, local_setup.maximal_order,
                    sflint::FmpzConstRef(local_bound), local_options, 128) ||
            local_class_group.certification_status() !=
                    silex::CertificationMode::grh ||
            local_units.certification_status() !=
                    silex::CertificationMode::grh ||
            !local_class_group.order(
                    sflint::FmpzRef(local_class_order)) ||
            sflint::fmpz_cmp_ui(
                    sflint::FmpzConstRef(local_class_order),
                    expected_class_order) != 0 ||
            local_units.free_rank() != expected_unit_rank) {
            std::cerr << name << ": result mismatch\n";
            return false;
        }
        return true;
    };

    const slong real_quadratic_5[] = {-5, 0};
    const slong cubic_disc_minus_23[] = {-1, -1, 0};
    if (!check_grh_pair(
                "paired GRH real quadratic", real_quadratic_5, 2, 1, 1) ||
        !check_grh_pair(
                "paired GRH cubic", cubic_disc_minus_23, 3, 1, 1)) {
        return 1;
    }
    return 0;
}

int test_public_paired_cross_order_failure_preserves_existing_outputs() {
    const slong real_quadratic_5[] = {-5, 0};
    FieldSetup populated_setup =
            setup_from_coefficients(real_quadratic_5, 2);

    sflint::Fmpz populated_factor_base_bound;
    silex::ClassGroupComputeOptions populated_options;
    if (!configure_matrix_options(
                populated_options, populated_factor_base_bound,
                populated_setup.maximal_order,
                silex::CertificationMode::proven)) {
        std::cerr << "paired cross-order atomicity: populated options "
                     "unavailable\n";
        return 1;
    }

    silex::ClassGroupContext populated_class_group;
    silex::OrderUnitGroup populated_units;
    if (!populated_units.compute_with_class_group(
                populated_class_group, populated_setup.maximal_order,
                sflint::FmpzConstRef(populated_factor_base_bound),
                populated_options, 128) ||
        populated_class_group.certification_status() !=
                silex::CertificationMode::proven ||
        populated_units.certification_status() !=
                silex::CertificationMode::proven) {
        std::cerr << "paired cross-order atomicity: populated proven pair "
                     "unavailable\n";
        return 1;
    }

    const slong populated_free_rank = populated_units.free_rank();
    const slong populated_relation_kernel_count =
            populated_class_group.relation_kernel_unit_count();
    if (populated_free_rank <= 0 ||
        populated_relation_kernel_count < populated_free_rank) {
        std::cerr << "paired cross-order atomicity: relation-kernel "
                     "witnesses unavailable\n";
        return 1;
    }

    auto capture_relation_kernel_units =
            [&](std::vector<FactoredElementSnapshot>& out) noexcept {
        const slong count =
                populated_class_group.relation_kernel_unit_count();
        if (count < 0) {
            return false;
        }
        out.clear();
        out.reserve(static_cast<std::size_t>(count));
        for (slong i = 0; i < count; ++i) {
            silex::FactoredElement witness(populated_setup.field);
            FactoredElementSnapshot snapshot;
            if (!populated_class_group.relation_kernel_unit(witness, i) ||
                !capture_factored_element(snapshot, witness)) {
                return false;
            }
            out.emplace_back(std::move(snapshot));
        }
        return true;
    };

    std::vector<FactoredElementSnapshot> relation_kernel_units_before;
    if (!capture_relation_kernel_units(relation_kernel_units_before)) {
        std::cerr << "paired cross-order atomicity: witness baseline "
                     "unavailable\n";
        return 1;
    }

    PairPublicSnapshot populated_before;
    if (!capture_pair_public_snapshot(
                populated_before, populated_class_group, populated_units)) {
        std::cerr << "paired cross-order atomicity: populated snapshot "
                     "unavailable\n";
        return 1;
    }

    const slong cubic_disc81[] = {1, -3, 0};
    FieldSetup attempted_setup = setup_from_coefficients(cubic_disc81, 3);
    const silex::Order* class_parent = populated_class_group.parent();
    const silex::Order* units_parent = populated_units.parent();
    if (class_parent == nullptr || units_parent == nullptr ||
        !class_parent->has_same_data(populated_setup.maximal_order) ||
        !units_parent->has_same_data(populated_setup.maximal_order) ||
        class_parent->has_same_data(attempted_setup.maximal_order) ||
        units_parent->has_same_data(attempted_setup.maximal_order)) {
        std::cerr << "paired cross-order atomicity: order identities "
                     "not distinct\n";
        return 1;
    }

    sflint::Fmpz attempted_factor_base_bound;
    silex::ClassGroupComputeOptions attempted_options;
    if (!configure_matrix_options(
                attempted_options, attempted_factor_base_bound,
                attempted_setup.maximal_order,
                silex::CertificationMode::proven)) {
        std::cerr << "paired cross-order atomicity: attempted options "
                     "unavailable\n";
        return 1;
    }
    attempted_options.max_candidates = 0;
    attempted_options.max_relations = 0;

    std::size_t verbose_callback_count = 0;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_verbose(
            diagnostics, silex::VerboseLevel::trace,
            silex::diagnostics_all_modules, count_verbose_callback,
            &verbose_callback_count);
    attempted_options.diagnostics = &diagnostics;

    verbose_callback_count = 0;
    const bool attempted_result = populated_units.compute_with_class_group(
            populated_class_group, attempted_setup.maximal_order,
            sflint::FmpzConstRef(attempted_factor_base_bound),
            attempted_options, 128);
    const std::size_t attempted_callback_count = verbose_callback_count;
    if (attempted_result || attempted_callback_count == 0) {
        std::cerr << "paired cross-order atomicity: attempted call did not "
                     "fail after observable work\n";
        return 1;
    }

    PairPublicSnapshot populated_after;
    if (!capture_pair_public_snapshot(
                populated_after, populated_class_group, populated_units)) {
        std::cerr << "paired cross-order atomicity: post-failure snapshot "
                     "unavailable\n";
        return 1;
    }

    std::vector<FactoredElementSnapshot> relation_kernel_units_after;
    if (!capture_relation_kernel_units(relation_kernel_units_after) ||
        !same_pair_public_snapshot(populated_before, populated_after) ||
        relation_kernel_units_before.size() !=
                relation_kernel_units_after.size()) {
        std::cerr << "paired cross-order atomicity: populated state "
                     "changed\n";
        return 1;
    }
    for (std::size_t i = 0; i < relation_kernel_units_before.size(); ++i) {
        if (!same_factored_element(relation_kernel_units_before[i],
                                   relation_kernel_units_after[i])) {
            std::cerr << "paired cross-order atomicity: relation-kernel "
                         "witness changed\n";
            return 1;
        }
    }
    return 0;
}

int test_stable_class_group_rows() {
    const slong quartic[] = {2, 0, -1, -1};
    const slong sextic[] = {-1, -1, 0, 0, 0, 0};

    FieldSetup quartic_setup = setup_from_coefficients(quartic, 4);
    if (!check_class_group_candidate("quartic class group", quartic_setup, 1,
                                     9)) {
        return 1;
    }

    FieldSetup sextic_setup = setup_from_coefficients(sextic, 6);
    if (!check_class_group_candidate("sextic class group", sextic_setup, 1,
                                     10)) {
        return 1;
    }
    return 0;
}

int test_release_target_proven_transactions() {
    const slong degree_one[] = {0};
    const slong real_quadratic_5[] = {-5, 0};
    const slong real_quadratic_210[] = {-210, 0};
    const slong imaginary_quadratic_14[] = {14, 0};
    const slong imaginary_quadratic_47[] = {47, 0};
    const slong real_quadratic_210_invariants[] = {2, 2};
    const slong imaginary_quadratic_14_invariants[] = {4};
    const slong imaginary_quadratic_47_invariants[] = {5};

    FieldSetup degree_one_setup =
            setup_from_coefficients(degree_one, 1);
    if (!check_class_unit_pair(
                "degree-one release-target proven transaction",
                degree_one_setup, silex::CertificationMode::proven, 1, 0,
                nullptr, 0)) {
        return 1;
    }

    FieldSetup real_quadratic_5_setup =
            setup_from_coefficients(real_quadratic_5, 2);
    if (!check_class_unit_pair(
                "real quadratic 5 release-target proven transaction",
                real_quadratic_5_setup, silex::CertificationMode::proven, 1,
                1, nullptr, 0)) {
        return 1;
    }

    FieldSetup real_quadratic_210_setup =
            setup_from_coefficients(real_quadratic_210, 2);
    if (!check_class_unit_pair(
                "real quadratic 210 release-target proven transaction",
                real_quadratic_210_setup,
                silex::CertificationMode::proven, 4, 1,
                real_quadratic_210_invariants, 2)) {
        return 1;
    }

    FieldSetup imaginary_quadratic_47_setup =
            setup_from_coefficients(imaginary_quadratic_47, 2);
    if (!check_class_unit_pair(
                "imaginary quadratic 47 release-target proven transaction",
                imaginary_quadratic_47_setup,
                silex::CertificationMode::proven, 5, 0,
                imaginary_quadratic_47_invariants, 1, 2, true)) {
        return 1;
    }

    FieldSetup imaginary_quadratic_14_setup =
            setup_from_coefficients(imaginary_quadratic_14, 2);
    if (!check_class_unit_pair(
                "imaginary quadratic 14 release-target proven transaction",
                imaginary_quadratic_14_setup,
                silex::CertificationMode::proven, 4, 0,
                imaginary_quadratic_14_invariants, 1, 2, true)) {
        return 1;
    }

    return 0;
}

int test_exact_edge_transactions() {
    const slong degree_one[] = {0};
    const slong canonical_discriminant_minus_three[] = {3, 0};
    const slong shifted_discriminant_minus_three[] = {1, -1};
    const slong discriminant_minus_four[] = {1, 0};

    FieldSetup degree_one_setup =
            setup_from_coefficients(degree_one, 1);
    if (!check_class_unit_pair(
                "degree-one exact GRH transaction", degree_one_setup,
                silex::CertificationMode::grh, 1, 0, nullptr, 0, 2, true) ||
        !check_exact_edge_failure_atomicity(
                "degree-one exact GRH transaction", degree_one_setup,
                silex::CertificationMode::grh)) {
        return 1;
    }

    FieldSetup canonical_minus_three_setup = setup_from_coefficients(
            canonical_discriminant_minus_three, 2);
    if (!check_class_unit_pair(
                "canonical discriminant -3 proven transaction",
                canonical_minus_three_setup,
                silex::CertificationMode::proven, 1, 0, nullptr, 0, 6,
                true) ||
        !check_class_unit_pair(
                "canonical discriminant -3 GRH preservation",
                canonical_minus_three_setup,
                silex::CertificationMode::grh, 1, 0, nullptr, 0, 6,
                false) ||
        !check_exact_edge_failure_atomicity(
                "canonical discriminant -3 proven transaction",
                canonical_minus_three_setup,
                silex::CertificationMode::proven)) {
        return 1;
    }

    FieldSetup shifted_minus_three_setup = setup_from_coefficients(
            shifted_discriminant_minus_three, 2);
    if (!check_class_unit_pair(
                "shifted discriminant -3 proven preservation",
                shifted_minus_three_setup,
                silex::CertificationMode::proven, 1, 0, nullptr, 0, 6,
                true)) {
        return 1;
    }

    FieldSetup minus_four_setup = setup_from_coefficients(
            discriminant_minus_four, 2);
    return check_class_unit_pair(
                   "discriminant -4 proven preservation", minus_four_setup,
                   silex::CertificationMode::proven, 1, 0, nullptr, 0, 4,
                   true)
            ? 0
            : 1;
}

int test_higher_degree_completion_boundaries() {
    // These fields exercise completion paths that previously failed closed
    // despite having valid, certifiable class/unit results.
    const slong cubic[] = {-1, 2, -1};
    const slong cubic_nontrivial[] = {-5, -2, 0};
    const slong cubic_seeded_h4[] = {2, -1, 4};
    const slong cubic_disc81[] = {1, -3, 0};
    const slong cubic_disc2213[] = {-1, -8, -4};
    const slong quartic[] = {2, 0, -1, -1};
    const slong quartic_disc70640[] = {-3, 2, 1, 4};
    const slong quartic_disc223479[] = {3, -3, -3, -4};
    const slong quartic_disc35019[] = {-5, 3, -2, 0};
    const slong quartic_disc1412343[] = {4, 1, 3, -8};
    const slong quartic_disc6067408[] = {7, 2, -3, -8};
    const slong quartic_x4_minus_x_minus_1[] = {-1, -1, 0, 0};
    const slong quintic[] = {-1, -1, 0, 0, 0};
    const slong quintic_disc11119[] = {1, -3, 1, 1, -2};
    const slong quintic_disc401370255[] = {3, -3, -3, -6, -7};
    const slong cubic_nontrivial_invariants[] = {2};
    const slong quartic_disc1412343_invariants[] = {4};

    FieldSetup cubic_setup = setup_from_coefficients(cubic, 3);
    if (!check_class_unit_pair("random cubic proven completion diagnostic",
                               cubic_setup,
                               silex::CertificationMode::proven, 1, 1)) {
        return 1;
    }

    FieldSetup cubic_nontrivial_setup =
            setup_from_coefficients(cubic_nontrivial, 3);
    if (!check_class_unit_pair(
                "nontrivial cubic proven completion regression",
                cubic_nontrivial_setup, silex::CertificationMode::proven, 2,
                1, cubic_nontrivial_invariants, 1, 0, false, true)) {
        return 1;
    }

    FieldSetup cubic_seeded_h4_setup =
            setup_from_coefficients(cubic_seeded_h4, 3);
    if (!check_class_unit_pair(
                "deterministic random cubic proven completion regression",
                cubic_seeded_h4_setup, silex::CertificationMode::proven, 1,
                1, nullptr, -1, 0, false, true)) {
        return 1;
    }

    FieldSetup cubic_disc81_setup = setup_from_coefficients(cubic_disc81, 3);
    if (!check_class_unit_pair(
                "cubic discriminant 81 proven completion diagnostic",
                cubic_disc81_setup, silex::CertificationMode::proven, 1, 2)) {
        return 1;
    }

    FieldSetup cubic_disc2213_setup =
            setup_from_coefficients(cubic_disc2213, 3);
    if (!check_class_unit_pair(
                "cubic discriminant 2213 proven completion diagnostic",
                cubic_disc2213_setup, silex::CertificationMode::proven, 1,
                2)) {
        return 1;
    }

    FieldSetup quartic_setup = setup_from_coefficients(quartic, 4);
    if (!check_candidate_boundary_and_proven_pair(
                "quartic proven completion regression", quartic_setup, 1,
                1, 1)) {
        return 1;
    }

    FieldSetup quartic_disc70640_setup =
            setup_from_coefficients(quartic_disc70640, 4);
    if (!check_candidate_boundary_and_proven_pair(
                "quartic discriminant -70640 proven completion regression",
                quartic_disc70640_setup, 3, 1, 2)) {
        return 1;
    }

    FieldSetup quartic_disc223479_setup =
            setup_from_coefficients(quartic_disc223479, 4);
    if (!check_candidate_boundary_and_proven_pair(
                "quartic discriminant -223479 proven completion regression",
                quartic_disc223479_setup, 1, 1, 2)) {
        return 1;
    }

    FieldSetup quartic_disc35019_setup =
            setup_from_coefficients(quartic_disc35019, 4);
    if (!check_candidate_boundary_and_proven_pair(
                "quartic discriminant -35019 proven completion regression",
                quartic_disc35019_setup, 1, 1, 2)) {
        return 1;
    }

    FieldSetup quartic_disc1412343_setup =
            setup_from_coefficients(quartic_disc1412343, 4);
    if (!check_candidate_boundary_and_proven_pair(
                "quartic discriminant -1412343 proven completion regression",
                quartic_disc1412343_setup, 4, 4, 2,
                quartic_disc1412343_invariants, 1)) {
        return 1;
    }

    FieldSetup quartic_disc6067408_setup =
            setup_from_coefficients(quartic_disc6067408, 4);
    if (!check_zero_resource_failure_atomicity(
                "quartic discriminant 6067408 zero-resource regression",
                quartic_disc6067408_setup,
                silex::CertificationMode::proven)) {
        return 1;
    }

    FieldSetup quartic_x4_minus_x_minus_1_setup =
            setup_from_coefficients(quartic_x4_minus_x_minus_1, 4);
    if (!check_class_unit_pair(
                "quartic x^4 - x - 1 proven completion diagnostic",
                quartic_x4_minus_x_minus_1_setup,
                silex::CertificationMode::proven, 1, 2)) {
        return 1;
    }

    FieldSetup quintic_setup = setup_from_coefficients(quintic, 5);
    if (!check_class_unit_pair("quintic proven completion diagnostic",
                               quintic_setup,
                               silex::CertificationMode::proven, 1, 2)) {
        return 1;
    }

    FieldSetup quintic_disc11119_setup =
            setup_from_coefficients(quintic_disc11119, 5);
    if (!check_class_unit_pair(
                "quintic discriminant -11119 proven completion diagnostic",
                quintic_disc11119_setup, silex::CertificationMode::proven, 1,
                3)) {
        return 1;
    }

    FieldSetup quintic_disc401370255_setup =
            setup_from_coefficients(quintic_disc401370255, 5);
    if (!check_candidate_boundary_and_proven_pair(
                "quintic discriminant -401370255 proven completion regression",
                quintic_disc401370255_setup, 1, 1, 3)) {
        return 1;
    }
    return 0;
}

int test_expanded_random_class_group_candidates() {
    const slong cubic_h8[] = {-1, -8, -4};
    const slong quartic_h4[] = {2, -2, -1, 1};

    FieldSetup cubic_setup = setup_from_coefficients(cubic_h8, 3);
    if (!check_class_group_candidate("random cubic H8 class candidate",
                                     cubic_setup, 1, 10)) {
        return 1;
    }

    FieldSetup quartic_setup = setup_from_coefficients(quartic_h4, 4);
    if (!check_class_group_candidate("random quartic H4 class candidate",
                                     quartic_setup, 1, 10)) {
        return 1;
    }

    return 0;
}

int test_random_quadratic_h4_proven_pair() {
    const slong quadratic[] = {1, -3};
    const char* pair_name = "random quadratic H4 class/unit proven";

    FieldSetup setup = setup_from_coefficients(quadratic, 2);
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_matrix_options(options, factor_base_bound,
                                  setup.maximal_order,
                                  silex::CertificationMode::proven)) {
        std::cerr << pair_name << ": options unavailable\n";
        return 1;
    }

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    if (!units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128)) {
        std::cerr << pair_name << ": class/unit computation failed\n";
        return 1;
    }

    sflint::Fmpz class_order;
    if (!class_group.order(sflint::FmpzRef(class_order)) ||
        !sflint::fmpz_is_one(sflint::FmpzConstRef(class_order))) {
        std::cerr << pair_name << ": unexpected class order\n";
        return 1;
    }
    if (!units.is_set() || units.free_rank() != 1) {
        std::cerr << pair_name << ": unexpected unit rank\n";
        return 1;
    }
    if (class_group.certification_status() !=
        silex::CertificationMode::proven) {
        std::cerr << pair_name << ": class group not proven\n";
        return 1;
    }
    if (units.certification_status() != silex::CertificationMode::proven) {
        std::cerr << pair_name << ": unit group not proven\n";
        return 1;
    }
    if (class_group.unit_proof_status() != silex::ProofState::verified ||
        class_group.regulator_proof_status() != silex::ProofState::verified) {
        std::cerr << pair_name << ": missing verified unit proof metadata\n";
        return 1;
    }

    return 0;
}

}  // namespace

int main() {
    int status = 0;
    status |= test_public_paired_certification_request_validation();
    status |=
            test_public_paired_cross_order_failure_preserves_existing_outputs();
    status |= test_stable_class_group_rows();
    status |= test_release_target_proven_transactions();
    status |= test_exact_edge_transactions();
    status |= test_higher_degree_completion_boundaries();
    status |= test_expanded_random_class_group_candidates();
    status |= test_random_quadratic_h4_proven_pair();
    return status;
}
