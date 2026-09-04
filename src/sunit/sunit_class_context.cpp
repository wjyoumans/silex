#include "sunit_internal.hpp"

namespace silex::detail {
namespace {

void fail(SUnitClassBuildResult& result,
          SUnitClassBuildStage stage,
          slong selected_index = -1) noexcept {
    result.success = false;
    result.stage = stage;
    result.selected_index = selected_index;
}

bool copy_row(flint::FmpzMatRef out,
              slong out_row,
              flint::FmpzMatConstRef in,
              slong in_row) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    for (slong col = 0; col < flint::fmpz_mat_ncols(in); ++col) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, out_row, col),
                        flint::fmpz_mat_entry(in, in_row, col));
    }
    return true;
}

bool copy_selected_primes(std::vector<PrimeIdeal>& out,
                          const std::vector<PrimeIdeal>& input,
                          const Order& order) noexcept {
    out.clear();
    out.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (!input[i].is_defined() || !input[i].has_prime_data() ||
            input[i].parent() == nullptr ||
            !input[i].parent()->has_same_data(order)) {
            return false;
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (input[i].equal(input[j])) {
                return false;
            }
        }
        PrimeIdeal copy(order);
        if (!copy.is_defined() || !copy.set(input[i])) {
            return false;
        }
        out.push_back(std::move(copy));
    }
    return true;
}

bool selected_prime_row_ideal(
        FractionalIdeal& out,
        const std::vector<PrimeIdeal>& selected_primes,
        flint::FmpzMatConstRef row) noexcept {
    const Order* order = out.parent();
    if (order == nullptr || flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) !=
                static_cast<slong>(selected_primes.size())) {
        return false;
    }

    FractionalIdeal accumulator(*order);
    FractionalIdeal prime_ideal(*order);
    FractionalIdeal power(*order);
    if (!accumulator.is_defined() || !prime_ideal.is_defined() ||
        !power.is_defined() || !accumulator.one()) {
        return false;
    }
    for (slong i = 0; i < static_cast<slong>(selected_primes.size()); ++i) {
        flint::FmpzConstRef exponent = flint::fmpz_mat_entry(row, 0, i);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!prime_to_fractional_ideal(
                    prime_ideal,
                    selected_primes[static_cast<std::size_t>(i)]) ||
            !power.pow_fmpz(prime_ideal, exponent) ||
            !accumulator.multiply(accumulator, power)) {
            return false;
        }
    }
    out.swap(accumulator);
    return true;
}

bool compose_augmented_witness(
        FactoredElement& out,
        const WitnessedClassRelationHnfBasis& class_hnf,
        const std::vector<FactoredElement>& selected_multipliers,
        flint::FmpzMatConstRef coefficients,
        slong row) noexcept {
    if (out.parent() == nullptr || row < 0 ||
        row >= flint::fmpz_mat_nrows(coefficients) ||
        flint::fmpz_mat_ncols(coefficients) !=
                static_cast<slong>(class_hnf.witnesses.size() +
                                   selected_multipliers.size())) {
        return false;
    }

    FactoredElement candidate(*out.parent());
    if (!candidate.is_defined() || !candidate.one()) {
        return false;
    }
    slong col = 0;
    for (const FactoredElement& witness : class_hnf.witnesses) {
        if (!multiply_factored_element_power_fmpz(
                    candidate, witness,
                    flint::fmpz_mat_entry(coefficients, row, col))) {
            return false;
        }
        ++col;
    }
    for (const FactoredElement& multiplier : selected_multipliers) {
        if (!multiply_factored_element_power_fmpz(
                    candidate, multiplier,
                    flint::fmpz_mat_entry(coefficients, row, col))) {
            return false;
        }
        ++col;
    }
    candidate.normalize();
    out.swap(candidate);
    return true;
}

bool factored_principal_ideal(FractionalIdeal& out,
                              const FactoredElement& element,
                              const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = out.parent();
    if (order == nullptr || order->parent() == nullptr ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*order->parent())) {
        return false;
    }
    Element value(*order->parent());
    return value.is_defined() && element.evaluate(value) &&
           out.set_principal(value, diagnostics);
}

bool verify_s_class_invariant_witness(
        const SUnitClassContext& context,
        const FractionalIdeal& invariant_ideal,
        flint::FmpzConstRef invariant,
        const FactoredElement& witness,
        flint::FmpzMatConstRef selected_exponents,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = invariant_ideal.parent();
    if (order == nullptr || !order->has_same_data(context.order)) {
        return false;
    }
    FractionalIdeal ideal_power(*order);
    FractionalIdeal selected_product(*order);
    FractionalIdeal expected(*order);
    FractionalIdeal principal(*order);
    return ideal_power.is_defined() && selected_product.is_defined() &&
           expected.is_defined() && principal.is_defined() &&
           ideal_power.pow_fmpz(invariant_ideal, invariant) &&
           selected_prime_row_ideal(selected_product, context.selected_primes,
                                    selected_exponents) &&
           expected.multiply(ideal_power, selected_product) &&
           factored_principal_ideal(principal, witness, diagnostics) &&
           principal.equal(expected);
}

bool verify_sunit_generator(
        const SUnitClassContext& context,
        const FactoredElement& generator,
        flint::FmpzMatConstRef valuation_row,
        const DiagnosticsContext* diagnostics) noexcept {
    FractionalIdeal expected(context.order);
    FractionalIdeal principal(context.order);
    if (!expected.is_defined() || !principal.is_defined() ||
        !selected_prime_row_ideal(expected, context.selected_primes,
                                  valuation_row) ||
        !factored_principal_ideal(principal, generator, diagnostics) ||
        !principal.equal(expected)) {
        return false;
    }

    for (slong i = 0; i < static_cast<slong>(context.selected_primes.size());
         ++i) {
        flint::FmpzConstRef expected_valuation =
                flint::fmpz_mat_entry(valuation_row, 0, i);
        slong actual = 0;
        if (!flint::fmpz_fits_si(expected_valuation) ||
            !context.selected_primes[static_cast<std::size_t>(i)].valuation(
                    actual, generator, diagnostics) ||
            actual != flint::fmpz_get_si(expected_valuation)) {
            return false;
        }
    }
    return true;
}

bool matrix_is_zero(flint::FmpzMatConstRef matrix) noexcept {
    return ::fmpz_mat_is_zero(matrix.raw()) != 0;
}

bool build_s_class_invariants(
        SUnitClassContext& context,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong invariant_count = context.s_class_group.invariant_count();
    const slong generator_count = context.s_class_group.generator_count();
    const slong relation_count = context.s_class_group.relation_count();
    const slong selected_count =
            static_cast<slong>(context.selected_primes.size());

    flint::FmpzMat generator_rows(invariant_count, generator_count);
    flint::FmpzMat relation_coefficients(invariant_count, relation_count);
    if (!context.s_class_group.invariant_generator_matrix(
                flint::FmpzMatRef(generator_rows)) ||
        !context.s_class_group.invariant_generator_relation_matrix(
                flint::FmpzMatRef(relation_coefficients))) {
        return false;
    }

    flint::FmpzMat relation_product(invariant_count, generator_count);
    flint::FmpzMat powered_generators(invariant_count, generator_count);
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(relation_product),
            flint::FmpzMatConstRef(relation_coefficients),
            flint::FmpzMatConstRef(context.augmented_relations));
    for (slong i = 0; i < invariant_count; ++i) {
        flint::Fmpz invariant;
        if (!context.s_class_group.invariant(flint::FmpzRef(invariant), i)) {
            return false;
        }
        for (slong j = 0; j < generator_count; ++j) {
            flint::fmpz_mul(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(powered_generators), i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(generator_rows), i, j),
                    flint::FmpzConstRef(invariant));
        }
    }
    if (!flint::fmpz_mat_equal(
                flint::FmpzMatConstRef(relation_product),
                flint::FmpzMatConstRef(powered_generators))) {
        return false;
    }

    context.s_class_power_selected_exponents =
            flint::FmpzMat(invariant_count, selected_count);
    context.s_class_invariant_ideals.reserve(
            static_cast<std::size_t>(invariant_count));
    context.s_class_power_witnesses.reserve(
            static_cast<std::size_t>(invariant_count));
    flint::FmpzMat generator_row(1, generator_count);
    flint::FmpzMat selected_row(1, selected_count);
    for (slong i = 0; i < invariant_count; ++i) {
        if (!copy_row(flint::FmpzMatRef(generator_row), 0,
                      flint::FmpzMatConstRef(generator_rows), i)) {
            return false;
        }
        FractionalIdeal invariant_ideal(context.order);
        FactoredElement witness(*context.order.parent());
        flint::Fmpz invariant;
        if (!invariant_ideal.is_defined() || !witness.is_defined() ||
            !factor_base_row_ideal(
                    invariant_ideal, context.factor_base,
                    flint::FmpzMatConstRef(generator_row)) ||
            !context.s_class_group.invariant(flint::FmpzRef(invariant), i) ||
            !compose_augmented_witness(
                    witness, context.class_hnf,
                    context.selected_multipliers,
                    flint::FmpzMatConstRef(relation_coefficients), i)) {
            return false;
        }
        for (slong j = 0; j < selected_count; ++j) {
            flint::FmpzRef exponent = flint::fmpz_mat_entry(
                    flint::FmpzMatRef(selected_row), 0, j);
            flint::fmpz_neg(
                    exponent,
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(relation_coefficients), i,
                            generator_count + j));
            flint::fmpz_set(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(
                                    context.s_class_power_selected_exponents),
                            i, j),
                    flint::FmpzConstRef(exponent.raw()));
        }
        if (!verify_s_class_invariant_witness(
                    context, invariant_ideal,
                    flint::FmpzConstRef(invariant), witness,
                    flint::FmpzMatConstRef(selected_row), diagnostics)) {
            return false;
        }
        context.s_class_invariant_ideals.push_back(
                std::move(invariant_ideal));
        context.s_class_power_witnesses.push_back(std::move(witness));
    }
    return true;
}

bool build_sunit_kernel(SUnitClassContext& context,
                        SUnitClassBuildStage& failure_stage,
                        const DiagnosticsContext* diagnostics) noexcept {
    failure_stage = SUnitClassBuildStage::relation_kernel;
    const slong generator_count = context.s_class_group.generator_count();
    const slong relation_count = context.s_class_group.relation_count();
    const slong selected_count =
            static_cast<slong>(context.selected_primes.size());
    if (context.s_class_group.relation_kernel_count() != selected_count) {
        return false;
    }

    context.relation_kernel =
            flint::FmpzMat(selected_count, relation_count);
    if (!context.s_class_group.relation_kernel_matrix(
                flint::FmpzMatRef(context.relation_kernel))) {
        return false;
    }
    flint::FmpzMat kernel_product(selected_count, generator_count);
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(kernel_product),
            flint::FmpzMatConstRef(context.relation_kernel),
            flint::FmpzMatConstRef(context.augmented_relations));
    if (!matrix_is_zero(flint::FmpzMatConstRef(kernel_product))) {
        return false;
    }

    if (selected_count == 0) {
        context.generator_coefficients = flint::FmpzMat(0, relation_count);
        context.valuation_rows = flint::FmpzMat(0, 0);
        return true;
    }

    failure_stage = SUnitClassBuildStage::valuation_hnf;
    flint::FmpzMat raw_valuations(selected_count, selected_count);
    for (slong i = 0; i < selected_count; ++i) {
        for (slong j = 0; j < selected_count; ++j) {
            flint::fmpz_neg(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(raw_valuations), i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(context.relation_kernel), i,
                            generator_count + j));
        }
    }
    if (flint::fmpz_mat_rank(flint::FmpzMatConstRef(raw_valuations)) !=
        selected_count) {
        return false;
    }

    flint::FmpzMat transform(selected_count, selected_count);
    context.valuation_rows =
            flint::FmpzMat(selected_count, selected_count);
    ::fmpz_mat_hnf_transform(context.valuation_rows.raw(), transform.raw(),
                             raw_valuations.raw());
    flint::FmpzMat transformed_valuations(selected_count, selected_count);
    flint::fmpz_mat_mul(flint::FmpzMatRef(transformed_valuations),
                        flint::FmpzMatConstRef(transform),
                        flint::FmpzMatConstRef(raw_valuations));
    if (!flint::fmpz_mat_equal(
                flint::FmpzMatConstRef(transformed_valuations),
                flint::FmpzMatConstRef(context.valuation_rows))) {
        return false;
    }

    context.generator_coefficients =
            flint::FmpzMat(selected_count, relation_count);
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(context.generator_coefficients),
            flint::FmpzMatConstRef(transform),
            flint::FmpzMatConstRef(context.relation_kernel));
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(kernel_product),
            flint::FmpzMatConstRef(context.generator_coefficients),
            flint::FmpzMatConstRef(context.augmented_relations));
    if (!matrix_is_zero(flint::FmpzMatConstRef(kernel_product))) {
        return false;
    }

    failure_stage = SUnitClassBuildStage::generator_composition;
    context.generators_mod_units.reserve(
            static_cast<std::size_t>(selected_count));
    flint::FmpzMat valuation_row(1, selected_count);
    for (slong i = 0; i < selected_count; ++i) {
        FactoredElement generator(*context.order.parent());
        if (!generator.is_defined() ||
            !compose_augmented_witness(
                    generator, context.class_hnf,
                    context.selected_multipliers,
                    flint::FmpzMatConstRef(context.generator_coefficients), i) ||
            !copy_row(flint::FmpzMatRef(valuation_row), 0,
                      flint::FmpzMatConstRef(context.valuation_rows), i) ||
            !verify_sunit_generator(context, generator,
                                    flint::FmpzMatConstRef(valuation_row),
                                    diagnostics)) {
            return false;
        }
        context.generators_mod_units.push_back(std::move(generator));
    }
    return true;
}

}  // namespace

const char* sunit_class_build_stage_name(SUnitClassBuildStage stage) noexcept {
    switch (stage) {
        case SUnitClassBuildStage::none:
            return "none";
        case SUnitClassBuildStage::input_validation:
            return "input_validation";
        case SUnitClassBuildStage::class_hnf_basis:
            return "class_hnf_basis";
        case SUnitClassBuildStage::selected_ideal_relation:
            return "selected_ideal_relation";
        case SUnitClassBuildStage::s_class_presentation:
            return "s_class_presentation";
        case SUnitClassBuildStage::s_class_invariant_witnesses:
            return "s_class_invariant_witnesses";
        case SUnitClassBuildStage::relation_kernel:
            return "relation_kernel";
        case SUnitClassBuildStage::valuation_hnf:
            return "valuation_hnf";
        case SUnitClassBuildStage::generator_composition:
            return "generator_composition";
    }
    return "unknown";
}

bool build_sunit_class_context(
        SUnitClassBuildResult& result,
        SUnitClassContext& out,
        const ClassGroupContext& class_group,
        const std::vector<PrimeIdeal>& selected_primes,
        const ClassGroupIdealRelationWitnessOptions& witness_options) noexcept {
    result = {};
    const Order* order = class_group.parent();
    const FactorBase* factor_base = class_group.factor_base();
    if (order == nullptr || order->parent() == nullptr ||
        !order->is_maximal() || factor_base == nullptr ||
        !class_group.has_presentation() ||
        class_group.certification_status() != CertificationMode::proven ||
        selected_primes.size() > static_cast<std::size_t>(WORD_MAX)) {
        fail(result, SUnitClassBuildStage::input_validation);
        return false;
    }

    SUnitClassContext candidate;
    candidate.order = *order;
    if (!candidate.factor_base.set(*factor_base) ||
        !copy_selected_primes(candidate.selected_primes, selected_primes,
                              candidate.order)) {
        fail(result, SUnitClassBuildStage::input_validation);
        return false;
    }
    if (!class_relation_witnessed_hnf_basis(candidate.class_hnf,
                                             class_group)) {
        fail(result, SUnitClassBuildStage::class_hnf_basis);
        return false;
    }

    const slong generator_count = class_group.generator_count();
    const slong selected_count =
            static_cast<slong>(candidate.selected_primes.size());
    candidate.selected_relation_rows =
            flint::FmpzMat(selected_count, generator_count);
    candidate.selected_multipliers.reserve(
            static_cast<std::size_t>(selected_count));
    candidate.selected_relation_results.reserve(
            static_cast<std::size_t>(selected_count));
    flint::FmpzMat selected_row(1, generator_count);
    for (slong i = 0; i < selected_count; ++i) {
        Ideal ideal(candidate.order);
        FactoredElement multiplier(*candidate.order.parent());
        ClassGroupIdealRelationWitnessResult witness_result;
        if (!ideal.is_defined() || !multiplier.is_defined() ||
            !candidate.selected_primes[static_cast<std::size_t>(i)].get_ideal(
                    ideal) ||
            !class_group_ideal_relation_witness(
                    witness_result, multiplier,
                    flint::FmpzMatRef(selected_row), class_group, ideal,
                    witness_options) ||
            !copy_row(flint::FmpzMatRef(candidate.selected_relation_rows), i,
                      flint::FmpzMatConstRef(selected_row), 0)) {
            result.selected_relation = witness_result;
            fail(result, SUnitClassBuildStage::selected_ideal_relation, i);
            return false;
        }
        candidate.selected_relation_results.push_back(witness_result);
        candidate.selected_multipliers.push_back(std::move(multiplier));
    }

    candidate.augmented_relations = flint::FmpzMat(
            generator_count + selected_count, generator_count);
    for (slong i = 0; i < generator_count; ++i) {
        if (!copy_row(flint::FmpzMatRef(candidate.augmented_relations), i,
                      flint::FmpzMatConstRef(candidate.class_hnf.rows), i)) {
            fail(result, SUnitClassBuildStage::s_class_presentation);
            return false;
        }
    }
    for (slong i = 0; i < selected_count; ++i) {
        if (!copy_row(
                    flint::FmpzMatRef(candidate.augmented_relations),
                    generator_count + i,
                    flint::FmpzMatConstRef(candidate.selected_relation_rows),
                    i)) {
            fail(result, SUnitClassBuildStage::s_class_presentation);
            return false;
        }
    }
    if (!candidate.s_class_group.set_relation_matrix(
                flint::FmpzMatConstRef(candidate.augmented_relations))) {
        fail(result, SUnitClassBuildStage::s_class_presentation);
        return false;
    }
    if (!build_s_class_invariants(candidate, class_group.diagnostics())) {
        fail(result, SUnitClassBuildStage::s_class_invariant_witnesses);
        return false;
    }
    SUnitClassBuildStage kernel_failure_stage =
            SUnitClassBuildStage::relation_kernel;
    if (!build_sunit_kernel(candidate, kernel_failure_stage,
                            class_group.diagnostics())) {
        fail(result, kernel_failure_stage);
        return false;
    }

    candidate.source_class_certification =
            class_group.certification_status();
    candidate.s_class_proof_status = ProofState::verified;
    candidate.s_unit_mod_units_proof_status = ProofState::verified;
    candidate.defined = true;
    out = std::move(candidate);
    result.success = true;
    result.stage = SUnitClassBuildStage::none;
    return true;
}

}  // namespace silex::detail
