#include "class_group_internal.hpp"

#include <silex/ideal_factorization.hpp>

#include "lll_relation_search_internal.hpp"
#include "ideal_lattice_reduction_internal.hpp"
#include "relation_candidate_internal.hpp"

namespace silex::detail {
namespace {

enum class FactorAttempt {
    nonsmooth,
    success,
    failure,
};

enum class SearchAttempt {
    no_match,
    success,
    failure,
};

bool valid_witness_inputs(
        const ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzMatConstRef factor_base_row,
        const ClassGroupIdealRelationWitnessOptions& options) noexcept {
    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    return order != nullptr && order->parent() != nullptr &&
           ideal.parent() != nullptr &&
           ideal.parent()->has_same_data(*order) && base != nullptr &&
           base->parent() != nullptr && base->parent()->has_same_data(*order) &&
           flint::fmpz_mat_nrows(factor_base_row) == 1 &&
           flint::fmpz_mat_ncols(factor_base_row) == base->length() &&
           options.reduction_precision > 0 && options.max_candidates >= 0 &&
           options.max_candidates_per_ideal > 0 &&
           options.max_random_products >= 0;
}

void set_failure(ClassGroupIdealRelationWitnessResult& result,
                 ClassGroupIdealRelationWitnessStatus status,
                 ClassGroupIdealRelationWitnessStage stage) noexcept {
    result.status = status;
    result.stage = stage;
}

FactorAttempt factor_if_smooth(
        flint::FmpzMat& row,
        const Ideal& ideal,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics) noexcept {
    bool smooth = false;
    if (!ideal_is_smooth(smooth, ideal, base)) {
        return FactorAttempt::failure;
    }
    if (!smooth) {
        return FactorAttempt::nonsmooth;
    }
    return ideal_factor_over_base(flint::FmpzMatRef(row), ideal, base,
                                  diagnostics)
            ? FactorAttempt::success
            : FactorAttempt::failure;
}

bool reconstruct_factor_base_ideal(
        FractionalIdeal& out,
        const FactorBase& base,
        flint::FmpzMatConstRef row) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*order) ||
        flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != base.length()) {
        return false;
    }

    FractionalIdeal accumulator(*order);
    FractionalIdeal prime_ideal(*order);
    FractionalIdeal power(*order);
    PrimeIdeal prime(*order);
    if (!accumulator.is_defined() || !prime_ideal.is_defined() ||
        !power.is_defined() || !prime.is_defined() || !accumulator.one()) {
        return false;
    }

    for (slong i = 0; i < base.length(); ++i) {
        const flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(row, 0, i);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!base.prime(prime, i) ||
            !prime_to_fractional_ideal(prime_ideal, prime) ||
            !power.pow_fmpz(prime_ideal, exponent) ||
            !accumulator.multiply(accumulator, power)) {
            return false;
        }
    }
    out.swap(accumulator);
    return true;
}

bool publish_verified_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& original_ideal,
        FactoredElement& candidate_multiplier,
        flint::FmpzMat& candidate_row,
        ClassGroupIdealRelationWitnessStage success_stage) noexcept {
    if (!verify_class_group_ideal_relation_witness(
                context, original_ideal, candidate_multiplier,
                flint::FmpzMatConstRef(candidate_row))) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::verification_failure,
                    ClassGroupIdealRelationWitnessStage::exact_verification);
        return false;
    }

    multiplier.swap(candidate_multiplier);
    flint::fmpz_mat_set(factor_base_row,
                        flint::FmpzMatConstRef(candidate_row));
    result.status = ClassGroupIdealRelationWitnessStatus::success;
    result.stage = success_stage;
    return true;
}

SearchAttempt search_ideal_for_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& original_ideal,
        const Ideal& search_ideal,
        const Element& reduction_multiplier,
        flint::FmpzMatConstRef random_product_row,
        bool random_product,
        const ClassGroupIdealRelationWitnessOptions& options) noexcept {
    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    const DiagnosticsContext* diagnostics = context.diagnostics();
    relation_search::IdealLatticeLllData t2_context;
    relation_search::OrderCoordinateElementConversion conversion;
    if (!relation_search::build_ideal_lattice_lll_data(
                t2_context, search_ideal, options.reduction_precision) ||
        !conversion.reset(*order)) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_setup);
        return SearchAttempt::failure;
    }

    std::vector<slong> selected;
    if (!relation_search::select_small_lll_rows(
                selected, t2_context, search_ideal, *order)) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_setup);
        return SearchAttempt::failure;
    }

    flint::FmpzMat coordinates(1, order->degree());
    flint::Fmpz one;
    flint::fmpz_one(flint::FmpzRef(one));
    slong candidates_on_ideal = 0;
    slong count = 0;
    while (result.candidates_tried < options.max_candidates &&
           candidates_on_ideal < options.max_candidates_per_ideal &&
           relation_search::set_small_lll_next_candidate(
                   coordinates, count, t2_context, selected)) {
        Element alpha(*order->parent());
        ++result.candidates_tried;
        ++candidates_on_ideal;
        if (!alpha.is_defined() ||
            !conversion.set(alpha, flint::FmpzMatConstRef(coordinates),
                            flint::FmpzConstRef(one))) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_enumeration);
            return SearchAttempt::failure;
        }

        FractionalIdeal alpha_ideal(*order);
        FractionalIdeal search_fractional(*order);
        FractionalIdeal quotient(*order);
        Ideal quotient_integral(*order);
        flint::Fmpz quotient_denominator;
        if (!alpha_ideal.is_defined() || !search_fractional.is_defined() ||
            !quotient.is_defined() || !quotient_integral.is_defined() ||
            !alpha_ideal.set_principal(alpha, diagnostics) ||
            !search_fractional.set_integral(search_ideal) ||
            !quotient.colon(alpha_ideal, search_fractional) ||
            !quotient.get_integral_den(
                    quotient_integral,
                    flint::FmpzRef(quotient_denominator)) ||
            !flint::fmpz_is_one(quotient_denominator)) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_factorization);
            return SearchAttempt::failure;
        }

        flint::FmpzMat quotient_row(1, base->length());
        const FactorAttempt factor_attempt = factor_if_smooth(
                quotient_row, quotient_integral, *base, diagnostics);
        if (factor_attempt == FactorAttempt::nonsmooth) {
            continue;
        }
        if (factor_attempt == FactorAttempt::failure) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::factorization_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_factorization);
            return SearchAttempt::failure;
        }

        flint::FmpzMat candidate_row(1, base->length());
        for (slong i = 0; i < base->length(); ++i) {
            flint::FmpzRef entry = flint::fmpz_mat_entry(
                    flint::FmpzMatRef(candidate_row), 0, i);
            flint::fmpz_add(
                    entry,
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(quotient_row), 0, i),
                    flint::fmpz_mat_entry(random_product_row, 0, i));
            flint::fmpz_neg(entry, flint::FmpzConstRef(entry.raw()));
        }

        FactoredElement candidate_multiplier(*order->parent());
        if (!candidate_multiplier.is_defined() ||
            !candidate_multiplier.one() ||
            !candidate_multiplier.push(reduction_multiplier, 1) ||
            !candidate_multiplier.push(alpha, -1)) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::continuation_enumeration);
            return SearchAttempt::failure;
        }
        candidate_multiplier.normalize();
        result.used_random_product = random_product;
        const ClassGroupIdealRelationWitnessStage success_stage =
                random_product
                ? ClassGroupIdealRelationWitnessStage::random_product_enumeration
                : ClassGroupIdealRelationWitnessStage::continuation_enumeration;
        return publish_verified_witness(
                       result, multiplier, factor_base_row, context,
                       original_ideal, candidate_multiplier, candidate_row,
                       success_stage)
                ? SearchAttempt::success
                : SearchAttempt::failure;
    }
    return SearchAttempt::no_match;
}

}  // namespace

const char* class_group_ideal_relation_witness_status_name(
        ClassGroupIdealRelationWitnessStatus status) noexcept {
    switch (status) {
        case ClassGroupIdealRelationWitnessStatus::not_started:
            return "not_started";
        case ClassGroupIdealRelationWitnessStatus::success:
            return "success";
        case ClassGroupIdealRelationWitnessStatus::invalid_input:
            return "invalid_input";
        case ClassGroupIdealRelationWitnessStatus::arithmetic_failure:
            return "arithmetic_failure";
        case ClassGroupIdealRelationWitnessStatus::factorization_failure:
            return "factorization_failure";
        case ClassGroupIdealRelationWitnessStatus::verification_failure:
            return "verification_failure";
        case ClassGroupIdealRelationWitnessStatus::exhausted:
            return "exhausted";
    }
    return "unknown";
}

const char* class_group_ideal_relation_witness_stage_name(
        ClassGroupIdealRelationWitnessStage stage) noexcept {
    switch (stage) {
        case ClassGroupIdealRelationWitnessStage::none:
            return "none";
        case ClassGroupIdealRelationWitnessStage::input_validation:
            return "input_validation";
        case ClassGroupIdealRelationWitnessStage::direct_factorization:
            return "direct_factorization";
        case ClassGroupIdealRelationWitnessStage::ideal_reduction:
            return "ideal_reduction";
        case ClassGroupIdealRelationWitnessStage::reduced_factorization:
            return "reduced_factorization";
        case ClassGroupIdealRelationWitnessStage::continuation_setup:
            return "continuation_setup";
        case ClassGroupIdealRelationWitnessStage::continuation_enumeration:
            return "continuation_enumeration";
        case ClassGroupIdealRelationWitnessStage::continuation_factorization:
            return "continuation_factorization";
        case ClassGroupIdealRelationWitnessStage::random_product_construction:
            return "random_product_construction";
        case ClassGroupIdealRelationWitnessStage::random_product_enumeration:
            return "random_product_enumeration";
        case ClassGroupIdealRelationWitnessStage::exact_verification:
            return "exact_verification";
        case ClassGroupIdealRelationWitnessStage::search_exhaustion:
            return "search_exhaustion";
    }
    return "unknown";
}

bool verify_class_group_ideal_relation_witness(
        const ClassGroupContext& context,
        const Ideal& ideal,
        const FactoredElement& multiplier,
        flint::FmpzMatConstRef factor_base_row) noexcept {
    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    if (order == nullptr || order->parent() == nullptr || base == nullptr ||
        ideal.parent() == nullptr ||
        !ideal.parent()->has_same_data(*order) || multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*order->parent()) ||
        flint::fmpz_mat_nrows(factor_base_row) != 1 ||
        flint::fmpz_mat_ncols(factor_base_row) != base->length()) {
        return false;
    }

    Element expanded_multiplier(*order->parent());
    FractionalIdeal multiplier_ideal(*order);
    FractionalIdeal input_ideal(*order);
    FractionalIdeal left(*order);
    FractionalIdeal right(*order);
    return expanded_multiplier.is_defined() && multiplier_ideal.is_defined() &&
           input_ideal.is_defined() && left.is_defined() && right.is_defined() &&
           multiplier.evaluate(expanded_multiplier) &&
           multiplier_ideal.set_principal(expanded_multiplier,
                                          context.diagnostics()) &&
           input_ideal.set_integral(ideal) &&
           left.multiply(multiplier_ideal, input_ideal) &&
           factor_base_row_ideal(right, *base, factor_base_row) &&
           left.equal(right);
}

bool factor_base_row_ideal(
        FractionalIdeal& out,
        const FactorBase& base,
        flint::FmpzMatConstRef factor_base_row) noexcept {
    return reconstruct_factor_base_ideal(out, base, factor_base_row);
}

bool class_group_ideal_relation_continuation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& original_ideal,
        const Ideal& reduced_ideal,
        const Element& reduction_multiplier,
        const ClassGroupIdealRelationWitnessOptions& options) noexcept {
    result = {};
    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    if (!valid_witness_inputs(
                context, original_ideal,
                flint::FmpzMatConstRef(factor_base_row.raw()),
                              options) ||
        reduced_ideal.parent() == nullptr ||
        !reduced_ideal.parent()->has_same_data(*order) ||
        reduction_multiplier.parent() == nullptr ||
        !reduction_multiplier.parent()->has_same_data(*order->parent())) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::invalid_input,
                    ClassGroupIdealRelationWitnessStage::input_validation);
        return false;
    }
    result.used_reduction = true;

    if (options.max_candidates == 0) {
        set_failure(result, ClassGroupIdealRelationWitnessStatus::exhausted,
                    ClassGroupIdealRelationWitnessStage::search_exhaustion);
        return false;
    }

    flint::FmpzMat random_product_row(1, base->length());
    flint::fmpz_mat_zero(flint::FmpzMatRef(random_product_row));
    SearchAttempt attempt = search_ideal_for_witness(
            result, multiplier, factor_base_row, context, original_ideal,
            reduced_ideal, reduction_multiplier,
            flint::FmpzMatConstRef(random_product_row), false, options);
    if (attempt == SearchAttempt::success) {
        return true;
    }
    if (attempt == SearchAttempt::failure) {
        return false;
    }

    Ideal random_product(*order);
    Ideal search_ideal(*order);
    if (!random_product.is_defined() || !search_ideal.is_defined()) {
        set_failure(
                result,
                ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                ClassGroupIdealRelationWitnessStage::random_product_construction);
        return false;
    }
    for (slong i = 0;
         i < options.max_random_products &&
         result.candidates_tried < options.max_candidates;
         ++i) {
        if (base->length() <= 0) {
            break;
        }
        if (!relation_search::build_random_factor_base_product(
                    random_product, flint::FmpzMatRef(random_product_row),
                    *base, i, options.random_seed, context.diagnostics()) ||
            !search_ideal.multiply(reduced_ideal, random_product)) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::random_product_construction);
            return false;
        }
        ++result.random_products_tried;
        attempt = search_ideal_for_witness(
                result, multiplier, factor_base_row, context, original_ideal,
                search_ideal, reduction_multiplier,
                flint::FmpzMatConstRef(random_product_row), true, options);
        if (attempt == SearchAttempt::success) {
            return true;
        }
        if (attempt == SearchAttempt::failure) {
            return false;
        }
    }

    set_failure(result, ClassGroupIdealRelationWitnessStatus::exhausted,
                ClassGroupIdealRelationWitnessStage::search_exhaustion);
    return false;
}

bool class_group_ideal_relation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& ideal,
        const ClassGroupIdealRelationWitnessOptions& options) noexcept {
    result = {};
    if (!valid_witness_inputs(
                context, ideal,
                flint::FmpzMatConstRef(factor_base_row.raw()), options)) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::invalid_input,
                    ClassGroupIdealRelationWitnessStage::input_validation);
        return false;
    }

    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    flint::FmpzMat direct_row(1, base->length());
    const FactorAttempt direct = factor_if_smooth(
            direct_row, ideal, *base, context.diagnostics());
    if (direct == FactorAttempt::failure) {
        set_failure(
                result,
                ClassGroupIdealRelationWitnessStatus::factorization_failure,
                ClassGroupIdealRelationWitnessStage::direct_factorization);
        return false;
    }
    if (direct == FactorAttempt::success) {
        FactoredElement one(*order->parent());
        if (!one.is_defined() || !one.one()) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::direct_factorization);
            return false;
        }
        return publish_verified_witness(
                result, multiplier, factor_base_row, context, ideal, one,
                direct_row,
                ClassGroupIdealRelationWitnessStage::direct_factorization);
    }

    Ideal reduced(*order);
    Element reduction_multiplier(*order->parent());
    if (!reduced.is_defined() || !reduction_multiplier.is_defined() ||
        !reduce_ideal_lattice(
                reduced, reduction_multiplier, ideal,
                options.reduction_precision, context.diagnostics())) {
        set_failure(result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::ideal_reduction);
        return false;
    }
    result.used_reduction = true;

    flint::FmpzMat reduced_row(1, base->length());
    const FactorAttempt reduced_factor = factor_if_smooth(
            reduced_row, reduced, *base, context.diagnostics());
    if (reduced_factor == FactorAttempt::failure) {
        set_failure(
                result,
                ClassGroupIdealRelationWitnessStatus::factorization_failure,
                ClassGroupIdealRelationWitnessStage::reduced_factorization);
        return false;
    }
    if (reduced_factor == FactorAttempt::success) {
        FactoredElement reduced_multiplier(*order->parent());
        if (!reduced_multiplier.is_defined() ||
            !reduced_multiplier.set_element(reduction_multiplier)) {
            set_failure(
                    result,
                    ClassGroupIdealRelationWitnessStatus::arithmetic_failure,
                    ClassGroupIdealRelationWitnessStage::reduced_factorization);
            return false;
        }
        return publish_verified_witness(
                result, multiplier, factor_base_row, context, ideal,
                reduced_multiplier, reduced_row,
                ClassGroupIdealRelationWitnessStage::reduced_factorization);
    }

    return class_group_ideal_relation_continuation_witness(
            result, multiplier, factor_base_row, context, ideal, reduced,
            reduction_multiplier, options);
}

}  // namespace silex::detail
