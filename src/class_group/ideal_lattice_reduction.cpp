#include <silex/class_group.hpp>

#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/order_element.hpp>

#include "ideal_lattice_reduction_internal.hpp"
#include "ideal_lattice_lll_internal.hpp"

#include <limits>

namespace silex {

using namespace detail::relation_search;

namespace {

bool element_scalar_div_fmpz(
        Element& out,
        const Element& input,
        flint::FmpzConstRef denominator) noexcept {
    if (!out.has_same_parent(input) ||
        flint::fmpz_sgn(denominator) <= 0) {
        return false;
    }

    flint::FmpqPoly polynomial;
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }

    flint::fmpq_poly_scalar_div_fmpz(
            polynomial, polynomial, denominator);
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool principal_times_integral_ideal(
        Ideal& out,
        const Element& multiplier,
        const Ideal& ideal,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*order) ||
        !multiplier.has_parent(*field) || !ideal.has_hnf()) {
        return false;
    }

    return detail::multiply_integral_ideal_by_element(out, ideal, multiplier,
                                                      diagnostics);
}

bool abs_order_discriminant(
        flint::Fmpz& out,
        const Order& order) noexcept {
    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }

    flint::fmpz_abs(flint::FmpzRef(out),
                    flint::FmpzConstRef(discriminant));
    return !flint::fmpz_is_zero(flint::FmpzConstRef(out));
}

bool ideal_norm_gt(
        bool& out,
        const Ideal& ideal,
        flint::FmpzConstRef bound) noexcept {
    flint::Fmpz norm;
    if (!ideal.norm(flint::FmpzRef(norm))) {
        return false;
    }

    out = flint::fmpz_cmp(flint::FmpzConstRef(norm), bound) > 0;
    return true;
}

bool ideal_norm_square_gt(
        bool& out,
        const Ideal& ideal,
        flint::FmpzConstRef bound) noexcept {
    flint::Fmpz norm;
    flint::Fmpz square;
    if (!ideal.norm(flint::FmpzRef(norm))) {
        return false;
    }

    flint::fmpz_mul(flint::FmpzRef(square),
                    flint::FmpzConstRef(norm),
                    flint::FmpzConstRef(norm));
    out = flint::fmpz_cmp(flint::FmpzConstRef(square), bound) > 0;
    return true;
}

bool ideal_norm_product_gt(
        bool& out,
        const Ideal& left,
        const Ideal& right,
        flint::FmpzConstRef bound) noexcept {
    flint::Fmpz left_norm;
    flint::Fmpz right_norm;
    flint::Fmpz product;
    if (!left.norm(flint::FmpzRef(left_norm)) ||
        !right.norm(flint::FmpzRef(right_norm))) {
        return false;
    }

    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(left_norm),
                    flint::FmpzConstRef(right_norm));
    out = flint::fmpz_cmp(flint::FmpzConstRef(product), bound) > 0;
    return true;
}

bool factored_multiply_in_place(
        FactoredElement& target,
        const FactoredElement& factor) noexcept {
    const NumberField* field = target.parent();
    if (field == nullptr || factor.parent() == nullptr ||
        !field->has_same_data(*factor.parent())) {
        return false;
    }

    FactoredElement product(*field);
    if (!product.multiply(target, factor)) {
        return false;
    }

    target.swap(product);
    return true;
}

bool factored_push_field_integer(
        FactoredElement& target,
        const NumberField& field,
        flint::FmpzConstRef value,
        slong exponent) noexcept {
    if (target.parent() == nullptr ||
        !target.parent()->has_same_data(field)) {
        return false;
    }
    if (exponent == 0) {
        return true;
    }

    Element scalar(field);
    return scalar.is_defined() && scalar.set_fmpz(value) &&
           target.push(scalar, exponent);
}

}  // namespace

namespace detail {

bool ideal_lattice_short_element(
        Element& out,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.short_element.unweighted");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !out.has_parent(*field) || !ideal.has_hnf()) {
        return false;
    }

    IdealLatticeLllData context;
    if (!build_ideal_lattice_lll_data(context, ideal, precision)) {
        return false;
    }

    const slong degree = order->degree();
    if (degree <= 0 || flint::fmpz_mat_nrows(context.basis) != degree ||
        flint::fmpz_mat_ncols(context.basis) != degree) {
        return false;
    }

    flint::FmpzMat coordinates(1, degree);
    OrderElement order_element(*order);
    Element candidate(*field);
    if (!order_element.is_defined() || !candidate.is_defined() ||
        !copy_ideal_basis_row_coordinates(
                coordinates, flint::FmpzMatConstRef(context.basis), 0) ||
        !order_element.set_coordinates(flint::FmpzMatConstRef(coordinates)) ||
        !order_element.get_element(candidate)) {
        return false;
    }

    return out.set(candidate);
}

bool weighted_ideal_lattice_short_element(
        Element& out,
        const Ideal& ideal,
        flint::FmpzMatConstRef weights,
        slong precision,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* cache) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.short_element.weighted");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !out.has_parent(*field) || !ideal.has_hnf() ||
        precision <= 0) {
        return false;
    }

    const slong degree = order->degree();
    if (!ideal_lattice_weight_vector_shape_is_valid(weights, degree)) {
        return false;
    }
    if (ideal_lattice_weight_vector_is_zero(weights)) {
        return ideal_lattice_short_element(out, ideal, precision,
                                           diagnostics);
    }

    flint::FmpzMat basis(degree, degree);
    flint::FmpzMat transform(degree, degree);
    flint::FmpzMat scaled_gram(degree, degree);
    flint::Fmpz gram_denominator;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.short_element.weighted.reduction");
        if (!build_weighted_ideal_lattice_reduction(
                    basis, transform, scaled_gram, gram_denominator, ideal,
                    weights, precision, diagnostics, cache)) {
            return false;
        }
    }

    flint::FmpzMat coordinates(1, degree);
    OrderElement order_element(*order);
    Element candidate(*field);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.short_element.weighted.extract");
        if (!order_element.is_defined() || !candidate.is_defined() ||
            !copy_ideal_basis_row_coordinates(
                    coordinates, flint::FmpzMatConstRef(basis), 0) ||
            !order_element.set_coordinates(
                    flint::FmpzMatConstRef(coordinates)) ||
            !order_element.get_element(candidate)) {
            return false;
        }
    }

    return out.set(candidate);
}

bool weighted_ideal_lattice_short_element(
        Element& out,
        const FractionalIdeal& ideal,
        flint::FmpzMatConstRef weights,
        slong precision,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* cache) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.short_element.weighted.fractional");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !out.has_parent(*field) || precision <= 0) {
        return false;
    }

    Ideal numerator(*order);
    flint::Fmpz denominator;
    Element numerator_short(*field);
    if (!numerator.is_defined() || !numerator_short.is_defined()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.short_element.weighted.fractional.integral_den");
        if (!ideal.get_integral_den(numerator, flint::FmpzRef(denominator))) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.short_element.weighted.fractional.numerator");
        if (!weighted_ideal_lattice_short_element(
                    numerator_short, numerator, weights, precision,
                    diagnostics, cache)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.ideal_lattice_reduction.short_element.weighted.fractional.divide");
        if (!element_scalar_div_fmpz(out, numerator_short,
                                     flint::FmpzConstRef(denominator))) {
            return false;
        }
    }
    return true;
}

bool reduce_ideal_lattice(
        Ideal& reduced,
        Element& multiplier,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.reduce.single");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || reduced.parent() == nullptr ||
        !reduced.parent()->has_same_data(*order) ||
        !multiplier.has_parent(*field) || !ideal.has_hnf() ||
        precision <= 0) {
        return false;
    }

    FractionalIdeal integral(*order);
    FractionalIdeal inverse(*order);
    Ideal inverse_numerator(*order);
    flint::Fmpz inverse_denominator;
    Element short_numerator(*field);
    Element source_multiplier(*field);
    Ideal candidate(*order);
    if (!integral.is_defined() || !inverse.is_defined() ||
        !inverse_numerator.is_defined() || !short_numerator.is_defined() ||
        !source_multiplier.is_defined() || !candidate.is_defined() ||
        !integral.set_integral(ideal) || !inverse.invert(integral) ||
        !inverse.get_integral_den(inverse_numerator,
                                  flint::FmpzRef(inverse_denominator)) ||
        !ideal_lattice_short_element(short_numerator, inverse_numerator,
                                     precision, diagnostics) ||
        !element_scalar_div_fmpz(source_multiplier, short_numerator,
                                 flint::FmpzConstRef(inverse_denominator)) ||
        !principal_times_integral_ideal(candidate, source_multiplier, ideal,
                                        diagnostics)) {
        return false;
    }

    return reduced.set(candidate) && multiplier.set(source_multiplier);
}

bool reduce_ideal_product(
        Ideal& reduced,
        Element& multiplier,
        const Ideal& left,
        const Ideal& right,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.reduce.product");
    const Order* order = left.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || right.parent() == nullptr ||
        reduced.parent() == nullptr || !right.parent()->has_same_data(*order) ||
        !reduced.parent()->has_same_data(*order) ||
        !multiplier.has_parent(*field) || !left.has_hnf() ||
        !right.has_hnf()) {
        return false;
    }

    Ideal product(*order);
    if (!product.is_defined() || !product.multiply(left, right)) {
        return false;
    }

    // Preserve exact product materialization before ordinary ideal reduction.
    // A specialized product-basis fast path remains a deferred optimization.
    return reduce_ideal_lattice(reduced, multiplier, product, precision,
                                diagnostics);
}

bool reduce_ideal_signed_power(
        Ideal& reduced,
        FactoredElement& multiplier,
        const Ideal& ideal,
        slong exponent,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_lattice_reduction.reduce.signed_power");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || reduced.parent() == nullptr ||
        !reduced.parent()->has_same_data(*order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) || !ideal.has_hnf() ||
        precision <= 0 ||
        exponent == std::numeric_limits<slong>::min()) {
        return false;
    }

    if (exponent == 0) {
        Ideal one(*order);
        FactoredElement one_multiplier(*field);
        return one.is_defined() && one_multiplier.is_defined() &&
               one.one() && one_multiplier.one() &&
               reduced.set(one) && multiplier.set(one_multiplier);
    }

    flint::Fmpz abs_discriminant;
    if (!abs_order_discriminant(abs_discriminant, *order)) {
        return false;
    }

    Ideal active(*order);
    FactoredElement active_multiplier(*field);
    if (!active.is_defined() || !active_multiplier.is_defined() ||
        !active.set(ideal) || !active_multiplier.one()) {
        return false;
    }

    slong active_exponent = exponent;
    bool norm_exceeds_discriminant = false;
    if (!ideal_norm_gt(norm_exceeds_discriminant, active,
                       flint::FmpzConstRef(abs_discriminant))) {
        return false;
    }
    if (norm_exceeds_discriminant) {
        Ideal reduced_active(*order);
        Element reduction_multiplier(*field);
        if (!reduced_active.is_defined() ||
            !reduction_multiplier.is_defined() ||
            !reduce_ideal_lattice(reduced_active, reduction_multiplier,
                                active, precision, diagnostics) ||
            !active.set(reduced_active) ||
            !active_multiplier.push(reduction_multiplier,
                                    -active_exponent)) {
            return false;
        }
    }

    if (active_exponent < 0) {
        FractionalIdeal fractional_active(*order);
        FractionalIdeal inverse_active(*order);
        Ideal inverse_numerator(*order);
        flint::Fmpz inverse_denominator;
        if (!fractional_active.is_defined() ||
            !inverse_active.is_defined() || !inverse_numerator.is_defined() ||
            !fractional_active.set_integral(active) ||
            !inverse_active.invert(fractional_active) ||
            !inverse_active.get_integral_den(
                    inverse_numerator, flint::FmpzRef(inverse_denominator)) ||
            !active.set(inverse_numerator) ||
            !factored_push_field_integer(
                    active_multiplier, *field,
                    flint::FmpzConstRef(inverse_denominator),
                    active_exponent)) {
            return false;
        }
        active_exponent = -active_exponent;
    }

    if (active_exponent == 1) {
        active_multiplier.normalize();
        return reduced.set(active) && multiplier.set(active_multiplier);
    }

    Ideal half_power(*order);
    FactoredElement half_multiplier(*field);
    if (!half_power.is_defined() || !half_multiplier.is_defined() ||
        !reduce_ideal_signed_power(half_power, half_multiplier, active,
                            active_exponent / 2, precision, diagnostics)) {
        return false;
    }

    FactoredElement half_multiplier_square(*field);
    if (!half_multiplier_square.is_defined() ||
        !half_multiplier_square.pow_si(half_multiplier, 2) ||
        !factored_multiply_in_place(active_multiplier,
                                    half_multiplier_square)) {
        return false;
    }

    Ideal squared(*order);
    if (!squared.is_defined()) {
        return false;
    }
    bool square_norm_exceeds_discriminant = false;
    if (!ideal_norm_square_gt(square_norm_exceeds_discriminant, half_power,
                              flint::FmpzConstRef(abs_discriminant))) {
        return false;
    }
    if (square_norm_exceeds_discriminant) {
        Element reduction_multiplier(*field);
        if (!reduction_multiplier.is_defined() ||
            !reduce_ideal_product(squared, reduction_multiplier, half_power,
                                  half_power, precision, diagnostics) ||
            !active_multiplier.push(reduction_multiplier, -1)) {
            return false;
        }
    } else if (!squared.multiply(half_power, half_power)) {
        return false;
    }

    Ideal result(*order);
    if (!result.is_defined()) {
        return false;
    }
    if (active_exponent % 2 != 0) {
        bool product_norm_exceeds_discriminant = false;
        if (!ideal_norm_product_gt(product_norm_exceeds_discriminant,
                                   squared, active,
                                   flint::FmpzConstRef(abs_discriminant))) {
            return false;
        }
        if (product_norm_exceeds_discriminant) {
            Element reduction_multiplier(*field);
            if (!reduction_multiplier.is_defined() ||
                !reduce_ideal_product(result, reduction_multiplier, squared,
                                      active, precision, diagnostics) ||
                !active_multiplier.push(reduction_multiplier, -1)) {
                return false;
            }
        } else if (!result.multiply(squared, active)) {
            return false;
        }
    } else if (!result.set(squared)) {
        return false;
    }

    active_multiplier.normalize();
    return reduced.set(result) && multiplier.set(active_multiplier);
}

}  // namespace detail

}  // namespace silex
