#include <silex/class_group.hpp>

#include <silex/order_element.hpp>
#include <silex/relation.hpp>

#include "relation_admission_cache_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"

#include <utility>

namespace silex {
namespace detail::relation_search {

bool OrderCoordinateElementConversion::reset(
        const Order& next_order) noexcept {
    if (next_order.parent() == nullptr || next_order.degree() <= 0) {
        return false;
    }

    const slong degree = next_order.degree();
    flint::FmpqMat next_basis(degree, degree);
    flint::FmpqMat next_power_row(1, degree);
    if (!next_order.get_basis(flint::FmpqMatRef(next_basis))) {
        return false;
    }

    order = &next_order;
    basis = std::move(next_basis);
    power_row = std::move(next_power_row);
    return true;
}

bool OrderCoordinateElementConversion::set(
        Element& out,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den) noexcept {
    if (order == nullptr || order->parent() == nullptr ||
        out.parent() == nullptr ||
        !out.parent()->has_same_data(*order->parent()) ||
        fmpz_sgn(den.raw()) <= 0 ||
        flint::fmpz_mat_nrows(coordinates) != 1 ||
        flint::fmpz_mat_ncols(coordinates) != order->degree()) {
        return false;
    }

    ::fmpq_mat_mul_r_fmpz_mat(power_row.raw(), coordinates.raw(),
                              basis.raw());
    flint::fmpq_poly_zero(flint::FmpqPolyRef(polynomial));
    for (slong col = 0; col < order->degree(); ++col) {
        flint::fmpq_poly_set_coeff_fmpq(
                flint::FmpqPolyRef(polynomial), col,
                flint::fmpq_mat_entry(
                        flint::FmpqMatConstRef(power_row), 0, col));
    }
    if (fmpz_is_one(den.raw()) == 0) {
        flint::fmpq_poly_scalar_div_fmpz(
                flint::FmpqPolyRef(polynomial),
                flint::FmpqPolyConstRef(polynomial), den);
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

void coordinates_from_lattice_combination(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef coefficients,
        flint::FmpzMatConstRef basis) noexcept {
    flint::fmpz_mat_zero(flint::FmpzMatRef(out));

    for (slong i = 0; i < flint::fmpz_mat_nrows(basis); ++i) {
        for (slong j = 0; j < flint::fmpz_mat_ncols(basis); ++j) {
            flint::fmpz_addmul(flint::fmpz_mat_entry(out, 0, j),
                               flint::fmpz_mat_entry(coefficients, 0, i),
                               flint::fmpz_mat_entry(basis, i, j));
        }
    }
}

bool fmpq_poly_to_monic_integral_fmpz_poly(
        flint::FmpzPoly& out,
        const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree <= 0) {
        return false;
    }

    flint::Fmpq coefficient;
    fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial, degree);
    if (fmpq_is_one(coefficient.raw()) == 0) {
        return false;
    }

    fmpz_poly_zero(out.raw());
    flint::Fmpz integer_coefficient;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coefficient.raw())) == 0) {
            fmpz_poly_zero(out.raw());
            return false;
        }
        fmpz_set(integer_coefficient.raw(), fmpq_numref(coefficient.raw()));
        fmpz_poly_set_coeff_fmpz(out.raw(), i, integer_coefficient.raw());
    }
    return true;
}

void integral_coordinates_to_polynomial(
        flint::FmpzPoly& out,
        flint::FmpzMatConstRef coordinates) noexcept {
    fmpz_poly_zero(out.raw());
    const slong cols = flint::fmpz_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        fmpz_poly_set_coeff_fmpz(
                out.raw(), j,
                flint::fmpz_mat_entry(coordinates, 0, j).raw());
    }
}

bool build_norm_prefilter(NormPrefilter& out,
                          const FactorBase& base,
                          flint::FmpzConstRef factor_base_bound,
                          bool allow_large_prime) noexcept {
    out.valid = false;
    out.allow_large_prime = allow_large_prime;
    out.rejected_count = 0;
    out.has_equation_order_norm_polynomial = false;
    if (!base.is_defined() || base.length() <= 0) {
        return true;
    }

    const Order* order = base.parent();
    if (order == nullptr) {
        return false;
    }

    flint::fmpz_one(flint::FmpzRef(out.rational_prime_product));
    flint::fmpz_mul(flint::FmpzRef(out.large_prime_norm_bound),
                    factor_base_bound, factor_base_bound);

    flint::Fmpz rational_prime;
    for (slong block = 0; block < base.rational_prime_block_count();
         ++block) {
        slong block_length = 0;
        if (!base.rational_prime_block_data(
                    flint::FmpzRef(rational_prime), block_length, block)) {
            return false;
        }
        if (block_length <= 0) {
            continue;
        }
        flint::fmpz_mul(flint::FmpzRef(out.rational_prime_product),
                        flint::FmpzConstRef(out.rational_prime_product),
                        flint::FmpzConstRef(rational_prime));
    }

    const NumberField* field = order->parent();
    const nf_struct* raw_field = field == nullptr
                                         ? nullptr
                                         : field->raw_flint_field();
    out.has_equation_order_norm_polynomial =
            order->is_equation_order() && raw_field != nullptr &&
            fmpq_poly_to_monic_integral_fmpz_poly(out.field_polynomial,
                                                  raw_field->pol);

    out.valid = true;
    return true;
}

bool norm_prefilter_accepts_current_norm(
        bool& accepts,
        NormPrefilter& prefilter,
        const DiagnosticsContext* diagnostics) noexcept {
    accepts = true;
    if (!flint::fmpz_is_one(flint::fmpq_den_ref(
                flint::FmpqConstRef(prefilter.norm_scratch)))) {
        return true;
    }

    flint::fmpz_abs(flint::FmpzRef(prefilter.norm_remainder),
                    flint::fmpq_num_ref(
                            flint::FmpqConstRef(prefilter.norm_scratch)));
    if (flint::fmpz_cmp_ui(
                flint::FmpzConstRef(prefilter.norm_remainder), 1) <= 0) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.factor_base_gcd");
        flint::fmpz_gcd(flint::FmpzRef(prefilter.norm_gcd),
                        flint::FmpzConstRef(prefilter.rational_prime_product),
                        flint::FmpzConstRef(prefilter.norm_remainder));
        while (!flint::fmpz_is_one(flint::FmpzConstRef(prefilter.norm_gcd))) {
            flint::fmpz_divexact(
                    flint::FmpzRef(prefilter.norm_remainder),
                    flint::FmpzConstRef(prefilter.norm_remainder),
                    flint::FmpzConstRef(prefilter.norm_gcd));
            flint::fmpz_gcd(flint::FmpzRef(prefilter.norm_gcd),
                            flint::FmpzConstRef(prefilter.norm_gcd),
                            flint::FmpzConstRef(prefilter.norm_remainder));
        }
    }

    if (flint::fmpz_is_one(flint::FmpzConstRef(prefilter.norm_remainder))) {
        return true;
    }
    if (prefilter.allow_large_prime) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.large_prime_screen");
        if (flint::fmpz_cmp(flint::FmpzConstRef(prefilter.norm_remainder),
                            flint::FmpzConstRef(
                                    prefilter.large_prime_norm_bound)) <= 0 &&
            flint::fmpz_is_prime(
                    flint::FmpzConstRef(prefilter.norm_remainder))) {
            return true;
        }
    }

    accepts = false;
    ++prefilter.rejected_count;
    return true;
}

bool norm_prefilter_accepts(bool& accepts,
                            NormPrefilter* prefilter,
                            const Element& alpha,
                            const DiagnosticsContext* diagnostics) noexcept {
    accepts = true;
    if (prefilter == nullptr || !prefilter->valid) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_prefilter.norm");
        if (!alpha.norm(flint::FmpqRef(prefilter->norm_scratch))) {
            return false;
        }
    }
    return norm_prefilter_accepts_current_norm(accepts, *prefilter,
                                               diagnostics);
}

bool norm_prefilter_accepts_integral_equation_coordinates(
        bool& handled,
        bool& accepts,
        NormPrefilter* prefilter,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den,
        const DiagnosticsContext* diagnostics) noexcept {
    handled = false;
    accepts = true;
    if (prefilter == nullptr || !prefilter->valid ||
        !prefilter->has_equation_order_norm_polynomial ||
        fmpz_is_one(den.raw()) == 0 ||
        flint::fmpz_mat_nrows(coordinates) != 1) {
        return true;
    }
    if (flint::fmpz_mat_ncols(coordinates) !=
        fmpz_poly_degree(prefilter->field_polynomial.raw())) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.coordinate_candidate_prefilter.coordinate_norm");
        integral_coordinates_to_polynomial(prefilter->coordinate_polynomial,
                                           coordinates);
        fmpz_poly_resultant(prefilter->norm_resultant_scratch.raw(),
                            prefilter->field_polynomial.raw(),
                            prefilter->coordinate_polynomial.raw());
        flint::fmpq_set_fmpz(
                flint::FmpqRef(prefilter->norm_scratch),
                flint::FmpzConstRef(prefilter->norm_resultant_scratch));
    }
    handled = true;
    return norm_prefilter_accepts_current_norm(accepts, *prefilter,
                                               diagnostics);
}

bool try_generator_relation(ClassGroupContext& context,
                            const Element& alpha,
                            const fmpq* known_norm,
                            const fmpz_mat_struct* known_integral_coordinates,
                            const flint::FmpzPoly*
                                    known_integral_coordinate_polynomial,
                            ClassGroupRelationSource source,
                            slong& accepted_relations,
                            bool& partial_throttle_exit,
                            bool& relation_appended) noexcept {
    const slong relation_count_before = context.relation_count();
    bool candidate_partial_throttle = false;
    relation_appended = false;
    bool ok = false;
    if (known_norm != nullptr && known_integral_coordinates != nullptr) {
        // reference Rel_add.jl and reference factorgen keep the exact norm/factor data
        // produced during admission.  Keep this prefilter data on the same
        // private path instead of reconstructing it for relation materialization.
        ok = detail::ClassGroupRelationSearchAccess::
                try_append_integral_generator_relation(
                        context, candidate_partial_throttle, alpha,
                        flint::FmpzMatConstRef(known_integral_coordinates),
                        flint::FmpqConstRef(known_norm),
                        known_integral_coordinate_polynomial, source);
    } else if (known_norm != nullptr) {
        ok = context.try_append_generator_relation_with_norm(
                candidate_partial_throttle, alpha,
                flint::FmpqConstRef(known_norm), source);
    } else {
        ok = context.try_append_generator_relation(
                candidate_partial_throttle, alpha, source);
    }
    if (!ok) {
        return false;
    }
    if (context.relation_count() > relation_count_before) {
        ++accepted_relations;
        relation_appended = true;
    }
    if (candidate_partial_throttle &&
        source == ClassGroupRelationSource::Search) {
        partial_throttle_exit = true;
    }
    return true;
}

bool element_from_order_coordinates_den(
        Element& out,
        const Order& order,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef den) noexcept {
    if (out.parent() == nullptr || order.parent() == nullptr ||
        !out.parent()->has_same_data(*order.parent()) ||
        fmpz_sgn(den.raw()) <= 0) {
        return false;
    }

    OrderElement order_element(order);
    Element candidate(*order.parent());
    if (!order_element.is_defined() || !candidate.is_defined() ||
        !order_element.set_coordinates(coordinates) ||
        !order_element.get_element(candidate)) {
        return false;
    }

    if (fmpz_is_one(den.raw()) == 0) {
        flint::FmpqPoly polynomial;
        if (!candidate.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
            return false;
        }
        flint::fmpq_poly_scalar_div_fmpz(
                polynomial, polynomial, den);
        if (!candidate.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial))) {
            return false;
        }
    }

    out.swap(candidate);
    return true;
}

bool element_is_scalar_rational(bool& scalar, const Element& element) noexcept {
    scalar = false;
    flint::FmpqPoly polynomial;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }

    scalar = flint::fmpq_poly_degree(
                     flint::FmpqPolyConstRef(polynomial)) <= 0;
    return true;
}

bool reduced_basis_first_row_is_scalar_rational(
        bool& scalar,
        const Order& order,
        flint::FmpzMatConstRef basis) noexcept {
    scalar = false;
    const NumberField* field = order.parent();
    const slong degree = order.degree();
    if (field == nullptr || degree <= 0 ||
        flint::fmpz_mat_nrows(basis) != degree ||
        flint::fmpz_mat_ncols(basis) != degree) {
        return false;
    }

    flint::FmpzMat row(1, degree);
    for (slong j = 0; j < degree; ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(row, 0, j),
                        flint::fmpz_mat_entry(basis, 0, j));
    }

    flint::Fmpz denominator;
    flint::fmpz_one(flint::FmpzRef(denominator));
    Element element(*field);
    return element.is_defined() &&
           element_from_order_coordinates_den(
                   element, order, flint::FmpzMatConstRef(row),
                   flint::FmpzConstRef(denominator)) &&
           element_is_scalar_rational(scalar, element);
}

bool try_coordinate_candidate_den(ClassGroupContext& context,
                                  const Order& order,
                                  flint::FmpzMat& coordinates,
                                  flint::FmpzConstRef den,
                                  NormPrefilter* norm_prefilter,
                                  ClassGroupRelationSource source,
                                  slong target_relation_kernel_units,
                                  slong max_candidates,
                                  slong max_relations,
                                  slong& candidates_tried,
                                  slong& accepted_relations,
                                  CertificationMode requested_certification,
                                  bool& goal_reached,
                                  bool& partial_throttle_exit,
                                  detail::RelationAdmissionCache* admission_cache,
                                  bool random_relation,
                                  slong* factor_attempts,
                                  slong max_factor_attempts,
                                  bool* factor_attempt_limit_reached,
                                  Relation* admission_scratch,
                                  OrderCoordinateElementConversion*
                                          coordinate_conversion,
                                  bool admission_scratch_base_verified)
        noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate");
    if (candidates_tried >= max_candidates ||
        accepted_relations >= max_relations) {
        return true;
    }

    bool scalar_check_done = false;
    if (admission_cache != nullptr && order.is_equation_order() &&
        flint::fmpz_mat_nrows(flint::FmpzMatConstRef(coordinates)) == 1 &&
        flint::fmpz_mat_ncols(flint::FmpzMatConstRef(coordinates)) ==
                order.degree()) {
        bool scalar = true;
        for (slong col = 1; col < order.degree(); ++col) {
            if (!flint::fmpz_is_zero(
                        flint::fmpz_mat_entry(coordinates, 0, col))) {
                scalar = false;
                break;
            }
        }
        scalar_check_done = true;
        if (scalar) {
            return true;
        }
    }

    Element alpha(*order.parent());
    bool alpha_ready = false;
    auto ensure_alpha = [&]() noexcept -> bool {
        if (alpha_ready) {
            return true;
        }
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_element");
        const bool converted =
                alpha.is_defined() &&
                (coordinate_conversion == nullptr
                         ? element_from_order_coordinates_den(
                                   alpha, order,
                                   flint::FmpzMatConstRef(coordinates), den)
                         : coordinate_conversion->set(
                                   alpha,
                                   flint::FmpzMatConstRef(coordinates), den));
        if (!converted) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate conversion failed");
            return false;
        }
        alpha_ready = true;
        return true;
    };

    if (admission_cache != nullptr && !scalar_check_done) {
        bool scalar = false;
        if (!ensure_alpha() || !element_is_scalar_rational(scalar, alpha)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate scalar check failed");
            return false;
        }
        if (scalar) {
            return true;
        }
    }

    if (admission_cache != nullptr && factor_attempts != nullptr &&
        max_factor_attempts > 0) {
        // The source finite quadratic-form ideal search increments try_factor before
        // factorgen, so candidates rejected by the norm/factorability gate
        // still consume the per-ideal factor-attempt budget.
        if (*factor_attempts >= max_factor_attempts) {
            if (factor_attempt_limit_reached != nullptr) {
                *factor_attempt_limit_reached = true;
            }
            return true;
        }
        ++*factor_attempts;
    }

    bool prefilter_accepts = true;
    const fmpq* known_norm = nullptr;
    const fmpz_mat_struct* known_integral_coordinates = nullptr;
    const flint::FmpzPoly* known_integral_coordinate_polynomial = nullptr;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_prefilter");
        bool norm_prefilter_handled = false;
        if (!norm_prefilter_accepts_integral_equation_coordinates(
                    norm_prefilter_handled, prefilter_accepts,
                    norm_prefilter,
                    flint::FmpzMatConstRef(coordinates), den,
                    diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate coordinate norm prefilter failed");
            return false;
        }
        if (!norm_prefilter_handled &&
            (!ensure_alpha() ||
             !norm_prefilter_accepts(prefilter_accepts, norm_prefilter,
                                     alpha, diagnostics))) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "coordinate candidate norm prefilter failed");
            return false;
        }
        if (norm_prefilter != nullptr && norm_prefilter->valid) {
            known_norm = norm_prefilter->norm_scratch.raw();
            if (flint::fmpz_is_one(den)) {
                known_integral_coordinates =
                        flint::FmpzMatConstRef(coordinates).raw();
            }
            if (norm_prefilter_handled) {
                known_integral_coordinate_polynomial =
                        &norm_prefilter->coordinate_polynomial;
            }
        }
    }
    if (!prefilter_accepts) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                            "class_group.candidate_prefilter_reject");
        return true;
    }
    const bool can_defer_integral_alpha =
            admission_cache != nullptr && known_norm != nullptr &&
            flint::fmpz_is_one(den) && order.is_equation_order() &&
            flint::fmpz_mat_nrows(flint::FmpzMatConstRef(coordinates)) == 1 &&
            flint::fmpz_mat_ncols(flint::FmpzMatConstRef(coordinates)) ==
                    order.degree();
    if (!can_defer_integral_alpha && !ensure_alpha()) {
        return false;
    }

    ++candidates_tried;
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.candidate_relation");
    bool candidate_partial_throttle = false;
    bool relation_appended = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_try_relation");
        if (admission_cache != nullptr) {
            const FactorBase* base = context.factor_base();
            if (base == nullptr) {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "relation admission missing factor base");
                return false;
            }
            Relation local_relation;
            Relation* relation = admission_scratch;
            if (relation == nullptr) {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.admission_scratch_define");
                if (!local_relation.define(*base)) {
                    SILEX_LOG(
                            diagnostics, DiagnosticsModule::class_group,
                            LogLevel::detail,
                            "relation admission scratch allocation failed");
                    return false;
                }
                relation = &local_relation;
            } else if (!relation->is_defined() ||
                       relation->factor_base() == nullptr ||
                       (admission_scratch_base_verified
                                ? relation->length() != base->length()
                                : !relation->factor_base()->equal(*base))) {
                SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "relation admission scratch mismatch");
                return false;
            }
            bool relation_set = false;
            bool deferred_integral_relation_row = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate_set_relation");
                if (can_defer_integral_alpha) {
                    bool handled = false;
                    bool smooth = false;
                    if (!detail::
                                factor_relation_row_from_integral_coordinates_and_norm(
                                        *relation, handled, smooth,
                                        flint::FmpzMatConstRef(coordinates),
                                        flint::FmpqConstRef(known_norm),
                                        known_integral_coordinate_polynomial,
                                        diagnostics)) {
                        return false;
                    }
                    if (handled) {
                        relation_set = smooth;
                        deferred_integral_relation_row = smooth;
                    } else {
                        if (!ensure_alpha()) {
                            return false;
                        }
                        relation_set =
                                detail::set_relation_from_integral_coordinates_and_norm(
                                        *relation, alpha,
                                        flint::FmpzMatConstRef(coordinates),
                                        flint::FmpqConstRef(known_norm),
                                        diagnostics);
                    }
                } else if (known_norm != nullptr &&
                           flint::fmpz_is_one(den)) {
                    relation_set =
                            detail::set_relation_from_integral_coordinates_and_norm(
                                    *relation, alpha,
                                    flint::FmpzMatConstRef(coordinates),
                                    flint::FmpqConstRef(known_norm),
                                    diagnostics);
                } else {
                    relation_set = known_norm == nullptr
                            ? relation->set_generator(alpha, diagnostics)
                            : relation->set_generator_with_norm(
                                      alpha,
                                      flint::FmpqConstRef(known_norm),
                                      diagnostics);
                }
            }
            if (!relation_set) {
                return true;
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.coordinate_candidate_admit_relation");
                bool admitted = false;
                if (deferred_integral_relation_row) {
                    admitted = detail::try_admit_deferred_integral_relation(
                            context, relation_appended, *admission_cache,
                            *relation,
                            flint::FmpzMatConstRef(coordinates),
                            random_relation,
                            admission_scratch_base_verified);
                } else {
                    admitted = detail::try_admit_relation(
                            context, relation_appended, *admission_cache,
                            *relation, random_relation,
                            admission_scratch_base_verified);
                }
                if (!admitted) {
                    SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                              LogLevel::detail,
                              "relation admission failed");
                    return false;
                }
            }
            if (relation_appended) {
                ++accepted_relations;
            }
        } else {
            if (!try_generator_relation(
                        context, alpha, known_norm,
                        known_integral_coordinates,
                        known_integral_coordinate_polynomial, source,
                        accepted_relations, candidate_partial_throttle,
                        relation_appended)) {
                return false;
            }
        }
    }

    if (relation_appended) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.coordinate_candidate_publish");
        if (admission_cache != nullptr) {
            goal_reached = context.relation_count() >= max_relations;
        } else {
            goal_reached = publish_and_check_compute_goal(
                    context, target_relation_kernel_units,
                    requested_certification);
        }
    }
    if (admission_cache == nullptr && candidate_partial_throttle &&
        !goal_reached) {
        partial_throttle_exit = true;
    }
    return true;
}

bool try_coordinate_candidate(ClassGroupContext& context,
                              const Order& order,
                              flint::FmpzMat& coordinates,
                              NormPrefilter* norm_prefilter,
                              ClassGroupRelationSource source,
                              slong target_relation_kernel_units,
                              slong max_candidates,
                              slong max_relations,
                              slong& candidates_tried,
                              slong& accepted_relations,
                              CertificationMode requested_certification,
                              bool& goal_reached,
                              bool& partial_throttle_exit,
                              OrderCoordinateElementConversion*
                                      coordinate_conversion) noexcept {
    flint::Fmpz den;
    fmpz_one(den.raw());
    return try_coordinate_candidate_den(
            context, order, coordinates, flint::FmpzConstRef(den),
            norm_prefilter, source, target_relation_kernel_units,
            max_candidates, max_relations, candidates_tried,
            accepted_relations, requested_certification, goal_reached,
            partial_throttle_exit, nullptr, false, nullptr, 0, nullptr,
            nullptr, coordinate_conversion);
}

}  // namespace detail::relation_search
}  // namespace silex
