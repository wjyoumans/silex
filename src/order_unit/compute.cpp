#include <silex/order_unit.hpp>

#include "compute_internal.hpp"
#include "class_unit_transaction_internal.hpp"
#include "compact_reconstruction_bound_internal.hpp"
#include "relation_unit_internal.hpp"
#include "order_unit_internal.hpp"
#include "../class_group/class_group_internal.hpp"
#include "../class_group/relation_saturation_internal.hpp"
#include "../factored_element/compact_reconstruction_internal.hpp"
#include "../lat/flatter_backend_internal.hpp"
#include "../lat/fplll_backend_internal.hpp"

#include <silex/class_group.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/signature.hpp>
#include <silex/unit.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace silex {

namespace detail {

struct RelationUnitCandidateSortKey {
    flint::Arf max_log_upper;
    std::size_t original_index = 0;
};

struct IncrementalUnitContext {
    const NumberField* field = nullptr;
    slong target_rank = 0;
    std::vector<FactoredElement> units;
    bool full_rank = false;

    bool define(const NumberField& parent, slong rank) noexcept {
        if (rank <= 0) {
            return false;
        }
        field = &parent;
        target_rank = rank;
        units.clear();
        units.reserve(static_cast<std::size_t>(rank));
        full_rank = false;
        return true;
    }

    bool add_independent_unit(bool& changed,
                              const FactoredElement& candidate,
                              EmbeddingContext& embeddings,
                              slong precision) noexcept {
        changed = false;
        if (field == nullptr || target_rank <= 0 || full_rank ||
            precision <= 0 || candidate.parent() == nullptr ||
            !candidate.parent()->has_same_data(*field)) {
            return false;
        }

        units.emplace_back(*field);
        if (!units.back().is_defined() || !units.back().set(candidate)) {
            units.pop_back();
            return false;
        }

        bool independent = false;
        if (!compact_independent(independent, embeddings,
                                 FactoredElementSpan(units.data(),
                                                     units.size()),
                                 precision) ||
            !independent) {
            units.pop_back();
            return true;
        }

        changed = true;
        full_rank = static_cast<slong>(units.size()) == target_rank;
        return true;
    }
};

slong order_unit_candidate_used_count(
        const RelationUnitCandidateBatchState& state) noexcept {
    return static_cast<slong>(state.relations_used.size());
}

bool unit_candidate_has_unused_relations(
        const RelationUnitCandidateBatchState& state) noexcept {
    return order_unit_candidate_used_count(state) <
           state.extra_relation_count;
}

slong unit_not_larger_bound(
        const RelationUnitCandidateBatchState& state,
        slong rank) noexcept {
    return min_slong(min_slong(WORD(20), state.extra_relation_count), rank);
}

slong gcd_slong(slong a, slong b) noexcept {
    while (b != 0) {
        const slong t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

ulong validation_bf_max_cutoff_for_options(
        const Order& order,
        const ClassGroupComputeOptions& options) noexcept {
    return options.requested_certification == CertificationMode::proven &&
                   order.degree() > 2
            ? options.zeta_bf_max_cutoff
            : 0;
}

bool seed_validation_cache_from_finish_product(
        AnalyticClassRegulatorCache& analytic_cache,
        const ClassGroupContext& class_group,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    if (precision <= 0) {
        return false;
    }

    flint::Arb analytic_hR;
    slong product_precision = 0;
    if (!class_group.analytic_finish_product(
                flint::ArbRef(analytic_hR), product_precision) ||
        product_precision < precision ||
        !flint::arb_is_finite(analytic_hR) ||
        !flint::arb_is_positive(analytic_hR)) {
        return false;
    }

    // The source pipeline computes invhr from the inverse residue before the relation
    // loop and accepts it through reconstruct_regulator/bad_check.  The reference-contract
    // route publishes that same accepted product as class-group route state;
    // reuse it for the following validation estimate instead of recomputing
    // the generic zeta product for the same order.
    analytic_cache.seed(order, flint::ArbConstRef(analytic_hR),
                        product_precision);
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.validation_analytic_index_bound.seed_finish_product");
    return true;
}

bool rank_zero_torsion(flint::FmpzRef best_order,
                       OrderElement& best_generator,
                       const Order& order) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr ||
        !same_order_parent(best_generator.parent(), &order)) {
        return false;
    }

    flint::Fmpz discriminant;
    if (order.is_maximal() && order.degree() == 2 &&
        order.discriminant(flint::FmpzRef(discriminant)) &&
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) < 0 &&
        !flint::fmpz_equal_si(flint::FmpzConstRef(discriminant), -3) &&
        !flint::fmpz_equal_si(flint::FmpzConstRef(discriminant), -4)) {
        // reference quad.c:quadunitindex and Silex C's quadratic roots-of-unity
        // path isolate discriminants -3 and -4 as the only exceptional
        // imaginary quadratic fields.  A maximal order in every other case
        // has exactly {1, -1} as its torsion subgroup.
        flint::fmpz_set_ui(best_order, 2);
        return best_generator.set_si(-1);
    }

    flint::Fmpz field_order;
    Element root(*field);
    Element power(*field);
    OrderElement candidate(order);
    if (!root.is_defined() || !power.is_defined() || !candidate.is_defined() ||
        !roots_of_unity(flint::FmpzRef(field_order), root, *field) ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(field_order))) {
        return false;
    }

    const slong n = flint::fmpz_get_si(flint::FmpzConstRef(field_order));
    if (n <= 0 || !best_generator.one()) {
        return false;
    }

    slong best = 1;
    for (slong k = 0; k < n; ++k) {
        if (!compute_power(power, root, k)) {
            return false;
        }
        if (candidate.set_element(power)) {
            const slong ord = n / gcd_slong(n, k);
            if (ord > best) {
                best = ord;
                if (!best_generator.set(candidate)) {
                    return false;
                }
            }
        }
    }

    flint::fmpz_set_si(best_order, best);
    return true;
}

bool unit_candidate_log_key(RelationUnitCandidateSortKey& out,
                                  const FactoredElement& candidate,
                                  EmbeddingContext& embeddings,
                                  std::size_t original_index) noexcept {
    slong places = 0;
    if (!compact_places(places, embeddings) || places <= 0) {
        return false;
    }

    flint::ArbVec logs(places);
    if (!candidate.logarithmic_embedding(
                flint::ArbVecRef(logs), embeddings,
                LogEmbeddingMode::product,
                kUnitCandidateSortPrecision)) {
        return false;
    }

    flint::Arb absolute;
    flint::Arf upper;
    bool have_key = false;
    for (slong i = 0; i < places; ++i) {
        flint::arb_abs(flint::ArbRef(absolute),
                       flint::ArbConstRef(logs.data() + i));
        flint::arb_get_ubound_arf(
                upper, absolute, kUnitCandidateSortPrecision);
        if (!flint::arf_is_finite(upper)) {
            return false;
        }
        if (!have_key || flint::arf_cmp(upper, out.max_log_upper) > 0) {
            flint::arf_set(out.max_log_upper, upper.raw());
            have_key = true;
        }
    }

    out.original_index = original_index;
    return have_key;
}

bool sort_unit_candidates_by_log_height(
        std::vector<FactoredElement>& candidates,
        EmbeddingContext& embeddings) noexcept {
    if (candidates.size() < 2) {
        return true;
    }

    std::vector<RelationUnitCandidateSortKey> keys(candidates.size());
    std::vector<std::size_t> order(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (!unit_candidate_log_key(keys[i], candidates[i], embeddings,
                                          i)) {
            return false;
        }
        order[i] = i;
    }

    std::sort(order.begin(), order.end(),
              [&](std::size_t left, std::size_t right) noexcept {
                  const int cmp = flint::arf_cmp(keys[left].max_log_upper,
                                                 keys[right].max_log_upper);
                  if (cmp != 0) {
                      return cmp < 0;
                  }
                  return keys[left].original_index <
                         keys[right].original_index;
              });

    std::vector<FactoredElement> sorted;
    sorted.reserve(candidates.size());
    for (std::size_t index : order) {
        sorted.push_back(std::move(candidates[index]));
    }
    candidates.swap(sorted);
    return true;
}

bool torsion_log_bound(flint::Arb& out,
                             slong& log_precision,
                             slong degree) noexcept {
    if (degree <= 0 || log_precision <= 0) {
        return false;
    }

    flint::Arb log_degree;
    flint::Arb term;
    flint::Arb one_plus_term;
    flint::Arb log_bound;
    flint::Arb neg_log_bound;
    flint::Arb log_two;
    flint::Arb ratio;
    flint::Arf upper;
    const ulong degree_ui = static_cast<ulong>(degree);
    const slong input_precision =
            max_slong(log_precision, kTorsionStartPrecision);

    flint::arb_log_ui(log_degree, degree_ui, input_precision);
    flint::arb_div_ui(term, log_degree, 6, input_precision);
    flint::arb_div_ui(term, term, degree_ui, input_precision);
    flint::arb_div_ui(term, term, degree_ui, input_precision);
    flint::arb_add_ui(one_plus_term, term, 1, input_precision);
    flint::arb_log(out, one_plus_term, input_precision);
    if (!flint::arb_is_positive(out)) {
        return false;
    }

    flint::arb_log(log_bound, out, input_precision);
    flint::arb_neg(neg_log_bound, log_bound);
    flint::arb_log_ui(log_two, 2, input_precision);
    flint::arb_div(ratio, neg_log_bound, log_two, input_precision);
    if (!flint::arb_is_positive(ratio)) {
        return false;
    }

    flint::arb_get_ubound_arf(upper, ratio, input_precision);
    if (!flint::arf_is_finite(upper)) {
        return false;
    }
    const double precision_estimate = flint::arf_get_d(upper, ARF_RND_CEIL);
    if (!std::isfinite(precision_estimate) || precision_estimate < 0.0 ||
        precision_estimate >
                static_cast<double>(std::numeric_limits<slong>::max() - 2)) {
        return false;
    }
    log_precision = static_cast<slong>(precision_estimate) + 2;
    if (log_precision < kTorsionStartPrecision) {
        log_precision = kTorsionStartPrecision;
    }
    return true;
}

bool unit_candidate_torsion_status(
        RelationTorsionStatus& status,
        RelationUnitExtractionState& extraction_state,
        const FactoredElement& candidate,
        EmbeddingContext& embeddings) noexcept {
    status = RelationTorsionStatus::inconclusive;

    const NumberField* field = embeddings.parent();
    if (field == nullptr || candidate.parent() == nullptr ||
        !candidate.parent()->has_same_data(*field)) {
        return false;
    }

    Signature sig;
    if (embeddings.is_set()) {
        sig = embeddings.signature();
    } else if (!signature(sig, *field)) {
        return false;
    }

    const slong degree = sig.degree();
    const slong places = sig.r1() + sig.r2();
    slong log_precision = extraction_state.torsion_precision;
    flint::Arb bound;
    if (places <= 0 ||
        !torsion_log_bound(bound, log_precision, degree)) {
        return false;
    }

    flint::ArbVec logs(places);
    if (!candidate.logarithmic_embedding(flint::ArbVecRef(logs), embeddings,
                                         LogEmbeddingMode::product,
                                         log_precision)) {
        return false;
    }
    extraction_state.torsion_precision =
            max_slong(extraction_state.torsion_precision, log_precision);

    flint::Arb value;
    flint::Arb difference;
    slong small_places = 0;
    for (slong i = 0; i < places; ++i) {
        flint::arb_abs(flint::ArbRef(value),
                       flint::ArbConstRef(logs.data() + i));
        if (i >= sig.r1()) {
            flint::arb_div_ui(value, value, 2, log_precision);
        }
        if (flint::arb_is_positive(value)) {
            status = RelationTorsionStatus::non_torsion;
            return true;
        }
        flint::arb_sub(difference, bound, value, log_precision);
        if (flint::arb_is_nonnegative(difference.raw())) {
            ++small_places;
        } else {
            status = RelationTorsionStatus::inconclusive;
            return true;
        }
    }

    status = small_places == places ? RelationTorsionStatus::torsion
                                    : RelationTorsionStatus::inconclusive;
    return true;
}

bool normalized_log_cutoff_row(flint::ArbMat& out,
                                     slong row,
                                     const FactoredElement& value,
                                     EmbeddingContext& embeddings,
                                     slong precision,
                                     const DiagnosticsContext* diagnostics)
        noexcept {
    Signature sig;
    if (embeddings.is_set()) {
        sig = embeddings.signature();
    } else if (embeddings.parent() == nullptr ||
               !signature(sig, *embeddings.parent())) {
        return false;
    }

    const slong places = sig.r1() + sig.r2();
    const slong rank = places - 1;
    if (precision <= 0 || row < 0 ||
        row >= flint::arb_mat_nrows_value(out) ||
        flint::arb_mat_ncols_value(out) != rank || rank <= 0 ||
        value.parent() == nullptr || embeddings.parent() == nullptr ||
        !value.parent()->has_same_data(*embeddings.parent())) {
        return false;
    }

    flint::ArbVec logs(places);
    if (!value.logarithmic_embedding(flint::ArbVecRef(logs), embeddings,
                                     LogEmbeddingMode::product, precision,
                                     diagnostics)) {
        return false;
    }

    flint::Arb sum;
    flint::Arb average;
    flint::Arb correction;
    flint::arb_zero(sum);
    for (slong j = 0; j < places; ++j) {
        flint::arb_add(sum, sum, logs.data() + j, precision);
    }
    flint::arb_div_ui(average, sum, static_cast<ulong>(sig.degree()),
                      precision);

    for (slong j = 0; j < rank; ++j) {
        if (j < sig.r1()) {
            flint::arb_set(correction, average);
        } else {
            flint::arb_mul_ui(correction, average, 2, precision);
        }
        ::arb_sub(arb_mat_entry(out.raw(), row, j),
                  logs.data() + j, correction.raw(), precision);
    }
    return true;
}

bool normalized_log_cutoff_matrix(
        flint::ArbMat& out,
        const std::vector<FactoredElement>& values,
        EmbeddingContext& embeddings,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    if (flint::arb_mat_nrows_value(out) !=
                static_cast<slong>(values.size()) ||
        precision <= 0) {
        return false;
    }
    for (slong i = 0; i < static_cast<slong>(values.size()); ++i) {
        if (!normalized_log_cutoff_row(
                    out, i, values[static_cast<std::size_t>(i)],
                    embeddings, precision, diagnostics)) {
            return false;
        }
    }
    return true;
}

bool unit_log_cutoff_matrix(
        flint::ArbMat& out,
        const std::vector<FactoredElement>& units,
        EmbeddingContext& embeddings,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong rank = static_cast<slong>(units.size());
    slong places = 0;
    if (rank <= 0 || precision <= 0 ||
        !compact_places(places, embeddings) || places != rank + 1 ||
        flint::arb_mat_nrows_value(out) != rank ||
        flint::arb_mat_ncols_value(out) != rank) {
        return false;
    }

    // reference reduce_mod_units inverts _conj_log_mat_cutoff(U.units, p), which
    // uses the first r raw conjugate-log columns.  Only the value rows are
    // normalized by _conj_arb_log_matrix_normalise_cutoff.
    flint::ArbVec logs(places);
    for (slong i = 0; i < rank; ++i) {
        if (!units[static_cast<std::size_t>(i)].logarithmic_embedding(
                    flint::ArbVecRef(logs), embeddings,
                    LogEmbeddingMode::product, precision, diagnostics)) {
            return false;
        }
        for (slong j = 0; j < rank; ++j) {
            arb_set(arb_mat_entry(out.raw(), i, j), logs.data() + j);
        }
    }
    return true;
}

bool copy_unit_generators(std::vector<FactoredElement>& out,
                                const OrderUnitGroup& group) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || rank <= 0) {
        return false;
    }

    out.clear();
    out.reserve(static_cast<std::size_t>(rank));
    for (slong i = 0; i < rank; ++i) {
        out.emplace_back(*field);
        if (!group.free_generator(out.back(), i)) {
            return false;
        }
    }
    return true;
}

bool round_arb_to_nearest_fmpz(flint::FmpzRef out,
                                     const arb_struct* value,
                                     slong precision) noexcept {
    if (value == nullptr || precision <= 0) {
        return false;
    }

    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value));
    if (!flint::arf_is_finite(midpoint)) {
        return false;
    }
    flint::arf_get_fmpz(out, midpoint, ARF_RND_NEAR);

    flint::Arb rounded;
    flint::Arb difference;
    flint::Arb absolute;
    flint::Arb half;
    flint::arb_set_fmpz(rounded, flint::FmpzConstRef(out.raw()));
    flint::arb_sub(difference, value, rounded, precision);
    flint::arb_abs(absolute, difference);
    flint::arb_set_ui(half, 1);
    flint::arb_div_ui(half, half, 2, precision);
    return flint::arb_lt(absolute, half);
}

bool round_arb_approx_fmpz(flint::FmpzRef out,
                                 const arb_struct* value,
                                 slong precision) noexcept {
    if (value == nullptr || precision <= 0 || !arb_is_finite(value)) {
        return false;
    }

    flint::Arf upper_bound;
    flint::Arf lower_bound;
    flint::Fmpz upper;
    flint::Fmpz lower;
    arb_get_ubound_arf(upper_bound.raw(), value, precision);
    arb_get_lbound_arf(lower_bound.raw(), value, precision);
    if (!flint::arf_is_finite(upper_bound) ||
        !flint::arf_is_finite(lower_bound)) {
        return false;
    }
    flint::arf_get_fmpz(upper, upper_bound, ARF_RND_CEIL);
    flint::arf_get_fmpz(lower, lower_bound, ARF_RND_FLOOR);

    flint::Fmpz source_width;
    fmpz_sub(source_width.raw(), lower.raw(), upper.raw());
    flint::Arb width;
    flint::Arb absolute;
    flint::Arb allowed;
    arb_set_fmpz(width.raw(), source_width.raw());
    arb_abs(absolute.raw(), value);
    arb_sqrt(allowed.raw(), absolute.raw(), precision);
    if (!arb_is_finite(allowed.raw()) ||
        arb_gt(width.raw(), allowed.raw())) {
        return false;
    }

    flint::Fmpz sum;
    fmpz_add(sum.raw(), upper.raw(), lower.raw());
    fmpz_tdiv_q_2exp(out.raw(), sum.raw(), 1);
    return true;
}

slong next_precision_bucket(slong precision) noexcept {
    if (precision <= 0) {
        return precision;
    }

    slong bucket = 1;
    while (bucket <= precision) {
        if (bucket > std::numeric_limits<slong>::max() / 2) {
            return bucket;
        }
        bucket *= 2;
    }
    return bucket;
}

enum class UnitReductionPrecisionEventKind {
    start,
    selected,
    inverse_work,
    inverse_hit,
};

void profile_reduce_precision(
        const DiagnosticsContext* diagnostics,
        UnitReductionPrecisionEventKind kind,
        slong precision) noexcept {
    if (kind == UnitReductionPrecisionEventKind::start) {
        if (precision <= 32) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.le32");
        } else if (precision <= 64) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.le64");
        } else if (precision <= 128) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.le128");
        } else if (precision <= 256) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.le256");
        } else if (precision <= 512) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.le512");
        } else {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.start.gt512");
        }
        return;
    }

    if (kind == UnitReductionPrecisionEventKind::selected) {
        if (precision <= 32) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.le32");
        } else if (precision <= 64) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.le64");
        } else if (precision <= 128) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.le128");
        } else if (precision <= 256) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.le256");
        } else if (precision <= 512) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.le512");
        } else {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.selected.gt512");
        }
        return;
    }

    if (kind == UnitReductionPrecisionEventKind::inverse_work) {
        if (precision <= 32) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.le32");
        } else if (precision <= 64) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.le64");
        } else if (precision <= 128) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.le128");
        } else if (precision <= 256) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.le256");
        } else if (precision <= 512) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.le512");
        } else {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.reduce_mod_units.precision.inverse_work.gt512");
        }
        return;
    }

    if (precision <= 32) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.le32");
    } else if (precision <= 64) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.le64");
    } else if (precision <= 128) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.le128");
    } else if (precision <= 256) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.le256");
    } else if (precision <= 512) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.le512");
    } else {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.precision.inverse_hit.gt512");
    }
}

void profile_reduce_stored_precision(
        const DiagnosticsContext* diagnostics,
        slong precision) noexcept {
    if (precision <= 16) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.le16");
    } else if (precision <= 32) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.le32");
    } else if (precision <= 64) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.le64");
    } else if (precision <= 128) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.le128");
    } else if (precision <= 256) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.le256");
    } else {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_stored_units.precision.gt256");
    }
}

bool reduce_mod_units_start_precision(
        slong& out,
        const std::vector<FactoredElement>& values,
        const OrderUnitGroup& group,
        EmbeddingContext& embeddings,
        slong start_precision) noexcept {
    const DiagnosticsContext* diagnostics = group.diagnostics();
    out = start_precision;
    profile_reduce_precision(
            diagnostics, UnitReductionPrecisionEventKind::start,
            start_precision);
    const slong rank = group.free_rank();
    if (values.empty() || rank <= 0) {
        return true;
    }

    flint::Arb regulator;
    if (!group.regulator(flint::ArbRef(regulator)) ||
        !flint::arb_is_finite(regulator) ||
        !flint::arb_is_positive(regulator)) {
        return true;
    }

    flint::ArbMat value_logs(static_cast<slong>(values.size()), rank);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.reduce_mod_units.start_value_log_matrix");
        if (!normalized_log_cutoff_matrix(
                    value_logs, values, embeddings, start_precision,
                    diagnostics)) {
            return false;
        }
    }

    flint::Arb max_norm;
    flint::Arb row_sum;
    flint::Arb square;
    flint::Arb row_norm;
    bool have_norm = false;
    for (slong i = 0; i < static_cast<slong>(values.size()); ++i) {
        flint::arb_zero(row_sum);
        for (slong j = 0; j < rank; ++j) {
            arb_mul(square.raw(), arb_mat_entry(value_logs.raw(), i, j),
                    arb_mat_entry(value_logs.raw(), i, j), start_precision);
            flint::arb_add(row_sum, row_sum, square, start_precision);
        }
        arb_sqrt(row_norm.raw(), row_sum.raw(), start_precision);
        if (!have_norm) {
            flint::arb_set(max_norm, row_norm);
            have_norm = true;
        } else {
            arb_max(max_norm.raw(), max_norm.raw(), row_norm.raw(),
                    start_precision);
        }
    }
    if (!have_norm) {
        return true;
    }

    flint::Arb regulator_root;
    flint::Arb bound;
    arb_root_ui(regulator_root.raw(), regulator.raw(),
                static_cast<ulong>(rank), start_precision);
    arb_div(bound.raw(), max_norm.raw(), regulator_root.raw(),
            start_precision);
    if (!arb_is_finite(bound.raw())) {
        return true;
    }

    const slong bits = arb_bits(bound.raw());
    if (bits > out) {
        out = bits;
    }
    out = next_precision_bucket(out);
    profile_reduce_precision(
            diagnostics, UnitReductionPrecisionEventKind::selected, out);
    return out > 0;
}

bool reduction_coordinates(flint::FmpzMat& out,
                                 const std::vector<FactoredElement>& values,
                                 const OrderUnitGroup& group,
                                 EmbeddingContext& embeddings,
                                 bool& exact,
                                 slong precision,
                                 const flint::ArbMat* unit_log_inverse)
        noexcept {
    exact = false;
    const slong len = static_cast<slong>(values.size());
    const slong rank = group.free_rank();
    if (len < 0 || rank <= 0 ||
        flint::fmpz_mat_nrows(flint::FmpzMatRef(out)) != len ||
        flint::fmpz_mat_ncols(flint::FmpzMatRef(out)) != rank ||
        precision <= 0) {
        return false;
    }

    flint::ArbMat value_logs(len, rank);
    const DiagnosticsContext* diagnostics = group.diagnostics();
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.reduce_mod_units.value_log_matrix");
        if (!normalized_log_cutoff_matrix(
                    value_logs, values, embeddings, precision, diagnostics)) {
            return false;
        }
    }

    flint::ArbMat coordinate_rows(len, rank);
    if (unit_log_inverse != nullptr) {
        if (flint::arb_mat_nrows_value(*unit_log_inverse) != rank ||
            flint::arb_mat_ncols_value(*unit_log_inverse) != rank) {
            return false;
        }
        arb_mat_mul(coordinate_rows.raw(), value_logs.raw(),
                    unit_log_inverse->raw(), precision);
    } else {
        std::vector<FactoredElement> units;
        if (!copy_unit_generators(units, group)) {
            return false;
        }

        flint::ArbMat unit_logs(rank, rank);
        flint::ArbMat system(rank, rank);
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_mod_units.fallback_unit_log_matrix");
            if (!unit_log_cutoff_matrix(
                        unit_logs, units, embeddings, precision,
                        diagnostics)) {
                return false;
            }
        }
        arb_mat_transpose(system.raw(), unit_logs.raw());

        flint::ArbMat rhs(rank, 1);
        flint::ArbMat coordinates(rank, 1);
        for (slong i = 0; i < len; ++i) {
            for (slong j = 0; j < rank; ++j) {
                arb_set(arb_mat_entry(rhs.raw(), j, 0),
                        arb_mat_entry(value_logs.raw(), i, j));
            }
            if (!arb_mat_solve(coordinates.raw(), system.raw(), rhs.raw(),
                               precision)) {
                return false;
            }
            for (slong j = 0; j < rank; ++j) {
                arb_set(arb_mat_entry(coordinate_rows.raw(), i, j),
                        arb_mat_entry(coordinates.raw(), j, 0));
            }
        }
    }

    flint::FmpzMat candidate(len, rank);
    bool exact_rounding = true;
    for (slong i = 0; exact_rounding && i < len; ++i) {
        for (slong j = 0; j < rank; ++j) {
            exact_rounding = round_arb_to_nearest_fmpz(
                    flint::fmpz_mat_entry(candidate, i, j),
                    arb_mat_entry(coordinate_rows.raw(), i, j), precision);
            if (!exact_rounding) {
                break;
            }
        }
    }
    if (!exact_rounding) {
        for (slong i = 0; i < len; ++i) {
            for (slong j = 0; j < rank; ++j) {
                if (!round_arb_approx_fmpz(
                            flint::fmpz_mat_entry(candidate, i, j),
                            arb_mat_entry(coordinate_rows.raw(), i, j),
                            precision)) {
                    return false;
                }
            }
        }
    }

    exact = exact_rounding;
    out = std::move(candidate);
    return true;
}

bool fmpz_matrix_is_zero(const flint::FmpzMat& matrix) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            if (!flint::fmpz_is_zero(flint::fmpz_mat_entry(matrix, i, j))) {
                return false;
            }
        }
    }
    return true;
}

bool apply_unit_reduction(std::vector<FactoredElement>& values,
                                const OrderUnitGroup& group,
                                const flint::FmpzMat& coefficients) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong len = static_cast<slong>(values.size());
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || rank <= 0 ||
        flint::fmpz_mat_nrows(coefficients) != len ||
        flint::fmpz_mat_ncols(coefficients) != rank) {
        return false;
    }

    std::vector<FactoredElement> reduced;
    reduced.reserve(values.size());
    for (slong i = 0; i < len; ++i) {
        reduced.emplace_back(*field);
        if (!reduced.back().set(values[static_cast<std::size_t>(i)])) {
            return false;
        }
        for (slong j = 0; j < rank; ++j) {
            const flint::FmpzConstRef coefficient =
                    flint::fmpz_mat_entry(coefficients, i, j);
            if (flint::fmpz_is_zero(coefficient)) {
                continue;
            }

            flint::Fmpz exponent;
            flint::fmpz_neg(flint::FmpzRef(exponent), coefficient);
            if (!flint::fmpz_fits_si(flint::FmpzConstRef(exponent))) {
                return true;
            }

            FactoredElement unit(*field);
            FactoredElement power(*field);
            FactoredElement product(*field);
            if (!group.free_generator(unit, j) ||
                !power.pow_si(unit,
                              flint::fmpz_get_si(
                                      flint::FmpzConstRef(exponent))) ||
                !product.multiply(reduced.back(), power) ||
                !reduced.back().set(product)) {
                return false;
            }
        }
    }

    values = std::move(reduced);
    return true;
}

bool scaled_log_matrix(flint::FmpzMat& out,
                             const std::vector<FactoredElement>& units,
                             EmbeddingContext& embeddings,
                             slong& precision) noexcept {
    if (precision <= 0 || units.empty() ||
        embeddings.parent() == nullptr) {
        return false;
    }

    slong places = 0;
    if (!compact_places(places, embeddings) ||
        flint::fmpz_mat_nrows(out) != static_cast<slong>(units.size()) ||
        flint::fmpz_mat_ncols(out) != places) {
        return false;
    }

    for (;;) {
        bool too_low = false;
        for (slong i = 0; i < static_cast<slong>(units.size()); ++i) {
            flint::ArbVec logs(places);
            if (!units[static_cast<std::size_t>(i)].logarithmic_embedding(
                        flint::ArbVecRef(logs), embeddings,
                        LogEmbeddingMode::product, precision)) {
                return false;
            }
            for (slong j = 0; j < places; ++j) {
                if (mag_get_d(arb_radref(logs.data() + j)) > 1.0e-9) {
                    too_low = true;
                    break;
                }
                arf_get_fmpz_fixed_si(
                        flint::fmpz_mat_entry(out, i, j).raw(),
                        arb_midref(logs.data() + j), -precision);
            }
            if (too_low) {
                break;
            }
        }
        if (!too_low) {
            return true;
        }
        if (precision > (WORD(1) << 29)) {
            return false;
        }
        precision *= 2;
    }
}

bool row_norm_product_bits(flint_bitcnt_t& out,
                                 const flint::FmpzMat& matrix) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    flint::Fmpz product;
    flint::Fmpz norm;
    flint::fmpz_one(flint::FmpzRef(product));
    for (slong i = 0; i < rows; ++i) {
        flint::fmpz_zero(flint::FmpzRef(norm));
        for (slong j = 0; j < cols; ++j) {
            const flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(matrix, i, j);
            flint::fmpz_addmul(flint::FmpzRef(norm), entry, entry);
        }
        flint::fmpz_mul(flint::FmpzRef(product),
                        flint::FmpzConstRef(product),
                        flint::FmpzConstRef(norm));
    }
    out = fmpz_bits(product.raw());
    return true;
}

void profile_unit_lll_backend(const DiagnosticsContext* diagnostics,
                                    const char* label) noexcept {
    (void)label;
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group, label);
}

bool unit_transform_is_unimodular(
        const flint::FmpzMat& transform) noexcept {
    if (flint::fmpz_mat_nrows(transform) !=
        flint::fmpz_mat_ncols(transform)) {
        return false;
    }
    flint::Fmpz det;
    flint::fmpz_mat_det(flint::FmpzRef(det),
                        flint::FmpzMatConstRef(transform));
    return flint::fmpz_is_pm1(flint::FmpzConstRef(det));
}

bool unit_transform_matches_scaled(
        const flint::FmpzMat& transform,
        const flint::FmpzMat& scaled,
        const flint::FmpzMat& reduced) noexcept {
    if (flint::fmpz_mat_nrows(transform) !=
                flint::fmpz_mat_nrows(scaled) ||
        flint::fmpz_mat_ncols(transform) !=
                flint::fmpz_mat_nrows(scaled) ||
        flint::fmpz_mat_nrows(reduced) !=
                flint::fmpz_mat_nrows(scaled) ||
        flint::fmpz_mat_ncols(reduced) !=
                flint::fmpz_mat_ncols(scaled)) {
        return false;
    }
    flint::FmpzMat expected(flint::fmpz_mat_nrows(scaled),
                            flint::fmpz_mat_ncols(scaled));
    flint::fmpz_mat_mul(flint::FmpzMatRef(expected),
                        flint::FmpzMatConstRef(transform),
                        flint::FmpzMatConstRef(scaled));
    return flint::fmpz_mat_equal(flint::FmpzMatConstRef(expected),
                                 flint::FmpzMatConstRef(reduced));
}

bool unit_lll_result_is_valid(
        const flint::FmpzMat& reduced,
        const flint::FmpzMat& transform,
        const flint::FmpzMat& scaled) noexcept {
    return unit_transform_is_unimodular(transform) &&
           unit_transform_matches_scaled(transform, scaled, reduced);
}

bool reduce_unit_log_rows_with_flatter(
        flint::FmpzMat& reduced,
        flint::FmpzMat& transform,
        const flint::FmpzMat& scaled,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong rows = flint::fmpz_mat_nrows(scaled);
    const slong cols = flint::fmpz_mat_ncols(scaled);
    if (rows <= 0 || cols <= 0 || flint::fmpz_mat_rank(scaled) != rows) {
        return false;
    }

    // reference's FLATTER backend is a column-image reducer.  Stored units need row
    // LLL, so transpose into that contract and transpose the square transform
    // back before applying it to the unit basis.
    flint::FmpzMat transposed(cols, rows);
    flint::FmpzMat reduced_transposed(cols, rows);
    flint::FmpzMat column_transform(rows, rows);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transposed),
                              flint::FmpzMatConstRef(scaled));
    const lat::detail::FlatterLllResult result =
            lat::detail::flatter_column_lll_transform(
                    flint::FmpzMatRef(reduced_transposed),
                    flint::FmpzMatRef(column_transform),
                    flint::FmpzMatConstRef(transposed), 1.02, 1);
    if (result.status != lat::detail::FlatterBackendStatus::success ||
        result.rank != rows) {
        profile_unit_lll_backend(
                diagnostics,
                result.status ==
                                lat::detail::FlatterBackendStatus::unavailable
                        ? "unit_group.reduce_stored_units.lll_backend."
                          "flatter_unavailable"
                        : "unit_group.reduce_stored_units.lll_backend."
                          "flatter_failed");
        return false;
    }

    flint::fmpz_mat_transpose(flint::FmpzMatRef(reduced),
                              flint::FmpzMatConstRef(reduced_transposed));
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transform),
                              flint::FmpzMatConstRef(column_transform));
    if (!unit_lll_result_is_valid(reduced, transform, scaled)) {
        profile_unit_lll_backend(
                diagnostics,
                "unit_group.reduce_stored_units.lll_backend."
                "flatter_rejected");
        return false;
    }

    profile_unit_lll_backend(
            diagnostics,
            "unit_group.reduce_stored_units.lll_backend.flatter");
    return true;
}

bool reduce_unit_log_rows_with_fplll(
        flint::FmpzMat& reduced,
        flint::FmpzMat& transform,
        const flint::FmpzMat& scaled,
        const DiagnosticsContext* diagnostics) noexcept {
    const lat::detail::FplllLllResult result =
            lat::detail::fplll_row_lll_transform(
                    flint::FmpzMatRef(reduced), flint::FmpzMatRef(transform),
                    flint::FmpzMatConstRef(scaled), 0.99);
    if (result.status != lat::detail::FplllBackendStatus::success) {
        profile_unit_lll_backend(
                diagnostics,
                result.status == lat::detail::FplllBackendStatus::unavailable
                        ? "unit_group.reduce_stored_units.lll_backend."
                          "fplll_unavailable"
                        : "unit_group.reduce_stored_units.lll_backend."
                          "fplll_failed");
        return false;
    }
    if (!unit_lll_result_is_valid(reduced, transform, scaled)) {
        profile_unit_lll_backend(
                diagnostics,
                "unit_group.reduce_stored_units.lll_backend.fplll_rejected");
        return false;
    }

    profile_unit_lll_backend(
            diagnostics, "unit_group.reduce_stored_units.lll_backend.fplll");
    return true;
}

void reduce_unit_log_rows_with_flint(flint::FmpzMat& reduced,
                                           flint::FmpzMat& transform,
                                           const flint::FmpzMat& scaled,
                                           const DiagnosticsContext* diagnostics)
        noexcept {
    flint::fmpz_mat_set(flint::FmpzMatRef(reduced),
                        flint::FmpzMatConstRef(scaled));
    flint::fmpz_mat_one(flint::FmpzMatRef(transform));
    flint::FmpzLll lll;
    fmpz_lll(reduced.raw(), transform.raw(), lll.raw());
    profile_unit_lll_backend(
            diagnostics, "unit_group.reduce_stored_units.lll_backend.flint");
}

bool reduce_unit_log_rows(flint::FmpzMat& reduced,
                                flint::FmpzMat& transform,
                                const flint::FmpzMat& scaled,
                                const DiagnosticsContext* diagnostics)
        noexcept {
#if defined(SILEX_WITH_FLATTER) && SILEX_WITH_FLATTER
    if (reduce_unit_log_rows_with_flatter(reduced, transform, scaled,
                                               diagnostics)) {
        return true;
    }
#endif
#if defined(SILEX_WITH_FPLLL) && SILEX_WITH_FPLLL
    if (reduce_unit_log_rows_with_fplll(reduced, transform, scaled,
                                             diagnostics)) {
        return true;
    }
#endif

    reduce_unit_log_rows_with_flint(reduced, transform, scaled,
                                          diagnostics);
    return unit_lll_result_is_valid(reduced, transform, scaled);
}

bool transform_factored_units(std::vector<FactoredElement>& units,
                                    const flint::FmpzMat& transform) noexcept {
    const slong rank = static_cast<slong>(units.size());
    if (rank <= 0 || flint::fmpz_mat_nrows(transform) != rank ||
        flint::fmpz_mat_ncols(transform) != rank ||
        units.front().parent() == nullptr) {
        return false;
    }

    const NumberField& field = *units.front().parent();
    std::vector<FactoredElement> transformed;
    transformed.reserve(units.size());
    for (slong i = 0; i < rank; ++i) {
        transformed.emplace_back(field);
        if (!transformed.back().is_defined() || !transformed.back().one()) {
            return false;
        }
        for (slong j = 0; j < rank; ++j) {
            const flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(transform, i, j);
            if (flint::fmpz_is_zero(entry)) {
                continue;
            }
            if (!flint::fmpz_fits_si(entry)) {
                return false;
            }

            FactoredElement power(field);
            FactoredElement product(field);
            if (!power.pow_si(units[static_cast<std::size_t>(j)],
                              flint::fmpz_get_si(entry)) ||
                !product.multiply(transformed.back(), power) ||
                !transformed.back().set(product)) {
                return false;
            }
        }
        transformed.back().normalize();
    }

    units = std::move(transformed);
    return true;
}

bool reduce_stored_units_vector(bool& changed,
                                      std::vector<FactoredElement>& units,
                                      EmbeddingContext& embeddings,
                                      slong precision,
                                      const DiagnosticsContext* diagnostics)
        noexcept {
    changed = false;
    if (units.empty()) {
        return true;
    }
    if (precision <= 0 || embeddings.parent() == nullptr) {
        return false;
    }

    slong work_precision = precision;
    for (;;) {
        profile_reduce_stored_precision(diagnostics, work_precision);
        slong places = 0;
        if (!compact_places(places, embeddings)) {
            return false;
        }
        flint::FmpzMat scaled(static_cast<slong>(units.size()), places);
        if (!scaled_log_matrix(scaled, units, embeddings,
                                     work_precision)) {
            return false;
        }

        flint::FmpzMat reduced(flint::fmpz_mat_nrows(scaled),
                               flint::fmpz_mat_ncols(scaled));
        flint::FmpzMat transform(flint::fmpz_mat_nrows(scaled),
                                 flint::fmpz_mat_nrows(scaled));
        if (!reduce_unit_log_rows(reduced, transform, scaled,
                                        diagnostics)) {
            return false;
        }
        if (fmpz_mat_is_one(transform.raw()) != 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_stored_units.transform_identity");
            return true;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.reduce_stored_units.transform_nonidentity");

        flint_bitcnt_t before_bits = 0;
        flint_bitcnt_t after_bits = 0;
        if (!row_norm_product_bits(before_bits, scaled) ||
            !row_norm_product_bits(after_bits, reduced) ||
            !transform_factored_units(units, transform)) {
            return false;
        }
        changed = true;
        if (after_bits >= before_bits) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_stored_units.norm_bits_not_improved");
            return true;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.reduce_stored_units.norm_bits_improved");
    }
}

bool reduce_stored_relation_units(bool& changed,
                               OrderUnitGroup& group,
                               EmbeddingContext& embeddings,
                               slong precision) noexcept {
    changed = false;
    if (!group.is_set() || group.free_rank() <= 0) {
        return true;
    }

    const Order* order = group.parent();
    if (order == nullptr || order->parent() == nullptr) {
        return false;
    }

    std::vector<FactoredElement> units;
    if (!copy_unit_generators(units, group) ||
        !reduce_stored_units_vector(changed, units, embeddings,
                                          precision, group.diagnostics())) {
        return false;
    }
    if (!changed) {
        return true;
    }

    OrderUnitGroup reduced(*order);
    reduced.set_diagnostics(group.diagnostics());
    if (!reduced.is_defined() ||
        !order_unit_group_set_units_internal(
                reduced, *order,
                FactoredElementSpan(units.data(), units.size()), embeddings,
                // reference assigns `U.units = reduce(U.units, p)` directly.  The
                // reduction is a unit-basis transform of already stored units,
                // so keep public validation out of this internal path.
                precision, true)) {
        return false;
    }
    group.swap(reduced);
    return true;
}

bool reconstruct_native_unit_log_generators(
        std::vector<FactoredElement>& out,
        const Order& order,
        std::span<const FactoredElement> generators,
        EmbeddingContext& embeddings,
        const DiagnosticsContext* diagnostics) noexcept {
    // reference vec_chinese_units derives one batch-wide coordinate bound before
    // reconstructing the complete compact fundamental-unit basis.
    flint::Fmpz coordinate_bound;
    CompactCoordinateBoundReport bound_report;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.reduce_mod_units.expand_unit_generators."
                "coordinate_bound");
        if (!compact_unit_coordinate_bound(
                    bound_report, flint::FmpzRef(coordinate_bound), order,
                    FactoredElementSpan(generators.data(), generators.size()),
                    embeddings)) {
            return false;
        }
    }

    std::vector<Element> reconstructed;
    BoundedCompactReconstructionReport reconstruction_report;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.reduce_mod_units.expand_unit_generators."
                "bounded_crt");
        if (!bounded_compact_reconstruct(
                    reconstruction_report, reconstructed, order, generators,
                    flint::FmpzConstRef(coordinate_bound))) {
            return false;
        }
    }

    const NumberField* field = order.parent();
    if (field == nullptr || reconstructed.size() != generators.size()) {
        return false;
    }
    std::vector<FactoredElement> candidate;
    candidate.reserve(reconstructed.size());
    for (const Element& value : reconstructed) {
        candidate.emplace_back(*field);
        if (!candidate.back().is_defined() ||
            !candidate.back().set_element(value)) {
            return false;
        }
    }
    out.swap(candidate);
    return true;
}

bool reduce_mod_units_cutoff_inverse(
        const flint::ArbMat*& out,
        RelationUnitExtractionState& extraction_state,
        const OrderUnitGroup& group,
        EmbeddingContext& embeddings,
        slong& precision,
        slong max_precision) noexcept {
    out = nullptr;
    const slong rank = group.free_rank();
    if (!group.is_set() || rank <= 0 || precision <= 0 ||
        max_precision < precision || embeddings.parent() == nullptr) {
        return false;
    }

    const DiagnosticsContext* diagnostics = group.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.reduce_mod_units.cutoff_inverse");

    for (const auto& entry : extraction_state.reduce_mod_units_inverses) {
        if (entry.precision == precision && entry.rank == rank) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_mod_units.cutoff_inverse_cache_hit");
            profile_reduce_precision(
                    diagnostics, UnitReductionPrecisionEventKind::inverse_hit,
                    precision);
            out = &entry.inverse;
            return true;
        }
    }

    std::vector<FactoredElement> uncached_units;
    std::vector<FactoredElement>* units = &uncached_units;
    if (extraction_state.expand_reduce_mod_units_log_generators) {
        units = &extraction_state.reduce_mod_units_log_generators;
    }
    if (units->empty()) {
        if (!copy_unit_generators(*units, group)) {
            return false;
        }
    }
    if (extraction_state.expand_reduce_mod_units_log_generators &&
        units == &extraction_state.reduce_mod_units_log_generators &&
        !extraction_state
                 .reduce_mod_units_log_generator_expansion_attempted) {
        extraction_state.reduce_mod_units_log_generator_expansion_attempted =
                true;
        bool all_single_factor = true;
        for (const FactoredElement& unit : *units) {
            if (unit.length() != 1) {
                all_single_factor = false;
                break;
            }
        }
        if (!all_single_factor) {
            const Order* order = group.parent();
            if (order == nullptr) {
                return false;
            }
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_mod_units.expand_unit_generators");
            std::vector<FactoredElement> expanded_units;
            if (!reconstruct_native_unit_log_generators(
                        expanded_units, *order,
                        std::span<const FactoredElement>(units->data(),
                                                         units->size()),
                        embeddings, diagnostics)) {
                return false;
            }
            units->swap(expanded_units);
            extraction_state.reduce_mod_units_log_generators_expanded = true;
        }
    }
    if (static_cast<slong>(units->size()) != rank) {
        units->clear();
        extraction_state.reduce_mod_units_log_generator_expansion_attempted =
                false;
        extraction_state.reduce_mod_units_log_generators_expanded = false;
        return false;
    }

    slong work_precision = precision;
    for (;;) {
        for (const auto& entry :
             extraction_state.reduce_mod_units_inverses) {
            if (entry.precision == work_precision && entry.rank == rank) {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.reduce_mod_units.cutoff_inverse_cache_hit");
                profile_reduce_precision(
                        diagnostics,
                        UnitReductionPrecisionEventKind::inverse_hit,
                        work_precision);
                precision = work_precision;
                out = &entry.inverse;
                return true;
            }
        }

        flint::ArbMat unit_logs(rank, rank);
        {
            profile_reduce_precision(
                    diagnostics,
                    UnitReductionPrecisionEventKind::inverse_work,
                    work_precision);
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_mod_units.unit_log_matrix");
            if (!unit_log_cutoff_matrix(
                        unit_logs, *units, embeddings, work_precision,
                        diagnostics)) {
                if (!extraction_state
                             .reduce_mod_units_log_generators_expanded) {
                    return false;
                }

                // Exact expansion can make a low-precision embedding
                // inconclusive through coefficient cancellation.  Preserve
                // the source representation and retry the same log matrix
                // through the original factored generators.
                units->clear();
                extraction_state.reduce_mod_units_log_generators_expanded =
                        false;
                if (!copy_unit_generators(*units, group) ||
                    !unit_log_cutoff_matrix(
                            unit_logs, *units, embeddings, work_precision,
                            diagnostics)) {
                    units->clear();
                    return false;
                }
            }
        }

        flint::ArbMat inverse(rank, rank);
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.reduce_mod_units.unit_log_inverse");
            if (arb_mat_inv(inverse.raw(), unit_logs.raw(),
                            work_precision) != 0) {
                extraction_state.reduce_mod_units_inverses.emplace_back(
                        work_precision, rank);
                arb_mat_set(
                        extraction_state.reduce_mod_units_inverses.back()
                                .inverse.raw(),
                        inverse.raw());
                precision = work_precision;
                out = &extraction_state.reduce_mod_units_inverses.back()
                               .inverse;
                return true;
            }
        }

        if (work_precision > max_precision / 2) {
            return true;
        }
        work_precision *= 2;
    }
}

bool reduce_relation_units_modulo(std::vector<FactoredElement>& values,
                            const OrderUnitGroup& group,
                            EmbeddingContext& embeddings,
                            slong start_precision,
                            RelationUnitExtractionState* extraction_state)
        noexcept {
    if (values.empty() || group.free_rank() <= 0) {
        return true;
    }
    if (!group.is_set() || embeddings.parent() == nullptr ||
        start_precision <= 0) {
        return false;
    }

    slong precision = start_precision;
    if (!reduce_mod_units_start_precision(
                precision, values, group, embeddings, start_precision)) {
        return false;
    }

    const slong max_precision =
            max_slong(precision, kReduceModUnitsMaxPrecision);
    for (; precision <= max_precision;) {
        const flint::ArbMat* unit_log_inverse = nullptr;
        if (extraction_state != nullptr) {
            if (!reduce_mod_units_cutoff_inverse(
                        unit_log_inverse, *extraction_state, group,
                        embeddings, precision, max_precision)) {
                return false;
            }
        }

        flint::FmpzMat coefficients(
                static_cast<slong>(values.size()), group.free_rank());
        bool exact = false;
        if (reduction_coordinates(coefficients, values, group,
                                        embeddings, exact, precision,
                                        unit_log_inverse)) {
            if (fmpz_matrix_is_zero(coefficients)) {
                return true;
            }
            if (!apply_unit_reduction(values, group, coefficients)) {
                return false;
            }
            if (exact) {
                return true;
            }
            continue;
        }
        if (precision > max_precision / 2) {
            break;
        }
        precision *= 2;
    }

    return true;
}

bool reduce_relation_units_modulo(FactoredElement& value,
                            const OrderUnitGroup& group,
                            EmbeddingContext& embeddings,
                            slong start_precision,
                            RelationUnitExtractionState* extraction_state =
                                    nullptr) noexcept {
    std::vector<FactoredElement> values;
    const NumberField* field = value.parent();
    if (field == nullptr) {
        return false;
    }
    values.emplace_back(*field);
    if (!values.back().set(value) ||
        !reduce_relation_units_modulo(values, group, embeddings, start_precision,
                                extraction_state)) {
        return false;
    }
    return value.set(values.front());
}

bool conj_log_cutoff_inverse(const flint::ArbMat*& out,
                                   RelationUnitExtractionState& extraction_state,
                                   const OrderUnitGroup& group,
                                   EmbeddingContext& embeddings,
                                   slong& precision) noexcept {
    out = nullptr;
    const slong rank = group.free_rank();
    if (!group.is_set() || rank <= 0 || precision <= 0 ||
        embeddings.parent() == nullptr) {
        return false;
    }

    if (extraction_state.cutoff_inverse.has_value() &&
        extraction_state.cutoff_inverse_precision == precision &&
        extraction_state.cutoff_inverse_rank == rank) {
        out = &*extraction_state.cutoff_inverse;
        return true;
    }

    std::vector<FactoredElement> units;
    if (!copy_unit_generators(units, group)) {
        return false;
    }

    slong work_precision = precision;
    for (;;) {
        slong places = 0;
        if (!compact_places(places, embeddings) || places != rank + 1) {
            return false;
        }

        flint::ArbMat logs(rank, places);
        if (!compact_log_matrix(
                    logs, embeddings,
                    FactoredElementSpan(units.data(), units.size()),
                    work_precision)) {
            return false;
        }

        flint::ArbMat cutoff(rank, rank);
        for (slong i = 0; i < rank; ++i) {
            for (slong j = 0; j < rank; ++j) {
                arb_set(arb_mat_entry(cutoff.raw(), i, j),
                        arb_mat_entry(logs.raw(), i, j));
            }
        }

        flint::ArbMat inverse(rank, rank);
        if (arb_mat_inv(inverse.raw(), cutoff.raw(), work_precision) != 0) {
            extraction_state.cutoff_inverse.emplace(rank, rank);
            arb_mat_set(extraction_state.cutoff_inverse->raw(),
                        inverse.raw());
            extraction_state.cutoff_inverse_precision = work_precision;
            extraction_state.cutoff_inverse_rank = rank;
            precision = work_precision;
            out = &*extraction_state.cutoff_inverse;
            return true;
        }

        if (work_precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        work_precision *= 2;
    }
}

bool search_dependent_relation(
        bool& recovered,
        FactoredElement& root,
        flint::FmpzMat& relation,
        flint::Fmpz& torsion_exp,
        RelationUnitExtractionState& extraction_state,
        const OrderUnitGroup& group,
        const FactoredElement& candidate,
        EmbeddingContext& embeddings,
        flint::FmpzConstRef denominator_bound,
        slong& precision) noexcept {
    recovered = false;
    if (precision <= 0) {
        return false;
    }

    slong work_precision = precision;
    for (;;) {
        const flint::ArbMat* inverse = nullptr;
        if (!conj_log_cutoff_inverse(inverse, extraction_state, group,
                                           embeddings, work_precision) ||
            inverse == nullptr) {
            return false;
        }

        if (!dependent_relation_bounded_with_inverse(
                    recovered, root, relation, torsion_exp, group, candidate,
                    embeddings, *inverse, denominator_bound, 1,
                    work_precision, true, false)) {
            return false;
        }
        if (recovered) {
            precision = work_precision;
            return true;
        }

        if (work_precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        work_precision *= 2;
    }
}

bool add_dependent_unit(bool& changed,
                              OrderUnitGroup& group,
                              const FactoredElement& candidate,
                              EmbeddingContext& embeddings,
                              RelationUnitExtractionState& extraction_state,
                              slong precision) noexcept {
    changed = false;
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || rank <= 0 ||
        candidate.parent() == nullptr ||
        !candidate.parent()->has_same_data(*field) ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) || precision <= 0) {
        return false;
    }

    flint::Fmpz denominator_bound;
    if (!group.regulator_index_bound(flint::FmpzRef(denominator_bound),
                                     precision) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(denominator_bound), 1) <= 0 ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(denominator_bound))) {
        return true;
    }

    FactoredElement root(*field);
    flint::FmpzMat relation(1, rank + 1);
    flint::Fmpz torsion_exp;
    bool recovered = false;
    slong relation_precision = extraction_state.rel_add_precision;
    if (!search_dependent_relation(
                recovered, root, relation, torsion_exp, extraction_state,
                group, candidate, embeddings,
                flint::FmpzConstRef(denominator_bound),
                relation_precision)) {
        return true;
    }
    if (!recovered) {
        return true;
    }
    extraction_state.rel_add_precision = relation_precision;

    OrderUnitGroup refined(*order);
    refined.set_diagnostics(group.diagnostics());
    if (!refined.is_defined() ||
        !adjoin_verified_dependent_relation(
                changed, refined, group, root,
                flint::FmpzMatConstRef(relation), 0, embeddings,
                precision)) {
        changed = false;
        return true;
    }
    if (changed) {
        group.swap(refined);
        bool stored_units_reduced = false;
        if (!reduce_stored_relation_units(stored_units_reduced, group,
                                       embeddings, relation_precision)) {
            changed = false;
            return false;
        }
        extraction_state.clear_dependent_unit_cache();
    }
    return true;
}

void validation_stop_after_improvement(
        bool& stop,
        RelationUnitExtractionState& extraction_state,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision,
        ulong validation_bf_max_cutoff) noexcept {
    stop = false;
    if (!extraction_state.has_expected_regulator) {
        if (validation_bf_max_cutoff != 0 && order.degree() > 2) {
            extraction_state.has_expected_regulator =
                    class_unit_bf_validation_estimate(
                            extraction_state.validation_index_bound,
                            extraction_state.expected_regulator, units,
                            class_group, analytic_cache, order,
                            validation_bf_max_cutoff, precision,
                            units.diagnostics());
        } else {
            extraction_state.has_expected_regulator =
                    class_unit_validation_estimate(
                            extraction_state.validation_index_bound,
                            extraction_state.expected_regulator, units,
                            class_group, analytic_cache, order, precision);
        }
    }
    if (extraction_state.has_expected_regulator &&
        expected_regulator_stop(
                units, extraction_state.expected_regulator, precision)) {
        stop = true;
    }
}

bool unknown_candidate_validated(
        bool& validated,
        RelationUnitExtractionState& extraction_state,
        ValidateRefineSummary& summary,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        const Order& order,
        AnalyticClassRegulatorCache& analytic_cache,
        slong expected_rank,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(units.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.unknown_validation");
    validated = false;
    summary.reset();
    if (precision <= 0 || expected_rank < 0 || !units.is_set() ||
        !same_order_parent(units.parent(), &order) ||
        !class_group.has_presentation() ||
        !same_order_parent(class_group.parent(), &order) ||
        units.free_rank() != expected_rank) {
        summary.outcome = ValidateRefineOutcome::local_rank_unavailable;
        return true;
    }

    flint::Fmpz class_order;
    if (!class_group.order(flint::FmpzRef(class_order))) {
        summary.outcome = ValidateRefineOutcome::local_rank_unavailable;
        return true;
    }

    if ((order.degree() == 1 || expected_rank == 0) &&
        flint::fmpz_is_one(flint::FmpzConstRef(class_order))) {
        flint::Arb regulator;
        if (!units.regulator(flint::ArbRef(regulator))) {
            summary.outcome = ValidateRefineOutcome::local_rank_unavailable;
            return true;
        }
        if (flint::arb_contains_si(regulator, 1)) {
            flint::fmpz_one(flint::FmpzRef(summary.last_index_bound));
            flint::fmpz_one(
                    flint::FmpzRef(extraction_state.validation_index_bound));
            flint::arb_abs(extraction_state.expected_regulator, regulator);
            extraction_state.has_expected_regulator = true;
            summary.outcome = ValidateRefineOutcome::proven;
            validated = true;
            return true;
        }
        if (order.degree() == 1) {
            summary.outcome = ValidateRefineOutcome::no_progress;
            return true;
        }
    }

    extraction_state.has_expected_regulator =
            class_unit_validation_estimate(
                    extraction_state.validation_index_bound,
                    extraction_state.expected_regulator, units, class_group,
                    analytic_cache, order, precision);
    if (!extraction_state.has_expected_regulator) {
        summary.outcome = ValidateRefineOutcome::analytic_unavailable;
        return true;
    }

    flint::fmpz_set(flint::FmpzRef(summary.last_index_bound),
                    flint::FmpzConstRef(
                            extraction_state.validation_index_bound));
    if (flint::fmpz_is_one(
                flint::FmpzConstRef(
                        extraction_state.validation_index_bound))) {
        summary.outcome = ValidateRefineOutcome::proven;
        validated = true;
    } else {
        summary.outcome = ValidateRefineOutcome::no_progress;
    }
    return true;
}

bool class_unit_index_bound_for_improvement(
        flint::Fmpz& out,
        const OrderUnitGroup& units,
        const ClassGroupContext& class_group,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong precision) noexcept {
    flint::Arb expected_regulator;
    return class_unit_validation_estimate(
            out, expected_regulator, units, class_group, analytic_cache, order,
            precision);
}

bool unit_improvement_ratio(slong& out,
                                  flint::FmpzConstRef starting_index,
                                  const OrderUnitGroup& units,
                                  const ClassGroupContext& class_group,
                                  AnalyticClassRegulatorCache& analytic_cache,
                                  const Order& order,
                                  slong precision) noexcept {
    out = 0;
    flint::Fmpz ending_index;
    if (!class_unit_index_bound_for_improvement(
                ending_index, units, class_group, analytic_cache, order,
                precision) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(ending_index))) {
        return false;
    }

    flint::Fmpz ratio;
    flint::fmpz_fdiv_q(flint::FmpzRef(ratio), starting_index,
                       flint::FmpzConstRef(ending_index));
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(ratio))) {
        return false;
    }
    out = flint::fmpz_get_si(flint::FmpzConstRef(ratio));
    return out >= 0;
}

bool next_prime_slong(slong& next, slong current) noexcept {
    if (current < 2) {
        next = 2;
        return true;
    }

    flint::Fmpz value;
    flint::fmpz_set_si(flint::FmpzRef(value), current);
    flint::fmpz_nextprime(flint::FmpzRef(value),
                          flint::FmpzConstRef(value));
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(value))) {
        return false;
    }
    next = flint::fmpz_get_si(flint::FmpzConstRef(value));
    return next > current;
}

bool unknown_saturation_at_two(
        bool& validated,
        RelationUnitExtractionState& extraction_state,
        ValidateRefineSummary& summary,
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong expected_rank,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.saturation_at_two");
    validated = false;
    if (extraction_state.saturation_at_two_done || precision <= 0 ||
        expected_rank < 0 || !units.is_set() ||
        units.free_rank() != expected_rank ||
        !same_order_parent(units.parent(), &order) ||
        !class_group.has_presentation() ||
        !same_order_parent(class_group.parent(), &order)) {
        return true;
    }

    constexpr double kSaturationAtTwoStable = 3.5;
    for (;;) {
        bool saturated = false;
        SaturationCandidateProcessingResult last_result;
        if (!saturate_class_unit_context(
                    saturated, last_result, class_group, units, embeddings,
                    extraction_state, 2, kSaturationAtTwoStable,
                    precision)) {
            return false;
        }
        if (!saturated) {
            break;
        }
    }

    extraction_state.saturation_at_two_done = true;
    return unknown_candidate_validated(
            validated, extraction_state, summary, units, class_group, order,
            analytic_cache, expected_rank, precision);
}

bool unknown_small_index_saturation(
        bool& validated,
        bool& attempted,
        RelationUnitExtractionState& extraction_state,
        ValidateRefineSummary& summary,
        ClassGroupContext& class_group,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        AnalyticClassRegulatorCache& analytic_cache,
        const Order& order,
        slong expected_rank,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.small_index_saturation");
    validated = false;
    attempted = false;
    if (precision <= 0 || expected_rank < 0 || !units.is_set() ||
        units.free_rank() != expected_rank ||
        !same_order_parent(units.parent(), &order) ||
        !class_group.has_presentation() ||
        !same_order_parent(class_group.parent(), &order)) {
        return true;
    }
    if (!flint::fmpz_fits_si(
                flint::FmpzConstRef(summary.last_index_bound))) {
        return true;
    }

    slong index_bound = flint::fmpz_get_si(
            flint::FmpzConstRef(summary.last_index_bound));
    constexpr double kSmallIndexStable = 3.5;
    while (index_bound > 1 && index_bound < 20) {
        attempted = true;
        bool saturated = false;
        SaturationCandidateProcessingResult last_result;
        for (slong p = 2; !saturated && p < 2 * index_bound;) {
            if (!saturate_class_unit_context(
                        saturated, last_result, class_group, units,
                        embeddings, extraction_state, p,
                        kSmallIndexStable, precision)) {
                return false;
            }
            if (!saturated && !next_prime_slong(p, p)) {
                return false;
            }
        }

        if (!saturated) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "small-index saturation source assertion failed");
            return false;
        }

        bool revalidated = false;
        if (!unknown_candidate_validated(
                    revalidated, extraction_state, summary, units,
                    class_group, order, analytic_cache, expected_rank,
                    precision)) {
            return false;
        }
        if (revalidated) {
            validated = true;
            return true;
        }
        if (!flint::fmpz_fits_si(
                    flint::FmpzConstRef(summary.last_index_bound))) {
            return false;
        }

        const slong next_index_bound = flint::fmpz_get_si(
                flint::FmpzConstRef(summary.last_index_bound));
        if (next_index_bound <= 1) {
            validated = true;
            return true;
        }
        if (next_index_bound >= index_bound) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "small-index saturation index did not decrease");
            return false;
        }
        index_bound = next_index_bound;
    }

    return true;
}

bool improve_full_rank_class_relation_units(
        bool& improved,
        OrderUnitGroup& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        AnalyticClassRegulatorCache& analytic_cache,
        RelationUnitExtractionState& extraction_state,
        slong precision,
        slong rank,
        slong outer_add,
        ulong validation_bf_max_cutoff) noexcept {
    improved = false;
    if (!out.is_set() || !same_order_parent(out.parent(), &order) ||
        out.free_rank() != rank || rank <= 0 || outer_add < 0 ||
        outer_add > WORD_MAX - 2) {
        return false;
    }

    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    RelationUnitCandidateBatchState& batch_state =
            extraction_state.unit_candidate_batch;
    slong add = outer_add + 2;
    slong new_add = 0;
    slong not_larger = 0;
    for (;;) {
        if (new_add > 0 && add > WORD_MAX - new_add) {
            return false;
        }
        add += new_add;
        new_add = 2;

        std::vector<FactoredElement> candidates;
        if (!unit_candidate_witnesses(candidates, batch_state,
                                            class_group, add)) {
            return false;
        }
        const slong not_larger_bound =
                unit_not_larger_bound(batch_state, rank);
        if (not_larger_bound <= 0) {
            return true;
        }
        if (!reduce_relation_units_modulo(candidates, out, embeddings,
                                    extraction_state.torsion_precision,
                                    &extraction_state)) {
            return false;
        }
        if (!sort_unit_candidates_by_log_height(candidates,
                                                      embeddings)) {
            return false;
        }
        if (candidates.empty()) {
            ++not_larger;
        }

        bool finished = false;
        for (FactoredElement& candidate : candidates) {
            // reference `_unit_group_find_units` only verifies these kernel
            // elements are units under `@hassert`; production follows the
            // relation-kernel construction and goes straight to reduction and
            // torsion testing.
            if (!reduce_relation_units_modulo(
                        candidate, out, embeddings,
                        extraction_state.torsion_precision,
                        &extraction_state)) {
                return false;
            }

            RelationTorsionStatus torsion_status =
                    RelationTorsionStatus::inconclusive;
            if (!unit_candidate_torsion_status(
                        torsion_status, extraction_state, candidate,
                        embeddings)) {
                return false;
            }
            if (torsion_status != RelationTorsionStatus::non_torsion) {
                ++not_larger;
                continue;
            }

            bool changed = false;
            if (!add_dependent_unit(changed, out, candidate,
                                          embeddings, extraction_state,
                                          precision)) {
                return false;
            }
            if (changed) {
                improved = true;
                not_larger = 0;
                bool stop = false;
                validation_stop_after_improvement(
                        stop, extraction_state, out, class_group,
                        analytic_cache, order, precision,
                        validation_bf_max_cutoff);
                if (stop) {
                    finished = true;
                    break;
                }
            } else {
                ++not_larger;
            }
        }

        bool stored_units_reduced = false;
        if (!reduce_stored_relation_units(
                    stored_units_reduced, out, embeddings,
                    extraction_state.torsion_precision)) {
            return false;
        }
        if (stored_units_reduced) {
            extraction_state.clear_dependent_unit_cache();
        }
        if (finished) {
            return true;
        }

        if (not_larger > not_larger_bound &&
            unit_candidate_has_unused_relations(batch_state)) {
            not_larger = 0;
        }
        if (not_larger >= not_larger_bound &&
            !unit_candidate_has_unused_relations(batch_state)) {
            break;
        }
    }

    return true;
}

bool set_class_relation_units(bool& ready,
                                    slong& improved,
                                    OrderUnitGroup& out,
                                    const Order& order,
                                    const ClassGroupContext& class_group,
                                    EmbeddingContext& embeddings,
                                    AnalyticClassRegulatorCache& analytic_cache,
                                    RelationUnitExtractionState& extraction_state,
                                    slong precision,
                                    slong rank,
                                    slong outer_add,
                                    ulong validation_bf_max_cutoff) noexcept {
    ready = false;
    improved = 0;
    const NumberField* field = order.parent();
    if (field == nullptr || rank < 0 || outer_add < 0) {
        return false;
    }
    if (rank == 0) {
        ready = out.compute(order);
        return ready;
    }

    if (out.is_set() && same_order_parent(out.parent(), &order) &&
        out.free_rank() == rank) {
        flint::Fmpz starting_index;
        if (!class_unit_index_bound_for_improvement(
                    starting_index, out, class_group, analytic_cache, order,
                    precision)) {
            return false;
        }

        bool changed = false;
        if (!improve_full_rank_class_relation_units(
                    changed, out, order, class_group, embeddings,
                    analytic_cache, extraction_state, precision, rank,
                    outer_add, validation_bf_max_cutoff)) {
            return false;
        }
        (void)changed;
        if (!unit_improvement_ratio(
                    improved, flint::FmpzConstRef(starting_index), out,
                    class_group, analytic_cache, order, precision)) {
            return false;
        }
        ready = true;
        return true;
    }

    RelationUnitCandidateBatchState& batch_state =
            extraction_state.unit_candidate_batch;
    IncrementalUnitContext partial_units;
    if (!partial_units.define(*field, rank)) {
        return false;
    }
    slong add = outer_add;
    slong new_add = 0;
    slong not_larger = 0;
    for (;;) {
        if (new_add > 0 && add > WORD_MAX - new_add) {
            return false;
        }
        add += new_add;
        new_add = 2;

        std::vector<FactoredElement> candidates;
        if (!unit_candidate_witnesses(candidates, batch_state,
                                            class_group, add)) {
            return false;
        }
        const slong not_larger_bound =
                unit_not_larger_bound(batch_state, rank);
        if (not_larger_bound <= 0) {
            return true;
        }
        if (!sort_unit_candidates_by_log_height(candidates,
                                                      embeddings)) {
            return false;
        }
        std::vector<char> done(candidates.size(), 0);
        if (candidates.empty()) {
            ++not_larger;
        }

        bool finished = false;
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            FactoredElement& candidate = candidates[i];
            // reference only checks `_isunit(candidate)` through `@hassert`;
            // normal unit search trusts the relation-kernel witnesses.
            if (out.is_set() && same_order_parent(out.parent(), &order) &&
                out.free_rank() == rank &&
                !reduce_relation_units_modulo(
                        candidate, out, embeddings,
                        extraction_state.torsion_precision,
                        &extraction_state)) {
                return false;
            }

            RelationTorsionStatus torsion_status =
                    RelationTorsionStatus::inconclusive;
            if (!unit_candidate_torsion_status(
                        torsion_status, extraction_state, candidate,
                        embeddings)) {
                return false;
            }
            if (torsion_status != RelationTorsionStatus::non_torsion) {
                done[i] = 1;
                ++not_larger;
                continue;
            }

            if (out.is_set() && same_order_parent(out.parent(), &order) &&
                out.free_rank() == rank) {
                bool changed = false;
                if (!add_dependent_unit(
                            changed, out, candidate, embeddings,
                            extraction_state, precision)) {
                    return false;
                }
                done[i] = 1;
                if (changed) {
                    not_larger = 0;
                    bool stop = false;
                    validation_stop_after_improvement(
                            stop, extraction_state, out, class_group,
                            analytic_cache, order, precision,
                            validation_bf_max_cutoff);
                    if (stop) {
                        finished = true;
                        break;
                    }
                } else {
                    ++not_larger;
                }
                continue;
            }

            bool added = false;
            if (!partial_units.add_independent_unit(added, candidate,
                                                    embeddings, precision)) {
                return false;
            }
            if (!added) {
                ++not_larger;
                continue;
            }
            done[i] = 1;
            not_larger = 0;
            if (partial_units.full_rank) {
                OrderUnitGroup full_rank_units(order);
                full_rank_units.set_diagnostics(out.diagnostics());
                bool stored_units_reduced = false;
                if (!full_rank_units.is_defined() ||
                    !order_unit_group_set_units_internal(
                            full_rank_units, order,
                            FactoredElementSpan(partial_units.units.data(),
                                                partial_units.units.size()),
                            // reference trusts relation-kernel candidates here
                            // after `add_unit!` has accepted independence; its
                            // `_isunit` check is only an assertion.
                            embeddings, precision, true) ||
                    !reduce_stored_relation_units(
                            stored_units_reduced, full_rank_units, embeddings,
                            extraction_state.torsion_precision)) {
                    return false;
                }
                if (stored_units_reduced) {
                    extraction_state.clear_dependent_unit_cache();
                }
                out.swap(full_rank_units);
                bool stop = false;
                validation_stop_after_improvement(
                        stop, extraction_state, out, class_group,
                        analytic_cache, order, precision,
                        validation_bf_max_cutoff);
                if (stop) {
                    finished = true;
                    break;
                }
            }
        }

        if (out.is_set() && same_order_parent(out.parent(), &order) &&
            out.free_rank() == rank) {
            for (std::size_t i = 0; !finished && i < candidates.size(); ++i) {
                if (done[i] != 0) {
                    continue;
                }

                if (!reduce_relation_units_modulo(
                            candidates[i], out, embeddings,
                            extraction_state.torsion_precision,
                            &extraction_state)) {
                    return false;
                }

                RelationTorsionStatus torsion_status =
                        RelationTorsionStatus::inconclusive;
                if (!unit_candidate_torsion_status(
                            torsion_status, extraction_state, candidates[i],
                            embeddings)) {
                    return false;
                }
                if (torsion_status != RelationTorsionStatus::non_torsion) {
                    continue;
                }

                bool changed = false;
                if (!add_dependent_unit(
                            changed, out, candidates[i], embeddings,
                            extraction_state, precision)) {
                    return false;
                }
                if (changed) {
                    bool stop = false;
                    validation_stop_after_improvement(
                            stop, extraction_state, out, class_group,
                            analytic_cache, order, precision,
                            validation_bf_max_cutoff);
                    if (stop) {
                        finished = true;
                        break;
                    }
                }
            }

            bool stored_units_reduced = false;
            if (!reduce_stored_relation_units(
                        stored_units_reduced, out, embeddings,
                        extraction_state.torsion_precision)) {
                return false;
            }
            if (stored_units_reduced) {
                extraction_state.clear_dependent_unit_cache();
            }
            if (finished) {
                break;
            }
        }

        if (not_larger > not_larger_bound &&
            unit_candidate_has_unused_relations(batch_state)) {
            not_larger = 0;
        }
        if (not_larger >= not_larger_bound &&
            !unit_candidate_has_unused_relations(batch_state)) {
            break;
        }
    }

    ready = out.is_set() && same_order_parent(out.parent(), &order) &&
            out.free_rank() == rank;
    return true;
}

}  // namespace detail

using detail::AnalyticClassRegulatorCache;
using detail::ValidateRefineSummary;
using detail::class_group_order_is_one;
using detail::initial_relation_kernel_target;
using detail::kComputeAnalyticBumpMaxPasses;
using detail::kComputePostFinitePhaseBudget;
using detail::kDefaultComputePrecision;
using detail::kRankOneAdaptiveSurplus;
using detail::max_slong;
using detail::missing_unit_target_bump;
using detail::rank_zero_torsion;
using detail::relation_kernel_independent_unit_count;
using detail::recompute_units_from_class_context;
using detail::set_class_relation_units;
using detail::set_hnf_units;
using detail::try_validate_refine_extra_pass;
using detail::try_validate_refine_loop;
using detail::validation_target_bump;

bool OrderUnitGroup::compute(const Order& order) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.compute");
    if (!order.is_defined() || order.parent() == nullptr || !order.has_basis()) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }

    if (rank != 0) {
        // `quadratic_fundamental_unit` accepts both the explicit quadratic
        // backend and a canonical polynomial-defined x^2-d field.
        if (rank != 1 || !order.is_maximal() ||
            order.degree() != 2) {
            return false;
        }

        Element epsilon(*order.parent());
        FactoredElement generator(*order.parent());
        EmbeddingContext embeddings(*order.parent());
        OrderUnitGroup candidate(order);
        candidate.set_diagnostics(diagnostics_);
        if (!candidate.is_defined() ||
            !quadratic_fundamental_unit(epsilon, *order.parent()) ||
            !generator.set_element(epsilon) ||
            !candidate.set_units(order, FactoredElementSpan(&generator, 1),
                                 embeddings, kDefaultComputePrecision)) {
            return false;
        }
        candidate.certification_ = CertificationMode::proven;
        swap(candidate);
        return true;
    }

    OrderUnitGroup candidate(order);
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined() ||
        !rank_zero_torsion(flint::FmpzRef(candidate.torsion_order_),
                           candidate.torsion_generator_, order)) {
        return false;
    }

    candidate.clear_free_generators_();
    flint::arb_one(candidate.regulator_);
    candidate.has_regulator_ = true;
    candidate.certification_ = CertificationMode::proven;
    candidate.is_set_ = true;

    swap(candidate);
    return true;
}

bool OrderUnitGroup::compute_with_relation_class_group_(
        ClassGroupContext& class_group,
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupComputeOptions& options,
        slong precision,
        slong rank,
        const DiagnosticsContext* active_diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.compute_with_relation_class_group");

    const bool transaction =
            detail::uses_class_unit_kernel(
                    class_group.class_unit_transaction_context_);
    if (options.requested_certification == CertificationMode::grh &&
        !transaction) {
        return false;
    }

    detail::ClassGroupRelationOptions local_options =
            detail::class_group_relation_options(order, options);
    local_options.diagnostics = active_diagnostics;
    // reference `class_group_ctx` is tentative in every degree; unconditional
    // proof belongs to the enclosing `_class_unit_group` computation.  Keep
    // the legacy selector's specialized quadratic request, but do not require
    // quadratic-backend metadata at this strict source-route boundary.
    local_options.requested_certification =
            order.degree() == 2 && !transaction
            ? options.requested_certification
            : CertificationMode::unknown;

    ClassGroupContext candidate_class_group;
    candidate_class_group.class_unit_transaction_context_ =
            class_group.class_unit_transaction_context_;
    candidate_class_group.set_diagnostics(active_diagnostics);
    if (!candidate_class_group.compute_tentative_candidate_(
                order, factor_base_bound, local_options)) {
        return false;
    }
    if (options.requested_certification == CertificationMode::proven ||
        (options.requested_certification == CertificationMode::grh &&
         candidate_class_group.factor_base_generation_status() ==
                 ProofState::verified)) {
        flint::Fmpz required_bound;
        if (!candidate_class_group.factor_base_generation_bound(
                    flint::FmpzRef(required_bound)) ||
            !candidate_class_group.check_factor_base_generation_bound(
                    flint::FmpzConstRef(required_bound))) {
            return false;
        }
    }

    OrderUnitGroup candidate_units(order);
    OrderUnitGroup scratch_units(order);
    candidate_units.set_diagnostics(active_diagnostics);
    scratch_units.set_diagnostics(active_diagnostics);
    EmbeddingContext embeddings(*order.parent());
    if (!candidate_units.is_defined() || !scratch_units.is_defined()) {
        return false;
    }

    AnalyticClassRegulatorCache analytic_cache;
    detail::RelationUnitExtractionState relation_unit_state;
    slong unit_add = 0;
    slong unit_improved = 0;
    detail::ValidationUnitRefreshOptions refresh_options;
    if (rank > 1 || (transaction && rank > 0)) {
        refresh_options.relation_unit_state = &relation_unit_state;
        refresh_options.unit_add = &unit_add;
        refresh_options.unit_improved = &unit_improved;
        refresh_options.rank = rank;
        refresh_options.validation_bf_max_cutoff =
                detail::validation_bf_max_cutoff_for_options(order, options);
        refresh_options.use_class_relation_units = true;
        refresh_options.require_source_units =
                transaction;
    }
    bool do_units = rank == 0 ||
            candidate_class_group.relation_kernel_unit_count() >= rank;
    for (;;) {
        bool candidate_units_ready = false;
        bool certified = false;
        bool rank_complete = false;
        ValidateRefineSummary validate_summary;
        slong revalidate_passes_done = 0;

        if (rank == 0) {
            SILEX_PROFILE_EVENT(active_diagnostics,
                                DiagnosticsModule::unit_group,
                                "unit_group.rank_zero_units");
            if (!candidate_units.compute(order)) {
                return false;
            }
            candidate_units_ready = true;
        } else if (order.degree() == 2 && !transaction &&
                   candidate_units.compute(order)) {
            candidate_units_ready = true;
        } else if (options.requested_certification ==
                           CertificationMode::proven &&
                   refresh_options.use_class_relation_units) {
            SILEX_PROFILE_EVENT(active_diagnostics,
                                DiagnosticsModule::unit_group,
                                "unit_group.class_relation_units");
            if (!set_class_relation_units(
                        candidate_units_ready, unit_improved,
                        candidate_units, order, candidate_class_group,
                        embeddings, analytic_cache, relation_unit_state,
                        precision, rank, unit_add,
                        refresh_options.validation_bf_max_cutoff)) {
                return false;
            }
            if (!candidate_units_ready) {
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "class/unit unit extraction did not yet "
                          "produce a full-rank unit group");
            }
            if (unit_add == WORD_MAX) {
                return false;
            }
            ++unit_add;
        } else if (options.requested_certification ==
                   CertificationMode::proven) {
            SILEX_PROFILE_EVENT(active_diagnostics,
                                DiagnosticsModule::unit_group,
                                "unit_group.recompute_units");
            if (!recompute_units_from_class_context(
                        candidate_units, order, candidate_class_group,
                        embeddings, analytic_cache, true,
                        detail::ValidationRecomputeCause::
                                initial_requested_proven,
                        precision)) {
                return false;
            }
            candidate_units_ready = true;
        } else if (!do_units) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "class/unit loop skipped unit extraction");
            if (candidate_units.is_set() &&
                same_order_parent(candidate_units.parent(), &order) &&
                candidate_units.free_rank() == rank) {
                // reference's do_units flag only skips _unit_group_find_units;
                // an already full-rank unit context is still validated against
                // the extended class context in the same loop iteration.
                candidate_units_ready = true;
            }
        } else {
            SILEX_PROFILE_EVENT(active_diagnostics,
                                DiagnosticsModule::unit_group,
                                "unit_group.class_relation_units");
            if (!set_class_relation_units(
                        candidate_units_ready, unit_improved,
                        candidate_units, order, candidate_class_group,
                        embeddings, analytic_cache, relation_unit_state,
                        precision, rank, unit_add,
                        refresh_options.validation_bf_max_cutoff)) {
                return false;
            }
            if (!candidate_units_ready) {
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "class/unit unit extraction did not yet "
                          "produce a full-rank unit group");
            }
            if (unit_add == WORD_MAX) {
                return false;
            }
            ++unit_add;
        }
        auto validate_candidate = [&]() noexcept {
            rank_complete = candidate_units.free_rank() == rank;
            if (options.requested_certification == CertificationMode::grh) {
                if (try_validate_refine_loop(
                            candidate_class_group, candidate_units,
                            scratch_units, order, options, embeddings,
                            analytic_cache, validate_summary, precision,
                            false, refresh_options)) {
                    return true;
                }
                return try_validate_refine_extra_pass(
                        revalidate_passes_done, candidate_class_group,
                        candidate_units, scratch_units, order, options,
                        embeddings, analytic_cache, validate_summary,
                        precision, false, refresh_options);
            }
            if (options.requested_certification != CertificationMode::proven) {
                bool validated = false;
                if (!detail::unknown_candidate_validated(
                            validated, relation_unit_state,
                            validate_summary, candidate_units,
                            candidate_class_group, order, analytic_cache,
                            rank, precision)) {
                    return false;
                }
                return validated;
            }
            if (try_validate_refine_loop(
                        candidate_class_group, candidate_units, scratch_units,
                        order, options, embeddings, analytic_cache,
                        validate_summary, precision, false,
                        refresh_options)) {
                return true;
            }
            return try_validate_refine_extra_pass(
                    revalidate_passes_done, candidate_class_group,
                    candidate_units, scratch_units, order, options, embeddings,
                    analytic_cache, validate_summary, precision, false,
                    refresh_options);
        };

        if (candidate_units_ready && validate_candidate()) {
            certified = true;
        }

        if (!certified && candidate_units_ready && rank_complete &&
            options.requested_certification == CertificationMode::unknown) {
            if (unit_improved == 1 &&
                !relation_unit_state.saturation_at_two_done) {
                bool saturation_at_two_validated = false;
                if (!detail::unknown_saturation_at_two(
                            saturation_at_two_validated, relation_unit_state,
                            validate_summary, candidate_class_group,
                            candidate_units, embeddings, analytic_cache,
                            order, rank, precision, active_diagnostics)) {
                    return false;
                }
                if (saturation_at_two_validated) {
                    certified = true;
                }
            }
        }

        if (!certified && candidate_units_ready && rank_complete &&
            options.requested_certification == CertificationMode::unknown) {
            bool saturation_validated = false;
            bool saturation_attempted = false;
            if (!detail::unknown_small_index_saturation(
                        saturation_validated, saturation_attempted,
                        relation_unit_state, validate_summary,
                        candidate_class_group, candidate_units, embeddings,
                        analytic_cache, order, rank, precision,
                        active_diagnostics)) {
                return false;
            }
            if (saturation_validated) {
                certified = true;
            } else if (saturation_attempted) {
                return false;
            }
        }

        if (certified) {
            if (!detail::complete_requested_proven_relation_saturation(
                        candidate_class_group, candidate_units, order,
                        options, analytic_cache, precision)) {
                certified = false;
            }
        }

        if (certified &&
            options.requested_certification == CertificationMode::grh) {
            // reference Clgp.jl:_class_unit_group publishes this completed pair
            // after analytic index-one validation and skips only the later
            // unconditional class- and unit-group proof passes.
            candidate_class_group.certification_ = CertificationMode::grh;
            candidate_units.certification_ = CertificationMode::grh;
        }

        if (certified) {
            class_group.swap(candidate_class_group);
            class_group.set_diagnostics(active_diagnostics);
            swap(candidate_units);
            set_diagnostics(active_diagnostics);
            SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::unit_group,
                          VerboseLevel::progress,
                          "class/unit computation finished");
            return true;
        }

        const slong relation_count_before =
                candidate_class_group.relation_count();
        flint::Fmpz class_order_before;
        if (!candidate_class_group.order(flint::FmpzRef(class_order_before))) {
            return false;
        }
        detail::ClassGroupRelationOptions extend_options = local_options;
        if (!candidate_class_group.extend_tentative_relations_(
                    order, factor_base_bound, extend_options)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "class/unit continuation did not add a relation "
                      "slice");
            return false;
        }
        if (candidate_class_group.relation_count() <= relation_count_before) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "class/unit continuation made no relation "
                      "progress");
            return false;
        }
        flint::Fmpz class_order_after;
        if (!candidate_class_group.order(flint::FmpzRef(class_order_after))) {
            return false;
        }
        const bool class_order_unchanged =
                flint::fmpz_equal(flint::FmpzConstRef(class_order_before),
                                  flint::FmpzConstRef(class_order_after));
        if (class_order_unchanged) {
            do_units =
                    candidate_class_group.relation_kernel_unit_count() >= rank;
        } else {
            if (unit_add > WORD_MAX - 2) {
                return false;
            }
            unit_add += 2;
            do_units = false;
        }
    }
}

bool OrderUnitGroup::compute_with_class_group(
        ClassGroupContext& class_group,
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const ClassGroupComputeOptions& options,
        slong precision) noexcept {
    if (!detail::valid_paired_certification_request(
                options.requested_certification)) {
        return false;
    }

    if (class_group.class_unit_transaction_context_ == nullptr) {
        detail::ClassUnitTransactionReport audit;
        return detail::compute_class_unit_transaction(
                *this, class_group, order, factor_base_bound, options,
                precision, audit);
    }

    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.compute_with_class_group");
    SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::unit_group,
                  VerboseLevel::progress,
                  "class/unit computation started");
    if (!order.is_defined() || order.parent() == nullptr ||
        !order.has_basis() || !order.is_maximal() || precision <= 0 ||
        (class_group.is_defined() &&
         !same_order_parent(class_group.parent(), &order))) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }

    const bool paired_transaction =
            class_group.class_unit_transaction_context_ != nullptr;
    const detail::ClassUnitExecutionPolicy* policy = paired_transaction
            ? &class_group.class_unit_transaction_context_->audit.policy
            : nullptr;
    detail::NativeUnitStrategy expected_native_units =
            detail::NativeUnitStrategy::none;
    detail::NativeContinuationStrategy expected_native_continuation =
            detail::NativeContinuationStrategy::none;
    const bool expected_native_hnf =
            paired_transaction && rank > 0 &&
            policy->units == detail::NativeUnitStrategy::hnf;
    if (paired_transaction &&
        options.requested_certification == CertificationMode::grh &&
        order.degree() != 1) {
        expected_native_units = rank == 0
                ? detail::NativeUnitStrategy::rank_zero
                : detail::NativeUnitStrategy::class_relation;
        expected_native_continuation =
                detail::NativeContinuationStrategy::lll;
    } else if (expected_native_hnf) {
        expected_native_units = detail::NativeUnitStrategy::hnf;
        expected_native_continuation =
                detail::NativeContinuationStrategy::regulator_multiple_reconstruction;
    } else if (rank == 0) {
        expected_native_units = detail::NativeUnitStrategy::rank_zero;
        expected_native_continuation =
                detail::NativeContinuationStrategy::native_target_growth;
    } else if (order.degree() == 2 && rank == 1) {
        expected_native_units =
                detail::NativeUnitStrategy::quadratic_rank_one;
        expected_native_continuation =
                detail::NativeContinuationStrategy::native_target_growth;
    } else {
        expected_native_units =
                detail::NativeUnitStrategy::class_relation;
        expected_native_continuation =
                detail::NativeContinuationStrategy::lll;
    }
    if (paired_transaction &&
        (!policy->selected || policy->degree != order.degree() ||
         policy->unit_rank != rank ||
         policy->units != expected_native_units ||
         policy->continuation != expected_native_continuation ||
         policy->validation !=
                 detail::NativeValidationStrategy::shared_exact)) {
        return false;
    }

    if (paired_transaction && order.degree() == 1 &&
        options.requested_certification == CertificationMode::grh) {
        ClassGroupComputeOptions exact_options = options;
        exact_options.requested_certification = CertificationMode::proven;
        if (!compute_with_class_group(
                    class_group, order, factor_base_bound, exact_options,
                    precision) ||
            class_group.certification_ != CertificationMode::proven ||
            certification_ != CertificationMode::proven) {
            return false;
        }

        // The caller requested the weaker conditional contract.  Retain the
        // exact proof receipts produced above, but keep both public coarse
        // labels equal to the requested certification mode.
        class_group.certification_ = CertificationMode::grh;
        certification_ = CertificationMode::grh;
        return true;
    }
    if (policy->relations == detail::NativeRelationStrategy::lll) {
        return compute_with_relation_class_group_(
                class_group, order, factor_base_bound, options, precision,
                rank, active_diagnostics);
    }

    const slong lower_target = rank;
    const slong rank_one_surplus_target = rank == 1
            ? max_slong(lower_target, rank + kRankOneAdaptiveSurplus)
            : lower_target;
    slong target =
            initial_relation_kernel_target(rank, order.degree(), lower_target);
    bool rank_one_surplus_started =
            rank == 1 && target == rank_one_surplus_target &&
            target > lower_target;
    slong post_finite_refinement_phase_budget = 0;
    slong analytic_bumps_done = 0;
    AnalyticClassRegulatorCache analytic_cache;
    if (paired_transaction) {
        analytic_cache.configure_validation(
                detail::validation_bf_max_cutoff_for_options(order, options));
    }
    bool validation_cache_seeded_from_relation_factor_base_plan = false;
    detail::RelationUnitExtractionState relation_unit_state;
    relation_unit_state.expand_reduce_mod_units_log_generators =
            paired_transaction;
    slong unit_add = 0;
    slong unit_improved = 0;
    const bool native_hnf_unit_candidates =
            rank > 0 &&
            policy->units == detail::NativeUnitStrategy::hnf;
    const bool native_unit_candidates =
            rank > 1 && policy->units ==
                    detail::NativeUnitStrategy::class_relation;
    const bool hnf_unit_candidates =
            native_hnf_unit_candidates;
    detail::ValidationUnitRefreshOptions refresh_options;
    if (native_unit_candidates) {
        refresh_options.relation_unit_state = &relation_unit_state;
        refresh_options.unit_add = &unit_add;
        refresh_options.unit_improved = &unit_improved;
        refresh_options.rank = rank;
        refresh_options.validation_bf_max_cutoff =
                detail::validation_bf_max_cutoff_for_options(order, options);
        refresh_options.use_class_relation_units = true;
    }
    EmbeddingContext embeddings(*order.parent());

    if (order.degree() == 2 && rank == 1) {
        detail::ClassGroupRelationOptions local_options =
                detail::class_group_relation_options(order, options);
        local_options.diagnostics = active_diagnostics;
        const bool native_relation_completion =
                policy->relations ==
                        detail::NativeRelationStrategy::
                                relation_completion_table;
        // The independently proven quadratic unit is the published unit-rank
        // witness.  Do not require relation-kernel *surplus* here.  The
        // relation_completion-style analytic finish still gathers enough kernel witnesses
        // for its own unit-rank check before it certifies the class lattice.
        local_options.target_relation_kernel_units =
                native_relation_completion ? 0 : rank_one_surplus_target;
        // reference validates class/unit pairs after tentative class-group work;
        // do not require standalone quadratic certification before units exist.
        local_options.requested_certification = CertificationMode::unknown;

        ClassGroupContext candidate_class_group;
        OrderUnitGroup candidate_units(order);
        OrderUnitGroup scratch_units(order);
        candidate_class_group.class_unit_transaction_context_ =
                class_group.class_unit_transaction_context_;
        candidate_class_group.set_diagnostics(active_diagnostics);
        candidate_units.set_diagnostics(active_diagnostics);
        scratch_units.set_diagnostics(active_diagnostics);
        ValidateRefineSummary validate_summary;

        if (!candidate_units.is_defined() || !scratch_units.is_defined()) {
            return false;
        }
        if (candidate_class_group.compute_relation_candidate_(
                    order, factor_base_bound, local_options) &&
            candidate_units.compute(order) &&
            try_validate_refine_loop(
                    candidate_class_group, candidate_units, scratch_units,
                    order, options, embeddings, analytic_cache,
                    validate_summary, precision, false)) {
            class_group.swap(candidate_class_group);
            class_group.set_diagnostics(active_diagnostics);
            swap(candidate_units);
            set_diagnostics(active_diagnostics);
            SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::unit_group,
                          VerboseLevel::progress,
                          "class/unit computation finished");
            return true;
        }
    }

    for (;;) {
        bool compute_ok = false;
        bool certified = false;
        ValidateRefineSummary validate_summary;
        bool rank_complete = false;
        slong independent_unit_count = -1;

        detail::ClassGroupRelationOptions local_options =
                detail::class_group_relation_options(order, options);
        local_options.diagnostics = active_diagnostics;
        local_options.target_relation_kernel_units = target;
        local_options.post_finite_refinement_phase_budget =
                post_finite_refinement_phase_budget;
        // The selected native reference package completes reconstruct_regulator before its
        // reference-HNF unit publisher consumes the relation table.
        const bool native_low_degree_finish =
                order.degree() <= 4 && rank > 0 &&
                policy->relations ==
                        detail::NativeRelationStrategy::
                                relation_completion_prepass;
        const bool native_finish =
                native_hnf_unit_candidates ||
                native_low_degree_finish;
        class_group.class_unit_transaction_context_
                ->defer_relation_saturation_until_units =
                options.requested_certification == CertificationMode::proven &&
                ((options.zeta_bf_max_cutoff != 0 && order.degree() > 2 &&
                  rank == 1) || native_finish);
        // Combined class/unit proof is handled by validation with units.
        local_options.requested_certification = CertificationMode::unknown;

        ClassGroupContext candidate_class_group;
        OrderUnitGroup candidate_units(order);
        candidate_class_group.class_unit_transaction_context_ =
                class_group.class_unit_transaction_context_;
        candidate_class_group.set_diagnostics(active_diagnostics);
        candidate_units.set_diagnostics(active_diagnostics);
        bool candidate_units_from_current_class_context = false;
        bool native_do_units = true;
        if (!candidate_units.is_defined()) {
            return false;
        }
        SILEX_PROFILE_EVENT(active_diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.compute_class_group");
        compute_ok = candidate_class_group.compute_relation_candidate_(
                order, factor_base_bound, local_options);

        while (!native_unit_candidates && compute_ok && rank > 1 &&
               options.requested_certification != CertificationMode::proven) {
            slong current_unit_count = -1;
            slong missing_unit_bump = 0;
            if (!relation_kernel_independent_unit_count(
                        current_unit_count, order, candidate_class_group,
                        embeddings, precision) ||
                !missing_unit_target_bump(
                        missing_unit_bump, rank, current_unit_count, target,
                        options.max_relations)) {
                break;
            }

            target += missing_unit_bump;
            local_options.target_relation_kernel_units = target;
            SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "extending relation-kernel target for missing unit rank");
            compute_ok = candidate_class_group.extend_relation_kernel_units_(
                    order, factor_base_bound, local_options);
        }

        auto try_current_candidate = [&]() noexcept {
            if (!compute_ok) {
                return false;
            }
            certified = false;
            validate_summary = ValidateRefineSummary{};
            rank_complete = false;
            independent_unit_count = -1;
            candidate_units_from_current_class_context = false;
            bool candidate_units_ready = false;
            bool used_candidate_units = false;
            bool used_hnf_candidate_units = false;
            const bool requested_proven =
                    options.requested_certification ==
                    CertificationMode::proven;
            const bool cache_bf_residue_degrees =
                    requested_proven && options.zeta_bf_max_cutoff != 0 &&
                    order.degree() > 2;
            if (cache_bf_residue_degrees &&
                !validation_cache_seeded_from_relation_factor_base_plan) {
                const detail::ZetaBfResidueDegreeCache* const producer_cache =
                        detail::ClassUnitTransactionAccess::
                                relation_factor_base_plan_residue_degrees(
                                        candidate_class_group, order);
                if (producer_cache != nullptr &&
                    analytic_cache.seed_bf_residue_degree_cache(
                            order, *producer_cache, active_diagnostics)) {
                    validation_cache_seeded_from_relation_factor_base_plan =
                            true;
                }
            }
            bool try_proven_unit_refinement =
                    requested_proven && rank > 1;
            auto try_hnf_units = [&]() noexcept {
                SILEX_PROFILE_SCOPE(active_diagnostics,
                                    DiagnosticsModule::unit_group,
                                    "unit_group.hnf_units");
                detail::HnfFinishWorkspace* const finish_workspace =
                        native_hnf_unit_candidates
                        ? detail::class_group_finish_workspace(
                                  candidate_class_group)
                        : nullptr;
                if (finish_workspace == nullptr) {
                    return false;
                }
                slong work_precision = max_slong(precision, 160);
                const slong finish_precision =
                        candidate_class_group.analytic_finish_precision();
                if (finish_precision > 0) {
                    work_precision =
                            max_slong(work_precision, finish_precision);
                }
                auto build_units_with_product =
                        [&](OrderUnitGroup& out,
                            flint::ArbConstRef analytic_hR) noexcept {
                            out.set_diagnostics(active_diagnostics);
                            return out.is_defined() &&
                                   embeddings.refine(work_precision) &&
                                   set_hnf_units(
                                           out, order, candidate_class_group,
                                           embeddings, analytic_hR,
                                           finish_workspace,
                                           work_precision);
                        };

                flint::Arb finish_hR;
                slong finish_product_precision = 0;
                if (candidate_class_group.analytic_finish_product(
                            flint::ArbRef(finish_hR),
                            finish_product_precision) &&
                    finish_product_precision >= work_precision) {
                    OrderUnitGroup route_product_units(order);
                    flint::Arb candidate_hR;
                    if (build_units_with_product(
                                route_product_units,
                                flint::ArbConstRef(finish_hR)) &&
                        route_product_units.class_regulator_product(
                                flint::ArbRef(candidate_hR),
                                candidate_class_group,
                                work_precision) &&
                        detail::
                                class_regulator_index_is_one_from_candidate_product(
                                        flint::ArbConstRef(candidate_hR),
                                        flint::ArbConstRef(finish_hR),
                                        work_precision,
                                        active_diagnostics)) {
                        candidate_units.swap(route_product_units);
                        return true;
                    }
                    SILEX_PROFILE_EVENT(
                            active_diagnostics,
                            DiagnosticsModule::unit_group,
                            "unit_group.hnf_units.route_product_index_gt_one");
                }

                if (!analytic_cache.ensure(
                            order, work_precision, active_diagnostics,
                            candidate_class_group.factor_base(),
                            cache_bf_residue_degrees)) {
                    return false;
                }
                OrderUnitGroup analytic_units(order);
                if (!build_units_with_product(analytic_units,
                                              analytic_cache.value())) {
                    return false;
                }
                candidate_units.swap(analytic_units);
                return true;
            };
            auto try_candidate_units = [&]() noexcept {
                used_candidate_units = true;
                if (!native_do_units) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "native class/unit loop skipped "
                              "unit extraction");
                    candidate_units_ready = candidate_units.is_set() &&
                            same_order_parent(candidate_units.parent(),
                                              &order) &&
                            candidate_units.free_rank() == rank;
                    return true;
                }

                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.class_relation_units");
                if (!set_class_relation_units(
                            candidate_units_ready, unit_improved,
                            candidate_units, order, candidate_class_group,
                            embeddings, analytic_cache, relation_unit_state,
                            precision, rank, unit_add,
                            refresh_options.validation_bf_max_cutoff)) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "native unit extraction unavailable");
                    return false;
                }
                if (!candidate_units_ready) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "native unit extraction did not "
                              "yet produce a full-rank unit group");
                }
                if (unit_add == WORD_MAX) {
                    return false;
                }
                ++unit_add;
                if (candidate_units_ready) {
                    candidate_units_from_current_class_context = true;
                }
                return true;
            };
            if (rank == 0) {
                SILEX_PROFILE_EVENT(active_diagnostics,
                                    DiagnosticsModule::unit_group,
                                    "unit_group.rank_zero_units");
                if (!candidate_units.compute(order)) {
                    return false;
                }
                candidate_units_ready = true;
            } else if (order.degree() == 2 && candidate_units.compute(order)) {
                candidate_units_ready = true;
                /* Prefer the proven quadratic unit path, matching the C driver. */
            } else if (hnf_unit_candidates && requested_proven) {
                if (try_hnf_units()) {
                    used_hnf_candidate_units = true;
                    candidate_units_ready = true;
                    candidate_units_from_current_class_context = true;
                } else {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "HNF unit extraction unavailable");
                    return false;
                }
            } else if (native_unit_candidates) {
                if (!try_candidate_units()) {
                    return false;
                }
            } else if (try_proven_unit_refinement) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.relation_kernel_units");
                candidate_units_ready = candidate_units.set_relation_kernel_units(
                        order, candidate_class_group, embeddings, precision);
            } else if (options.requested_certification ==
                       CertificationMode::proven) {
                SILEX_PROFILE_EVENT(active_diagnostics,
                                    DiagnosticsModule::unit_group,
                                    "unit_group.recompute_units");
                const bool defer_initial_unit_proof =
                        options.zeta_bf_max_cutoff != 0 &&
                        !class_group_order_is_one(candidate_class_group);
                // reference validates the tentative class/unit pair before final
                // proof.  For nontrivial quotients with an analytic/BF
                // validation route available, keep the initial
                // relation-kernel unit extraction unproved; validation-loop
                // refinements still request proof when the candidate survives
                // to that stage.  Class-order-one rows keep the immediate
                // proof so the trivial-quotient route remains non-analytic.
                if (!recompute_units_from_class_context(
                            candidate_units, order, candidate_class_group,
                            embeddings, analytic_cache,
                            !defer_initial_unit_proof,
                            detail::ValidationRecomputeCause::
                                    post_compute_initial_candidate,
                            precision)) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "relation-kernel unit refinement unavailable");
                } else {
                    candidate_units_ready = true;
                    candidate_units_from_current_class_context = true;
                }
            } else if (rank > 1 ||
                       !class_group_order_is_one(candidate_class_group)) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.relation_kernel_units");
                candidate_units_ready = candidate_units.set_relation_kernel_units(
                        order, candidate_class_group, embeddings, precision);
            }
            OrderUnitGroup scratch_units(order);
            scratch_units.set_diagnostics(active_diagnostics);
            slong revalidate_passes_done = 0;
            if (!scratch_units.is_defined()) {
                return false;
            }
            auto ensure_factor_base_generation = [&]() noexcept {
                flint::Fmpz required_bound;
                if (candidate_class_group.factor_base_generation_bound(
                            flint::FmpzRef(required_bound)) &&
                    candidate_class_group.check_factor_base_generation_bound(
                            flint::FmpzConstRef(required_bound))) {
                    return true;
                }
                return false;
            };
            auto validate_candidate = [&]() noexcept {
                if (!ensure_factor_base_generation()) {
                    return false;
                }
                rank_complete = candidate_units.free_rank() == rank;
                if (used_candidate_units && requested_proven) {
                    bool validated = false;
                    if (!detail::unknown_candidate_validated(
                                validated, relation_unit_state,
                                validate_summary, candidate_units,
                                candidate_class_group, order, analytic_cache,
                                rank, precision)) {
                        return false;
                    }
                    return validated;
                }
                if (try_validate_refine_loop(
                            candidate_class_group, candidate_units,
                            scratch_units, order, options, embeddings,
                            analytic_cache, validate_summary, precision,
                            candidate_units_from_current_class_context,
                            refresh_options)) {
                    return true;
                }
                return try_validate_refine_extra_pass(
                        revalidate_passes_done, candidate_class_group,
                        candidate_units, scratch_units, order, options,
                        embeddings, analytic_cache, validate_summary, precision,
                        candidate_units_from_current_class_context,
                        refresh_options);
            };
            auto try_hnf_analytic_certification = [&]() noexcept {
                if (!used_hnf_candidate_units ||
                    !native_hnf_unit_candidates || !requested_proven) {
                    return false;
                }
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.hnf_units.analytic_certify");
                if (!ensure_factor_base_generation() ||
                    !analytic_cache.ensure(
                            order, precision, active_diagnostics,
                            candidate_class_group.factor_base(),
                            order.degree() > 2 &&
                                    options.zeta_bf_max_cutoff != 0)) {
                    return false;
                }
                return candidate_class_group
                        .try_certify_class_unit_with_units(
                                candidate_units, analytic_cache.value(),
                                precision);
            };

            const bool try_quadratic_rank_one_unit_proof =
                    !requested_proven && order.degree() == 2 && rank == 1 &&
                    class_group_order_is_one(candidate_class_group);
            const bool candidate_validated = candidate_units_ready &&
                    (try_hnf_analytic_certification() ||
                     validate_candidate());

            if (candidate_validated) {
                certified = true;
            } else if (!candidate_units_ready &&
                       !used_candidate_units) {
                if (!recompute_units_from_class_context(
                            candidate_units, order, candidate_class_group,
                            embeddings, analytic_cache,
                            requested_proven ||
                                    try_quadratic_rank_one_unit_proof,
                            detail::ValidationRecomputeCause::
                                    missing_candidate_units,
                            precision, refresh_options)) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group, LogLevel::detail,
                              "relation-kernel unit refinement unavailable");
                } else {
                    candidate_units_from_current_class_context = true;
                    if (validate_candidate()) {
                        certified = true;
                    }
                }
            } else if (!certified && try_proven_unit_refinement &&
                       !native_unit_candidates) {
                if (!recompute_units_from_class_context(
                            candidate_units, order, candidate_class_group,
                            embeddings, analytic_cache, true,
                            detail::ValidationRecomputeCause::
                                    proven_refinement_retry,
                            precision, refresh_options)) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group, LogLevel::detail,
                              "relation-kernel unit proof refinement unavailable");
                } else {
                    candidate_units_from_current_class_context = true;
                    if (validate_candidate()) {
                        certified = true;
                    }
                }
            }
            if (!certified && candidate_units_ready && rank_complete &&
                used_candidate_units && requested_proven) {
                if (unit_improved == 1 &&
                    !relation_unit_state.saturation_at_two_done) {
                    bool saturation_at_two_validated = false;
                    if (!detail::unknown_saturation_at_two(
                                saturation_at_two_validated, relation_unit_state,
                                validate_summary, candidate_class_group,
                                candidate_units, embeddings, analytic_cache,
                                order, rank, precision, active_diagnostics)) {
                        return false;
                    }
                    if (saturation_at_two_validated) {
                        certified = true;
                    }
                }
            }
            if (!certified && candidate_units_ready && rank_complete &&
                used_candidate_units && requested_proven) {
                bool saturation_validated = false;
                bool saturation_attempted = false;
                if (!detail::unknown_small_index_saturation(
                            saturation_validated, saturation_attempted,
                            relation_unit_state, validate_summary,
                            candidate_class_group, candidate_units, embeddings,
                            analytic_cache, order, rank, precision,
                            active_diagnostics)) {
                    return false;
                }
                if (saturation_validated) {
                    certified = true;
                } else if (saturation_attempted) {
                    return false;
                }
            }
            if (certified && used_candidate_units &&
                requested_proven) {
                validate_summary = ValidateRefineSummary{};
                revalidate_passes_done = 0;
                if (!try_validate_refine_loop(
                            candidate_class_group, candidate_units,
                            scratch_units, order, options, embeddings,
                            analytic_cache, validate_summary, precision,
                            false, refresh_options) &&
                    !try_validate_refine_extra_pass(
                            revalidate_passes_done, candidate_class_group,
                            candidate_units, scratch_units, order, options,
                            embeddings, analytic_cache, validate_summary,
                            precision, false, refresh_options)) {
                    certified = false;
                }
            }
            if (certified) {
                if (!detail::complete_requested_proven_relation_saturation(
                            candidate_class_group, candidate_units, order,
                            options, analytic_cache, precision)) {
                    certified = false;
                }
            }

            if (certified) {
                class_group.swap(candidate_class_group);
                class_group.set_diagnostics(active_diagnostics);
                swap(candidate_units);
                set_diagnostics(active_diagnostics);
                SILEX_VERBOSE(active_diagnostics, DiagnosticsModule::unit_group,
                              VerboseLevel::progress,
                              "class/unit computation finished");
                return true;
            }

            if (!rank_complete && rank > 1) {
                independent_unit_count = candidate_units.is_set()
                        ? candidate_units.free_rank()
                        : -1;
                if (independent_unit_count < 0 &&
                    used_candidate_units) {
                    independent_unit_count = 0;
                } else if (independent_unit_count < 0 &&
                    !relation_kernel_independent_unit_count(
                            independent_unit_count, order,
                            candidate_class_group, embeddings, precision)) {
                    independent_unit_count = -1;
                }
            }
            return false;
        };

        for (;;) {
            if (try_current_candidate()) {
                return true;
            }
            SILEX_PROFILE_EVENT(
                    active_diagnostics, DiagnosticsModule::unit_group,
                    detail::validate_refine_outcome_profile_label(
                            validate_summary.outcome));

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
            if (log_enabled(active_diagnostics, DiagnosticsModule::unit_group,
                            LogLevel::detail)) {
                std::string message = "validation summary target=";
                message += std::to_string(static_cast<long long>(target));
                message += " outcome=";
                message += validate_refine_outcome_name(
                        validate_summary.outcome);
                message += " rank_complete=";
                message += rank_complete ? "1" : "0";
                if (flint::fmpz_fits_si(
                            flint::FmpzConstRef(validate_summary
                                                        .last_index_bound))) {
                    message += " last_index_bound=";
                    message += std::to_string(static_cast<long long>(
                            flint::fmpz_get_si(flint::FmpzConstRef(
                                    validate_summary.last_index_bound))));
                }
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail, message.c_str());
            }
#endif

            const bool requested_proven =
                    options.requested_certification ==
                    CertificationMode::proven;
            if (compute_ok && native_unit_candidates &&
                requested_proven && rank_complete &&
                validate_summary.outcome ==
                        detail::ValidateRefineOutcome::no_progress &&
                flint::fmpz_cmp_ui(
                        flint::FmpzConstRef(
                                validate_summary.last_index_bound),
                        1) > 0) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.lll_relation_slice");
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "extending current class-group candidate with "
                          "LLL relation slice");

                const slong relation_count_before =
                        candidate_class_group.relation_count();
                flint::Fmpz class_order_before;
                const bool have_class_order_before =
                        candidate_class_group.order(
                                flint::FmpzRef(class_order_before));
                bool slice_added_without_retry = false;
                if (policy->continuation !=
                    detail::NativeContinuationStrategy::lll) {
                    return false;
                }
                const bool relation_slice_extended =
                        candidate_class_group
                                .extend_lll_relation_slice_(
                                        order, factor_base_bound,
                                        local_options);
                if (relation_slice_extended &&
                    candidate_class_group.relation_count() >
                            relation_count_before) {
                    bool retry_current_candidate = true;
                    if (have_class_order_before) {
                        flint::Fmpz class_order_after;
                        if (!candidate_class_group.order(
                                    flint::FmpzRef(class_order_after))) {
                            return false;
                        }
                        const bool class_order_unchanged =
                                flint::fmpz_equal(
                                        flint::FmpzConstRef(
                                                class_order_before),
                                        flint::FmpzConstRef(
                                                class_order_after));
                        if (class_order_unchanged) {
                            native_do_units =
                                    candidate_class_group
                                            .relation_kernel_unit_count() >=
                                    rank;
                            retry_current_candidate = false;
                        } else {
                            if (unit_add > WORD_MAX - 2) {
                                return false;
                            }
                            unit_add += 2;
                            native_do_units = false;
                        }
                    } else {
                        native_do_units =
                                candidate_class_group
                                        .relation_kernel_unit_count() >= rank;
                    }
                    if (!retry_current_candidate) {
                        SILEX_LOG(active_diagnostics,
                                  DiagnosticsModule::unit_group,
                                  LogLevel::detail,
                                  "LLL relation slice preserved the "
                                  "class-group order; trying validation "
                                  "target bump");
                        slice_added_without_retry = true;
                    } else {
                        SILEX_PROFILE_EVENT(
                                active_diagnostics,
                                DiagnosticsModule::unit_group,
                                "unit_group.candidate_retry."
                                "extended_current_candidate_with_lll");
                        continue;
                    }
                }
                if (!slice_added_without_retry) {
                    SILEX_LOG(active_diagnostics,
                              DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "LLL relation slice did not progress");
                }
            }

            slong missing_unit_bump = 0;
            if (compute_ok &&
                missing_unit_target_bump(
                        missing_unit_bump, rank, independent_unit_count, target,
                        options.max_relations)) {
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "raising relation-kernel target for missing unit rank");
            }

            slong analytic_bump = 0;
            if (compute_ok &&
                analytic_bumps_done < kComputeAnalyticBumpMaxPasses &&
                validation_target_bump(
                        analytic_bump, validate_summary, rank, target,
                        options.max_relations,
                        options.requested_certification ==
                                CertificationMode::proven,
                        rank_complete)) {
                post_finite_refinement_phase_budget = max_slong(
                        post_finite_refinement_phase_budget,
                        kComputePostFinitePhaseBudget);
                ++analytic_bumps_done;
            }

            if (missing_unit_bump > 0) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.missing_unit_target_bump");
                target += missing_unit_bump;
                local_options.target_relation_kernel_units = target;
                local_options.post_finite_refinement_phase_budget =
                        post_finite_refinement_phase_budget;

                const slong relation_count_before =
                        candidate_class_group.relation_count();
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "extending current class-group candidate after "
                          "validation target bump");
                bool extended_candidate = false;
                if (candidate_class_group.extend_relation_kernel_units_(
                            order, factor_base_bound, local_options) &&
                    candidate_class_group.relation_count() >
                            relation_count_before) {
                    extended_candidate = true;
                }
                if (extended_candidate) {
                    SILEX_PROFILE_EVENT(
                            active_diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.candidate_retry.extended_current_candidate");
                    continue;
                }

                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.restart_after_missing_unit_target_bump");
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "current class-group candidate extension did not "
                          "make progress; restarting candidate");
            } else if (analytic_bump > 0) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.analytic_target_bump");
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "raising relation-kernel target after validation");
                target += analytic_bump;
                local_options.target_relation_kernel_units = target;
                local_options.post_finite_refinement_phase_budget =
                        post_finite_refinement_phase_budget;

                const slong relation_count_before =
                        candidate_class_group.relation_count();
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "extending current class-group candidate after "
                          "analytic validation target bump");
                bool extended_candidate = false;
                if (candidate_class_group.extend_relation_kernel_units_(
                            order, factor_base_bound, local_options) &&
                    candidate_class_group.relation_count() >
                            relation_count_before) {
                    extended_candidate = true;
                }
                if (extended_candidate) {
                    SILEX_PROFILE_EVENT(
                            active_diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.candidate_retry.extended_current_candidate_after_analytic_bump");
                    continue;
                }

                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.restart_after_analytic_bump");
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "current class-group candidate extension did not "
                          "make progress; restarting candidate");
            } else if (target == lower_target) {
                if (rank == 1 && !rank_one_surplus_started &&
                    rank_one_surplus_target > lower_target) {
                    SILEX_PROFILE_EVENT(
                            active_diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.candidate_retry.restart_after_rank_one_surplus");
                    SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "raising relation-kernel target");
                    target = rank_one_surplus_target;
                    rank_one_surplus_started = true;
                } else {
                    SILEX_PROFILE_EVENT(
                            active_diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.candidate_retry.exhausted_lower_target");
                    SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                              LogLevel::detail,
                              "class/unit computation did not certify");
                    return false;
                }
            } else if (rank == 1 && rank_one_surplus_started &&
                       target == lower_target + 1) {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.exhausted_rank_one_fallback");
                SILEX_LOG(active_diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "class/unit target fallback exhausted");
                return false;
            } else {
                SILEX_PROFILE_EVENT(
                        active_diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.candidate_retry.restart_after_target_decrement");
                --target;
            }
            break;
        }
    }
}

}  // namespace silex
