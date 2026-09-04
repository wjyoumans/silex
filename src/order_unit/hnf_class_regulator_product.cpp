#include "order_unit_internal.hpp"

#include <vector>

#include <silex/flint/arf.hpp>
#include <silex/signature.hpp>
#include <silex/unit.hpp>

namespace silex::detail {

namespace {

bool arb_abs_upper_lt_2exp(const arb_struct* value,
                           slong exponent,
                           slong precision) noexcept {
    if (precision <= 0 || ::arb_is_finite(value) == 0) {
        return false;
    }

    flint::Arf upper;
    ::arb_get_abs_ubound_arf(upper.raw(), value, precision);
    return flint::arf_is_finite(upper) &&
           ::arf_cmpabs_2exp_si(upper.raw(), exponent) < 0;
}

}  // namespace

bool hnf_regulator_multiple_from_basis(
        flint::ArbRef out,
        const flint::ArbMat& basis,
        slong degree,
        slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(basis);
    const slong cols = flint::arb_mat_ncols_value(basis);
    if (precision <= 0 || degree <= 0 || rows < 2 || rows != cols) {
        return false;
    }

    flint::Arb determinant;
    flint::Arb regulator_multiple;
    arb_mat_det(determinant.raw(), basis.raw(), precision);
    arb_abs(determinant.raw(), determinant.raw());
    flint::arb_div_ui(regulator_multiple, determinant,
                      static_cast<ulong>(degree), precision);

    flint::Arb minimum;
    flint::arb_one(minimum);
    ::arb_mul_2exp_si(minimum.raw(), minimum.raw(), -3);
    if (!flint::arb_is_finite(regulator_multiple) ||
        !flint::arb_is_positive(regulator_multiple) ||
        ::arb_ge(regulator_multiple.raw(), minimum.raw()) == 0) {
        return false;
    }

    const slong rank = rows - 1;
    flint::ArbMat lower_right(rank, rank);
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < rank; ++j) {
            flint::arb_set(
                    flint::arb_mat_entry_ref(lower_right, i, j),
                    flint::arb_mat_entry_ref(
                            flint::ArbMatConstRef(basis), i + 1, j + 1));
        }
    }

    flint::Arb minor;
    arb_mat_det(minor.raw(), lower_right.raw(), precision);
    arb_abs(minor.raw(), minor.raw());
    if (!flint::arb_is_finite(minor) || !flint::arb_is_positive(minor)) {
        return false;
    }

    flint::Arb difference;
    flint::arb_sub(difference, minor, regulator_multiple, precision);
    if (!flint::arb_is_finite(difference)) {
        return false;
    }
    if (!flint::arb_is_zero(difference)) {
        flint::Arf difference_upper;
        flint::Arf minor_lower;
        ::arb_get_abs_ubound_arf(
                difference_upper.raw(), difference.raw(), precision);
        ::arb_get_abs_lbound_arf(minor_lower.raw(), minor.raw(), precision);
        if (!flint::arf_is_finite(difference_upper) ||
            ::arf_is_zero(difference_upper.raw()) != 0 ||
            !flint::arf_is_finite(minor_lower) ||
            ::arf_is_zero(minor_lower.raw()) != 0) {
            return false;
        }

        flint::Fmpz difference_exponent;
        flint::Fmpz minor_exponent;
        flint::Fmpz relative_exponent;
        ::arf_abs_bound_lt_2exp_fmpz(
                difference_exponent.raw(), difference_upper.raw());
        ::fmpz_sub_ui(
                difference_exponent.raw(), difference_exponent.raw(), 1);
        ::arf_abs_bound_lt_2exp_fmpz(
                minor_exponent.raw(), minor_lower.raw());
        ::fmpz_sub_ui(minor_exponent.raw(), minor_exponent.raw(), 1);
        ::fmpz_sub(relative_exponent.raw(), difference_exponent.raw(),
                   minor_exponent.raw());
        if (::fmpz_cmp_si(relative_exponent.raw(), -20) > 0) {
            return false;
        }
    }

    flint::arb_set(out, flint::ArbConstRef(regulator_multiple));
    return true;
}

namespace {

bool inverse_identity_residual_bits_ok(const flint::ArbMat& inverse,
                                            const flint::ArbMat& basis,
                                            slong precision) noexcept {
    const slong rows = flint::arb_mat_nrows_value(inverse);
    const slong cols = flint::arb_mat_ncols_value(inverse);
    if (precision <= 0 || rows <= 0 || rows != cols ||
        flint::arb_mat_nrows_value(basis) != rows ||
        flint::arb_mat_ncols_value(basis) != cols) {
        return false;
    }

    flint::ArbMat product(rows, cols);
    arb_mat_mul(product.raw(), inverse.raw(), basis.raw(), precision);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            flint::Arb residual;
            arb_set(residual.raw(), arb_mat_entry(product.raw(), i, j));
            if (i == j) {
                arb_sub_ui(residual.raw(), residual.raw(), 1, precision);
            }
            if (!arb_abs_upper_lt_2exp(residual.raw(), -16, precision)) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

RegulatorPivotOutcome
hnf_class_regulator_product_and_unit_matrix(
        flint::ArbRef out,
        slong* independent_count_out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        flint::FmpzMat* unit_matrix_out,
        flint::Arb* regulator_out,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept {
    if (independent_count_out != nullptr) {
        *independent_count_out = 0;
    }
    const DiagnosticsContext* const diagnostics = class_group.diagnostics();
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision) ||
        !flint::arb_is_finite(analytic_hR) ||
        !flint::arb_is_positive(analytic_hR)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: invalid inputs");
        return RegulatorPivotOutcome::invalid;
    }

    slong rank = -1;
    slong places = 0;
    Signature sig;
    if (!unit_rank(rank, *order.parent()) || rank <= 0 ||
        !compact_places(places, embeddings) || places != rank + 1 ||
        !signature(sig, *order.parent())) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: rank/signature setup");
        return RegulatorPivotOutcome::invalid;
    }

    flint::ArbMat unit_logs(
            class_group.relation_count() - class_group.relation_rank(),
            places);
    std::vector<slong> selected;
    if (!hnf_unit_log_matrix(unit_logs, class_group, embeddings,
                                  workspace, precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: HNF unit log matrix");
        // This legacy bool conflates structural, exact, and numeric failures.
        // The five-path pivot correction cannot safely relabel that opaque
        // upstream result as certified precision ambiguity: only the Arb
        // selector below owns the adjudicated tri-state distinction.
        return RegulatorPivotOutcome::invalid;
    }
    const RegulatorPivotOutcome pivot_outcome =
            regulator_pivot_unit_indices(
                    selected, unit_logs, sig, rank, places, precision);
    if (pivot_outcome != RegulatorPivotOutcome::success) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: pivot elimination");
        return pivot_outcome;
    }
    const slong independent_count = static_cast<slong>(selected.size());
    if (static_cast<slong>(selected.size()) != rank) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: unit rank defect");
        if (independent_count_out != nullptr) {
            *independent_count_out = independent_count;
        }
        return RegulatorPivotOutcome::success;
    }

    flint::ArbMat basis(places, places);
    for (slong i = 0; i < sig.r1(); ++i) {
        flint::arb_set_si(flint::arb_mat_entry_ref(basis, i, 0), 1);
    }
    for (slong i = sig.r1(); i < places; ++i) {
        flint::arb_set_si(flint::arb_mat_entry_ref(basis, i, 0), 2);
    }
    for (slong col = 0; col < rank; ++col) {
        for (slong row = 0; row < places; ++row) {
            arb_set(arb_mat_entry(basis.raw(), row, col + 1),
                    arb_mat_entry(unit_logs.raw(), selected[col], row));
        }
    }

    flint::Arb regulator_multiple;
    if (!hnf_regulator_multiple_from_basis(
                flint::ArbRef(regulator_multiple), basis, order.degree(),
                precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: regulator multiple");
        return RegulatorPivotOutcome::precision_inconclusive;
    }

    flint::ArbMat inverse(places, places);
    if (arb_mat_inv(inverse.raw(), basis.raw(), precision) == 0) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: inverse");
        return RegulatorPivotOutcome::precision_inconclusive;
    }
    if (!inverse_identity_residual_bits_ok(inverse, basis, precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "certify_regulator_multiple failed: inverse residual bits");
        return RegulatorPivotOutcome::precision_inconclusive;
    }

    const slong unit_columns = flint::arb_mat_nrows_value(unit_logs);
    flint::ArbMat coordinates(rank, unit_columns);
    for (slong j = 0; j < unit_columns; ++j) {
        for (slong coord = 0; coord < rank; ++coord) {
            flint::arb_zero(flint::arb_mat_entry_ref(coordinates, coord, j));
            for (slong i = 0; i < places; ++i) {
                flint::Arb product;
                flint::arb_mul(
                        product,
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(inverse), coord + 1, i),
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(unit_logs), j, i),
                        precision);
                arb_add(arb_mat_entry(coordinates.raw(), coord, j),
                        arb_mat_entry(coordinates.raw(), coord, j),
                        product.raw(), precision);
            }
        }
    }

    flint::Fmpz class_order;
    flint::Arb class_order_arb;
    flint::Arb z;
    flint::Arb regulator;
    flint::Arb product;
    if (!class_group.order(flint::FmpzRef(class_order))) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: class order unavailable");
        return RegulatorPivotOutcome::invalid;
    }
    flint::arb_set_fmpz(class_order_arb, flint::FmpzConstRef(class_order));
    flint::arb_div(z, class_order_arb, analytic_hR, precision);
    if (!reduced_regulator_from_coordinates(
                flint::ArbRef(regulator),
                coordinates, flint::ArbConstRef(regulator_multiple),
                flint::ArbConstRef(z), unit_matrix_out, precision,
                diagnostics)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: reduced regulator reconstruction");
        return RegulatorPivotOutcome::precision_inconclusive;
    }
    flint::arb_mul_fmpz(product, regulator,
                        flint::FmpzConstRef(class_order), precision);
    flint::arb_set(out, flint::ArbConstRef(product));
    if (regulator_out != nullptr) {
        flint::arb_set(*regulator_out, flint::ArbConstRef(regulator));
    }
    if (independent_count_out != nullptr) {
        *independent_count_out = independent_count;
    }
    return RegulatorPivotOutcome::success;
}

RegulatorPivotOutcome hnf_class_regulator_product(
        flint::ArbRef out,
        slong& independent_count,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept {
    independent_count = 0;
    slong rank = -1;
    if (order.parent() == nullptr ||
        !unit_rank(rank, *order.parent())) {
        return RegulatorPivotOutcome::invalid;
    }
    if (rank == 0) {
        // The source regulator-multiple certification and reconstruction both return
        // regulator 1 when RU = R1 + R2 = 1.  The candidate hR is therefore
        // the exact class order, with no HNF unit matrix to reconstruct.
        flint::Fmpz class_order;
        if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                             precision) ||
            !flint::arb_is_finite(analytic_hR) ||
            !flint::arb_is_positive(analytic_hR) ||
            !class_group.order(flint::FmpzRef(class_order))) {
            return RegulatorPivotOutcome::invalid;
        }
        flint::Arb candidate_hR;
        flint::arb_set_fmpz(candidate_hR,
                            flint::FmpzConstRef(class_order));
        flint::arb_set(out, flint::ArbConstRef(candidate_hR));
        return RegulatorPivotOutcome::success;
    }
    return hnf_class_regulator_product_and_unit_matrix(
            out, &independent_count, order, class_group, embeddings,
            analytic_hR, nullptr, nullptr, workspace, precision);
}

}  // namespace silex::detail
