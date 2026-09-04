#include "order_unit_internal.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include <silex/flint/arf.hpp>
#include <silex/signature.hpp>

namespace silex::detail {

struct RegulatorPivotCandidate {
    flint::Fmpz exponent;
    flint::Fmpz upper_exponent;
    bool exact_zero = false;
    bool exponent_certified = false;
};

static RegulatorPivotOutcome pivot_candidate(
        RegulatorPivotCandidate& out,
        const arb_struct* value,
        slong precision) noexcept {
    out.exact_zero = false;
    out.exponent_certified = false;
    if (precision <= 0 || ::arb_is_finite(value) == 0) {
        return RegulatorPivotOutcome::invalid;
    }
    if (::arb_is_zero(value) != 0) {
        out.exact_zero = true;
        return RegulatorPivotOutcome::success;
    }

    flint::Arf upper;
    ::arb_get_abs_ubound_arf(upper.raw(), value, precision);
    if (!flint::arf_is_finite(upper) || ::arf_is_zero(upper.raw()) != 0) {
        return RegulatorPivotOutcome::invalid;
    }
    ::arf_abs_bound_lt_2exp_fmpz(out.upper_exponent.raw(), upper.raw());
    ::fmpz_sub_ui(out.upper_exponent.raw(), out.upper_exponent.raw(), 1);

    // Exact zero is distinct from a nonpoint ball containing zero.  The
    // latter has no certified nonzero lower exponent, but its upper exponent
    // can still prove that it cannot affect reference's ordered winner.
    if (::arb_contains_zero(value) != 0) {
        return RegulatorPivotOutcome::success;
    }

    flint::Arf lower;
    ::arb_get_abs_lbound_arf(lower.raw(), value, precision);
    if (!flint::arf_is_finite(lower) || ::arf_is_zero(lower.raw()) != 0) {
        return RegulatorPivotOutcome::invalid;
    }
    ::arf_abs_bound_lt_2exp_fmpz(out.exponent.raw(), lower.raw());
    ::fmpz_sub_ui(out.exponent.raw(), out.exponent.raw(), 1);
    out.exponent_certified =
            ::fmpz_equal(out.exponent.raw(), out.upper_exponent.raw()) != 0;
    return RegulatorPivotOutcome::success;
}

static RegulatorPivotOutcome regulator_pivot_rows(
        std::vector<slong>& out,
        flint::ArbMat& working,
        slong precision) noexcept {
    out.clear();
    const slong rows = flint::arb_mat_nrows_value(working);
    const slong columns = flint::arb_mat_ncols_value(working);
    if (rows <= 0 || columns < 0 || precision <= 0) {
        return RegulatorPivotOutcome::invalid;
    }
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < columns; ++col) {
            if (::arb_is_finite(
                        arb_mat_entry(working.raw(), row, col)) == 0) {
                return RegulatorPivotOutcome::invalid;
            }
        }
    }

    std::vector<slong> pivot_rows(static_cast<std::size_t>(columns), -1);
    std::vector<slong> pivot_column_for_row(
            static_cast<std::size_t>(rows), -1);
    for (slong col = 0; col < columns; ++col) {
        std::vector<RegulatorPivotCandidate> candidates(
                static_cast<std::size_t>(rows));
        slong pivot_row = -1;
        for (slong row = 0; row < rows; ++row) {
            if (pivot_column_for_row[static_cast<std::size_t>(row)] >= 0) {
                continue;
            }

            RegulatorPivotCandidate& candidate =
                    candidates[static_cast<std::size_t>(row)];
            const RegulatorPivotOutcome candidate_outcome =
                    pivot_candidate(
                            candidate,
                            arb_mat_entry(working.raw(), row, col),
                            precision);
            if (candidate_outcome != RegulatorPivotOutcome::success) {
                return candidate_outcome;
            }
            if (!candidate.exponent_certified ||
                ::fmpz_cmp_si(candidate.exponent.raw(), -32) <= 0) {
                continue;
            }
            if (pivot_row < 0 ||
                ::fmpz_cmp(
                        candidate.exponent.raw(),
                        candidates[static_cast<std::size_t>(pivot_row)]
                                .exponent.raw()) > 0) {
                pivot_row = row;
            }
        }

        // reference scans unused rows in order and replaces only on a strictly
        // greater gexpo.  An uncertain earlier row can therefore change the
        // winner on equality; a later row can change it only by exceeding it.
        if (pivot_row < 0) {
            for (slong row = 0; row < rows; ++row) {
                if (pivot_column_for_row[static_cast<std::size_t>(row)] >= 0) {
                    continue;
                }
                const RegulatorPivotCandidate& candidate =
                        candidates[static_cast<std::size_t>(row)];
                if (!candidate.exact_zero &&
                    !candidate.exponent_certified &&
                    ::fmpz_cmp_si(candidate.upper_exponent.raw(), -32) > 0) {
                    return RegulatorPivotOutcome::precision_inconclusive;
                }
            }
            continue;
        }

        const flint::Fmpz& best_exponent =
                candidates[static_cast<std::size_t>(pivot_row)].exponent;
        for (slong row = 0; row < rows; ++row) {
            if (pivot_column_for_row[static_cast<std::size_t>(row)] >= 0) {
                continue;
            }
            const RegulatorPivotCandidate& candidate =
                    candidates[static_cast<std::size_t>(row)];
            if (candidate.exact_zero || candidate.exponent_certified) {
                continue;
            }
            const int comparison = ::fmpz_cmp(
                    candidate.upper_exponent.raw(), best_exponent.raw());
            if ((row < pivot_row && comparison >= 0) ||
                (row > pivot_row && comparison > 0)) {
                return RegulatorPivotOutcome::precision_inconclusive;
            }
        }

        pivot_rows[static_cast<std::size_t>(col)] = pivot_row;
        pivot_column_for_row[static_cast<std::size_t>(pivot_row)] = col;
        const arb_struct* pivot =
                arb_mat_entry(working.raw(), pivot_row, col);

        flint::Arb negative_inverse;
        arb_inv(negative_inverse.raw(), pivot, precision);
        arb_neg(negative_inverse.raw(), negative_inverse.raw());
        for (slong later = col + 1; later < columns; ++later) {
            arb_mul(arb_mat_entry(working.raw(), pivot_row, later),
                    negative_inverse.raw(),
                    arb_mat_entry(working.raw(), pivot_row, later),
                    precision);
        }

        for (slong row = 0; row < rows; ++row) {
            if (pivot_column_for_row[static_cast<std::size_t>(row)] >= 0) {
                continue;
            }

            flint::Arb multiplier;
            arb_set(multiplier.raw(),
                    arb_mat_entry(working.raw(), row, col));
            arb_zero(arb_mat_entry(working.raw(), row, col));
            for (slong later = col + 1; later < columns; ++later) {
                flint::Arb term;
                arb_mul(term.raw(), multiplier.raw(),
                        arb_mat_entry(working.raw(), pivot_row, later),
                        precision);
                arb_add(arb_mat_entry(working.raw(), row, later),
                        arb_mat_entry(working.raw(), row, later),
                        term.raw(), precision);
            }
        }

        for (slong later = col; later < columns; ++later) {
            arb_zero(arb_mat_entry(working.raw(), pivot_row, later));
        }
    }

    out = std::move(pivot_rows);
    return RegulatorPivotOutcome::success;
}

RegulatorPivotOutcome regulator_pivot_rows_for_testing(
        std::vector<slong>& out,
        flint::ArbMatConstRef matrix,
        slong precision) noexcept {
    out.clear();
    const slong rows = flint::arb_mat_nrows_value(matrix);
    const slong columns = flint::arb_mat_ncols_value(matrix);
    if (rows <= 0 || columns < 0 || precision <= 0) {
        return RegulatorPivotOutcome::invalid;
    }

    flint::ArbMat working(rows, columns);
    flint::arb_mat_set(flint::ArbMatRef(working), matrix);
    std::vector<slong> pivot_rows;
    const RegulatorPivotOutcome outcome =
            regulator_pivot_rows(pivot_rows, working, precision);
    if (outcome == RegulatorPivotOutcome::success) {
        out = std::move(pivot_rows);
    }
    return outcome;
}

static RegulatorPivotOutcome clean_regulator_unit_columns(
        std::vector<slong>& retained,
        const flint::ArbMat& unit_logs,
        slong places,
        slong precision) noexcept {
    retained.clear();
    const slong unit_columns = flint::arb_mat_nrows_value(unit_logs);
    if (unit_columns < 0 || places <= 0 || precision <= 0 ||
        flint::arb_mat_ncols_value(unit_logs) != places) {
        return RegulatorPivotOutcome::invalid;
    }
    for (slong col = 0; col < unit_columns; ++col) {
        for (slong row = 0; row < places; ++row) {
            if (::arb_is_finite(
                        arb_mat_entry(unit_logs.raw(), col, row)) == 0) {
                return RegulatorPivotOutcome::invalid;
            }
        }
    }

    std::vector<slong> candidate;
    candidate.reserve(static_cast<std::size_t>(unit_columns));
    for (slong col = 0; col < unit_columns; ++col) {
        bool certified_retain = false;
        bool possible_retain = false;
        // Scan the complete source column.  A later certified entry resolves
        // earlier retention uncertainty, while any later nonfinite entry
        // still makes the whole input invalid.
        for (slong row = 0; row < places; ++row) {
            const arb_struct* value =
                    arb_mat_entry(unit_logs.raw(), col, row);
            if (::arb_is_zero(value) != 0) {
                continue;
            }

            flint::Arf upper;
            flint::Arf lower;
            ::arb_get_abs_ubound_arf(upper.raw(), value, precision);
            ::arb_get_abs_lbound_arf(lower.raw(), value, precision);
            if (!flint::arf_is_finite(upper) ||
                !flint::arf_is_finite(lower)) {
                return RegulatorPivotOutcome::invalid;
            }
            if (::arf_cmpabs_2exp_si(lower.raw(), -2) >= 0) {
                certified_retain = true;
            } else if (::arf_cmpabs_2exp_si(upper.raw(), -2) >= 0) {
                possible_retain = true;
            }
        }

        if (certified_retain) {
            candidate.push_back(col);
        } else if (possible_retain) {
            return RegulatorPivotOutcome::precision_inconclusive;
        }
    }

    retained = std::move(candidate);
    return RegulatorPivotOutcome::success;
}

RegulatorPivotOutcome regulator_pivot_unit_indices(
        std::vector<slong>& out,
        const flint::ArbMat& unit_logs,
        const Signature& sig,
        slong rank,
        slong places,
        slong precision,
        std::vector<slong>* pivot_rows_out) noexcept {
    out.clear();
    if (pivot_rows_out != nullptr) {
        pivot_rows_out->clear();
    }
    if (places <= 0 || rank < 0 || rank != places - 1 || precision <= 0 ||
        flint::arb_mat_ncols_value(unit_logs) != places ||
        sig.r1() < 0 || sig.r2() < 0 ||
        sig.r1() > places || sig.r2() != places - sig.r1()) {
        return RegulatorPivotOutcome::invalid;
    }
    if (rank == 0) {
        // reference certify_regulator_multiple's RU == 1 branch bypasses filter_unit_log_columns
        // and pivoting.  Do not classify irrelevant unit-log intervals here.
        return RegulatorPivotOutcome::success;
    }

    std::vector<slong> retained_unit_columns;
    const RegulatorPivotOutcome clean_outcome =
            clean_regulator_unit_columns(
                    retained_unit_columns, unit_logs, places, precision);
    if (clean_outcome != RegulatorPivotOutcome::success) {
        return clean_outcome;
    }

    // Source regulator-multiple certification first filters unit-log columns, then
    // calls RgM_pivots on shallowconcat(T, Ar).  Its callback compares integer
    // gexpo bins, preserves the first row on ties, and accepts only gexpo >
    // -32.  The retained-to-original map keeps consumers on uncleaned Ar.
    const slong columns =
            static_cast<slong>(retained_unit_columns.size()) + 1;
    flint::ArbMat working(places, columns);
    for (slong row = 0; row < places; ++row) {
        flint::arb_set_si(flint::arb_mat_entry_ref(working, row, 0),
                          row < sig.r1() ? 1 : 2);
    }
    for (slong col = 0;
         col < static_cast<slong>(retained_unit_columns.size()); ++col) {
        const slong source_col =
                retained_unit_columns[static_cast<std::size_t>(col)];
        for (slong row = 0; row < places; ++row) {
            arb_set(arb_mat_entry(working.raw(), row, col + 1),
                    arb_mat_entry(unit_logs.raw(), source_col, row));
        }
    }

    std::vector<slong> pivot_rows;
    const RegulatorPivotOutcome pivot_outcome =
            regulator_pivot_rows(pivot_rows, working, precision);
    if (pivot_outcome != RegulatorPivotOutcome::success) {
        return pivot_outcome;
    }

    std::vector<slong> selected;
    for (slong col = 1; col < columns; ++col) {
        if (pivot_rows[static_cast<std::size_t>(col)] >= 0) {
            selected.push_back(retained_unit_columns[
                    static_cast<std::size_t>(col - 1)]);
        }
    }
    out = std::move(selected);
    if (pivot_rows_out != nullptr) {
        *pivot_rows_out = std::move(pivot_rows);
    }
    return RegulatorPivotOutcome::success;
}

RegulatorPivotOutcome regulator_pivot_unit_indices_for_testing(
        std::vector<slong>& out,
        std::vector<slong>& pivot_rows,
        const flint::ArbMat& unit_logs,
        const Signature& sig,
        slong rank,
        slong places,
        slong precision) noexcept {
    return regulator_pivot_unit_indices(
            out, unit_logs, sig, rank, places, precision, &pivot_rows);
}

}  // namespace silex::detail
