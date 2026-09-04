#include "order_unit_internal.hpp"

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>

#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_vec.hpp>

namespace silex::detail {

namespace {

slong trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 && ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool fmpq_abs_error_less(flint::FmpqConstRef left,
                         flint::FmpqConstRef right,
                         flint::FmpqConstRef target) noexcept {
    flint::Fmpq left_error;
    flint::Fmpq right_error;
    fmpq_sub(left_error.raw(), left.raw(), target.raw());
    fmpq_abs(left_error.raw(), left_error.raw());
    fmpq_sub(right_error.raw(), right.raw(), target.raw());
    fmpq_abs(right_error.raw(), right_error.raw());
    return fmpq_cmp(left_error.raw(), right_error.raw()) < 0;
}

bool fmpq_set_best_denominator_bounded(
        flint::Fmpq& out,
        flint::FmpqConstRef value,
        flint::FmpzConstRef max_denominator) noexcept {
    if (flint::fmpz_sgn(max_denominator) <= 0) {
        return false;
    }

    const slong cf_bound = fmpq_cfrac_bound(value.raw());
    if (cf_bound <= 0) {
        return false;
    }

    flint::FmpzVec coefficients(cf_bound);
    flint::Fmpq remainder;
    const slong len =
            fmpq_get_cfrac(coefficients.data(), remainder.raw(),
                           value.raw(), cf_bound);
    if (len <= 0) {
        return false;
    }

    flint::Fmpz p_prev2;
    flint::Fmpz p_prev1;
    flint::Fmpz q_prev2;
    flint::Fmpz q_prev1;
    flint::fmpz_zero(flint::FmpzRef(p_prev2));
    flint::fmpz_one(flint::FmpzRef(p_prev1));
    flint::fmpz_one(flint::FmpzRef(q_prev2));
    flint::fmpz_zero(flint::FmpzRef(q_prev1));

    bool have = false;
    for (slong i = 0; i < len; ++i) {
        if (i > 0 &&
            flint::fmpz_cmp(
                    flint::FmpzConstRef(coefficients.data() + i),
                    max_denominator) > 0) {
            flint::Fmpz remaining;
            flint::Fmpz t;
            flint::fmpz_sub(flint::FmpzRef(remaining),
                            max_denominator, flint::FmpzConstRef(q_prev2));
            if (flint::fmpz_sgn(flint::FmpzConstRef(remaining)) > 0) {
                fmpz_fdiv_q(t.raw(), remaining.raw(), q_prev1.raw());
                if (flint::fmpz_sgn(flint::FmpzConstRef(t)) > 0) {
                    flint::Fmpz semi_p;
                    flint::Fmpz semi_q;
                    fmpz_mul(semi_p.raw(), t.raw(), p_prev1.raw());
                    fmpz_add(semi_p.raw(), semi_p.raw(), p_prev2.raw());
                    fmpz_mul(semi_q.raw(), t.raw(), q_prev1.raw());
                    fmpz_add(semi_q.raw(), semi_q.raw(), q_prev2.raw());

                    flint::Fmpq semi;
                    flint::fmpq_set_fmpz_frac(
                            semi, flint::FmpzConstRef(semi_p),
                            flint::FmpzConstRef(semi_q));
                    if (!have ||
                        fmpq_abs_error_less(flint::FmpqConstRef(semi),
                                            flint::FmpqConstRef(out),
                                            value)) {
                        flint::fmpq_set(flint::FmpqRef(out),
                                        flint::FmpqConstRef(semi));
                        have = true;
                    }
                }
            }
            break;
        }

        flint::Fmpz p;
        flint::Fmpz q;
        fmpz_mul(p.raw(), coefficients.data() + i, p_prev1.raw());
        fmpz_add(p.raw(), p.raw(), p_prev2.raw());
        fmpz_mul(q.raw(), coefficients.data() + i, q_prev1.raw());
        fmpz_add(q.raw(), q.raw(), q_prev2.raw());

        if (flint::fmpz_cmp(flint::FmpzConstRef(q),
                            max_denominator) <= 0) {
            flint::fmpq_set_fmpz_frac(
                    out, flint::FmpzConstRef(p), flint::FmpzConstRef(q));
            have = true;
            p_prev2.swap(p_prev1);
            p_prev1.swap(p);
            q_prev2.swap(q_prev1);
            q_prev1.swap(q);
            continue;
        }
        break;
    }

    return have;
}

bool arb_best_denominator_bounded_rational(
        flint::Fmpq& out,
        flint::ArbConstRef value,
        flint::FmpzConstRef max_denominator) noexcept {
    if (flint::fmpz_sgn(max_denominator) <= 0 ||
        !flint::arb_is_finite(value)) {
        return false;
    }

    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value.raw()));
    if (!flint::arf_is_finite(midpoint)) {
        return false;
    }

    flint::Fmpq midpoint_q;
    arf_get_fmpq(midpoint_q.raw(), midpoint.raw());
    if (!fmpq_set_best_denominator_bounded(
                out, flint::FmpqConstRef(midpoint_q), max_denominator) ||
        !flint::arb_contains_fmpq(value, out)) {
        return false;
    }

    return true;
}

bool arb_floor_positive(flint::Fmpz& out,
                        flint::ArbConstRef value,
                        slong precision) noexcept {
    if (precision <= 0 || !flint::arb_is_finite(value) ||
        !flint::arb_is_positive(value)) {
        return false;
    }

    flint::Arf lower;
    ::arb_get_lbound_arf(lower.raw(), value.raw(), precision);
    if (!flint::arf_is_finite(lower)) {
        return false;
    }

    flint::arf_get_fmpz(out, lower, ARF_RND_FLOOR);
    return flint::fmpz_sgn(flint::FmpzConstRef(out)) > 0;
}

bool column_lattice_hnf_determinant(flint::Fmpz& out,
                                    const flint::FmpzMat& matrix) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    if (rows <= 0 || cols <= 0) {
        return false;
    }

    flint::FmpzMat transpose(cols, rows);
    flint::FmpzMat hnf(cols, rows);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transpose),
                              flint::FmpzMatConstRef(matrix));
    ::fmpz_mat_hnf(hnf.raw(), transpose.raw());

    const slong rank = trimmed_nonzero_rows(flint::FmpzMatConstRef(hnf));
    if (rank != rows) {
        return false;
    }

    flint::FmpzMat square(rows, rows);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < rows; ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(square, i, j),
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(hnf), i, j));
        }
    }

    flint::fmpz_mat_det(flint::FmpzRef(out),
                        flint::FmpzMatConstRef(square));
    fmpz_abs(out.raw(), out.raw());
    return flint::fmpz_sgn(flint::FmpzConstRef(out)) > 0;
}

bool arb_abs_upper_exponent_bound(slong& out,
                                  const arb_struct* value,
                                  slong precision) noexcept {
    if (precision <= 0 || ::arb_is_finite(value) == 0) {
        return false;
    }

    flint::Arf upper;
    ::arb_get_abs_ubound_arf(upper.raw(), value, precision);
    if (!flint::arf_is_finite(upper)) {
        return false;
    }
    if (arf_is_zero(upper.raw())) {
        out = std::numeric_limits<slong>::min();
        return true;
    }

    out = ::arf_abs_bound_lt_2exp_si(upper.raw());
    return true;
}

slong fmpz_exponent_bound(const fmpz* value) noexcept {
    if (fmpz_is_zero(value) != 0) {
        return std::numeric_limits<slong>::min();
    }
    const flint_bitcnt_t bits = fmpz_bits(value);
    if (bits == 0) {
        return std::numeric_limits<slong>::min();
    }
    if (bits > static_cast<flint_bitcnt_t>(
                       std::numeric_limits<slong>::max())) {
        return std::numeric_limits<slong>::max();
    }
    return static_cast<slong>(bits) - 1;
}

bool rational_approximation_error_bits_ok(
        const flint::FmpzMat& integer_coordinates,
        const std::vector<flint::Fmpq>& approximations,
        const flint::ArbMat& coordinates,
        flint::FmpzConstRef common_denominator,
        slong precision) noexcept {
    const slong rank = flint::arb_mat_nrows_value(coordinates);
    const slong columns = flint::arb_mat_ncols_value(coordinates);
    if (precision <= 0 || rank <= 0 || columns <= 0 ||
        flint::fmpz_sgn(common_denominator) <= 0 ||
        flint::fmpz_mat_nrows(integer_coordinates) != rank ||
        flint::fmpz_mat_ncols(integer_coordinates) != columns ||
        approximations.size() !=
                static_cast<std::size_t>(rank * columns)) {
        return false;
    }

    slong max_error_exponent = std::numeric_limits<slong>::min();
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < columns; ++j) {
            const flint::Fmpq& approximation =
                    approximations[static_cast<std::size_t>(i * columns + j)];
            flint::Arb approximation_arb;
            flint::Arb difference;
            flint::arb_set_fmpq(approximation_arb, approximation, precision);
            arb_sub(difference.raw(), approximation_arb.raw(),
                    arb_mat_entry(coordinates.raw(), i, j), precision);

            slong exponent = std::numeric_limits<slong>::min();
            if (!arb_abs_upper_exponent_bound(exponent, difference.raw(),
                                              precision)) {
                return false;
            }
            if (exponent > max_error_exponent) {
                max_error_exponent = exponent;
            }
        }
    }

    if (max_error_exponent == std::numeric_limits<slong>::min()) {
        return true;
    }

    slong max_integer_exponent = std::numeric_limits<slong>::min();
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < columns; ++j) {
            const slong exponent = fmpz_exponent_bound(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(integer_coordinates),
                            i, j).raw());
            if (exponent > max_integer_exponent) {
                max_integer_exponent = exponent;
            }
        }
    }
    const slong denominator_exponent =
            fmpz_exponent_bound(common_denominator.raw());
    if (max_integer_exponent == std::numeric_limits<slong>::min() ||
        denominator_exponent == std::numeric_limits<slong>::min()) {
        return true;
    }

    // The source regulator reconstruction rejects a dubious approximation when
    // gexpo(L) + expi(den) > bit - 32, with bit = -gexpo(rational_approximation-lambda).
    if (max_error_exponent > -32) {
        return false;
    }
    const slong budget = -32 - max_error_exponent;
    return max_integer_exponent <= budget &&
           denominator_exponent <= budget - max_integer_exponent;
}

}  // namespace

bool reduced_regulator_from_coordinates(
        flint::ArbRef out,
        const flint::ArbMat& coordinates,
        flint::ArbConstRef regulator_multiple,
        flint::ArbConstRef z,
        flint::FmpzMat* integer_coordinates_out,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong rank = flint::arb_mat_nrows_value(coordinates);
    const slong columns = flint::arb_mat_ncols_value(coordinates);
    if (precision <= 0 ||
        rank <= 0 || columns <= 0 ||
        !flint::arb_is_finite(regulator_multiple) ||
        !flint::arb_is_positive(regulator_multiple) ||
        !flint::arb_is_finite(z) || !flint::arb_is_positive(z)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: invalid rational reconstruction inputs");
        return false;
    }

    flint::Arb denominator_bound_arb;
    flint::arb_mul(denominator_bound_arb, regulator_multiple, z, precision);
    flint::arb_mul_ui(denominator_bound_arb, denominator_bound_arb, 2,
                      precision);

    flint::Fmpz denominator_bound;
    if (!arb_floor_positive(denominator_bound,
                            flint::ArbConstRef(denominator_bound_arb),
                            precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: denominator bound below one");
        return false;
    }

    flint::Fmpz common_denominator;
    flint::fmpz_one(flint::FmpzRef(common_denominator));

    std::vector<flint::Fmpq> approximations;
    approximations.reserve(static_cast<std::size_t>(rank * columns));
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < columns; ++j) {
            approximations.emplace_back();
            if (!arb_best_denominator_bounded_rational(
                        approximations.back(),
                        flint::arb_mat_entry_ref(coordinates, i, j),
                        flint::FmpzConstRef(denominator_bound))) {
                SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "reconstruct_regulator failed: rational_approximation rational reconstruction");
                return false;
            }
            flint::fmpz_lcm(
                    flint::FmpzRef(common_denominator),
                    flint::FmpzConstRef(common_denominator),
                    flint::fmpq_den_ref(approximations.back()));
        }
    }

    if (flint::fmpz_cmp(flint::FmpzConstRef(common_denominator),
                        flint::FmpzConstRef(denominator_bound)) > 0) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: common denominator exceeds bound");
        return false;
    }

    flint::FmpzMat integer_coordinates(rank, columns);
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < columns; ++j) {
            const flint::Fmpq& approximation =
                    approximations[static_cast<std::size_t>(i * columns + j)];
            flint::Fmpz scaled;
            fmpz_divexact(scaled.raw(), common_denominator.raw(),
                          flint::fmpq_den_ref(approximation).raw());
            flint::fmpz_mul(
                    flint::fmpz_mat_entry(integer_coordinates, i, j),
                    flint::FmpzConstRef(scaled),
                    flint::fmpq_num_ref(approximation));
        }
    }

    if (!rational_approximation_error_bits_ok(
                integer_coordinates, approximations, coordinates,
                flint::FmpzConstRef(common_denominator), precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: dubious rational_approximation precision");
        return false;
    }

    flint::Fmpz hnf_determinant;
    if (!column_lattice_hnf_determinant(hnf_determinant,
                                        integer_coordinates)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: integer coordinate HNF determinant");
        return false;
    }

    flint::Fmpz denominator_power;
    flint::fmpz_pow_ui(flint::FmpzRef(denominator_power),
                       flint::FmpzConstRef(common_denominator),
                       static_cast<ulong>(rank));

    flint::Arb reduced;
    ::arb_mul_fmpz(reduced.raw(), regulator_multiple.raw(),
                   hnf_determinant.raw(), precision);
    flint::arb_div_fmpz(reduced, reduced,
                        flint::FmpzConstRef(denominator_power), precision);
    if (!flint::arb_is_finite(reduced) ||
        !flint::arb_is_positive(reduced)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: reduced regulator invalid");
        return false;
    }

    flint::Arb minimum;
    flint::arb_set_ui(minimum, 1);
    flint::arb_div_ui(minimum, minimum, 8, precision);
    if (::arb_ge(reduced.raw(), minimum.raw()) == 0) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "reconstruct_regulator failed: tentative regulator too small");
        return false;
    }

    flint::arb_set(out, flint::ArbConstRef(reduced));
    if (integer_coordinates_out != nullptr) {
        *integer_coordinates_out = std::move(integer_coordinates);
    }
    return true;
}

}  // namespace silex::detail
