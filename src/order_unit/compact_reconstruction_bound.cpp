#include "compact_reconstruction_bound_internal.hpp"

#include "order_unit_internal.hpp"

#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/signature.hpp>

#include <optional>
#include <utility>

namespace silex::detail {
namespace {

using OptionalOrderBasisTransform = std::optional<flint::FmpzMat>;

bool order_basis_begins_with_one(const Order& order) noexcept {
    const slong degree = order.degree();
    flint::FmpqMat basis(degree, degree);
    if (degree <= 0 || !order.get_basis(flint::FmpqMatRef(basis)) ||
        fmpq_mat_nrows(basis.raw()) != degree ||
        fmpq_mat_ncols(basis.raw()) != degree) {
        return false;
    }

    for (slong column = 0; column < degree; ++column) {
        const slong expected = column == 0 ? 1 : 0;
        if (!flint::fmpq_equal_si(
                    flint::fmpq_mat_entry(flint::FmpqMatConstRef(basis), 0,
                                          column),
                    expected)) {
            return false;
        }
    }
    return true;
}

bool normalize_order_basis(flint::FmpzMat& normalized_to_original,
                           ulong& coordinate_scale_bits,
                           const Order& order) noexcept {
    const slong degree = order.degree();
    const NumberField* field = order.parent();
    if (degree <= 0 || field == nullptr ||
        flint::fmpz_mat_nrows(normalized_to_original) != degree ||
        flint::fmpz_mat_ncols(normalized_to_original) != degree) {
        return false;
    }

    Element one(*field);
    flint::FmpqMat rational_coordinates(1, degree);
    flint::FmpzMat one_coordinates(1, degree);
    if (!one.one() ||
        !order.coordinates(flint::FmpqMatRef(rational_coordinates), one) ||
        fmpq_mat_get_fmpz_mat(one_coordinates.raw(),
                              rational_coordinates.raw()) == 0) {
        return false;
    }

    // Construct V by exact column operations so that c V = e_1, where c is
    // the coordinate row of 1 in the supplied order basis. Then U = V^-1 is
    // unimodular, its first row is c, and U times the supplied basis begins
    // with the element 1.
    flint::FmpzMat current(1, degree);
    flint::FmpzMat original_coordinates(1, degree);
    flint::FmpzMat reduction(degree, degree);
    flint::fmpz_mat_set(flint::FmpzMatRef(current),
                        flint::FmpzMatConstRef(one_coordinates));
    flint::fmpz_mat_set(flint::FmpzMatRef(original_coordinates),
                        flint::FmpzMatConstRef(one_coordinates));
    flint::fmpz_mat_one(flint::FmpzMatRef(reduction));

    flint::Fmpz a;
    flint::Fmpz b;
    flint::Fmpz gcd;
    flint::Fmpz bezout_a;
    flint::Fmpz bezout_b;
    flint::Fmpz negative_b_over_gcd;
    flint::Fmpz a_over_gcd;
    flint::Fmpz old_first;
    flint::Fmpz old_other;
    flint::Fmpz new_first;
    flint::Fmpz new_other;
    for (slong column = 1; column < degree; ++column) {
        flint::fmpz_set(
                flint::FmpzRef(b),
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(current), 0,
                                      column));
        if (flint::fmpz_is_zero(flint::FmpzConstRef(b))) {
            continue;
        }
        flint::fmpz_set(
                flint::FmpzRef(a),
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(current), 0, 0));
        ::fmpz_xgcd_canonical_bezout(gcd.raw(), bezout_a.raw(),
                                      bezout_b.raw(), a.raw(), b.raw());
        if (flint::fmpz_sgn(flint::FmpzConstRef(gcd)) <= 0) {
            return false;
        }
        flint::fmpz_divexact(flint::FmpzRef(negative_b_over_gcd),
                             flint::FmpzConstRef(b),
                             flint::FmpzConstRef(gcd));
        flint::fmpz_neg(flint::FmpzRef(negative_b_over_gcd),
                        flint::FmpzConstRef(negative_b_over_gcd));
        flint::fmpz_divexact(flint::FmpzRef(a_over_gcd),
                             flint::FmpzConstRef(a),
                             flint::FmpzConstRef(gcd));

        for (slong row = 0; row < degree; ++row) {
            flint::fmpz_set(
                    flint::FmpzRef(old_first),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(reduction), row, 0));
            flint::fmpz_set(
                    flint::FmpzRef(old_other),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(reduction), row, column));

            flint::fmpz_mul(flint::FmpzRef(new_first),
                            flint::FmpzConstRef(old_first),
                            flint::FmpzConstRef(bezout_a));
            flint::fmpz_addmul(flint::FmpzRef(new_first),
                               flint::FmpzConstRef(old_other),
                               flint::FmpzConstRef(bezout_b));
            flint::fmpz_mul(flint::FmpzRef(new_other),
                            flint::FmpzConstRef(old_first),
                            flint::FmpzConstRef(negative_b_over_gcd));
            flint::fmpz_addmul(flint::FmpzRef(new_other),
                               flint::FmpzConstRef(old_other),
                               flint::FmpzConstRef(a_over_gcd));
            flint::fmpz_set(
                    flint::fmpz_mat_entry(flint::FmpzMatRef(reduction), row,
                                          0),
                    flint::FmpzConstRef(new_first));
            flint::fmpz_set(
                    flint::fmpz_mat_entry(flint::FmpzMatRef(reduction), row,
                                          column),
                    flint::FmpzConstRef(new_other));
        }

        flint::fmpz_set(
                flint::fmpz_mat_entry(flint::FmpzMatRef(current), 0, 0),
                flint::FmpzConstRef(gcd));
        flint::fmpz_zero(
                flint::fmpz_mat_entry(flint::FmpzMatRef(current), 0,
                                      column));
    }

    const auto reduced_first = flint::fmpz_mat_entry(
            flint::FmpzMatConstRef(current), 0, 0);
    if (flint::fmpz_equal_si(reduced_first, -1)) {
        for (slong row = 0; row < degree; ++row) {
            auto entry = flint::fmpz_mat_entry(
                    flint::FmpzMatRef(reduction), row, 0);
            ::fmpz_neg(entry.raw(), entry.raw());
        }
    } else if (!flint::fmpz_is_one(reduced_first)) {
        return false;
    }

    flint::Fmpz inverse_denominator;
    if (!flint::fmpz_mat_inv(flint::FmpzMatRef(normalized_to_original),
                             flint::FmpzRef(inverse_denominator),
                             flint::FmpzMatConstRef(reduction))) {
        return false;
    }
    if (flint::fmpz_equal_si(flint::FmpzConstRef(inverse_denominator), -1)) {
        flint::fmpz_mat_neg(flint::FmpzMatRef(normalized_to_original),
                            flint::FmpzMatConstRef(normalized_to_original));
    } else if (!flint::fmpz_is_one(
                       flint::FmpzConstRef(inverse_denominator))) {
        return false;
    }

    for (slong column = 0; column < degree; ++column) {
        if (!flint::fmpz_equal(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(normalized_to_original), 0,
                            column),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(original_coordinates), 0,
                            column))) {
            return false;
        }
    }

    // If normalized coordinates satisfy |x'_i| <= C, then original
    // coordinates x = x' U satisfy |x_j| <= C sum_i |U_ij|. Round the
    // largest column sum up to a power of two so the public bound remains a
    // power of two and reconstruction in the original basis is certified.
    flint::Fmpz column_sum;
    flint::Fmpz maximum_column_sum;
    flint::Fmpz absolute_entry;
    for (slong column = 0; column < degree; ++column) {
        flint::fmpz_zero(flint::FmpzRef(column_sum));
        for (slong row = 0; row < degree; ++row) {
            flint::fmpz_abs(
                    flint::FmpzRef(absolute_entry),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(normalized_to_original), row,
                            column));
            flint::fmpz_add(flint::FmpzRef(column_sum),
                            flint::FmpzConstRef(column_sum),
                            flint::FmpzConstRef(absolute_entry));
        }
        if (flint::fmpz_cmp(flint::FmpzConstRef(column_sum),
                            flint::FmpzConstRef(maximum_column_sum)) > 0) {
            flint::fmpz_set(flint::FmpzRef(maximum_column_sum),
                            flint::FmpzConstRef(column_sum));
        }
    }
    if (flint::fmpz_sgn(flint::FmpzConstRef(maximum_column_sum)) <= 0) {
        return false;
    }
    coordinate_scale_bits = ::fmpz_bits(maximum_column_sum.raw());
    if (::fmpz_popcnt(maximum_column_sum.raw()) == 1) {
        --coordinate_scale_bits;
    }
    return true;
}

bool validate_inputs(const Order& order, FactoredElementSpan generators,
                     EmbeddingContext& embeddings,
                     const CompactCoordinateBoundOptions& options,
                     Signature& sig, flint::Fmpz& absolute_discriminant,
                     OptionalOrderBasisTransform& normalized_to_original,
                     ulong& coordinate_scale_bits,
                     bool& basis_is_unchanged,
                     CompactCoordinateBoundStatus& status) noexcept {
    status = CompactCoordinateBoundStatus::invalid_input;
    const NumberField* field = order.parent();
    if (!order.is_defined() || field == nullptr || !order.has_basis() ||
        !order.maximality_known() || !order.is_maximal() ||
        options.start_precision <= 0 ||
        options.max_precision < options.start_precision ||
        !signature(sig, *field) || sig.degree() != order.degree()) {
        return false;
    }

    basis_is_unchanged = order_basis_begins_with_one(order);
    if (basis_is_unchanged) {
        coordinate_scale_bits = 0;
    } else {
        normalized_to_original.emplace(order.degree(), order.degree());
        if (!normalize_order_basis(*normalized_to_original,
                                   coordinate_scale_bits, order)) {
            status = CompactCoordinateBoundStatus::unsupported_order_basis;
            return false;
        }
    }

    const slong rank = sig.r1() + sig.r2() - 1;
    if (rank < 0 || generators.size() != static_cast<std::size_t>(rank)) {
        return false;
    }
    for (const FactoredElement& generator : generators) {
        if (!generator.is_defined() || generator.parent() == nullptr ||
            !generator.parent()->has_same_data(*field)) {
            return false;
        }
    }

    if (embeddings.is_defined()) {
        if (embeddings.parent() == nullptr ||
            !embeddings.parent()->has_same_data(*field)) {
            return false;
        }
    } else if (!embeddings.define(*field)) {
        return false;
    }

    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(discriminant)) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(absolute_discriminant),
                    flint::FmpzConstRef(discriminant));
    return true;
}

bool compute_log2_minkowski_bound(flint::Arb& out, const Order& order,
                                  const Signature& sig,
                                  flint::FmpzConstRef absolute_discriminant,
                                  const OptionalOrderBasisTransform&
                                          normalized_to_original,
                                  bool basis_is_unchanged,
                                  const flint::Arb& log_two,
                                  slong precision) noexcept {
    const slong degree = order.degree();
    flint::ArbMat rows(degree, degree);
    if (basis_is_unchanged) {
        if (!order_minkowski_embedding_rows(flint::ArbMatRef(rows), order,
                                            precision, nullptr)) {
            return false;
        }
    } else {
        if (!normalized_to_original.has_value() ||
            flint::fmpz_mat_nrows(*normalized_to_original) != degree ||
            flint::fmpz_mat_ncols(*normalized_to_original) != degree) {
            return false;
        }
        const flint::FmpzMatConstRef normalization(
                *normalized_to_original);
        flint::ArbMat original_rows(degree, degree);
        if (!order_minkowski_embedding_rows(flint::ArbMatRef(original_rows),
                                            order, precision, nullptr)) {
            return false;
        }
        for (slong row = 0; row < degree; ++row) {
            for (slong column = 0; column < degree; ++column) {
                flint::ArbRef entry =
                        flint::arb_mat_entry_ref(rows, row, column);
                flint::arb_zero(entry);
                for (slong k = 0; k < degree; ++k) {
                    flint::arb_addmul_fmpz(
                            entry,
                            flint::arb_mat_entry_ref(
                                    flint::ArbMatConstRef(original_rows), k,
                                    column),
                            flint::fmpz_mat_entry(normalization, row, k),
                            precision);
                }
            }
        }
    }
    flint::Arb e;
    flint::Arb norm_squared;
    flint::Arb entry;
    flint::Arb term;
    flint::Arb d;
    flint::arb_log_ui(e, static_cast<ulong>(degree), precision);
    flint::arb_div(e, e, log_two, precision);

    // The source compact-bound loop skips G's first column, which is the basis
    // element 1, and applies Hadamard to the remaining squared T2 norms.
    for (slong row = 1; row < degree; ++row) {
        flint::arb_zero(norm_squared);
        for (slong column = 0; column < degree; ++column) {
            flint::arb_set(entry, arb_mat_entry(rows.raw(), row, column));
            flint::arb_sqr(term, entry, precision);
            flint::arb_add(norm_squared, norm_squared, term, precision);
        }
        if (!flint::arb_is_finite(norm_squared) ||
            !flint::arb_is_positive(norm_squared)) {
            return false;
        }
        flint::arb_log(term, norm_squared, precision);
        flint::arb_div(term, term, log_two, precision);
        flint::arb_add(e, e, term, precision);
    }
    flint::arb_div_ui(e, e, 2, precision);

    flint::arb_log_fmpz(d, absolute_discriminant, precision);
    flint::arb_div(d, d, log_two, precision);
    flint::arb_div_ui(d, d, 2, precision);

    // This intentionally reproduces the literal reference implementation expression
    // dbllog2(D)/2 - r2*M_LN2. It is not independently renormalized to r2.
    flint::arb_mul_ui(term, log_two, static_cast<ulong>(sig.r2()), precision);
    flint::arb_sub(d, d, term, precision);
    flint::arb_sub(out, e, d, precision);
    return flint::arb_is_finite(out);
}

bool compute_log2_unit_bound(flint::Arb& out, EmbeddingContext& embeddings,
                             FactoredElementSpan generators,
                             const Signature& sig, const flint::Arb& log_two,
                             slong precision) noexcept {
    const slong places = sig.r1() + sig.r2();
    flint::ArbMat logs(static_cast<slong>(generators.size()), places);
    if (!compact_log_matrix(logs, embeddings, generators, precision)) {
        return false;
    }

    flint::Arb maximum;
    flint::Arb entry;
    flint::arb_zero(maximum);
    for (slong row = 0; row < static_cast<slong>(generators.size()); ++row) {
        for (slong column = 0; column < places; ++column) {
            flint::arb_set(entry, arb_mat_entry(logs.raw(), row, column));
            // Product-mode complex logs are 2*log|sigma(u)|, matching reference's
            // division by two for the complex places in log2fubound.
            if (column >= sig.r1()) {
                flint::arb_div_ui(entry, entry, 2, precision);
            }
            flint::arb_max(maximum, maximum, entry, precision);
        }
    }
    flint::arb_div(out, maximum, log_two, precision);
    return flint::arb_is_finite(out);
}

bool compute_bound_interval(CompactCoordinateBoundReport& report,
                            const Order& order, FactoredElementSpan generators,
                            EmbeddingContext& embeddings, const Signature& sig,
                            flint::FmpzConstRef absolute_discriminant,
                            const OptionalOrderBasisTransform&
                                    normalized_to_original,
                            ulong coordinate_scale_bits,
                            bool basis_is_unchanged,
                            slong precision) noexcept {
    flint::Arb log_two;
    flint::arb_log_ui(log_two, 2, precision);
    if (!flint::arb_is_finite(log_two) || !flint::arb_is_positive(log_two) ||
        !compute_log2_minkowski_bound(report.log2_minkowski_bound, order, sig,
                                      absolute_discriminant,
                                      normalized_to_original,
                                      basis_is_unchanged, log_two, precision) ||
        !compute_log2_unit_bound(report.log2_unit_bound, embeddings, generators,
                                 sig, log_two, precision)) {
        return false;
    }

    flint::arb_add(report.log2_coordinate_bound, report.log2_minkowski_bound,
                   report.log2_unit_bound, precision);
    if (coordinate_scale_bits != 0) {
        flint::arb_add_ui(report.log2_coordinate_bound,
                          report.log2_coordinate_bound,
                          coordinate_scale_bits, precision);
    }
    return flint::arb_is_finite(report.log2_coordinate_bound);
}

bool stable_ceiling(flint::Fmpz& ceiling, const flint::Arb& value,
                    slong precision) noexcept {
    flint::Arf lower;
    flint::Arf upper;
    flint::Fmpz lower_ceiling;
    flint::Fmpz upper_ceiling;
    flint::arb_get_lbound_arf(lower, value, precision);
    flint::arb_get_ubound_arf(upper, value, precision);
    if (!flint::arf_is_finite(lower) || !flint::arf_is_finite(upper)) {
        return false;
    }
    flint::arf_get_fmpz(lower_ceiling, lower, ARF_RND_CEIL);
    flint::arf_get_fmpz(upper_ceiling, upper, ARF_RND_CEIL);
    if (!flint::fmpz_equal(lower_ceiling, upper_ceiling)) {
        return false;
    }
    flint::fmpz_set(flint::FmpzRef(ceiling),
                    flint::FmpzConstRef(lower_ceiling));
    return true;
}

void publish_status(CompactCoordinateBoundReport& report,
                    CompactCoordinateBoundStatus status) noexcept {
    CompactCoordinateBoundReport result;
    result.status = status;
    report = std::move(result);
}

}  // namespace

bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options) noexcept {
    Signature sig;
    flint::Fmpz absolute_discriminant;
    OptionalOrderBasisTransform normalized_to_original;
    ulong coordinate_scale_bits = 0;
    bool basis_is_unchanged = false;
    CompactCoordinateBoundStatus validation_status;
    if (!validate_inputs(order, generators, embeddings, options, sig,
                         absolute_discriminant, normalized_to_original,
                         coordinate_scale_bits, basis_is_unchanged,
                         validation_status)) {
        publish_status(report, validation_status);
        return false;
    }

    slong precision = options.start_precision;
    slong attempts = 0;
    for (;;) {
        ++attempts;
        CompactCoordinateBoundReport candidate;
        candidate.precision = precision;
        candidate.attempts = attempts;
        flint::Fmpz bit_bound;
        if (compute_bound_interval(
                    candidate, order, generators, embeddings, sig,
                    flint::FmpzConstRef(absolute_discriminant),
                    normalized_to_original, coordinate_scale_bits,
                    basis_is_unchanged, precision) &&
            stable_ceiling(bit_bound, candidate.log2_coordinate_bound,
                           precision)) {
            if (flint::fmpz_sgn(flint::FmpzConstRef(bit_bound)) < 0 ||
                !flint::fmpz_abs_fits_ui(flint::FmpzConstRef(bit_bound))) {
                candidate.status =
                        CompactCoordinateBoundStatus::bound_overflow;
                report = std::move(candidate);
                return false;
            }

            candidate.bit_bound =
                    flint::fmpz_get_ui(flint::FmpzConstRef(bit_bound));
            flint::fmpz_one(flint::FmpzRef(candidate.coordinate_bound));
            flint::fmpz_mul_2exp(
                    flint::FmpzRef(candidate.coordinate_bound),
                    flint::FmpzConstRef(candidate.coordinate_bound),
                    candidate.bit_bound);
            candidate.status = CompactCoordinateBoundStatus::success;
            flint::fmpz_set(out,
                            flint::FmpzConstRef(candidate.coordinate_bound));
            report = std::move(candidate);
            return true;
        }

        if (precision >= options.max_precision) {
            candidate.status =
                    CompactCoordinateBoundStatus::precision_exhausted;
            report = std::move(candidate);
            return false;
        }
        if (precision > options.max_precision / 2) {
            precision = options.max_precision;
        } else {
            precision *= 2;
        }
    }
}

}  // namespace silex::detail
