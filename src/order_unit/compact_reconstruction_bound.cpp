#include "compact_reconstruction_bound_internal.hpp"

#include "order_unit_internal.hpp"
#include "../factored_element/wide_factored_internal.hpp"

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

// The units whose logarithms feed log2fubound: word-size factored
// generators, or a batch with arbitrary-size exponents.
struct UnitLogSource {
    FactoredElementSpan generators;
    const WideFactoredUnits* wide = nullptr;

    std::size_t count() const noexcept {
        return wide == nullptr
                       ? generators.size()
                       : static_cast<std::size_t>(
                                 flint::fmpz_mat_nrows(wide->exponents));
    }
};

bool unit_source_has_field(const UnitLogSource& source,
                           const NumberField& field) noexcept {
    if (source.wide == nullptr) {
        for (const FactoredElement& generator : source.generators) {
            if (!generator.is_defined() || generator.parent() == nullptr ||
                !generator.parent()->has_same_data(field)) {
                return false;
            }
        }
        return true;
    }
    if (flint::fmpz_mat_ncols(source.wide->exponents) !=
        static_cast<slong>(source.wide->bases.size())) {
        return false;
    }
    for (const Element& base : source.wide->bases) {
        if (!base.is_defined() || !base.has_parent(field)) {
            return false;
        }
    }
    return true;
}

// Product-mode log rows of the source units.  A wide unit's row is
// sum_k e_k * log(base_k); Arb tracks the cancellation in that sum, and the
// caller raises the precision until the bound is stable.
bool unit_source_log_matrix(flint::ArbMat& out, EmbeddingContext& embeddings,
                            const UnitLogSource& source,
                            slong precision) noexcept {
    if (source.wide == nullptr) {
        return compact_log_matrix(out, embeddings, source.generators,
                                  precision);
    }

    const WideFactoredUnits& wide = *source.wide;
    const NumberField* field = embeddings.parent();
    slong places = 0;
    if (field == nullptr || !compact_places(places, embeddings) ||
        flint::arb_mat_nrows_value(out) !=
                flint::fmpz_mat_nrows(wide.exponents) ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    std::vector<FactoredElement> bases;
    bases.reserve(wide.bases.size());
    for (const Element& base : wide.bases) {
        bases.emplace_back(*field);
        if (!bases.back().is_defined() || !bases.back().set_element(base)) {
            return false;
        }
    }
    flint::ArbMat base_logs(static_cast<slong>(bases.size()), places);
    if (!compact_log_matrix(base_logs, embeddings,
                            FactoredElementSpan(bases.data(), bases.size()),
                            precision)) {
        return false;
    }

    const slong rows = flint::fmpz_mat_nrows(wide.exponents);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < places; ++j) {
            arb_struct* entry = arb_mat_entry(out.raw(), i, j);
            arb_zero(entry);
            for (slong k = 0; k < static_cast<slong>(bases.size()); ++k) {
                arb_addmul_fmpz(entry, arb_mat_entry(base_logs.raw(), k, j),
                                fmpz_mat_entry(wide.exponents.raw(), i, k),
                                precision);
            }
        }
    }
    return true;
}

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

bool validate_inputs(const Order& order, const UnitLogSource& generators,
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
    if (rank < 0 || generators.count() != static_cast<std::size_t>(rank) ||
        !unit_source_has_field(generators, *field)) {
        return false;
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

// Rows of the weighted Minkowski embedding (the reference's T2 matrix G, one
// row per basis element) of the basis that begins with 1.
bool normalized_minkowski_rows(flint::ArbMat& rows, const Order& order,
                               const OptionalOrderBasisTransform&
                                       normalized_to_original,
                               bool basis_is_unchanged,
                               slong precision) noexcept {
    const slong degree = order.degree();
    if (flint::arb_mat_nrows_value(rows) != degree ||
        flint::arb_mat_ncols_value(rows) != degree) {
        return false;
    }
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
    return true;
}

bool compute_log2_minkowski_bound(flint::Arb& out, const flint::ArbMat& rows,
                                  const Signature& sig,
                                  flint::FmpzConstRef absolute_discriminant,
                                  const flint::Arb& log_two,
                                  slong precision) noexcept {
    const slong degree = flint::arb_mat_nrows_value(rows);
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

// log2 of the quantity the reference's log2Mbound bounds (second reference
// buch2.c:1201): the largest ratio |x|_oo / |y|_oo over the solutions x of
// the split real/imaginary system S x = y (RgM_solve_realimag), where S has
// the real parts of all r1 + r2 embeddings and the imaginary parts of the r2
// complex ones, unweighted.  That ratio is at most the induced max-row-sum
// norm ||S^-1||_oo, and log2Mbound is the Hadamard and Cramer estimate of
// it; this evaluates the norm itself from a certified Arb inverse.
//
// `rows` is G^T: row j holds sqrt(w_k) times the embedding coordinate k of
// basis element j, with w_k = 1 for the r1 real columns and w_k = 2 for the
// interleaved real and imaginary columns of the complex places.  If c is
// the coordinate row of an element with split embedding row y, then
// c rows = y W^(1/2), so c = y W^(1/2) rows^-1 and
//   |c_j| <= |y|_oo sum_k sqrt(w_k) |(rows^-1)_kj|,
// a column sum of rows^-1 that equals the row sum of S^-1.  Every Arb
// operation encloses its exact value, so the returned ball encloses the
// exact norm and its upper end is an outward-rounded bound.
bool compute_log2_inverse_norm_bound(flint::Arb& out,
                                     const flint::ArbMat& rows,
                                     const Signature& sig,
                                     const flint::Arb& log_two,
                                     slong precision) noexcept {
    const slong degree = flint::arb_mat_nrows_value(rows);
    flint::ArbMat inverse(degree, degree);
    if (degree <= 0 || flint::arb_mat_ncols_value(rows) != degree ||
        arb_mat_inv(inverse.raw(), rows.raw(), precision) == 0) {
        return false;
    }

    flint::Arb sqrt_two;
    flint::Arb column_sum;
    flint::Arb entry;
    flint::Arb maximum;
    flint::arb_sqrt_ui(sqrt_two, 2, precision);
    flint::arb_zero(maximum);
    for (slong column = 0; column < degree; ++column) {
        flint::arb_zero(column_sum);
        for (slong k = 0; k < degree; ++k) {
            flint::arb_set(entry, arb_mat_entry(inverse.raw(), k, column));
            flint::arb_abs(entry, entry);
            if (k >= sig.r1()) {
                flint::arb_mul(entry, entry, sqrt_two, precision);
            }
            flint::arb_add(column_sum, column_sum, entry, precision);
        }
        flint::arb_max(maximum, maximum, column_sum, precision);
    }
    if (!flint::arb_is_finite(maximum) || !flint::arb_is_positive(maximum)) {
        return false;
    }
    flint::arb_log(out, maximum, precision);
    flint::arb_div(out, out, log_two, precision);
    return flint::arb_is_finite(out);
}

bool compute_log2_unit_bound(flint::Arb& out, EmbeddingContext& embeddings,
                             const UnitLogSource& generators,
                             const Signature& sig, const flint::Arb& log_two,
                             slong precision) noexcept {
    const slong places = sig.r1() + sig.r2();
    flint::ArbMat logs(static_cast<slong>(generators.count()), places);
    if (!unit_source_log_matrix(logs, embeddings, generators, precision)) {
        return false;
    }

    flint::Arb maximum;
    flint::Arb entry;
    flint::arb_zero(maximum);
    for (slong row = 0; row < static_cast<slong>(generators.count()); ++row) {
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
                            const Order& order,
                            const UnitLogSource& generators,
                            EmbeddingContext& embeddings, const Signature& sig,
                            flint::FmpzConstRef absolute_discriminant,
                            const OptionalOrderBasisTransform&
                                    normalized_to_original,
                            ulong coordinate_scale_bits,
                            bool basis_is_unchanged,
                            slong precision) noexcept {
    flint::Arb log_two;
    flint::arb_log_ui(log_two, 2, precision);
    const slong degree = order.degree();
    flint::ArbMat rows(degree, degree);
    if (!flint::arb_is_finite(log_two) || !flint::arb_is_positive(log_two) ||
        !normalized_minkowski_rows(rows, order, normalized_to_original,
                                   basis_is_unchanged, precision) ||
        !compute_log2_minkowski_bound(report.log2_minkowski_bound, rows, sig,
                                      absolute_discriminant, log_two,
                                      precision) ||
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
    if (!flint::arb_is_finite(report.log2_coordinate_bound)) {
        return false;
    }

    // An inverse that is not certified at this precision leaves the
    // operator-norm bound indeterminate, and the caller then publishes the
    // reference expression; that expression alone decides success and the
    // precision ladder.
    if (compute_log2_inverse_norm_bound(report.log2_inverse_norm_bound, rows,
                                        sig, log_two, precision)) {
        flint::arb_add(report.log2_inverse_norm_coordinate_bound,
                       report.log2_inverse_norm_bound, report.log2_unit_bound,
                       precision);
        if (coordinate_scale_bits != 0) {
            flint::arb_add_ui(report.log2_inverse_norm_coordinate_bound,
                              report.log2_inverse_norm_coordinate_bound,
                              coordinate_scale_bits, precision);
        }
    } else {
        arb_indeterminate(report.log2_inverse_norm_bound.raw());
        arb_indeterminate(report.log2_inverse_norm_coordinate_bound.raw());
    }
    return true;
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

namespace {

bool compact_unit_coordinate_bound_from_source(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, const UnitLogSource& generators,
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

            candidate.hadamard_bit_bound =
                    flint::fmpz_get_ui(flint::FmpzConstRef(bit_bound));
            candidate.bit_bound = candidate.hadamard_bit_bound;
            // Publish the operator-norm bound whenever it is certified: the
            // Arb inverse succeeded at this precision and its total has a
            // stable, nonnegative ceiling that fits a word.  It bounds the
            // coordinates directly, for every signature.  Only otherwise
            // publish the reference expression.  The operator norm is not
            // capped by the reference expression: for r2 = 0 it never
            // exceeds it (Cramer's rule, Cauchy-Schwarz and Hadamard's
            // inequality, with the column of 1 of smallest norm,
            // sqrt(n) <= sqrt(T2(w)) for a nonzero integer w), but for
            // r2 > 0 the reference subtracts r2 * ln 2 for
            // log2 |det S| = log2 |D| / 2 - r2 and that argument does not
            // cover it, so there the operator norm may be the larger.  The
            // reference expression alone still decides success, the
            // precision ladder and failure.
            flint::Fmpz inverse_norm_bit_bound;
            if (flint::arb_is_finite(
                        candidate.log2_inverse_norm_coordinate_bound) &&
                stable_ceiling(inverse_norm_bit_bound,
                               candidate.log2_inverse_norm_coordinate_bound,
                               precision) &&
                flint::fmpz_sgn(
                        flint::FmpzConstRef(inverse_norm_bit_bound)) >= 0 &&
                flint::fmpz_abs_fits_ui(
                        flint::FmpzConstRef(inverse_norm_bit_bound))) {
                candidate.bit_bound = flint::fmpz_get_ui(
                        flint::FmpzConstRef(inverse_norm_bit_bound));
                candidate.inverse_norm_bound_used = true;
            }
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

}  // namespace

bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options) noexcept {
    return compact_unit_coordinate_bound_from_source(
            report, out, order, UnitLogSource{generators, nullptr},
            embeddings, options);
}

bool compact_unit_coordinate_bound(
        CompactCoordinateBoundReport& report, flint::FmpzRef out,
        const Order& order, const WideFactoredUnits& generators,
        EmbeddingContext& embeddings,
        const CompactCoordinateBoundOptions& options) noexcept {
    return compact_unit_coordinate_bound_from_source(
            report, out, order,
            UnitLogSource{FactoredElementSpan(), &generators}, embeddings,
            options);
}

}  // namespace silex::detail
