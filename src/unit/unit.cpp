#include <silex/unit.hpp>

#include "../element/element_internal.hpp"
#include "unit_internal.hpp"

#include <flint/arb_mat.h>

#include <silex/flint/arf.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/signature.hpp>

namespace silex {
namespace {

bool divide_exact_checked(flint::FmpzRef out,
                          flint::FmpzConstRef numerator,
                          flint::FmpzConstRef denominator) noexcept {
    if (flint::fmpz_sgn(denominator) <= 0) {
        return false;
    }
    return ::fmpz_divides(out.raw(), numerator.raw(), denominator.raw()) != 0;
}

bool set_quadratic_coefficients(Element& out,
                                flint::FmpzConstRef constant,
                                flint::FmpzConstRef linear,
                                bool half_integral) noexcept {
    flint::FmpqPoly polynomial;
    flint::Fmpq coeff;

    flint::fmpq_poly_zero(polynomial);
    if (half_integral) {
        flint::fmpq_set_fmpz(coeff, constant);
        flint::fmpq_div_2exp(coeff, coeff, 1);
        flint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coeff);
        flint::fmpq_set_fmpz(coeff, linear);
        flint::fmpq_div_2exp(coeff, coeff, 1);
        flint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coeff);
    } else {
        flint::fmpq_set_fmpz(coeff, constant);
        flint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coeff);
        flint::fmpq_set_fmpz(coeff, linear);
        flint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coeff);
    }

    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool is_positive_squarefree_nonsquare(
        flint::FmpzConstRef value) noexcept {
    if (flint::fmpz_sgn(value) <= 0 || flint::fmpz_is_square(value)) {
        return false;
    }

    flint::FmpzFactor factorization;
    flint::fmpz_factor(flint::FmpzFactorRef(factorization), value);
    for (slong i = 0;
         i < flint::fmpz_factor_num(
                     flint::FmpzFactorConstRef(factorization));
         ++i) {
        if (flint::fmpz_factor_exp(
                    flint::FmpzFactorConstRef(factorization), i) != 1) {
            return false;
        }
    }
    return true;
}

bool pure_quadratic_defining_radicand(
        flint::FmpzRef out,
        const NumberField& field) noexcept {
    // reference `quadunit` is presentation-independent.  For the exact canonical
    // presentation x^2-d, the FLINT field generator is already sqrt(d), so
    // the existing quadratic coefficients need no generator-change map.
    const nf_struct* raw_field = field.raw_flint_field();
    if (raw_field == nullptr || field.degree() != 2) {
        return false;
    }

    const flint::FmpqPolyConstRef polynomial(raw_field->pol);
    flint::Fmpq constant;
    flint::Fmpq linear;
    flint::Fmpq leading;
    flint::fmpq_poly_get_coeff_fmpq(
            flint::FmpqRef(constant), polynomial, 0);
    flint::fmpq_poly_get_coeff_fmpq(
            flint::FmpqRef(linear), polynomial, 1);
    flint::fmpq_poly_get_coeff_fmpq(
            flint::FmpqRef(leading), polynomial, 2);
    if (!flint::fmpq_equal_si(leading, 1) ||
        !flint::fmpq_equal_si(linear, 0) ||
        !flint::fmpz_is_one(flint::fmpq_den_ref(constant)) ||
        flint::fmpz_sgn(flint::fmpq_num_ref(constant)) >= 0) {
        return false;
    }

    flint::Fmpz radicand;
    flint::fmpz_neg(flint::FmpzRef(radicand),
                    flint::fmpq_num_ref(constant));
    if (!is_positive_squarefree_nonsquare(
                flint::FmpzConstRef(radicand))) {
        return false;
    }

    flint::fmpz_set(out, flint::FmpzConstRef(radicand));
    return true;
}

bool fundamental_unit_radicand(flint::FmpzRef out,
                               const NumberField& field) noexcept {
    flint::Fmpz radicand;
    if (!field.quadratic_radicand(flint::FmpzRef(radicand)) &&
        !pure_quadratic_defining_radicand(
                flint::FmpzRef(radicand), field)) {
        return false;
    }
    if (flint::fmpz_sgn(flint::FmpzConstRef(radicand)) <= 0) {
        return false;
    }
    flint::fmpz_set(out, flint::FmpzConstRef(radicand));
    return true;
}

bool signature_from_embeddings(Signature& out,
                               const EmbeddingContext& embeddings) noexcept {
    if (embeddings.parent() == nullptr) {
        return false;
    }
    out = embeddings.signature();
    if (out.degree() == embeddings.degree()) {
        return true;
    }
    return out.compute(*embeddings.parent());
}

bool valid_unit_span(EmbeddingContext& embeddings,
                     ElementSpan units) noexcept {
    const NumberField* parent = embeddings.parent();
    if (parent == nullptr) {
        return false;
    }
    if (units.size() > 0 && units.data() == nullptr) {
        return false;
    }
    for (const Element& unit : units) {
        if (!unit.has_parent(*parent) || unit.equal_si(0)) {
            return false;
        }
    }
    return true;
}

}  // namespace

bool unit_rank(slong& rank, const NumberField& field) noexcept {
    Signature sig;
    if (!sig.compute(field)) {
        return false;
    }
    rank = sig.r1() + sig.r2() - 1;
    return true;
}

bool unit_lower_regulator_bound(flint::ArbRef out,
                                const NumberField& field,
                                slong precision) noexcept {
    if (!field.is_defined() || precision <= 0) {
        return false;
    }

    Signature sig;
    if (!sig.compute(field)) {
        return false;
    }

    // Every bound used here increases with w, so it is valid for any
    // w <= w_K.  For an order O in K, O^x is a finite-index subgroup of
    // O_K^x, so Reg(O) >= R_K and a lower bound for R_K bounds Reg(O); w is
    // therefore the root-of-unity count of the field, not of the order.  If
    // it cannot be computed, w = 2 (every number field contains -1) keeps
    // the bound valid.
    flint::Fmpz roots_order;
    if (!root_of_unity_order(flint::FmpzRef(roots_order), field)) {
        fmpz_set_ui(roots_order.raw(), 2);
    }

    flint::Arb bound;
    if (!detail::regulator_lower_bound_from_signature(
                flint::ArbRef(bound), sig.r1(), sig.r2(),
                flint::FmpzConstRef(roots_order.raw()), precision)) {
        return false;
    }
    arb_set(out.raw(), bound.raw());
    return true;
}

bool quadratic_fundamental_unit(Element& out,
                                const NumberField& field) noexcept {
    if (!field.is_defined() || !detail::ensure_parent(out, field)) {
        return false;
    }

    flint::Fmpz radicand;
    if (!fundamental_unit_radicand(flint::FmpzRef(radicand), field)) {
        return false;
    }

    const bool odd_discriminant =
            flint::fmpz_fdiv_ui(flint::FmpzConstRef(radicand), 4) == 1;
    flint::Fmpz discriminant;
    if (odd_discriminant) {
        flint::fmpz_set(flint::FmpzRef(discriminant),
                        flint::FmpzConstRef(radicand));
    } else {
        flint::fmpz_mul_ui(flint::FmpzRef(discriminant),
                           flint::FmpzConstRef(radicand), 4);
    }

    flint::Fmpz root_floor;
    flint::Fmpz root_remainder;
    flint::fmpz_sqrtrem(flint::FmpzRef(root_floor),
                        flint::FmpzRef(root_remainder),
                        flint::FmpzConstRef(discriminant));
    if (flint::fmpz_sgn(flint::FmpzConstRef(root_floor)) <= 0 ||
        flint::fmpz_sgn(flint::FmpzConstRef(root_remainder)) <= 0) {
        return false;
    }

    flint::Fmpz reduced_numerator;
    flint::Fmpz previous_denominator;
    flint::Fmpz denominator;
    flint::Fmpz two;
    flint::fmpz_set(flint::FmpzRef(reduced_numerator),
                    flint::FmpzConstRef(root_floor));
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(previous_denominator),
                            flint::FmpzConstRef(root_remainder), 1);
    flint::fmpz_set_ui(flint::FmpzRef(two), 2);
    flint::fmpz_set(flint::FmpzRef(denominator),
                    flint::FmpzConstRef(two));

    const bool root_floor_is_odd =
            flint::fmpz_fdiv_ui(flint::FmpzConstRef(root_floor), 2) == 1;
    if (root_floor_is_odd != odd_discriminant) {
        flint::fmpz_sub_ui(flint::FmpzRef(reduced_numerator),
                           flint::FmpzConstRef(root_floor), 1);
        flint::fmpz_add(flint::FmpzRef(previous_denominator),
                        flint::FmpzConstRef(previous_denominator),
                        flint::FmpzConstRef(root_floor));
    }
    if (flint::fmpz_sgn(flint::FmpzConstRef(previous_denominator)) <= 0) {
        return false;
    }

    flint::Fmpz convergent_constant_previous;
    flint::Fmpz convergent_constant_current;
    flint::Fmpz convergent_linear_previous;
    flint::Fmpz convergent_linear_current;
    flint::fmpz_set_ui(flint::FmpzRef(convergent_constant_previous), 2);
    flint::fmpz_set(flint::FmpzRef(convergent_constant_current),
                    flint::FmpzConstRef(reduced_numerator));
    flint::fmpz_zero(flint::FmpzRef(convergent_linear_previous));
    flint::fmpz_one(flint::FmpzRef(convergent_linear_current));

    flint::Fmpz constant_product;
    flint::Fmpz linear_product;
    flint::Fmpz sum_product;
    flint::Fmpz denominator_snapshot;
    flint::Fmpz dividend;
    flint::Fmpz quotient;
    flint::Fmpz remainder;
    flint::Fmpz old_reduced_numerator;
    flint::Fmpz next_constant_convergent;
    flint::Fmpz next_linear_convergent;
    if (flint::fmpz_equal(flint::FmpzConstRef(previous_denominator),
                          flint::FmpzConstRef(denominator))) {
        flint::fmpz_mul(
                flint::FmpzRef(constant_product),
                flint::FmpzConstRef(convergent_constant_previous),
                flint::FmpzConstRef(convergent_constant_current));
        flint::fmpz_mul(
                flint::FmpzRef(linear_product),
                flint::FmpzConstRef(convergent_linear_previous),
                flint::FmpzConstRef(convergent_linear_current));
        flint::fmpz_add(
                flint::FmpzRef(dividend),
                flint::FmpzConstRef(convergent_constant_previous),
                flint::FmpzConstRef(convergent_linear_previous));
        flint::fmpz_add(
                flint::FmpzRef(remainder),
                flint::FmpzConstRef(convergent_constant_current),
                flint::FmpzConstRef(convergent_linear_current));
        flint::fmpz_mul(flint::FmpzRef(sum_product),
                         flint::FmpzConstRef(dividend),
                         flint::FmpzConstRef(remainder));
    } else {
        denominator.swap(previous_denominator);
        for (;;) {
            denominator.swap(denominator_snapshot);
            if (flint::fmpz_sgn(
                        flint::FmpzConstRef(denominator_snapshot)) <= 0) {
                return false;
            }

            flint::fmpz_add(flint::FmpzRef(dividend),
                            flint::FmpzConstRef(reduced_numerator),
                            flint::FmpzConstRef(root_floor));
            if (flint::fmpz_sgn(flint::FmpzConstRef(dividend)) < 0) {
                return false;
            }
            ::fmpz_fdiv_qr(quotient.raw(), remainder.raw(), dividend.raw(),
                           denominator_snapshot.raw());
            if (flint::fmpz_sgn(flint::FmpzConstRef(quotient)) < 0 ||
                flint::fmpz_sgn(flint::FmpzConstRef(remainder)) < 0 ||
                flint::fmpz_cmp(
                        flint::FmpzConstRef(remainder),
                        flint::FmpzConstRef(denominator_snapshot)) >= 0) {
                return false;
            }

            old_reduced_numerator.swap(reduced_numerator);
            flint::fmpz_sub(flint::FmpzRef(reduced_numerator),
                            flint::FmpzConstRef(root_floor),
                            flint::FmpzConstRef(remainder));
            if (flint::fmpz_equal(
                        flint::FmpzConstRef(old_reduced_numerator),
                        flint::FmpzConstRef(reduced_numerator))) {
                flint::fmpz_mul(
                        flint::FmpzRef(constant_product),
                        flint::FmpzConstRef(convergent_constant_current),
                        flint::FmpzConstRef(convergent_constant_current));
                flint::fmpz_mul(
                        flint::FmpzRef(linear_product),
                        flint::FmpzConstRef(convergent_linear_current),
                        flint::FmpzConstRef(convergent_linear_current));
                flint::fmpz_add(
                        flint::FmpzRef(dividend),
                        flint::FmpzConstRef(convergent_constant_current),
                        flint::FmpzConstRef(convergent_linear_current));
                flint::fmpz_mul(flint::FmpzRef(sum_product),
                                 flint::FmpzConstRef(dividend),
                                 flint::FmpzConstRef(dividend));
                denominator.swap(denominator_snapshot);
                break;
            }

            flint::fmpz_set(
                    flint::FmpzRef(next_constant_convergent),
                    flint::FmpzConstRef(convergent_constant_previous));
            flint::fmpz_addmul(
                    flint::FmpzRef(next_constant_convergent),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(convergent_constant_current));
            convergent_constant_previous.swap(
                    convergent_constant_current);
            convergent_constant_current.swap(next_constant_convergent);

            flint::fmpz_set(
                    flint::FmpzRef(next_linear_convergent),
                    flint::FmpzConstRef(convergent_linear_previous));
            flint::fmpz_addmul(
                    flint::FmpzRef(next_linear_convergent),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(convergent_linear_current));
            convergent_linear_previous.swap(convergent_linear_current);
            convergent_linear_current.swap(next_linear_convergent);

            flint::fmpz_sub(flint::FmpzRef(remainder),
                            flint::FmpzConstRef(reduced_numerator),
                            flint::FmpzConstRef(old_reduced_numerator));
            flint::fmpz_set(flint::FmpzRef(denominator),
                            flint::FmpzConstRef(previous_denominator));
            flint::fmpz_submul(flint::FmpzRef(denominator),
                               flint::FmpzConstRef(quotient),
                               flint::FmpzConstRef(remainder));
            previous_denominator.swap(denominator_snapshot);
            if (flint::fmpz_sgn(
                        flint::FmpzConstRef(denominator)) <= 0) {
                return false;
            }
            if (flint::fmpz_equal(
                        flint::FmpzConstRef(denominator),
                        flint::FmpzConstRef(previous_denominator))) {
                flint::fmpz_mul(
                        flint::FmpzRef(constant_product),
                        flint::FmpzConstRef(convergent_constant_previous),
                        flint::FmpzConstRef(convergent_constant_current));
                flint::fmpz_mul(
                        flint::FmpzRef(linear_product),
                        flint::FmpzConstRef(convergent_linear_previous),
                        flint::FmpzConstRef(convergent_linear_current));
                flint::fmpz_add(
                        flint::FmpzRef(dividend),
                        flint::FmpzConstRef(convergent_constant_previous),
                        flint::FmpzConstRef(convergent_linear_previous));
                flint::fmpz_add(
                        flint::FmpzRef(remainder),
                        flint::FmpzConstRef(convergent_constant_current),
                        flint::FmpzConstRef(convergent_linear_current));
                flint::fmpz_mul(flint::FmpzRef(sum_product),
                                 flint::FmpzConstRef(dividend),
                                 flint::FmpzConstRef(remainder));
                break;
            }
        }
    }

    flint::Fmpz constant_numerator;
    flint::Fmpz linear_numerator;
    flint::Fmpz constant;
    flint::Fmpz linear;
    flint::fmpz_set(flint::FmpzRef(constant_numerator),
                    flint::FmpzConstRef(constant_product));
    flint::fmpz_addmul(flint::FmpzRef(constant_numerator),
                       flint::FmpzConstRef(discriminant),
                       flint::FmpzConstRef(linear_product));
    if (!divide_exact_checked(flint::FmpzRef(constant),
                              flint::FmpzConstRef(constant_numerator),
                              flint::FmpzConstRef(denominator))) {
        return false;
    }

    flint::fmpz_sub(flint::FmpzRef(linear_numerator),
                    flint::FmpzConstRef(sum_product),
                    flint::FmpzConstRef(constant_product));
    flint::fmpz_sub(flint::FmpzRef(linear_numerator),
                    flint::FmpzConstRef(linear_numerator),
                    flint::FmpzConstRef(linear_product));
    if (!divide_exact_checked(flint::FmpzRef(linear),
                              flint::FmpzConstRef(linear_numerator),
                              flint::FmpzConstRef(denominator))) {
        return false;
    }
    if (odd_discriminant) {
        flint::fmpz_sub(flint::FmpzRef(constant),
                        flint::FmpzConstRef(constant),
                        flint::FmpzConstRef(linear));
    }

    flint::Fmpz power_basis_constant;
    if (!divide_exact_checked(flint::FmpzRef(power_basis_constant),
                              flint::FmpzConstRef(constant),
                              flint::FmpzConstRef(two))) {
        return false;
    }

    flint::Fmpz published_constant;
    if (odd_discriminant) {
        flint::fmpz_mul_ui(flint::FmpzRef(published_constant),
                           flint::FmpzConstRef(power_basis_constant), 2);
        flint::fmpz_add(flint::FmpzRef(published_constant),
                        flint::FmpzConstRef(published_constant),
                        flint::FmpzConstRef(linear));
    } else {
        flint::fmpz_set(flint::FmpzRef(published_constant),
                        flint::FmpzConstRef(power_basis_constant));
    }

    flint::fmpz_mul(flint::FmpzRef(constant_product),
                    flint::FmpzConstRef(published_constant),
                    flint::FmpzConstRef(published_constant));
    flint::fmpz_mul(flint::FmpzRef(linear_product),
                    flint::FmpzConstRef(linear),
                    flint::FmpzConstRef(linear));
    flint::fmpz_submul(flint::FmpzRef(constant_product),
                       flint::FmpzConstRef(radicand),
                       flint::FmpzConstRef(linear_product));
    const slong expected_norm_numerator = odd_discriminant ? 4 : 1;
    if (!flint::fmpz_equal_si(flint::FmpzConstRef(constant_product),
                              expected_norm_numerator) &&
        !flint::fmpz_equal_si(flint::FmpzConstRef(constant_product),
                              -expected_norm_numerator)) {
        return false;
    }
    return set_quadratic_coefficients(
            out, flint::FmpzConstRef(published_constant),
            flint::FmpzConstRef(linear), odd_discriminant);
}

bool unit_log_matrix(flint::ArbMatRef out,
                     EmbeddingContext& embeddings,
                     ElementSpan units,
                     LogEmbeddingMode mode,
                     slong precision) noexcept {
    if (precision <= 0 ||
        (mode != LogEmbeddingMode::plain &&
         mode != LogEmbeddingMode::product) ||
        !valid_unit_span(embeddings, units)) {
        return false;
    }

    Signature sig;
    if (!signature_from_embeddings(sig, embeddings)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    if (flint::arb_mat_nrows_value(out) !=
                static_cast<slong>(units.size()) ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    flint::ArbMat tmp(static_cast<slong>(units.size()), places);
    flint::ArbVec row(places);
    for (std::size_t i = 0; i < units.size(); ++i) {
        if (!logarithmic_embedding(flint::ArbVecRef(row), embeddings,
                                   units[i], mode, precision)) {
            return false;
        }
        for (slong j = 0; j < places; ++j) {
            arb_set(arb_mat_entry(tmp.raw(), static_cast<slong>(i), j),
                    row.data() + j);
        }
    }

    arb_mat_set(out.raw(), tmp.raw());
    return true;
}

bool unit_regulator(flint::ArbRef out,
                    EmbeddingContext& embeddings,
                    ElementSpan units,
                    slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }

    Signature sig;
    if (!signature_from_embeddings(sig, embeddings)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    const slong rank = places - 1;
    if (static_cast<slong>(units.size()) != rank ||
        !valid_unit_span(embeddings, units)) {
        return false;
    }

    if (rank == 0) {
        arb_set_ui(out.raw(), 1);
        return true;
    }

    flint::ArbMat logs(rank, places);
    flint::ArbMat minor(rank, rank);
    flint::Arb determinant;
    if (!unit_log_matrix(flint::ArbMatRef(logs), embeddings, units,
                         LogEmbeddingMode::product, precision)) {
        return false;
    }

    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < rank; ++j) {
            arb_set(arb_mat_entry(minor.raw(), i, j),
                    arb_mat_entry(logs.raw(), i, j));
        }
    }

    arb_mat_det(determinant.raw(), minor.raw(), precision);
    arb_abs(determinant.raw(), determinant.raw());
    if (arb_is_finite(determinant.raw()) == 0 ||
        arb_contains_zero(determinant.raw()) != 0) {
        return false;
    }

    arb_set(out.raw(), determinant.raw());
    return true;
}

bool units_independent(bool& independent,
                       EmbeddingContext& embeddings,
                       ElementSpan units,
                       slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }
    if (units.empty()) {
        independent = true;
        return true;
    }
    if (!valid_unit_span(embeddings, units)) {
        return false;
    }

    Signature sig;
    if (!signature_from_embeddings(sig, embeddings)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    const slong rank = places - 1;
    if (static_cast<slong>(units.size()) > rank) {
        independent = false;
        return true;
    }

    const slong len = static_cast<slong>(units.size());
    flint::ArbMat logs(len, places);
    flint::ArbMat transpose(places, len);
    flint::ArbMat gram(len, len);
    flint::Arb determinant;

    if (!unit_log_matrix(flint::ArbMatRef(logs), embeddings, units,
                         LogEmbeddingMode::product, precision)) {
        return false;
    }
    arb_mat_transpose(transpose.raw(), logs.raw());
    arb_mat_mul(gram.raw(), logs.raw(), transpose.raw(), precision);
    arb_mat_det(determinant.raw(), gram.raw(), precision);

    if (flint::arb_is_positive(determinant)) {
        independent = true;
        return true;
    }
    if (flint::arb_is_zero(determinant)) {
        independent = false;
        return true;
    }
    return false;
}

namespace detail {
namespace {

// Working precision for the signature bounds.  The bounds are lower bounds
// only, so a lower endpoint computed at 64 bits loses at most a relative
// 2^-60 or so against a tighter evaluation and stays rigorous at any
// precision; capping it keeps the Satz 3 special-function evaluations cheap
// for callers that ask for several hundred bits.
constexpr slong kRegulatorBoundPrecision = 64;

bool valid_signature_bound_input(slong r1,
                                 slong r2,
                                 flint::FmpzConstRef w,
                                 slong precision) noexcept {
    return r1 >= 0 && r2 >= 0 && r1 + r2 >= 1 && precision > 0 &&
           fmpz_sgn(w.raw()) > 0;
}

bool finish_bound(flint::ArbRef out, const flint::Arb& value) noexcept {
    if (arb_is_finite(value.raw()) == 0 || arb_is_positive(value.raw()) == 0) {
        return false;
    }
    arb_set(out.raw(), value.raw());
    return true;
}

// Friedman 1989, Table 6 "Lower bounds for the regulator" (p. 621), valid
// for all fields of the given signature.  Entries are in units of 10^-4 and
// are the rounded-down values printed in the table.  The table also lists
// (r1, r2) = (1, 0) and (0, 1) with value 1, the conventional regulator of Q
// and of imaginary quadratic fields; those have unit rank zero and are not
// used.  The (0, 3) entry is 0.27 "with the three exceptions
// D_K = -10051, -10571 and -12167" (note b), whose regulators are 0.2052,
// 0.2132 and 0.2372 (Theorem B, p. 599), so it is capped at 0.2052 here.
struct SignatureRegulatorBound {
    slong r1;
    slong r2;
    ulong ten_thousandths;
};

constexpr SignatureRegulatorBound kFriedmanTable6[] = {
        {0, 2, 3000},  {0, 3, 2052},  {0, 4, 2960},  {0, 5, 2500},
        {0, 6, 3000},  {0, 7, 4800},  {0, 8, 8200},  {0, 9, 4700},
        {1, 1, 2800},  {1, 2, 2600},  {1, 3, 3700},  {1, 4, 2600},
        {1, 5, 3900},  {1, 6, 6200},  {1, 7, 10500}, {1, 8, 6100},
        {2, 0, 4800},  {2, 1, 3600},  {2, 2, 4700},  {2, 3, 3600},
        {2, 4, 5200},  {2, 5, 3000},  {3, 0, 5200},  {3, 1, 6200},
        {3, 2, 8600},  {3, 3, 6500},  {3, 4, 3900},  {4, 0, 8200},
        {4, 1, 12300}, {4, 2, 3100},  {5, 0, 16000}, {5, 1, 4000},
        {6, 0, 32200}, {7, 0, 10800},
};

// Rational gamma values at which Zimmert Satz 3 is evaluated.  Every
// gamma > 0 gives a valid bound; this fixed set is within about 6% of the
// optimum over gamma for every signature of degree at most 20 (GP scan with
// step 0.02 on [0.05, 5]).
struct RationalGamma {
    slong numerator;
    slong denominator;
};

constexpr RationalGamma kZimmertSatz3Gammas[] = {
        {1, 10}, {3, 20}, {1, 5}, {3, 10}, {2, 5}, {1, 2}, {3, 5}, {4, 5},
        {1, 1},  {6, 5},  {7, 5}, {8, 5},  {2, 1}, {5, 2}, {3, 1},
};

void raise_lower_bound(flint::Arf& best,
                       bool& have_best,
                       const flint::Arb& candidate,
                       slong precision) noexcept {
    flint::Arf lower;
    arb_get_lbound_arf(lower.raw(), candidate.raw(), precision);
    if (arf_is_finite(lower.raw()) == 0 || arf_sgn(lower.raw()) <= 0) {
        return;
    }
    if (!have_best || arf_cmp(lower.raw(), best.raw()) > 0) {
        arf_set(best.raw(), lower.raw());
        have_best = true;
    }
}

}  // namespace

bool regulator_lower_bound_friedman_minimum(flint::ArbRef out,
                                            slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }
    // Theorem B: the smallest regulator of all number fields is that of the
    // totally complex sextic field of discriminant -10051, 0.2052 to four
    // places (GP 2.17: 0.20521646...), so R >= 0.2052 for every field.
    flint::Arb value;
    arb_set_ui(value.raw(), 2052);
    arb_div_ui(value.raw(), value.raw(), 10000, precision);
    return finish_bound(out, value);
}

bool regulator_lower_bound_zimmert_corollary(flint::ArbRef out,
                                             slong r1,
                                             slong r2,
                                             flint::FmpzConstRef w,
                                             slong precision) noexcept {
    if (!valid_signature_bound_input(r1, r2, w, precision)) {
        return false;
    }
    flint::Arb exponent;
    flint::Arb term;
    arb_set_si(exponent.raw(), r1);
    arb_mul_ui(exponent.raw(), exponent.raw(), 46, precision);
    arb_set_si(term.raw(), r2);
    arb_mul_ui(term.raw(), term.raw(), 10, precision);
    arb_add(exponent.raw(), exponent.raw(), term.raw(), precision);
    arb_div_ui(exponent.raw(), exponent.raw(), 100, precision);

    flint::Arb value;
    arb_exp(value.raw(), exponent.raw(), precision);
    arb_mul_fmpz(value.raw(), value.raw(), w.raw(), precision);
    arb_div_ui(value.raw(), value.raw(), 50, precision);
    return finish_bound(out, value);
}

bool regulator_lower_bound_friedman_corollary(flint::ArbRef out,
                                              slong r1,
                                              slong r2,
                                              flint::FmpzConstRef w,
                                              slong precision) noexcept {
    if (!valid_signature_bound_input(r1, r2, w, precision)) {
        return false;
    }
    flint::Arb exponent;
    flint::Arb term;
    arb_set_si(exponent.raw(), r1);
    arb_mul_ui(exponent.raw(), exponent.raw(), 241 + 497, precision);
    arb_set_si(term.raw(), r2);
    arb_mul_ui(term.raw(), term.raw(), 2 * 241, precision);
    arb_add(exponent.raw(), exponent.raw(), term.raw(), precision);
    arb_div_ui(exponent.raw(), exponent.raw(), 1000, precision);

    flint::Arb value;
    arb_exp(value.raw(), exponent.raw(), precision);
    arb_mul_fmpz(value.raw(), value.raw(), w.raw(), precision);
    arb_mul_ui(value.raw(), value.raw(), 31, precision);
    arb_div_ui(value.raw(), value.raw(), 10000, precision);
    return finish_bound(out, value);
}

bool regulator_lower_bound_zimmert_satz3(flint::ArbRef out,
                                         slong r1,
                                         slong r2,
                                         flint::FmpzConstRef w,
                                         slong gamma_numerator,
                                         slong gamma_denominator,
                                         slong precision) noexcept {
    if (!valid_signature_bound_input(r1, r2, w, precision) ||
        gamma_numerator <= 0 || gamma_denominator <= 0) {
        return false;
    }

    // R/w >= (1+g)(1+2g)/2 * Gamma(1+g)^(r1+r2) * Gamma(3/2+g)^r2
    //        * 2^(-r1-r2) * pi^(-r2/2)
    //        * exp{(-1-g)[(r1+r2) psi((1+g)/2) + r2 psi(1+g/2)
    //                     + 2/g + 1/(1+g)]},
    // evaluated as the exponential of its logarithm.
    flint::Arb g;
    flint::Arb one_plus_g;
    flint::Arb a;
    flint::Arb b;
    flint::Arb x;
    flint::Arb y;
    flint::Arb log_bound;
    flint::Arb bracket;

    arb_set_si(g.raw(), gamma_numerator);
    arb_div_si(g.raw(), g.raw(), gamma_denominator, precision);
    arb_add_ui(one_plus_g.raw(), g.raw(), 1, precision);
    arb_set_si(a.raw(), r1 + r2);
    arb_set_si(b.raw(), r2);

    // log((1+g)(1+2g)/2)
    arb_mul_2exp_si(x.raw(), g.raw(), 1);
    arb_add_ui(x.raw(), x.raw(), 1, precision);
    arb_mul(x.raw(), x.raw(), one_plus_g.raw(), precision);
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_log(log_bound.raw(), x.raw(), precision);

    // (r1+r2) log Gamma(1+g)
    arb_lgamma(x.raw(), one_plus_g.raw(), precision);
    arb_addmul(log_bound.raw(), a.raw(), x.raw(), precision);

    // r2 log Gamma(3/2+g)
    arb_set_ui(x.raw(), 3);
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_add(x.raw(), x.raw(), g.raw(), precision);
    arb_lgamma(y.raw(), x.raw(), precision);
    arb_addmul(log_bound.raw(), b.raw(), y.raw(), precision);

    // -(r1+r2) log 2
    arb_const_log2(x.raw(), precision);
    arb_submul(log_bound.raw(), a.raw(), x.raw(), precision);

    // -(r2/2) log pi
    arb_const_pi(x.raw(), precision);
    arb_log(x.raw(), x.raw(), precision);
    arb_mul_2exp_si(x.raw(), x.raw(), -1);
    arb_submul(log_bound.raw(), b.raw(), x.raw(), precision);

    // (r1+r2) psi((1+g)/2)
    arb_mul_2exp_si(x.raw(), one_plus_g.raw(), -1);
    arb_digamma(y.raw(), x.raw(), precision);
    arb_mul(bracket.raw(), a.raw(), y.raw(), precision);

    // + r2 psi(1+g/2)
    arb_mul_2exp_si(x.raw(), g.raw(), -1);
    arb_add_ui(x.raw(), x.raw(), 1, precision);
    arb_digamma(y.raw(), x.raw(), precision);
    arb_addmul(bracket.raw(), b.raw(), y.raw(), precision);

    // + 2/g + 1/(1+g)
    arb_set_ui(x.raw(), 2);
    arb_div(x.raw(), x.raw(), g.raw(), precision);
    arb_add(bracket.raw(), bracket.raw(), x.raw(), precision);
    arb_inv(x.raw(), one_plus_g.raw(), precision);
    arb_add(bracket.raw(), bracket.raw(), x.raw(), precision);

    // -(1+g) [ ... ]
    arb_submul(log_bound.raw(), one_plus_g.raw(), bracket.raw(), precision);

    flint::Arb value;
    arb_exp(value.raw(), log_bound.raw(), precision);
    arb_mul_fmpz(value.raw(), value.raw(), w.raw(), precision);
    return finish_bound(out, value);
}

bool regulator_lower_bound_friedman_table6(flint::ArbRef out,
                                           slong r1,
                                           slong r2,
                                           slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }
    for (const SignatureRegulatorBound& entry : kFriedmanTable6) {
        if (entry.r1 == r1 && entry.r2 == r2) {
            flint::Arb value;
            arb_set_ui(value.raw(), entry.ten_thousandths);
            arb_div_ui(value.raw(), value.raw(), 10000, precision);
            return finish_bound(out, value);
        }
    }
    return false;
}

bool regulator_lower_bound_from_signature(flint::ArbRef out,
                                          slong r1,
                                          slong r2,
                                          flint::FmpzConstRef w,
                                          slong precision) noexcept {
    if (!valid_signature_bound_input(r1, r2, w, precision)) {
        return false;
    }
    const slong prec = precision < kRegulatorBoundPrecision
                               ? precision
                               : kRegulatorBoundPrecision;

    flint::Arf best;
    bool have_best = false;
    flint::Arb term;

    if (regulator_lower_bound_friedman_minimum(flint::ArbRef(term), prec)) {
        raise_lower_bound(best, have_best, term, prec);
    }
    if (regulator_lower_bound_zimmert_corollary(flint::ArbRef(term), r1, r2,
                                                w, prec)) {
        raise_lower_bound(best, have_best, term, prec);
    }
    if (regulator_lower_bound_friedman_corollary(flint::ArbRef(term), r1, r2,
                                                 w, prec)) {
        raise_lower_bound(best, have_best, term, prec);
    }
    for (const RationalGamma& gamma : kZimmertSatz3Gammas) {
        if (regulator_lower_bound_zimmert_satz3(
                    flint::ArbRef(term), r1, r2, w, gamma.numerator,
                    gamma.denominator, prec)) {
            raise_lower_bound(best, have_best, term, prec);
        }
    }
    // Table 6 covers only fields of positive unit rank.
    if (r1 + r2 >= 2 &&
        regulator_lower_bound_friedman_table6(flint::ArbRef(term), r1, r2,
                                              prec)) {
        raise_lower_bound(best, have_best, term, prec);
    }

    if (!have_best) {
        return false;
    }
    arb_set_arf(out.raw(), best.raw());
    return true;
}

}  // namespace detail

}  // namespace silex
