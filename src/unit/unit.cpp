#include <silex/unit.hpp>

#include "../element/element_internal.hpp"

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

    flint::Fmpz roots_order;
    if (!root_of_unity_order(flint::FmpzRef(roots_order), field)) {
        fmpz_set_ui(roots_order.raw(), 2);
    }

    flint::Arb t;
    flint::Arb u;
    flint::Arb zimmert;
    flint::Arb floor;
    flint::Arb candidate;
    flint::Arf lower;

    arb_set_si(t.raw(), sig.r1());
    arb_mul_ui(t.raw(), t.raw(), 46, precision);
    arb_set_si(u.raw(), sig.r2());
    arb_add(t.raw(), t.raw(), u.raw(), precision);
    arb_div_ui(t.raw(), t.raw(), 100, precision);
    arb_exp(t.raw(), t.raw(), precision);

    arb_set_fmpz(zimmert.raw(), roots_order.raw());
    arb_mul(zimmert.raw(), zimmert.raw(), t.raw(), precision);
    arb_mul_ui(zimmert.raw(), zimmert.raw(), 4, precision);
    arb_div_ui(zimmert.raw(), zimmert.raw(), 100, precision);

    arb_set_ui(floor.raw(), 54);
    arb_div_ui(floor.raw(), floor.raw(), 1000, precision);
    arb_max(candidate.raw(), zimmert.raw(), floor.raw(), precision);
    if (arb_is_finite(candidate.raw()) == 0) {
        return false;
    }

    arb_get_lbound_arf(lower.raw(), candidate.raw(), precision);
    if (arf_is_finite(lower.raw()) == 0 || arf_sgn(lower.raw()) <= 0) {
        return false;
    }

    arb_set_arf(out.raw(), lower.raw());
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

}  // namespace silex
