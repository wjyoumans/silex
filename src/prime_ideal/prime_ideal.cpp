#include <silex/prime_ideal.hpp>

#include "prime_ideal_internal.hpp"

#include <silex/diagnostics.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/lat.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_poly.h>

#include <algorithm>
#include <new>
#include <utility>
#include <vector>

namespace silex {
namespace {

#define SILEX_PRIME_IDEAL_PROFILE_EVENT(diagnostics, label) \
    SILEX_PROFILE_EVENT((diagnostics), DiagnosticsModule::prime_ideal, (label))

bool fmpq_poly_is_monic_integral(const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coeff.raw())) == 0) {
            return false;
        }
    }

    fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, degree);
    return fmpz_is_one(fmpq_numref(coeff.raw())) != 0;
}

void fmpq_poly_get_fmpz_mod_poly(flint::FmpzModPoly& out,
                                 const fmpq_poly_t polynomial,
                                 const flint::FmpzModCtx& ctx) noexcept {
    flint::Fmpz coeff;
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpz(coeff.raw(), polynomial, i);
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff.raw(), ctx.raw());
    }
}

bool fmpq_poly_to_integral_fmpz_poly(flint::FmpzPoly& out,
                                     const fmpq_poly_t input) noexcept {
    if (!fmpq_poly_is_monic_integral(input)) {
        return false;
    }

    fmpz_poly_zero(out.raw());
    flint::Fmpq coeff;
    const slong degree = fmpq_poly_degree(input);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), input, i);
        fmpz_poly_set_coeff_fmpz(out.raw(), i, fmpq_numref(coeff.raw()));
    }
    return true;
}

void quadratic_integral_generator_minimal_polynomial(
        flint::FmpzPoly& out, flint::FmpzConstRef radicand) noexcept {
    flint::Fmpz constant;
    fmpz_poly_zero(out.raw());
    fmpz_poly_set_coeff_si(out.raw(), 2, 1);
    if (flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        fmpz_poly_set_coeff_si(out.raw(), 1, -1);
        fmpz_sub_ui(constant.raw(), radicand.raw(), 1);
        fmpz_neg(constant.raw(), constant.raw());
        fmpz_divexact_ui(constant.raw(), constant.raw(), 4);
    } else {
        fmpz_neg(constant.raw(), radicand.raw());
    }
    fmpz_poly_set_coeff_fmpz(out.raw(), 0, constant.raw());
}

bool set_quadratic_integral_generator(
        Element& out,
        const NumberField& field,
        flint::FmpzConstRef radicand) noexcept {
    Element theta(field);
    Element shifted(field);
    if (!out.is_defined() || !theta.is_defined() || !shifted.is_defined() ||
        !theta.gen()) {
        return false;
    }
    if (flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        return shifted.add_si(theta, 1) && out.scalar_div_si(shifted, 2);
    }
    return out.set(theta);
}

bool fmpq_poly_reduce_modulus(flint::FmpzModPoly& out,
                              const fmpq_poly_t polynomial,
                              const flint::FmpzModCtx& ctx,
                              const fmpz_t modulus) noexcept {
    flint::Fmpq coeff;
    flint::Fmpz coeff_mod;
    flint::Fmpz inverse_denominator;

    fmpz_mod_poly_zero(out.raw(), ctx.raw());
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_invmod(inverse_denominator.raw(),
                        fmpq_denref(coeff.raw()), modulus) == 0) {
            return false;
        }
        fmpz_mod_set_fmpz(coeff_mod.raw(), fmpq_numref(coeff.raw()),
                          ctx.raw());
        fmpz_mod_mul(coeff_mod.raw(), coeff_mod.raw(),
                     inverse_denominator.raw(), ctx.raw());
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff_mod.raw(),
                                     ctx.raw());
    }
    return true;
}

bool hensel_lift_factor_to_precision(
        flint::FmpzModPoly& out,
        const fmpq_poly_t field_polynomial,
        const flint::FmpzModPoly& field_mod_prime,
        const flint::FmpzModPoly& factor_mod_prime,
        const flint::FmpzModCtx& prime_ctx,
        const fmpz_t prime,
        slong target_exponent,
        const fmpz_t target_modulus,
        const flint::FmpzModCtx& target_ctx) noexcept {
    const slong degree =
            fmpz_mod_poly_degree(field_mod_prime.raw(), prime_ctx.raw());
    if (fmpz_mod_poly_degree(factor_mod_prime.raw(), prime_ctx.raw()) ==
        degree) {
        return fmpq_poly_reduce_modulus(out, field_polynomial, target_ctx,
                                        target_modulus);
    }
    if (target_exponent <= 0) {
        return false;
    }

    flint::FmpzPoly field_integral;
    if (!fmpq_poly_to_integral_fmpz_poly(field_integral, field_polynomial)) {
        return false;
    }

    flint::FmpzModPoly quotient(prime_ctx);
    flint::FmpzModPoly remainder(prime_ctx);
    flint::FmpzModPoly gcd(prime_ctx);
    flint::FmpzModPoly a_mod(prime_ctx);
    flint::FmpzModPoly b_mod(prime_ctx);
    if (!quotient.is_initialized() || !remainder.is_initialized() ||
        !gcd.is_initialized() || !a_mod.is_initialized() ||
        !b_mod.is_initialized()) {
        return false;
    }

    fmpz_mod_poly_divrem(quotient.raw(), remainder.raw(),
                         field_mod_prime.raw(), factor_mod_prime.raw(),
                         prime_ctx.raw());
    if (fmpz_mod_poly_is_zero(remainder.raw(), prime_ctx.raw()) == 0) {
        return false;
    }
    fmpz_mod_poly_xgcd(gcd.raw(), a_mod.raw(), b_mod.raw(),
                       factor_mod_prime.raw(), quotient.raw(),
                       prime_ctx.raw());
    if (fmpz_mod_poly_is_one(gcd.raw(), prime_ctx.raw()) == 0) {
        return false;
    }

    flint::FmpzPoly g;
    flint::FmpzPoly h;
    flint::FmpzPoly a;
    flint::FmpzPoly b;
    fmpz_mod_poly_get_fmpz_poly(g.raw(), factor_mod_prime.raw(),
                                prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(h.raw(), quotient.raw(), prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(a.raw(), a_mod.raw(), prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(b.raw(), b_mod.raw(), prime_ctx.raw());

    flint::Fmpz current_modulus;
    fmpz_set(current_modulus.raw(), prime);
    slong current_exponent = 1;
    while (current_exponent < target_exponent) {
        const slong step = std::min(current_exponent,
                                    target_exponent - current_exponent);
        flint::Fmpz next_factor;
        fmpz_pow_ui(next_factor.raw(), prime, static_cast<ulong>(step));

        flint::FmpzPoly g_out;
        flint::FmpzPoly h_out;
        flint::FmpzPoly a_out;
        flint::FmpzPoly b_out;
        fmpz_poly_hensel_lift(g_out.raw(), h_out.raw(), a_out.raw(),
                              b_out.raw(), field_integral.raw(), g.raw(),
                              h.raw(), a.raw(), b.raw(),
                              current_modulus.raw(), next_factor.raw());
        g = std::move(g_out);
        h = std::move(h_out);
        a = std::move(a_out);
        b = std::move(b_out);
        fmpz_mul(current_modulus.raw(), current_modulus.raw(),
                 next_factor.raw());
        current_exponent += step;
    }

    fmpz_poly_scalar_mod_fmpz(g.raw(), g.raw(), target_modulus);
    fmpz_mod_poly_set_fmpz_poly(out.raw(), g.raw(), target_ctx.raw());
    return true;
}

slong fmpz_valuation_capped(const fmpz_t input,
                            const fmpz_t prime,
                            slong cap) noexcept {
    if (fmpz_is_zero(input) != 0) {
        return cap;
    }

    slong valuation = 0;
    flint::Fmpz remaining;
    fmpz_set(remaining.raw(), input);
    while (valuation < cap &&
           fmpz_divisible(remaining.raw(), prime) != 0) {
        fmpz_divexact(remaining.raw(), remaining.raw(), prime);
        ++valuation;
    }
    return valuation;
}

bool fmpz_mod_poly_content_valuation_capped(
        slong& out,
        const fmpz_mod_poly_t polynomial,
        const flint::FmpzModCtx& ctx,
        const fmpz_t prime,
        slong cap) noexcept {
    if (cap < 0) {
        return false;
    }

    const slong degree = fmpz_mod_poly_degree(polynomial, ctx.raw());
    if (degree < 0) {
        out = cap;
        return true;
    }

    flint::Fmpz coeff;
    slong valuation = cap;
    for (slong i = 0; i <= degree; ++i) {
        fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), polynomial, i, ctx.raw());
        valuation = std::min(
                valuation,
                fmpz_valuation_capped(coeff.raw(), prime, cap));
        if (valuation == 0) {
            break;
        }
    }

    out = valuation;
    return true;
}

bool factorization_is_squarefree(
        const flint::FmpzModPolyFactor& factorization) noexcept {
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (factorization.raw()->exp[i] != 1) {
            return false;
        }
    }
    return true;
}

void fmpz_mod_poly_get_fmpz_poly(flint::FmpzPoly& out,
                                 const fmpz_mod_poly_t polynomial,
                                 const flint::FmpzModCtx& ctx) noexcept {
    fmpz_mod_poly_get_fmpz_poly(out.raw(), polynomial, ctx.raw());
}

bool element_set_mod_poly(Element& out,
                          const fmpz_mod_poly_t polynomial,
                          const flint::FmpzModCtx& ctx) noexcept {
    flint::FmpqPoly rational;
    flint::Fmpz coeff;
    const slong degree = fmpz_mod_poly_degree(polynomial, ctx.raw());
    for (slong i = 0; i <= degree; ++i) {
        fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), polynomial, i, ctx.raw());
        fmpq_poly_set_coeff_fmpz(rational.raw(), i, coeff.raw());
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(rational));
}

bool order_element_set_fmpz(OrderElement& out,
                            const Order& order,
                            flint::FmpzConstRef value) noexcept {
    if (order.parent() == nullptr) {
        return false;
    }

    Element element(*order.parent());
    flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, value.raw());
    return element.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial)) &&
           out.set_element(element);
}

bool fmpq_mat_entries_are_integral(const fmpq_mat_t matrix) noexcept {
    for (slong i = 0; i < fmpq_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < fmpq_mat_ncols(matrix); ++j) {
            if (fmpz_is_one(fmpq_mat_entry_den(matrix, i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

void integral_coords_to_poly(flint::FmpzPoly& out,
                             const fmpq_mat_t coordinates) noexcept {
    fmpz_poly_zero(out.raw());
    const slong cols = fmpq_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        fmpz_poly_set_coeff_fmpz(out.raw(), j,
                                 fmpq_mat_entry_num(coordinates, 0, j));
    }
}

void integral_coords_to_poly(flint::FmpzPoly& out,
                             const fmpz_mat_t coordinates) noexcept {
    fmpz_poly_zero(out.raw());
    const slong cols = fmpz_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        fmpz_poly_set_coeff_fmpz(out.raw(), j,
                                 fmpz_mat_entry(coordinates, 0, j));
    }
}

bool reduce_poly_mod_residue(flint::FmpzPoly& out,
                             const flint::FmpzPoly& input,
                             const flint::FmpzPoly& residue_polynomial,
                             flint::FmpzConstRef p) noexcept {
    if (fmpz_poly_degree(residue_polynomial.raw()) <= 0) {
        return false;
    }

    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPoly modulus(ctx);
    flint::FmpzModPoly remainder(ctx);
    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), input.raw(), ctx.raw());
    fmpz_mod_poly_set_fmpz_poly(modulus.raw(), residue_polynomial.raw(),
                                ctx.raw());
    if (fmpz_mod_poly_is_zero(modulus.raw(), ctx.raw()) != 0) {
        return false;
    }

    fmpz_mod_poly_rem(remainder.raw(), reduced.raw(), modulus.raw(), ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(out.raw(), remainder.raw(), ctx.raw());
    return true;
}

bool linear_root_from_residue_polynomial(
        flint::Fmpz& out,
        flint::FmpzPolyConstRef residue_polynomial,
        flint::FmpzConstRef p) noexcept {
    if (fmpz_poly_degree(residue_polynomial.raw()) != 1) {
        return false;
    }

    flint::Fmpz constant;
    flint::Fmpz leading;
    flint::Fmpz inverse_leading;
    fmpz_poly_get_coeff_fmpz(constant.raw(), residue_polynomial.raw(), 0);
    fmpz_poly_get_coeff_fmpz(leading.raw(), residue_polynomial.raw(), 1);
    fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }
    fmpz_neg(out.raw(), constant.raw());
    fmpz_mod(out.raw(), out.raw(), p.raw());
    fmpz_mul(out.raw(), out.raw(), inverse_leading.raw());
    fmpz_mod(out.raw(), out.raw(), p.raw());
    return true;
}

bool fmpz_mat_shape(flint::FmpzMatConstRef matrix,
                    slong rows,
                    slong cols) noexcept {
    return matrix.raw() != nullptr && fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

void fmpz_vec_zero(flint::FmpzVec& vector) noexcept {
    for (slong i = 0; i < vector.length(); ++i) {
        fmpz_zero(vector.data() + i);
    }
}

void fmpz_vec_set(flint::FmpzVec& out,
                  const flint::FmpzVec& in) noexcept {
    for (slong i = 0; i < out.length(); ++i) {
        fmpz_set(out.data() + i, in.data() + i);
    }
}

bool fmpz_vec_is_zero(const flint::FmpzVec& value) noexcept {
    for (slong i = 0; i < value.length(); ++i) {
        if (fmpz_is_zero(value.data() + i) == 0) {
            return false;
        }
    }
    return true;
}

bool fmpz_vec_equal(const flint::FmpzVec& left,
                    const flint::FmpzVec& right) noexcept {
    if (left.length() != right.length()) {
        return false;
    }
    for (slong i = 0; i < left.length(); ++i) {
        if (fmpz_equal(left.data() + i, right.data() + i) == 0) {
            return false;
        }
    }
    return true;
}

void fmpz_mat_set_row(flint::FmpzMat& matrix,
                      slong row,
                      const flint::FmpzVec& vector) noexcept {
    for (slong j = 0; j < vector.length(); ++j) {
        fmpz_set(fmpz_mat_entry(matrix.raw(), row, j), vector.data() + j);
    }
}

void fmpz_mat_center_row_mod_prime(flint::FmpzMat& matrix,
                                   slong row,
                                   flint::FmpzConstRef p) noexcept {
    flint::Fmpz twice;
    for (slong j = 0; j < fmpz_mat_ncols(matrix.raw()); ++j) {
        fmpz* entry = fmpz_mat_entry(matrix.raw(), row, j);
        fmpz_mod(entry, entry, p.raw());
        fmpz_mul_2exp(twice.raw(), entry, 1);
        if (fmpz_cmp(twice.raw(), p.raw()) > 0) {
            fmpz_sub(entry, entry, p.raw());
        }
    }
}

void fmpz_center_mod_prime(fmpz_t out,
                           const fmpz_t input,
                           const fmpz_t p) noexcept {
    flint::Fmpz twice;
    fmpz_mod(out, input, p);
    fmpz_mul_2exp(twice.raw(), out, 1);
    if (fmpz_cmp(twice.raw(), p) > 0) {
        fmpz_sub(out, out, p);
    }
}

bool order_one_coords(flint::FmpzVec& out, const Order& order) noexcept {
    if (out.length() != order.degree()) {
        return false;
    }

    OrderElement one(order);
    flint::FmpzMat row(1, order.degree());
    if (!one.is_defined() || !one.one() ||
        !one.get_coordinates(flint::FmpzMatRef(row))) {
        return false;
    }

    for (slong j = 0; j < order.degree(); ++j) {
        fmpz_set(out.data() + j, fmpz_mat_entry(row.raw(), 0, j));
    }
    return true;
}

bool order_element_from_coords(OrderElement& out,
                               const Order& order,
                               const flint::FmpzVec& coords) noexcept {
    if (coords.length() != order.degree()) {
        return false;
    }

    flint::FmpzMat row(1, order.degree());
    for (slong j = 0; j < order.degree(); ++j) {
        fmpz_set(fmpz_mat_entry(row.raw(), 0, j), coords.data() + j);
    }
    return out.set_coordinates(flint::FmpzMatConstRef(row));
}

void order_coords_mul(flint::FmpzVec& out,
                      const flint::FmpzVec& left,
                      const flint::FmpzVec& right,
                      const flint::FmpzMat& multiplication_table) noexcept;

struct FiniteAlgebraComponent {
    explicit FiniteAlgebraComponent(slong degree) noexcept
        : idempotent(degree), basis(0, degree) {
    }

    flint::FmpzVec idempotent;
    flint::FmpzMat basis;
    slong dimension = 0;
};

struct IndexPrimeData {
    explicit IndexPrimeData(const Order& order) noexcept : ideal(order) {}

    Ideal ideal;
    slong ramification_index = 1;
    slong residue_degree = 0;
};

bool valuation_by_ideal_power_containment(slong& out,
                                          const Ideal& prime_ideal,
                                          const OrderElement& element,
                                          slong bound,
                                          const DiagnosticsContext* diagnostics)
        noexcept {
    // reference documents prime-ideal element valuation as the largest i with
    // a in P^i.  This exact fallback uses that definition directly, avoiding
    // construction of the whole principal ideal aO when fast element kernels
    // do not apply.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_containment");
    if (!prime_ideal.has_hnf() || bound < 0) {
        return false;
    }
    if (bound == 0) {
        out = 0;
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                            "prime_ideal.valuation_containment_initial_contains");
        if (!prime_ideal.contains(element)) {
            out = 0;
            return true;
        }
    }

    slong value = 1;
    if (value < bound) {
        const Order* order = prime_ideal.parent();
        if (order == nullptr) {
            return false;
        }

        Ideal power(*order);
        Ideal next(*order);
        if (!power.is_defined() || !next.is_defined() ||
            !power.set(prime_ideal)) {
            return false;
        }

        while (value < bound) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_containment_power_multiply");
                if (!next.multiply(power, prime_ideal)) {
                    return false;
                }
            }
            bool contained = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_containment_power_contains");
                contained = next.contains(element);
            }
            if (!contained) {
                break;
            }
            ++value;
            power.swap(next);
        }
    }

    out = value;
    return true;
}

bool quotient_multiplier_coordinates(flint::FmpzVec& out,
                                     const Order& order,
                                     const fmpq_poly_t field_polynomial,
                                     const flint::FmpzPoly& residue_polynomial,
                                     flint::FmpzConstRef p) noexcept {
    // Source trace: reference `idealprimedec_kummer` stores `tau` as the
    // multiplication table by `T/u mod p` for Kummer primes `(p, u)`;
    // `ZC_nfvalrem` then multiplies coordinates by that `tau` to test and
    // divide p-adically.  `poltobasis` expresses the quotient in reference's
    // integral basis; do the same through OrderElement for a maximal order.
    const slong n = order.degree();
    if (out.length() != n || n <= 0 ||
        (!order.is_equation_order() && !order.is_maximal()) ||
        order.parent() == nullptr ||
        !fmpq_poly_is_monic_integral(field_polynomial) ||
        fmpz_poly_degree(residue_polynomial.raw()) <= 0) {
        return false;
    }

    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzModPoly field_mod_prime(ctx);
    flint::FmpzModPoly factor_mod_prime(ctx);
    flint::FmpzModPoly quotient(ctx);
    flint::FmpzModPoly remainder(ctx);
    if (ctx.raw() == nullptr || !field_mod_prime.is_initialized() ||
        !factor_mod_prime.is_initialized() || !quotient.is_initialized() ||
        !remainder.is_initialized()) {
        return false;
    }

    fmpq_poly_get_fmpz_mod_poly(field_mod_prime, field_polynomial, ctx);
    fmpz_mod_poly_set_fmpz_poly(factor_mod_prime.raw(),
                                residue_polynomial.raw(), ctx.raw());
    if (fmpz_mod_poly_degree(factor_mod_prime.raw(), ctx.raw()) <= 0) {
        return false;
    }
    fmpz_mod_poly_divrem(quotient.raw(), remainder.raw(),
                         field_mod_prime.raw(), factor_mod_prime.raw(),
                         ctx.raw());
    if (fmpz_mod_poly_is_zero(remainder.raw(), ctx.raw()) == 0 ||
        fmpz_mod_poly_degree(quotient.raw(), ctx.raw()) >= n) {
        return false;
    }

    if (order.is_equation_order()) {
        fmpz_vec_zero(out);
        flint::Fmpz coeff;
        const slong quotient_degree =
                fmpz_mod_poly_degree(quotient.raw(), ctx.raw());
        for (slong i = 0; i <= quotient_degree; ++i) {
            fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), quotient.raw(), i,
                                         ctx.raw());
            fmpz_center_mod_prime(out.data() + i, coeff.raw(), p.raw());
        }
        return true;
    }

    Element multiplier(*order.parent());
    OrderElement order_multiplier(order);
    flint::FmpzMat coordinates(1, n);
    if (!element_set_mod_poly(multiplier, quotient.raw(), ctx) ||
        !order_multiplier.set_element(multiplier) ||
        !order_multiplier.get_coordinates(flint::FmpzMatRef(coordinates))) {
        return false;
    }
    for (slong j = 0; j < n; ++j) {
        fmpz_center_mod_prime(out.data() + j,
                              fmpz_mat_entry(coordinates.raw(), 0, j),
                              p.raw());
    }
    return true;
}

bool order_multiplier_matrix_from_coordinates(
        flint::FmpzMat& out,
        const Order& order,
        const flint::FmpzVec& multiplier) noexcept {
    const slong n = order.degree();
    if (!fmpz_mat_shape(flint::FmpzMatConstRef(out), n, n) ||
        multiplier.length() != n) {
        return false;
    }

    flint::FmpzMat multiplication_table(n * n, n);
    Order table_order = order;
    if (!table_order.multiplication_table(
                flint::FmpzMatRef(multiplication_table))) {
        return false;
    }

    flint::fmpz_mat_zero(flint::FmpzMatRef(out));
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (fmpz_is_zero(multiplier.data() + j) != 0) {
                continue;
            }
            for (slong k = 0; k < n; ++k) {
                fmpz_addmul(
                        fmpz_mat_entry(out.raw(), i, k),
                        multiplier.data() + j,
                        fmpz_mat_entry(multiplication_table.raw(), i * n + j,
                                       k));
            }
        }
    }
    return true;
}

bool quadratic_quotient_multiplier_matrix(
        flint::FmpzMat& out,
        const Order& order,
        flint::FmpzConstRef radicand,
        const flint::FmpzPoly& residue_polynomial,
        flint::FmpzConstRef p) noexcept {
    if (order.degree() != 2 ||
        fmpz_poly_degree(residue_polynomial.raw()) <= 0) {
        return false;
    }

    flint::FmpzPoly minimal_polynomial;
    quadratic_integral_generator_minimal_polynomial(minimal_polynomial,
                                                     radicand);
    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzModPoly field_mod_prime(ctx);
    flint::FmpzModPoly factor_mod_prime(ctx);
    flint::FmpzModPoly quotient(ctx);
    flint::FmpzModPoly remainder(ctx);
    if (ctx.raw() == nullptr || !field_mod_prime.is_initialized() ||
            !factor_mod_prime.is_initialized() || !quotient.is_initialized() ||
            !remainder.is_initialized()) {
        return false;
    }

    fmpz_mod_poly_set_fmpz_poly(
            field_mod_prime.raw(), minimal_polynomial.raw(), ctx.raw());
    fmpz_mod_poly_set_fmpz_poly(
            factor_mod_prime.raw(), residue_polynomial.raw(), ctx.raw());
    fmpz_mod_poly_divrem(quotient.raw(),
            remainder.raw(),
            field_mod_prime.raw(),
            factor_mod_prime.raw(),
            ctx.raw());
    if (fmpz_mod_poly_is_zero(remainder.raw(), ctx.raw()) == 0 ||
            fmpz_mod_poly_degree(quotient.raw(), ctx.raw()) >= 2) {
        return false;
    }

    // The direct maximal quadratic basis is exactly [1, omega], so the
    // quotient polynomial coefficients are already order coordinates.
    flint::FmpzVec multiplier(2);
    fmpz_vec_zero(multiplier);
    flint::Fmpz coefficient;
    const slong quotient_degree =
            fmpz_mod_poly_degree(quotient.raw(), ctx.raw());
    for (slong i = 0; i <= quotient_degree; ++i) {
        fmpz_mod_poly_get_coeff_fmpz(coefficient.raw(), quotient.raw(), i,
                                     ctx.raw());
        fmpz_center_mod_prime(multiplier.data() + i, coefficient.raw(),
                              p.raw());
    }
    return order_multiplier_matrix_from_coordinates(out, order, multiplier);
}

bool quotient_multiplier_matrix(flint::FmpzMat& out,
        const Order& order,
        const fmpq_poly_t field_polynomial,
        const flint::FmpzPoly& residue_polynomial,
        flint::FmpzConstRef p) noexcept {
    // reference `idealprimedec_kummer` stores `zk_multable(nf, tau)` in the prime
    // record.  Build that same fixed linear map once and leave repeated
    // products to FLINT.
    flint::FmpzVec multiplier(order.degree());
    return quotient_multiplier_coordinates(multiplier, order,
                                            field_polynomial,
                                            residue_polynomial, p) &&
           order_multiplier_matrix_from_coordinates(out, order, multiplier);
}

bool valuation_by_coordinate_divisibility(
        slong& out,
        const Order& order,
        const OrderElement& element,
        flint::FmpzMatConstRef multiplier_matrix,
        flint::FmpzConstRef p,
        slong bound,
        const DiagnosticsContext* diagnostics) noexcept {
    // Source trace: reference `base3.c:ZC_nfvalrem` and reference `val_func_index`
    // repeat `x = (tau * x) / p` while every order coordinate of `tau * x`
    // is divisible by p.  The cached map is reference's T/u map expressed in the
    // active equation or maximal-order basis.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_coordinate_divisibility");
    if (bound < 0 ||
        (!order.is_equation_order() && !order.is_maximal()) ||
        !same_order_parent(element.parent(), &order)) {
        return false;
    }
    if (bound == 0) {
        out = 0;
        return true;
    }

    const slong n = order.degree();
    if (!fmpz_mat_shape(multiplier_matrix, n, n)) {
        return false;
    }

    flint::FmpzMat current(1, n);
    flint::FmpzMat product(1, n);
    if (!element.get_coordinates(flint::FmpzMatRef(current))) {
        return false;
    }

    slong value = 0;
    while (value < bound) {
        flint::fmpz_mat_mul(flint::FmpzMatRef(product),
                            flint::FmpzMatConstRef(current),
                            multiplier_matrix);
        for (slong j = 0; j < n; ++j) {
            if (fmpz_divisible(fmpz_mat_entry(product.raw(), 0, j),
                               p.raw()) == 0) {
                out = value;
                return true;
            }
        }
        fmpz_mat_scalar_divexact_fmpz(current.raw(), product.raw(), p.raw());
        ++value;
    }

    out = value;
    return true;
}

bool element_from_order_coords(Element& out,
                               const Order& order,
                               const flint::FmpzVec& coords) noexcept {
    OrderElement order_element(order);
    return order_element.is_defined() &&
           order_element_from_coords(order_element, order, coords) &&
           order_element.get_element(out);
}

bool residue_degree_allowed(slong residue_degree,
                            slong max_residue_degree) noexcept {
    return max_residue_degree <= 0 || residue_degree <= max_residue_degree;
}

bool maximal_order_has_nonindex_defining_polynomial_prime(
        const Order& order,
        flint::FmpzConstRef p) noexcept {
    // Source trace: reference
    // `NumFieldOrd/NfOrd/Ideal/Prime.jl:prime_decomposition` routes maximal
    // orders with a nice defining polynomial and non-index prime through
    // `prime_dec_nonindex`, factoring the field polynomial modulo p instead
    // of first searching for a global monogenic generator of the order.
    const NumberField* field = order.parent();
    if (field == nullptr || !order.is_maximal()) {
        return false;
    }

    const nf_struct* raw_field = field->raw_flint_field();
    if (raw_field == nullptr ||
        !fmpq_poly_is_monic_integral(raw_field->pol)) {
        return false;
    }

    Order equation = Order::equation_order(*field);
    flint::Fmpz index;
    if (!equation.is_defined() ||
        !order_index(flint::FmpzRef(index), equation, order)) {
        return false;
    }

    return fmpz_divisible(index.raw(), p.raw()) == 0;
}

void order_coords_mul(flint::FmpzVec& out,
                      const flint::FmpzVec& left,
                      const flint::FmpzVec& right,
                      const flint::FmpzMat& multiplication_table) noexcept {
    const slong n = out.length();
    flint::FmpzVec result(n);
    flint::Fmpz product;

    for (slong i = 0; i < n; ++i) {
        if (fmpz_is_zero(left.data() + i) != 0) {
            continue;
        }
        for (slong j = 0; j < n; ++j) {
            if (fmpz_is_zero(right.data() + j) != 0) {
                continue;
            }

            fmpz_mul(product.raw(), left.data() + i, right.data() + j);
            for (slong k = 0; k < n; ++k) {
                fmpz_addmul(
                        result.data() + k, product.raw(),
                        fmpz_mat_entry(multiplication_table.raw(), i * n + j,
                                       k));
            }
        }
    }

    fmpz_vec_set(out, result);
}

void set_identity_coords(flint::FmpzVec& out, slong index) noexcept {
    fmpz_vec_zero(out);
    fmpz_one(out.data() + index);
}

void set_small_trinary_coords(flint::FmpzVec& out, slong value) noexcept {
    slong u = value;
    for (slong i = 0; i < out.length(); ++i) {
        switch (u % 3) {
            case 0:
                fmpz_set_si(out.data() + i, -1);
                break;
            case 1:
                fmpz_zero(out.data() + i);
                break;
            default:
                fmpz_one(out.data() + i);
                break;
        }
        u /= 3;
    }
}

bool build_multiplication_table(flint::FmpzMat& out,
                                const Order& order) noexcept {
    const slong n = order.degree();
    if (!fmpz_mat_shape(flint::FmpzMatConstRef(out), n * n, n)) {
        return false;
    }

    OrderElement left(order);
    OrderElement right(order);
    OrderElement product(order);
    flint::FmpzVec left_coords(n);
    flint::FmpzVec right_coords(n);
    flint::FmpzMat row(1, n);
    flint::FmpzMat coords(1, n);
    if (!left.is_defined() || !right.is_defined() || !product.is_defined()) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        set_identity_coords(left_coords, i);
        if (!order_element_from_coords(left, order, left_coords)) {
            return false;
        }

        for (slong j = 0; j < n; ++j) {
            set_identity_coords(right_coords, j);
            if (!order_element_from_coords(right, order, right_coords) ||
                !product.multiply(left, right) ||
                !product.get_coordinates(flint::FmpzMatRef(coords))) {
                return false;
            }
            for (slong k = 0; k < n; ++k) {
                fmpz_set(fmpz_mat_entry(out.raw(), i * n + j, k),
                         fmpz_mat_entry(coords.raw(), 0, k));
            }
        }
    }

    return true;
}

void set_identity_coords_mod(flint::FmpzVec& out,
                             slong index,
                             const flint::FmpzModCtx& ctx) noexcept {
    fmpz_vec_zero(out);
    fmpz_mod_set_ui(out.data() + index, 1, ctx.raw());
}

void fmpz_vec_set_mod(flint::FmpzVec& out,
                      const flint::FmpzVec& input,
                      const flint::FmpzModCtx& ctx) noexcept {
    for (slong i = 0; i < out.length(); ++i) {
        fmpz_mod_set_fmpz(out.data() + i, input.data() + i, ctx.raw());
    }
}

void fmpz_vec_scalar_addmul_mod(flint::FmpzVec& out,
                                const flint::FmpzVec& input,
                                const fmpz_t scalar,
                                const flint::FmpzModCtx& ctx) noexcept {
    flint::Fmpz term;
    for (slong i = 0; i < out.length(); ++i) {
        fmpz_mod_mul(term.raw(), input.data() + i, scalar, ctx.raw());
        fmpz_add(out.data() + i, out.data() + i, term.raw());
        fmpz_mod_set_fmpz(out.data() + i, out.data() + i, ctx.raw());
    }
}

void order_coords_mul_mod(flint::FmpzVec& out,
                          const flint::FmpzVec& left,
                          const flint::FmpzVec& right,
                          const flint::FmpzMat& multiplication_table,
                          const flint::FmpzModCtx& ctx) noexcept {
    const slong n = out.length();
    flint::FmpzVec result(n);
    flint::Fmpz ab;
    flint::Fmpz coeff;
    flint::Fmpz term;

    for (slong i = 0; i < n; ++i) {
        if (fmpz_is_zero(left.data() + i) != 0) {
            continue;
        }
        for (slong j = 0; j < n; ++j) {
            if (fmpz_is_zero(right.data() + j) != 0) {
                continue;
            }

            fmpz_mod_mul(ab.raw(), left.data() + i, right.data() + j,
                         ctx.raw());
            for (slong k = 0; k < n; ++k) {
                fmpz_mod_set_fmpz(
                        coeff.raw(),
                        fmpz_mat_entry(multiplication_table.raw(), i * n + j,
                                       k),
                        ctx.raw());
                fmpz_mod_mul(term.raw(), ab.raw(), coeff.raw(), ctx.raw());
                fmpz_add(result.data() + k, result.data() + k, term.raw());
                fmpz_mod_set_fmpz(result.data() + k, result.data() + k,
                                  ctx.raw());
            }
        }
    }

    fmpz_vec_set(out, result);
}

void order_coords_pow_mod(flint::FmpzVec& out,
                          const flint::FmpzVec& base,
                          flint::FmpzConstRef exponent,
                          const flint::FmpzVec& one,
                          const flint::FmpzMat& multiplication_table,
                          const flint::FmpzModCtx& ctx) noexcept {
    const slong n = out.length();
    flint::FmpzVec acc(n);
    flint::FmpzVec power(n);
    flint::Fmpz e;

    fmpz_set(e.raw(), exponent.raw());
    fmpz_vec_set(acc, one);
    fmpz_vec_set(power, base);

    while (fmpz_is_zero(e.raw()) == 0) {
        if (fmpz_is_odd(e.raw()) != 0) {
            order_coords_mul_mod(acc, acc, power, multiplication_table, ctx);
        }
        fmpz_fdiv_q_2exp(e.raw(), e.raw(), 1);
        if (fmpz_is_zero(e.raw()) == 0) {
            order_coords_mul_mod(power, power, power, multiplication_table,
                                 ctx);
        }
    }

    fmpz_vec_set(out, acc);
}

void algebra_poly_eval_mod(flint::FmpzVec& out,
                           const fmpz_mod_poly_t polynomial,
                           const flint::FmpzVec& element,
                           const flint::FmpzVec& one,
                           const flint::FmpzMat& multiplication_table,
                           const flint::FmpzModCtx& ctx) noexcept {
    const slong n = out.length();
    flint::FmpzVec result(n);
    flint::FmpzVec product(n);
    flint::Fmpz coeff;
    const slong degree = fmpz_mod_poly_degree(polynomial, ctx.raw());

    for (slong i = degree; i >= 0; --i) {
        if (i != degree) {
            order_coords_mul_mod(product, result, element,
                                 multiplication_table, ctx);
            fmpz_vec_set(result, product);
        }
        fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), polynomial, i, ctx.raw());
        fmpz_vec_scalar_addmul_mod(result, one, coeff.raw(), ctx);
    }

    for (slong i = 0; i < n; ++i) {
        fmpz_mod_set_fmpz(out.data() + i, result.data() + i, ctx.raw());
    }
}

bool fixed_frobenius_kernel(flint::FmpzModMat& out,
                            slong& nullity,
                            const flint::FmpzMat& multiplication_table,
                            const flint::FmpzVec& one,
                            flint::FmpzConstRef p,
                            const flint::FmpzModCtx& ctx) noexcept {
    const slong n = one.length();
    flint::FmpzModMat frobenius(n, n, ctx);
    flint::FmpzVec base(n);
    flint::FmpzVec image(n);

    if (out.context() != ctx.raw()) {
        return false;
    }

    for (slong j = 0; j < n; ++j) {
        set_identity_coords_mod(base, j, ctx);
        order_coords_pow_mod(image, base, p, one, multiplication_table, ctx);
        fmpz_sub(image.data() + j, image.data() + j, base.data() + j);
        fmpz_mod_set_fmpz(image.data() + j, image.data() + j, ctx.raw());

        for (slong i = 0; i < n; ++i) {
            fmpz_set(fmpz_mod_mat_entry(frobenius.raw(), i, j),
                     image.data() + i);
        }
    }

    nullity = fmpz_mod_mat_nullspace(out.raw(), frobenius.raw(), ctx.raw());
    return nullity >= 0;
}

bool pradical_frobenius_kernel(flint::FmpzModMat& out,
                               slong& nullity,
                               const flint::FmpzMat& multiplication_table,
                               const flint::FmpzVec& one,
                               flint::FmpzConstRef p,
                               const flint::FmpzModCtx& ctx) noexcept {
    // Source trace: reference `Ideal.jl:pradical_frobenius` computes
    // Ker(x -> x^(p^j)) in O/pO for p^j >= degree(O); reference
    // `base2.c:pradical` builds the same Frobenius-power kernel before
    // decomposing the etale quotient.
    const slong n = one.length();
    flint::Fmpz p_power;
    flint::FmpzModMat frobenius_power(n, n, ctx);
    flint::FmpzVec base(n);
    flint::FmpzVec image(n);

    if (n <= 0 || out.context() != ctx.raw()) {
        return false;
    }

    fmpz_one(p_power.raw());
    while (fmpz_cmp_ui(p_power.raw(), static_cast<ulong>(n)) < 0) {
        fmpz_mul(p_power.raw(), p_power.raw(), p.raw());
    }

    for (slong j = 0; j < n; ++j) {
        set_identity_coords_mod(base, j, ctx);
        order_coords_pow_mod(image, base, flint::FmpzConstRef(p_power),
                             one, multiplication_table, ctx);

        for (slong i = 0; i < n; ++i) {
            fmpz_set(fmpz_mod_mat_entry(frobenius_power.raw(), i, j),
                     image.data() + i);
        }
    }

    nullity = fmpz_mod_mat_nullspace(out.raw(), frobenius_power.raw(),
                                     ctx.raw());
    return nullity >= 0;
}

bool set_kernel_combination(flint::FmpzVec& out,
                            const flint::FmpzModMat& kernel,
                            slong nullity,
                            slong combination_index,
                            slong p_si,
                            const flint::FmpzModCtx& ctx) noexcept {
    if (out.length() <= 0 || nullity <= 0 || p_si <= 1) {
        return false;
    }

    fmpz_vec_zero(out);
    flint::Fmpz coeff;
    flint::Fmpz entry;
    for (slong col = 0; col < nullity; ++col) {
        const slong digit = combination_index % p_si;
        combination_index /= p_si;
        if (digit == 0) {
            continue;
        }

        fmpz_set_si(coeff.raw(), digit);
        for (slong row = 0; row < out.length(); ++row) {
            fmpz_mod_mat_get_entry(entry.raw(), kernel.raw(), row, col,
                                   ctx.raw());
            fmpz_mod_mul(entry.raw(), entry.raw(), coeff.raw(), ctx.raw());
            fmpz_add(out.data() + row, out.data() + row, entry.raw());
            fmpz_mod_set_fmpz(out.data() + row, out.data() + row,
                              ctx.raw());
        }
    }
    return true;
}

bool multiplication_matrix_mod(flint::FmpzModMat& out,
                               const flint::FmpzVec& element,
                               const flint::FmpzMat& multiplication_table,
                               const flint::FmpzModCtx& ctx) noexcept {
    const slong n = element.length();
    flint::FmpzVec basis(n);
    flint::FmpzVec image(n);
    if (out.context() != ctx.raw()) {
        return false;
    }

    for (slong j = 0; j < n; ++j) {
        set_identity_coords_mod(basis, j, ctx);
        order_coords_mul_mod(image, basis, element, multiplication_table, ctx);
        for (slong i = 0; i < n; ++i) {
            fmpz_set(fmpz_mod_mat_entry(out.raw(), i, j), image.data() + i);
        }
    }
    return true;
}

bool find_separating_fixed_element(
        flint::FmpzVec& out,
        flint::FmpzModPolyFactor& out_factorization,
        flint::FmpzModPoly& out_minpoly,
        const flint::FmpzModMat& kernel,
        slong nullity,
        const flint::FmpzMat& multiplication_table,
        flint::FmpzConstRef p,
        const flint::FmpzModCtx& ctx) noexcept {
    // Source trace: reference `AlgAss/Decompose.jl:_dec_com_finite` chooses a
    // random element of `ker(x -> x^q - x)`, factors its minimal polynomial,
    // and uses CRT to obtain idempotents.  This branch uses bounded
    // deterministic combinations instead so regression tests are stable.
    if (!fmpz_fits_si(p.raw()) || nullity <= 1) {
        return false;
    }
    const slong p_si = fmpz_get_si(p.raw());
    if (p_si <= 1) {
        return false;
    }

    slong attempts = 1;
    for (slong i = 0; i < nullity; ++i) {
        if (attempts > 4096 / p_si) {
            return false;
        }
        attempts *= p_si;
    }

    flint::FmpzVec candidate(out.length());
    flint::FmpzModMat representation(out.length(), out.length(), ctx);
    flint::FmpzModPoly minpoly(ctx);
    flint::FmpzModPolyFactor factorization(ctx);

    for (slong attempt = 1; attempt < attempts; ++attempt) {
        if (!set_kernel_combination(candidate, kernel, nullity, attempt, p_si,
                                    ctx) ||
            !multiplication_matrix_mod(representation, candidate,
                                       multiplication_table, ctx)) {
            return false;
        }

        fmpz_mod_mat_minpoly(minpoly.raw(), representation.raw(), ctx.raw());
        if (fmpz_mod_poly_degree(minpoly.raw(), ctx.raw()) < 2) {
            continue;
        }

        fmpz_mod_poly_factor(factorization.raw(), minpoly.raw(), ctx.raw());
        if (factorization.raw()->num != nullity) {
            continue;
        }

        bool squarefree_linear = true;
        for (slong i = 0; i < factorization.raw()->num; ++i) {
            if (factorization.raw()->exp[i] != 1 ||
                fmpz_mod_poly_degree(factorization.raw()->poly + i,
                                     ctx.raw()) != 1) {
                squarefree_linear = false;
                break;
            }
        }
        if (!squarefree_linear) {
            continue;
        }

        fmpz_vec_set(out, candidate);
        fmpz_mod_poly_set(out_minpoly.raw(), minpoly.raw(), ctx.raw());
        fmpz_mod_poly_factor_set(out_factorization.raw(),
                                 factorization.raw(), ctx.raw());
        return true;
    }

    return false;
}

bool crt_idempotent_polynomial(flint::FmpzModPoly& out,
                               const fmpz_mod_poly_t minpoly,
                               const fmpz_mod_poly_t factor,
                               const flint::FmpzModCtx& ctx) noexcept {
    flint::FmpzModPoly quotient(ctx);
    flint::FmpzModPoly remainder(ctx);
    flint::FmpzModPoly inverse(ctx);
    flint::FmpzModPoly product(ctx);
    if (!quotient.is_initialized() || !remainder.is_initialized() ||
        !inverse.is_initialized() || !product.is_initialized()) {
        return false;
    }

    fmpz_mod_poly_divrem(quotient.raw(), remainder.raw(), minpoly, factor,
                         ctx.raw());
    if (fmpz_mod_poly_is_zero(remainder.raw(), ctx.raw()) == 0 ||
        fmpz_mod_poly_invmod(inverse.raw(), quotient.raw(), factor,
                             ctx.raw()) == 0) {
        return false;
    }

    fmpz_mod_poly_mul(product.raw(), quotient.raw(), inverse.raw(), ctx.raw());
    fmpz_mod_poly_rem(out.raw(), product.raw(), minpoly, ctx.raw());
    return true;
}

bool component_basis_from_idempotent(
        flint::FmpzMat& out,
        slong& dimension,
        const flint::FmpzVec& idempotent,
        const flint::FmpzMat& multiplication_table,
        const flint::FmpzModCtx& ctx) noexcept {
    const slong n = idempotent.length();
    flint::FmpzModMat products(n, n, ctx);
    flint::FmpzModMat rref(n, n, ctx);
    flint::FmpzVec basis(n);
    flint::FmpzVec image(n);
    flint::Fmpz entry;

    for (slong row = 0; row < n; ++row) {
        set_identity_coords_mod(basis, row, ctx);
        order_coords_mul_mod(image, idempotent, basis,
                             multiplication_table, ctx);
        for (slong col = 0; col < n; ++col) {
            fmpz_set(fmpz_mod_mat_entry(products.raw(), row, col),
                     image.data() + col);
        }
    }

    dimension = fmpz_mod_mat_rref(rref.raw(), products.raw(), ctx.raw());
    if (dimension <= 0 || dimension > n) {
        return false;
    }

    flint::FmpzMat basis_rows(dimension, n);
    for (slong i = 0; i < dimension; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_mod_mat_get_entry(entry.raw(), rref.raw(), i, j, ctx.raw());
            fmpz_set(fmpz_mat_entry(basis_rows.raw(), i, j), entry.raw());
        }
    }

    out = std::move(basis_rows);
    return true;
}

bool components_from_fixed_idempotents(
        std::vector<FiniteAlgebraComponent>& out,
        const flint::FmpzModMat& kernel,
        slong nullity,
        const flint::FmpzMat& multiplication_table,
        flint::FmpzConstRef p,
        const flint::FmpzModCtx& ctx) noexcept {
    // Source trace: reference `AlgAss/Decompose.jl:_dec_com_finite` decomposes
    // a finite etale algebra through idempotents.  A single separating
    // element cannot split more components than the base field has elements;
    // for small fixed spaces, enumerate fixed idempotents directly and keep
    // the primitive ones.
    if (!fmpz_fits_si(p.raw()) || nullity <= 1) {
        return false;
    }
    const slong p_si = fmpz_get_si(p.raw());
    if (p_si <= 1) {
        return false;
    }

    slong attempts = 1;
    for (slong i = 0; i < nullity; ++i) {
        if (attempts > 4096 / p_si) {
            return false;
        }
        attempts *= p_si;
    }

    const slong n = fmpz_mod_mat_nrows(kernel.raw(), ctx.raw());
    std::vector<FiniteAlgebraComponent> idempotents;
    flint::FmpzVec candidate(n);
    flint::FmpzVec square(n);

    for (slong attempt = 1; attempt < attempts; ++attempt) {
        if (!set_kernel_combination(candidate, kernel, nullity, attempt,
                                    p_si, ctx)) {
            return false;
        }
        if (fmpz_vec_is_zero(candidate)) {
            continue;
        }
        order_coords_mul_mod(square, candidate, candidate,
                             multiplication_table, ctx);
        if (!fmpz_vec_equal(square, candidate)) {
            continue;
        }

        idempotents.emplace_back(n);
        FiniteAlgebraComponent& component = idempotents.back();
        fmpz_vec_set(component.idempotent, candidate);
        if (!component_basis_from_idempotent(
                    component.basis, component.dimension,
                    component.idempotent, multiplication_table, ctx)) {
            return false;
        }
    }

    if (idempotents.empty()) {
        return false;
    }

    std::vector<FiniteAlgebraComponent> primitive;
    flint::FmpzVec product(n);
    for (std::size_t i = 0; i < idempotents.size(); ++i) {
        bool is_primitive = true;
        for (std::size_t j = 0; j < idempotents.size(); ++j) {
            if (i == j ||
                idempotents[j].dimension >= idempotents[i].dimension) {
                continue;
            }
            order_coords_mul_mod(product, idempotents[j].idempotent,
                                 idempotents[i].idempotent,
                                 multiplication_table, ctx);
            if (fmpz_vec_equal(product, idempotents[j].idempotent)) {
                is_primitive = false;
                break;
            }
        }
        if (!is_primitive) {
            continue;
        }

        primitive.emplace_back(n);
        FiniteAlgebraComponent& component = primitive.back();
        fmpz_vec_set(component.idempotent, idempotents[i].idempotent);
        component.dimension = idempotents[i].dimension;
        component.basis = flint::FmpzMat(component.dimension, n);
        for (slong row = 0; row < component.dimension; ++row) {
            for (slong col = 0; col < n; ++col) {
                fmpz_set(fmpz_mat_entry(component.basis.raw(), row, col),
                         fmpz_mat_entry(idempotents[i].basis.raw(), row,
                                        col));
            }
        }
    }

    slong dimension_sum = 0;
    for (const FiniteAlgebraComponent& component : primitive) {
        dimension_sum += component.dimension;
    }
    if (dimension_sum != n) {
        return false;
    }

    out.swap(primitive);
    return true;
}

bool components_from_fixed_subspace(
        std::vector<FiniteAlgebraComponent>& out,
        const flint::FmpzModMat& kernel,
        slong nullity,
        const flint::FmpzVec& one,
        const flint::FmpzMat& multiplication_table,
        flint::FmpzConstRef p,
        const flint::FmpzModCtx& ctx) noexcept {
    const slong n = one.length();
    flint::FmpzVec separator(n);
    flint::FmpzModPoly minpoly(ctx);
    flint::FmpzModPolyFactor factorization(ctx);
    if (!find_separating_fixed_element(separator, factorization, minpoly,
                                       kernel, nullity, multiplication_table,
                                       p, ctx)) {
        return components_from_fixed_idempotents(
                out, kernel, nullity, multiplication_table, p, ctx);
    }

    std::vector<FiniteAlgebraComponent> components;
    components.reserve(static_cast<std::size_t>(nullity));
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        flint::FmpzModPoly idempotent_polynomial(ctx);
        if (!crt_idempotent_polynomial(idempotent_polynomial, minpoly.raw(),
                                       factorization.raw()->poly + i, ctx)) {
            return false;
        }

        components.emplace_back(n);
        FiniteAlgebraComponent& component = components.back();
        algebra_poly_eval_mod(component.idempotent,
                              idempotent_polynomial.raw(), separator, one,
                              multiplication_table, ctx);
        if (!component_basis_from_idempotent(
                    component.basis, component.dimension,
                    component.idempotent, multiplication_table, ctx)) {
            return false;
        }
    }

    slong dimension_sum = 0;
    for (const FiniteAlgebraComponent& component : components) {
        dimension_sum += component.dimension;
    }
    if (dimension_sum != n) {
        return false;
    }

    out.swap(components);
    return true;
}

bool ideal_from_mod_rows(Ideal& out,
                         const Order& order,
                         flint::FmpzConstRef p,
                         const std::vector<FiniteAlgebraComponent>& components,
                         slong omitted_component) noexcept {
    const slong n = order.degree();
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        omitted_component < 0 ||
        omitted_component >= static_cast<slong>(components.size())) {
        return false;
    }

    slong extra_rows = 0;
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        if (i != omitted_component) {
            extra_rows += components[static_cast<std::size_t>(i)].dimension;
        }
    }

    flint::FmpzMat rows(n + extra_rows, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i), p.raw());
    }

    slong out_row = n;
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        if (i == omitted_component) {
            continue;
        }
        const FiniteAlgebraComponent& component =
                components[static_cast<std::size_t>(i)];
        for (slong r = 0; r < component.dimension; ++r) {
            for (slong c = 0; c < n; ++c) {
                fmpz_set(fmpz_mat_entry(rows.raw(), out_row, c),
                         fmpz_mat_entry(component.basis.raw(), r, c));
            }
            ++out_row;
        }
    }

    lat::Lat row_lattice(n);
    lat::Lat hnf_lattice(n);
    if (!row_lattice.set_basis(flint::FmpzMatConstRef(rows)) ||
        !row_lattice.hnf(hnf_lattice) || hnf_lattice.nrows() != n) {
        return false;
    }

    return out.set_hnf(hnf_lattice.basis_ref());
}

bool ideal_from_mod_kernel(Ideal& out,
                           const Order& order,
                           flint::FmpzConstRef p,
                           const flint::FmpzModMat& kernel,
                           slong nullity,
                           const flint::FmpzModCtx& ctx) noexcept {
    const slong n = order.degree();
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        nullity < 0 || nullity > n ||
        fmpz_mod_mat_nrows(kernel.raw(), ctx.raw()) != n ||
        fmpz_mod_mat_ncols(kernel.raw(), ctx.raw()) < nullity) {
        return false;
    }

    flint::FmpzMat rows(n + nullity, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i), p.raw());
    }

    flint::Fmpz entry;
    for (slong col = 0; col < nullity; ++col) {
        for (slong row = 0; row < n; ++row) {
            fmpz_mod_mat_get_entry(entry.raw(), kernel.raw(), row, col,
                                   ctx.raw());
            fmpz_set(fmpz_mat_entry(rows.raw(), n + col, row), entry.raw());
        }
    }

    lat::Lat row_lattice(n);
    lat::Lat hnf_lattice(n);
    if (!row_lattice.set_basis(flint::FmpzMatConstRef(rows)) ||
        !row_lattice.hnf(hnf_lattice) || hnf_lattice.nrows() != n) {
        return false;
    }

    return out.set_hnf(hnf_lattice.basis_ref());
}

void fmpz_mod_mat_column_to_vec(flint::FmpzVec& out,
                                const flint::FmpzModMat& matrix,
                                slong column,
                                const flint::FmpzModCtx& ctx) noexcept {
    for (slong row = 0; row < out.length(); ++row) {
        fmpz_mod_mat_get_entry(out.data() + row, matrix.raw(), row, column,
                               ctx.raw());
    }
}

void fmpz_mod_mat_set_column_from_kernel(flint::FmpzModMat& out,
                                         slong out_column,
                                         const flint::FmpzModMat& kernel,
                                         slong kernel_column,
                                         const flint::FmpzModCtx& ctx)
        noexcept {
    flint::Fmpz entry;
    const slong n = fmpz_mod_mat_nrows(out.raw(), ctx.raw());
    for (slong row = 0; row < n; ++row) {
        fmpz_mod_mat_get_entry(entry.raw(), kernel.raw(), row,
                               kernel_column, ctx.raw());
        fmpz_set(fmpz_mod_mat_entry(out.raw(), row, out_column),
                 entry.raw());
    }
}

void fmpz_mod_mat_set_identity_column(flint::FmpzModMat& out,
                                      slong column,
                                      slong index,
                                      const flint::FmpzModCtx& ctx)
        noexcept {
    const slong n = fmpz_mod_mat_nrows(out.raw(), ctx.raw());
    for (slong row = 0; row < n; ++row) {
        if (row == index) {
            fmpz_mod_set_ui(fmpz_mod_mat_entry(out.raw(), row, column), 1,
                            ctx.raw());
        } else {
            fmpz_zero(fmpz_mod_mat_entry(out.raw(), row, column));
        }
    }
}

void fmpz_mod_mat_copy_prefix_columns(flint::FmpzModMat& out,
                                      const flint::FmpzModMat& input,
                                      slong columns,
                                      const flint::FmpzModCtx& ctx)
        noexcept {
    flint::Fmpz entry;
    const slong rows = fmpz_mod_mat_nrows(input.raw(), ctx.raw());
    for (slong col = 0; col < columns; ++col) {
        for (slong row = 0; row < rows; ++row) {
            fmpz_mod_mat_get_entry(entry.raw(), input.raw(), row, col,
                                   ctx.raw());
            fmpz_set(fmpz_mod_mat_entry(out.raw(), row, col), entry.raw());
        }
    }
}

bool complete_quotient_basis(flint::FmpzModMat& basis,
                             flint::FmpzModMat& inverse,
                             const flint::FmpzModMat& radical_kernel,
                             slong radical_dimension,
                             const flint::FmpzModCtx& ctx) noexcept {
    const slong n = fmpz_mod_mat_nrows(radical_kernel.raw(), ctx.raw());
    if (radical_dimension < 0 || radical_dimension > n ||
        fmpz_mod_mat_nrows(basis.raw(), ctx.raw()) != n ||
        fmpz_mod_mat_ncols(basis.raw(), ctx.raw()) != n ||
        fmpz_mod_mat_nrows(inverse.raw(), ctx.raw()) != n ||
        fmpz_mod_mat_ncols(inverse.raw(), ctx.raw()) != n) {
        return false;
    }

    for (slong col = 0; col < radical_dimension; ++col) {
        fmpz_mod_mat_set_column_from_kernel(basis, col, radical_kernel, col,
                                            ctx);
    }

    slong columns = radical_dimension;
    for (slong standard = 0; standard < n && columns < n; ++standard) {
        fmpz_mod_mat_set_identity_column(basis, columns, standard, ctx);

        flint::FmpzModMat trial(n, columns + 1, ctx);
        fmpz_mod_mat_copy_prefix_columns(trial, basis, columns + 1, ctx);
        const slong rank = fmpz_mod_mat_rank(trial.raw(), ctx.raw());
        if (rank == columns + 1) {
            ++columns;
        }
    }

    if (columns != n) {
        return false;
    }

    return fmpz_mod_mat_inv(inverse.raw(), basis.raw(), ctx.raw()) != 0;
}

void quotient_coords_from_lift(flint::FmpzVec& out,
                               const flint::FmpzVec& lift,
                               const flint::FmpzModMat& inverse_basis,
                               slong radical_dimension,
                               const flint::FmpzModCtx& ctx) noexcept {
    const slong quotient_dimension = out.length();
    const slong n = quotient_dimension + radical_dimension;
    flint::Fmpz entry;
    flint::Fmpz term;
    for (slong row = 0; row < quotient_dimension; ++row) {
        fmpz_zero(out.data() + row);
        for (slong col = 0; col < n; ++col) {
            fmpz_mod_mat_get_entry(entry.raw(), inverse_basis.raw(),
                                   radical_dimension + row, col,
                                   ctx.raw());
            fmpz_mod_mul(term.raw(), entry.raw(), lift.data() + col,
                         ctx.raw());
            fmpz_add(out.data() + row, out.data() + row, term.raw());
            fmpz_mod_set_fmpz(out.data() + row, out.data() + row,
                              ctx.raw());
        }
    }
}

bool quotient_algebra_table(flint::FmpzMat& quotient_table,
                            flint::FmpzVec& quotient_one,
                            const flint::FmpzMat& multiplication_table,
                            const flint::FmpzVec& one,
                            const flint::FmpzModMat& basis,
                            const flint::FmpzModMat& inverse_basis,
                            slong radical_dimension,
                            const flint::FmpzModCtx& ctx) noexcept {
    const slong quotient_dimension = quotient_one.length();
    if (!fmpz_mat_shape(flint::FmpzMatConstRef(quotient_table),
                        quotient_dimension * quotient_dimension,
                        quotient_dimension)) {
        return false;
    }

    quotient_coords_from_lift(quotient_one, one, inverse_basis,
                              radical_dimension, ctx);

    flint::FmpzVec left(quotient_dimension + radical_dimension);
    flint::FmpzVec right(quotient_dimension + radical_dimension);
    flint::FmpzVec product(quotient_dimension + radical_dimension);
    flint::FmpzVec projected(quotient_dimension);
    for (slong i = 0; i < quotient_dimension; ++i) {
        fmpz_mod_mat_column_to_vec(left, basis, radical_dimension + i, ctx);
        for (slong j = 0; j < quotient_dimension; ++j) {
            fmpz_mod_mat_column_to_vec(right, basis, radical_dimension + j,
                                       ctx);
            order_coords_mul_mod(product, left, right, multiplication_table,
                                 ctx);
            quotient_coords_from_lift(projected, product, inverse_basis,
                                      radical_dimension, ctx);
            for (slong k = 0; k < quotient_dimension; ++k) {
                fmpz_set(fmpz_mat_entry(quotient_table.raw(),
                                        i * quotient_dimension + j, k),
                         projected.data() + k);
            }
        }
    }
    return true;
}

void lift_quotient_component_row(flint::FmpzVec& out,
                                 const FiniteAlgebraComponent& component,
                                 slong component_row,
                                 const flint::FmpzModMat& basis,
                                 slong radical_dimension,
                                 const flint::FmpzModCtx& ctx) noexcept {
    const slong n = out.length();
    const slong quotient_dimension = n - radical_dimension;
    flint::Fmpz coeff;
    flint::Fmpz entry;
    flint::Fmpz term;
    fmpz_vec_zero(out);
    for (slong q = 0; q < quotient_dimension; ++q) {
        fmpz_mod_set_fmpz(
                coeff.raw(),
                fmpz_mat_entry(component.basis.raw(), component_row, q),
                ctx.raw());
        if (fmpz_is_zero(coeff.raw()) != 0) {
            continue;
        }
        for (slong row = 0; row < n; ++row) {
            fmpz_mod_mat_get_entry(entry.raw(), basis.raw(), row,
                                   radical_dimension + q, ctx.raw());
            fmpz_mod_mul(term.raw(), coeff.raw(), entry.raw(), ctx.raw());
            fmpz_add(out.data() + row, out.data() + row, term.raw());
            fmpz_mod_set_fmpz(out.data() + row, out.data() + row,
                              ctx.raw());
        }
    }
}

bool ideal_from_pradical_quotient_component(
        Ideal& out,
        const Order& order,
        flint::FmpzConstRef p,
        const flint::FmpzModMat& radical_kernel,
        slong radical_dimension,
        const flint::FmpzModMat& quotient_basis,
        const std::vector<FiniteAlgebraComponent>& components,
        slong omitted_component,
        const flint::FmpzModCtx& ctx) noexcept {
    const slong n = order.degree();
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        radical_dimension < 0 || radical_dimension > n ||
        omitted_component < 0 ||
        omitted_component >= static_cast<slong>(components.size())) {
        return false;
    }

    slong extra_rows = radical_dimension;
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        if (i != omitted_component) {
            extra_rows += components[static_cast<std::size_t>(i)].dimension;
        }
    }

    flint::FmpzMat rows(n + extra_rows, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i), p.raw());
    }

    slong out_row = n;
    flint::Fmpz entry;
    for (slong col = 0; col < radical_dimension; ++col) {
        for (slong row = 0; row < n; ++row) {
            fmpz_mod_mat_get_entry(entry.raw(), radical_kernel.raw(), row,
                                   col, ctx.raw());
            fmpz_set(fmpz_mat_entry(rows.raw(), out_row, row), entry.raw());
        }
        ++out_row;
    }

    flint::FmpzVec lifted(n);
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        if (i == omitted_component) {
            continue;
        }
        const FiniteAlgebraComponent& component =
                components[static_cast<std::size_t>(i)];
        for (slong r = 0; r < component.dimension; ++r) {
            lift_quotient_component_row(lifted, component, r,
                                        quotient_basis, radical_dimension,
                                        ctx);
            for (slong c = 0; c < n; ++c) {
                fmpz_set(fmpz_mat_entry(rows.raw(), out_row, c),
                         lifted.data() + c);
            }
            ++out_row;
        }
    }

    lat::Lat row_lattice(n);
    lat::Lat hnf_lattice(n);
    if (!row_lattice.set_basis(flint::FmpzMatConstRef(rows)) ||
        !row_lattice.hnf(hnf_lattice) || hnf_lattice.nrows() != n) {
        return false;
    }

    return out.set_hnf(hnf_lattice.basis_ref());
}

bool ramification_index_from_containment(slong& out,
                                         const Ideal& prime_ideal,
                                         const Order& order,
                                         flint::FmpzConstRef p) noexcept {
    OrderElement p_element(order);
    if (!p_element.is_defined() ||
        !order_element_set_fmpz(p_element, order, p) ||
        !valuation_by_ideal_power_containment(
                out, prime_ideal, p_element, order.degree(), nullptr) ||
        out <= 0) {
        return false;
    }
    return true;
}

bool unramified_maximal_index_divisor_prime_data(
        std::vector<IndexPrimeData>& out,
        const Order& order,
        flint::FmpzConstRef p) noexcept {
    // Source trace: reference `Prime.jl:prime_dec_gen` forms `O / pradical(O,p)`,
    // `Polygons.jl:_decomposition` splits that finite algebra, and for
    // p not dividing the maximal-order discriminant sets every e to 1.  In
    // this unramified case the p-radical is pO, so the quotient is O/pO.
    if (!order.is_maximal() || order.is_equation_order()) {
        return false;
    }

    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        fmpz_divisible(discriminant.raw(), p.raw()) != 0) {
        return false;
    }

    const slong n = order.degree();
    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzMat multiplication_table(n * n, n);
    flint::FmpzVec one(n);
    flint::FmpzVec one_mod(n);
    if (!build_multiplication_table(multiplication_table, order) ||
        !order_one_coords(one, order)) {
        return false;
    }
    fmpz_vec_set_mod(one_mod, one, ctx);

    flint::FmpzModMat fixed_kernel(n, n, ctx);
    slong nullity = 0;
    if (!fixed_frobenius_kernel(fixed_kernel, nullity,
                                multiplication_table, one_mod, p, ctx) ||
        nullity <= 0) {
        return false;
    }

    if (nullity == 1) {
        std::vector<IndexPrimeData> data;
        data.emplace_back(order);
        IndexPrimeData& prime = data.back();
        OrderElement p_element(order);
        if (!prime.ideal.is_defined() || !p_element.is_defined() ||
            !order_element_set_fmpz(p_element, order, p) ||
            !prime.ideal.set_principal(p_element)) {
            return false;
        }
        prime.residue_degree = n;

        out.swap(data);
        return true;
    }

    std::vector<FiniteAlgebraComponent> components;
    if (!components_from_fixed_subspace(components, fixed_kernel, nullity,
                                        one_mod, multiplication_table, p,
                                        ctx)) {
        return false;
    }

    std::vector<IndexPrimeData> data;
    data.reserve(components.size());
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        const FiniteAlgebraComponent& component =
                components[static_cast<std::size_t>(i)];
        data.emplace_back(order);
        IndexPrimeData& prime = data.back();
        if (!prime.ideal.is_defined() ||
            !ideal_from_mod_rows(prime.ideal, order, p, components, i)) {
            return false;
        }
        prime.residue_degree = component.dimension;
    }

    out.swap(data);
    return true;
}

bool ramified_maximal_index_divisor_prime_data(
        std::vector<IndexPrimeData>& out,
        const Order& order,
        flint::FmpzConstRef p) noexcept {
    // Source trace: reference `Prime.jl:prime_dec_gen` forms the p-radical and
    // decomposes O/(p, pradical); reference `base2.c:primedec_aux` builds the
    // same p-radical, splits the etale quotient, and recovers ramification
    // indices from valuations.  Native uses the existing finite algebra
    // decomposition helpers for the quotient and exact containment valuation
    // for v_P(p).
    if (!order.is_maximal() || order.is_equation_order()) {
        return false;
    }

    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        fmpz_divisible(discriminant.raw(), p.raw()) == 0) {
        return false;
    }

    const slong n = order.degree();
    flint::FmpzModCtx ctx(p.raw());

    flint::FmpzMat multiplication_table(n * n, n);
    flint::FmpzVec one(n);
    flint::FmpzVec one_mod(n);
    if (!build_multiplication_table(multiplication_table, order) ||
        !order_one_coords(one, order)) {
        return false;
    }
    fmpz_vec_set_mod(one_mod, one, ctx);

    flint::FmpzModMat radical_kernel(n, n, ctx);
    slong radical_nullity = 0;
    if (!pradical_frobenius_kernel(radical_kernel, radical_nullity,
                                   multiplication_table, one_mod, p, ctx)) {
        return false;
    }

    const slong residue_degree = n - radical_nullity;
    if (radical_nullity <= 0 || residue_degree <= 0) {
        return false;
    }

    flint::FmpzModMat quotient_basis(n, n, ctx);
    flint::FmpzModMat inverse_basis(n, n, ctx);
    if (!complete_quotient_basis(quotient_basis, inverse_basis,
                                 radical_kernel, radical_nullity, ctx)) {
        return false;
    }

    flint::FmpzMat quotient_multiplication_table(
            residue_degree * residue_degree, residue_degree);
    flint::FmpzVec quotient_one(residue_degree);
    if (!quotient_algebra_table(quotient_multiplication_table, quotient_one,
                                multiplication_table, one_mod,
                                quotient_basis, inverse_basis,
                                radical_nullity, ctx)) {
        return false;
    }

    flint::FmpzModMat fixed_kernel(residue_degree, residue_degree, ctx);
    slong fixed_nullity = 0;
    if (!fixed_frobenius_kernel(fixed_kernel, fixed_nullity,
                                quotient_multiplication_table,
                                quotient_one, p, ctx) ||
        fixed_nullity <= 0) {
        return false;
    }

    if (fixed_nullity == 1) {
        std::vector<IndexPrimeData> data;
        data.emplace_back(order);
        IndexPrimeData& prime = data.back();
        if (!prime.ideal.is_defined() ||
            !ideal_from_mod_kernel(prime.ideal, order, p, radical_kernel,
                                   radical_nullity, ctx) ||
            !ramification_index_from_containment(
                    prime.ramification_index, prime.ideal, order, p)) {
            return false;
        }
        prime.residue_degree = residue_degree;

        out.swap(data);
        return true;
    }

    std::vector<FiniteAlgebraComponent> components;
    if (!components_from_fixed_subspace(
                components, fixed_kernel, fixed_nullity, quotient_one,
                quotient_multiplication_table, p, ctx)) {
        return false;
    }

    std::vector<IndexPrimeData> data;
    data.reserve(components.size());
    for (slong i = 0; i < static_cast<slong>(components.size()); ++i) {
        const FiniteAlgebraComponent& component =
                components[static_cast<std::size_t>(i)];
        data.emplace_back(order);
        IndexPrimeData& prime = data.back();
        if (!prime.ideal.is_defined() ||
            !ideal_from_pradical_quotient_component(
                    prime.ideal, order, p, radical_kernel,
                    radical_nullity, quotient_basis, components, i, ctx) ||
            !ramification_index_from_containment(
                    prime.ramification_index, prime.ideal, order, p)) {
            return false;
        }
        prime.residue_degree = component.dimension;
    }

    out.swap(data);
    return true;
}

bool monogenic_candidate(flint::FmpzPoly& chi,
                         Element& beta,
                         const Order& order,
                         const flint::FmpzVec& candidate,
                         const flint::FmpzMat& multiplication_table) noexcept {
    const slong n = order.degree();
    flint::FmpzMat powers(n, n);
    flint::FmpzMat multiplication_by_beta(n, n);
    flint::FmpzVec basis_vec(n);
    flint::FmpzVec image(n);
    flint::FmpzVec power(n);
    flint::FmpzVec next(n);
    flint::Fmpz det;

    if (!order_one_coords(power, order)) {
        return false;
    }

    fmpz_mat_set_row(powers, 0, power);
    for (slong i = 1; i < n; ++i) {
        order_coords_mul(next, power, candidate, multiplication_table);
        fmpz_mat_set_row(powers, i, next);
        fmpz_vec_set(power, next);
    }

    fmpz_mat_det(det.raw(), powers.raw());
    fmpz_abs(det.raw(), det.raw());
    if (fmpz_is_one(det.raw()) == 0) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        set_identity_coords(basis_vec, i);
        order_coords_mul(image, basis_vec, candidate, multiplication_table);
        fmpz_mat_set_row(multiplication_by_beta, i, image);
    }

    fmpz_mat_charpoly(chi.raw(), multiplication_by_beta.raw());
    return element_from_order_coords(beta, order, candidate);
}

bool find_monogenic_generator(flint::FmpzPoly& chi,
                              Element& beta,
                              const Order& order) noexcept {
    const slong n = order.degree();
    if (n <= 1 || order.parent() == nullptr ||
        !beta.has_parent(*order.parent())) {
        return false;
    }

    flint::FmpzMat multiplication_table(n * n, n);
    flint::FmpzVec candidate(n);
    if (!build_multiplication_table(multiplication_table, order)) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        set_identity_coords(candidate, i);
        if (monogenic_candidate(chi, beta, order, candidate,
                                multiplication_table)) {
            return true;
        }
    }

    if (n > 7) {
        return false;
    }

    slong total = 1;
    for (slong i = 0; i < n; ++i) {
        total *= 3;
    }

    for (slong t = 1; t < total; ++t) {
        set_small_trinary_coords(candidate, t);
        if (monogenic_candidate(chi, beta, order, candidate,
                                multiplication_table)) {
            return true;
        }
    }

    return false;
}

bool element_evaluate_mod_poly_at(Element& out,
                                  const fmpz_mod_poly_t polynomial,
                                  const Element& beta,
                                  const flint::FmpzModCtx& ctx) noexcept {
    if (out.parent() == nullptr || beta.parent() == nullptr ||
        !out.parent()->has_same_data(*beta.parent())) {
        return false;
    }

    Element result(*beta.parent());
    Element product(*beta.parent());
    Element coefficient_element(*beta.parent());
    flint::FmpqPoly coefficient_polynomial;
    flint::Fmpz coefficient;
    if (!result.is_defined() || !product.is_defined() ||
        !coefficient_element.is_defined() || !result.zero()) {
        return false;
    }

    const slong degree = fmpz_mod_poly_degree(polynomial, ctx.raw());
    for (slong i = degree; i >= 0; --i) {
        if (i != degree) {
            if (!product.multiply(result, beta) || !result.set(product)) {
                return false;
            }
        }

        fmpq_poly_zero(coefficient_polynomial.raw());
        fmpz_mod_poly_get_coeff_fmpz(coefficient.raw(), polynomial, i,
                                     ctx.raw());
        fmpq_poly_set_coeff_fmpz(coefficient_polynomial.raw(), 0,
                                 coefficient.raw());
        if (!coefficient_element.set_fmpq_poly(
                    flint::FmpqPolyConstRef(coefficient_polynomial)) ||
            !result.add(result, coefficient_element)) {
            return false;
        }
    }

    return out.set(result);
}

}  // namespace

struct PrimeIdeal::HenselValuationCache {
    HenselValuationCache(const fmpz_t modulus, slong exponent_value) noexcept
        : target_ctx(modulus),
          lifted_factor(target_ctx),
          exponent(exponent_value) {
        fmpz_set(target_modulus.raw(), modulus);
    }

    bool is_initialized() const noexcept {
        return target_ctx.raw() != nullptr && lifted_factor.is_initialized();
    }

    flint::Fmpz target_modulus;
    flint::FmpzModCtx target_ctx;
    flint::FmpzModPoly lifted_factor;
    slong exponent = 0;
};

struct PrimeIdeal::ContainmentPowerCache {
    std::vector<Ideal> powers;
};

struct PrimeIdeal::CoordinateValuationCache {
    explicit CoordinateValuationCache(slong degree) noexcept
        : multiplier_matrix(degree, degree) {}

    flint::FmpzMat multiplier_matrix;
};

PrimeIdeal::PrimeIdeal() noexcept = default;

PrimeIdeal::PrimeIdeal(const Order& parent) noexcept {
    define(parent);
}

PrimeIdeal::~PrimeIdeal() noexcept {
    clear();
}

PrimeIdeal::PrimeIdeal(PrimeIdeal&& other) noexcept {
    swap(other);
}

PrimeIdeal& PrimeIdeal::operator=(PrimeIdeal&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void PrimeIdeal::swap(PrimeIdeal& other) noexcept {
    parent_.swap(other.parent_);
    p_.swap(other.p_);
    ideal_.swap(other.ideal_);
    kummer_generator_.swap(other.kummer_generator_);
    residue_poly_.swap(other.residue_poly_);
    linear_residue_root_.swap(other.linear_residue_root_);
    hensel_valuation_cache_.swap(other.hensel_valuation_cache_);
    containment_power_cache_.swap(other.containment_power_cache_);
    coordinate_valuation_cache_.swap(
            other.coordinate_valuation_cache_);
    std::swap(e_, other.e_);
    std::swap(f_, other.f_);
    std::swap(has_prime_, other.has_prime_);
    std::swap(has_kummer_generator_, other.has_kummer_generator_);
    std::swap(has_residue_poly_, other.has_residue_poly_);
    std::swap(has_linear_residue_root_,
              other.has_linear_residue_root_);
}

void PrimeIdeal::clear() noexcept {
    parent_.clear();
    fmpz_zero(p_.raw());
    ideal_.clear();
    kummer_generator_ = flint::FmpzMat(0, 0);
    fmpz_poly_zero(residue_poly_.raw());
    fmpz_zero(linear_residue_root_.raw());
    hensel_valuation_cache_.reset();
    containment_power_cache_.reset();
    coordinate_valuation_cache_.reset();
    e_ = 0;
    f_ = 0;
    has_prime_ = false;
    has_kummer_generator_ = false;
    has_residue_poly_ = false;
    has_linear_residue_root_ = false;
}

bool PrimeIdeal::define(const Order& parent) noexcept {
    if (!parent.has_basis()) {
        return false;
    }

    PrimeIdeal next;
    next.parent_ = parent;
    next.ideal_ = Ideal(parent);
    if (!next.ideal_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool PrimeIdeal::set(const PrimeIdeal& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    PrimeIdeal copy(other.parent_);
    if (!copy.is_defined()) {
        return false;
    }

    fmpz_set(copy.p_.raw(), other.p_.raw());
    if (!copy.ideal_.set(other.ideal_)) {
        return false;
    }
    copy.kummer_generator_ = flint::FmpzMat(
            fmpz_mat_nrows(other.kummer_generator_.raw()),
            fmpz_mat_ncols(other.kummer_generator_.raw()));
    fmpz_mat_set(copy.kummer_generator_.raw(), other.kummer_generator_.raw());
    fmpz_poly_set(copy.residue_poly_.raw(),
                  other.residue_poly_.raw());
    fmpz_set(copy.linear_residue_root_.raw(),
             other.linear_residue_root_.raw());
    copy.e_ = other.e_;
    copy.f_ = other.f_;
    copy.has_prime_ = other.has_prime_;
    copy.has_kummer_generator_ = other.has_kummer_generator_;
    copy.has_residue_poly_ = other.has_residue_poly_;
    copy.has_linear_residue_root_ = other.has_linear_residue_root_;

    swap(copy);
    return true;
}

bool PrimeIdeal::is_defined() const noexcept {
    return parent_.is_defined() && ideal_.is_defined();
}

bool PrimeIdeal::has_prime_data() const noexcept {
    return is_defined() && has_prime_ && ideal_.has_hnf();
}

const flint::FmpzMat*
PrimeIdeal::coordinate_valuation_matrix_cached() const noexcept {
    if (coordinate_valuation_cache_ != nullptr) {
        return &coordinate_valuation_cache_->multiplier_matrix;
    }
    if (!has_prime_data() || parent_.parent() == nullptr) {
        return nullptr;
    }

    std::unique_ptr<CoordinateValuationCache> next(
            new (std::nothrow)
                    CoordinateValuationCache(parent_.degree()));
    if (next == nullptr) {
        return nullptr;
    }

    bool initialized = false;
    if (has_residue_poly_ &&
        (parent_.is_equation_order() || parent_.is_maximal())) {
        flint::Fmpz radicand;
        flint::Fmpz conductor;
        if (parent_.is_maximal() && parent_.degree() == 2 &&
            parent_.parent()->backend_kind() ==
                    NumberFieldBackendKind::quadratic &&
            parent_.parent()->quadratic_radicand(flint::FmpzRef(radicand)) &&
            parent_.quadratic_conductor(flint::FmpzRef(conductor)) &&
            flint::fmpz_is_one(conductor)) {
            initialized = quadratic_quotient_multiplier_matrix(
                    next->multiplier_matrix, parent_,
                    flint::FmpzConstRef(radicand), residue_poly_,
                    flint::FmpzConstRef(p_));
        } else {
            const nf_struct* raw_field = parent_.parent()->raw_flint_field();
            initialized = raw_field != nullptr &&
                          quotient_multiplier_matrix(
                                  next->multiplier_matrix, parent_,
                                  raw_field->pol, residue_poly_,
                                  flint::FmpzConstRef(p_));
        }
    }
    if (!initialized) {
        return nullptr;
    }

    coordinate_valuation_cache_ = std::move(next);
    return &coordinate_valuation_cache_->multiplier_matrix;
}

const Ideal* PrimeIdeal::containment_power_cached(
        slong exponent,
        const DiagnosticsContext* diagnostics) const noexcept {
    // Source trace: reference attaches reusable valuation state to a prime ideal
    // (`assure_valuation_function`).  This keeps the existing exact
    // containment contract, but reuses the exact P^k ideals that the fallback
    // would otherwise rebuild for repeated valuations at the same prime.
    if (!has_prime_data() || exponent <= 0) {
        return nullptr;
    }
    if (exponent == 1) {
        return &ideal_;
    }

    const std::size_t target =
            static_cast<std::size_t>(exponent - static_cast<slong>(2));
    if (containment_power_cache_ != nullptr &&
        target < containment_power_cache_->powers.size() &&
        containment_power_cache_->powers[target].has_hnf()) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_containment_power_cache_hit");
        return &containment_power_cache_->powers[target];
    }

    if (containment_power_cache_ == nullptr) {
        containment_power_cache_.reset(
                new (std::nothrow) ContainmentPowerCache());
        if (containment_power_cache_ == nullptr) {
            return nullptr;
        }
    }
    ContainmentPowerCache& cache = *containment_power_cache_;

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::prime_ideal,
            "prime_ideal.valuation_containment_power_cache_build");
    while (cache.powers.size() <= target) {
        Ideal next(parent_);
        if (!next.is_defined()) {
            return nullptr;
        }

        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::prime_ideal,
                    "prime_ideal.valuation_containment_power_cache_multiply");
            if (cache.powers.empty()) {
                if (!next.multiply(ideal_, ideal_)) {
                    return nullptr;
                }
            } else if (!next.multiply(cache.powers.back(), ideal_)) {
                return nullptr;
            }
        }

        cache.powers.push_back(std::move(next));
    }

    return &cache.powers[target];
}

const Order* PrimeIdeal::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong PrimeIdeal::degree() const noexcept {
    return is_defined() ? parent_.degree() : 0;
}

bool PrimeIdeal::rational_prime(flint::FmpzRef out) const noexcept {
    if (!has_prime_data()) {
        return false;
    }
    fmpz_set(out.raw(), p_.raw());
    return true;
}

bool PrimeIdeal::get_ideal(Ideal& out) const noexcept {
    return has_prime_data() && same_order_parent(out.parent(), &parent_) &&
           out.set(ideal_);
}

bool PrimeIdeal::kummer_generator_coordinates(
        flint::FmpzMatRef out) const noexcept {
    if (!has_prime_data() || !has_kummer_generator_ ||
        fmpz_mat_nrows(out.raw()) != 1 ||
        fmpz_mat_ncols(out.raw()) != parent_.degree()) {
        return false;
    }
    fmpz_mat_set(out.raw(), kummer_generator_.raw());
    return true;
}

bool PrimeIdeal::residue_polynomial(flint::FmpzPolyRef out) const noexcept {
    if (!has_prime_data() || !has_residue_poly_) {
        return false;
    }
    fmpz_poly_set(out.raw(), residue_poly_.raw());
    return true;
}

bool PrimeIdeal::norm(flint::FmpzRef out) const noexcept {
    if (!has_prime_data()) {
        return false;
    }
    fmpz_pow_ui(out.raw(), p_.raw(), static_cast<ulong>(f_));
    return true;
}

slong PrimeIdeal::ramification_index() const noexcept {
    return has_prime_data() ? e_ : 0;
}

slong PrimeIdeal::residue_degree() const noexcept {
    return has_prime_data() ? f_ : 0;
}

bool PrimeIdeal::equal(const PrimeIdeal& other) const noexcept {
    if (!has_prime_data() || !other.has_prime_data() ||
        !parent_.has_same_data(other.parent_) || e_ != other.e_ ||
        f_ != other.f_ ||
        fmpz_equal(p_.raw(), other.p_.raw()) == 0 ||
        !ideal_.equal(other.ideal_) ||
        has_residue_poly_ != other.has_residue_poly_) {
        return false;
    }
    return !has_residue_poly_ ||
           fmpz_poly_equal(residue_poly_.raw(),
                           other.residue_poly_.raw()) != 0;
}

bool PrimeIdeal::reduce(flint::FmpzPolyRef out,
                        const OrderElement& element) const noexcept {
    if (!has_prime_data() || !has_residue_poly_ ||
        !same_order_parent(element.parent(), &parent_)) {
        return false;
    }

    Element ambient(*parent_.parent());
    flint::FmpqMat coordinates(1, parent_.degree());
    flint::FmpzPoly input;
    flint::FmpzPoly result;
    if (!element.get_element(ambient) ||
        !parent_.coordinates(flint::FmpqMatRef(coordinates), ambient) ||
        !fmpq_mat_entries_are_integral(coordinates.raw())) {
        return false;
    }

    integral_coords_to_poly(input, coordinates.raw());
    if (!reduce_poly_mod_residue(result, input, residue_poly_,
                                 flint::FmpzConstRef(p_))) {
        return false;
    }
    fmpz_poly_set(out.raw(), result.raw());
    return true;
}

bool PrimeIdeal::reduce(flint::FmpzPolyRef out,
                        const Element& element) const noexcept {
    if (!is_defined() || parent_.parent() == nullptr ||
        !element.has_parent(*parent_.parent())) {
        return false;
    }

    OrderElement order_element(parent_);
    return order_element.set_element(element) && reduce(out, order_element);
}

bool PrimeIdeal::valuation_by_power_containment(
        slong& out,
        const OrderElement& element,
        slong bound,
        const DiagnosticsContext* diagnostics) const noexcept {
    // Exact definition fallback: largest i with a in P^i.  The only
    // deviation from the free helper is reusing the same exact P^k ideals
    // across repeated valuations for this prime.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_containment");
    if (!has_prime_data() || !same_order_parent(element.parent(), &parent_) ||
        bound < 0) {
        return false;
    }
    if (bound == 0) {
        out = 0;
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                            "prime_ideal.valuation_containment_initial_contains");
        if (!ideal_.contains(element)) {
            out = 0;
            return true;
        }
    }

    slong value = 1;
    while (value < bound) {
        const Ideal* power = containment_power_cached(value + 1, diagnostics);
        if (power == nullptr) {
            return false;
        }

        bool contained = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::prime_ideal,
                    "prime_ideal.valuation_containment_power_contains");
            contained = power->contains(element);
        }
        if (!contained) {
            break;
        }
        ++value;
    }

    out = value;
    return true;
}

bool PrimeIdeal::valuation(slong& out,
                           const OrderElement& element) const noexcept {
    return valuation(out, element, nullptr);
}

bool PrimeIdeal::valuation(
        slong& out,
        const OrderElement& element,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_order_element");
    if (!has_prime_data() ||
        !same_order_parent(element.parent(), &parent_) ||
        parent_.parent() == nullptr || !parent_.is_maximal() || e_ <= 0 ||
        f_ <= 0 || element.equal_si(0)) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics, "prime_ideal.valuation_order_element.error_input");
        return false;
    }

    Element ambient(*parent_.parent());
    flint::Fmpq norm;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                            "prime_ideal.valuation_order_element_norm");
        if (!element.get_element(ambient) ||
            !ambient.norm(flint::FmpqRef(norm)) ||
            fmpz_is_one(fmpq_denref(norm.raw())) == 0) {
            SILEX_PRIME_IDEAL_PROFILE_EVENT(
                    diagnostics,
                    "prime_ideal.valuation_order_element.error_norm");
            return false;
        }
    }

    flint::Fmpz norm_abs;
    slong norm_vp = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                            "prime_ideal.valuation_order_element_norm_vp");
        fmpz_abs(norm_abs.raw(), fmpq_numref(norm.raw()));
        while (fmpz_divisible(norm_abs.raw(), p_.raw()) != 0) {
            if (norm_vp == WORD_MAX) {
                SILEX_PRIME_IDEAL_PROFILE_EVENT(
                        diagnostics,
                        "prime_ideal.valuation_order_element.error_norm_vp_overflow");
                return false;
            }
            fmpz_divexact(norm_abs.raw(), norm_abs.raw(), p_.raw());
            ++norm_vp;
        }
    }

    return valuation_with_norm_vp_impl(out, element, norm_vp, diagnostics);
}

namespace detail {

bool prime_ideal_valuation_with_norm_vp(
        slong& out,
        const PrimeIdeal& prime,
        const OrderElement& element,
        slong norm_vp,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_order_element");
    if (!prime.has_prime_data() ||
        !same_order_parent(element.parent(), &prime.parent_) ||
        prime.parent_.parent() == nullptr || !prime.parent_.is_maximal() ||
        prime.e_ <= 0 || prime.f_ <= 0 || element.equal_si(0) ||
        norm_vp < 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics, "prime_ideal.valuation_order_element.error_input");
        return false;
    }

    return prime.valuation_with_norm_vp_impl(out, element, norm_vp,
                                             diagnostics);
}

}  // namespace detail

bool PrimeIdeal::valuation_with_norm_vp_impl(
        slong& out,
        const OrderElement& element,
        slong norm_vp,
        const DiagnosticsContext* diagnostics) const noexcept {
    // Source trace: reference
    // `NumFieldOrd/NfOrd/Ideal/Valuation.jl:val_func_no_index_small`
    // computes element valuations in maximal no-index, unramified equation
    // orders by Hensel-lifting the residue factor and taking the p-adic
    // content of the reduced coordinate polynomial.  reference
    // `assure_valuation_function` also uses the norm shortcut
    // `valuation(norm(a), p) / f` when a single prime ideal lies over the
    // rational prime (`e*f == degree(O)`) in small-degree orders.  reference
    // `base3.c:ZC_nfvalrem` exposes the same integral-coordinate valuation
    // contract, with a broader anti-uniformizer route that is not ported here.
    if (norm_vp == 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_order_element.path_norm_zero");
        out = 0;
        return true;
    }

    const slong degree = parent_.degree();
    if (degree > 0 && degree < 40 && degree % f_ == 0 &&
        e_ == degree / f_) {
        if (norm_vp % f_ != 0) {
            SILEX_PRIME_IDEAL_PROFILE_EVENT(
                    diagnostics,
                    "prime_ideal.valuation_order_element.error_norm_shortcut");
            return false;
        }
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_order_element.path_norm_shortcut");
        out = norm_vp / f_;
        return true;
    }

    const slong bound = norm_vp / f_;
    if (bound <= 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_order_element.path_bound_zero");
        out = 0;
        return true;
    }

    if (has_residue_poly_ && parent_.is_equation_order() && e_ == 1) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                            "prime_ideal.valuation_no_index_hensel");
        flint::FmpzMat coordinates(1, parent_.degree());
        flint::FmpzPoly input_polynomial;
        if (!element.get_coordinates(flint::FmpzMatRef(coordinates))) {
            SILEX_PRIME_IDEAL_PROFILE_EVENT(
                    diagnostics,
                    "prime_ideal.valuation_order_element.fallback_hensel_coordinates");
        } else {
            integral_coords_to_poly(input_polynomial, coordinates.raw());

            auto cache_for_bound = [&]() noexcept
                    -> const HenselValuationCache* {
                const HenselValuationCache* cache =
                        hensel_valuation_cache_.get();
                if (cache != nullptr &&
                    cache->is_initialized() &&
                    cache->exponent >= bound) {
                    return cache;
                }
                return nullptr;
            };
            auto try_cached_hensel = [&](const HenselValuationCache* cache,
                                         slong& valuation_out) noexcept {
                if (cache == nullptr) {
                    return false;
                }
                flint::FmpzModPoly input_mod(cache->target_ctx);
                flint::FmpzModPoly remainder(cache->target_ctx);
                if (!input_mod.is_initialized() ||
                    !remainder.is_initialized()) {
                    return false;
                }
                fmpz_mod_poly_set_fmpz_poly(input_mod.raw(),
                                            input_polynomial.raw(),
                                            cache->target_ctx.raw());
                fmpz_mod_poly_rem(remainder.raw(), input_mod.raw(),
                                  cache->lifted_factor.raw(),
                                  cache->target_ctx.raw());

                return fmpz_mod_poly_content_valuation_capped(
                        valuation_out, remainder.raw(), cache->target_ctx,
                        p_.raw(), bound);
            };
            slong cached_valuation = 0;
            if (try_cached_hensel(cache_for_bound(), cached_valuation)) {
                SILEX_PRIME_IDEAL_PROFILE_EVENT(
                        diagnostics,
                        "prime_ideal.valuation_order_element.path_no_index_hensel");
                out = cached_valuation;
                return true;
            }

            const nf_struct* raw_field = parent_.parent()->raw_flint_field();
            if (raw_field != nullptr &&
                fmpq_poly_is_monic_integral(raw_field->pol)) {
                flint::FmpzModCtx prime_ctx(p_.raw());
                flint::FmpzModPoly field_mod_prime(prime_ctx);
                flint::FmpzModPoly factor_mod_prime(prime_ctx);
                if (prime_ctx.raw() != nullptr &&
                    field_mod_prime.is_initialized() &&
                    factor_mod_prime.is_initialized()) {
                    fmpq_poly_get_fmpz_mod_poly(field_mod_prime,
                                                raw_field->pol, prime_ctx);
                    fmpz_mod_poly_set_fmpz_poly(
                            factor_mod_prime.raw(), residue_poly_.raw(),
                            prime_ctx.raw());
                    if (fmpz_mod_poly_degree(factor_mod_prime.raw(),
                                             prime_ctx.raw()) > 0) {
                        const HenselValuationCache* cache =
                                cache_for_bound();
                        if (cache == nullptr) {
                            flint::Fmpz target_modulus;
                            fmpz_pow_ui(target_modulus.raw(), p_.raw(),
                                        static_cast<ulong>(bound));
                            std::unique_ptr<HenselValuationCache> next(
                                    new (std::nothrow) HenselValuationCache(
                                            target_modulus.raw(), bound));
                            if (next != nullptr && next->is_initialized() &&
                                hensel_lift_factor_to_precision(
                                        next->lifted_factor, raw_field->pol,
                                        field_mod_prime, factor_mod_prime,
                                        prime_ctx, p_.raw(), bound,
                                        next->target_modulus.raw(),
                                        next->target_ctx) &&
                                fmpz_mod_poly_degree(
                                        next->lifted_factor.raw(),
                                        next->target_ctx.raw()) > 0) {
                                hensel_valuation_cache_ = std::move(next);
                                cache = hensel_valuation_cache_.get();
                            }
                        }

                        slong valuation = 0;
                        if (try_cached_hensel(cache, valuation)) {
                            SILEX_PRIME_IDEAL_PROFILE_EVENT(
                                    diagnostics,
                                    "prime_ideal.valuation_order_element.path_no_index_hensel");
                            out = valuation;
                            return true;
                        }
                    }
                }
            }
        }
    }

    const flint::FmpzMat* multiplier =
            coordinate_valuation_matrix_cached();
    if (multiplier != nullptr &&
        valuation_by_coordinate_divisibility(
                out, parent_, element,
                flint::FmpzMatConstRef(*multiplier),
                flint::FmpzConstRef(p_), bound, diagnostics)) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_order_element.path_coordinate_divisibility");
        return true;
    }

    SILEX_PRIME_IDEAL_PROFILE_EVENT(
            diagnostics,
            "prime_ideal.valuation_order_element.path_containment_fallback");
    return valuation_by_power_containment(out, element, bound, diagnostics);
}

bool PrimeIdeal::valuation(slong& out, const Element& element) const noexcept {
    return valuation(out, element, nullptr);
}

bool PrimeIdeal::valuation(
        slong& out,
        const Element& element,
        const DiagnosticsContext* diagnostics) const noexcept {
    // reference `valuation(a, p)` clears an order denominator and subtracts its
    // prime valuation. reference `nfval` makes the same primitive-part/content
    // split. Reuse the exact principal fractional-ideal representation here.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_element");
    if (!has_prime_data() || parent_.parent() == nullptr ||
        !element.has_parent(*parent_.parent()) || element.equal_si(0)) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics, "prime_ideal.valuation_element.error_input");
        return false;
    }

    FractionalIdeal principal(parent_);
    slong candidate = 0;
    if (!principal.is_defined() ||
        !principal.set_principal(element, diagnostics) ||
        !valuation(candidate, principal, diagnostics)) {
        return false;
    }

    out = candidate;
    return true;
}

bool PrimeIdeal::valuation(slong& out, const Ideal& ideal) const noexcept {
    return valuation(out, ideal, nullptr);
}

bool PrimeIdeal::valuation(slong& out,
                           const Ideal& ideal,
                           const DiagnosticsContext* diagnostics) const noexcept {
    // Source trace: reference `NfOrd/Ideal/Valuation.jl:valuation(A, p)` first
    // returns valuation 0 when the ideal norm/minimum is prime-to-p, then uses
    // a known principal generator when available before falling back to ideal
    // containment/power tests for arbitrary ideals.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_ideal");
    if (!has_prime_data() || !same_order_parent(ideal.parent(), &parent_) ||
        !ideal.has_hnf() || !parent_.is_maximal() ||
        ideal_.is_one() || f_ <= 0) {
        return false;
    }

    flint::Fmpz norm;
    if (!ideal.norm(flint::FmpzRef(norm))) {
        return false;
    }

    slong bound = 0;
    while (fmpz_divisible(norm.raw(), p_.raw()) != 0) {
        if (bound == WORD_MAX) {
            return false;
        }
        fmpz_divexact(norm.raw(), norm.raw(), p_.raw());
        ++bound;
    }

    bound /= f_;
    if (bound == 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics, "prime_ideal.valuation_ideal.path_norm_zero");
        out = 0;
        return true;
    }

    if (const OrderElement* generator = ideal.known_principal_generator();
        generator != nullptr) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_ideal.path_principal_generator");
        return valuation(out, *generator, diagnostics);
    }

    SILEX_PRIME_IDEAL_PROFILE_EVENT(
            diagnostics, "prime_ideal.valuation_ideal.path_containment");
    slong value = 0;
    if (ideal_.contains(ideal)) {
        value = 1;

        if (value < bound) {
            Ideal power(parent_);
            Ideal next(parent_);
            if (!power.is_defined() || !next.is_defined() ||
                !power.set(ideal_)) {
                return false;
            }

            while (value < bound) {
                if (!next.multiply(power, ideal_) ||
                    !next.contains(ideal)) {
                    break;
                }

                ++value;
                power.swap(next);
            }
        }
    }

    out = value;
    return true;
}

bool PrimeIdeal::valuation(slong& out,
                           const FractionalIdeal& ideal) const noexcept {
    return valuation(out, ideal, nullptr);
}

bool PrimeIdeal::valuation(
        slong& out,
        const FractionalIdeal& ideal,
        const DiagnosticsContext* diagnostics) const noexcept {
    // reference `valuation(A::FractionalIdeal, p)` is numerator valuation minus
    // denominator valuation. For the stored positive integer denominator d,
    // v_P(d O) = e(P/p) * v_p(d), matching reference `nfval` denominator content.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_fractional_ideal");
    if (!has_prime_data() || !same_order_parent(ideal.parent(), &parent_) ||
        !ideal.has_integral_denominator() || !parent_.is_maximal() || e_ <= 0 ||
        f_ <= 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_fractional_ideal.error_input");
        return false;
    }

    Ideal numerator(parent_);
    flint::Fmpz denominator;
    if (!numerator.is_defined() ||
        !ideal.get_integral_den(numerator, flint::FmpzRef(denominator))) {
        return false;
    }

    slong numerator_valuation = 0;
    if (!valuation(numerator_valuation, numerator, diagnostics)) {
        return false;
    }

    flint::Fmpz remaining_denominator;
    flint::fmpz_set(flint::FmpzRef(remaining_denominator),
                    flint::FmpzConstRef(denominator));
    slong denominator_vp = 0;
    while (fmpz_divisible(remaining_denominator.raw(), p_.raw()) != 0) {
        if (denominator_vp == WORD_MAX) {
            SILEX_PRIME_IDEAL_PROFILE_EVENT(
                    diagnostics,
                    "prime_ideal.valuation_fractional_ideal.error_denominator_vp_overflow");
            return false;
        }
        fmpz_divexact(remaining_denominator.raw(), remaining_denominator.raw(),
                      p_.raw());
        ++denominator_vp;
    }

    flint::Fmpz candidate;
    flint::Fmpz ramification;
    flint::Fmpz denominator_exponent;
    flint::fmpz_set_si(flint::FmpzRef(candidate), numerator_valuation);
    flint::fmpz_set_si(flint::FmpzRef(ramification), e_);
    flint::fmpz_set_si(flint::FmpzRef(denominator_exponent), denominator_vp);
    flint::fmpz_submul(flint::FmpzRef(candidate),
                       flint::FmpzConstRef(ramification),
                       flint::FmpzConstRef(denominator_exponent));
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(candidate))) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_fractional_ideal.error_result_overflow");
        return false;
    }

    out = flint::fmpz_get_si(flint::FmpzConstRef(candidate));
    return true;
}

bool PrimeIdeal::valuation(slong& out,
                           const FactoredElement& element) const noexcept {
    return valuation(out, element, nullptr);
}

bool PrimeIdeal::valuation(
        slong& out,
        const FactoredElement& element,
        const DiagnosticsContext* diagnostics) const noexcept {
    // reference `famat_nfvalrem` and reference `valuation(::FacElem, P)` both sum
    // factor valuations times signed exponents without expanding the product.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.valuation_factored_element");
    if (!has_prime_data() || parent_.parent() == nullptr ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*parent_.parent()) ||
        !parent_.is_maximal() || e_ <= 0 || f_ <= 0) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_factored_element.error_input");
        return false;
    }

    flint::Fmpz candidate;
    flint::fmpz_zero(flint::FmpzRef(candidate));
    for (const auto& entry : element.factors()) {
        if (entry.exponent == 0) {
            continue;
        }

        slong factor_valuation = 0;
        if (!valuation(factor_valuation, entry.factor, diagnostics)) {
            return false;
        }

        flint::Fmpz factor_valuation_z;
        flint::Fmpz exponent_z;
        flint::fmpz_set_si(flint::FmpzRef(factor_valuation_z),
                           factor_valuation);
        flint::fmpz_set_si(flint::FmpzRef(exponent_z), entry.exponent);
        flint::fmpz_addmul(flint::FmpzRef(candidate),
                           flint::FmpzConstRef(factor_valuation_z),
                           flint::FmpzConstRef(exponent_z));
    }

    if (!flint::fmpz_fits_si(flint::FmpzConstRef(candidate))) {
        SILEX_PRIME_IDEAL_PROFILE_EVENT(
                diagnostics,
                "prime_ideal.valuation_factored_element.error_result_overflow");
        return false;
    }

    out = flint::fmpz_get_si(flint::FmpzConstRef(candidate));
    return true;
}

bool PrimeIdeal::set_data(
        flint::FmpzConstRef p,
        slong ramification_index,
        slong residue_degree,
        const Ideal& ideal,
        flint::FmpzPolyConstRef residue_polynomial) noexcept {
    if (!is_defined() || !same_order_parent(ideal.parent(), &parent_) ||
        !ideal.has_hnf() || ramification_index <= 0 || residue_degree <= 0 ||
        fmpz_poly_degree(residue_polynomial.raw()) <= 0) {
        return false;
    }

    Ideal ideal_copy(parent_);
    if (!ideal_copy.is_defined() || !ideal_copy.set(ideal)) {
        return false;
    }

    fmpz_set(p_.raw(), p.raw());
    ideal_.swap(ideal_copy);
    kummer_generator_ = flint::FmpzMat(0, 0);
    fmpz_poly_set(residue_poly_.raw(), residue_polynomial.raw());
    has_linear_residue_root_ = linear_root_from_residue_polynomial(
            linear_residue_root_, flint::FmpzPolyConstRef(residue_poly_),
            p);
    hensel_valuation_cache_.reset();
    containment_power_cache_.reset();
    coordinate_valuation_cache_.reset();
    e_ = ramification_index;
    f_ = residue_degree;
    has_prime_ = true;
    has_kummer_generator_ = false;
    has_residue_poly_ = true;
    return true;
}

bool PrimeIdeal::set_data_no_residue(flint::FmpzConstRef p,
                                     slong ramification_index,
                                     slong residue_degree,
                                     const Ideal& ideal) noexcept {
    if (!is_defined() || !same_order_parent(ideal.parent(), &parent_) ||
        !ideal.has_hnf() || ramification_index <= 0 || residue_degree <= 0) {
        return false;
    }

    Ideal ideal_copy(parent_);
    if (!ideal_copy.is_defined() || !ideal_copy.set(ideal)) {
        return false;
    }

    fmpz_set(p_.raw(), p.raw());
    ideal_.swap(ideal_copy);
    kummer_generator_ = flint::FmpzMat(0, 0);
    fmpz_poly_zero(residue_poly_.raw());
    fmpz_zero(linear_residue_root_.raw());
    hensel_valuation_cache_.reset();
    containment_power_cache_.reset();
    coordinate_valuation_cache_.reset();
    e_ = ramification_index;
    f_ = residue_degree;
    has_prime_ = true;
    has_kummer_generator_ = false;
    has_residue_poly_ = false;
    has_linear_residue_root_ = false;
    return true;
}

bool PrimeIdeal::set_kummer_generator(
        flint::FmpzMatConstRef generator) noexcept {
    if (!has_prime_data() || fmpz_mat_nrows(generator.raw()) != 1 ||
        fmpz_mat_ncols(generator.raw()) != parent_.degree()) {
        return false;
    }

    flint::FmpzMat copy(1, parent_.degree());
    fmpz_mat_set(copy.raw(), generator.raw());
    kummer_generator_.swap(copy);
    has_kummer_generator_ = true;
    return true;
}

namespace detail {

const flint::FmpzPoly* residue_polynomial_ptr(
        const PrimeIdeal& prime) noexcept {
    return prime.has_prime_data() && prime.has_residue_poly_
            ? &prime.residue_poly_
            : nullptr;
}

const flint::Fmpz* linear_residue_root_ptr(
        const PrimeIdeal& prime) noexcept {
    return prime.has_prime_data() && prime.has_linear_residue_root_
            ? &prime.linear_residue_root_
            : nullptr;
}

bool MaximalQuadraticPrimeAccess::set_from_integral_generator_factor(
        PrimeIdeal& out,
        const Order& order,
        flint::FmpzConstRef p,
        flint::FmpzModPolyConstRef factor,
        slong ramification_index,
        const Element& integral_generator,
        const flint::FmpzModCtx& context) noexcept {
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        order.parent() == nullptr || factor.raw() == nullptr ||
        context.raw() == nullptr || ramification_index <= 0) {
        return false;
    }

    const slong residue_degree =
            fmpz_mod_poly_degree(factor.raw(), context.raw());
    Element generator_element(*order.parent());
    OrderElement generator(order);
    Ideal candidate_ideal(order);
    flint::FmpzMat kummer_generator(1, order.degree());
    flint::FmpzPoly residue_polynomial;
    if (residue_degree <= 0 || !generator_element.is_defined() ||
        !generator.is_defined() || !candidate_ideal.is_defined() ||
        !element_evaluate_mod_poly_at(generator_element, factor.raw(),
                                      integral_generator, context) ||
        !generator.set_element(generator_element) ||
        !generator.get_coordinates(flint::FmpzMatRef(kummer_generator)) ||
        !detail::set_known_two_generator_ideal(candidate_ideal, p,
                                               generator)) {
        return false;
    }

    fmpz_mat_center_row_mod_prime(kummer_generator, 0, p);
    fmpz_mod_poly_get_fmpz_poly(residue_polynomial, factor.raw(), context);
    return out.set_data(p, ramification_index, residue_degree,
                        candidate_ideal,
                        flint::FmpzPolyConstRef(residue_polynomial)) &&
           out.set_kummer_generator(
                   flint::FmpzMatConstRef(kummer_generator));
}

bool MaximalQuadraticPrimeAccess::set_first_degree_one_prime(
        PrimeIdeal& out,
        RetainedQuadraticPrimeKind& kind,
        const Order& order,
        flint::FmpzConstRef p,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::prime_ideal,
            "prime_ideal.maximal_quadratic.first_degree_one");
    kind = RetainedQuadraticPrimeKind::inert;
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        !order.has_basis() || order.parent() == nullptr ||
        order.degree() != 2 || !order.is_maximal() ||
        order.parent()->backend_kind() !=
                NumberFieldBackendKind::quadratic ||
        fmpz_is_prime(p.raw()) == 0 ||
        !fmpq_poly_is_monic_integral(order.parent()->raw_flint_field()->pol)) {
        return false;
    }

    flint::Fmpz radicand;
    flint::Fmpz conductor;
    Element integral_generator(*order.parent());
    flint::FmpzPoly minimal_polynomial;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::prime_ideal,
                "prime_ideal.maximal_quadratic.prepare");
        if (!order.parent()->quadratic_radicand(flint::FmpzRef(radicand)) ||
            !order.quadratic_conductor(flint::FmpzRef(conductor)) ||
            !flint::fmpz_is_one(flint::FmpzConstRef(conductor)) ||
            !integral_generator.is_defined() ||
            !set_quadratic_integral_generator(
                    integral_generator, *order.parent(),
                    flint::FmpzConstRef(radicand))) {
            return false;
        }
        quadratic_integral_generator_minimal_polynomial(
                minimal_polynomial, flint::FmpzConstRef(radicand));
    }

    flint::FmpzModCtx context(p.raw());
    if (context.raw() == nullptr) {
        return false;
    }
    flint::FmpzModPoly reduced(context);
    flint::FmpzModPolyFactor factorization(context);
    if (!reduced.is_initialized()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::prime_ideal,
                "prime_ideal.maximal_quadratic.factor");
        fmpz_mod_poly_set_fmpz_poly(reduced.raw(), minimal_polynomial.raw(),
                                    context.raw());
        fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), context.raw());
    }

    const slong count = factorization.raw()->num;
    const auto degree = [&](slong index) noexcept {
        return fmpz_mod_poly_degree(factorization.raw()->poly + index,
                                    context.raw());
    };
    if (count == 1 && degree(0) == 2 &&
        factorization.raw()->exp[0] == 1) {
        kind = RetainedQuadraticPrimeKind::inert;
        return true;
    }

    slong retained_ramification_index = 0;
    if (count == 1 && degree(0) == 1 &&
        factorization.raw()->exp[0] == 2) {
        kind = RetainedQuadraticPrimeKind::ramified;
        retained_ramification_index = 2;
    } else if (count == 2 && degree(0) == 1 && degree(1) == 1 &&
               factorization.raw()->exp[0] == 1 &&
               factorization.raw()->exp[1] == 1) {
        kind = RetainedQuadraticPrimeKind::split;
        retained_ramification_index = 1;
    } else {
        return false;
    }

    PrimeIdeal candidate(order);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::prime_ideal,
                "prime_ideal.maximal_quadratic.materialize");
        if (!candidate.is_defined() ||
            !set_from_integral_generator_factor(
                    candidate, order, p,
                    flint::FmpzModPolyConstRef(factorization.raw()->poly),
                    retained_ramification_index, integral_generator,
                    context)) {
            return false;
        }
    }
    out.swap(candidate);
    return true;
}

bool set_degree_one_prime_ideal_from_root(PrimeIdeal& out,
                                          const Order& order,
                                          flint::FmpzConstRef p,
                                          flint::FmpzConstRef root) noexcept {
    if (!out.is_defined() || !same_order_parent(out.parent(), &order) ||
        !order.has_basis() || order.parent() == nullptr ||
        (!order.is_equation_order() && !order.is_maximal()) ||
        fmpz_is_prime(p.raw()) == 0 ||
        !fmpq_poly_is_monic_integral(order.parent()->raw_flint_field()->pol)) {
        return false;
    }

    const NumberField* field = order.parent();
    flint::Fmpz radicand;
    flint::Fmpz conductor;
    const bool use_quadratic_integral_generator =
            order.degree() == 2 && order.is_maximal() &&
            field->backend_kind() == NumberFieldBackendKind::quadratic &&
            field->quadratic_radicand(flint::FmpzRef(radicand)) &&
            order.quadratic_conductor(flint::FmpzRef(conductor)) &&
            flint::fmpz_is_one(conductor);

    flint::Fmpz residue_root;
    flint::fmpz_set(flint::FmpzRef(residue_root), root);
    if (use_quadratic_integral_generator &&
        flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        flint::Fmpz two;
        flint::Fmpz inverse_two;
        flint::fmpz_set_ui(flint::FmpzRef(two), 2);
        if (!flint::fmpz_invmod(
                    flint::FmpzRef(inverse_two),
                    flint::FmpzConstRef(two), p)) {
            return false;
        }
        flint::fmpz_add_ui(
                flint::FmpzRef(residue_root), root, 1);
        flint::fmpz_mul(
                flint::FmpzRef(residue_root),
                flint::FmpzConstRef(residue_root),
                flint::FmpzConstRef(inverse_two));
        fmpz_mod(residue_root.raw(), residue_root.raw(), p.raw());
    }

    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzModPoly residue_mod(ctx);
    if (!residue_mod.is_initialized()) {
        return false;
    }

    flint::Fmpz neg_root;
    fmpz_mod(residue_root.raw(), residue_root.raw(), p.raw());
    fmpz_neg(neg_root.raw(), residue_root.raw());
    fmpz_mod(neg_root.raw(), neg_root.raw(), p.raw());
    fmpz_mod_poly_set_coeff_fmpz(residue_mod.raw(), 0, neg_root.raw(),
                                 ctx.raw());
    fmpz_mod_poly_set_coeff_ui(residue_mod.raw(), 1, 1, ctx.raw());

    flint::FmpzPoly residue_polynomial;
    fmpz_mod_poly_get_fmpz_poly(residue_polynomial, residue_mod.raw(), ctx);

    Ideal p_ideal(order);
    OrderElement p_element(order);
    Element integral_generator(*field);
    Element generator_element(*field);
    OrderElement generator(order);
    Ideal generator_ideal(order);
    Ideal candidate_ideal(order);
    flint::FmpzMat kummer_generator(1, order.degree());
    if (!p_ideal.is_defined() || !p_element.is_defined() ||
        !integral_generator.is_defined() || !generator_element.is_defined() ||
        !generator.is_defined() || !generator_ideal.is_defined() ||
        !candidate_ideal.is_defined()) {
        return false;
    }

    bool generator_set = false;
    if (use_quadratic_integral_generator) {
        generator_set = set_quadratic_integral_generator(
                                integral_generator, *field,
                                flint::FmpzConstRef(radicand)) &&
                        element_evaluate_mod_poly_at(
                                generator_element, residue_mod.raw(),
                                integral_generator, ctx);
    } else {
        generator_set = element_set_mod_poly(
                generator_element, residue_mod.raw(), ctx);
    }
    if (!generator_set || !order_element_set_fmpz(p_element, order, p) ||
        !p_ideal.set_principal(p_element) ||
        !generator.set_element(generator_element) ||
        !generator.get_coordinates(flint::FmpzMatRef(kummer_generator)) ||
        !generator_ideal.set_principal(generator) ||
        !candidate_ideal.add(p_ideal, generator_ideal)) {
        return false;
    }

    fmpz_mat_center_row_mod_prime(kummer_generator, 0, p);
    return out.set_data(p, 1, 1, candidate_ideal,
                        flint::FmpzPolyConstRef(residue_polynomial)) &&
           out.set_kummer_generator(flint::FmpzMatConstRef(kummer_generator));
}

}  // namespace detail

PrimeIdealList::PrimeIdealList(const Order& parent, slong length) noexcept {
    define(parent, length);
}

PrimeIdealList::~PrimeIdealList() noexcept {
    clear();
}

PrimeIdealList::PrimeIdealList(PrimeIdealList&& other) noexcept {
    swap(other);
}

PrimeIdealList& PrimeIdealList::operator=(PrimeIdealList&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void PrimeIdealList::swap(PrimeIdealList& other) noexcept {
    entries_.swap(other.entries_);
    std::swap(size_, other.size_);
    std::swap(capacity_, other.capacity_);
    std::swap(defined_, other.defined_);
}

void PrimeIdealList::clear() noexcept {
    entries_.reset();
    size_ = 0;
    capacity_ = 0;
    defined_ = false;
}

bool PrimeIdealList::define(const Order& parent, slong length) noexcept {
    if (!parent.has_basis() || length < 0) {
        return false;
    }

    PrimeIdealList next;
    next.defined_ = true;

    if (!next.reserve(length)) {
        return false;
    }
    for (slong i = 0; i < length; ++i) {
        if (!next.append(parent)) {
            return false;
        }
    }

    swap(next);
    return true;
}

bool PrimeIdealList::is_defined() const noexcept {
    return defined_;
}

slong PrimeIdealList::size() const noexcept {
    return defined_ ? size_ : 0;
}

PrimeIdeal* PrimeIdealList::at(slong index) noexcept {
    if (!defined_ || index < 0 || index >= size_) {
        return nullptr;
    }
    return &entries_[index];
}

const PrimeIdeal* PrimeIdealList::at(slong index) const noexcept {
    if (!defined_ || index < 0 || index >= size_) {
        return nullptr;
    }
    return &entries_[index];
}

bool PrimeIdealList::append(const Order& parent) noexcept {
    if (!defined_ || !parent.has_basis()) {
        return false;
    }
    if (size_ == capacity_) {
        const slong next_capacity = capacity_ == 0 ? 1 : 2 * capacity_;
        if (!reserve(next_capacity)) {
            return false;
        }
    }
    entries_[size_] = PrimeIdeal(parent);
    if (!entries_[size_].is_defined()) {
        return false;
    }
    ++size_;
    return true;
}

PrimeIdeal* PrimeIdealList::back() noexcept {
    if (!defined_ || size_ == 0) {
        return nullptr;
    }
    return &entries_[size_ - 1];
}

bool PrimeIdealList::reserve(slong capacity) noexcept {
    if (capacity < 0) {
        return false;
    }
    if (capacity <= capacity_) {
        return true;
    }

    std::unique_ptr<PrimeIdeal[]> next(new (std::nothrow) PrimeIdeal[capacity]);
    if (!next) {
        return false;
    }
    for (slong i = 0; i < size_; ++i) {
        next[i] = std::move(entries_[i]);
    }
    entries_.swap(next);
    capacity_ = capacity;
    return true;
}

bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p) noexcept {
    return decompose_prime(out, order, p, 0, nullptr);
}

bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p,
                     slong max_residue_degree) noexcept {
    return decompose_prime(out, order, p, max_residue_degree, nullptr);
}

bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p,
                     slong max_residue_degree,
                     const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.decompose_prime");
    if (!order.has_basis() || order.parent() == nullptr ||
        max_residue_degree < 0 ||
        fmpz_is_prime(p.raw()) == 0 ||
        !fmpq_poly_is_monic_integral(order.parent()->raw_flint_field()->pol)) {
        return false;
    }

    PrimeIdealList candidates;
    candidates.defined_ = true;

    if (order.degree() == 1) {
        Ideal p_ideal(order);
        OrderElement p_element(order);
        flint::FmpzPoly residue_polynomial;
        fmpz_poly_set_coeff_si(residue_polynomial.raw(), 1, 1);

        if (!p_ideal.is_defined() || !p_element.is_defined() ||
            !order_element_set_fmpz(p_element, order, p) ||
            !p_ideal.set_principal(p_element)) {
            return false;
        }

        if (residue_degree_allowed(1, max_residue_degree)) {
            if (!candidates.append(order) ||
                !candidates.back()->set_data(
                        p, 1, 1, p_ideal,
                        flint::FmpzPolyConstRef(residue_polynomial))) {
                return false;
            }
        }

        out.swap(candidates);
        return true;
    }

    flint::FmpzModCtx ctx(p.raw());

    flint::Fmpz quadratic_radicand;
    flint::Fmpz quadratic_conductor;
    if (order.degree() == 2 && order.is_maximal() &&
            order.parent()->backend_kind() ==
                    NumberFieldBackendKind::quadratic &&
            order.parent()->quadratic_radicand(
                    flint::FmpzRef(quadratic_radicand)) &&
            order.quadratic_conductor(flint::FmpzRef(quadratic_conductor)) &&
            flint::fmpz_is_one(quadratic_conductor)) {
        // reference `quadgen`/`quadpoly` use the integral generator omega and its
        // minimal polynomial.  Factoring in that generator keeps the stored
        // residue polynomial aligned with this order's natural basis.
        Element omega(*order.parent());
        flint::FmpzPoly minimal_polynomial;
        if (!set_quadratic_integral_generator(
                    omega, *order.parent(),
                    flint::FmpzConstRef(quadratic_radicand))) {
            return false;
        }
        quadratic_integral_generator_minimal_polynomial(
                minimal_polynomial, flint::FmpzConstRef(quadratic_radicand));

        flint::FmpzModPoly reduced_quadratic(ctx);
        flint::FmpzModPolyFactor quadratic_factorization(ctx);
        fmpz_mod_poly_set_fmpz_poly(
                reduced_quadratic.raw(), minimal_polynomial.raw(), ctx.raw());
        fmpz_mod_poly_factor(quadratic_factorization.raw(),
                reduced_quadratic.raw(),
                ctx.raw());

        const slong factor_count = quadratic_factorization.raw()->num;
        if (factor_count <= 0 || !candidates.reserve(factor_count)) {
            return false;
        }
        for (slong i = 0; i < factor_count; ++i) {
            const slong residue_degree = fmpz_mod_poly_degree(
                    quadratic_factorization.raw()->poly + i, ctx.raw());
            if (!residue_degree_allowed(residue_degree, max_residue_degree)) {
                continue;
            }
            if (!candidates.append(order) ||
                !detail::MaximalQuadraticPrimeAccess::
                        set_from_integral_generator_factor(
                                *candidates.back(), order, p,
                                flint::FmpzModPolyConstRef(
                                        quadratic_factorization.raw()->poly +
                                        i),
                                quadratic_factorization.raw()->exp[i], omega,
                                ctx)) {
                return false;
            }
        }

        out.swap(candidates);
        return true;
    }

    flint::FmpzModPoly reduced_field_polynomial(ctx);
    fmpq_poly_get_fmpz_mod_poly(reduced_field_polynomial,
                                order.parent()->raw_flint_field()->pol, ctx);
    auto append_defining_polynomial_factor =
            [&](const flint::FmpzModPolyFactor& factorization,
                slong factor_index) noexcept -> bool {
        const slong residue_degree =
                fmpz_mod_poly_degree(factorization.raw()->poly + factor_index,
                                     ctx.raw());
        if (!residue_degree_allowed(residue_degree, max_residue_degree)) {
            return true;
        }

        Element generator_element(*order.parent());
        OrderElement generator(order);
        Ideal candidate_ideal(order);
        flint::FmpzMat kummer_generator(1, order.degree());
        flint::FmpzPoly residue_polynomial;

        if (!generator_element.is_defined() || !generator.is_defined() ||
            !candidate_ideal.is_defined()) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                                "prime_ideal.decompose_prime.factor_to_generator");
            if (!element_set_mod_poly(
                        generator_element,
                        factorization.raw()->poly + factor_index,
                        ctx) ||
                !generator.set_element(generator_element) ||
                !generator.get_coordinates(
                        flint::FmpzMatRef(kummer_generator))) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::prime_ideal,
                                "prime_ideal.decompose_prime.two_generator_ideal");
            if (!detail::set_known_two_generator_ideal(
                        candidate_ideal, p, generator)) {
                return false;
            }
        }
        fmpz_mat_center_row_mod_prime(kummer_generator, 0, p);

        fmpz_mod_poly_get_fmpz_poly(residue_polynomial,
                                    factorization.raw()->poly + factor_index,
                                    ctx);
        return candidates.append(order) &&
               candidates.back()->set_data(
                       p, factorization.raw()->exp[factor_index],
                       residue_degree, candidate_ideal,
                       flint::FmpzPolyConstRef(residue_polynomial)) &&
               candidates.back()->set_kummer_generator(
                       flint::FmpzMatConstRef(kummer_generator));
    };

    if (max_residue_degree == 1 &&
        (order.is_equation_order() || order.is_maximal()) &&
        fmpz_mod_poly_is_squarefree(reduced_field_polynomial.raw(),
                                    ctx.raw()) != 0) {
        flint::FmpzModPolyFactor roots(ctx);
        fmpz_mod_poly_roots(roots.raw(), reduced_field_polynomial.raw(), 0,
                            ctx.raw());
        if (!candidates.reserve(roots.raw()->num)) {
            return false;
        }
        for (slong i = 0; i < roots.raw()->num; ++i) {
            const slong residue_degree =
                    fmpz_mod_poly_degree(roots.raw()->poly + i, ctx.raw());
            if (residue_degree != 1) {
                continue;
            }

            Element generator_element(*order.parent());
            OrderElement generator(order);
            Ideal candidate_ideal(order);
            flint::FmpzMat kummer_generator(1, order.degree());
            flint::FmpzPoly residue_polynomial;

            if (!generator_element.is_defined() ||
                !generator.is_defined() ||
                !candidate_ideal.is_defined() ||
                !element_set_mod_poly(generator_element,
                                      roots.raw()->poly + i, ctx) ||
                !generator.set_element(generator_element) ||
                !generator.get_coordinates(
                        flint::FmpzMatRef(kummer_generator)) ||
                !detail::set_known_two_generator_ideal(
                        candidate_ideal, p, generator)) {
                return false;
            }
            fmpz_mat_center_row_mod_prime(kummer_generator, 0, p);
            fmpz_mod_poly_get_fmpz_poly(residue_polynomial,
                                        roots.raw()->poly + i, ctx);
            if (!candidates.append(order) ||
                !candidates.back()->set_data(
                        p, roots.raw()->exp[i], residue_degree,
                        candidate_ideal,
                        flint::FmpzPolyConstRef(residue_polynomial)) ||
                !candidates.back()->set_kummer_generator(
                        flint::FmpzMatConstRef(kummer_generator))) {
                return false;
            }
        }

        out.swap(candidates);
        return true;
    }

    if (!order.is_equation_order()) {
        if (!order.is_maximal()) {
            return false;
        }

        if (maximal_order_has_nonindex_defining_polynomial_prime(order, p)) {
            flint::FmpzModPoly reduced(ctx);
            flint::FmpzModPolyFactor factorization(ctx);

            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::prime_ideal,
                        "prime_ideal.decompose_prime.nonindex_factor_field_polynomial");
                fmpz_mod_poly_set(reduced.raw(), reduced_field_polynomial.raw(),
                                  ctx.raw());
                fmpz_mod_poly_factor(factorization.raw(), reduced.raw(),
                                     ctx.raw());
            }

            const slong num_factors = factorization.raw()->num;
            if (num_factors <= 0 || !candidates.reserve(num_factors)) {
                return false;
            }
            for (slong i = 0; i < num_factors; ++i) {
                if (!append_defining_polynomial_factor(factorization, i)) {
                    return false;
                }
            }

            out.swap(candidates);
            return true;
        }

        std::vector<IndexPrimeData> unramified_primes;
        if (unramified_maximal_index_divisor_prime_data(
                    unramified_primes, order, p)) {
            if (!candidates.reserve(
                        static_cast<slong>(unramified_primes.size()))) {
                return false;
            }
            for (IndexPrimeData& prime : unramified_primes) {
                if (!residue_degree_allowed(prime.residue_degree,
                                            max_residue_degree)) {
                    continue;
                }
                if (!candidates.append(order) ||
                    !candidates.back()->set_data_no_residue(
                            p, prime.ramification_index,
                            prime.residue_degree, prime.ideal)) {
                    return false;
                }
            }
            out.swap(candidates);
            return true;
        }

        std::vector<IndexPrimeData> pradical_primes;
        if (ramified_maximal_index_divisor_prime_data(
                    pradical_primes, order, p)) {
            if (!candidates.reserve(
                        static_cast<slong>(pradical_primes.size()))) {
                return false;
            }
            for (IndexPrimeData& prime : pradical_primes) {
                if (!residue_degree_allowed(prime.residue_degree,
                                            max_residue_degree)) {
                    continue;
                }
                if (!candidates.append(order) ||
                    !candidates.back()->set_data_no_residue(
                            p, prime.ramification_index,
                            prime.residue_degree, prime.ideal)) {
                    return false;
                }
            }
            out.swap(candidates);
            return true;
        }

        Element beta(*order.parent());
        flint::FmpzPoly chi;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::prime_ideal,
                    "prime_ideal.decompose_prime.find_monogenic_generator");
            if (!beta.is_defined() ||
                !find_monogenic_generator(chi, beta, order)) {
                return false;
            }
        }

        flint::FmpzModPoly reduced(ctx);
        flint::FmpzModPolyFactor factorization(ctx);
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::prime_ideal,
                    "prime_ideal.decompose_prime.monogenic_factor_polynomial");
            fmpz_mod_poly_set_fmpz_poly(reduced.raw(), chi.raw(), ctx.raw());
            fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), ctx.raw());
        }
        if (factorization.raw()->num <= 0) {
            return false;
        }

        if (!candidates.reserve(factorization.raw()->num)) {
            return false;
        }
        for (slong i = 0; i < factorization.raw()->num; ++i) {
            const slong residue_degree =
                    fmpz_mod_poly_degree(factorization.raw()->poly + i,
                                         ctx.raw());
            if (!residue_degree_allowed(residue_degree,
                                        max_residue_degree)) {
                continue;
            }

            Element generator_element(*order.parent());
            OrderElement generator(order);
            Ideal candidate_ideal(order);
            if (!generator_element.is_defined() || !generator.is_defined() ||
                !candidate_ideal.is_defined() ||
                !element_evaluate_mod_poly_at(generator_element,
                                              factorization.raw()->poly + i,
                                              beta, ctx) ||
                !generator.set_element(generator_element) ||
                !detail::set_known_two_generator_ideal(
                        candidate_ideal, p, generator) ||
                !candidates.append(order) ||
                !candidates.back()->set_data_no_residue(
                        p, factorization.raw()->exp[i],
                        residue_degree, candidate_ideal)) {
                return false;
            }
        }

        out.swap(candidates);
        return true;
    }

    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPolyFactor factorization(ctx);

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::prime_ideal,
                "prime_ideal.decompose_prime.equation_factor_field_polynomial");
        fmpz_mod_poly_set(reduced.raw(), reduced_field_polynomial.raw(),
                          ctx.raw());
        fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), ctx.raw());
    }
    if (!order.is_maximal() && !factorization_is_squarefree(factorization)) {
        return false;
    }

    const slong num_factors = factorization.raw()->num;
    if (!candidates.reserve(num_factors)) {
        return false;
    }
    for (slong i = 0; i < num_factors; ++i) {
        if (!append_defining_polynomial_factor(factorization, i)) {
            return false;
        }
    }

    out.swap(candidates);
    return true;
}

}  // namespace silex
