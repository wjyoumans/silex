#include <silex/class_group.hpp>

#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/qfb.hpp>
#include <silex/signature.hpp>

#include "class_group_internal.hpp"
#include "class_group_certification_internal.hpp"
#include "factor_base_honesty_internal.hpp"
#include "../factor_base/factor_base_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"
#include "relation_search_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace silex {
namespace {

constexpr ulong kQuadraticRelationPhase = UWORD(11);
constexpr ulong kQuadraticHonestyPhase = UWORD(7);
constexpr ulong kQuadraticRelationModulePrime = UWORD(27449);
constexpr slong kQuadraticRandomExponentMax = 15;
constexpr slong kQuadraticPrincipalReductionPrecision = 200;

struct QuadraticFormPrime {
    flint::Fmpz rational_prime;
    flint::Qfb form;
    Ideal conjugate_ideal;
    std::vector<flint::Qfb> powers;
    slong oriented_factor_base_index = -1;
    slong conjugate_factor_base_index = -1;
    slong factor_base_orientation = 1;
    bool ramified = false;
};

struct QuadraticFormFactor {
    slong form_prime_index = -1;
    slong exponent = 0;
};

struct QuadraticFormFactorization {
    std::vector<QuadraticFormFactor> factors;
    ulong residual = 0;
    bool smooth_or_partial = false;
};

struct QuadraticLargePrimePartial {
    ulong residual = 0;
    std::vector<slong> exponents;
    slong pivot = -1;
};

bool element_scalar_div_fmpz(Element& out,
                             const Element& input,
                             flint::FmpzConstRef denominator) noexcept {
    if (!out.has_same_parent(input) || flint::fmpz_sgn(denominator) <= 0) {
        return false;
    }
    flint::FmpqPoly polynomial;
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    flint::fmpq_poly_scalar_div_fmpz(polynomial, polynomial, denominator);
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

void quadratic_norm_form_value(flint::Fmpz& out,
                               flint::FmpzConstRef x,
                               flint::FmpzConstRef y,
                               flint::FmpzConstRef omega_trace,
                               flint::FmpzConstRef omega_norm) noexcept {
    flint::Fmpz product;
    flint::fmpz_mul(flint::FmpzRef(out), x, x);
    flint::fmpz_mul(flint::FmpzRef(product), x, y);
    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(product), omega_trace);
    flint::fmpz_add(flint::FmpzRef(out), flint::FmpzConstRef(out),
                    flint::FmpzConstRef(product));
    flint::fmpz_mul(flint::FmpzRef(product), y, y);
    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(product), omega_norm);
    flint::fmpz_add(flint::FmpzRef(out), flint::FmpzConstRef(out),
                    flint::FmpzConstRef(product));
}

void quadratic_norm_form_cross_coefficient(
        flint::Fmpz& out,
        flint::FmpzConstRef x1,
        flint::FmpzConstRef y1,
        flint::FmpzConstRef x2,
        flint::FmpzConstRef y2,
        flint::FmpzConstRef omega_trace,
        flint::FmpzConstRef omega_norm) noexcept {
    flint::Fmpz product;
    flint::fmpz_mul(flint::FmpzRef(out), x1, x2);
    flint::fmpz_mul_2exp(flint::FmpzRef(out),
                         flint::FmpzConstRef(out), 1);

    flint::fmpz_mul(flint::FmpzRef(product), x1, y2);
    flint::fmpz_addmul(flint::FmpzRef(product), x2, y1);
    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(product), omega_trace);
    flint::fmpz_add(flint::FmpzRef(out), flint::FmpzConstRef(out),
                    flint::FmpzConstRef(product));

    flint::fmpz_mul(flint::FmpzRef(product), y1, y2);
    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(product), omega_norm);
    flint::fmpz_mul_2exp(flint::FmpzRef(product),
                         flint::FmpzConstRef(product), 1);
    flint::fmpz_add(flint::FmpzRef(out), flint::FmpzConstRef(out),
                    flint::FmpzConstRef(product));
}

bool quadratic_order_trace_and_norm(flint::Fmpz& omega_trace,
                                    flint::Fmpz& omega_norm,
                                    flint::FmpzConstRef discriminant) noexcept {
    if (flint::fmpz_sgn(discriminant) >= 0) {
        return false;
    }
    if (::fmpz_is_odd(discriminant.raw()) != 0) {
        flint::Fmpz numerator;
        flint::fmpz_one(flint::FmpzRef(omega_trace));
        flint::fmpz_neg(flint::FmpzRef(numerator), discriminant);
        flint::fmpz_add_ui(flint::FmpzRef(numerator),
                           flint::FmpzConstRef(numerator), 1);
        flint::fmpz_fdiv_q_2exp(flint::FmpzRef(omega_norm),
                                flint::FmpzConstRef(numerator), 2);
        return true;
    }
    flint::fmpz_zero(flint::FmpzRef(omega_trace));
    flint::fmpz_neg(flint::FmpzRef(omega_norm), discriminant);
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(omega_norm),
                            flint::FmpzConstRef(omega_norm), 2);
    return true;
}

void quadratic_form_reduce_basis_step(flint::Fmpz& b,
                                      flint::Fmpz& c,
                                      const flint::Fmpz& a,
                                      const flint::Fmpz& u1,
                                      flint::Fmpz& u2) noexcept {
    flint::Fmpz quotient;
    flint::Fmpz remainder;
    flint::Fmpz twice_a;
    flint::Fmpz half_sum;
    flint::Fmpz change;
    flint::fmpz_mul_2exp(flint::FmpzRef(twice_a),
                         flint::FmpzConstRef(a), 1);
    flint::fmpz_tdiv_qr(
            flint::FmpzRef(quotient), flint::FmpzRef(remainder),
            flint::FmpzConstRef(b), flint::FmpzConstRef(twice_a));
    const int remainder_comparison = flint::fmpz_cmpabs(
            flint::FmpzConstRef(remainder), flint::FmpzConstRef(a));
    if (flint::fmpz_sgn(flint::FmpzConstRef(b)) >= 0) {
        if (remainder_comparison > 0) {
            flint::fmpz_add_ui(flint::FmpzRef(quotient),
                               flint::FmpzConstRef(quotient), 1);
            flint::fmpz_sub(flint::FmpzRef(remainder),
                            flint::FmpzConstRef(remainder),
                            flint::FmpzConstRef(twice_a));
        }
    } else if (remainder_comparison >= 0) {
        flint::fmpz_sub_ui(flint::FmpzRef(quotient),
                           flint::FmpzConstRef(quotient), 1);
        flint::fmpz_add(flint::FmpzRef(remainder),
                        flint::FmpzConstRef(remainder),
                        flint::FmpzConstRef(twice_a));
    }
    flint::fmpz_add(flint::FmpzRef(half_sum), flint::FmpzConstRef(b),
                    flint::FmpzConstRef(remainder));
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(half_sum),
                            flint::FmpzConstRef(half_sum), 1);
    flint::fmpz_mul(flint::FmpzRef(change),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(half_sum));
    flint::fmpz_sub(flint::FmpzRef(c), flint::FmpzConstRef(c),
                    flint::FmpzConstRef(change));
    flint::fmpz_set(flint::FmpzRef(b), flint::FmpzConstRef(remainder));
    flint::fmpz_mul(flint::FmpzRef(change),
                    flint::FmpzConstRef(quotient),
                    flint::FmpzConstRef(u1));
    flint::fmpz_sub(flint::FmpzRef(u2), flint::FmpzConstRef(u2),
                    flint::FmpzConstRef(change));
}

bool quadratic_form_shortest_vector(flint::Fmpz& u1,
                                    flint::Fmpz& v1,
                                    flint::Fmpz& shortest_value,
                                    flint::FmpzConstRef input_a,
                                    flint::FmpzConstRef input_b,
                                    flint::FmpzConstRef input_c) noexcept {
    if (flint::fmpz_sgn(input_a) <= 0 || flint::fmpz_sgn(input_c) <= 0) {
        return false;
    }

    // Source: reference implementation `Qfb.c:qfbredsl2_imag_basecase`.  The first
    // column of its returned SL(2,Z) transform is a shortest vector for the
    // positive definite input form.
    flint::Fmpz a;
    flint::Fmpz b;
    flint::Fmpz c;
    flint::Fmpz u2;
    flint::fmpz_set(flint::FmpzRef(a), input_a);
    flint::fmpz_set(flint::FmpzRef(b), input_b);
    flint::fmpz_set(flint::FmpzRef(c), input_c);
    flint::fmpz_one(flint::FmpzRef(u1));
    flint::fmpz_zero(flint::FmpzRef(u2));

    int comparison = flint::fmpz_cmpabs(
            flint::FmpzConstRef(a), flint::FmpzConstRef(b));
    if (comparison < 0) {
        quadratic_form_reduce_basis_step(b, c, a, u1, u2);
    } else if (comparison == 0 && flint::fmpz_sgn(flint::FmpzConstRef(b)) < 0) {
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::fmpz_one(flint::FmpzRef(u2));
    }

    for (;;) {
        comparison = flint::fmpz_cmpabs(
                flint::FmpzConstRef(a), flint::FmpzConstRef(c));
        if (comparison <= 0) {
            break;
        }
        a.swap(c);
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::Fmpz old_u1;
        flint::fmpz_set(flint::FmpzRef(old_u1), flint::FmpzConstRef(u1));
        flint::fmpz_set(flint::FmpzRef(u1), flint::FmpzConstRef(u2));
        flint::fmpz_neg(flint::FmpzRef(u2), flint::FmpzConstRef(old_u1));
        quadratic_form_reduce_basis_step(b, c, a, u1, u2);
    }
    if (comparison == 0 && flint::fmpz_sgn(flint::FmpzConstRef(b)) < 0) {
        flint::fmpz_neg(flint::FmpzRef(b), flint::FmpzConstRef(b));
        flint::Fmpz old_u1;
        flint::fmpz_set(flint::FmpzRef(old_u1), flint::FmpzConstRef(u1));
        flint::fmpz_set(flint::FmpzRef(u1), flint::FmpzConstRef(u2));
        flint::fmpz_neg(flint::FmpzRef(u2), flint::FmpzConstRef(old_u1));
    }

    flint::Fmpz z;
    flint::Fmpz numerator;
    flint::Fmpz product;
    flint::fmpz_sub(flint::FmpzRef(z), flint::FmpzConstRef(b), input_b);
    if (flint::fmpz_fdiv_ui(flint::FmpzConstRef(z), 2) != 0) {
        return false;
    }
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(z), flint::FmpzConstRef(z), 1);
    flint::fmpz_mul(flint::FmpzRef(numerator), flint::FmpzConstRef(z),
                    flint::FmpzConstRef(u1));
    flint::fmpz_mul(flint::FmpzRef(product), flint::FmpzConstRef(a),
                    flint::FmpzConstRef(u2));
    flint::fmpz_sub(flint::FmpzRef(numerator),
                    flint::FmpzConstRef(numerator),
                    flint::FmpzConstRef(product));
    if (!flint::fmpz_divisible(flint::FmpzConstRef(numerator), input_c)) {
        return false;
    }
    flint::fmpz_divexact(flint::FmpzRef(v1),
                         flint::FmpzConstRef(numerator), input_c);
    flint::fmpz_set(flint::FmpzRef(shortest_value),
                    flint::FmpzConstRef(a));
    return true;
}

bool integral_ideal_pow_ui(Ideal& out,
                           const Ideal& input,
                           ulong exponent) noexcept {
    if (out.parent() == nullptr || input.parent() == nullptr ||
        !out.parent()->has_same_data(*input.parent())) {
        return false;
    }
    if (exponent == 1) {
        return out.set(input);
    }

    Ideal result(*input.parent());
    Ideal power(*input.parent());
    if (!result.is_defined() || !power.is_defined() || !result.one() ||
        !power.set(input)) {
        return false;
    }
    while (exponent != 0) {
        if ((exponent & 1U) != 0) {
            Ideal next(*input.parent());
            if (!next.is_defined() || !next.multiply(result, power)) {
                return false;
            }
            result.swap(next);
        }
        exponent >>= 1U;
        if (exponent != 0) {
            Ideal square(*input.parent());
            if (!square.is_defined() || !square.multiply(power, power)) {
                return false;
            }
            power.swap(square);
        }
    }
    out.swap(result);
    return true;
}

bool quadratic_factor_base_row_integral_den(
        Ideal& numerator,
        flint::Fmpz& denominator,
        const FactorBase& base,
        const std::vector<QuadraticFormPrime>& form_primes,
        flint::FmpzMatConstRef row,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || numerator.parent() == nullptr ||
        !numerator.parent()->has_same_data(*order) ||
        flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != base.length()) {
        return false;
    }

    Ideal accumulator(*order);
    PrimeIdeal prime(*order);
    Ideal prime_ideal(*order);
    Ideal prime_power(*order);
    if (!accumulator.is_defined() || !prime.is_defined() ||
        !prime_ideal.is_defined() || !prime_power.is_defined() ||
        !accumulator.one()) {
        return false;
    }
    flint::fmpz_one(flint::FmpzRef(denominator));

    auto multiply_prime_power =
            [&](slong index, slong conjugate_index,
                const Ideal* conjugate_ideal,
                flint::FmpzConstRef rational_prime) noexcept {
        const flint::FmpzConstRef exponent_ref =
                flint::fmpz_mat_entry(row, 0, index);
        if (flint::fmpz_is_zero(exponent_ref)) {
            return true;
        }
        if (!flint::fmpz_fits_si(exponent_ref)) {
            return false;
        }
        const slong signed_exponent = flint::fmpz_get_si(exponent_ref);
        if (signed_exponent == std::numeric_limits<slong>::min()) {
            return false;
        }
        const ulong exponent = static_cast<ulong>(
                signed_exponent < 0 ? -signed_exponent : signed_exponent);
        const slong selected_index =
                signed_exponent < 0 ? conjugate_index : index;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.maximal_imaginary_quadratic_relations."
                    "witness.row_ideal.prime");
            const bool selected_ideal_ok = signed_exponent < 0 &&
                                                   selected_index < 0
                    ? conjugate_ideal != nullptr &&
                              conjugate_ideal->is_defined() &&
                              prime_ideal.set(*conjugate_ideal)
                    : selected_index >= 0 &&
                              base.prime(prime, selected_index) &&
                              prime.get_ideal(prime_ideal);
            if (!selected_ideal_ok) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.maximal_imaginary_quadratic_relations."
                    "witness.row_ideal.power");
            if (!integral_ideal_pow_ui(
                        prime_power, prime_ideal, exponent)) {
                return false;
            }
        }
        Ideal next(*order);
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.maximal_imaginary_quadratic_relations."
                    "witness.row_ideal.multiply");
            if (!next.is_defined() ||
                !next.multiply(accumulator, prime_power)) {
                return false;
            }
        }
        accumulator.swap(next);
        if (signed_exponent < 0) {
            flint::Fmpz denominator_power;
            flint::fmpz_pow_ui(flint::FmpzRef(denominator_power),
                               rational_prime, exponent);
            flint::fmpz_mul(flint::FmpzRef(denominator),
                            flint::FmpzConstRef(denominator),
                            flint::FmpzConstRef(denominator_power));
        }
        return true;
    };

    // For split primes P * conjugate(P) = (p); for ramified primes P^2 =
    // (p).  Thus a negative row exponent is represented integrally by the
    // conjugate prime and a rational p denominator.
    for (const QuadraticFormPrime& form_prime : form_primes) {
        const slong conjugate = form_prime.ramified
                ? form_prime.oriented_factor_base_index
                : form_prime.conjugate_factor_base_index;
        if (!multiply_prime_power(
                    form_prime.oriented_factor_base_index, conjugate,
                    &form_prime.conjugate_ideal,
                    flint::FmpzConstRef(form_prime.rational_prime)) ||
            (!form_prime.ramified &&
             form_prime.conjugate_factor_base_index >= 0 &&
             !multiply_prime_power(
                     form_prime.conjugate_factor_base_index,
                     form_prime.oriented_factor_base_index,
                     nullptr,
                     flint::FmpzConstRef(form_prime.rational_prime)))) {
            return false;
        }
    }
    numerator.swap(accumulator);
    return true;
}

bool quadratic_integral_generator(Element& out,
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

bool quadratic_form_prime_shift(
        flint::Fmpz& out,
        flint::FmpzConstRef discriminant,
        flint::QfbConstRef form) noexcept {
    flint::fmpz_set(flint::FmpzRef(out), flint::qfb_b(form));
    if (::fmpz_is_odd(discriminant.raw()) != 0) {
        flint::fmpz_add_ui(flint::FmpzRef(out),
                           flint::FmpzConstRef(out), 1);
    }
    if (::fmpz_is_even(out.raw()) == 0) {
        return false;
    }
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(out),
                            flint::FmpzConstRef(out), 1);
    return true;
}

bool quadratic_form_prime_ideal(
        Ideal& out,
        const Order& order,
        const Element& omega,
        flint::FmpzConstRef discriminant,
        flint::FmpzConstRef rational_prime,
        flint::QfbConstRef form) noexcept {
    if (out.parent() == nullptr || !out.parent()->has_same_data(order) ||
        order.parent() == nullptr || !omega.has_parent(*order.parent())) {
        return false;
    }

    flint::Fmpz shift;
    if (!quadratic_form_prime_shift(shift, discriminant, form)) {
        return false;
    }

    Element shift_element(*order.parent());
    Element generator_element(*order.parent());
    OrderElement generator(order);
    if (!shift_element.is_defined() || !generator_element.is_defined() ||
        !generator.is_defined() ||
        !shift_element.set_fmpz(flint::FmpzConstRef(shift)) ||
        !generator_element.subtract(omega, shift_element) ||
        !generator.set_element(generator_element) ||
        !detail::set_known_two_generator_ideal(
                out, rational_prime, generator)) {
        return false;
    }
    return true;
}

bool quadratic_form_prime_orientation_matches(
        bool& matches,
        const PrimeIdeal& prime,
        flint::FmpzConstRef discriminant,
        flint::FmpzConstRef rational_prime,
        flint::QfbConstRef form) noexcept {
    matches = false;
    const flint::Fmpz* root = detail::linear_residue_root_ptr(prime);
    if (prime.residue_degree() != 1 || root == nullptr) {
        return false;
    }

    // reference `idealhnf_shallow` maps (p, b, c) to
    // p Z + (-b + sqrt(D))/2 Z.  With Silex's natural integral generator
    // omega this is (p, omega - shift), where shift is (b + 1)/2 for odd D
    // and b/2 for even D.  The quadratic decomposition path factors the
    // minimal polynomial of that same omega and caches omega mod p, so the
    // ideal orientations agree exactly when the cached root is shift mod p.
    flint::Fmpz shift;
    flint::Fmpz reduced_shift;
    if (!quadratic_form_prime_shift(shift, discriminant, form)) {
        return false;
    }
    ::fmpz_mod(reduced_shift.raw(), shift.raw(), rational_prime.raw());
    matches = flint::fmpz_equal(flint::FmpzConstRef(*root),
                                flint::FmpzConstRef(reduced_shift));
    return true;
}

bool identify_quadratic_form_primes(
        std::vector<QuadraticFormPrime>& out,
        const Order& order,
        const FactorBase& base,
        flint::FmpzConstRef discriminant,
        const DiagnosticsContext* diagnostics) noexcept {
    out.clear();
    if (order.parent() == nullptr || base.parent() == nullptr ||
        !base.parent()->has_same_data(order)) {
        return false;
    }

    flint::Fmpz radicand;
    Element omega(*order.parent());
    if (!order.parent()->quadratic_radicand(flint::FmpzRef(radicand)) ||
        !quadratic_integral_generator(
                omega, *order.parent(), flint::FmpzConstRef(radicand))) {
        return false;
    }

    out.reserve(static_cast<std::size_t>(base.rational_prime_block_count()));
    for (slong block = 0; block < base.rational_prime_block_count(); ++block) {
        flint::Fmpz rational_prime;
        slong start = 0;
        slong length = 0;
        bool complete = false;
        if (!detail::FactorBaseBlockAccess::
                    rational_prime_block_is_complete(
                            complete, base, block) ||
            !base.rational_prime_block(
                    flint::FmpzRef(rational_prime), start, length, block) ||
            (length != 1 && length != 2) ||
            !flint::fmpz_abs_fits_ui(
                    flint::FmpzConstRef(rational_prime))) {
            return false;
        }

        QuadraticFormPrime entry;
        flint::fmpz_set(flint::FmpzRef(entry.rational_prime),
                        flint::FmpzConstRef(rational_prime));
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.maximal_imaginary_quadratic_relations."
                    "form_setup.prime_form");
            flint::qfb_prime_form(
                    flint::QfbRef(entry.form), discriminant,
                    flint::FmpzConstRef(entry.rational_prime));
        }
        const int splitting = flint::fmpz_kronecker(
                discriminant,
                flint::FmpzConstRef(entry.rational_prime));
        entry.ramified = splitting == 0;
        if (splitting < 0 || (entry.ramified && length != 1) ||
            (!entry.ramified && !complete && length != 1) ||
            (!entry.ramified && complete && length != 2)) {
            return false;
        }

        for (slong offset = 0; offset < length; ++offset) {
            slong factor_base_index = -1;
            bool orientation_matches = false;
            if (!base.rational_prime_block_index(
                        factor_base_index, block, offset)) {
                return false;
            }
            const PrimeIdeal* prime = base.prime_at(factor_base_index);
            if (prime == nullptr ||
                !quadratic_form_prime_orientation_matches(
                        orientation_matches, *prime, discriminant,
                        flint::FmpzConstRef(entry.rational_prime),
                        flint::QfbConstRef(entry.form))) {
                return false;
            }
            if (orientation_matches) {
                entry.oriented_factor_base_index = factor_base_index;
            } else {
                entry.conjugate_factor_base_index = factor_base_index;
            }
        }
        if (!entry.ramified && length == 1 &&
            entry.oriented_factor_base_index < 0) {
            if (entry.conjugate_factor_base_index != start) {
                return false;
            }
            entry.oriented_factor_base_index = start;
            entry.conjugate_factor_base_index = -1;
            entry.factor_base_orientation = -1;
        }
        if (!entry.ramified && length == 1) {
            flint::Qfb inverse;
            if (entry.factor_base_orientation > 0) {
                flint::qfb_inverse(flint::QfbRef(inverse),
                                   flint::QfbConstRef(entry.form));
            } else {
                flint::qfb_set(flint::QfbRef(inverse),
                               flint::QfbConstRef(entry.form));
            }
            bool conjugate_ideal_ok = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.maximal_imaginary_quadratic_relations."
                        "form_setup.missing_conjugate_hnf");
                entry.conjugate_ideal = Ideal(order);
                conjugate_ideal_ok = entry.conjugate_ideal.is_defined() &&
                        quadratic_form_prime_ideal(
                                entry.conjugate_ideal, order, omega,
                                discriminant,
                                flint::FmpzConstRef(entry.rational_prime),
                                flint::QfbConstRef(inverse));
            }
            if (!conjugate_ideal_ok) {
                return false;
            }
        }
        if (entry.oriented_factor_base_index < 0 ||
            (!entry.ramified && length == 2 &&
             entry.conjugate_factor_base_index < 0)) {
            return false;
        }
        out.emplace_back(std::move(entry));
    }
    return !out.empty();
}

bool build_quadratic_subfactor_base(
        std::vector<slong>& subfactor_base,
        std::vector<QuadraticFormPrime>& form_primes,
        flint::FmpzConstRef discriminant,
        flint::FmpzConstRef fourth_root,
        const DiagnosticsContext* diagnostics) noexcept {
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "form_setup.subfactor_base");
        subfactor_base.clear();
        flint::Fmpz abs_discriminant;
        flint::fmpz_abs(flint::FmpzRef(abs_discriminant), discriminant);
        const double reduction_bound =
                std::sqrt(flint::fmpz_get_d(
                                  flint::FmpzConstRef(abs_discriminant)) /
                          3.0) +
                0.5;
        if (!std::isfinite(reduction_bound)) {
            return false;
        }

        const slong minimum = ::fmpz_bits(discriminant.raw()) > 16 ? 3 : 2;
        double product = 1.0;
        for (slong i = 0; i < static_cast<slong>(form_primes.size()); ++i) {
            QuadraticFormPrime& prime =
                    form_primes[static_cast<std::size_t>(i)];
            if (prime.ramified) {
                continue;
            }
            subfactor_base.push_back(i);
            product *= flint::fmpz_get_d(
                    flint::FmpzConstRef(prime.rational_prime));
            if (static_cast<slong>(subfactor_base.size()) >= minimum &&
                product > reduction_bound) {
                break;
            }
        }
        if (static_cast<slong>(subfactor_base.size()) < minimum) {
            return false;
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "form_setup.successive_powers");
        for (slong index : subfactor_base) {
            QuadraticFormPrime& prime =
                    form_primes[static_cast<std::size_t>(index)];
            prime.powers.clear();
            prime.powers.reserve(kQuadraticRandomExponentMax);
            // reference buch1.c:powsubFBquad builds F^j successively from F^(j-1)
            // and F instead of running an independent powering for every j.
            for (slong exponent = 1; exponent <= kQuadraticRandomExponentMax;
                 ++exponent) {
                flint::Qfb power;
                if (exponent == 1) {
                    flint::qfb_set(flint::QfbRef(power),
                                   flint::QfbConstRef(prime.form));
                } else {
                    flint::qfb_nucomp(
                            flint::QfbRef(power),
                            flint::QfbConstRef(prime.powers.back()),
                            flint::QfbConstRef(prime.form), discriminant,
                            fourth_root);
                    flint::qfb_reduce(flint::QfbRef(power),
                                      flint::QfbConstRef(power), discriminant);
                }
                prime.powers.emplace_back(std::move(power));
            }
        }
    }
    return true;
}

void compose_quadratic_forms(flint::Qfb& out,
                             const flint::Qfb& left,
                             const flint::Qfb& right,
                             flint::FmpzConstRef discriminant,
                             flint::FmpzConstRef fourth_root) noexcept {
    flint::Qfb composed;
    flint::qfb_nucomp(
            flint::QfbRef(composed), flint::QfbConstRef(left),
            flint::QfbConstRef(right), discriminant, fourth_root);
    flint::qfb_reduce(flint::QfbRef(out), flint::QfbConstRef(composed),
                      discriminant);
}

bool random_quadratic_form(
        flint::Qfb& out,
        std::vector<slong>& exponents,
        ulong& random_state,
        const std::vector<slong>& subfactor_base,
        const std::vector<QuadraticFormPrime>& form_primes,
        flint::FmpzConstRef discriminant,
        flint::FmpzConstRef fourth_root) noexcept {
    exponents.assign(subfactor_base.size(), 0);
    for (;;) {
        bool set = false;
        flint::Qfb product;
        for (std::size_t i = 0; i < subfactor_base.size(); ++i) {
            const slong exponent = static_cast<slong>(
                    detail::relation_search::next_relation_random_exponent(
                            random_state));
            exponents[i] = exponent;
            if (exponent == 0) {
                continue;
            }
            const QuadraticFormPrime& prime = form_primes[
                    static_cast<std::size_t>(subfactor_base[i])];
            const flint::Qfb& power =
                    prime.powers[static_cast<std::size_t>(exponent - 1)];
            if (!set) {
                flint::qfb_set(flint::QfbRef(product),
                               flint::QfbConstRef(power));
                set = true;
            } else {
                flint::Qfb next;
                compose_quadratic_forms(next, product, power, discriminant,
                                        fourth_root);
                product.swap(next);
            }
        }
        if (set) {
            out.swap(product);
            return true;
        }
    }
}

bool factor_quadratic_form(
        QuadraticFormFactorization& out,
        const flint::Qfb& form,
        const std::vector<QuadraticFormPrime>& form_primes,
        flint::FmpzConstRef large_prime_limit) noexcept {
    out = {};
    flint::Fmpz remaining;
    flint::fmpz_abs(flint::FmpzRef(remaining),
                    flint::qfb_a(flint::QfbConstRef(form)));
    if (flint::fmpz_is_zero(flint::FmpzConstRef(remaining))) {
        return false;
    }

    for (slong i = 0; i < static_cast<slong>(form_primes.size()); ++i) {
        const QuadraticFormPrime& prime =
                form_primes[static_cast<std::size_t>(i)];
        slong exponent = 0;
        while (flint::fmpz_divisible(
                flint::FmpzConstRef(remaining),
                flint::FmpzConstRef(prime.rational_prime))) {
            flint::fmpz_divexact(
                    flint::FmpzRef(remaining),
                    flint::FmpzConstRef(remaining),
                    flint::FmpzConstRef(prime.rational_prime));
            ++exponent;
        }
        if (exponent > 0) {
            out.factors.push_back({i, exponent});
        }
    }

    if (flint::fmpz_is_one(flint::FmpzConstRef(remaining))) {
        out.residual = 1;
        out.smooth_or_partial = true;
        return true;
    }
    if (flint::fmpz_cmp(flint::FmpzConstRef(remaining),
                        large_prime_limit) > 0 ||
        !flint::fmpz_abs_fits_ui(flint::FmpzConstRef(remaining))) {
        return true;
    }
    out.residual = flint::fmpz_get_ui(flint::FmpzConstRef(remaining));
    out.smooth_or_partial = out.residual > 1;
    return true;
}

slong oriented_factor_exponent(
        const flint::Qfb& form,
        const QuadraticFormPrime& prime,
        slong exponent) noexcept {
    const ulong p = flint::fmpz_get_ui(
            flint::FmpzConstRef(prime.rational_prime));
    if (p > std::numeric_limits<ulong>::max() / 2) {
        return 0;
    }
    const ulong residue = ::fmpz_fdiv_ui(
            flint::qfb_b(flint::QfbConstRef(form)).raw(), p << 1);
    return residue > p ? -exponent : exponent;
}

void add_form_factorization_row(
        flint::FmpzMat& row,
        const flint::Qfb& form,
        const QuadraticFormFactorization& factorization,
        const std::vector<QuadraticFormPrime>& form_primes,
        slong multiplier) noexcept {
    for (const QuadraticFormFactor& factor : factorization.factors) {
        const QuadraticFormPrime& prime = form_primes[
                static_cast<std::size_t>(factor.form_prime_index)];
        const slong exponent = oriented_factor_exponent(
                form, prime, factor.exponent);
        const slong change =
                multiplier * exponent * prime.factor_base_orientation;
        flint::FmpzRef entry = flint::fmpz_mat_entry(
                row, 0, prime.oriented_factor_base_index);
        if (change >= 0) {
            flint::fmpz_add_ui(entry, flint::FmpzConstRef(entry.raw()),
                               static_cast<ulong>(change));
        } else {
            flint::fmpz_sub_ui(entry, flint::FmpzConstRef(entry.raw()),
                               static_cast<ulong>(-change));
        }
    }
}

bool principal_generator_from_factor_base_row(
        Element& out,
        const Order& order,
        const FactorBase& base,
        const std::vector<QuadraticFormPrime>& form_primes,
        flint::FmpzMatConstRef row,
        flint::FmpzConstRef discriminant,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.maximal_imaginary_quadratic_relations.witness");
    if (order.parent() == nullptr || !out.has_parent(*order.parent())) {
        return false;
    }

    Ideal numerator(order);
    flint::Fmpz denominator;
    flint::Fmpz omega_trace;
    flint::Fmpz omega_norm;
    flint::FmpzMat hnf(2, 2);
    flint::Fmpz form_a;
    flint::Fmpz form_b;
    flint::Fmpz form_c;
    flint::Fmpz u1;
    flint::Fmpz v1;
    flint::Fmpz shortest_value;
    flint::Fmpz numerator_norm;
    flint::Fmpz product;
    flint::FmpzMat coordinates(1, 2);
    OrderElement integral_generator(order);
    Element generator(*order.parent());
    if (!numerator.is_defined() || !integral_generator.is_defined() ||
        !generator.is_defined()) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "quadratic witness allocation failed");
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "witness.row_ideal");
        if (!quadratic_factor_base_row_integral_den(
                    numerator, denominator, base, form_primes, row,
                    diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "quadratic witness row ideal reconstruction failed");
            return false;
        }
    }
    if (!quadratic_order_trace_and_norm(
                omega_trace, omega_norm, discriminant) ||
        !numerator.get_hnf(flint::FmpzMatRef(hnf))) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "quadratic witness norm-form setup failed");
        return false;
    }

    const flint::FmpzConstRef x1 = flint::fmpz_mat_entry(
            flint::FmpzMatConstRef(hnf), 0, 0);
    const flint::FmpzConstRef y1 = flint::fmpz_mat_entry(
            flint::FmpzMatConstRef(hnf), 0, 1);
    const flint::FmpzConstRef x2 = flint::fmpz_mat_entry(
            flint::FmpzMatConstRef(hnf), 1, 0);
    const flint::FmpzConstRef y2 = flint::fmpz_mat_entry(
            flint::FmpzMatConstRef(hnf), 1, 1);
    quadratic_norm_form_value(
            form_a, x1, y1, flint::FmpzConstRef(omega_trace),
            flint::FmpzConstRef(omega_norm));
    quadratic_norm_form_cross_coefficient(
            form_b, x1, y1, x2, y2,
            flint::FmpzConstRef(omega_trace),
            flint::FmpzConstRef(omega_norm));
    quadratic_norm_form_value(
            form_c, x2, y2, flint::FmpzConstRef(omega_trace),
            flint::FmpzConstRef(omega_norm));
    if (!quadratic_form_shortest_vector(
                u1, v1, shortest_value, flint::FmpzConstRef(form_a),
                flint::FmpzConstRef(form_b),
                flint::FmpzConstRef(form_c))) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "quadratic witness binary-form reduction failed");
        return false;
    }

    flint::FmpzRef coordinate_x = flint::fmpz_mat_entry(coordinates, 0, 0);
    flint::FmpzRef coordinate_y = flint::fmpz_mat_entry(coordinates, 0, 1);
    flint::fmpz_mul(coordinate_x, flint::FmpzConstRef(u1), x1);
    flint::fmpz_mul(flint::FmpzRef(product), flint::FmpzConstRef(v1), x2);
    flint::fmpz_add(coordinate_x, flint::FmpzConstRef(coordinate_x.raw()),
                    flint::FmpzConstRef(product));
    flint::fmpz_mul(coordinate_y, flint::FmpzConstRef(u1), y1);
    flint::fmpz_mul(flint::FmpzRef(product), flint::FmpzConstRef(v1), y2);
    flint::fmpz_add(coordinate_y, flint::FmpzConstRef(coordinate_y.raw()),
                    flint::FmpzConstRef(product));

    Element integral_element(*order.parent());
    // The coordinates prove membership in numerator.  Equal element and
    // ideal norms then make (integral_element) a contained ideal of index one.
    if (!integral_element.is_defined() ||
        !integral_generator.set_coordinates(
                flint::FmpzMatConstRef(coordinates)) ||
        !integral_generator.get_element(integral_element) ||
        !numerator.norm(flint::FmpzRef(numerator_norm)) ||
        !flint::fmpz_equal(flint::FmpzConstRef(shortest_value),
                           flint::FmpzConstRef(numerator_norm)) ||
        !element_scalar_div_fmpz(
                generator, integral_element,
                flint::FmpzConstRef(denominator))) {
        SILEX_LOG(diagnostics, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "quadratic witness exact norm/index verification failed");
        return false;
    }
    return out.set(generator);
}

bool fmpz_mat_row_add_si(flint::FmpzMat& row,
                         slong column,
                         slong value) noexcept {
    if (column < 0 || column >= flint::fmpz_mat_ncols(row)) {
        return false;
    }
    flint::FmpzRef entry = flint::fmpz_mat_entry(row, 0, column);
    if (value >= 0) {
        flint::fmpz_add_ui(entry, flint::FmpzConstRef(entry.raw()),
                           static_cast<ulong>(value));
    } else {
        flint::fmpz_sub_ui(entry, flint::FmpzConstRef(entry.raw()),
                           static_cast<ulong>(-value));
    }
    return true;
}

const QuadraticLargePrimePartial* matching_large_prime_partial(
        const std::vector<QuadraticLargePrimePartial>& partials,
        ulong residual,
        const std::vector<slong>& exponents,
        slong pivot) noexcept {
    for (const QuadraticLargePrimePartial& partial : partials) {
        if (partial.residual != residual) {
            continue;
        }
        if (partial.exponents != exponents || partial.pivot != pivot) {
            return &partial;
        }
    }
    return nullptr;
}

bool same_or_opposite_middle_coefficient(bool& same,
                                         const flint::Qfb& left,
                                         const flint::Qfb& right,
                                         ulong residual) noexcept {
    if (residual == 0 ||
        residual > std::numeric_limits<ulong>::max() / 2) {
        return false;
    }
    const ulong modulus = residual << 1;
    const ulong left_b = ::fmpz_fdiv_ui(
            flint::qfb_b(flint::QfbConstRef(left)).raw(), modulus);
    const ulong right_b = ::fmpz_fdiv_ui(
            flint::qfb_b(flint::QfbConstRef(right)).raw(), modulus);
    same = left_b == right_b;
    return same || left_b + right_b == modulus;
}

}  // namespace

ClassGroupContext::SpecializedRelationBackendStatus
ClassGroupContext::run_maximal_imaginary_quadratic_relation_backend_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options) noexcept {
    const DiagnosticsContext* active_diagnostics =
            options.diagnostics != nullptr ? options.diagnostics : diagnostics_;
    SILEX_PROFILE_SCOPE(
            active_diagnostics, DiagnosticsModule::class_group,
            "class_group.maximal_imaginary_quadratic_relations");

    Signature signature_value;
    flint::Fmpz discriminant;
    flint::Fmpz conductor;
    if (order.parent() == nullptr || order.degree() != 2 ||
        !order.is_maximal() ||
        order.parent()->backend_kind() != NumberFieldBackendKind::quadratic ||
        !order.quadratic_conductor(flint::FmpzRef(conductor)) ||
        !flint::fmpz_is_one(flint::FmpzConstRef(conductor)) ||
        !signature(signature_value, *order.parent()) ||
        signature_value.r1() != 0 || signature_value.r2() != 1 ||
        !order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0 ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(discriminant)) ||
        !flint::fmpz_abs_fits_ui(factor_base_bound) ||
        options.target_relation_kernel_units != 0 ||
        options.max_candidates <= 0 || options.max_relations <= 0 ||
        base_.parent() == nullptr || !base_.parent()->has_same_data(order)) {
        return SpecializedRelationBackendStatus::unavailable;
    }

    flint::Fmpz abs_discriminant;
    flint::Fmpz fourth_root;
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    ::fmpz_root(fourth_root.raw(), abs_discriminant.raw(), 4);

    std::vector<QuadraticFormPrime> form_primes;
    std::vector<slong> subfactor_base;
    flint::Fmpz exact_exceptional_class_order;
    bool exact_discriminant_minus_three = false;
    bool form_setup_ok = false;
    {
        SILEX_PROFILE_SCOPE(
                active_diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "form_setup");
        const bool form_primes_ready = identify_quadratic_form_primes(
                form_primes, order, base_,
                flint::FmpzConstRef(discriminant), active_diagnostics);
        exact_discriminant_minus_three =
                form_primes_ready &&
                flint::fmpz_equal_si(
                        flint::FmpzConstRef(discriminant), -3) &&
                detail::ClassGroupCertificationAccess::
                        exact_imaginary_quadratic_class_order_for_run(
                                flint::FmpzRef(
                                        exact_exceptional_class_order),
                                *this,
                                flint::FmpzConstRef(discriminant)) &&
                flint::fmpz_is_one(
                        flint::FmpzConstRef(
                                exact_exceptional_class_order));
        form_setup_ok = form_primes_ready &&
                (exact_discriminant_minus_three ||
                 build_quadratic_subfactor_base(
                         subfactor_base, form_primes,
                         flint::FmpzConstRef(discriminant),
                         flint::FmpzConstRef(fourth_root),
                         active_diagnostics));
    }
    if (!form_setup_ok) {
        return SpecializedRelationBackendStatus::unavailable;
    }

    if (base_.length() > options.max_relations) {
        return SpecializedRelationBackendStatus::failed;
    }

    flint::Fmpz required_generation_bound;
    bool compact_base_honest = true;
    if (!factor_base_class_group_bound(
                flint::FmpzRef(required_generation_bound), order)) {
        return SpecializedRelationBackendStatus::failed;
    }
    const bool compact_base_needs_honesty =
            flint::fmpz_cmp(
                    factor_base_bound,
                    flint::FmpzConstRef(required_generation_bound)) < 0;
    if (compact_base_needs_honesty) {
        SILEX_PROFILE_SCOPE(
                active_diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "factor_base_honesty");
        const ulong honesty_seed =
                detail::relation_search::relation_search_phase_seed(
                        *this, base_.length(), kQuadraticHonestyPhase, 0);
        if (!detail::relation_search::factor_base_honesty_check(
                    compact_base_honest, base_, factor_base_bound,
                    flint::FmpzConstRef(required_generation_bound), nullptr,
                    honesty_seed, false,
                    kQuadraticPrincipalReductionPrecision,
                    active_diagnostics)) {
            return SpecializedRelationBackendStatus::failed;
        }
        if (!compact_base_honest) {
            return SpecializedRelationBackendStatus::unavailable;
        }
    }

    if (exact_discriminant_minus_three) {
        // The reduced-form enumeration above proves class number one.  Build
        // the identity relation lattice directly, reconstructing and exactly
        // verifying a principal generator for every factor-base ideal.  This
        // is the deterministic exceptional-field analogue of the witnessed
        // HNF publication below; it does not relax the random-search restart
        // envelope used by nontrivial class groups.
        if (base_.length() > options.max_candidates) {
            return SpecializedRelationBackendStatus::failed;
        }
        for (slong i = 0; i < base_.length(); ++i) {
            flint::FmpzMat row(1, base_.length());
            flint::fmpz_one(flint::fmpz_mat_entry(row, 0, i));
            Element generator(*order.parent());
            Relation relation(base_);
            RelationAppendOutcome outcome = RelationAppendOutcome::none;
            if (!generator.is_defined() || !relation.is_defined() ||
                !principal_generator_from_factor_base_row(
                        generator, order, base_, form_primes,
                        flint::FmpzMatConstRef(row),
                        flint::FmpzConstRef(discriminant),
                        active_diagnostics) ||
                !detail::set_relation_from_known_row(
                        relation, base_, generator,
                        flint::FmpzMatConstRef(row)) ||
                !append_relation_with_outcome_(
                        outcome, relation, ClassGroupRelationSource::Search,
                        DependentRelationPolicy::keep)) {
                return SpecializedRelationBackendStatus::failed;
            }
        }

        if (!publish_presentation()) {
            return SpecializedRelationBackendStatus::failed;
        }
        const bool generation_verified = compact_base_needs_honesty
                ? detail::ClassGroupCertificationAccess::
                          record_factor_base_honesty_proof(
                                  *this,
                                  flint::FmpzConstRef(
                                          required_generation_bound))
                : check_factor_base_generation_bound(
                          flint::FmpzConstRef(required_generation_bound));
        if (!generation_verified) {
            return SpecializedRelationBackendStatus::failed;
        }
        if (options.requested_certification == CertificationMode::proven &&
            !detail::ClassGroupCertificationAccess::
                    try_certify_imaginary_quadratic_from_exact_order(
                            *this, options.requested_certification,
                            flint::FmpzConstRef(discriminant),
                            flint::FmpzConstRef(
                                    exact_exceptional_class_order))) {
            return SpecializedRelationBackendStatus::failed;
        }
        return SpecializedRelationBackendStatus::succeeded;
    }

    flint::Fmpz large_prime_limit;
    flint::fmpz_mul(flint::FmpzRef(large_prime_limit), factor_base_bound,
                    factor_base_bound);

    relation_kernel_units_target_ = 0;
    configure_partial_relations_(options);

    fmpz_smat::HnfContext staged_relation_module;
    staged_relation_module.set_diagnostics(active_diagnostics);
    // Candidate form rows stay private.  Only the canonical HNF basis below
    // is admitted after reconstructing an exact public element witness.
    std::vector<flint::FmpzMat> staged_relation_rows;
    std::vector<ClassGroupRelationSource> staged_relation_sources;
    if (!staged_relation_module.reset(
                base_.length(), kQuadraticRelationModulePrime)) {
        return SpecializedRelationBackendStatus::failed;
    }

    auto stage_relation =
            [&](flint::FmpzMatConstRef row,
                ClassGroupRelationSource source) noexcept {
        bool independent = false;
        bool index_refined = false;
        if (!staged_relation_module.add_fmpz_mat_row(
                    &independent, &index_refined, row, 0)) {
            return false;
        }
        flint::FmpzMat stored(1, base_.length());
        flint::fmpz_mat_set(flint::FmpzMatRef(stored), row);
        staged_relation_rows.emplace_back(std::move(stored));
        staged_relation_sources.push_back(source);
        return true;
    };

    auto append_verified_relation =
            [&](flint::FmpzMatConstRef row,
                const Element& generator,
                ClassGroupRelationSource source) noexcept {
        Relation relation(base_);
        RelationAppendOutcome outcome = RelationAppendOutcome::none;
        return relation.is_defined() &&
               detail::set_relation_from_known_row(
                       relation, base_, generator, row) &&
               append_relation_with_outcome_(
                       outcome, relation, source, DependentRelationPolicy::keep);
    };

    {
        SILEX_PROFILE_SCOPE(
                active_diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "rational_relations");
        for (const QuadraticFormPrime& prime : form_primes) {
            if (!prime.ramified &&
                prime.conjugate_factor_base_index < 0) {
                continue;
            }
            flint::FmpzMat row(1, base_.length());
            if (!fmpz_mat_row_add_si(
                        row, prime.oriented_factor_base_index,
                                     prime.ramified ? 2 : 1) ||
                (!prime.ramified &&
                 !fmpz_mat_row_add_si(
                         row, prime.conjugate_factor_base_index, 1)) ||
                !stage_relation(
                        flint::FmpzMatConstRef(row),
                        ClassGroupRelationSource::Search)) {
                SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                          LogLevel::detail,
                          "quadratic rational-prime relation failed");
                return SpecializedRelationBackendStatus::failed;
            }
        }
    }

    slong target_relations = base_.length();
    slong candidates_tried = 0;
    slong current = 0;
    ulong random_state = detail::relation_search::relation_search_phase_seed(
            *this, base_.length(), kQuadraticRelationPhase, 0);
    std::vector<QuadraticLargePrimePartial> partials;
    std::vector<slong> exponents;
    flint::Fmpz exact_class_order;
    bool exact_class_order_ready = false;

    while (true) {
        const slong staged_relation_count = static_cast<slong>(
                staged_relation_rows.size());
        if (staged_relation_count >= target_relations) {
            const bool full_rank =
                    staged_relation_module.rank() >= generator_count();
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      full_rank
                              ? "quadratic relation checkpoint has full rank"
                              : "quadratic relation checkpoint lacks full rank");
            bool exact_relation_lattice = false;
            if (full_rank) {
                flint::Fmpz candidate_order;
                bool candidate_order_ready = false;
                {
                    SILEX_PROFILE_SCOPE(
                            active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group."
                            "maximal_imaginary_quadratic_relations."
                            "checkpoint.index");
                    candidate_order_ready =
                            staged_relation_module.full_rank_index(
                                    flint::FmpzRef(candidate_order));
                }
                if (candidate_order_ready) {
                    if (!exact_class_order_ready) {
                        SILEX_PROFILE_SCOPE(
                                active_diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group."
                                "maximal_imaginary_quadratic_relations."
                                "checkpoint.exact_order");
                        exact_class_order_ready = detail::
                                ClassGroupCertificationAccess::
                                        exact_imaginary_quadratic_class_order_for_run(
                                                flint::FmpzRef(
                                                        exact_class_order),
                                                *this,
                                                flint::FmpzConstRef(
                                                        discriminant));
                    }
                    exact_relation_lattice =
                            exact_class_order_ready &&
                            flint::fmpz_equal(
                                    flint::FmpzConstRef(candidate_order),
                                    flint::FmpzConstRef(exact_class_order));
                }
            }
            bool published = false;
            if (exact_relation_lattice) {
                SILEX_PROFILE_SCOPE(
                        active_diagnostics, DiagnosticsModule::class_group,
                        "class_group.maximal_imaginary_quadratic_relations."
                        "checkpoint.witness_hnf");
                flint::FmpzMat hnf_basis(
                        generator_count(), generator_count());
                if (!staged_relation_module.get_hnf_rows(
                            flint::FmpzMatRef(hnf_basis))) {
                    return SpecializedRelationBackendStatus::failed;
                }
                for (slong row_index = 0;
                     row_index < generator_count(); ++row_index) {
                    flint::FmpzMatConstWindow row_window(
                            hnf_basis, row_index, 0, row_index + 1,
                            generator_count());
                    ClassGroupRelationSource source =
                            ClassGroupRelationSource::RandomProduct;
                    for (std::size_t staged_index = 0;
                         staged_index < staged_relation_rows.size();
                         ++staged_index) {
                        if (flint::fmpz_mat_equal(
                                    row_window.const_ref(),
                                    flint::FmpzMatConstRef(
                                            staged_relation_rows[
                                                    staged_index]))) {
                            source = staged_relation_sources[staged_index];
                            break;
                        }
                    }
                    Element generator(*order.parent());
                    if (!generator.is_defined() ||
                        !principal_generator_from_factor_base_row(
                                generator, order, base_, form_primes,
                                row_window.const_ref(),
                                flint::FmpzConstRef(discriminant),
                                active_diagnostics) ||
                        !append_verified_relation(
                                row_window.const_ref(), generator, source)) {
                        return SpecializedRelationBackendStatus::failed;
                    }
                }
                published = publish_presentation();
            }
            if (published) {
                bool generation_verified = false;
                {
                    SILEX_PROFILE_SCOPE(
                            active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group."
                            "maximal_imaginary_quadratic_relations."
                            "checkpoint.generation");
                    generation_verified = compact_base_needs_honesty
                            ? detail::ClassGroupCertificationAccess::
                                      record_factor_base_honesty_proof(
                                              *this,
                                              flint::FmpzConstRef(
                                                      required_generation_bound))
                            : check_factor_base_generation_bound(
                                      flint::FmpzConstRef(
                                              required_generation_bound));
                }
                bool exact_presentation = false;
                if (generation_verified &&
                    options.requested_certification ==
                            CertificationMode::proven) {
                    SILEX_PROFILE_SCOPE(
                            active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group."
                            "maximal_imaginary_quadratic_relations."
                            "checkpoint.certification");
                    exact_presentation =
                            detail::ClassGroupCertificationAccess::
                                    try_certify_imaginary_quadratic_from_exact_order(
                                            *this,
                                            options.requested_certification,
                                            flint::FmpzConstRef(
                                                    discriminant),
                                            flint::FmpzConstRef(
                                                    exact_class_order));
                } else if (generation_verified) {
                    exact_presentation = true;
                }
                if (exact_presentation) {
                    SILEX_PROFILE_EVENT(
                            active_diagnostics,
                            DiagnosticsModule::class_group,
                            "class_group."
                            "maximal_imaginary_quadratic_relations.success");
                    return SpecializedRelationBackendStatus::succeeded;
                }
            }
            if (target_relations >= options.max_relations) {
                return SpecializedRelationBackendStatus::failed;
            }
            ++target_relations;
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "quadratic relation checkpoint needs more relations");
        }

        SILEX_PROFILE_SCOPE(
                active_diagnostics, DiagnosticsModule::class_group,
                "class_group.maximal_imaginary_quadratic_relations."
                "candidate");

        if (candidates_tried >= options.max_candidates ||
            staged_relation_count >= options.max_relations) {
            return SpecializedRelationBackendStatus::failed;
        }

        flint::Qfb random_form;
        if (!random_quadratic_form(
                    random_form, exponents, random_state, subfactor_base,
                    form_primes, flint::FmpzConstRef(discriminant),
                    flint::FmpzConstRef(fourth_root))) {
            return SpecializedRelationBackendStatus::failed;
        }
        flint::Qfb form;
        compose_quadratic_forms(
                form, random_form,
                form_primes[static_cast<std::size_t>(current)].form,
                flint::FmpzConstRef(discriminant),
                flint::FmpzConstRef(fourth_root));
        ++candidates_tried;

        QuadraticFormFactorization factorization;
        if (!factor_quadratic_form(
                    factorization, form, form_primes,
                    flint::FmpzConstRef(large_prime_limit))) {
            return SpecializedRelationBackendStatus::failed;
        }
        if (!factorization.smooth_or_partial) {
            if ((candidates_tried & 0xff) == 0) {
                current = (current + 1) %
                        static_cast<slong>(form_primes.size());
            }
            continue;
        }

        flint::FmpzMat row(1, base_.length());
        ClassGroupRelationSource source =
                ClassGroupRelationSource::RandomProduct;
        bool have_relation = factorization.residual == 1;
        if (have_relation) {
            for (std::size_t i = 0; i < subfactor_base.size(); ++i) {
                const QuadraticFormPrime& prime = form_primes[
                        static_cast<std::size_t>(subfactor_base[i])];
                if (!fmpz_mat_row_add_si(
                            row, prime.oriented_factor_base_index,
                            -exponents[i] *
                                    prime.factor_base_orientation)) {
                    return SpecializedRelationBackendStatus::failed;
                }
            }
            add_form_factorization_row(
                    row, form, factorization, form_primes, 1);
        } else {
            const QuadraticLargePrimePartial* match =
                    matching_large_prime_partial(
                            partials, factorization.residual, exponents,
                            current);
            if (match == nullptr) {
                partials.push_back(
                        {factorization.residual, exponents, current});
                continue;
            }

            flint::Qfb saved_random_form;
            std::vector<slong> saved_exponents = match->exponents;
            bool set = false;
            for (std::size_t i = 0; i < subfactor_base.size(); ++i) {
                const slong exponent = saved_exponents[i];
                if (exponent == 0) {
                    continue;
                }
                const QuadraticFormPrime& prime = form_primes[
                        static_cast<std::size_t>(subfactor_base[i])];
                const flint::Qfb& power =
                        prime.powers[static_cast<std::size_t>(exponent - 1)];
                if (!set) {
                    flint::qfb_set(
                            flint::QfbRef(saved_random_form),
                            flint::QfbConstRef(power));
                    set = true;
                } else {
                    flint::Qfb next;
                    compose_quadratic_forms(
                            next, saved_random_form, power,
                            flint::FmpzConstRef(discriminant),
                            flint::FmpzConstRef(fourth_root));
                    saved_random_form.swap(next);
                }
            }
            if (!set) {
                return SpecializedRelationBackendStatus::failed;
            }
            flint::Qfb form2;
            compose_quadratic_forms(
                    form2, saved_random_form,
                    form_primes[static_cast<std::size_t>(match->pivot)].form,
                    flint::FmpzConstRef(discriminant),
                    flint::FmpzConstRef(fourth_root));
            bool same = false;
            if (!same_or_opposite_middle_coefficient(
                        same, form2, form, factorization.residual)) {
                continue;
            }

            QuadraticFormFactorization factorization2;
            if (!factor_quadratic_form(
                        factorization2, form2, form_primes,
                        flint::FmpzConstRef(large_prime_limit)) ||
                !factorization2.smooth_or_partial ||
                factorization2.residual != factorization.residual) {
                return SpecializedRelationBackendStatus::failed;
            }

            add_form_factorization_row(
                    row, form, factorization, form_primes, 1);
            for (std::size_t i = 0; i < subfactor_base.size(); ++i) {
                const QuadraticFormPrime& prime = form_primes[
                        static_cast<std::size_t>(subfactor_base[i])];
                const slong exponent = same
                        ? saved_exponents[i] - exponents[i]
                        : -saved_exponents[i] - exponents[i];
                if (!fmpz_mat_row_add_si(
                            row, prime.oriented_factor_base_index,
                            exponent * prime.factor_base_orientation)) {
                    return SpecializedRelationBackendStatus::failed;
                }
            }
            add_form_factorization_row(
                    row, form2, factorization2, form_primes,
                    same ? -1 : 1);
            if (!fmpz_mat_row_add_si(
                        row,
                        form_primes[static_cast<std::size_t>(match->pivot)]
                                .oriented_factor_base_index,
                        (same ? 1 : -1) *
                                form_primes[static_cast<std::size_t>(
                                        match->pivot)]
                                        .factor_base_orientation)) {
                return SpecializedRelationBackendStatus::failed;
            }
            source = ClassGroupRelationSource::LargePrimeMatch;
            have_relation = true;
        }

        if (!have_relation ||
            !fmpz_mat_row_add_si(
                    row,
                    form_primes[static_cast<std::size_t>(current)]
                            .oriented_factor_base_index,
                    -form_primes[static_cast<std::size_t>(current)]
                             .factor_base_orientation)) {
            return SpecializedRelationBackendStatus::failed;
        }

        if (!stage_relation(flint::FmpzMatConstRef(row), source)) {
            SILEX_LOG(active_diagnostics, DiagnosticsModule::class_group,
                      LogLevel::detail,
                      "quadratic random relation staging failed");
            return SpecializedRelationBackendStatus::failed;
        }
        current = (current + 1) % static_cast<slong>(form_primes.size());
    }
}

}  // namespace silex
