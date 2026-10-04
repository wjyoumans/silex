#pragma once

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/signature.hpp>

#include <cassert>

namespace silex::test {

inline NumberField field_by_polynomial(
        flint::FmpqPolyConstRef polynomial) noexcept {
    NumberField field = NumberField::by_polynomial(polynomial);
    assert(field.is_defined());
    return field;
}

inline NumberField field_by_polynomial(
        flint::FmpzPolyConstRef polynomial) noexcept {
    NumberField field = NumberField::by_polynomial(polynomial);
    assert(field.is_defined());
    return field;
}

// Degree-4 field with two conjugate pairs whose imaginary parts agree to
// about 303 bits:
//
//     f = (x^2 + 2x + 2)(x^2 - 8x + 17) - 2^-300 x
//       = x^4 - 6x^3 + 3x^2 + (18 - 2^-300)x + 34.
//
// The unperturbed roots are A = -1 +/- i and B = 4 +/- i.  To first order the
// perturbation moves a root r by 2^-300 r / f0'(r), which gives
// |Im A| - |Im B| ~ +4.9e-92 (about 2^-303.3), so the pair with the smaller
// real part has the larger modulus of its imaginary part.
//
// FLINT's _acb_vec_sort_pretty (acb/vec_sort_pretty.c, acb_cmp_pretty) orders
// by |Im| and falls back to the real part only when the |Im| difference ball
// contains zero.  With roots accurate to about 64 bits the difference is
// unresolved and A sorts first; with roots accurate to 512 bits it is resolved
// and B sorts first.  A refine path that re-sorts the refined roots with
// that comparator lets a 64 -> 512 refine swap the two complex places;
// the place-order tests in t-embedding.cpp and t-archimedean.cpp failed on
// that code.
inline NumberField close_imaginary_pairs_field() noexcept {
    flint::FmpqPoly polynomial;
    flint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    flint::fmpq_poly_set_coeff_si(polynomial, 3, -6);
    flint::fmpq_poly_set_coeff_si(polynomial, 2, 3);
    flint::fmpq_poly_set_coeff_si(polynomial, 0, 34);

    fmpq_t linear;
    fmpq_init(linear);
    fmpz_one(fmpq_denref(linear));
    fmpz_mul_2exp(fmpq_denref(linear), fmpq_denref(linear), 300);
    fmpz_mul_si(fmpq_numref(linear), fmpq_denref(linear), 18);
    fmpz_sub_ui(fmpq_numref(linear), fmpq_numref(linear), 1);
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, linear);
    fmpq_clear(linear);

    return field_by_polynomial(flint::FmpqPolyConstRef(polynomial));
}

inline NumberField quadratic_field(slong radicand) noexcept {
    flint::Fmpz d;
    flint::fmpz_set_si(flint::FmpzRef(d), radicand);
    NumberField field = NumberField::quadratic(flint::FmpzConstRef(d));
    assert(field.is_defined());
    return field;
}

// White-box helper: returns a defined field whose FLINT context has been
// replaced by one for the non-squarefree polynomial x^2.  NumberField
// construction rejects x^2, so this is the only way to keep the internal
// non-squarefree failure paths covered.  The field starts as the generic
// (non-quadratic-backend) field defined by x^2 - x - 1 so that no quadratic
// backend data disagrees with the installed polynomial.
inline NumberField nonsquarefree_field() noexcept {
    flint::FmpqPoly polynomial;
    flint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    flint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    flint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    NumberField field =
            field_by_polynomial(flint::FmpqPolyConstRef(polynomial));
    assert(field.backend_kind() == NumberFieldBackendKind::generic);

    flint::FmpqPoly square;
    flint::fmpq_poly_set_coeff_si(square, 2, 1);
    nf_struct* raw = field.raw_flint_field();
    assert(raw != nullptr);
    nf_clear(raw);
    nf_init(raw, square.raw());
    assert(field.degree() == 2);
    return field;
}

inline void check_square_rejected() noexcept {
    flint::FmpqPoly square;
    flint::fmpq_poly_set_coeff_si(square, 2, 1);
    assert(!NumberField::by_polynomial(flint::FmpqPolyConstRef(square))
                    .is_defined());
}

inline Order equation_order(const NumberField& parent) noexcept {
    Order order = Order::equation_order(parent);
    assert(order.is_defined());
    return order;
}

// Returns `order` with computed maximality.  Order has no public maximality
// setter, so a fixture that already names the maximal order (for example in a
// permuted basis) runs it through Order::maximal_order.  The assertions check
// that `order` was maximal and that the result keeps its exact basis, so the
// fixture is unchanged apart from the computed maximality flag.
inline Order verified_maximal_order(const Order& order) noexcept {
    assert(order.is_defined() && order.parent() != nullptr);
    Order maximal(*order.parent());
    assert(maximal.maximal_order(order));
    assert(maximal.maximality_known() && maximal.is_maximal());
    const auto input_basis = order.basis();
    const auto output_basis = maximal.basis();
    assert(input_basis.has_value() && output_basis.has_value());
    assert(flint::fmpq_mat_equal(*input_basis, *output_basis));
    return maximal;
}

// Minkowski's bound with 2^r2 in place of (4/pi)^r2 in exact integer form,
// max(1, ceil(n! ceil(sqrt|d|) 2^r2 / n^n)): the proven generation bound that
// factor_base_class_group_bound returned in degree n >= 3 before it adopted
// Zimmert's bound.
inline void former_generic_generation_bound(flint::Fmpz& out,
                                            const flint::Fmpz& abs_discriminant,
                                            slong degree,
                                            slong complex_pairs) noexcept {
    flint::Fmpz root;
    flint::Fmpz remainder;
    flint::Fmpz denominator;
    fmpz_sqrtrem(root.raw(), remainder.raw(), abs_discriminant.raw());
    if (!fmpz_is_zero(remainder.raw())) {
        fmpz_add_ui(root.raw(), root.raw(), 1);
    }
    fmpz_fac_ui(out.raw(), static_cast<ulong>(degree));
    fmpz_mul(out.raw(), out.raw(), root.raw());
    fmpz_mul_2exp(out.raw(), out.raw(), static_cast<ulong>(complex_pairs));
    fmpz_set_ui(denominator.raw(), static_cast<ulong>(degree));
    fmpz_pow_ui(denominator.raw(), denominator.raw(),
                static_cast<ulong>(degree));
    fmpz_cdiv_q(out.raw(), out.raw(), denominator.raw());
    if (fmpz_is_zero(out.raw())) {
        fmpz_one(out.raw());
    }
}

// The proven generation bound that factor_base_class_group_bound returned
// before it adopted Zimmert's bound: 1 in degree one, floor(sqrt(|d|/3)) and
// floor(sqrt(d))/2 for imaginary and real quadratic fields, and otherwise
// Minkowski's bound with 2^r2 in place of (4/pi)^r2,
// max(1, ceil(n! ceil(sqrt|d|) 2^r2 / n^n)).  The bound tests compare the
// current bound with it, and fixtures that were built on a factor base of
// this size use it so that they keep exercising the same scenario.
inline bool former_generation_bound(flint::Fmpz& out,
                                    const Order& order) noexcept {
    if (!order.is_maximal() || order.parent() == nullptr) {
        return false;
    }
    const slong degree = order.degree();
    if (degree == 1) {
        fmpz_one(out.raw());
        return true;
    }
    flint::Fmpz discriminant;
    Signature sig;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        !sig.compute(*order.parent())) {
        return false;
    }
    flint::Fmpz abs_discriminant;
    fmpz_abs(abs_discriminant.raw(), discriminant.raw());
    if (degree == 2) {
        if (fmpz_sgn(discriminant.raw()) < 0) {
            fmpz_fdiv_q_ui(out.raw(), abs_discriminant.raw(), 3);
            fmpz_sqrt(out.raw(), out.raw());
        } else {
            fmpz_sqrt(out.raw(), abs_discriminant.raw());
            fmpz_fdiv_q_2exp(out.raw(), out.raw(), 1);
        }
        return true;
    }
    former_generic_generation_bound(out, abs_discriminant, degree, sig.r2());
    return true;
}

}  // namespace silex::test
