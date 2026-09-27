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

}  // namespace silex::test
