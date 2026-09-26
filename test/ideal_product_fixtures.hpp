#pragma once

#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/lat.hpp>

namespace silex::test::ideal_product {

// A deliberately independent oracle: form every pair of basis products.
inline bool exhaustive_product(Ideal& out, const Ideal& left,
                               const Ideal& right) noexcept {
    const Order* order = left.parent();
    if (order == nullptr || right.parent() == nullptr ||
        !order->has_same_data(*right.parent())) {
        return false;
    }
    const slong n = order->degree();
    flint::FmpzMat h(n, n), k(n, n), row(1, n), rows(n * n, n);
    OrderElement a(*order), b(*order), ab(*order);
    if (!left.get_hnf(flint::FmpzMatRef(h)) ||
        !right.get_hnf(flint::FmpzMatRef(k))) {
        return false;
    }
    for (slong i = 0; i < n; ++i) {
        for (slong c = 0; c < n; ++c) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, c),
                     fmpz_mat_entry(h.raw(), i, c));
        }
        if (!a.set_coordinates(flint::FmpzMatConstRef(row))) return false;
        for (slong j = 0; j < n; ++j) {
            for (slong c = 0; c < n; ++c) {
                fmpz_set(fmpz_mat_entry(row.raw(), 0, c),
                         fmpz_mat_entry(k.raw(), j, c));
            }
            if (!b.set_coordinates(flint::FmpzMatConstRef(row)) ||
                !ab.multiply(a, b) ||
                !ab.get_coordinates(flint::FmpzMatRef(row))) return false;
            for (slong c = 0; c < n; ++c) {
                fmpz_set(fmpz_mat_entry(rows.raw(), i * n + j, c),
                         fmpz_mat_entry(row.raw(), 0, c));
            }
        }
    }
    lat::Lat lattice(n), reduced(n);
    return lattice.set_basis(flint::FmpzMatConstRef(rows)) &&
           lattice.hnf(reduced) && out.set_hnf(reduced.basis_ref());
}

inline bool polynomial(flint::FmpqPoly& f, slong degree,
                       bool nonmaximal = false) noexcept {
    fmpq_poly_zero(f.raw());
    fmpq_poly_set_coeff_si(f.raw(), degree, 1);
    if (degree == 2) {
        fmpq_poly_set_coeff_si(f.raw(), 0, nonmaximal ? -5 : 47);
    } else if (degree == 3) {
        fmpq_poly_set_coeff_si(f.raw(), 0, -2);
    } else if (degree == 4) {
        fmpq_poly_set_coeff_si(f.raw(), 0, -5);
        fmpq_poly_set_coeff_si(f.raw(), 1, 3);
        fmpq_poly_set_coeff_si(f.raw(), 2, -2);
    } else {
        fmpq_poly_set_coeff_si(f.raw(), 0, -1);
        fmpq_poly_set_coeff_si(f.raw(), 1, -1);
    }
    return true;
}

// (f(k), theta-k) is an exact integral ideal; no principal metadata is set.
inline bool evaluation_ideal(Ideal& out, const flint::FmpqPoly& f,
                             slong k) noexcept {
    const Order* order = out.parent();
    if (order == nullptr || order->parent() == nullptr) return false;
    Element theta(*order->parent()), value(*order->parent());
    OrderElement beta(*order);
    flint::Fmpq x, fx;
    flint::Fmpz a;
    fmpq_set_si(x.raw(), k, 1);
    fmpq_poly_evaluate_fmpq(fx.raw(), f.raw(), x.raw());
    fmpz_abs(a.raw(), fmpq_numref(fx.raw()));
    return theta.gen() && value.add_si(theta, -k) &&
           beta.set_element(value) &&
           detail::set_known_two_generator_ideal(
                   out, flint::FmpzConstRef(a), beta);
}

// kinds: scalar, principal, mixed, general, square, nonmaximal.
inline bool fixture(Order& order, Ideal& left, Ideal& right,
                    slong degree, slong kind) noexcept {
    flint::FmpqPoly f;
    polynomial(f, degree, kind == 5);
    NumberField field = NumberField::by_polynomial(flint::FmpqPolyConstRef(f));
    Order equation = Order::equation_order(field);
    if (!field.is_defined() || !equation.is_defined() ||
        !order.define(field)) return false;
    if (kind == 5) {
        if (!order.set(equation)) return false;
    } else if (!order.maximal_order(equation)) {
        return false;
    }
    if (!left.define(order) || !right.define(order)) return false;
    OrderElement a(order), b(order);
    if (kind == 0) {
        return a.set_si(-2) && b.set_si(3) &&
               left.set_principal(a) && right.set_principal(b);
    }
    Element theta(field), alpha(field), beta(field);
    if (!theta.gen() || !alpha.add_si(theta, 2) ||
        !beta.add_si(theta, 3) || !a.set_element(alpha) ||
        !b.set_element(beta)) return false;
    if (kind == 1) return left.set_principal(a) && right.set_principal(b);
    if (!evaluation_ideal(right, f, kind == 5 ? 4 : 3)) return false;
    if (kind == 2) return left.set_principal(a);
    if (!evaluation_ideal(left, f, kind == 5 ? 3 : 2)) return false;
    if (kind == 4) return right.set(left);
    return !left.is_one() && !right.is_one();
}

inline bool random_search_fixture(Order& order, Ideal& ideal) noexcept {
    Order original;
    Ideal left, right;
    if (!fixture(original, left, right, 3, 3)) return false;
    flint::FmpqMat basis(3, 3), u(3, 3), changed(3, 3);
    if (!original.get_basis(flint::FmpqMatRef(basis))) return false;
    const slong entries[3][3] = {{7, 5, -21}, {-3, -2, 9}, {-4, -4, 13}};
    for (slong i = 0; i < 3; ++i)
        for (slong j = 0; j < 3; ++j)
            fmpq_set_si(fmpq_mat_entry(u.raw(), i, j), entries[i][j], 1);
    fmpq_mat_mul(changed.raw(), u.raw(), basis.raw());
    // det(u)=1: the order is unchanged by this basis replacement.  Order has
    // no public maximality setter, so maximality is recomputed; for an
    // already maximal input the result keeps the replacement basis.
    Order replaced = Order::from_basis(*original.parent(), flint::FmpqMatConstRef(changed));
    if (!replaced.is_defined() || !order.define(*original.parent()) ||
        !order.maximal_order(replaced) || !order.is_maximal()) return false;
    flint::FmpqMat kept(3, 3);
    if (!order.get_basis(flint::FmpqMatRef(kept)) ||
        fmpq_mat_equal(kept.raw(), changed.raw()) == 0) return false;
    flint::FmpqPoly f;
    polynomial(f, 3);
    return ideal.define(order) && evaluation_ideal(ideal, f, 2);
}

}  // namespace silex::test::ideal_product
