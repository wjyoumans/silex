// Conductor of an order in an overorder and contraction of an ideal to a
// suborder.  Internal routines; see order_internal.hpp for the contracts.
//
// Source lineage:
// - reference v0.38.6, src/NumFieldOrd/NfOrd/NfOrd.jl:1488,
//   `conductor(R::AbsSimpleNumFieldOrder, S::AbsSimpleNumFieldOrder)`.
// - reference v0.38.6, src/NumFieldOrd/NfOrd/Ideal/Arithmetic.jl:952,
//   `contract(A::AbsNumFieldOrderIdeal, O::AbsNumFieldOrder)`.
// Both routines keep the reference's row conventions: an order basis matrix
// has the basis elements as rows in power-basis coordinates, and an ideal
// basis has its generators as rows in coordinates of its order.
//
// Local mechanics that differ from the reference without changing the
// lattice:
// - The final ideal bases are normalized with FLINT's `fmpz_mat_hnf` through
//   `Ideal::set_hnf`, Silex's existing ideal normalization, instead of the
//   reference's `ideal(R, M)` and
//   `_hnf_modular_eldiv(M, minimum(A), :lowerleft)`.  The modulus there only
//   bounds intermediate entries; the lattice spanned is the same.
// - The integer left kernel in `contract` is read from the transform of
//   `fmpz_mat_hnf_transform` instead of Nemo's `kernel(H, side = :left)`.
//   Both give a Z-basis of the same kernel lattice.
// - The reference's zero-ideal branch in `contract` has no counterpart: a
//   Silex ideal always has a full-rank HNF.  Its `is_prime` flag propagation
//   has no counterpart either, since Silex ideals carry no such flag.

#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>

#include "order_internal.hpp"

namespace silex::detail {

namespace {

bool orders_share_field(const Order& left, const Order& right) noexcept {
    return left.is_defined() && right.is_defined() && left.has_basis() &&
           right.has_basis() && left.parent() != nullptr &&
           right.parent() != nullptr &&
           left.parent()->has_same_data(*right.parent()) &&
           left.degree() == right.degree() && left.degree() > 0;
}

// `basis_matrix(sub) * basis_mat_inv(over)`: the basis of `sub` in
// coordinates of `over`.  Fails unless the result is integral, which is the
// reference's containment check `isone(t.den)` / `@assert M.den == 1`.
bool suborder_basis_in_overorder(flint::FmpzMat& out,
                                 const Order& sub,
                                 const Order& over) noexcept {
    const slong n = sub.degree();
    flint::FmpqMat sub_basis(n, n);
    flint::FmpqMat over_basis(n, n);
    flint::FmpqMat over_basis_inv(n, n);
    flint::FmpqMat product(n, n);
    flint::FmpzMat integral(n, n);
    if (!sub.get_basis(flint::FmpqMatRef(sub_basis)) ||
        !over.get_basis(flint::FmpqMatRef(over_basis)) ||
        ::fmpq_mat_inv(over_basis_inv.raw(), over_basis.raw()) == 0) {
        return false;
    }
    ::fmpq_mat_mul(product.raw(), sub_basis.raw(), over_basis_inv.raw());
    if (::fmpq_mat_get_fmpz_mat(integral.raw(), product.raw()) == 0) {
        return false;
    }
    out.swap(integral);
    return true;
}

}  // namespace

bool order_conductor(Ideal& out,
                     const Order& order,
                     const Order& overorder) noexcept {
    if (!out.is_defined() || out.parent() == nullptr ||
        !out.parent()->has_same_data(order) ||
        !orders_share_field(order, overorder)) {
        return false;
    }

    const slong n = order.degree();

    // t = basis_matrix(R) * basis_mat_inv(S); R must lie in S.
    flint::FmpzMat t(n, n);
    if (!suborder_basis_in_overorder(t, order, overorder)) {
        return false;
    }

    // basis_mat_R_in_S_inv_num, d = pseudo_inv(t.num)
    flint::FmpzMat t_inv_num(n, n);
    flint::Fmpz d;
    if (::fmpz_mat_inv(t_inv_num.raw(), d.raw(), t.raw()) == 0) {
        return false;
    }

    // Row k*n + i of the multiplication table holds the S-coordinates of
    // b_k * b_i, so its block k is representation_matrix(b_k) in S.
    Order over(overorder);
    flint::FmpzMat table(n * n, n);
    if (!over.multiplication_table(flint::FmpzMatRef(table))) {
        return false;
    }

    // M[(k - 1)*n + i, j] = N[j, i] with
    // N = representation_matrix(b_k) * basis_mat_R_in_S_inv_num.
    flint::FmpzMat m(n * n, n);
    flint::Fmpz entry;
    for (slong k = 0; k < n; ++k) {
        for (slong row = 0; row < n; ++row) {
            for (slong col = 0; col < n; ++col) {
                // N[row, col]
                ::fmpz_zero(entry.raw());
                for (slong l = 0; l < n; ++l) {
                    ::fmpz_addmul(entry.raw(),
                                  fmpz_mat_entry(table.raw(), k * n + row, l),
                                  fmpz_mat_entry(t_inv_num.raw(), l, col));
                }
                ::fmpz_set(fmpz_mat_entry(m.raw(), k * n + col, row),
                           entry.raw());
            }
        }
    }

    // hnf!(M); H = sub(M, 1:n, 1:n)
    flint::FmpzMat m_hnf(n * n, n);
    ::fmpz_mat_hnf(m_hnf.raw(), m.raw());
    flint::FmpzMat h_transpose(n, n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            ::fmpz_set(fmpz_mat_entry(h_transpose.raw(), j, i),
                       fmpz_mat_entry(m_hnf.raw(), i, j));
        }
    }

    // Hinv, new_den = pseudo_inv(transpose(H)); Hinv = Hinv * t_inv_num
    flint::FmpzMat h_inv(n, n);
    flint::Fmpz new_den;
    if (::fmpz_mat_inv(h_inv.raw(), new_den.raw(), h_transpose.raw()) == 0) {
        return false;
    }
    flint::FmpzMat conductor_rows(n, n);
    ::fmpz_mat_mul(conductor_rows.raw(), h_inv.raw(), t_inv_num.raw());

    // divexact(Hinv, new_den): the conductor lies in R, so this is exact.
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz* value = fmpz_mat_entry(conductor_rows.raw(), i, j);
            if (::fmpz_divisible(value, new_den.raw()) == 0) {
                return false;
            }
            ::fmpz_divexact(value, value, new_den.raw());
        }
    }

    // ideal(R, ...)
    Ideal candidate(*out.parent());
    if (!candidate.set_hnf(flint::FmpzMatConstRef(conductor_rows))) {
        return false;
    }
    out.swap(candidate);
    return true;
}

bool ideal_contract(Ideal& out,
                    const Ideal& ideal,
                    const Order& suborder) noexcept {
    if (!out.is_defined() || out.parent() == nullptr ||
        !out.parent()->has_same_data(suborder) || !ideal.is_defined() ||
        ideal.parent() == nullptr || !ideal.has_hnf() ||
        !orders_share_field(suborder, *ideal.parent())) {
        return false;
    }

    // if order(A) === O: return A
    if (ideal.parent()->has_same_data(suborder)) {
        return out.set(ideal);
    }

    const Order& overorder = *ideal.parent();
    const slong n = suborder.degree();

    // M = basis_matrix(O) * basis_mat_inv(order(A)); @assert M.den == 1
    flint::FmpzMat t(n, n);
    flint::FmpzMat ideal_basis(n, n);
    if (!suborder_basis_in_overorder(t, suborder, overorder) ||
        !ideal.get_hnf(flint::FmpzMatRef(ideal_basis))) {
        return false;
    }

    // H = vcat(basis_matrix(A), M.num); K = kernel(H, side = :left).
    // With U * H = HNF(H) and U unimodular, the last n rows of U are a
    // Z-basis of the left kernel because H has rank n.
    flint::FmpzMat stacked(2 * n, n);
    flint::FmpzMat stacked_hnf(2 * n, n);
    flint::FmpzMat transform(2 * n, 2 * n);
    ::fmpz_mat_concat_vertical(stacked.raw(), ideal_basis.raw(), t.raw());
    ::fmpz_mat_hnf_transform(stacked_hnf.raw(), transform.raw(),
                             stacked.raw());
    for (slong i = n; i < 2 * n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (::fmpz_is_zero(fmpz_mat_entry(stacked_hnf.raw(), i, j)) == 0) {
                return false;
            }
        }
    }

    // M = view(K, 1:d, 1:d) * basis_matrix(A)
    flint::FmpzMat kernel_left(n, n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            ::fmpz_set(fmpz_mat_entry(kernel_left.raw(), i, j),
                       fmpz_mat_entry(transform.raw(), n + i, j));
        }
    }
    flint::FmpzMat in_overorder(n, n);
    ::fmpz_mat_mul(in_overorder.raw(), kernel_left.raw(), ideal_basis.raw());

    // M = M * basis_matrix(order(A)) * basis_mat_inv(O); @assert M.den == 1
    flint::FmpqMat over_basis(n, n);
    flint::FmpqMat sub_basis(n, n);
    flint::FmpqMat sub_basis_inv(n, n);
    flint::FmpqMat rows(n, n);
    flint::FmpqMat scratch(n, n);
    flint::FmpzMat contracted(n, n);
    if (!overorder.get_basis(flint::FmpqMatRef(over_basis)) ||
        !suborder.get_basis(flint::FmpqMatRef(sub_basis)) ||
        ::fmpq_mat_inv(sub_basis_inv.raw(), sub_basis.raw()) == 0) {
        return false;
    }
    ::fmpq_mat_set_fmpz_mat(rows.raw(), in_overorder.raw());
    ::fmpq_mat_mul(scratch.raw(), rows.raw(), over_basis.raw());
    ::fmpq_mat_mul(rows.raw(), scratch.raw(), sub_basis_inv.raw());
    if (::fmpq_mat_get_fmpz_mat(contracted.raw(), rows.raw()) == 0) {
        return false;
    }

    // ideal(O, HNF(M))
    Ideal candidate(*out.parent());
    if (!candidate.set_hnf(flint::FmpzMatConstRef(contracted))) {
        return false;
    }
    out.swap(candidate);
    return true;
}

}  // namespace silex::detail
