#include <silex/ideal.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/lat.hpp>

#include "../order/order_internal.hpp"
#include "../order_element/order_element_internal.hpp"
#include "hnf_membership_internal.hpp"
#include "ideal_access_internal.hpp"
#include "prime_inverse_internal.hpp"
#include "two_generator_internal.hpp"

#include <cassert>
#include <utility>
#include <vector>

namespace silex {
namespace detail {

bool is_full_rank_upper_triangular(const fmpz_mat_t basis) noexcept {
    const slong n = fmpz_mat_nrows(basis);
    if (n <= 0 || fmpz_mat_ncols(basis) != n) {
        return false;
    }
    for (slong i = 0; i < n; ++i) {
        if (fmpz_is_zero(fmpz_mat_entry(basis, i, i)) != 0) {
            return false;
        }
        for (slong j = 0; j < i; ++j) {
            if (fmpz_is_zero(fmpz_mat_entry(basis, i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

namespace {

// The pivots B_jj of the triangular B, stored as their two's-complement words
// in pivots[0, n), and d = |det B| = prod_j |B_jj| in `determinant`, when
// every pivot is a single signed word and d fits in a word.
bool triangular_pivots_word(ulong& determinant,
                            ulong* pivots,
                            const fmpz_mat_t basis) noexcept {
    const slong n = fmpz_mat_ncols(basis);
    ulong product = 1;
    for (slong j = 0; j < n; ++j) {
        const fmpz* pivot = fmpz_mat_entry(basis, j, j);
        if (fmpz_fits_si(pivot) == 0) {
            return false;
        }
        const slong value = fmpz_get_si(pivot);
        if (value == WORD_MIN) {
            return false;
        }
        const ulong magnitude =
                static_cast<ulong>(value < 0 ? -value : value);
        if (n_mul_checked(&product, product, magnitude) != 0) {
            return false;
        }
        pivots[j] = static_cast<ulong>(value);
    }
    determinant = product;
    return true;
}

// The forward substitution of `triangular_basis_contains` on residues modulo
// d = |det B| > 1.  Since adj(B) B = det(B) I, d Z^n lies in the row lattice
// L of B (the precondition of the HNF modulo D, Cohen, GTM 138,
// Section 2.4.2, the remark before Algorithm 2.4.6, which is also the
// precondition of Algorithm 2.4.8), so v and v mod d have the same membership.
// Each step subtracts an integer multiple of row j, which keeps membership,
// and clears column j modulo d.  Divisibility of the residue r_j by
// b_j = |B_jj| is well defined modulo d because b_j divides d; when it fails,
// lifting the residual to an integer vector that is zero in columns < j shows
// that the exact substitution would stop at column j as well.  When every
// column clears, the residual lies in d Z^n, hence in L.  `residue` has n
// entries of scratch space.
bool triangular_basis_contains_mod_det(const fmpz_mat_t basis,
                                       const fmpz* vector,
                                       ulong determinant,
                                       const ulong* pivots,
                                       ulong* residue) noexcept {
    const slong n = fmpz_mat_ncols(basis);
    nmod_t mod{};
    nmod_init(&mod, determinant);
    for (slong l = 0; l < n; ++l) {
        residue[l] = fmpz_fdiv_ui(vector + l, determinant);
    }
    for (slong j = 0; j < n; ++j) {
        const ulong r = residue[j];
        if (r == 0) {
            continue;
        }
        const slong pivot = static_cast<slong>(pivots[j]);
        const ulong b = static_cast<ulong>(pivot < 0 ? -pivot : pivot);
        if (r % b != 0) {
            return false;
        }
        // x_j = r_j / B_jj as a residue modulo d.
        ulong quotient = r / b;
        if (pivot < 0) {
            quotient = nmod_neg(quotient, mod);
        }
        for (slong l = j + 1; l < n; ++l) {
            const fmpz* entry = fmpz_mat_entry(basis, j, l);
            if (fmpz_is_zero(entry) != 0) {
                continue;
            }
            residue[l] = nmod_sub(
                    residue[l],
                    nmod_mul(quotient, fmpz_fdiv_ui(entry, determinant), mod),
                    mod);
        }
    }
    return true;
}

enum class WordMembership { member, nonmember, not_word };

// The word-size half of `triangular_basis_contains`: decides membership when
// every pivot is a single word and |det B| fits in a word, reading `vector`
// only.  Bases up to kStackDegree use stack scratch instead of a heap vector.
WordMembership triangular_basis_contains_word(const fmpz_mat_t basis,
                                              const fmpz* vector) noexcept {
    constexpr slong kStackDegree = 32;
    const slong n = fmpz_mat_ncols(basis);
    ulong stack_scratch[2 * kStackDegree];
    std::vector<ulong> heap_scratch;
    ulong* scratch = stack_scratch;
    if (n > kStackDegree) {
        heap_scratch.resize(static_cast<std::size_t>(2 * n));
        scratch = heap_scratch.data();
    }
    // pivots in [0, n), residues in [n, 2 n).
    ulong determinant = 0;
    if (!triangular_pivots_word(determinant, scratch, basis)) {
        return WordMembership::not_word;
    }
    // A unimodular B spans all of Z^n.
    return determinant == 1 ||
                   triangular_basis_contains_mod_det(basis, vector,
                                                     determinant, scratch,
                                                     scratch + n)
           ? WordMembership::member
           : WordMembership::nonmember;
}

// The exact forward substitution, using `vector` as the residual.
bool triangular_basis_contains_exact(const fmpz_mat_t basis,
                                     fmpz* vector) noexcept {
    const slong n = fmpz_mat_ncols(basis);
    flint::Fmpz quotient;
    for (slong j = 0; j < n; ++j) {
        if (fmpz_is_zero(vector + j) != 0) {
            continue;
        }
        const fmpz* pivot = fmpz_mat_entry(basis, j, j);
        if (fmpz_divisible(vector + j, pivot) == 0) {
            return false;
        }
        fmpz_divexact(quotient.raw(), vector + j, pivot);
        for (slong l = j; l < n; ++l) {
            fmpz_submul(vector + l, quotient.raw(),
                        fmpz_mat_entry(basis, j, l));
        }
    }
    return true;
}

}  // namespace

// B is invertible, so v lies in the row lattice of B iff the unique rational
// solution x of x B = v is integral.  Because B is upper triangular, column j
// of x B = v reads v_j = sum_{i <= j} x_i B_ij, so x_0, x_1, ... are
// determined one at a time by dividing the reduced v_j by B_jj; v is in the
// lattice iff every such division is exact.  This is the coordinate test of
// reference `NumFieldOrd/NfOrd/Ideal/Ideal.jl:containment_by_matrices` (x in I
// iff coordinates(x) * basis_mat_inv(I) is integral), evaluated by forward
// substitution instead of a stored inverse.  When |det B| fits in a word the
// same substitution runs on residues modulo |det B| (see
// `triangular_basis_contains_mod_det`); the answer is the same.
bool triangular_basis_contains(const fmpz_mat_t basis, fmpz* vector) noexcept {
    const WordMembership word =
            triangular_basis_contains_word(basis, vector);
    if (word != WordMembership::not_word) {
        return word == WordMembership::member;
    }
    return triangular_basis_contains_exact(basis, vector);
}

bool triangular_basis_contains_row(const fmpz_mat_t basis,
                                   const fmpz* vector) noexcept {
    const WordMembership word =
            triangular_basis_contains_word(basis, vector);
    if (word != WordMembership::not_word) {
        return word == WordMembership::member;
    }
    const slong n = fmpz_mat_ncols(basis);
    flint::FmpzMat residual(1, n);
    _fmpz_vec_set(fmpz_mat_row(residual.raw(), 0), vector, n);
    return triangular_basis_contains_exact(basis,
                                           fmpz_mat_row(residual.raw(), 0));
}

}  // namespace detail

namespace {

bool fmpz_mat_shape(flint::FmpzMatConstRef matrix,
                    slong rows,
                    slong cols) noexcept {
    return matrix.raw() != nullptr && fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

bool fmpz_mat_shape(flint::FmpzMatRef matrix,
                    slong rows,
                    slong cols) noexcept {
    return matrix.raw() != nullptr && fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

void copy_row(flint::FmpzMatRef out,
              slong out_row,
              flint::FmpzMatConstRef in,
              slong in_row) noexcept {
    for (slong j = 0; j < fmpz_mat_ncols(in.raw()); ++j) {
        fmpz_set(fmpz_mat_entry(out.raw(), out_row, j),
                 fmpz_mat_entry(in.raw(), in_row, j));
    }
}

void multiplication_from_coordinates(
        flint::FmpzMatRef out, flint::FmpzMatConstRef element_coords,
        flint::FmpzMatConstRef multiplication_table) noexcept {
    const slong n = fmpz_mat_ncols(out.raw());
    fmpz_mat_zero(out.raw());
    for (slong basis_index = 0; basis_index < n; ++basis_index) {
        for (slong coeff_index = 0; coeff_index < n; ++coeff_index) {
            const fmpz* coeff = fmpz_mat_entry(
                    element_coords.raw(), 0, coeff_index);
            if (fmpz_is_zero(coeff) != 0) {
                continue;
            }

            const slong table_row = coeff_index * n + basis_index;
            for (slong out_col = 0; out_col < n; ++out_col) {
                fmpz_addmul(fmpz_mat_entry(out.raw(), basis_index, out_col),
                            coeff,
                            fmpz_mat_entry(multiplication_table.raw(),
                                           table_row, out_col));
            }
        }
    }
}

bool order_element_multiplication_matrix(
        flint::FmpzMatRef out,
        Order& order,
        const OrderElement& element) noexcept {
    const Order* element_parent = element.parent();
    const slong n = order.degree();
    if (element_parent == nullptr ||
        !element_parent->has_same_data(order) ||
        !fmpz_mat_shape(out, n, n)) return false;
    flint::FmpzMat element_coords(1, n);
    // A view of the order's cached table; `order` is not reassigned below.
    const flint::FmpzMat* multiplication_table =
            detail::OrderAccess::multiplication_table(order);
    if (!element.get_coordinates(flint::FmpzMatRef(element_coords)) ||
        multiplication_table == nullptr) return false;
    multiplication_from_coordinates(out, flint::FmpzMatConstRef(element_coords),
                                    flint::FmpzMatConstRef(*multiplication_table));
    return true;
}

void set_identity_row(flint::FmpzMatRef row, slong index) noexcept {
    fmpz_mat_zero(row.raw());
    fmpz_one(fmpz_mat_entry(row.raw(), 0, index));
}

// True when `hnf` is an n x n full-rank upper triangular matrix for the
// ideal's degree n, the shape on which the forward substitution of
// `detail::triangular_basis_contains` reads and writes only n entries of an
// n-entry row.
bool is_triangular_basis_of_degree(const flint::FmpzMat& hnf,
                                   slong degree) noexcept {
    return fmpz_mat_ncols(hnf.raw()) == degree &&
           detail::is_full_rank_upper_triangular(hnf.raw());
}

bool hnf_lattice_contains(flint::FmpzMatConstRef hnf,
                          flint::FmpzMatConstRef row) noexcept {
    lat::Lat lattice(fmpz_mat_ncols(hnf.raw()));
    if (!lattice.set_basis(hnf)) {
        return false;
    }
    return lattice.contains_row(row, 0);
}

bool hnf_is_unit_sum(flint::FmpzMatConstRef hnf, slong degree) noexcept {
    if (hnf.raw() == nullptr ||
        fmpz_mat_nrows(hnf.raw()) != 2 * degree ||
        fmpz_mat_ncols(hnf.raw()) != degree) {
        return false;
    }

    for (slong i = 0; i < fmpz_mat_nrows(hnf.raw()); ++i) {
        for (slong j = 0; j < degree; ++j) {
            if (i < degree) {
                if (fmpz_equal_si(fmpz_mat_entry(hnf.raw(), i, j),
                                  i == j ? 1 : 0) == 0) {
                    return false;
                }
            } else if (fmpz_is_zero(fmpz_mat_entry(hnf.raw(), i, j)) == 0) {
                return false;
            }
        }
    }

    return true;
}

bool normalize_full_rank_hnf(flint::FmpzMatRef out,
                             flint::FmpzMatConstRef input) noexcept {
    if (input.raw() == nullptr ||
        fmpz_mat_nrows(input.raw()) != fmpz_mat_ncols(input.raw()) ||
        !fmpz_mat_shape(out, fmpz_mat_nrows(input.raw()),
                        fmpz_mat_ncols(input.raw()))) {
        return false;
    }

    fmpz_mat_hnf(out.raw(), input.raw());

    flint::Fmpz det;
    fmpz_mat_det(det.raw(), out.raw());
    return !flint::fmpz_is_zero(det) &&
           fmpz_mat_is_in_hnf(out.raw()) != 0;
}

bool hnf_is_ideal(const Order& order, flint::FmpzMatConstRef hnf) noexcept {
    const slong n = order.degree();
    if (!fmpz_mat_shape(hnf, n, n)) {
        return false;
    }

    OrderElement gen(order);
    OrderElement omega(order);
    OrderElement product(order);
    flint::FmpzMat row(1, n);
    flint::FmpzMat omega_row(1, n);
    lat::Lat lattice(n);
    if (!gen.is_defined() || !omega.is_defined() || !product.is_defined() ||
        !lattice.set_basis(hnf)) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        copy_row(flint::FmpzMatRef(row), 0, hnf, i);
        if (!gen.set_coordinates(flint::FmpzMatConstRef(row))) {
            return false;
        }

        for (slong j = 0; j < n; ++j) {
            set_identity_row(flint::FmpzMatRef(omega_row), j);
            if (!omega.set_coordinates(flint::FmpzMatConstRef(omega_row)) ||
                !product.multiply(gen, omega) ||
                !product.get_coordinates(flint::FmpzMatRef(row)) ||
                !lattice.contains_row(flint::FmpzMatConstRef(row), 0)) {
                return false;
            }
        }
    }

    return true;
}

void common_content(flint::FmpzRef out,
                    flint::FmpzMatConstRef hnf,
                    flint::FmpzConstRef den) noexcept {
    flint::fmpz_set(out, den);
    flint::Fmpz abs_entry;
    for (slong i = 0; i < fmpz_mat_nrows(hnf.raw()); ++i) {
        for (slong j = 0; j < fmpz_mat_ncols(hnf.raw()); ++j) {
            flint::fmpz_abs(flint::FmpzRef(abs_entry),
                            flint::FmpzConstRef(
                                    fmpz_mat_entry(hnf.raw(), i, j)));
            flint::fmpz_gcd(out, flint::FmpzConstRef(out.raw()),
                            flint::FmpzConstRef(abs_entry));
        }
    }
}

bool multiply_element_by_fmpz(Element& out,
                              const Element& input,
                              flint::FmpzConstRef scalar) noexcept {
    if (out.parent() == nullptr || input.parent() == nullptr ||
        !out.parent()->has_same_data(*input.parent())) {
        return false;
    }

    flint::FmpqPoly polynomial;
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    fmpq_poly_scalar_mul_fmpz(polynomial.raw(), polynomial.raw(), scalar.raw());
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool scale_integral_ideal(Ideal& out,
                          const Ideal& input,
                          flint::FmpzConstRef scalar) noexcept {
    if (out.parent() == nullptr || input.parent() == nullptr ||
        !out.parent()->has_same_data(*input.parent()) || !input.has_hnf() ||
        flint::fmpz_sgn(scalar) <= 0) {
        return false;
    }

    const slong n = input.degree();
    flint::FmpzMat hnf(n, n);
    if (!input.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }
    fmpz_mat_scalar_mul_fmpz(hnf.raw(), hnf.raw(), scalar.raw());
    // `input` is an Ideal, so its lattice I is an ideal, and sI is closed
    // under the order: w(sx) = s(wx) lies in sI for x in I and w in O.  The
    // hnf_is_ideal product test of set_hnf is therefore skipped.
    return detail::IdealAccess::set_known_ideal_hnf(
            out, flint::FmpzMatConstRef(hnf));
}

bool normalize_lattice_num_den(flint::FmpzMat& numerator,
                               flint::Fmpz& den) noexcept {
    if (flint::fmpz_sgn(flint::FmpzConstRef(den)) <= 0) {
        return false;
    }

    const slong n = fmpz_mat_ncols(numerator.raw());
    lat::Lat lattice(n);
    lat::Lat hnf_lattice(n);
    if (!lattice.set_basis(flint::FmpzMatConstRef(numerator)) ||
        !lattice.hnf(hnf_lattice) || hnf_lattice.nrows() != n) {
        return false;
    }

    flint::FmpzMat normalized(n, n);
    flint::Fmpz content;
    fmpz_mat_set(normalized.raw(), hnf_lattice.basis_ref().raw());
    common_content(flint::FmpzRef(content),
                   flint::FmpzMatConstRef(normalized),
                   flint::FmpzConstRef(den));
    if (!flint::fmpz_is_one(content)) {
        fmpz_mat_scalar_divexact_fmpz(normalized.raw(), normalized.raw(),
                                      content.raw());
        flint::fmpz_divexact(flint::FmpzRef(den),
                             flint::FmpzConstRef(den),
                             flint::FmpzConstRef(content));
    }

    numerator.swap(normalized);
    return true;
}

bool rational_mat_to_lattice(flint::FmpzMat& numerator,
                             flint::Fmpz& den,
                             flint::FmpqMatConstRef matrix) noexcept {
    flint::FmpzMat integral(fmpq_mat_nrows(matrix.raw()),
                            fmpq_mat_ncols(matrix.raw()));
    fmpq_mat_get_fmpz_mat_matwise(integral.raw(), den.raw(), matrix.raw());
    if (!normalize_lattice_num_den(integral, den)) {
        return false;
    }
    numerator = std::move(integral);
    return true;
}

bool intersect_rational_lattices(flint::FmpzMat& numerator,
                                 flint::Fmpz& den,
                                 flint::FmpzMatConstRef left_num,
                                 flint::FmpzConstRef left_den,
                                 flint::FmpzMatConstRef right_num,
                                 flint::FmpzConstRef right_den) noexcept {
    const slong n = fmpz_mat_ncols(left_num.raw());
    if (n != fmpz_mat_ncols(right_num.raw()) ||
        flint::fmpz_sgn(left_den) <= 0 ||
        flint::fmpz_sgn(right_den) <= 0) {
        return false;
    }

    flint::Fmpz common_den;
    flint::Fmpz scale;
    flint::fmpz_lcm(flint::FmpzRef(common_den), left_den, right_den);

    flint::FmpzMat left_scaled(fmpz_mat_nrows(left_num.raw()), n);
    flint::FmpzMat right_scaled(fmpz_mat_nrows(right_num.raw()), n);
    flint::fmpz_divexact(flint::FmpzRef(scale),
                         flint::FmpzConstRef(common_den), left_den);
    fmpz_mat_scalar_mul_fmpz(left_scaled.raw(), left_num.raw(), scale.raw());
    flint::fmpz_divexact(flint::FmpzRef(scale),
                         flint::FmpzConstRef(common_den), right_den);
    fmpz_mat_scalar_mul_fmpz(right_scaled.raw(), right_num.raw(), scale.raw());

    lat::Lat left_lattice(n);
    lat::Lat right_lattice(n);
    lat::Lat intersection(n);
    if (!left_lattice.set_basis(flint::FmpzMatConstRef(left_scaled)) ||
        !right_lattice.set_basis(flint::FmpzMatConstRef(right_scaled)) ||
        !left_lattice.intersection(intersection, right_lattice)) {
        return false;
    }

    flint::FmpzMat out(intersection.nrows(), n);
    fmpz_mat_set(out.raw(), intersection.basis_ref().raw());
    if (!normalize_lattice_num_den(out, common_den)) {
        return false;
    }

    numerator = std::move(out);
    flint::fmpz_set(flint::FmpzRef(den), flint::FmpzConstRef(common_den));
    return true;
}

bool multiplication_matrix_order_element(flint::FmpqMatRef out,
                                         const Order& order,
                                         const OrderElement& multiplier) noexcept {
    const slong n = order.degree();
    if (out.raw() == nullptr || fmpq_mat_nrows(out.raw()) != n ||
        fmpq_mat_ncols(out.raw()) != n ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(order)) {
        return false;
    }

    OrderElement basis_elem(order);
    OrderElement product(order);
    flint::FmpzMat row(1, n);
    flint::FmpzMat coords(1, n);
    if (!basis_elem.is_defined() || !product.is_defined()) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        set_identity_row(flint::FmpzMatRef(row), i);
        if (!basis_elem.set_coordinates(flint::FmpzMatConstRef(row)) ||
            !product.multiply(basis_elem, multiplier) ||
            !product.get_coordinates(flint::FmpzMatRef(coords))) {
            return false;
        }

        for (slong j = 0; j < n; ++j) {
            fmpq_set_fmpz(fmpq_mat_entry(out.raw(), i, j),
                          fmpz_mat_entry(coords.raw(), 0, j));
        }
    }

    return true;
}

bool preimage_under_mul(flint::FmpzMat& numerator,
                        flint::Fmpz& den,
                        flint::FmpzMatConstRef hnf,
                        const Order& order,
                        const OrderElement& multiplier) noexcept {
    const slong n = order.degree();
    flint::FmpqMat matrix(n, n);
    flint::FmpqMat inverse(n, n);
    flint::FmpqMat hnf_q(n, n);
    flint::FmpqMat preimage_q(n, n);

    if (!multiplication_matrix_order_element(flint::FmpqMatRef(matrix),
                                             order, multiplier) ||
        fmpq_mat_inv(inverse.raw(), matrix.raw()) == 0) {
        return false;
    }

    fmpq_mat_set_fmpz_mat(hnf_q.raw(), hnf.raw());
    fmpq_mat_mul(preimage_q.raw(), hnf_q.raw(), inverse.raw());
    return rational_mat_to_lattice(numerator, den,
                                   flint::FmpqMatConstRef(preimage_q));
}

bool colon_one_integral_hnf_inverse(FractionalIdeal& out,
                                          const Ideal& denominator) noexcept {
    const Order* order = out.parent();
    if (order == nullptr || denominator.parent() == nullptr ||
        !order->has_same_data(*denominator.parent()) ||
        !denominator.has_hnf()) {
        return false;
    }

    const slong n = order->degree();
    flint::FmpzMat denominator_hnf(n, n);
    flint::FmpzMat row(1, n);
    flint::FmpqMat multiplication(n, n);
    flint::FmpzMat constraints(n * n, n);
    OrderElement basis_elem(*order);
    if (!basis_elem.is_defined() ||
        !denominator.get_hnf(flint::FmpzMatRef(denominator_hnf))) {
        return false;
    }

    // Source trace: reference NfOrd/Ideal/Ideal.jl colon(a, b) builds
    // representation_matrix(B[i]) for b's basis and returns the dual lattice.
    // For a = O in order coordinates, basis_mat_inv(a) is the identity.
    for (slong i = 0; i < n; ++i) {
        copy_row(flint::FmpzMatRef(row), 0,
                 flint::FmpzMatConstRef(denominator_hnf), i);
        if (!basis_elem.set_coordinates(flint::FmpzMatConstRef(row)) ||
            !multiplication_matrix_order_element(
                    flint::FmpqMatRef(multiplication), *order, basis_elem)) {
            return false;
        }

        for (slong col = 0; col < n; ++col) {
            for (slong r = 0; r < n; ++r) {
                const fmpq* entry =
                        fmpq_mat_entry(multiplication.raw(), r, col);
                if (!fmpz_is_one(fmpq_denref(entry))) {
                    return false;
                }
                fmpz_set(fmpz_mat_entry(constraints.raw(), i * n + col, r),
                         fmpq_numref(entry));
            }
        }
    }

    lat::Lat constraint_lattice(n);
    lat::Lat hnf_lattice(n);
    if (!constraint_lattice.set_basis(
                flint::FmpzMatConstRef(constraints)) ||
        !constraint_lattice.hnf(hnf_lattice) ||
        hnf_lattice.nrows() != n) {
        return false;
    }

    flint::FmpqMat hnf_q(n, n);
    flint::FmpqMat hnf_transpose(n, n);
    flint::FmpqMat dual(n, n);
    fmpq_mat_set_fmpz_mat(hnf_q.raw(), hnf_lattice.basis_ref().raw());
    fmpq_mat_transpose(hnf_transpose.raw(), hnf_q.raw());
    if (fmpq_mat_inv(dual.raw(), hnf_transpose.raw()) == 0) {
        return false;
    }

    flint::FmpzMat numerator_hnf(n, n);
    flint::Fmpz denominator_value;
    Ideal numerator(*order);
    FractionalIdeal candidate(*order);
    // The dual lattice is {x : <x, c> in Z for every constraint row c}, that
    // is {x : x b_i in O for every basis element b_i of the denominator},
    // the colon (O : I).  It is closed under the order because O is a ring
    // (w x b_i lies in wO, inside O, for w in O), so numerator_hnf, which
    // spans denominator_value * (O : I) inside O, is an ideal, and the
    // hnf_is_ideal product test of set_hnf is skipped.
    if (!rational_mat_to_lattice(numerator_hnf, denominator_value,
                                 flint::FmpqMatConstRef(dual)) ||
        !detail::IdealAccess::set_known_ideal_hnf(
                numerator, flint::FmpzMatConstRef(numerator_hnf)) ||
        !candidate.set_integral_den(numerator,
                                    flint::FmpzConstRef(denominator_value))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

}  // namespace

Ideal::Ideal(const Order& parent) noexcept {
    define(parent);
}

Ideal::~Ideal() noexcept {
    clear();
}

Ideal::Ideal(Ideal&& other) noexcept {
    swap(other);
}

Ideal& Ideal::operator=(Ideal&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void Ideal::swap(Ideal& other) noexcept {
    parent_.swap(other.parent_);
    hnf_.swap(other.hnf_);
    principal_generator_.swap(other.principal_generator_);
    norm_.swap(other.norm_);
    std::swap(has_hnf_, other.has_hnf_);
    std::swap(has_principal_generator_,
              other.has_principal_generator_);
    std::swap(has_norm_, other.has_norm_);
}

void Ideal::clear() noexcept {
    parent_.clear();
    hnf_ = flint::FmpzMat(0, 0);
    principal_generator_.clear();
    flint::fmpz_zero(flint::FmpzRef(norm_));
    has_hnf_ = false;
    has_principal_generator_ = false;
    has_norm_ = false;
}

bool Ideal::define(const Order& parent) noexcept {
    if (!parent.is_defined() || !parent.has_basis()) {
        return false;
    }

    Ideal next;
    next.parent_ = parent;
    next.hnf_ = flint::FmpzMat(0, parent.degree());
    swap(next);
    return true;
}

bool Ideal::set(const Ideal& other) noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_) || !other.has_hnf_) {
        return false;
    }
    if (this == &other) {
        return true;
    }

    Ideal copy(parent_);
    if (!copy.set_hnf_direct(flint::FmpzMatConstRef(other.hnf_))) {
        return false;
    }
    if (other.has_principal_generator_) {
        if (!copy.principal_generator_.define(parent_) ||
            !copy.principal_generator_.set(other.principal_generator_)) {
            return false;
        }
        copy.has_principal_generator_ = true;
    }
    flint::fmpz_set(flint::FmpzRef(copy.norm_),
                    flint::FmpzConstRef(other.norm_));
    copy.has_norm_ = other.has_norm_;
    swap(copy);
    return true;
}

bool Ideal::is_defined() const noexcept {
    return parent_.is_defined();
}

const Order* Ideal::parent() const noexcept {
    return !is_defined() ? nullptr : &parent_;
}

slong Ideal::degree() const noexcept {
    return !is_defined() ? 0 : parent_.degree();
}

bool Ideal::has_hnf() const noexcept {
    return is_defined() && has_hnf_;
}

bool Ideal::one() noexcept {
    if (!is_defined()) {
        return false;
    }

    const slong n = parent_.degree();
    Ideal candidate(parent_);
    flint::FmpzMat identity(n, n);
    fmpz_mat_one(identity.raw());
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(identity)) ||
        !candidate.principal_generator_.define(parent_) ||
        !candidate.principal_generator_.one()) {
        return false;
    }
    candidate.has_principal_generator_ = true;
    flint::fmpz_one(flint::FmpzRef(candidate.norm_));
    candidate.has_norm_ = true;
    swap(candidate);
    return true;
}

bool Ideal::is_one() const noexcept {
    return is_defined() && has_hnf_ &&
           fmpz_mat_is_one(hnf_.raw()) != 0;
}

bool Ideal::set_hnf(flint::FmpzMatConstRef hnf) noexcept {
    return detail::IdealAccess::store_hnf(
        *this, hnf, detail::IdealAccess::HnfInput::unnormalized,
        detail::IdealAccess::IdealTest::always);
}

bool Ideal::get_hnf(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !has_hnf_ ||
        !fmpz_mat_shape(out, parent_.degree(), parent_.degree())) {
        return false;
    }
    fmpz_mat_set(out.raw(), hnf_.raw());
    return true;
}

bool Ideal::set_principal(const OrderElement& generator) noexcept {
    return detail::IdealAccess::set_principal(*this, generator, nullptr);
}

const OrderElement* Ideal::known_principal_generator() const noexcept {
    return is_defined() && has_hnf_ && has_principal_generator_ &&
                           principal_generator_.is_defined()
                   ? &principal_generator_
                   : nullptr;
}

bool Ideal::contains(const OrderElement& element) const noexcept {
    if (!is_defined() || !has_hnf_ ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(parent_)) {
        return false;
    }

    // The element's stored coordinates are read in place; a copy is made
    // only for the fallbacks that need a residual or a matrix argument.
    const fmpz_mat_struct* coordinates =
            detail::OrderElementAccess::coordinates(element);
    if (coordinates == nullptr || fmpz_mat_nrows(coordinates) != 1 ||
        fmpz_mat_ncols(coordinates) != parent_.degree()) {
        return false;
    }
    // Every construction route is expected to store a full-rank HNF, which
    // is upper triangular with a nonzero diagonal; membership is then decided
    // by forward substitution against the stored matrix, with no lattice
    // copy. Both routes decide membership in the Z-row span of `hnf_`. The
    // lattice test below is a defensive fallback for any other stored shape;
    // it is not reachable through the API today.  Tests reach it through the
    // noninstalled detail::IdealAccess hook.
    if (is_triangular_basis_of_degree(hnf_, parent_.degree())) {
        return detail::triangular_basis_contains_row(
                hnf_.raw(), fmpz_mat_row(coordinates, 0));
    }
    return hnf_lattice_contains(flint::FmpzMatConstRef(hnf_),
                                flint::FmpzMatConstRef(coordinates));
}

bool Ideal::contains(const Ideal& ideal) const noexcept {
    if (!is_defined() || !ideal.is_defined() ||
        !parent_.has_same_data(ideal.parent_) || !has_hnf_ ||
        !ideal.has_hnf_) {
        return false;
    }

    const slong n = degree();
    // As in contains(OrderElement): I is contained in this ideal iff each of
    // the n rows of I's stored HNF lies in the Z-row span of `hnf_`, decided
    // by forward substitution against `hnf_` with one reused residual row
    // instead of a lattice copy.  The lattice test below decides the same
    // span membership and is a defensive fallback for any other stored
    // shape; it is not expected to be reachable through the API.  Tests
    // reach it through the noninstalled detail::IdealAccess hook.
    if (is_triangular_basis_of_degree(hnf_, n) &&
        fmpz_mat_nrows(ideal.hnf_.raw()) == n &&
        fmpz_mat_ncols(ideal.hnf_.raw()) == n) {
        flint::FmpzMat residual(1, n);
        fmpz* row = fmpz_mat_row(residual.raw(), 0);
        for (slong i = 0; i < n; ++i) {
            _fmpz_vec_set(row, fmpz_mat_row(ideal.hnf_.raw(), i), n);
            if (!detail::triangular_basis_contains(hnf_.raw(), row)) {
                return false;
            }
        }
        return true;
    }

    lat::Lat lattice(n);
    if (!lattice.set_basis(flint::FmpzMatConstRef(hnf_))) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        if (!lattice.contains_row(flint::FmpzMatConstRef(ideal.hnf_), i)) {
            return false;
        }
    }
    return true;
}

bool Ideal::norm(flint::FmpzRef out) const noexcept {
    if (!is_defined() || !has_hnf_) {
        return false;
    }

    if (!has_norm_) {
        fmpz_mat_det(norm_.raw(), hnf_.raw());
        flint::fmpz_abs(flint::FmpzRef(norm_),
                        flint::FmpzConstRef(norm_));
        has_norm_ = true;
    }
    flint::fmpz_set(out, flint::FmpzConstRef(norm_));
    return true;
}

bool Ideal::equal(const Ideal& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.has_same_data(other.parent_) && has_hnf_ &&
           other.has_hnf_ &&
           fmpz_mat_equal(hnf_.raw(), other.hnf_.raw()) != 0;
}

bool Ideal::is_coprime(bool& result, const Ideal& other) const noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_) || !has_hnf_ ||
        !other.has_hnf_) {
        return false;
    }

    Ideal sum(parent_);
    if (!sum.add(*this, other)) {
        return false;
    }

    result = sum.is_one();
    return true;
}

bool Ideal::add_to_one(OrderElement& left_witness,
                       OrderElement& right_witness,
                       const Ideal& other) const noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_) || !has_hnf_ ||
        !other.has_hnf_ || left_witness.parent() == nullptr ||
        !left_witness.parent()->has_same_data(parent_) ||
        right_witness.parent() == nullptr ||
        !right_witness.parent()->has_same_data(parent_) ||
        &left_witness == &right_witness) {
        return false;
    }

    const slong n = parent_.degree();
    flint::FmpzMat stacked(2 * n, n);
    flint::FmpzMat hnf(2 * n, n);
    flint::FmpzMat transform(2 * n, 2 * n);
    flint::FmpzMat one_row(1, n);
    flint::FmpzMat left_coeff(1, n);
    flint::FmpzMat right_coeff(1, n);
    flint::FmpzMat left_row(1, n);
    flint::FmpzMat right_row(1, n);
    OrderElement one(parent_);
    OrderElement left_candidate(parent_);
    OrderElement right_candidate(parent_);
    OrderElement sum(parent_);
    OrderElement publish_left(parent_);
    OrderElement publish_right(parent_);

    if (!one.is_defined() || !left_candidate.is_defined() ||
        !right_candidate.is_defined() || !sum.is_defined() ||
        !publish_left.is_defined() || !publish_right.is_defined()) {
        return false;
    }

    fmpz_mat_concat_vertical(stacked.raw(), hnf_.raw(),
                             other.hnf_.raw());
    fmpz_mat_hnf_transform(hnf.raw(), transform.raw(), stacked.raw());
    if (!hnf_is_unit_sum(flint::FmpzMatConstRef(hnf), n) ||
        !one.one() || !one.get_coordinates(flint::FmpzMatRef(one_row))) {
        return false;
    }

    for (slong j = 0; j < n; ++j) {
        for (slong i = 0; i < n; ++i) {
            fmpz_addmul(fmpz_mat_entry(left_coeff.raw(), 0, j),
                        fmpz_mat_entry(one_row.raw(), 0, i),
                        fmpz_mat_entry(transform.raw(), i, j));
            fmpz_addmul(fmpz_mat_entry(right_coeff.raw(), 0, j),
                        fmpz_mat_entry(one_row.raw(), 0, i),
                        fmpz_mat_entry(transform.raw(), i, n + j));
        }
    }

    fmpz_mat_mul(left_row.raw(), left_coeff.raw(), hnf_.raw());
    fmpz_mat_mul(right_row.raw(), right_coeff.raw(), other.hnf_.raw());

    if (!left_candidate.set_coordinates(flint::FmpzMatConstRef(left_row)) ||
        !right_candidate.set_coordinates(flint::FmpzMatConstRef(right_row)) ||
        !sum.add(left_candidate, right_candidate) || !sum.equal_si(1) ||
        !contains(left_candidate) || !other.contains(right_candidate) ||
        !publish_left.set(left_candidate) ||
        !publish_right.set(right_candidate)) {
        return false;
    }

    left_witness.swap(publish_left);
    right_witness.swap(publish_right);
    return true;
}

bool Ideal::add(const Ideal& left, const Ideal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_hnf_ ||
        !right.has_hnf_) {
        return false;
    }

    const slong n = parent_.degree();
    Ideal candidate(parent_);
    flint::FmpzMat rows(2 * n, n);
    for (slong i = 0; i < n; ++i) {
        copy_row(flint::FmpzMatRef(rows), i,
                 flint::FmpzMatConstRef(left.hnf_), i);
        copy_row(flint::FmpzMatRef(rows), n + i,
                 flint::FmpzMatConstRef(right.hnf_), i);
    }

    // `left` and `right` are Ideals of this order holding an HNF, so their
    // lattices I and J are ideals: every installed route that stores an HNF
    // either tests closure or builds an ideal by construction, and only the
    // test hook IdealAccess::set_hnf_direct stores an unchecked matrix.  The
    // rows span I + J, and for x in I, y in J and w in O, w(x + y) = wx + wy
    // lies in I + J, so the sum is an ideal and the hnf_is_ideal product test
    // of set_rows is skipped (assertion builds still run it).
    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool Ideal::intersect(const Ideal& left, const Ideal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_hnf_ ||
        !right.has_hnf_) {
        return false;
    }

    lat::Lat left_lattice(parent_.degree());
    lat::Lat right_lattice(parent_.degree());
    lat::Lat intersection(parent_.degree());
    if (!left_lattice.set_basis(flint::FmpzMatConstRef(left.hnf_)) ||
        !right_lattice.set_basis(flint::FmpzMatConstRef(right.hnf_)) ||
        !left_lattice.intersection(intersection, right_lattice) ||
        intersection.nrows() != parent_.degree()) {
        return false;
    }

    Ideal candidate(parent_);
    if (!candidate.set_rows(intersection.basis_ref())) {
        return false;
    }
    swap(candidate);
    return true;
}

bool Ideal::multiply(const Ideal& left, const Ideal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_hnf_ ||
        !right.has_hnf_) {
        return false;
    }
    if (left.is_one()) {
        return set(right);
    }
    if (right.is_one()) {
        return set(left);
    }

    if (!parent_.is_maximal() && !left.has_principal_generator_ &&
        !right.has_principal_generator_) return multiply_generic(left, right);
    flint::FmpzMat fast_hnf(parent_.degree(), parent_.degree());
    if (detail::try_ideal_product_hnf(flint::FmpzMatRef(fast_hnf), left, right)) {
        Ideal candidate(parent_);
        if (candidate.set_hnf_direct(flint::FmpzMatConstRef(fast_hnf))) {
            swap(candidate);
            return true;
        }
    }
    // The generic two-generator search can fail (for example, small-norm
    // ideals above a prime that splits completely).  On a maximal order, take
    // the HNF of the basis products modulo the product of the minima of the
    // primitive parts (contents restored by scaling; the reference mul_gen
    // uses the minima of the ideals themselves) instead of the exact HNF that
    // multiply_generic computes; the HNF is the same.
    if (parent_.is_maximal() &&
        detail::basis_product_hnf(flint::FmpzMatRef(fast_hnf), left, right)) {
        Ideal candidate(parent_);
        if (candidate.set_hnf_direct(flint::FmpzMatConstRef(fast_hnf))) {
            swap(candidate);
            return true;
        }
    }
    return multiply_generic(left, right);
}

bool Ideal::multiply_generic(const Ideal& left, const Ideal& right) noexcept {
    const slong n = parent_.degree();
    Ideal candidate(parent_);
    OrderElement product(parent_);
    flint::FmpzMat rows(n * n, n);
    flint::FmpzMat row(1, n);
    std::vector<OrderElement> left_basis;
    std::vector<OrderElement> right_basis;
    left_basis.reserve(static_cast<std::size_t>(n));
    right_basis.reserve(static_cast<std::size_t>(n));
    if (!product.is_defined()) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        left_basis.emplace_back(parent_);
        copy_row(flint::FmpzMatRef(row), 0,
                 flint::FmpzMatConstRef(left.hnf_), i);
        if (!left_basis.back().is_defined() ||
            !left_basis.back().set_coordinates(flint::FmpzMatConstRef(row))) {
            return false;
        }

        right_basis.emplace_back(parent_);
        copy_row(flint::FmpzMatRef(row), 0,
                 flint::FmpzMatConstRef(right.hnf_), i);
        if (!right_basis.back().is_defined() ||
            !right_basis.back().set_coordinates(flint::FmpzMatConstRef(row))) {
            return false;
        }
    }

    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (!product.multiply(left_basis[static_cast<std::size_t>(i)],
                                  right_basis[static_cast<std::size_t>(j)]) ||
                !product.get_coordinates(flint::FmpzMatRef(row))) {
                return false;
            }
            copy_row(flint::FmpzMatRef(rows), i * n + j,
                     flint::FmpzMatConstRef(row), 0);
        }
    }

    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool Ideal::multiplier_ring(Order& out) const noexcept {
    if (!is_defined() || !has_hnf_ || out.parent() == nullptr ||
        parent_.parent() == nullptr ||
        !out.parent()->has_same_data(*parent_.parent()) ||
        out.has_same_data(parent_)) {
        return false;
    }

    const slong n = parent_.degree();
    FractionalIdeal integral(parent_);
    FractionalIdeal colon_result(parent_);
    Ideal numerator(parent_);
    flint::Fmpz den;
    flint::FmpzMat hnf(n, n);
    flint::FmpqMat hnf_q(n, n);
    flint::FmpqMat parent_basis(n, n);
    flint::FmpqMat basis(n, n);
    Order candidate(*parent_.parent());
    flint::Fmpz index;

    if (!integral.set_integral(*this) ||
        !colon_result.colon(integral, integral) ||
        !colon_result.get_integral_den(numerator, flint::FmpzRef(den)) ||
        !numerator.get_hnf(flint::FmpzMatRef(hnf)) ||
        !parent_.get_basis(flint::FmpqMatRef(parent_basis))) {
        return false;
    }

    fmpq_mat_set_fmpz_mat(hnf_q.raw(), hnf.raw());
    fmpq_mat_mul(basis.raw(), hnf_q.raw(), parent_basis.raw());
    fmpq_mat_scalar_div_fmpz(basis.raw(), basis.raw(), den.raw());

    if (!candidate.set_basis(flint::FmpqMatConstRef(basis)) ||
        !order_index(flint::FmpzRef(index), parent_, candidate)) {
        return false;
    }

    candidate.clear_maximality();
    out.swap(candidate);
    return true;
}

bool Ideal::set_hnf_direct(flint::FmpzMatConstRef hnf) noexcept {
    if (!is_defined() ||
        !fmpz_mat_shape(hnf, parent_.degree(), parent_.degree())) {
        return false;
    }

    flint::FmpzMat copy(parent_.degree(), parent_.degree());
    fmpz_mat_set(copy.raw(), hnf.raw());
    hnf_.swap(copy);
    has_hnf_ = true;
    has_principal_generator_ = false;
    has_norm_ = false;
    flint::fmpz_zero(flint::FmpzRef(norm_));
    return true;
}

bool Ideal::set_rows(flint::FmpzMatConstRef rows) noexcept {
    if (!is_defined() || rows.raw() == nullptr ||
        fmpz_mat_ncols(rows.raw()) != parent_.degree()) {
        return false;
    }

    lat::Lat lattice(parent_.degree());
    lat::Lat hnf_lattice(parent_.degree());
    if (!lattice.set_basis(rows) || !lattice.hnf(hnf_lattice) ||
        hnf_lattice.nrows() != parent_.degree()) {
        return false;
    }
    return set_hnf(hnf_lattice.basis_ref());
}

bool Ideal::set_known_ideal_rows(flint::FmpzMatConstRef rows) noexcept {
    if (!is_defined() || rows.raw() == nullptr ||
        fmpz_mat_ncols(rows.raw()) != parent_.degree()) {
        return false;
    }

    lat::Lat lattice(parent_.degree());
    lat::Lat hnf_lattice(parent_.degree());
    if (!lattice.set_basis(rows) || !lattice.hnf(hnf_lattice) ||
        hnf_lattice.nrows() != parent_.degree()) {
        return false;
    }

    // The lattice HNF is already normalized; the caller vouches for closure.
    return detail::IdealAccess::store_hnf(
        *this, hnf_lattice.basis_ref(),
        detail::IdealAccess::HnfInput::normalized,
        detail::IdealAccess::IdealTest::assert_only);
}

namespace detail {

bool IdealAccess::store_hnf(Ideal& ideal, flint::FmpzMatConstRef hnf,
                            HnfInput input, IdealTest test) noexcept {
    if (!ideal.is_defined()) {
        return false;
    }

    const slong n = ideal.parent_.degree();
    flint::FmpzMat normalized(n, n);
    if (input == HnfInput::unnormalized) {
        if (!normalize_full_rank_hnf(flint::FmpzMatRef(normalized), hnf)) {
            return false;
        }
    } else if (!fmpz_mat_shape(hnf, n, n)) {
        return false;
    } else {
        fmpz_mat_set(normalized.raw(), hnf.raw());
    }

    if (test == IdealTest::always) {
        if (!hnf_is_ideal(ideal.parent_, flint::FmpzMatConstRef(normalized))) {
            return false;
        }
    } else {
        assert(hnf_is_ideal(ideal.parent_, flint::FmpzMatConstRef(normalized)));
    }

    Ideal candidate(ideal.parent_);
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(normalized))) {
        return false;
    }
    ideal.swap(candidate);
    return true;
}

namespace {

// True when the diagonal of the n x n matrix `hnf` multiplies to `modulus`.
bool hnf_diagonal_product_equals(flint::FmpzMatConstRef hnf,
                                 const fmpz_t modulus) noexcept {
    flint::Fmpz product;
    fmpz_one(product.raw());
    for (slong i = 0; i < fmpz_mat_nrows(hnf.raw()); ++i) {
        fmpz_mul(product.raw(), product.raw(), fmpz_mat_entry(hnf.raw(), i, i));
    }
    return fmpz_equal(product.raw(), modulus) != 0;
}

// Debug cross-check of the determinant-modular route: `hnf` equals the
// canonical HNF that fmpz_mat_hnf returns for `rows`.
[[maybe_unused]] bool equals_full_hnf(flint::FmpzMatConstRef rows,
                                      flint::FmpzMatConstRef hnf) noexcept {
    flint::FmpzMat full(fmpz_mat_nrows(rows.raw()), fmpz_mat_ncols(rows.raw()));
    fmpz_mat_hnf(full.raw(), rows.raw());
    return fmpz_mat_equal(full.raw(), hnf.raw()) != 0;
}

}  // namespace

bool IdealAccess::set_principal(Ideal& ideal, const OrderElement& generator,
                                const fmpq* norm) noexcept {
    if (!ideal.is_defined() || generator.parent() == nullptr ||
        !generator.parent()->has_same_data(ideal.parent_) ||
        generator.equal_si(0)) {
        return false;
    }

    // Source trace: reference `base4.c:idealmulelt` and `idealhnf_two` form HNF
    // from multiplication-table rows.  For a principal ideal, those rows are
    // the multiplication matrix of the generator in the order basis.
    const slong n = ideal.parent_.degree();
    flint::FmpzMat rows(n, n);
    if (!order_element_multiplication_matrix(flint::FmpzMatRef(rows),
                                             ideal.parent_, generator)) {
        return false;
    }

    Ideal candidate(ideal.parent_);
    bool stored = false;
    if (norm != nullptr && fmpz_is_one(fmpq_denref(norm)) != 0 &&
        fmpz_is_zero(fmpq_numref(norm)) == 0) {
        // The rows are the matrix M of multiplication by the generator g in
        // the order's Z-basis, so |det M| = |N_{K/Q}(g)| = [O : gO], the
        // known norm.  adj(M) M = det(M) I puts det(M) Z^n inside the row
        // lattice, so |N(g)| is a positive multiple of the determinant of the
        // full-rank row lattice: the precondition of FLINT 3.6.0
        // fmpz_mat_hnf_modular (src/fmpz_mat/hnf_modular.c), the HNF modulo D
        // of Domich, Kannan and Trotter (Math. Oper. Res. 12, 1987) and Cohen,
        // "A Course in Computational Algebraic Number Theory", Algorithm
        // 2.4.8.  FLINT documents its output as the unique row HNF of the
        // input, the canonical form fmpz_mat_hnf returns, so the stored HNF is
        // the one the full-HNF route below stores.  A diagonal product other
        // than the modulus means the norm was not |det M|; that case takes
        // the full-HNF route.
        flint::Fmpz modulus;
        fmpz_abs(modulus.raw(), fmpq_numref(norm));
        flint::FmpzMat hnf(n, n);
        ::fmpz_mat_hnf_modular(hnf.raw(), rows.raw(), modulus.raw());
        if (hnf_diagonal_product_equals(flint::FmpzMatConstRef(hnf),
                                        modulus.raw())) {
            assert(equals_full_hnf(flint::FmpzMatConstRef(rows),
                                   flint::FmpzMatConstRef(hnf)));
            // The caller vouches for closure, as with set_known_ideal_rows.
            if (!store_hnf(candidate, flint::FmpzMatConstRef(hnf),
                           HnfInput::normalized, IdealTest::assert_only)) {
                return false;
            }
            stored = true;
        }
    }
    if ((!stored &&
         !candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) ||
        !candidate.principal_generator_.define(ideal.parent_) ||
        !candidate.principal_generator_.set(generator)) {
        return false;
    }
    candidate.has_principal_generator_ = true;
    ideal.swap(candidate);
    return true;
}

bool IdealAccess::set_known_ideal_hnf(Ideal& ideal,
                                      flint::FmpzMatConstRef hnf) noexcept {
    // As Ideal::set_hnf, minus the release-build hnf_is_ideal test: the caller
    // vouches that the row span is an ideal.
    return store_hnf(ideal, hnf, HnfInput::unnormalized,
                     IdealTest::assert_only);
}

namespace {

class IdealSearchRandom {
public:
    IdealSearchRandom() noexcept {
        flint_rand_init(state_);
        flint_rand_set_seed(state_, 0x53494c45UL, 0x58UL);
    }
    ~IdealSearchRandom() noexcept { flint_rand_clear(state_); }
    IdealSearchRandom(const IdealSearchRandom&) = delete;
    IdealSearchRandom& operator=(const IdealSearchRandom&) = delete;
    void coefficient(fmpz_t out, const fmpz_t modulus) noexcept {
        fmpz_randm(out, state_, modulus);
    }
private:
    flint_rand_t state_;
};

bool scalar_hnf(flint::FmpzMatConstRef h) noexcept {
    const slong n = fmpz_mat_nrows(h.raw());
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (i == j) {
                if (!fmpz_equal(fmpz_mat_entry(h.raw(), i, j),
                                fmpz_mat_entry(h.raw(), 0, 0))) return false;
            } else if (!fmpz_is_zero(fmpz_mat_entry(h.raw(), i, j))) {
                return false;
            }
        }
    }
    return true;
}

void primitive_hnf(flint::FmpzMat& h, flint::Fmpz& content) noexcept {
    fmpz_zero(content.raw());
    for (slong i = 0; i < fmpz_mat_nrows(h.raw()); ++i) {
        for (slong j = 0; j < fmpz_mat_ncols(h.raw()); ++j) {
            fmpz_gcd(content.raw(), content.raw(), fmpz_mat_entry(h.raw(), i, j));
        }
    }
    fmpz_mat_scalar_divexact_fmpz(h.raw(), h.raw(), content.raw());
}

}  // namespace

bool ideal_minimum_generic(flint::Fmpz& minimum,
                           flint::FmpzMatConstRef coordinates,
                           flint::FmpzMatConstRef h) noexcept {
    // Source: assure_has_minimum solves in the ideal basis when 1 is not the
    // first order-basis element. The rows of H are the ideal basis, so solve
    // H^t x = [1]_O^t; this is ideal_minimum's route for nontriangular H.
    const slong n = fmpz_mat_ncols(h.raw());
    flint::FmpzMat transposed(n, n), rhs(n, 1), solution(n, 1);
    flint::Fmpz denominator, common;
    fmpz_mat_transpose(transposed.raw(), h.raw());
    fmpz_mat_transpose(rhs.raw(), coordinates.raw());
    if (!fmpz_mat_solve(solution.raw(), denominator.raw(), transposed.raw(), rhs.raw())) {
        return false;
    }
    fmpz_set(common.raw(), denominator.raw());
    for (slong i = 0; i < n; ++i) {
        fmpz_gcd(common.raw(), common.raw(), fmpz_mat_entry(solution.raw(), i, 0));
    }
    fmpz_divexact(minimum.raw(), denominator.raw(), common.raw());
    fmpz_abs(minimum.raw(), minimum.raw());
    return fmpz_sgn(minimum.raw()) > 0;
}

bool ideal_minimum(flint::Fmpz& minimum, flint::FmpzMatConstRef coordinates,
                   flint::FmpzMatConstRef h) noexcept {
    // Source: assure_has_minimum solves in the ideal basis when 1 is not the
    // first order-basis element. The minimum is the least m > 0 with m x
    // integral, where x H = [1]_O, that is H^t x^t = [1]_O^t.
    //
    // Every caller passes a stored ideal HNF (fmpz_mat_hnf or its modular
    // variants, or that matrix divided by its content): FLINT's row HNF is
    // upper triangular with a positive diagonal, so H^t is lower triangular
    // and x follows by forward substitution, column j giving
    //   x_j = (c_j - sum_{i<j} H[i][j] x_i) / H[j][j].
    // Any other square input keeps the generic rational solve.
    //
    // Exactness: x is kept as y / d with integer y and d > 0. With
    // s = c_j d - sum_{i<j} H[i][j] y_i, x_j = s / (d H[j][j]). Put
    // g = gcd(s, H[j][j]) > 0, a = H[j][j] / g and b = s / g; both
    // divisions are exact because g divides each dividend. After a sign
    // change making a > 0, x_j = b / (d a), so d becomes d a, the earlier
    // y_i are multiplied by a and y_j = b. Then m x is integral exactly when
    // d / gcd(d, y) divides m, the value the generic solve returns as
    // den / gcd(den, solution): x is unique because H is nonsingular.
    const slong n = fmpz_mat_ncols(h.raw());
    if (fmpz_mat_nrows(h.raw()) != n || fmpz_mat_nrows(coordinates.raw()) != 1 ||
        fmpz_mat_ncols(coordinates.raw()) != n) return false;
    for (slong j = 0; j < n; ++j) {
        bool triangular = !fmpz_is_zero(fmpz_mat_entry(h.raw(), j, j));
        for (slong i = j + 1; triangular && i < n; ++i) {
            triangular = fmpz_is_zero(fmpz_mat_entry(h.raw(), i, j)) != 0;
        }
        if (!triangular) return ideal_minimum_generic(minimum, coordinates, h);
    }
    flint::FmpzMat numerators(1, n);
    flint::Fmpz denominator, sum, divisor, scale;
    fmpz_one(denominator.raw());
    for (slong j = 0; j < n; ++j) {
        const fmpz* diagonal = fmpz_mat_entry(h.raw(), j, j);
        fmpz_mul(sum.raw(), fmpz_mat_entry(coordinates.raw(), 0, j), denominator.raw());
        for (slong i = 0; i < j; ++i) {
            fmpz_submul(sum.raw(), fmpz_mat_entry(h.raw(), i, j),
                        fmpz_mat_entry(numerators.raw(), 0, i));
        }
        fmpz_gcd(divisor.raw(), sum.raw(), diagonal);
        fmpz_divexact(scale.raw(), diagonal, divisor.raw());
        fmpz* numerator = fmpz_mat_entry(numerators.raw(), 0, j);
        fmpz_divexact(numerator, sum.raw(), divisor.raw());
        if (fmpz_sgn(scale.raw()) < 0) {
            fmpz_neg(scale.raw(), scale.raw());
            fmpz_neg(numerator, numerator);
        }
        if (fmpz_is_one(scale.raw())) continue;
        fmpz_mul(denominator.raw(), denominator.raw(), scale.raw());
        for (slong i = 0; i < j; ++i) {
            fmpz_mul(fmpz_mat_entry(numerators.raw(), 0, i),
                     fmpz_mat_entry(numerators.raw(), 0, i), scale.raw());
        }
    }
    fmpz_set(divisor.raw(), denominator.raw());
    for (slong i = 0; i < n && !fmpz_is_one(divisor.raw()); ++i) {
        fmpz_gcd(divisor.raw(), divisor.raw(), fmpz_mat_entry(numerators.raw(), 0, i));
    }
    fmpz_divexact(minimum.raw(), denominator.raw(), divisor.raw());
    return true;
}

void zero_skipping_mul(flint::FmpzMatRef out, flint::FmpzMatConstRef left,
                       flint::FmpzMatConstRef right) noexcept {
    // Not a source algorithm: the classical row-by-row product, omitting the
    // terms whose left factor is zero.  Each output row is accumulated as
    // sum_j left[i, j] right[j, *] over the nonzero left[i, j] only, so the
    // result is the integer product fmpz_mat_mul returns.
    const slong rows = fmpz_mat_nrows(left.raw());
    const slong inner = fmpz_mat_ncols(left.raw());
    const slong columns = fmpz_mat_ncols(right.raw());
    assert(fmpz_mat_nrows(right.raw()) == inner);
    assert(fmpz_mat_nrows(out.raw()) == rows && fmpz_mat_ncols(out.raw()) == columns);
    for (slong i = 0; i < rows; ++i) {
        fmpz* target = fmpz_mat_entry(out.raw(), i, 0);
        _fmpz_vec_zero(target, columns);
        for (slong j = 0; j < inner; ++j) {
            const fmpz* coefficient = fmpz_mat_entry(left.raw(), i, j);
            if (fmpz_is_zero(coefficient)) continue;
            const fmpz* source = fmpz_mat_entry(right.raw(), j, 0);
            if (fmpz_is_one(coefficient)) {
                _fmpz_vec_add(target, target, source, columns);
            } else {
                _fmpz_vec_scalar_addmul_fmpz(target, source, columns, coefficient);
            }
        }
    }
}

namespace {

void two_generator_rows(flint::FmpzMat& rows, flint::FmpzMatConstRef h,
                        flint::FmpzConstRef scalar,
                        flint::FmpzMatConstRef multiplication) noexcept {
    const slong n = fmpz_mat_ncols(h.raw());
    flint::FmpzMatWindow scalar_rows(rows, 0, 0, n, n);
    flint::FmpzMatWindow element_rows(rows, n, 0, 2 * n, n);
    fmpz_mat_scalar_mul_fmpz(scalar_rows.raw(), h.raw(), scalar.raw());
    zero_skipping_mul(flint::FmpzMatRef(element_rows.raw()), h, multiplication);
}

bool extract_two_generator(flint::FmpzMat& coordinates,
                           flint::FmpzMat& multiplication,
                           Order& order, flint::FmpzMatConstRef h,
                           const flint::Fmpz& minimum,
                           const TwoGeneratorSearchOptions& options,
                           TwoGeneratorSearchReport& report) noexcept {
    // Source: base4.c:get_random_a / ok_elt. Reduce multiplication matrices
    // modulo a, and accept only HNF([a I; M_beta]) == H. Unlike the source's
    // unbounded loop, exhaustion returns to the exact basis-product algorithm.
    const slong n = order.degree();
    flint::FmpzMat row(1, n), matrix(n, n), rows(2 * n, n);
    // A view of the order's cached table; `order` is not reassigned below.
    const flint::FmpzMat* table_view =
            detail::OrderAccess::multiplication_table(order);
    if (table_view == nullptr) return false;
    const flint::FmpzMat& table = *table_view;
    std::vector<flint::FmpzMat> basis_rows;
    std::vector<flint::FmpzMat> basis_matrices;
    // Source: base4.c ok_elt takes ZM_hnfmodid(M_beta, a), the HNF modulo the
    // minimum a.  The rows a I and M_beta span a O + beta O, which contains
    // a Z^n, so every elementary divisor divides a and the rank is n: the
    // contract of FLINT fmpz_mat_hnf_modular_eldiv with D = a, whose result
    // is the unique HNF that fmpz_mat_hnf would return.
    const auto verifies = [&]() noexcept {
        // Shortcut, not a correctness requirement: for a = 1, a O + beta O = O,
        // so compare H with I directly.
        if (fmpz_is_one(minimum.raw())) return fmpz_mat_is_one(h.raw()) != 0;
        fmpz_mat_zero(rows.raw());
        for (slong i = 0; i < n; ++i) {
            fmpz_set(fmpz_mat_entry(rows.raw(), i, i), minimum.raw());
            copy_row(flint::FmpzMatRef(rows), n + i,
                     flint::FmpzMatConstRef(matrix), i);
        }
        fmpz_mat_hnf_modular_eldiv(rows.raw(), minimum.raw());
        // The window is taken after the HNF, which may permute row storage.
        flint::FmpzMatConstWindow top(rows, 0, 0, n, n);
        return fmpz_mat_equal(top.raw(), h.raw()) != 0;
    };
    const auto publish = [&]() noexcept {
        // Modular matrices certify the presentation only. Entrywise reduction
        // does not preserve the action on the other ideal in a product.
        multiplication_from_coordinates(flint::FmpzMatRef(matrix),
                flint::FmpzMatConstRef(row), flint::FmpzMatConstRef(table));
        coordinates.swap(row);
        multiplication.swap(matrix);
        report.verified = true;
        return true;
    };
    for (slong i = 0; i < n; ++i) {
        copy_row(flint::FmpzMatRef(row), 0, h, i);
        multiplication_from_coordinates(flint::FmpzMatRef(matrix),
                flint::FmpzMatConstRef(row), flint::FmpzMatConstRef(table));
        fmpz_mat_scalar_mod_fmpz(matrix.raw(), matrix.raw(), minimum.raw());
        if (fmpz_mat_is_zero(matrix.raw())) continue;
        if (options.basis_candidates) {
            ++report.basis_trials;
            if (verifies()) return publish();
        }
        basis_rows.emplace_back(1, n);
        basis_matrices.emplace_back(n, n);
        fmpz_mat_set(basis_rows.back().raw(), row.raw());
        fmpz_mat_set(basis_matrices.back().raw(), matrix.raw());
    }
    if (basis_rows.empty()) {
        // The scalar ideal (a) has the valid presentation (a, 0).
        fmpz_mat_zero(row.raw());
        fmpz_mat_zero(matrix.raw());
        return verifies() && publish();
    }
    IdealSearchRandom random;
    flint::Fmpz coefficient;
    for (slong trial = 0; trial < options.random_trials; ++trial) {
        ++report.random_trials;
        fmpz_mat_zero(row.raw());
        fmpz_mat_zero(matrix.raw());
        for (std::size_t k = 0; k < basis_rows.size(); ++k) {
            random.coefficient(coefficient.raw(), minimum.raw());
            for (slong j = 0; j < n; ++j) {
                fmpz_addmul(fmpz_mat_entry(row.raw(), 0, j), coefficient.raw(),
                            fmpz_mat_entry(basis_rows[k].raw(), 0, j));
                for (slong i = 0; i < n; ++i) {
                    fmpz_addmul(fmpz_mat_entry(matrix.raw(), i, j), coefficient.raw(),
                                fmpz_mat_entry(basis_matrices[k].raw(), i, j));
                }
            }
        }
        fmpz_mat_scalar_mod_fmpz(row.raw(), row.raw(), minimum.raw());
        fmpz_mat_scalar_mod_fmpz(matrix.raw(), matrix.raw(), minimum.raw());
        if (!fmpz_mat_is_zero(matrix.raw()) && verifies()) return publish();
    }
    return false;
}

}  // namespace

bool ideal_two_generator(flint::FmpzRef scalar, OrderElement& element,
                         const Ideal& ideal,
                         const TwoGeneratorSearchOptions& options,
                         TwoGeneratorSearchReport* output_report) noexcept {
    TwoGeneratorSearchReport report;
    if (output_report != nullptr) *output_report = report;
    const Order* parent = ideal.parent();
    if (scalar.raw() == nullptr || parent == nullptr || !ideal.has_hnf() ||
        !parent->is_maximal() || element.parent() == nullptr ||
        !element.parent()->has_same_data(*parent) || options.random_trials < 0 ||
        options.random_trials > 32) return false;
    Order order(*parent);
    const slong n = order.degree();
    flint::FmpzMat h(n, n), row(1, n), matrix(n, n);
    flint::Fmpz minimum;
    OrderElement one(order);
    if (!ideal.get_hnf(flint::FmpzMatRef(h)) ||
        !one.one() || !one.get_coordinates(flint::FmpzMatRef(row)) ||
        !ideal_minimum(minimum, flint::FmpzMatConstRef(row), flint::FmpzMatConstRef(h))) return false;
    const bool found = extract_two_generator(row, matrix, order,
            flint::FmpzMatConstRef(h), minimum, options, report);
    if (output_report != nullptr) *output_report = report;
    OrderElement candidate(order);
    if (!found || !candidate.set_coordinates(flint::FmpzMatConstRef(row))) return false;
    fmpz_set(scalar.raw(), minimum.raw());
    element.swap(candidate);
    return true;
}

bool try_ideal_product_hnf(flint::FmpzMatRef out, const Ideal& left,
                           const Ideal& right,
                           const TwoGeneratorSearchOptions& options,
                           TwoGeneratorSearchReport* output_report) noexcept {
    TwoGeneratorSearchReport report;
    if (output_report != nullptr) *output_report = report;
    const Order* parent = left.parent();
    if (parent == nullptr || right.parent() == nullptr ||
        !parent->has_same_data(*right.parent()) || !left.has_hnf() || !right.has_hnf() ||
        !fmpz_mat_shape(out, parent->degree(), parent->degree()) ||
        options.random_trials < 0 || options.random_trials > 32) return false;
    const slong n = parent->degree();
    flint::FmpzMat h(n, n), k(n, n), result(n, n);
    if (!left.get_hnf(flint::FmpzMatRef(h)) || !right.get_hnf(flint::FmpzMatRef(k))) return false;
    if (scalar_hnf(flint::FmpzMatConstRef(h))) {
        fmpz_mat_scalar_mul_fmpz(result.raw(), k.raw(), fmpz_mat_entry(h.raw(), 0, 0));
    } else if (scalar_hnf(flint::FmpzMatConstRef(k))) {
        fmpz_mat_scalar_mul_fmpz(result.raw(), h.raw(), fmpz_mat_entry(k.raw(), 0, 0));
    } else if (left.known_principal_generator() != nullptr ||
               right.known_principal_generator() != nullptr) {
        const bool use_left = left.known_principal_generator() != nullptr;
        const OrderElement* generator = use_left ? left.known_principal_generator()
                                                : right.known_principal_generator();
        Element multiplier(*parent->parent());
        Ideal product(*parent);
        if (!generator->get_element(multiplier) ||
            !multiply_integral_ideal_by_element(product, use_left ? right : left, multiplier) ||
            !product.get_hnf(flint::FmpzMatRef(result))) return false;
    } else {
        if (!parent->is_maximal()) return false;
        Order order(*parent);
        flint::Fmpz ch, ck, mh, mk;
        primitive_hnf(h, ch);
        primitive_hnf(k, ck);
        flint::FmpzMat one_coordinates(1, n);
        OrderElement one(order);
        if (!one.one() || !one.get_coordinates(flint::FmpzMatRef(one_coordinates)) ||
            !ideal_minimum(mh, flint::FmpzMatConstRef(one_coordinates), flint::FmpzMatConstRef(h)) ||
            !ideal_minimum(mk, flint::FmpzMatConstRef(one_coordinates), flint::FmpzMatConstRef(k))) return false;
        if (fmpz_cmp(mh.raw(), mk.raw()) < 0) {
            h.swap(k);
            mh.swap(mk);
        }
        // Extract from k (the smaller minimum; right wins ties).
        flint::FmpzMat row(1, n), matrix(n, n), rows(2 * n, n);
        const bool found = extract_two_generator(row, matrix, order,
                flint::FmpzMatConstRef(k), mk, options, report);
        if (output_report != nullptr) *output_report = report;
        if (!found) return false;
        two_generator_rows(rows, flint::FmpzMatConstRef(h),
                           flint::FmpzConstRef(mk), flint::FmpzMatConstRef(matrix));
        // Source trace: reference `base4.c:idealHNF_mul_two` returns
        // ZM_hnfmodid(m, a * x[1,1]) for the rows m = [x a; M_alpha x], with
        // a the scalar generator and x[1,1] the minimum of x.  Here the rows
        // mk H and H M_beta span mk h + beta h = h (mk O + beta O) = h k,
        // because extract_two_generator verified k = mk O + beta O.  Since
        // mh is in h and mk is in k, mh mk lies in the ideal h k, so
        // mh mk Z^n is inside the row lattice and every elementary divisor
        // divides mh mk.  That is the contract of FLINT
        // fmpz_mat_hnf_modular_eldiv with D = mh mk, whose result is the
        // unique HNF that fmpz_mat_hnf returns.  The argument uses only that
        // h k is an ideal of the order, not that the order is maximal.
        // mh mk > 1: a minimum of one would put 1 in the primitive h (or k),
        // making it O and the original ideal scalar, which the scalar
        // branches above have already handled.
        flint::Fmpz exponent;
        fmpz_mul(exponent.raw(), mh.raw(), mk.raw());
        fmpz_mat_hnf_modular_eldiv(rows.raw(), exponent.raw());
        // The window is taken after the HNF, which may permute row storage.
        flint::FmpzMatConstWindow top(rows, 0, 0, n, n);
        fmpz_mul(ch.raw(), ch.raw(), ck.raw());
        fmpz_mat_scalar_mul_fmpz(result.raw(), top.raw(), ch.raw());
    }
    fmpz_mat_set(out.raw(), result.raw());
    return true;
}

bool ideal_two_generator_form(IdealTwoGeneratorForm& out, const Ideal& ideal,
                              const TwoGeneratorSearchOptions& options,
                              TwoGeneratorSearchReport* output_report) noexcept {
    TwoGeneratorSearchReport report;
    if (output_report != nullptr) *output_report = report;
    const Order* parent = ideal.parent();
    if (parent == nullptr || !ideal.has_hnf() || !parent->is_maximal() ||
        options.random_trials < 0 || options.random_trials > 32) return false;
    Order order(*parent);
    const slong n = order.degree();
    flint::FmpzMat h(n, n), k(n, n), one_coordinates(1, n), row(1, n), matrix(n, n);
    flint::Fmpz content, minimum;
    OrderElement one(order);
    if (!ideal.get_hnf(flint::FmpzMatRef(h)) || scalar_hnf(flint::FmpzMatConstRef(h)) ||
        !one.one() || !one.get_coordinates(flint::FmpzMatRef(one_coordinates))) return false;
    // The extraction try_ideal_product_hnf runs on its factor k: the primitive
    // part of the HNF and its minimum.
    fmpz_mat_set(k.raw(), h.raw());
    primitive_hnf(k, content);
    if (!ideal_minimum(minimum, flint::FmpzMatConstRef(one_coordinates),
                       flint::FmpzMatConstRef(k))) return false;
    const bool found = extract_two_generator(row, matrix, order,
            flint::FmpzMatConstRef(k), minimum, options, report);
    if (output_report != nullptr) *output_report = report;
    if (!found) return false;
    out.hnf.swap(h);
    out.content.swap(content);
    out.minimum.swap(minimum);
    out.multiplication.swap(matrix);
    return true;
}

bool try_ideal_product_hnf_with_form(flint::FmpzMatRef out, const Ideal& left,
                                     const Ideal& right,
                                     const IdealTwoGeneratorForm& right_form) noexcept {
    const Order* parent = left.parent();
    if (parent == nullptr || right.parent() == nullptr ||
        !parent->has_same_data(*right.parent()) || !left.has_hnf() || !right.has_hnf() ||
        !fmpz_mat_shape(out, parent->degree(), parent->degree())) return false;
    const slong n = parent->degree();
    flint::FmpzMat h(n, n), k(n, n);
    if (!left.get_hnf(flint::FmpzMatRef(h)) || !right.get_hnf(flint::FmpzMatRef(k)) ||
        !fmpz_mat_shape(flint::FmpzMatConstRef(right_form.hnf), n, n) ||
        !fmpz_mat_shape(flint::FmpzMatConstRef(right_form.multiplication), n, n) ||
        !fmpz_mat_equal(k.raw(), right_form.hnf.raw())) return false;
    // The branches of try_ideal_product_hnf that do not search for a
    // two-element form are taken as they are.
    if (scalar_hnf(flint::FmpzMatConstRef(h)) || scalar_hnf(flint::FmpzMatConstRef(k)) ||
        left.known_principal_generator() != nullptr ||
        right.known_principal_generator() != nullptr) {
        return try_ideal_product_hnf(out, left, right);
    }
    if (!parent->is_maximal()) return false;
    flint::Fmpz ch, mh;
    primitive_hnf(h, ch);
    flint::FmpzMat one_coordinates(1, n);
    OrderElement one(*parent);
    if (!one.one() || !one.get_coordinates(flint::FmpzMatRef(one_coordinates)) ||
        !ideal_minimum(mh, flint::FmpzMatConstRef(one_coordinates), flint::FmpzMatConstRef(h))) return false;
    // The rows mk H and H M_beta span h (mk O + beta O) = h k for the
    // primitive parts h of left and k of right, because
    // ideal_two_generator_form verified k = mk O + beta O, and mh mk lies in
    // h k: the argument of the generic branch of try_ideal_product_hnf, with
    // the form of k found once instead of on every product.  The HNF of an
    // ideal is unique, so this is the HNF that try_ideal_product_hnf returns,
    // whichever factor it would have searched.
    flint::FmpzMat rows(2 * n, n);
    two_generator_rows(rows, flint::FmpzMatConstRef(h), flint::FmpzConstRef(right_form.minimum),
                       flint::FmpzMatConstRef(right_form.multiplication));
    flint::Fmpz exponent;
    fmpz_mul(exponent.raw(), mh.raw(), right_form.minimum.raw());
    fmpz_mat_hnf_modular_eldiv(rows.raw(), exponent.raw());
    // The window is taken after the HNF, which may permute row storage.
    flint::FmpzMatConstWindow top(rows, 0, 0, n, n);
    fmpz_mul(ch.raw(), ch.raw(), right_form.content.raw());
    fmpz_mat_scalar_mul_fmpz(out.raw(), top.raw(), ch.raw());
    return true;
}

bool basis_product_hnf(flint::FmpzMatRef out, const Ideal& left,
                       const Ideal& right) noexcept {
    const Order* parent = left.parent();
    if (parent == nullptr || right.parent() == nullptr ||
        !parent->has_same_data(*right.parent()) || !left.has_hnf() || !right.has_hnf() ||
        !fmpz_mat_shape(out, parent->degree(), parent->degree())) return false;
    const slong n = parent->degree();
    Order order(*parent);
    // A view of the order's cached table; `order` is not reassigned below.
    const flint::FmpzMat* table = detail::OrderAccess::multiplication_table(order);
    flint::FmpzMat h(n, n), k(n, n), one_coordinates(1, n), row(1, n), matrix(n, n),
            rows(n * n, n), result(n, n);
    flint::Fmpz ch, ck, mh, mk, modulus;
    OrderElement one(order);
    if (table == nullptr || !left.get_hnf(flint::FmpzMatRef(h)) ||
        !right.get_hnf(flint::FmpzMatRef(k))) return false;
    primitive_hnf(h, ch);
    primitive_hnf(k, ck);
    if (!one.one() || !one.get_coordinates(flint::FmpzMatRef(one_coordinates)) ||
        !ideal_minimum(mh, flint::FmpzMatConstRef(one_coordinates), flint::FmpzMatConstRef(h)) ||
        !ideal_minimum(mk, flint::FmpzMatConstRef(one_coordinates), flint::FmpzMatConstRef(k))) return false;
    // Rows j n .. j n + n - 1 are H M_{k_j}: row i is the coordinate vector of
    // h_i k_j, since row b of M_{k_j} (multiplication_from_coordinates) is
    // e_b k_j.  Together they are the n^2 basis products that span h k.
    for (slong j = 0; j < n; ++j) {
        copy_row(flint::FmpzMatRef(row), 0, flint::FmpzMatConstRef(k), j);
        multiplication_from_coordinates(flint::FmpzMatRef(matrix),
                flint::FmpzMatConstRef(row), flint::FmpzMatConstRef(*table));
        flint::FmpzMatWindow block(rows, j * n, 0, j * n + n, n);
        zero_skipping_mul(flint::FmpzMatRef(block.raw()), flint::FmpzMatConstRef(h),
                          flint::FmpzMatConstRef(matrix));
    }
    // Source trace: reference `NumFieldOrd/NfOrd/Ideal/Arithmetic.jl:mul_gen`
    // takes hnf_modular_eldiv of the basis products modulo
    // l = minimum(x) minimum(y), on any order.  The modulus is valid because
    // mh is in h and mk is in k, so mh mk lies in the ideal h k; then
    // mh mk Z^n is inside the row lattice, which has rank n, and its largest
    // elementary divisor divides mh mk.  That is the contract of FLINT
    // fmpz_mat_hnf_modular_eldiv with D = mh mk (an m x n matrix with
    // m >= n), whose result is the unique HNF that fmpz_mat_hnf (and the
    // exact Lat::hnf route of multiply_generic) returns.  The argument uses
    // only that h and k are ideals of the order, not maximality.
    fmpz_mul(modulus.raw(), mh.raw(), mk.raw());
    if (fmpz_is_one(modulus.raw())) {
        // Both primitive parts contain 1, so both are O and so is h k.
        fmpz_mat_one(result.raw());
    } else {
        fmpz_mat_hnf_modular_eldiv(rows.raw(), modulus.raw());
        // The window is taken after the HNF, which may permute row storage.
        flint::FmpzMatConstWindow top(rows, 0, 0, n, n);
        fmpz_mat_set(result.raw(), top.raw());
    }
    fmpz_mul(ch.raw(), ch.raw(), ck.raw());
    fmpz_mat_scalar_mul_fmpz(out.raw(), result.raw(), ch.raw());
    return true;
}

bool IdealAccess::set_known_two_generator(
        Ideal& out,
        flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator) noexcept {
    const Order* element_parent = element_generator.parent();
    if (!out.is_defined() || element_parent == nullptr ||
        !element_parent->has_same_data(out.parent_) ||
        flint::fmpz_sgn(scalar_generator) <= 0) {
        return false;
    }

    // Source trace: reference `Ideal.jl` (lower-left convention, lines 422 and
    // 444-455) builds a two-element ideal basis as
    // `hnf_modular_eldiv!(representation_matrix(gen_two), gen_one, :lowerleft)`,
    // which is FLINT `fmpz_mat_hnf_modular_eldiv`, and reference
    // `base4.c:idealhnf_two` computes the HNF of the same lattice (in the
    // reference's column convention) as `ZM_hnfmodid` of the multiplication
    // table of beta modulo the scalar generator. Silex keeps FLINT's
    // upper-right HNF convention.
    //
    // The rows a I and M_beta span a O + beta O, which contains a Z^n, so
    // every elementary divisor of the 2n x n row matrix divides a and its rank
    // is n.  That is the documented contract of `fmpz_mat_hnf_modular_eldiv`
    // with D = a (FLINT `fmpz_mat/hnf_modular_eldiv.c`, Fieker-Hofmann 2014).
    // The HNF of a full-rank lattice is unique, so the result is the same
    // matrix as the generic `fmpz_mat_hnf` of these rows.
    const slong n = out.parent_.degree();
    flint::FmpzMat element_multiplication(n, n);
    flint::FmpzMat rows(2 * n, n);
    if (!order_element_multiplication_matrix(
                flint::FmpzMatRef(element_multiplication),
                out.parent_,
                element_generator)) {
        return false;
    }

    fmpz_mat_zero(rows.raw());
    for (slong i = 0; i < n; ++i) {
        fmpz_set(fmpz_mat_entry(rows.raw(), i, i),
                 scalar_generator.raw());
    }

    for (slong row = 0; row < n; ++row) {
        copy_row(flint::FmpzMatRef(rows), n + row,
                 flint::FmpzMatConstRef(element_multiplication), row);
    }

    Ideal candidate(out.parent_);
    if (flint::fmpz_is_one(scalar_generator)) {
        // a O + beta O = O.  FLINT's word path would reduce modulo 1.
        fmpz_mat_zero(rows.raw());
        for (slong i = 0; i < n; ++i) {
            fmpz_one(fmpz_mat_entry(rows.raw(), i, i));
        }
    } else {
        fmpz_mat_hnf_modular_eldiv(rows.raw(), scalar_generator.raw());
    }
    // The window is taken after the HNF, which may permute row storage.
    flint::FmpzMatConstWindow top(rows, 0, 0, n, n);
    if (!candidate.set_hnf_direct(top.const_ref())) {
        return false;
    }

    out.swap(candidate);
    return true;
}

namespace {

// Debug cross-check of the closed form below: `hnf` equals the HNF that
// set_known_two_generator stores for a O + (theta - r) O.
[[maybe_unused]] bool linear_hnf_matches_modular(const Order& order,
                                                 flint::FmpzConstRef scalar,
                                                 flint::FmpzConstRef root,
                                                 flint::FmpzMatConstRef hnf) noexcept {
    // theta - r in the basis 1, theta, ..., theta^(n-1); for n = 1, theta is
    // the integer -F_0 of the monic F = x + F_0.
    const slong n = order.degree();
    flint::FmpzMat coordinates(1, n), expected(n, n);
    fmpz* constant = fmpz_mat_entry(coordinates.raw(), 0, 0);
    if (n >= 2) {
        fmpz_one(fmpz_mat_entry(coordinates.raw(), 0, 1));
    } else {
        fmpz_neg(constant, fmpq_poly_numref(order.parent()->raw_flint_field()->pol));
    }
    fmpz_sub(constant, constant, root.raw());
    OrderElement beta(order);
    Ideal reference(order);
    return beta.set_coordinates(flint::FmpzMatConstRef(coordinates)) &&
           IdealAccess::set_known_two_generator(reference, scalar, beta) &&
           reference.get_hnf(flint::FmpzMatRef(expected)) &&
           fmpz_mat_equal(expected.raw(), hnf.raw()) != 0;
}

}  // namespace

bool IdealAccess::set_known_linear_two_generator(
        Ideal& out,
        flint::FmpzConstRef scalar_generator,
        flint::FmpzConstRef root) noexcept {
    if (!out.is_defined() || !out.parent_.is_equation_order() ||
        out.parent_.parent() == nullptr ||
        flint::fmpz_sgn(scalar_generator) <= 0) {
        return false;
    }
    const fmpq_poly_struct* polynomial =
            out.parent_.parent()->raw_flint_field()->pol;
    const slong n = out.parent_.degree();
    // FLINT stores an fmpq_poly as a numerator vector over one common
    // denominator, so F is monic integral iff both are one at degree n.
    if (n < 1 || fmpq_poly_degree(polynomial) != n ||
        !fmpz_is_one(fmpq_poly_denref(polynomial)) ||
        !fmpz_is_one(fmpq_poly_numref(polynomial) + n)) {
        return false;
    }

    // Source trace: reference v0.39.19
    // `src/NumFieldOrd/NfOrd/Ideal/Ideal.jl:basis_mat_prime_deg_1` (lines
    // 459-477) writes the basis matrix of a degree-one prime P directly: in
    // its lower-left convention row i is e_i - (w_i mod P) e_1, with the
    // minimum p in the corner, because O / P = F_p.  This routine is the
    // analogous construction in Silex's upper-right FLINT convention, for the
    // basis w_i = theta^i of the equation order Z[theta] = Z[x] / (F): the
    // pivot a sits in the last column, the other entries are -r^j / r^(n-1)
    // mod a, and r = 0 (theta^(n-1) in the ideal) needs its own matrix.
    //
    // When F(r) = 0 mod a, f(theta) -> f(r) mod a is a well-defined surjective
    // ring map Z[theta] -> Z / a.  Its kernel is a O + (theta - r) O, since
    // f(theta) = f(r) + (theta - r) q(theta) for f(x) = f(r) + (x - r) q(x).
    // In coordinates the ideal is the index-a lattice of c in Z^n with
    // sum_i c_i r^i = 0 mod a.  Its upper triangular HNF (pivots positive,
    // entries above a pivot reduced into [0, pivot)) is
    //   - r = 0 mod a: diag(a, 1, ..., 1) (the rows a e_0, e_1, ..., e_(n-1));
    //   - r a unit mod a, s = r^-1 mod a: the rows e_j + c_j e_(n-1) for
    //     j < n - 1 with c_j = -r^j s^(n-1) = -s^(n-1-j) mod a, and a e_(n-1).
    // Each row lies in the lattice and the rows have determinant a, its index,
    // so they are a basis, and the matrix is in HNF, hence it is the unique
    // HNF that set_known_two_generator computes for a and theta - r.  Any
    // other r (a nonzero zero divisor modulo a composite a) returns false.
    flint::Fmpz r, value, unit, power;
    fmpz_mod(r.raw(), root.raw(), scalar_generator.raw());
    const fmpz* coefficients = fmpq_poly_numref(polynomial);
    for (slong i = n; i >= 0; --i) {
        fmpz_mul(value.raw(), value.raw(), r.raw());
        fmpz_add(value.raw(), value.raw(), coefficients + i);
        fmpz_mod(value.raw(), value.raw(), scalar_generator.raw());
    }
    if (!fmpz_is_zero(value.raw())) {
        return false;
    }

    flint::FmpzMat hnf(n, n);
    for (slong i = 0; i < n; ++i) {
        fmpz_one(fmpz_mat_entry(hnf.raw(), i, i));
    }
    if (fmpz_is_zero(r.raw())) {
        fmpz_set(fmpz_mat_entry(hnf.raw(), 0, 0), scalar_generator.raw());
    } else {
        if (fmpz_invmod(unit.raw(), r.raw(), scalar_generator.raw()) == 0) {
            return false;
        }
        fmpz_one(power.raw());
        for (slong j = n - 2; j >= 0; --j) {
            fmpz_mul(power.raw(), power.raw(), unit.raw());
            fmpz_mod(power.raw(), power.raw(), scalar_generator.raw());
            fmpz_sub(fmpz_mat_entry(hnf.raw(), j, n - 1),
                     scalar_generator.raw(), power.raw());
        }
        fmpz_set(fmpz_mat_entry(hnf.raw(), n - 1, n - 1),
                 scalar_generator.raw());
    }
    // Builds with assertions compare the matrix with the general route,
    // which also shows that its row span is the ideal.
    assert(linear_hnf_matches_modular(out.parent_, scalar_generator,
                                      flint::FmpzConstRef(r),
                                      flint::FmpzMatConstRef(hnf)));
    Ideal candidate(out.parent_);
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(hnf))) {
        return false;
    }
    out.swap(candidate);
    return true;
}

namespace {

// Debug cross-check of the closed form below: `hnf` equals the HNF that
// set_known_two_generator stores for a O + g(theta) O, where `reduced` holds
// the f + 1 coefficients of the monic g reduced modulo a.
[[maybe_unused]] bool factor_hnf_matches_modular(const Order& order,
                                                 flint::FmpzConstRef scalar,
                                                 const fmpz* reduced,
                                                 slong factor_degree,
                                                 flint::FmpzMatConstRef hnf) noexcept {
    // g(theta) in the basis 1, theta, ..., theta^(n-1); for f = n, g(theta)
    // = (g - F)(theta) with deg(g - F) < n.
    const slong n = order.degree();
    const fmpz* polynomial = fmpq_poly_numref(order.parent()->raw_flint_field()->pol);
    flint::FmpzMat coordinates(1, n), expected(n, n);
    for (slong i = 0; i < n && i <= factor_degree; ++i) {
        fmpz_set(fmpz_mat_entry(coordinates.raw(), 0, i), reduced + i);
        if (factor_degree == n) {
            fmpz_sub(fmpz_mat_entry(coordinates.raw(), 0, i),
                     fmpz_mat_entry(coordinates.raw(), 0, i), polynomial + i);
        }
    }
    OrderElement beta(order);
    Ideal reference(order);
    return beta.set_coordinates(flint::FmpzMatConstRef(coordinates)) &&
           IdealAccess::set_known_two_generator(reference, scalar, beta) &&
           reference.get_hnf(flint::FmpzMatRef(expected)) &&
           fmpz_mat_equal(expected.raw(), hnf.raw()) != 0;
}

}  // namespace

bool IdealAccess::set_known_factor_two_generator(
        Ideal& out,
        flint::FmpzConstRef scalar_generator,
        flint::FmpzPolyConstRef factor) noexcept {
    if (!out.is_defined() || !out.parent_.is_equation_order() ||
        out.parent_.parent() == nullptr ||
        flint::fmpz_sgn(scalar_generator) <= 0) {
        return false;
    }
    const fmpq_poly_struct* polynomial =
            out.parent_.parent()->raw_flint_field()->pol;
    const slong n = out.parent_.degree();
    const slong f = fmpz_poly_degree(factor.raw());
    // F must be monic integral (see set_known_linear_two_generator) and g
    // monic of degree 1 <= f <= n.
    if (n < 1 || fmpq_poly_degree(polynomial) != n ||
        !fmpz_is_one(fmpq_poly_denref(polynomial)) ||
        !fmpz_is_one(fmpq_poly_numref(polynomial) + n) || f < 1 || f > n ||
        !fmpz_is_one(factor.raw()->coeffs + f)) {
        return false;
    }

    // Source trace: the degree-f analogue of set_known_linear_two_generator
    // above (reference v0.39.19 `Ideal.jl:basis_mat_prime_deg_1`, lines
    // 459-477, writes the basis of a degree-one prime from O / P = F_p).
    // That routine covers degree one only; no upstream routine computes this
    // HNF for f > 1, so the closed form below is derived here.  Dedekind-
    // Kummer primes P = (p, g(theta)) of residue degree f arise in the
    // setting of Cohen, GTM 138, Theorem 4.8.13, but the formula does not
    // use that theorem or its hypothesis that p does not divide the index:
    // the identification of Z[theta] / (a, g(theta)) with (Z / a)[x] / (g)
    // is elementary and is shown below for every a > 0.
    //
    // When g | F modulo a, f(theta) -> f mod (a, g) is a well-defined
    // surjective ring map Z[theta] -> R = (Z / a)[x] / (g), and its kernel is
    // a O + g(theta) O (f lies in the kernel iff f is in (a, g) Z[x]).  As g
    // is monic, R is free over Z / a on 1, x, ..., x^(f-1), so the kernel has
    // index a^f.  When g(0) is a unit modulo a, x is a unit of R (x h = 1 for
    // h = -(g - g(0)) / (x g(0))), so x^(n-f), ..., x^(n-1) are also a basis
    // of R, with u = x^-1: x^j = x^(n-f) u^(n-f-j) in R.  In Silex's upper
    // triangular convention (pivots positive, entries above a pivot reduced
    // into [0, pivot)) the HNF of the kernel is therefore
    //   - rows j < n - f: e_j - sum_k w_(j,k) e_(n-f+k) with entries reduced
    //     modulo a, where u^(n-f-j) = sum_k w_(j,k) x^k in R;
    //   - rows n - f + k, k < f: a e_(n-f+k).
    // Each row lies in the kernel and the determinant is a^f, its index, so
    // the rows are a basis in HNF: the unique HNF set_known_two_generator
    // computes for a and g(theta).  For f = 1 and g = x - r this is the
    // matrix of set_known_linear_two_generator with r a unit.  Any g that is
    // not a factor of F modulo a, or has g(0) not a unit, returns false.
    const fmpz* a = scalar_generator.raw();
    const fmpz* coefficients = fmpq_poly_numref(polynomial);
    flint::FmpzVec reduced_storage(f + 1);
    flint::FmpzVec remainder_storage(n + 1);
    flint::FmpzVec power_storage(f);
    fmpz* reduced = reduced_storage.data();
    fmpz* remainder = remainder_storage.data();
    fmpz* power = power_storage.data();
    flint::Fmpz inverse, lead;
    bool usable = true;
    _fmpz_vec_scalar_mod_fmpz(reduced, factor.raw()->coeffs, f + 1, a);
    _fmpz_vec_scalar_mod_fmpz(remainder, coefficients, n + 1, a);
    // F mod g over Z / a by schoolbook division with the monic g.
    for (slong k = n; k >= f; --k) {
        fmpz_set(lead.raw(), remainder + k);
        if (!fmpz_is_zero(lead.raw())) {
            for (slong i = 0; i < f; ++i) {
                fmpz_submul(remainder + k - f + i, lead.raw(), reduced + i);
                fmpz_mod(remainder + k - f + i, remainder + k - f + i, a);
            }
            fmpz_zero(remainder + k);
        }
    }
    usable = _fmpz_vec_is_zero(remainder, f) != 0;
    if (usable) {
        if (fmpz_is_one(a)) {
            fmpz_zero(inverse.raw());
        } else {
            usable = fmpz_invmod(inverse.raw(), reduced, a) != 0;
        }
    }

    flint::FmpzMat hnf(n, n);
    if (usable) {
        // power holds u^(n-f-j) in the basis 1, ..., x^(f-1); u * w is
        // (w - (w_0 / g_0) g) / x.
        fmpz_one(power);
        for (slong j = n - f - 1; j >= 0; --j) {
            fmpz_mul(lead.raw(), power, inverse.raw());
            fmpz_mod(lead.raw(), lead.raw(), a);
            for (slong k = 0; k + 1 < f; ++k) {
                fmpz_set(power + k, power + k + 1);
                fmpz_submul(power + k, lead.raw(), reduced + k + 1);
                fmpz_mod(power + k, power + k, a);
            }
            fmpz_neg(power + f - 1, lead.raw());
            fmpz_mod(power + f - 1, power + f - 1, a);
            fmpz_one(fmpz_mat_entry(hnf.raw(), j, j));
            for (slong k = 0; k < f; ++k) {
                fmpz* entry = fmpz_mat_entry(hnf.raw(), j, n - f + k);
                fmpz_neg(entry, power + k);
                fmpz_mod(entry, entry, a);
            }
        }
        for (slong k = n - f; k < n; ++k) {
            fmpz_set(fmpz_mat_entry(hnf.raw(), k, k), a);
        }
        // Builds with assertions compare the matrix with the general route,
        // which also shows that its row span is the ideal.
        assert(factor_hnf_matches_modular(out.parent_, scalar_generator,
                                          reduced, f,
                                          flint::FmpzMatConstRef(hnf)));
    }
    if (!usable) {
        return false;
    }
    Ideal candidate(out.parent_);
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(hnf))) {
        return false;
    }
    out.swap(candidate);
    return true;
}

// Trust boundary of the six silex::detail ideal shortcuts
// (set_known_two_generator_ideal, set_known_linear_two_generator_ideal and
// set_known_factor_two_generator_ideal above, this routine,
// multiply_integral_ideal_by_element below, and
// set_integral_ideal_known_hnf near FractionalIdeal::set_integral_den).  "Known" means that they skip the
// hnf_is_ideal closure test of Ideal::set_hnf, not that they store a claim
// the caller makes.  For every input that passes their checks they store
// only the HNF of a lattice that is an ideal by construction.  The first
// five also store no norm and no principal generator; the sixth does not
// make that claim (see its bullet):
//   - set_known_two_generator_ideal(a, beta): rows a I_n and M_beta span
//     a Z^n + beta O = a O + beta O, an ideal for every a > 0 and beta in O,
//     and every elementary divisor divides a, so the modular-HNF contract
//     holds whatever beta is.
//   - set_known_linear_two_generator_ideal(a, r): it checks F(r) = 0 mod a
//     itself, so its closed form is the HNF of a O + (theta - r) O, the
//     ideal set_known_two_generator_ideal(a, theta - r) stores.  It refuses
//     every input it cannot confirm, so it does not rely on the class-group
//     test named below, which does not call it.
//   - set_known_factor_two_generator_ideal(a, g): likewise, it checks that
//     the monic g divides F modulo a and that g(0) is a unit modulo a, so
//     its closed form is the HNF of a O + g(theta) O, the ideal
//     set_known_two_generator_ideal(a, g(theta)) stores, and it refuses
//     every input it cannot confirm.
//   - multiply_integral_ideal_by_two_generator(I, a, beta): rows a H and
//     H M_beta span a I + beta I = I (a O + beta O), an ideal whenever I is.
//   - multiply_integral_ideal_by_element(I, x): the rows are x times the
//     basis of I, so x I, kept only when it is integral and in HNF.
//   - set_integral_ideal_known_hnf(I): I / 1, a copy of the Ideal I, which
//     carries over whatever cached data (norm, principal generator) I
//     already holds.
// The inputs are Ideal, OrderElement and Element objects whose own
// invariants hold on every installed route (an Ideal with an HNF stores an
// ideal; the only unchecked store is the noninstalled test hook
// IdealAccess::set_hnf_direct).  So a "false" generator, one that does not
// generate the ideal the caller intends, gives a true description of the
// ideal it does generate, the same ideal Ideal::set_hnf, Ideal::add of
// principal ideals, or Ideal::multiply would produce.  None of the six
// reads or writes a ClassGroupContext, OrderUnitGroup or PrimeIdeal, or any
// certification state: they are friends of Ideal and FractionalIdeal only.
// Their results can reach labelled computations, for example
// ClassGroupContext::ideal_class_coordinates on a proven group, but only as
// ideals the installed API can also build, so no caller-supplied claim
// reaches a proven or grh label through them.  Library callers that attach
// meaning to a result (a prime power or a prime decomposition from
// PrimeIdeal data) rely on their own source-traced claims about the
// generators, not on these routines.  The class-group test
// test_detail_ideal_shortcuts_reach_labelled_dlog_as_true_ideals pins this.
bool IdealAccess::multiply_two_generator(
        Ideal& out,
        const Ideal& ideal,
        flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* element_parent = element_generator.parent();
    if (!out.is_defined() || !ideal.is_defined() || !ideal.has_hnf_ ||
        element_parent == nullptr ||
        !out.parent_.has_same_data(ideal.parent_) ||
        !element_parent->has_same_data(out.parent_) ||
        flint::fmpz_sgn(scalar_generator) <= 0) {
        return false;
    }

    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.multiply_two_generator");

    // Source trace: reference `base4.c:idealHNF_mul_two` forms the generators
    // [a * H; H * M_beta], then HNF-reduces the resulting module.
    const slong n = out.parent_.degree();
    flint::FmpzMat element_multiplication(n, n);
    flint::FmpzMat rows(2 * n, n);
    if (!order_element_multiplication_matrix(
                flint::FmpzMatRef(element_multiplication),
                out.parent_,
                element_generator)) {
        return false;
    }

    two_generator_rows(rows, flint::FmpzMatConstRef(ideal.hnf_),
                       scalar_generator, flint::FmpzMatConstRef(element_multiplication));

    Ideal candidate(out.parent_);
    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

// Ideal::multiply(left, right) on a maximal order takes the HNF from
// try_ideal_product_hnf and stores it with set_hnf_direct.  This takes the
// same HNF from try_ideal_product_hnf_with_form, which skips the two-element
// search, and stores it the same way.  The HNF of an ideal is unique, so the
// stored ideal and its cached state (no norm, no principal generator) are
// those Ideal::multiply would store.  The O factors, nonmaximal orders and
// any failure of the form route go to Ideal::multiply itself.
bool IdealAccess::multiply_two_generator_form(
        Ideal& out,
        const Ideal& left,
        const Ideal& right,
        const IdealTwoGeneratorForm& right_form) noexcept {
    if (out.is_defined() && left.is_defined() && right.is_defined() &&
        out.parent_.has_same_data(left.parent_) &&
        left.parent_.has_same_data(right.parent_) && left.has_hnf_ &&
        right.has_hnf_ && out.parent_.is_maximal() && !left.is_one() &&
        !right.is_one()) {
        const slong n = out.parent_.degree();
        flint::FmpzMat hnf(n, n);
        if (try_ideal_product_hnf_with_form(flint::FmpzMatRef(hnf), left,
                                            right, right_form)) {
            Ideal candidate(out.parent_);
            if (candidate.set_hnf_direct(flint::FmpzMatConstRef(hnf))) {
                out.swap(candidate);
                return true;
            }
        }
    }
    return out.multiply(left, right);
}

bool IdealAccess::multiply_by_element(
        Ideal& out,
        const Ideal& ideal,
        const Element& multiplier,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::fractional_ideal,
                        "fractional_ideal.principal_times_integral_direct");
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!out.is_defined() || order == nullptr || field == nullptr ||
        !out.parent_.has_same_data(*order) || !ideal.has_hnf() ||
        !multiplier.has_parent(*field) || multiplier.equal_si(0)) {
        return false;
    }

    // Source trace: reference `base4.c:idealmulelt` handles principal-times-matrix
    // products by multiplying the element through the ideal basis and reducing
    // that matrix directly.  reference's compact presentation forms
    // `simplify(b*eA)` and then requires denominator one.
    const slong n = order->degree();
    flint::FmpzMat ideal_hnf(n, n);
    flint::FmpzMat row(1, n);
    flint::FmpqMat product_rows(n, n);
    flint::FmpqMat product_coords(1, n);
    OrderElement basis_element(*order);
    Element basis_field_element(*field);
    Element product(*field);
    if (!basis_element.is_defined() || !basis_field_element.is_defined() ||
        !product.is_defined() ||
        !ideal.get_hnf(flint::FmpzMatRef(ideal_hnf))) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::fractional_ideal,
                "fractional_ideal.principal_times_integral_rows");
        for (slong i = 0; i < n; ++i) {
            copy_row(flint::FmpzMatRef(row), 0,
                     flint::FmpzMatConstRef(ideal_hnf), i);
            if (!basis_element.set_coordinates(flint::FmpzMatConstRef(row)) ||
                !basis_element.get_element(basis_field_element) ||
                !product.multiply(multiplier, basis_field_element) ||
                !order->coordinates(flint::FmpqMatRef(product_coords),
                                    product)) {
                return false;
            }
            for (slong j = 0; j < n; ++j) {
                fmpq_set(fmpq_mat_entry(product_rows.raw(), i, j),
                         fmpq_mat_entry(product_coords.raw(), 0, j));
            }
        }
    }

    flint::FmpzMat numerator(0, n);
    flint::Fmpz denominator;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::fractional_ideal,
                "fractional_ideal.principal_times_integral_hnf");
        if (!rational_mat_to_lattice(numerator, denominator,
                                     flint::FmpqMatConstRef(product_rows)) ||
            !flint::fmpz_is_one(flint::FmpzConstRef(denominator)) ||
            fmpz_mat_is_in_hnf(numerator.raw()) == 0) {
            return false;
        }
    }

    Ideal candidate(*order);
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(numerator))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

// Source trace: reference `base4.c:pr_inv_p` returns p P^{-1} as
// `ZM_hnfmodid(tau, p)`, the HNF of p Z_K + t Z_K where tau is the
// multiplication table of t, and `base2.c:get_pr` takes
// t = `FpM_deplin(zk_multable(nf, u), p)`, a kernel vector modulo p of
// multiplication by u for P = p Z_K + u Z_K ("Z_K + t/p Z_K = P^(-1)").
// The reference `anti_uniformizer` keeps the same
// kernel of the uniformizer's representation matrix modulo p as its
// alternative route.  This routine takes the kernel over the HNF rows of I
// instead of one two-element generator, so it needs no stored generator.
//
// For O maximal and p O contained in I, p I^{-1} = {y in O : y I in p O}.
// Its image in O/pO is the common kernel modulo p of y -> y b over the HNF
// rows b of I, and has F_p-dimension f = log_p N(I) because
// [p I^{-1} : p O] = N(I).  Every partial kernel contains that image, so the
// running intersection equals it as soon as its dimension is f.  For t in
// the kernel, p O + t O lies in p I^{-1}; it is equal when its index in O is
// N(p I^{-1}) = p^(n - f), which is checked.  Since I != O, p I^{-1} is not
// contained in p O, so (p I^{-1}, p) is already the normalized
// numerator/denominator pair.
bool prime_ideal_inverse_num_den(Ideal& numerator,
                                 flint::FmpzRef den,
                                 const Ideal& ideal) noexcept {
    const Order* parent = ideal.parent();
    if (parent == nullptr || !ideal.has_hnf() || !parent->is_maximal() ||
        numerator.parent() == nullptr ||
        !numerator.parent()->has_same_data(*parent) || den.raw() == nullptr) {
        return false;
    }

    const slong n = parent->degree();
    flint::FmpzMat hnf(n, n);
    if (n <= 0 || !ideal.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }

    // The diagonal of the HNF of an ideal containing p O is in {1, p}, and
    // N(I) = p^f with f the number of entries p.
    flint::Fmpz p;
    slong f = 0;
    for (slong i = 0; i < n; ++i) {
        const fmpz* diagonal = fmpz_mat_entry(hnf.raw(), i, i);
        if (fmpz_is_one(diagonal) != 0) {
            continue;
        }
        if (f == 0) {
            fmpz_set(p.raw(), diagonal);
        } else if (fmpz_equal(p.raw(), diagonal) == 0) {
            return false;
        }
        ++f;
    }
    if (f == 0 || fmpz_cmp_si(p.raw(), 1) <= 0 ||
        fmpz_is_prime(p.raw()) != 1) {
        return false;
    }

    // Check p O in I: reduce each p e_j by the upper-triangular HNF rows.
    {
        flint::FmpzMat residual(1, n);
        flint::Fmpz quotient;
        for (slong j = 0; j < n; ++j) {
            fmpz_mat_zero(residual.raw());
            fmpz_set(fmpz_mat_entry(residual.raw(), 0, j), p.raw());
            for (slong i = j; i < n; ++i) {
                fmpz* entry = fmpz_mat_entry(residual.raw(), 0, i);
                if (fmpz_is_zero(entry) != 0) {
                    continue;
                }
                const fmpz* pivot = fmpz_mat_entry(hnf.raw(), i, i);
                if (fmpz_divisible(entry, pivot) == 0) {
                    return false;
                }
                fmpz_divexact(quotient.raw(), entry, pivot);
                for (slong c = i; c < n; ++c) {
                    fmpz_submul(fmpz_mat_entry(residual.raw(), 0, c),
                                quotient.raw(),
                                fmpz_mat_entry(hnf.raw(), i, c));
                }
            }
        }
    }

    Order order(*parent);
    // A view of the order's cached table; `order` is not reassigned below.
    const flint::FmpzMat* table_view =
            detail::OrderAccess::multiplication_table(order);
    if (table_view == nullptr) {
        return false;
    }
    const flint::FmpzMat& table = *table_view;

    // Rows of `kernel` span the current kernel in O/pO (order coordinates).
    flint::FmpzModCtx context(p.raw());
    flint::FmpzModMat kernel(n, n, context);
    fmpz_mod_mat_one(kernel.raw(), context.raw());
    bool kernel_is_identity = true;
    slong k = n;
    flint::FmpzMat row(1, n);
    flint::FmpzMat multiplication(n, n);
    for (slong i = 0; i < n && k > f; ++i) {
        bool row_is_zero = true;
        for (slong c = 0; c < n; ++c) {
            fmpz* entry = fmpz_mat_entry(row.raw(), 0, c);
            fmpz_mod(entry, fmpz_mat_entry(hnf.raw(), i, c), p.raw());
            row_is_zero = row_is_zero && fmpz_is_zero(entry) != 0;
        }
        if (row_is_zero) {
            continue;
        }

        // Row r of `multiplication` is the coordinate row of e_r * b.
        multiplication_from_coordinates(flint::FmpzMatRef(multiplication),
                                        flint::FmpzMatConstRef(row),
                                        flint::FmpzMatConstRef(table));
        flint::FmpzModMat image(k, n, context);
        if (kernel_is_identity) {
            fmpz_mod_mat_set_fmpz_mat(image.raw(), multiplication.raw(),
                                      context.raw());
        } else {
            flint::FmpzModMat multiplication_mod(n, n, context);
            fmpz_mod_mat_set_fmpz_mat(multiplication_mod.raw(),
                                      multiplication.raw(), context.raw());
            fmpz_mod_mat_mul(image.raw(), kernel.raw(),
                             multiplication_mod.raw(), context.raw());
        }

        // Left kernel of image: columns of null span {z : z image = 0}.
        flint::FmpzModMat image_transpose(n, k, context);
        flint::FmpzModMat null(k, k, context);
        fmpz_mod_mat_transpose(image_transpose.raw(), image.raw(),
                               context.raw());
        const slong nullity = fmpz_mod_mat_nullspace(
                null.raw(), image_transpose.raw(), context.raw());
        if (nullity < f) {
            return false;
        }
        if (nullity == k) {
            continue;
        }

        flint::FmpzModMat combination(nullity, k, context);
        for (slong a = 0; a < nullity; ++a) {
            for (slong b = 0; b < k; ++b) {
                fmpz_set(fmpz_mod_mat_entry(combination.raw(), a, b),
                         fmpz_mod_mat_entry(null.raw(), b, a));
            }
        }
        flint::FmpzModMat next(nullity, n, context);
        fmpz_mod_mat_mul(next.raw(), combination.raw(), kernel.raw(),
                         context.raw());
        fmpz_mod_mat_swap(kernel.raw(), next.raw(), context.raw());
        kernel_is_identity = false;
        k = nullity;
    }
    if (k != f) {
        return false;
    }

    flint::FmpzMat generator_coordinates(1, n);
    for (slong c = 0; c < n; ++c) {
        fmpz_set(fmpz_mat_entry(generator_coordinates.raw(), 0, c),
                 fmpz_mod_mat_entry(kernel.raw(), 0, c));
    }
    OrderElement generator(order);
    Ideal candidate(*parent);
    flint::FmpzMat candidate_hnf(n, n);
    if (!generator.is_defined() || !candidate.is_defined() ||
        !generator.set_coordinates(
                flint::FmpzMatConstRef(generator_coordinates)) ||
        !IdealAccess::set_known_two_generator(candidate, flint::FmpzConstRef(p),
                                              generator) ||
        !candidate.get_hnf(flint::FmpzMatRef(candidate_hnf))) {
        return false;
    }

    // [O : p O + t O] = p^(n - f) gives p O + t O = p I^{-1}.  Every
    // diagonal entry divides p because p O is contained in p O + t O.
    slong p_count = 0;
    for (slong i = 0; i < n; ++i) {
        const fmpz* diagonal = fmpz_mat_entry(candidate_hnf.raw(), i, i);
        if (fmpz_equal(diagonal, p.raw()) != 0) {
            ++p_count;
        } else if (fmpz_is_one(diagonal) == 0) {
            return false;
        }
    }
    if (p_count != n - f) {
        return false;
    }

    numerator.swap(candidate);
    fmpz_set(den.raw(), p.raw());
    return true;
}

}  // namespace detail

FractionalIdeal::FractionalIdeal(const Order& parent) noexcept {
    define(parent);
}

FractionalIdeal::~FractionalIdeal() noexcept {
    clear();
}

FractionalIdeal::FractionalIdeal(FractionalIdeal&& other) noexcept {
    swap(other);
}

FractionalIdeal& FractionalIdeal::operator=(FractionalIdeal&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FractionalIdeal::swap(FractionalIdeal& other) noexcept {
    parent_.swap(other.parent_);
    numerator_.swap(other.numerator_);
    den_.swap(other.den_);
    norm_.swap(other.norm_);
    std::swap(has_num_, other.has_num_);
    std::swap(has_norm_, other.has_norm_);
}

void FractionalIdeal::clear() noexcept {
    parent_.clear();
    numerator_.clear();
    flint::fmpz_zero(flint::FmpzRef(den_));
    flint::fmpq_zero(flint::FmpqRef(norm_));
    has_num_ = false;
    has_norm_ = false;
}

bool FractionalIdeal::define(const Order& parent) noexcept {
    if (!parent.is_defined() || !parent.has_basis()) {
        return false;
    }

    FractionalIdeal next;
    next.parent_ = parent;
    next.numerator_ = Ideal(parent);
    if (!next.numerator_.is_defined()) {
        return false;
    }
    flint::fmpz_one(flint::FmpzRef(next.den_));

    swap(next);
    return true;
}

bool FractionalIdeal::set(const FractionalIdeal& other) noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_) || !other.has_num_) {
        return false;
    }
    if (this == &other) {
        return true;
    }

    FractionalIdeal copy(parent_);
    if (!copy.is_defined() ||
        !copy.numerator_.set(other.numerator_)) {
        return false;
    }
    flint::fmpz_set(flint::FmpzRef(copy.den_),
                    flint::FmpzConstRef(other.den_));
    flint::fmpq_set(flint::FmpqRef(copy.norm_),
                    flint::FmpqConstRef(other.norm_));
    copy.has_num_ = true;
    copy.has_norm_ = other.has_norm_;
    swap(copy);
    return true;
}

bool FractionalIdeal::is_defined() const noexcept {
    return parent_.is_defined() && numerator_.is_defined();
}

const Order* FractionalIdeal::parent() const noexcept {
    return !is_defined() ? nullptr : &parent_;
}

slong FractionalIdeal::degree() const noexcept {
    return !is_defined() ? 0 : parent_.degree();
}

bool FractionalIdeal::has_integral_denominator() const noexcept {
    return is_defined() && has_num_;
}

bool FractionalIdeal::one() noexcept {
    if (!is_defined() || !numerator_.one()) {
        return false;
    }
    flint::fmpz_one(flint::FmpzRef(den_));
    flint::fmpq_one(flint::FmpqRef(norm_));
    has_num_ = true;
    has_norm_ = true;
    return true;
}

bool FractionalIdeal::set_principal(
        const Element& generator,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::fractional_ideal,
                        "fractional_ideal.set_principal");
    if (!is_defined() || parent_.parent() == nullptr ||
        !generator.has_parent(*parent_.parent()) ||
        generator.equal_si(0)) {
        return false;
    }

    const slong n = parent_.degree();
    flint::FmpqMat coords(1, n);
    flint::FmpzMat row(1, n);
    flint::Fmpz den;
    {
        SILEX_PROFILE_SCOPE(diagnostics,
                            DiagnosticsModule::fractional_ideal,
                            "fractional_ideal.principal_coordinates");
        if (!parent_.coordinates(flint::FmpqMatRef(coords), generator)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics,
                            DiagnosticsModule::fractional_ideal,
                            "fractional_ideal.principal_denominator");
        fmpq_mat_get_fmpz_mat_matwise(row.raw(), den.raw(), coords.raw());
    }

    OrderElement integral_generator(parent_);
    Ideal numerator(parent_);
    FractionalIdeal candidate(parent_);
    if (!integral_generator.is_defined() || !numerator.is_defined() ||
        !candidate.is_defined()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics,
                            DiagnosticsModule::fractional_ideal,
                            "fractional_ideal.principal_integral_generator");
        if (!integral_generator.set_coordinates(
                    flint::FmpzMatConstRef(row))) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics,
                            DiagnosticsModule::fractional_ideal,
                            "fractional_ideal.principal_integral_ideal");
        if (!numerator.set_principal(integral_generator)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics,
                            DiagnosticsModule::fractional_ideal,
                            "fractional_ideal.principal_set_integral_den");
        if (!candidate.set_integral_den(numerator, flint::FmpzConstRef(den))) {
            return false;
        }
    }

    swap(candidate);
    return true;
}

bool FractionalIdeal::set_integral(const Ideal& numerator) noexcept {
    flint::Fmpz one;
    flint::fmpz_one(flint::FmpzRef(one));
    return set_integral_den(numerator, flint::FmpzConstRef(one));
}

bool FractionalIdeal::set_integral_den(const Ideal& numerator,
                                       flint::FmpzConstRef den) noexcept {
    if (!is_defined() || numerator.parent() == nullptr ||
        !numerator.parent()->has_same_data(parent_) ||
        !numerator.has_hnf() || flint::fmpz_sgn(den) <= 0) {
        return false;
    }

    const slong n = parent_.degree();
    flint::FmpzMat hnf(n, n);
    flint::Fmpz g;
    if (!numerator.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }
    common_content(flint::FmpzRef(g), flint::FmpzMatConstRef(hnf), den);

    FractionalIdeal candidate(parent_);
    if (!candidate.is_defined()) {
        return false;
    }
    flint::fmpz_set(flint::FmpzRef(candidate.den_), den);
    if (!flint::fmpz_is_one(g)) {
        fmpz_mat_scalar_divexact_fmpz(hnf.raw(), hnf.raw(), g.raw());
        flint::fmpz_divexact(flint::FmpzRef(candidate.den_),
                             flint::FmpzConstRef(candidate.den_),
                             flint::FmpzConstRef(g));
    }

    // `numerator` is an Ideal holding an HNF, so its lattice I is an ideal:
    // Ideal::set_hnf tests closure, the other routes that store an HNF build
    // an ideal by construction (the private Ideal::set_hnf_direct is reached
    // only from those routes: copy, the unit ideal, the product fast path,
    // aO + bO, and an integral ideal times an element), and outside the
    // construction routes only the test hook IdealAccess::set_hnf_direct
    // stores an unchecked matrix.  Dividing the
    // content g out keeps an ideal: I = gJ as lattices, and for x in J and w
    // in O, g(wx) = w(gx) lies in I, so wx lies in J.  The hnf_is_ideal
    // product test of set_hnf is therefore skipped.
    if (!detail::IdealAccess::set_known_ideal_hnf(
                candidate.numerator_, flint::FmpzMatConstRef(hnf))) {
        return false;
    }

    candidate.has_num_ = true;
    candidate.has_norm_ = false;
    flint::fmpq_zero(flint::FmpqRef(candidate.norm_));
    swap(candidate);
    return true;
}

namespace detail {

bool IdealAccess::set_integral_known_hnf(FractionalIdeal& out,
                                         const Ideal& numerator) noexcept {
    if (!out.is_defined() || numerator.parent() == nullptr ||
        !numerator.parent()->has_same_data(out.parent_) ||
        !numerator.has_hnf()) {
        return false;
    }

    FractionalIdeal candidate(out.parent_);
    if (!candidate.is_defined() || !candidate.numerator_.set(numerator)) {
        return false;
    }

    flint::fmpz_one(flint::FmpzRef(candidate.den_));
    candidate.has_num_ = true;
    candidate.has_norm_ = false;
    flint::fmpq_zero(flint::FmpqRef(candidate.norm_));
    out.swap(candidate);
    return true;
}

}  // namespace detail

bool FractionalIdeal::get_integral_den(Ideal& numerator,
                                       flint::FmpzRef den) const noexcept {
    if (!is_defined() || !has_num_ ||
        numerator.parent() == nullptr ||
        !numerator.parent()->has_same_data(parent_) ||
        !numerator.set(numerator_)) {
        return false;
    }
    flint::fmpz_set(den, flint::FmpzConstRef(den_));
    return true;
}

bool FractionalIdeal::contains(const Element& element) const noexcept {
    if (!is_defined() || !has_num_ ||
        parent_.parent() == nullptr ||
        !element.has_parent(*parent_.parent())) {
        return false;
    }

    Element scaled(*parent_.parent());
    OrderElement order_element(parent_);
    return multiply_element_by_fmpz(scaled, element,
                                    flint::FmpzConstRef(den_)) &&
           order_element.set_element(scaled) &&
           numerator_.contains(order_element);
}

bool FractionalIdeal::norm(flint::FmpqRef out) const noexcept {
    if (!is_defined() || !has_num_) {
        return false;
    }

    if (!has_norm_) {
        flint::Fmpz num_norm;
        flint::Fmpz den_pow;
        if (!numerator_.norm(flint::FmpzRef(num_norm))) {
            return false;
        }
        flint::fmpz_pow_ui(flint::FmpzRef(den_pow),
                           flint::FmpzConstRef(den_),
                           static_cast<ulong>(parent_.degree()));
        flint::fmpq_set_fmpz_frac(flint::FmpqRef(norm_),
                                  flint::FmpzConstRef(num_norm),
                                  flint::FmpzConstRef(den_pow));
        has_norm_ = true;
    }
    flint::fmpq_set(out, flint::FmpqConstRef(norm_));
    return true;
}

bool FractionalIdeal::equal(const FractionalIdeal& other) const noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_) || !has_num_ ||
        !other.has_num_) {
        return false;
    }

    if (fmpz_equal(den_.raw(), other.den_.raw()) != 0) {
        return numerator_.equal(other.numerator_);
    }

    const slong n = parent_.degree();
    flint::FmpzMat left(n, n);
    flint::FmpzMat right(n, n);
    if (!numerator_.get_hnf(flint::FmpzMatRef(left)) ||
        !other.numerator_.get_hnf(flint::FmpzMatRef(right))) {
        return false;
    }

    fmpz_mat_scalar_mul_fmpz(left.raw(), left.raw(), other.den_.raw());
    fmpz_mat_scalar_mul_fmpz(right.raw(), right.raw(), den_.raw());
    fmpz_mat_hnf(left.raw(), left.raw());
    fmpz_mat_hnf(right.raw(), right.raw());
    return fmpz_mat_equal(left.raw(), right.raw()) != 0;
}

bool FractionalIdeal::add(const FractionalIdeal& left,
                          const FractionalIdeal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_num_ ||
        !right.has_num_) {
        return false;
    }

    FractionalIdeal candidate(parent_);
    Ideal left_num(parent_);
    Ideal right_num(parent_);
    Ideal sum(parent_);
    flint::Fmpz den;
    flint::Fmpz left_scale;
    flint::Fmpz right_scale;
    if (!candidate.is_defined() || !left_num.is_defined() ||
        !right_num.is_defined() || !sum.is_defined()) {
        return false;
    }

    flint::fmpz_lcm(flint::FmpzRef(den),
                    flint::FmpzConstRef(left.den_),
                    flint::FmpzConstRef(right.den_));
    flint::fmpz_divexact(flint::FmpzRef(left_scale),
                         flint::FmpzConstRef(den),
                         flint::FmpzConstRef(left.den_));
    flint::fmpz_divexact(flint::FmpzRef(right_scale),
                         flint::FmpzConstRef(den),
                         flint::FmpzConstRef(right.den_));

    if (!scale_integral_ideal(left_num, left.numerator_,
                              flint::FmpzConstRef(left_scale)) ||
        !scale_integral_ideal(right_num, right.numerator_,
                              flint::FmpzConstRef(right_scale)) ||
        !sum.add(left_num, right_num) ||
        !candidate.set_integral_den(sum, flint::FmpzConstRef(den))) {
        return false;
    }

    swap(candidate);
    return true;
}

bool FractionalIdeal::intersect(const FractionalIdeal& left,
                                const FractionalIdeal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_num_ ||
        !right.has_num_) {
        return false;
    }

    FractionalIdeal candidate(parent_);
    Ideal left_num(parent_);
    Ideal right_num(parent_);
    Ideal intersection(parent_);
    flint::Fmpz den;
    flint::Fmpz left_scale;
    flint::Fmpz right_scale;
    if (!candidate.is_defined() || !left_num.is_defined() ||
        !right_num.is_defined() || !intersection.is_defined()) {
        return false;
    }

    flint::fmpz_lcm(flint::FmpzRef(den),
                    flint::FmpzConstRef(left.den_),
                    flint::FmpzConstRef(right.den_));
    flint::fmpz_divexact(flint::FmpzRef(left_scale),
                         flint::FmpzConstRef(den),
                         flint::FmpzConstRef(left.den_));
    flint::fmpz_divexact(flint::FmpzRef(right_scale),
                         flint::FmpzConstRef(den),
                         flint::FmpzConstRef(right.den_));

    if (!scale_integral_ideal(left_num, left.numerator_,
                              flint::FmpzConstRef(left_scale)) ||
        !scale_integral_ideal(right_num, right.numerator_,
                              flint::FmpzConstRef(right_scale)) ||
        !intersection.intersect(left_num, right_num) ||
        !candidate.set_integral_den(intersection, flint::FmpzConstRef(den))) {
        return false;
    }

    swap(candidate);
    return true;
}

bool FractionalIdeal::multiply(const FractionalIdeal& left,
                               const FractionalIdeal& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_) || !left.has_num_ ||
        !right.has_num_) {
        return false;
    }

    FractionalIdeal candidate(parent_);
    Ideal product(parent_);
    flint::Fmpz den;
    if (!candidate.is_defined() || !product.is_defined()) {
        return false;
    }

    flint::fmpz_mul(flint::FmpzRef(den),
                    flint::FmpzConstRef(left.den_),
                    flint::FmpzConstRef(right.den_));
    if (!product.multiply(left.numerator_, right.numerator_) ||
        !candidate.set_integral_den(product, flint::FmpzConstRef(den))) {
        return false;
    }

    swap(candidate);
    return true;
}

bool FractionalIdeal::pow_fmpz(const FractionalIdeal& input,
                               flint::FmpzConstRef exponent) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.has_same_data(input.parent_) || !input.has_num_) {
        return false;
    }

    FractionalIdeal result(parent_);
    FractionalIdeal base(parent_);
    flint::Fmpz n;
    if (!result.one()) {
        return false;
    }

    if (flint::fmpz_sgn(exponent) < 0) {
        if (!base.invert(input)) {
            return false;
        }
        flint::fmpz_neg(flint::FmpzRef(n), exponent);
    } else {
        if (!base.set(input)) {
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(n), exponent);
    }

    while (!flint::fmpz_is_zero(n)) {
        if (flint::fmpz_tstbit(flint::FmpzConstRef(n), 0) &&
            !result.multiply(result, base)) {
            return false;
        }

        flint::fmpz_fdiv_q_2exp(flint::FmpzRef(n),
                                flint::FmpzConstRef(n), 1);
        if (!flint::fmpz_is_zero(n) && !base.multiply(base, base)) {
            return false;
        }
    }

    swap(result);
    return true;
}

bool FractionalIdeal::colon(const FractionalIdeal& numerator,
                            const FractionalIdeal& denominator) noexcept {
    if (!is_defined() || !numerator.is_defined() ||
        !denominator.is_defined() ||
        !parent_.has_same_data(numerator.parent_) ||
        !numerator.parent_.has_same_data(denominator.parent_) ||
        !numerator.has_num_ || !denominator.has_num_) {
        return false;
    }

    const slong n = parent_.degree();
    FractionalIdeal candidate(parent_);
    Ideal integral_num(parent_);
    flint::FmpzMat numerator_hnf(n, n);
    flint::FmpzMat denominator_hnf(n, n);
    flint::FmpzMat current(0, n);
    flint::FmpzMat preimage(0, n);
    flint::Fmpz current_den;
    flint::Fmpz preimage_den;
    flint::Fmpz den;
    OrderElement basis_elem(parent_);
    flint::FmpzMat row(1, n);
    if (!candidate.is_defined() || !integral_num.is_defined() ||
        !basis_elem.is_defined() ||
        !numerator.numerator_.get_hnf(
                flint::FmpzMatRef(numerator_hnf)) ||
        !denominator.numerator_.get_hnf(
                flint::FmpzMatRef(denominator_hnf))) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        copy_row(flint::FmpzMatRef(row), 0,
                 flint::FmpzMatConstRef(denominator_hnf), i);
        if (!basis_elem.set_coordinates(flint::FmpzMatConstRef(row)) ||
            !preimage_under_mul(preimage, preimage_den,
                                flint::FmpzMatConstRef(numerator_hnf),
                                parent_, basis_elem)) {
            return false;
        }

        if (i == 0) {
            current = std::move(preimage);
            flint::fmpz_set(flint::FmpzRef(current_den),
                            flint::FmpzConstRef(preimage_den));
        } else if (!intersect_rational_lattices(
                           current, current_den,
                           flint::FmpzMatConstRef(current),
                           flint::FmpzConstRef(current_den),
                           flint::FmpzMatConstRef(preimage),
                           flint::FmpzConstRef(preimage_den))) {
            return false;
        }
    }

    fmpz_mat_scalar_mul_fmpz(current.raw(), current.raw(),
                             denominator.den_.raw());
    flint::fmpz_mul(flint::FmpzRef(den),
                    flint::FmpzConstRef(current_den),
                    flint::FmpzConstRef(numerator.den_));
    if (!normalize_lattice_num_den(current, den) ||
        !integral_num.set_hnf(flint::FmpzMatConstRef(current)) ||
        !candidate.set_integral_den(integral_num, flint::FmpzConstRef(den))) {
        return false;
    }

    swap(candidate);
    return true;
}

bool FractionalIdeal::invert(const FractionalIdeal& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.has_same_data(input.parent_) || !input.has_num_) {
        return false;
    }

    FractionalIdeal one(parent_);
    FractionalIdeal colon_result(parent_);
    if (!one.one()) {
        return false;
    }

    // reference's generic colon(O, I) route computes the inverse over known
    // maximal orders directly from representation matrices.  Keep the older
    // checked colon path for nonmaximal/unknown orders, where colon need not
    // be an inverse.
    if (parent_.is_maximal()) {
        if (!colon_one_integral_hnf_inverse(colon_result,
                                                  input.numerator_)) {
            return false;
        }

        if (!flint::fmpz_is_one(input.den_)) {
            Ideal scaled(parent_);
            FractionalIdeal scaled_fractional(parent_);
            flint::Fmpz direct_den;
            if (!scaled.is_defined() || !scaled_fractional.is_defined() ||
                !colon_result.get_integral_den(scaled,
                                               flint::FmpzRef(direct_den)) ||
                !scale_integral_ideal(scaled, scaled,
                                      flint::FmpzConstRef(input.den_)) ||
                !scaled_fractional.set_integral_den(
                        scaled, flint::FmpzConstRef(direct_den))) {
                return false;
            }
            colon_result.swap(scaled_fractional);
        }
    } else {
        if (!colon_result.colon(one, input)) {
            return false;
        }

        FractionalIdeal product(parent_);
        if (!product.multiply(input, colon_result) || !product.equal(one)) {
            return false;
        }
    }

    swap(colon_result);
    return true;
}

}  // namespace silex
