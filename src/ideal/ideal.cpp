#include <silex/ideal.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/lat.hpp>

#include <utility>
#include <vector>

namespace silex {
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

bool order_element_multiplication_matrix(
        flint::FmpzMatRef out,
        Order& order,
        const OrderElement& element) noexcept {
    const Order* element_parent = element.parent();
    const slong n = order.degree();
    if (element_parent == nullptr ||
        !element_parent->has_same_data(order) ||
        !fmpz_mat_shape(out, n, n)) {
        return false;
    }

    flint::FmpzMat element_coords(1, n);
    flint::FmpzMat multiplication_table(n * n, n);
    if (!element.get_coordinates(flint::FmpzMatRef(element_coords)) ||
        !order.multiplication_table(
                flint::FmpzMatRef(multiplication_table))) {
        return false;
    }

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
    return true;
}

void set_identity_row(flint::FmpzMatRef row, slong index) noexcept {
    fmpz_mat_zero(row.raw());
    fmpz_one(fmpz_mat_entry(row.raw(), 0, index));
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
    if (!gen.is_defined() || !omega.is_defined() || !product.is_defined()) {
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
                !hnf_lattice_contains(hnf, flint::FmpzMatConstRef(row))) {
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
    return out.set_hnf(flint::FmpzMatConstRef(hnf));
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
    if (!rational_mat_to_lattice(numerator_hnf, denominator_value,
                                 flint::FmpqMatConstRef(dual)) ||
        !numerator.set_hnf(flint::FmpzMatConstRef(numerator_hnf)) ||
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
    if (!is_defined()) {
        return false;
    }

    const slong n = parent_.degree();
    flint::FmpzMat normalized(n, n);
    if (!normalize_full_rank_hnf(flint::FmpzMatRef(normalized), hnf) ||
        !hnf_is_ideal(parent_, flint::FmpzMatConstRef(normalized))) {
        return false;
    }

    Ideal candidate(parent_);
    if (!candidate.set_hnf_direct(flint::FmpzMatConstRef(normalized))) {
        return false;
    }
    swap(candidate);
    return true;
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
    if (!is_defined() || generator.parent() == nullptr ||
        !generator.parent()->has_same_data(parent_) ||
        generator.equal_si(0)) {
        return false;
    }

    // Source trace: reference `base4.c:idealmulelt` and `idealhnf_two` form HNF
    // from multiplication-table rows.  For a principal ideal, those rows are
    // the multiplication matrix of the generator in the order basis.
    const slong n = parent_.degree();
    flint::FmpzMat rows(n, n);
    if (!order_element_multiplication_matrix(flint::FmpzMatRef(rows), parent_,
                                             generator)) {
        return false;
    }

    Ideal candidate(parent_);
    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows)) ||
        !candidate.principal_generator_.define(parent_) ||
        !candidate.principal_generator_.set(generator)) {
        return false;
    }
    candidate.has_principal_generator_ = true;
    swap(candidate);
    return true;
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

    flint::FmpzMat row(1, parent_.degree());
    return element.get_coordinates(flint::FmpzMatRef(row)) &&
           hnf_lattice_contains(flint::FmpzMatConstRef(hnf_),
                                flint::FmpzMatConstRef(row));
}

bool Ideal::contains(const Ideal& ideal) const noexcept {
    if (!is_defined() || !ideal.is_defined() ||
        !parent_.has_same_data(ideal.parent_) || !has_hnf_ ||
        !ideal.has_hnf_) {
        return false;
    }

    const slong n = degree();
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

    if (!candidate.set_rows(flint::FmpzMatConstRef(rows))) {
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

    return set_hnf_direct(hnf_lattice.basis_ref());
}

namespace detail {

bool set_known_two_generator_ideal(
        Ideal& out,
        flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator) noexcept {
    const Order* element_parent = element_generator.parent();
    if (!out.is_defined() || element_parent == nullptr ||
        !element_parent->has_same_data(out.parent_) ||
        flint::fmpz_sgn(scalar_generator) <= 0) {
        return false;
    }

    // Source trace: reference `Ideal.jl:assure_has_basis_matrix` builds a
    // two-element ideal basis from `representation_matrix_mod(gen_two, m)`,
    // and reference `base4.c:idealhnf_two` computes HNF from the multiplication
    // table of beta together with the scalar generator.
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
    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

bool multiply_integral_ideal_by_two_generator(
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

    flint::FmpzMatWindow scalar_rows(rows, 0, 0, n, n);
    flint::FmpzMatWindow element_rows(rows, n, 0, 2 * n, n);
    fmpz_mat_scalar_mul_fmpz(scalar_rows.raw(), ideal.hnf_.raw(),
                             scalar_generator.raw());
    fmpz_mat_mul(element_rows.raw(), ideal.hnf_.raw(),
                 element_multiplication.raw());

    Ideal candidate(out.parent_);
    if (!candidate.set_known_ideal_rows(flint::FmpzMatConstRef(rows))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

bool multiply_integral_ideal_by_element(
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

    if (!candidate.numerator_.set_hnf(flint::FmpzMatConstRef(hnf))) {
        return false;
    }

    candidate.has_num_ = true;
    candidate.has_norm_ = false;
    flint::fmpq_zero(flint::FmpqRef(candidate.norm_));
    swap(candidate);
    return true;
}

namespace detail {

bool set_integral_ideal_known_hnf(FractionalIdeal& out,
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
