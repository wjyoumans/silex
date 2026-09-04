#include <silex/residue_ring.hpp>

#include <flint/fmpz_mat.h>

#include <utility>

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

bool set_row_si(flint::FmpzMatRef row, slong value) noexcept {
    if (row.raw() == nullptr || fmpz_mat_nrows(row.raw()) != 1 ||
        fmpz_mat_ncols(row.raw()) <= 0) {
        return false;
    }
    fmpz_mat_zero(row.raw());
    fmpz_set_si(fmpz_mat_entry(row.raw(), 0, 0), value);
    return true;
}

void set_identity_row(flint::FmpzMatRef row, slong index) noexcept {
    fmpz_mat_zero(row.raw());
    fmpz_one(fmpz_mat_entry(row.raw(), 0, index));
}

bool build_multiplication_table(flint::FmpzMatRef out,
                                const Order& order) noexcept {
    const slong n = order.degree();
    if (!fmpz_mat_shape(out, n * n, n)) {
        return false;
    }

    OrderElement left(order);
    OrderElement right(order);
    OrderElement product(order);
    flint::FmpzMat row(1, n);
    flint::FmpzMat coords(1, n);
    if (!left.is_defined() || !right.is_defined() || !product.is_defined()) {
        return false;
    }

    for (slong i = 0; i < n; ++i) {
        set_identity_row(flint::FmpzMatRef(row), i);
        if (!left.set_coordinates(flint::FmpzMatConstRef(row))) {
            return false;
        }

        for (slong j = 0; j < n; ++j) {
            set_identity_row(flint::FmpzMatRef(row), j);
            if (!right.set_coordinates(flint::FmpzMatConstRef(row)) ||
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

}  // namespace

ResidueRing::ResidueRing(const Ideal& modulus) noexcept {
    define(modulus);
}

ResidueRing::~ResidueRing() noexcept = default;

ResidueRing::ResidueRing(ResidueRing&& other) noexcept {
    swap(other);
}

ResidueRing& ResidueRing::operator=(ResidueRing&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueRing::swap(ResidueRing& other) noexcept {
    modulus_.swap(other.modulus_);
    multiplication_table_.swap(other.multiplication_table_);
}

void ResidueRing::clear() noexcept {
    modulus_.clear();
    multiplication_table_ = flint::FmpzMat(0, 0);
}

bool ResidueRing::define(const Ideal& modulus) noexcept {
    if (modulus.parent() == nullptr || !modulus.has_hnf()) {
        return false;
    }

    ResidueRing candidate;
    candidate.modulus_ = Ideal(*modulus.parent());
    candidate.multiplication_table_ =
            flint::FmpzMat(modulus.degree() * modulus.degree(),
                           modulus.degree());
    if (!candidate.modulus_.is_defined() ||
        !candidate.modulus_.set(modulus) ||
        !build_multiplication_table(
                flint::FmpzMatRef(candidate.multiplication_table_),
                *modulus.parent())) {
        return false;
    }

    swap(candidate);
    return true;
}

bool ResidueRing::set(const ResidueRing& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined() || other.modulus_.parent() == nullptr) {
        clear();
        return true;
    }

    ResidueRing candidate;
    candidate.modulus_ = Ideal(*other.modulus_.parent());
    candidate.multiplication_table_ =
            flint::FmpzMat(
                    flint::fmpz_mat_nrows(other.multiplication_table_),
                    flint::fmpz_mat_ncols(other.multiplication_table_));
    if (!candidate.modulus_.is_defined() ||
        !candidate.modulus_.set(other.modulus_)) {
        return false;
    }
    fmpz_mat_set(candidate.multiplication_table_.raw(),
                 other.multiplication_table_.raw());

    swap(candidate);
    return true;
}

bool ResidueRing::is_defined() const noexcept {
    return modulus_.has_hnf();
}

const Order* ResidueRing::parent_order() const noexcept {
    return modulus_.parent();
}

const Ideal* ResidueRing::modulus() const noexcept {
    return is_defined() ? &modulus_ : nullptr;
}

slong ResidueRing::degree() const noexcept {
    return is_defined() ? modulus_.degree() : 0;
}

bool ResidueRing::get_modulus(Ideal& out) const noexcept {
    return is_defined() && same_order_parent(out.parent(), modulus_.parent()) &&
           out.set(modulus_);
}

bool ResidueRing::cardinality(flint::FmpzRef out) const noexcept {
    return is_defined() && modulus_.norm(out);
}

std::optional<flint::Fmpz> ResidueRing::cardinality() const noexcept {
    flint::Fmpz out;
    if (!cardinality(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ResidueRing::equal(const ResidueRing& other) const noexcept {
    return is_defined() && other.is_defined() && modulus_.equal(other.modulus_);
}

bool ResidueRing::reduce_row(flint::FmpzMatRef row) const noexcept {
    if (!is_defined() || !fmpz_mat_shape(row, 1, modulus_.degree())) {
        return false;
    }

    const slong n = modulus_.degree();
    flint::FmpzMat hnf(n, n);
    if (!modulus_.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }

    flint::Fmpz q;
    for (slong i = 0; i < n; ++i) {
        const fmpz* pivot = fmpz_mat_entry(hnf.raw(), i, i);
        fmpz_fdiv_q(q.raw(), fmpz_mat_entry(row.raw(), 0, i), pivot);

        if (fmpz_is_zero(q.raw()) == 0) {
            for (slong j = i; j < n; ++j) {
                fmpz_submul(fmpz_mat_entry(row.raw(), 0, j), q.raw(),
                            fmpz_mat_entry(hnf.raw(), i, j));
            }
        }
    }

    return true;
}

bool ResidueRing::multiply_rows(flint::FmpzMatRef out,
                                flint::FmpzMatConstRef left,
                                flint::FmpzMatConstRef right) const noexcept {
    if (!is_defined() || !fmpz_mat_shape(out, 1, modulus_.degree()) ||
        !fmpz_mat_shape(left, 1, modulus_.degree()) ||
        !fmpz_mat_shape(right, 1, modulus_.degree())) {
        return false;
    }

    const slong n = modulus_.degree();
    flint::Fmpz ab;
    fmpz_mat_zero(out.raw());
    for (slong i = 0; i < n; ++i) {
        if (fmpz_is_zero(fmpz_mat_entry(left.raw(), 0, i)) != 0) {
            continue;
        }
        for (slong j = 0; j < n; ++j) {
            if (fmpz_is_zero(fmpz_mat_entry(right.raw(), 0, j)) != 0) {
                continue;
            }

            fmpz_mul(ab.raw(), fmpz_mat_entry(left.raw(), 0, i),
                     fmpz_mat_entry(right.raw(), 0, j));
            for (slong k = 0; k < n; ++k) {
                fmpz_addmul(
                        fmpz_mat_entry(out.raw(), 0, k), ab.raw(),
                        fmpz_mat_entry(multiplication_table_.raw(),
                                       i * n + j, k));
            }
        }
    }

    return true;
}

ResidueElement::ResidueElement(const ResidueRing& parent) noexcept {
    define(parent);
}

ResidueElement::~ResidueElement() noexcept = default;

ResidueElement::ResidueElement(ResidueElement&& other) noexcept {
    swap(other);
}

ResidueElement& ResidueElement::operator=(ResidueElement&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueElement::swap(ResidueElement& other) noexcept {
    parent_.swap(other.parent_);
    coordinates_.swap(other.coordinates_);
}

void ResidueElement::clear() noexcept {
    parent_.clear();
    coordinates_ = flint::FmpzMat(0, 0);
}

bool ResidueElement::define(const ResidueRing& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    ResidueElement candidate;
    if (!candidate.parent_.set(parent)) {
        return false;
    }
    candidate.coordinates_ = flint::FmpzMat(1, parent.degree());
    swap(candidate);
    return true;
}

bool ResidueElement::set(const ResidueElement& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!is_defined() || !other.is_defined() ||
        !parent_.equal(other.parent_)) {
        return false;
    }

    fmpz_mat_set(coordinates_.raw(), other.coordinates_.raw());
    return true;
}

bool ResidueElement::is_defined() const noexcept {
    return parent_.is_defined();
}

const ResidueRing* ResidueElement::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong ResidueElement::degree() const noexcept {
    return parent_.degree();
}

bool ResidueElement::zero() noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_mat_zero(coordinates_.raw());
    return true;
}

bool ResidueElement::one() noexcept {
    if (!is_defined()) {
        return false;
    }

    flint::FmpzMat row(1, parent_.degree());
    ResidueElement candidate(parent_);
    if (!candidate.is_defined() || !set_row_si(flint::FmpzMatRef(row), 1) ||
        !parent_.reduce_row(flint::FmpzMatRef(row))) {
        return false;
    }
    fmpz_mat_set(candidate.coordinates_.raw(), row.raw());
    swap(candidate);
    return true;
}

bool ResidueElement::set_order_element(const OrderElement& element) noexcept {
    if (!is_defined() ||
        !same_order_parent(element.parent(), parent_.parent_order())) {
        return false;
    }

    flint::FmpzMat row(1, parent_.degree());
    ResidueElement candidate(parent_);
    if (!candidate.is_defined() ||
        !element.get_coordinates(flint::FmpzMatRef(row)) ||
        !parent_.reduce_row(flint::FmpzMatRef(row))) {
        return false;
    }

    fmpz_mat_set(candidate.coordinates_.raw(), row.raw());
    swap(candidate);
    return true;
}

bool ResidueElement::set_element(const Element& element) noexcept {
    if (!is_defined() ||
        parent_.parent_order() == nullptr ||
        parent_.parent_order()->parent() == nullptr ||
        !element.has_parent(*parent_.parent_order()->parent())) {
        return false;
    }

    OrderElement order_element(*parent_.parent_order());
    return order_element.set_element(element) &&
           set_order_element(order_element);
}

bool ResidueElement::get_coordinates(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !fmpz_mat_shape(out, 1, parent_.degree())) {
        return false;
    }
    fmpz_mat_set(out.raw(), coordinates_.raw());
    return true;
}

std::optional<flint::FmpzMat> ResidueElement::coordinates() const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(1, parent_.degree());
    if (!get_coordinates(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ResidueElement::lift(OrderElement& out) const noexcept {
    return is_defined() &&
           same_order_parent(out.parent(), parent_.parent_order()) &&
           out.set_coordinates(flint::FmpzMatConstRef(coordinates_));
}

bool ResidueElement::lift(Element& out) const noexcept {
    if (!is_defined() ||
        parent_.parent_order() == nullptr ||
        parent_.parent_order()->parent() == nullptr ||
        !out.has_parent(*parent_.parent_order()->parent())) {
        return false;
    }

    OrderElement order_element(*parent_.parent_order());
    return lift(order_element) && order_element.get_element(out);
}

bool ResidueElement::equal(const ResidueElement& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.equal(other.parent_) &&
           fmpz_mat_equal(coordinates_.raw(), other.coordinates_.raw()) != 0;
}

bool ResidueElement::add(const ResidueElement& left,
                         const ResidueElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }

    ResidueElement candidate(parent_);
    if (!candidate.is_defined()) {
        return false;
    }
    fmpz_mat_add(candidate.coordinates_.raw(), left.coordinates_.raw(),
                 right.coordinates_.raw());
    if (!parent_.reduce_row(flint::FmpzMatRef(candidate.coordinates_))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ResidueElement::negate(const ResidueElement& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.equal(input.parent_)) {
        return false;
    }

    ResidueElement candidate(parent_);
    if (!candidate.is_defined()) {
        return false;
    }
    fmpz_mat_neg(candidate.coordinates_.raw(), input.coordinates_.raw());
    if (!parent_.reduce_row(flint::FmpzMatRef(candidate.coordinates_))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ResidueElement::subtract(const ResidueElement& left,
                              const ResidueElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }

    ResidueElement candidate(parent_);
    if (!candidate.is_defined()) {
        return false;
    }
    fmpz_mat_sub(candidate.coordinates_.raw(), left.coordinates_.raw(),
                 right.coordinates_.raw());
    if (!parent_.reduce_row(flint::FmpzMatRef(candidate.coordinates_))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ResidueElement::multiply(const ResidueElement& left,
                              const ResidueElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }

    ResidueElement candidate(parent_);
    if (!candidate.is_defined() ||
        !parent_.multiply_rows(
                flint::FmpzMatRef(candidate.coordinates_),
                flint::FmpzMatConstRef(left.coordinates_),
                flint::FmpzMatConstRef(right.coordinates_)) ||
        !parent_.reduce_row(flint::FmpzMatRef(candidate.coordinates_))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ResidueElement::crt(const ResidueElement& left,
                         const ResidueElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !same_order_parent(left.parent_.parent_order(),
                           right.parent_.parent_order()) ||
        !same_order_parent(parent_.parent_order(),
                           left.parent_.parent_order())) {
        return false;
    }

    Ideal product(*parent_.parent_order());
    if (!product.multiply(*left.parent_.modulus(),
                          *right.parent_.modulus()) ||
        !product.equal(*parent_.modulus())) {
        return false;
    }

    OrderElement u(*parent_.parent_order());
    OrderElement v(*parent_.parent_order());
    if (!left.parent_.modulus()->add_to_one(u, v,
                                            *right.parent_.modulus())) {
        return false;
    }

    OrderElement left_lift(*parent_.parent_order());
    OrderElement right_lift(*parent_.parent_order());
    OrderElement left_part(*parent_.parent_order());
    OrderElement right_part(*parent_.parent_order());
    OrderElement combined(*parent_.parent_order());
    ResidueElement candidate(parent_);
    if (!left.lift(left_lift) || !right.lift(right_lift) ||
        !left_part.multiply(left_lift, v) ||
        !right_part.multiply(right_lift, u) ||
        !combined.add(left_part, right_part) ||
        !candidate.set_order_element(combined)) {
        return false;
    }

    swap(candidate);
    return true;
}

}  // namespace silex
