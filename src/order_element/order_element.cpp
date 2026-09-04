#include <silex/order_element.hpp>

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>

#include <utility>

namespace silex {
namespace {

bool fmpz_mat_shape(flint::FmpzMatConstRef matrix,
                    slong rows,
                    slong cols) noexcept {
    return fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

bool fmpz_mat_shape(flint::FmpzMatRef matrix, slong rows, slong cols) noexcept {
    return fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

void set_element_from_power_row(Element& out,
                                flint::FmpqMatConstRef power_row) noexcept {
    flint::FmpqPoly polynomial;
    for (slong j = 0; j < fmpq_mat_ncols(power_row.raw()); ++j) {
        fmpq_poly_set_coeff_fmpq(polynomial.raw(), j,
                                 fmpq_mat_entry(power_row.raw(), 0, j));
    }
    out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

}  // namespace

OrderElement::OrderElement(const Order& parent) noexcept {
    define(parent);
}

OrderElement::~OrderElement() noexcept {
    clear();
}

OrderElement::OrderElement(OrderElement&& other) noexcept {
    swap(other);
}

OrderElement& OrderElement::operator=(OrderElement&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void OrderElement::swap(OrderElement& other) noexcept {
    parent_.swap(other.parent_);
    value_.swap(other.value_);
    coords_.swap(other.coords_);
}

void OrderElement::clear() noexcept {
    parent_.clear();
    value_.clear();
    coords_ = flint::FmpzMat(0, 0);
}

bool OrderElement::define(const Order& parent) noexcept {
    if (!parent.is_defined() || !parent.has_basis() || parent.parent() == nullptr) {
        return false;
    }

    OrderElement next;
    next.parent_ = parent;
    next.value_ = Element(*parent.parent());
    if (!next.value_.is_defined()) {
        return false;
    }
    next.coords_ = flint::FmpzMat(1, parent.degree());
    if (!next.zero()) {
        return false;
    }

    swap(next);
    return true;
}

bool OrderElement::set(const OrderElement& other) noexcept {
    if (!is_defined() || !other.is_defined() ||
        !parent_.has_same_data(other.parent_)) {
        return false;
    }
    if (this == &other) {
        return true;
    }

    if (!value_.set(other.value_)) {
        return false;
    }
    fmpz_mat_set(coords_.raw(), other.coords_.raw());
    return true;
}

bool OrderElement::is_defined() const noexcept {
    return parent_.is_defined() && value_.is_defined();
}

const Order* OrderElement::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong OrderElement::degree() const noexcept {
    return parent_.is_defined() ? parent_.degree() : 0;
}

bool OrderElement::zero() noexcept {
    if (!is_defined()) {
        return false;
    }
    if (!value_.zero()) {
        return false;
    }
    fmpz_mat_zero(coords_.raw());
    return true;
}

bool OrderElement::one() noexcept {
    if (!is_defined()) {
        return false;
    }

    Element one(*parent_.parent());
    if (!one.one()) {
        return false;
    }
    return set_element(one);
}

bool OrderElement::set_si(slong value) noexcept {
    if (!is_defined()) {
        return false;
    }

    Element element(*parent_.parent());
    if (!element.set_si(value)) {
        return false;
    }
    return set_element(element);
}

bool OrderElement::set_element(const Element& element) noexcept {
    if (!is_defined() || parent_.parent() == nullptr ||
        !element.has_parent(*parent_.parent())) {
        return false;
    }

    flint::FmpqMat rational_coords(1, parent_.degree());
    flint::FmpzMat integral_coords(1, parent_.degree());
    if (!parent_.coordinates(flint::FmpqMatRef(rational_coords), element) ||
        fmpq_mat_get_fmpz_mat(integral_coords.raw(), rational_coords.raw()) == 0) {
        return false;
    }

    Element next_value(*parent_.parent());
    if (!next_value.set(element)) {
        return false;
    }

    value_.swap(next_value);
    fmpz_mat_set(coords_.raw(), integral_coords.raw());
    return true;
}

bool OrderElement::set_coordinates(flint::FmpzMatConstRef coordinates) noexcept {
    if (!is_defined() || !fmpz_mat_shape(coordinates, 1, parent_.degree())) {
        return false;
    }

    flint::FmpqMat basis(parent_.degree(), parent_.degree());
    flint::FmpqMat power_row(1, parent_.degree());
    if (!parent_.get_basis(flint::FmpqMatRef(basis))) {
        return false;
    }

    fmpq_mat_mul_r_fmpz_mat(power_row.raw(), coordinates.raw(), basis.raw());

    Element next_value(*parent_.parent());
    set_element_from_power_row(next_value, flint::FmpqMatConstRef(power_row));

    value_.swap(next_value);
    fmpz_mat_set(coords_.raw(), coordinates.raw());
    return true;
}

bool OrderElement::get_element(Element& out) const noexcept {
    if (!is_defined() || parent_.parent() == nullptr ||
        !out.has_parent(*parent_.parent())) {
        return false;
    }
    return out.set(value_);
}

bool OrderElement::get_coordinates(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !fmpz_mat_shape(out, 1, parent_.degree())) {
        return false;
    }
    fmpz_mat_set(out.raw(), coords_.raw());
    return true;
}

std::optional<flint::FmpzMat> OrderElement::coordinates() const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(1, parent_.degree());
    if (!get_coordinates(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool OrderElement::negate(const OrderElement& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.has_same_data(input.parent_)) {
        return false;
    }

    Element next_value(*parent_.parent());
    if (!next_value.negate(input.value_)) {
        return false;
    }
    value_.swap(next_value);
    fmpz_mat_neg(coords_.raw(), input.coords_.raw());
    return true;
}

bool OrderElement::add(const OrderElement& left,
                       const OrderElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_)) {
        return false;
    }

    Element next_value(*parent_.parent());
    flint::FmpzMat next_coords(1, parent_.degree());
    if (!next_value.add(left.value_, right.value_)) {
        return false;
    }
    fmpz_mat_add(next_coords.raw(), left.coords_.raw(), right.coords_.raw());

    value_.swap(next_value);
    fmpz_mat_set(coords_.raw(), next_coords.raw());
    return true;
}

bool OrderElement::subtract(const OrderElement& left,
                            const OrderElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_)) {
        return false;
    }

    Element next_value(*parent_.parent());
    flint::FmpzMat next_coords(1, parent_.degree());
    if (!next_value.subtract(left.value_, right.value_)) {
        return false;
    }
    fmpz_mat_sub(next_coords.raw(), left.coords_.raw(), right.coords_.raw());

    value_.swap(next_value);
    fmpz_mat_set(coords_.raw(), next_coords.raw());
    return true;
}

bool OrderElement::multiply(const OrderElement& left,
                            const OrderElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_)) {
        return false;
    }

    Element product(*parent_.parent());
    if (!product.multiply(left.value_, right.value_)) {
        return false;
    }
    return set_element(product);
}

bool OrderElement::equal(const OrderElement& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.has_same_data(other.parent_) &&
           fmpz_mat_equal(coords_.raw(), other.coords_.raw()) != 0;
}

bool OrderElement::equal_si(slong value) const noexcept {
    return is_defined() && value_.equal_si(value);
}

}  // namespace silex
