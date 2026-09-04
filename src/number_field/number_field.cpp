#include <silex/number_field.hpp>

#include <flint/fmpz_factor.h>

#include <silex/flint/fmpz_factor.hpp>

#include <utility>

namespace silex {
namespace detail {

struct NumberFieldData {
    flint::Nf field;
    NumberFieldBackendKind backend_kind = NumberFieldBackendKind::generic;
    flint::Fmpz quadratic_radicand;
    bool has_quadratic_radicand = false;
};

}  // namespace detail

namespace {

bool degree_is_positive(flint::FmpqPolyConstRef polynomial) noexcept {
    return fmpq_poly_degree(polynomial.raw()) >= 1;
}

bool is_squarefree_nonzero(flint::FmpzConstRef value) noexcept {
    if (fmpz_is_zero(value.raw()) != 0) {
        return false;
    }

    flint::FmpzFactor factorization;
    fmpz_factor(factorization.raw(), value.raw());
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (factorization.raw()->exp[i] > 1) {
            return false;
        }
    }
    return true;
}

bool is_valid_quadratic_radicand(flint::FmpzConstRef radicand) noexcept {
    return fmpz_is_square(radicand.raw()) == 0 &&
           is_squarefree_nonzero(radicand);
}

bool canonical_quadratic_radicand(
        flint::FmpzRef out, flint::FmpqPolyConstRef polynomial) noexcept {
    if (flint::fmpq_poly_degree(polynomial) != 2) {
        return false;
    }

    flint::Fmpq coefficient;
    flint::fmpq_poly_get_coeff_fmpq(flint::FmpqRef(coefficient), polynomial, 2);
    if (!flint::fmpq_equal_si(coefficient, 1)) {
        return false;
    }

    flint::fmpq_poly_get_coeff_fmpq(flint::FmpqRef(coefficient), polynomial, 1);
    if (!flint::fmpq_equal_si(coefficient, 0)) {
        return false;
    }

    flint::fmpq_poly_get_coeff_fmpq(flint::FmpqRef(coefficient), polynomial, 0);
    if (!flint::fmpz_is_one(flint::fmpq_den_ref(coefficient))) {
        return false;
    }

    flint::fmpz_neg(out, flint::fmpq_num_ref(coefficient));
    return is_valid_quadratic_radicand(flint::FmpzConstRef(out.raw()));
}

std::shared_ptr<detail::NumberFieldData> make_field_data(
        flint::FmpqPolyConstRef polynomial) noexcept {
    auto next = std::make_shared<detail::NumberFieldData>();
    if (!next->field.define(polynomial.raw())) {
        return {};
    }
    return next;
}

void set_quadratic_backend(
        detail::NumberFieldData& data, flint::FmpzConstRef radicand) noexcept {
    data.backend_kind = NumberFieldBackendKind::quadratic;
    data.has_quadratic_radicand = true;
    flint::fmpz_set(flint::FmpzRef(data.quadratic_radicand), radicand);
}

}  // namespace

NumberField::~NumberField() noexcept = default;

NumberField::NumberField(const NumberField& other) noexcept = default;

NumberField& NumberField::operator=(const NumberField& other) noexcept = default;

NumberField::NumberField(NumberField&& other) noexcept {
    swap(other);
}

NumberField& NumberField::operator=(NumberField&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void NumberField::swap(NumberField& other) noexcept {
    data_.swap(other.data_);
}

void NumberField::clear() noexcept {
    data_.reset();
}

bool NumberField::set(const NumberField& other) noexcept {
    data_ = other.data_;
    return true;
}

NumberField NumberField::by_polynomial(
        flint::FmpqPolyConstRef polynomial) noexcept {
    NumberField field;
    if (!field.define_by_polynomial(polynomial)) {
        return {};
    }
    return field;
}

NumberField NumberField::by_polynomial(
        flint::FmpzPolyConstRef polynomial) noexcept {
    NumberField field;
    if (!field.define_by_polynomial(polynomial)) {
        return {};
    }
    return field;
}

NumberField NumberField::quadratic(flint::FmpzConstRef radicand) noexcept {
    NumberField field;
    if (!field.define_quadratic(radicand)) {
        return {};
    }
    return field;
}

bool NumberField::define_by_polynomial(flint::FmpqPolyConstRef polynomial) noexcept {
    if (!degree_is_positive(polynomial)) {
        return false;
    }

    auto next = make_field_data(polynomial);
    if (next == nullptr) {
        return false;
    }

    flint::Fmpz radicand;
    if (canonical_quadratic_radicand(flint::FmpzRef(radicand), polynomial)) {
        set_quadratic_backend(*next, flint::FmpzConstRef(radicand));
    }

    data_ = std::move(next);
    return true;
}

bool NumberField::define_by_polynomial(flint::FmpzPolyConstRef polynomial) noexcept {
    flint::FmpqPoly rational_polynomial;
    fmpq_poly_set_fmpz_poly(rational_polynomial.raw(), polynomial.raw());
    return define_by_polynomial(flint::FmpqPolyConstRef(rational_polynomial));
}

bool NumberField::define_quadratic(flint::FmpzConstRef radicand) noexcept {
    if (!is_valid_quadratic_radicand(radicand)) {
        return false;
    }

    flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_si(polynomial.raw(), 2, 1);
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, radicand.raw());
    fmpq_poly_neg(polynomial.raw(), polynomial.raw());
    fmpq_poly_set_coeff_si(polynomial.raw(), 2, 1);

    auto next = make_field_data(flint::FmpqPolyConstRef(polynomial));
    if (next == nullptr) {
        return false;
    }
    set_quadratic_backend(*next, radicand);

    data_ = std::move(next);
    return true;
}

bool NumberField::is_defined() const noexcept {
    return data_ != nullptr && data_->field.is_defined();
}

bool NumberField::has_same_data(const NumberField& other) const noexcept {
    return is_defined() && other.is_defined() && data_ == other.data_;
}

slong NumberField::degree() const noexcept {
    if (!is_defined()) {
        return 0;
    }
    return fmpq_poly_degree(data_->field.raw()->pol);
}

NumberFieldBackendKind NumberField::backend_kind() const noexcept {
    if (!is_defined()) {
        return NumberFieldBackendKind::generic;
    }
    return data_->backend_kind;
}

bool NumberField::quadratic_radicand(flint::FmpzRef out) const noexcept {
    if (!is_defined() || !data_->has_quadratic_radicand) {
        return false;
    }
    fmpz_set(out.raw(), data_->quadratic_radicand.raw());
    return true;
}

flint::NfRef NumberField::flint_field_ref() noexcept {
    return flint::NfRef(data_->field);
}

flint::NfConstRef NumberField::flint_field_ref() const noexcept {
    return flint::NfConstRef(data_->field);
}

nf_struct* NumberField::raw_flint_field() noexcept {
    return is_defined() ? data_->field.raw() : nullptr;
}

const nf_struct* NumberField::raw_flint_field() const noexcept {
    return is_defined() ? data_->field.raw() : nullptr;
}

}  // namespace silex
