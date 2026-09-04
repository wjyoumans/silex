#include <silex/aut.hpp>

#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>

#include <flint/nf.h>

#include <utility>

namespace silex {
namespace {

bool quadratic_conjugate_image(Element& out,
                               const NumberField& parent) noexcept {
    if (parent.degree() != 2) {
        return false;
    }

    const nf_struct* field = parent.raw_flint_field();
    if (field == nullptr) {
        return false;
    }

    flint::Fmpq leading;
    flint::Fmpq linear;
    flint::Fmpq constant;
    flint::fmpq_poly_get_coeff_fmpq(
            flint::FmpqRef(leading),
            flint::FmpqPolyConstRef(field->pol),
            2);
    flint::fmpq_poly_get_coeff_fmpq(
            flint::FmpqRef(linear),
            flint::FmpqPolyConstRef(field->pol),
            1);
    flint::fmpq_div(constant, linear, leading);
    flint::fmpq_neg(constant, constant);

    flint::FmpqPoly image;
    flint::fmpq_poly_set_coeff_fmpq(image, 0, constant);
    flint::fmpq_poly_set_coeff_si(image, 1, -1);
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(image));
}

}  // namespace

FieldAutomorphism::FieldAutomorphism(const NumberField& parent) noexcept {
    define(parent);
}

FieldAutomorphism::FieldAutomorphism(FieldAutomorphism&& other) noexcept {
    swap(other);
}

FieldAutomorphism& FieldAutomorphism::operator=(
        FieldAutomorphism&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FieldAutomorphism::swap(FieldAutomorphism& other) noexcept {
    parent_.swap(other.parent_);
    hom_.swap(other.hom_);
    std::swap(has_homomorphism_, other.has_homomorphism_);
}

void FieldAutomorphism::clear() noexcept {
    parent_.clear();
    hom_.clear();
    has_homomorphism_ = false;
}

bool FieldAutomorphism::define(const NumberField& parent) noexcept {
    FieldHom next_hom(parent, parent);
    if (!next_hom.is_defined()) {
        return false;
    }

    parent_ = parent;
    hom_.swap(next_hom);
    has_homomorphism_ = false;
    return true;
}

bool FieldAutomorphism::set(const FieldAutomorphism& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    FieldAutomorphism copy(other.parent_);
    if (other.has_homomorphism_) {
        if (!copy.hom_.set(other.hom_)) {
            return false;
        }
        copy.has_homomorphism_ = true;
    }

    swap(copy);
    return true;
}

bool FieldAutomorphism::is_defined() const noexcept {
    return parent_.is_defined() && hom_.is_defined();
}

const NumberField* FieldAutomorphism::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

bool FieldAutomorphism::has_homomorphism() const noexcept {
    return has_homomorphism_;
}

bool FieldAutomorphism::set_identity() noexcept {
    if (!is_defined()) {
        return false;
    }

    FieldAutomorphism candidate(parent_);
    Element generator(parent_);
    if (!candidate.is_defined() || !generator.gen() ||
        !candidate.hom_.set_generator_image(generator)) {
        return false;
    }

    candidate.has_homomorphism_ = true;
    swap(candidate);
    return true;
}

bool FieldAutomorphism::set_quadratic_conjugation() noexcept {
    if (!is_defined() || parent_.degree() != 2) {
        return false;
    }

    FieldAutomorphism candidate(parent_);
    Element image(parent_);
    if (!candidate.is_defined() ||
        !quadratic_conjugate_image(image, parent_) ||
        !candidate.hom_.set_generator_image(image)) {
        return false;
    }

    candidate.has_homomorphism_ = true;
    swap(candidate);
    return true;
}

bool FieldAutomorphism::apply(Element& out,
                              const Element& input) const noexcept {
    return has_homomorphism_ && hom_.apply(out, input);
}

bool FieldAutomorphism::homomorphism(FieldHom& out) const noexcept {
    if (!has_homomorphism_) {
        return false;
    }
    return out.set(hom_);
}

bool FieldAutomorphism::is_identity() const noexcept {
    return has_homomorphism_ && hom_.is_identity();
}

}  // namespace silex
