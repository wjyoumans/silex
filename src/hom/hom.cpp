#include <silex/hom.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nf_elem.h>

#include <silex/flint/fmpq_mat.hpp>

#include <utility>

namespace silex {
namespace {

bool evaluate_polynomial_at(
        Element& out,
        flint::FmpqPolyConstRef polynomial,
        const Element& beta) noexcept {
    if (out.parent() == nullptr || beta.parent() == nullptr ||
        !out.parent()->has_same_data(*beta.parent())) {
        return false;
    }

    const nf_struct* field = beta.parent()->raw_flint_field();
    if (field == nullptr) {
        return false;
    }

    flint::NfElem accumulator(field);
    flint::NfElem product(field);
    flint::Fmpq coefficient;
    nf_elem_zero(accumulator.raw(), field);

    for (slong i = fmpq_poly_length(polynomial.raw()) - 1; i >= 0; --i) {
        nf_elem_mul(product.raw(), accumulator.raw(),
                    beta.raw_flint_element(), field);
        fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial.raw(), i);
        nf_elem_add_fmpq(accumulator.raw(), product.raw(),
                         coefficient.raw(), field);
    }

    nf_elem_set(out.raw_flint_element(), accumulator.raw(), field);
    return true;
}

bool fmpq_mat_entries_are_integral(flint::FmpqMatConstRef matrix) noexcept {
    for (slong i = 0; i < fmpq_mat_nrows(matrix.raw()); ++i) {
        for (slong j = 0; j < fmpq_mat_ncols(matrix.raw()); ++j) {
            if (fmpz_is_one(fmpq_mat_entry_den(matrix.raw(), i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

bool output_shape(flint::FmpzMatRef matrix,
                  slong rows,
                  slong cols) noexcept {
    return fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

bool set_element_from_power_row(Element& out,
                                flint::FmpqMatConstRef power_row) noexcept {
    if (out.parent() == nullptr ||
        fmpq_mat_nrows(power_row.raw()) != 1 ||
        fmpq_mat_ncols(power_row.raw()) != out.parent()->degree()) {
        return false;
    }

    flint::FmpqPoly polynomial;
    for (slong j = 0; j < fmpq_mat_ncols(power_row.raw()); ++j) {
        fmpq_poly_set_coeff_fmpq(polynomial.raw(), j,
                                 fmpq_mat_entry(power_row.raw(), 0, j));
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

void set_fmpz_row_from_integral_fmpq_row(flint::FmpzMat& matrix,
                                         slong row,
                                         flint::FmpqMatConstRef coords) noexcept {
    for (slong j = 0; j < fmpq_mat_ncols(coords.raw()); ++j) {
        fmpz_set(fmpz_mat_entry(matrix.raw(), row, j),
                 fmpq_mat_entry_num(coords.raw(), 0, j));
    }
}

bool set_order_image_matrix(flint::FmpzMat& out,
                            const Order& source,
                            const Order& target,
                            const FieldHom& field_hom) noexcept {
    if (!source.has_basis() || !target.has_basis() ||
        source.parent() == nullptr || target.parent() == nullptr ||
        field_hom.domain() == nullptr || field_hom.codomain() == nullptr ||
        !source.parent()->has_same_data(*field_hom.domain()) ||
        !target.parent()->has_same_data(*field_hom.codomain())) {
        return false;
    }

    const slong source_degree = source.degree();
    const slong target_degree = target.degree();
    flint::FmpzMat matrix(source_degree, target_degree);
    flint::FmpqMat source_basis(source_degree, source_degree);
    flint::FmpqMat basis_row(1, source_degree);
    flint::FmpqMat coords(1, target_degree);
    Element x(*source.parent());
    Element y(*target.parent());

    if (!source.get_basis(flint::FmpqMatRef(source_basis))) {
        return false;
    }

    for (slong i = 0; i < source_degree; ++i) {
        for (slong j = 0; j < source_degree; ++j) {
            fmpq_set(fmpq_mat_entry(basis_row.raw(), 0, j),
                     fmpq_mat_entry(source_basis.raw(), i, j));
        }
        if (!set_element_from_power_row(
                    x, flint::FmpqMatConstRef(basis_row)) ||
            !field_hom.apply(y, x) ||
            !target.coordinates(flint::FmpqMatRef(coords), y) ||
            !fmpq_mat_entries_are_integral(flint::FmpqMatConstRef(coords))) {
            return false;
        }
        set_fmpz_row_from_integral_fmpq_row(
                matrix, i, flint::FmpqMatConstRef(coords));
    }

    out = std::move(matrix);
    return true;
}

}  // namespace

FieldHom::FieldHom(const NumberField& domain,
                   const NumberField& codomain) noexcept {
    define(domain, codomain);
}

FieldHom::~FieldHom() noexcept = default;

FieldHom::FieldHom(FieldHom&& other) noexcept {
    swap(other);
}

FieldHom& FieldHom::operator=(FieldHom&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FieldHom::swap(FieldHom& other) noexcept {
    domain_.swap(other.domain_);
    codomain_.swap(other.codomain_);
    generator_image_.swap(other.generator_image_);
    std::swap(has_generator_image_, other.has_generator_image_);
}

void FieldHom::clear() noexcept {
    generator_image_.clear();
    domain_.clear();
    codomain_.clear();
    has_generator_image_ = false;
}

bool FieldHom::define(const NumberField& domain,
                      const NumberField& codomain) noexcept {
    if (!domain.is_defined() || !codomain.is_defined()) {
        return false;
    }

    FieldHom next;
    next.domain_ = domain;
    next.codomain_ = codomain;
    next.generator_image_ = Element(codomain);
    if (!next.generator_image_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool FieldHom::set(const FieldHom& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    FieldHom copy(other.domain_, other.codomain_);
    if (other.has_generator_image_) {
        if (!copy.generator_image_.set(other.generator_image_)) {
            return false;
        }
        copy.has_generator_image_ = true;
    }

    swap(copy);
    return true;
}

bool FieldHom::is_defined() const noexcept {
    return domain_.is_defined() && codomain_.is_defined() &&
           generator_image_.is_defined();
}

const NumberField* FieldHom::domain() const noexcept {
    return is_defined() ? &domain_ : nullptr;
}

const NumberField* FieldHom::codomain() const noexcept {
    return is_defined() ? &codomain_ : nullptr;
}

bool FieldHom::has_generator_image() const noexcept {
    return is_defined() && has_generator_image_;
}

bool FieldHom::set_generator_image(const Element& image) noexcept {
    if (!is_defined() || !image.has_parent(codomain_)) {
        return false;
    }

    const nf_struct* source = domain_.raw_flint_field();
    if (source == nullptr) {
        return false;
    }

    Element value(codomain_);
    if (!evaluate_polynomial_at(
                value, flint::FmpqPolyConstRef(source->pol), image) ||
        !value.equal_si(0)) {
        return false;
    }

    Element candidate(codomain_);
    if (!candidate.set(image)) {
        return false;
    }

    generator_image_.swap(candidate);
    has_generator_image_ = true;
    return true;
}

bool FieldHom::generator_image(Element& out) const noexcept {
    if (!is_defined() || !has_generator_image_ || !out.has_parent(codomain_)) {
        return false;
    }
    return out.set(generator_image_);
}

bool FieldHom::apply(Element& out, const Element& input) const noexcept {
    if (!is_defined() || !has_generator_image_ ||
        !input.has_parent(domain_) || !out.has_parent(codomain_)) {
        return false;
    }

    flint::FmpqPoly polynomial;
    Element candidate(codomain_);
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) ||
        !evaluate_polynomial_at(
                candidate, flint::FmpqPolyConstRef(polynomial),
                generator_image_)) {
        return false;
    }

    return out.set(candidate);
}

bool FieldHom::is_identity() const noexcept {
    if (!is_defined() || !has_generator_image_ ||
        !domain_.has_same_data(codomain_)) {
        return false;
    }

    Element generator(codomain_);
    return generator.gen() && generator_image_.equal(generator);
}

bool FieldHom::is_isomorphism() const noexcept {
    return is_defined() && has_generator_image_ &&
           domain_.degree() == codomain_.degree();
}

OrderHom::OrderHom(const Order& source, const Order& target) noexcept {
    define(source, target);
}

OrderHom::~OrderHom() noexcept = default;

OrderHom::OrderHom(OrderHom&& other) noexcept {
    swap(other);
}

OrderHom& OrderHom::operator=(OrderHom&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void OrderHom::swap(OrderHom& other) noexcept {
    source_.swap(other.source_);
    target_.swap(other.target_);
    field_hom_.swap(other.field_hom_);
    image_matrix_.swap(other.image_matrix_);
    std::swap(has_field_hom_, other.has_field_hom_);
}

void OrderHom::clear() noexcept {
    field_hom_.clear();
    image_matrix_ = flint::FmpzMat(0, 0);
    source_.clear();
    target_.clear();
    has_field_hom_ = false;
}

bool OrderHom::define(const Order& source, const Order& target) noexcept {
    if (!source.is_defined() || !source.has_basis() ||
        !target.is_defined() || !target.has_basis() ||
        source.parent() == nullptr || target.parent() == nullptr) {
        return false;
    }

    OrderHom next;
    next.source_ = source;
    next.target_ = target;
    next.field_hom_ = FieldHom(*source.parent(), *target.parent());
    if (!next.field_hom_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool OrderHom::set(const OrderHom& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    OrderHom copy(other.source_, other.target_);
    if (other.has_field_hom_) {
        if (!copy.field_hom_.set(other.field_hom_)) {
            return false;
        }
        copy.image_matrix_ =
                flint::FmpzMat(
                        fmpz_mat_nrows(other.image_matrix_.raw()),
                        fmpz_mat_ncols(other.image_matrix_.raw()));
        fmpz_mat_set(copy.image_matrix_.raw(), other.image_matrix_.raw());
        copy.has_field_hom_ = true;
    }

    swap(copy);
    return true;
}

bool OrderHom::is_defined() const noexcept {
    return source_.is_defined() && target_.is_defined() &&
           field_hom_.is_defined();
}

const Order* OrderHom::source_order() const noexcept {
    return is_defined() ? &source_ : nullptr;
}

const Order* OrderHom::target_order() const noexcept {
    return is_defined() ? &target_ : nullptr;
}

bool OrderHom::has_field_homomorphism() const noexcept {
    return is_defined() && has_field_hom_;
}

bool OrderHom::set_field_homomorphism(const FieldHom& field_hom) noexcept {
    if (!is_defined() || !field_hom.has_generator_image() ||
        field_hom.domain() == nullptr || source_.parent() == nullptr ||
        field_hom.codomain() == nullptr || target_.parent() == nullptr ||
        !field_hom.domain()->has_same_data(*source_.parent()) ||
        !field_hom.codomain()->has_same_data(*target_.parent())) {
        return false;
    }

    OrderHom candidate(source_, target_);
    if (!set_order_image_matrix(candidate.image_matrix_,
                                source_,
                                target_,
                                field_hom) ||
        !candidate.field_hom_.set(field_hom)) {
        return false;
    }

    candidate.has_field_hom_ = true;
    swap(candidate);
    return true;
}

bool OrderHom::apply(Element& out, const Element& input) const noexcept {
    return is_defined() && has_field_hom_ && field_hom_.apply(out, input);
}

bool OrderHom::field_homomorphism(FieldHom& out) const noexcept {
    if (!is_defined() || !has_field_hom_) {
        return false;
    }
    return out.set(field_hom_);
}

bool OrderHom::image_matrix(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !has_field_hom_ ||
        !output_shape(out,
                      fmpz_mat_nrows(image_matrix_.raw()),
                      fmpz_mat_ncols(image_matrix_.raw()))) {
        return false;
    }
    fmpz_mat_set(out.raw(), image_matrix_.raw());
    return true;
}

bool OrderHom::image_matrix(flint::FmpzMat& out) const noexcept {
    if (!is_defined() || !has_field_hom_) {
        return false;
    }

    flint::FmpzMat copy(fmpz_mat_nrows(image_matrix_.raw()),
                        fmpz_mat_ncols(image_matrix_.raw()));
    fmpz_mat_set(copy.raw(), image_matrix_.raw());
    out = std::move(copy);
    return true;
}

}  // namespace silex
