#include <silex/element.hpp>

#include "element_internal.hpp"

#include <silex/diagnostics.hpp>

namespace silex {

namespace detail {

bool ensure_parent(Element& out, const NumberField& field) noexcept {
    if (out.has_parent(field)) {
        return true;
    }
    if (out.parent() != nullptr) {
        return false;
    }
    out = Element(field);
    return out.is_defined();
}

}  // namespace detail

Element::Element(const NumberField& parent) noexcept {
    define(parent);
}

Element::~Element() noexcept = default;

Element::Element(Element&& other) noexcept {
    swap(other);
}

Element& Element::operator=(Element&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void Element::swap(Element& other) noexcept {
    parent_.swap(other.parent_);
    value_.swap(other.value_);
}

void Element::clear() noexcept {
    value_.clear();
    parent_.clear();
}

bool Element::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    Element next;
    next.parent_ = parent;
    next.value_ = flint::NfElem(parent.raw_flint_field());
    if (!next.value_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool Element::set(const Element& other) noexcept {
    if (!has_same_parent(other)) {
        return false;
    }
    if (this == &other) {
        return true;
    }
    nf_elem_set(value_.raw(), other.value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::is_defined() const noexcept {
    return parent_.is_defined() && value_.is_defined();
}

const NumberField* Element::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

bool Element::has_parent(const NumberField& parent) const noexcept {
    return is_defined() && parent_.has_same_data(parent);
}

bool Element::has_same_parent(const Element& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.has_same_data(other.parent_);
}

bool Element::zero() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_zero(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::one() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_one(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::gen() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_gen(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::set_si(slong value) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_si(value_.raw(), value, parent_.raw_flint_field());
    return true;
}

bool Element::set_fmpz(flint::FmpzConstRef value) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_fmpz(value_.raw(), value.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::set_si_over_si(slong numerator, slong denominator) noexcept {
    if (!is_defined() || denominator == 0) {
        return false;
    }

    Element tmp(parent_);
    nf_elem_set_si(tmp.value_.raw(), numerator, parent_.raw_flint_field());
    nf_elem_scalar_div_si(tmp.value_.raw(), tmp.value_.raw(), denominator,
                          parent_.raw_flint_field());
    swap(tmp);
    return true;
}

bool Element::set_fmpq_poly(flint::FmpqPolyConstRef polynomial) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_fmpq_poly(value_.raw(), polynomial.raw(),
                          parent_.raw_flint_field());
    return true;
}

bool Element::get_fmpq_poly(flint::FmpqPolyRef polynomial) const noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(),
                          parent_.raw_flint_field());
    return true;
}

std::optional<flint::FmpqPoly> Element::to_fmpq_poly() const noexcept {
    flint::FmpqPoly polynomial;
    if (!get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return std::nullopt;
    }
    return polynomial;
}

bool Element::equal(const Element& other) const noexcept {
    if (!has_same_parent(other)) {
        return false;
    }
    return nf_elem_equal(value_.raw(), other.value_.raw(),
                         parent_.raw_flint_field()) != 0;
}

bool Element::equal_si(slong value) const noexcept {
    if (!is_defined()) {
        return false;
    }
    return nf_elem_equal_si(value_.raw(), value,
                            parent_.raw_flint_field()) != 0;
}

bool Element::negate(const Element& input) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    nf_elem_neg(value_.raw(), input.value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::add(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_add(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::add_si(const Element& input, slong value) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    nf_elem_add_si(value_.raw(), input.value_.raw(), value,
                   parent_.raw_flint_field());
    return true;
}

bool Element::subtract(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_sub(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::multiply(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_mul(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::scalar_div_si(const Element& input,
                            slong denominator) noexcept {
    if (!has_same_parent(input) || denominator == 0) {
        return false;
    }

    Element tmp(parent_);
    nf_elem_scalar_div_si(tmp.value_.raw(), input.value_.raw(), denominator,
                          parent_.raw_flint_field());
    swap(tmp);
    return true;
}

bool Element::invert(const Element& input) noexcept {
    if (!has_same_parent(input) || input.equal_si(0)) {
        return false;
    }

    // NumberField construction only accepts irreducible defining
    // polynomials, so every nonzero element is invertible.  A zero divisor is
    // still rejected here, exactly and per FLINT representation, rather than
    // letting nf_elem_inv return a non-inverse.  Linear fields are Q, where
    // nonzero already implies invertible.
    const nf_struct* field = parent_.raw_flint_field();
    Element tmp(parent_);
    if ((field->flag & NF_LINEAR) != 0) {
        nf_elem_inv(tmp.value_.raw(), input.value_.raw(), field);
    } else if ((field->flag & NF_QUADRATIC) != 0) {
        // A quadratic element is a zero divisor exactly when its norm, the
        // resultant with the defining polynomial, vanishes.  The quadratic
        // norm is a constant number of fmpz operations.
        flint::Fmpq norm;
        nf_elem_norm(norm.raw(), input.value_.raw(), field);
        if (fmpq_is_zero(norm.raw()) != 0) {
            return false;
        }
        nf_elem_inv(tmp.value_.raw(), input.value_.raw(), field);
    } else {
        // The generic branch of FLINT's _nf_elem_inv (nf_elem/inv.c) computes
        // the Bezout cofactor of fmpq_poly_xgcd and discards the monic gcd.
        // Perform the same call and keep the gcd, which is 1 exactly when the
        // element is invertible.
        flint::FmpqPoly gcd;
        flint::FmpqPoly cofactor;
        fmpq_poly_xgcd(gcd.raw(), NF_ELEM(tmp.value_.raw()), cofactor.raw(),
                       NF_ELEM(input.value_.raw()), field->pol);
        if (fmpq_poly_is_one(gcd.raw()) == 0) {
            return false;
        }
    }
    swap(tmp);
    return true;
}

bool Element::pow_fmpz(const Element& input,
                       flint::FmpzConstRef exponent) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    const int exponent_sign = fmpz_sgn(exponent.raw());
    if (exponent_sign < 0 && input.equal_si(0)) {
        return false;
    }

    const nf_struct* raw_field = parent_.raw_flint_field();
    if (flint::fmpz_abs_fits_ui(exponent)) {
        if (exponent_sign >= 0) {
            Element tmp(parent_);
            nf_elem_pow(tmp.value_.raw(), input.value_.raw(),
                        flint::fmpz_get_ui(exponent), raw_field);
            swap(tmp);
            return true;
        }

        Element inverse(parent_);
        Element tmp(parent_);
        flint::Fmpz abs_exponent;
        if (!inverse.invert(input)) {
            return false;
        }
        flint::fmpz_neg(flint::FmpzRef(abs_exponent), exponent);
        nf_elem_pow(tmp.value_.raw(), inverse.value_.raw(),
                    flint::fmpz_get_ui(flint::FmpzConstRef(abs_exponent)),
                    raw_field);
        swap(tmp);
        return true;
    }

    Element accumulator(parent_);
    Element base(parent_);
    Element product(parent_);
    flint::Fmpz k;

    accumulator.one();
    if (exponent_sign < 0) {
        if (!base.invert(input)) {
            return false;
        }
        fmpz_neg(k.raw(), exponent.raw());
    } else {
        if (!base.set(input)) {
            return false;
        }
        fmpz_set(k.raw(), exponent.raw());
    }

    while (fmpz_is_zero(k.raw()) == 0) {
        if (fmpz_is_odd(k.raw()) != 0) {
            if (!product.multiply(accumulator, base) ||
                !accumulator.set(product)) {
                return false;
            }
        }

        fmpz_fdiv_q_2exp(k.raw(), k.raw(), 1);
        if (fmpz_is_zero(k.raw()) == 0) {
            if (!product.multiply(base, base) || !base.set(product)) {
                return false;
            }
        }
    }

    swap(accumulator);
    return true;
}

bool Element::is_square(
        bool& is_square,
        Element& root,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.is_square");
    if (!has_same_parent(root)) {
        return false;
    }

    Element candidate(parent_);
    bool flag = false;
    bool known = detail::is_square_rational_constant(flag, candidate, *this);

    if (known && !flag && parent_.degree() != 1) {
        known = false;
    }

    if (!known && parent_.degree() != 1 &&
        parent_.backend_kind() == NumberFieldBackendKind::quadratic) {
        flint::Fmpz radicand;
        if (parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
            known = detail::is_square_quadratic(flag, candidate, *this, radicand.raw());
        }
    }

    if (!known &&
        detail::pure_square_hensel_root(flag, candidate, *this, diagnostics)) {
        known = true;
    }

    if (!known &&
        detail::pure_square_residue_disproves(flag, *this, diagnostics)) {
        known = true;
    }

    if (!known) {
        return false;
    }

    is_square = flag;
    if (flag) {
        root.swap(candidate);
    }
    return true;
}

bool Element::is_power(bool& is_power,
                       Element& root,
                       flint::FmpzConstRef exponent,
                       const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.is_power");
    if (!has_same_parent(root) || fmpz_cmp_ui(exponent.raw(), 1) < 0) {
        return false;
    }

    if (fmpz_is_one(exponent.raw()) != 0) {
        Element candidate(parent_);
        if (!candidate.set(*this)) {
            return false;
        }
        is_power = true;
        root.swap(candidate);
        return true;
    }

    Element candidate(parent_);
    bool flag = false;
    if (detail::is_power_rational_constant(flag, candidate, *this, exponent)) {
        if (flag) {
            is_power = true;
            root.swap(candidate);
            return true;
        }
        if (parent_.degree() == 1) {
            is_power = false;
            return true;
        }
    }

    if (fmpz_equal_ui(exponent.raw(), 2) != 0) {
        return is_square(is_power, root, diagnostics);
    }

    if (fmpz_fits_si(exponent.raw()) != 0) {
        const slong exponent_si = fmpz_get_si(exponent.raw());
        if (exponent_si > 1) {
            if (detail::pure_power_hensel_root(
                        flag, candidate, *this, exponent_si, diagnostics)) {
                is_power = flag;
                if (flag) {
                    root.swap(candidate);
                }
                return true;
            }
            if (detail::pure_power_residue_disproves(
                        flag, *this, exponent_si, diagnostics)) {
                is_power = flag;
                return true;
            }
        }
    }

    return false;
}

bool Element::trace(flint::FmpqRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    if (parent_.backend_kind() == NumberFieldBackendKind::quadratic) {
        flint::FmpqPoly polynomial;
        flint::Fmpq constant;

        nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
        fmpq_poly_get_coeff_fmpq(constant.raw(), polynomial.raw(), 0);
        fmpq_mul_2exp(out.raw(), constant.raw(), 1);
        return true;
    }

    detail::generic_trace(out, field, parent_.degree(), value_.raw());
    return true;
}

bool Element::norm(flint::FmpqRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    flint::Fmpz radicand;
    if (parent_.backend_kind() == NumberFieldBackendKind::quadratic &&
        parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
        flint::FmpqPoly polynomial;
        flint::Fmpq constant;
        flint::Fmpq linear;
        flint::Fmpq term;

        nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
        fmpq_poly_get_coeff_fmpq(constant.raw(), polynomial.raw(), 0);
        fmpq_poly_get_coeff_fmpq(linear.raw(), polynomial.raw(), 1);

        fmpq_mul(out.raw(), constant.raw(), constant.raw());
        fmpq_mul(term.raw(), linear.raw(), linear.raw());
        fmpq_mul_fmpz(term.raw(), term.raw(), radicand.raw());
        fmpq_sub(out.raw(), out.raw(), term.raw());
        return true;
    }

    detail::generic_norm(out, field, parent_.degree(), value_.raw());
    return true;
}

bool Element::conjugate(Element& out) const noexcept {
    if (!has_same_parent(out) ||
        parent_.backend_kind() != NumberFieldBackendKind::quadratic) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    flint::FmpqPoly polynomial;
    flint::Fmpq linear;

    nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
    fmpq_poly_get_coeff_fmpq(linear.raw(), polynomial.raw(), 1);
    fmpq_neg(linear.raw(), linear.raw());
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, linear.raw());
    nf_elem_set_fmpq_poly(out.value_.raw(), polynomial.raw(), field);
    return true;
}

flint::NfElemRef Element::flint_element_ref() noexcept {
    return flint::NfElemRef(value_);
}

flint::NfElemConstRef Element::flint_element_ref() const noexcept {
    return flint::NfElemConstRef(value_);
}

nf_elem_struct* Element::raw_flint_element() noexcept {
    return value_.raw();
}

const nf_elem_struct* Element::raw_flint_element() const noexcept {
    return value_.raw();
}

}  // namespace silex
