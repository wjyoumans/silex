#include <silex/element.hpp>

#include "element_internal.hpp"

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include <flint/fmpz_vec.h>

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

namespace {

// Private shape test for Element::norm: true when the FLINT field polynomial
// is x^n + c with integer c != 0 and n even, n >= 4.  Degrees 1 and 2 use
// FLINT's linear and quadratic element representations and keep their
// existing paths.
bool is_even_binomial_field(const nf_struct* field) noexcept {
    if ((field->flag & (NF_LINEAR | NF_QUADRATIC)) != 0) {
        return false;
    }
    const fmpq_poly_struct* polynomial = field->pol;
    const slong degree = polynomial->length - 1;
    if (degree < 4 || degree % 2 != 0 ||
        fmpz_is_one(polynomial->den) == 0 ||
        fmpz_is_one(polynomial->coeffs + degree) == 0 ||
        fmpz_is_zero(polynomial->coeffs) != 0) {
        return false;
    }
    for (slong i = degree - 1; i > 0; --i) {
        if (fmpz_is_zero(polynomial->coeffs + i) == 0) {
            return false;
        }
    }
    return true;
}

// Exact norm of `element` in a field defined by f = x^n + c with n even, as
// checked by is_even_binomial_field.
//
// FLINT 3.6 `src/nf_elem/norm.c` (`_nf_elem_norm`, generic branch) computes
// the norm of a = Z/d (Z integral, d its denominator) as the resultant
// Res(f, a); for monic integral f this is Res(f, Z) / d^n, the product of
// Z(alpha) over the roots alpha of f.  Write f = h(x^2) with
// h = y^(n/2) + c, and Z(x) = E(x^2) + x O(x^2).  The roots of f are the two
// square roots +-sqrt(beta) of each root beta of h, so
//
//     Res(f, Z) = prod_beta Z(sqrt(beta)) Z(-sqrt(beta))
//               = prod_beta G(beta) = Res(h, G),
//     G(y) = E(y)^2 - y O(y)^2,
//
// where G is the Graeffe transform of Z, G(x^2) = Z(x) Z(-x) (source anchor
// in docs/reference/algorithms_and_sources.rst, "Number fields, orders, and
// ideals").
// Since h is monic, G may be replaced by its remainder modulo h, which
// folds y^k to -c y^(k - n/2).  This is the norm transitivity
// N_{K/Q} = N_{Q(alpha^2)/Q} o N_{K/Q(alpha^2)} for alpha a root of f.  The
// step repeats while the degree is even (12 -> 6 -> 3); the remaining odd
// degree uses FLINT's integer resultant.  Every step is exact integer
// arithmetic, and fmpq_set_fmpz_frac gives the same canonical rational as
// nf_elem_norm.
void even_binomial_norm(flint::FmpqRef out,
                        const nf_struct* field,
                        const nf_elem_t element) noexcept {
    const fmpq_poly_struct* value = NF_ELEM(element);
    const slong length = value->length;
    if (length == 0) {
        fmpq_zero(out.raw());
        return;
    }

    const fmpz* constant = field->pol->coeffs;
    const slong field_degree = field->pol->length - 1;
    slong degree = field_degree;

    flint::FmpzPoly current;
    flint::FmpzPoly even;
    flint::FmpzPoly odd;
    fmpz_poly_fit_length(current.raw(), length);
    _fmpz_vec_set(current.raw()->coeffs, value->coeffs, length);
    _fmpz_poly_set_length(current.raw(), length);

    while (degree % 2 == 0) {
        const slong half = degree / 2;
        const slong current_length = current.raw()->length;
        const slong even_length = (current_length + 1) / 2;
        const slong odd_length = current_length / 2;

        fmpz_poly_fit_length(even.raw(), even_length);
        fmpz_poly_fit_length(odd.raw(), odd_length);
        for (slong k = 0; k < even_length; ++k) {
            fmpz_set(even.raw()->coeffs + k, current.raw()->coeffs + 2 * k);
        }
        for (slong k = 0; k < odd_length; ++k) {
            fmpz_set(odd.raw()->coeffs + k,
                     current.raw()->coeffs + 2 * k + 1);
        }
        _fmpz_poly_set_length(even.raw(), even_length);
        _fmpz_poly_set_length(odd.raw(), odd_length);
        _fmpz_poly_normalise(even.raw());
        _fmpz_poly_normalise(odd.raw());

        fmpz_poly_sqr(even.raw(), even.raw());
        fmpz_poly_sqr(odd.raw(), odd.raw());
        fmpz_poly_shift_left(odd.raw(), odd.raw(), 1);
        fmpz_poly_sub(current.raw(), even.raw(), odd.raw());

        fmpz* coeffs = current.raw()->coeffs;
        const slong product_length = current.raw()->length;
        for (slong k = product_length - 1; k >= half; --k) {
            if (fmpz_is_zero(coeffs + k) == 0) {
                fmpz_submul(coeffs + (k - half), coeffs + k, constant);
                fmpz_zero(coeffs + k);
            }
        }
        if (product_length > half) {
            _fmpz_poly_set_length(current.raw(), half);
        }
        _fmpz_poly_normalise(current.raw());
        degree = half;
    }

    flint::FmpzPoly modulus;
    fmpz_poly_set_coeff_ui(modulus.raw(), static_cast<ulong>(degree), 1);
    fmpz_poly_set_coeff_fmpz(modulus.raw(), 0, constant);

    flint::Fmpz resultant;
    flint::Fmpz denominator;
    fmpz_poly_resultant(resultant.raw(), modulus.raw(), current.raw());
    fmpz_pow_ui(denominator.raw(), value->den,
                static_cast<ulong>(field_degree));
    fmpq_set_fmpz_frac(out.raw(), resultant.raw(), denominator.raw());
}

}  // namespace

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

    if (is_even_binomial_field(field)) {
        even_binomial_norm(out, field, value_.raw());
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
