#include <silex/residue_field.hpp>

#include "residue_field_internal.hpp"

#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpq_mat.hpp>

#include <flint/fmpz_factor.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>

#include <utility>

namespace silex {
namespace {

bool reduce_polynomial(flint::FmpzModPoly& out,
                       flint::FmpzPolyConstRef polynomial,
                       const flint::FmpzModCtx& ctx,
                       const flint::FmpzModPoly& modulus) noexcept {
    flint::FmpzModPoly input(ctx);
    fmpz_mod_poly_set_fmpz_poly(input.raw(), polynomial.raw(), ctx.raw());
    fmpz_mod_poly_rem(out.raw(), input.raw(), modulus.raw(), ctx.raw());
    return true;
}

void row_denominator_lcm(flint::Fmpz& out,
                         const flint::FmpqMat& row) noexcept {
    fmpz_one(out.raw());
    for (slong j = 0; j < fmpq_mat_ncols(row.raw()); ++j) {
        fmpz_lcm(out.raw(), out.raw(),
                 fmpq_mat_entry_den(row.raw(), 0, j));
    }
}

void row_scale_to_polynomial(flint::FmpzPoly& out,
                             const flint::FmpqMat& row,
                             const flint::Fmpz& denominator) noexcept {
    flint::Fmpz quotient;
    flint::Fmpz coefficient;

    fmpz_poly_zero(out.raw());
    for (slong j = 0; j < fmpq_mat_ncols(row.raw()); ++j) {
        fmpz_divexact(quotient.raw(), denominator.raw(),
                      fmpq_mat_entry_den(row.raw(), 0, j));
        fmpz_mul(coefficient.raw(), fmpq_mat_entry_num(row.raw(), 0, j),
                 quotient.raw());
        fmpz_poly_set_coeff_fmpz(out.raw(), j, coefficient.raw());
    }
}

bool set_residue_field_element_from_index(
        ResidueFieldElement& out,
        const ResidueField& parent,
        flint::FmpzConstRef index) noexcept {
    if (!parent.is_defined() || !out.is_defined() ||
        out.parent() == nullptr || !out.parent()->equal(parent)) {
        return false;
    }

    flint::Fmpz p;
    flint::Fmpz t;
    flint::Fmpz quotient;
    flint::Fmpz remainder;
    flint::FmpzPoly polynomial;
    const slong degree = parent.degree();
    if (degree <= 0 || !parent.characteristic(flint::FmpzRef(p))) {
        return false;
    }

    fmpz_set(t.raw(), index.raw());
    fmpz_poly_zero(polynomial.raw());
    for (slong i = 0; i < degree; ++i) {
        fmpz_fdiv_qr(quotient.raw(), remainder.raw(), t.raw(), p.raw());
        fmpz_poly_set_coeff_fmpz(polynomial.raw(), i, remainder.raw());
        fmpz_swap(t.raw(), quotient.raw());
    }

    return out.set_polynomial(flint::FmpzPolyConstRef(polynomial));
}

bool quotient_log_mod_prime_apply_ui(flint::FmpzRef out,
                                     flint::FmpzConstRef value,
                                     flint::FmpzConstRef cofactor,
                                     flint::FmpzConstRef quotient_generator,
                                     flint::FmpzConstRef p,
                                     flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_abs_fits_ui(value) ||
        !flint::fmpz_abs_fits_ui(cofactor) ||
        !flint::fmpz_abs_fits_ui(quotient_generator) ||
        !flint::fmpz_abs_fits_ui(p) ||
        !flint::fmpz_abs_fits_ui(ell)) {
        return false;
    }

    const ulong p_ui = flint::fmpz_get_ui(p);
    const ulong ell_ui = flint::fmpz_get_ui(ell);
    if (p_ui < 2 || ell_ui < 2) {
        return false;
    }

    const ulong value_ui = flint::fmpz_get_ui(value) % p_ui;
    const ulong quotient_generator_ui =
            flint::fmpz_get_ui(quotient_generator) % p_ui;
    if (value_ui == 0 || quotient_generator_ui == 0) {
        return false;
    }

    const ulong p_inverse = n_preinvert_limb(p_ui);
    const ulong image = n_powmod2_ui_preinv(
            value_ui, flint::fmpz_get_ui(cofactor), p_ui, p_inverse);
    nmod_t mod;
    nmod_init(&mod, p_ui);
    ulong current = 1 % p_ui;
    for (ulong i = 0; i < ell_ui; ++i) {
        if (current == image) {
            flint::fmpz_set_ui(out, i);
            return true;
        }
        current = nmod_mul(current, quotient_generator_ui, mod);
    }

    return false;
}

}  // namespace

ResidueField::ResidueField(flint::FmpzConstRef characteristic) noexcept
    : ctx_(characteristic.raw()),
      modulus_(ctx_) {
}

ResidueField::ResidueField(const PrimeIdeal& prime) noexcept {
    set_prime(prime);
}

ResidueField::~ResidueField() noexcept = default;

ResidueField::ResidueField(ResidueField&& other) noexcept {
    swap(other);
}

ResidueField& ResidueField::operator=(ResidueField&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueField::swap(ResidueField& other) noexcept {
    prime_.swap(other.prime_);
    p_.swap(other.p_);
    ctx_.swap(other.ctx_);
    modulus_.swap(other.modulus_);
    std::swap(defined_, other.defined_);
}

void ResidueField::clear() noexcept {
    ResidueField empty;
    swap(empty);
}

bool ResidueField::set(const ResidueField& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }
    return set_prime(other.prime_);
}

bool ResidueField::set_prime(const PrimeIdeal& prime) noexcept {
    if (!prime.has_prime_data() || prime.parent() == nullptr) {
        return false;
    }

    flint::Fmpz p;
    flint::FmpzPoly modulus_z;
    if (!prime.rational_prime(flint::FmpzRef(p)) ||
        !prime.residue_polynomial(flint::FmpzPolyRef(modulus_z)) ||
        fmpz_is_prime(p.raw()) == 0) {
        return false;
    }

    ResidueField next{flint::FmpzConstRef(p)};
    next.prime_ = PrimeIdeal(*prime.parent());
    if (!next.prime_.is_defined() || !next.prime_.set(prime)) {
        return false;
    }

    fmpz_set(next.p_.raw(), p.raw());
    fmpz_mod_poly_set_fmpz_poly(next.modulus_.raw(), modulus_z.raw(),
                                next.ctx_.raw());
    if (fmpz_mod_poly_degree(next.modulus_.raw(), next.ctx_.raw()) <= 0 ||
        fmpz_mod_poly_is_monic(next.modulus_.raw(), next.ctx_.raw()) == 0 ||
        fmpz_mod_poly_is_irreducible(next.modulus_.raw(),
                                     next.ctx_.raw()) == 0) {
        return false;
    }

    next.defined_ = true;
    swap(next);
    return true;
}

bool ResidueField::is_defined() const noexcept {
    return defined_ && prime_.has_prime_data() && modulus_.is_initialized();
}

const Order* ResidueField::parent_order() const noexcept {
    return is_defined() ? prime_.parent() : nullptr;
}

const PrimeIdeal* ResidueField::prime() const noexcept {
    return is_defined() ? &prime_ : nullptr;
}

slong ResidueField::degree() const noexcept {
    return is_defined() ? fmpz_mod_poly_degree(modulus_.raw(), ctx_.raw()) : 0;
}

bool ResidueField::get_prime(PrimeIdeal& out) const noexcept {
    return is_defined() && same_order_parent(out.parent(), prime_.parent()) &&
           out.set(prime_);
}

bool ResidueField::characteristic(flint::FmpzRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_set(out.raw(), p_.raw());
    return true;
}

bool ResidueField::cardinality(flint::FmpzRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_pow_ui(out.raw(), p_.raw(), static_cast<ulong>(degree()));
    return true;
}

bool ResidueField::modulus(flint::FmpzPolyRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_mod_poly_get_fmpz_poly(out.raw(), modulus_.raw(), ctx_.raw());
    return true;
}

std::optional<flint::Fmpz> ResidueField::characteristic() const noexcept {
    flint::Fmpz out;
    if (!characteristic(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

std::optional<flint::Fmpz> ResidueField::cardinality() const noexcept {
    flint::Fmpz out;
    if (!cardinality(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

std::optional<flint::FmpzPoly> ResidueField::modulus() const noexcept {
    flint::FmpzPoly out;
    if (!modulus(flint::FmpzPolyRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ResidueField::equal(const ResidueField& other) const noexcept {
    return is_defined() && other.is_defined() && prime_.equal(other.prime_);
}

bool ResidueField::multiplicative_generator(
        ResidueFieldElement& out) const noexcept {
    if (!is_defined() ||
        (out.is_defined() &&
         (out.parent() == nullptr || !out.parent()->equal(*this)))) {
        return false;
    }

    ResidueFieldElement candidate(*this);
    flint::Fmpz q;
    flint::Fmpz target;
    flint::Fmpz order;
    flint::Fmpz k;

    cardinality(flint::FmpzRef(q));
    fmpz_sub_ui(target.raw(), q.raw(), 1);

    for (fmpz_one(k.raw()); fmpz_cmp(k.raw(), q.raw()) < 0;
         fmpz_add_ui(k.raw(), k.raw(), 1)) {
        if (set_residue_field_element_from_index(
                    candidate, *this, flint::FmpzConstRef(k)) &&
            candidate.multiplicative_order(flint::FmpzRef(order)) &&
            fmpz_equal(order.raw(), target.raw()) != 0) {
            if (!out.is_defined()) {
                out = ResidueFieldElement(*this);
                if (!out.is_defined()) {
                    return false;
                }
            }
            return out.set(candidate);
        }
    }

    return false;
}

ResidueFieldElement::ResidueFieldElement(const ResidueField& parent) noexcept {
    define(parent);
}

ResidueFieldElement::~ResidueFieldElement() noexcept = default;

ResidueFieldElement::ResidueFieldElement(ResidueFieldElement&& other) noexcept {
    swap(other);
}

ResidueFieldElement& ResidueFieldElement::operator=(
        ResidueFieldElement&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueFieldElement::swap(ResidueFieldElement& other) noexcept {
    parent_.swap(other.parent_);
    rep_.swap(other.rep_);
}

void ResidueFieldElement::clear() noexcept {
    rep_ = flint::FmpzModPoly();
    parent_.clear();
}

bool ResidueFieldElement::define(const ResidueField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    ResidueFieldElement candidate;
    if (!candidate.parent_.set(parent)) {
        return false;
    }
    flint::FmpzModPoly next(candidate.parent_.ctx_);
    if (!next.is_initialized()) {
        return false;
    }
    candidate.rep_ = std::move(next);

    swap(candidate);
    return true;
}

bool ResidueFieldElement::set(const ResidueFieldElement& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!is_defined() || !other.is_defined() ||
        !parent_.equal(other.parent_)) {
        return false;
    }
    fmpz_mod_poly_set(rep_.raw(), other.rep_.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::is_defined() const noexcept {
    return parent_.is_defined() && rep_.is_initialized();
}

const ResidueField* ResidueFieldElement::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

bool ResidueFieldElement::zero() noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_mod_poly_zero(rep_.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::one() noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_mod_poly_one(rep_.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::is_zero() const noexcept {
    return is_defined() &&
           fmpz_mod_poly_is_zero(rep_.raw(), parent_.ctx_.raw()) != 0;
}

bool ResidueFieldElement::set_polynomial(
        flint::FmpzPolyConstRef polynomial) noexcept {
    if (!is_defined()) {
        return false;
    }

    flint::FmpzModPoly reduced(parent_.ctx_);
    reduce_polynomial(reduced, polynomial, parent_.ctx_,
                      parent_.modulus_);
    fmpz_mod_poly_swap(rep_.raw(), reduced.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::set_order_element(
        const OrderElement& element) noexcept {
    if (!is_defined() ||
        !same_order_parent(element.parent(), parent_.parent_order())) {
        return false;
    }

    flint::FmpzPoly reduced;
    if (!parent_.prime_.reduce(flint::FmpzPolyRef(reduced), element)) {
        return false;
    }
    return set_polynomial(flint::FmpzPolyConstRef(reduced));
}

bool ResidueFieldElement::set_element(const Element& element) noexcept {
    if (!is_defined() || parent_.parent_order() == nullptr ||
        parent_.parent_order()->parent() == nullptr ||
        !element.has_parent(*parent_.parent_order()->parent())) {
        return false;
    }

    const Order* order = parent_.parent_order();
    flint::FmpqMat coordinates(1, order->degree());
    flint::Fmpz denominator;
    flint::Fmpz denominator_inverse;
    flint::FmpzPoly scaled;
    flint::FmpzModPoly reduced(parent_.ctx_);

    if (!order->coordinates(flint::FmpqMatRef(coordinates), element)) {
        return false;
    }

    row_denominator_lcm(denominator, coordinates);
    if (fmpz_divisible(denominator.raw(), parent_.p_.raw()) != 0 ||
        fmpz_invmod(denominator_inverse.raw(), denominator.raw(),
                    parent_.p_.raw()) == 0) {
        return false;
    }

    row_scale_to_polynomial(scaled, coordinates, denominator);
    reduce_polynomial(reduced, flint::FmpzPolyConstRef(scaled),
                      parent_.ctx_, parent_.modulus_);
    fmpz_mod_poly_scalar_mul_fmpz(reduced.raw(), reduced.raw(),
                                  denominator_inverse.raw(),
                                  parent_.ctx_.raw());
    fmpz_mod_poly_swap(rep_.raw(), reduced.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::set_factored_element(
        const FactoredElement& element) noexcept {
    if (!is_defined() || !element.is_defined() ||
        parent_.parent_order() == nullptr ||
        parent_.parent_order()->parent() == nullptr ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*parent_.parent_order()->parent())) {
        return false;
    }

    ResidueFieldElement product(parent_);
    ResidueFieldElement factor(parent_);
    ResidueFieldElement power(parent_);
    flint::Fmpz exponent_value;
    if (!product.one()) {
        return false;
    }

    for (slong i = 0; i < element.length(); ++i) {
        const Element* factor_element = element.factor(i);
        slong exponent = 0;
        if (factor_element == nullptr ||
            !element.exponent(exponent, i) ||
            !factor.set_element(*factor_element)) {
            return false;
        }
        fmpz_set_si(exponent_value.raw(), exponent);
        if (!power.pow_fmpz(factor, flint::FmpzConstRef(exponent_value)) ||
            !product.multiply(product, power)) {
            return false;
        }
    }

    return set(product);
}

bool ResidueFieldElement::get_polynomial(
        flint::FmpzPolyRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }
    fmpz_mod_poly_get_fmpz_poly(out.raw(), rep_.raw(), parent_.ctx_.raw());
    return true;
}

std::optional<flint::FmpzPoly> ResidueFieldElement::polynomial() const noexcept {
    flint::FmpzPoly out;
    if (!get_polynomial(flint::FmpzPolyRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ResidueFieldElement::degree_one_scalar(
        flint::FmpzRef out) const noexcept {
    if (!is_defined() || out.raw() == nullptr || parent_.degree() != 1 ||
        fmpz_mod_poly_degree(rep_.raw(), parent_.ctx_.raw()) > 0) {
        return false;
    }

    fmpz_mod_poly_get_coeff_fmpz(out.raw(), rep_.raw(), 0,
                                 parent_.ctx_.raw());
    fmpz_mod(out.raw(), out.raw(), parent_.p_.raw());
    return true;
}

bool ResidueFieldElement::equal(
        const ResidueFieldElement& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.equal(other.parent_) &&
           fmpz_mod_poly_equal(rep_.raw(), other.rep_.raw(),
                               parent_.ctx_.raw()) != 0;
}

bool ResidueFieldElement::add(const ResidueFieldElement& left,
                              const ResidueFieldElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }
    fmpz_mod_poly_add(rep_.raw(), left.rep_.raw(), right.rep_.raw(),
                      parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::negate(const ResidueFieldElement& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.equal(input.parent_)) {
        return false;
    }
    fmpz_mod_poly_neg(rep_.raw(), input.rep_.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::subtract(const ResidueFieldElement& left,
                                   const ResidueFieldElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }
    fmpz_mod_poly_sub(rep_.raw(), left.rep_.raw(), right.rep_.raw(),
                      parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::multiply(const ResidueFieldElement& left,
                                   const ResidueFieldElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.equal(left.parent_) ||
        !left.parent_.equal(right.parent_)) {
        return false;
    }
    fmpz_mod_poly_mulmod(rep_.raw(), left.rep_.raw(), right.rep_.raw(),
                         parent_.modulus_.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::invert(
        const ResidueFieldElement& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.equal(input.parent_) ||
        fmpz_mod_poly_is_zero(input.rep_.raw(), parent_.ctx_.raw()) != 0) {
        return false;
    }

    flint::FmpzModPoly inverse(parent_.ctx_);
    const int ok = fmpz_mod_poly_invmod(
            inverse.raw(), input.rep_.raw(),
            parent_.modulus_.raw(), parent_.ctx_.raw());
    if (ok == 0) {
        return false;
    }
    fmpz_mod_poly_swap(rep_.raw(), inverse.raw(), parent_.ctx_.raw());
    return true;
}

bool ResidueFieldElement::pow_fmpz(
        const ResidueFieldElement& input,
        flint::FmpzConstRef exponent) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.equal(input.parent_)) {
        return false;
    }

    const flint::FmpzModCtx& ctx = parent_.ctx_;
    const flint::FmpzModPoly& modulus = parent_.modulus_;
    const fmpz_mod_poly_struct* base = input.rep_.raw();
    flint::Fmpz abs_exponent;
    flint::FmpzModPoly inverse(ctx);
    flint::FmpzModPoly result(ctx);

    if (fmpz_sgn(exponent.raw()) < 0) {
        if (fmpz_mod_poly_is_zero(input.rep_.raw(), ctx.raw()) != 0) {
            return false;
        }
        if (fmpz_mod_poly_invmod(inverse.raw(), input.rep_.raw(),
                                 modulus.raw(), ctx.raw()) == 0) {
            return false;
        }
        fmpz_neg(abs_exponent.raw(), exponent.raw());
        base = inverse.raw();
    } else {
        fmpz_set(abs_exponent.raw(), exponent.raw());
    }

    fmpz_mod_poly_powmod_fmpz_binexp(result.raw(), base, abs_exponent.raw(),
                                    modulus.raw(), ctx.raw());
    fmpz_mod_poly_swap(rep_.raw(), result.raw(), ctx.raw());
    return true;
}

bool ResidueFieldElement::multiplicative_order(
        flint::FmpzRef out) const noexcept {
    if (!is_defined() ||
        fmpz_mod_poly_is_zero(rep_.raw(), parent_.ctx_.raw()) != 0) {
        return false;
    }

    ResidueFieldElement power(parent_);
    flint::FmpzFactor factorization;
    flint::Fmpz q;
    flint::Fmpz candidate;
    flint::Fmpz test;
    ResidueFieldElement one(parent_);
    one.one();

    parent_.cardinality(flint::FmpzRef(q));
    fmpz_sub_ui(candidate.raw(), q.raw(), 1);
    fmpz_factor(factorization.raw(), candidate.raw());

    for (slong i = 0; i < factorization.raw()->num; ++i) {
        const fmpz* prime = factorization.raw()->p + i;
        while (fmpz_divisible(candidate.raw(), prime) != 0) {
            fmpz_divexact(test.raw(), candidate.raw(), prime);
            if (!power.pow_fmpz(*this, flint::FmpzConstRef(test)) ||
                !power.equal(one)) {
                break;
            }
            fmpz_set(candidate.raw(), test.raw());
        }
    }

    fmpz_set(out.raw(), candidate.raw());
    return true;
}

bool ResidueFieldElement::discrete_log(
        flint::FmpzRef out,
        const ResidueFieldElement& base) const noexcept {
    if (!is_defined() || !base.is_defined() ||
        !parent_.equal(base.parent_) ||
        fmpz_mod_poly_is_zero(rep_.raw(), parent_.ctx_.raw()) != 0 ||
        fmpz_mod_poly_is_zero(base.rep_.raw(), parent_.ctx_.raw()) != 0) {
        return false;
    }

    ResidueFieldElement current(parent_);
    flint::Fmpz i;
    flint::Fmpz order;
    if (!base.multiplicative_order(flint::FmpzRef(order))) {
        return false;
    }

    current.one();
    for (fmpz_zero(i.raw()); fmpz_cmp(i.raw(), order.raw()) < 0;
         fmpz_add_ui(i.raw(), i.raw(), 1)) {
        if (current.equal(*this)) {
            fmpz_set(out.raw(), i.raw());
            return true;
        }
        if (!current.multiply(current, base)) {
            return false;
        }
    }

    return false;
}

bool ResidueFieldElement::quotient_log_mod_prime(
        flint::FmpzRef out,
        flint::FmpzConstRef ell) const noexcept {
    if (!is_defined() ||
        fmpz_mod_poly_is_zero(rep_.raw(), parent_.ctx_.raw()) != 0 ||
        fmpz_sgn(ell.raw()) <= 0 ||
        fmpz_is_prime(ell.raw()) == 0) {
        return false;
    }

    flint::Fmpz q;
    flint::Fmpz group_order;
    flint::Fmpz cofactor;
    flint::Fmpz local_log;
    ResidueFieldElement generator(parent_);
    ResidueFieldElement quotient_generator(parent_);
    ResidueFieldElement image(parent_);

    parent_.cardinality(flint::FmpzRef(q));
    fmpz_sub_ui(group_order.raw(), q.raw(), 1);
    if (fmpz_divisible(group_order.raw(), ell.raw()) == 0) {
        return false;
    }
    fmpz_divexact(cofactor.raw(), group_order.raw(), ell.raw());

    if (!parent_.multiplicative_generator(generator) ||
        !quotient_generator.pow_fmpz(generator,
                                     flint::FmpzConstRef(cofactor)) ||
        !image.pow_fmpz(*this, flint::FmpzConstRef(cofactor)) ||
        !image.discrete_log(flint::FmpzRef(local_log),
                            quotient_generator)) {
        return false;
    }

    fmpz_set(out.raw(), local_log.raw());
    return true;
}

ResidueFieldQuotientLog::ResidueFieldQuotientLog(
        const ResidueField& parent) noexcept {
    define(parent);
}

ResidueFieldQuotientLog::~ResidueFieldQuotientLog() noexcept = default;

ResidueFieldQuotientLog::ResidueFieldQuotientLog(
        ResidueFieldQuotientLog&& other) noexcept {
    swap(other);
}

ResidueFieldQuotientLog& ResidueFieldQuotientLog::operator=(
        ResidueFieldQuotientLog&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueFieldQuotientLog::swap(
        ResidueFieldQuotientLog& other) noexcept {
    parent_.swap(other.parent_);
    ell_.swap(other.ell_);
    order_.swap(other.order_);
    cofactor_.swap(other.cofactor_);
    generator_.swap(other.generator_);
    quotient_generator_.swap(other.quotient_generator_);
    std::swap(is_set_, other.is_set_);
}

void ResidueFieldQuotientLog::clear() noexcept {
    ResidueFieldQuotientLog empty;
    swap(empty);
}

bool ResidueFieldQuotientLog::define(
        const ResidueField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    ResidueFieldQuotientLog next;
    if (!next.parent_.set(parent)) {
        return false;
    }
    next.generator_ = ResidueFieldElement(next.parent_);
    next.quotient_generator_ = ResidueFieldElement(next.parent_);
    if (!next.generator_.is_defined() ||
        !next.quotient_generator_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool ResidueFieldQuotientLog::set(
        const ResidueFieldQuotientLog& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    ResidueFieldQuotientLog next(other.parent_);
    if (!next.is_defined()) {
        return false;
    }
    fmpz_set(next.ell_.raw(), other.ell_.raw());
    fmpz_set(next.order_.raw(), other.order_.raw());
    fmpz_set(next.cofactor_.raw(), other.cofactor_.raw());
    if (!next.generator_.set(other.generator_) ||
        !next.quotient_generator_.set(other.quotient_generator_)) {
        return false;
    }
    next.is_set_ = other.is_set_;

    swap(next);
    return true;
}

bool ResidueFieldQuotientLog::is_defined() const noexcept {
    return parent_.is_defined() &&
           generator_.is_defined() && quotient_generator_.is_defined();
}

bool ResidueFieldQuotientLog::is_set() const noexcept {
    return is_defined() && is_set_;
}

const ResidueField* ResidueFieldQuotientLog::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

bool ResidueFieldQuotientLog::degree_one_scalar_(
        flint::Fmpz& out,
        const ResidueFieldElement& element) noexcept {
    return element.degree_one_scalar(flint::FmpzRef(out));
}

bool ResidueFieldQuotientLog::set_ell(
        flint::FmpzConstRef ell_in) noexcept {
    if (!is_defined() || fmpz_sgn(ell_in.raw()) <= 0 ||
        fmpz_is_prime(ell_in.raw()) == 0) {
        return false;
    }

    ResidueFieldQuotientLog next(parent_);
    flint::Fmpz q;
    flint::Fmpz ell_free_cofactor;
    flint::Fmpz ell_primary_order;
    flint::Fmpz quotient_exponent;
    flint::Fmpz k;
    ResidueFieldElement candidate(parent_);
    ResidueFieldElement primary_part(parent_);
    ResidueFieldElement test(parent_);
    ResidueFieldElement one(parent_);
    if (!next.is_defined() || !parent_.cardinality(flint::FmpzRef(q))) {
        return false;
    }

    fmpz_sub_ui(next.order_.raw(), q.raw(), 1);
    if (fmpz_divisible(next.order_.raw(), ell_in.raw()) == 0) {
        return false;
    }

    fmpz_set(next.ell_.raw(), ell_in.raw());
    fmpz_divexact(next.cofactor_.raw(), next.order_.raw(), ell_in.raw());
    const slong ell_valuation = fmpz_remove(
            ell_free_cofactor.raw(), next.order_.raw(), ell_in.raw());
    if (ell_valuation <= 0 ||
        !candidate.is_defined() || !primary_part.is_defined() ||
        !test.is_defined() || !one.is_defined() || !one.one()) {
        return false;
    }
    fmpz_divexact(ell_primary_order.raw(), next.order_.raw(),
                  ell_free_cofactor.raw());
    fmpz_divexact(quotient_exponent.raw(), ell_primary_order.raw(),
                  ell_in.raw());

    // reference `compute_candidates_for_saturate` only needs the ell-primary
    // quotient of the residue-field unit group: it raises a random element to
    // the ell-free cofactor and rejects it if the result dies modulo ell.  The
    // deterministic scan here is the same test over a fixed candidate order.
    bool found = false;
    for (fmpz_one(k.raw()); fmpz_cmp(k.raw(), q.raw()) < 0;
         fmpz_add_ui(k.raw(), k.raw(), 1)) {
        if (!set_residue_field_element_from_index(
                    candidate, parent_, flint::FmpzConstRef(k)) ||
            !primary_part.pow_fmpz(
                    candidate, flint::FmpzConstRef(ell_free_cofactor)) ||
            !test.pow_fmpz(primary_part,
                           flint::FmpzConstRef(quotient_exponent))) {
            return false;
        }
        if (!test.equal(one)) {
            if (!next.generator_.set(primary_part) ||
                !next.quotient_generator_.set(test)) {
                return false;
            }
            found = true;
            break;
        }
    }
    if (!found) {
        return false;
    }

    next.is_set_ = true;
    swap(next);
    return true;
}

bool ResidueFieldQuotientLog::ell(flint::FmpzRef out) const noexcept {
    if (!is_set()) {
        return false;
    }
    fmpz_set(out.raw(), ell_.raw());
    return true;
}

bool ResidueFieldQuotientLog::generator(
        ResidueFieldElement& out) const noexcept {
    if (!is_set()) {
        return false;
    }
    if (!out.is_defined()) {
        out = ResidueFieldElement(parent_);
        if (!out.is_defined()) {
            return false;
        }
    }
    return out.parent() != nullptr && out.parent()->equal(parent_) &&
           out.set(generator_);
}

bool ResidueFieldQuotientLog::quotient_generator(
        ResidueFieldElement& out) const noexcept {
    if (!is_set()) {
        return false;
    }
    if (!out.is_defined()) {
        out = ResidueFieldElement(parent_);
        if (!out.is_defined()) {
            return false;
        }
    }
    return out.parent() != nullptr && out.parent()->equal(parent_) &&
           out.set(quotient_generator_);
}

bool ResidueFieldQuotientLog::apply(
        flint::FmpzRef out,
        const ResidueFieldElement& element) const noexcept {
    if (!is_set() || element.parent() == nullptr ||
        !element.parent()->equal(parent_) || element.is_zero()) {
        return false;
    }

    if (parent_.degree() == 1) {
        flint::Fmpz value;
        flint::Fmpz quotient_generator;
        if (degree_one_scalar_(value, element) &&
            degree_one_scalar_(quotient_generator, quotient_generator_)) {
            if (quotient_log_mod_prime_apply_ui(
                        out, flint::FmpzConstRef(value),
                        flint::FmpzConstRef(cofactor_),
                        flint::FmpzConstRef(quotient_generator),
                        flint::FmpzConstRef(parent_.p_),
                        flint::FmpzConstRef(ell_))) {
                return true;
            }
            return detail::quotient_log_mod_prime_apply(
                    out, flint::FmpzConstRef(value),
                    flint::FmpzConstRef(cofactor_),
                    flint::FmpzConstRef(quotient_generator),
                    flint::FmpzConstRef(parent_.p_),
                    flint::FmpzConstRef(ell_));
        }
    }

    ResidueFieldElement image(parent_);
    ResidueFieldElement current(parent_);
    flint::Fmpz i;
    if (!image.pow_fmpz(element, flint::FmpzConstRef(cofactor_)) ||
        !current.one()) {
        return false;
    }

    for (fmpz_zero(i.raw()); fmpz_cmp(i.raw(), ell_.raw()) < 0;
         fmpz_add_ui(i.raw(), i.raw(), 1)) {
        if (current.equal(image)) {
            fmpz_set(out.raw(), i.raw());
            return true;
        }
        if (!current.multiply(current, quotient_generator_)) {
            return false;
        }
    }

    return false;
}

}  // namespace silex
