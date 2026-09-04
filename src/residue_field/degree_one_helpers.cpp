#include "residue_field_internal.hpp"

#include <limits>

#include <flint/fmpq_poly.h>
#include <flint/fmpz_factor.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>

#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpq.hpp>

namespace silex::detail {
namespace {

class NmodPoly {
public:
    explicit NmodPoly(ulong modulus) noexcept {
        nmod_poly_init(value_, modulus);
    }

    ~NmodPoly() noexcept {
        nmod_poly_clear(value_);
    }

    NmodPoly(const NmodPoly&) = delete;
    NmodPoly& operator=(const NmodPoly&) = delete;

    nmod_poly_t& raw() noexcept {
        return value_;
    }

private:
    nmod_poly_t value_;
};

bool fmpq_poly_get_nmod_poly(NmodPoly& out,
                             flint::FmpqPolyConstRef polynomial,
                             ulong modulus,
                             nmod_t mod) noexcept {
    flint::Fmpq coefficient;
    const slong degree = flint::fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        flint::fmpq_poly_get_coeff_fmpq(flint::FmpqRef(coefficient),
                                        polynomial, i);
        const ulong denominator =
                fmpz_fdiv_ui(fmpq_denref(coefficient.raw()), modulus);
        if (denominator == 0) {
            return false;
        }
        const ulong numerator =
                fmpz_fdiv_ui(fmpq_numref(coefficient.raw()), modulus);
        const ulong coeff = nmod_mul(numerator, n_invmod(denominator, modulus),
                                     mod);
        nmod_poly_set_coeff_ui(out.raw(), i, coeff);
    }
    return true;
}

bool pow_mod_abs_slong(flint::Fmpz& out,
                       const flint::Fmpz& base,
                       slong exponent,
                       flint::FmpzConstRef p) noexcept {
    if (exponent < 0) {
        if (exponent == std::numeric_limits<slong>::min()) {
            return false;
        }
        exponent = -exponent;
    }
    ::fmpz_powm_ui(out.raw(), base.raw(), static_cast<ulong>(exponent),
                   p.raw());
    return true;
}

bool pow_mod_abs_slong_ui(ulong& out,
                          ulong base,
                          slong exponent,
                          ulong modulus,
                          ulong modulus_inverse) noexcept {
    if (exponent < 0) {
        if (exponent == std::numeric_limits<slong>::min()) {
            return false;
        }
        exponent = -exponent;
    }
    out = n_powmod2_ui_preinv(base, static_cast<ulong>(exponent), modulus,
                              modulus_inverse);
    return true;
}

bool multiplicative_order_mod_prime(flint::Fmpz& out,
                                    flint::FmpzConstRef base,
                                    flint::FmpzConstRef p) noexcept {
    flint::Fmpz reduced_base;
    ::fmpz_mod(reduced_base.raw(), base.raw(), p.raw());
    if (::fmpz_is_zero(reduced_base.raw()) != 0) {
        return false;
    }

    flint::FmpzFactor factorization;
    flint::Fmpz candidate;
    flint::Fmpz test;
    flint::Fmpz power;
    ::fmpz_sub_ui(candidate.raw(), p.raw(), 1);
    ::fmpz_factor(factorization.raw(), candidate.raw());

    for (slong i = 0; i < factorization.raw()->num; ++i) {
        const fmpz* prime = factorization.raw()->p + i;
        while (::fmpz_divisible(candidate.raw(), prime) != 0) {
            ::fmpz_divexact(test.raw(), candidate.raw(), prime);
            ::fmpz_powm(power.raw(), reduced_base.raw(), test.raw(), p.raw());
            if (!flint::fmpz_is_one(flint::FmpzConstRef(power))) {
                break;
            }
            ::fmpz_set(candidate.raw(), test.raw());
        }
    }

    ::fmpz_set(out.raw(), candidate.raw());
    return true;
}

bool primitive_root_mod_prime(flint::Fmpz& out,
                              flint::FmpzConstRef p) noexcept {
    if (!flint::fmpz_is_prime(p)) {
        return false;
    }

    flint::Fmpz group_order;
    flint::Fmpz candidate;
    flint::Fmpz order;
    ::fmpz_sub_ui(group_order.raw(), p.raw(), 1);
    for (::fmpz_one(candidate.raw());
         ::fmpz_cmp(candidate.raw(), p.raw()) < 0;
         ::fmpz_add_ui(candidate.raw(), candidate.raw(), 1)) {
        if (multiplicative_order_mod_prime(order, flint::FmpzConstRef(candidate),
                                           p) &&
            ::fmpz_equal(order.raw(), group_order.raw()) != 0) {
            ::fmpz_set(out.raw(), candidate.raw());
            return true;
        }
    }

    return false;
}

}  // namespace

bool degree_one_prime_root_mod_p(flint::Fmpz& root,
                                 flint::Fmpz& p,
                                 const PrimeIdeal& prime) noexcept {
    flint::FmpzPoly residue;
    flint::Fmpz leading;
    flint::Fmpz constant;
    flint::Fmpz inverse_leading;
    if (prime.residue_degree() != 1 ||
        !prime.rational_prime(flint::FmpzRef(p)) ||
        !prime.residue_polynomial(flint::FmpzPolyRef(residue)) ||
        ::fmpz_poly_degree(residue.raw()) != 1) {
        return false;
    }

    ::fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    ::fmpz_poly_get_coeff_fmpz(leading.raw(), residue.raw(), 1);
    ::fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (::fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }

    ::fmpz_neg(root.raw(), constant.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());
    ::fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());

    const Order* order = prime.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    flint::Fmpz radicand;
    flint::Fmpz conductor;
    if (field != nullptr && order->degree() == 2 && order->is_maximal() &&
        field->backend_kind() == NumberFieldBackendKind::quadratic &&
        field->quadratic_radicand(flint::FmpzRef(radicand)) &&
        order->quadratic_conductor(flint::FmpzRef(conductor)) &&
        flint::fmpz_is_one(conductor) &&
        flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        // Direct maximal-quadratic residue polynomials use reference's integral
        // generator omega=(1+theta)/2.  The degree-one fast paths evaluate
        // ambient Element polynomials in the supplied generator theta, so
        // return theta=2*omega-1 here.  Order-coordinate residue scans keep
        // using PrimeIdeal's stored omega root directly.
        ::fmpz_mul_ui(root.raw(), root.raw(), 2);
        ::fmpz_sub_ui(root.raw(), root.raw(), 1);
        ::fmpz_mod(root.raw(), root.raw(), p.raw());
    }
    return true;
}

bool evaluate_fmpq_poly_at_degree_one_prime(
        flint::Fmpz& out,
        flint::FmpqPolyConstRef polynomial,
        flint::FmpzConstRef root,
        flint::FmpzConstRef p) noexcept {
    flint::Fmpq coefficient;
    flint::Fmpz numerator;
    flint::Fmpz denominator;
    flint::Fmpz inverse_denominator;
    flint::Fmpz coefficient_mod;

    ::fmpz_zero(out.raw());
    const slong degree = ::fmpq_poly_degree(polynomial.raw());
    for (slong i = degree; i >= 0; --i) {
        ::fmpz_mul(out.raw(), out.raw(), root.raw());
        ::fmpz_mod(out.raw(), out.raw(), p.raw());

        ::fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial.raw(), i);
        ::fmpz_set(numerator.raw(), fmpq_numref(coefficient.raw()));
        ::fmpz_set(denominator.raw(), fmpq_denref(coefficient.raw()));
        ::fmpz_mod(denominator.raw(), denominator.raw(), p.raw());
        if (::fmpz_invmod(inverse_denominator.raw(), denominator.raw(),
                          p.raw()) == 0) {
            return false;
        }
        ::fmpz_mod(numerator.raw(), numerator.raw(), p.raw());
        ::fmpz_mul(coefficient_mod.raw(), numerator.raw(),
                   inverse_denominator.raw());
        ::fmpz_mod(coefficient_mod.raw(), coefficient_mod.raw(), p.raw());
        ::fmpz_add(out.raw(), out.raw(), coefficient_mod.raw());
        ::fmpz_mod(out.raw(), out.raw(), p.raw());
    }

    return true;
}

bool quotient_log_mod_prime_setup(flint::Fmpz& cofactor,
                                  flint::Fmpz& quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_is_prime(p) || !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::Fmpz group_order;
    flint::Fmpz generator;
    ::fmpz_sub_ui(group_order.raw(), p.raw(), 1);
    if (::fmpz_divisible(group_order.raw(), ell.raw()) == 0 ||
        !primitive_root_mod_prime(generator, p)) {
        return false;
    }

    ::fmpz_divexact(cofactor.raw(), group_order.raw(), ell.raw());
    ::fmpz_powm(quotient_generator.raw(), generator.raw(), cofactor.raw(),
                p.raw());
    return !flint::fmpz_is_one(flint::FmpzConstRef(quotient_generator));
}

bool quotient_log_mod_prime_apply(flint::FmpzRef out,
                                  flint::FmpzConstRef value,
                                  flint::FmpzConstRef cofactor,
                                  flint::FmpzConstRef quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept {
    flint::Fmpz reduced;
    flint::Fmpz image;
    flint::Fmpz current;
    flint::Fmpz i;
    ::fmpz_mod(reduced.raw(), value.raw(), p.raw());
    if (::fmpz_is_zero(reduced.raw()) != 0) {
        return false;
    }

    ::fmpz_powm(image.raw(), reduced.raw(), cofactor.raw(), p.raw());
    ::fmpz_one(current.raw());
    for (::fmpz_zero(i.raw()); ::fmpz_cmp(i.raw(), ell.raw()) < 0;
         ::fmpz_add_ui(i.raw(), i.raw(), 1)) {
        if (::fmpz_equal(current.raw(), image.raw()) != 0) {
            ::fmpz_set(out.raw(), i.raw());
            return true;
        }
        ::fmpz_mul(current.raw(), current.raw(), quotient_generator.raw());
        ::fmpz_mod(current.raw(), current.raw(), p.raw());
    }

    return false;
}

bool factored_value_at_degree_one_root(flint::Fmpz& out,
                                       const FactoredElement& element,
                                       flint::FmpzConstRef p,
                                       flint::FmpzConstRef root) noexcept {
    if (!element.is_defined() || !flint::fmpz_is_prime(p)) {
        return false;
    }

    flint::Fmpz base;
    flint::Fmpz inverse_base;
    flint::Fmpz power;
    flint::fmpz_set_ui(flint::FmpzRef(out), 1);

    for (const auto& entry : element.factors()) {
        if (entry.exponent == 0) {
            continue;
        }

        flint::FmpqPoly polynomial;
        if (!entry.factor.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) ||
            !evaluate_fmpq_poly_at_degree_one_prime(
                    base, flint::FmpqPolyConstRef(polynomial), root, p)) {
            return false;
        }
        if (::fmpz_is_zero(base.raw()) != 0) {
            return false;
        }

        const flint::Fmpz* exponent_base = &base;
        if (entry.exponent < 0) {
            if (::fmpz_invmod(inverse_base.raw(), base.raw(), p.raw()) == 0) {
                return false;
            }
            exponent_base = &inverse_base;
        }
        if (!pow_mod_abs_slong(power, *exponent_base, entry.exponent, p)) {
            return false;
        }
        ::fmpz_mul(out.raw(), out.raw(), power.raw());
        ::fmpz_mod(out.raw(), out.raw(), p.raw());
    }

    return ::fmpz_is_zero(out.raw()) == 0;
}

bool factored_value_at_degree_one_root_nmod(ulong& out,
                                           const FactoredElement& element,
                                           ulong p,
                                           ulong root) noexcept {
    if (!element.is_defined() || p < 2) {
        return false;
    }

    nmod_t mod;
    nmod_init(&mod, p);
    const ulong p_inverse = n_preinvert_limb(p);
    ulong product = 1 % p;

    for (const auto& entry : element.factors()) {
        if (entry.exponent == 0) {
            continue;
        }

        flint::FmpqPoly polynomial;
        NmodPoly reduced(p);
        if (!entry.factor.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) ||
            !fmpq_poly_get_nmod_poly(
                    reduced, flint::FmpqPolyConstRef(polynomial), p, mod)) {
            return false;
        }

        ulong base = nmod_poly_evaluate_nmod(reduced.raw(), root);
        if (base == 0) {
            return false;
        }
        if (entry.exponent < 0) {
            base = n_invmod(base, p);
        }

        ulong power = 1;
        if (!pow_mod_abs_slong_ui(power, base, entry.exponent, p, p_inverse)) {
            return false;
        }
        product = nmod_mul(product, power, mod);
    }

    out = product;
    return product != 0;
}

bool element_value_at_degree_one_root(flint::Fmpz& out,
                                      const Element& element,
                                      flint::FmpzConstRef p,
                                      flint::FmpzConstRef root) noexcept {
    flint::FmpqPoly polynomial;
    return element.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) &&
           evaluate_fmpq_poly_at_degree_one_prime(
                   out, flint::FmpqPolyConstRef(polynomial), root, p) &&
           ::fmpz_is_zero(out.raw()) == 0;
}

}  // namespace silex::detail
