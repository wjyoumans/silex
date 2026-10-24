#include "residue_field_internal.hpp"

#include "../prime_ideal/prime_ideal_internal.hpp"

#include <limits>
#include <utility>

#include <flint/fmpq_poly.h>
#include <flint/fmpz_factor.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
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

// The stored form of `value` in `field` is the integer polynomial num(x) over
// the positive integer den that nf_elem_get_fmpq_poly copies out: one
// coefficient for NF_LINEAR, two for NF_QUADRATIC, and the fmpq_poly
// otherwise.  Sets `numerator` = num(root) mod p and `denominator` = den mod
// p, reading the coefficients in place.  `root` must be reduced mod p.
// Horner's rule on the reduced coefficients gives the residue that
// nmod_poly_evaluate_nmod gives on their reduction.
void stored_form_value_mod_p(ulong& numerator,
                             ulong& denominator,
                             const nf_elem_struct* value,
                             const nf_struct* field,
                             ulong root,
                             nmod_t mod) noexcept {
    const fmpz* coefficients = nullptr;
    const fmpz* stored_denominator = nullptr;
    slong length = 0;
    if ((field->flag & NF_LINEAR) != 0) {
        coefficients = LNF_ELEM_NUMREF(value);
        stored_denominator = LNF_ELEM_DENREF(value);
        length = 1;
    } else if ((field->flag & NF_QUADRATIC) != 0) {
        coefficients = QNF_ELEM_NUMREF(value);
        stored_denominator = QNF_ELEM_DENREF(value);
        length = 2;
    } else {
        coefficients = NF_ELEM_NUMREF(value);
        stored_denominator = NF_ELEM_DENREF(value);
        length = NF_ELEM(value)->length;
    }

    ulong accumulator = 0;
    for (slong i = length - 1; i >= 0; --i) {
        accumulator = nmod_add(nmod_mul(accumulator, root, mod),
                               ::fmpz_get_nmod(coefficients + i, mod), mod);
    }
    numerator = accumulator;
    denominator = ::fmpz_get_nmod(stored_denominator, mod);
}

// The factor's value at `root` from its canonical fmpq_poly, inverting each
// coefficient's own denominator: fails exactly when one of them is divisible
// by p.  On success `numerator` is the value and `denominator` is 1.  This is
// the route factored_value_at_degree_one_root_nmod takes when the stored
// common denominator is divisible by p, so its failures stay those of the
// per-coefficient evaluation (and of evaluate_fmpq_poly_at_degree_one_prime).
bool fmpq_poly_value_at_root_nmod(ulong& numerator,
                                  ulong& denominator,
                                  const Element& factor,
                                  ulong p,
                                  ulong root,
                                  nmod_t mod) noexcept {
    flint::FmpqPoly polynomial;
    NmodPoly reduced(p);
    if (!factor.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) ||
        !fmpq_poly_get_nmod_poly(reduced, flint::FmpqPolyConstRef(polynomial),
                                 p, mod)) {
        return false;
    }
    numerator = nmod_poly_evaluate_nmod(reduced.raw(), root);
    denominator = 1;
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

}  // namespace

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

ulong primitive_root_mod_prime_ui(ulong p) noexcept {
    // Same smallest primitive root as primitive_root_mod_prime: candidate c
    // has order p - 1 exactly when c^((p - 1) / q) != 1 for every prime q
    // dividing p - 1, so p - 1 is factored once instead of once per
    // candidate (the test of FLINT ulong_extras/primitive_root_prime.c,
    // n_primitive_root_prime_prefactor, run from the smallest candidate).
    // For p = 2 the group is trivial and the answer is 1.
    if (p < 2) {
        return 0;
    }
    if (p == 2) {
        return 1;
    }

    const ulong group_order = p - 1;
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, group_order, 1);
    const ulong p_inverse = n_preinvert_limb(p);
    for (ulong candidate = 2; candidate < p; ++candidate) {
        bool generates = true;
        for (int i = 0; i < factors.num; ++i) {
            if (n_powmod2_ui_preinv(candidate, group_order / factors.p[i], p,
                                    p_inverse) == 1) {
                generates = false;
                break;
            }
        }
        if (generates) {
            return candidate;
        }
    }

    return 0;
}

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

    // Read the residue variable recorded in the prime when it was built, not
    // the order's current maximality record (see
    // PrimeIdeal::residue_polynomial).
    const Order* order = prime.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    flint::Fmpz radicand;
    if (field != nullptr &&
        PrimeIdealAccess::residue_uses_integral_generator(prime) &&
        field->quadratic_radicand(flint::FmpzRef(radicand)) &&
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

    if (flint::fmpz_abs_fits_ui(p)) {
        // Word-size route, with the results of the fmpz route below.  A prime
        // ell above the word size cannot divide p - 1.
        const ulong p_ui = flint::fmpz_get_ui(p);
        const ulong group_order_ui = p_ui - 1;
        if (!flint::fmpz_abs_fits_ui(ell)) {
            return false;
        }
        const ulong ell_ui = flint::fmpz_get_ui(ell);
        if (group_order_ui % ell_ui != 0) {
            return false;
        }
        const ulong generator_ui = primitive_root_mod_prime_ui(p_ui);
        if (generator_ui == 0) {
            return false;
        }
        const ulong cofactor_ui = group_order_ui / ell_ui;
        const ulong quotient_generator_ui = n_powmod2_ui_preinv(
                generator_ui, cofactor_ui, p_ui, n_preinvert_limb(p_ui));
        ::fmpz_set_ui(cofactor.raw(), cofactor_ui);
        ::fmpz_set_ui(quotient_generator.raw(), quotient_generator_ui);
        return quotient_generator_ui != 1;
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
    const NumberField* parent = element.parent();
    const nf_struct* field =
            parent == nullptr ? nullptr : parent->raw_flint_field();
    if (field == nullptr) {
        return false;
    }

    nmod_t mod;
    nmod_init(&mod, p);
    const ulong p_inverse = n_preinvert_limb(p);
    const ulong reduced_root = root % p;
    // The product is kept as numerator / denominator, so one inversion at
    // the end replaces one per factor.
    ulong numerator_product = 1;
    ulong denominator_product = 1;

    for (const auto& entry : element.factors()) {
        if (entry.exponent == 0) {
            continue;
        }

        // base = numerator / denominator is the factor's value at the root.
        ulong numerator = 0;
        ulong denominator = 0;
        stored_form_value_mod_p(numerator, denominator,
                                entry.factor.raw_flint_element(), field,
                                reduced_root, mod);
        // The fallback below succeeds only for a non-canonical stored form
        // (p divides den but no reduced coefficient denominator); canonical
        // FLINT elements fail on both routes.
        if (denominator == 0 &&
            !fmpq_poly_value_at_root_nmod(numerator, denominator,
                                          entry.factor, p, reduced_root,
                                          mod)) {
            return false;
        }
        if (numerator == 0) {
            return false;
        }

        ulong numerator_power = 1;
        ulong denominator_power = 1;
        if (!pow_mod_abs_slong_ui(numerator_power, numerator, entry.exponent,
                                  p, p_inverse) ||
            !pow_mod_abs_slong_ui(denominator_power, denominator,
                                  entry.exponent, p, p_inverse)) {
            return false;
        }
        if (entry.exponent < 0) {
            std::swap(numerator_power, denominator_power);
        }
        numerator_product = nmod_mul(numerator_product, numerator_power, mod);
        denominator_product =
                nmod_mul(denominator_product, denominator_power, mod);
    }

    // For prime p every power above is a unit, so the inverse exists; a
    // composite p (outside the contract) whose denominator is not a unit
    // fails here instead of in n_invmod.
    ulong denominator_inverse = 0;
    if (denominator_product == 0 ||
        n_gcdinv(&denominator_inverse, denominator_product, p) != 1) {
        return false;
    }
    const ulong product =
            nmod_mul(numerator_product, denominator_inverse, mod);
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
