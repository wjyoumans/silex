#include <silex/unit.hpp>

#include "../element/element_internal.hpp"

#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>

#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/signature.hpp>

namespace silex {
namespace {

constexpr slong kTrivialGoodPrimeLimit = 64;

bool set_rational_half(Element& out, slong constant, slong linear) noexcept {
    flint::FmpqPoly polynomial;
    flint::Fmpq coeff;

    fmpq_set_si(coeff.raw(), constant, 2);
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 0, coeff.raw());
    fmpq_set_si(coeff.raw(), linear, 2);
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, coeff.raw());
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool defining_fmpz_poly(flint::FmpzPoly& out,
                        const NumberField& field) noexcept {
    const nf_struct* raw = field.raw_flint_field();
    if (raw == nullptr) {
        return false;
    }

    flint::Fmpz den;
    flint::Fmpz lc;
    fmpq_poly_get_denominator(den.raw(), raw->pol);
    if (!flint::fmpz_is_one(den)) {
        return false;
    }

    fmpq_poly_get_numerator(out.raw(), raw->pol);
    fmpz_poly_get_coeff_fmpz(lc.raw(), out.raw(),
                             fmpz_poly_degree(out.raw()));
    return flint::fmpz_is_one(lc);
}

bool c_parity_defining_polynomial_is_x4_plus_one(
        const NumberField& field) noexcept {
    const nf_struct* raw = field.raw_flint_field();
    if (raw == nullptr || fmpq_poly_degree(raw->pol) != 4) {
        return false;
    }

    flint::Fmpq coeff;
    fmpq_poly_get_coeff_fmpq(coeff.raw(), raw->pol, 0);
    bool ok = flint::fmpq_equal_si(coeff, 1);
    fmpq_poly_get_coeff_fmpq(coeff.raw(), raw->pol, 1);
    ok = ok && fmpq_is_zero(coeff.raw()) != 0;
    fmpq_poly_get_coeff_fmpq(coeff.raw(), raw->pol, 2);
    ok = ok && fmpq_is_zero(coeff.raw()) != 0;
    fmpq_poly_get_coeff_fmpq(coeff.raw(), raw->pol, 3);
    ok = ok && fmpq_is_zero(coeff.raw()) != 0;
    fmpq_poly_get_coeff_fmpq(coeff.raw(), raw->pol, 4);
    return ok && flint::fmpq_equal_si(coeff, 1);
}

bool c_parity_roots_x4_plus_one(Element& generator,
                                const NumberField& field) noexcept {
    if (!c_parity_defining_polynomial_is_x4_plus_one(field)) {
        return false;
    }

    Element theta(field);
    Element theta2(field);
    Element theta4(field);
    Element theta8(field);
    if (!theta.is_defined() || !theta2.is_defined() ||
        !theta4.is_defined() || !theta8.is_defined() ||
        !theta.gen() ||
        !theta2.multiply(theta, theta) ||
        !theta4.multiply(theta2, theta2) ||
        !theta8.multiply(theta4, theta4) ||
        !theta4.equal_si(-1) || !theta8.equal_si(1)) {
        return false;
    }

    return generator.set(theta);
}

bool primitive_fourth_root(Element& generator,
                           const NumberField& field) noexcept {
    Element minus_one(field);
    Element root(field);
    bool is_square = false;
    if (!minus_one.is_defined() || !root.is_defined() ||
        !minus_one.set_si(-1) ||
        !minus_one.is_square(is_square, root) ||
        !is_square) {
        return false;
    }
    return generator.set(root);
}

bool primitive_sixth_root(Element& generator,
                          const NumberField& field) noexcept {
    Element minus_three(field);
    Element root(field);
    Element numerator(field);
    bool is_square = false;
    if (!minus_three.is_defined() || !root.is_defined() ||
        !numerator.is_defined() ||
        !minus_three.set_si(-3) ||
        !minus_three.is_square(is_square, root) ||
        !is_square ||
        !numerator.add_si(root, 1)) {
        return false;
    }
    return generator.scalar_div_si(numerator, 2);
}

bool roots_cyclotomic_square_factors(flint::FmpzRef order,
                                     Element& generator,
                                     const NumberField& field) noexcept {
    // Source trace: reference `TorsionUnits.jl:_torsion_units_gen` factors a
    // torsion-order divisor and searches for roots of cyclotomic prime-power
    // factors with `_roots_hensel`.  Native reuses the existing reference-shaped
    // exact square-root/Hensel primitive for the quadratic cyclotomic factors
    // Phi_4 (sqrt(-1)) and Phi_6 (sqrt(-3)) after the good-prime divisor has
    // failed to prove that the torsion order is two.
    Element zeta4(field);
    Element zeta6(field);
    Element zeta3(field);
    Element zeta12(field);
    const bool have_zeta4 = zeta4.is_defined() &&
                            primitive_fourth_root(zeta4, field);
    const bool have_zeta6 = zeta6.is_defined() &&
                            primitive_sixth_root(zeta6, field);

    if (have_zeta4 && have_zeta6 &&
        zeta3.is_defined() && zeta12.is_defined() &&
        zeta3.multiply(zeta6, zeta6) &&
        zeta12.multiply(zeta4, zeta3)) {
        fmpz_set_ui(order.raw(), 12);
        return generator.set(zeta12);
    }
    if (have_zeta6) {
        fmpz_set_ui(order.raw(), 6);
        return generator.set(zeta6);
    }
    if (have_zeta4) {
        fmpz_set_ui(order.raw(), 4);
        return generator.set(zeta4);
    }
    return false;
}

ulong factor_residue_degree_gcd(const flint::FmpzPoly& polynomial,
                                ulong p) noexcept {
    flint::FmpzModCtx ctx(p);
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPolyFactor factorization(ctx);

    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), polynomial.raw(), ctx.raw());
    if (fmpz_mod_poly_is_squarefree(reduced.raw(), ctx.raw()) == 0) {
        return 0;
    }
    fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), ctx.raw());

    ulong degree_gcd = 0;
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        const slong d =
                fmpz_mod_poly_degree(factorization.raw()->poly + i,
                                     ctx.raw());
        if (d <= 0) {
            return 0;
        }
        degree_gcd = degree_gcd == 0
                ? static_cast<ulong>(d)
                : n_gcd(degree_gcd, static_cast<ulong>(d));
    }
    return degree_gcd;
}

bool roots_trivial_by_good_primes(flint::FmpzRef order,
                                  Element& generator,
                                  const NumberField& field) noexcept {
    flint::FmpzPoly polynomial;
    if (!defining_fmpz_poly(polynomial, field)) {
        return false;
    }

    flint::Fmpz discriminant;
    flint::Fmpz local_bound;
    flint::Fmpz bound;
    fmpz_poly_discriminant(discriminant.raw(), polynomial.raw());
    fmpz_zero(bound.raw());

    slong good_primes = 0;
    for (ulong p = 3; good_primes < kTrivialGoodPrimeLimit;) {
        if (fmpz_fdiv_ui(discriminant.raw(), p) == 0) {
            p = n_nextprime(p, 1);
            continue;
        }

        const ulong degree_gcd = factor_residue_degree_gcd(polynomial, p);
        if (degree_gcd == 0) {
            p = n_nextprime(p, 1);
            continue;
        }

        fmpz_ui_pow_ui(local_bound.raw(), p, degree_gcd);
        fmpz_sub_ui(local_bound.raw(), local_bound.raw(), 1);
        if (fmpz_is_zero(bound.raw()) != 0) {
            fmpz_set(bound.raw(), local_bound.raw());
        } else {
            fmpz_gcd(bound.raw(), bound.raw(), local_bound.raw());
        }
        ++good_primes;

        if (flint::fmpz_equal_si(bound, 2)) {
            fmpz_set_ui(order.raw(), 2);
            return generator.set_si(-1);
        }

        p = n_nextprime(p, 1);
    }

    return false;
}

bool compute_roots_of_unity(flint::FmpzRef order,
                            Element& generator,
                            const NumberField& field) noexcept {
    Signature sig;
    if (!sig.compute(field)) {
        return false;
    }

    if (sig.r1() > 0) {
        fmpz_set_ui(order.raw(), 2);
        return generator.set_si(-1);
    }

    flint::Fmpz radicand;
    if (field.quadratic_radicand(flint::FmpzRef(radicand))) {
        if (flint::fmpz_equal_si(radicand, -1)) {
            fmpz_set_ui(order.raw(), 4);
            return generator.gen();
        }
        if (flint::fmpz_equal_si(radicand, -3)) {
            fmpz_set_ui(order.raw(), 6);
            return set_rational_half(generator, 1, 1);
        }
        fmpz_set_ui(order.raw(), 2);
        return generator.set_si(-1);
    }

    if (c_parity_roots_x4_plus_one(generator, field)) {
        fmpz_set_ui(order.raw(), 8);
        return true;
    }
    // reference `nfrootsof1` and reference `_torsion_units_gen` both compute a proven
    // good-prime order multiple before attempting cyclotomic root searches.
    if (roots_trivial_by_good_primes(order, generator, field)) {
        return true;
    }
    if (roots_cyclotomic_square_factors(order, generator, field)) {
        return true;
    }
    return false;
}

}  // namespace

bool roots_of_unity(flint::FmpzRef order,
                    Element& generator,
                    const NumberField& field) noexcept {
    if (!field.is_defined() || !detail::ensure_parent(generator, field)) {
        return false;
    }

    flint::Fmpz tmp_order;
    Element tmp_generator(field);
    if (!tmp_generator.is_defined() ||
        !compute_roots_of_unity(flint::FmpzRef(tmp_order), tmp_generator,
                                field)) {
        return false;
    }

    fmpz_set(order.raw(), tmp_order.raw());
    return generator.set(tmp_generator);
}

bool root_of_unity_order(flint::FmpzRef order,
                         const NumberField& field) noexcept {
    Element generator(field);
    if (!field.is_defined() || !generator.is_defined()) {
        return false;
    }

    flint::Fmpz tmp_order;
    if (!compute_roots_of_unity(flint::FmpzRef(tmp_order), generator,
                                field)) {
        return false;
    }
    fmpz_set(order.raw(), tmp_order.raw());
    return true;
}

bool root_of_unity_generator(Element& generator,
                             const NumberField& field) noexcept {
    if (!field.is_defined() || !detail::ensure_parent(generator, field)) {
        return false;
    }

    flint::Fmpz tmp_order;
    Element tmp_generator(field);
    if (!tmp_generator.is_defined() ||
        !compute_roots_of_unity(flint::FmpzRef(tmp_order), tmp_generator,
                                field)) {
        return false;
    }
    return generator.set(tmp_generator);
}

}  // namespace silex
