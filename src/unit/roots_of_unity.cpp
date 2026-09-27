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

#include "unit_internal.hpp"

namespace silex {
namespace {

bool set_rational_half(Element& out, slong constant, slong linear) noexcept {
    flint::FmpqPoly polynomial;
    flint::Fmpq coeff;

    fmpq_set_si(coeff.raw(), constant, 2);
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 0, coeff.raw());
    fmpq_set_si(coeff.raw(), linear, 2);
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, coeff.raw());
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
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

// Monic integral defining polynomial of theta' = scale * theta, where theta
// is the field generator.  With T = sum a_i x^i the integral numerator of the
// defining polynomial and a_n its leading coefficient, theta' = a_n theta is
// a root of sum a_i a_n^(n-1-i) y^i, which is monic with integer
// coefficients.  scale = 1 when T is already monic.
bool monic_integral_defining_poly(flint::FmpzPoly& out,
                                  flint::Fmpz& scale,
                                  const NumberField& field) noexcept {
    const nf_struct* raw = field.raw_flint_field();
    if (raw == nullptr || fmpq_poly_degree(raw->pol) < 1) {
        return false;
    }

    flint::FmpzPoly numerator;
    fmpq_poly_get_numerator(numerator.raw(), raw->pol);
    const slong degree = fmpz_poly_degree(numerator.raw());
    fmpz_poly_get_coeff_fmpz(scale.raw(), numerator.raw(), degree);
    if (fmpz_is_zero(scale.raw()) != 0) {
        return false;
    }

    flint::Fmpz power;
    flint::Fmpz coeff;
    fmpz_one(power.raw());
    fmpz_poly_zero(out.raw());
    fmpz_poly_set_coeff_ui(out.raw(), degree, 1);
    for (slong i = degree - 1; i >= 0; --i) {
        fmpz_poly_get_coeff_fmpz(coeff.raw(), numerator.raw(), i);
        fmpz_mul(coeff.raw(), coeff.raw(), power.raw());
        fmpz_poly_set_coeff_fmpz(out.raw(), i, coeff.raw());
        fmpz_mul(power.raw(), power.raw(), scale.raw());
    }
    return true;
}

// Resource bound only: the search below stops by the stationarity rule long
// before this in practice.  Reaching it leaves the caller without a bound,
// which fails closed.
constexpr slong kGoodPrimeSearchLimit = slong{1} << 16;

// Proven multiple of the number w of roots of unity, from good primes.
//
// Source trace: reference 2.17.3 src/basemath/nffactor.c `guess_roots`.  For a
// prime p not dividing disc(T), every prime P | p has N(P) = p^f(P) and
// mu(K) injects into (O_K/P)^*, so w | p^f(P) - 1 for every P | p, hence
// w | gcd_P (p^f(P) - 1) = p^(gcd f(P)) - 1.  The gcd over primes p >= 3 is a
// proven multiple of w.  As in the reference, the loop stops once the gcd fits in a
// word and has been unchanged for more than B = n + 20 consecutive good
// primes.  Following reference
// `TorsionUnits.jl:_torsion_group_order_divisor`, it also stops as soon as the gcd is
// 2, which is exact since w is even and every odd p^f - 1 is even.
bool good_prime_root_bound(ulong& bound,
                           const flint::FmpzPoly& polynomial,
                           const flint::Fmpz& discriminant) noexcept {
    const slong degree = fmpz_poly_degree(polynomial.raw());
    const slong stable_limit = degree + 20;

    flint::Fmpz gcd;
    flint::Fmpz local;
    flint::Fmpz previous;
    bool have = false;
    slong stable = 0;
    ulong p = 3;
    for (slong attempt = 0; attempt < kGoodPrimeSearchLimit;
         ++attempt, p = n_nextprime(p, 1)) {
        if (fmpz_fdiv_ui(discriminant.raw(), p) == 0) {
            continue;
        }
        const ulong degree_gcd = factor_residue_degree_gcd(polynomial, p);
        if (degree_gcd == 0) {
            continue;
        }

        fmpz_ui_pow_ui(local.raw(), p, degree_gcd);
        fmpz_sub_ui(local.raw(), local.raw(), 1);
        fmpz_set(previous.raw(), gcd.raw());
        if (have) {
            fmpz_gcd(gcd.raw(), gcd.raw(), local.raw());
        } else {
            fmpz_set(gcd.raw(), local.raw());
        }

        if (flint::fmpz_equal_si(gcd, 2)) {
            bound = 2;
            return true;
        }

        if (have && fmpz_equal(previous.raw(), gcd.raw()) != 0) {
            if (fmpz_abs_fits_ui(gcd.raw()) != 0 && ++stable > stable_limit) {
                bound = fmpz_get_ui(gcd.raw());
                return true;
            }
        } else {
            stable = 0;
        }
        have = true;
    }
    return false;
}

// Degree and ramification reduction of the good-prime bound.
//
// Source trace: reference 2.17.3 src/basemath/nffactor.c `nfrootsof1`, step 1
// after `guess_roots`.  Q(zeta_(p^k)) has degree (p-1) p^(k-1) and
// v_p(disc) = ((p-1)k - 1) p^(k-1), so Q(zeta_(p^k)) subset K needs
// (p-1) | n and v_p(disc_K) >= kn - n/(p-1) (for p = 2: v_2 >= n(k-1)), and
// p^(k-1) | n/(p-1).  The reference polynomial branch uses disc(T) for disc_K;
// v_p(disc(T)) >= v_p(disc_K), so the test is weaker but the reduced bound is
// still a proven multiple of w.
ulong reduce_root_bound(ulong bound,
                        slong degree,
                        const flint::Fmpz& discriminant) noexcept {
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, bound, 1);

    flint::Fmpz prime;
    flint::Fmpz cofactor;
    for (int i = 0; i < factors.num; ++i) {
        const ulong p = factors.p[i];
        const slong e = static_cast<slong>(factors.exp[i]);
        fmpz_set_ui(prime.raw(), p);
        const slong v = static_cast<slong>(
                fmpz_remove(cofactor.raw(), discriminant.raw(), prime.raw()));
        if (p == 2) {
            if (e == 1) {
                continue;
            }
            const slong vnf = static_cast<slong>(
                    flint_ctz(static_cast<ulong>(degree)));
            for (slong k = e < vnf + 1 ? e : vnf + 1; k >= 1; --k) {
                if (v >= degree * (k - 1)) {
                    bound >>= e - k;
                    break;
                }
            }
            continue;
        }

        if (static_cast<ulong>(degree) % (p - 1) != 0) {
            bound /= n_pow(p, static_cast<ulong>(e));
            continue;
        }
        const slong q = degree / static_cast<slong>(p - 1);
        ulong q_rest = static_cast<ulong>(q);
        const slong vnf = static_cast<slong>(n_remove(&q_rest, p));
        for (slong k = e < vnf + 1 ? e : vnf + 1; k >= 0; --k) {
            if (v >= degree * k - q) {
                bound /= n_pow(p, static_cast<ulong>(e - k));
                break;
            }
        }
    }
    return bound;
}

// True when z has exact multiplicative order m: z^m = 1 and z^(m/l) != 1 for
// every prime l | m.  This is the certificate that w is at least m.
bool has_exact_order(const Element& z, ulong m) noexcept {
    const NumberField* parent = z.parent();
    if (parent == nullptr || m == 0) {
        return false;
    }
    Element power(*parent);
    flint::Fmpz exponent;
    fmpz_set_ui(exponent.raw(), m);
    if (!power.is_defined() ||
        !power.pow_fmpz(z, flint::FmpzConstRef(exponent)) ||
        !power.equal_si(1)) {
        return false;
    }

    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, m, 1);
    for (int i = 0; i < factors.num; ++i) {
        fmpz_set_ui(exponent.raw(), m / factors.p[i]);
        if (!power.pow_fmpz(z, flint::FmpzConstRef(exponent)) ||
            power.equal_si(1)) {
            return false;
        }
    }
    return true;
}

// Candidate root of unity when T(y) = Phi_N(+-y + c) for the scaled
// generator y = theta'.
//
// Source trace: reference 2.17.3 src/basemath/nffactor.c
// `ZXirred_is_cyclo_translate` (step 1.5 of `nfrootsof1`, used when
// phi(bound) = n).  The coefficient of y^(n-1) in Phi_N(+-y + c) is
// +-(n c - mu(N)), so c is the quotient of that coefficient by n, with
// remainder +-1 when N is squarefree and 0 otherwise.  The reference then confirms the
// translate with a Graeffe comparison and returns +-(theta' + c).  Here the
// two candidates +-(theta' + c) are instead confirmed by has_exact_order, an
// exact certificate that does not depend on the guess.
bool cyclotomic_translate_root(Element& generator,
                               const NumberField& field,
                               const flint::FmpzPoly& polynomial,
                               const flint::Fmpz& scale,
                               ulong bound) noexcept {
    const slong degree = fmpz_poly_degree(polynomial.raw());
    flint::Fmpz c;
    flint::Fmpz r;
    flint::Fmpz trace_coeff;
    flint::Fmpz degree_fmpz;
    fmpz_poly_get_coeff_fmpz(trace_coeff.raw(), polynomial.raw(), degree - 1);
    fmpz_set_si(degree_fmpz.raw(), degree);
    fmpz_tdiv_qr(c.raw(), r.raw(), trace_coeff.raw(), degree_fmpz.raw());

    if (n_is_squarefree(bound) == 0) {
        if (fmpz_is_zero(r.raw()) == 0) {
            return false;
        }
    } else {
        if (fmpz_cmp_si(r.raw(), -1) < 0) {
            fmpz_add_si(r.raw(), r.raw(), degree);
            fmpz_sub_ui(c.raw(), c.raw(), 1);
        } else if (fmpz_equal_si(r.raw(), degree - 1) != 0) {
            fmpz_set_si(r.raw(), -1);
            fmpz_add_ui(c.raw(), c.raw(), 1);
        }
        if (!fmpz_is_pm1(r.raw())) {
            return false;
        }
    }

    Element candidate(field);
    Element term(field);
    if (!candidate.is_defined() || !term.is_defined() ||
        !candidate.gen() ||
        !term.set_fmpz(flint::FmpzConstRef(scale)) ||
        !candidate.multiply(candidate, term) ||
        !term.set_fmpz(flint::FmpzConstRef(c)) ||
        !candidate.add(candidate, term)) {
        return false;
    }
    if (has_exact_order(candidate, bound)) {
        return generator.set(candidate);
    }
    if (candidate.negate(candidate) && has_exact_order(candidate, bound)) {
        return generator.set(candidate);
    }
    return false;
}

// A primitive p^k-th root of unity for p = 2 (k >= 2) or p = 3 (k >= 1),
// built with the exact square- and power-root primitives: zeta_4 = sqrt(-1)
// and zeta_(2^(j+1)) = sqrt(zeta_(2^j)); zeta_3 = (-1 + sqrt(-3)) / 2 and
// zeta_(3^(j+1)) = zeta_(3^j)^(1/3).  Any p-th root of a primitive
// p^j-th root of unity (j >= 1) is a primitive p^(j+1)-th root, and if
// zeta_(p^(j+1)) is in K then every p^j-th root of unity is a p-th power in K,
// so the choice of root at each step does not matter.  This replaces the
// `nfisincl(polcyclo(p^k), T)` test of reference `nfrootsof1` step 2 (and
// the `_roots_hensel` search of reference `_torsion_units_gen`) for these primes;
// false when a root is not found or not decided.
bool prime_power_root_of_unity(Element& out,
                               const NumberField& field,
                               ulong p,
                               slong k) noexcept {
    Element z(field);
    Element next(field);
    if (!z.is_defined() || !next.is_defined()) {
        return false;
    }
    flint::Fmpz exponent;
    fmpz_set_ui(exponent.raw(), p);
    slong have = 0;
    if (p == 2 && k >= 2) {
        if (!primitive_fourth_root(z, field)) {
            return false;
        }
        have = 2;
    } else if (p == 3 && k >= 1) {
        if (!primitive_sixth_root(z, field) || !z.multiply(z, z)) {
            return false;
        }
        have = 1;
    } else {
        return false;
    }

    for (; have < k; ++have) {
        bool found = false;
        if (!z.is_power(found, next, flint::FmpzConstRef(exponent)) ||
            !found) {
            return false;
        }
        z.swap(next);
    }
    return out.set(z);
}

// The good-prime bound of `guess_roots` reduced by degree and ramification:
// a proven multiple of w for the field defined by the monic integral
// `polynomial`.  False when the discriminant vanishes or the good-prime
// search reaches its resource limit.
bool proven_root_bound(ulong& bound,
                       const flint::FmpzPoly& polynomial) noexcept {
    flint::Fmpz discriminant;
    fmpz_poly_discriminant(discriminant.raw(), polynomial.raw());
    if (fmpz_is_zero(discriminant.raw()) != 0) {
        return false;
    }
    ulong good_prime_bound = 0;
    if (!good_prime_root_bound(good_prime_bound, polynomial, discriminant)) {
        return false;
    }
    if (good_prime_bound != 2) {
        good_prime_bound = reduce_root_bound(
                good_prime_bound, fmpz_poly_degree(polynomial.raw()),
                discriminant);
    }
    bound = good_prime_bound;
    return true;
}

// Step 2 of reference `nfrootsof1` against a proven multiple `bound` of w,
// with a fail-closed rule: the search succeeds only when it certifies a root
// of unity of exact order `bound`, which with the upper bound proves
// w = bound.  The reference instead accepts a smaller p-power after a "wrong
// guess" warning; Silex fails, because w is then not certified by the bound.
// This also covers the unported case p >= 5 outside cyclotomic-translate
// presentations.
bool search_roots_for_bound(flint::FmpzRef order,
                            Element& generator,
                            const NumberField& field,
                            const flint::FmpzPoly& polynomial,
                            const flint::Fmpz& scale,
                            ulong bound) noexcept {
    if (bound < 2 || bound % 2 != 0) {
        return false;
    }
    Element result(field);
    if (!result.is_defined()) {
        return false;
    }
    if (bound == 2) {
        fmpz_set_ui(order.raw(), 2);
        return result.set_si(-1) && generator.set(result);
    }

    const slong degree = fmpz_poly_degree(polynomial.raw());
    if (degree >= 1 &&
        n_euler_phi(bound) == static_cast<ulong>(degree) &&
        cyclotomic_translate_root(result, field, polynomial, scale, bound)) {
        fmpz_set_ui(order.raw(), bound);
        return generator.set(result);
    }

    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, bound, 1);
    if (!result.set_si(1)) {
        return false;
    }
    Element component(field);
    if (!component.is_defined()) {
        return false;
    }
    for (int i = 0; i < factors.num; ++i) {
        const ulong p = factors.p[i];
        const slong e = static_cast<slong>(factors.exp[i]);
        if (p == 2 && e == 1) {
            if (!component.set_si(-1)) {
                return false;
            }
        } else if (!prime_power_root_of_unity(component, field, p, e)) {
            return false;
        }
        if (!result.multiply(result, component)) {
            return false;
        }
    }
    if (!has_exact_order(result, bound)) {
        return false;
    }
    fmpz_set_ui(order.raw(), bound);
    return generator.set(result);
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

    // Source trace: reference 2.17.3 src/basemath/nffactor.c `nfrootsof1`
    // (polynomial branch) and reference `_torsion_units_gen`: a proven good-prime
    // multiple of w, reduced by degree and ramification, then a search for a
    // root of unity of that order.
    flint::FmpzPoly polynomial;
    flint::Fmpz scale;
    if (!monic_integral_defining_poly(polynomial, scale, field)) {
        return false;
    }
    ulong bound = 0;
    if (!proven_root_bound(bound, polynomial)) {
        return false;
    }
    return search_roots_for_bound(order, generator, field, polynomial, scale,
                                  bound);
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

    // Publish the generator first: it is the only write that can fail, so a
    // failure leaves both outputs unchanged.
    if (!generator.set(tmp_generator)) {
        return false;
    }
    fmpz_set(order.raw(), tmp_order.raw());
    return true;
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


namespace detail {

bool roots_of_unity_for_bound(flint::FmpzRef order,
                              Element& generator,
                              const NumberField& field,
                              ulong bound) noexcept {
    if (!field.is_defined() || !detail::ensure_parent(generator, field)) {
        return false;
    }
    flint::FmpzPoly polynomial;
    flint::Fmpz scale;
    if (!monic_integral_defining_poly(polynomial, scale, field)) {
        return false;
    }

    // Guard the precondition: `bound` must be a multiple of the proven
    // reduced good-prime bound, hence of w.  Without it a proper divisor of w
    // (bound 6 for Q(zeta_9)) would pass the exact-order check and be
    // published as w.
    ulong proven = 0;
    if (bound == 0 || !proven_root_bound(proven, polynomial) ||
        bound % proven != 0) {
        return false;
    }

    flint::Fmpz tmp_order;
    Element tmp_generator(field);
    if (!tmp_generator.is_defined() ||
        !search_roots_for_bound(flint::FmpzRef(tmp_order), tmp_generator,
                                field, polynomial, scale, bound)) {
        return false;
    }
    if (!generator.set(tmp_generator)) {
        return false;
    }
    fmpz_set(order.raw(), tmp_order.raw());
    return true;
}

}  // namespace detail

}  // namespace silex
