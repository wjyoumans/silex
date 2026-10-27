#include <silex/element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/unit.hpp>

#include "unit_internal.hpp"

#include <flint/arb.h>
#include <flint/arb_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_lll.h>
#include <flint/fmpz_mat.h>
#include <flint/ulong_extras.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// Cyclotomic units of a field given by a cyclotomic polynomial Phi_m,
// following source nfcyclotomicunits (src/basemath/buch2.c), which
// Buchall_param seeds into the class/unit relation search before it looks
// for any relation (add_cyclotomic_units, buch2.c:3626-3636, called at
// :3873).  See docs/reference/algorithms_and_sources.rst, "Class groups
// and order units".

namespace silex::detail {
namespace {

bool is_prime_power(ulong q) noexcept {
    if (q < 2) {
        return false;
    }
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, q, 1);
    return factors.num == 1;
}

bool push_copy(std::vector<Element>& out, const Element& value) noexcept {
    Element copy(*value.parent());
    if (!copy.is_defined() || !copy.set(value)) {
        return false;
    }
    out.push_back(std::move(copy));
    return true;
}

// An integral element whose norm is +-1 is a unit of the maximal order.
// Integrality is checked exactly: the element is a polynomial in the
// generator with integer coefficients, and the generator is integral
// because the defining polynomial is monic with integer coefficients (the
// callers only accept Phi_m and psi_n).
bool is_integral_unit(const Element& value) noexcept {
    flint::FmpqPoly polynomial;
    if (!value.get_fmpq_poly(flint::FmpqPolyRef(polynomial)) ||
        fmpz_is_one(fmpq_poly_denref(polynomial.raw())) == 0) {
        return false;
    }
    flint::Fmpq norm;
    if (!value.norm(flint::FmpqRef(norm))) {
        return false;
    }
    return fmpz_is_one(fmpq_denref(norm.raw())) != 0 &&
           fmpz_is_pm1(fmpq_numref(norm.raw())) != 0;
}

}  // namespace

bool cyclotomic_field_units(std::vector<Element>& out,
                            const NumberField& field) noexcept {
    out.clear();
    const ulong m = cyclotomic_defining_polynomial_index(field);
    if (m == 0) {
        return true;
    }
    // The field has w = lcm(2, m) roots of unity.  Source: no units when
    // w <= 6; otherwise z is a primitive n-th root of unity with n = w,
    // or n = w / 2 when w = 2 (mod 4), so that n != 2 (mod 4).  The
    // generator x is a primitive m-th root: n = m and z = x unless
    // m = 2 (mod 4), where n = m / 2 and z = -x.
    const ulong w = (m % 2 == 0) ? m : 2 * m;
    if (w <= 6) {
        return true;
    }
    const bool negate = m % 4 == 2;
    const ulong n = negate ? m / 2 : m;
    const ulong half = n / 2;

    // powers[i] = z^i for 0 <= i < n / 2.
    std::vector<Element> powers;
    powers.reserve(half);
    {
        Element one(field);
        if (!one.is_defined() || !one.one()) {
            return false;
        }
        powers.push_back(std::move(one));
    }
    Element z(field);
    if (!z.is_defined() || !z.gen()) {
        return false;
    }
    if (negate) {
        Element negated(field);
        if (!negated.is_defined() || !negated.negate(z)) {
            return false;
        }
        z.swap(negated);
    }
    for (ulong i = 1; i < half; ++i) {
        Element next(field);
        if (!next.is_defined() || !next.multiply(powers.back(), z)) {
            return false;
        }
        powers.push_back(std::move(next));
    }

    std::vector<Element> units;
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, n, 1);
    // Second kind: (Z^a - 1) / (Z - 1) = 1 + Z + ... + Z^(a - 1) for each
    // p^k || n, Z = z^(n / p^k), 1 < a <= (p^k - 1) / 2, p not dividing a.
    for (slong f = 0; f < factors.num; ++f) {
        const ulong p = factors.p[f];
        const ulong pk = n_pow(p, factors.exp[f]);
        const ulong step = n / pk;
        const ulong top = (pk - 1) / 2;
        Element partial(field);
        if (!partial.is_defined() || !partial.one()) {
            return false;
        }
        for (ulong a = 2; a <= top; ++a) {
            const ulong exponent = step * (a - 1);
            if (exponent >= powers.size()) {
                return false;
            }
            Element next(field);
            if (!next.is_defined() ||
                !next.add(partial, powers[exponent])) {
                return false;
            }
            partial.swap(next);
            if (a % p != 0 && !push_copy(units, partial)) {
                return false;
            }
        }
    }
    // First kind, when n is not a prime power: z^a - 1 for 1 <= a < n / 2,
    // skipping a > 1 that divides n or has n / gcd(a, n) a prime power.
    if (factors.num > 1) {
        for (ulong a = 1; a < half; ++a) {
            if (a > 1 && (n % a == 0 || is_prime_power(n / n_gcd(a, n)))) {
                continue;
            }
            Element value(field);
            if (!value.is_defined() || !value.add_si(powers[a], -1)) {
                return false;
            }
            units.push_back(std::move(value));
        }
    }

    // Each listed element is checked exactly before use.
    for (const Element& unit : units) {
        if (!is_integral_unit(unit)) {
            return false;
        }
    }
    out.swap(units);
    return true;
}

namespace {

// values[k] = T_k(y) for 0 <= k < count, where T_0 = t0, T_1 = y and
// T_(k+1) = y T_k - T_(k-1).  With t0 = 2 these are the Dickson polynomials
// D_k(y), D_k(z + 1/z) = z^k + z^-k; with t0 = 1 they are S_k(y) = U_k(y / 2)
// (U_k the Chebyshev polynomial of the second kind),
// S_k(z + 1/z) = (z^(k+1) - z^-(k+1)) / (z - z^-1).
bool recurrence_values(std::vector<Element>& values, const Element& y,
                       slong t0, ulong count) noexcept {
    values.clear();
    values.reserve(count);
    for (ulong k = 0; k < count; ++k) {
        Element value(*y.parent());
        if (!value.is_defined()) {
            return false;
        }
        if (k == 0) {
            if (!value.set_si(t0)) {
                return false;
            }
        } else if (k == 1) {
            if (!value.set(y)) {
                return false;
            }
        } else {
            Element product(*y.parent());
            if (!product.is_defined() ||
                !product.multiply(y, values[k - 1]) ||
                !value.subtract(product, values[k - 2])) {
                return false;
            }
        }
        values.push_back(std::move(value));
    }
    return true;
}

// Reference _cyclotomic_units_totally_real_prime_conductor
// (src/NumField/NfAbs/Cyclotomic.jl:47-87, after Schoof, "Class numbers of
// real cyclotomic fields of prime conductor", Math. Comp. 72 (2003), Section
// 2): for p prime and g a primitive root modulo p, eta = U_(g-1)(y / 2) =
// S_(g-1)(y) and its images sigma_a(eta) = S_(ga-1)(y) / S_(a-1)(y) for
// 2 <= a < p / 2 - 1; (p - 3) / 2 elements, the unit rank, after the reference's
// leading -1.  S_(ga-1)(y) is replaced by S_(c-1)(y), c = ga mod p, which
// has the same value at every root of psi_p because S_(k-1)(2 cos t) =
// sin(k t) / sin t and t is a multiple of 2 pi / p.  The reference takes the
// generator of (Z/pZ)^* that its unit_group returns; this takes FLINT's
// n_primitive_root_prime.
bool prime_conductor_units(std::vector<Element>& out, const Element& y,
                           ulong p) noexcept {
    std::vector<Element> s;
    if (!recurrence_values(s, y, 1, p - 1)) {
        return false;
    }
    const ulong g = n_primitive_root_prime(p);
    if (g < 2 || g >= p) {
        return false;
    }
    Element eta(*y.parent());
    if (!eta.is_defined() || !eta.set(s[g - 1])) {
        return false;
    }
    out.push_back(std::move(eta));
    for (ulong a = 2; 2 * a + 2 < p; ++a) {
        const ulong c = n_mulmod2(g, a, p);
        Element inverse(*y.parent());
        Element image(*y.parent());
        if (c == 0 || !inverse.is_defined() || !image.is_defined() ||
            !inverse.invert(s[a - 1]) ||
            !image.multiply(s[c - 1], inverse)) {
            return false;
        }
        out.push_back(std::move(image));
    }
    return true;
}

// Reference _cyclotomic_units_totally_real_prime_power_conductor
// (src/NumField/NfAbs/Cyclotomic.jl:89-136, after Washington, Introduction
// to Cyclotomic Fields, Lemma 8.1): for q a prime power and z = zeta_q,
// xi_a = z^((1 - a) / 2) (1 - z^a) / (1 - z) for 1 < a < q / 2 prime to q,
// phi(q) / 2 - 1 elements.  The reference forms xi_a in Q(zeta_q) and maps
// it back; here it is written in y = z + 1/z directly.  Up to sign xi_a is
// sin(pi a / q) / sin(pi / q), and with a' the odd one of a and q - a,
// a' = 2 m + 1, that is sum_(|j| <= m) z^j = 1 + D_1(y) + ... + D_m(y).
// The sign is a torsion factor and does not change the unit group.
bool prime_power_conductor_units(std::vector<Element>& out,
                                 const Element& y, ulong q) noexcept {
    std::vector<Element> d;
    if (!recurrence_values(d, y, 2, q / 2 + 1)) {
        return false;
    }
    // partial[m] = 1 + D_1 + ... + D_m.
    std::vector<Element> partial;
    partial.reserve(d.size());
    for (std::size_t m = 0; m < d.size(); ++m) {
        Element value(*y.parent());
        if (!value.is_defined() ||
            (m == 0 ? !value.one() : !value.add(partial[m - 1], d[m]))) {
            return false;
        }
        partial.push_back(std::move(value));
    }
    for (ulong a = 2; 2 * a < q; ++a) {
        if (n_gcd(a, q) != 1) {
            continue;
        }
        const ulong odd = a % 2 == 1 ? a : q - a;
        const ulong m = (odd - 1) / 2;
        if (m >= partial.size() || !push_copy(out, partial[m])) {
            return false;
        }
    }
    return true;
}

// The real elements of the group D generated by +-z and the 1 - z^a,
// z = zeta_n, for n = 0 mod 4 or n odd, written in y = z + 1/z.  Neither
// the source nor the reference has a routine for composite n; these generators are
// Silex's own, from the identities below.
// For n = 4 m: D_j(y) = |1 - z^(2(m - j))| for 0 <= j < m, and
// D_((a-1)/2)(y) - D_((a+1)/2)(y) = |1 - z| |1 - z^a| for odd a < 2 m.  For
// n odd: S_(b-1)(y) = |1 - z^(2b)| / |1 - z^2| for 2 <= b <= (n - 1) / 2, and
// 2 - D_2(y) = |1 - z^2|^2.
bool composite_conductor_generators(std::vector<Element>& out,
                                    const Element& y, ulong n) noexcept {
    if (n % 4 == 0) {
        const ulong m = n / 4;
        std::vector<Element> d;
        if (!recurrence_values(d, y, 2, m + 1)) {
            return false;
        }
        for (ulong j = 0; j < m; ++j) {
            if (!push_copy(out, d[j])) {
                return false;
            }
        }
        for (ulong a = 1; a < 2 * m; a += 2) {
            Element difference(*y.parent());
            if (!difference.is_defined() ||
                !difference.subtract(d[(a - 1) / 2], d[(a + 1) / 2])) {
                return false;
            }
            out.push_back(std::move(difference));
        }
        return true;
    }
    std::vector<Element> s;
    if (!recurrence_values(s, y, 1, (n - 1) / 2)) {
        return false;
    }
    for (std::size_t k = 1; k < s.size(); ++k) {
        if (!push_copy(out, s[k])) {
            return false;
        }
    }
    // 2 - D_2(y) = 4 - y^2.
    Element square(*y.parent());
    Element negated(*y.parent());
    Element value(*y.parent());
    if (!square.is_defined() || !negated.is_defined() ||
        !value.is_defined() || !square.multiply(y, y) ||
        !negated.negate(square) || !value.add_si(negated, 4)) {
        return false;
    }
    out.push_back(std::move(value));
    return true;
}

// The units among products of `generators`: the integer kernel of the
// exponents of the primes dividing n in their norms (every generator norm
// must be supported there), LLL-reduced, each kernel vector giving the
// product of the generators to its exponents.  The norm exponents stand in
// for the valuations at the primes above each p | n; every product is then
// checked exactly to be an integral unit, so a product that is not one makes
// the construction fail rather than enter the list.
bool valuation_kernel_units(std::vector<Element>& out,
                            const std::vector<Element>& generators,
                            ulong n) noexcept {
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, n, 1);
    const slong count = static_cast<slong>(generators.size());
    const slong primes = factors.num;
    flint::FmpzMat valuations(count, primes);
    for (slong i = 0; i < count; ++i) {
        flint::Fmpq norm;
        if (!generators[static_cast<std::size_t>(i)].norm(
                    flint::FmpqRef(norm)) ||
            fmpz_is_one(fmpq_denref(norm.raw())) == 0 ||
            fmpz_is_zero(fmpq_numref(norm.raw())) != 0) {
            return false;
        }
        flint::Fmpz rest;
        fmpz_set(rest.raw(), fmpq_numref(norm.raw()));
        for (slong j = 0; j < primes; ++j) {
            flint::Fmpz prime;
            fmpz_set_ui(prime.raw(), factors.p[j]);
            const slong exponent = fmpz_remove(rest.raw(), rest.raw(),
                                               prime.raw());
            fmpz_set_si(fmpz_mat_entry(valuations.raw(), i, j), exponent);
        }
        if (fmpz_is_pm1(rest.raw()) == 0) {
            return false;
        }
    }
    // H = U V with U unimodular: the rows of U whose row of H is zero are a
    // basis of the integer left kernel of V.
    flint::FmpzMat hnf(count, primes);
    flint::FmpzMat transform(count, count);
    ::fmpz_mat_hnf_transform(hnf.raw(), transform.raw(), valuations.raw());
    slong kernel_rows = 0;
    for (slong i = 0; i < count; ++i) {
        if (fmpz_mat_is_zero_row(hnf.raw(), i) != 0) {
            ++kernel_rows;
        }
    }
    if (kernel_rows == 0) {
        return true;
    }
    flint::FmpzMat kernel(kernel_rows, count);
    for (slong i = 0, row = 0; i < count; ++i) {
        if (fmpz_mat_is_zero_row(hnf.raw(), i) == 0) {
            continue;
        }
        for (slong j = 0; j < count; ++j) {
            fmpz_set(fmpz_mat_entry(kernel.raw(), row, j),
                     fmpz_mat_entry(transform.raw(), i, j));
        }
        ++row;
    }
    flint::FmpzLll lll;
    ::fmpz_lll(kernel.raw(), nullptr, lll.raw());

    const NumberField& field = *generators.front().parent();
    for (slong i = 0; i < kernel_rows; ++i) {
        Element numerator(field);
        Element denominator(field);
        if (!numerator.is_defined() || !denominator.is_defined() ||
            !numerator.one() || !denominator.one()) {
            return false;
        }
        for (slong j = 0; j < count; ++j) {
            const fmpz* exponent = fmpz_mat_entry(kernel.raw(), i, j);
            if (fmpz_is_zero(exponent) != 0) {
                continue;
            }
            flint::Fmpz magnitude;
            fmpz_abs(magnitude.raw(), exponent);
            Element power(field);
            Element product(field);
            if (!power.is_defined() || !product.is_defined() ||
                !power.pow_fmpz(generators[static_cast<std::size_t>(j)],
                                flint::FmpzConstRef(magnitude))) {
                return false;
            }
            Element& side = fmpz_sgn(exponent) > 0 ? numerator : denominator;
            if (!product.multiply(side, power)) {
                return false;
            }
            side.swap(product);
        }
        Element inverse(field);
        Element unit(field);
        if (!inverse.is_defined() || !unit.is_defined() ||
            !inverse.invert(denominator) ||
            !unit.multiply(numerator, inverse) ||
            !is_integral_unit(unit)) {
            return false;
        }
        out.push_back(std::move(unit));
    }
    return true;
}

// Keeps `rank` candidates proved independent by the test of
// units_independent: the Gram determinant of the product-mode logarithmic
// embeddings of the kept units is a positive ball.  When there are more
// candidates than `rank`, a floating-point proposal picks them and the ball
// test decides; the in-order per-row test is the fallback.  The embeddings
// of all candidates are computed once.
bool select_independent_units(std::vector<Element>& out,
                              const std::vector<Element>& candidates,
                              EmbeddingContext& embeddings, slong rank,
                              slong precision) noexcept {
    out.clear();
    const slong count = static_cast<slong>(candidates.size());
    // psi_n is totally real: deg = rank + 1 real places.
    const slong places = rank + 1;
    flint::ArbMat logs(count, places);
    if (!unit_log_matrix(flint::ArbMatRef(logs), embeddings,
                         ElementSpan(candidates.data(), candidates.size()),
                         LogEmbeddingMode::product, precision)) {
        return true;
    }
    // Whether the rows `kept` of `logs` have a positive Gram determinant.
    auto positive_gram = [&](const std::vector<slong>& kept) noexcept {
        const slong size = static_cast<slong>(kept.size());
        flint::ArbMat rows(size, places);
        flint::ArbMat transpose(places, size);
        flint::ArbMat gram(size, size);
        flint::Arb determinant;
        for (slong a = 0; a < size; ++a) {
            for (slong b = 0; b < places; ++b) {
                arb_set(arb_mat_entry(rows.raw(), a, b),
                        arb_mat_entry(logs.raw(),
                                      kept[static_cast<std::size_t>(a)],
                                      b));
            }
        }
        arb_mat_transpose(transpose.raw(), rows.raw());
        arb_mat_mul(gram.raw(), rows.raw(), transpose.raw(), precision);
        arb_mat_det(determinant.raw(), gram.raw(), precision);
        return flint::arb_is_positive(determinant);
    };
    std::vector<slong> kept;
    if (count == rank) {
        // All of them at once (the prime and prime-power lists).
        for (slong i = 0; i < count; ++i) {
            kept.push_back(i);
        }
        if (!positive_gram(kept)) {
            kept.clear();
        }
    } else {
        // Gram-Schmidt in double precision on the ball midpoints proposes
        // `rank` rows (a row whose residual is below 1e-6 of its length is
        // passed over); one ball test then decides them.  Only when it is
        // not positive does each row get its own test, in order.
        std::vector<double> basis;
        for (slong i = 0;
             i < count && static_cast<slong>(kept.size()) < rank; ++i) {
            std::vector<double> row(static_cast<std::size_t>(places));
            double length = 0.0;
            for (slong b = 0; b < places; ++b) {
                row[static_cast<std::size_t>(b)] = arf_get_d(
                        arb_midref(arb_mat_entry(logs.raw(), i, b)),
                        ARF_RND_NEAR);
                length += row[static_cast<std::size_t>(b)] *
                          row[static_cast<std::size_t>(b)];
            }
            for (std::size_t k = 0; k < kept.size(); ++k) {
                const double* q = basis.data() + k * row.size();
                double dot = 0.0;
                for (std::size_t b = 0; b < row.size(); ++b) {
                    dot += q[b] * row[b];
                }
                for (std::size_t b = 0; b < row.size(); ++b) {
                    row[b] -= dot * q[b];
                }
            }
            double residual = 0.0;
            for (double value : row) {
                residual += value * value;
            }
            if (!(residual > 1e-12 * length) || !(residual > 0.0)) {
                continue;
            }
            const double scale = 1.0 / std::sqrt(residual);
            for (double value : row) {
                basis.push_back(value * scale);
            }
            kept.push_back(i);
        }
        if (static_cast<slong>(kept.size()) != rank || !positive_gram(kept)) {
            kept.clear();
            for (slong i = 0; i < count; ++i) {
                if (static_cast<slong>(kept.size()) == rank) {
                    break;
                }
                kept.push_back(i);
                if (!positive_gram(kept)) {
                    kept.pop_back();
                }
            }
        }
    }
    for (slong i : kept) {
        if (!push_copy(out, candidates[static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    return true;
}

}  // namespace

bool real_cyclotomic_field_units(std::vector<Element>& out,
                                 const NumberField& field,
                                 EmbeddingContext& embeddings) noexcept {
    out.clear();
    const ulong index = real_cyclotomic_defining_polynomial_index(field);
    if (index == 0) {
        return true;
    }
    // psi_n has degree phi(n) / 2 and only real roots, so the unit rank is
    // the degree minus one.
    const slong rank = field.degree() - 1;
    if (rank < 1) {
        return true;
    }
    // y = x is a root of psi_n.  For n = 2 m, m odd, the field is also
    // Q(zeta_m)^+ and y = -x = 2 cos(pi - pi / m) = 2 cos(2 pi k / m),
    // k = (m - 1) / 2 prime to m, is a root of psi_m.
    const bool negate = index % 4 == 2;
    const ulong n = negate ? index / 2 : index;
    Element y(field);
    if (!y.is_defined() || !y.gen()) {
        return false;
    }
    if (negate) {
        Element negated(field);
        if (!negated.is_defined() || !negated.negate(y)) {
            return false;
        }
        y.swap(negated);
    }

    std::vector<Element> candidates;
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, n, 1);
    bool constructed = false;
    if (factors.num == 1 && factors.exp[0] == 1) {
        constructed = prime_conductor_units(candidates, y, n);
    } else if (factors.num == 1) {
        constructed = prime_power_conductor_units(candidates, y, n);
    } else {
        std::vector<Element> generators;
        constructed =
                composite_conductor_generators(generators, y, n) &&
                !generators.empty() &&
                valuation_kernel_units(candidates, generators, n);
    }
    if (!constructed) {
        return false;
    }
    // Every candidate is checked exactly to be an integral unit.
    for (const Element& unit : candidates) {
        if (!is_integral_unit(unit)) {
            return false;
        }
    }
    // Exactly `rank` independent units, or none.
    std::vector<Element> units;
    for (slong precision = 128; precision <= 512; precision *= 4) {
        if (!select_independent_units(units, candidates, embeddings, rank,
                                      precision)) {
            return false;
        }
        if (static_cast<slong>(units.size()) == rank) {
            out.swap(units);
            return true;
        }
    }
    return false;
}

}  // namespace silex::detail
