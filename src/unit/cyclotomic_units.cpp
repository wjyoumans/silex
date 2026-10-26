#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/number_field.hpp>

#include "unit_internal.hpp"

#include <flint/fmpq.h>
#include <flint/ulong_extras.h>

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

// An integral element (it is a polynomial in the generator with integer
// coefficients) whose norm is +-1 is a unit of the maximal order.
bool has_unit_norm(const Element& value) noexcept {
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
        if (!has_unit_norm(unit)) {
            return false;
        }
    }
    out.swap(units);
    return true;
}

}  // namespace silex::detail
