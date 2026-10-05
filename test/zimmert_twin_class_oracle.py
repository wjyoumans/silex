#!/usr/bin/env python3
"""Independent mpmath oracle for Zimmert's twin-class generation bound.

H. Zimmert, "Ideale kleiner Norm in Idealklassen und eine
Regulatorabschaetzung", Invent. Math. 62 (1981), 367-380.  Korollar 1
(p. 376) follows from Satz 4, whose proof (pp. 379-380) applies the
Hilfssatz (p. 378) to H = -Q_a^a Q_b^b U V for fixed natural numbers m and
j >= 1/2 + gamma.  For every such m the Hilfssatz already gives

    G_{m,j}(gamma) <= 2 log A + f'/f(1 + gamma) + g'/g(1 + gamma),

and Satz 4 is its limit m -> infinity.  This script evaluates the left side
G_{m,j} in closed form, checks that closed form against numerical
log-derivatives of Zimmert's definitions of Q_a, Q_b, U and V, evaluates the
series F_{a,b} of p. 376 independently, and prints:

* the Tabelle 3 (p. 377) rows with 3 <= n <= 20, with
  exp((G_{m,j} + n log pi) / 2) and the printed, rounded-down value;
* the gamma table (hundredths, integers excluded) that maximises G_{m,j} for
  each signature with 3 <= n <= 20;
* the bound ceil(10^15 exp(-(G_{m,j} + n log pi) / 2)) at |d| = 10^30 for
  the signatures pinned in test/t-generation-bound.cpp.

Its values are mpmath floating point, not interval enclosures; they are a
differential oracle for the Arb implementation in
src/factor_base/factor_base.cpp, which test/t-generation-bound.cpp pins.

Usage: python3 test/zimmert_twin_class_oracle.py [--self-check] [--table]
       [--gammas] [--bounds]   (no option: all of them)
"""

import argparse
import sys

from mpmath import (ceil, diff, exp, inf, log, loggamma, mp, mpf, nsum, pi,
                    psi)

mp.dps = 40

# Must match kZimmertTwinClassTerms in src/factor_base/factor_base.cpp.
M_TERMS = 128


def digamma(x):
    return psi(0, x)


def j_of(gamma):
    """Smallest natural number j with j >= 1/2 + gamma."""
    j = int(ceil(mpf(1) / 2 + gamma))
    return max(j, 1)


def g_parts(gamma, m, j=None):
    """The real-place, complex-place and U V parts of g_closed."""
    g = mpf(gamma)
    if j is None:
        j = j_of(g)
    d = 2 * (1 + 2 * g)
    half = mpf(1) / 2

    def delta(y):
        return digamma(y + half) - digamma(y)

    sa = mpf(0)
    sb = mpf(0)
    for l in range(1, m + j + 1):
        y = (2 * l - 1 + g) / d
        yp = (2 * l + g) / d
        if l <= m:
            sa += 4 / d * delta(y) - 1 / (2 * l - 2 - g) - 1 / (2 * l - 1 + g)
            sb += 4 / d * delta(yp) - 1 / (2 * l - 1 - g) - 1 / (2 * l + g)
        else:
            sa += 2 / d * delta(y) - 1 / (2 * l - 2 - g)
            sb += 2 / d * delta(yp) - 1 / (2 * l - 1 - g)
    u = 4 / d * (digamma(g / d) + digamma((1 + g) / d)
                 - digamma((1 + 3 * g) / d) - digamma((2 + 3 * g) / d))
    ca = -(digamma(half + g / 2) + digamma(-g / 2)) / 2
    cb = -(digamma(1 + g / 2) + digamma(half - g / 2)) / 2
    return sa + ca, sb + cb, u


def g_closed(a, b, gamma, m, j=None):
    """Closed form of the left side of the Hilfssatz for Zimmert's H_m."""
    pa, pb, u = g_parts(gamma, m, j)
    return a * pa + b * pb + u


def log_h(a, b, gamma, m, j, s):
    """log |H(s)| from Zimmert's definitions (p. 379), for real s."""
    g = mpf(gamma)
    d = 2 * (1 + 2 * g)

    def lg(x):
        return loggamma(x).real

    qa = mpf(0)
    qb = mpf(0)
    for l in range(1, m + j + 1):
        qa += lg((s + 2 * l - 2) / d) - lg((-s + 2 * l + 2 * g) / d)
        qb += lg((s + 2 * l - 1) / d) - lg((-s + 2 * l + 1 + 2 * g) / d)
    for l in range(1, m + 1):
        qa += lg((-s + 2 * l + 1 + 4 * g) / d) - lg((s + 2 * l - 1 + 2 * g) / d)
        qb += lg((-s + 2 * l + 2 + 4 * g) / d) - lg((s + 2 * l + 2 * g) / d)

    def log_u(t):
        return (lg((-t + 2 + 2 * g) / d) - lg(t / d) + lg((t + 2 * g) / d)
                - lg((-t + 2 + 4 * g) / d))

    return a * qa + b * qb + log_u(s) + log_u(1 - s)


def g_from_definitions(a, b, gamma, m, j=None):
    """-H'/H(1+gamma) + H'/H(-gamma) minus the Hilfssatz digamma terms."""
    g = mpf(gamma)
    if j is None:
        j = j_of(g)
    half = mpf(1) / 2
    right = diff(lambda s: log_h(a, b, g, m, j, s), 1 + g)
    left = diff(lambda s: log_h(a, b, g, m, j, s), -g)
    return (-right + left
            - a / 2 * (digamma(half + g / 2) + digamma(-g / 2))
            - b / 2 * (digamma(1 + g / 2) + digamma(half - g / 2)))


def f_series(a, b, gamma):
    """Zimmert's F_{a,b}(gamma), p. 376, summed by mpmath."""
    g = mpf(gamma)
    d = 2 * (1 + 2 * g)
    half = mpf(1) / 2
    sa = nsum(lambda l: 2 / (1 + 2 * g) * (digamma((2 * l + 3 * g) / d)
                                           - digamma((2 * l - 1 + g) / d))
              - 1 / (2 * l - 2 - g) - 1 / (2 * l - 1 + g), [1, inf])
    sb = nsum(lambda l: 2 / (1 + 2 * g) * (digamma((2 * l + 1 + 3 * g) / d)
                                           - digamma((2 * l + g) / d))
              - 1 / (2 * l - 1 - g) - 1 / (2 * l + g), [1, inf])
    ca = -(digamma(half + g / 2) + digamma(-g / 2)) / 2
    cb = -(digamma(1 + g / 2) + digamma(half - g / 2)) / 2
    c0 = (-4 / g + 2 / (1 + 2 * g)
          * (digamma((1 + g) / d) - digamma((1 + 3 * g) / d)
             + digamma((2 + 5 * g) / d) - digamma((2 + 3 * g) / d)))
    return a * (sa + ca) + b * (sb + cb) + c0


def constant(r1, r2, gamma, m=M_TERMS):
    """exp((G_{m,j} + n log pi) / 2), the constant c with N <= sqrt|d| / c."""
    n = r1 + 2 * r2
    return exp((g_closed(r1 + r2, r2, gamma, m) + n * log(pi)) / 2)


# Tabelle 3 (p. 377), rows with 3 <= n <= 20: r1, r2, printed gamma, printed
# value (rounded down) and one unit in its last printed digit.
TABLE3 = [
    (3, 0, "1.53", "5.56", "0.01"),
    (1, 1, "1.82", "3.95", "0.01"),
    (4, 0, "1.15", "18.5", "0.1"),
    (2, 1, "1.31", "12.1", "0.1"),
    (0, 2, "1.53", "8.32", "0.01"),
    (5, 0, "0.95", "68.5", "0.1"),
    (3, 1, "1.06", "42.7", "0.1"),
    (1, 2, "1.19", "27.4", "0.1"),
    (6, 0, "0.82", "273", "1"),
    (0, 3, "1.09", "62.9", "0.1"),
    (8, 0, "0.66", "5110", "10"),
    (0, 4, "0.87", "570", "1"),
    (10, 0, "0.56", "110000", "1000"),
    (0, 5, "0.74", "5770", "10"),
    (20, 0, "0.36", "1630000000000", "10000000000"),
    (0, 10, "0.47", "1480000000", "10000000"),
]

# Signatures whose bound at |d| = 10^30 test/t-generation-bound.cpp pins.
BOUND_SIGNATURES = [(3, 0), (1, 1), (4, 0), (2, 1), (0, 2), (6, 0), (0, 3),
                    (8, 0), (0, 4), (10, 0), (0, 10), (16, 0), (20, 0),
                    (1, 9)]


def self_check():
    ok = True
    for a, b, g, m in [(3, 0, "1.53", 3), (2, 1, "1.82", 4), (1, 1, "1.31", 2),
                       (5, 2, "0.95", 3), (4, 3, "1.37", 2), (6, 3, "0.47", 5),
                       (2, 2, "2.30", 2)]:
        c = g_closed(a, b, mpf(g), m)
        h = g_from_definitions(a, b, mpf(g), m)
        err = abs(c - h)
        print(f"self-check a={a} b={b} gamma={g} m={m}: closed={mp.nstr(c, 20)}"
              f" definitions={mp.nstr(h, 20)} |diff|={mp.nstr(err, 3)}")
        ok = ok and err < mpf(10) ** -15
    for a, b, g in [(3, 0, "1.53"), (1, 1, "1.82"), (20, 0, "0.36")]:
        gm = g_closed(a, b, mpf(g), M_TERMS)
        f = f_series(a, b, mpf(g))
        print(f"F - G_m a={a} b={b} gamma={g} m={M_TERMS}:"
              f" {mp.nstr(f - gm, 6)}")
        ok = ok and 0 < f - gm < mpf(10) ** -3
    return ok


def table():
    ok = True
    for r1, r2, g, printed, unit in TABLE3:
        v = constant(r1, r2, mpf(g))
        f = exp((f_series(r1 + r2, r2, mpf(g)) + (r1 + 2 * r2) * log(pi)) / 2)
        p = mpf(printed)
        good = p <= v < p + mpf(unit)
        ok = ok and good
        print(f"table3 r1={r1} r2={r2} gamma={g} printed={printed}"
              f" G_m={mp.nstr(v, 10)} F={mp.nstr(f, 10)}"
              f" margin={mp.nstr(v / p - 1, 3)} {'ok' if good else 'FAIL'}")
    return ok


_PARTS = {}


def best_gamma(r1, r2):
    """The gamma in hundredths, 0 < gamma < 3 and not an integer, that
    maximises G_{m,j} for the signature (first one on ties)."""
    if not _PARTS:
        for h in range(1, 300):
            if h % 100 != 0:
                _PARTS[h] = g_parts(mpf(h) / 100, M_TERMS)
    a = r1 + r2
    best = None
    for h, (pa, pb, u) in sorted(_PARTS.items()):
        v = a * pa + r2 * pb + u
        if best is None or v > best[1]:
            best = (h, v)
    return best[0]


def gammas():
    rows = {}
    for n in range(3, 21):
        row = [best_gamma(n - 2 * r2, r2) for r2 in range(n // 2 + 1)]
        rows[n] = row
        print(f"    {{{', '.join(str(x) for x in row)}}},  // n = {n}")
    return rows


def bounds():
    for r1, r2 in BOUND_SIGNATURES:
        g = mpf(best_gamma(r1, r2)) / 100
        c = constant(r1, r2, g)
        print(f"bound r1={r1} r2={r2} gamma={mp.nstr(g, 3)}"
              f" c={mp.nstr(c, 12)} ceil(1e15/c)={int(ceil(mpf(10)**15 / c))}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--table", action="store_true")
    parser.add_argument("--gammas", action="store_true")
    parser.add_argument("--bounds", action="store_true")
    args = parser.parse_args()
    run_all = not (args.self_check or args.table or args.gammas or args.bounds)
    ok = True
    if run_all or args.self_check:
        ok = self_check() and ok
    if run_all or args.table:
        ok = table() and ok
    if run_all or args.gammas:
        gammas()
    if run_all or args.bounds:
        bounds()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
