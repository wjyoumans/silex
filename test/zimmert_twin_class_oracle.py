#!/usr/bin/env python3
"""Independent mpmath oracle for Zimmert's twin-class generation bound.

H. Zimmert, "Ideale kleiner Norm in Idealklassen und eine
Regulatorabschaetzung", Invent. Math. 62 (1981), 367-380.  Korollar 1
(p. 376) bounds N(a) N(b) <= d exp(-(F_{a,b}(gamma) + n log pi)) for
a = r1 + r2, b = r2 and Zimmert's series F_{a,b} of p. 376.  This script
sums that series directly with mpmath's nsum (Richardson/Shanks, checked
against Levin's transformation by --self-check), with no tail estimate, and
prints:

* the Tabelle 3 (p. 377) rows with 3 <= n <= 20, with
  exp((F + n log pi) / 2) and the printed, rounded-down value;
* the gamma table (hundredths, integers excluded, where the series has
  poles) that maximises F for each signature with 3 <= n <= 20;
* the bound ceil(10^15 exp(-(F + n log pi) / 2)) at |d| = 10^30 for the
  signatures pinned in test/t-generation-bound.cpp.

Its values are mpmath floating point, not interval enclosures; they are a
differential oracle for the Arb enclosure of F in
src/factor_base/factor_base.cpp (GB-zimmert-F-tail-bound), which
test/t-generation-bound.cpp pins.

Usage: python3 test/zimmert_twin_class_oracle.py [--self-check] [--table]
       [--gammas] [--bounds]   (no option: all of them)
"""

import argparse
import sys

from mpmath import ceil, exp, inf, log, mp, mpf, nsum, pi, psi

mp.dps = 40


def digamma(x):
    return psi(0, x)


def f_parts(gamma, method="r+s"):
    """The a-part, b-part and constant of Zimmert's F_{a,b}(gamma), p. 376:
    F_{a,b} = a * pa + b * pb + c0.  gamma must not be an integer."""
    g = mpf(gamma)
    d = 2 * (1 + 2 * g)
    half = mpf(1) / 2
    sa = nsum(lambda l: 2 / (1 + 2 * g) * (digamma((2 * l + 3 * g) / d)
                                           - digamma((2 * l - 1 + g) / d))
              - 1 / (2 * l - 2 - g) - 1 / (2 * l - 1 + g), [1, inf],
              method=method)
    sb = nsum(lambda l: 2 / (1 + 2 * g) * (digamma((2 * l + 1 + 3 * g) / d)
                                           - digamma((2 * l + g) / d))
              - 1 / (2 * l - 1 - g) - 1 / (2 * l + g), [1, inf],
              method=method)
    ca = -(digamma(half + g / 2) + digamma(-g / 2)) / 2
    cb = -(digamma(1 + g / 2) + digamma(half - g / 2)) / 2
    c0 = (-4 / g + 2 / (1 + 2 * g)
          * (digamma((1 + g) / d) - digamma((1 + 3 * g) / d)
             + digamma((2 + 5 * g) / d) - digamma((2 + 3 * g) / d)))
    return sa + ca, sb + cb, c0


def f_series(a, b, gamma, method="r+s"):
    """Zimmert's F_{a,b}(gamma), p. 376, summed by mpmath."""
    pa, pb, c0 = f_parts(gamma, method)
    return a * pa + b * pb + c0


def constant(r1, r2, gamma):
    """exp((F + n log pi) / 2), the constant c with N <= sqrt|d| / c."""
    n = r1 + 2 * r2
    return exp((f_series(r1 + r2, r2, gamma) + n * log(pi)) / 2)


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
    """Checks the nsum values against Levin's transformation."""
    ok = True
    for a, b, g in [(3, 0, "1.53"), (1, 1, "1.82"), (20, 0, "0.36"),
                    (2, 2, "2.30")]:
        f = f_series(a, b, mpf(g))
        h = f_series(a, b, mpf(g), method="levin")
        err = abs(f - h)
        print(f"self-check a={a} b={b} gamma={g}: F={mp.nstr(f, 30)}"
              f" |nsum - levin|={mp.nstr(err, 3)}")
        ok = ok and err < mpf(10) ** -30
    return ok


def table():
    ok = True
    for r1, r2, g, printed, unit in TABLE3:
        v = constant(r1, r2, mpf(g))
        p = mpf(printed)
        good = p <= v < p + mpf(unit)
        ok = ok and good
        print(f"table3 r1={r1} r2={r2} gamma={g} printed={printed}"
              f" value={mp.nstr(v, 10)}"
              f" margin={mp.nstr(v / p - 1, 3)} {'ok' if good else 'FAIL'}")
    return ok


_PARTS = {}


def best_gamma(r1, r2):
    """The gamma in hundredths, 0 < gamma < 3 and not an integer, that
    maximises F_{r1+r2,r2} for the signature (first one on ties)."""
    if not _PARTS:
        for h in range(1, 300):
            if h % 100 != 0:
                _PARTS[h] = f_parts(mpf(h) / 100)
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
