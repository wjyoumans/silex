#!/usr/bin/env python3
"""Independent mpmath oracle for Zimmert's Satz 2 generation bound.

H. Zimmert, "Ideale kleiner Norm in Idealklassen und eine
Regulatorabschaetzung", Invent. Math. 62 (1981), 367-380.  Satz 2 (p. 372):
every ideal class of a number field of degree n = r1 + 2 r2 and
discriminant d contains an integral ideal a with log(sqrt|d| / N(a)) >= R,

  R = r1 (-psi((1+g)/2) - log Gamma(1/2+g) + log Gamma(1+g) + (1/2) log pi)
    + r2 (-2 psi(1+g) + 2 log 2 + log(1/2+g) + log pi)
    - 2/(g-a) - log[(1+1/a) (1+1/g)^(-2) (1+1/(2g-a))^(-1)],

for every g > a > 0, with a = g - g(g+1)/sqrt(1+3g+3g^2) from Bemerkung 1
(p. 373).  This script evaluates R with mpmath and prints:

* --table1: the Tabelle 1 (p. 368) rows with n = 20 and n = 100, with
  exp(R) at the printed gamma and the printed, rounded-down value;
* --gammas N: for every signature with 3 <= n <= N, the gamma = k/100,
  1 <= k <= 300, maximising R; with --check FILE it compares the
  21 <= n <= 64 rows with src/factor_base/zimmert_gamma_table.inc;
* --bounds: gamma and ceil(10^(e/2) exp(-R)) at |d| = 10^e for the
  signatures pinned in test/t-generation-bound.cpp (kZimmertHighRows);
* --fields: the combined generation bound (smallest of Satz 2, Minkowski
  with (4/pi)^r2 and the integer Minkowski form with 2^r2) for the fields
  of degree above 20 pinned in kFieldRows.

Its values are mpmath floating point at 80 digits, not interval enclosures;
they are a differential oracle for the Arb evaluation in
src/factor_base/factor_base.cpp, which test/t-generation-bound.cpp pins.

Usage: python3 test/zimmert_satz2_oracle.py [--table1] [--gammas N
       [--check FILE]] [--bounds] [--fields]
"""

import argparse
import re
import sys

from mpmath import (ceil, digamma, exp, factorial, log, loggamma, mp, mpf,
                    pi, sqrt)

mp.dps = 80


def satz2(r1, r2, gamma):
    g = mpf(gamma)
    a = g - g * (g + 1) / sqrt(1 + 3 * g + 3 * g * g)
    assert 0 < a < g
    real = (-digamma((1 + g) / 2) - loggamma(mpf(1) / 2 + g)
            + loggamma(1 + g) + log(pi) / 2)
    cplx = -2 * digamma(1 + g) + 2 * log(2) + log(mpf(1) / 2 + g) + log(pi)
    return (r1 * real + r2 * cplx - 2 / (g - a)
            - log((1 + 1 / a) * (1 + 1 / g) ** -2 * (1 + 1 / (2 * g - a)) ** -1))


def best_gamma(r1, r2):
    best = None
    for k in range(1, 301):
        value = satz2(r1, r2, mpf(k) / 100)
        if best is None or value > best[1]:
            best = (k, value)
    return best[0]


def table1():
    rows = [(20, 0, 36, "4.332e11"), (0, 10, 46, "5.736e8"),
            (100, 0, 14, "8.448e72"), (0, 50, 18, "2.417e55")]
    for r1, r2, g, printed in rows:
        value = exp(satz2(r1, r2, mpf(g) / 100))
        print(f"r1={r1} r2={r2} gamma={g} exp(R)={mp.nstr(value, 10)} "
              f"printed={printed}")


def gammas(nmax, check):
    table = None
    if check:
        text = open(check).read()
        blocks = re.findall(r"// n = (\d+)\n([\d,\s]+)\}", text)
        table = {int(n): [int(x) for x in re.findall(r"\d+", b)]
                 for n, b in blocks}
    mismatches = 0
    for n in range(3, nmax + 1):
        row = [best_gamma(n - 2 * r2, r2) for r2 in range(n // 2 + 1)]
        print(f"n={n}", " ".join(map(str, row)))
        if table is not None and 21 <= n <= 64:
            if table.get(n) != row:
                mismatches += 1
                print(f"MISMATCH n={n} table={table.get(n)}")
    if table is not None:
        checked = sum(1 for n in table if 21 <= n <= min(nmax, 64))
        print(f"checked rows {checked} mismatches {mismatches}")


# Signatures pinned in kZimmertHighRows.  |d| = 10^e with e = 2 ceil(R /
# log 10 + 15), so that each bound is about 10^15.
BOUND_ROWS = [(21, 0), (1, 10), (0, 11), (0, 12), (32, 0), (0, 16), (0, 20),
              (64, 0), (0, 32), (65, 0), (1, 32), (81, 0), (100, 0), (0, 50)]


def bounds():
    for r1, r2 in BOUND_ROWS:
        g = best_gamma(r1, r2)
        r = satz2(r1, r2, mpf(g) / 100)
        e = 2 * int(ceil(r / log(10) + 15))
        value = mpf(10) ** (e // 2) * exp(-r)
        print(f"{{{r1}, {r2}, {e}, {g}, \"{int(ceil(value))}\"}},")


def fields(path):
    # Signatures and discriminants of the kFieldRows fields of degree above
    # 20, read from the pinned rows (computed with an independent computer
    # algebra system).
    text = open(path).read()
    pattern = r'\{"([^"]+)", (\d+), (\d+),\s*((?:"\d+"\s*)+),'
    for m in re.finditer(pattern, text):
        name, r1, r2 = m.group(1), int(m.group(2)), int(m.group(3))
        disc = int("".join(re.findall(r'"(\d+)"', m.group(4))))
        n = r1 + 2 * r2
        if n <= 20:
            continue
        root = sqrt(mpf(disc))
        g = best_gamma(r1, r2)
        zimmert = int(ceil(root * exp(-satz2(r1, r2, mpf(g) / 100))))
        mink = int(ceil(factorial(n) / mpf(n) ** n * (4 / pi) ** r2 * root))
        isqrt = int(ceil(root))
        while isqrt * isqrt < disc:
            isqrt += 1
        while (isqrt - 1) ** 2 >= disc:
            isqrt -= 1
        num = int(factorial(n)) * isqrt * 2 ** r2
        generic = max(1, -(-num // n ** n))
        print(f"{name} gamma={g} zimmert={zimmert} minkowski={mink} "
              f"generic={generic} bound={min(zimmert, mink, generic)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table1", action="store_true")
    parser.add_argument("--gammas", type=int)
    parser.add_argument("--check")
    parser.add_argument("--bounds", action="store_true")
    parser.add_argument("--fields", metavar="T_GENERATION_BOUND_CPP")
    args = parser.parse_args()
    if args.table1:
        table1()
    if args.gammas:
        gammas(args.gammas, args.check)
    if args.bounds:
        bounds()
    if args.fields:
        fields(args.fields)
    return 0


if __name__ == "__main__":
    sys.exit(main())
