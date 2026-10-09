# Offline verification tool; not part of the build or CTest.
#
# What it checks: every entry E = s/10^k of
# src/factor_base/zimmert_twin_class_table.inc (117 signatures, 3 <= n <= 20)
# is at least c_up(r1, r2, gamma), the certified upper bound on
# k(r1, r2; gamma) = exp(-(F + n log pi)/2) computed by silex-research
# (src/silex_research/zimmert.py, Arb at 128 bits), using exact rational
# comparison.  E >= c_up implies ceil(E sqrt|d|) >= B(gamma), which is what
# lets Silex's twin-class bound replace B(gamma) (see
# docs/reference/algorithms_and_sources.rst).
#
# Requires: Sage (python module sage.all) and a silex-research checkout at
# commit 8926120.
#
# Run (from the repository root; the sage launcher runs .py files):
#   sage tools/check_zimmert_twin_class_table.py \
#       <silex-research>/src src/factor_base/zimmert_twin_class_table.inc
#
# Expected: one "ok" line per signature and the final line
#   entries 117 failures 0 min E/c_up 1.0 at r1,r2,gamma ('5', '5', '50')
# (minimum ratio 1.0 means equality to the printed precision at gamma = 1/2).
import re, sys
sys.path.insert(0, sys.argv[1])
from fractions import Fraction
from sage.all import QQ
from silex_research.zimmert import c_up
text = open(sys.argv[2]).read()
rows = re.findall(r'\{(\d+), "(\d+)", (\d+)\},\s*// r1 = (\d+), r2 = (\d+)', text)
bad = 0
worst = None
for g, s, k, r1, r2 in rows:
    E = QQ(int(s)) / QQ(10) ** int(k)
    c = QQ(c_up(int(r1), int(r2), QQ(int(g)) / 100, 128))
    ratio = E / c
    ok = E >= c
    if not ok:
        bad += 1
    if worst is None or ratio < worst[0]:
        worst = (ratio, r1, r2, g)
    print(r1, r2, g, "ok" if ok else "FAIL", float(ratio))
print("entries", len(rows), "failures", bad, "min E/c_up", float(worst[0]), "at r1,r2,gamma", worst[1:])
sys.exit(1 if bad else 0)
