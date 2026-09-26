# Tests

Use CTest. Add focused tests for public behavior, source-backed algorithms,
failure preservation, parent lifetimes, and benchmark-facing workflows.

Every mathematical behavior change needs parity tests and documented source
lineage.

The fixed class/unit adapter inventory in `data/class_unit_fields.json` has
seventeen proven replays and five provisional GRH replays. The proven rows use
twenty-second adapter budgets. The five unchanged GRH fixtures retain their
ten-second manifest defaults and the protocol test's existing twenty-second
override. A timeout is a failure; these budgets must not be increased after a
failed observation. The manifest test checks the exact population and the five
added rows' coefficients, invariants, ranks, signatures and maximal-order
discriminants.

The additions replay existing native transactions: degree one, discriminants
-3 (both presentations) and -4, and canonical real quadratic 210, exercising
the exact degree-one/imaginary-quadratic and canonical Dirichlet routes. These
are coverage additions, not new support promises. Existing higher-degree proven
labels remain regression expectations; they do not establish general
unconditional applicability of BF-based completion. Five GRH fixtures retain
paired `grh` labels and regulator proof `not_checked`. Exact proven routes may
also retain analytic/zeta `not_checked`, while their factor-base, saturation,
unit and regulator components must remain `verified`.
The canonical real quadratic 210 route instead verifies the analytic
class-regulator component and leaves relation saturation and BF `not_checked`;
its source S-unit receipts preserve that distinction. These exact component
expectations follow `ClassGroupContext::try_certify_analytic_class_regulator_`
and `try_promote_proven_certification_` rather than a blanket interpretation
of the coarse `proven` label.

`t-class-unit-matrix.cpp` verifies invariant ideal powers against their compact
principal witnesses, ideal-class coordinates, free-unit principal ideals,
mathematical parents, torsion orders and finite positive regulators (exactly
one in rank zero). It preserves the incomplete-candidate, cross-parent,
zero-resource and paired failure snapshots.

`test_sunit_instance.py` executes all six rows of `data/sunit_fields.json` with
twenty-second observation budgets. It checks manifest-order `(p, beta)` prime
descriptors, exact invariants/ranks/valuation indices, component proof states,
and the adapter's membership results. `t-sunit.cpp` additionally checks:

- quotient identities `A_i^d_i * product(P_j^e_ij) = (alpha_i)`;
- exact compact/expanded valuations and nonunit principal-ideal identities,
  which exclude support at every prime outside S;
- compact and expanded coordinate round trips, signed exponents and torsion
  normalization, including empty-S projections of the exact edge fields;
- the empty-S regulator `R_K` and nonempty-S convention
  `R_K * h_S * product(log N(P))`, where `h_S` is the S-class number, not the
  valuation-lattice index;
- verified membership, mathematical nonmembership (including the index-three
  lattice obstruction), unknown results for invalid input, and unchanged
  coordinate or paired output objects on rejection/failure.

The expanded-element preimage overload rejects zero during input validation
because `FactoredElement::set_element` excludes zero. The tests preserve its
`unknown` failure and unchanged coordinates, separately from mathematical
nonmembership for valid nonzero inputs.

Native prime selection matches the exact two-generator ideals, independently
of decomposition order. Existing split-prime reordering and parent-lifetime
regressions remain. The source contracts are documented in
`docs/reference/class_units_compact.rst`, `docs/reference/sunit_groups.rst` and
`docs/reference/algorithms_and_sources.rst`: PARI 2.17.3 `buch2.c`/`quad.c`,
Hecke v0.38.6 class/unit validation, and the S-unit baselines PARI
`bnfunits.c:119-235` and Hecke v0.39.19 `Sunits.jl:1-265`. No new upstream code
is translated here. Oracle agreement has not been established for the full
inventory.
