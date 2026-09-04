S-Class Groups and S-Unit Groups
================================

The public S-class/S-unit layer is declared in ``silex/sunit.hpp``.  It
publishes owned, move-only ``SClassGroup`` and ``SUnitGroup`` results from an
already proven ordinary ``ClassGroupContext`` and ``OrderUnitGroup``.  The
entry point is ``compute_sunit_groups``.

Construction Contract
---------------------

Construction requires:

- a maximal order;
- a proven class-group context with a finite presentation;
- a proven ordinary unit group for the same order; and
- a list of distinct, defined prime ideals of that order.

Both output objects are replaced atomically only after the S-class
presentation, compact witnesses, nonunit generators, valuation basis, and
S-regulator have been verified and copied into owned storage.  On failure the
previous outputs are preserved.  ``SUnitComputeResult`` records the failed
stage and, when applicable, the selected-prime index.

This layer does not select or expose a class/unit computation strategy.  It
consumes the already certified ordinary groups and has no alternate
class/unit dispatch path.

S-Class Publication
-------------------

``SClassGroup`` owns the quotient of the ordinary ideal class group by the
classes of the selected primes.  It exposes:

- finite-group invariants and order;
- invariant ideal representatives;
- compact principal witnesses for invariant powers; and
- the corresponding selected-prime exponent rows.

For invariant representative ``A_i``, invariant ``d_i``, witness ``alpha_i``,
and selected exponent row ``e_i``, the published data has been checked
exactly as

.. math::

   A_i^{d_i} \prod_j P_j^{e_{ij}} = (\alpha_i).

The object also records its proof state and source class-group certification.

S-Unit Publication
------------------

``SUnitGroup`` owns the ordinary unit group together with one compact
nonunit generator for each selected prime and the square valuation matrix of
those generators at the selected primes.  Its stable coordinate order is:

1. torsion exponent;
2. ordinary free-unit exponents; and
3. nonunit exponents.

``SUnitCoordinates`` stores those three components explicitly.  This order
matches the domain-oriented public contract; callers do not need to depend on
external-engine container layouts.  Within the nonunit component, selected
primes retain the exact order supplied to ``compute_sunit_groups``; the same
order indexes selected-prime accessors, valuation columns, and nonunit
coordinates.

``image`` accepts exact coordinates and publishes either a compact
``FactoredElement`` or an expanded ``Element``.  ``preimage`` first computes
the selected valuations, solves the nonunit valuation lattice exactly,
removes that component, and invokes the ordinary-unit coordinate kernel on
the residual.  It then verifies the complete image/preimage round trip as an
exact field-element equality.

Membership has three outcomes:

- ``verified``: exact coordinates were published;
- ``not_sunit``: a selected valuation was outside the nonunit lattice or the
  residual was not an ordinary unit; and
- ``unknown``: input, arithmetic, precision, or exact verification failed.

Coordinates are replaced only for ``verified`` membership.  In particular,
rejection of an element with support outside the selected primes preserves the
caller's previous coordinate object.

S-Regulator
-----------

For nonempty ``S``, the regulator follows PARI ``bnfsunit``:

.. math::

   R_S = R_K h_S \prod_{P \in S} \log N(P),

where ``R_K`` is the proven ordinary regulator and ``h_S`` is the published
S-class number.  Arb carries the result and its requested working precision.

PARI has a distinct empty-``S`` branch that returns the ordinary regulator
directly.  Silex preserves that behavior, so an empty selected-prime list
reproduces the ordinary unit group and ordinary regulator.

Proof Metadata
--------------

The public objects separate source and derived proof facts:

- source class and ordinary-unit certification;
- source relation-saturation, unit, and regulator proof states;
- S-class, S-unit, and S-regulator proof states;

Successful publication requires proven source groups and exact private
S-class/mod-units proof records.

Source Lineage
--------------

PARI 2.17.3 ``src/basemath/bnfunits.c:119-235`` and Hecke
``src/NumFieldOrd/NfOrd/Clgp/Sunits.jl:1-265`` provide the selected-prime
quotient, valuation-HNF, compact-generator, coordinate, and regulator
baselines.  Silex implements its owned result and publication contracts in
``src/sunit``.  The focused regression target is ``t-silex-sunit``; external
comparison commands and schema policy are documented in
:doc:`../benchmarking`.  See :doc:`algorithms_and_sources` for the complete
public source map.
