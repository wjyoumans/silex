Source Fidelity
===============

Mathematical behavior must be traceable to an approved requirement, an
established algorithm, or an existing tested Silex contract.  The public
:doc:`../algorithms_and_sources` page explains the policy and
:doc:`../legal_and_provenance` fixes the external source identities.

Authority order
---------------

When sources differ, apply this order:

#. the approved requirement for the current change;
#. cited PARI/GP or Hecke behavior for the algorithm being implemented;
#. current Silex behavior and regression tests; and
#. contracts of the FLINT primitives in use.

SageMath, OSCAR, and Nemo may provide additional comparison evidence, but they
do not silently replace a selected primary source.  Independent reasoning is
appropriate for glue code, tests, build code, diagnostics, documentation, and
clearly identified local mechanics; it must not invent a mathematical proof
rule or stopping condition.

Required change record
----------------------

Before changing mathematical behavior, record:

.. code-block:: text

   Source lineage:
   - Task and intended contract:
   - Existing Silex files and functions:
   - Upstream versions, files, and routines:
   - FLINT primitives involved:
   - Behavior and edge cases preserved:
   - Intentional deviations:
   - Tests added or updated:
   - Benchmarks added or updated:

Put durable conclusions in :doc:`../reference/algorithms_and_sources` or the
relevant conceptual API page.  Keep chronological investigation logs and
machine-local paths out of the public repository.

Guardrails
----------

Do not introduce a new algorithm, normalization, cache, backend default,
mandatory dependency, public API, or performance shortcut without an approved
need and evidence.  In particular, a failing fixture or benchmark is not by
itself permission to change relation selection, saturation, certification,
precision growth, or failure semantics.

If sources disagree, record the conflict, identify the applicable authority,
preserve the existing tested behavior until a choice is approved, and add a
test for the chosen contract.  If no source exists for nontrivial mathematics,
stop at research, scaffolding, or tests for known examples rather than
inventing production behavior.

Experiments must be isolated, labeled, and given an acceptance or stop
criterion.  Benchmark improvement alone is insufficient to replace
source-backed mathematics.

Recorded intentional deviations
-------------------------------

Unit reduction applies approximate coefficients.  Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Map.jl:reduce_mod_units`` (lines 638-658) computes
``round_approx`` when exact rounding raises ``InexactError``, then doubles the
precision and takes an unconditional ``continue``, so the approximate
coefficients are never applied.  Silex applies them and retries at the same
precision (``reduction_coordinates`` and ``reduce_relation_units_modulo`` in
``src/order_unit/compute.cpp``).  The function also reduces class-group
relation generators and unit candidates before the torsion test, so the
argument is not about a unit basis only: each output is its input times an
exact product of the current unit generators (exact ``FactoredElement``
arithmetic), so it generates the same principal ideal and, with the current
units, the same group.  Correctness does not depend on which integer
coefficients are applied; only which representatives are kept differs, so the
published generators and basis can differ from Hecke's.  The behavior is kept
so that published outputs stay stable; the user decided to keep it and
document it; outputs do not change.  An approximate rounding that gives all
zeros ends the loop at the current precision, where Hecke would double the
precision and might round exactly to a nonzero value; there the output can be
less reduced than Hecke's.  See
:doc:`../reference/algorithms_and_sources`.

The analytic imaginary-quadratic class number first evaluates
``L(1, chi_D)`` to an absolute tolerance sized for a ball of radius about
1/16 around ``h(D)``, instead of FLINT's mapping from a target precision to
a tolerance (FLINT 3.6.0 ``acb_dirichlet_l_fmpq_afe``).  FLINT's tail bound
and per-term working precisions are kept, so the enclosure stays rigorous,
and ``h(D)`` is still taken only when the ball holds one integer.  When that
step does not isolate ``h(D)``, the previous precision steps run unchanged,
so the result is the same and no input fails that succeeded before.  The
change was approved as a change of precision rule; see "Class groups and
order units" in :doc:`../reference/algorithms_and_sources`.
