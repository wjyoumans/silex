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
