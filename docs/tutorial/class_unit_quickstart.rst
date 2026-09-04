Class and Unit Computation
==========================

The class/unit pipeline computes a finite class-group presentation and an
order-unit subgroup over the same order.  The native C++ API makes
certification status explicit: usable standalone class-group candidate data
and certified paired data are different outcomes.  The paired entry point is
one native transaction and accepts only ``grh`` or ``proven`` certification.
It computes class and unit candidates privately and publishes both only after
the requested proof contract succeeds.  ``unknown``, ``heuristic``, and
invalid paired requests fail before computation.  The main compiled example is
:download:`class_unit_computation.cpp <../../examples/class_unit_computation.cpp>`.
The focused fragments below emphasize API contracts; the complete compiled
example is included at the end.

Compute With Explicit Options
-----------------------------

Construct the field and order first, then pass a bound and options object to
the paired class/unit compute path:

.. code-block:: cpp

   namespace sflint = silex::flint;

   sflint::Fmpz d;
   sflint::fmpz_set_si(sflint::FmpzRef(d), 2);

   silex::NumberField K =
           silex::NumberField::quadratic(sflint::FmpzConstRef(d));
   silex::Order O = silex::Order::equation_order(K);
   if (!K.is_defined() || !O.is_defined()) {
       return false;
   }

   sflint::Fmpz factor_base_bound;
   sflint::fmpz_set_si(sflint::FmpzRef(factor_base_bound), 2);

   silex::ClassGroupComputeOptions options;
   options.max_candidates = 256;
   options.max_relations = 48;
   options.requested_certification = silex::CertificationMode::proven;

   silex::ClassGroupContext class_group;
   silex::OrderUnitGroup units;
   if (!units.compute_with_class_group(
           class_group, O, sflint::FmpzConstRef(factor_base_bound),
           options, 160)) {
       return false;
   }

The paired options expose resource ceilings and proof resources.  Relation
search radii, unit-witness targets, and post-finite schedules are private,
deterministic policy; callers cannot select alternate implementations through
this options object.

For an explicitly unproven standalone class-relation candidate, use the
separate resource-only options type:

.. code-block:: cpp

   silex::ClassGroupCandidateOptions candidate_options;
   candidate_options.max_candidates = 256;
   candidate_options.max_relations = 48;

   silex::ClassGroupContext candidate;
   if (!candidate.compute_candidate(
           O, sflint::FmpzConstRef(factor_base_bound),
           candidate_options)) {
       return false;
   }
   assert(candidate.certification_status() ==
          silex::CertificationMode::unknown);

Read Published Data And Proof Metadata
--------------------------------------

Standalone class-group candidate outputs are exact mathematical data for the
computed presentation, but they are not automatically proof of the full class
group.  Paired outputs are published only with their requested certification.
Query the coarse certification mode and component proof records separately:

.. code-block:: cpp

   auto class_order = class_group.order();
   auto torsion_order = units.torsion_order();
   auto regulator = units.regulator();

   silex::CertificationMode class_cert =
           class_group.certification_status();
   silex::CertificationMode unit_cert =
           units.certification_status();
   silex::ProofState saturation =
           class_group.relation_saturation_status();
   silex::ProofState zeta =
           class_group.zeta_bf_proof_status();

Use ``CertificationMode::proven`` only when a documented proof gate has
verified the published state.  ``CertificationMode::unknown`` remains useful
for standalone ``ClassGroupContext::compute_candidate`` calls; it is not a
paired transaction mode.  ``CertificationMode::grh`` means the paired
transaction completed and its result is conditional on GRH; it does not imply
unconditional component proof records.  GRH certification semantics and
coverage remain work in progress, and successful results retain the requested
``grh`` label rather than being silently promoted or downgraded.

When requested-proven output is required, set the requested certification
field in ``ClassGroupComputeOptions`` and branch on the returned status and
proof metadata.  The native transaction publishes both class and unit objects
atomically only after its proof contract passes.  Increasing proof
budgets is appropriate only when metadata shows that a budgeted proof
component was the limiting factor.

Next Steps
----------

- Read :doc:`../reference/class_units_compact` for class-group,
  finite-presentation, order-unit, certification, and proof-record APIs.
- Read :doc:`../support_matrix` for the qualified platform and exact paired
  class/unit presentations.
- Read :doc:`../reference/factored_zeta` for compact/factored elements and
  Belabas-Friedman zeta audit records used by proof drivers.
- Run ``example-class-unit-computation`` for the paired class/unit workflow.
- Run ``example-relation-kernel-order-units`` and
  ``example-supplied-order-units`` for lower-level unit-publication paths.

Complete Compiled Example
-------------------------

.. literalinclude:: ../../examples/class_unit_computation.cpp
   :language: cpp
   :caption: examples/class_unit_computation.cpp
