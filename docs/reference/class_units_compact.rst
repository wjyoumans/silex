Class Groups and Unit Groups
============================

The higher-level number-field stack exposes native C++ APIs for finite
presentations, candidate class groups, order-unit subgroups, field-level unit
helpers, proof metadata, and saturation/proof-driver workflows.  Factored
elements, compact presentations, and zeta/class-regulator product audit APIs
are documented separately in :doc:`factored_zeta`.

Class Groups
------------

``silex::FiniteAbelianGroup`` is the finite row-relation quotient used by
relation matrices and class-group publication.  It exposes Smith-normal-form
invariants, class order, invariant-generator matrices, invariant-generator
relation combinations, row reduction, invariant coordinates, and Smith
left-kernel relation-dependency rows.  It is documented here as part of the
class/unit presentation layer rather than as a separate algebra hierarchy.

``silex::ClassGroupContext`` owns the native class-group candidate state and
stores its ``Order`` parent handle by value.  It manages copied factor bases,
relation rows, retained
relation generators, row-module HNF state, finite quotient presentations,
class-generator representatives, ideal-class coordinates, compact invariant
power witnesses, relation-kernel unit witnesses, and proof/certification
metadata.
Simple presentation exports such as relation rows, finite presentation,
invariants, class order, and invariant-generator matrices have owned-return
convenience helpers in addition to the explicit output-buffer APIs.

``ClassGroupCandidateOptions`` is the resource-only options object for
``compute_candidate``.  It exposes hard candidate/relation ceilings and an
optional diagnostics context.  Search radii, relation-kernel targets, and
post-finite schedules are deterministic library policy rather than public
algorithm selectors.  A successful call publishes an exact presentation of
the relation lattice that was found, always with
``CertificationMode::unknown``; it does not claim that this lattice is the
full relation lattice.  Invalid input or resource exhaustion leaves the
existing context unchanged.

``ClassGroupComputeOptions`` belongs to the paired class/unit transaction.  It
exposes hard resource ceilings, bounded relation-saturation resources, the BF
zeta cutoff, requested certification, and diagnostics.  Its proof request is
not accepted by ``compute_candidate``.  Proven status is published only by
the paired transaction or by a separate explicit source-backed proof gate.

Paired Class/Unit Transaction
-----------------------------

``OrderUnitGroup::compute_with_class_group`` is the single paired entry point.
It accepts requested ``grh`` or ``proven`` certification, fixes its private
source-backed execution policy before relation production, and computes into
temporary class/unit objects.  It publishes both objects atomically only after
the requested class, unit, and regulator contract succeeds.  Unsupported
``unknown`` or ``heuristic`` paired requests fail without modifying either
caller output.  Standalone class-group candidate computation remains available
through ``ClassGroupContext::compute_candidate``.  Completion is qualified
only for the presentations in :doc:`../support_matrix`; other quadratic and
higher-degree routes are fail-closed best effort.

Certification Metadata
----------------------

Class-group and order-unit objects expose certification as data, not as
implicit success.  A successful compute call means the object was constructed
for the supported path; it does not by itself mean the class group or full
unit group is proven.  Callers should read the coarse mode and then inspect
the relevant proof records when they need to explain or audit the result.

``CertificationMode`` records the strength of the published result:

- ``unknown``: a candidate was constructed, but no accepted proof gate has
  certified it;
- ``heuristic``: a future heuristic path may publish this mode, but current
  proof-oriented callers should not treat it as proven;
- ``grh``: the paired class/unit transaction completed under the current GRH
  route and retained the requested public label; its semantics and coverage
  are work in progress, and it is conditional and provisional rather than an
  unconditional component proof;
- ``proven``: a source-backed proof gate has certified the published state.

Detailed proof fields use ``ProofState``:

- ``not_checked``: the proof path was not attempted for the current object;
- ``unavailable``: the proof path was attempted or requested but its
  prerequisites were not available;
- ``verified``: the corresponding proof record verified the relevant
  condition.

These enums are intentionally small and serializable.  They are suitable for
Python/Sage bindings and for proof-driver logs because they do not encode
hidden global state or depend on raw parent pointer identity.

For ``ClassGroupContext``, use ``certification_status()`` for the coarse
published mode, and inspect the specific proof metadata when deciding why a
result is or is not proven.  The current API exposes checked factor-base
generation status, relation-saturation status and records, analytic
class-regulator proof status, Belabas-Friedman zeta proof records, and the
unit/regulator proof statuses used by paired class/unit certification.
Record-shaped metadata can be read through named result structs such as
``ClassGroupFactorBaseGenerationRecord``,
``ClassGroupRelationSaturationRecord``, and
``ClassGroupZetaBfProofRecord``.  The older explicit output-parameter methods
remain available for failure-preserving callers.

The class-group proof metadata is deliberately componentized:

- factor-base generation records audit rational primes against the stored
  build/generation bounds;
- relation-saturation records audit the ``ell``-local saturation checks;
- unit and regulator proof statuses record the paired class/unit inputs used
  by the certification gate;
- analytic class-regulator status records whether the analytic index-bound
  gate verified the finite presentation;
- BF zeta proof records keep the cutoff, maximum cutoff, requested precision,
  working precision, and error bound used by the Belabas-Friedman proof path.

The coarse class-group mode is promoted to ``proven`` only when the relevant
source-backed gate has verified the components it needs.  For example, a
bounded relation search may produce a usable finite presentation while still
reporting ``CertificationMode::unknown`` if the proof records are
``not_checked`` or ``unavailable``.

For ``OrderUnitGroup``, use ``certification_status()`` for the coarse mode and
``unit_proof_record_count()``, ``unit_proof_record()``,
``unit_proof_verified()``, ``regulator_index_bound()``,
``class_regulator_product()``, and ``class_regulator_index_bound()`` for the
metadata produced by local saturation and analytic index-bound proof paths.
Single-value queries such as torsion order and regulator have owned-return
convenience helpers.  Multi-output proof records remain explicit output
parameters and also have the named ``OrderUnitProofRecord`` convenience form.

Order-unit proof records are per-prime records.  Each
``OrderUnitProofRecord`` stores the proof prime ``ell``, the proof state, the
auxiliary-prime bound used for local evidence, the number of local primes
selected, and whether the proof/saturation pass changed the subgroup.  A
positive regulator or a full-rank supplied subgroup is not the same thing as a
proven full unit group; the certification mode and proof records are the
public source of truth.

The examples ``class_unit_computation.cpp``,
``relation_kernel_order_units.cpp``, and ``supplied_order_units.cpp`` show the
expected public pattern: compute with explicit options, query group data
through domain objects, and then read certification/proof metadata explicitly.

Order Units and Saturation
--------------------------

``silex::OrderUnitGroup`` stores its ``Order`` parent handle by value together
with torsion data, compact free generators, subgroup regulator metadata,
certification status, saturation status, and unit-proof records.  It covers
proven rank-zero torsion groups, caller-supplied compact full-rank subgroups,
relation-kernel unit subgroups, bounded dependent-relation refinement, local
saturation, index-bounded saturation, proof-prime selection, local proof
publication, restartable index-bound proof, and paired class/unit transactions.

Public ``parent()`` accessors on these objects are borrowed compatibility
views of the stored handles.  They are suitable for local access, but pointer
addresses should not be used as mathematical identity tokens.

The generic positive-rank paths remain bounded and source-backed.  They do
not claim full unit-group certification unless a documented proof slice
records the required proof status.

The standalone constructors have distinct fixed contracts:

- ``compute`` publishes a proven full group only when its exact standalone
  route succeeds; the presentations qualified for 0.1.1 are listed in
  :doc:`../support_matrix`;
- ``set_units`` publishes the exact full-rank subgroup generated by the
  supplied units with unknown certification;
- ``set_relation_kernel_units`` delegates rank zero to the proven ``compute``
  path; in positive rank it publishes the exact full-rank subgroup found in
  the supplied relation kernel with unknown certification;
- ``set_relation_kernel_units_bounded`` applies its fixed bounded dependent-
  relation refinement and retains the last valid subgroup if no later
  refinement verifies;
- ``set_relation_kernel_units_index_bounded`` derives one regulator index
  bound, applies the bounded refinement when available, and otherwise retains
  the initial relation-kernel subgroup;
- ``set_relation_kernel_units_index_bounded_saturated`` applies the fixed
  index-bounded construction followed by the caller-bounded saturation pass;
  if that optional pass cannot complete, it publishes the index-bounded
  subgroup with ``changed=false`` and ``stable=false``.

All input/construction failures preserve the caller's existing group.  A
successful subgroup constructor is not proof of the full unit group unless
its returned certification status explicitly says so.

Benchmark Coverage
------------------

``bench/nf_clgp/b-class_group_pipeline.cpp`` contains the class/unit benchmark
surface.  The target includes class-group candidate rows, the paired
class/unit matrix, and public random-sweep/random-boundary rows.

Private diagnostic rows and counters are not part of the native public API.
Do not add private benchmark rows unless they can be expressed through the
native diagnostics/profiling surface without expanding the public API.

Source Lineage
--------------

PARI 2.17.3 ``src/basemath/buch2.c`` is the primary source for relation
completion, analytic finishing, saturation, HNF/regulator reconstruction, and
unit-log validation.  Hecke's class-group implementation supplies selected
LLL relation-search, unit-extraction, validation, and saturation control flow.
Silex owns the native presentation, options, proof records, request validation,
and transactional publication.  See :doc:`algorithms_and_sources` for the
file-level public map; update that map whenever mathematical lineage changes.

Examples
--------

The following examples exercise this layer:

- ``examples/class_unit_computation.cpp``
- ``examples/rank_zero_order_units.cpp``
- ``examples/real_quadratic_order_units.cpp``
- ``examples/relation_kernel_order_units.cpp``
- ``examples/supplied_order_units.cpp``

See :doc:`factored_zeta` for compact/factored examples and zeta/BF audit
examples used by the proof-driver layer.
