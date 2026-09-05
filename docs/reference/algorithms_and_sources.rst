Algorithms and Sources
======================

Silex is a native C++ library built on FLINT.  Its public API, storage model,
diagnostics, and failure contracts are Silex-owned.  Several mathematical
algorithms retain documented lineage from established systems; that lineage
is an engineering and attribution obligation, not a runtime dependency.

Authority and provenance
------------------------

Production mathematical behavior follows the authority order in
:doc:`../algorithms_and_sources`: an approved Silex requirement, cited
upstream behavior, current Silex behavior and tests, then the contracts of the
FLINT primitives in use.  A source translation does not automatically define
the public C++ API.

The 0.1.1 foundation was imported from a historical Silex native-library tree.
SHA ``7fdce40f30abbb08e024253347efba0d4e5fcc5a`` is a private archival anchor,
not an object in the public repository and not a public reproducibility claim.
Retained translations and adaptations require ``GPL-3.0-or-later``
distribution.  See :doc:`../legal_and_provenance` for public source anchors
and legal attribution.

External source anchors are fixed for reproducibility:

* PARI/GP 2.17.3 is the primary source for the class-group, unit, ideal,
  relation-search, enumeration, zeta, and S-unit behavior identified below.
* Hecke.jl v0.38.6, commit
  ``74215ba3d34f296e6f709e415e8007d225524287`` and root tree
  ``5758221d5c6c176b4781dbafe267bb056099d56b``, supplies selected order,
  ideal, LLL relation-search, class/unit, compact-element, and saturation
  behavior.
* Hecke.jl v0.39.19, commit
  ``122658620f5ac3c8260785c06d9ce7062f037498`` and root tree
  ``b538e80d8c6d6eb783b99b0a217c864458cc7b99``, separately supplies the
  S-class/S-unit source and comparison baseline.
* FLINT contracts govern exact arithmetic, storage, HNF/SNF, polynomial,
  number-field, and modular-linear-algebra primitives used directly by Silex.

Sparse matrices and lattices
----------------------------

``src/fmpz_smat`` implements canonical sparse integer rows and matrices,
dense/modular conversion, transpose, multiplication, one-shot and persistent
modular rank, and incremental exact HNF/index contexts.  These operations use
FLINT ``fmpz``, ``fmpz_mat``, and ``nmod_mat`` contracts; sparse rows are
strictly column-sorted and never store zero values.

``src/lat`` implements exact row lattices, HNF and transformation matrices,
containment, sum, intersection, index, saturation, LLL reduction, and bounded
short-vector enumeration.  FLINT HNF/SNF and exact solve contracts are the
baseline for canonical operations; optional fplll and flatter integrations
remain backend comparisons rather than mathematical authorities.

The LLL reference checks compare FLINT's optional transformation argument in
three forms: a zero matrix, an identity-initialized matrix, and a null pointer.
They check the unimodular row identity, canonical HNF, independent reducedness,
and identical reduced bases within the same build.  This is a regression gate,
not a promise of a canonical reduced basis across backend versions.  PARI
2.17.3 ``src/basemath/lll.c:ZM_lll_norms`` omits transformation tracking for
``LLL_INPLACE``, as used by ``src/basemath/buch2.c:SPLIT``; this supplies the
basis-only source precedent.  FLINT's ``src/fmpz_lll/lll.c``,
``lll_with_removal_ulll.c``, and ``test/t-lll.c`` establish the optional
transformation boundary, including retention of internal truncation
transformations.  Native tests cover both sides of the 250-bit truncation
threshold and the honesty-search benchmarks check every searched prime-ideal
lattice in their fixed quartic and quintic fixtures.

Number fields, orders, and ideals
---------------------------------

``src/number_field` and ``src/element`` retain established Silex behavior for
validated field construction, the polynomial generator, exact arithmetic,
trace, norm, conjugation, roots, and failure preservation.  PARI
``src/basemath/quad.c``, ``base1.c``, and ``base2.c`` are source anchors for
quadratic normalization, generic field data, and local data.  The natural
integral generator of a quadratic maximal order is kept distinct from the
field's polynomial generator.

``src/order`` implements validated order bases, trace forms, discriminants,
multiplication tables, indices, p-maximal overorders, and global maximal-order
construction.  PARI ``src/basemath/base2.c`` (including the ``nfmaxord``
family) and Hecke
``src/NumFieldOrd/NfOrd/MaxOrd/MaxOrd.jl`` and
``MaxOrd/DedekindCriterion.jl`` are the main external anchors.  The standard
quadratic identity ``disc(O_f) = disc(O_K) f^2`` governs direct conductor
orders.

``src/ideal`` implements integral and fractional canonical row-HNF ideal
lattices, principal construction, containment, norm, sum, intersection,
product, inverse, colon, and normalized numerator/denominator storage.  PARI
``src/basemath/base4.c`` and Hecke
``src/NumFieldOrd/NfOrd/Ideal/Ideal.jl`` are the principal ideal-arithmetic
anchors.  Output is published only after validation, so documented failures
preserve prior caller state.

General ideal multiplication uses PARI 2.17.3
``src/basemath/base4.c:idealmul_aux``, ``idealmulelt``,
``get_random_a``, ``ok_elt``, ``idealHNF_mul`` and ``idealHNF_mul_two``
for scalar/principal handling, content removal, verified two-generator
extraction and multiplication.  Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Ideal/Ideal.jl:assure_has_minimum`` supplies the
coordinate-solve treatment of order bases whose first element is not one.
Silex uses row HNF throughout, so upstream column indices are not copied.

The extraction checks equality of the complete ideal lattice, not merely
norm equality.  Reduced multiplication matrices are used only during search;
the exact multiplication matrix is reconstructed before forming the product.
Search tries basis generators and at most 32 reproducibly seeded random
combinations.  This finite budget is a Silex resource-policy deviation from
the upstream unbounded search, not a completeness claim.  Exhaustion falls
back to the original exhaustive basis-product algorithm.  General extraction
requires a known maximal order; scalar and known-principal arithmetic retain
their existing order scope.  Witnesses are operation-local, with no new
public cache, normalization, certification or failure-publication contract.
The fixed general-product lattices in ``test/data/ideal_product_lattices.json``
were checked in a common polynomial power basis against PARI/GP 2.17.4;
representative quadratic, cubic and quartic products were also checked with
Hecke v0.39.19.  These executable comparison versions are distinct from the
algorithm-source versions above.

Local algebra and relations
---------------------------

``src/prime_ideal``, ``src/residue_ring``, ``src/residue_field``,
``src/factor_base``, ``src/ideal_factorization``, and ``src/relation`` cover
prime decomposition, reduction, valuation, residue arithmetic, factor bases,
factorization over a base, and relation rows.  Hecke
``src/NumFieldOrd/NfOrd/Ideal/Prime.jl`` and ``ResidueField.jl`` and PARI's
``idealprimedec`` behavior are the external decomposition and residue-field
anchors.  Dedekind--Kummer factor data is interpreted in the correct integral
generator, not silently relabeled as polynomial-generator data.

Class groups and order units
----------------------------

``src/class_group``, ``src/order_unit``, and ``src/zeta`` implement
factor-base relation collection, finite abelian presentations, relation and
unit saturation, analytic class-regulator validation, proof metadata, and
transactional paired publication.  PARI 2.17.3
``src/basemath/buch2.c`` supplies primary relation-completion, analytic
finishing, HNF/regulator reconstruction, saturation, and unit-log validation
logic.  ``base1.c``, ``base2.c``, ``base3.c``, and ``bibli1.c`` supply
supporting ideal, relation-search, enumeration, and zeta behavior.

Hecke ``src/NumFieldOrd/NfOrd/Clgp.jl`` and its class-group implementation
supply selected LLL relation-search, relation-unit reduction, validation, and
saturation control flow.  Silex keeps candidate, ``grh``, and proven outcomes
distinct and never upgrades certification without the corresponding proof
record.  ``test/t-class-group.cpp`` and ``test/t-order-unit.cpp`` are the
focused regression surfaces.  Cross-engine campaign orchestration is outside
this repository's scope.

The legacy factor-base honesty search tests principal witnesses directly in
order coordinates.  PARI 2.17.3
``src/basemath/buch2.c:divide_p_elt`` and ``can_factor`` supply the exact
norm-accounting criterion: the required prime has valuation one, and
residue-degree-weighted valuations over that prime and the factor base
account for the entire absolute norm.  Missing primes, including other
primes above the same rational prime, cannot be discarded.  A completed
negative classification is final; unsupported or failed direct evaluation
falls back to the existing full ideal factorization predicate.

This routing change preserves Silex's scalar/lattice candidate search,
random draws, bounds, proof-target selection, and certification labels.
It does not select the separate T2 witness search.  The noninstalled
reference selector and scan audit support differential tests without a
public option or persistent cache.  Frozen quartic and quintic benchmark
batches contain 672 and 3937 candidates respectively; each classification
was checked against PARI/GP 2.17.4 full ideal factorization using exact
prime-ideal lattices in a common polynomial power basis.  The executable
oracle version is distinct from the algorithm-source version above.

The noninstalled integral-scalar predicate additionally follows
``src/basemath/base4.c:Q_nffactor`` and ``prV_e_muls`` from PARI 2.17.3:
``v_P(m) = e_P v_p(m)`` for a nonzero integer ``m``.  Factoring ``abs(m)``
suffices; at each rational prime dividing it, the sum of ``e_P f_P`` over
the retained ideals and the required ideal must equal the field degree.
The required ideal must have valuation exactly one.  This is exact support
accounting, not a norm-only comparison of arbitrary ideals.  The predicate
requires a known maximal order and reports unsupported evaluations separately
from negative classifications.  The legacy search uses it only for candidates
already known to be integral scalars.  Failed evaluations fall back through
the general order-element predicate and its full-factorization fallback;
completed negative classifications do not fall back.  The scalar helper's
arbitrary-precision input does not enlarge the search: the initial rational
prime retains its machine-integer eligibility check, and the scalar loop,
candidate order, bounds, random draws, and non-scalar path are unchanged.

Two exact edge routes have narrower routine-level anchors.  PARI 2.17.3
``src/basemath/buch2.c:Buchall_deg1`` and the degree-at-most-one branch in
``Buchall_param`` publish the trivial class group, regulator one, torsion
order two, and no free units for ``QQ``.  Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp.jl:_validate_class_unit_group`` likewise
short-circuits degree-one validation only when the tentative class number and
regulator are both one.  Silex uses its exact degree-one route for either
paired request, while retaining the caller's requested public ``grh`` label
when that was the request.

For negative fundamental discriminants of absolute value at most 12, PARI
2.17.3 ``src/basemath/quad.c:classno`` and ``classno2`` return exact class
number one before their general algorithms.  Silex's discriminant ``-3``
route independently enumerates reduced forms through FLINT
``qfb_reduced_forms``, requires exact class order one, reconstructs and
verifies a principal generator for every factor-base ideal, and only then
publishes the identity relation lattice.  The exceptional route does not
change the generic restart policy, certification labels, or failure-atomic
publication contract.

Factored elements and analytic products
---------------------------------------

``src/factored_element`` implements formal products and compact elements,
structural roots, explicit evaluation, logarithmic embeddings, and compact
base representations.  PARI ``famat`` behavior and Hecke
``src/Misc/FactoredElem.jl`` and ``CompactRepresentation.jl`` supply the
formal-product and compact-power baselines.  Exact expansion is explicit and
is not a hidden side effect of structural operations.

``src/zeta`` follows the PARI analytic class-regulator and
Belabas--Friedman-style bounds used by the class/unit proof consumers.  Arb
balls preserve explicit precision and error information; an inconclusive
interval does not become a proof by heuristic rounding.

S-class and S-unit publication
------------------------------

``src/sunit`` publishes proven S-class and S-unit data from proven ordinary
class/unit inputs.  PARI 2.17.3 ``src/basemath/bnfunits.c:119-235`` supplies
selected-prime class rows, HNF/kernel construction, valuation bases, nonunit
generators, and the nonempty-S regulator formula.  Hecke 0.39.19
``src/NumFieldOrd/NfOrd/Clgp/Sunits.jl:1-265`` supplies the corresponding
class-group quotient, valuation-HNF, compact-generator, and combined-group
maps.  Silex deliberately orders coordinates as torsion, ordinary free, then
nonunit and exposes owned result objects rather than a reference-system
container layout.

Maintenance requirements
------------------------

A mathematically significant change must identify the upstream version and
file or routine when one exists, state the behavior being preserved or
changed, add focused invariant or differential tests, and update this page
when lineage changes.  Temporary range-level investigation belongs with task
evidence, not in the public repository.  Do not introduce a proof rule,
normalization, stopping condition, certification label, or performance
heuristic merely to make a fixture pass.

The tracked-tree source-name gate allows source-project names in documentation
and legal attribution.  It rejects retired project identities and upstream
project-name tokens in active source paths and contents, public headers,
identifiers, tests, fixtures, examples, native benchmarks, tools, diagnostics,
protocol fields, and build labels.  Cross-engine corpora and campaign
orchestrators are not part of the Silex source tree.  See
:doc:`../development/source_fidelity` for the contributor-facing procedure.
