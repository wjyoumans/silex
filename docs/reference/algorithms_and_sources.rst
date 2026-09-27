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

The independent test-only rational Gram--Schmidt oracle follows the
size-reduction and Lovasz inequalities in section 1.2.2 of Damien Stehle's
`Floating-Point LLL: Theoretical and Practical Aspects
<https://perso.ens-lyon.fr/damien.stehle/downloads/LLL25.pdf>`_
(author's chapter manuscript).  Parameters are converted from binary doubles
exactly, not replaced by decimal fractions.  Certifier benchmarks separately
measure ``fmpz_lll_is_reduced`` and ``fmpz_mat_is_reduced`` on prebuilt reduced
bases.  FLINT's ``src/fmpz_mat/is_reduced.c`` implements the latter using
rigorous ball arithmetic with rational fallback; this is distinct from the
``src/fmpz_lll/is_reduced_mpfr.c`` checker.

``Lat::lll_reduce`` first normalizes to full-row-rank HNF.  For 2--14 rows,
at most 14 columns, and coefficients of at most 250 bits in absolute value,
it tries ``fmpz_lll_d_with_removal_knapsack`` with no removal bound and accepts
only after ``fmpz_mat_is_reduced`` certifies the result with the unchanged
default parameters.  The reducer's status is not a reducedness certificate
or a rank.  A failed reduction or certificate restores the original HNF
before calling the original ``fmpz_lll`` path; ineligible inputs use that path
directly.  The existing zero-initialized transformation argument is retained.
This is Silex-specific routing, not FLINT's default checker order.  The source
anchors are FLINT ``src/fmpz_lll/wrapper_with_removal_knapsack.c`` (first phase
and mandatory certification) and ``lll_with_removal_ulll.c`` (the 250-bit
truncation boundary), inspected at revision
``b574c48f46ca505d6d48dd843c0fe95e1199ecdc``.  PARI 2.17.3
``src/basemath/lll.c:ZM_lll_norms`` likewise retains provable reduction after
its fast phases unless certification is explicitly disabled; Silex does not
disable it.  Tests force both rejection routes after scratch mutation and
require exact agreement with the original path.

Number fields, orders, and ideals
---------------------------------

``src/number_field`` and ``src/element`` retain established Silex behavior for
validated field construction, the polynomial generator, exact arithmetic,
trace, norm, conjugation, roots, and failure preservation.  PARI
``src/basemath/quad.c``, ``base1.c``, and ``base2.c`` are source anchors for
quadratic normalization, generic field data, and local data.  The natural
integral generator of a quadratic maximal order is kept distinct from the
field's polynomial generator.

Field construction requires a defining polynomial that is irreducible over
``Q``.  The check factors the primitive integral multiple of the polynomial
with FLINT ``fmpz_poly_factor`` (``src/fmpz_poly_factor``), which returns the
content separately and the remaining factors as primitive irreducible
polynomials with multiplicities; by Gauss's lemma the polynomial is
irreducible over ``Q`` exactly when that list is one factor of multiplicity
one.  ``Element::invert`` also rejects zero divisors exactly: through the
zero norm (FLINT ``nf_elem_norm``) in the quadratic representation, and in the
generic representation through the gcd of the same ``fmpq_poly_xgcd`` call
that FLINT's generic ``_nf_elem_inv`` branch (``src/nf_elem/inv.c``) makes.

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

Element powers and roots
------------------------

``Element::is_power`` and ``Element::is_square`` answer rational constants
exactly, route quadratic squares through the radicand, and otherwise use a
private port of the pure-power Hensel root finder for monic integral fields of
degree less than 10.  A definite ``true`` is published only after the lifted
candidate ``c`` satisfies ``c^n == a`` exactly; every unported or unverified
case returns failure (unsupported) and leaves the caller's root unchanged.

* Hensel lifting follows Hecke v0.38.6
  ``src/NumFieldOrd/NfOrd/Hensel.jl:_roots_hensel`` (lines 57-232, with
  ``ispure = true`` and ``is_normal = true``), its lifting loop and exact
  candidate check (lines 612-622), and the Friedrich-Fieker lifting exponent
  ``_lifting_expo`` (lines 674-756), with the equation order in place of
  ``any_order(K)``.
* Hecke skips a prime when the reduction of ``y^n - a`` is not squarefree
  (``Hensel.jl`` lines 153-155).  Since ``y^n - a`` is inseparable modulo
  ``p`` whenever ``p | n``, Silex skips every prime ``p | n`` for lifting
  directly.  The residue-field membership disproof below is valid at every
  good prime, including ``p | n``, so it is not skipped there.
* The residue test follows PARI 2.17.3
  ``src/basemath/FpX.c:Fq_ispower`` (lines 2509-2523): ``a`` is an ``n``-th
  power in ``F_q^*`` exactly when ``a^((q - 1)/gcd(q - 1, n)) = 1``.  A failed
  test at a good prime proves that ``a`` is not an ``n``-th power in ``K``.
* Residue roots follow PARI 2.17.3
  ``src/basemath/bb_group.c:gen_Shanks_sqrtn`` (lines 899-953): the Bezout
  exponent of ``n`` modulo ``q - 1`` gives the unique root when
  ``gcd(n, q - 1) = 1``, and otherwise ``y^n = a`` is reduced to
  ``y^g = a^u`` with ``g = gcd(n, q - 1)`` and ``n u + (q - 1) v = g``.  For
  an ``a`` that passed the residue test above, the two equations have the
  same roots in ``F_q``.  Without that test they can differ: for ``q = 11``,
  ``n = 6``, ``g = 2`` and ``u = 2``, every ``a`` solves ``y^2 = a^2``, but
  not every ``a`` is a sixth power.

Two Silex pre-filters run before the exact check ``c^n == a`` and only reject
candidates, so a rejection leaves the query unsupported and never changes a
definite answer.  Neither is in the upstream sources: Hecke's polynomial
degree ``n`` is a machine integer bounded by the polynomial it builds, while
Silex accepts ``n`` up to ``2^63 - 1`` and must not form ``c^n`` blindly.

* Norm filter: ``N(c)^n = N(a)`` by multiplicativity of the norm.  When
  ``N(c)`` is not ``0`` or ``+-1``, the comparison needs ``n`` below the bit
  size of ``N(a)``, so ``c^n`` is formed only for inputs at least that large.
* Height filter, for candidates of norm ``+-1``: if ``c^n = a`` then
  ``h(a) = n h(c)`` for the absolute logarithmic Weil height, so
  ``n log M(P_c) = log M(P_a)``, where ``P_x`` is the primitive integral
  characteristic polynomial of ``x`` over ``Q`` and ``M`` is its Mahler
  measure, ``log M(P_x) = log |lc(P_x)| + sum_i log max(1, |sigma_i(x)|)``
  over the ``[K : Q]`` complex embeddings (Bombieri and Gubler, *Heights in
  Diophantine Geometry*, sections 1.5 and 1.6).  The leading coefficient
  accounts for non-integral candidates.  A candidate is rejected only when
  Arb enclosures prove ``n log M(P_c) > log M(P_a) + 1``.  A passing
  candidate therefore has ``h(c^n) <= h(a) + (1 + r)/d``, where
  ``d = [K : Q]`` and ``r`` is the sum of the widths (twice the radii) of
  the two enclosures ``n log M(P_c)`` and ``log M(P_a)``.
  A root of unity (``M(P_c) = 1``) always passes, and binary powering keeps
  its coefficients bounded, so ``c^n`` stays cheap.  The working precision
  is 128 bits plus the largest coefficient bit size of ``c``, ``a`` and the
  defining polynomial.  If the enclosures are too wide to decide, the
  candidate goes on to the exact check.

Embedding root contexts
-----------------------

``src/embedding`` stores certified complex roots of the defining polynomial
for archimedean evaluation.  Hecke v0.38.6 ``src/Misc/acb_root_ctx.jl``
(``_roots!``) is the source: it refines by calling FLINT
``acb_poly_find_roots`` with the previous roots as initial approximations,
doubling the working precision, and accepts the result once all roots are
isolated, accurate enough, and pass ``acb_poly_validate_real_roots``.  The
first isolation, and the fallback, use FLINT
``arb_fmpz_poly_complex_roots``.  The signature comes from FLINT
``fmpz_poly_signature`` and is computed when the context is defined.

Silex deviates from the source in three intentional ways:

* No re-sort after refinement.  Hecke re-sorts with FLINT
  ``_acb_vec_sort_pretty``, whose ``acb_cmp_pretty`` comparator orders by
  ``|Im|`` only when the difference is resolved and otherwise by real part,
  so place indices can change with precision.  Silex matches each refined
  ball to the previous ball it overlaps and requires the match to be
  one-to-one.  Because both vectors isolate the roots with pairwise disjoint
  balls, a new ball that overlaps exactly one previous ball encloses the same
  root.  Place order is therefore fixed once roots are set, which cached
  logarithmic embeddings rely on.
* A precision cap on refinement from previous roots.  Hecke doubles up to
  ``2^22`` bits and then raises an error.  Silex stops that loop at eight times
  the requested precision and falls back to FLINT's certified isolation.  The
  factor is a tuning constant, not a correctness condition.
* A retrying fallback.  When a full isolation does not yet match the previous
  balls one-to-one, Silex isolates again at double the precision.  This
  terminates because each root is at positive distance from the closed
  previous balls that do not contain it.  Refinement fails only when the
  working precision would overflow, or for a non-squarefree polynomial.

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

Saturation-backed class-group proof follows Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Proof.jl:_class_group_proof``: verify generation
up to the Minkowski-type bound, then saturate at every prime dividing the
candidate class order.  Once the factor base generates the class group, the
full relation lattice contains the computed one with index ``h_cand / h``,
which divides ``h_cand``; the relations are complete exactly when they are
saturated at every prime ``p | h_cand``.  Silex therefore publishes a
saturation-backed ``proven`` only when all of the following hold: verified
generation up to the bound, with every prime ideal up to the bound in the
factor base; proven units and regulator; and a verified ``ell``-local proof
for every prime ``p`` dividing ``h_cand``.  The required primes are derived
from the published presentation; the internal index-bound gate proves the
union of the primes up to its bound and the prime divisors of ``h_cand``.
Gates that accept a supplied index bound or analytic class-regulator product
are not installed; public certification uses only values Silex computes or
verifies itself, and proven unit publication through the class/unit gates
re-derives the torsion subgroup of the order rather than trusting a stored
one.

An analytic ``hR`` proves a class/unit pair only when it is unconditional.
Silex has two unconditional routes: degree one, where the residue is exactly
one, and maximal orders of explicit quadratic-backend fields whose
discriminant ``|D|`` fits in a machine word, where the residue comes from
FLINT Dirichlet ``L(1, chi_D)``.  When that ``L(1, chi_D)`` evaluation fails
(for example the Dirichlet group cannot be initialized, or the ``L``-value
ball is not finite with a positive real part), the quadratic route falls back to Belabas--Friedman, and
the value is then GRH-conditional like any other; certification uses the
route that actually produced the value, not the field type.  Every other
``hR`` comes from the Belabas--Friedman evaluation (K. Belabas and E. Friedman,
"Computing the residue of the Dedekind zeta function", *Math. Comp.* 84
(2015), 357--369, Theorem 1), whose truncation-error bound assumes GRH.  Such
an ``hR`` is therefore at most ``grh`` evidence: the analytic check is
recorded with ``analytic_class_regulator_certification() == grh``, it never
marks units, the regulator, or relation saturation as proven, and it never
promotes a result to ``proven``.  For degree three and higher a ``proven``
request succeeds only through the saturation route above (generation, proven
units, and saturation at every ``p | h_cand``, following Hecke's
``_class_group_proof``) and otherwise fails closed.  Silex deviates from
Hecke's routine in two documented, conservative ways: it proves the unit
group first (as Hecke's ``_unit_group_proof`` does) and then runs the class
saturation against the proven units, instead of saturating relations and
unproven units jointly; and when an ``ell``-local test finds a missing root it
fails closed instead of enlarging the relations as Hecke's ``saturate!``
does.  Inside the paired transaction a Belabas--Friedman index-one check only
selects the candidate pair to prove; the unit group is then proven from the
regulator lower bound and unit saturation, and the class group by
saturation.  A failed saturation proof leaves the candidate units and class
group unchanged.

A Belabas--Friedman check recorded after an unconditional one (for example a
BF audit of a quadratic field already proven through ``L(1, chi_D)``) keeps
the unconditional label, and a Belabas--Friedman check is accepted only
against units that are already proven.  Because the Belabas--Friedman value
no longer contributes to a degree-three-or-higher ``proven`` result, such
results report ``analytic_class_regulator_status()`` and
``zeta_bf_proof_status()`` as ``not_checked`` unless an audit is requested
separately; these statuses are audit records, not proof components.

The ``ell``-local test takes the Smith-invariant ``beta`` rows (compact power
witnesses of invariant generators whose invariant ``ell`` divides), the free
units, and the torsion generator when ``ell`` divides the torsion order, and
requires their ``ell``-th-power characters at degree-one primes ``q`` with
``ell | q - 1`` to reach full rank modulo ``ell``.  This row choice follows
PARI 2.17.3 ``src/basemath/buch3.c:check_prime`` and ``primecertify`` (the
``bnfcertify`` local check); the torsion row matches Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Saturate.jl:compute_candidates_for_saturate``,
which adds the torsion generator when ``ell`` divides the torsion order.
Where Hecke would enlarge the relations with ``saturate!`` and repeat, Silex
fails closed: an ``ell`` whose local test does not reach full rank within the
auxiliary-prime bound is recorded as ``unavailable`` and does not promote.
Each certification entry point is transactional; a failed call restores all
class-group certification metadata, including the unit and regulator proof
statuses.  A successful ``try_prove_relation_saturation_with_units`` call
means only that its ``ell`` verified; the relation-saturation status becomes
``verified`` only once every ``p | h_cand`` is covered.

The order-unit saturation proofs build degree-one primes ``(q, theta - r)``
from roots ``r`` of the defining polynomial ``f`` modulo ``q`` with an
internal constructor.  It accepts a root only when ``f(r) = 0 mod q`` and
``f`` is squarefree modulo ``q``, the Dedekind-Kummer hypothesis under which
``(q, theta - r)`` is a prime of residue degree and ramification index one
(Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138,
Theorem 4.8.13); this is the same condition as the degree-one fast path of
``decompose_prime``.

The unconditional unit proof ``OrderUnitGroup::prove_index_bound`` bounds the
index of the supplied full-rank subgroup from a regulator lower bound, then
proves ``ell``-saturation for every prime ``ell`` up to that index bound.  An
``ell`` is verified only when the ``ell``-th-power characters at degree-one
primes ``q = 1 mod ell``, not dividing the order discriminant, have an empty
common kernel.  The auxiliary-prime search follows Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Saturate.jl``: ``compute_candidates_for_saturate``
iterates the primes ``q = 1 mod ell`` (``PrimesSet`` with no upper end) and
stops when the kernel is empty or its dimension has stayed the same for more
than ``stable`` times the number of input elements; ``saturate!`` starts at
``stable = 3.5``, adjoins candidates that are ``ell``-th powers and repeats,
and doubles ``stable`` after a candidate that is not.  Silex multiplies
``stable`` by the current kernel dimension instead of the number of input
elements; the stopping rule decides only when candidates are tested, never
whether ``ell`` is verified.  Silex first scans the primes
``q <= aux_bound`` (the caller's bound; 1000 in class/unit validation) and
continues with the Hecke scan when that pre-scan selects no usable prime,
when it leaves a nonempty kernel without adjoining a root, or when
``ell >= aux_bound``.  The auxiliary bound is therefore not a proof
parameter; a unit proof record stores it as the pre-scan bound of the call,
not as the largest prime used.  Hecke starts its scan at
``next_prime(2^60)``; Silex starts at ``q = ell + 1``, which changes only
which valid characters are used.  Roots are adjoined only after the exact
order-unit check described below.

The index bound comes from ``unit_lower_regulator_bound``.  For a subgroup
``G`` of full rank in ``O^x``, with ``O`` an order of the field ``K``,
``[O^x : G] = R_G / Reg(O) <= R_G / R_K``.  So any proven lower bound
``R_lower <= R_K`` gives ``[O^x : G] <= floor(R_G / R_lower)``.
``OrderUnitGroup::regulator_index_bound`` takes the floor of the rigorous
upper endpoint of the Arb quotient, following Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Proof.jl:_unit_group_proof``
(``floor(tent_reg / low_reg)``).  If that upper endpoint is below one, the
inputs contradict ``R_K >= R_lower``, and the bound fails closed instead of
reporting index one.  ``R_lower`` is the largest lower endpoint of the
following proven bounds.  ``r1`` and ``r2`` are the numbers of real and
complex places, ``n = r1 + 2 r2``, and ``w`` is the number of roots of unity
of ``K``:

* E. Friedman, "Analytic formulas for the regulator of a number field",
  *Invent. Math.* 98 (1989), 599--622, Theorem B (p. 599): ``R >= 0.2052``
  for every number field.  The minimum is the regulator 0.20521... of the
  sextic field of discriminant ``-10051``.
* H. Zimmert, "Ideale kleiner Norm in Idealklassen und eine
  Regulatorabschätzung", *Invent. Math.* 62 (1981), 367--380, Korollar (i)
  (p. 375): ``2R/w >= 0.04 exp(0.46 r1 + 0.1 r2)``, used as
  ``R >= 0.02 w exp(0.46 r1 + 0.1 r2)``.
* Friedman 1989, Corollary on p. 620:
  ``R/w > 0.0031 exp(0.241 n + 0.497 r1)``.
* Zimmert 1981, Satz 3 (p. 374), which holds for every ``gamma > 0``:
  ``R/w >= (1+gamma)(1+2 gamma)/2 * Gamma(1+gamma)^(r1+r2) *
  Gamma(3/2+gamma)^r2 * 2^(-r1-r2) * pi^(-r2/2) *
  exp{(-1-gamma)[(r1+r2) psi((1+gamma)/2) + r2 psi(1+gamma/2) + 2/gamma +
  1/(1+gamma)]}``.  Silex evaluates it with Arb at 15 fixed rational values
  of ``gamma`` between 1/10 and 3.  For every signature of degree at most 20,
  the best of these is within 6% of the optimum over ``gamma``.
* Friedman 1989, Table 6 (p. 621): per-signature lower bounds for the
  signatures of positive unit rank that it lists, transcribed exactly.  The
  ``(r1, r2) = (0, 3)`` entry 0.27 excludes three fields (note b:
  ``D_K = -10051, -10571, -12167``).  Theorem B gives their regulators
  0.2052, 0.2132 and 0.2372, so Silex uses 0.2052 for ``(0, 3)``.

All of these bounds increase with ``w``.  For an order ``O``, ``O^x`` has
finite index in ``O_K^x``, so ``Reg(O) >= R_K``, and the field's ``w`` is the
right value to use.  If ``w`` cannot be computed, Silex uses ``w = 2``.  The
terms are evaluated at a working precision of at most 64 bits; the result is
an exact lower endpoint, so this only makes the bound slightly weaker.

Hecke's ``lower_regulator_bound``
(``src/NumFieldOrd/NfOrd/Unit/Regulator.jl``) uses
``max(0.054, 0.04 w exp(0.46 r1 + 0.01 r2))``.  That does not match
Zimmert's Korollar: it keeps Zimmert's 0.04 for ``2R/w`` and also multiplies
by ``w``, and it uses ``0.01 r2`` where Zimmert has ``0.1 r2``.  No source
proves the Hecke value, so Silex does not use it.

Termination needs only that the supplied subgroup ``U`` has full rank, so
that ``[O^x : U]`` is finite; it does not depend on the index bound being
correct.  For ``ell`` not dividing the torsion order ``w``, a kernel element
that is not an ``ell``-th power in ``K`` is cut by a positive density of the
primes ``q`` by the Chebotarev density theorem, and each adjoined root
multiplies the index of the subgroup by ``ell``, which can happen only
finitely often.  For a non-maximal order the scan may instead stop at a root
outside the order, as described below.  When ``ell`` divides ``w`` (always
for ``ell = 2``, since ``w`` is even), the characters also cover the torsion
generator, so the kernel has one column more than the free rank.  The root
step accepts only kernel rows over the free generators, so it cannot adjoin a
root of a torsion-twisted row: when such a scan stops at its stability
threshold with a nonempty kernel, the ``ell``-local proof fails closed as
``unavailable`` rather than continuing.

Silex adds one resource guard: a single scan examines at most ``2^16``
rational primes ``q``, and reaching it records the ``ell``-local proof as
``unavailable``.  The guard now also bounds the scans for
``ell >= aux_bound``, which previously ran until the machine-word range of
``q`` was exhausted.  A scan counts at most ``d`` degree-one primes per
``q`` in a field of degree ``d``, so once the doubled stability threshold
exceeds ``d * 2^16`` a scan can end only with an empty kernel or at the
guard.  One ``ell`` therefore runs at most ``17 + log2 d`` scans that end
without adjoining a root, plus one scan per adjoined root, and adjoined roots
are bounded by the ``ell``-part of the finite index ``[O^x : U]`` (and by the
caller's restart limit).  In total one ``ell`` examines at most about
``(17 + log2 d + log_ell [O^x : U]) * 2^16`` primes ``q``.

Residue-character ``ell``-saturation of unit subgroups of a non-maximal order
``O`` is a Silex extension; neither primary source applies this proof to
non-maximal orders.  Its soundness argument: for a maximal ideal ``P`` of
``O`` with ``O/P = F_q`` and ``q = 1 mod ell``, reduction ``O -> O/P`` is a
ring homomorphism, so an ``ell``-th power in ``O^x`` maps to an ``ell``-th
power in ``F_q^x``.  An empty character kernel for a subgroup ``U`` whose
torsion lies in ``U`` therefore gives ``U ∩ (O^x)^ell = U^ell``.  An
``ell``-th root found for a kernel row is computed in the field and is
adjoined only after an exact ``O``-unit check.  When the root lies in
``O_K \ O``, the step fails closed with status ``unavailable``: the root still
reduces into ``O/P`` at the primes the proof uses, so such a row never leaves
the kernel, and the characters cannot distinguish it from a row whose root
lies in ``O``.  The source-backed alternative is the exact sequence
``1 -> O^x -> O_K^x -> (O_K/f)^x / (O/f)^x`` for the conductor ``f``, so
that ``O^x`` is the kernel of the map from ``O_K^x``, used by Hecke v0.38.6
``_unit_group_non_maximal`` (``PicardGroup.jl`` under
``src/NumFieldOrd/NfOrd``); it is a planned follow-up.

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
Belabas--Friedman-style bounds used by the class/unit proof consumers.  The
Belabas--Friedman error bound assumes GRH; see "Class groups and order units"
for how certification consumes it.  Arb
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
