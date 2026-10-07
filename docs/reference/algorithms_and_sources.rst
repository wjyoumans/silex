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
baseline for canonical operations.  The coordinate bound that routes
enumeration is the Fincke--Pohst bound: U. Fincke and M. Pohst, "Improved
methods for calculating vectors of short length in a lattice, including a
complexity analysis", *Math. Comp.* 44 (1985), 463--471, and H. Cohen, *A
Course in Computational Algebraic Number Theory*, GTM 138, section 2.7.3 (the
citations given at ``enum_double_coordinates_bounded`` in ``src/lat/lat.cpp``).
Silex uses the bound only to choose the route for uncapped calls
(``max_coord < 0``); the enumeration itself is a witness search as described
in :doc:`lat`.  Optional fplll and flatter integrations remain backend
comparisons rather than mathematical authorities.

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
with FLINT ``fmpz_poly_factor`` (``src/fmpz_poly_factor``; contract from FLINT
3.6 ``src/fmpz_poly_factor/factor.c``), which returns the
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

Roots of unity
--------------

``roots_of_unity``, ``root_of_unity_order`` and ``root_of_unity_generator``
in ``src/unit/roots_of_unity.cpp`` compute the number ``w`` of roots of unity
of a field and a generator.  Fields with a real place, quadratic fields and
``x^4 + 1`` have direct answers.  Every other field follows the polynomial
branch of PARI/GP 2.17.3 ``src/basemath/nffactor.c:nfrootsof1``:

1. ``guess_roots``: for primes ``p >= 3`` not dividing ``disc(T)`` of the
   monic integral defining polynomial ``T`` of ``theta' = L theta``, ``w``
   divides ``p^g - 1`` with ``g`` the gcd of the residue degrees above ``p``.
   ``T`` and the positive integer ``L`` come from PARI's
   ``src/basemath/base1.c:ZX_primitive_to_monic`` (lines 390-446), which
   ``nfmaxord`` applies to non-monic polynomials: for each ``p^e`` exactly
   dividing the leading coefficient of the primitive integral numerator
   ``a(x)`` of the defining polynomial, ``L`` takes the factor ``p^k`` for
   the least ``k`` that makes ``p^(k n - e) a(x / p^k)`` integral.
   ``L = 1`` and ``theta' = theta`` when ``a`` is monic.
   The gcd of these values is a proven multiple of ``w``.  As in PARI, the
   loop stops when the gcd fits in a word and has not changed for more than
   ``n + 20`` consecutive good primes.  As in Hecke v0.38.6
   ``src/NumFieldOrd/NfOrd/TorsionUnits.jl:_torsion_group_order_divisor``
   (lines 274-345), it also stops at once when the gcd is 2, and it resets
   the stability count while ``phi(gcd)`` does not divide ``n = [K : Q]``:
   ``Q(zeta_w)`` is a subfield of ``K``, so ``phi(w) | n``, and a gcd with
   ``phi(gcd)`` not dividing ``n`` is a strict multiple of ``w``, so the
   reset only delays the stop, so the final gcd still divides the one
   the reference ``guess_roots`` rule would return.

   By the Chebotarev density theorem the gcd over all good primes equals
   ``w``: for each prime ``l`` dividing ``w``, let
   ``l^a = l^(v_l(w) + 1)``, so ``zeta_(l^a)`` is not in ``K``, and let
   ``M`` be the Galois closure of ``K(zeta_(l^a))`` over ``Q``.  The
   subgroup of ``Gal(M/Q)`` fixing ``K`` does not fix ``zeta_(l^a)``, so it
   contains some ``sigma`` that moves it.  The good primes whose Frobenius
   is conjugate to ``sigma`` have positive density; each has a
   residue-degree-one prime in ``K``, so ``gcd f = 1``, and satisfies
   ``p != 1 (mod l^a)``, so ``l^a`` does not divide ``p^(gcd f) - 1 = p - 1``.
   Only finitely many ``l`` divide the first term, so the gcd over all good
   primes is ``w``.  So a bound larger than ``w`` comes only from stopping
   early, never from a wrong gcd value.
2. The degree and ramification conditions of ``nfrootsof1`` lower each prime
   power of the bound: ``Q(zeta_(p^k))`` has degree ``(p - 1) p^(k-1)`` and a
   known ``p``-adic discriminant valuation.  Like PARI's polynomial branch,
   Silex uses ``disc(T)``, which gives a weaker but still valid reduction.
3. When ``phi(bound) = n``, ``ZXirred_is_cyclo_translate`` tests whether
   ``T(y) = Phi_N(+-y + c)``, reading ``c`` from the trace coefficient.  Silex
   takes the two candidates ``+-(theta' + c)`` and certifies them by exact order
   instead of PARI's Graeffe comparison.
4. Otherwise Silex looks for a primitive ``p^e``-th root of unity for each
   prime power ``p^e > 2`` of the bound and multiplies them.  PARI uses
   ``nfisincl(polcyclo(p^e), T)``; Silex follows Hecke v0.38.6
   ``TorsionUnits.jl:_torsion_units_gen`` (lines 414-456).  The search runs
   in a monic integral model: the field itself when its defining polynomial
   is monic and integral, and otherwise ``Q(theta')`` defined by ``T``.  A
   root ``g(theta')`` found there is the element ``g(L theta)`` of the
   field.

   a. Silex first tries its exact square and power roots, which need degree
      less than 10: ``zeta_4 = sqrt(-1)``, then square roots up to
      ``zeta_(2^e)``; ``zeta_3 = (-1 + sqrt(-3))/2``, then cube roots up to
      ``zeta_(3^e)``.  Any ``p``-th root of a primitive ``p^j``-th root of
      unity is a primitive ``p^(j+1)``-th root.  If ``zeta_(p^(j+1))`` lies
      in ``K``, the ``p``-th power map from the cyclic group
      ``mu_(p^(j+1))`` onto ``mu_(p^j)`` is surjective, so every ``p^j``-th
      root of unity is a ``p``-th power in ``K``.  Together these make the
      choice of root at each step irrelevant.
   b. When that fails or is unsupported (``p >= 5``, or degree 10 or more),
      Silex finds a root of ``Phi_(p^e)`` in ``K`` as Hecke's
      ``_roots_hensel(Phi_(p^e), max_roots = 1, is_normal = true,
      root_bound = ones)`` does (``src/NumFieldOrd/NfOrd/Hensel.jl`` lines
      57-232 and ``_hensel``, lines 313-637).  It chooses a prime ``P`` of
      the equation order where ``T`` and ``Phi_(p^e)`` are squarefree
      (lines 83-181), takes a root of ``Phi_(p^e)`` in the residue field,
      and lifts it and the inverse of the derivative by Newton steps
      (lines 593-596).  It multiplies by ``T'(theta)`` and reconstructs the
      element by rounding against an LLL-reduced basis of ``P^k``
      (lines 599-611).  The exponent ``k`` is the Friedrich-Fieker bound
      ``_lifting_expo`` (lines 674-721) with root bound 1 at every place,
      since every conjugate of a root of unity has absolute value 1.  Hecke
      stops when a residue field has fewer than ``deg Phi_(p^e)`` roots,
      which proves that ``K`` has none, because ``Q(zeta_(p^e))`` is normal.
      Silex reports failure there.  Hecke lifts every residue root.  Silex
      lifts one, which loses nothing: a normal ``Phi_(p^e)`` with a root in
      ``K`` has all its roots in ``K``, and their reductions are all the
      residue roots.

Silex publishes ``w`` only when it certifies a root of unity whose exact
order is the reduced bound.  The certificate is ``z^w = 1`` and
``z^(w/l) != 1`` for every prime ``l | w``.  With the upper bound, this
proves ``w``.  The certificate is checked in the field itself, after mapping
back from the monic model, so neither step 4a nor step 4b has to be
trusted.  PARI accepts a smaller prime power after a "wrong guess" warning.
Silex fails closed instead.  This happens when the good-prime loop stops
with a bound above ``w``, including ``w = 2`` whenever the reduced bound is
not 2.  It also happens when the Hensel search finds no root, for example
at its limit of 256 primes or when its lifting bound cannot be evaluated.
Fields with a real place and quadratic fields, including those with a
non-monic polynomial, take the direct cases above and do not reach these
steps.
``test/t-unit.cpp`` checks ``Q(zeta_n)`` for ``n = 3..30`` against
``nfrootsof1`` from PARI/GP 2.17.4.  It also checks non-cyclotomic
presentations of ``Q(zeta_n)`` for ``n = 5, 7, 9, 11, 12, 16, 27`` and 32,
where degrees 16 and 18 use step 4b for 2- and 3-power roots.  It checks
non-monic presentations with ``w = 2``, a non-monic cyclotomic translate,
and non-monic presentations with ``w = 10``, 12, 18, 32 and 50.  It checks
the Hensel search on its own, and fail-closed behavior against bounds that
are too large or not multiples of the proven bound.
Earlier versions used only ``sqrt(-1)`` and ``sqrt(-3)`` without comparing
the result to the bound, and returned ``w = 6`` for ``Q(zeta_9)``, where the
true value is 18.

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
* For a non-integral input ``a`` (equivalently ``lc(P_a) > 1``, with ``P_x``
  as in the height filter below), Silex follows Hecke v0.38.6
  ``src/NumField/NfAbs/Elem.jl:is_power`` (lines 674-711): it lifts a root of
  ``y^n = a d^n``, where ``d`` is the power-basis denominator of ``a`` (the
  denominator of ``a`` in ``Z[theta]``, not in a maximal order), and returns
  ``y / d``.  This is needed because the Hensel reconstruction above recovers
  ``f'(theta) c`` as an element of ``Z[theta]``, which holds for every
  integral ``c`` (since ``f'(theta) O_K`` is contained in ``Z[theta]``) but
  fails for most non-integral ``c``; ``d c`` is integral by construction, so
  the rescaled root is reconstructible.  A residue disproof on ``a d^n`` is
  valid for ``a``, since ``a d^n = (d c)^n`` whenever ``a = c^n``.
  **Deviation** (documented at the call site): Silex rescales only when ``a``
  is not an algebraic integer.  Hecke rescales whenever ``d != 1``.  Every
  integral root is already reconstructible from the unscaled input, so
  behavior for integral ``a`` is unchanged either way; Hecke's choice of
  ``d`` from a maximal order, when known, only gives a smaller denominator
  and is not a correctness difference.

Three Silex pre-filters run before the lift or the exact check ``c^n == a``
and only reject, so a rejection leaves the query unsupported and never
changes a definite answer.  None is in the upstream sources: Hecke's
polynomial degree ``n`` is a machine integer bounded by the polynomial it
builds, while Silex accepts ``n`` up to ``2^63 - 1`` and must not form
``c^n``, or the rescaled radicand ``a d^n`` above, blindly.

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
* Denominator (leading-coefficient) filter, ahead of the rescaling above: by
  the same finite-place identity used for the height filter,
  ``log lc(P_x) = sum`` over finite places ``v`` of
  ``d_v log max(1, |x|_v)``, so ``lc(P_{c^n}) = lc(P_c)^n`` exactly whenever
  ``c^n = a``.  A non-integral ``a`` has only non-integral roots, and a
  non-integral ``c`` has ``lc(P_c) >= 2``, so ``lc(P_a)`` must be an exact
  ``n``-th power of an integer ``>= 2``; this also bounds
  ``n < bits(lc(P_a))``, and hence the size of ``d^n``, before it is formed.
  The query is left unsupported when ``lc(P_a)`` is not such a power or when
  ``n >= bits(lc(P_a))``.  This is a stronger, denominator-ideal-norm form of
  the norm filter's reasoning (``lc(P_x)`` is the norm of the denominator
  ideal of ``x``), and it is sound as a disproof: per project decision, it is
  not promoted to a definite ``false``, and the existing residue disproof
  above still runs on the original, unscaled ``a``.

Embedding root contexts
-----------------------

``src/embedding`` stores certified complex roots of the defining polynomial
for archimedean evaluation.  Hecke v0.38.6 ``src/Misc/acb_root_ctx.jl``
(``_roots!``) is the source: it refines by calling FLINT
``acb_poly_find_roots`` with the previous roots as initial approximations,
doubling the working precision, and accepts the result once all roots are
isolated, accurate enough, and pass ``acb_poly_validate_real_roots``.  Silex's
first isolation, and its fallback, use FLINT
``arb_fmpz_poly_complex_roots``.  The signature comes from FLINT
``fmpz_poly_signature`` and is computed when the context is defined.

Silex deviates from the source in five intentional ways:

* No re-sort after refinement.  Hecke re-sorts with FLINT
  ``_acb_vec_sort_pretty``, whose ``acb_cmp_pretty`` comparator orders by
  ``|Im|`` only when the difference is resolved and otherwise by real part,
  so place indices can change with precision.  Silex matches each refined
  ball to the previous ball it overlaps and requires the match to be
  one-to-one.  Because both vectors isolate the roots with pairwise disjoint
  balls, a new ball that overlaps exactly one previous ball encloses the same
  root.  Place order is therefore fixed once roots are set, which cached
  logarithmic embeddings rely on.
* A guarded starting precision.  Hecke's first refinement attempt works at
  the requested precision.  Root finding there returns balls a few bits short
  of the requested accuracy, so that attempt fails and is repeated at double
  the precision.  Silex starts 32 guard bits above the requested precision
  and keeps the same acceptance checks and the doubling.  This changes only
  the precision schedule, not what is verified.
* A precision cap on refinement from previous roots.  Hecke doubles up to
  ``2^22`` bits and then raises an error.  Silex stops that loop at eight times
  its starting precision and falls back to FLINT's certified isolation.  The
  factor is a tuning constant, not a correctness condition.
* No refinement of roots that already meet the request.  When every stored
  ball already has the requested relative accuracy, Silex keeps the stored
  balls: they passed the isolation, real-root validation, and place-order
  checks when they were stored, and the accuracy check depends only on the
  balls.  Refining them anyway is slow with FLINT 3.6.0, whose
  ``acb_poly_find_roots`` then runs its full iteration budget because every
  correction is zero and its early-stop test never fires.
* A retrying fallback.  When a full isolation does not yet match the previous
  balls one-to-one, Silex isolates again at double the precision.  This
  terminates because each root is at positive distance from the closed
  previous balls that do not contain it.  Refinement fails only when the
  working precision would overflow, for a non-squarefree polynomial, or for
  invalid input.

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

Residue reduction follows PARI 2.17.3 ``src/basemath/base2.c:modprinit``
in the branch where ``p`` does not divide the index: an order element is
mapped to its ``alpha``-polynomial through the order basis (PARI
``nf_get_zkprimpart`` and ``nf_get_zkden``) and reduced modulo
``(p, g)``.  Silex does not implement PARI's anti-uniformizer handling of
denominators divisible by ``p`` (``Rg_to_ff``); such inputs fail.  PARI
takes that polynomial branch only for residue degree ``f > 1``; for
``f = 1`` it projects with ``dim1proj`` on the prime's HNF (``base2.c``
lines 2520--2526), whereas Silex uses the polynomial reduction for every
``f``.  PARI's ``nf`` is always the maximal order, so applying the same
reduction on an equation order is Silex's own extension, justified by
``Z[alpha]/(p, g(alpha)) = F_p[x]/(g)`` for an irreducible factor ``g`` of
the defining polynomial modulo ``p``.  On a
maximal quadratic-backend order the residue polynomial is instead taken in
the integral generator ``omega`` of PARI ``quadgen``/``quadpoly``, whose
``[1, omega]`` coordinates are already the ``omega``-polynomial.  The
generator in use is chosen by one shared predicate when a prime is built and
stored in the prime, so it cannot drift from the stored residue polynomial.

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

The ideal-lattice LLL step of relation search follows Hecke 0.39.19
``src/NumFieldOrd/NfOrd/LLL.jl:_lll``, which rounds the scaled Minkowski Gram
matrix from ``NfOrd.jl:minkowski_gram_mat_scaled``.  Hecke computes that
matrix from ``minkowski_matrix(O, prec + 2)``, whose entries meet the
absolute tolerance ``2^-(prec + 2)``.  Silex embeds at the relative precision
``prec + 64`` and checks that every scaled entry rounds with radius at most
``0.1``.  On high-degree power bases (for example the real cyclotomic fields
of conductor 37 and 41) the large entries miss that check at every
precision, so doubling ``prec`` never succeeds.  As an intentional deviation
in the precision schedule only, when the rounding check fails Silex keeps
``prec`` and raises the embedding precision by the measured shortfall,
again while the shortfall at least halves, before it falls back to
doubling.  A raise is kept when ``prec`` doubles.  The same schedule applies
to the weighted reduction.  Hecke's weighted branch of ``_lll``
(``LLL.jl:104-128``) differs from its unweighted one: it rounds
``minkowski_matrix(K, prec)`` (``src/conjugates.jl:105-158``), BigFloat
values at relative precision ``prec``, with ``round_scale!``
(``Misc/Matrix.jl:405``) and no radius check.  Silex's radius check on that
path is therefore already stricter than Hecke's, and the raise moves Silex
toward Hecke's behavior there of rounding at the first precision.
Attempts that pass the rounding check are unchanged, and every relation is
still verified exactly.

Adjoining a dependent unit ``y`` to a full-rank unit group follows Hecke
v0.38.6 ``src/NumFieldOrd/NfOrd/Unit/UnitGrpCtx.jl:_add_dependent_unit!`` and
``Unit/Relation.jl:_find_rational_relation!``.  Silex computes the
coordinates of ``y`` against the generators' log matrix in Arb, replaces each
coordinate ball by the simplest rational inside it (FLINT
``_fmpq_simplest_between`` on the exact ball endpoints, as Nemo
``simplest_rational_inside``), and takes the least common multiple ``m`` of
the denominators.  A ball of radius above one, or a denominator or ``m``
above the regulator index bound ``B``, makes the attempt inconclusive.
Silex accepts each denominator ``q`` in the closed range ``q <= B``.
Hecke's ``simplest_inside`` requires each denominator below the bound, while
its lcm test accepts ``m`` equal to it.  The closed range is the one the
bound proves: with ``G`` the group generated so far and ``U`` the unit group,
every denominator divides the least ``m`` with ``y^m`` in ``G`` modulo
torsion, and ``m <= [U:G] <= B``, so ``m = B`` can occur.  The candidate
relation is only a guess until it passes the exact power-root or torsion
check.  When adjoining ``y``, Silex decides the candidate as Hecke does, with
the certified log screen of ``Unit/Relation.jl:_check_relation_mod_torsion``
and adjoins ``y`` itself.  When the screen is inconclusive, Hecke raises a
precision error; Silex instead adds its own exact-expansion fallback to the
check.  It extracts an ``m``-th root only when that check cannot complete.  An inconclusive attempt or a failed check doubles the
precision.  Hecke's loop is unbounded; Silex stops at a fixed precision cap
and then leaves the group unchanged.  The loop terminates below the cap in
the usual case: when the true coordinate is ``a/b`` with ``b`` at most the
bound ``B`` and the ball radius is below ``1/(2B^2)``, the ball has diameter
below ``1/B^2``, while every other fraction ``p/q`` with ``q <= b`` lies at
distance at least ``1/(qb) >= 1/B^2`` from ``a/b``, so the simplest rational
in the ball is ``a/b``.

Reducing relation units modulo the unit group follows Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Map.jl:reduce_mod_units``: round the Arb
coordinates of each unit against the generators' log matrix and divide by
the rounded product, doubling the precision while the rounding fails.  The
rounding differs at a half-integer.  Hecke's ``round(ZZMatrix, C)`` calls
Nemo ``src/HeckeMoreStuff.jl:Base.round(::Type{ZZRingElem}, ::ArbFieldElem)``,
which accepts any coordinate ball of radius at most 0.1 and rounds its
midpoint.  Silex needs a certified nearest integer, so an exact half-integer
coordinate never rounds.  Silex stops doubling and leaves the unit unreduced
once some coordinate ball has radius at most 0.1 and contains a
half-integer, and every coordinate ball that defeats the approximate
rounding (for example one with 0 strictly inside it, see below) also has
radius at most 0.1.
When those balls hold an exact tie and exact zeros, this is the result the
precision ladder reached at the cap.  It can differ from the cap result in a
near-tie, or beside a small nonzero coordinate whose narrow ball still
contains 0, where more precision would have let the rounding succeed.

The approximate rounding itself ports Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Map.jl:round_approx(ZZMatrix, ArbMatrix)``.  Its
width check forms ``lower - upper`` (the floor of the lower bound minus the
ceiling of the upper bound) and compares it with ``sqrt(abs(x))``.  That
difference is never positive, so the check never rejects; Hecke has the same
sign and Silex keeps it.  A wide ball rounds to ``div(upper + lower, 2)``.

**Deviation** (documented in the code): Silex also fails the rounding
when the ball or ``sqrt(abs(x))`` is not finite; Hecke has no such test.
FLINT ``arb_abs`` keeps the radius, and ``arb_sqrt`` returns an indeterminate
ball when that ball contains negative numbers.  So Silex fails for a ball
with 0 strictly inside it, such as ``0 +/- 1/2``.  Hecke's ``>`` reaches
FLINT ``arb_lt`` through Nemo, which is false for an indeterminate ball, so
Hecke returns 0 there.  An exact 0, or a ball such as ``1/2 +/- 1/2`` whose endpoint is 0,
rounds in both.  The half-integer stop rule above relies on this failure.

Saturation-backed class-group proof follows Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Proof.jl:_class_group_proof``: verify generation
up to the proven generation requirement (the Minkowski-type bound, or
Zimmert's twin-class bound where it applies; see below), then saturate at every prime dividing the
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
A prime ``ell`` whose ``ell``-local target rank is zero (``ell`` divides
neither ``h_cand`` nor the computed torsion order ``w``, and the unit rank is
zero) needs no discrete-logarithm check: the index of the computed relation
lattice in the full one divides ``h_cand``, and the unit index divides ``w``,
so the ``ell``-part of the S-unit index is trivial.  Silex records a verified
``ell``-local proof with rank 0 and target 0 for such a prime.  This rests on
that index argument, not on a matching upstream routine.
Gates that accept a supplied index bound or analytic class-regulator product
are not installed; public certification uses only values Silex computes or
verifies itself, and proven unit publication through the class/unit gates
re-derives the torsion subgroup of the order rather than trusting a stored
one.

``factor_base_class_group_bound`` returns the generation bound for a
maximal order: every ideal class contains an integral ideal of norm at most the
bound, so verifying that the prime ideals up to it lie in the factor base
proves generation.  It is 1 in degree one, ``floor(sqrt(|d|/3))`` for
imaginary quadratic and ``floor(sqrt(d)/2)`` for real quadratic fields.  In
degree ``n = r1 + 2 r2 >= 3`` it is the smallest of three proven bounds,
which is valid because each of them is:

* H. Zimmert, "Ideale kleiner Norm in Idealklassen und eine
  Regulatorabschätzung", *Invent. Math.* 62 (1981), 367--380, Satz 2
  (p. 372), used for ``n <= 20``: every ideal class contains an integral
  ideal ``a`` with ``log(sqrt|d| / N(a)) >= L``, where
  ``L = r1 (-psi((1+gamma)/2) - log Gamma(1/2+gamma) + log Gamma(1+gamma) +
  (1/2) log pi) + r2 (-2 psi(1+gamma) + 2 log 2 + log(1/2+gamma) + log pi) -
  2/(gamma-alpha) - log[(1+1/alpha) (1+1/gamma)^(-2)
  (1+1/(2 gamma-alpha))^(-1)]`` for every ``gamma > alpha > 0``.  Silex
  takes ``alpha = gamma - gamma(gamma+1)/sqrt(1+3 gamma+3 gamma^2)`` from
  Bemerkung 1 (p. 373) and, for each signature, a fixed ``gamma`` in
  hundredths that maximizes ``L`` over that grid; where Tabelle 1 (p. 368)
  prints a ``gamma``, the tabulated value equals it, and ``exp(L)``
  reproduces the printed constants.  The choice of ``gamma`` affects only
  the size of the bound, not its validity.  Silex evaluates ``L`` with Arb
  at 128 bits and uses the ceiling of the upper endpoint of
  ``sqrt|d| exp(-L)``.
* Minkowski's bound ``n! n^(-n) (4/pi)^r2 sqrt|d|`` (S. Lang, *Algebraic
  Number Theory*, p. 119, Theorem 4, as quoted by Zimmert 1981, p. 367),
  evaluated the same way.
* The exact integer form ``max(1, ceil(n! ceil(sqrt|d|) 2^r2 / n^n))``, which
  replaces ``4/pi`` with the larger ``2``.  This was the whole bound before
  Zimmert's bound was adopted; it remains the fallback if an Arb evaluation
  does not give a finite enclosure, and it guarantees that the bound never
  exceeds its former value.

Neither of the first two dominates: for ``(r1, r2) = (1, 1)`` Minkowski's
constant ``27/(6 (4/pi))`` exceeds Zimmert's ``exp(L)``, so Minkowski's bound
is the smaller one there.  For larger degrees Zimmert's bound is much smaller;
for ``x^19 - x - 1`` it is 6135 against the former 43837247.
``test/t-generation-bound.cpp`` pins Zimmert's bound against an independent
evaluation of Satz 2 and the combined bound for a corpus of fields.
``factor_base_class_group_bound`` promises a small integral ideal in every
class; it does not use Zimmert's twin-class bound below, which promises less.

A ``proven`` class/unit computation verifies generation up to the *proven
generation requirement*, which is ``factor_base_class_group_bound`` unless
Zimmert's twin-class bound applies and is smaller.  Zimmert 1981, Korollar 1
(p. 376), a consequence of Satz 4, says: for every ideal class ``K`` there
are integral ideals ``a`` in ``K`` and ``b`` in ``[D] K^(-1)``, where ``[D]``
is the class of the different, with
``log(|d| / (N(a) N(b))) >= F_{r1+r2, r2}(gamma) + n log pi`` for every
``gamma > 0``.  Hence ``K`` or ``[D] K^(-1)`` contains an integral ideal of
norm at most ``sqrt|d| exp(-(F + n log pi)/2)``.  Silex does not evaluate
the series ``F_{a,b}`` of p. 376.  The proof of Satz 4 (pp. 379--380) applies
the Hilfssatz (p. 378) to ``H = -Q_a^a Q_b^b U V``, where ``Q_a`` and ``Q_b``
are products of ``m`` factors, and for every ``m >= 1`` it establishes the
inequality with the left side ``G_{m,j}(gamma)`` of the Hilfssatz in place of
``F``, for ``j`` the smallest integer with ``j >= 1/2 + gamma``; ``F`` is the
limit as ``m`` tends to infinity.  Silex evaluates that finite-``m``
inequality itself, with ``m = 128``, so it needs no estimate of a series
tail.  With ``a = r1 + r2``, ``b = r2``, ``D = 2 (1 + 2 gamma)``,
``y_l = (2l - 1 + gamma)/D``, ``y'_l = (2l + gamma)/D`` and
``delta(y) = psi(y + 1/2) - psi(y)``, the logarithmic derivatives of
Zimmert's definitions give::

    G = a [sum_{l=1}^{m} ((4/D) delta(y_l) - 1/(2l-2-gamma) - 1/(2l-1+gamma))
           + sum_{l=m+1}^{m+j} ((2/D) delta(y_l) - 1/(2l-2-gamma))]
      + b [sum_{l=1}^{m} ((4/D) delta(y'_l) - 1/(2l-1-gamma) - 1/(2l+gamma))
           + sum_{l=m+1}^{m+j} ((2/D) delta(y'_l) - 1/(2l-1-gamma))]
      + (4/D) [psi(gamma/D) + psi((1+gamma)/D)
               - psi((1+3 gamma)/D) - psi((2+3 gamma)/D)]
      - (a/2) [psi(1/2 + gamma/2) + psi(-gamma/2)]
      - (b/2) [psi(1 + gamma/2) + psi(1/2 - gamma/2)]

Silex evaluates this sum of digamma values with Arb at 128 bits, using
``psi(x) = psi(x + k) - sum_{i<k} 1/(x + i)`` at negative arguments, but
only offline: the constant depends on the signature alone, so Silex stores a
precomputed, test-verified table.  For each signature with ``3 <= n <= 20``,
the generator ``tools/zimmert_twin_class_table.cpp`` (target
``silex-zimmert-twin-class-table``) takes the upper endpoint of the Arb
enclosure of ``exp(-(G + n log pi)/2)`` and rounds it up to an exact decimal
rational ``E`` with 30 significant digits, written to
``src/factor_base/zimmert_twin_class_table.inc``; to regenerate, build that
target and redirect its output to that file.  At run time the twin-class
bound ``B`` is ``ceil(E sqrt|d|)``, computed exactly with integers as
``ceil(ceil(s sqrt|d|) / 10^k)`` for ``E = s / 10^k``.
``test/t-generation-bound.cpp`` recomputes every entry with the Arb
evaluation and checks that it is on the safe side (``E`` at least the upper
endpoint ``U``, so ``B`` is at least Zimmert's ``sqrt|d| exp(-c)``) and
tight (``E <= U (1 + 10^-25)``).  PARI's ``zimmertbound`` (below) likewise
hard-codes a constant ``c`` per ``(n, r2)``, but as a C ``double`` given to
about nine digits with no stated rounding direction, and evaluates
``exp(-c) sqrt|d|`` in floating point.  Zimmert's terms are singular at integer
``gamma``; every other ``gamma > 0`` gives a valid bound, so ``gamma``
affects only its size.  For each signature with ``3 <= n <= 20`` Silex fixes
the ``gamma`` in hundredths, ``0 < gamma < 3`` and not an integer, that
maximizes ``G``.  ``G_{m,j}`` converges to ``F`` as ``m`` grows; at
``m = 128``, ``F - G < 3e-4`` at these ``gamma``.  The validity of the
bound does not depend on how close ``G`` is to ``F``.
``test/t-generation-bound.cpp`` checks the constant
``exp((G + n log pi)/2)`` against all 16 rows of Tabelle 3 (p. 377) with
``3 <= n <= 20`` (the printed value is a lower bound within one unit of its
last digit) and pins ``B`` at ``|d| = 10^30``;
``test/zimmert_twin_class_oracle.py`` is the independent mpmath evaluation
those values come from, and it also checks the closed form above against
Zimmert's definitions of ``Q_a``, ``Q_b``, ``U`` and ``V``.

The twin-class bound proves generation only together with the class of the
different: if every prime ideal of norm at most ``R >= B`` lies in the
subgroup ``H`` generated by the factor base and ``[D]`` lies in ``H``, then
every class lies in ``H``, since ``K`` or ``[D] K^(-1)`` does.  Silex settles
``[D]`` exactly, only for a monic defining polynomial ``f`` with integer
coefficients and root ``theta``, in one of two ways:

* ``O_K = Z[theta]`` (index one): the conductor ``F`` below is ``O_K``
  (W. Narkiewicz, *Elementary and Analytic Theory of Algebraic Numbers*,
  3rd ed., Springer, 2004, Corollary to Proposition 4.18, p. 153), so the
  different is the principal ideal ``f'(theta) O_K``, ``[D]`` is trivial
  and the requirement is ``B``.
* ``[O_K : Z[theta]] > 1``: ``f'(theta) O_K = F D`` for the conductor
  ``F = {x in K : x O_K in Z[theta]}``, so ``[D] = [F]^(-1)``.  The index
  lies in ``F``, so every prime ideal dividing ``F`` lies above a prime
  ``p`` dividing the index.  The requirement is the larger of ``B`` and the
  largest norm of a prime ideal of ``O_K`` above such a ``p``, which puts
  ``[F]``, and so ``[D]``, in ``H``.  For ``f'(theta) O_K = F D`` Silex cites
  Narkiewicz 2004, Proposition 4.18(i) (p. 152), which states
  ``f_A = delta_{L/K}(a) D_{L/K}^(-1)`` for the conductor ``f_A`` of
  ``A = R[a]``, where ``a`` is an integral element generating ``L/K`` and
  ``delta_{L/K}(a) = f'(a)``, with ``K = Q``, ``R = Z``, ``L`` the number
  field and ``a = theta``; its proof rests on Proposition 4.17 (p. 151),
  ``R[a]^* = f'(a)^(-1) R[a]``.  H. Cohen, *Advanced Topics in
  Computational Number Theory*, GTM 193, Lemma 10.1.20 and the proof of
  Corollary 10.1.21 (pp. 484--485), give the same argument in a relative
  setting.

Silex uses the twin-class route only when its requirement is strictly below
``factor_base_class_group_bound``.  Otherwise, and whenever the degree is
outside ``3..20``, the defining polynomial is not monic with integer
coefficients, or the index or a prime decomposition cannot be computed, the
requirement is ``factor_base_class_group_bound`` (Satz 2, Minkowski or the
integer form above) and the generation basis is ``minkowski_type``.  Silex
does not settle ``[D]`` from relations, so a field whose index primes have
large norm keeps that bound.  With the twin-class route the generation basis
is ``zimmert_twin_class`` and the generation certification bound is the
requirement; the detail log names the route (``index_one`` or
``conductor_primes``), and ``silex-class-unit-instance`` prints it as
``generation_different_route`` in its ``class_group`` object, null for any
other basis.  The same requirement applies wherever unconditional
generation is recorded, including by a ``grh`` run whose base reaches it,
and in the relation-search honesty check.

PARI 2.17.3 ``src/basemath/buch3.c:zimmertbound``, the bound of
``bnfcertify``, uses a table of constants that are larger than the optimum
of Satz 2 for every signature.  The table follows Zimmert's Korollar 1,
and ``bnfcertify0`` accordingly expresses the class of the different over
the computed class group, through ``bnfisprincipal0(bnf, D, nf_FORCE)``,
whenever the polynomial order is not maximal; it discards the result, so
this is not a principality test.  Earlier releases make the same call (for
example 2.7.6 ``src/basemath/buch2.c:testprimes``); it moved to
``bnfcertify0`` in 2.9.  Silex does not copy the table: for every
signature with ``5 <= n <= 20`` except ``(n, r2) = (5, 1)``, ``(5, 2)`` and
``(6, 3)``, its constant exceeds the optimum of Korollar 1 over ``gamma``,
so PARI's bound is smaller than any bound Korollar 1 proves, by up to about
0.33% of the bound at ``n = 20``.  This comparison was computed with
non-interval mpmath at 30 digits, not with Arb, so the small excesses for
``5 <= n <= 7`` are not certified; from ``n = 8`` up the excess is clearly
beyond floating-point error.

The analytic index bound is the ceiling of the rigorous upper endpoint of the
quotient ``h_cand R_cand / hR``.  When the analytic ``hR`` is correct and the
factor base generates the class group, this quotient is the combined
class/unit index, a positive integer.  If the upper endpoint is below one, the
enclosure contains no positive integer: the analytic value or the candidate is
wrong, or the factor base does not generate the class group (in general the
quotient is ``[Lambda_S : L] / [Cl : <S>]``, which can lie below one when the
factor base ``S`` does not generate).  In every such case the bound fails
closed, and the check reports the analytic value as unavailable instead of
reporting index one.  The public accessor
``OrderUnitGroup::class_regulator_index_bound`` does not itself check
generation, so it fails closed for any upper endpoint below one.  Only once
factor-base generation is proven is the quotient a positive integer.

This is a deliberate deviation from Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp.jl:_validate_class_unit_group``, which in this
case returns index one (through ``abs_upper_bound`` of a value in ``(0, 1)``)
and so accepts.  Silex fails closed instead (a user decision, 2026-09-26).
PARI 2.17.3 ``src/basemath/buch2.c:bad_check`` also rejects a quotient below
0.75, treating it as a precision failure, though it then retries with more
precision where Silex fails.

The exact imaginary-quadratic ``proven`` route uses the same index argument
with the exact class number ``h = h(D)``, counted from the reduced forms
returned by FLINT ``qfb_reduced_forms``, in place of ``ell``-local tests.
Here ``D`` is the fundamental discriminant of the maximal order (the route
applies to maximal orders only), so every form of discriminant ``D`` is
primitive and the count of reduced forms is the class number of the order.
With generation verified up to the Minkowski-type bound and ``h_cand = h``,
the computed relation lattice has index one in the full one, so it is
saturated at every prime.  The route stores a verified relation-saturation
record for every prime ``ell`` dividing ``h`` (none when ``h = 1``), backed
internally by a proof record that names the exact class number, not an
``ell``-local test, as its basis.  A prime that already carries a verified
``ell``-local proof keeps it.  Later gates, such as the Belabas--Friedman
audit, therefore see the same per-prime saturation records as for a
saturation-backed proof and can promote from them.

An analytic ``hR`` proves a class/unit pair only when it is unconditional.
Silex has two unconditional routes: degree one, where the residue is exactly
one, and maximal orders of explicit quadratic-backend fields whose
discriminant satisfies ``|D| < 2^44`` (about ``1.76e13``), where the residue
comes from ``L(1, chi_D)``.  The size cap bounds the cost of the
evaluation, whose number of terms grows like ``sqrt(|D| prec)``.  The cap
itself does not depend on the precision, and fields with ``|D| >= 2^44`` use
the Belabas--Friedman fallback.  Silex evaluates ``L(1, chi_D)`` for the
primitive real character ``chi_D(n) = (D/n)`` (Kronecker symbol) of the
fundamental discriminant ``D`` with the approximate functional equation of FLINT
``acb_dirichlet_l_fmpq_afe`` (FLINT 3.6.0 ``src/acb_dirichlet/l_fmpq_afe.c``,
copyright (C) 2021 Fredrik Johansson, LGPL-3.0-or-later; the adapted code is
in ``src/zeta/zeta.cpp``).  It keeps FLINT's two sums, truncation rule, tail
bound, tolerances, and working precisions, and the result is a rigorous ball
enclosure.  It deviates from FLINT in three places: ``chi_D(n)`` is the
Kronecker symbol instead of a Dirichlet-group discrete logarithm, so no
Dirichlet group is built and the work grows like ``sqrt(|D| prec)``; the root
number is ``W(chi_D) = 1``, because the Gauss sum of a primitive real
character is ``sqrt(q)`` for even and ``i sqrt(q)`` for odd characters
(Davenport, *Multiplicative Number Theory*, Chapters 2 and 9); and the
incomplete gamma values are the closed forms ``Gamma(1, z) = exp(-z)``,
``Gamma(1/2, z) = sqrt(pi) erfc(sqrt(z))`` and ``Gamma(0, z) = E_1(z)``
(DLMF 6.2 and 8.4), evaluated with arb ball functions at block anchors and
obtained by Taylor stepping between them, where FLINT evaluates each term
from scratch.

The stepping is Silex's own evaluator and uses no FLINT code.  Its closest
upstream analogue is in the same FLINT file: the term loop there has a
disabled branch that would continue each incomplete gamma value from the
previous term with ``_arb_gamma_upper_fmpq_step_bsplit`` (a source comment
says evaluating from scratch "seems to be faster"), and FLINT uses that
routine for the bit-burst refinement within each term.  That routine takes a
single step between two points, with its own number of terms, and sums the
Taylor series by binary splitting.  Silex instead fixes one polynomial per
block of consecutive indices: the coefficients are computed once at the
block anchor, each value in the block is a Horner evaluation, and one
truncation bound at the largest step covers the whole block.  Write
``a = (s + parity) / 2`` in ``{0, 1/2, 1}``, ``z_n = pi n^2 / q``, and for an
anchor ``z = z_{n0}`` put ``K = exp(-z) z^(a-1)`` and ``v = z_n - z =
(pi / q) (n - n0) (n + n0) >= 0``.  Substituting ``u = z + w`` in
``Gamma(a, z) - Gamma(a, z + v) = int_z^(z+v) u^(a-1) exp(-u) du`` and
expanding ``exp(-w) (1 + w / z)^(a-1)`` with the binomial series (DLMF 4.6.7,
absolutely convergent for ``0 <= w < z``) gives, for ``0 <= v < z`` (any
``v >= 0`` when ``a = 1``)::

    Gamma(a, z + v) = Gamma(a, z) + K sum_{m >= 0} p_m (-v)^(m+1) / (m + 1),
    p_m = sum_{j=0..m} |binom(a - 1, j)| z^(-j) / (m - j)!,

because ``binom(a - 1, j)`` is ``(-1)^j`` for ``a = 0``, ``(-1)^j binom(2j, j)
/ 4^j`` for ``a = 1/2``, and zero for ``j >= 1`` when ``a = 1``.  In all three
cases ``|binom(a - 1, j)| <= 1``, so with ``rho = max(1, 1 / z)`` (``rho = 1``
for ``a = 1``) one has ``0 <= p_m < e rho^m``, and for ``rho v < 1`` (which
gives ``v < z`` when ``a < 1``) the terms with ``m >= M`` sum to at most
``|K| e v (rho v)^M / (1 - rho v)`` (using ``1 / (m + 1) <= 1`` and the
geometric series).  A block anchored at ``n0`` covers ``n0
<= n <= n0 + B``, so ``v <= v_max = (pi / q) B (2 n0 + B)``.  Within it each
value is the anchor ball plus the degree-``M`` polynomial in ``-v`` evaluated
in ball arithmetic, with the truncation bound at ``v_max`` added to the
radius; that bound is computed in upward-rounded ``mag`` arithmetic from an
upper bound of ``|K|``, of ``v_max`` and of ``rho``, and a block is used only
if these bounds give ``rho v_max < 1``.  The first anchor is a direct
evaluation and each later anchor is the stepped value at the end of the
previous block, so by induction every value is a ball containing the exact
``Gamma(a, z_n)``.  ``B`` (so that ``rho v_max`` is about ``2^-10``), ``M``,
and the extra precision carried by the anchors are chosen in floating point
and affect only the width and the cost, not the validity, of the enclosure;
where no block of at least four steps fits, or a block would need more than
32 terms (as at 512-bit precision, where stepping was measured not to pay),
the value is evaluated directly.

When the ``L``-value ball is not finite and positive, the quadratic route
falls back to Belabas--Friedman, and the value is then GRH-conditional like
any other; certification uses the route that actually produced the value, not
the field type.  Every other ``hR`` comes from the Belabas--Friedman
evaluation (K. Belabas and E. Friedman, "Computing the residue of the
Dedekind zeta function", *Math. Comp.* 84 (2015), 357--369, Theorem 1), whose
truncation-error bound assumes GRH.  Such an ``hR`` is therefore at most
``grh`` evidence: the analytic check is recorded with
``analytic_class_regulator_certification() == grh``, it never
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
the unconditional label.  A Belabas--Friedman evaluation never serves as a
proof component; in degree one the Belabas--Friedman entry points use the exact
value ``hR = 1`` and evaluate no Belabas--Friedman series.  Outside degree one
a Belabas--Friedman check is recorded against unproven units only as the
``grh``-mode acceptance record described below.  Because the Belabas--Friedman
value no longer contributes to a degree-three-or-higher ``proven`` result, such
results report ``analytic_class_regulator_status()`` and
``zeta_bf_proof_status()`` as ``not_checked`` unless an audit is requested
separately; these statuses are audit records, not proof components.

A ``grh`` request for a paired class/unit transaction sizes the factor base
with the GRH bound, the minimum of the Belabas--Diaz y Diaz--Friedman
criterion and Bach's bound, as Hecke v0.38.6
``src/NumFieldOrd/NfOrd/FactorBaseBound.jl:factor_base_bound_grh`` does.
The Belabas--Diaz y Diaz--Friedman bound is the smallest integer cutoff at
which an Arb comparison proves the criterion of Belabas--Diaz y
Diaz--Friedman 2008, Section 3.  Silex finds it as PARI 2.17.3
``src/basemath/buch2.c:Buchall_param`` (lines 3814-3821) finds its
``LIMC2``: doubling the cutoff from 2 until the criterion holds, then an
exact bisection on the integers.  Unlike PARI, which compares doubles with a
small slack, every comparison is an interval comparison, and a cutoff where
the bisection cannot decide the comparison counts as one where the criterion
is not proven, so the result is always a proven cutoff.  Hecke's search
instead doubles from 100 and stops its dichotomy at a step below 100, so its
bound is never below about 50.  When the Belabas--Diaz y Diaz--Friedman
criterion cannot be computed or decided, the GRH bound is Bach's bound
alone, and the generation record names Bach's theorem.  This happens for a
non-monic defining polynomial, since the criterion's residue-degree route
needs a monic integral one, for a residue-degree failure, when an interval
comparison in the doubling cannot be decided, and when the doubling passes
``2^52``.  The factor base is then built to that GRH bound, doubled until
the base has more than 20 prime ideals.  Bach's theorem
does not depend on the criterion, and Belabas--Diaz y Diaz--Friedman 2008,
Section 3, itself takes Bach's bound once its search passes it.  The search also
stops, as an intentional deviation from that rule, once the criterion fails
at a cutoff at or above Bach's bound; since the search statistic is
non-decreasing in the cutoff, this leaves the minimum unchanged.  A ``grh``
request on a non-monic presentation still fails closed today, because prime
decomposition, and so the factor base, also needs a monic integral defining
polynomial.
In degree three and higher and for real quadratic fields, factor-base
generation may rest on GRH.  The grh base contains every prime ideal of norm
at most its build bound, inert primes and both primes of a split pair
included, and that build bound is at least the undoubled GRH bound.  Under
GRH such a base generates the class group: Belabas--Diaz y Diaz--Friedman
2008, Corollary 2.2, assumes the Riemann hypothesis for the ``L``-functions
of the nontrivial class-group characters; Bach 1990, Theorem 4 (with
conductor one, ``12 log^2 |D|``) and the remark following it (p. 376,
``6 log^2 |D|`` in degree two), assumes it for ``zeta_K`` as well.  Silex
records this as a separate GRH generation record on the factor base, written
only when the grh policy builds the base.  It never sets
``factor_base_generation_status()`` to ``verified`` and is not the
relation-search honesty receipt, so public generation readers and every
``proven`` route still see only unconditional generation, and a ``proven``
promotion still requires generation verified up to the proven generation
requirement (``factor_base_class_group_bound``, or Zimmert's twin-class
bound where it applies).  When the GRH-sized base reaches that
requirement anyway, generation is recorded unconditionally as before.  If the
GRH generation record cannot be written, the ``grh`` transaction fails
closed.  As in PARI 2.17.3
``src/basemath/buch2.c:Buchall_param``, where the primes up to the
``GRHchk`` bound generate the class group under GRH, and Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp.jl:_class_unit_group``, fields such as
``x^3 + x + 200`` and ``x^2 - 100003``, whose Minkowski-type bound is far
above the GRH bound, are therefore accepted in ``grh`` mode.  Relation and
unit completeness then rest on the analytic index-one test: an enclosure of
``h_cand R_cand / hR`` whose upper endpoint is below two, as in Hecke
``_validate_class_unit_group`` and PARI's ``bad_check``.  With a
Belabas--Friedman ``hR`` that test holds under GRH (Belabas--Friedman 2015,
Theorem 1).  With the quadratic ``L(1, chi_D)`` value the test itself is
unconditional, but the result still rests on GRH through generation unless
the base reached the Minkowski-type bound.  Such an acceptance is not a
relation-saturation proof: ``relation_saturation_status()`` is not set to
``verified`` by it.
The pair is published with both labels ``grh``, and the analytic check that
accepted it is recorded on the published class group with its own
conditionality.  ``analytic_class_regulator_certification()`` is ``grh`` for
a Belabas--Friedman ``hR`` and ``proven`` for the unconditional quadratic
``L(1, chi)`` value.  With a Belabas--Friedman ``hR`` (in degree three and
higher, the default zeta route's Belabas--Friedman evaluation) the BF audit
data of that evaluation is recorded as well: its error bound, cutoff,
maximum cutoff, and precisions (``zeta_bf_proof_record()``).  The record is
informational.  The ``grh`` transaction never changes or promotes the
class-group or unit labels, and the record does not change the public zeta
and BF gates, which still fail closed in degree three and higher unless
saturation is proven.  There is one exception to "the record never promotes
labels": a ``proven`` record stored by a real-quadratic ``grh`` run can later
let an explicit ``try_certify_with_units(units, proven)`` with proven units
promote the class group to ``proven``.  This is sound, because the index-one
test against the unconditional ``L(1, chi_D)`` value, with generation verified
up to the Minkowski-type bound, gives ``h_cand = h`` and ``R_cand = R_K``
unconditionally.  A ``grh`` run whose base covers only the GRH bound has no
unconditional generation, so that promotion fails for it.  The exact
imaginary-quadratic ``grh`` route uses an exact class number instead of an
analytic value.  Its GRH dependence is in factor-base generation, which is
checked only up to the GRH bound, and it records no analytic check.  The
degree-one ``grh`` route uses the exact degree-one route described below and
records none either.

Each published class group reports where GRH entered, component by
component, through ``factor_base_generation_certification()``,
``factor_base_generation_basis()``,
``factor_base_generation_certification_bound()``,
``class_unit_completeness_certification()`` and
``class_unit_completeness_basis()``; ``silex-class-unit-instance`` prints the
same values in its ``class_group`` object.  Generation is ``proven``, with
basis ``minkowski_type`` or ``zimmert_twin_class`` and the proven generation
requirement, when ``factor_base_generation_status()`` is ``verified``; it is ``grh`` when only
the GRH generation record above holds, with basis ``bdf`` or ``bach`` for the
theorem that gave the GRH bound and that bound before the factor-base size
doubling.  A ``grh`` generation basis ``bdf`` assumes the Riemann hypothesis
for the ``L``-functions of the nontrivial class-group characters
(Belabas--Diaz y Diaz--Friedman 2008, Corollary 2.2); ``bach`` assumes it for
``zeta_K`` as well (Bach 1990, Theorem 4).  Completeness of the relations and
units is stated relative to generation, because the index argument needs a
generating factor base.  It is ``grh``, basis ``belabas_friedman``, when the
pair was accepted against a Belabas--Friedman ``hR``, whose error bound
assumes the Riemann hypothesis for ``zeta_K`` and for ``zeta_Q``
(Belabas--Friedman 2015, Theorem 1).  It is ``proven`` when the pair was
accepted against the unconditional quadratic ``L(1, chi_D)`` value
(``unconditional_analytic``) or the exact imaginary-quadratic class number
(``exact_class_number``), and for a ``proven`` class group or the exact
degree-one route (``unconditional_certification``).
``certification_status()`` stays the requested label: a ``grh`` result whose
two components are both ``proven`` still reads ``grh``.  Per route:

.. list-table::
   :header-rows: 1

   * - ``grh`` route
     - Generation
     - Completeness
   * - Degree one
     - ``proven`` (``minkowski_type``)
     - ``proven`` (``unconditional_certification``)
   * - Imaginary quadratic, exact class number
     - ``grh``, or ``proven`` when the base reaches the proven bound
     - ``proven`` (``exact_class_number``)
   * - Real quadratic, ``L(1, chi_D)``
     - ``grh`` or ``proven``
     - ``proven`` (``unconditional_analytic``)
   * - Real quadratic, ``L(1, chi_D)`` failed or ``|D| >= 2^44``,
       Belabas--Friedman fallback
     - ``grh`` or ``proven``
     - ``grh`` (``belabas_friedman``)
   * - Degree two outside the explicit quadratic backend
     - ``grh`` or ``proven``
     - ``grh`` (``belabas_friedman``)
   * - Degree three and higher
     - ``grh`` or ``proven``
     - ``grh`` (``belabas_friedman``)

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
internal constructor.  It accepts a root when ``f(r) = 0 mod q`` and ``f``
is squarefree modulo ``q``.  This is a sufficient condition for
``(q, theta - r)`` to be a prime of residue degree and ramification index one:
``q`` then does not divide the index of ``Z[theta]``, so the Dedekind-Kummer
theorem applies (Cohen, *A Course in Computational Algebraic Number Theory*,
GTM 138, Theorem 4.8.13).  It is not a necessary condition, and the
constructor refuses every root that fails it, even when ``(q, theta - r)`` is
still such a prime.  This is the same condition as the degree-one fast path of
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
elements.  This is an accepted deviation: the stopping rule decides only when
candidates are tested, never whether ``ell`` is verified.  The same constant
``3.5`` is used for the initial ``ell = 2`` pass and for the small-index
saturation loop in ``src/order_unit/compute.cpp``.  Silex first scans the
primes ``q <= aux_bound`` (the caller's bound; 1000 in class/unit
validation) and continues with the Hecke scan when that pre-scan selects no usable prime,
when it leaves a nonempty kernel without adjoining a root, or when
``ell >= aux_bound``.  The auxiliary bound is therefore not a proof
parameter; a unit proof record stores it as the pre-scan bound of the call,
not as the largest prime used.  Hecke starts its scan at
``next_prime(2^60)``; Silex starts at ``q = ell + 1``, which changes only
which valid characters are used.  Roots are adjoined only after the exact
order-unit check described below.

When ``ell`` divides the torsion order, the candidates include the torsion
generator ``zeta``, as in Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp/Saturate.jl``: ``compute_candidates_for_saturate``
(lines 130--137) appends ``zeta`` to the input units, and ``saturate!``
(lines 380--395) tests ``prod u_j^(e_j) * zeta^(e_t)`` for every candidate
kernel vector.  Silex tests the same products from the torsion-extended
character kernel, so a unit that is an ``ell``-th power only up to torsion,
such as ``-(3 + 2 sqrt2) = -(1 + sqrt2)^2`` at ``ell = 2``, is enlarged
instead of failing closed.  An exact root ``r`` with
``r^ell = prod u_j^(e_j) * zeta^(e_t)`` satisfies ``r^ell = prod u_j^(e_j)``
modulo torsion, so the new free basis uses only the free exponents.
Silex's ``zeta`` is ``group.torsion_generator``, the generator of the
*order's* unit group's torsion subgroup ``mu(O)``, not of ``mu(K)`` as in
Hecke's ``torsion_units_gen_order(K)``; ``mu(O)`` is the group that matters
for ``[O^x : G]`` when ``O`` is non-maximal, and this is a pre-existing Silex
choice, not one introduced with the torsion-extended candidates.  The
pre-scan tests the free-exponent kernel first and then, only if that leaves
the group unchanged, the kernel vectors with ``e_t != 0 mod ell``; this
ordering has no Hecke counterpart (``saturate!`` tests every candidate from
one combined kernel) and only affects which valid root a pre-scan pass finds
first, never whether ``ell`` is verified, since verification still requires
an empty torsion-extended kernel.

The non-proof saturation routines use the same torsion-extended candidates.
``OrderUnitGroup::residue_dlog_kernel`` returns the kernel of the
``ell``-th-power residue characters on the free generators and, when ``ell``
divides the torsion order ``w``, on ``zeta``: its rows have
``free_rank() + 1`` columns, the last one the exponent of ``zeta``, when
``ell | w`` and ``free_rank()`` columns otherwise, one row per candidate
exponent vector of ``compute_candidates_for_saturate`` (which returns the
candidates as the columns of its matrix).  ``saturate_row`` takes a row of
that width and tests ``prod u_j^(e_j) * zeta^(e_t)``, as in ``saturate!``;
``saturate_local_once``, ``saturate_bounded``, and the index-bounded passes
built on them test every row of that kernel.  As in the proof routes, a
root is adjoined only after the exact order-unit check, and these routines
never record a verified proof status.

A ``proven`` unit computation from a class-group context runs its
relation-kernel units through the index-bounded saturation pass
(``set_relation_kernel_units_index_bounded_saturated``) before
``prove_index_bound``, so that pass now also considers torsion-twisted
candidates when ``ell | w``.  It can therefore hand the proof a larger
subgroup, or a different generator of the same subgroup, than before.  The
proof routine and its rule are unchanged: every adjoined root is an exact
unit of the order, and ``ell`` is verified only when the torsion-extended
kernel is empty.

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

All of these bounds are non-decreasing in ``w``.  For an order ``O``, ``O^x``
has finite index in ``O_K^x``, so ``Reg(O) >= R_K``, and the field's ``w`` is
the right value to use.  If ``w`` cannot be computed, Silex uses ``w = 2``.  The
terms are evaluated at a working precision of at most 64 bits; the result is
an exact lower endpoint, so this only weakens the bound.  The floor in
``regulator_index_bound`` is taken of a rigorous upper endpoint, whereas Hecke
floors a floating-point quotient, so Silex follows Hecke's rule but,
unlike a floating-point floor, can never report a bound below the true
index.  The Zimmert Korollar and the Friedman Corollary are
dominated by the Satz 3 values at ``gamma = 1`` and ``gamma = 3/5``, which are
already evaluated; they are kept as documented cross-checks.

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

The honesty witness search for each required prime starts at fixed caps:
lattice radius 8, 16 twists, 50 random sub-factor-base tries (PARI 2.17.3
``src/basemath/buch2.c:be_honest``, ``maxtry_HONEST``), and, on the T2 route,
500 factored candidates and a 1e6-node enumeration step limit per ideal.
PARI's ``be_honest`` gives up at those caps and ``bnfinit`` enlarges the
factor base.  Silex keeps that flow wherever a factor-base restart or the
imaginary quadratic fallback route is available.  Where an unwitnessed prime
is a final failure, Silex instead doubles every one of these caps at each of
at most three further stages (the step limit doubles with the candidate cap
rather than growing quadratically, and the T2 enumeration bound targets four
times the candidate cap as before), and only then reports the prime
unwitnessed.  The doubling, the caps it applies to, and the three stages are
Silex's own bounded choices with no upstream source.  The nearest reference
schedule, ``SPLIT`` (used by ``bnftestprimes`` and ``isprincipal`` once the
class group is known), keeps a fixed set of twisting directions and instead
widens the number of ``Vbase`` primes in each random product, doubling its
try limit as that number grows and then removing the limit; the Silex random
tries keep the sub-factor base fixed.  The class-group detail log records the
caps of the last stage searched, and on the T2 route how many ideals stopped
at the candidate cap or the step limit.  When such a final failure ends the
last candidate attempt of a class/unit transaction, the transaction fails
with reason ``factor_base_honesty_unwitnessed`` instead of
``class_unit_computation_failed``, as ``silex-class-unit-instance`` prints in
``failure_reason``; the failure itself and every label are unchanged.  Every
candidate at every stage passes the same exact witness check, and stage 0 is
the previous search, so results that succeeded at stage 0 are unchanged.

In the ``silex-class-unit-instance`` JSON, the top-level ``failure_detail``
object is the one place for structured failure detail; its keys depend on
``failure_reason``, and it is null for every reason other than
``factor_base_honesty_unwitnessed``.  For that reason it gives the
unwitnessed prime as ``p``, a decimal string such as ``"3"`` like every big
integer in this JSON, and ``residue_degree``, the route as ``search``
(``lattice`` or ``t2``), the last stage searched and the last stage allowed
as ``stage`` and ``max_stage``, and that stage's caps under the detail log's
names: ``radius``, ``twists`` and ``random_tries`` on the lattice route, or
``random_tries``, ``factor_attempts`` and ``element_steps`` on the T2 route.
The residue degree, stages and caps are JSON numbers.

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
is not a hidden side effect of structural operations.  Expansion follows
Hecke's ``evaluate(FacElem)`` (``_ev`` in ``src/Misc/FactoredElem.jl``),
which halves exponents of absolute value at least 10 and squares the
recursive part; the small powers, odd residual factors, and that square at
each level are multiplied with a balanced product tree, as PARI's
``gen_product`` (``src/basemath/bb_group.c``) does for ``gen_factorback`` and
``nffactorback``, rather than into one running accumulator.  Field products
are exact and canonical, so the expanded value does not depend on the
multiplication order.  The reduced power and square tests expand only the
leftover factors whose exponents are not multiples of the root degree.

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
selected-prime class rows, HNF/kernel construction, valuation bases, and
nonunit generators.  The S-regulator is the S-unit lattice regulator
``R h^(S) prod log N(P)`` of Bernard and Roux-Langlois (ASIACRYPT 2020,
Prop. 2.2), with ``h^(S) = h / h_S`` taken exactly as the index of the
nonunit valuation HNF; PARI's ``bnfsunit`` reports ``R h_S prod log N(P)``
instead (``bnfunits.c:217-233``), which is ``h_S^2 / h`` times it.  See
:doc:`sunit_groups`.  Hecke 0.39.19
``src/NumFieldOrd/NfOrd/Clgp/Sunits.jl:1-265`` supplies the corresponding
class-group quotient, valuation-HNF, compact-generator, and combined-group
maps.  Silex deliberately orders coordinates as torsion, ordinary free, then
nonunit and exposes owned result objects rather than a reference-system
container layout.

The published witnesses are factored products whose exponents come from HNF
transforms and can be very large, so Silex does not verify them by expanding
the product.  For each class-relation HNF row, the witness is
``w_i = prod_j g_j^C[i,j]``, a product of stored relation generators.
For every generator ``g_j`` that some witness uses, Silex proves once
that its principal ideal equals the ideal of its stored factor-base row,
``(g_j) = prod_P P^R[j,P]``.  This is an exact fractional-ideal
equality, so it also rejects a generator that carries a prime outside the
factor base.  This check is required because some relation rows are installed
from known data without being factored again.  The ideal map
``v -> prod_P P^v[P]`` is a homomorphism, so the exact matrix
identities already checked give the claimed principal ideals:

- the HNF transform identity ``C * stored == rows``;
- the S-class identity ``coef * augmented == d * generator_rows``;
- the two relation-kernel products that must be zero;
- the valuation transform identity ``transform * raw == valuation_rows``.

The ideal relation of each selected prime is still checked by expansion when
it is found; its multiplier has only a few small factors.  Each S-class
invariant witness and each S-unit generator is also checked against its
factored valuations at every selected prime.  Where a selected prime equals a
factor-base prime, the expected valuation adds the contributions from both.
Correctness therefore rests on three things: the per-generator ideal proofs,
the checked matrix identities, and exact factored-element arithmetic.
Hecke's FacElem S-unit code places the same trust in factored arithmetic.
Silex's check is stricter than the upstream S-unit routine: Hecke 0.39.19
``Sunits.jl`` leaves both its ideal-relation check and its final generator
check commented out and trusts the relation rows.  A composition coefficient
that does not fit in ``slong`` is not refused: ``compact_multiply_power_fmpz``
evaluates the factored base and raises it to the ``fmpz`` power, which is exact
but may be expensive.  Factored exponent products that overflow ``slong`` in
``FactoredElement::pow_si``, and factored valuations that overflow ``slong``,
fail the build at the failing stage instead of wrapping.  A failed proof never
falls back to another check.

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
