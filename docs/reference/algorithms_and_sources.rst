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
strictly column-sorted and never store zero values.  When ``HnfContext``
refreshes a retained HNF of less than full rank without a transformation, it
folds the new rows into that HNF with the row-addition step of the
Pernet--Stein HNF algorithm as FLINT 3.6.0 implements it
(``src/fmpz_mat/hnf_pernet_stein.c``, ``add_rows``; C. Pernet and W. Stein,
"Fast computation of Hermite normal forms of random integer matrices",
*J. Number Theory* 130 (2010), 1675--1683), as the relation-module
checkpoint already does; other refreshes call ``fmpz_mat_hnf`` or
``fmpz_mat_hnf_transform``.  The canonical row HNF is unique, so the result
is the one ``fmpz_mat_hnf`` returns for all rows.

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
it takes the certified first LLL phase described next; larger inputs call
``fmpz_lll`` directly.  The existing zero-initialized transformation argument
is retained.

The certified first LLL phase is one internal helper shared by
``Lat::lll_reduce``, the unit-candidate kernel and the relation-saturation
kernel basis of the class-group computation, the rounded-embedding
integer LLL and the finite-ideal T2 enumeration setup (below).  For a ``Z_BASIS``, ``APPROX`` context, at least two rows
and coefficients of at most 250 bits it runs
``fmpz_lll_d_with_removal_knapsack`` with no removal bound and accepts the
result only when an exact reducedness certificate holds with the unchanged
parameters.  The reducer's status is not a reducedness certificate or a rank.
A failed reduction or certificate restores the original matrix and
transformation before calling ``fmpz_lll``; ineligible inputs call it
directly.  The certificate decides the predicate of ``fmpz_mat_is_reduced``
(every Gram--Schmidt norm nonzero, ``|mu_ij| <= eta`` and
``(delta - mu_{i,i-1}^2) |b*_{i-1}|^2 <= |b*_i|^2``, with ``delta`` and ``eta``
read as the exact rationals of the doubles).  It equals that predicate for
FLINT's parameter range ``eta^2 < delta``; all callers use
``eta = 0.51`` with ``delta = 0.99`` or ``delta = 0.5``.
When the reducer's output has
entries of at most 16 bits it is decided in integers with the integral
Gram--Schmidt recurrence of Cohen, *A Course in Computational Algebraic
Number Theory*, GTM 138, Algorithm 2.6.7, step 2; otherwise by
``fmpz_mat_is_reduced`` itself.  The cutoff only chooses the cheaper exact
test: on a busy machine the integer test was faster on every matrix the
class/unit rows produced (at most 14 bits after the reducer) and on dense
reduced matrices up to about 16 bits with 81 rows, and slower above about 20
bits as the row count grows.  The result is the matrix and transformation
``fmpz_lll`` returns, for three reasons: a passing certificate excludes zero
Gram--Schmidt vectors, so FLINT's own check (``fmpz_mat_is_reduced`` after
stripping initial zero rows) accepts the same output; a failed certificate
reruns ``fmpz_lll`` from the restored input; and the double phase is
deterministic.  Only "certificate true implies FLINT accepts" is used.  This
is Silex-specific routing, not FLINT's default checker order.  The source
anchors are FLINT ``src/fmpz_lll/lll.c`` (``fmpz_lll`` is
``fmpz_lll_with_removal_ulll`` with a 250-bit truncation),
``wrapper_with_removal_knapsack.c`` (first phase and mandatory
certification), ``lll_with_removal_ulll.c`` (the 250-bit truncation
boundary) and ``is_reduced.c``, inspected at revision
``b574c48f46ca505d6d48dd843c0fe95e1199ecdc`` (FLINT 3.6.0).  PARI 2.17.3
``src/basemath/lll.c:ZM_lll_norms`` likewise retains provable reduction after
its fast phases unless certification is explicitly disabled; Silex does not
disable it.  Tests force both rejection routes after scratch mutation at
every call site and require exact agreement with ``fmpz_lll``.

The same helper serves the two-LLL ideal-lattice route (grh class group,
``ideal_lattice_lll.cpp``).  Its two pre-reduction ``Z_BASIS`` sites use
parameters (0.5, 0.51) with ``eta^2 < delta`` and the certificate above.  Its
two ``GRAM`` sites are eligible when the context is ``APPROX`` and the matrix
is square with at least two rows, with no bit bound (FLINT's
``lll_with_removal_ulll.c`` makes one wrapper call for ``GRAM`` at every
size).  The GRAM certificate is the integer test above (Cohen, GTM 138,
Algorithm 2.6.7, step 2, applied to the Gram matrix) while rows times maximal
entry bits is at most 3000, and FLINT's own ``fmpz_lll_is_reduced(B, fl,
120)`` above that, the very call FLINT's wrapper makes on the same matrix.
``fmpz_mat_is_reduced_gram`` is not used: it divides by zero on a zero
Gram--Schmidt norm, whereas a passing integer test excludes zero norms, so
the exact quantities agree and FLINT's check accepts too.  The cutoff is
measured: on corpus site matrices the integer test took at most 0.79 of
FLINT's check up to about 2500, was even near 3400 (degree 32) and was far
slower at degree 81.  The transformation argument does not steer the
``GRAM`` d-pass: at the pinned revision ``lll_d.c`` reads it only for its row
count (452-456) and in the ``update_b`` test (458), which needs more columns
than rows and never holds for a square Gram matrix; all other uses write row
operations into it.  The original matrix and transformation are restored
before any ``fmpz_lll`` fallback.

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

``Element::norm`` uses FLINT ``nf_elem_norm`` except in two field shapes.
Quadratic fields use ``a^2 - d b^2``.  When the defining polynomial is
``f = x^n + c`` with integer ``c`` and even ``n >= 4``, the norm is computed
by halving the degree.  FLINT 3.6 ``src/nf_elem/norm.c`` (``_nf_elem_norm``,
generic branch) computes the norm of ``Z/d`` as a resultant: for monic
integral ``f`` it is ``Res(f, Z) / d^n`` (Cohen, *A Course in Computational
Algebraic Number Theory*, GTM 138, Proposition 4.3.4: ``N(a) = d^(-n)
R(T, A)`` for ``a = A(theta)/d``), the product of ``Z(alpha)`` over
the roots ``alpha`` of ``f``.  Write ``f = h(x^2)`` with ``h = y^(n/2) + c``
and ``Z(x) = E(x^2) + x O(x^2)``.  The roots of ``f`` are the two square
roots of each root ``beta`` of ``h``, so ``Res(f, Z) = Res(h, G)``, where
``G(y) = E(y)^2 - y O(y)^2`` is the Graeffe transform of ``Z``,
``G(x^2) = Z(x) Z(-x)`` (PARI 2.17.3 ``src/basemath/rootpol.c``,
``ZX_graeffe``; GP ``polgraeffe``).  ``G`` is first reduced modulo ``h``,
and the step repeats while the degree is even; for ``x^12 + 7`` it goes
from degree 12 to 6 to 3, and the last step is FLINT's integer resultant.
The identity ``Res(f, Z) = Res(h, G)`` is the proof; it does not need ``h``
to be irreducible.  (When ``f`` is irreducible, ``h`` is too, because
``[K : Q(alpha^2)] <= 2``, and the step is transitivity of the norm through
``Q(alpha^2)``.)  The arithmetic is
exact and the result is the same canonical rational that ``nf_elem_norm``
returns, which a randomized differential test checks.

``src/order`` implements validated order bases, trace forms, discriminants,
multiplication tables, indices, p-maximal overorders, and global maximal-order
construction.  PARI ``src/basemath/base2.c`` (including the ``nfmaxord``
family) and Hecke
``src/NumFieldOrd/NfOrd/MaxOrd/MaxOrd.jl`` and
``MaxOrd/DedekindCriterion.jl`` are the main external anchors.  The standard
quadratic identity ``disc(O_f) = disc(O_K) f^2`` governs direct conductor
orders.

Before the Round 2 multiplier-ring loop, ``maximal_order`` and
``pmaximal_overorder`` apply Dedekind's criterion to an equation-order input
``Z[theta]`` whose defining polynomial ``f`` is monic integral (Cohen, GTM 138,
Theorem 6.1.4 (2) and Algorithm 6.1.8 steps 4-5, the ``m = 0`` case; PARI 2.17.3
``src/basemath/base2.c:ZX_Dedekind`` and ``maxord``, where a gcd of degree 0
means ``p``-maximal).  When ``Z[theta]`` is ``p``-maximal the input order is
kept as is, which is the order the multiplier ring would return (Cohen,
Theorem 6.1.3).  In ``maximal_order`` the test is also valid after
enlargements at other primes: the running order contains ``Z[theta]``, and
the criterion holds for any order containing ``Z[theta]`` with index prime to
``p`` (Cohen, GTM 138, p. 306).  The criterion only replaces a computation that would
return index 1; nonmaximal primes still go through the multiplier-ring loop.

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

A principal ideal ``g O`` is the row lattice of the matrix ``M`` of
multiplication by ``g`` in the order basis, and
``|det M| = |N(g)| = [O : g O]``.  When the relation
search already knows N(g) exactly (the fallback that factors an integral
relation generator as an ideal), it computes the HNF of M with FLINT 3.6.0
``fmpz_mat_hnf_modular`` (``src/fmpz_mat/hnf_modular.c``), the HNF modulo D
of Domich, Kannan and Trotter (1987; H. Cohen, *A Course in Computational
Algebraic Number Theory*, Algorithm 2.4.8), with ``D = |N(g)|``.  Its
precondition, D a positive multiple of the lattice determinant, holds since
``adj(M) M = det(M) I``.  FLINT documents the output as the unique row HNF, the
one ``fmpz_mat_hnf`` returns, so the stored HNF is unchanged.  A diagonal
product other than D, or a norm that is not a nonzero integer, takes the
``fmpz_mat_hnf`` route.  The norm is trusted, as the residue screen already
trusts it for valuations: the diagonal check catches only some wrong norms (one
that properly divides ``det M`` can pass) and is not a proof of the norm.  Other principal-ideal constructions use
``fmpz_mat_hnf``.

The relation search inverts each search-base prime ideal P of a maximal order
with ``detail::prime_ideal_inverse_num_den`` instead of the general inverse.
It follows PARI 2.17.3 ``src/basemath/base4.c:pr_inv_p`` and ``pr_inv``
(p P^(-1) = p Z_K + t Z_K) with ``src/basemath/base2.c:get_pr`` (t a kernel
vector modulo p, "Z_K + t/p Z_K = P^(-1)"), and Hecke
``src/NumFieldOrd/NfOrd/Ideal/Prime.jl:anti_uniformizer``.  Silex takes the
common kernel modulo p of y -> y b over the HNF rows b of P, stopping at
dimension log_p N(P), and accepts p Z_K + t Z_K only if its index is
p^(n - f); otherwise the general inverse runs.  Silex's implementation is its
own; tests compare the numerator HNF and denominator with the general inverse.

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
candidate ``c`` satisfies ``c^n == a`` exactly, and a definite ``false`` only
after a residue disproof (below); every unported or unverified case returns
failure (unsupported) and leaves the caller's root unchanged.  In degree 10
and higher there is no lift, so only the residue disproof can answer, and
only with ``false``.

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
* **Silex extension** (not in Hecke): pre-lift residue disproof.  Hecke tests
  only the first factor of the first good prime before lifting (in degree
  below 10), and when no lifted candidate verifies it answers ``false``;
  Silex instead fails closed there.  So that such inputs still get an exact
  answer without the lift, Silex first applies the residue test to every
  factor of the chosen prime and of up to 40 further admissible primes, and
  answers ``false`` at the first residue that is not an ``n``-th power; if
  there is none, the lift runs exactly as before.  A rational prime ``p`` is
  admissible when the monic integral defining polynomial ``f`` is
  squarefree modulo ``p`` and the power-basis coordinates of ``a`` are
  ``p``-integral; only factors ``g`` of ``f mod p`` with ``a mod g`` nonzero
  are tested.  Argument: ``p`` does not divide ``disc(f)``, so ``p`` is
  unramified and prime to ``[O_K : Z[theta]]``, and by Dedekind-Kummer
  (Cohen, *A Course in Computational Algebraic Number Theory*, Theorem
  4.8.13) ``P = (p, g(theta))`` is a prime of ``O_K`` with
  ``O_K / P = F_p[x] / (g)``.  ``a`` lies in ``Z_(p)[theta]`` and is a
  ``P``-unit, so if ``a = b^n`` with ``b`` in ``K`` then ``b`` is a
  ``P``-unit and ``a mod P = (b mod P)^n``.  A residue that fails the test
  above therefore proves that ``a`` is not an ``n``-th power; ``p | n`` is
  allowed.  The test runs on each block of the distinct-degree factorisation
  of ``f mod p`` (with the factors where ``a`` vanishes divided out), which by
  the Chinese remainder theorem is equivalent to testing each factor of that
  degree.  The bound 40 is a resource limit only, and any value is sound: for
  ``n = 2`` a rational number is a square in every ``F_(p^f)`` with ``f``
  even, so ``c b^2`` with ``c`` in ``Q`` can be disproved only at primes of
  odd residue degree, and the class-group saturation rows
  ``multiquadratic_m3_m5_m7`` and ``lmfdb_8_8_73116160000_2`` need up to 21
  and 36 further primes.
* **Silex extension**: residue disproof in degree 10 and higher.  When the
  monic integral defining polynomial has degree at least 10 and no earlier
  path decided, Silex applies the same residue test (and the same
  admissibility argument) to every factor of each of the first 40
  admissible primes from 5 upward, and answers ``false`` at the first
  residue that is not an ``n``-th power.  This corresponds to Hecke v0.38.6
  ``_roots_hensel`` (``Hensel.jl`` lines 159-161), which returns no roots as soon as
  ``y^n - a`` has no root modulo a factor of a good prime, and which in
  degree 10 and higher examines further factors before lifting
  (``STABILIZED = 2``, lines 102-106).  Hecke stops after a few factors and
  then lifts; Silex has no lift here and scans more primes.  A true power is
  never disproved, so it and every input with no disproving residue within
  the prime bound still return failure (unsupported).  The bound 40 (the
  same as the pre-lift disproof) is a resource limit only, and any value is
  sound.  Every true power runs to the bound.  A
  residue can fail to be an ``n``-th power only in a residue field of order
  ``q`` with ``gcd(n, q - 1) > 1`` (``F_q^*`` is cyclic of order ``q - 1``;
  Lidl and Niederreiter, *Finite Fields*, Theorem 2.8), and every residue
  degree is at most ``N``.  So a prime ``p`` with ``gcd(n, p^f - 1) = 1``
  for every ``f <= N`` cannot disprove the input; Silex skips it exactly,
  without factoring, and does not count it towards the bound.  The scan also stops after 256
  attempts, skipped primes included, and for a large ``n`` that attempt cap
  is the bound that ends it.  Without that
  skip, the 9 non-powers each of the compact-root kernel rows of the
  ``proven`` runs on ``Q(zeta_25)`` and ``Q(zeta_44)`` needed 9 to 218
  admissible primes; with it, all 18 are disproved within the bound.
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

When ``p`` divides the index of the maximal order, ``decompose_prime`` splits
the etale algebra ``O/(p, I_p)``, where ``I_p`` is the ``p``-radical (equal
to ``pO`` when ``p`` is unramified), into its fields through Frobenius-fixed
idempotents.  It first tries two bounded enumerations of the fixed
subalgebra, a separating element as in Hecke
``AlgAss/Decompose.jl:_dec_com_finite`` but taken from deterministic
combinations, then the fixed idempotents themselves, and both stop when
``p`` raised to the number of primes above ``p`` exceeds 4096.  Beyond that
bound it next tries the Dedekind--Kummer route through a monogenic generator
of the maximal order, and only if that search finds none does it split
recursively as PARI 2.17.3 ``src/basemath/base2.c:primedec_aux`` does
(lines 2316--2350, after Cohen, GTM 138, Algorithm 6.2.9): a piece that is not a field has a non-scalar fixed
element whose minimal polynomial splits into distinct linear factors, and
each root cuts off a smaller piece, with no bound on ``p`` and no
randomness.  Silex keeps a piece as an idempotent ``e`` rather than PARI's
subspace, takes the first ``e k_j`` (``k_j`` the fixed-kernel basis) that is
not a multiple of ``e``, and returns the fields in the order the idempotent
enumeration lists them (by the base-``p`` combination index of their
fixed-kernel coordinates).  The primitive idempotents are unique, so neither
the prime ideals nor their order depend on the element chosen.  Because the
recursive split runs last, decompositions the bounded enumerations or a
monogenic generator already produced are unchanged, prime order included.
Unlike PARI, Silex does not split off Dedekind--Kummer factors before the
radical when ``p`` divides the index.

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

The one-large-prime step of relation factoring takes valuations of an element
``alpha`` of the maximal order at the primes above each rational prime ``p``
of its norm that the factor base does not account for.  When
``decompose_prime`` would build those primes by Dedekind--Kummer (``p`` prime
to ``[O : Z[theta]]``, off the quadratic integral-generator route), it builds
only the ``P_i = (p, g_i(theta))`` with ``g_i`` dividing the reduction of
``A`` modulo ``p``, where ``alpha = A(theta)/d`` with ``d`` prime to ``p``.
``P_i`` meets ``Z[theta]`` in ``p Z[theta] + g_i(theta) Z[theta]`` (Cohen,
GTM 138, Theorem 4.8.13 and its proof), so ``alpha`` lies in ``P_i`` exactly
when ``g_i`` divides ``A`` modulo ``p``; the primes left out have valuation 0
and were skipped before.  The primes that are built, and their order, are
those of ``decompose_prime``.

A relation generator ``alpha`` in a maximal order is first factored by a
residue screen.  For each rational prime ``p`` dividing ``N(alpha)``, it
tests ``alpha`` against the residue polynomial of every factor-base prime
over ``p``.  If the residue degrees of the primes containing ``alpha`` add up
to ``v_p(N(alpha))``, each of those primes has valuation one.  Otherwise
their exact valuations are computed, and the row is published only when
those valuations account for all of ``v_p(N(alpha))``, so a published row is
the exact valuation vector of ``(alpha)``.  Any other case goes to the exact
factorization of the principal ideal.  On the equation order the integral
coordinates are the ``theta``-polynomial coefficients.  On a maximal order
that is not the equation order, the screen tests a prime over ``p`` not
dividing ``[O : Z[theta]]`` on the ``theta``-polynomial numerator of the
element, by the ``modprinit`` branch above (Cohen, GTM 138, Theorem 4.8.13).
That numerator's denominator divides the order basis denominator, so it is
prime to ``p``.  The index is the order's cached ``[O : Z[theta]]``, which
exists when ``Z[theta]`` is contained in ``O``.  An element whose norm is
divisible by an index divisor goes to the principal-ideal factorization.  The
exception is a prime stored in the quadratic generator ``omega``, which is
tested on the element's ``[1, omega]`` coordinates, its ``omega``-polynomial,
for every ``p``, since the order is ``Z[omega]``.

Relation factoring over the factor base (the integral residue screen of
``Relation`` and the direct one-large-prime partial-relation row of the
class-group search) needs, for each rational prime ``p`` of ``N(alpha)``,
the valuations of ``alpha`` at the factor-base primes above ``p``.  On a
maximal equation order the residue test above finds the primes above ``p``
that contain ``alpha`` (the hits), and exact valuations are computed only at
hits.  The last hit's valuation is not computed when the factor-base block
of ``p`` holds every prime above ``p``, that is, when the sum of ``e f`` over
the block equals ``[K : Q]`` (Cohen, GTM 138, Section 4.8; the factor base
keeps distinct primes, one block per rational prime).  Then
``|N(alpha)| = N(alpha O_K)`` and ``N(P) = p^f`` give
``v_p(N(alpha)) = sum over P | p of f(P) v_P(alpha)``, every prime above
``p`` with positive valuation is a hit, and the last hit's valuation is what
is left of ``v_p(N(alpha))`` after the other hits, divided by its residue
degree.  PARI 2.17.3 ``src/basemath/buch2.c:divide_p_elt`` does the same
norm accounting, subtracting ``v f`` per prime, but computes every
valuation.  The deduction is used only when what is left is a positive
multiple of the residue degree; otherwise, and for blocks that miss a prime
above ``p`` (a norm-bounded base, for example), the valuation is computed,
so the exponent rows are unchanged.

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

The ``proven`` analytic finish needs a small basis of the integer kernel
of the relation matrix ``A`` (the unit exponent rows).  Silex reduces the
rows of ``[c*A | I]`` with ``fmpz_lll`` for a power-of-two scale ``c`` of
``2^40``, and takes the identity parts of the rows whose ``A`` part is zero
as the kernel basis.  This follows the orthogonal-lattice construction of
J. Chen, D. Stehle, and G. Villard, "Computing an LLL-reduced basis of the
orthogonal lattice" (arXiv:1805.03418), Proposition 4.1, adapted from
Nguyen and Stern (1997), whose sufficient scale is
``K > 2^((n-1)/2) * lambda_{n-k}``, where ``k`` is the rank of the matrix.
The proposition assumes a matrix of full column rank, which the relation
matrix generally is not, so Silex does not rely on it.  Silex does not rely
on the scale being large enough: the number of candidate rows must equal ``n`` minus
the rank, and every row is checked exactly against ``v*A = 0``.  When
either check fails, the helper uses the rows of the HNF transform that
belong to the zero rows of the HNF, which are an exact kernel basis; the
HNF transform is computed only in that case.  As
in the reducer paragraph below, the basis claim relies on FLINT's guarantee
that ``fmpz_lll`` applies only unimodular row operations; Silex adds no
determinant check, and keeps FLINT's final reducedness check.  The rank
passed in must be the exact rank of ``A``; the full-lattice claim depends on
it.  The finish checks that the relation rank it passes equals
``fmpz_mat_rank`` of ``A``, the number of nonzero rows of its HNF.  The
two routes return different bases of the same kernel lattice.  The class
group, the regulator, and the unit lattice do not depend on the
route; later LLL reductions start from the basis they receive, so the
specific unit generators are not shown to be route-independent.  On the
tested proven corpus rows the published outputs were identical.

Continued-fraction fundamental unit of a real quadratic field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``quadratic_fundamental_unit`` in ``src/unit/unit.cpp`` computes the
fundamental unit of the maximal order of a real quadratic field by the
classical continued-fraction period theorem.  It accepts the explicit
quadratic backend and any field with a monic integral defining polynomial
``x^2 + b*x + c`` of positive nonsquare discriminant.  The basis is
unconditional: no GRH or heuristic hypothesis enters.

Source.  The core follows PARI 2.17.3 ``src/basemath/quad.c``.  It is the
sequential convergent recurrence of ``quadunit_uv_basecase`` (lines 280-316):
the ``sqrtremi`` start with its parity correction, the even-period and
odd-period exits, ``u = (a + D*b)/q`` and ``v = (c - a - b)/q``, the
``u -= v`` adjustment for odd discriminant, and the final halving.
``quadunit_uv`` (lines 428-474) is the same algorithm with a product tree of
the partial-quotient matrices; PARI's ``quadunit`` (lines 475-518) calls it
for discriminants of at least ``2000000`` and the basecase below that.  The
nonmaximal-conductor branch of ``quadunit`` (lines 492-516), which raises the
maximal-order unit to the index from ``quadunitindex``, is not ported: Silex
rejects nonmaximal orders instead (see Runtime checks).

Presentation map.  With ``D = b^2 - 4*c = f^2*d`` and ``d`` squarefree, a root
``theta`` satisfies ``(2*theta + b)^2 = D``, so
``sqrt(d) = (2*theta + b)/f`` is the square root of ``d`` in the field that is
positive at the real embedding sending ``theta`` to the larger root
(``monic_integral_quadratic_radicand_map``, ``src/unit/unit.cpp:72-146``,
which factors ``D`` to split off ``f``).  For the explicit quadratic backend
the generator is ``sqrt(d)`` itself and the map is ``b = 0``, ``f = 2``
(``fundamental_unit_radicand_map``, ``src/unit/unit.cpp:148-163``); it relies
on the backend radicand being squarefree, which
``is_valid_quadratic_radicand`` (``src/number_field/number_field.cpp:66-69``)
enforces whenever that backend is set.  The squarefree ``d`` matters: on a
non-squarefree radicand the core would in general return a unit that is not
fundamental for the maximal order.  The continued-fraction core runs on
``d``, so the unit does not depend on the presentation.  Its coefficients are rewritten through this map: the ``theta``
coefficient is ``2*B/f`` and the constant is ``A + B*b/f``, where the unit is
``A + B*sqrt(d)`` (``set_quadratic_coefficients``,
``src/unit/unit.cpp:37-70``).  The ``x^2 - d`` presentation with ``d``
squarefree is the case ``b = 0``, ``f = 2``.  Here ``f`` is not the conductor
of ``Z[theta]``, which is ``f/2`` when ``d = 2, 3 (mod 4)``.  For a
nonmaximal ``Z[theta]`` the unit is generally not in ``Z[theta]``; only the
maximal order is claimed.

Runtime checks.  Together these establish that the result is a non-torsion
unit of the maximal order; none of them certifies that it is fundamental.

* The order is maximal: ``OrderUnitGroup::compute`` returns failure unless
  ``order.is_maximal()`` (``src/order_unit/compute.cpp:3496``); the equation
  order of a nonmaximal ``Z[theta]`` is rejected.
* The divisions in the coefficient reconstruction are checked to be exact
  (``src/unit/unit.cpp:431-459``).
* The unit has norm ``+1`` or ``-1``, tested exactly in integers as a norm
  numerator of ``+-4`` for odd discriminant and ``+-1`` otherwise
  (``src/unit/unit.cpp:482-487``).
* The unit is then published through ``set_units``
  (``src/order_unit/compute.cpp:3509``), whose validating path
  (``src/order_unit/order_unit.cpp:1174-1319``) computes the torsion exactly
  (``rank_zero_torsion``, line 1230), checks that the generator is a unit of
  the order (``evaluated_is_order_unit``, line 530), and checks log
  independence, so the generator is not torsion (``compact_independent``,
  line 1286).

That the unit is fundamental rests on the period theorem and the reference
``quadunit``, the same basis the ``x^2 - d`` route already had.

Certification label.  On success of ``quadratic_fundamental_unit`` the caller
sets ``CertificationMode::proven``
(``src/order_unit/compute.cpp:3507-3513``) with no unit proof record
(``unit_proof_record_count`` is 0).  This does not contradict the sentence
above that Silex never upgrades certification without the corresponding proof
record: nothing is upgraded, because the label is ``proven`` from the start by
the unconditional theorem rather than by a record, and no candidate or ``grh``
result is promoted.

The class group depends on this route.  For degree-2, rank-1 orders the
paired class and unit computations take this unit directly
(``src/order_unit/compute.cpp:3699``, ``:4161``, ``:4410``).  Validation then
treats a ``proven`` quadratic unit as the full unit group and attributes any
remaining ``h*R`` index to missing class relations
(``independently_proven_quadratic_unit``,
``src/order_unit/validation.cpp:1478-1491``).  A class number certified there
from ``h*R`` uses this unit's ``proven`` label and its regulator, so a proven
class group of a real quadratic field also rests on the period theorem.

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

The small-LLL phase of the LLL relation search reduces each ideal with one
integer LLL, as PARI 2.17.3 ``buch2.c:Fincke_Pohst_ideal`` (2461-2462) does
on the rounded embedding matrix ``nf_get_roundG``, which
``base1.c`` (1851) stores as ``RM_round_maxrank(G)`` (``base4.c``
2631-2640): the embeddings scaled by ``2^e`` for ``e`` = 4, 8, 16, ...
and rounded, until the result has full rank.  Silex rounds the order's
Minkowski rows this way once per order and per relation route, and reduces
``B R`` for the ideal's HNF basis ``B``.  PARI rounds the embeddings of a
T2-reduced integral basis; Silex rounds those of its order basis, which is
not reduced.  The rounded matrix only steers the reduction: Silex checks
exactly that the transformed basis reproduces the LLL output, forms the
scaled T2 Gram of the reduced basis in integers from the order Gram of the
Hecke precision schedule above, and applies only the two diagonal bounds of
the Hecke route's checks, on ``|b_1|^2`` and on the product of the
``|b_i|^2``.  It skips the ``LLL.jl:_lll`` check on the size of the Gram
transform, as this route has no Gram transform.  When a check fails, or the
order data cannot be formed, the ideal is reduced by the Hecke route above.
Because Silex's order basis is not reduced, its rounded rows can be too
coarse to steer the LLL: on a badly conditioned basis the reduced basis
passes the two diagonal bounds yet is far longer than the Hecke route's.
Silex therefore adds a check that is not in either source: for the first
ideal whose one-LLL basis is accepted in an order (once per order in a
class-group context, not once per relation slice), it also reduces the ideal
by the Hecke route, and when the product of the one-LLL basis's T2 diagonal
entries exceeds ``2^n`` times the Hecke route's (``n`` the degree, so on
average a factor 2 per entry), that ideal and every later ideal of the order
are reduced by the Hecke route.  When the Hecke reduction of that ideal
fails, the one-LLL basis is used and the check waits for the next ideal;
when the comparison inputs are malformed, the order is switched to the
Hecke route.

The one integer LLL of this route takes the certified first LLL phase
described with ``Lat::lll_reduce`` (above): for at least two rows and
entries of at most 250 bits it runs ``fmpz_lll_d_with_removal_knapsack`` and
accepts the result only when the exact certificate holds, and otherwise
restores the input and calls ``fmpz_lll``.  The FLINT files and functions
named here are those of the revision pinned in that paragraph (b574c48).
The matrix and transform are those ``fmpz_lll`` returns, for the three
reasons given there.  Only the certificate is cheaper.
The two-LLL grh route uses the same helper at its ``Z_BASIS`` and ``GRAM``
integer LLL sites, with the ``GRAM`` certificate and size cutoff described
there.

The finite-ideal T2 enumeration reduces the ideal lattice from its embeddings
rounded at an integer scale, starting at 4 as PARI ``src/basemath/base4.c:
RM_round_maxrank`` does.  PARI rounds the embeddings of a T2-reduced integral
basis; Silex's order basis is not reduced, so the rounding error times the LLL
transform can swamp the reduced vectors (for the degree-one primes above 89 in
the real cyclotomic field of conductor 11 and degree 10, the last
Gram-Schmidt pivot was about ``-1e-8`` at every precision and the ten
precision retries repeated the same input).  As an intentional deviation,
when the double-precision Gram-Schmidt test fails, or the reduced basis
misses Hecke's ideal-LLL acceptance bounds (``LLL.jl:_lll``: the bound on
``|b_1|^2`` and the bound on the product of the ``|b_i|^2``, with ``d`` the
determinant of the Gram matrix), Silex reduces the same lattice again from
the next finer full-rank scale (4, 8, ... up to the working precision).
When a finer reduction fails, or no finer scale is left, the first reduced
basis that had quadratic-form data is kept.  The lattice and the bound rule
are unchanged; the bound value and the enumeration order follow the new
reduced basis.  Each of these integer LLLs (``Z_BASIS``, ``APPROX``,
``delta = 0.99``, ``eta = 0.51``) takes the certified first LLL phase
described with ``Lat::lll_reduce`` (above): for rounded rows with entries of
at most 250 bits it runs ``fmpz_lll_d_with_removal_knapsack`` and accepts the
result only when the exact certificate holds, and otherwise restores the
rows and transform and calls ``fmpz_lll``; rows above 250 bits call
``fmpz_lll`` directly.  The reduced rows and transform are those
``fmpz_lll`` returns, for the three reasons given there, so the reduced
basis, its quadratic-form data and the bound are unchanged.

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
check.  It extracts an ``m``-th root only when that check cannot complete.
The bounded relation-kernel refinement
(``OrderUnitGroup::set_relation_kernel_units_bounded``) uses a different,
Silex-specific rule that has no cited source.  For ``m >= 2`` it first looks
for any exact ``m``-th root of ``prod u_i^{a_i} zeta^i``, where ``zeta``
generates the torsion.  It adjoins the first exact root found without
comparing it with ``y``.  For that root the relation certifies that its
``m``-th power lies in the current group modulo torsion; it certifies nothing
about ``y``.  When the candidate is a true relation for ``y``, any such root
equals ``y`` up to torsion.  When the candidate is a wrong guess, the refined
group may not contain ``y``.  When that search does not find a root or cannot
complete, Silex checks ``y`` with ``_check_relation_mod_torsion`` and, if the
relation holds, adjoins ``y`` itself, as Hecke does.  The adjoined element is
always the one whose relation was verified.  Hecke's ``_add_dependent_unit!`` always decides
with ``_check_relation_mod_torsion`` on ``y`` (``UnitGrpCtx.jl`` line 97 in
v0.38.6) and adjoins ``y`` itself through ``_transform(vcat(U.units, y), m)``
(line 131).  The root is an exact field element, and it is adjoined only
after an exact check that it lies in the order.  The refined group is
therefore still a subgroup of the order's units.  An inconclusive attempt or
a failed check doubles the precision.  Hecke's loop is unbounded; Silex stops at a fixed precision cap
and then leaves the group unchanged.  The loop terminates below the cap in
the usual case: when the true coordinate is ``a/b`` with ``b`` at most the
bound ``B`` and the ball radius is below ``1/(2B^2)``, the ball has diameter
below ``1/B^2``, while every other fraction ``p/q`` with ``q <= b`` lies at
distance at least ``1/(qb) >= 1/B^2`` from ``a/b``, so the simplest rational
in the ball is ``a/b``.

Hecke's ``B`` is the a priori bound ``R_G / R_low``, with ``R_G`` the
regulator of ``G`` and ``R_low`` the unconditional lower regulator bound,
and it may have any size.  Silex's search takes a word-size ``B``, so when
the a priori bound does not fit in ``slong`` the attempt leaves the group
unchanged, except in a ``grh`` class/unit request.  There, while the class
group relations improve the units, the bound is replaced by the class/unit
index bound, the outward-rounded ceiling of the upper endpoint of the ball
quotient ``h_cand R_G / hR``, with ``hR`` the GRH-conditional
Belabas-Friedman value.  The fallback always uses the class/unit validation
estimate, whereas validation may use the Belabas-Friedman validation
estimate instead.  The replacement is made provided the factor base generates
the class group unconditionally or under GRH and the replacement fits in
``slong``.  Under GRH ``h`` divides ``h_cand``, so ``[U:G] = R_G / R`` is at
most ``h_cand R_G / hR`` and the closed range above still holds.  The bound
is recomputed for the current group after each adjoined unit.  It is a
Silex choice, not Hecke's; candidate relations are still verified exactly,
so a bound that is too small only misses relations.  ``proven`` requests and
the other callers keep the a priori bound.

Extracting independent units for a ``grh`` class/unit computation decides, for
each candidate unit, whether its logarithmic row is independent of the rows
already kept.  Hecke ``src/NumFieldOrd/NfOrd/Unit/Independent.jl:_isindependent``
tests the Gram determinant of the kept rows plus the candidate against a lower
bound.  Silex keeps the rows, the inverse of their Gram matrix ``G`` and ``det G``
(rebuilding ``G`` from the rows when the kept count changes) and evaluates
that determinant incrementally.  With ``L`` the kept rows, ``v`` the
candidate row, ``g = L v`` and ``c = v . v``, the Schur complement identity
(Horn and Johnson, *Matrix Analysis*, 0.8.5) gives
``det [[G, g], [g^T, c]] = det(G) (c - g^T G^{-1} g)``.  Every factor is a ball
enclosure, so the result contains the true determinant, and the decision rule
(dependent below the lower bound or exactly zero, independent when positive) is
the one applied to the full determinant.  If the matrix inverse cannot be
proved, or the ball is undecided, the precision is doubled up to the same cap
as before, and the full Gram determinant is used when the test stays
undecided.  Whenever the full test decides, the incremental test gives the
same decision.  Two cases differ, both sound.  A candidate that the full test
dropped at the precision cap may now be decided, since the incremental ball
can decide at the same level, and no cap event is then emitted.  And when the
logarithmic embedding failed at some precision ``p``, the full test dropped the
candidate without doubling, whereas the incremental loop retries at ``2p``,
where the embedding's own cap is larger, and may decide.  Both cases need
about ``2^20`` bits of precision.

The unit-basis transform that follows a reduction applies the integer matrix
to the stored factored units.  Hecke keeps factored elements as ``FacElem``
dictionaries with ``ZZRingElem`` exponents (``FactoredElem.jl:_transform``,
339-365; ``Reduction.jl``, 88-112), so the transform never overflows.
Silex's public ``FactoredElement`` has word-size exponents.  The transform
first runs in word-size form; when an exponent overflows, Silex repeats it
with integer exponents of any size, merging equal bases as a ``FacElem``
dictionary does.  When the merged exponents fit in ``slong``, those products
are the new units.  Otherwise the group keeps the wide form in private
storage, and ``free_generator`` returns the equal unit rebuilt in word-size
form by the bounded compact reconstruction, which ports PARI 2.17.3
``buch2.c`` ``vec_chinese_units`` (1214-1240), ``chinese_unit`` and
``FlxqX_chinese_unit`` (1026-1100), the bounds ``log2fubound`` and
``log2Mbound`` (1180-1212), and ``ZM_remove_unused`` (1104-1118).  Returning
the rebuilt unit instead of the wide factored form is a deviation from
Hecke; the unit is equal, so no published value changes.  Equality rests on
the coordinate bound, not on a separate check afterwards: PARI evaluates
``log2Mbound`` and ``log2fubound`` in doubles, while Silex evaluates them
with Arb balls and raises the precision until the ceiling of the bound is
stable, and the CRT runs until the modulus exceeds twice that bound, so the
centred coordinates are exact.

PARI's comment on ``log2Mbound`` (``buch2.c``:1201) states the quantity it
bounds: :math:`\log_2(|x|_\infty / |y|_\infty)` for the solution ``x`` of
``RgM_solve_realimag(M, y)``, that is :math:`\log_2 \|S^{-1}\|_\infty` (the
largest absolute row sum) for the matrix ``S`` of real parts of the
:math:`r_1 + r_2` embeddings and imaginary parts of the :math:`r_2` complex
ones.  ``log2Mbound`` itself is a Hadamard and Cramer estimate of that norm
on PARI's T2-LLL-reduced integral basis.  Silex evaluates the formula on its
own order basis (with 1 moved first), often a power basis, where the
estimate also carries the basis's orthogonality defect.  As a Silex
evaluation choice, not PARI's literal expression, Silex evaluates
:math:`\|S^{-1}\|_\infty` itself from a certified Arb inverse of the same
weighted embedding matrix (the upper end of the ball is an outward-rounded
bound) and publishes that bound whenever it is certified: the inverse
succeeds at the working precision and the ceiling of the total is stable.
Only when it is not certified does Silex publish the bound from PARI's
expression.  The PARI expression alone still decides success, the precision
ladder and failure.  The norm never exceeds
``log2Mbound`` on any signature, so the published bit bound never exceeds
the bound from PARI's expression (Cramer's rule, Cauchy-Schwarz and
Hadamard's inequality, with the basis element 1 of smallest T2 norm).  For
:math:`r_2 > 0` the proof inverts the weighted matrix, whose determinant is
:math:`\sqrt{|D|}` and whose weights sum to :math:`n + 2 r_2`; this pays for
PARI's :math:`r_2 \ln 2`, with margin at least
:math:`r_2 \ln 2 - \tfrac12 \log_2 (1 + 2 r_2 / n)`, sharp: attained for
some bases of imaginary quadratic fields (for example
:math:`\mathbb{Z}[i]` with basis :math:`1, 1+i`).  The CRT therefore never uses more primes than before.
For the totally real degree-32 subfield of the 136th
cyclotomic field the Minkowski term falls from 442.6 to 18.9 bits, so the
wide-unit bound falls from 454 to 30 bits and the CRT uses 1 prime instead
of 8.

Modulo each CRT prime ``p``, the reconstruction evaluates the power products
in ``R = (Z/pZ)[x]/(f)`` for the reduced defining polynomial ``f``, with
interleaved sliding-window exponentiation (Moller, "Algorithms for
multi-exponentiation", SAC 2001, Section 3.1; Menezes, van Oorschot and
Vanstone, *Handbook of Applied Cryptography*, Algorithm 14.85).  PARI's
``FlxqV_factorback`` uses the exponents as given.  As a Silex evaluation
choice, when ``f`` is squarefree modulo ``p`` and some exponent exceeds
``p - 1`` in magnitude, Silex first reduces the large exponents.  By the
Chinese remainder theorem, ``R`` is then the product of the fields
``F_p[x]/(f_i) = F_(p^d_i)`` over the distinct irreducible factors ``f_i`` of
``f`` (Lidl and Niederreiter, *Finite Fields*, Theorem 1.61), and every
``a`` in ``F_q`` satisfies ``a^q = a``, zero included (Lemma 2.3).  So with
``lambda = lcm(p^d_i - 1)``, the exponent of the unit group of ``R``, every
element satisfies ``a^(lambda + 1) = a``.  An exponent ``e > lambda`` is
replaced by ``((e - 1) mod lambda) + 1``, which gives the same power of every
element of ``R``.  Plain ``e mod lambda`` would not: it turns ``a^lambda``
into ``1`` for a zero divisor ``a``.  An exponent ``e < -lambda`` is replaced
in the same way only for a base that is a unit of ``R``.  The prime is
rejected before this point when a base with a negative exponent is not a
unit, as before.  The degrees ``d_i`` come from FLINT's
``nmod_poly_factor_distinct_deg``.  When ``f`` is not squarefree modulo
``p``, the exponents are used as given.  The residues are the same elements
of ``R`` either way, so no reconstructed unit changes.  The reduction saves
work only when an exponent exceeds ``lambda``.

At the first prime above ``2^63`` the 41st cyclotomic polynomial splits into
8 quintics, so ``lambda = p^5 - 1`` has 316 bits and the wide unit exponents
of ``cyclotomic_41_phi40`` (at most 136 bits) are not reduced.  Silex
therefore chooses the CRT primes differently from PARI's
``init_modular_big`` sequence when some exponent has more than ``FLINT_BITS``
bits: for each prime the CRT needs, it scans at most ``4 * degree``
candidates (a Silex limit, not a sourced constant) and takes the first prime
``p`` with ``x^p = x`` modulo ``f``.  Such a prime has ``f`` split into
distinct linear factors, so
``lambda = p - 1`` (Lidl and Niederreiter, Theorem 1.61 and Lemma 2.3); by
the Chebotarev density theorem a positive proportion of primes splits a
Galois polynomial.  If the scan finds no such prime, the next prime of the
PARI sequence is used and the preference stops.  The skipped primes are
counted in the internal report field ``nonsplit_prime_skips``.  For a split
prime the unit group exponent is at most 64 bits, so a negative
exponent is also made positive (``e`` plus a multiple of ``lambda``) when the
sliding-window cost model finds that cheaper; a small negative exponent is
otherwise kept, since the positive form would need about 63 bits.  This
applies to any row with an inverse, at any prime where ``lambda`` is
available (non-split squarefree primes included), once some exponent in the
batch exceeds ``p - 1``; it affects cost only.  Rows with no exponent above ``FLINT_BITS`` bits use the
PARI primes and no scan.  The residues are the same elements of ``R``, so no
reconstructed unit changes.  For ``cyclotomic_41_phi40`` 11 primes are
skipped, ``lambda`` has 64 bits, and the exponent bits fall from 373k to
179k; for the degree-32 real subfield of the 136th cyclotomic field, 40
primes are skipped and the exponent bits fall from 159k to 150k.

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

**Deviation** (decided, kept): when exact rounding fails, Hecke v0.38.6
``reduce_mod_units`` computes ``round_approx``, sets ``exact = false``, doubles
the precision and takes an unconditional ``continue``, so the approximate
coefficients are never applied (checked against the v0.38.6 source,
``Map.jl`` lines 638-658).  Silex applies the approximate coefficients and
retries at the same precision.  The loop ends when an exact rounding succeeds,
when the coefficients are all zero (whether from exact or approximate
rounding), or, through the branch where the approximate rounding fails, when a
half-integer tie stops it or the cap is reached; only that failing branch
changes the precision.  The all-zero exit is where Silex can be less reduced
than Hecke: an approximate rounding of a ball such as 0.5 +- 0.45 gives zero,
where Hecke would double the precision and might round exactly to a nonzero
value.  The function also reduces class-group relation generators and unit
candidates before the torsion test, so the argument is not about a unit basis
only: each output is its input times an exact product of the current unit
generators (exact ``FactoredElement`` arithmetic), so it generates the same
principal ideal and, with the current units, the same group.  Correctness does
not depend on which integer coefficients are applied; only which
representatives are kept differs, so the published generators and basis can
differ from Hecke's.  Silex keeps this
behavior on purpose.

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
factor base or discharged by the principal ideal ``pO_K``; proven units and
regulator; and a verified ``ell``-local proof for every prime ``p`` dividing
``h_cand``.  When every prime ideal above ``p`` has norm at most the bound
and the last one in the decomposition is unramified, that one is
discharged: ``pO_K = P_1^e_1 ... P_g^e_g`` with ``e_g = 1`` puts ``[P_g]``
in the subgroup generated by the others.  This includes an inert ``p``
(``g = 1``), whose only prime ideal is ``pO_K`` itself and so is principal;
the relation-completion base omits inert primes.  PARI 2.17.3
``src/basemath/buch2.c:be_honest`` and ``buch2.c:bnftestprimes`` skip the
same last unramified ideal, also when the list of primes is cut off at the
norm bound (Silex skips it only for a complete decomposition), and ``buch2.c:FBgen`` keeps inert primes out of
the factor base; see also Cohen, *A Course in Computational Algebraic
Number Theory*, GTM 138, Section 4.8, for ``sum e_i f_i = n``.  A totally
ramified ``p`` with a single prime ideal stays a target, because ``pO_K``
gives only ``e [P] = 0``.  The required primes are derived
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
norm at most ``sqrt|d| exp(-(F + n log pi)/2)``.  Silex evaluates ``F_{a,b}`` with a certified enclosure, using the theorem
GB-zimmert-F-tail-bound of the silex-research repository,
``research/proofs/generation/GB-zimmert-F-tail-bound.md`` at commit
``06191ff``.  With ``a = r1 + r2``, ``b = r2``, ``c = 1 + 2 gamma`` and
``h = 2c``, its part (iii) gives, for an integer ``L > gamma/2`` and every
``gamma`` in ``(0, 2L)`` that is not an integer, ``F = S_L + a sum_{l>L}
T(l) + b sum_{l>L} T(l + 1/2)``, where ``S_L`` is a finite sum of digamma
values at positive arguments and ``T`` is the term of Zimmert's series.  Its
part (i) bounds each tail ``sum_{l>L} T(l + delta)``, ``delta`` in
``{0, 1/2}``, by an explicit main term ``mu(N, K)`` built from exact
rational coefficients ``Phi_k`` (polynomials in ``gamma`` given by the
theorem's recurrences) and an explicit error ``E(N, K)``, for any admissible
``(L, delta, K)``; the bound holds for every real ``gamma > 0``, integers
included, and part (iii) shows that ``F`` extends continuously to integer
``gamma``.  The evaluation follows the theorem's procedure "Evaluation with a
certified enclosure" with Arb at 128 bits, every digamma argument being an
exact positive rational, and ``L`` and ``K`` chosen as in "Choosing L and
K" (``L`` between 190 and 294 and ``K`` between 32 and 33 for the table
below), with the admissibility conditions rechecked in exact arithmetic
and the size of ``E`` rechecked through the upper endpoint of an Arb ball.  The enclosure is used only
offline: the constant depends on the signature alone, so Silex stores a precomputed, test-verified table.
For each signature with ``3 <= n <= 20``, the generator
``tools/zimmert_twin_class_table.cpp`` (target
``silex-zimmert-twin-class-table``) takes the upper endpoint of the Arb
enclosure of ``exp(-(F + n log pi)/2)`` and rounds it up to an exact decimal
rational ``E`` with 30 significant digits, written to
``src/factor_base/zimmert_twin_class_table.inc``; to regenerate, build that
target and redirect its output to that file.  At run time the twin-class
bound ``B`` is ``ceil(E sqrt|d|)``, computed exactly with integers as
``ceil(ceil(s sqrt|d|) / 10^k)`` for ``E = s / 10^k``.
``test/t-generation-bound.cpp`` recomputes every entry with the Arb
evaluation and checks that it is on the safe side (``E`` at least the upper
endpoint ``U``, so ``B`` is at least Zimmert's ``sqrt|d| exp(-(F + n log pi)/2)``) and tight (``E <= U (1 + 10^-25)``).  PARI's
``zimmertbound`` (below) likewise hard-codes a constant ``c`` per
``(n, r2)``, but as a C ``double`` given to about nine digits with no stated
rounding direction, and evaluates ``exp(-c) sqrt|d|`` in floating point.
Every ``gamma > 0`` gives a valid bound, so ``gamma`` affects only its size.
For each signature with ``3 <= n <= 20`` Silex fixes the ``gamma`` in
hundredths, ``0 < gamma < 3``, that maximizes ``F``
(``test/zimmert_twin_class_oracle.py --gammas``); integers are excluded
from that search in the oracle, which cannot evaluate integer ``gamma``.
The integer case was checked once, outside the repository, with
silex-research's ``zimmert_F``: for all 117 signatures, ``F`` at
``gamma = 1`` and ``2`` never exceeds the tabulated value.  Relative to the earlier
finite-``m`` evaluation, no table entry is larger and every entry with ``gamma != 1/2``
is smaller, by at most about ``1.3e-4`` relative.
``test/t-generation-bound.cpp`` checks the constant
``exp((F + n log pi)/2)`` against all 16 rows of Tabelle 3 (p. 377) with
``3 <= n <= 20`` (the printed value is a lower bound within one unit of its
last digit) and pins ``B`` at ``|d| = 10^30``;
``test/zimmert_twin_class_oracle.py`` is the independent mpmath evaluation
those values come from.

The twin-class bound is a generation bound with no hypothesis on the class
``[D]`` of the different.  Zimmert's Korollar 1 alone gives generation only
when ``[D]`` lies in the subgroup ``H`` generated by the prime ideals of
norm at most the bound, since ``K`` or ``[D] K^(-1)`` contains such an
ideal.  The theorem GB-kor1-generation of the silex-research repository,
``research/proofs/generation/GB-kor1-generation.md`` at commit
``8926120``, removes that side condition: for every number field of degree
``n = r1 + 2 r2`` with absolute discriminant ``d`` and every real
``gamma > 0``, the classes of the prime ideals of norm at most
``B(gamma) = sqrt(d) exp(-(F_{r1+r2, r2}(gamma) + n log pi)/2)`` generate
the class group (non-strict inequality).  Its proof combines Korollar 1
with a twin-class generation lemma resting on Hecke's theorem that the
class of the different is a square.  The proof file also states that any
certified upper bound ``B' >= B(gamma)`` at a fixed ``gamma`` may replace
``B(gamma)``.  Silex's bound is such a ``B'``: every entry ``E`` is at
least the upper endpoint of the certified enclosure of
``exp(-(F + n log pi)/2)`` at the tabulated ``gamma`` described above, so
``ceil(E sqrt|d|) >= B(gamma)``.  As an independent cross-check,
``tools/check_zimmert_twin_class_table.py`` (an offline tool that needs Sage
and silex-research commit ``06191ff``) compares every entry with
``c_up(r1, r2, gamma)``, the certified upper bound computed by
silex-research's ``src/silex_research/zimmert.py`` with Arb at 128 bits,
using exact rational comparisons; for all 117 signatures with
``3 <= n <= 20`` it finds ``E >= c_up``, with equality to the printed
precision at ``gamma = 1/2``.

Silex uses the twin-class route only when its bound is strictly below
``factor_base_class_group_bound``, for any defining polynomial: monic or
not, and whatever the index ``[O_K : Z[theta]]``.  Otherwise, and whenever
the degree is outside ``3..20`` or the discriminant or signature cannot be
computed, the requirement is ``factor_base_class_group_bound`` (Satz 2,
Minkowski or the integer form above) and the generation basis is
``minkowski_type``.  With the twin-class route the generation basis is
``zimmert_twin_class`` and the generation certification bound is the
twin-class bound; the detail log names the route ``kor1_generation`` (no
check of the different), and ``silex-class-unit-instance`` prints it as
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
with the exact class number ``h = h(D)`` in place of ``ell``-local tests.
Here ``D`` is the fundamental discriminant of the maximal order (the route
applies to maximal orders only).  For ``|D| < 2^44`` ``h(D)`` is counted from
the reduced forms returned by FLINT ``qfb_reduced_forms``; every form of
discriminant ``D`` is primitive, so the count of reduced forms is the class
number of the order.  The listing takes time and memory growing like
``sqrt(|D|)`` (about 1.4 GB at ``|D| = 1e14``), so for ``|D| >= 2^44``
``h(D)`` comes instead from Dirichlet's class number formula
``h(D) = w sqrt(|D|) L(1, chi_D) / (2 pi)``, with ``w = 2`` for ``D < -4``
(Davenport, *Multiplicative Number Theory*, 3rd ed., Chapter 6), evaluated
in ball arithmetic from the rigorous ``L(1, chi_D)`` enclosure described
below.  The ball contains ``h(D)``, and ``h(D)`` is taken only when the ball
contains exactly one integer, which must be positive.  The enclosure is
computed in up to three steps, each tried only when the previous one gives
no such integer:

#. ``L(1, chi_D)`` to an absolute tolerance sized for ``h(D)`` (below);
#. FLINT's precision mapping at the bit length of ``|D|`` plus 32 bits;
#. the same at twice that precision.

Any failure of the first step moves on to the second.  From the second step
on, a non-finite or non-positive enclosure fails closed, and so does a ball
that still holds more than one integer after the third step.  So the
computation succeeds wherever the last two steps alone succeed, with the
same ``h(D)``.

The first step sizes the work for ``rad(h) <= 1/16`` instead of for a target
precision of ``L(1, chi_D)``.  Since ``h = sqrt(q) L(1, chi_D) / pi`` with
``q = |D|``, this needs ``rad L <= pi / (16 sqrt(q))``.  The tolerance of the
first sum is ``2^t`` with ``t = floor(log2(pi / (16 sqrt(q))) - log2(N / 16 +
2) - 2)``.  The term ``N / 16 + 2`` allows for the per-term errors of the
``N`` terms (each about ``tol / 32`` under FLINT's per-term precision) and for
the two tails, and the last 2 bits are a margin.  ``N = sqrt(q z / pi) + 2``
estimates the truncation index, where ``z = max(1, (1 - t) ln 2 + ln(1 + q /
pi))`` is about where FLINT's tail bound falls below ``2^t``; ``t`` and ``N``
are iterated four times from ``N = sqrt(60 q / pi)``.  The second sum keeps
FLINT's tolerance relative to the first, about ``2^t sqrt(q / pi)``.  Both
sums are accumulated at ``-t + log2 sqrt(q) + 8`` bits, at least 32.  FLINT's
tail bound and per-term working precisions are unchanged, so the enclosure
stays rigorous; the tolerance and precisions are cost parameters only.  This
step deviates from FLINT's mapping of a target precision ``p`` to the
tolerance ``2^(mmax - p')`` with ``p' = 1.001 p + 2 bits(q)`` (FLINT 3.6.0
``src/acb_dirichlet/l_fmpq_afe.c``, ``acb_dirichlet_l_fmpq_afe``), which the
last two steps keep.  At ``|D| = 1e14 + 319`` that mapping targets
``2^-174``, while the first step needs about ``2^-49``; the estimated
truncation index falls from about ``6.9e7`` to about ``4.6e7``, and each term
is evaluated at about one limb of precision instead of three.  The sources
consulted select no such tolerance: PARI 2.17.4 ``classno2``
(``qfbclassno(D, 1)``) truncates the series of Cohen, *A Course in
Computational Algebraic Number Theory*, Prop. 5.3.14, at the index of
Cor. 5.3.16, works at 64 bits without a tail bound, and rounds.  The rule
is therefore a Silex cost parameter, adopted as an approved change of the
precision rule, and not a source's precision rule.

At ``|D| = 1e14 + 31`` a direct call of this function gives
``h = 14414435``, PARI's value.  A full proven run reaches this step at such
sizes only through the structure certificate route described below.
With generation verified up to the Minkowski-type bound and ``h_cand = h``,
the computed relation lattice has index one in the full one, so it is
saturated at every prime.  The route stores a verified relation-saturation
record for every prime ``ell`` dividing ``h`` (none when ``h = 1``), backed
internally by a proof record that names the exact class number, not an
``ell``-local test, as its basis.  A prime that already carries a verified
``ell``-local proof keeps it.  Later gates, such as the Belabas--Friedman
audit, therefore see the same per-prime saturation records as for a
saturation-backed proof and can promote from them.

The relation modules of a class-group context keep dense ``n x n``
matrices modulo a word-size prime for a factor base of ``n`` generators
(three per context, 2 GiB each at ``n = 2^14``), and FLINT aborts the
process when such an allocation fails.  A factor base with more than
``2^14`` generators is therefore refused before anything is allocated:
building or installing it fails, and a class/unit transaction fails with
reason ``factor_base_size_limit_exceeded``, as ``silex-class-unit-instance``
prints in ``failure_reason``.  The limit is a fixed Silex choice with no
upstream source, so the same fields fail the same way on every machine.  The
proven base of the exact imaginary-quadratic route reaches the Minkowski-type
bound ``sqrt(|D| / 3)`` and passes the limit between ``|D| = 4e11`` and
``5e11``.  There a ``proven`` request keeps the compact relation base and is
certified by the structure certificate instead (next paragraphs); the
refusal reason stays for every other route that needs a base over the limit.

*Structure certificate (imaginary quadratic, ``proven``).*  When the base up
to the generation bound would exceed the limit, the route computes the class
group from the compact base (the staged relation rows with their exact
witnesses are the published relations; no relation outside the stage is
added) and then checks, in exact reduced binary quadratic form arithmetic,
that the published invariant generators ``g_j`` (products of factor-base
primes, with invariants ``d_j``) satisfy: the product of the ``d_j`` is the
exact class number ``h`` described above; each ``g_j^(d_j)`` is principal;
and for every prime ``ell`` dividing ``h`` no nonzero
``F_ell``-combination of the ``g_j^(d_j / ell)`` with ``ell | d_j`` is
principal.  The induced map from the direct sum of the ``Z/d_j`` to the
class group is then injective between groups of order ``h``, hence bijective,
so the ``g_j``, and with them the factor base, generate the class group.
The ideal-to-form isomorphism is Cohen, *A Course in Computational
Algebraic Number Theory*, GTM 138, Section 5.2.  The stopping rule, that
classes generating a subgroup of order exactly ``h`` generate the class
group, is Cohen, GTM 138, Algorithm 5.4.1 and its proof (pp. 240--242), and
Mosunov and Jacobson, *Math. Comp.* 85 (2016), Section 5; the class number
is the exact one above.  A relation lattice of index ``h`` alone
would not suffice, because the factor base may fail to generate the class
group (the ``Cl(-199)`` example of the route's design), which is why the
``ell``-step test is required.  The route publishes
``factor_base_generation_basis() == structure_certificate`` with
certification ``proven`` and no generation bound
(``factor_base_generation_certification_bound`` is false).  For
``x^2 + 1000000000039`` a ``proven`` request now succeeds in under a second
with ``h = 1113261`` (the factor base would have 23592 generators).  The
largest tested proven run on this route is at ``|D|`` of order ``1e12``;
``h(D)`` alone takes about 100 s at ``|D| = 1e14`` and its cost grows like
``sqrt(|D|)``; ``fmpz_fits_si`` limits ``D`` to ``|D| < 2^63``.  When the
structure certificate fails, or the relation index is below ``h``, the run
fails closed with the generic reason ``class_unit_computation_failed``, with
no fallback.

The proven imaginary quadratic witness step needs, for each factor-base row,
the Hermite normal form of the row's ideal in the order basis ``(1, omega)``
of the maximal quadratic order.  Silex builds it directly instead of
multiplying ideals, following Cohen (GTM 138): a prime above ``p`` is
``p Z + (omega - r) Z`` with its conjugate (Prop. 5.1.4); the product of
ideals of coprime norm is the primitive ideal ``A Z + (omega - s) Z`` (Lemma
5.4.5 with ``d = 1`` and Lemma 5.5.1), and an ideal ``m`` times a primitive
ideal of norm ``A`` has norm ``m^2 A`` (Prop. 5.2.1).  The root ``s`` is the
Chinese remainder combination of the roots of the individual prime powers,
each obtained from the stored root modulo ``p`` by the quadratic Hensel lift
(Thm. 3.5.3, Alg. 3.5.5 and 3.5.6).  The builder relies on one contract:
the factor-base prime of each column is ``p Z + (omega - r) Z`` for its
stored root ``r``, and the conjugate column is ``p Z + (omega - (t - r)) Z``.
Once per prime the code checks this by comparing the prime's Hermite normal
form with that of ``p Z + (omega - r) Z``; a row that needs the root of a
prime failing the check is built by the ideal route.  Before the result is
used the code checks that ``A`` divides ``N(omega - s)`` and that ``s`` is
the selected root modulo every ``p`` dividing ``A``; with the norm ``A`` this
identifies the lattice with the product.  If the builder declines or a check
fails, the witness step uses the earlier ideal-arithmetic route, so the
proven statement is unchanged.

The reduction of the very large binary quadratic forms that the proven
imaginary quadratic witness step produces (forms of more than 9000 bits)
is a port of the recursive fast reduction of PARI 2.17.3
``src/basemath/Qfb.c`` (``pqfbred_1`` through ``pqfbred_rec``, lines
655--843, which implement Schönhage, "Fast reduction and composition of
binary quadratic forms", ISSAC 1991, pp. 128--133), with the switch to it
taken from ``qfbred_imag_av`` (lines 991--1002) and ``qfbredsl2_imag``
(lines 886--907).  Below the switch the earlier basecase reduction runs
unchanged.  The port deviates from PARI in one place: for a form with
``b < 0`` PARI negates the second row of the transformation matrix, which
gives determinant ``-1`` and the reduced form of the wrong class; Silex
conjugates the matrix by ``diag(1, -1)`` instead, so the determinant stays
``1``.  If a step of the fast path meets the condition where PARI raises
its internal error, the basecase reduction of a step fails, or the input
form does not take the reduced minimum at the computed vector, the
reduction falls back to the full basecase reduction.  Reduced forms, and so
class-relation rows, norms and labels, are those of the basecase; the
generator of a relation can differ from it by sign above the switch (by a
unit when ``D = -3`` or ``-4``).

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
from scratch.  The analytic imaginary-quadratic class number above also
evaluates the same sums with a tolerance sized for ``h(D)`` instead of
FLINT's precision mapping (first step of its ladder).

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
The Belabas--Diaz y Diaz--Friedman bound is a cutoff at which an Arb
comparison proves the criterion of Belabas--Diaz y Diaz--Friedman 2008,
Section 3; it is the smallest such integer cutoff unless a near-tie that
Arb cannot decide at a bisection midpoint leaves a larger proven one.  Silex finds it as PARI 2.17.3
``src/basemath/buch2.c:Buchall_param`` (lines 3810-3817) finds its
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
Theorem 1).  In ``grh`` mode the test first evaluates ``hR`` coarsely: the
Euler-product cutoff is the smallest one whose Theorem 7 error bound
(Belabas--Friedman 2015, Theorem 7, with ``beta`` from equation (16) and
``sigma = 1 + 1/sqrt(log D)``) is below 1/4 - 2^-20, an absolute error on
``log hR`` that is enough to decide an index below two (footnote 1 of the
paper; PARI ``buch2.c:primeneeded`` uses the same 0.25).  The margin 2^-20
is reserved for numerical error.  The coarse value is used only when it
shows index one.  If the scan finds no cutoff up to 20007, or the coarse
value does not show index one, it is discarded and the default zeta route
computes the value (so such fields pay up to part of one extra analytic
evaluation).  A field whose default route fails can therefore be accepted
under GRH when the coarse value shows index one.
A ``grh`` request takes the coarse value, and the default route's
Belabas--Friedman ``hR`` when the coarse value does not show index one, in
every degree above one, including a quadratic field written as ``x^2 - d``,
which also has the unconditional ``L(1, chi_D)`` route; ``proven`` requests
and the public zeta gates keep that route first.  Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp.jl:_class_unit_group`` takes the
Belabas--Friedman residue (``zeta_log_residue``, ``Zeta.jl``) for every
degree including two.  PARI 2.17.3 ``quadclassunit``, which is
GRH-conditional, approximates ``hR`` by a truncated Euler product with no
explicit error bound (``src/basemath/buch1.c:compute_invresquad``).
Neither has an unconditional evaluation on this path.  A ``grh`` result
does not use the unconditionality of ``L(1, chi_D)``, whose cost grows like
``sqrt(|D|)``, while the Belabas--Friedman cutoff is set by the precision:
on the real quadratic benchmark rows from ``|D| = 840`` to ``4e5`` the
coarse value took under about 1 ms, against 35 to 120 ms for
``L(1, chi_D)`` at ``|D| = 4e4`` to ``4e5``.  Silex has no ``|D|``
crossover, so the recorded completeness basis of a ``grh`` result does not
depend on ``|D|`` or on how the field is presented.  Such an acceptance is not a
relation-saturation proof: ``relation_saturation_status()`` is not set to
``verified`` by it.
The pair is published with both labels ``grh``, and the analytic check that
accepted it is recorded on the published class group with its own
conditionality.  ``analytic_class_regulator_certification()`` is ``grh``,
since the accepting ``hR`` is a Belabas--Friedman value (the default zeta
route's Belabas--Friedman evaluation), and the BF audit data of that
evaluation is recorded as well: its error bound, cutoff, maximum cutoff, and
precisions (``zeta_bf_proof_record()``).  The record is
informational.  The ``grh`` transaction never changes or promotes the
class-group or unit labels, and the record does not change the public zeta
and BF gates, which still fail closed in degree three and higher unless
saturation is proven.  A later explicit
``try_certify_with_units(units, proven)`` with proven units can still promote
a real-quadratic ``grh`` class group to ``proven``: it evaluates the
unconditional ``L(1, chi_D)`` value itself, and the index-one test against
it, with generation verified up to the Minkowski-type bound, gives
``h_cand = h`` and ``R_cand = R_K`` unconditionally.  A ``grh`` run whose
base covers only the GRH bound has no unconditional generation, so that
promotion fails for it.  For an
imaginary quadratic field with ``|D| < 2^34`` the ``grh`` route uses an exact
class number, counted from the reduced forms returned by FLINT
``qfb_reduced_forms``, instead of an analytic value.  Its GRH dependence is in
factor-base generation, which is checked only up to the GRH bound, and it
records no analytic check.  Listing the reduced forms takes time and memory
growing like ``sqrt(|D|)`` (about 1.4 GB at ``|D| = 1e14``), so from
``|D| >= 2^34`` an imaginary-quadratic ``grh`` request takes ``h`` from the
GRH relation lattice and the analytic index-one test, like every other degree;
there the test uses a Belabas--Friedman ``hR``.  The cutoff ``2^34`` was
measured on an idle machine over ``|D|`` from about ``2^30`` to ``2^44``
(:doc:`../development/class_number_route_cutoff`).  The
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
requirement, when ``factor_base_generation_status()`` is ``verified``, or
with basis ``structure_certificate`` and no generation bound on the
imaginary quadratic structure certificate route above, where the status is
not ``verified``; it is ``grh`` when only
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
accepted against an unconditional ``hR`` (``unconditional_analytic``, which no
public transaction reaches now that a ``grh`` request takes the
Belabas--Friedman value) or the exact imaginary-quadratic
class number
(``exact_class_number``; reported only by a ``grh`` request on an imaginary quadratic field with ``|D| < 2^34``), and for a ``proven`` class group or
the exact
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
   * - Imaginary quadratic, ``|D| < 2^34``, exact class number
     - ``grh``, or ``proven`` when the base reaches the proven bound
     - ``proven`` (``exact_class_number``)
   * - Imaginary quadratic, ``|D| >= 2^34``, Belabas--Friedman
     - ``grh`` or ``proven``
     - ``grh`` (``belabas_friedman``)
   * - Real quadratic, Belabas--Friedman (any presentation and ``|D|``)
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
A caller of ``try_prove_relation_saturation_with_units`` chooses that bound.
Inside a ``proven`` paired class/unit transaction, the saturation proof run at
analytic index one has no fixed bound: as
in ``primecertify``, which walks every prime ``q = 1 mod p`` up to
``ULONG_MAX`` until the matrix has full rank (2.17.3 ``buch3.c`` lines
1086--1146), the test runs over the whole machine-word range of ``q``.  The
pair it proves passed the Belabas--Friedman index-one test, so under GRH its
relations are saturated and the scan stops at full rank.  Every auxiliary
prime it uses is checked as before.  Silex has no time limit of its own: a
request whose scan runs longer is stopped only by the caller's time or memory
limit, and a stopped request publishes nothing, so no result is labelled
``proven``.  For example, ``Q(sqrt(1575029))`` (``h = 89``) needs
``q = 2671``: the only ``q = 1 mod 89`` below 1000 is 179, which is inert.
Each certification entry point is transactional; a failed call restores all
class-group certification metadata, including the unit and regulator proof
statuses.  A successful ``try_prove_relation_saturation_with_units`` call
means only that its ``ell`` verified; the relation-saturation status becomes
``verified`` only once every ``p | h_cand`` is covered.

The ``beta`` row of an invariant generator ``g_i`` is formed from its
relation coefficients: a vector ``c_i`` for which ``c_i R``, with ``R`` the
relation matrix, is the factor-base exponent vector of ``g_i^{d_i}``, so the
product of the relation generators to the powers ``c_i`` generates
``g_i^{d_i}``.  The presentation's ``c_i`` comes
from the Smith transform and can exceed the reduced size by orders of
magnitude (figures below), which makes the
compact witness, and its images at the auxiliary primes, too large to use.
Measured on commit 49dfe1d with a temporary, uncommitted print of the
maximum coefficient bit size, on a busy machine: over the 11 beta rows
reduced in ``t-silex-class-unit-matrix``, the raw coefficients are at most
414 bits (286 to 414 bits for the rows with 42 relations) and the reduced
ones at most 12 bits (2 to 12 bits for those rows).  Larger presentations
were not measured.
Silex size-reduces ``c_i`` against an LLL-reduced basis of the relation
kernel before forming the witness.  If the basis that yields the unit
candidates is already built, it is reused.  Otherwise, for a kernel with at
least two rows and entries of at most 250 bits, the basis comes from the
certified first LLL phase described with ``Lat::lll_reduce`` above:
``fmpz_lll_d_with_removal_knapsack``, accepted only when the exact
certificate holds with FLINT's default parameters, and ``fmpz_lll`` (whose
result is also kept for the unit candidates) when the reducer fails or the
certificate fails.  For a kernel with at least two rows and entries of more
than 250 bits, where ``fmpz_lll`` of
the kernel itself took most of a 90 s ``proven`` run of
``real_quadratic_d1000000000009`` (38 rows of up to 151685 bits, busy
machine), the basis comes from the
small kernel basis of the ``proven`` analytic finish above: the reduction
of ``[c*A | I]`` for the presentation's relation matrix ``A``, with the
same two exact checks, accepted only when ``fmpz_mat_is_reduced`` certifies
the kernel rows with FLINT's default parameters; ``fmpz_lll`` is used when
the reducer fails or a check or the certificate fails.  A kernel with one
row takes ``fmpz_lll``.  The reducer's row operations are unimodular, and
the checked rows of ``[c*A | I]`` are a kernel basis, so the result is a
basis of the kernel on every route.  With FLINT's default routing a certified first phase
is exactly the basis ``fmpz_lll`` returns for that kernel (``fmpz_lll`` keeps
it when ``fmpz_lll_is_reduced`` holds), so for kernels of at most 250 bits
the basis does not depend on the route.  The ``[c*A | I]`` route returns an
LLL-reduced basis of the same lattice that need not equal that one, so the
reduced ``c_i`` below can differ by a kernel vector from those of the
``fmpz_lll`` basis.  As shown below, this does not change the rank modulo
``ell`` at any set of auxiliary primes usable on both routes.  Whether an
auxiliary prime ``q`` is usable also depends on which relation generators
occur in ``beta_i`` (a generator that vanishes at the prime, or whose
denominator is divisible by ``q``, makes Silex skip it), so the number of
primes the proof uses is not shown to be route-independent.  Within the
auxiliary-prime bound, the proof's status (verified, or unavailable within
the bound) and the relations appended after it can then also differ.  The
``[c*A | I]`` route runs only when the ``fmpz_lll`` basis has not been built
yet; when the unit candidates were requested first, that basis is reused.
So on kernels with entries above 250 bits the proof record can depend on
whether kernel units were requested before the saturation proof.  Soundness
and the published labels do not depend on the route.  On
``real_quadratic_d1000000000009`` the proof records were identical to those
of the ``fmpz_lll`` route.  With the integral Gram--Schmidt data of that
basis, Silex subtracts
``round(mu)`` times each basis row, last row first, as in sub-algorithm RED
of Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138,
Algorithm 2.6.7, so that every Gram--Schmidt coefficient of the result has
absolute value at most ``1/2``.  This changes ``c_i`` by a kernel vector,
which multiplies ``beta_i`` by a unit ``u`` and keeps
``(beta_i) = g_i^{d_i}``.  PARI 2.17.3 ``buch3.c:check_prime``
(lines 1159-1180) fixes its ``beta`` list (``cycgen``) only by that ideal
identity, so ``beta_i`` is determined up to a unit there as well.  The
reduction is sound only when ``u`` lies in the group generated by the rows
in the local system together with ``ell``-th powers.  With proven units,
``u = zeta^a prod eps_j^{b_j}`` for the free-unit rows ``eps_j``.  Then
``zeta^a`` is either covered by the torsion row (``ell | w``) or is an
``ell``-th power with character zero, so the character row of ``beta_i``
changes by a unimodular row operation modulo ``ell``, and the rank modulo
``ell`` is unchanged.  Silex therefore applies the reduction only after the
units are proven, and refuses it otherwise.

The order-unit saturation proofs build degree-one primes ``(q, theta - r)``
from roots ``r`` of the defining polynomial ``f`` modulo ``q`` with an
internal constructor.  It accepts a root when ``f(r) = 0 mod q`` and ``f``
is squarefree modulo ``q``.  This is a sufficient condition for
``(q, theta - r)`` to be a prime of residue degree and ramification index one:
``q`` then does not divide the index of ``Z[theta]``, so the Dedekind-Kummer
theorem applies (Cohen, *A Course in Computational Algebraic Number Theory*,
GTM 138, Theorem 4.8.13).  It is not a necessary condition, and the
constructor refuses every root that fails it, even when ``(q, theta - r)`` is
still such a prime.  For example, ``f = x^3 - 3x^2 + 2x - 8`` has
discriminant ``-1724``, field discriminant ``-431`` and
``[O_K : Z[theta]] = 2``.  Modulo 2, ``f = x^2 (x + 1)`` is not squarefree,
yet ``r = 1`` is a simple root mod 2 (2 does not divide ``f'(1)``), which
by itself makes ``(2, theta - 1) O_K`` a prime with ``e = f = 1``, one that
the constructor refuses.  (Hensel's lemma lifts ``r`` to a root of ``f`` in
``Z_2``, so ``f`` has a linear factor over ``Q_2``, which gives a prime above
2 with ``e = f = 1``.)  The constructor's squarefree condition is the same
condition as the degree-one fast path of ``decompose_prime`` uses.

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
``ell >= aux_bound``.  The pre-scan, its root attempts and the bound are a
Silex addition; ``saturate!`` has no pre-scan and tests candidates only after
its own scan.  In ``prove_index_bound``, when the pre-scan leaves a nonempty
kernel that it has not certified, Silex first runs the first scan of the Hecke
loop (``stable = 3.5``) on the same group, before testing any pre-scan row for
a root.  If that kernel is empty, ``ell`` is verified with that scan's local
primes: the group is then ``ell``-saturated, so no pre-scan row has a root to
adjoin, and the path below would end in the same scan.  Otherwise the
pre-scan rows are tested and the Hecke loop runs as before.  A kernel row whose ``ell``-th root Silex cannot
decide (the compact root and the factored root both leave it undecided, for
example because ``Element::is_power`` decides in degree 10 and higher only
by residue disproof, never by a root) is treated like a row
that is not an ``ell``-th power: the scan doubles ``stable`` and continues.
This follows ``saturate!`` (Hecke v0.38.6 ``Saturate.jl``, lines 396-411),
whose ``is_power`` reports ``false`` whenever it finds no root, and it does
not change the certification rule: ``ell`` is verified only by an empty
character kernel, and the group grows only by an exactly verified root.  A
row that is an ``ell``-th power whose root Silex cannot compute stays in
every kernel, so its proof fails closed as ``unavailable`` only at the
resource guard below.  For ``ell >= 3`` in degree 10 and higher,
``Element::is_power`` finds a root only of a rational constant, so a root
attempt can end only in a non-power, undecided (the scan doubles
``stable`` after either), or with the root of a true ``ell``-th power.
Until the first doubling of ``stable`` for such an ``ell``, the scan
therefore does not attempt a row that is not a formal ``ell``-th power of
its factored generators; it treats the row as undecided and rescans with
the doubled ``stable`` at once.  Rows of later scans are attempted as
usual.  The rows that can be rooted are the same: Hecke tests each
candidate before doubling, but a true ``ell``-th power stays in every later
kernel, so it still reaches the root test at a later scan.  Such a root is
then accepted at a larger ``stable``, possibly from another kernel row, so
the published unit generators and the ``local_primes`` of the unit proof
record may differ from those of an attempt in the first scan; the
certification rule is unchanged.  When the scan has used the whole prime
range, rows are attempted at once, as before.  The auxiliary bound is therefore not a proof
parameter; a unit proof record stores it as the pre-scan bound of the call,
not as the largest prime used.  Hecke starts its scan at
``next_prime(2^60)``; Silex starts at ``q = ell + 1``, which changes only
which valid characters are used.  Roots are adjoined only after the exact
order-unit check described below.

Each auxiliary prime contributes one column of ``ell``-th-power residue
characters, the quotient logarithms of the unit images, following Hecke
``Saturate.jl:mod_p`` (lines 6-46).  Hecke fills a table of all generator
powers when ``nbits(pp) < 19`` (line 32) and otherwise takes each logarithm by
``baby_step_giant_step`` (``ResidueRingMultGrp.jl``, lines 716-738, line 43 of
``mod_p``).  Silex works in the order-``ell`` quotient on word-size
degree-one primes, where its shared linear walk costs up to ``ell`` steps
against Hecke's ``pp`` table entries.  As an intentional deviation from the
19-bit cutoff (``Saturate.jl`` line 32), Silex takes the logarithms by
baby-step giant-step once ``ell`` has 7 or more bits and keeps the walk only
for smaller ``ell``: on the measured proven unit saturation rows, where
``ell`` ran up to 16 bits, baby-step giant-step was cheaper from 7 bits
on (busy-machine probe data, not a controlled measurement).  This is a cost choice
only: baby-step giant-step runs only when ``ell`` is prime and the quotient
generator has order ``ell``, so both routes return the unique logarithm in
``[0, ell)`` and the columns are the same.

``ResidueFieldQuotientLog::apply`` returns, for ``x`` in ``F_q^*`` with
``q = p^f``, the ``i`` in ``[0, ell)`` with ``zeta^i = x^((q - 1)/ell)``,
where ``zeta`` is the stored element of order ``ell``.  At residue degree
``f > 1`` with ``ell | p - 1`` it computes the same ``i`` through the norm
``N`` from ``F_q`` to ``F_p``: the conjugates of ``x`` are ``x^(p^j)`` for
``0 <= j < f``, so ``N(x) = x^((q - 1)/(p - 1))`` and
``N(x)^((p - 1)/ell) = x^((q - 1)/ell)``; the ``ell``-th roots of unity lie
in ``F_p``, so ``zeta`` is a constant and the degree-one walk applies.  For
``x = a(theta)`` with monic irreducible modulus ``m``, ``N(x) = Res(m, a)``,
taken with FLINT ``fmpz_mod_poly_resultant``.  The route follows PARI 2.17.4
``src/basemath/FpX.c``: ``FpXQ_issquare`` (lines 2490-2497) reduces the
``ell = 2`` test to a Kronecker symbol of ``FpXQ_norm(x, T, p)``, with the
comment ``Ng = g^((q-1)/(p-1))``, and ``FpXQ_norm`` (lines 2770-2778) is
``FpX_resultant(T, x, p)`` divided by ``lc(T)^deg x``.  Cohen, GTM 138,
Proposition 4.3.4 gives the resultant formula for the norm (stated over
``Z``; its proof uses only the conjugates of ``a(theta)`` and so carries
over to ``F_p``).  This is a cost choice only: it
replaces a power in ``F_q`` by a resultant and a power in ``F_p``, and
returns the same logarithm.  Other ``ell`` keep the power in ``F_q``.

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

The class-relation saturation with units
(``ClassGroupContext::saturate_relations_bounded_with_units`` and the
index-bounded pass of class/unit validation) scans the residue-character
kernels at one prime at a time over the relation generators, the free units
and, when ``ell`` divides ``w``, the torsion generator, and appends the first
kernel row whose product is an ``ell``-th power.  In each append round
it keeps a stack of the ``ell``-th-power residue characters of the same
inputs at degree-one primes ``q = 1 mod ell`` from ``next_prime(2^60)`` up,
and grows it lazily: a candidate that passes the divisibility check and that
the current stack does not rule out first adds one more prime, and is tested
only if the common kernel including that prime does not rule it out either.
The stack therefore never holds more primes than the unfiltered scan would
run power tests, and no clock enters.  The stack starts only after the
round has run as many power tests as there are inputs; until then every
candidate that passes the divisibility check is power tested.  This count
gate is a Silex heuristic, not in Hecke, whose ``saturate!`` stacks every
prime before its first power test: on small fields a failed power test costs
less than the degree-one prime the stack adds for it, and the count was
chosen by measurement on the native class-group benchmark rows.  Like the
stopping rule, it changes only how many power tests are skipped.  Stacking
stops by the rule of Hecke
v0.38.6 ``src/NumFieldOrd/NfOrd/Clgp/Saturate.jl``,
``compute_candidates_for_saturate`` (lines 130--203): when the common kernel
is empty or its dimension has stayed the same for more than ``stable`` times
the number of inputs.  ``stable`` starts at 1.5, the default of that
routine; ``saturate!`` (lines 419--509) calls it with ``stable = 3.5``, so
1.5 is a Silex choice.  ``stable`` doubles after a failed power test, as in
``saturate!`` (line 503), but only after a failed test on a row that the
stack let through once it had stopped with a nonempty common kernel; a
failure while the stack is still growing leaves ``stable`` unchanged.  The
doubling lets the stopped stack grow again.  Each append round builds a new
stack with ``stable`` back at 1.5, whereas the doubled ``stable`` of
``saturate!`` carries over to its later rounds.  These departures change
only how many power tests are skipped.  A row outside the common
kernel modulo ``ell`` is skipped without a power test, and an empty common
kernel ends the round.  Every character vanishes on ``ell``-th powers, so a
skipped row would have failed its power test, and a stack that stops at any
point is sound: the appended relations, and with them the published class
group, generators and units, are those of the unfiltered scan.  Hecke
instead tests the basis of the stacked kernel; Silex keeps its own candidate
order, so the stopping rule decides only how many power tests are skipped.
The filter is not used for non-maximal orders,
for a defining polynomial that is not monic with integer coefficients (Hecke
applies the index check that keeps the auxiliary primes regular only to such
polynomials), or for an ``ell`` that does not fit a machine word, and the
scan tests every row when the filter cannot be computed.

Within one round the per-prime kernels often repeat a row (the same
exponent vector of the relation generators, free units and torsion
generator).  Whether a row gives an appended relation depends only on that
vector and on the relations, units and factor base, which do not change
before the round appends: the divisibility check, the product, its power
test (deterministic) and the refinement check.  The filter only removes
tests of non-powers.  So the round records each row it has run without an
appended relation, by its exact entries, and skips that row at later primes
(a Silex memo, not in Hecke, whose candidates are the basis of the stacked
kernel).  This applies to undecided power tests too, such as the square
tests of degree-16 leftovers: from degree 10
(``kElementPowerRouteDegreeLimit``) ``Element::is_square`` can only disprove
a power, so a true power there fails the test and is recorded as settled.
The appended relations are those of the full scan.  A skipped row is not
counted as a power test, so the filter's count gate opens later or not at
all; like the stopping rule, that changes only which non-powers are tested.
The scan still visits every auxiliary prime up to the bound; it does not
stop early when only settled rows recur.

The class-relation saturation candidate search reads the degree-one primes
over an auxiliary prime ``q`` from the roots of the defining polynomial ``f``
modulo ``q``.  For a monic integral ``f`` it rejects index and discriminant
primes, so each simple root ``r`` gives the prime ``(q, theta - r)``.  For any
other ``f`` that check does not apply, and a root modulo ``q | lc(F) *
disc(F)`` (``F`` the integer numerator of ``f``, which keeps its content) need not give a
degree-one prime, so its character could rule out a true ``ell``-th power.
Hecke ``compute_candidates_for_saturate`` takes
``prime_decomposition(OK, q, 1)`` in that case.  Silex skips every such ``q``:
for ``q`` not dividing ``lc(F) * disc(F)``, ``theta`` is ``q``-integral and
``Z_(q)[theta]`` is ``q``-maximal, so the Dedekind-Kummer theorem applies.
Skipping an auxiliary prime only leaves the candidate space larger.  Monic
integral polynomials are unaffected.  When ``f`` has non-integral
coefficients the search fails (closed), because the roots of ``f`` modulo
``q`` cannot be read.

Unit candidates from the relation kernel are formed from an LLL-reduced
basis of the saturated kernel, following Hecke.jl ``find_candidates``
(commit 44b5409, ``FindUnits.jl:53-58``: ``solve_dixon_sf``, then
``saturate``, then ``lll!``).  Silex deviates in the input to the reduction:
it reduces the full saturated Smith kernel, whereas Hecke reduces the
saturated kernel of a random subset of the relations.

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

In prime degree ``n``, ``OrderUnitGroup::regulator_index_bound`` also uses a
bound that depends on the field discriminant ``d_K``.  It raises ``R_lower`` to
the larger of the bounds above and this one:

* M. Pohst and H. Zassenhaus, *Algorithmic Algebraic Number Theory*,
  Cambridge University Press, 1989, Chapter 5, (6.22) (p. 366), excerpted
  there from M. Pohst, *J. Number Theory* 9 (1977), 459--492, and
  *Abh. Math. Sem. Univ. Hamburg* 47 (1978), 221--231, after R. Remak, *J.
  Reine Angew. Math.* 167 (1932), 360--378: if ``K`` has no proper subfield
  other than ``Q`` and ``|d_K| > c``, then with unit rank ``r = r1 + r2 - 1``
  ``R >= ((3 log^2(|d_K|/c) / ((n-1) n (n+1) - 6 r2))^r 2^r2 /
  (n gamma_r^r))^(1/2)``.  Here ``c = n^n``, replaced by ``4^floor(n/2)`` for
  totally real fields with ``n <= 11`` (the remark after (6.22)).
  ``gamma_r^r`` is Hermite's constant to the power ``r``: the exact values
  ``1, 4/3, 2, 4, 8, 64/3, 64, 256`` for ``r <= 8`` (3.35a, p. 198), and
  Blichfeldt's upper bound ``(2/pi)^r Gamma(2 + r/2)^2`` (3.35b, p. 199)
  beyond, which only lowers the bound.

Silex certifies that ``K`` has no proper subfield other than ``Q`` only by
``n`` being prime, so the bound is never used in composite degree, even for
primitive fields.  ``d_K`` is the discriminant of the maximal order, so the
bound is used only when the unit group's order is known to be maximal; for
any other order, and whenever ``r = 0`` or ``|d_K| <= c``, the signature
bounds above are used unchanged.  The bound is evaluated with Arb from the
exact integers ``|d_K|`` and ``c`` and its lower endpoint is used, at the same
precision cap as the other terms.  The public ``unit_lower_regulator_bound``,
which takes only the field, does not use it.  PARI
(``src/basemath/buch3.c``, ``regulatorbound``) evaluates the same bound in
floating point and also applies it in composite degree after a primitivity test on a
partial factorization of ``d_K``; Silex does neither.  On the prime-degree
fields of ``test_field_lower_regulator_bound`` in ``test/t-unit.cpp`` (14 fields)
the bound never exceeds the regulator, which the test checks against certified
PARI regulators.  For
example, the index bound for the quintic of discriminant ``-401370255`` falls
from 1184 to 236.

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
The same bound covers kernel rows with an undecided root, which end a scan
without adjoining a root.  If a scan ever ends because no prime ``q`` is left
in the machine-word range, a row that is not an ``ell``-th power and a row
with an undecided root both fail closed at once, since a rescan would return
the same kernel.

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
primes above the required prime's rational prime, cannot be discarded, with
one exception.  A witness may also carry a positive valuation at a prime
ideal ``Q`` above a rational prime ``q`` that differs from the required
prime's, when ``Q`` is unramified (``e = 1``), the factor base holds at
least one prime ideal above ``q``, and ``Q`` is the only prime ideal above
``q`` outside the factor base.  The argument is the step used above for the
proof targets: ``qO_K = Q_1^e_1 ... Q_g^e_g`` with ``e_g = 1`` puts
``[Q_g]`` in the subgroup generated by the other ideals above ``q``, so the
witness still places the required class in the subgroup generated by the
factor base (Silex ``select_factor_base_proof_targets``; PARI 2.17.3
``buch2.c:be_honest`` and ``buch2.c:bnftestprimes`` skip the last
unramified ideal above a rational prime; Cohen, GTM 138, Section 4.8, for
``sum e_i f_i = n``).  Accepting such a witness is a Silex extension: the
cited PARI step only skips testing the last unramified ideal, PARI's factor
base keeps every prime ideal above its rational primes, and
``buch2.c:can_factor`` does not accept a witness through a prime ideal
outside the factor base.  An inert ``q``, with no prime ideal above it in the
factor base, is not accepted, as in ``buch2.c:can_factor``.  A completed
negative classification is final; unsupported or failed direct evaluation
falls back to the existing full ideal factorization predicate.

On the maximal order of a quadratic field the scan classifies each
rational prime ``p`` by the Kronecker symbol ``(D/p)`` of the field
discriminant and skips an inert ``p`` (symbol ``-1``), whose only prime ideal
``pO_K`` is never a proof target (Cohen, GTM 138, Proposition 5.1.4;
PARI 2.17.3 ``buch1.c:nthidealquad``).

This routing change preserves Silex's scalar/lattice candidate search,
random draws, bounds, proof-target selection, and certification labels.
It does not select the separate T2 witness search.  The relation-search
honesty checkpoint, however, searches witnesses with that T2 route, as the
reference does: ``buch2.c:Fincke_Pohst_ideal`` enumerates small elements
of the required prime ideal by T2 norm and stops after ``maxtry_FACT`` =
500 factored candidates, and ``be_honest`` retries with at most
``maxtry_HONEST`` = 50 random sub-factor-base twists.  The scalar/lattice
search has no candidate cap: in high degree it could run on one prime with
no smooth element until the caller's time limit.  Every T2 candidate passes
the same exact norm-accounting check, so a witness proves the same thing on
either route.  The T2 route also keeps the scalar/lattice route's fallbacks:
a candidate the direct test cannot evaluate is classified by the full ideal
factorization, and a missing sub-factor-base schedule, or a schedule entry
outside the factor base, ends the random twists without a witness rather
than failing the scan.  The noninstalled
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

For degrees 2 to 4 the direct honesty path rejects most candidates by the
norm-accounting criterion above, so the lattice searches apply that rejection
from the order-basis coordinates before any element is built.  The absolute
norm of ``a = sum c_i w_i`` is the determinant of the multiplication-by-``a``
matrix on the order basis (Cohen, *A Course in Computational Algebraic Number
Theory*, GTM 138, section 4.3, ``N(a) = det(x -> a x)``, and the determinant
is independent of the Q-basis, so the order basis may be used), assembled from
the precomputed table of ``w_i w_j`` coordinates; the layout follows PARI 2.17.3
``src/basemath/base3.c:zk_multable``.  The smoothness test is the one of
``buch2.c:can_factor`` and ``base4.c:Z_ppo``: dividing out the part of the
norm supported on the product of the factor-base rational primes and the
required prime's rational prime must leave a unit, and a nonzero norm with a
larger cofactor cannot be a smooth witness.  In degree 2, when the table
entries and coordinates are at most ``2^31 - 1`` in absolute value (on 64-bit
words) and the norm is below ``2^64``, the same determinant and cofactor loop
run exactly in machine words, using ``gcd(N, P) = gcd(N, P mod N)`` for the
first step; any other input takes the multiprecision path, so the decision is
the same.  The prefilter only rejects a
candidate that the direct path rejects by that same cofactor test; every
other candidate, every unsupported or unreadable case, and degree 5 and above
continue through the direct path unchanged, so accepted witnesses, search
order, bounds, random draws and labels are unchanged.

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
``failure_reason``; the failure itself and every label are unchanged.  A
factor base over the fixed size limit of ``2^14`` generators is refused
before any matrix is allocated, and the transaction then reports
``factor_base_size_limit_exceeded``.  A ``grh`` computation that fails
because a search cap stopped its LLL relation search reports
``relation_search_budget_exhausted``: the tentative search spent its
candidate cap (``ClassGroupComputeOptions::max_candidates``, counted per call
of the relation route) or reached the relation cap (``max_relations``,
counted over all relations of the candidate) before it had a presentation,
the relation cap refused a continuation slice, or a continuation slice spent
its candidate cap without adding a relation.  A ``proven`` computation
reports the same reason when it gives up after a cap was reached by one of
its relation searches (the relation-production prepass, the relation-search
driver, the imaginary-quadratic relation backend, or the LLL relation route
and its continuation slices, which the relation cap refuses once the
relation count has reached it), counting a search that ended at its goal
with a cap reached, or after the relation cap alone refused it a larger
relation-kernel target.  The rule is broad on purpose: a cap can leave a
search with too few relations for a later step, which then fails for its
own reason (``x^4 - 10x^2 + 1`` with ``--max-relations 12`` ends in a unit
rank defect and succeeds with 13), so a ``proven`` run that reached a cap
names it even when its final failure has another cause, and
``failure_detail`` then keeps that cause next to the cap.  Both caps are unlimited by
default; ``silex-class-unit-instance`` sets them.  The caps and these stops
are Silex's own bounds: Hecke v0.38.6
``src/NumFieldOrd/NfOrd/Clgp.jl:_class_unit_group`` and
``src/NumFieldOrd/NfOrd/Clgp/Main_LLL.jl:class_group_new_relations_via_lll``
keep searching at these points, and PARI 2.17.3
``src/basemath/buch2.c:Buchall_param`` enlarges its factor base instead of
failing.  Only the reason is new; the failure itself is unchanged.  Any other
computation failure reports ``class_unit_computation_failed``.  Every
candidate at every stage passes the same exact witness check.  For the
relation-search checkpoint, stage 0 is the capped T2 search, not the previous
search; the unchanged results of the 65 succeeding ``proven`` benchmark rows
were established by measurement, not by construction.

Within one class/unit transaction, Silex keeps a passed honesty scan and
reuses it, with the same receipt, for a later candidate whose factor base is
equal (``FactorBase::equal``) and whose active and required bounds are equal;
only passed scans are kept.  PARI's ``be_honest`` has no such reuse: it runs
once per ``bnfinit`` pass.

In the ``silex-class-unit-instance`` JSON, the top-level ``failure_detail``
object is the one place for structured failure detail; its keys depend on
``failure_reason``, and it is null for every reason other than
``factor_base_honesty_unwitnessed``, ``class_unit_computation_failed`` and
``relation_search_budget_exhausted``.  For ``class_unit_computation_failed``
it is an object with one key,
``cause``, free text naming the last unit-group validation outcome the
computation reached (for example ``index_one_publication_failed``); it is
null when the computation failed before validation ran.  For
``relation_search_budget_exhausted`` it is an object with the key ``cap``,
the option name of the cap that stopped the search (``max_candidates`` or
``max_relations``), preceded by the same ``cause`` key when a validation
ran.  For ``factor_base_honesty_unwitnessed`` it gives the
unwitnessed prime as ``p``, a decimal string such as ``"3"`` like every big
integer in this JSON, and ``residue_degree``, the route as ``search``
(``lattice`` or ``t2``), the last stage searched and the last stage allowed
as ``stage`` and ``max_stage``, and that stage's caps under the detail log's
names: ``radius``, ``twists`` and ``random_tries`` on the lattice route, or
``random_tries``, ``factor_attempts`` and ``element_steps`` on the T2 route.
The residue degree, stages and caps are JSON numbers.

After each small-norm pass, the relation-production prepass drops the
outstanding relation surplus.  PARI 2.17.3
``src/basemath/buch2.c:Buchall_param`` (lines 3938-3946) runs ``small_norm``
and then sets ``cache.end = cache.last``, so the surplus
(``RELSUP + RU-1``, passed to ``init_rel`` at line 3871 and included in the
initial request recorded at line 3872) is not demanded again; the next
request is the recomputed ``need`` from the current rank and unit deficits
(lines 4047-4061).  Silex sets
``cache.target_relation_count`` to the current relation count after each of
its two small-norm passes in ``src/class_group/relation_production.cpp``.

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
leftover factors whose exponents are not multiples of the root degree.  In
degree 10 and higher, the reduced power test for ``n >= 3`` first runs the
first admissible prime of the residue disproof of degree 10 and higher
(above) on the leftover ``r = prod g_i^(e_i)`` without expanding it.  It
reduces each factor modulo ``p`` and multiplies the residues modulo
``f mod p``.  This is valid because reduction from ``Z_(p)[x]`` is a ring
homomorphism and ``f`` is monic integral.  The screen walks the same primes
with the same admissibility as the scan of the expanded ``r``.  At a prime
where some factor is not ``p``-integral, it stops without deciding, since
the admissibility of ``r`` there is unknown.  It therefore answers ``false``
only where the scan of the expanded leftover would.  There, no other step of
``Element::is_power`` can answer "not a power", so every answer is
unchanged, except that a non-power the screen decides now returns "not a
power" where expanding the leftover would have failed for lack of
resources.  When it decides nothing, the leftover is expanded and tested as
before, with the full scan.  The one-prime limit is a resource bound only:
each prime costs a pass over all factors, which a true power pays without
benefit.

The compact power test of a factored unit or relation-saturation candidate
follows Hecke's ``is_power(::FacElem, n)`` (Hecke.jl commit 44b5409,
``src/NumField/NfAbs/CompactRepresentation.jl``, lines 370-432).  It removes
the formal exponent quotient, writing the input as ``a = s^n * r``.  When the
conjugate-log bound of ``r`` is at most the integer square root of that of
``a``, Hecke tests ``r`` with ``_ispower(r, n)``; that call drops the ideal
decomposition, so ``compact_presentation`` recomputes the factorization of
``r`` with ``factor_coprime`` (lines 18-20).  Otherwise Hecke tests ``a`` with
its decomposition.  Silex has no port of ``factor_coprime``: its test of
``r`` reduces from the unit ideal, which is valid only when ``r`` is a unit.
Silex therefore takes that branch only for a unit ``r``.  A saturation
candidate with a nonzero factor-base row always takes the branch for ``a``
with that row.  For a unit input, ``r`` is a unit when ``s`` is empty
(``r = a``).  Otherwise an ``r`` whose factored norm, the product of
``N(f)^e`` over its factors, is not ``+-1`` is not a unit, since units have
norm ``+-1``, and is rejected without being expanded; any other ``r`` is
expanded and checked to generate the unit ideal of the order.  When the
check fails, Silex tests ``a``.  This is a recorded deviation.  Both
branches decide the same question, since ``a = s^n * r``, so only the cost
differs.
Whether to port ``factor_coprime`` is a separate decision.

The split of each exponent ``v`` into ``v = q*n + e`` also differs.  Hecke
uses Julia's ``divrem`` (line 390), which truncates, so on the path with the
conjugate-log comparison (``easy = false``) ``e`` has the sign of ``v`` and
lies in ``(-n, n)``.  Silex floors (``compact_floor_divrem_positive`` in
``src/class_group/class_group.cpp``, used by ``compact_split_power``), so
``e`` lies in ``[0, n)``.  For a negative ``v`` that ``n`` does not divide,
the two leftovers ``r`` differ, and so can the conjugate-log comparison that
chooses between the branches.  The answer does not change, since both
branches decide whether ``a = s^n * r`` is an ``n``-th power.  This is an
existing deviation.

The split also applies to different entries.  Hecke iterates over the
``FacElem`` dictionary (lines 383-389), whose equal bases are already merged
into one entry with the summed exponent, and skips bases equal to one.
``compact_split_power`` splits each stored entry of the factored element on
its own, merging no equal bases and keeping bases equal to one.  For
example, ``f^1 * f^1`` with ``n = 2`` leaves ``r = 1`` in Hecke and
``r = f * f`` in Silex.  The leftover and its conjugate-log bound can
therefore differ, and so can the branch chosen.  When both sides decide, the
answer is the same, for the same reason as above.  Not merging has two
further consequences, neither of which affects a label.  In the ``f * f``
case Silex skips the structural-power shortcut that Hecke takes at once, so
it can fail where ``Element::is_power`` is unsupported while Hecke returns.
The root returned can also differ from Hecke's by an ``n``-th root of unity.
This is an existing deviation.

The logarithmic embedding follows Hecke's ``conjugates_arb_log`` for
``FacElem`` (``src/NumFieldOrd/NfOrd/FacElem.jl``): equal bases are merged
into one dictionary entry with summed exponents, each base is evaluated with
guard bits, and the sum is widened by Hecke's ``expand!``
(``src/Misc/acb_root_ctx.jl``) as the final rounding.  The published radius
of each logarithm is at most ``2^-precision``, but it is not monotone in the
working precision, so it may be larger or smaller than before merging while
remaining a rigorous enclosure.

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
