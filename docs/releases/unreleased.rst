Unreleased
==========

This page was built for Silex |silex-release| on the |docs-channel|
documentation channel.

It records user-visible changes since the development workflow moved to task
branches; earlier post-0.1.1 changes are not yet listed here.  Listed entries
are available on the ``dev`` channel and become the next release's notes when
that release is prepared.

Added
-----

- ``ClassGroupContext`` reports where GRH entered a class/unit result:
  ``factor_base_generation_certification()``,
  ``factor_base_generation_basis()`` and
  ``factor_base_generation_certification_bound()`` say whether factor-base
  generation is ``proven`` (the proven generation bound) or ``grh`` (the
  Belabas-Diaz y Diaz-Friedman or Bach bound, whichever was used), and
  ``class_unit_completeness_certification()`` and
  ``class_unit_completeness_basis()`` say whether relation and unit
  completeness is ``proven`` (exact class number, unconditional
  ``L(1, chi)``, or a ``proven`` certification) or ``grh``
  (Belabas-Friedman ``hR``).  New enums ``ClassGroupGenerationBasis`` and
  ``ClassUnitCompletenessBasis`` name the theorems.  Labels are unchanged.
  ``silex-class-unit-instance`` adds the five values to its ``class_group``
  JSON object.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.
- ``silex-class-unit-instance`` reports ``failure_reason``
  ``factor_base_honesty_unwitnessed`` (``failure_stage`` ``total``) when a
  class/unit computation fails because the last candidate attempt stopped at
  a required prime that the factor-base honesty witness search could not
  witness at its final caps, with no factor-base restart left.  Other
  computation failures still report ``class_unit_computation_failed``, and
  results, labels and successful runs are unchanged.  A new top-level JSON
  field ``failure_detail`` holds structured detail about a failure, with keys
  that depend on ``failure_reason``.  It is null except for
  ``factor_base_honesty_unwitnessed``, where it names the unwitnessed prime
  (``p``, a decimal string, and ``residue_degree``), the witness search
  route (``search``, ``lattice`` or ``t2``), the stage of its last search
  (``stage`` and ``max_stage``), and that stage's caps: ``radius``,
  ``twists`` and ``random_tries`` on the lattice route, or
  ``random_tries``, ``factor_attempts`` and ``element_steps`` on the T2
  route.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

Changed
-------

- **Breaking:** ``SUnitGroup::regulator`` now returns the regulator of the
  S-unit lattice, ``R_S = R h/h_S prod log N(P)`` for nonempty ``S``: the
  absolute value of a maximal minor of the log-S embedding of the free
  generators, so ``h_S R_S = h R prod log N(P)``.  It previously returned
  PARI ``bnfsunit``'s number ``R h_S prod log N(P)``, which equals
  ``R_S h_S^2 / h`` and agrees with ``R_S`` only when ``h = h_S^2``.  For
  example, ``Q(sqrt(-5))`` with ``S = {P2}`` now gives ``2 log 2`` (PARI and
  0.1.1: ``log 2``), and ``Q(sqrt(10))`` with ``S = {P3}`` gives
  ``3.99553...`` (PARI and 0.1.1: ``1.99776...``).  The index ``h / h_S`` is
  exact: the absolute determinant of the nonunit valuation HNF.  To compare
  with PARI, multiply by ``h_S^2 / h``.  Empty ``S`` still gives the ordinary
  regulator, and the value is still an Arb enclosure with the same proof
  status.  ``silex-class-unit-instance`` prints the new value as the S-unit
  ``regulator_midpoint``.  Sources: Bernard and Roux-Langlois (ASIACRYPT
  2020), Prop. 2.2.  See :doc:`../reference/sunit_groups`.

- The GRH generation bound of a ``grh`` class/unit request falls back to
  Bach's bound alone when the Belabas-Diaz y Diaz-Friedman criterion cannot
  be computed or decided (for example for a non-monic defining polynomial,
  or an interval comparison in its search that cannot be decided), and the
  generation record then names Bach's theorem
  (``ClassGroupGenerationBasis::bach``).  A ``grh`` request on a monic
  presentation whose criterion search cannot be decided used to fail closed
  at policy selection and now uses Bach's bound.  Proven results and the
  failure rules are unchanged; a ``grh`` request on a non-monic presentation
  still fails closed at policy selection, because the factor base needs a
  monic integral defining polynomial.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- **Breaking:** ``FiniteAbelianGroup::set_relation_matrix_with_hnf_basis``,
  public in 0.1.1, is no longer public API; it is now internal.
  ``FiniteAbelianGroup`` also has a new private member, so its layout changed
  and code compiled against an earlier header must be rebuilt.

- ``factor_base_class_group_bound`` now returns a smaller proven bound in
  degree 3 and above: the smallest of Zimmert's bound (Zimmert 1981,
  Satz 2, for degree at most 20), Minkowski's bound with ``(4/pi)^r2``, and
  the former bound, which used ``2^r2`` in place of ``(4/pi)^r2``.  The
  bound never exceeds its former value; for example it drops from 463 to
  295 for ``x^3 + x + 200`` and from 43837247 to 6135 for ``x^19 - x - 1``.
  Proven class/unit computations therefore verify factor-base generation
  over fewer prime ideals; their results and certification labels are
  unchanged.  Degree one and quadratic bounds are unchanged.  See
  :doc:`../reference/algorithms_and_sources`.

- Proven class/unit computations in degree 3 to 20 can verify factor-base
  generation up to a smaller bound: Zimmert's twin-class bound (Zimmert
  1981, Korollar 1, evaluated in the finite form of its proof), when the
  defining polynomial is monic with integer coefficients and the class of
  the different is settled exactly.  That holds when ``O_K = Z[theta]``,
  where the different is principal, or when the bound is raised to cover
  every prime ideal above the primes dividing ``[O_K : Z[theta]]``.  The
  bound is used only when it is below ``factor_base_class_group_bound``,
  which is unchanged.  Such results report the new generation basis
  ``ClassGroupGenerationBasis::zimmert_twin_class`` (printed as
  ``zimmert_twin_class`` by ``silex-class-unit-instance``) with that bound
  as ``factor_base_generation_certification_bound()``; other results still
  report ``minkowski_type``.  For example the bound drops from 295 to 263
  for ``x^3 + x + 200``.  A ``grh`` request whose factor base reaches the
  smaller bound now records unconditional (``proven``) generation where it
  used to record ``grh`` generation, as for the quartic field of
  ``x^4 - 8 x^3 + 3 x^2 + x + 4`` (bound 98 against 122, base 100).
  Results and overall certification labels are unchanged.  See "Class groups and
  order units" in :doc:`../reference/algorithms_and_sources`.

- The index-bounded relation-saturation pass of class/unit validation now
  reports ``relation_saturation_status()`` as ``verified`` only when its
  index bound holds without GRH or verified records cover every prime
  dividing the candidate class number (in this pass coverage holds only once
  the candidate class number has been reduced to 1, so it never rescues a
  GRH-conditional bound); otherwise it reports ``unavailable``.  This
  changes the status after a ``grh`` imaginary quadratic run whose candidate
  class number exceeded the exact one, and with it
  ``SUnitGroup::source_relation_saturation_status()``, which copies it.  The
  unconditional analytic index-one route keeps ``verified``.  Labels are
  unaffected.

- A relation-saturation proof at a prime ``ell`` whose local target rank is
  zero (``ell`` divides neither the candidate class number nor the torsion
  order, and the unit rank is zero) now succeeds and records a verified
  proof with rank 0 and target 0, instead of failing.  Such a prime needs no
  check: the saturation index is prime to ``ell``.

- The exact imaginary-quadratic index route of class/unit validation now
  rolls back its factor-base generation check when it fails.  A rank-zero
  quadratic pair whose input units under-claim the torsion is then
  recomputed with the computed torsion and can be published ``proven``;
  the under-claimed units are never published ``proven``.

- Bounded saturation (``saturate_bounded``, index-bounded saturation) may now
  return different, equally valid unit generators.  Its prime selector now
  skips rational primes dividing the order discriminant, like the public
  ``select_saturation_primes``, so the two pick the same primes.  Proof
  status is unaffected.

- A ``grh`` paired class/unit transaction that is accepted by the analytic
  index-one test now records that check on the published class group.
  ``analytic_class_regulator_status()`` reads ``verified``, and
  ``analytic_class_regulator_certification()`` gives the check's own
  conditionality: ``grh`` for a Belabas-Friedman ``hR`` and ``proven`` for
  the quadratic ``L(1, chi)`` value.  With a Belabas-Friedman ``hR`` the
  evaluation's BF audit data is recorded too, so ``zeta_bf_proof_status()``
  reads ``verified`` and ``zeta_bf_proof_record()`` gives its error bound,
  cutoff, and precisions.  The record is informational.  The transaction
  keeps the class-group and unit labels at ``grh``, and that change leaves
  the public zeta and BF gates unchanged.  The one exception: a ``proven``
  record stored by a real-quadratic ``grh`` run lets a later explicit
  ``try_certify_with_units(units, proven)`` with proven units promote the
  class group to ``proven``.  The exact imaginary-quadratic and degree-one
  ``grh`` routes use no analytic check and record none.  See "Class groups
  and order units" in :doc:`../reference/algorithms_and_sources`.

- A ``grh`` paired class/unit transaction in degree three and higher or for
  a real quadratic field may now take factor-base generation from GRH: its
  factor base contains every prime ideal up to the GRH bound (the minimum of
  the Belabas-Diaz y Diaz-Friedman and Bach bounds), which generates the
  class group under GRH.  Fields whose Minkowski-type bound lies above the
  GRH bound, such as ``x^3 + x + 200``, ``x^2 - 100003`` and many quartic
  and quintic fields, used to fail closed in ``grh`` mode and are now
  published with both labels ``grh``.  GRH generation is kept apart from
  unconditional generation: ``factor_base_generation_status()`` and
  ``relation_saturation_status()`` never read ``verified`` because of it,
  a ``proven`` request is unchanged, and a later ``proven`` promotion still
  requires generation verified up to the proven generation requirement (the
  Minkowski-type bound, or Zimmert's twin-class bound where it applies).
  See "Class
  groups and order units" in :doc:`../reference/algorithms_and_sources`.

- **Breaking:** ``OrderUnitGroup::residue_dlog_kernel`` and
  ``saturate_row`` now include the torsion generator when ``ell`` divides the
  torsion order ``w``, as Hecke's ``saturate!`` does.  When ``ell | w``,
  kernel rows have ``free_rank() + 1`` columns, the last one the exponent of
  the torsion generator, and ``saturate_row`` rejects rows of width
  ``free_rank()``.  When ``ell`` does not divide ``w`` the shape is
  unchanged.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- ``try_certify_with_units`` no longer evaluates a Belabas-Friedman
  ``hR``, and its ``zeta_bf_max_cutoff`` argument is ignored.  A ``proven``
  request still needs proven units.  It uses an analytic ``hR`` only when
  that value is unconditional (degree one or the quadratic ``L(1, chi)``
  route); otherwise, and always in degree three and higher, it runs the
  relation-saturation proof at every prime dividing the candidate class
  number itself.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The public zeta gates ``try_certify_class_unit_with_zeta`` and
  ``try_certify_class_unit_with_zeta_bf`` are stricter outside degree one.
  When the analytic ``hR`` is a Belabas-Friedman value, which assumes GRH, they
  now succeed only if the units are already ``proven`` and the class group
  already has an unconditional proof component for this presentation: relation
  saturation proven at every prime dividing the candidate class number, or an
  unconditional analytic check recorded earlier (degree one or the quadratic
  ``L(1, chi)`` route, including one stored by a real-quadratic ``grh`` run),
  and then they only record the GRH-conditional check; otherwise they return
  ``false`` and leave both objects unchanged.
  ``try_certify_class_unit_with_zeta_bf`` uses a Belabas-Friedman value for
  quadratic fields too, and in degree one the exact value ``hR = 1``.  Degree
  one and the quadratic ``L(1, chi)`` route of
  ``try_certify_class_unit_with_zeta``, whose ``hR`` is unconditional, are
  unchanged.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- ``OrderUnitGroup::prove_index_bound``, ``saturate_index_bounded``,
  ``saturate_index_bounded_adaptive`` and ``prove_local_saturated`` now check
  on entry that the unit group's torsion is the torsion Silex computes for
  the order, and return ``false`` with the output unchanged otherwise.  A
  ``grh`` paired transaction makes the same check before it publishes the
  ``grh`` labels.  Every unit group built through the installed API already
  carries the computed torsion, so no result from the installed API changes.
  See :doc:`../reference/class_units_compact`.

- The unconditional quadratic ``L(1, chi)`` route now evaluates
  ``L(1, chi_D)`` with the Kronecker character ``(D/.)`` and an approximate
  functional equation adapted from FLINT ``acb_dirichlet_l_fmpq_afe``,
  instead of FLINT's Dirichlet-character ``L``-function.  It no longer builds
  a Dirichlet group, so its cost grows like ``sqrt(|D|)`` rather than
  ``|D|``.  Quadratic fields whose ``L(1, chi)`` evaluation previously failed
  (the Dirichlet group could not be initialized) or did not finish, and which
  therefore fell back to a GRH-conditional Belabas-Friedman ``hR``, now get an
  unconditional ``hR`` when ``|D| < 2^44`` (about ``1.76e13``).  The route
  takes only those fields, which bounds its running time; larger quadratic
  fields use the Belabas-Friedman fallback.  The class/unit routes that use
  it, such as ``try_certify_class_unit_with_zeta`` and the analytic
  index-one check, then record a ``proven`` analytic check where they
  previously recorded ``grh`` or failed.  The certification rule is
  unchanged: a finite, positive ``L``-value ball is unconditional, and
  anything else still falls back to Belabas-Friedman.  The incomplete gamma
  values of the evaluation are obtained by Taylor stepping over blocks of
  consecutive terms instead of being evaluated term by term, which lowers
  the cost for large ``|D|``.  At ``|D|`` near ``10^12`` and 64-bit
  precision, an evaluation now takes 10.0 s instead of 123 s for a real field
  (12.4x) and 12.5 s instead of 34.7 s for an imaginary field (2.8x), on one
  core of a Ryzen AI 9 HX PRO 370.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

Fixed
-----

- ``FiniteAbelianGroup`` invariant generators now stay stable across witness
  requests.  A signed diagonal relation matrix such as ``diag(-2, 3)``
  previously could yield generators that changed between calls.

- A class group is promoted to ``proven`` only when its verified
  factor-base generation check covers the required generation bound
  (``factor_base_generation_checked_bound()`` at least
  ``factor_base_generation_bound()``).  Previously a public
  ``check_factor_base_generation_bound`` call at a smaller bound, such as 1,
  could stand in for the full check and let a ``proven`` real-quadratic
  analytic record promote the class group through
  ``try_certify_with_units``.  Routes that run the full check themselves,
  such as the relation-saturation proof, are unchanged.

- ``try_certify_quadratic`` on a real quadratic field with class number one
  no longer sets ``unit_proof_status()`` and ``regulator_proof_status()`` to
  ``verified``.  That route proves ``h = 1`` from factor-base generation and
  the candidate class number alone and never proves a unit group, so both
  statuses now keep their previous value (``not_checked`` on a fresh
  context).  The class group is still ``proven``.

- A ``proven`` class/unit result certified by the unconditional analytic
  index-one test (degree one or the quadratic ``L(1, chi)`` value, through
  ``try_certify_with_units``, including its promotion from an
  unconditional check stored earlier by a real-quadratic ``grh`` run,
  ``try_certify_class_unit_with_zeta``, or, in degree one, where
  ``hR = 1`` exactly, ``try_certify_class_unit_with_zeta_bf``) now
  reports ``relation_saturation_status()`` as ``verified``, as documented:
  index one with verified generation means the relations are saturated at
  every prime.  It was ``not_checked`` on that route, for example on a real
  quadratic field with class number greater than one such as
  ``x^2 - 40001``, and on an imaginary quadratic field certified through
  ``try_certify_with_units``.  These routes clear the per-prime records
  unless verified records already cover every prime dividing the candidate
  class number, in which case they are kept; this includes the internal
  fallback that ``try_certify_with_units`` runs when the direct index-one
  check fails.  Labels are unchanged.  The internal Belabas-Friedman audit
  route that class/unit validation uses when ``zeta_bf_max_cutoff`` is nonzero
  follows the same rule in degree one.  Validation only takes that route for
  orders of degree greater than two, so computations through it see no change
  from this entry.

- ``try_certify_class_unit_with_zeta_bf`` (and the ``--zeta-bf-audit`` option
  of the class/unit instance tool) no longer fails on an imaginary quadratic
  field with class number greater than one that the exact
  imaginary-quadratic route proved, such as ``x^2 + 5``.  That route now
  stores a verified relation-saturation record for every prime dividing the
  exact class number, so ``relation_saturation_record_count()`` is now
  positive for such fields (one record for each distinct prime of ``h``; it
  was zero), and the GRH-conditional audit keeps ``proven`` by promoting from
  those records.  Fields with class number one are unchanged.  See "Class
  groups and order units" in :doc:`../reference/algorithms_and_sources`.

- A ``grh`` paired class/unit transaction of positive unit rank whose
  GRH-sized factor base does not reach the Minkowski-type bound (for example
  ``x^3 + x + 200`` or ``x^2 - 100003``) now fails closed once the factor
  base and its initial relations are built.  It previously kept extending
  relations, which could never lead to acceptance, and did not finish in
  practice.  Which requests succeed is unchanged.

- ``roots_of_unity`` and ``root_of_unity_order`` no longer return a
  wrong ``w``.  ``Q(zeta_9)`` previously gave ``w = 6``, where the true value
  is 18, which made ``zeta_class_regulator_product`` three times too small.
  The search now builds prime-power roots of unity up to the good-prime bound
  and returns ``w`` only when it certifies a root of unity of exactly that
  order.  In every other case it fails.  Cyclotomic fields defined by
  (translates of) cyclotomic polynomials, such as ``Q(zeta_7)``, which
  previously failed, are now supported.  See "Roots of unity" in
  :doc:`../reference/algorithms_and_sources`.

- ``roots_of_unity``, ``root_of_unity_order`` and
  ``root_of_unity_generator`` now find ``w`` in cases that previously failed
  closed.  These include ``zeta_p`` for ``p >= 5`` when the defining
  polynomial is not a cyclotomic translate, such as ``Q(zeta_5)`` given by
  ``x^4 + 3x^3 + 9x^2 + 7x + 11``.  They also include 2- and 3-power roots
  of unity in fields of degree 10 or more, and defining polynomials that are
  not monic and integral.  A primitive ``p^e``-th root of unity that the
  exact square and cube roots cannot supply is now found by Hensel lifting a
  root of ``Phi_(p^e)`` from a good prime, following Hecke.  A non-monic
  field is searched in a monic integral model.  ``w`` is still returned only
  when a root of unity of exactly the proven bound's order is certified.

- ``roots_of_unity``, ``root_of_unity_order`` and
  ``root_of_unity_generator`` fail closed less often.  The good-prime bound
  search no longer stops while the Euler phi of its gcd does not divide the
  field degree, as in Hecke, since such a gcd cannot be ``w``.  For example,
  a presentation of ``Q(zeta_40)`` whose bound previously stalled at 120 now
  gives ``w = 40``.  ``w`` is still returned only when a root of unity of
  exactly the proven bound's order is certified.

- ``Element::is_power`` and ``Element::is_square`` now find non-integral roots
  that were previously reported as unsupported, for example
  ``((3 + 4i) / 5)^3`` in ``Q(i)``.  A non-integral input is rescaled by its
  power-basis denominator before the root is lifted, and the root is scaled
  back and verified exactly.  Answers for algebraic integers are unchanged.
  See "Element powers and roots" in :doc:`../reference/algorithms_and_sources`.
- ``OrderUnitGroup::class_regulator_index_bound`` now returns ``false`` when
  the upper endpoint of ``h_cand R_cand / hR`` is below one, instead of
  reporting an index of one.  No positive integer index is consistent with
  such an analytic value.  See :doc:`../reference/class_units_compact`.
- Unit saturation proofs now enlarge a unit that is an ``ell``-th power only
  up to torsion, such as ``-(3 + 2 sqrt2) = -(1 + sqrt2)^2`` at ``ell = 2``
  in ``Q(sqrt2)``.  Such cases previously failed closed.  Candidate roots now
  include the torsion generator, and verification still requires an empty
  torsion-extended kernel.
- ``saturate_local_once``, ``saturate_bounded``, ``saturate_index_bounded``,
  ``saturate_index_bounded_adaptive``, and
  ``set_relation_kernel_units_index_bounded_saturated`` now also find a unit
  that is an ``ell``-th power only up to torsion, such as ``-(3 + 2 sqrt2)``
  at ``ell = 2``, when ``ell`` divides the torsion order ``w``.  They test
  the torsion-extended candidates described under "Changed".  A ``proven``
  unit computation from a class-group context passes through
  ``set_relation_kernel_units_index_bounded_saturated`` before
  ``prove_index_bound``, so that step now considers torsion too and can hand
  the proof a different generator; ``prove_index_bound`` itself is
  unchanged, and ``ell`` is still verified only when the torsion-extended
  kernel is empty.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.
- ``zeta_residue_bf_audit`` and ``zeta_class_regulator_product_bf_audit`` now
  leave every output unchanged when they fail.  Every BF audit API now writes
  the value last, so a call that passes the same ``Arb`` for the value and the
  error bound receives the value.  A
  ``max_cutoff`` near ``UWORD_MAX`` is now accepted, where it previously
  wrapped and made every call fail.  See :doc:`../reference/factored_zeta`.
