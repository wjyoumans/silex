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

- ``compute_sunit_groups`` accepts ``grh`` class-group and ordinary-unit
  inputs.  The new option ``SUnitComputeOptions::requested_certification``
  (default ``proven``) states the accepted label: a ``proven`` request fails
  at ``SUnitComputeStage::input_validation`` for any ``grh`` input, a ``grh``
  request accepts ``grh`` or ``proven`` inputs, and any other request fails at
  ``input_validation``.  ``SClassGroup::certification_status()`` is the class
  label; ``SUnitGroup::certification_status()`` is the weaker of the class and
  unit labels; ``proven`` inputs give ``proven`` outputs.  ``proof_status()``
  and ``regulator_proof_status()`` are ``verified`` for every published group
  and mean exact verification relative to the inputs, not an unconditional
  result; the ``grh`` label is in ``certification_status()``.  With a ``grh``
  class group, ``SUnitGroup::preimage`` returns ``unknown`` (``success``
  false, stage ``valuation_solve``) where a ``proven`` class group gives
  ``not_sunit``; a residual that is not a unit stays ``not_sunit``, and
  ``grh`` ordinary units are accepted in the residual step.
  ``silex-class-unit-instance --compute-sunit`` now follows ``--mode``, so
  ``grh`` rows publish ``grh`` S-groups.  Default requests are unchanged.  See
  :doc:`../reference/sunit_groups`.

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
- ``silex-class-unit-instance`` adds ``generation_different_route`` to its
  ``class_group`` JSON object.  It names how the class of the different is
  handled for the ``zimmert_twin_class`` generation basis,
  ``kor1_generation`` (the spelling of the detail log: no hypothesis on the
  different is needed), and is null for every other basis.  Additive; the
  format is unchanged.
- ``silex-class-unit-instance`` raises its default ``max_candidates`` for
  degree 3 and up from 5000 to ``max(5000, n^2 + 64)``, so degrees above 70
  are no longer stopped before the small-LLL phase finishes its ``n^2``
  candidates.  The slack of 64 is headroom after that phase for the later
  phases; it was chosen, not measured.  ``--max-candidates`` still overrides;
  smaller degrees are unchanged.
- ``silex-class-unit-instance`` reports ``failure_reason``
  ``factor_base_honesty_unwitnessed`` (``failure_stage`` ``total``) when a
  class/unit computation fails because the last candidate attempt stopped at
  a required prime that the factor-base honesty witness search could not
  witness at its final caps, with no factor-base restart left.  It reports
  ``factor_base_size_limit_exceeded`` when the factor base would exceed the
  size limit described below.  Other computation failures still report
  ``class_unit_computation_failed``, and
  results, labels and successful runs are unchanged.  A new top-level JSON
  field ``failure_detail`` holds structured detail about a failure, with keys
  that depend on ``failure_reason``.  For ``class_unit_computation_failed``
  it holds ``cause``, free text naming the last unit-group validation
  outcome reached (null when validation never ran).  For
  ``relation_search_budget_exhausted`` it holds ``cap`` (below).  It is null
  except for those reasons and ``factor_base_honesty_unwitnessed``, where it
  names the unwitnessed prime (``p``, a decimal string, and
  ``residue_degree``), the witness search route (``search``, ``lattice`` or
  ``t2``), the stage of its last search
  (``stage`` and ``max_stage``), and that stage's caps: ``radius``,
  ``twists`` and ``random_tries`` on the lattice route, or
  ``random_tries``, ``factor_attempts`` and ``element_steps`` on the T2
  route.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.
- ``silex-class-unit-instance`` reports the new ``failure_reason``
  ``relation_search_budget_exhausted`` (``failure_stage`` ``total``) when a
  ``grh`` class/unit computation fails because a search cap stopped its LLL
  relation search: the candidate cap (``--max-candidates``, option
  ``ClassGroupComputeOptions::max_candidates``) was spent, or the relation cap
  (``--max-relations``, option ``max_relations``) was reached, before the
  search found the relations it needed.  The tool's default caps for
  degree 3 and higher count too, so a ``grh`` run that exhausts them with no
  flag set also reports the new reason.  Such runs previously reported
  ``class_unit_computation_failed``.  Its ``failure_detail`` holds ``cap``,
  the option name of the cap that stopped the search (``max_candidates`` or
  ``max_relations``), preceded by the same ``cause`` key as for
  ``class_unit_computation_failed`` when a unit-group validation ran.
  ``proven`` runs report the same reason and detail when they fail after a
  cap was reached by one of their relation searches (even one that reached
  its goal) or after the relation cap alone refused a larger relation-kernel
  target; the detail keeps the final ``cause`` next to the cap, so a later
  failure with another cause stays visible.  For example
  ``x^3 - 6x^2 - 9x - 1`` with ``--max-relations 12`` now reports
  ``relation_search_budget_exhausted`` with ``cause`` and ``cap``
  ``max_relations``, and with ``--max-relations 1`` only ``cap``.  The tool's
  default caps count for ``proven`` runs as well.  Only failing runs
  change; successful runs, results and labels are unchanged.  See "Class
  groups and order units" in :doc:`../reference/algorithms_and_sources`.

- A ``proven`` request for a maximal imaginary quadratic order whose base up
  to the generation bound would exceed the factor-base size limit of 2^14
  generators no longer fails with ``factor_base_size_limit_exceeded`` (for
  example ``x^2 + 1000000000039``, ``h = 1113261``).  It keeps its compact
  base and is certified by a structure certificate, checked in exact reduced
  form arithmetic against the exact class number; the result is ``proven``
  with ``ClassGroupContext::factor_base_generation_basis()`` equal to the new
  value ``ClassGroupGenerationBasis::structure_certificate`` (printed as
  ``structure_certificate`` by ``silex-class-unit-instance``), and no
  generation bound is published
  (``factor_base_generation_certification_bound`` returns false).  On this
  route the published relations are the staged relation rows with exact
  witnesses.  When the certificate fails, or the relation index is below
  ``h``, the run fails closed with the generic reason
  ``class_unit_computation_failed``, with no fallback.  Results on other
  routes are unchanged; the refusal reason
  remains for other fields.  See
  :doc:`../reference/algorithms_and_sources`.

Changed
-------

- When the relation search factors an integral relation generator ``g`` through
  its principal ideal ``g O``, it computes the HNF of ``g O`` modulo the norm
  ``|N(g)|`` it already knows, with FLINT's ``fmpz_mat_hnf_modular``, instead of a full
  HNF of the multiplication matrix.  The result is the same canonical HNF,
  so published outputs are unchanged on the tested ``grh`` and ``proven``
  corpus rows.  See "Number fields, orders, and ideals" in
  :doc:`../reference/algorithms_and_sources`.

- The compact-representation reconstruction of ``grh`` units, when an
  exponent has more than 64 bits, prefers CRT primes at which the defining
  polynomial splits into linear factors and reduces the exponents of the
  modular power products modulo the unit group exponent ``p - 1`` (``lcm(p^d - 1)`` at squarefree
  primes that do not split), making
  negative exponents positive where that is cheaper.  The primes differ from
  PARI's ``init_modular_big`` sequence only for such rows; every reconstructed
  unit and published output is unchanged.  See
  :doc:`../reference/algorithms_and_sources`.

- ``Lat::lll_reduce``, the rounded-embedding LLL of the LLL relation search,
  and the unit-candidate and relation-saturation kernel bases of the
  class-group computation share one certified first LLL phase.  Its exact
  reducedness certificate is an integer Gram-Schmidt test when the reduced
  entries have at most 16 bits and ``fmpz_mat_is_reduced`` otherwise; both
  decide the same predicate for FLINT's parameter range ``eta^2 < delta`` (all
  callers use eta 0.51 with delta 0.99 or 0.5), and the integer test measured faster at those
  sizes.  ``Lat::lll_reduce``, the rounded-embedding LLL and the saturation
  kernel previously always used ``fmpz_mat_is_reduced``.  Every route still
  returns the basis and transform ``fmpz_lll`` returns, so every published
  result is unchanged.  See ``Lat::lll_reduce`` in
  :doc:`../reference/algorithms_and_sources`.
- The two integer LLL calls of the grh two-LLL ideal-lattice route (the
  ``Z_BASIS`` pre-reduction and the ``GRAM`` reduction) take the same
  certified first LLL phase.  The ``GRAM`` certificate is the integer
  Gram-Schmidt test up to rows times bits of 3000 and FLINT's
  ``fmpz_lll_is_reduced`` above that; a failed certificate restores the input
  and calls ``fmpz_lll``, so every published result is unchanged.  See
  :doc:`../reference/algorithms_and_sources`.
- The integer LLL of the finite-ideal T2 enumeration setup (the reduction of
  a prime ideal's rounded embedding rows, run for each ideal of the
  factor-base honesty check, and for each ideal of the finite-ideal relation
  search in both ``grh`` and ``proven`` mode) takes the same certified first
  LLL phase, which skips FLINT's multiprecision check of the phase output.
  A failed certificate restores the input and calls ``fmpz_lll``, and rows
  with entries above 250 bits call ``fmpz_lll`` directly, so the reduced bases and
  every published result are unchanged.  On a busy machine the ``proven``
  class/unit computation of the field defined by ``x^12 + 7`` took about 16%
  less CPU time (median of seven alternating runs); no controlled measurement
  has been made.  See :doc:`../reference/algorithms_and_sources`.

- Publishing a class-group presentation no longer recomputes the Smith form
  (``fmpz_mat_snf_transform``) of the HNF basis of the relation lattice when
  that basis is equal to the one of the previous publish, as it is when new
  relations leave the lattice unchanged; the Smith form and right transform
  of that basis are reused.  Published outputs are unchanged on the tested
  ``grh`` and ``proven`` corpus rows.  On a busy machine, the change in the
  ``grh`` class/unit computation of the field defined by ``x^32 + 1`` was
  small and within run-to-run noise; no controlled measurement has been made.

- The ``grh`` and ``proven`` class saturation with units no longer repeats the
  power test of a kernel row whose exact integer vector it has already
  settled in the same round of the prime scan: within a round no relation,
  unit or factor-base entry is appended, so the outcome of such a row cannot
  change at a later prime.  Published outputs are unchanged on the tested
  ``grh`` and ``proven`` corpus rows.  The scan runs on many ``grh`` corpus
  rows at the default settings (for example cyclotomic, pure and quadratic
  fields), and the memo skips rows on most of them; no controlled
  measurement has been made.  See :doc:`../reference/algorithms_and_sources`.

- The rounded-embedding LLL of the LLL relation search no longer pays
  FLINT's multiprecision reducedness pre-check on every call.  For at least
  two rows and entries of at most 250 bits it runs FLINT's double-precision
  first phase (``fmpz_lll_d_with_removal_knapsack``) and accepts it only when
  an exact certificate shows it LLL-reduced, as ``Lat::lll_reduce`` does;
  otherwise it restores the input and calls ``fmpz_lll``.  The reduced
  basis and transform are those ``fmpz_lll`` returns, so every ideal's basis
  and every published result is unchanged.  See "Class groups and order
  units" in :doc:`../reference/algorithms_and_sources`.

- The one-large-prime step of relation factoring builds only the prime
  ideals above the large rational prime that can contain the element, found
  by a Dedekind--Kummer divisibility test, instead of every prime above it,
  when the prime does not divide the index of the equation order.  The
  primes it builds are those built before, in the same order, so published
  outputs are unchanged.  See "Local algebra and relations" in
  :doc:`../reference/algorithms_and_sources`.

- The ``proven`` analytic finish no longer computes the HNF transform of
  the relation matrix at every check.  It checks the relation rank with
  ``fmpz_mat_rank`` and computes the transform only when the LLL kernel
  basis fails its exact checks.  The kernel basis and every published
  output are unchanged on the tested ``proven`` corpus rows.  On a busy
  machine, the 17 analytic finish checks of the field defined by
  ``x^2 - 10000000005`` (``proven``) took 1.6 to 3.1 s instead of 71 to 143 s over three runs
  each.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The relation search residue screen now also runs on maximal orders
  that are not the equation order, for primes that do not divide the
  equation-order index: it tests the theta-numerator of an element modulo
  the prime ideal's polynomial instead of falling back to a principal-ideal
  factorization.  Published results are unchanged: 185 corpus rows in
  ``grh`` and 185 in ``proven`` mode gave identical published output before
  and after (``grh``: 184 finished, 1 timeout on both; ``proven``: 168
  finished, 17 timeouts on both), and the screen agreed with the
  principal-ideal fallback on every decision checked.  See "Class groups and
  order units" in :doc:`../reference/algorithms_and_sources`.

- The ``proven`` analytic finish now obtains a small kernel basis for
  the unit exponent rows by LLL reduction of ``[c*A | I]`` with exact
  checks and a fallback to the HNF-transform rows, instead of using the
  HNF transform rows directly.  Published outputs are unchanged.  On a
  busy machine, ``real_quadratic_d1000000000009`` now finishes ``proven``
  within a 300 s limit where it timed out before.  See "Class groups and
  order units" in :doc:`../reference/algorithms_and_sources`.

- ``FactoredElement::is_power_si`` (the reduced test, ``n >= 3``) in degree
  10 and higher first runs the residue disproof at the first admissible
  prime on the factors of the leftover, before forming the leftover's
  product.  The disproof answers
  "not a power" only where the residue disproof on the expanded leftover
  would.  Otherwise the leftover is expanded and tested as before, so every
  answer is unchanged, with one accepted exception: a non-power that the
  disproof decides now returns "not a power" where expanding the leftover
  would have run out of resources and returned a failure.  ``grh`` saturation on ``Q(zeta_52)`` (degree 24)
  no longer expands most of its non-cube leftovers.  See "Factored elements and
  analytic products" in
  :doc:`../reference/algorithms_and_sources`.

- The ``proven`` relation-saturation proof now takes the relation-kernel
  basis that its ``beta`` rows are size-reduced against from the same
  ``[c*A | I]`` reduction, with the same exact checks and an
  ``fmpz_mat_is_reduced`` certificate, when the kernel has entries above
  250 bits, instead of running ``fmpz_lll`` on that kernel.  The
  ``fmpz_lll`` basis remains the fallback.
  ``real_quadratic_d1000000000009`` ``proven`` takes about 14 s instead of
  about 90 s (release build, three runs each on a busy machine; the
  controlled comparison is pending), with identical published outputs and
  saturation proof records; the proven corpus rows, whose kernels are
  smaller, take the route they took before.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The exact class number of an imaginary quadratic field with
  ``|D| >= 2^44`` (``proven`` class/unit computations) now first evaluates
  ``L(1, chi_D)`` to an absolute tolerance sized for the class number,
  instead of FLINT's precision mapping at the bit length of ``|D|`` plus 32
  bits.  The enclosure is still rigorous and the class number is still taken
  only when its ball holds exactly one positive integer.  When that first
  step does not isolate it, the previous precisions are tried as before, so
  the class number and every published output are unchanged and no field
  fails that succeeded before.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The ``proven`` imaginary quadratic witness step builds each factor-base
  row's ideal Hermite normal form by Hensel lifting the prime roots and
  combining them with the Chinese remainder theorem, instead of multiplying
  ideals.  Published results are unchanged.

- The stable auxiliary-prime scan of the ``proven`` unit saturation proof
  (proof rank two or more) builds each ``ell``-th-power character column by
  evaluating the units at the degree-one roots of the defining polynomial
  modulo ``q``, as the proof-selector scan already does, instead of
  decomposing ``q`` into prime ideals and reducing into each residue field.
  The roots are taken in the order ``decompose_prime`` lists its primes, and a
  root whose direct evaluation declines falls back to the prime ideal over
  it.  The two columns differ only by a nonzero scalar modulo ``ell``, so the
  proof kernel is unchanged.  The direct evaluation also gains a word-size
  route for word-size ``q``.  Published results are unchanged on the 70
  ``proven`` corpus rows.

- The word-size route of that column evaluation reads each unit factor's
  stored numerator polynomial and denominator in place, reducing them modulo
  ``q`` without copying the factor to a rational polynomial or inverting
  each coefficient's denominator, and inverts one accumulated denominator per
  unit.  A factor whose denominator is divisible by ``q`` takes the previous
  per-coefficient route, so the route fails where it failed before and the
  columns are unchanged.  Published results are unchanged on the 155
  ``proven`` corpus rows that finish within 60 seconds.

- That stable scan now builds the columns over each ``q`` one at a time and
  stops at the column that ends the scan (full rank, the stability threshold,
  or a failure), instead of building every column over ``q`` first.  The
  columns it no longer builds were discarded unused: whether ``q`` is skipped
  is decided before its first column, and an unavailable column is left out
  without affecting the others, so the proof kernel and the proof records'
  ``local_primes`` are unchanged.  On ``pure_x12_plus_7`` this avoids 5596 of
  14892 column evaluations.  Published results are unchanged on the 155
  ``proven`` corpus rows that finish within 60 seconds.

- The rational reconstruction of the reduced regulator in the ``proven``
  analytic finish stores each continued-fraction convergent of a unit
  coordinate as it is, instead of reducing it with a gcd.  Convergents and
  semiconvergents are already in lowest terms with a positive denominator, so
  the reconstructed rationals, the regulator and every published output are
  unchanged.  At high working precision (up to 40960 bits) that gcd was most
  of the finish time on large real quadratic fields.

- Unit extraction's memo of factored logarithmic embeddings finds each stored
  factor through a hash of its coefficients and a hash table, instead of
  copying the factor to a rational polynomial and comparing against every
  stored factor.  Every hash match is still confirmed by element equality, so
  the memo returns the same values and published results are unchanged on
  the 184 ``grh`` corpus rows that finish within 120 seconds.

- Each factor stored in a ``FactoredElement`` keeps the hash of its
  coefficients and an origin shared by the copies the factored operations
  make of it.  The factored logarithmic embedding and unit extraction's log
  memo reuse the stored hash instead of hashing every factor again, and two
  factors of the same origin are known to be equal without comparing their
  coefficients.  Equal factors pushed separately are still compared, so the
  merged factors, the log embeddings and every published result are
  unchanged.

- Unit extraction's memo of factored logarithmic embeddings remembers the
  result of each comparison between two stored factors of different origins,
  true or false, and answers a later comparison of the same two origins from
  it.  This covers matching a candidate against the memo's entries and, when
  the memo misses, merging the candidate's equal factors.  An origin is never
  reused and a stored factor never changes, so each answer is the comparison
  it replaces, and the log embeddings and every published result are
  unchanged.  On ``cyclotomic_real_31_deg15`` in ``grh`` mode, 135,174 comparisons
  went through the memo; a further 69,643 come from callers that use
  ``FactoredElement::logarithmic_embedding`` directly, and those are not
  cached.  The row took about 8% less CPU time (busy-machine comparison).

- In fields of degree 10 and higher, the residue test that can prove an
  element is not a square or an ``n``-th power keeps its answer per defining
  polynomial, element and exponent, so asking it again does not scan the
  primes again.  Class-group saturation asks it of the same few elements many
  times; on ``Q(zeta_19)`` one of them is a true square that no residue can
  disprove, so every repeat used to scan all 40 primes.  The test reads
  nothing else, so its answers are unchanged: such a square is still reported
  as undecided, and published results are unchanged on the ``grh`` and
  ``proven`` corpus rows that finish within 120 seconds.

- Unit extraction from class-group relations reduces its kernel basis
  without FLINT's high-precision reducedness certificate.  ``fmpz_lll``
  first runs a double-precision pass and then certifies the result with an
  ``mpfr`` check before returning it; on these small kernels that check cost
  far more than the reduction.  Silex now runs the same first pass and
  accepts its output when an exact integer Gram-Schmidt test shows it meets
  FLINT's reducedness definition, and otherwise calls ``fmpz_lll`` as
  before.  The returned basis is the one ``fmpz_lll`` returns, so units and
  every published result are unchanged on the ``grh`` and ``proven`` corpus
  rows that finish within 120 seconds.

- The Belabas-Friedman ``hR`` estimate counts the roots of the defining
  polynomial modulo the large cutoff primes with its word-size kernel for
  defining polynomials of degree up to 16, instead of up to 6.  Degrees 7 to
  16 previously took the per-prime FLINT route.  The kernel covers primes
  below ``2^30`` up to degree 8 and below ``2^29`` for degrees 9 to 16,
  where its lazily reduced sums stay below ``2^64``; other primes keep the
  per-prime route.  The counts are the same, so the Belabas-Friedman
  values and published results are unchanged on the 184 grh corpus rows.

- In ``OrderUnitGroup::prove_index_bound`` (the ``proven`` unit proof),
  when the scan of the auxiliary primes ``q <= aux_bound`` at a prime ``ell``
  leaves a nonempty kernel that it has not certified, the first scan of the
  unbounded stable auxiliary-prime loop now runs before the ``ell``-th-root
  attempts on that kernel's rows.  If that scan's kernel is empty, ``ell`` is
  verified at once with the record the loop would have published, since no
  row can then have a root in the order; otherwise the root attempts and the
  loop run as before.  This skips root attempts that cannot succeed (on
  ``random_d9_h2_s8`` all 30 such primes ``ell``, on ``pure_x12_plus_7`` all
  15).  Published results, including the proof records' ``aux_prime_bound``
  and ``local_primes``, are unchanged on the 70 ``proven`` corpus rows.  See
  "Class groups and order units" in :doc:`../reference/algorithms_and_sources`.

- ``quadratic_fundamental_unit`` and the ``proven`` unit route now use the
  continued-fraction fundamental unit for every monic integral ``x^2+bx+c``
  with positive nonsquare discriminant, normalised as the existing quadratic
  route publishes it.  The published unit may differ from earlier releases on
  such polynomials.  On those rows ``regulator_radius``,
  ``unit_proof_record_count``, ``relation_count``,
  ``class_group.relation_sources.search``, ``relation_kernel_unit_count``,
  ``generator_count`` and ``relation_rank`` change, and the documented
  contract of ``quadratic_fundamental_unit`` widens from explicit quadratic
  fields and ``x^2-d`` to monic integral ``x^2+bx+c``.  Per-row cost moves
  both ways: ``real_quadratic_d100000005`` ``proven`` gets faster, while
  ``real_quadratic_d100000001`` ``proven`` gets slower.

- The ``proven`` factor-base honesty check accepts a principal witness that
  also carries a prime ideal ``Q`` outside the factor base, when ``Q`` lies
  above a rational prime ``q`` other than the required prime's, ``Q`` is
  unramified, and every other prime ideal above ``q`` is in the factor base.
  The relation ``qO_K = Q_1^e_1 ... Q_g^e_g`` puts the class of ``Q`` in the
  subgroup of the factor base, the step the proof-target selection already
  uses (PARI 2.17.3 ``buch2.c:be_honest`` and ``bnftestprimes``).  Accepting
  such a witness is a Silex extension; PARI's ``buch2.c:can_factor`` does
  not accept it.  Results are unchanged on the 65 ``proven`` benchmark rows
  that passed before (measured, not by construction); a computation that
  failed with ``factor_base_honesty_unwitnessed`` or timed out searching for
  a witness may now succeed.  For example, the ``proven`` computation for the imaginary
  quadratic field of discriminant ``-100000000000319`` now finishes, with
  class group invariants ``[2, 2, 3953768]``.  See
  :doc:`../reference/algorithms_and_sources`.

- When a unit group keeps units with exponents beyond ``slong`` and rebuilds
  them by bounded CRT reconstruction, and when compact unit generators are
  expanded to reduce relation units modulo the unit group, the CRT
  coordinate bound now evaluates the quantity PARI's ``log2Mbound``
  estimates, the norm of the inverse embedding matrix, with a certified Arb
  inverse, and uses it whenever it is certified.  PARI's expression is used
  only when the inverse cannot be certified at the working precision; it
  still decides success, the precision ladder and failure.  On every
  signature the new bound is never larger than before, so the
  reconstruction never uses more CRT primes.  For the ``grh``
  computation on the totally real degree-32 subfield of the 136th
  cyclotomic field the bound falls from 454 to 30 bits and the
  reconstruction uses 1 CRT prime instead of 8.  Published results are
  unchanged on the test corpus.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The relation-production prepass drops its outstanding relation surplus
  after each small-norm pass, as PARI's ``Buchall_param`` does
  (``cache.end = cache.last``), so a later request asks only for the missing
  rank.  It previously kept asking for the surplus and stopped without
  collecting more relations.  The ``proven`` computation for
  ``x^5 - 2x^4 + x^3 + x^2 - 3x + 1`` (discriminant ``-11119``) now publishes
  the fundamental units ``2 - x + x^3 - x^4``,
  ``3 - x - x^2 + 2x^3 - x^4`` and ``-2 + 2x - 2x^3 + x^4`` instead of
  ``1 - x``, ``-3 + x + x^2 - 2x^3 + x^4`` and ``1 - x + x^2``.  The class
  group, regulator and ``proven`` labels are unchanged.

- In prime degree, ``OrderUnitGroup::regulator_index_bound`` on an order
  known to be maximal now also uses the discriminant-dependent regulator lower
  bound of Pohst-Zassenhaus (6.22), taking the larger of it and the
  signature bounds.  The unit-index scan of a ``proven`` unit computation
  then checks fewer primes: for the quintic of discriminant ``-401370255`` the
  index bound falls from 1184 to 236, for the cubic of discriminant
  ``-1080004`` from 1063 to 84.  Published groups, regulators and labels are
  unchanged on the test corpus; ``OrderUnitGroup::regulator_index_bound`` and
  the unit proof records (``unit_proof_record_count`` and
  ``unit_proof_record``) change: the records cover fewer primes.  For direct
  callers of ``OrderUnitGroup::set_relation_kernel_units_index_bounded`` and
  ``OrderUnitGroup::saturate_index_bounded``, which label a group ``proven``
  when the index bound reaches 1, the smaller bound can make that happen more
  often.  Composite degrees, non-maximal orders and
  ``unit_lower_regulator_bound`` are unchanged.  See
  "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- ``grh`` class/unit computations on an imaginary quadratic field with
  ``|D| >= 2^44`` no longer list the reduced forms of discriminant ``D``,
  whose time and memory grow like ``sqrt(|D|)``.  They take ``h`` from the
  GRH relation lattice and the analytic index-one test, as in every other
  degree; for ``x^2 + 10^14 + 319`` this took 0.12 s and 22 MB instead of
  3.2 s and 1.4 GB on a busy machine.  ``h``, the invariants and the labels
  are unchanged.  For these fields completeness is now ``grh`` with basis
  ``belabas_friedman`` (previously ``proven``, ``exact_class_number``), and
  ``analytic_class_regulator_status()`` is ``verified`` with certification
  ``grh`` (previously ``not_checked``).  Fields with ``|D| < 2^44`` keep the
  exact class number and their output is unchanged.  See "Class groups and
  order units" in :doc:`../reference/algorithms_and_sources` and
  :doc:`../development/class_number_route_cutoff`.

- ``grh`` class/unit computations on an imaginary quadratic field with
  ``2^34 <= |D| < 2^44`` now also take ``h`` from the GRH relation lattice
  instead of the exact class number (the cutoff of the previous entry moves
  from ``2^44`` to ``2^34``).  ``grh`` is faster for these fields: on an idle
  machine the sum over sampled fields is about 0.83 times the exact route's
  time in the bin from ``n = 2^34``, and about 0.26 times near ``2^44`` (the
  sampled fields are ``x^2 + n``).  ``h``, the invariants and the invariant
  generator matrix are unchanged (compared on ``D = -(2^34 - 1)``,
  ``D = -(2^34 + 3)`` and ``D = -(2^40 + 3)``, which agree between the two
  routes), but the completeness record of these fields changes: completeness
  is now ``grh`` with basis ``belabas_friedman`` (previously ``proven``,
  ``exact_class_number``), ``analytic_class_regulator_status()`` is
  ``verified`` (previously ``not_checked``),
  ``analytic_class_regulator_certification()`` is ``grh``, and
  ``zeta_bf_proof_record()`` is now populated (previously none).  Fields with ``|D| < 2^34`` keep the
  exact class number and their output is unchanged; ``proven`` requests are
  unaffected (their cutoff stays ``2^44``).  See
  :doc:`../development/class_number_route_cutoff`.

- ``proven`` class/unit computations on an imaginary quadratic field with
  ``|D| >= 2^44`` take the exact class number ``h(D)`` from Dirichlet's
  class number formula with a rigorous ``L(1, chi_D)`` enclosure in ball
  arithmetic, instead of listing the reduced forms of discriminant ``D``
  (time and memory growing like ``sqrt(|D|)``).  The value is taken only
  when the ball contains exactly one positive integer, and otherwise the
  computation fails closed.  At ``|D| = 1e14 + 31`` a direct call of the
  function gives ``h = 14414435``, PARI's value, in about 100 s and 14 MB.
  No full proven run reaches this step today: ``x^2 + n`` fields above about
  ``4e11`` fail first at the factor-base size limit (next entry), and the
  ``x^2 + x + c`` presentation at ``1e14`` did not finish.  The formula is
  checked directly by the slow test ``t-zeta-quadratic-l1`` (run with
  ``SILEX_TEST_SLOW``).  Fields with ``|D| < 2^44`` and all labels are
  unchanged.

- A class-group factor base with more than ``2^14`` generators is refused
  before any matrix is allocated: ``ClassGroupContext::build_factor_base``,
  ``set_factor_base`` and ``set_relation_matrix`` return false, and a
  class/unit transaction fails with ``failure_reason`` ``factor_base_size_limit_exceeded`` in
  ``silex-class-unit-instance``.  The limit is a fixed Silex choice, so the
  same fields fail the same way on every machine.  Previously the dense
  relation matrices (2 GiB each at ``2^14``) were allocated, and FLINT
  aborted the process when the allocation failed.  For ``proven`` requests
  this affects imaginary quadratic fields with ``|D|`` above about
  ``4e11``, whose factor base reaches ``sqrt(|D| / 3)``: for
  ``x^2 + 10^14 + 31`` the run now fails with that reason in about 3 s
  instead of aborting.  Largest factor-base size seen in the test suite is
  380 and on the succeeding benchmark corpus rows 34, so no existing
  fixture changes.

- The ``proven`` factor-base honesty check of class-group relation search
  looks for each required prime's witness with the capped T2 enumeration of
  PARI's ``be_honest`` (at most 500 factored candidates per ideal and 50
  random sub-factor-base twists, ``maxtry_FACT`` and ``maxtry_HONEST`` in
  ``buch2.c``) instead of the uncapped scalar/lattice search.  A prime with
  no witness within the caps now ends in a factor-base restart or, without
  one, doubles the caps at up to three further stages (to 4000 candidates
  and 400 twists) and then fails closed with ``failure_reason``
  ``factor_base_honesty_unwitnessed`` (``failure_detail`` ``search``
  ``t2``) instead of searching until the caller's time limit.  Every
  witness passes the same exact check as before, and results and labels of
  the benchmark corpus's 65 succeeding ``proven`` rows are unchanged.  Seven
  benchmark rows now finish ``proven``, including three fields that ran past
  a 600 s limit, ``Q(zeta_31)^+`` and
  ``Q(sqrt(-3), sqrt(-5), sqrt(-7))`` and ``Q(sqrt(5), sqrt(-7), sqrt(-11))``, with the class groups and regulators of PARI's
  certified ``bnfinit``.

- The ``proven`` factor-base honesty check on the maximal order of a
  quadratic field reads each rational prime's splitting from the Kronecker
  symbol ``(D/p)`` of the field discriminant and skips an inert ``p``
  without computing its prime decomposition; an inert ``p`` has no proof
  target (its only prime ideal is ``pO_K``), so results are unchanged.  Split
  and ramified primes are decomposed as before.  Published outputs, labels
  and counters are identical on the 65 ``proven`` and the 183 completed
  ``grh`` benchmark rows.

- The class-relation saturation with units keeps a stack of ``ell``-th-power
  residue characters of the relation generators and units at degree-one
  primes, grown one prime at a time for each candidate that survives the
  divisibility check and the current stack, and skips the candidates outside
  the common kernel, which cannot be ``ell``-th powers.  The stack starts
  only after the round has run as many power tests as there are inputs, a
  Silex heuristic that Hecke does not have.  Stacking stops by
  Hecke's ``compute_candidates_for_saturate`` rule with ``stable`` starting
  at 1.5 and doubling after a failed power test on a candidate that the
  stopped stack let through.  The appended relations,
  class groups, generators and units are unchanged.  See "Class groups and
  order units" in :doc:`../reference/algorithms_and_sources`.

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

- The Belabas-Diaz y Diaz-Friedman bound of a ``grh`` class/unit request is
  now a cutoff at which an interval comparison proves the criterion, the
  smallest one unless an undecidable near-tie leaves a larger one
  (doubling from 2, then an exact bisection, as PARI 2.17.3
  ``Buchall_param``), instead of the result of a search that doubled from
  100 and stopped at a step below 100 and so never returned less than about
  50.  The reported GRH generation bound
  (``factor_base_generation_certification_bound()`` with basis ``bdf``) and
  the factor base built from it are smaller: for ``x^3 + x + 200`` the bound
  is 19 instead of 50 and the base is built to 76 instead of 100.  The base
  is still doubled until it has more than 20 prime ideals.  A ``grh`` base
  can now stay below the Minkowski-type bound where it used to reach it, so
  factor-base generation of such a result is recorded as ``grh`` rather than
  ``proven`` in ``factor_base_generation_certification()``, and a later zeta
  promotion of it fails closed, as documented for any base that covers only the GRH bound.  Labels and failure rules are
  unchanged.  See "Class groups and order units" in
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
  1981, Korollar 1, with the series F enclosed by the certified tail bound
  GB-zimmert-F-tail-bound of silex-research, commit 06191ff), for every
  defining polynomial.  The bound generates the class group with no
  hypothesis on the class of the different, by the theorem
  GB-kor1-generation of silex-research (commit 8926120), so non-monogenic
  fields and non-monic defining polynomials get it too.  The
  bound is used only when it is below ``factor_base_class_group_bound``,
  which is unchanged.  Such results report the new generation basis
  ``ClassGroupGenerationBasis::zimmert_twin_class`` (printed as
  ``zimmert_twin_class`` by ``silex-class-unit-instance``) with that bound
  as ``factor_base_generation_certification_bound()``; other results still
  report ``minkowski_type``.  For example the bound drops from 295 to 263
  for ``x^3 + x + 200``.  Its constants come from a precomputed table of
  exact upper bounds, generated offline from an Arb enclosure and checked
  against it by the tests, so a run does not evaluate them.  A ``grh`` request whose factor base reaches the
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
  GRH-conditional bound); otherwise it reports ``unavailable``.  The exact
  imaginary-quadratic index bound holds without GRH when factor-base
  generation has been verified up to the required bound, under a ``grh``
  request as under a ``proven`` one.  This changes the status after a
  ``grh`` imaginary quadratic run whose candidate class number exceeded the
  exact one and whose generation was not checked to that bound, and with it
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

- The analytic index-one test of a ``grh`` paired class/unit transaction in
  degree three and higher now first evaluates ``hR`` to an absolute error
  below 1/4 on ``log hR`` with the Belabas-Friedman Theorem 7 error bound,
  instead of the full default zeta route.  When that value shows index one it
  is used, and the recorded ``zeta_bf_proof_record()`` values (cutoff and
  error bound) are those of the coarse evaluation.  Otherwise (no cutoff up
  to 20007, or the coarse value does not show index one) the coarse value is
  discarded and the default zeta route is used, so such fields pay up to
  part of one extra analytic evaluation.  Accept/reject results were
  unchanged on the checked corpus, but a field whose default route fails can
  now be accepted under GRH when the coarse value shows index one.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- A ``grh`` paired class/unit transaction that is accepted by the analytic
  index-one test now records that check on the published class group.
  ``analytic_class_regulator_status()`` reads ``verified``, and
  ``analytic_class_regulator_certification()`` gives the check's own
  conditionality: ``grh`` for a Belabas-Friedman ``hR`` and ``proven`` for
  the quadratic ``L(1, chi)`` value (no longer reached by a ``grh`` request;
  see the x^2 - d entry under Changed).  With a Belabas-Friedman ``hR`` the
  evaluation's BF audit data is recorded too, so ``zeta_bf_proof_status()``
  reads ``verified`` and ``zeta_bf_proof_record()`` gives its error bound,
  cutoff, and precisions.  The record is informational.  The transaction
  keeps the class-group and unit labels at ``grh``, and that change leaves
  the public zeta and BF gates unchanged.  The one exception: a ``proven``
  record stored by another route lets a later explicit
  ``try_certify_with_units(units, proven)`` with proven units promote the
  class group to ``proven``; a ``grh`` run on ``x^2 - d`` stores no such
  record.  The exact imaginary-quadratic and degree-one
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
  ``L(1, chi)`` route, which a ``grh`` run on ``x^2 - d`` no longer stores),
  and then they only record the GRH-conditional check; otherwise they return
  ``false`` and leave both objects unchanged.
  ``try_certify_class_unit_with_zeta_bf`` uses a Belabas-Friedman value for
  quadratic fields too, and in degree one the exact value ``hR = 1``.  Degree
  one and the quadratic ``L(1, chi)`` route of
  ``try_certify_class_unit_with_zeta``, whose ``hR`` is unconditional, are
  unchanged.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The unconditional unit proof (``OrderUnitGroup::prove_index_bound`` and
  ``proven`` class/unit validation) no longer fails at once when the
  ``ell``-th root of a residue-character kernel row cannot be decided, which
  happens mostly in degree 10 and higher.  As in Hecke's ``saturate!``, it treats
  the row like a non-power and rescans with more auxiliary primes until the
  kernel is empty.  ``proven`` class/unit computations that previously failed
  this way now succeed, for example ``Q(zeta_13)`` at ``ell = 17``,
  ``Q(zeta_17)``, ``Q(zeta_19)`` and ``Q(sqrt(-1), sqrt(-3), sqrt(5))``,
  with the class number and regulator of PARI's ``bnfinit`` and
  ``bnfcertify``.  The certification rule is unchanged: ``ell`` is verified
  only by an empty character kernel, and roots are adjoined only after
  exact verification.  The failure behavior changes: a row that really is an
  ``ell``-th power but whose root cannot be computed now fails as
  ``unavailable`` only after the ``2^16``-prime scan guard instead of at
  once, so such a failure takes longer.  See "Class groups and order units"
  in :doc:`../reference/algorithms_and_sources`.

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

- A ``grh`` class/unit request on a quadratic field written as ``x^2 - d``
  now takes a Belabas-Friedman ``hR`` for the analytic index-one test, first
  the coarse index-one value, as for every other presentation, instead of
  the unconditional ``L(1, chi_D)`` value, whose cost grows like
  ``sqrt(|D|)``.  Hecke (``Clgp.jl``, ``zeta_log_residue``) takes the
  Belabas-Friedman residue in every degree, and PARI's GRH-conditional
  ``quadclassunit`` a truncated Euler product approximation
  (``buch1.c``, ``compute_invresquad``).  On the real quadratic benchmark
  rows the analytic step fell from 38 to 138 ms to about 1 ms (``|D|`` of
  4e4 to 4e5; one run, not a timing claim).  The class number,
  regulator and generators are unchanged.  Three published informational
  fields of such a ``grh`` result change: ``class_unit_completeness_basis()``
  is ``belabas_friedman`` instead of ``unconditional_analytic``,
  ``class_unit_completeness_certification()`` and
  ``analytic_class_regulator_certification()`` are ``grh`` instead of
  ``proven``, and ``zeta_bf_proof_status()`` is ``verified`` instead of
  ``not_checked``, with the BF audit record stored.  Consequently a ``grh``
  run no longer stores a ``proven`` analytic record, a later
  ``try_certify_with_units(units, proven)`` computes ``L(1, chi)`` itself
  (at a usable precision) to promote, and ``unconditional_analytic`` is no
  longer reachable through public transactions.  ``proven`` requests and the
  public zeta gates are unchanged.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- ``Element::is_power`` and ``Element::is_square`` (and the factored-element
  root tests built on them) now answer ``false`` for many inputs they
  previously reported as unsupported, such as ``2^(1/3)`` as a cube in
  ``Q(2^(1/3))`` or ``2`` as a square in ``Q(sqrt(-3), sqrt(-5), sqrt(-7))``.
  Before the Hensel lift, the residue of the input is tested at every prime
  above the chosen prime and above up to 40 further unramified primes; a
  residue that is not an ``n``-th power proves the input is not an ``n``-th
  power.  In degree 10 and higher, where there is no Hensel lift, the
  same residue test runs at up to 40 unramified primes, so ``is_power`` and
  ``is_square`` there also return a definite ``false`` when a residue
  disproves the input.  Inputs with no such residue are handled as before,
  and ``true`` answers and their roots are unchanged.  See "Element powers
  and roots" in :doc:`../reference/algorithms_and_sources`.

- ``FactoredElement::logarithmic_embedding`` merges equal bases before
  summing, as Hecke's ``conjugates_arb_log`` sums over a ``FacElem``
  dictionary (``src/NumFieldOrd/NfOrd/FacElem.jl``, 149-205): each distinct
  factor contributes one term with the exact sum of its exponents (exponents
  that cancel contribute nothing), and the guard bits come from the merged
  exponents and the number of terms.  The result is still a rigorous Arb
  enclosure within the requested ``2^-precision``, but its bits differ from
  the entry-by-entry sum.  Unit candidates built from products of relation
  units repeat a few dozen bases thousands of times, so this shortens unit
  computations on such fields (for example ``grh`` on the maximal real
  subfield of ``Q(zeta_31)``).  The published regulator enclosure changes in
  its last bits: on the class/unit benchmark corpus, 36 of 237 ``grh`` and
  ``proven`` results report a different ``regulator_radius`` (between 0.71
  and 1.11 times the previous value) and some a regulator midpoint that
  differs in its last bits; every new enclosure overlaps the previous one.
  Units, class groups and labels are unchanged on that corpus.

- The relation-kernel basis that the class-relation saturation proof
  size-reduces its ``beta`` rows against no longer runs a full ``fmpz_lll``
  when the basis for the unit candidates is not already built.  For a kernel
  with at least two rows and entries of at most 250 bits, it takes FLINT's
  double-precision first phase (``fmpz_lll_d_with_removal_knapsack``) and
  accepts it only when an exact certificate shows it LLL-reduced with the
  default parameters, as ``Lat::lll_reduce`` does; otherwise it falls
  back to ``fmpz_lll``.  ``fmpz_lll``'s own certification of that phase
  dominated the proof on a 71 x 94 kernel for ``x^2 - x - 25000000``.
  With FLINT's default routing a certified first phase is exactly what
  ``fmpz_lll`` returns, so the basis, and with it every published result,
  is unchanged; this was also checked on the ``proven`` corpus and the
  class/unit test matrix.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- The compact power test of a unit or saturation candidate embeds each base
  once.  It bounds the conjugate logarithms of the candidate and of its
  formal residual, which have the same bases, and then reduces the chosen
  one to a compact presentation, which embeds it again.  The three steps
  previously used separate embedding contexts, so a base whose logarithm
  needs a high working precision (a base with large coefficients and a tiny
  conjugate) was evaluated up the whole precision ladder three times.  They
  now share one context per candidate, and the later steps take each base's
  logarithm from its cache, stored at the highest precision the first step
  reached.  The cached values are the same Arb enclosures, so every bound
  stays rigorous; published results are unchanged on the class/unit test
  corpus.

- Copying a factored value (``FactoredElement::set``) and multiplying or
  dividing a value into itself now share the factor storage and copy it only
  when a shared value is modified, instead of deep-copying every factor.
  Published results are unchanged.

- ``Order::maximal_order`` and ``Order::pmaximal_overorder`` apply Dedekind's
  criterion before Round 2 when the input is the equation order of a monic
  integral defining polynomial, and return the input order at primes where it
  is ``p``-maximal.  Results are unchanged.

- ``decompose_prime`` on a maximal order keeps the decomposition of each
  prime dividing the index ``[O : Z[theta]]`` with the order, before any
  residue-degree bound, and later calls for that prime on the same order
  rebuild the primes from it instead of splitting ``O/pO`` again.  When a
  prime has more prime ideals above it than ``p`` itself, the split no longer
  searches for a separating element, which cannot exist.  Decompositions,
  prime order included, are unchanged.

- The Bach-Friedman Euler factors of a ``grh`` class/unit computation read
  the residue degrees of a prime dividing the index from the decomposition
  the maximal order keeps, instead of rebuilding its prime ideals.  When the
  decomposition is not yet kept, the quotient ``O/pO`` of an unramified
  prime is built once for the bounded and the recursive split, and the
  decomposition kept is the one ``decompose_prime`` would keep.
  ``decompose_prime`` and published results are unchanged (identical on the
  ``grh`` corpus and the ``proven`` matrix compared).  On
  ``cyclotomic_real_q10_c11_deg10`` (``grh``) the class/unit CPU time went
  from a median of 354 ms to 331 ms over 11 interleaved runs (Release build
  on a busy machine, not a controlled measurement).

- The ``grh`` factor-base probe no longer decomposes a rational prime that has
  no prime ideal of norm at most the probe bound.  It reads the residue
  degrees of the prime from the defining polynomial modulo the prime, and when
  ``p`` raised to the smallest residue degree exceeds the bound it keeps a
  placeholder and decomposes the prime only if a later, larger bound admits
  it.  Factor bases and outputs are unchanged (identical on the 13 ``grh``
  rows compared).  On ``cyclotomic_real_163_deg81`` the probe went from 48.5 s
  to 8.9 s (busy machine, one sample each, not a controlled measurement).
  Behavior change: a skipped prime is never decomposed, so an internal error
  that decomposing it would have raised no longer fails the build; the prime
  contributes nothing to the factor base either way.

- The validation loop of a ``grh`` class/unit computation no longer repeats
  a bounded class-relation saturation pass that appended no relation while
  the candidate class number, the number of factor-base generators and the
  candidate regulator (within a certified factor two) are unchanged; it
  keeps the result of the earlier pass (``relation_saturation_status()``
  ``unavailable``) instead.  The skip assumes that the repeated pass would
  again append nothing, which is not proven; published results are
  unchanged on the class/unit benchmark corpus.  If a skipped pass would
  itself have failed, the computation now keeps the earlier ``unavailable``
  result instead of failing.  The skip never adds a proof label, and
  ``proven`` requests are unchanged.

- Unit-basis reduction of the stored units no longer lets their factored
  form grow without bound.  When a transform would give the stored units more
  than 20000 factors in total, or they already exceed that total (as after
  an adjoin), on an order known to be maximal the units are rebuilt as single
  field elements by the existing wide path (``rebuild_wide_unit_generators``,
  the port of PARI 2.17.3 ``buch2.c`` ``vec_chinese_units``, 1214-1240).  The
  cutoff is a Silex cost limit, not a source constant.  If the rebuild fails
  the factored units are kept.  Units, class groups, regulator values and
  labels are unchanged.  The published ``unit_group.regulator_radius``
  changes (still a valid enclosure; some radii shrink, some grow) because the
  regulator is computed from the rebuilt generators: on the class/unit
  benchmark corpus (degree at most 40, ``grh``) 10 fields differ in it, among
  them ``cyclotomic_37_phi36`` and ``cyclotomic_64_phi32``.  No other published field differs; in ``proven`` mode
  no compared field differs.

- A ``proven`` class/unit computation over an imaginary quadratic field
  reduces the very large binary forms of its witness step (more than 9000
  bits) with a port of the recursive fast reduction of PARI 2.17.3
  (Schönhage, ISSAC 1991; see the source map).  For ``x^2 + 100000007``
  the witness step falls from about 25--31 s to 9--10 s and the total from
  27--35 s to 11--13 s (release build on a busy machine; the controlled
  comparison is pending).  Relation rows, norms, labels and instance JSON
  are unchanged on the test corpus.  Above the switch the reduced form is
  the same, but ``Relation::generator`` can differ in sign (171 of the 394
  witnesses at ``x^2 + 100000007``), or by a unit when ``D = -3`` or
  ``-4``; below it nothing changes.  The port
  does not copy PARI's treatment of forms with ``b < 0``, which yields a
  transformation of determinant ``-1``.  If a step of the fast path fails,
  or the input form does not take the reduced minimum at the computed
  vector, the reduction falls back to the earlier basecase reduction, so
  published results and failures are as before.

- When ``decompose_prime`` splits ``O/pO`` for an unramified prime dividing
  the index, it no longer tests that each new prime ideal's lattice is an
  ideal in builds with ``NDEBUG`` defined (for example Release).  The
  lattice is ``pO`` plus sums of ideals ``e(O/pO)`` of ``O/pO``, so it is an
  ideal by construction.  Builds with assertions enabled still assert the
  test.  Decompositions are unchanged.

- ``FractionalIdeal::set_integral_den`` and ``set_integral``, the
  maximal-order ``FractionalIdeal::invert`` and the scaling step of
  fractional ``add`` and ``intersect`` no longer test that a lattice is an
  ideal in builds with ``NDEBUG`` defined (for example Release) when the
  lattice is an ideal by construction: an ``Ideal`` argument divided by its
  common content with the denominator, a positive integer multiple of an
  ``Ideal``, or the colon ``(O : I)`` scaled into ``O``.  ``Ideal::set_hnf``
  still rejects a matrix whose lattice is not an ideal, and builds with
  assertions enabled still assert the test.  Results are unchanged.

- The LLL relation search reduces each ideal once, with integer LLL on the
  ideal basis mapped through a rounded embedding of the order that is
  computed once per order, as PARI 2.17.3 ``Fincke_Pohst_ideal``
  (``buch2.c`` 2461-2462) does with the rounded ``M`` of ``RM_round_maxrank``
  (``base4.c`` 2631-2640).  It previously built a floating-point T2 Gram
  matrix for every ideal.  When the reduced basis fails the T2 bounds, or
  the per-order data cannot be formed, the ideal is reduced by the earlier
  route instead.  The reduced bases differ, so the relation search
  tries other candidates and ``grh`` results change in how they were found:
  class groups, regulator values and labels are unchanged on the test corpus,
  while ``relation_count``, ``relation_kernel_unit_count``,
  ``relation_sources``, ``skipped_dependent_relation_count``, the chosen
  fundamental units, ``unit_group.regulator_radius`` (still a valid
  enclosure) and the ``zeta_bf_proof_record()`` cutoff and error bound
  (still GRH-conditional) can differ.  For the degree-12 field
  ``x^12 - x^11 - x^9 + ... + 1`` the coarse index-one check now succeeds,
  so the record has the coarse cutoff 2547 instead of the default route's.
  ``cyclotomic_56_phi24`` reports
  ``relation_saturation_status`` ``not_checked`` instead of ``unavailable``.
  Of 40 compared ``grh`` runs (the test-suite class/unit field manifest and
  eight further fields up to degree 40), 22 differ in these fields, 8 of them
  in ``regulator_radius``; the 19 ``proven`` runs are unchanged.
  One row is slower: for ``x^3 + x + 200`` in ``grh`` mode
  (``cubic_x3_x_200_grh``) the first unit candidate fails the coarse
  Belabas-Friedman check, so the run refines its validation and saturates
  the class group, and takes about 26 ms instead of 9 ms (release build on a
  busy machine).
- The analytic finish check of class-group computations, which runs at each
  relation checkpoint, keeps its embedding context between checks at the
  same precision instead of rebuilding it, and the regulator product reuses
  the unit-log matrix and pivot selection that the independent-unit check
  just built from the same relations.  The kept context is released when
  the check finishes.  This is a speed-up for computations that run the
  check many times; outputs are unchanged.

- The one-LLL reduction of the LLL relation search now checks, for the first
  ideal in each order whose one-LLL basis is accepted and whose two-LLL
  reduction also succeeds, that the reduced basis is about as short
  in T2 as the basis of the earlier two-LLL route: the product of its T2
  diagonal entries may exceed the two-LLL product by at most ``2^n`` for
  degree ``n``.  When the check fails, that ideal and every later ideal of
  the order use the two-LLL route.  Silex rounds the embeddings of its order
  basis, which, unlike PARI's T2-reduced integral basis, can be badly
  conditioned; on the power basis of the real subfield of ``Q(zeta_136)``
  (``cyclotomic_real_136_deg32``) the rounded embedding steered the LLL to
  bases that found almost no relations, so ``grh`` ran for more than 1200 s
  instead of about 8 s.  It now takes about 7.9 s (release build on a busy
  machine).  Where the check passes, the reduction and all results are
  unchanged; it costs one two-LLL reduction per order and class-group
  context, or more when that reduction itself fails and the check repeats on
  the next ideal.  The check runs once per order in a class-group context, on the
  first accepted ideal, and its result is shared by all relation slices of
  that context; each slice does not repeat it on its own first ideal.
  Of the 172 ``grh`` corpus rows (the comparison did not cover the fixture
  rows), three more take the two-LLL
  route for some orders: ``cyclotomic_real_37_deg18``,
  ``cyclotomic_real_41_deg20`` and ``multiquadratic_p5_m7_m11``.  Their class groups, regulator values and
  labels are unchanged, while ``relation_count``,
  ``relation_kernel_unit_count``, ``relation_sources``,
  ``skipped_dependent_relation_count``, ``unit_group.regulator_radius``
  (still a valid enclosure) and the chosen fundamental units can differ.
  ``cyclotomic_real_41_deg20`` falls from about 1.00 s to 0.44 s and
  ``multiquadratic_p5_m7_m11`` from 0.12 s to 0.05 s.  ``proven`` runs are unchanged on the same rows.  Because the check runs
  once per context, ``cyclotomic_real_31_deg15`` (whose first ideal passes
  the check, while a later slice's first ideal failed it narrowly) keeps the
  one-LLL route and returns to the outputs and time it had before this
  check.

- ``roots_of_unity``, ``root_of_unity_order`` and ``root_of_unity_generator``
  remember their result per thread for the last few defining polynomials,
  so repeated queries on the same field, such as those made during a
  class-group and unit computation, no longer repeat the good-prime search
  for ``w``.  Results and failures are unchanged: an entry is keyed by the
  defining polynomial and backend, the computation reads nothing else, and
  failures are not stored.  ``grh`` on ``Q(zeta_45)`` falls from about 0.70 s
  to 0.61 s and on ``Q(zeta_72)`` from about 1.23 s to 1.16 s (release
  build, three runs each on a busy machine; the controlled comparison is
  pending).  Instance JSON is unchanged on the ``grh`` corpus and the 65
  ``proven`` benchmark rows.

- The factored-element root tests and the compact power roots of the
  class-group and unit computation verify a Hensel-lifted root against a norm
  computed from the small factors of the input, as the product of the factor
  norms to their exponents, instead of the norm of the expanded input, which
  on ``proven`` runs of degree 6 to 9 can have tens of thousands of bits.
  The value is the same norm, so every answer, root and failure is
  unchanged; only the time spent rejecting candidate roots falls.

- The large-prime partial relations of the class-group relation search test
  which factor-base primes divide an element with word-size arithmetic
  modulo the rational prime when it fits in a machine word, instead of
  multi-precision modular polynomials.  The test is the same, so outputs are
  unchanged: instance JSON is identical on the ``grh`` corpus and the 65
  ``proven`` benchmark rows.  ``grh`` CPU time on ``Q(zeta_41)`` falls by
  about 4% (7.73 s to 7.40 s) and on ``Q(zeta_37)`` by about 3% (2.63 s to
  2.55 s) (release build, medians of nine paired runs on a busy machine; the
  controlled comparison is pending).

- The exact HNF context ``fmpz_smat::HnfContext``, which the class-group
  relation search uses, folds newly added rows into the HNF it already holds
  when that HNF has less than full rank, instead of recomputing the HNF of
  all rows with ``fmpz_mat_hnf``.  The row HNF is unique, so the
  lattice, rank, index and every dependent-row decision are unchanged:
  instance JSON is identical on the 65 ``proven`` and the 183 completed
  ``grh`` benchmark rows.  ``grh`` CPU time on ``Q(zeta_64)`` falls by about 4%
  (1.43 s to 1.37 s) and on ``Q(zeta_41)`` by about 2% (4.98 s to 4.88 s)
  (release build, medians of nine paired runs on a busy machine; the
  controlled comparison is pending).

- ``FiniteAbelianGroup`` computes the inverse of its Smith right transform
  on the first call of ``invariant_generator_matrix`` or ``reduce`` instead
  of on every ``set_relation_matrix``, so class-group candidates that publish
  a presentation several times during relation collection no longer invert
  that matrix each time (four inversions fewer on ``grh`` runs of
  ``Q(zeta_64)`` and ``Q(zeta_41)``).  Outputs are unchanged: instance JSON
  is identical on the ``grh`` corpus and the 70 ``proven`` benchmark rows.
  One failure point moves: the right transform from FLINT's Smith form is
  unimodular, but if its inverse were not integral, ``set_relation_matrix``
  would now succeed and the first ``invariant_generator_matrix`` or
  ``reduce`` call would return ``false`` instead, before writing any output.
  ``FiniteAbelianGroup`` has a new private member.  As with the existing
  left-transform cache, ``const`` calls may fill internal caches, so one
  object must not be used from several threads at once.

- ``ResidueField::set`` copies a defined residue field's prime,
  characteristic and modulus instead of re-running the primality and
  irreducibility checks, which the field already passed when
  ``set_prime`` defined it.  Every ``ResidueFieldElement`` constructor
  copies its parent this way, so residue-field arithmetic and
  ``ResidueFieldQuotientLog::apply`` no longer repeat those checks per
  call.  ``set_prime`` and the ``PrimeIdeal`` constructor still validate
  their input, and results are unchanged on the ``grh`` class/unit corpus
  and the 70 ``proven`` benchmark rows.  On ``cyclotomic_19_phi18`` in
  ``grh`` mode the computation takes about 13% less time (release build on
  a busy machine, five runs each).

- Class-relation saturation (with and without units) installs the evaluated
  ``ell``-th root with the row it already knows, the kernel combination of
  the stored rows divided by ``ell``, instead of factoring the root over the
  factor base again.  The root is exact and valuations are additive, so the
  row is the one factoring would give; the saturation pass of the class/unit
  pipeline and the reference ``saturate!`` install their roots the same way.
  One consequence: a resource or factoring failure in that re-factoring can no
  longer fail the saturation pass.  Instance JSON is identical on the 184
  ``grh`` corpus rows and the 70 ``proven`` benchmark rows.  ``grh`` CPU time
  on ``real_quadratic_d100000000000009`` falls from about 0.53 s to 0.42 s
  (about -20%; release build on a busy machine, measured after the rebase
  with ``cyclotomic_real_41_deg20`` as a control row that stayed flat; the
  controlled comparison is pending).
- ``PrimeIdeal::valuation`` on an ``Ideal`` with a known principal generator
  ``a`` reuses the power of ``p`` in the ideal norm, which equals
  ``|N(a)|``, instead of computing ``N(a)`` again for every prime.  Class-group
  relation candidates in maximal orders that are not the equation order take
  this route for every factor-base prime they are tested against.
  Valuations are unchanged: instance JSON is identical on the 184 ``grh``
  corpus rows and the 70 ``proven`` benchmark rows.  ``grh`` CPU time on
  ``cyclotomic_real_q10_c11_deg10`` falls by about 9% (0.64 s to 0.58 s)
  (release build, medians of 15 paired runs on a busy machine; the
  controlled comparison is pending).

- ``ResidueFieldQuotientLog::apply`` at a prime of residue degree
  ``f > 1`` takes the logarithm through the norm to ``F_p`` when ``ell``
  divides ``p - 1``: ``N(x)^((p - 1)/ell) = x^((p^f - 1)/ell)``, with the
  norm computed as a resultant, instead of raising ``x`` to
  ``(p^f - 1)/ell`` in ``F_(p^f)``.  It returns the same logarithm, and
  results are unchanged on the ``grh`` class/unit corpus and the 70
  ``proven`` benchmark rows.  On ``cyclotomic_19_phi18`` in ``grh`` mode,
  where nearly all of these logs are at residue degree 3 to 18 with
  ``ell = 2``, the computation takes about 28% less time (release build on
  a busy machine, eleven interleaved runs each).

- Unit extraction for ``grh`` class/unit computations tests each candidate
  unit for independence incrementally, from the Gram matrix and inverse of the
  units already kept (Schur complement identity), instead of rebuilding the
  whole logarithmic Gram determinant per candidate.  Wherever the previous test
  decided, the decision is the same.  A candidate it dropped at the precision
  cap, or after a failed logarithmic embedding, may now be decided at higher
  precision (both need about ``2^20`` bits).  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- A ``grh`` class/unit computation no longer rebuilds the units in its
  validation loop when the units were just built from the same class
  context and that context has not changed since.  Hecke finds units once
  per outer iteration and does not rebuild them during validation
  (``Clgp.jl``, ``_class_unit_group`` and ``_validate_class_unit_group``).
  The rebuild drew a fresh random subset of relations, so skipping it
  shifts later random draws.  Class groups, class generators, regulator
  midpoints and certification labels are unchanged on every ``grh`` corpus
  row that completed (``cyclotomic_real_163_deg81`` timed out on both builds
  and was not compared).  The published unit generators change on five rows,
  and on each the new generators generate the same unit group as the old:
  ``cyclotomic_35_phi24``, ``cyclotomic_52_phi24``,
  ``cyclotomic_real_41_deg20``, ``cyclotomic_real_q10_c11_deg10`` and
  ``simplest_cubic_a12``.  On all of these except ``cyclotomic_52_phi24``
  the regulator radius changes too.

Fixed
-----

- A ``grh`` class/unit computation whose relation-search factor base exceeds
  the size limit of ``2^14`` generators now fails with ``failure_reason``
  ``factor_base_size_limit_exceeded``, as documented, instead of
  ``class_unit_computation_failed``.  The relation search of fields of degree
  at least 2 built its base on a context without the run's transaction, so
  the refusal was not recorded (for example ``x^2 - (10^140 + 1399)``, whose
  GRH base has more than ``2^14`` generators).  The same inputs fail as before; results,
  labels and successful runs are unchanged.

- ``FactorBase::index`` and ``FactorBase::contains`` compare prime ideals as
  ideals.  A prime equal to a base prime but recorded with other residue data
  (another residue polynomial or variable convention, or none) was reported
  absent; it is now found at the base prime's index.  Factor-base builds
  refuse to list a prime ideal twice, which the block-completeness checks of
  factorization over a base assume; no builder produced a repeat.
  Test-suite outputs are unchanged.

- The class-relation saturation proof with proven units size-reduces the
  relation coefficients of each invariant-generator ``beta`` row against an
  LLL-reduced basis of the relation kernel before forming its compact
  witness.  The
  presentation's coefficients can exceed the reduced size by orders of
  magnitude: in
  ``t-silex-class-unit-matrix`` the largest raw coefficients are 414 bits
  and the reduced ones 12 bits or fewer (measured at commit 49dfe1d on a
  busy machine).  The reduced row differs from the original by a relation-kernel unit, so the ``ell``-th
  power rank modulo ``ell`` is unchanged.
  Published class groups, units and labels are unchanged on the test
  manifest.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- Proven mode no longer fails in the finite-ideal T2 enumeration on fields
  such as the real cyclotomic field of conductor 11 and degree 10, where the
  ideal-lattice reduction ran on embeddings rounded too coarsely and every
  precision retry repeated the same input.  The reduction now repeats from a
  finer rounding scale when its result fails a Gram-Schmidt check or Hecke's
  ideal-LLL acceptance bounds.  The lattice and the bound rule are unchanged;
  the bound value and the enumeration order follow the new reduced basis.  On
  rows whose ideal bases are now reduced better, the published
  ``regulator_radius`` can change (still a valid enclosure; the regulator,
  units, class group and labels are the same).

- The class-relation saturation candidate search no longer uses an auxiliary
  prime ``q`` that divides ``lc(F) * disc(F)`` when the defining polynomial is
  not monic with integer coefficients (``F`` is its integer numerator, which
  keeps the content).  A root of the polynomial modulo such a ``q`` need not give a
  degree-one prime, and its ``ell``-th-power character could rule out a true
  ``ell``-th power.  Skipping the prime only leaves the candidate space
  larger.  Monic integral polynomials, including every field in the test
  manifests, are unaffected.

- A ``grh`` class/unit computation now improves its candidate units when the
  a priori unit index bound ``R / R_low`` does not fit in a machine word, as
  for ``Q(zeta_41)`` (degree 40).  Such a computation previously never
  adjoined a dependent unit and kept extending relations; the grh run of
  ``Q(zeta_41)`` did not finish and now succeeds with class group
  ``[11, 11]`` and unit rank 19, agreeing with PARI.  The search then uses
  the GRH-conditional class/unit index bound that validation computes.
  ``proven`` requests, and ``grh`` computations whose a priori bound fits,
  are unchanged.  See
  "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- A class/unit computation no longer exhausts memory while extracting units
  from the relation kernel when the raw kernel coefficients do not fit in a
  machine word, as for ``x^2 + 2x - 104002`` and ``x^2 - x - 25000000``, which
  now finish ``proven`` with ``h = 1`` and ``h = 720``.  Unit candidates now
  come from an LLL-reduced basis of the relation kernel, as in Hecke's
  ``find_candidates``; the reduced basis spans the same kernel lattice.  An
  intermediate unit-group candidate may differ from before, and is often
  closer to a fundamental system; no published unit basis changed on the
  test or benchmark corpora.  In the unit-extraction power routine, a power of
  a factored element with an exponent outside a machine word now fails
  instead of being expanded unless the base is a root of unity, which is
  raised to the exponent modulo the order of the roots of unity; the
  separate copy used by the relation-HNF witness still expands such powers.
  The setters ``OrderUnitGroup::set_relation_kernel_units*`` may return
  different generators and stronger (sound) labels than before: in
  Q(sqrt 2), ``set_relation_kernel_units_index_bounded`` now publishes
  ``proven`` where it published ``unknown``.

- The unit computation no longer fails with a class-unit computation error
  when the unit-basis transform needs exponents that do not fit in a machine
  word, as in the grh run of the real cyclotomic field of conductor 136
  (degree 32).  The transform now uses exponents of any size, and
  ``OrderUnitGroup::free_generator`` returns an equal unit rebuilt in
  word-size form, a documented deviation from Hecke.

- ``decompose_prime`` on a maximal order no longer fails at a prime ``p``
  dividing the index when ``p`` raised to the number of primes above ``p``
  exceeds 4096, as for ``p = 47`` in the sextic field of LMFDB 6.6.46411625.1
  (six primes) or ``p = 11`` in ``Q(sqrt(-3), sqrt(-5), sqrt(-7))`` (four
  primes).  When the bounded enumerations and the monogenic-generator route
  both fail, it now splits the residue algebra recursively, as PARI's
  ``idealprimedec`` does; decompositions that succeeded before are
  unchanged, prime order included.  The ``grh``
  class/unit computations of seven corpus fields that failed when choosing
  the factor-base bound now succeed with the expected class numbers and unit
  ranks.  See "Local algebra and relations" in
  :doc:`../reference/algorithms_and_sources`.

- A ``proven`` class/unit computation no longer fails when an inert rational
  prime lies below the proven generation requirement.  Such a prime has a
  single prime ideal ``pO_K``, which is principal; the generation proof now
  treats it as covered, as it already did for the last unramified prime
  ideal of a fully listed decomposition, instead of requiring it in a
  relation-completion base that omits inert primes.  The simplest cubic fields
  ``x^3 - a x^2 - (a + 3) x - 1`` for ``a = 6, 8, 10`` now succeed with
  ``proven`` class numbers 3, 1, 1 and unit rank 2, matching PARI.  A totally
  ramified prime with one prime ideal is still checked.  See "Class groups
  and order units" in :doc:`../reference/algorithms_and_sources`.

- A ``proven`` class/unit computation no longer fails closed when its
  relation-saturation proof at a prime ``ell`` needs an auxiliary prime
  ``q = 1 mod ell`` above 1000.  The proof now scans ``q`` with no fixed
  bound until it reaches full rank, as PARI's ``bnfcertify`` does, and checks
  every ``q`` it uses as before.  Only the caller's time or memory limit stops
  a longer scan, and a stopped request publishes no result.  For example,
  ``Q(sqrt(1575029))`` now succeeds with ``proven`` class group ``Z/89`` and
  the regulator PARI reports; it needs ``q = 2671``.  Results that succeeded
  before are unchanged.  See "Class groups and order units" in
  :doc:`../reference/algorithms_and_sources`.

- Class/unit computations no longer fail in relation search when the order's
  Minkowski embedding has large entries, as for the power bases of the real
  cyclotomic fields of conductor 37 and 41 (degree 18 and 20), which now
  succeed with ``grh``.  The ideal-lattice reduction raised its working
  precision without ever meeting its rounding check; it now raises the
  embedding precision by the measured shortfall.  Computations whose
  reductions passed that check are unchanged.  See "Class groups and order
  units" in :doc:`../reference/algorithms_and_sources`.

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
  unconditional check stored earlier by another route,
  ``try_certify_class_unit_with_zeta``, or, in degree one, where
  ``hR = 1`` exactly, ``try_certify_class_unit_with_zeta_bf``) now
  reports ``relation_saturation_status()`` as ``verified``, as documented:
  index one with verified generation means the relations are saturated at
  every prime.  It was ``not_checked`` on that route, for example on a real
  quadratic field with class number greater than one such as
  ``x^2 - 12905``, and on an imaginary quadratic field certified through
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
