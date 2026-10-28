Class number route cutoff
=========================

A ``grh`` class/unit request on an imaginary quadratic field takes the class
number from one of two routes, chosen by the bit length of the discriminant
(``grh_imaginary_quadratic_lattice_min_abs_discriminant_bits`` in
``src/class_group/class_group.cpp``, value 35):

- ``|D| < 2^34``: the exact class number, counted from the reduced forms of
  discriminant ``D``, in time growing like ``sqrt(|D|)``.  The completeness
  basis is ``exact_class_number``.
- ``|D| >= 2^34``: the GRH relation lattice and the analytic index-one test,
  which never lists forms.  The completeness basis is ``belabas_friedman``.

The cutoff was ``2^44`` before it was lowered to ``2^34`` (see "Lowering the
cutoff to 2^34" below).  The exact class number ``h(D)`` that ``proven``
requests use has its own cutoff, ``2^54``
(``imaginary_quadratic_count_max_abs_discriminant_bits``, raised from ``2^44``;
see "Raising the proven cutoff to 2^54" below): below it ``h(D)`` is counted
from the reduced forms, and from ``2^54`` on it comes from the class number
formula with a rigorous ``L(1, chi_D)`` enclosure (see "Class groups and order
units" in :doc:`../reference/algorithms_and_sources`).

Counting instead of listing
---------------------------

Both cutoffs were measured when the reduced forms were listed by FLINT
``qfb_reduced_forms``, whose time and memory grow like ``sqrt(|D|)`` (about
1.4 GB at ``|D| = 1e14``).  For a fundamental ``D``, which is every
discriminant of a maximal order and so every discriminant these routes see,
Silex now counts the reduced forms without building them (the source chain is
in "Class groups and order units" in
:doc:`../reference/algorithms_and_sources`); a non-fundamental ``D`` keeps the
listing.  The count needs about 5 bytes per ``a <= sqrt(|D| / 3)``, about
12 MB just below ``2^44``, so the memory reason for the ``2^44`` cutoff no
longer applies, and it is faster: on ``imaginary_quadratic_d10000000799``
(``|D|`` about ``2^33.2``) callgrind counts 298.9M instructions for ``h(D)``
by listing and 39.3M by counting.  The time comparisons below, against the
class number formula and against the GRH relation lattice, were taken with
the listing.  The ``grh`` cutoff stands unchanged; the ``proven`` cutoff was
re-measured with the count and raised (last section).

The 2^44 cutoff was set by a first measurement (a busy machine; the lattice
route's index-one test then used a different ``hR`` below ``2^44``) and was
re-measured with the faster listing of T-193 (below) before it was lowered to
``2^34`` on an idle machine (last section).  The sections "Setup" and
"Results" describe the re-measurement; it was taken with the cutoff at
``2^44``, so its "lattice" columns are forced routes, not the shipped
behavior.

Coupling to the zeta cap
------------------------

The shipped state is class-group cutoff ``2^34`` with zeta cap ``2^44``
(``quadratic_route_max_abs_discriminant_bits`` in ``src/zeta/zeta.cpp``,
value 44, line 905).  Whether the lattice route's index-one test can take
``hR`` from the unconditional analytic evaluation (the AFE of ``L(1, chi_D)``)
below the zeta cap is the question here.  The forced-route re-measurement below
found that, with the cutoff at ``2^44`` and the zeta cap at 44 in scratch trees,
only the index-2 presentation ``x^2 + n`` showed a slow lattice route, 23 to
79 s where listing took 0.3 to 1.2 s; the AFE was inferred from the
completeness basis field and the cause was not investigated.

A later check on main (T-312, busy machine, ``taskset -c 0-3``, median of 3)
did not reproduce this.  In ``class_group.cpp`` a ``grh`` request with
``|D| >= 2^34`` leaves the exact class-number route (line 5421), and
``try_certify_class_unit_with_bf_audit`` (lines 5144-5156) sets
``hr_unconditional`` to ``(degree == 1)``, false for a quadratic field, so a
``grh`` request on the lattice route gets no ``unconditional_analytic`` record
and takes the Belabas-Friedman ``hR``.  The AFE (``quadratic_residue``,
``zeta.cpp`` line 941, gated at line 932) is reached from the ``proven``
class-unit path and the zeta residue path.  Index-1 and index-2 fields of
about ``2^36``, ``2^40`` and ``2^43`` all ran in under 0.35 s with basis
``belabas_friedman``, and T-242's idle sweep on ``x^2 + n`` agrees (1.9 s at
bin 43.9).  So the shipped build does not pay the AFE cost in this range.
Why the scratch trees of the re-measurement did take 23 to 79 s is unexplained;
something in their ``hR`` selection may have differed from the shipped code,
but this was not confirmed.  Lowering the zeta cap would still remove the
unconditional ``L(1, chi)`` route for ``proven`` requests between the new cap
and ``2^44``, which is a scope decision.

Setup
-----

Each field is ``Q(sqrt(-n))`` with ``n`` squarefree and ``n = 3 mod 4``, so
``D = -n``.  Each run is one ``silex-class-unit-instance --mode grh``
process, whose wall time is the median over 7 runs (5 in the sweep below
``1e12``), interleaved per repetition with the order rotating.  The source was
main at 4e7433a (which includes T-193).  Because no option forces the route,
three scratch Release builds differing in one or two constants were used:

- "listing": the class-group cutoff set to 99, so forms are listed at every
  ``|D|`` (equal to the shipped build below ``2^44`` at 4e7433a);
- "lattice": the cutoff set to 0, so the lattice is used at every ``|D|``,
  with the shipped ``hR`` choice (the AFE below the zeta cap, Belabas-Friedman
  at or above it; equal to the shipped build above ``2^44``);
- "lattice BF": the cutoff and the zeta cap both set to 0, so the lattice is
  used with the Belabas-Friedman ``hR`` at every ``|D|``.

Each field was run in four presentations: ``x^2 + x + (n+1)/4`` (equation
order index 1, the maximal order, "lucky"), ``x^2 + n`` (index 2),
``x^2 + x + (9n+1)/4`` (index 3) and ``x^2 + 9n`` (index 6).  The tables give
the lucky (index 1) and the index-2 presentation.  ``h`` and the class-group
invariants agreed in all trees and presentations for every field.  The
machine was busy (``machine-idle.sh`` reported BUSY at the start and IDLE at
the end, after all runs; load average 5 to 8 on 24 cores), so these are
busy-machine timings, not a controlled benchmark.  Memory was not measured.

Results
-------

Fields with ``|D|`` below ``2^44``, wall time in seconds (index 1 / index 2).
The AFE cost appears only in the index-2 "lattice" column.

.. list-table::
   :header-rows: 1

   * - ``n``
     - bits of ``|D|``
     - ``h``
     - listing
     - lattice
     - lattice BF
   * - 1000000000003
     - 40
     - 124568
     - 0.60 / 0.43
     - 0.51 / 22.9
     - 0.50 / 0.36
   * - 3000000000003
     - 42
     - 412512
     - 0.39 / 0.33
     - 0.21 / 27.7
     - 0.21 / 0.14
   * - 10000000000003
     - 44
     - 497056
     - 1.03 / 1.16
     - 0.74 / 64.6
     - 0.72 / 0.88
   * - 15000000000007
     - 44
     - 1331780
     - 0.80 / 0.82
     - 0.29 / 78.6
     - 0.29 / 0.31
   * - 17592186044415
     - 44
     - 3384192
     - 0.99 / 1.03
     - 0.12 / 52.3
     - 0.12 / 0.17

Fields with ``|D|`` at or above ``2^44``: the shipped lattice route (equal to
"lattice") against listing.

.. list-table::
   :header-rows: 1

   * - ``n``
     - bits of ``|D|``
     - ``h``
     - listing
     - lattice
   * - 18000000000003
     - 45
     - 888096
     - 0.68 / 0.75
     - 0.22 / 0.28
   * - 30000000000003
     - 45
     - 1290536
     - 0.93 / 0.94
     - 0.32 / 0.34
   * - 100000000000003
     - 47
     - 1425472
     - 1.64 / 2.12
     - 0.72 / 1.18
   * - 100000000000319
     - 47
     - 15815072
     - 2.94 / 2.87
     - 0.18 / 0.12

At ``n = 10^14 + 319``, ``h = 15815072`` agrees with PARI/GP 2.17
``quadclassunit``.  Below ``1e12``, listing against "lattice BF" (index 1 /
index 2, seconds) shows a tie within the noise at ``1e8`` to ``1e10`` and a
5 to 10% gain for the lattice at ``1e11``:

.. list-table::
   :header-rows: 1

   * - ``n``
     - bits of ``|D|``
     - ``h``
     - listing
     - lattice BF
   * - 100000003
     - 27
     - 1702
     - 0.019 / 0.022
     - 0.019 / 0.025
   * - 1000000003
     - 30
     - 3680
     - 0.047 / 0.070
     - 0.054 / 0.073
   * - 10000000003
     - 34
     - 10538
     - 0.079 / 0.111
     - 0.073 / 0.108
   * - 100000000003
     - 37
     - 31057
     - 0.262 / 0.375
     - 0.233 / 0.351

Listing now takes 0.3 to 1.2 s for ``|D|`` from ``1e12`` to ``2^44`` in every
presentation, and 0.7 to 2.9 s up to ``1e14``.  The earlier table had 27 to
91 s there, which came from the listing before T-193.

The lattice with the shipped ``hR`` choice is slower than listing below
``2^44`` only on the presentations that reach the AFE, here the index-2
presentation ``x^2 + n`` (23 to 79 s against 0.3 to 1.2 s); the index 1, 3 and
6 presentations take the Belabas-Friedman route even there, and are as fast
as or faster than listing from ``1e12``.  Why only index 2 reaches the AFE
was not investigated.  The lattice with the Belabas-Friedman ``hR`` is at
least as fast as listing at every presentation from ``1e12``, by about 1.2 to 1.9x
at ``1e12`` and ``1e13`` and up to 8x at ``2^44 - 1``.  At and above ``2^44``
the lattice is faster than listing in every row, by 1.5 to 25x.

Lowering the cutoff to 2^34
---------------------------

The measurement on an idle machine (``machine-idle.sh`` reported IDLE before
the build, before the sweep and after it) used two Release builds that differ
only in the constant (45 and 34), one process per run, GCC, FLINT 3.6.0.  It
covered 11 bins of ``n`` for ``x^2 + n`` (the bins are labelled by ``n``;
for ``n = 3 mod 4`` this is ``|D|``, otherwise ``|D| = 4n``) (bit length 30, 32, 33, 34, 35, 36,
38, 40, 42, 43, and 43.9), 8 fields per bin, 5 interleaved repetitions per
field with the order alternating.  The class-group invariants agreed between
the routes on all 88 fields, and the completeness basis was
``exact_class_number`` on the exact route and ``belabas_friedman`` on the
lattice route.  The ratio is lattice time over exact time, as sums of
per-field medians:

.. list-table::
   :header-rows: 1

   * - bin (bits of ``n``)
     - exact: s
     - lattice: s
     - ratio
     - lattice wins
   * - 33
     - 0.562
     - 0.494
     - 0.878
     - 6 of 8
   * - 34
     - 0.701
     - 0.583
     - 0.832
     - 7 of 8
   * - 35
     - 0.724
     - 0.539
     - 0.745
     - 8 of 8
   * - 36
     - 0.875
     - 0.564
     - 0.645
     - 8 of 8
   * - 38
     - 1.808
     - 1.274
     - 0.705
     - 8 of 8
   * - 40
     - 2.669
     - 1.487
     - 0.557
     - 8 of 8
   * - 42
     - 4.913
     - 2.457
     - 0.500
     - 8 of 8
   * - 43
     - 5.310
     - 1.790
     - 0.337
     - 8 of 8
   * - 43.9
     - 7.460
     - 1.914
     - 0.257
     - 8 of 8

The lattice route is faster on the sum in every bin from 33 up and on every
field from bin 35 up (bins by ``n``).  In bin 34 one field is about 9% slower, within its
repetition spread.  Bin 33 has a per-field ratio up to 1.19, so the evidence
below 2^34 is weaker; the cutoff ``2^34`` is the conservative choice, and
``|D| < 2^34`` keeps the exact route.  The two routes agree on the invariants
at every field measured.  A test pins the route and the completeness basis at
``|D| = 2^34 - 1`` (exact) and ``|D| = 2^34 + 3`` (lattice)
(``test_imaginary_quadratic_grh_refine_uses_exact_index`` and
``test_imaginary_quadratic_grh_lowered_cutoff_boundary`` in
``test/t-class-unit-matrix.cpp``).

Raising the proven cutoff to 2^54
---------------------------------

With the count in place of the listing, the ``proven`` cutoff was
re-measured by forcing each route in copies of a profile-release build (the
compare immediate of the cutoff patched to 0 for the formula at every
``|D|`` and to 127 for the count at every ``|D|``), on ``x^2 - x + c`` with
``D = 1 - 4c`` fundamental.  The instruction counts are ``perf stat -e
instructions:u`` of the whole ``silex-class-unit-instance --mode proven``
process (at ``|D| = 1e14 + 319`` they agree with callgrind of the class/unit
computation within 1%), and the peak memory is the process's maximum
resident set.  The machine was busy, so the wall times are indicative only.
Both routes gave the same full output (the instance JSON without timings) and
exit code on every row:

.. list-table::
   :header-rows: 1

   * - ``|D|``
     - bits
     - count: instructions
     - formula: instructions
     - ratio
     - peak memory, count / formula
   * - 100000000000319
     - 47
     - 55.7G
     - 362.8G
     - 6.5
     - 49 / 38 MB
   * - 100000000000015
     - 47
     - 59.9G
     - 315.0G
     - 5.3
     - 49 / 38 MB
   * - 1000000000000003
     - 50
     - 223.5G
     - 1293.8G
     - 5.8
     - 109 / 40 MB
   * - 1000000000000023
     - 50
     - 188.8G
     - 901.9G
     - 4.8
     - 108 / 38 MB
   * - 10000000000000003
     - 54
     - 734.9G
     - 4055.8G
     - 5.5
     - 297 / 40 MB
   * - 10000000000000023
     - 54
     - 609.7G
     - 3184.4G
     - 5.2
     - 297 / 40 MB
   * - 100000000000000003
     - 57
     - 2871.1G
     - 15584.6G
     - 5.4
     - 893 / 40 MB

The count is 4.8 to 6.5 times cheaper at every size measured, with no sign
of a crossover up to ``1e17``, so time does not bound the cutoff in this
range; memory does.  The count's peak grows like ``sqrt(|D|)`` (about 5 bytes
per ``a <= sqrt(|D| / 3)``, roughly three times per decade) while the
formula stays near 40 MB, and FLINT and the count abort on a failed
allocation.  The cutoff ``2^54`` keeps the count's peak under the peak of
the listing that the same function already allows below its own cutoff:
FLINT ``qfb_reduced_forms`` on the non-fundamental ``D = -(2^44 - 4)`` peaks
at 461 MB, and the whole ``proven`` computation on
``D = -18014398509481979``, the largest fundamental ``D = 1 mod 4`` below
``2^54``, peaks at 391 MB with the count.  ``2^57`` would have
been the edge of the measured range, at about 1.1 GB.

The listing cutoff of a non-fundamental ``D`` stays at ``2^44``
(``imaginary_quadratic_listing_max_abs_discriminant_bits``).  No maximal
order has a non-fundamental discriminant, so no class/unit computation
reaches it; such a ``D`` is listed below ``2^44`` and passed to the formula
from ``2^44`` on, as before.  The formula's contract is a fundamental ``D``,
and it gives no ``h(D)`` for ``D = -(2^44 + 23)`` (conductor 3), so the call
fails closed there, as before; the slow test below checks this.

On ``imaginary_quadratic_d100000000000319`` (``proven``), callgrind counts
363.19G Ir of ``compute_class_unit_transaction`` before the change and 55.29G
after it (-84.8%), and the full output is unchanged on that row, on
``|D| = 18000000000003`` (45 bits) and on ``D = -18014398509481979``.  The
slow test ``t-zeta-quadratic-l1`` (``SILEX_TEST_SLOW``) checks that the count
and the formula give the same ``h(D)`` at ``D = -(2^44 + 3)`` and
``D = -(2^50 + 3)``.
