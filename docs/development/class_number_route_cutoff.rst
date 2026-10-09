Class number route cutoff
=========================

A ``grh`` class/unit request on an imaginary quadratic field takes the class
number from one of two routes, chosen by the bit length of the discriminant
(``grh_imaginary_quadratic_lattice_min_abs_discriminant_bits`` in
``src/class_group/class_group.cpp``):

- ``|D| < 2^44``: the exact class number, counted from the reduced forms of
  discriminant ``D`` (FLINT ``qfb_reduced_forms``), whose time and memory grow
  like ``sqrt(|D|)``.
- ``|D| >= 2^44``: the GRH relation lattice and the analytic index-one test,
  which never lists forms.  At these sizes the test uses a Belabas--Friedman
  ``hR``, because ``2^44`` is also the cap of the unconditional quadratic
  ``L(1, chi_D)`` route (``quadratic_route_max_abs_discriminant_bits`` in
  ``src/zeta/zeta.cpp``).  Below the cap the test would use that
  unconditional value, whose approximate functional equation costs more
  than listing the forms.

This page records the measurement that set the cutoff.  The requirement was
that the combined route be no slower than the reduced-form listing alone at
any measured ``|D|`` from ``1e7`` to ``1e14``, with at least two
discriminants per decade.

Setup
-----

Each field is ``x^2 + n`` with ``n`` squarefree and ``n = 3 mod 4``, so
``D = -n``.  Each run is one ``silex-class-unit-instance --coeffs n,0,1
--mode grh`` process.  Wall time is the median over the listed runs.  Peak
memory is the process's maximum resident set size; the process uses about
21 MB before any class-group work.  Three separate Release trees were built
with identical settings (GCC 16.2.1, FLINT 3.6.0, AMD Ryzen AI 9 HX PRO
370): "main" lists forms at every ``|D|``, "lattice only" uses the lattice
at every ``|D|``, and "combined" uses the cutoff.  The machine was busy with
other work (load average 5 to 12 on 24 cores, ``machine-idle.sh`` reported
BUSY at start and end), so these are busy-machine timings, not a controlled
benchmark.  ``h`` and the class-group invariants agreed in all three trees
in every row; at ``n = 10^14 + 319`` ``h = 15815072`` agrees with PARI/GP
2.17 ``quadclassunit``.

Results
-------

.. list-table::
   :header-rows: 1

   * - ``n``
     - bits of ``|D|``
     - ``h``
     - runs
     - main: s, MB
     - lattice only: s, MB
     - combined: s, MB
   * - 10000003
     - 24
     - 706
     - 3
     - 0.011, 21
     - 0.183, 21
     - 0.011, 21
   * - 30000003
     - 25
     - 1380
     - 3
     - 0.011, 21
     - 0.214, 21
     - 0.012, 21
   * - 100000003
     - 27
     - 1702
     - 3
     - 0.016, 22
     - 0.661, 21
     - 0.017, 22
   * - 300000003
     - 29
     - 3072
     - 3
     - 0.031, 23
     - 0.823, 22
     - 0.030, 23
   * - 1000000003
     - 30
     - 3680
     - 1
     - 2.100, 25
     - 2.085, 22
     - 2.088, 25
   * - 3000000003
     - 32
     - 9728
     - 3
     - 0.080, 28
     - 2.232, 22
     - 0.078, 28
   * - 10000000003
     - 34
     - 10538
     - 1
     - 5.181, 32
     - 5.203, 22
     - 5.506, 32
   * - 30000000003
     - 35
     - 33216
     - 3
     - 0.084, 40
     - 5.701, 22
     - 0.084, 40
   * - 100000000003
     - 37
     - 31057
     - 1
     - 10.677, 55
     - 10.696, 24
     - 10.815, 55
   * - 300000000007
     - 39
     - 173988
     - 3
     - 0.245, 82
     - 16.946, 23
     - 0.246, 82
   * - 1000000000003
     - 40
     - 124568
     - 1
     - 26.825, 125
     - 27.516, 24
     - 28.299, 126
   * - 3000000000003
     - 42
     - 412512
     - 3
     - 0.383, 207
     - 34.430, 23
     - 0.373, 207
   * - 10000000000003
     - 44
     - 497056
     - 1
     - 77.044, 355
     - 78.642, 25
     - 78.804, 355
   * - 15000000000007
     - 44
     - 1331780
     - 1
     - 91.277, 446
     - 92.895, 24
     - 86.812, 446
   * - 18000000000003
     - 45
     - 888096
     - 3
     - 0.854, 473
     - 0.329, 24
     - 0.321, 24
   * - 20000000000003
     - 45
     - 1081476
     - 3
     - 0.932, 501
     - 0.406, 24
     - 0.396, 24
   * - 30000000000003
     - 45
     - 1290536
     - 3
     - 1.031, 607
     - 0.378, 24
     - 0.397, 23
   * - 50000000000003
     - 46
     - 1469496
     - 3
     - 1.820, 774
     - 1.017, 26
     - 0.982, 26
   * - 100000000000003
     - 47
     - 1425472
     - 3
     - 2.384, 1071
     - 1.345, 27
     - 1.359, 27
   * - 100000000000319
     - 47
     - 15815072
     - 3
     - 3.211, 1396
     - 0.118, 22
     - 0.118, 22

The combined route equals main below ``2^44`` and is faster above it.  The
lattice-only route is slower than main at every measured ``|D| < 2^44`` and
faster at every measured ``|D| >= 2^44``, so ``2^44`` is the cutoff.  The
fields ``n = 10^k + 3`` (``k = 9`` to ``13``) and ``n = 15000000000007``
take seconds to minutes in all three trees, by similar amounts in both
routes; that cost does not come from listing forms and is not studied here.
A second, interleaved run of main
and combined at ``n = 10^10 + 3`` and ``n = 10^12 + 3`` (three runs each)
gave medians of 5.07 s and 5.09 s, and 25.87 s and 25.92 s.
