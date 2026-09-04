Lat Benchmarks
==============

This directory contains lattice benchmark kernels through the native C++ API.

Current rows:

- ``b-silex-lat_hnf``
- ``b-silex-lat_lll``
- ``b-silex-lat_intersection``

The LLL target includes the established 4-by-4 kernel, a deterministic synthetic
degree-14 HNF, and a captured native degree-14 ideal HNF.  The captured row
reports first, maximum, and summed squared row norms plus log2 orthogonality
defect.  Smaller quality counters are better.

When configured with ``SILEX_WITH_FPLLL=ON`` and/or
``SILEX_WITH_FLATTER=ON``, the captured row also compares fplll LLL, bounded
one-tour fplll BKZ, and the flatter RHF-1.02 pipeline.  These functions remain
private benchmark plumbing and do not select a class/unit route backend.  For
the small captured matrix, the installed flatter pipeline delegates to its
fplll BKZ implementation.

These are fixed-size smoke/regression kernels.  They are not broad performance
claims and should be compared with a full input trajectory and complete route
before making optimization decisions.
