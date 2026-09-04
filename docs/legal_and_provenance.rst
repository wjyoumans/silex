Legal and Source Provenance
===========================

Silex is distributed under the GNU General Public License, version 3 or, at
your option, any later version (``GPL-3.0-or-later``).  Redistribution must
retain the :download:`license <../LICENSE>` and
:download:`third-party notices <../THIRD_PARTY_NOTICES.md>`.  The
:download:`repository notice <../NOTICE.md>` summarizes attribution and
:download:`citation metadata <../CITATION.cff>` supplies a machine-readable
citation record.  This page is a technical provenance index; the distributed
legal files are authoritative for license terms.

Development process disclosure
------------------------------

Silex and its companion repositories, `Silex Bench
<https://github.com/wjyoumans/silex-bench>`_ and `Silex Devtools
<https://github.com/wjyoumans/silex-devtools>`_, were built almost entirely
with OpenAI Codex, initially using GPT-5.5 and later GPT-5.6, under the
direction and review of William Youmans.  This development-process disclosure
is separate from the licenses, copyright notices, and upstream source
attribution recorded below.

Upstream source anchors
-----------------------

Several mathematical routines preserve behavior or algorithm lineage from
PARI/GP and Hecke.jl.  These projects are source and comparison baselines, not
runtime dependencies of the native Silex library.

``PARI/GP 2.17.3``
   Primary source for the class-group, unit-group, ideal, relation-search,
   enumeration, zeta, and S-unit families identified in
   :doc:`reference/algorithms_and_sources`.  The versioned release is available
   from the `PARI/GP distribution site
   <https://pari.math.u-bordeaux.fr/download.html>`_.

`Hecke.jl v0.38.6 <https://github.com/thofma/Hecke.jl/tree/v0.38.6>`_
   Source or behavioral baseline for the selected order, ideal, relation,
   class/unit, and compact-element families identified in the source map.
   The peeled tag commit is
   `74215ba3d34f296e6f709e415e8007d225524287
   <https://github.com/thofma/Hecke.jl/commit/74215ba3d34f296e6f709e415e8007d225524287>`_
   and its root Git tree is
   ``5758221d5c6c176b4781dbafe267bb056099d56b``.

`Hecke.jl v0.39.19 <https://github.com/thofma/Hecke.jl/tree/v0.39.19>`_
   Separate source and comparison baseline for the S-class/S-unit publication
   family.  The peeled tag commit is
   `122658620f5ac3c8260785c06d9ce7062f037498
   <https://github.com/thofma/Hecke.jl/commit/122658620f5ac3c8260785c06d9ce7062f037498>`_
   and its root Git tree is
   ``b538e80d8c6d6eb783b99b0a217c864458cc7b99``.

``FLINT``
   Supplies the exact arithmetic, storage, HNF/SNF, polynomial, number-field,
   and modular-linear-algebra contracts used by Silex.  The supported and
   qualified versions are recorded in :doc:`support_matrix`.

Historical Silex import
-----------------------

The 0.1.1 foundation was imported from a historical native-library tree.  SHA
``7fdce40f30abbb08e024253347efba0d4e5fcc5a`` is retained as a **private
archival anchor** for that import audit.  The object is not part of the public
repository history and is not a public reproducibility reference.  Publicly
reproducible behavior is defined by the selected Silex release, its tests, and
the upstream anchors above.

History rewriting or repository reorganization does not remove upstream
authorship, algorithm lineage, or license obligations.  Changes to adapted
algorithms must preserve the relevant attribution and update the public source
map when an anchor changes.
