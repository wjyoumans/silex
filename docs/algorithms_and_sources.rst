Algorithms and Source Policy
============================

Silex is a native C++ library over FLINT.  Its implementation owns its C++
object model, storage, diagnostics, and failure contracts, while several
mathematical algorithms retain documented lineage from established systems.
Lineage is both an engineering aid and an attribution obligation; it does not
create a runtime dependency on those systems.

Source-first rule
-----------------

For production mathematical behavior, apply this order of authority:

#. the current approved Silex requirement;
#. the cited upstream algorithm and behavior;
#. the current Silex behavior and regression tests;
#. the contracts of the FLINT operations used by the implementation.

Do not invent a proof rule, normalization, stopping condition, certification
label, or performance heuristic when these sources disagree or leave a gap.
Resolve a material public choice explicitly and document it.

Lineage map
-----------

``FLINT``
   Supplies the exact arithmetic and storage primitives used by the native
   implementation.  Silex wraps owned FLINT objects with RAII types and uses
   borrowed views at explicit interoperation boundaries.

``PARI/GP 2.17.3``
   Supplies primary source logic for relation completion, analytic
   class-number finishing, HNF and regulator reconstruction, saturation,
   unit-log validation, and supporting ideal, enumeration, and zeta routines.
   The detailed map identifies the relevant files and algorithm families.

``Hecke.jl v0.38.6``
   At commit ``74215ba3d34f296e6f709e415e8007d225524287`` and root tree
   ``5758221d5c6c176b4781dbafe267bb056099d56b``, supplies source or behavioral
   baselines for selected LLL relation search, relation-unit reduction, and
   class/unit workflows.

``Hecke.jl v0.39.19``
   At commit ``122658620f5ac3c8260785c06d9ce7062f037498`` and root tree
   ``b538e80d8c6d6eb783b99b0a217c864458cc7b99``, separately supplies the
   S-class/S-unit source and comparison baseline.

``Silex``
   Owns the final public API, mathematical parent/lifetime model, options and
   result records, diagnostics, and explicit failure behavior.  A source
   translation is not automatically a public-API specification.

The detailed per-family index is
:doc:`reference/algorithms_and_sources`.  See :doc:`legal_and_provenance` for
public source links and the distributed license and attribution files.

Maintaining source fidelity
---------------------------

A mathematically significant change should:

* identify the upstream file, routine or range, and source version where one
  exists;
* state the exact behavior being preserved or intentionally changed;
* preserve normalization, parent identity, exactness, failure, and
  certification contracts unless an approved decision changes them;
* add focused invariant, regression, or differential tests;
* separate correctness evidence from performance evidence; and
* update attribution when newly translated or adapted material is introduced.

Experiments must remain isolated and carry an acceptance or stop criterion.
Benchmark improvement alone is not sufficient to replace established
mathematical behavior.  The contributor workflow is documented in
:doc:`development/source_fidelity`.
