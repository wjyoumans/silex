0.1.0 Support Matrix
====================

This page is the authoritative support and limitations contract for the tagged
Silex 0.1.0 release.  "Qualified" below means that the exact configuration or
mathematical presentation was part of the release-candidate validation.  It
does not turn a bounded algorithm into a total one, and callers must still
handle ``false`` returns without reading unpublished output.

Qualified native stack
----------------------

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Component
     - Qualified for 0.1.0
     - Boundary
   * - Operating system and architecture
     - Linux x86_64
     - Other operating systems and architectures are unqualified.
   * - C++ compiler
     - GCC 16 and Clang 22
     - C++20 is required.  Other compiler families and major versions are
       unqualified.
   * - FLINT
     - 3.6
     - 3.0.0 is the declared package minimum, not a claim that every 3.x
       release has been qualified.
   * - CMake
     - 4.4
     - 3.20 is the declared project minimum, not a qualified version range.
   * - Package configuration
     - Default FLINT-backed build
     - The installed package exports ``Silex::silex`` and discovers FLINT
       through ``pkg-config``.
   * - Optional backends
     - None
     - fplll and flatter remain source-tree development options.  Their
       installed-package configurations are not supported for 0.1.0.

The 0.1.0 qualification used GCC 16.1.1, Clang 22.1.8, FLINT 3.6.0, and CMake
4.4.2.  The table deliberately states major/minor support while
:doc:`releases/0.1.0` records the exact patch versions, avoiding any
implication that untested intervening versions were qualified.

Public mathematical boundary
----------------------------

The default package provides the native C++20 domain types and exact,
FLINT-backed operations documented in this manual for absolute simple number
fields over ``QQ``.  This includes fields and elements, orders, ideals, local
and residue arithmetic, lattices, finite abelian groups, factor bases,
relations, compact elements, class/unit infrastructure, and S-class/S-unit
publication.  Exact data and a completion guarantee are different things:
resource-bounded searches and proof drivers may fail closed even when their
inputs are mathematically valid.

The following contracts apply across the public surface:

* long-lived domain objects keep mathematical parents by value-like handles;
  pointer-returning parent accessors are borrowed views, not identity tokens;
* canonical HNF/SNF data, normalized ideal numerator/denominator storage, and
  documented coordinate orders are preserved;
* mutating operations publish only after validation and preserve documented
  outputs on failure; the paired class/unit and S-class/S-unit entry points
  are failure-atomic across both output objects; and
* certification is explicit data.  Successful construction is not an
  unconditional proof unless the result and required component records say
  ``proven`` and ``verified``.

Paired class/unit qualification
-------------------------------

``OrderUnitGroup::compute_with_class_group`` accepts only explicit ``proven``
or ``grh`` requests and requires a maximal absolute order.  The following
presentations are the conservative 0.1.0 completion matrix under the shipped
test resource settings:

.. list-table::
   :header-rows: 1
   :widths: 22 36 42

   * - Request
     - Qualified presentations
     - Published contract
   * - ``proven``
     - ``QQ``; maximal orders for ``x^2 - 2``, ``x^2 - 5``,
       ``x^2 - 210``, ``x^2 + 1``, ``x^2 + 3``, ``x^2 - x + 1``,
       ``x^2 + 14``, and ``x^2 + 47``
     - Both coarse labels are ``proven`` and the route-specific required
       class, unit, regulator, generation, and saturation records are
       verified.
   * - ``grh``
     - ``QQ``; maximal orders for ``x^2 + 3``, ``x^2 + 34``,
       ``x^2 + 46``, ``x^2 + 47``, ``x^2 + 66``, and ``x^2 + 185``
     - Both coarse labels remain ``grh``.  The result is conditional and
       provisional; it is not an unconditional component proof.

Equivalent-looking polynomials are not automatically covered: the order must
be the same maximal order established by the public construction path.  The
two discriminant ``-3`` presentations above are both regression fixtures for
that parent/presentation boundary.

Higher-degree paired computation and unlisted quadratic presentations remain
available as fail-closed best effort, but they have no 0.1.0 completion or
compatibility guarantee.  Higher-degree native fixtures remain correctness
regressions for the implemented routes; they do not establish blanket cubic,
quartic, quintic, or general-degree support.  Nonmaximal paired orders and
paired ``unknown`` or ``heuristic`` requests are rejected without replacing
existing outputs.

GRH work-in-progress status
---------------------------

GRH certification semantics and coverage are work in progress.  A successful
``grh`` transaction preserves the requested public
``CertificationMode::grh`` label for both objects.  Callers must treat it as
conditional and provisional, inspect the component metadata they need, and
must not reinterpret it as ``proven``.  The degree-one route may retain exact
proof receipts internally while still publishing the caller's requested
``grh`` coarse label.

Other qualified group routes
----------------------------

The standalone ``OrderUnitGroup::compute`` route is release-qualified for its
explicit regression presentations: ``QQ``; maximal real-quadratic orders for
``x^2 - 2``, ``x^2 - 5``, and ``x^2 - 17345``; the maximal order for
``x^2 + 1``; and the equation and maximal orders exercised for ``x^2 + 3``.
Other rank-zero and quadratic inputs may use the same exact source-backed
implementation, but are not blanket completion guarantees for 0.1.0.

S-class and S-unit publication is supported for a maximal absolute order from
already ``proven`` ordinary class- and unit-group inputs for that same order,
with a validated selected-prime list.  It does not independently broaden the
paired class/unit completion matrix.

Unsupported or deferred
-----------------------

The following are outside the 0.1.0 support contract:

* source compatibility across 0.x versions and binary ABI compatibility;
* relative fields, ray class groups, class fields, and general Galois groups;
* a C API, Python or Sage bindings, runtime brokers, or asynchronous APIs;
* installed packages with optional fplll or flatter backends;
* general higher-degree or unlisted quadratic paired completion;
* unconditional interpretation of ``grh`` results;
* cross-engine performance targets or multi-engine campaigns; and
* production-readiness or long-term-support claims.

Consumers should rebuild against the exact Silex version they use.  Public
interfaces and object layouts may change throughout 0.x.
