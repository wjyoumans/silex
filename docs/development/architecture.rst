Native C++ Architecture
=======================

Silex uses a deliberately restricted C++20 profile: C++ domain objects expose
the public API, Silex RAII values own FLINT storage, and narrow borrowed views
cross into the FLINT C API.

.. code-block:: text

   Silex public API
       uses
   Silex domain objects and algorithms
       own/use
   Silex FLINT RAII values and borrowed views
       call
   FLINT C API

This is not a C-with-classes layer and not an expression-template algebra
system.  The goal is explicit ownership, visible algorithmic cost, and code
whose mathematical behavior remains reviewable.

Public API and ownership
------------------------

Use Silex domain types, named expensive operations, explicit options, and
explicit status/result records.  Mathematical parents such as ``NumberField``
and ``Order`` are value-like handles backed by shared data so child objects can
keep their parent alive.  Borrowed accessors such as ``parent()`` are local
views, not mathematical identity tokens.

Factories and direct constructors are the ordinary construction style.
Mutating ``define`` methods remain compatibility or scratch-object helpers
where failure-preserving reuse is useful.  ``flint::*Ref`` and
``flint::*ConstRef`` parameters are borrowed wrapper views for explicit input
and output buffers; ``raw_flint_*`` accessors are secondary interoperation
escape hatches.

Silex 0.x makes no source-compatibility or binary-ABI promise.  Consumers must
rebuild against the exact release they use.  Opaque implementation storage
improves ownership boundaries but does not imply ABI stability.

Restricted language profile
---------------------------

Use C++20 namespaces, RAII classes, move-only owners, ``enum class``, small
pure helpers, explicit ``Options``/``Result``/``Status`` records, and visibly
borrowed ranges.  Small templates are acceptable when they remove obvious
duplication and remain easy to instantiate and review.

Avoid by default:

* template metaprogramming and expression-template algebra;
* deep, multiple, or virtual inheritance;
* RTTI, ``dynamic_cast``, and ``typeid``;
* implicit expensive copies or conversions between mathematical domains;
* operator overloads that hide substantial work or normalization;
* hidden global registries or mutable singletons;
* exception-only public APIs; and
* manually owned raw FLINT storage in ordinary module code.

The core target disables exceptions and RTTI by default.  A profile exception
requires a concrete need, an ownership and failure analysis, and focused
tests.  ``tools/check-cpp-profile.py`` is a textual guard, not a substitute for
review.

FLINT ownership
---------------

Owning wrappers initialize in constructors, clear in nonthrowing destructors,
and are non-copyable unless copying has an explicit mathematical meaning.
Moves must be safe and tested.  Borrowed types use visible names such as
``Ref``, ``ConstRef``, ``View``, or ``ConstView``.  Raw handles are limited to
the wrapper layer, direct FLINT call sites, borrowed views, interoperation
bridges, or a justified low-level kernel.

Error contract
--------------

Use ``bool`` for one clear success/failure result, ``std::optional<T>`` for a
single owned value that may be absent, named result records for related
outputs or proof metadata, and ``silex::Status`` when callers need a coarse
failure category.  Recoverable input, domain, resource, or backend failures
must not be reported only through assertions.

Mutating operations publish only after validation and preserve documented
outputs on recoverable failure.  Assertions and debug checks enforce
programmer invariants; they do not replace input validation.  Third-party
exception behavior should be translated to an explicit Silex boundary where
practical.

Optional backends
-----------------

FLINT is the default exact arithmetic and storage engine.  The implemented
optional lattice backends are selected with ``SILEX_WITH_FPLLL`` and
``SILEX_WITH_FLATTER``.  Backend types stay out of public APIs, conversion
overhead is measured, and the default FLINT path remains available.  An
explicit request for an unimplemented backend must fail configuration rather
than create a no-op build.

Installed optional-backend packages are not qualified for 0.1.1.  A source
tree or experimental package must provide the required dependency through
pkg-config, ``CMAKE_PREFIX_PATH``, or the documented backend root variables.

Diagnostics
-----------

Verbose progress, logging, debug checks, profiling, and benchmarking remain
separate mechanisms.  Their compile-time gates, runtime context, facade
requirements, and validation rules are documented in
:doc:`instrumentation`.
