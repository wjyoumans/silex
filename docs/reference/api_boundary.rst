Native API Boundary
===================

The native C++ API is the primary Silex surface.  The public headers also
contain a small amount of compatibility, scratch-object, and FLINT interop
surface.  This page records how to read those entry points.

Ordinary Native API
-------------------

Ordinary C++ callers should use:

- domain objects such as ``NumberField``, ``Element``, ``Order``, ``Ideal``,
  ``FactorBase``, ``Relation``, ``ResidueRing``, ``ResidueField``,
  ``ClassGroupContext``, and ``OrderUnitGroup``;
- factory construction where available, especially
  ``NumberField::by_polynomial``, ``NumberField::quadratic``,
  ``Order::equation_order``, ``Order::from_basis``, and
  ``Order::quadratic_order``;
- named operations such as ``multiply``, ``norm``, ``trace``, ``contains``,
  ``factor``, ``reduce``, ``discrete_log``, and certification/proof accessors;
- owned-return convenience helpers such as ``Lat::basis``,
  ``SparseMat::to_fmpz_mat``, ``Element::to_fmpq_poly``, ``Order::basis``,
  ``OrderElement::coordinates``, ``FiniteAbelianGroup::order``,
  ``ResidueFieldElement::polynomial``, and ``OrderUnitGroup::regulator`` when
  caller code does not need to reuse an output buffer;
- Silex FLINT RAII values such as ``flint::Fmpz``, ``flint::Fmpq``,
  ``flint::FmpzMat``, ``flint::Arb``, and their ``Ref``/``ConstRef`` views as
  explicit input/output buffers.

The ``flint::*Ref`` and ``flint::*ConstRef`` parameters are borrowed views of
caller-owned wrapper objects or FLINT objects.  They are not ownership
transfer, and they are not a request for users to manage raw FLINT lifetimes.
They remain the preferred API shape for hot paths, reusable output storage,
and operations where callers need failure-preserving mutation.
Multi-output proof/audit APIs should use named result structs when they gain
owned-return convenience forms; do not hide several status fields in ad hoc
tuples.  Existing output-parameter forms remain useful for failure-preserving
call sites and compatibility tests.

Error Policy
------------

Public functions use explicit success/failure reporting rather than
exceptions.  Mutating operations normally return ``bool`` and preserve
documented outputs on failure.  Owned-return helpers use ``std::optional`` for
single optional values and named result structs for multi-field records.
``silex::Status`` is reserved for APIs that need to distinguish broad failure
classes.  The core target builds with exceptions disabled by default.

Compatibility and Scratch APIs
------------------------------

Many domain classes still expose mutating ``define`` or ``define_*`` methods.
They are compatibility and scratch-object helpers, not the preferred public
construction style.  They remain useful where tests or algorithms explicitly
check failure preservation, reuse an existing output object, or need a
default-constructed object filled later.

The current audited compatibility constructors are:

- ``NumberField::define_by_polynomial`` and ``NumberField::define_quadratic``;
- ``Order::define`` and ``Order::define_equation_order``;
- ``Element::define`` and ``OrderElement::define``;
- ``Ideal::define`` and ``FractionalIdeal::define``;
- ``PrimeIdeal::define`` and ``PrimeIdealList::define``;
- ``FactorBase::define`` and ``IdealFactorization::define``;
- ``Relation::define`` and ``RelationMatrix::define``;
- ``ResidueRing::define`` and ``ResidueElement::define``;
- ``ResidueFieldElement::define`` and ``ResidueFieldQuotientLog::define``;
- ``EmbeddingContext::define``, ``FieldHom::define``,
  ``OrderHom::define``, ``FieldAutomorphism::define``,
  ``FactoredElement::define``, ``CompactElement::define``,
  ``ClassGroupContext::define``, and ``OrderUnitGroup::define``.

Do not add new public mutating construction APIs as the default shape for new
modules.  Prefer factories for parent/context objects and direct constructors
for simple parented values.  Keep existing ``define`` methods documented as
compatibility until callers have a better replacement and the failure behavior
is covered elsewhere.

Ordinary examples and module setup paths use factories or direct constructors.
Remaining direct
``.define(...)`` call sites in the repository are compatibility/failure tests,
while wrapper-internal ``define`` calls remain part of low-level FLINT RAII
construction.

FLINT Wrapper Views
-------------------

The wrapper layer under ``silex/flint/`` is public support infrastructure for
the current native C++ API.  It gives Silex RAII ownership of FLINT storage and
provides named borrowed views for input/output parameters.

This layer is intentionally low-level.  Its ``raw()`` and ``data()`` accessors
are escape hatches for direct FLINT calls, interop code, parity tests, and
small kernels.  They should be used locally at FLINT call sites, not stored as
long-lived application state.

When writing ordinary Silex code or examples, prefer wrapper objects and Silex
domain operations.  Add narrow wrapper helper functions when repeated direct
FLINT calls make module code hard to read or hard to audit.

Raw FLINT Escape Hatches
------------------------

The current non-wrapper public headers expose only a few true raw FLINT
escape hatches:

- ``NumberField::flint_field_ref`` and ``NumberField::raw_flint_field``;
- ``Element::flint_element_ref`` and ``Element::raw_flint_element``;
- ``lat::ShortVectorCallback`` and ``Lat::raw_basis``.

These are for bridge code, parity tests, and low-level kernels that must call
FLINT directly.  They are not ordinary user-facing algebra APIs and should not
be used for mathematical identity, ownership, or persistent parent storage.

Some raw or direct-storage helpers are private implementation hooks rather
than public API, including ``Lat::set_basis_direct``,
``fmpz_smat::HnfContext::replace_hnf``, and ``Ideal::set_hnf_direct``.  Keep
them private unless a concrete benchmarked kernel or interop boundary needs
them.

Future Binding Surface
----------------------

Future Python/Sage bindings should expose the native domain objects,
factories, named operations, and explicit result/certification data.  They
should not expose raw FLINT pointers, wrapper ``raw()`` accessors, or
compatibility ``define`` methods as the primary construction style.

The reference docs should therefore describe compatibility and raw
interop APIs as available but secondary.  They should not be used in examples
unless the example is specifically about FLINT interop or failure-preserving
scratch-object mutation.

Compiled examples should prefer owned-return helpers for ordinary scalar,
matrix, polynomial, and proof-record queries where those helpers exist.  Keep
explicit output buffers in examples that demonstrate reusable storage,
matrix/log kernels, proof-driver workflows, or APIs whose failure-preserving
output semantics are the point.
