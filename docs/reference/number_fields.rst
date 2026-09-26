Number Fields, Elements, and Embeddings
=======================================

The native number-field foundation lives in ``silex/number_field.hpp``,
``silex/element.hpp``, ``silex/signature.hpp``, ``silex/embedding.hpp``,
``silex/archimedean.hpp``, and ``silex/unit.hpp``.

Object Model
------------

``silex::NumberField`` owns the FLINT ``nf_t`` context and is the parent for
``Element`` and ``EmbeddingContext`` objects.  It is a cheap value-like handle
to shared field data; copying a field handle preserves mathematical parent
identity, while ``set`` remains available for code that needs an explicit
scratch-object assignment.

Construction
------------

User code should construct fields with factories:

``NumberField::by_polynomial(flint::FmpqPolyConstRef)``
    Construct a generic number field from a rational polynomial.

``NumberField::by_polynomial(flint::FmpzPolyConstRef)``
    Construct a generic number field from an integral polynomial.

``NumberField::quadratic(flint::FmpzConstRef)``
    Construct the source-backed quadratic field ``Q(sqrt(d))`` when the
    radicand is valid.

A defining polynomial must have positive degree and be irreducible over
``Q``.  It may be nonmonic and may have rational coefficients; the stored
polynomial and generator are the ones supplied.  Construction checks
irreducibility exactly, by factoring the primitive integral multiple of the
polynomial with FLINT ``fmpz_poly_factor``, so reducible input, including
input that is not squarefree such as ``x^3``, is rejected rather than
producing a ring with zero divisors.  This is a one-time cost at
construction.  ``NumberField::quadratic`` accepts only squarefree nonsquare
radicands ``d``, for which ``x^2 - d`` is irreducible.

Each factory returns an undefined ``NumberField`` on invalid input, including
a constant or reducible polynomial.  Call ``is_defined()`` before building
dependent objects when input validity is not already known.

The mutating ``define_by_polynomial`` and ``define_quadratic`` methods remain
available as compatibility and scratch-object helpers.  They apply the same
validity rules, return ``false`` on invalid input, and leave the object
unchanged, so a previously defined field stays defined and an undefined field
stays undefined.  They are not the preferred construction style for ordinary
public code.

.. code-block:: cpp

   #include <silex/flint/fmpz.hpp>
   #include <silex/number_field.hpp>

   silex::flint::Fmpz d;
   silex::flint::fmpz_set_si(silex::flint::FmpzRef(d), 5);

   silex::NumberField K =
           silex::NumberField::quadratic(silex::flint::FmpzConstRef(d));
   if (!K.is_defined()) {
       return false;
   }

``silex::Element`` owns one field element and stores a ``NumberField`` parent
handle, so returned elements keep their parent field alive.  Element
operations are explicit named methods such as ``add``, ``multiply``,
``invert``, ``pow_fmpz``,
``trace``, ``norm``, ``conjugate``, ``is_square``, and ``is_power``.  Methods
return ``false`` on domain mismatch, undefined parents, unsupported exact root
cases, or invalid inputs, and leave documented outputs unchanged on failure.
``invert`` fails on zero; because defining polynomials are irreducible, every
other element is invertible.  ``invert`` still verifies exactly that its input
is not a zero divisor and fails, leaving the output unchanged, if it is.
``Element::to_fmpq_poly`` is the owned-return convenience form of
``get_fmpq_poly`` and returns ``std::nullopt`` when the element is undefined.

``silex::Signature`` is a value object for real and complex place counts.
``EmbeddingContext`` stores certified complex roots and a ``NumberField``
parent handle, then refines the roots to requested precision.  Its
``signature()`` is computed from the defining polynomial when the context is
defined, so it is available before the first ``refine``.  Roots are stored with
the real places first, then each complex place as its root with positive
imaginary part followed by the conjugate.  Place indices are fixed once roots
are set: each refinement matches the new root balls to the previous ones by
overlap, isolating at higher precision until the match is one-to-one.
``refine`` fails, leaving the context unchanged, only for a non-squarefree
defining polynomial or when the working precision would overflow.  Archimedean
helpers compute absolute values, logarithmic embeddings, and Minkowski
embeddings from an
``EmbeddingContext`` and ``Element``.

The field-level unit helpers in ``silex/unit.hpp`` cover Dirichlet unit rank,
supported roots of unity, lower regulator bounds, supplied-unit logarithmic
matrices and regulators, supplied-unit independence checks, and exact
continued-fraction fundamental units for the explicit and canonical
polynomial-defined ``x^2-d`` real-quadratic routes.  See
:doc:`../support_matrix` for the presentations qualified for 0.1.1.

Homomorphisms and Automorphisms
-------------------------------

``silex::FieldHom`` stores value-like domain and codomain field handles, a
certified generator image, and explicit application methods.  Its field
accessors return borrowed pointers to those stored handles.  Parent-only
construction leaves the generator image unset; ``set_generator_image`` checks
the defining polynomial exactly before accepting the image.
The identity and isomorphism predicates are data queries; they do not imply a
full automorphism-group enumeration.

``silex::OrderHom`` combines a field homomorphism with source and target
``Order`` handles.  It certifies that source-basis images land in the target
order and stores the integer image matrix.  Its ``apply`` method applies the
certified field homomorphism.

``silex::FieldAutomorphism`` is the current finite automorphism wrapper.  It
supports identity and quadratic conjugation constructors, optional certified
endomorphism storage through ``FieldHom``, application to elements, and
homomorphism extraction.  General automorphism lists/groups and cyclotomic
metadata remain future source-backed work.

The maps retain their field/order handles, keeping parent data alive when the
original local handles go out of scope.  Pointers returned by their parent
accessors are borrowed views of the map's stored handles; do not retain them
across mutation or destruction of that map.  Compare field parent identity
with ``NumberField::has_same_data``, rather than accessor pointer addresses.
Two separate constructions from the same polynomial create distinct parents;
copying a field handle preserves its parent identity.

Low-level FLINT interop remains available through ``flint_field_ref()`` and
``raw_flint_field()`` for bridge code and parity tests.  Ordinary user code
should prefer the domain operations above and should not inspect the
underlying ``nf_t`` unless it is deliberately crossing into FLINT-level code.
On an undefined field both accessors expose a null ``nf_struct`` pointer
(``flint_field_ref().raw()`` and ``raw_flint_field()`` return ``nullptr``).
Modifying the ``nf_t`` behind a defined field bypasses the construction
checks and is not supported.

.. _field-maps-walkthrough:

Compiled Field-Map Walkthrough
------------------------------

The program :download:`field_maps.cpp <../../examples/field_maps.cpp>` follows
the existing cases in ``test/t-hom.cpp`` and ``test/t-aut.cpp``.  Every fallible
operation is checked before its output is used, with a nonzero exit on failure.
Build and run it with:

.. code-block:: console

   cmake --preset default
   cmake --build --preset default --target example-field-maps
   ./build/default/examples/example-field-maps
   ctest --preset default -R '^silex-example-field-maps$'

First construct two distinct parents for ``Q(sqrt(2))``.  A copied handle
shares the first parent's data:

.. literalinclude:: ../../examples/field_maps.cpp
   :language: cpp
   :start-after: // field-maps-parents-begin
   :end-before: // field-maps-parents-end
   :dedent: 4

Certify the image of the first generator in the second field.  The resulting
map is an isomorphism, but ``is_identity()`` is false because the parents are
distinct.  Rejecting the nonroot ``1`` preserves the previously certified
image.  Application takes an input in the domain and an output in the
codomain.  The same-parent map that fixes the generator is the identity.
The private ``set_linear`` helper used in these excerpts sets an element from
an exact polynomial with the given constant and linear coefficients; ``fail``
prints a diagnostic and returns a nonzero exit code.

.. literalinclude:: ../../examples/field_maps.cpp
   :language: cpp
   :start-after: // field-maps-hom-begin
   :end-before: // field-maps-hom-end
   :dedent: 4

The automorphism wrapper supplies identity and quadratic conjugation.  Here
conjugation sends ``theta_K`` to ``-theta_K`` and ``3 + 4*theta_K`` to
``3 - 4*theta_K``.  Applying it twice recovers the original element; extracting
its ``FieldHom`` gives the same image:

.. literalinclude:: ../../examples/field_maps.cpp
   :language: cpp
   :start-after: // field-maps-aut-begin
   :end-before: // field-maps-aut-end
   :dedent: 4

Finally, certify that conjugation maps the equation order to itself.  Each
row of the exported integer matrix gives the target-basis coordinates of a
source-basis image.  For the basis ``(1, theta_K)``, this is ``diag(1, -1)``.
The owned ``FmpzMat`` overload resizes its output on success:

.. literalinclude:: ../../examples/field_maps.cpp
   :language: cpp
   :start-after: // field-maps-order-begin
   :end-before: // field-maps-order-end
   :dedent: 4

Minimal Example
---------------

.. code-block:: cpp

   #include <silex/archimedean.hpp>
   #include <silex/element.hpp>
   #include <silex/embedding.hpp>
   #include <silex/flint/arb_vec.hpp>
   #include <silex/flint/fmpz.hpp>
   #include <silex/number_field.hpp>

   silex::flint::Fmpz radicand;
   silex::flint::fmpz_set_si(silex::flint::FmpzRef(radicand), 5);

   silex::NumberField K =
           silex::NumberField::quadratic(silex::flint::FmpzConstRef(radicand));
   const bool field_ok = K.is_defined();

   silex::Element theta(K);
   const bool gen_ok = theta.gen();

   silex::EmbeddingContext embeddings(K);
   const bool refined = embeddings.refine(128);

   silex::flint::ArbVec logs(K.degree());
   const bool logs_ok = silex::logarithmic_embedding(
           silex::flint::ArbVecRef(logs), embeddings, theta,
           silex::LogEmbeddingMode::plain, 128);

Field/Order/Ideal Walkthrough
-----------------------------

The example ``examples/field_order_ideal_basics.cpp`` shows the intended
C++-first path from a field to an order, an integral ideal, and a fractional
ideal:

- define ``Q(theta)`` by ``theta^2 - 2`` with
  ``NumberField::quadratic``;
- build the equation order with ``Order::equation_order``;
- convert the generator to an ``OrderElement`` and construct the principal
  integral ideal ``(theta)``;
- construct the fractional principal ideal ``(theta/2)`` from an ambient
  ``Element``;
- query discriminants, containment, and integral/fractional norms through
  explicit output parameters.

The example uses Silex domain objects and Silex FLINT RAII values at the call
site; it does not require user code to own raw FLINT handles.
The integral/fractional ideal distinction is intentional: integral principal
ideals take ``OrderElement`` generators, while fractional principal ideals take
ambient ``Element`` generators so denominator clearing remains explicit in the
``FractionalIdeal`` representation.

Exact Element Arithmetic
------------------------

The example ``examples/element_arithmetic_basics.cpp`` shows exact arithmetic
over ``Q(sqrt(5))`` using the native ``Element`` API:

- construct the generator ``theta`` and an element ``alpha = theta + 3`` with
  ``Element::gen`` and ``Element::add_si``;
- query exact ``trace`` and ``norm`` values through ``flint::Fmpq`` output
  parameters;
- multiply elements with ``Element::multiply``;
- compute inverses with ``Element::invert`` and check products with
  ``Element::equal_si``;
- compute signed integer powers with ``Element::pow_fmpz``.

The example intentionally uses named methods rather than operator overloads so
the cost and failure behavior of exact algebraic operations remain visible.

Benchmark Coverage
------------------

There are no standalone Google Benchmark targets for many field-layer
components.  Performance for these pieces is currently covered through the
benchmark targets that consume them, especially ``nf_fac_elt`` for
element/root helpers, ``nf_ord`` for order construction, and ``nf_clgp`` for
unit, embedding, zeta, and class/unit proof-driver workflows.

Do not add standalone benchmark targets for these modules unless a current
higher-level benchmark identifies a specific bottleneck.

Source Lineage
--------------

Field and element behavior retains established Silex contracts and uses the
cited PARI quadratic/generic field sources and FLINT ``nf`` primitives.
Signature, embedding, archimedean, unit, homomorphism, and automorphism layers
follow the same source-first and failure-preserving policy.  See
:doc:`algorithms_and_sources` for the public file-level map.

Examples
--------

See ``examples/element_arithmetic_basics.cpp`` for exact element arithmetic,
``examples/field_maps.cpp`` for certified field and order maps,
``examples/log_unit_lattice.cpp`` for embeddings and log unit lattices,
``examples/field_order_ideal_basics.cpp`` for the field/order/ideal path, and
``examples/real_quadratic_order_units.cpp`` for an exact continued-fraction
fundamental unit on a release-qualified maximal real-quadratic example through
the order-unit layer.
