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

Each factory returns an undefined ``NumberField`` on invalid input.  Call
``is_defined()`` before building dependent objects when input validity is not
already known.

The mutating ``define_by_polynomial`` and ``define_quadratic`` methods remain
available as compatibility and scratch-object helpers.  They are tested for
failure-preservation behavior, but they are not the preferred construction
style for ordinary public code.

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
``Element::to_fmpq_poly`` is the owned-return convenience form of
``get_fmpq_poly`` and returns ``std::nullopt`` when the element is undefined.

``silex::Signature`` is a value object for real and complex place counts.
``EmbeddingContext`` stores certified complex roots and a ``NumberField``
parent handle, then refines the roots to requested precision.  Archimedean
helpers compute absolute values, logarithmic embeddings, and Minkowski
embeddings from an
``EmbeddingContext`` and ``Element``.

The field-level unit helpers in ``silex/unit.hpp`` cover Dirichlet unit rank,
supported roots of unity, lower regulator bounds, supplied-unit logarithmic
matrices and regulators, supplied-unit independence checks, and exact
continued-fraction fundamental units for the explicit and canonical
polynomial-defined ``x^2-d`` real-quadratic routes.  See
:doc:`../support_matrix` for the presentations qualified for 0.1.0.

Homomorphisms and Automorphisms
-------------------------------

``silex::FieldHom`` stores value-like domain and codomain field handles, a
certified generator image, and explicit application methods.  Its field
accessors return borrowed pointers to those stored handles.  Construction
checks the defining polynomial exactly before accepting the generator image.
The identity and isomorphism predicates are data queries; they do not imply a
full automorphism-group enumeration.

``silex::OrderHom`` combines a field homomorphism with source and target
``Order`` handles.  It certifies that source-basis images land in the target
order and stores the integer image matrix used by order-homomorphism
application.

``silex::FieldAutomorphism`` is the current finite automorphism wrapper.  It
supports identity and quadratic conjugation constructors, optional certified
endomorphism storage through ``FieldHom``, application to elements, and
homomorphism extraction.  General automorphism lists/groups and cyclotomic
metadata remain future source-backed work.

Low-level FLINT interop remains available through ``flint_field_ref()`` and
``raw_flint_field()`` for bridge code and parity tests.  Ordinary user code
should prefer the domain operations above and should not inspect the
underlying ``nf_t`` unless it is deliberately crossing into FLINT-level code.

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
           logs, embeddings, theta, silex::LogEmbeddingMode::plain, 128);

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
``examples/log_unit_lattice.cpp`` for embeddings and log unit lattices,
``examples/field_order_ideal_basics.cpp`` for the field/order/ideal path, and
``examples/real_quadratic_order_units.cpp`` for an exact continued-fraction
fundamental unit on a release-qualified maximal real-quadratic example through
the order-unit layer.
