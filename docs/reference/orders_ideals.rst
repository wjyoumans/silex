Orders and Ideals
=================

The order and ideal layer covers orders, order elements, integral ideals, and
fractional ideals.
Prime ideals, factor bases, relations, residue rings, and residue fields are
documented separately in :doc:`local_algebra`.

Orders and Order Elements
-------------------------

``silex::Order`` is a value-like handle to shared order data and stores its
``NumberField`` parent handle by value.  It can construct equation orders,
explicit-basis orders, direct quadratic/maximal orders, generic p-maximal
overorders, and global maximal orders through the source-backed
maximalization path.  It exposes explicit methods for
coordinates, containment, trace matrices, discriminants, multiplication
tables, order indices, p-radicals, p-maximal overorders, and maximal orders.

User code should construct orders with factories:

``Order::equation_order(field)``
    Construct the equation order of a defined ``NumberField``.

``Order::from_basis(field, basis)``
    Construct an order from a rational basis matrix when the basis is valid.

``Order::quadratic_order(field, conductor)``
    Construct a quadratic order with conductor metadata when supported.

Each factory returns an undefined ``Order`` on invalid input.  Mutating
``define`` and ``define_equation_order`` remain compatibility and
scratch-object helpers; ordinary public examples and tests use factories
except when mutation/failure behavior is what the test covers.

``Order::basis()`` and ``Order::coordinates(element)`` are owned-return
convenience forms of the output-buffer APIs.  They return ``std::nullopt`` on
the same domain and definition failures as the buffer forms.

``silex::OrderElement`` owns one order element and stores an ``Order`` parent
handle.  It is the integral element type accepted by integral ideals, residue
rings, and prime-ideal residue reduction.

An ambient ``Element`` can be converted to an ``OrderElement`` only when it
lies in the order; failed conversions leave the target object unchanged.
``OrderElement::coordinates()`` returns an owned coordinate row, or
``std::nullopt`` for an undefined element.

Ideals
------

``silex::Ideal`` stores an integral ideal as validated row-HNF order
coordinates and keeps its ``Order`` parent alive with a value-like handle.  It
supports unit and principal construction, HNF import/export, containment,
norm, equality, sum, intersection, multiplication, coprimality, add-to-one
witnesses, and multiplier rings.

``silex::FractionalIdeal`` stores an integral numerator ideal and a positive
integer denominator and keeps its ``Order`` parent alive with a value-like
handle.  It supports principal fractional ideals from ambient field elements,
normalization, containment, norm, equality, sum, intersection, multiplication,
powers, colon ideals, and certified inverses.

Integral ideals consume ``OrderElement`` values because they live inside the
order lattice.  Fractional principal ideals consume ambient ``Element`` values
because denominators are part of the fractional-ideal representation.  The
stored representation is normalized as an integral numerator over a positive
integer denominator, and ``norm`` is a logically const cached query.

Use the integral and fractional APIs according to the mathematical object
being represented:

- use ``Ideal::one`` for the unit ideal and ``Ideal::set_principal`` with an
  ``OrderElement`` for an integral principal ideal;
- use ``Ideal::set_hnf`` only when the caller already has order-coordinate HNF
  data; invalid shape, rank, or ideal-closure checks fail without replacing
  the current ideal;
- use ``FractionalIdeal::set_integral`` or
  ``FractionalIdeal::set_integral_den`` when an integral numerator and
  denominator are already known;
- use ``FractionalIdeal::set_principal`` with an ambient ``Element`` when the
  generator may have denominators in the order basis;
- use ``FractionalIdeal::colon`` and ``FractionalIdeal::invert`` for
  inverse-like operations; ``invert`` is certified by multiplying back to the
  unit fractional ideal and returns ``false`` for noninvertible inputs while
  preserving the output.

Arithmetic methods such as ``add``, ``intersect``, ``multiply``,
``pow_fmpz``, ``colon``, and ``invert`` publish a candidate only after the
operation succeeds, so the tested APIs support output/input aliasing and
preserve the previous output value on documented failure paths.  Parent order
compatibility is semantic: parented objects keep their order handles alive,
and compatibility checks use the stored mathematical data rather than raw
pointer identity.

The native test suite covers degree-one arithmetic, quadratic and
nontrivial-order principal ideals, nongalois cubic arithmetic over
``theta^3 - 2``, nontrivial-order colon/inverse cases, and deterministic
seeded principal identity checks.  The remaining readiness work for this layer
is broader benchmark and consumer coverage rather than a known unmirrored core
arithmetic fixture.

Basic Walkthrough
-----------------

``examples/field_order_ideal_basics.cpp`` is the minimal compiled example for
this layer.  It constructs ``Q(sqrt(2))``, builds the equation order, forms the
principal integral ideal ``(theta)`` from an ``OrderElement``, and forms the
fractional principal ideal ``(theta/2)`` from an ambient ``Element``.

The example deliberately shows the generator-type distinction: the integral
ideal is generated by an element known to lie in the order, while the
fractional ideal is generated by an ambient field element whose denominator is
absorbed into the normalized numerator/denominator representation.

The same parent and failure rules apply in larger code:

- parented objects store value-like handles where they must keep mathematical
  parents alive;
- borrowed accessors such as ``parent()`` are compatibility views and should
  not be used as identity tokens;
- compare order identity with ``Order::has_same_data`` or
  ``same_order_parent``;
- check each boolean return value before consuming outputs;
- prefer domain operations such as ``set_principal``, ``contains``, ``norm``,
  ``add``, ``intersect``, ``multiply``, ``colon``, and ``invert`` instead of
  reaching into FLINT storage directly.

Benchmark Coverage
------------------

Benchmark coverage for this layer is split by workload:

- ``nf_ord`` maximal-order and p-maximal-overorder rows are mirrored by
  ``b-silex-nf_ord_maximal_order`` and
  ``b-silex-nf_ord_pmaximal_overorder``.
- factor-over-base ideal exponent rows are mirrored under the local algebra
  benchmark target ``b-silex-nf_idl``.

There is no standalone benchmark target today for general ideal or
fractional-ideal arithmetic beyond the factor-over-base surface.  Broader
ideal and fractional-ideal performance should be measured through local
targets or higher-level class/unit consumers until a standalone workload is
needed.

Source Lineage
--------------

The order layer follows PARI ``src/basemath/base2.c`` and Hecke's maximal-order
implementation for maximalization control, with the quadratic conductor and
discriminant identities stated in :doc:`algorithms_and_sources`.  Ideal
arithmetic follows PARI ``src/basemath/base4.c`` and Hecke
``Ideal/Ideal.jl`` where cited.  Current Silex behavior, exact failure
preservation, and native parent identity are covered by the order, element,
ideal, and fractional-ideal tests.

Examples
--------

See ``examples/maximal_order_zeta.cpp`` for maximal orders and zeta input
orders, ``examples/field_order_ideal_basics.cpp`` for basic integral and
fractional ideals, and :doc:`local_algebra` for the local and relation-layer
example.
