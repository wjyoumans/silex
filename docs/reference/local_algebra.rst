Prime Ideals, Factor Bases, Relations, and Residues
===================================================

This page documents the local algebra layer that sits between orders/ideals
and class/unit computations.  It covers prime ideals, factor bases, ideal
factorization over a factor base, relations, residue rings, and residue
fields.

The public model is C++-first: objects own copied mathematical data, store
value-like parent handles where long-lived parent relationships are needed,
and return ``false`` for unsupported cases instead of exposing raw buffers or
fabricating results.

Prime Ideals
------------

``silex::PrimeIdeal`` stores a copied integral ideal together with rational
prime metadata, ramification index, residue degree, and optional
residue-polynomial data.  ``PrimeIdealList`` is the owned output container
filled by ``decompose_prime``.  It exposes list-style access without requiring
callers to manage raw FLINT or C arrays.

``decompose_prime`` covers the currently supported prime-decomposition
contract:

* Dedekind-Kummer decomposition for supported equation orders;
* an explicit known-maximal monogenic path that can produce prime ideals
  without residue-polynomial data;
* repeated-factor failure behavior for uncertified nonmaximal orders;
* integral and fractional ideal valuations over certified maximal orders;
* nonzero field-element and compact factored-element valuations.

Check ``PrimeIdeal::has_prime_data`` and boolean return values before
consuming decomposition-dependent outputs.  ``get_ideal`` materializes the
copied integral ideal.  ``residue_polynomial`` and residue reduction require a
decomposition route that produced residue-polynomial data; the explicit
maximal monogenic path intentionally does not always provide it.
``valuation`` accepts ``OrderElement``, ``Element``, ``Ideal``,
``FractionalIdeal``, and ``FactoredElement`` inputs on the currently
certified/implemented paths. Fractional inputs subtract the exact rational
denominator valuation, and factored inputs accumulate signed factor valuations
without expanding the product. Zero elements, parent mismatches, unsupported
orders, and unrepresentable ``slong`` results return ``false`` without changing
the caller's output.

Factor Bases and Ideal Factorization
------------------------------------

``silex::FactorBase`` owns a deterministic list of copied prime ideals and
stores its ``Order`` parent handle by value.  The implemented builder uses
rational primes up to a bound.
``prime(out, i)`` copies one prime ideal out of the base, while
``prime_at(i)`` is the checked borrowed internal-style view for read-only
scans.  ``index`` and ``contains`` check compatible copied prime ideals.

Rational-prime block views expose the grouped layout needed by factorization
and class-group code without exposing raw owning arrays or making callers
depend on the internal storage.  The block API gives the rational prime, start
index, and length for each block.

``silex::IdealFactorization`` factors integral ideals over known maximal
orders, reconstructs them for verification, and produces dense factor-base
exponent rows for integral and fractional ideals.  The exponent rows are over
the copied factor-base ordering.  Failure paths, such as unsupported orders or
non-smooth ideals with respect to the supplied base, preserve the previous
output value.

Integral exponent rows are nonnegative valuation rows.  Fractional exponent
rows are signed: the numerator contributes positive exponents and the
denominator contributes negative exponents after both sides have been
certified over the same factor base.  The row is published only after the
numerator and denominator rows are both certified, so nonsmooth inputs,
wrong-sized rows, and parent mismatches leave the caller's row unchanged.

Relations
---------

``silex::Relation`` stores one certified principal relation and a copied
``FactorBase`` handle.  ``Relation::set_generator`` succeeds only when the
supplied generator factors smoothly over the copied factor base.  On success
it copies the generator and materializes the exponent row; on failure it
leaves the previous relation unchanged.

``silex::RelationMatrix`` owns rows, copied generators, and one copied
factor-base handle.  It appends same-base relations, materializes dense rows,
exposes copied or borrowed generators for supported access patterns, and
converts finite presentations to ``FiniteAbelianGroup``.
``RelationMatrix::append`` rejects relations from a different factor base
without mutating the matrix.

Relation search and collection live in ``ClassGroupContext`` rather than in
the standalone relation-storage API.  The relation layer is the storage and
presentation bridge used by class-group and order-unit code.

Residue Rings
-------------

``silex::ResidueRing`` owns its modulus ideal and works modulo that integral
ideal.  ``ResidueElement`` stores a ``ResidueRing`` parent handle.  The layer
supports canonical row reduction, element set/lift, arithmetic, and CRT
through ideal add-to-one witnesses.

``ResidueRing::cardinality()`` and ``ResidueElement::coordinates()`` are
owned-return convenience forms of the corresponding output-buffer exports.
This is the quotient by an arbitrary integral ideal; it does not promise field
operations unless the modulus has the required mathematical property.

Residue Fields
--------------

``silex::ResidueField`` owns its prime ideal data and works at prime ideals
with residue-polynomial data.  ``ResidueFieldElement`` and
``ResidueFieldQuotientLog`` store ``ResidueField`` parent handles.  The layer
exposes characteristic, degree, cardinality, modulus, canonical polynomial
representatives, order/ambient/factored element reduction, field arithmetic,
multiplicative orders, exact discrete logarithms, stateless quotient logs, and
cached quotient-log contexts.

The simple characteristic, cardinality, modulus, and polynomial
representative exports also have ``std::optional`` owned-return forms for
ordinary native C++ callers.  Reduction from order elements, ambient elements,
and factored elements is accepted when denominator and residue-data conditions
are satisfied.  Quotient-log contexts cache the ``ell``-local setup used by
repeated local proof computations.

Ambient and factored reductions are exact maps into the current residue field,
not local-repair routines.  An element whose denominator is not invertible at
the residue characteristic, a factored element with an unmappable factor, or a
negative power of a factor reducing to zero is rejected and preserves the
previous residue-field element.  Prime ideals produced through paths without
stored residue-polynomial data cannot be used to construct residue fields.

Lifetime and Failure Rules
--------------------------

Prime ideals, factor bases, relations, residue rings, and residue fields are
parented objects.  They keep mathematical parents alive with value-like
handles where needed.  Borrowed accessors such as ``parent()``,
``parent_order()``, ``factor_base()``, ``prime()``, and ``modulus()`` are
compatibility views for local access, not ownership or identity tokens.

Public methods return ``false`` for unsupported cases, invalid dimensions,
parent mismatches, noninvertible denominators, no-residue prime data,
non-smooth ideals, or insufficient certification.  Output objects are designed
to preserve their previous value on documented failure paths.

Benchmark Coverage
------------------

Current Google Benchmark coverage includes two local-layer surfaces:

* prime-ideal local reduction and valuation through
  ``bench/nf_prime_idl/b-local.cpp``;
* ideal factor-over-base exponent rows through
  ``bench/nf_idl/b-factor_over_base.cpp``.

There are no standalone benchmark targets today for ``FactorBase``,
``Relation``, ``ResidueRing``, or ``ResidueField``.  Performance for those
pieces should therefore be measured through local benchmarks or through
higher-level class/unit consumers until a real standalone workload is needed.

Local/Relation Walkthrough
--------------------------

``examples/local_relation_basics.cpp`` is the minimal compiled example for
this layer.  It works over the degree-one order ``Z`` so the mathematics is
easy to inspect while still exercising the same native C++ objects used by
higher-degree code:

* build a ``FactorBase`` up to a rational-prime bound;
* fetch a copied ``PrimeIdeal`` from the factor base and inspect its norm;
* set a ``Relation`` from a principal generator and read its exponent row;
* append the relation to a ``RelationMatrix``;
* construct a ``ResidueField`` at a prime ideal and multiply reduced order
  elements.

Higher-degree examples use the same object model, but decomposition and
residue-data availability depend on the order certification and supported
Dedekind-Kummer paths.

Source Lineage
--------------

Prime decomposition and residue behavior follow the cited PARI
``idealprimedec`` and Hecke prime-ideal/residue-field sources.  Factor-base,
factorization, and relation behavior also retain the established Silex
contracts and exact FLINT arithmetic.  See :doc:`algorithms_and_sources` for
the public file-level map.
