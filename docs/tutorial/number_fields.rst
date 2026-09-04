Number Fields
=============

This tutorial follows the native C++ path for field construction and exact
element arithmetic.  The corresponding compiled examples are
:download:`element_arithmetic_basics.cpp <../../examples/element_arithmetic_basics.cpp>`
and :download:`log_unit_lattice.cpp <../../examples/log_unit_lattice.cpp>`.
The focused fragments below emphasize API contracts; the complete compiled
element-arithmetic example is included at the end.

Construct Fields With Factories
-------------------------------

Use factory functions for ordinary public code:

.. code-block:: cpp

   namespace sflint = silex::flint;

   sflint::Fmpz d;
   sflint::fmpz_set_si(sflint::FmpzRef(d), 5);

   silex::NumberField K =
           silex::NumberField::quadratic(sflint::FmpzConstRef(d));
   if (!K.is_defined()) {
       return false;
   }

``NumberField`` is a value-like handle to shared field data.  Child objects
such as ``Element`` and ``EmbeddingContext`` store field handles where they
must keep the parent alive.  Borrowed raw FLINT accessors exist for bridge
code, but ordinary user code should not use raw ``nf_t`` storage.

Work With Elements Explicitly
-----------------------------

``Element`` operations are named methods with boolean status returns.  This
keeps expensive algebraic operations and failure behavior visible:

.. code-block:: cpp

   silex::Element theta(K);
   if (!theta.gen()) {
       return false;
   }

   silex::Element alpha(K);
   if (!alpha.add_si(theta, 3)) {
       return false;
   }

   sflint::Fmpq trace;
   sflint::Fmpq norm;
   if (!alpha.trace(sflint::FmpqRef(trace)) ||
       !alpha.norm(sflint::FmpqRef(norm))) {
       return false;
   }

The native API deliberately avoids operator-heavy exact algebra here.  Use
``add``, ``multiply``, ``invert``, ``pow_fmpz``, ``trace``, ``norm``,
``conjugate``, ``is_square``, and ``is_power`` so the call site shows what is
being computed and whether it can fail.

Embeddings And Unit Lattices
----------------------------

Use ``EmbeddingContext`` to refine certified complex roots before asking for
archimedean data:

.. code-block:: cpp

   silex::EmbeddingContext embeddings(K);
   if (!embeddings.refine(128)) {
       return false;
   }

   silex::flint::ArbVec logs(K.degree());
   if (!silex::logarithmic_embedding(
           logs, embeddings, theta,
           silex::LogEmbeddingMode::plain, 128)) {
       return false;
   }

Field-level unit helpers in ``silex/unit.hpp`` cover Dirichlet unit rank,
supported roots of unity, lower regulator bounds, supplied-unit log matrices,
supplied-unit subgroup regulators, supplied-unit independence checks, and
exact continued-fraction fundamental units for the explicit and canonical
polynomial-defined ``x^2-d`` real-quadratic routes.  See
:doc:`../support_matrix` for the presentations qualified for 0.1.0.

Next Steps
----------

- Read :doc:`../reference/number_fields` for the full field, element,
  embedding, unit, homomorphism, and automorphism API summary.
- Run ``example-element-arithmetic-basics`` for exact element arithmetic.
- Run ``example-log-unit-lattice`` for embeddings and logarithmic unit
  matrices.

Complete Compiled Example
-------------------------

.. literalinclude:: ../../examples/element_arithmetic_basics.cpp
   :language: cpp
   :caption: examples/element_arithmetic_basics.cpp
